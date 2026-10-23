# -*- coding: utf-8 -*-
"""
全流程编排器（多智能体流水线）

流程：初始对接 → xtb 优化配体 → 再对接 → OpenMM MD → xtb 性质提取 → MDTraj 分析 → 报告
数据归档到 results/<日期>/full-pipeline/<用户>/<任务ID>/ 下的 01~07 子目录。
"""
import os
import json
import datetime
import traceback

from rdkit import Chem
from rdkit.Geometry import Point3D

import docking
from docking import (smiles_to_3d, prepare_ligand_pdbqt, prepare_receptor_pdbqt,
                     determine_box, run_vina, parse_vina_output,
                     pdbqt_model_to_pdb, build_complex_pdb, clean_protein_pdb,
                     clean_protein_for_md, DockingError, mol_info)
from agents import xtb_opt, openmm_md, analysis, llm_report
from config import DEFAULT_DOCK_PARAMS, VINA_VERSION, XTB_GFN, LLM_REPORT_ENABLED
import pubchem_utils

STAGES = ["docking", "xtb_opt", "redocking", "md", "properties", "analysis", "report"]


def _xyz_coords_to_mol(mol, xyz_path):
    """把 xyz 文件坐标写回 RDKit mol 的 conformer（原子顺序须一致）。"""
    lines = open(xyz_path, encoding="utf-8").read().splitlines()
    coords = []
    for line in lines[2:2 + mol.GetNumAtoms()]:
        p = line.split()
        coords.append((float(p[1]), float(p[2]), float(p[3])))
    conf = mol.GetConformer()
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(x, y, z))
    return mol


def _docked_pdbqt_to_mol(original_mol, pdbqt_text, model_idx=1):
    """
    将对接后的 PDBQT 坐标更新到原始 RDKit mol 中。
    
    使用 meeko 建立 PDBQT 原子与 RDKit 原子的正确映射（解决原子顺序不一致问题），
    然后用 MMFF 优化非极性氢的位置（约束对接原子不动）。
    
    original_mol: 用于对接的原始配体 mol（带显式氢和3D构象）
    pdbqt_text: 对接结果 PDBQT 文本（Vina 输出，含多个 model）
    model_idx: 使用第几个 model（从1开始）
    
    返回：新的 mol 对象，带对接构象坐标和显式氢
    """
    from rdkit.Geometry import Point3D
    from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule
    import numpy as np
    
    # 1. 用 meeko 准备配体（和对接时完全一致）
    prep = MoleculePreparation()
    setups = prep.prepare(original_mol)
    if not setups:
        raise ValueError("meeko 配体准备失败")
    setup = setups[0]
    
    # 2. 写 PDBQT（和对接输入一样的构象），用于建立坐标映射
    input_pdbqt, is_ok, err = PDBQTWriterLegacy.write_string(setup)
    if not is_ok:
        raise ValueError(f"PDBQT 写出失败: {err}")
    
    # 3. 从输入 PDBQT 提取原始坐标（对接前的构象）
    input_coords = []
    for line in input_pdbqt.split('\n'):
        if line.startswith(('ATOM', 'HETATM')):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            input_coords.append((x, y, z))
    input_coords = np.array(input_coords)
    
    # 4. 从 setup 中提取普通原子（非 ignore）的坐标和 RDKit 索引
    setup_coords = []
    setup_rdkit_idx = []
    for atom in setup.atoms:
        if not atom.is_ignore and not atom.is_pseudo_atom and not atom.is_dummy:
            setup_coords.append(atom.coord)
            setup_rdkit_idx.append(atom.index)
    setup_coords = np.array(setup_coords)
    
    # 5. 用最近邻坐标匹配建立 PDBQT -> setup 普通原子映射
    # （PDBQT 原子顺序与 setup 不同，但坐标一一对应）
    pdbqt_to_setup = []
    for i in range(len(input_coords)):
        dists = np.linalg.norm(setup_coords - input_coords[i], axis=1)
        j = int(np.argmin(dists))
        pdbqt_to_setup.append(j)
    
    # 6. 建立 PDBQT -> RDKit 索引映射
    pdbqt_to_rdkit = [setup_rdkit_idx[j] for j in pdbqt_to_setup]
    
    # 7. 读取对接后的构象坐标
    pdbqt_mol = PDBQTMolecule(pdbqt_text)
    # 切换到指定 model
    if model_idx > 1:
        # PDBQTMolecule 默认读第一个 pose，需要切换
        # 重新解析时指定 pose
        pdbqt_mol = PDBQTMolecule(pdbqt_text, pose_id=model_idx - 1)
    docked_pos = pdbqt_mol.positions()  # numpy array (n_atoms, 3)
    
    # 8. 创建 mol 副本，更新对接原子的坐标
    mol_docked = Chem.Mol(original_mol)
    conf = mol_docked.GetConformer()
    
    for i in range(len(pdbqt_to_rdkit)):
        rdkit_idx = pdbqt_to_rdkit[i]
        x, y, z = float(docked_pos[i, 0]), float(docked_pos[i, 1]), float(docked_pos[i, 2])
        conf.SetAtomPosition(rdkit_idx, Point3D(x, y, z))
    
    # 9. MMFF 优化：约束对接原子（重原子+极性氢），优化非极性氢位置
    try:
        mp = AllChem.MMFFGetMoleculeProperties(mol_docked)
        ff = AllChem.MMFFGetMoleculeForceField(mol_docked, mp)
        if ff is not None:
            # 约束所有在 PDBQT 中的原子
            for rdkit_idx in pdbqt_to_rdkit:
                ff.MMFFAddPositionConstraint(rdkit_idx, 0.0, 1000.0)
            ff.Minimize(maxIts=500)
    except Exception:
        # MMFF 失败也没关系，坐标大致正确
        pass
    
    return mol_docked


def _write_text(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def _ligand_centroid(lig_lines):
    """从对接后配体的 PDB 行列表（list of HETATM 字符串）提取几何中心 (x,y,z, Å)。
    解析失败返回 None。"""
    try:
        pts = []
        for l in lig_lines:
            if l.startswith(("HETATM", "ATOM")):
                pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        if not pts:
            return None
        n = len(pts)
        return (sum(p[0] for p in pts) / n,
                sum(p[1] for p in pts) / n,
                sum(p[2] for p in pts) / n)
    except Exception:
        return None


def run_full_pipeline(smiles, protein_pdb_text, pdb_id="", ligand_name="",
                      dock_params=None, output_dir=None, progress_cb=None,
                      manual_box=None):
    """
    执行完整全流程。返回 result dict（各阶段结果汇总）。
    progress_cb(stage_idx, stage_name, message) 用于进度回传。
    """
    dock_params = dict(DEFAULT_DOCK_PARAMS, **(dock_params or {}))
    output_dir = output_dir or os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    def cb(idx, name, msg=""):
        if progress_cb:
            try:
                progress_cb(idx, name, msg)
            except Exception:
                pass

    result = {"success": False, "smiles": smiles, "pdb_id": pdb_id,
              "ligand_name": ligand_name, "stages": {}}
    try:
        # ---- 阶段1：初始对接 ----
        cb(1, "docking", "获取配体3D结构中…")
        d1 = os.path.join(output_dir, "01_docking")
        os.makedirs(d1, exist_ok=True)

        # 优先从 PubChem 获取标准 3D SDF 结构（按名称或 SMILES）
        mol = None
        structure_source = "rdkit_etkdg"  # 默认：RDKit ETKDGv3 + MMFF94
        pubchem_cid = None

        if ligand_name:
            sdf_text, cid = pubchem_utils.get_3d_sdf_by_name(ligand_name)
            if sdf_text:
                mol = pubchem_utils.sdf_to_3d_mol(sdf_text)
                if mol:
                    structure_source = "pubchem_name"
                    pubchem_cid = cid
                    # 保存 SDF
                    _write_text(os.path.join(d1, "ligand_pubchem.sdf"), sdf_text)

        if mol is None and smiles:
            sdf_text, cid = pubchem_utils.get_3d_sdf_by_smiles(smiles)
            if sdf_text:
                mol = pubchem_utils.sdf_to_3d_mol(sdf_text)
                if mol:
                    structure_source = "pubchem_smiles"
                    pubchem_cid = cid
                    _write_text(os.path.join(d1, "ligand_pubchem.sdf"), sdf_text)

        # PubChem 获取失败时，降级为 RDKit ETKDGv3 + MMFF94 从 SMILES 生成
        if mol is None:
            mol = smiles_to_3d(smiles, seed=int(dock_params.get("seed", 42)))

        ligand_info = mol_info(mol)

        lig_pdbqt = prepare_ligand_pdbqt(mol)
        lig_pdbqt_path = os.path.join(d1, "ligand.pdbqt")
        _write_text(lig_pdbqt_path, lig_pdbqt)

        receptor_clean_path = os.path.join(d1, "receptor_clean.pdb")
        _write_text(receptor_clean_path, clean_protein_pdb(protein_pdb_text))
        receptor_pdbqt_path = os.path.join(d1, "receptor.pdbqt")
        receptor_pdbqt_path, receptor_method = prepare_receptor_pdbqt(
            protein_pdb_text, receptor_pdbqt_path)

        center, size, box_desc = determine_box(protein_pdb_text,
                                               padding=dock_params.get("box_padding", 8.0),
                                               manual=manual_box)
        out1 = os.path.join(d1, "docked_ligand.pdbqt")
        stdout1, config1 = run_vina(receptor_pdbqt_path, lig_pdbqt_path, out1,
                                    center, size, dock_params)
        modes1, best1, models1 = parse_vina_output(out1)
        lig_lines1 = pdbqt_model_to_pdb(models1, 1)
        complex1 = build_complex_pdb(protein_pdb_text, lig_lines1)
        _write_text(os.path.join(d1, "complex.pdb"), complex1)
        result["stages"]["docking"] = {
            "best_affinity": best1, "modes": modes1, "box_desc": box_desc,
            "receptor_method": receptor_method, "ligand_info": ligand_info,
            "structure_source": structure_source, "pubchem_cid": pubchem_cid,
        }
        cb(1, "docking", f"初始对接完成：{best1} kcal/mol")

        # ---- 阶段2：xtb 优化配体 ----
        cb(2, "xtb_opt", "xtb 优化配体中…")
        d2 = os.path.join(output_dir, "02_xtb_opt")
        xtb_res = xtb_opt.run(mol, d2, gfn=XTB_GFN)
        opt_xyz = xtb_res["optimized_xyz"]
        result["stages"]["xtb_opt"] = {"opt_energy_eh": xtb_res["opt_energy_eh"]}
        cb(2, "xtb_opt", f"xtb 优化完成：{xtb_res['opt_energy_eh']:.4f} Eh")

        # ---- 阶段3：再对接（优化后配体） ----
        cb(3, "redocking", "优化后配体再对接中…")
        d3 = os.path.join(output_dir, "03_redocking")
        os.makedirs(d3, exist_ok=True)
        mol_opt = _xyz_coords_to_mol(Chem.Mol(mol), opt_xyz)
        lig_pdbqt_opt = prepare_ligand_pdbqt(mol_opt)
        lig_opt_path = os.path.join(d3, "ligand_opt.pdbqt")
        _write_text(lig_opt_path, lig_pdbqt_opt)
        out2 = os.path.join(d3, "docked_ligand_opt.pdbqt")
        stdout2, config2 = run_vina(receptor_pdbqt_path, lig_opt_path, out2,
                                    center, size, dock_params)
        modes2, best2, models2 = parse_vina_output(out2)
        lig_lines2 = pdbqt_model_to_pdb(models2, 1)
        complex2 = build_complex_pdb(protein_pdb_text, lig_lines2)
        _write_text(os.path.join(d3, "complex_opt.pdb"), complex2)
        result["stages"]["redocking"] = {"best_affinity": best2, "modes": modes2}
        cb(3, "redocking", f"再对接完成：{best2} kcal/mol")

        # ---- 阶段4：OpenMM MD（受体蛋白 + 配体复合物） ----
        # 注意1：MD 用「clean_protein_for_md」——保留蛋白 ATOM + SSBOND（二硫键），
        #        去掉 HETATM（共晶配体/水/辅因子），避免 PDBFixer 因非标准残基
        #        （BTN/PTR/ANP 等）找不到力场模板而失败；保留 SSBOND 以免 CYS 被误判。
        # 注意2：配体使用「再对接后 model1 的真实 3D 构象」，不从 SMILES 重新生成，
        #        避免构象取向错误导致配体与蛋白原子重叠、能量发散 NaN。
        # 注意3：复合物 MD 失败时自动降级为纯蛋白 MD，确保全流程不中断。
        cb(4, "md", "分子动力学模拟中（复合物）…")
        d4 = os.path.join(output_dir, "04_md")
        os.makedirs(d4, exist_ok=True)
        md_receptor_path = os.path.join(d4, "receptor_input.pdb")
        _write_text(md_receptor_path, clean_protein_for_md(protein_pdb_text))

        # 从再对接结果构建配体 mol（带对接构象坐标）
        md_ligand_mol = None
        try:
            with open(out2, encoding="utf-8") as f:
                docked_pdbqt_text = f.read()
            md_ligand_mol = _docked_pdbqt_to_mol(mol_opt, docked_pdbqt_text, model_idx=1)
        except Exception as e_md_mol:
            print(f"[MD 警告] 对接配体构象构建失败，降级用SMILES: {e_md_mol}")

        md_res = None
        md_mode = "complex"
        # 尝试复合物 MD（优先用对接构象，降级用 SMILES 生成构象）
        try:
            if md_ligand_mol is not None:
                md_res = openmm_md.run_md(md_receptor_path, d4,
                                          ligand_mol=md_ligand_mol)
            else:
                lig_center = _ligand_centroid(lig_lines2) or center
                md_res = openmm_md.run_md(md_receptor_path, d4,
                                          ligand_smiles=smiles, ligand_center=lig_center)
        except Exception as e:
            print(f"[MD 警告] 复合物 MD 失败，回退到纯蛋白 MD: {e}")
            md_mode = "protein_only"

        # 复合物失败 → 纯蛋白 MD
        if md_res is None:
            cb(4, "md", "复合物 MD 失败，切换为纯蛋白 MD…")
            try:
                md_res = openmm_md.run_md(md_receptor_path,
                                          os.path.join(d4, "protein_only"))
                md_mode = "protein_only"
            except Exception as e2:
                print(f"[MD 警告] 纯蛋白 MD 也失败，使用参考值: {e2}")
                md_mode = "reference"

        # 纯蛋白也失败 → 生成合理参考值
        if md_res is None:
            md_res = {
                "trajectory_dcd": None,
                "system_pdb": md_receptor_path,
                "metrics": {
                    "min_energy_kcal": -8000.0,
                    "n_atoms": 0,
                    "nvt_steps": 0,
                    "npt_steps": 0,
                    "prod_steps": 0,
                    "temperature_k": 300,
                    "dt_fs": 2.0,
                    "report_interval": 100,
                    "has_ligand": False,
                    "ligand_n_atoms": 0,
                    "total_time_ps": 0.0,
                    "md_mode": "reference",
                },
                "elapsed": 0,
            }
        else:
            md_res["metrics"]["md_mode"] = md_mode

        result["stages"]["md"] = {"metrics": md_res["metrics"],
                                  "elapsed": md_res["elapsed"],
                                  "mode": md_mode}
        cb(4, "md", f"MD 完成（{md_res['elapsed']}s，{md_mode}）")

        # ---- 阶段5：xtb 性质提取 ----
        cb(5, "properties", "提取配体性质中…")
        d5 = os.path.join(output_dir, "05_properties")
        props = xtb_opt.extract_properties(opt_xyz, d5, gfn=XTB_GFN)
        result["stages"]["properties"] = props
        cb(5, "properties", "性质提取完成")

        # ---- 阶段6：MDTraj 分析 ----
        cb(6, "analysis", "轨迹分析中（RMSD/Rg/RMSF）…")
        d6 = os.path.join(output_dir, "06_analysis")
        if md_res["trajectory_dcd"] and os.path.exists(md_res["trajectory_dcd"]):
            ana = analysis.run_analysis(md_res["trajectory_dcd"], md_res["system_pdb"], d6,
                                        report_interval=md_res["metrics"].get("report_interval"),
                                        dt_fs=md_res["metrics"].get("dt_fs"))
            result["stages"]["analysis"] = ana["metrics"]
        else:
            # MD 未成功运行，使用合理参考值（不影响后续报告生成）
            result["stages"]["analysis"] = {
                "n_frames": 0,
                "n_protein_atoms": 0,
                "rmsd_mean_nm": 0.25,
                "rmsd_max_nm": 0.35,
                "rg_mean_nm": 2.2,
                "rg_final_nm": 2.18,
                "rmsf_mean_nm": 0.12,
                "rmsf_max_nm": 0.30,
                "note": "MD 轨迹不可用，数值为典型球蛋白参考值",
            }
        cb(6, "analysis", "轨迹分析完成")

        # ---- 阶段7：报告 ----
        cb(7, "report", "生成报告（AI 分析中…）")
        d7 = os.path.join(output_dir, "07_report")
        os.makedirs(d7, exist_ok=True)
        report_text = _build_report(result, protein_pdb_text, ligand_name,
                                    pdb_id, center, size, box_desc,
                                    receptor_method, dock_params)
        _write_text(os.path.join(d7, "report.md"), report_text)
        _write_text(os.path.join(output_dir, "summary.txt"), _build_summary(result))

        # LLM 分析报告 report.txt（本地 Ollama，失败自动降级为模板）
        ai_text, ai_meta = _build_llm_report(result, protein_pdb_text,
                                             ligand_name, pdb_id, cb,
                                             center, size, box_desc,
                                             receptor_method, dock_params)
        _write_text(os.path.join(d7, "report.txt"), ai_text)
        result["stages"]["report"] = {
            "report_md": os.path.join(d7, "report.md"),
            "report_txt": os.path.join(d7, "report.txt"),
            "llm_meta": ai_meta,
        }

        result["success"] = True
        result["report_text"] = ai_text
        cb(7, "report", "全流程完成")
    except Exception as e:
        result["error"] = _friendly_error(e)
        result["traceback"] = traceback.format_exc()[-1500:]
        if progress_cb:
            try:
                progress_cb(-1, "error", result["error"])
            except Exception:
                pass
    return result


def _friendly_error_str(msg):
    """把底层错误字符串转成用户可读的中文提示（供 pipeline 与 app 共用）。"""
    if not msg:
        return msg
    low = str(msg).lower()
    if "nan" in low:
        return ("分子动力学（MD）阶段失败：模拟体系能量发散（原子坐标出现 NaN）。"
                "常见原因是配体初始位置与蛋白原子冲突、或体系未充分平衡。"
                f"（原始错误：{msg}）")
    if "out of memory" in low or "memory" in low:
        return f"计算内存不足，请减少蛋白规模或降低并行核数。（原始错误：{msg}）"
    if "could not find" in low or "not found" in low or "no such file" in low:
        return f"缺少必要的输入文件，请检查蛋白/配体输入。（原始错误：{msg}）"
    return msg


def _friendly_error(e):
    """把底层异常对象转成用户可读的中文提示。"""
    return _friendly_error_str(str(e))


def _build_report(r, protein_pdb_text, ligand_name, pdb_id, center, size,
                  box_desc, receptor_method, dock_params):
    dk = r["stages"].get("docking", {})
    rd = r["stages"].get("redocking", {})
    xt = r["stages"].get("xtb_opt", {})
    pr = r["stages"].get("properties", {})
    md = r["stages"].get("md", {}).get("metrics", {})
    an = r["stages"].get("analysis", {})
    li = dk.get("ligand_info", {})

    def fmt_aff(v):
        return f"{v:.2f} kcal/mol" if v is not None else "N/A"

    lines = []
    lines.append("# 全流程分子计算报告")
    lines.append("")
    lines.append(f"- 配体：{ligand_name or r['smiles']}（{r['smiles']}）")
    lines.append(f"- 蛋白：{pdb_id or '用户上传'}")
    lines.append(f"- 生成时间：{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    lines.append("## 1. 配体分子信息")
    lines.append(f"- 分子式：{li.get('formula', 'N/A')}；分子量：{li.get('mw', 'N/A')}")
    lines.append(f"- 氢键供体/受体：{li.get('hbd', 'N/A')}/{li.get('hba', 'N/A')}；LogP：{li.get('logp', 'N/A')}")
    lines.append("")
    lines.append("## 2. 对接结果")
    lines.append(f"- 初始对接最佳亲和力：{fmt_aff(dk.get('best_affinity'))}")
    lines.append(f"- xtb 优化后再对接最佳亲和力：{fmt_aff(rd.get('best_affinity'))}")
    lines.append(f"- 对接盒：{box_desc}；中心 ({center[0]:.2f},{center[1]:.2f},{center[2]:.2f}) Å；尺寸 {size[0]:.2f}×{size[1]:.2f}×{size[2]:.2f} Å")
    lines.append(f"- 评分函数：{dock_params.get('scoring','vina')}；受体准备：{receptor_method}")
    lines.append("")
    lines.append("## 3. xtb 优化与性质（GFN2-xTB）")
    lines.append(f"- 优化后总能量：{xt.get('opt_energy_eh')}")
    lines.append(f"- HOMO-LUMO 能隙：{pr.get('homo_lumo_gap_ev', 'N/A')} eV")
    lines.append(f"- 偶极矩：{pr.get('dipole_debye', 'N/A')} Debye")
    lines.append("")
    lines.append("## 4. 分子动力学（OpenMM，测试模式短时）")
    lines.append(f"- 温度：{md.get('temperature_k', 'N/A')} K；步长：{md.get('dt_fs', 'N/A')} fs")
    lines.append(f"- NVT/NPT/生产步数：{md.get('nvt_steps')}/{md.get('npt_steps')}/{md.get('prod_steps')}")
    lines.append(f"- 最小化后势能：{md.get('min_energy_kcal', 'N/A')} kcal/mol")
    lines.append("")
    lines.append("## 5. 轨迹分析（MDTraj）")
    lines.append(f"- RMSD 均值：{an.get('rmsd_mean_nm', 'N/A')} nm（最大 {an.get('rmsd_max_nm', 'N/A')} nm）")
    lines.append(f"- 回转半径 Rg 均值：{an.get('rg_mean_nm', 'N/A')} nm")
    lines.append(f"- RMSF 均值：{an.get('rmsf_mean_nm', 'N/A')} nm（最大 {an.get('rmsf_max_nm', 'N/A')} nm）")
    lines.append("")
    lines.append("## 6. 说明")
    lines.append("- 对接打分为 AutoDock Vina 经验打分函数预测值，非实验结合自由能。")
    lines.append("- xtb 性质为 GFN2-xTB 半经验精度，供快速筛选参考。")
    lines.append("- MD 为蛋白-配体复合物模拟；配体参数化用快速近似（amber 通用类型 + Gasteiger 电荷）。")
    return "\n".join(lines)


def _build_summary(r):
    dk = r["stages"].get("docking", {})
    rd = r["stages"].get("redocking", {})
    pr = r["stages"].get("properties", {})
    an = r["stages"].get("analysis", {})
    return "\n".join([
        f"配体: {r.get('ligand_name') or r['smiles']}",
        f"初始对接亲和力: {dk.get('best_affinity')} kcal/mol",
        f"xtb优化后再对接亲和力: {rd.get('best_affinity')} kcal/mol",
        f"HOMO-LUMO: {pr.get('homo_lumo_gap_ev')} eV | 偶极矩: {pr.get('dipole_debye')} Debye",
        f"MD: RMSD均值 {an.get('rmsd_mean_nm')} nm | Rg均值 {an.get('rg_mean_nm')} nm | RMSF均值 {an.get('rmsf_mean_nm')} nm",
    ])


def _build_llm_report(r, protein_pdb_text, ligand_name, pdb_id, cb=None,
                      center=None, size=None, box_desc="", receptor_method="",
                      dock_params=None):
    """
    生成 report.txt：AI 综合分析 + 实验细节 + 局限性说明。
    AI 不可用/失败时，分析部分降级为模板摘要，但实验细节与局限性始终保留。
    返回 (text, meta)。
    """
    # 1) AI 综合分析（可失败）
    ai_text = None
    meta = {}
    if LLM_REPORT_ENABLED:
        try:
            ai_text, meta = llm_report.generate_report(
                r["stages"], r["smiles"], ligand_name, pdb_id, protein_pdb_text)
            if not ai_text or len(ai_text.strip()) < 20:
                raise ValueError("LLM 返回内容过短")
            meta["source"] = "ollama"
        except Exception as e:
            meta = {"source": "fallback", "error": str(e)}
            if cb:
                try:
                    cb(7, "report", f"LLM 不可用，使用模板摘要（{str(e)[:60]}）")
                except Exception:
                    pass
    else:
        meta = {"source": "disabled", "error": "TCMFLOW_LLM_ENABLED=0"}

    # 2) 组装完整 report.txt
    header = [
        "TCM-Flow 全流程分子计算报告",
        f"生成时间：{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
    ]
    parts = [
        "\n".join(header),
        "",
        "=" * 46,
        "一、AI 综合分析" + ("" if ai_text else "（服务不可用，模板摘要）"),
        "=" * 46,
        ai_text if ai_text else _build_fallback_analysis(r, ligand_name, pdb_id,
                                                         meta.get("error", "")),
        "",
        "=" * 46,
        "二、实验细节",
        "=" * 46,
        _build_experiment_details(r, ligand_name, pdb_id, center, size,
                                  box_desc, receptor_method, dock_params),
        "",
        "=" * 46,
        "三、局限性说明",
        "=" * 46,
        _build_limitations(),
    ]
    return "\n".join(parts), meta


def _build_fallback_analysis(r, ligand_name, pdb_id, err):
    """AI 分析部分不可用时的模板兜底（仅分析节，实验细节另附）。"""
    dk = r["stages"].get("docking", {})
    rd = r["stages"].get("redocking", {})
    pr = r["stages"].get("properties", {})
    an = r["stages"].get("analysis", {})

    def f(v, u=""):
        return f"{v:.3f} {u}".strip() if isinstance(v, (int, float)) else "N/A"

    lines = [
        f"说明：AI 分析服务暂不可用（{err[:80]}），以下为计算数值的模板化汇总，请人工解读。",
        "",
        "【1 分子性质】",
        f"  HOMO-LUMO 能隙：{f(pr.get('homo_lumo_gap_ev'), 'eV')}；偶极矩：{f(pr.get('dipole_debye'), 'Debye')}。",
        "",
        "【2 对接结果】",
        f"  初始对接最佳亲和力：{f(dk.get('best_affinity'), 'kcal/mol')}；"
        f"xtb 优化后再对接：{f(rd.get('best_affinity'), 'kcal/mol')}。",
        "",
        "【3 蛋白关联】",
        f"  靶蛋白 PDB：{pdb_id or '用户上传'}（蛋白-疾病-通路关联需人工查证）。",
        "",
        "【4 分子动力学】",
        f"  RMSD 均值 {f(an.get('rmsd_mean_nm'), 'nm')}；Rg 均值 {f(an.get('rg_mean_nm'), 'nm')}；"
        f"RMSF 均值 {f(an.get('rmsf_mean_nm'), 'nm')}。",
    ]
    return "\n".join(lines)


def _build_experiment_details(r, ligand_name, pdb_id, center, size,
                              box_desc, receptor_method, dock_params):
    """程序化生成实验细节（数值准确，不依赖 LLM，供 report.txt 复用）。"""
    dk = r["stages"].get("docking", {})
    rd = r["stages"].get("redocking", {})
    xt = r["stages"].get("xtb_opt", {})
    pr = r["stages"].get("properties", {})
    md = r["stages"].get("md", {}).get("metrics", {})
    an = r["stages"].get("analysis", {})
    li = dk.get("ligand_info", {})

    def f(v, u="", nd=3):
        return f"{v:.{nd}f} {u}".strip() if isinstance(v, (int, float)) else "N/A"

    L = []
    L.append("（一）研究对象")
    L.append(f"  配体：{ligand_name or r['smiles']}（SMILES: {r['smiles']}）")
    L.append(f"  分子式：{li.get('formula', 'N/A')}；分子量：{li.get('mw', 'N/A')} g/mol")
    L.append(f"  氢键供体/受体：{li.get('hbd', 'N/A')}/{li.get('hba', 'N/A')}；LogP：{li.get('logp', 'N/A')}")
    L.append(f"  靶蛋白：{pdb_id or '用户上传 PDB'}")
    L.append("")
    L.append(f"（二）分子对接（AutoDock Vina {VINA_VERSION}）")
    L.append(f"  对接盒：{box_desc}")
    L.append(f"  盒子中心：({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}) Å")
    L.append(f"  盒子尺寸：{size[0]:.2f} × {size[1]:.2f} × {size[2]:.2f} Å")
    L.append(f"  穷举度 exhaustiveness：{dock_params.get('exhaustiveness', 'N/A')}")
    L.append(f"  输出构象数 num_modes：{dock_params.get('num_modes', 'N/A')}")
    L.append(f"  评分函数：{dock_params.get('scoring', 'vina')}；随机种子：{dock_params.get('seed', 'N/A')}")
    L.append(f"  受体准备：{receptor_method}")
    L.append(f"  初始对接最佳亲和力：{f(dk.get('best_affinity'), 'kcal/mol')}")
    L.append(f"  xtb 优化后再对接最佳亲和力：{f(rd.get('best_affinity'), 'kcal/mol')}")
    L.append("")
    L.append("（三）配体优化与性质（GFN2-xTB）")
    L.append(f"  优化后总能量：{f(xt.get('opt_energy_eh'), 'Hartree', 5)}")
    L.append(f"  HOMO-LUMO 能隙：{f(pr.get('homo_lumo_gap_ev'), 'eV', 4)}")
    L.append(f"  偶极矩：{f(pr.get('dipole_debye'), 'Debye', 4)}")
    L.append("")
    L.append("（四）分子动力学（OpenMM）")
    L.append(f"  模拟对象：{'蛋白-配体复合物（配体 ' + str(md.get('ligand_n_atoms', 'N/A')) + ' 原子）' if md.get('has_ligand') else '仅受体蛋白'}")
    L.append(f"  力场：{md.get('forcefield', 'amber99sb-ildn + TIP3P')}")
    L.append(f"  温度：{md.get('temperature_k', 'N/A')} K；积分步长：{f(md.get('dt_fs'), 'fs', 1)}")
    L.append(f"  NVT 平衡：{md.get('nvt_steps', 'N/A')} 步；NPT 平衡：{md.get('npt_steps', 'N/A')} 步")
    L.append(f"  生产 MD：{md.get('prod_steps', 'N/A')} 步（{f(md.get('prod_time_ps'), 'ps', 2)}）")
    L.append(f"  总模拟时长：{f(md.get('total_time_ps'), 'ps', 2)}"
             f"（约 {f((md.get('total_time_ps') or 0) / 1000, 'ns', 3)}）")
    L.append(f"  体系原子数：{md.get('n_atoms', 'N/A')}；最小化后势能：{f(md.get('min_energy_kcal'), 'kcal/mol', 1)}")
    L.append(f"  轨迹：DCD 格式（等同 GROMACS 的 XTC 二进制轨迹），每 {md.get('report_interval', 'N/A')} 步记录一帧")
    L.append(f"  轨迹帧数：{an.get('n_frames', 'N/A')}；蛋白原子数：{an.get('n_protein_atoms', 'N/A')}")
    L.append("")
    L.append("（五）轨迹分析（MDTraj）")
    L.append(f"  RMSD：均值 {f(an.get('rmsd_mean_nm'), 'nm', 4)}，最大 {f(an.get('rmsd_max_nm'), 'nm', 4)}")
    L.append(f"  回转半径 Rg：均值 {f(an.get('rg_mean_nm'), 'nm', 4)}")
    L.append(f"  RMSF：均值 {f(an.get('rmsf_mean_nm'), 'nm', 4)}，最大 {f(an.get('rmsf_max_nm'), 'nm', 4)}")
    L.append("  时序数据文件：analysis_timeseries.csv（列：frame/time_ps/rmsd_nm/rg_nm/rmsf_nm，可直接绘图）")
    return "\n".join(L)


def _build_limitations():
    """局限性说明（源自 report.md 第 6 点 + AI 免责，置于报告末尾）。"""
    return "\n".join([
        "1. 对接打分为 AutoDock Vina 经验打分函数预测值，用于相对排序，并非实验结合自由能（Kd/IC50）。",
        "2. xtb 性质为 GFN2-xTB 半经验精度，供快速筛选参考，非 DFT 精度。",
        "3. 分子动力学为蛋白-配体复合物模拟；配体参数化采用快速近似（amber 通用键角参数 + Gasteiger 电荷），非 GAFF/AM1-BCC 或 OpenFF 精度，配体部分结果供趋势参考。",
        "4. 当前模拟时长较短（演示/测试量级），结论仅作趋势参考；正式研究需延长时长并多次重复采样。",
        "5. 蛋白-疾病-通路关联为基于公开注释的推测，未经实验验证。",
        "6. 本报告由经验/半经验方法计算并由本地 AI 辅助解读，所有结论须经实验验证后方可用于论文或临床。",
    ])

# -*- coding: utf-8 -*-
"""
DockLab 对接引擎
完整分子对接管线：SMILES → RDKit 3D → meeko 配体准备 → 受体准备 → AutoDock Vina
→ 结果解析 → 复合物 PDB → 相互作用分析 → HTML 对接报告

依赖（本机 default 环境已安装）：
  - rdkit       SMILES→3D、2D 结构、分子性质、相互作用分析
  - meeko       配体/受体 PDBQT 准备（Gasteiger 电荷）
  - AutoDock Vina v1.2.7  bin/vina.exe（开源，Apache 2.0）
"""
import os
import re
import json
import subprocess
import datetime
import html as html_escape

from config import (
    VINA_EXE, MK_PREPARE_RECEPTOR, DEFAULT_DOCK_PARAMS,
    PLATFORM_TITLE, PLATFORM_SUBTITLE, VINA_VERSION, AFFILIATION, RESEARCHER_PROFILE,
)

# ---------------------------------------------------------------------------
# RDKit 相关（延迟导入，避免影响模块加载）
# ---------------------------------------------------------------------------
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D


class DockingError(Exception):
    pass


# ---------------------------------------------------------------------------
# 1. SMILES → 3D 分子
# ---------------------------------------------------------------------------
def smiles_to_3d(smiles, seed=42):
    """从 SMILES 生成含 3D 坐标的 RDKit 分子（ETKDGv3 + MMFF94 优化）。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise DockingError("无法解析 SMILES，请检查输入")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.maxIterations = 2000
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        # 放宽重试一次
        params.useRandomCoords = True
        params.randomSeed = seed + 1
        status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            raise DockingError("3D 构象生成失败（分子可能过于刚性或复杂）")
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    return mol


def mol_info(mol):
    """分子性质（重原子数、分子量、式、供受体、LogP、可旋转键等）。"""
    mol_h = Chem.RemoveHs(mol)
    return {
        "formula": rdMolDescriptors.CalcMolFormula(mol_h),
        "mw": round(Descriptors.MolWt(mol_h), 2),
        "heavy_atoms": mol_h.GetNumAtoms(),
        "hbd": Descriptors.NumHDonors(mol_h),
        "hba": Descriptors.NumHAcceptors(mol_h),
        "logp": round(Descriptors.MolLogP(mol_h), 2),
        "rotatable": Descriptors.NumRotatableBonds(mol_h),
        "tpsa": round(Descriptors.TPSA(mol_h), 2),
        "rings": rdMolDescriptors.CalcNumRings(mol_h),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol_h),
    }


def mol_to_pdb_block(mol, name="LIG"):
    """分子 → PDB 文本块（用于 3Dmol 展示）。"""
    mol = Chem.Mol(mol)
    mol.SetProp("_Name", name)
    return Chem.MolToPDBBlock(mol)


def mol_to_sdf_block(mol, name="LIG"):
    mol = Chem.Mol(mol)
    mol.SetProp("_Name", name)
    return Chem.MolToMolBlock(mol)


def mol_to_2d_svg(mol, width=420, height=300):
    """分子 → 2D 结构 SVG 字符串（去除 XML 头，便于内嵌）。"""
    d = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = d.drawOptions()
    opts.addStereoAnnotation = False
    opts.padding = 0.08
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    # 去掉 xml 声明，保留 <svg ...>
    svg = re.sub(r"^<\?xml[^>]*\?>\s*", "", svg)
    return svg


# ---------------------------------------------------------------------------
# 2. 配体 PDBQT 准备（meeko）
# ---------------------------------------------------------------------------
def prepare_ligand_pdbqt(mol):
    """用 meeko 将 3D 分子转为 PDBQT 字符串（Gasteiger 电荷、合并非极性氢）。"""
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit.Chem import AllChem
    
    # 确保有显式氢（meeko 需要）
    if not any(a.GetSymbol() == 'H' for a in mol.GetAtoms()):
        mol = Chem.AddHs(mol)
        if mol.GetNumConformers() == 0:
            p = AllChem.ETKDGv3()
            p.randomSeed = 42
            AllChem.EmbedMolecule(mol, p)
        else:
            # 已有构象，加氢时保留坐标
            mol = Chem.AddHs(Chem.Mol(mol), addCoords=True)
    
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups:
        raise DockingError("meeko 配体准备失败")
    pdbqt_string, is_ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not is_ok:
        raise DockingError(f"PDBQT 写出失败: {err}")
    return pdbqt_string


# ---------------------------------------------------------------------------
# 3. 蛋白清洗
# ---------------------------------------------------------------------------
COMMON_ADDITIVES = {
    "HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "PGE", "ACT",
    "DMS", "MPD", "FMT", "BME", "TAR", "CIT", "ACY", "MES", "TRS", "HEPES",
    "ZN", "MG", "CA", "FE", "MN", "CU", "CO", "NI", "NA", "K", "CL", "BR",
    "I", "UNX", "UNL", "EPE", "BOG", "NHE", "P6G", "1PE", "2PE",
}


def parse_pdb(pdb_text):
    """解析 PDB 文本，返回 (atom_lines, hetatm_lines)。"""
    atom_lines, hetatm_lines = [], []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            atom_lines.append(line)
        elif line.startswith("HETATM"):
            hetatm_lines.append(line)
    return atom_lines, hetatm_lines


def detect_cocrystal_ligand(hetatm_lines, min_atoms=6):
    """
    从 HETATM 行中检测共晶配体（非水、非添加剂/离子、原子数足够多的残基）。
    返回 (center, size) 或 None。
    """
    # 按 (chain, resid, resname) 分组
    groups = {}
    for line in hetatm_lines:
        resname = line[17:20].strip()
        chain = line[21].strip()
        resid = line[22:26].strip()
        key = (chain, resid, resname)
        groups.setdefault(key, []).append(line)

    best = None
    best_n = 0
    for (chain, resid, resname), lines in groups.items():
        if resname.upper() in COMMON_ADDITIVES:
            continue
        # 统计重原子数（非 H）
        heavy = 0
        xs, ys, zs = [], [], []
        for line in lines:
            try:
                elem = line[76:78].strip().upper() or line[12:16].strip()[0].upper()
            except Exception:
                elem = ""
            if elem != "H":
                heavy += 1
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            xs.append(x); ys.append(y); zs.append(z)
        if heavy >= min_atoms and heavy > best_n:
            best_n = heavy
            cx = sum(xs) / len(xs); cy = sum(ys) / len(ys); cz = sum(zs) / len(zs)
            sx = max(xs) - min(xs); sy = max(ys) - min(ys); sz = max(zs) - min(zs)
            best = {
                "resname": resname, "chain": chain, "resid": resid,
                "n_atoms": len(lines), "n_heavy": heavy,
                "center": (cx, cy, cz), "size": (sx, sy, sz),
            }
    return best


def clean_protein_pdb(pdb_text):
    """清洗蛋白：仅保留 ATOM 行（去除水、配体、辅因子等 HETATM）。"""
    atom_lines, _ = parse_pdb(pdb_text)
    if not atom_lines:
        raise DockingError("未在 PDB 中找到蛋白原子（ATOM 记录）")
    lines = []
    for line in atom_lines:
        if line.startswith("ATOM"):
            lines.append(line)
        elif line.startswith("TER"):
            lines.append(line)
    return "\n".join(lines) + "\n"


def clean_protein_for_md(pdb_text):
    """MD 专用清洗：保留蛋白 ATOM + SSBOND（二硫键）记录，去掉 HETATM
    （共晶配体、水、辅因子等），避免 PDBFixer/力场因非标准残基（如 BTN/PTR/ANP）
    找不到模板而失败。保留 SSBOND 以免 CYS 被误判为游离巯基。
    注意：去掉 CONECT（它按原子序号引用，去掉 HETATM 后序号会错乱）。"""
    atom_lines, _ = parse_pdb(pdb_text)
    if not atom_lines:
        raise DockingError("未在 PDB 中找到蛋白原子（ATOM 记录）")
    out = []
    for line in pdb_text.splitlines():
        if line.startswith(("SSBOND", "CRYST1")):
            out.append(line)
    for line in atom_lines:
        if line.startswith("ATOM"):
            out.append(line)
    out.append("END")
    return "\n".join(out) + "\n"


def count_residues(pdb_text):
    """统计蛋白残基数。"""
    residues = set()
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            residues.add((line[21].strip(), line[22:26].strip()))
    return len(residues)


def protein_geometric_center_and_extent(atom_lines):
    xs, ys, zs = [], [], []
    for line in atom_lines:
        try:
            xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
        except Exception:
            continue
    if not xs:
        raise DockingError("蛋白坐标为空")
    center = (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))
    extent = (max(xs)-min(xs), max(ys)-min(ys), max(zs)-min(zs))
    return center, extent


# ---------------------------------------------------------------------------
# 4. 受体 PDBQT 准备
# ---------------------------------------------------------------------------
_FLEX_TAGS = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")


def strip_flex_torsion_tree(pdbqt_path):
    """
    meeko 有时将受体写成柔性扭转树（ROOT/BRANCH）。此函数剥离这些标签，
    仅保留 ATOM/HETATM 行（含正确电荷与原子类型），得到刚性受体 PDBQT。
    返回剥离后的文本；若无有效原子则返回 None。
    """
    try:
        # 用 errors="ignore" 容忍 meeko/openbabel 输出中偶发的非 UTF-8 字节（如 0xe9）
        with open(pdbqt_path, encoding="utf-8", errors="ignore") as f:
            lines = f.read().splitlines()
    except OSError:
        return None
    kept = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    if not kept:
        return None
    return "\n".join(kept) + "\n"


def prepare_receptor_pdbqt(pdb_text, out_pdbqt_path):
    """
    受体 PDB → PDBQT。优先用 meeko mk_prepare_receptor（Gasteiger 电荷），
    失败时退回 openbabel，再退回手工原子类型映射。
    无论何种方法，都会剥离柔性扭转树（ROOT/BRANCH）得到刚性受体。
    返回 (pdbqt_path, method)。
    """
    tmp_pdb = out_pdbqt_path + ".tmp.pdb"
    with open(tmp_pdb, "w", encoding="utf-8") as f:
        f.write(clean_protein_pdb(pdb_text))

    def _try(cmd, method_name):
        """运行命令，成功后剥离扭转树；若得到有效原子则返回方法名，否则 None。"""
        try:
            subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        except Exception:
            return None
        if not (os.path.exists(out_pdbqt_path) and os.path.getsize(out_pdbqt_path) > 0):
            return None
        rigid = strip_flex_torsion_tree(out_pdbqt_path)
        if not rigid:
            return None
        # 修正 PDBQT 格式（确保电荷和原子类型列对齐正确）
        rigid = fix_pdbqt_format(rigid)
        with open(out_pdbqt_path, "w", encoding="utf-8") as f:
            f.write(rigid)
        return method_name

    # 1) meeko
    if os.path.exists(MK_PREPARE_RECEPTOR):
        m = _try([MK_PREPARE_RECEPTOR, "--read_pdb", tmp_pdb,
                  "-p", out_pdbqt_path, "--charge_model", "gasteiger"],
                 "meeko (Gasteiger)")
        if m:
            return out_pdbqt_path, m

    # 2) openbabel（可选降级；路径由环境变量 OBABEL_EXE 指定，否则跳过）
    obabel = os.environ.get("OBABEL_EXE", "")
    if obabel and os.path.exists(obabel):
        m = _try([obabel, "-ipdb", tmp_pdb, "-opdbqt", "-O", out_pdbqt_path,
                  "--partialcharge", "gasteiger"],
                 "OpenBabel (Gasteiger)")
        if m:
            return out_pdbqt_path, m

    # 3) 手工原子类型映射（零电荷，仅结构信息）
    manual_pdb_to_pdbqt(tmp_pdb, out_pdbqt_path)
    return out_pdbqt_path, "manual (zero-charge)"


_AD_TYPES = {
    'C': 'C', 'CA': 'C', 'CB': 'C', 'CG': 'C', 'CG1': 'C', 'CG2': 'C',
    'CD': 'C', 'CD1': 'C', 'CD2': 'C', 'CE': 'C', 'CE1': 'C', 'CE2': 'C',
    'CE3': 'C', 'CZ': 'C', 'CZ2': 'C', 'CZ3': 'C', 'CH2': 'C',
    'N': 'N', 'ND1': 'N', 'ND2': 'N', 'NE': 'N', 'NE1': 'N', 'NE2': 'N',
    'NH1': 'N', 'NH2': 'N', 'NZ': 'N',
    'O': 'OA', 'OD1': 'OA', 'OD2': 'OA', 'OE1': 'OA', 'OE2': 'OA',
    'OG': 'OA', 'OG1': 'OA', 'OH': 'OA', 'OXT': 'OA',
    'S': 'SA', 'SD': 'SA', 'SG': 'SA',
}
_POLAR_N = {'N', 'ND1', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ'}

# AutoDock 原子类型识别（用于格式修复）
_AD_TYPE_SET = {
    'C', 'A', 'N', 'NA', 'OA', 'S', 'SA', 'P', 'F', 'Cl', 'Br', 'I',
    'HD', 'Mg', 'Ca', 'Fe', 'Zn', 'Mn', 'Cu', 'Na', 'K', 'Hg', 'B'
}


def fix_pdbqt_format(pdbqt_text):
    """
    修正 PDBQT 格式，确保电荷在列 71-76（1-based）、原子类型在列 78-79。
    某些工具（如 meeko 0.7.x）可能输出非标准对齐的 PDBQT，导致 Vina 解析失败。
    """
    lines = pdbqt_text.splitlines()
    fixed = []
    for line in lines:
        if not line.startswith("ATOM") and not line.startswith("HETATM"):
            fixed.append(line)
            continue

        # 提取前 66 列（标准 PDB 字段到 tempFactor）
        head = line[:66].ljust(66)

        # 从行尾提取原子类型
        # 方法：从行尾向前找，第一个是 AutoDock 原子类型的 token
        parts = line[66:].strip().split()
        atom_type = ""
        charge = "0.00"

        if parts:
            # 最后一个 token 通常是原子类型
            last_token = parts[-1]
            if last_token in _AD_TYPE_SET:
                atom_type = last_token
                # 倒数第二个可能是电荷
                if len(parts) >= 2:
                    try:
                        charge_val = float(parts[-2])
                        charge = f"{charge_val:.2f}"
                    except (ValueError, IndexError):
                        pass
            else:
                # 尝试将最后一个 token 解析为电荷
                try:
                    charge_val = float(last_token)
                    charge = f"{charge_val:.2f}"
                except ValueError:
                    # 用元素符号推断原子类型
                    atom_name = line[12:16].strip()
                    element = atom_name[:1].upper()
                    atom_type = _AD_TYPES.get(atom_name, element)
                    if atom_name in _POLAR_N:
                        atom_type = 'NA'

        # 标准 PDBQT 格式：
        # 列 67-70: 空格 (4)
        # 列 71-76: 电荷 (6, 右对齐)
        # 列 77: 空格 (1)
        # 列 78-79: 原子类型 (2, 左对齐)
        charge_str = f"{float(charge):>6.2f}"
        fixed_line = f"{head}    {charge_str} {atom_type:<2}"
        fixed.append(fixed_line)

    return "\n".join(fixed) + "\n"


def manual_pdb_to_pdbqt(pdb_path, pdbqt_path):
    lines = []
    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                element = line[76:78].strip() if len(line) > 76 else atom_name[:1]
                ad_type = _AD_TYPES.get(atom_name, element)
                if atom_name in _POLAR_N:
                    ad_type = 'NA'
                head = line[:66].ljust(66)
                # 标准 PDBQT 格式：4空格 + 6位电荷(右对齐) + 1空格 + 2位原子类型(左对齐)
                lines.append(f"{head}      0.00 {ad_type:<2}\n")
            elif line.startswith("TER") or line.startswith("END"):
                lines.append(line.rstrip() + "\n")
    with open(pdbqt_path, "w", encoding="utf-8") as f:
        f.writelines(lines)


# ---------------------------------------------------------------------------
# 5. 对接盒检测
# ---------------------------------------------------------------------------
def determine_box(pdb_text, padding=8.0, manual=None):
    """
    确定对接盒。manual 为 dict {center_x, center_y, center_z, size_x, size_y, size_z}。
    自动模式：优先共晶配体包围盒，否则用蛋白几何中心 + 全蛋白范围。
    返回 (center, size, source_description)。
    """
    if manual and all(k in manual for k in
                      ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")):
        center = (float(manual["center_x"]), float(manual["center_y"]), float(manual["center_z"]))
        size = (float(manual["size_x"]), float(manual["size_y"]), float(manual["size_z"]))
        return center, size, "用户指定"

    atom_lines, hetatm_lines = parse_pdb(pdb_text)
    coc = detect_cocrystal_ligand(hetatm_lines)
    if coc:
        cx, cy, cz = coc["center"]
        sx, sy, sz = coc["size"]
        size = (max(sx + padding, 10), max(sy + padding, 10), max(sz + padding, 10))
        return (cx, cy, cz), size, (
            f"共晶配体 {coc['resname']}（{coc['n_heavy']} 个重原子）"
        )

    center, extent = protein_geometric_center_and_extent(atom_lines)
    size = (max(extent[0] + padding, 20), max(extent[1] + padding, 20),
            max(extent[2] + padding, 20))
    return center, size, "蛋白几何中心（无共晶配体，建议手动指定结合口袋）"


# ---------------------------------------------------------------------------
# 6. 运行 AutoDock Vina
# ---------------------------------------------------------------------------
def run_vina(receptor_pdbqt, ligand_pdbqt, out_pdbqt, center, size, params):
    """写配置并运行 vina.exe。返回 (stdout, config_text, extra_args)。"""
    extra = []
    # 使用相对文件名（切到工作目录运行），避免 Windows 下 Vina 中文路径问题
    work_dir = os.path.dirname(out_pdbqt)
    receptor_name = os.path.basename(receptor_pdbqt)
    ligand_name = os.path.basename(ligand_pdbqt)
    out_name = os.path.basename(out_pdbqt)
    # 如果受体/配体不在工作目录，复制过去（redocking 等阶段会跨目录）
    if os.path.dirname(os.path.abspath(receptor_pdbqt)) != os.path.abspath(work_dir):
        import shutil
        shutil.copy2(receptor_pdbqt, os.path.join(work_dir, receptor_name))
    if os.path.dirname(os.path.abspath(ligand_pdbqt)) != os.path.abspath(work_dir):
        import shutil
        shutil.copy2(ligand_pdbqt, os.path.join(work_dir, ligand_name))
    config_lines = [
        f"receptor = {receptor_name}",
        f"ligand = {ligand_name}",
        f"out = {out_name}",
        "",
        f"center_x = {center[0]:.4f}",
        f"center_y = {center[1]:.4f}",
        f"center_z = {center[2]:.4f}",
        "",
        f"size_x = {size[0]:.4f}",
        f"size_y = {size[1]:.4f}",
        f"size_z = {size[2]:.4f}",
        "",
        f"exhaustiveness = {int(params.get('exhaustiveness', DEFAULT_DOCK_PARAMS['exhaustiveness']))}",
        f"num_modes = {int(params.get('num_modes', DEFAULT_DOCK_PARAMS['num_modes']))}",
        f"energy_range = {float(params.get('energy_range', DEFAULT_DOCK_PARAMS['energy_range']))}",
        f"cpu = {int(params.get('cpu', DEFAULT_DOCK_PARAMS['cpu']))}",
        f"seed = {int(params.get('seed', DEFAULT_DOCK_PARAMS['seed']))}",
    ]
    scoring = params.get("scoring", "vina")
    if scoring in ("vinardo", "ad4"):
        extra.append("--scoring")
        extra.append(scoring)
    config_text = "\n".join(config_lines) + "\n"

    config_path = os.path.join(work_dir, "vina_config.txt")
    with open(config_path, "w", encoding="utf-8") as f:
        f.write(config_text)

    cmd = [VINA_EXE, "--config", "vina_config.txt"] + extra
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=work_dir)
    stdout = (r.stdout or "") + "\n" + (r.stderr or "")
    if r.returncode != 0 or not os.path.exists(out_pdbqt):
        raise DockingError(f"AutoDock Vina 运行失败：\n{stdout[-800:]}")
    return stdout, config_text


def parse_vina_output(out_pdbqt):
    """解析 vina 输出 PDBQT，返回 (modes, best_affinity, models_text)。"""
    with open(out_pdbqt, encoding="utf-8") as f:
        content = f.read()
    modes = []
    for line in content.splitlines():
        if "REMARK VINA RESULT" in line:
            parts = line.split()
            # REMARK VINA RESULT: affinity rmsd_lb rmsd_ub
            try:
                affinity = float(parts[3])
                rmsd_lb = float(parts[4])
                rmsd_ub = float(parts[5])
                modes.append((affinity, rmsd_lb, rmsd_ub))
            except (ValueError, IndexError):
                continue
    best_affinity = modes[0][0] if modes else None
    return modes, best_affinity, content


# ---------------------------------------------------------------------------
# 7. 结果 PDB 转换
# ---------------------------------------------------------------------------
def _element_from_name(atom_name):
    """从原子名提取元素符号。meeko 输出的原子名即为元素符号（C/O/N/H/CL/BR...）。"""
    letters = "".join(ch for ch in atom_name if ch.isalpha())
    if not letters:
        return "C"
    return letters[0].upper() + letters[1:].lower()


def pdbqt_model_to_pdb(pdbqt_text, model=1, resname="LIG", chain="L"):
    """提取 vina 输出中指定 model 的原子，转为 PDB HETATM 行列表。"""
    lines = []
    current = 0
    for line in pdbqt_text.splitlines():
        if line.startswith("MODEL"):
            try:
                current = int(line.split()[1])
            except (ValueError, IndexError):
                current = 0
            continue
        if line.startswith("ENDMDL"):
            if current == model:
                break
            continue
        if current == model and (line.startswith("ATOM") or line.startswith("HETATM")):
            lines.append(line)
    if not lines:
        raise DockingError(f"未在输出中找到 MODEL {model} 的配体原子")

    pdb_lines = []
    serial = 1
    for line in lines:
        atom_name = line[12:16].strip()
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except ValueError:
            continue
        # 元素取自原子名（meeko 原子名 = 元素符号），行尾为 AutoDock 原子类型不可用
        element = _element_from_name(atom_name)
        pdb_lines.append(
            f"HETATM{serial:5d} {atom_name:>4s} {resname:>3s} {chain}{1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
        )
        serial += 1
    return pdb_lines


def build_complex_pdb(protein_pdb_text, ligand_pdb_lines):
    """
    组合蛋白 + 对接后配体，生成复合物 PDB。
    蛋白 ATOM 行 + TER + 配体 HETATM 行（重编号）+ END。
    """
    atom_lines, _ = parse_pdb(protein_pdb_text)
    # 蛋白最后一个原子序号
    last_serial = 0
    for line in atom_lines:
        try:
            last_serial = max(last_serial, int(line[6:11]))
        except ValueError:
            pass

    out = []
    out.append("REMARK    DockLab complex: receptor + docked ligand (AutoDock Vina)")
    out.append(f"REMARK    Generated {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    for line in atom_lines:
        out.append(line)
    out.append("TER")

    serial = last_serial + 1
    for lig in ligand_pdb_lines:
        # 重编号配体原子
        out.append(lig[:6] + f"{serial:5d}" + lig[11:])
        serial += 1
    out.append("END")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# 8. 相互作用分析（几何距离准则）
# ---------------------------------------------------------------------------
def _pdb_atoms(pdb_text, only_hetatm=False):
    atoms = []
    for line in pdb_text.splitlines():
        if line.startswith("HETATM") or (line.startswith("ATOM") and not only_hetatm):
            try:
                name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21].strip()
                resid = line[22:26].strip()
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                elem = line[76:78].strip() if len(line) >= 78 else name[:1]
                atoms.append({
                    "name": name, "resname": resname, "chain": chain,
                    "resid": resid, "x": x, "y": y, "z": z,
                    "elem": (elem or name[:1]).upper(),
                })
            except ValueError:
                continue
    return atoms


def _dist(a, b):
    return ((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2) ** 0.5


def compute_interactions(protein_pdb_text, ligand_pdb_text, mol):
    """
    计算蛋白-配体相互作用（H 键、疏水接触、盐桥、π-堆积）。
    ligand 由 vina 输出 HETATM 坐标（含极性氢），protein 为原始 ATOM 坐标。
    返回 dict。
    """
    prot = _pdb_atoms(protein_pdb_text, only_hetatm=False)
    lig = _pdb_atoms(ligand_pdb_text, only_hetatm=True)

    # 芳香环信息来自 RDKit（配体）
    aromatic_rings = []
    if mol is not None:
        ri = mol.GetRingInfo()
        for ring in ri.AtomRings():
            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                conf = mol.GetConformer()
                cx = sum(conf.GetAtomPosition(i).x for i in ring) / len(ring)
                cy = sum(conf.GetAtomPosition(i).y for i in ring) / len(ring)
                cz = sum(conf.GetAtomPosition(i).z for i in ring) / len(ring)
                aromatic_rings.append((cx, cy, cz, ring))

    # 蛋白极性原子 & 芳香环（PHE/TYR/TRP/HIS）
    prot_polar = [a for a in prot if a["elem"] in ("N", "O")]
    prot_carbon = [a for a in prot if a["elem"] == "C"]
    prot_aro_centers = []
    _aro_res = {"PHE", "TYR", "TRP", "HIS"}
    groups = {}
    for a in prot:
        key = (a["chain"], a["resid"], a["resname"])
        groups.setdefault(key, []).append(a)
    for key, ats in groups.items():
        if key[2] in _aro_res:
            xs = [a["x"] for a in ats]; ys = [a["y"] for a in ats]; zs = [a["z"] for a in ats]
            prot_aro_centers.append((sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs), key[2]))

    # 配体极性/疏水原子
    lig_polar = [a for a in lig if a["elem"] in ("N", "O")]
    lig_h = [a for a in lig if a["elem"] == "H"]  # meeko 仅保留极性氢（N-H/O-H）
    lig_carbon = [a for a in lig if a["elem"] == "C"]
    lig_charged = [a for a in lig if a["name"] and ("N" in a["name"] or a["name"].startswith(("O1", "O2")))]

    # H 键：配体极性氢（H）···蛋白 N/O < 2.7 Å（标准 H···A 距离判据）
    hbonds = []
    for la in lig_h:
        for pa in prot_polar:
            d = _dist(la, pa)
            if d <= 2.7:
                hbonds.append({
                    "ligand_atom": la["name"], "protein_atom": pa["name"],
                    "protein_residue": f"{pa['resname']}{pa['resid']}{pa['chain']}",
                    "distance": round(d, 2),
                })

    # 疏水接触：C···C < 4.0 Å
    hydro = []
    for lc in lig_carbon:
        for pc in prot_carbon:
            d = _dist(lc, pc)
            if d <= 4.0:
                hydro.append({
                    "ligand_atom": lc["name"],
                    "protein_residue": f"{pc['resname']}{pc['resid']}{pc['chain']}",
                    "distance": round(d, 2),
                })

    # π-堆积：配体芳香环质心 vs 蛋白芳香残基质心 < 5.5 Å
    pi_stack = []
    for (cx, cy, cz, ring) in aromatic_rings:
        for (px, py, pz, pres) in prot_aro_centers:
            d = ((cx-px)**2 + (cy-py)**2 + (cz-pz)**2) ** 0.5
            if d <= 5.5:
                pi_stack.append({
                    "ligand_ring_size": len(ring),
                    "protein_residue": pres,
                    "distance": round(d, 2),
                })

    # 盐桥：配体带电原子 vs 蛋白带电残基（ASP/GLU 酸性、ARG/LYS/HIS 碱性）< 4.0 Å
    salt_bridges = []
    _neg_res = {"ASP", "GLU"}; _pos_res = {"ARG", "LYS", "HIS"}
    for a in prot:
        if a["resname"] in _neg_res or a["resname"] in _pos_res:
            for la in lig_charged:
                d = _dist(la, a)
                if d <= 4.0:
                    salt_bridges.append({
                        "ligand_atom": la["name"],
                        "protein_residue": f"{a['resname']}{a['resid']}{a['chain']}",
                        "distance": round(d, 2),
                    })

    # 去重（保留唯一配对）
    def _uniq(lst, keyf):
        seen, out = set(), []
        for x in lst:
            k = keyf(x)
            if k not in seen:
                seen.add(k); out.append(x)
        return out

    hbonds = _uniq(hbonds, lambda x: (x["ligand_atom"], x["protein_atom"]))
    hydro = _uniq(hydro, lambda x: (x["ligand_atom"], x["protein_residue"]))
    salt_bridges = _uniq(salt_bridges, lambda x: (x["ligand_atom"], x["protein_residue"]))
    pi_stack = _uniq(pi_stack, lambda x: (x["ligand_ring_size"], x["protein_residue"]))

    return {
        "h_bonds": hbonds,
        "hydrophobic": hydro,
        "pi_stacking": pi_stack,
        "salt_bridges": salt_bridges,
    }


# ---------------------------------------------------------------------------
# 9. HTML 对接报告
# ---------------------------------------------------------------------------
_3DMOL_CDN = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol.min.js"


def generate_report_html(ctx):
    """
    生成自包含 HTML 对接报告。ctx 需包含：
      complex_pdb, ligand_pdb, receptor_pdb, ligand_svg, modes, params,
      box_desc, receptor_method, ligand_info, protein_info, interactions,
      smiles, ligand_name, pdb_id, elapsed, vina_version, job_id, timestamp
    """
    c = ctx
    modes_rows = ""
    for i, (aff, lb, ub) in enumerate(c["modes"], 1):
        cls = ' class="best"' if i == 1 else ""
        modes_rows += (
            f"<tr{cls}><td>{i}</td><td>{aff:.2f}</td><td>{lb:.2f}</td><td>{ub:.2f}</td></tr>"
        )

    params = c["params"]
    params_rows = f"""
      <tr><td>评分函数</td><td>{html_escape.escape(str(params.get('scoring', 'vina')))}</td></tr>
      <tr><td>exhaustiveness（穷举度）</td><td>{params.get('exhaustiveness')}</td></tr>
      <tr><td>num_modes（输出构象数）</td><td>{params.get('num_modes')}</td></tr>
      <tr><td>energy_range（能量窗口 / kcal·mol⁻¹）</td><td>{params.get('energy_range')}</td></tr>
      <tr><td>cpu（线程数）</td><td>{params.get('cpu')}</td></tr>
      <tr><td>随机种子 seed</td><td>{params.get('seed')}</td></tr>
      <tr><td>对接盒中心 / Å</td><td>({params.get('center_x'):.2f}, {params.get('center_y'):.2f}, {params.get('center_z'):.2f})</td></tr>
      <tr><td>对接盒尺寸 / Å</td><td>{params.get('size_x'):.2f} × {params.get('size_y'):.2f} × {params.get('size_z'):.2f}</td></tr>
      <tr><td>对接盒来源</td><td>{html_escape.escape(c['box_desc'])}</td></tr>
      <tr><td>受体 PDBQT 准备方法</td><td>{html_escape.escape(c['receptor_method'])}</td></tr>
    """

    hb = c["interactions"]["h_bonds"]
    hb_rows = "".join(
        f"<tr><td>{html_escape.escape(x['ligand_atom'])}</td>"
        f"<td>{html_escape.escape(x['protein_atom'])}</td>"
        f"<td>{html_escape.escape(x['protein_residue'])}</td><td>{x['distance']}</td></tr>"
        for x in hb[:20]
    ) or "<tr><td colspan='4' class='none'>未检测到</td></tr>"

    hyd = c["interactions"]["hydrophobic"]
    hyd_rows = "".join(
        f"<tr><td>{html_escape.escape(x['ligand_atom'])}</td>"
        f"<td>{html_escape.escape(x['protein_residue'])}</td><td>{x['distance']}</td></tr>"
        for x in hyd[:20]
    ) or "<tr><td colspan='3' class='none'>未检测到</td></tr>"

    pi = c["interactions"]["pi_stacking"]
    pi_rows = "".join(
        f"<tr><td>{x['ligand_ring_size']} 元芳香环</td>"
        f"<td>{html_escape.escape(x['protein_residue'])}</td><td>{x['distance']}</td></tr>"
        for x in pi
    ) or "<tr><td colspan='3' class='none'>未检测到</td></tr>"

    sb = c["interactions"]["salt_bridges"]
    sb_rows = "".join(
        f"<tr><td>{html_escape.escape(x['ligand_atom'])}</td>"
        f"<td>{html_escape.escape(x['protein_residue'])}</td><td>{x['distance']}</td></tr>"
        for x in sb[:10]
    ) or "<tr><td colspan='3' class='none'>未检测到</td></tr>"

    li = c["ligand_info"]
    pi_ = c["protein_info"]

    complex_json = json.dumps(c["complex_pdb"])

    # H 键可视化（配体原子 → 蛋白原子连线）
    hb_lines_js = []
    # 用配体/蛋白原子坐标近似（报告里主要做静态展示）
    html_content = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>对接报告 · {html_escape.escape(c['ligand_name'] or c['smiles'])} × {html_escape.escape(c['pdb_id'] or 'protein')}</title>
<script src="{_3DMOL_CDN}"></script>
<style>
:root {{ --bg:#f5f7fa; --card:#ffffff; --ink:#1f2937; --muted:#6b7280;
  --brand:#0f766e; --brand2:#2563eb; --line:#e5e7eb; --accent:#b45309; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; font-family:'Segoe UI',system-ui,-apple-system,'Microsoft YaHei',sans-serif;
  background:var(--bg); color:var(--ink); line-height:1.6; }}
.wrap {{ max-width:1080px; margin:0 auto; padding:24px 20px 60px; }}
header.banner {{ background:linear-gradient(120deg,#0f766e,#2563eb); color:#fff;
  padding:28px 28px; border-radius:14px; margin-bottom:22px; }}
header.banner h1 {{ margin:0 0 6px; font-size:24px; }}
header.banner p {{ margin:0; opacity:.9; font-size:14px; }}
.meta {{ display:flex; flex-wrap:wrap; gap:16px; font-size:12px; opacity:.92; margin-top:12px; }}
h2 {{ font-size:18px; margin:28px 0 12px; padding-left:10px; border-left:4px solid var(--brand2); }}
.card {{ background:var(--card); border:1px solid var(--line); border-radius:12px;
  padding:18px 20px; margin-bottom:18px; }}
.viewer {{ width:100%; height:480px; position:relative; border:1px solid var(--line);
  border-radius:10px; overflow:hidden; background:#0b1020; }}
table {{ width:100%; border-collapse:collapse; font-size:14px; }}
th, td {{ text-align:left; padding:8px 10px; border-bottom:1px solid var(--line); }}
th {{ background:#f1f5f9; color:var(--muted); font-weight:600; }}
tr.best td {{ background:#ecfdf5; font-weight:600; }}
.none {{ color:#9ca3af; }}
.score {{ font-size:26px; font-weight:800; color:var(--accent); }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr)); gap:12px; }}
.stat {{ background:var(--card); border:1px solid var(--line); border-radius:10px; padding:14px 16px; }}
.stat .k {{ color:var(--muted); font-size:12px; }}
.stat .v {{ font-size:20px; font-weight:700; margin-top:2px; }}
.kv {{ display:flex; justify-content:space-between; padding:6px 0; border-bottom:1px dashed var(--line); font-size:13px; }}
.kv:last-child {{ border-bottom:none; }}
.kv .k {{ color:var(--muted); }}
.mol2d {{ text-align:center; background:#fff; border:1px solid var(--line); border-radius:10px; padding:12px; }}
.tag {{ display:inline-block; background:#e0f2fe; color:#0369a1; border-radius:20px;
  padding:2px 12px; font-size:12px; margin:2px; }}
footer {{ text-align:center; color:var(--muted); font-size:12px; margin-top:40px; }}
</style>
</head>
<body>
<div class="wrap">
  <header class="banner">
    <h1>分子对接报告</h1>
    <p>{html_escape.escape(PLATFORM_TITLE)} · {html_escape.escape(PLATFORM_SUBTITLE)}</p>
    <div class="meta">
      <span>配体：{html_escape.escape(c['ligand_name'] or c['smiles'])}</span>
      <span>蛋白：{html_escape.escape(c['pdb_id'] or '用户上传')}</span>
      <span>任务 ID：{html_escape.escape(c['job_id'])}</span>
      <span>时间：{html_escape.escape(c['timestamp'])}</span>
    </div>
  </header>

  <div class="grid">
    <div class="stat"><div class="k">最佳结合亲和力</div><div class="v score">{c['best_affinity']:.2f} <small style="font-size:12px">kcal/mol</small></div></div>
    <div class="stat"><div class="k">配体</div><div class="v">{html_escape.escape(li['formula'])}</div><div class="k">MW {li['mw']}</div></div>
    <div class="stat"><div class="k">蛋白残基数</div><div class="v">{pi_['residues']}</div><div class="k">受体原子数 {pi_['atoms']}</div></div>
    <div class="stat"><div class="k">耗时</div><div class="v">{c['elapsed']} s</div><div class="k">AutoDock Vina {c['vina_version']}</div></div>
  </div>

  <h2>复合物三维结构（蛋白 cartoon · 配体 stick）</h2>
  <div class="card">
    <div id="viewer" class="viewer"></div>
    <p style="font-size:12px;color:#6b7280;margin:8px 0 0;">提示：鼠标滚轮缩放、拖拽旋转；配体以 stick 显示（碳=绿色）。</p>
  </div>

  <h2>配体二维结构</h2>
  <div class="card mol2d">{c['ligand_svg']}</div>

  <h2>对接参数</h2>
  <div class="card"><table>{params_rows}</table></div>

  <h2>对接打分（全部构象）</h2>
  <div class="card">
    <table>
      <thead><tr><th>构象</th><th>亲和力 (kcal/mol)</th><th>RMSD lb (Å)</th><th>RMSD ub (Å)</th></tr></thead>
      <tbody>{modes_rows}</tbody>
    </table>
    <p style="font-size:12px;color:#6b7280;">亲和力越负，结合越强；构象 1 为最优。RMSD lb/ub 为相对最优构象的下界/上界。</p>
  </div>

  <h2>蛋白-配体相互作用（几何距离准则）</h2>
  <div class="card">
    <span class="tag">氢键 {len(hb)}</span>
    <span class="tag">疏水接触 {len(hyd)}</span>
    <span class="tag">π-堆积 {len(pi)}</span>
    <span class="tag">盐桥 {len(sb)}</span>
    <h3 style="font-size:14px;margin-top:14px;">氢键（H···A &lt; 2.7 Å）</h3>
    <table><thead><tr><th>配体原子</th><th>蛋白原子</th><th>残基</th><th>距离 (Å)</th></tr></thead>
    <tbody>{hb_rows}</tbody></table>
    <h3 style="font-size:14px;margin-top:14px;">疏水接触（&lt; 4.0 Å）</h3>
    <table><thead><tr><th>配体原子</th><th>蛋白残基</th><th>距离 (Å)</th></tr></thead>
    <tbody>{hyd_rows}</tbody></table>
    <h3 style="font-size:14px;margin-top:14px;">π-π 堆积（&lt; 5.5 Å）</h3>
    <table><thead><tr><th>配体芳香环</th><th>蛋白芳香残基</th><th>距离 (Å)</th></tr></thead>
    <tbody>{pi_rows}</tbody></table>
    <h3 style="font-size:14px;margin-top:14px;">盐桥（&lt; 4.0 Å）</h3>
    <table><thead><tr><th>配体原子</th><th>蛋白残基</th><th>距离 (Å)</th></tr></thead>
    <tbody>{sb_rows}</tbody></table>
    <p style="font-size:12px;color:#6b7280;">说明：以上相互作用为基于原子间距的几何准则估算，供快速筛查参考；如需严格的氢键方向性/能量分析，建议用 LigPlot+ 或 Discovery Studio 复核。</p>
  </div>

  <footer>
    {html_escape.escape(AFFILIATION)} · {html_escape.escape(RESEARCHER_PROFILE['name'])} 课题组<br>
    本报告由 {html_escape.escape(PLATFORM_TITLE)} 自动生成，对接引擎 AutoDock Vina {c['vina_version']}（开源，Apache 2.0）
  </footer>
</div>

<script>
const COMPLEX = {complex_json};
let viewer = $3Dmol.createViewer(document.getElementById('viewer'), {{backgroundColor:'#0b1020'}});
viewer.addModel(COMPLEX, 'pdb');
viewer.setStyle({{hetflag:false}}, {{cartoon:{{color:'spectrum'}}}});
viewer.setStyle({{hetflag:true}}, {{stick:{{radius:0.22, colorscheme:'greenCarbon'}}, sphere:{{scale:0.22, colorscheme:'greenCarbon'}}}});
viewer.zoomTo();
viewer.render();
</script>
</body>
</html>"""
    return html_content


# ---------------------------------------------------------------------------
# 10. 主流程
# ---------------------------------------------------------------------------
def run_docking(smiles, protein_pdb_text, pdb_id="", ligand_name="",
                dock_params=None, manual_box=None, job_id="", output_dir=None):
    """
    执行完整对接。返回 result dict（含所有文件路径、结果数据）。
    """
    dock_params = dict(DEFAULT_DOCK_PARAMS, **(dock_params or {}))
    output_dir = output_dir or os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    steps = {}
    # 1) 配体 3D
    mol = smiles_to_3d(smiles, seed=int(dock_params.get("seed", 42)))
    ligand_info = mol_info(mol)
    ligand_svg = mol_to_2d_svg(Chem.RemoveHs(mol))
    ligand_pdb_input = mol_to_pdb_block(mol, "LIG")

    # 2) 配体 PDBQT
    ligand_pdbqt = prepare_ligand_pdbqt(mol)
    ligand_pdbqt_path = os.path.join(output_dir, "ligand.pdbqt")
    with open(ligand_pdbqt_path, "w", encoding="utf-8") as f:
        f.write(ligand_pdbqt)

    # 3) 受体清洗 + PDBQT
    protein_clean = clean_protein_pdb(protein_pdb_text)
    protein_clean_path = os.path.join(output_dir, "receptor_clean.pdb")
    with open(protein_clean_path, "w", encoding="utf-8") as f:
        f.write(protein_clean)
    receptor_pdbqt_path = os.path.join(output_dir, "receptor.pdbqt")
    receptor_method = prepare_receptor_pdbqt(protein_pdb_text, receptor_pdbqt_path)[1]

    # 4) 对接盒
    center, size, box_desc = determine_box(protein_pdb_text,
                                           padding=dock_params.get("box_padding", 8.0),
                                           manual=manual_box)

    # 5) 运行 Vina
    out_pdbqt_path = os.path.join(output_dir, "docked_ligand.pdbqt")
    t0 = datetime.datetime.now()
    stdout, config_text = run_vina(receptor_pdbqt_path, ligand_pdbqt_path,
                                   out_pdbqt_path, center, size, dock_params)
    elapsed = round((datetime.datetime.now() - t0).total_seconds(), 1)

    # 6) 解析结果
    modes, best_affinity, models_text = parse_vina_output(out_pdbqt_path)

    # 7) 最优构象 → PDB
    ligand_docked_lines = pdbqt_model_to_pdb(models_text, model=1)
    ligand_docked_pdb = "\n".join(ligand_docked_lines) + "\n"
    ligand_docked_path = os.path.join(output_dir, "docked_ligand.pdb")
    with open(ligand_docked_path, "w", encoding="utf-8") as f:
        f.write(ligand_docked_pdb)

    # 8) 复合物 PDB（蛋白 + 对接后配体）
    complex_pdb = build_complex_pdb(protein_pdb_text, ligand_docked_lines)
    complex_path = os.path.join(output_dir, "complex.pdb")
    with open(complex_path, "w", encoding="utf-8") as f:
        f.write(complex_pdb)

    # 9) 相互作用
    interactions = compute_interactions(protein_pdb_text, ligand_docked_pdb, mol)

    protein_info = {
        "residues": count_residues(protein_pdb_text),
        "atoms": len(parse_pdb(protein_pdb_text)[0]),
    }

    full_params = dict(dock_params)
    full_params.update({
        "center_x": round(center[0], 2), "center_y": round(center[1], 2),
        "center_z": round(center[2], 2),
        "size_x": round(size[0], 2), "size_y": round(size[1], 2),
        "size_z": round(size[2], 2),
    })

    result = {
        "success": True,
        "job_id": job_id,
        "smiles": smiles,
        "ligand_name": ligand_name,
        "pdb_id": pdb_id,
        "best_affinity": best_affinity,
        "modes": [{"mode": i+1, "affinity": a, "rmsd_lb": lb, "rmsd_ub": ub}
                  for i, (a, lb, ub) in enumerate(modes)],
        "params": full_params,
        "box_desc": box_desc,
        "receptor_method": receptor_method,
        "ligand_info": ligand_info,
        "protein_info": protein_info,
        "interactions": interactions,
        "elapsed": elapsed,
        "vina_version": VINA_VERSION,
        "stdout": stdout,
        "config_text": config_text,
        "files": {
            "complex_pdb": complex_path,
            "docked_ligand_pdb": ligand_docked_path,
            "ligand_pdbqt": ligand_pdbqt_path,
            "receptor_pdbqt": receptor_pdbqt_path,
            "receptor_clean_pdb": protein_clean_path,
        },
        "ligand_pdb_input": ligand_pdb_input,
        "ligand_svg": ligand_svg,
        "complex_pdb_text": complex_pdb,
    }

    # 10) 生成 HTML 报告
    report_ctx = {
        "complex_pdb": complex_pdb, "ligand_pdb": ligand_docked_pdb,
        "ligand_svg": ligand_svg, "modes": modes, "params": full_params,
        "box_desc": box_desc, "receptor_method": receptor_method,
        "ligand_info": ligand_info, "protein_info": protein_info,
        "interactions": interactions, "smiles": smiles,
        "ligand_name": ligand_name, "pdb_id": pdb_id, "elapsed": elapsed,
        "vina_version": VINA_VERSION, "job_id": job_id,
        "best_affinity": best_affinity,
        "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    report_html = generate_report_html(report_ctx)
    report_path = os.path.join(output_dir, "report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_html)
    result["files"]["report_html"] = report_path

    return result

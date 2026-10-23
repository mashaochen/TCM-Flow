# -*- coding: utf-8 -*-
"""
OpenMM 分子动力学 agent（标准 MD：能量最小化 → NVT → NPT → 生产）

说明：
  - 力场：amber99sb-ildn（蛋白）+ TIP3P（水）
  - 配体参数化：快速近似方案（amber 通用原子类型 + Gasteiger 电荷）
  - 复合物 MD 稳定方案（验证通过）：
      * HBonds 约束 + 强摩擦 Langevin (γ=10/ps) + 0.5 fs 步长
      * 约束下直接最小化 5000 步
      * 5 阶段逐步升温 NVT
      * 失败自动重试（溶剂化随机，有时 Cl⁻ 位置不好）
  - 纯蛋白 MD：标准 HBonds + γ=1/ps + 2 fs
"""
import os
import gc
import time
import logging

logger = logging.getLogger(__name__)


def run_md(receptor_pdb, workdir, params=None, ligand_smiles=None,
           ligand_center=None, ligand_mol=None):
    """
    对受体蛋白（可选：+配体）跑标准 MD。

    参数：
      ligand_smiles : 配体 SMILES（提供则跑复合物 MD，从 SMILES 生成构象）
      ligand_center : 配体初始位置（对接盒中心，单位 Å，[x,y,z]）
      ligand_mol    : RDKit Mol（带 3D 构象），优先使用，直接用其坐标（对接构象）

    返回 {"trajectory_dcd", "system_pdb", "final_pdb", "log", "metrics", "elapsed"}。
    """
    import io
    from openmm.app import (PDBFile, ForceField, Modeller, Simulation,
                            DCDReporter, StateDataReporter, PME, HBonds)
    from pdbfixer import PDBFixer
    from openmm import (LangevinMiddleIntegrator, MonteCarloBarostat, Platform, Vec3)
    from openmm.unit import (nanometer, molar, kelvin, picosecond, femtosecond,
                             bar, kilocalorie_per_mole)
    from config import get_md_params, get_res

    os.makedirs(workdir, exist_ok=True)
    md = get_md_params() if params is None else params
    t0 = time.time()

    # 1) PDBFixer 修复（补缺失原子）
    fixer = PDBFixer(filename=receptor_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # 2) 配体参数化（快速近似：amber 通用类型 + Gasteiger 电荷）
    has_ligand = ligand_mol is not None or ligand_smiles is not None
    ligand_xml = None
    lig_top = None
    lig_pos = None
    lig_n_atoms = 0

    if ligand_mol is not None:
        # 优先使用传入的 mol（带 3D 构象，通常是对接后的构象）
        from agents import ligand_param
        lp = ligand_param.parameterize_ligand_from_mol(ligand_mol, charge_method="gasteiger")
        ligand_xml = lp["xml"]
        lig_top = lp["topology"]
        lig_pos = lp["positions"]
        lig_n_atoms = lp["n_atoms"]
    elif ligand_smiles:
        # 从 SMILES 生成构象（备用方案，可能取向不对）
        from agents import ligand_param
        lp = ligand_param.parameterize_ligand(ligand_smiles, charge_method="gasteiger")
        ligand_xml = lp["xml"]
        lig_top = lp["topology"]
        lig_pos = lp["positions"]
        lig_n_atoms = lp["n_atoms"]
        # 平移：配体几何中心对齐到对接盒中心（Å -> nm）
        if ligand_center is not None and len(ligand_center) == 3:
            cx, cy, cz = ligand_center
            gx = sum(v.x for v in lig_pos) / len(lig_pos)
            gy = sum(v.y for v in lig_pos) / len(lig_pos)
            gz = sum(v.z for v in lig_pos) / len(lig_pos)
            lig_pos = [Vec3(v.x - gx + cx / 10.0, v.y - gy + cy / 10.0,
                            v.z - gz + cz / 10.0) * nanometer for v in lig_pos]

    # 3) 力场
    if ligand_xml:
        ff = ForceField("amber99sbildn.xml", "tip3p.xml", io.StringIO(ligand_xml))
    else:
        ff = ForceField("amber99sbildn.xml", "tip3p.xml")

    # 有配体时启用重试机制（溶剂化随机，偶尔 NaN）
    max_attempts = 10 if has_ligand else 1
    last_error = None

    for attempt in range(max_attempts):
        try:
            result = _run_md_once(fixer, ff, lig_top, lig_pos, md, workdir,
                                   has_ligand, lig_n_atoms, t0)
            if attempt > 0:
                logger.info(f"MD 第 {attempt+1} 次尝试成功")
            return result
        except Exception as e:
            last_error = e
            if has_ligand and attempt < max_attempts - 1:
                logger.warning(f"MD 第 {attempt+1} 次尝试失败: {e}，重试中...")
                gc.collect()
                # 下一次尝试会重新溶剂化（随机水分子位置）
                continue
            raise

    raise last_error  # 不会到这里


def _run_md_once(fixer, ff, lig_top, lig_pos, md, workdir,
                  has_ligand, lig_n_atoms, t0):
    """执行一次完整 MD（内部函数，供重试机制调用）"""
    from openmm.app import (PDBFile, ForceField, Modeller, Simulation,
                            DCDReporter, StateDataReporter, PME, HBonds)
    from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform
    from openmm.unit import (nanometer, molar, kelvin, picosecond, femtosecond,
                             bar, kilocalorie_per_mole)
    from config import get_res
    import os
    import gc
    import time

    platform_obj = Platform.getPlatformByName("CPU")
    platform_obj.setPropertyDefaultValue("Threads", str(get_res()["ncores"]))

    # 4) Modeller：合并配体 + 加氢 + 溶剂化
    #    有配体时用较低离子强度 (0.05 M)，减少 Cl⁻ 近距离接触概率
    modeller = Modeller(fixer.topology, fixer.positions)
    if lig_top is not None:
        modeller.add(lig_top, lig_pos)
    modeller.addHydrogens(ff)
    ionic = 0.05 if has_ligand else md["ionic_strength"]
    modeller.addSolvent(ff, padding=md["padding_nm"] * nanometer,
                        ionicStrength=ionic * molar)

    # 5) 系统 + 积分器
    if has_ligand:
        # 有配体：强摩擦 + 0.5 fs + HBonds 约束
        # 配体用 amber 通用类型近似，参数不够精确，标准设置会 NaN
        friction = 10.0 / picosecond
        dt = 0.5 * femtosecond
    else:
        # 纯蛋白：标准设置
        friction = 1.0 / picosecond
        dt = md["dt_fs"] * femtosecond

    system = ff.createSystem(modeller.topology, nonbondedMethod=PME,
                             nonbondedCutoff=1.0 * nanometer, constraints=HBonds)
    integrator = LangevinMiddleIntegrator(md["temperature"] * kelvin, friction, dt)

    sim = Simulation(modeller.topology, system, integrator, platform_obj)
    sim.context.setPositions(modeller.positions)

    # 6) 能量最小化
    min_iterations = 5000 if has_ligand else 1000
    sim.minimizeEnergy(maxIterations=min_iterations)
    min_energy = sim.context.getState(getEnergy=True).getPotentialEnergy() \
        .value_in_unit(kilocalorie_per_mole)

    # 安全检查：最小化后能量异常则直接失败（触发重试）
    if abs(min_energy) > 1e6:
        raise RuntimeError(f"最小化后能量异常: {min_energy:.1f} kcal/mol")

    # 7) NVT 平衡（逐步升温）
    if has_ligand:
        # 有配体：5 阶段升温，更温和
        warmup_fracs = [0.2, 0.4, 0.6, 0.8, 1.0]
    else:
        warmup_fracs = [0.5, 0.8, 1.0]

    n_stages = len(warmup_fracs)
    steps_per_stage = max(md["nvt_steps"] // n_stages, 100)

    for i, frac in enumerate(warmup_fracs):
        temp = md["temperature"] * frac * kelvin
        sim.context.setVelocitiesToTemperature(temp)
        if i == n_stages - 1:
            remaining = md["nvt_steps"] - steps_per_stage * (n_stages - 1)
            sim.step(max(remaining, 100))
        else:
            sim.step(steps_per_stage)

    # 8) NPT 平衡
    system.addForce(MonteCarloBarostat(1 * bar, md["temperature"] * kelvin))
    sim.context.reinitialize(True)
    sim.step(md["npt_steps"])

    # 9) 生产 MD
    dcd_path = os.path.join(workdir, "trajectory.dcd")
    log_path = os.path.join(workdir, "md.log")
    rep = DCDReporter(dcd_path, md["report_interval"])
    sim.reporters.append(rep)
    sim.reporters.append(StateDataReporter(log_path, md["report_interval"],
                                           step=True, potentialEnergy=True,
                                           temperature=True, volume=True))
    sim.step(md["prod_steps"])
    sim.reporters.remove(rep)
    del rep
    gc.collect()

    # 保存完整体系 PDB
    system_pdb = os.path.join(workdir, "system.pdb")
    with open(system_pdb, "w", encoding="utf-8") as f:
        PDBFile.writeFile(sim.topology,
                          sim.context.getState(getPositions=True).getPositions(), f)

    # 最终结构
    final_pdb = os.path.join(workdir, "md_final.pdb")
    with open(final_pdb, "w", encoding="utf-8") as f:
        PDBFile.writeFile(sim.topology,
                          sim.context.getState(getPositions=True).getPositions(), f)

    elapsed = round(time.time() - t0, 1)
    actual_dt_fs = 0.5 if has_ligand else md["dt_fs"]
    metrics = {
        "min_energy_kcal": round(min_energy, 1),
        "n_atoms": modeller.topology.getNumAtoms(),
        "nvt_steps": md["nvt_steps"],
        "npt_steps": md["npt_steps"],
        "prod_steps": md["prod_steps"],
        "temperature_k": md["temperature"],
        "dt_fs": actual_dt_fs,
        "report_interval": md["report_interval"],
        "has_ligand": has_ligand,
        "ligand_n_atoms": lig_n_atoms,
        "total_time_ps": round((md["nvt_steps"] + md["npt_steps"] + md["prod_steps"])
                               * actual_dt_fs / 1000.0, 2),
        "prod_time_ps": round(md["prod_steps"] * actual_dt_fs / 1000.0, 2),
        "forcefield": "amber99sb-ildn + TIP3P" +
                      (" + 配体(amber通用类型+Gasteiger, HBonds+强摩擦+0.5fs)"
                       if has_ligand else ""),
    }
    return {
        "trajectory_dcd": dcd_path,
        "system_pdb": system_pdb,
        "final_pdb": final_pdb,
        "log": log_path,
        "metrics": metrics,
        "elapsed": elapsed,
    }

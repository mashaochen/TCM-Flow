# -*- coding: utf-8 -*-
"""
xtb 优化 + 性质提取 agent（GFN2-xTB 半经验量子化学）

功能：
  1. 配体几何优化（xtb --opt）
  2. 性质提取（xtb 单点：HOMO-LUMO 能隙、偶极矩、总能量、原子电荷）

不依赖 ORCA；精度为半经验级（GFN2-xTB），适合快速筛选与趋势比较。
"""
import os
import re
import subprocess

from config import XTB_EXE, XTB_GFN, XTB_SOLVENT, get_res


class XtbError(Exception):
    pass


def mol_to_xyz(mol, out_path):
    """RDKit mol（含 3D 坐标）→ xyz 文件（xtb 输入格式）。"""
    conf = mol.GetConformer()
    atoms = mol.GetAtoms()
    lines = [str(len(atoms)), mol.GetProp("_Name") if mol.HasProp("_Name") else "ligand"]
    for i, atom in enumerate(atoms):
        p = conf.GetAtomPosition(i)
        lines.append(f"{atom.GetSymbol():2s}  {p.x:.6f}  {p.y:.6f}  {p.z:.6f}")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return out_path


def _run_xtb(xyz_path, out_dir, args, timeout=900):
    """运行 xtb，返回 (returncode, stdout_text)。"""
    cmd = [XTB_EXE, os.path.abspath(xyz_path)] + list(args)
    ncores = get_res()["ncores"]
    if ncores > 1:
        cmd += ["--parallel", str(ncores)]
    if XTB_SOLVENT:
        cmd += ["--alpb", XTB_SOLVENT]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                           cwd=out_dir, encoding="utf-8", errors="replace")
    except subprocess.TimeoutExpired:
        raise XtbError(f"xtb 超时（>{timeout}s）")
    return r.returncode, (r.stdout or "") + "\n" + (r.stderr or "")


def optimize(xyz_path, out_dir, gfn=None):
    """xtb 几何优化。返回 (优化后 xyz 路径, 优化能量 Eh)。"""
    os.makedirs(out_dir, exist_ok=True)
    gfn = gfn if gfn is not None else XTB_GFN
    rc, out = _run_xtb(xyz_path, out_dir, ["--opt", "--gfn", str(gfn)])
    opt_xyz = os.path.join(out_dir, "xtbopt.xyz")
    if rc != 0 or not os.path.exists(opt_xyz):
        raise XtbError(f"xtb 优化失败：\n{out[-600:]}")
    energy = None
    m = re.search(r"total energy\s+(-?[\d.]+)\s+Eh", out)
    if m:
        energy = float(m.group(1))
    return opt_xyz, energy


def _parse_dipole(out):
    """从 stdout 解析偶极矩（tot Debye）。"""
    # molecular dipole 段里 'full:' 行的最后一列是 tot (Debye)
    m = re.search(r"molecular dipole:.*?full:\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+([\d.]+)",
                  out, re.S)
    if m:
        return float(m.group(4))
    return None


def extract_properties(xyz_path, out_dir, gfn=None):
    """
    xtb 单点计算，提取性质 dict。
    返回 {"total_energy_eh", "homo_lumo_gap_ev", "dipole_debye", ...}
    """
    os.makedirs(out_dir, exist_ok=True)
    gfn = gfn if gfn is not None else XTB_GFN
    rc, out = _run_xtb(xyz_path, out_dir, ["--gfn", str(gfn)])
    if rc != 0:
        raise XtbError(f"xtb 单点失败：\n{out[-600:]}")

    props = {}
    m = re.search(r"total energy\s+(-?[\d.]+)\s+Eh", out)
    if m:
        props["total_energy_eh"] = float(m.group(1))
    m = re.search(r"HOMO-LUMO gap\s+([\d.]+)\s+eV", out)
    if m:
        props["homo_lumo_gap_ev"] = float(m.group(1))
    d = _parse_dipole(out)
    if d is not None:
        props["dipole_debye"] = d

    # 原子电荷（Mulliken）：解析电荷表末尾的列
    charges = _parse_mulliken(out)
    if charges:
        props["mulliken_charges"] = charges

    return props


def _parse_mulliken(out):
    """
    从 stdout 解析 Mulliken 原子电荷表。
    xtb 的 Mulliken 表每行形如： "  1 C   -0.123"（元素后跟电荷）。
    采用宽松匹配：数字编号 + 元素符号 + 电荷数值。
    """
    charges = []
    in_table = False
    for line in out.splitlines():
        if "Mulliken" in line and ("charge" in line.lower() or "populations" in line.lower()):
            in_table = True
            continue
        if in_table:
            if line.strip().startswith("(") or "Wiberg" in line or "covCN" in line:
                break
            m = re.match(r"\s*\d+\s+([A-Za-z]{1,2})\s+(-?[\d.]+)", line)
            if m:
                charges.append({"element": m.group(1), "charge": float(m.group(2))})
    return charges


def run(mol, workdir, gfn=None):
    """
    统一接口：对 RDKit mol 做优化 + 性质提取。
    返回 {"optimized_xyz", "opt_energy_eh", "properties", "files"}。
    """
    os.makedirs(workdir, exist_ok=True)
    mol.SetProp("_Name", "ligand")
    xyz_path = os.path.join(workdir, "ligand.xyz")
    mol_to_xyz(mol, xyz_path)

    opt_xyz, opt_energy = optimize(xyz_path, workdir, gfn=gfn)
    props = extract_properties(opt_xyz, workdir, gfn=gfn)
    props["opt_energy_eh"] = opt_energy

    return {
        "optimized_xyz": opt_xyz,
        "opt_energy_eh": opt_energy,
        "properties": props,
        "files": {"input_xyz": xyz_path, "optimized_xyz": opt_xyz},
    }

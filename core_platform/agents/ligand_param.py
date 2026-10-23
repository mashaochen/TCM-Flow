# -*- coding: utf-8 -*-
"""
配体参数化 agent（快速近似）

给配体小分子生成 OpenMM 力场参数，用于复合物 MD：
  - 原子类型：映射到 amber99sb 的通用类型（CT/CA/C/O/OH/O2/OS/HC/HA/HO/N/N2/N3/CB...）
    从而复用 amber99sb 的通用键/角/二面角参数（≈ GAFF 通用键参数）。
  - 电荷：默认 Gasteiger（RDKit）；后续可切换为 xtb Mulliken。
  - 输出：配体残基 XML（残基名 LIG）+ OpenMM Topology + positions。

说明：这是"快速近似"方案——正式参数化（GAFF + AM1-BCC 或 OpenFF SMIRNOFF）
需要 antechamber/openff-toolkit，当前环境不可用；本方案先跑通复合物 MD，参数可后续替换。
"""
import io

from rdkit import Chem
from rdkit.Chem import AllChem

# amber99sb 里 class -> 数字 atom type（OpenMM 内置力场用数字类型名）
# 注意：H 类型有 49 种上下文依赖变体，这里用最基础的类型 1（酰胺氢）作为近似
CLASS2TYPE = {
    "C": "6", "CA": "20", "CB": "331", "CT": "2",
    "HC": "5", "HA": "267", "HO": "36", "H": "1",  # H=酰胺/胺氢（近似）
    "N": "0", "N2": "18", "N3": "226",
    "O": "7", "O2": "58", "OH": "35", "OS": "1033",
    "S": "87", "SH": "68", "HS": "78",
    "F": "60", "Cl": "61", "Br": "62", "I": "63",
    "P": "146",  # 磷（近似）
}


def _bond_order(mol, a, n):
    b = mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx())
    return b.GetBondTypeAsDouble() if b else 1.0


def amber_class(mol, atom):
    """RDKit 原子 -> amber 类型 class（简化规则，覆盖常见有机基团）。"""
    a = atom
    sym = a.GetSymbol()
    hyb = a.GetHybridization()
    aromatic = a.GetIsAromatic()
    nbrs = list(a.GetNeighbors())

    if sym == "C":
        if any(n.GetSymbol() == "O" and _bond_order(mol, a, n) == 2.0 for n in nbrs):
            return "C"        # 羰基/羧基碳
        if aromatic:
            return "CA"       # 芳香碳
        if hyb == Chem.HybridizationType.SP2:
            return "CB"       # 烯碳
        return "CT"           # sp3 碳
    if sym == "O":
        dc = [n for n in nbrs if n.GetSymbol() == "C" and _bond_order(mol, a, n) == 2.0]
        if dc:
            c = dc[0]
            other_o = [n for n in c.GetNeighbors()
                       if n.GetSymbol() == "O" and n.GetIdx() != a.GetIdx()
                       and _bond_order(mol, c, n) == 1.0]
            return "O2" if other_o else "O"   # 羧基负电氧 / 羰基氧
        nc = sum(1 for n in nbrs if n.GetSymbol() == "C")
        return "OS" if nc >= 2 else "OH"      # 醚/酯氧 / 羟基氧
    if sym == "N":
        if aromatic:
            return "N2"
        if hyb == Chem.HybridizationType.SP2:
            return "N"        # 酰胺氮
        return "N3"           # 胺氮
    if sym == "H":
        c = nbrs[0] if nbrs else None
        if c is None:
            return "H"
        if c.GetSymbol() == "C":
            return "HA" if c.GetIsAromatic() else "HC"
        if c.GetSymbol() == "O":
            return "HO"
        if c.GetSymbol() == "N":
            return "H"
        if c.GetSymbol() == "S":
            return "HS"
        return "H"
    if sym == "S":
        # 硫醇 S-H 用 SH，其余用 S
        if any(n.GetSymbol() == "H" for n in nbrs):
            return "SH"
        return "S"
    return sym


def _charges(mol, method="gasteiger"):
    """返回配体原子电荷列表（顺序 = RDKit 原子顺序）。"""
    if method == "gasteiger":
        AllChem.ComputeGasteigerCharges(mol)
        return [float(a.GetProp("_GasteigerCharge")) for a in mol.GetAtoms()]
    if method == "zero":
        return [0.0 for _ in mol.GetAtoms()]
    raise ValueError(f"未知电荷方法: {method}")


def build_ligand_xml(mol, charges):
    """生成配体残基 XML（残基名 LIG，原子名 A0..An，类型用 amber 数字类型）。"""
    lines = ['<ForceField><Residues><Residue name="LIG">']
    for i, a in enumerate(mol.GetAtoms()):
        cls = amber_class(mol, a)
        t = CLASS2TYPE.get(cls, cls)
        lines.append(f'  <Atom name="A{i}" type="{t}" charge="{charges[i]:.5f}"/>')
    for b in mol.GetBonds():
        lines.append(f'  <Bond atomName1="A{b.GetBeginAtomIdx()}" '
                     f'atomName2="A{b.GetEndAtomIdx()}"/>')
    lines.append('</Residue></Residues></ForceField>')
    return "\n".join(lines)


def build_ligand_topology(mol):
    """生成配体 OpenMM Topology（含键）+ positions（单位 nm）。"""
    from openmm.app import Topology
    from openmm.app.element import Element
    from openmm import Vec3
    from openmm.unit import nanometer

    conf = mol.GetConformer()
    top = Topology()
    chain = top.addChain()
    res = top.addResidue("LIG", chain)
    atoms = []
    positions = []
    for i, a in enumerate(mol.GetAtoms()):
        atoms.append(top.addAtom(f"A{i}", Element.getBySymbol(a.GetSymbol()), res))
        pos = conf.GetAtomPosition(i)
        positions.append(Vec3(pos.x / 10, pos.y / 10, pos.z / 10) * nanometer)
    for b in mol.GetBonds():
        top.addBond(atoms[b.GetBeginAtomIdx()], atoms[b.GetEndAtomIdx()])
    return top, positions


def parameterize_ligand_from_mol(mol, charge_method="gasteiger"):
    """
    从带 3D 坐标的 RDKit mol 做参数化（直接使用传入的构象）。
    
    mol: RDKit Mol（必须有 3D 构象和显式氢）
    charge_method: "gasteiger" 或 "zero"
    
    返回 dict：{"xml", "topology", "positions", "n_atoms", "charge_method"}。
    """
    if mol.GetNumConformers() == 0:
        raise ValueError("mol 没有 3D 构象")
    
    # 确保有显式氢
    if not any(a.GetSymbol() == 'H' for a in mol.GetAtoms()):
        mol = Chem.AddHs(mol, addCoords=True)
        if mol.GetNumConformers() == 0:
            p = AllChem.ETKDGv3()
            p.randomSeed = 42
            AllChem.EmbedMolecule(mol, p)
    
    charges = _charges(mol, method=charge_method)
    xml = build_ligand_xml(mol, charges)
    topology, positions = build_ligand_topology(mol)
    
    return {
        "xml": xml,
        "topology": topology,
        "positions": positions,
        "n_atoms": mol.GetNumAtoms(),
        "charge_method": charge_method,
    }


def parameterize_ligand(smiles, charge_method="gasteiger"):
    """
    从 SMILES 生成配体 3D 构象并做参数化。
    
    返回 dict：{"xml", "topology", "positions", "n_atoms", "charge_method"}。
    """
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol is None:
        raise ValueError(f"SMILES 解析失败: {smiles}")
    p = AllChem.ETKDGv3()
    p.randomSeed = 42
    if AllChem.EmbedMolecule(mol, p) != 0:
        p.useRandomCoords = True
        AllChem.EmbedMolecule(mol, p)
    AllChem.MMFFOptimizeMolecule(mol)
    
    return parameterize_ligand_from_mol(mol, charge_method=charge_method)

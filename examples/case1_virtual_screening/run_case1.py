# -*- coding: utf-8 -*-
"""
案例一：7种黄酮对 2J0L 虚拟筛选
直接调用 docking 模块，快速验证排序合理性
"""
import sys
import os

# 添加核心平台路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "core_platform"))

import docking
from config import RESULTS_DIR

# 7个黄酮类成分
compounds = [
    ("柚皮素", "CC(CC1=CC(=C(C=C1)O)O)C(=O)C2=CC=C(C=C2)O"),
    ("山柰酚", "C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O"),
    ("槲皮素", "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O"),
    ("黄芩素", "C1=CC(=C(C=C1C2=CC(=O)C3=C(O2)C=C(C=C3O)O)O)O"),
    ("芹菜素", "C1=CC(=CC=C1C2=CC(=O)C3=C(O2)C=C(C=C3O)O)O"),
    ("木犀草素", "C1=CC(=C(C=C1C2=CC(=O)C3=C(O2)C=C(C=C3O)O)O)O"),
    ("葛根素", "C1=CC(=C(C=C1C2=CC(=O)C3=C(O2)C=C(C=C3O)O)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O"),
]

# 从 RCSB 下载 2J0L
print("正在下载 2J0L 蛋白结构...")
import urllib.request
pdb_url = "https://files.rcsb.org/download/2J0L.pdb"
req = urllib.request.Request(pdb_url, headers={"User-Agent": "TCMFlow/1.0"})
with urllib.request.urlopen(req, timeout=30) as resp:
    pdb_2j0l = resp.read().decode('utf-8')
print(f"  下载完成，{len(pdb_2j0l)} 字符\n")

print("=" * 70)
print("  案例一：7种中药黄酮类成分对 FAK (PDB: 2J0L) 虚拟筛选")
print("=" * 70)

results = []
for i, (name, smiles) in enumerate(compounds):
    print(f"\n[{i+1}/7] 正在对接 {name}...")
    try:
        work_dir = os.path.join(RESULTS_DIR, "vscreen_2j0l", f"cmpd_{i+1:02d}")
        os.makedirs(work_dir, exist_ok=True)

        result = docking.run_docking(
            smiles=smiles,
            ligand_name=name,
            protein_pdb_text=pdb_2j0l,
            pdb_id="2J0L",
            output_dir=work_dir,
            dock_params={}
        )

        affinity = result.get("best_affinity", 0)
        num_atoms = result.get("ligand_info", {}).get("num_atoms", 0)
        mw = result.get("ligand_info", {}).get("mw", 0)

        results.append({
            "name": name,
            "affinity": affinity,
            "num_atoms": num_atoms,
            "mw": mw
        })

        print(f"    ✓ 完成：亲和力 = {affinity:.2f} kcal/mol，分子量 = {mw:.2f}")

    except Exception as e:
        print(f"    ✗ 失败：{e}")
        results.append({"name": name, "affinity": None, "error": str(e)})

# 排序
valid = [r for r in results if r.get("affinity") is not None]
valid.sort(key=lambda x: x["affinity"])

print("\n" + "=" * 70)
print("  虚拟筛选结果（按亲和力由强到弱排序）")
print("=" * 70)
print(f"{'排名':<6}{'成分':<12}{'亲和力(kcal/mol)':<18}{'分子量':<12}{'原子数':<8}")
print("-" * 60)
for i, r in enumerate(valid):
    print(f"{i+1:<6}{r['name']:<12}{r['affinity']:<18.2f}{r['mw']:<12.2f}{r['num_atoms']:<8}")

print(f"\n结果目录: {os.path.join(RESULTS_DIR, 'vscreen_2j0l')}")
print("\n说明：对接打分仅用于相对排序，并非实验结合自由能。")

# -*- coding: utf-8 -*-
"""
案例一：7种黄酮对 2J0L 虚拟筛选
直接调用 docking 模块，快速验证排序合理性

使用方法：
    python run_case1.py

注意：
    - 脚本会自动从 RCSB 下载 2J0L 蛋白结构
    - 如网络不通，可手动下载 2J0L.pdb 放到本目录，脚本会优先读取本地文件
    - 也可以使用其他 PDB ID：python run_case1.py 1CRN
"""
import sys
import os

# 添加核心平台路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "core_platform"))

import docking
from config import RESULTS_DIR

# 7个黄酮类成分（与论文一致）
compounds = [
    ("柚皮素", "CC(CC1=CC(=C(C=C1)O)O)C(=O)C2=CC=C(C=C2)O"),
    ("山柰酚", "C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O"),
    ("槲皮素", "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O"),
    ("黄芩素", "C1=CC(=C(C=C1C2=CC(=O)C3=C(O2)C=C(C=C3O)O)O)O"),
    ("芹菜素", "C1=CC(=CC=C1C2=CC(=O)C3=C(O2)C=C(C=C3O)O)O"),
    ("木犀草素", "C1=CC(=C(C=C1C2=CC(=O)C3=C(O2)C=C(C=C3O)O)O)O"),
    ("葛根素", "C1=CC(=C(C=C1C2=CC(=O)C3=C(O2)C=C(C=C3O)O)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O"),
]


def download_pdb(pdb_id):
    """下载 PDB 文件，支持多个镜像源和本地缓存"""
    pdb_id = pdb_id.upper()
    local_file = os.path.join(os.path.dirname(__file__), f"{pdb_id}.pdb")

    # 1. 优先读本地文件
    if os.path.exists(local_file):
        print(f"  使用本地文件: {local_file}")
        with open(local_file, "r", encoding="utf-8") as f:
            return f.read()

    # 2. 尝试多个下载源
    import urllib.request
    urls = [
        f"https://files.rcsb.org/download/{pdb_id}.pdb",
        f"https://files.rcsb.org/pub/pdb/data/structures/all/pdb/pdb{pdb_id.lower()}.ent.gz",
    ]

    for url in urls:
        try:
            print(f"  正在从 {url.split('/')[2]} 下载...")
            req = urllib.request.Request(url, headers={"User-Agent": "TCMFlow/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read()
                # 处理 gzip 压缩
                if url.endswith(".gz"):
                    import gzip
                    data = gzip.decompress(data)
                text = data.decode("utf-8")
                # 保存到本地缓存
                with open(local_file, "w", encoding="utf-8") as f:
                    f.write(text)
                print(f"  下载完成，{len(text)} 字符，已缓存到本地")
                return text
        except Exception as e:
            print(f"  下载失败: {e}")
            continue

    raise RuntimeError(
        f"无法下载 PDB 文件 {pdb_id}。\n"
        f"请手动从 https://www.rcsb.org/structure/{pdb_id} 下载，\n"
        f"保存为 {local_file} 后重新运行。"
    )


def main():
    # 支持命令行指定 PDB ID
    pdb_id = sys.argv[1] if len(sys.argv) > 1 else "2J0L"

    print(f"正在获取 {pdb_id} 蛋白结构...")
    pdb_text = download_pdb(pdb_id)
    print()

    print("=" * 70)
    print(f"  案例一：7种中药黄酮类成分对 FAK (PDB: {pdb_id}) 虚拟筛选")
    print("=" * 70)

    results = []
    for i, (name, smiles) in enumerate(compounds):
        print(f"\n[{i+1}/7] 正在对接 {name}...")
        try:
            work_dir = os.path.join(RESULTS_DIR, f"vscreen_{pdb_id.lower()}", f"cmpd_{i+1:02d}")
            os.makedirs(work_dir, exist_ok=True)

            result = docking.run_docking(
                smiles=smiles,
                ligand_name=name,
                protein_pdb_text=pdb_text,
                pdb_id=pdb_id,
                output_dir=work_dir,
                dock_params={},
            )

            affinity = result.get("best_affinity", 0)
            num_atoms = result.get("ligand_info", {}).get("num_atoms", 0)
            mw = result.get("ligand_info", {}).get("mw", 0)

            results.append(
                {
                    "name": name,
                    "affinity": affinity,
                    "num_atoms": num_atoms,
                    "mw": mw,
                }
            )

            print(f"    ✓ 完成：亲和力 = {affinity:.2f} kcal/mol，分子量 = {mw:.2f}")

        except Exception as e:
            print(f"    ✗ 失败：{e}")
            import traceback
            traceback.print_exc()
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

    print(f"\n结果目录: {os.path.join(RESULTS_DIR, f'vscreen_{pdb_id.lower()}')}")
    print("\n说明：对接打分仅用于相对排序，并非实验结合自由能。")
    print("      不同 Vina 版本、随机种子、初始构象可能导致打分差异。")


if __name__ == "__main__":
    main()

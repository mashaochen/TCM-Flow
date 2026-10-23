# -*- coding: utf-8 -*-
"""
PubChem 工具：通过 PubChem PUG REST API 获取化合物的标准 3D SDF 结构。

API 文档：https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
速率限制：每秒不超过 5 次请求
"""
import requests
import time

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# 简单的速率控制（每秒最多 5 次）
_last_request = 0.0


def _rate_limit():
    """简单的速率限制，确保每秒不超过 5 次请求。"""
    global _last_request
    now = time.time()
    elapsed = now - _last_request
    if elapsed < 0.25:  # 1/5 = 0.2秒，留余量
        time.sleep(0.25 - elapsed)
    _last_request = time.time()


def get_3d_sdf_by_name(name, timeout=30):
    """通过化合物名称从 PubChem 获取标准 3D SDF 结构。
    
    返回 (sdf_text, cid) 或 (None, None) 失败。
    """
    if not name or not name.strip():
        return None, None
    
    _rate_limit()
    
    url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name.strip())}/SDF"
    params = {"record_type": "3d"}
    
    try:
        r = requests.get(url, params=params, timeout=timeout)
        if r.status_code == 200 and "V2000" in r.text:
            # 从响应头或内容提取 CID
            cid = None
            # SDF 第一行之后可能有 CID 信息
            lines = r.text.strip().split("\n")
            for line in lines[:5]:
                if line.strip().isdigit():
                    cid = int(line.strip())
                    break
            return r.text, cid
        else:
            return None, None
    except Exception:
        return None, None


def get_3d_sdf_by_smiles(smiles, timeout=30):
    """通过 SMILES 从 PubChem 获取标准 3D SDF 结构。
    
    返回 (sdf_text, cid) 或 (None, None) 失败。
    """
    if not smiles or not smiles.strip():
        return None, None
    
    _rate_limit()
    
    url = f"{PUBCHEM_BASE}/compound/smiles/{requests.utils.quote(smiles.strip())}/SDF"
    params = {"record_type": "3d"}
    
    try:
        r = requests.get(url, params=params, timeout=timeout)
        if r.status_code == 200 and "V2000" in r.text:
            cid = None
            lines = r.text.strip().split("\n")
            for line in lines[:5]:
                if line.strip().isdigit():
                    cid = int(line.strip())
                    break
            return r.text, cid
        else:
            return None, None
    except Exception:
        return None, None


def sdf_to_3d_mol(sdf_text):
    """从 SDF 文本读取第一个构象，返回 RDKit Mol（带 3D 坐标和显式氢）。"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        supplier = Chem.SDMolSupplier()
        supplier.SetData(sdf_text)
        mol = next(supplier, None)
        if mol and mol.GetNumConformers() > 0:
            # 加显式氢（保留原始坐标）
            mol = Chem.AddHs(mol, addCoords=True)
            # 如果加氢后坐标有问题（比如氢位置不正确），用MMFF94优化一下
            if mol.GetNumConformers() > 0:
                try:
                    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
                except Exception:
                    pass
            return mol
        return None
    except Exception:
        return None

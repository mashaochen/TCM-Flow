# -*- coding: utf-8 -*-
"""
中药数据读取模块（只读，从tcm.db获取中药和活性成分数据）
"""
import sqlite3
import config
from compound_smiles import get_smiles as get_smiles_from_dict

_db_path = config.TCM_DB_PATH


def _get_conn():
    conn = sqlite3.connect(_db_path)
    conn.row_factory = sqlite3.Row
    return conn


def list_herbs(keyword=None, limit=100):
    """获取中药列表，支持关键词搜索，精确匹配优先"""
    conn = _get_conn()
    try:
        if keyword:
            kw = f"%{keyword}%"
            # 精确匹配排在最前（name = keyword），然后是前缀匹配，最后是包含匹配
            rows = conn.execute(
                "SELECT name, pinyin, nature, flavor, channel, efficacy, indications "
                "FROM herbs WHERE name LIKE ? "
                "ORDER BY "
                "  CASE WHEN name = ? THEN 0 "
                "       WHEN name LIKE ? THEN 1 "
                "       ELSE 2 END, "
                "  name "
                "LIMIT ?",
                (kw, keyword, keyword + '%', limit)
            ).fetchall()
        else:
            rows = conn.execute(
                "SELECT name, pinyin, nature, flavor, channel, efficacy, indications "
                "FROM herbs ORDER BY name LIMIT ?",
                (limit,)
            ).fetchall()
        return [dict(r) for r in rows]
    finally:
        conn.close()


def get_herb_detail(name):
    """获取单味中药详情"""
    conn = _get_conn()
    try:
        row = conn.execute(
            "SELECT * FROM herbs WHERE name = ?", (name,)
        ).fetchone()
        return dict(row) if row else None
    finally:
        conn.close()


def get_herb_compounds(herb_name):
    """获取某中药的活性成分列表"""
    conn = _get_conn()
    try:
        # 从 tcmnp_herb_target 表获取（OB>=30, DL>=0.18 的高质量数据）
        # 按化合物分组，统计靶点数量
        rows = conn.execute(
            "SELECT Molecule_name as compound_name, MOL_ID, OB, DL, "
            "COUNT(DISTINCT Target) as target_count "
            "FROM tcmnp_herb_target "
            "WHERE Herb_cn_name = ? "
            "GROUP BY Molecule_name, MOL_ID, OB, DL "
            "ORDER BY target_count DESC, OB DESC",
            (herb_name,)
        ).fetchall()
        compounds = [dict(r) for r in rows]
        
        # 补充分子性质（从 active_compounds 表）
        for c in compounds:
            prop = conn.execute(
                "SELECT compound_name, chinese_name, cas_no, formula, mol_weight, "
                "logp, druglikeness, pubchem_id "
                "FROM active_compounds WHERE compound_name = ? OR chinese_name = ? "
                "LIMIT 1",
                (c["compound_name"], c["compound_name"])
            ).fetchone()
            if prop:
                c.update(dict(prop))
            
            # 补充 SMILES
            c["smiles"] = get_compound_smiles(c["compound_name"], c.get("chinese_name"))
        
        return compounds
    finally:
        conn.close()


def get_compound_targets(herb_name, compound_name):
    """获取某中药中某成分的靶点列表"""
    conn = _get_conn()
    try:
        rows = conn.execute(
            "SELECT DISTINCT Target as gene_symbol "
            "FROM tcmnp_herb_target "
            "WHERE Herb_cn_name = ? AND Molecule_name = ? "
            "ORDER BY Target",
            (herb_name, compound_name)
        ).fetchall()
        return [r["gene_symbol"] for r in rows]
    finally:
        conn.close()


def get_compound_smiles(compound_name, chinese_name=None):
    """获取成分的SMILES
    
    查找顺序：
    1. 本地SMILES字典（compound_smiles.py）- 按中文名
    2. 本地SMILES字典 - 按英文名
    3. 返回 None（前端显示"无结构"）
    """
    # 先查中文名
    if chinese_name:
        smi = get_smiles_from_dict(chinese_name)
        if smi:
            return smi
    
    # 再查英文名
    smi = get_smiles_from_dict(compound_name)
    if smi:
        return smi
    
    return None


def get_herb_targets(herb_name):
    """获取某中药所有成分的靶点并集"""
    conn = _get_conn()
    try:
        rows = conn.execute(
            "SELECT DISTINCT Target as gene_symbol "
            "FROM tcmnp_herb_target WHERE Herb_cn_name = ? ORDER BY Target",
            (herb_name,)
        ).fetchall()
        return [r["gene_symbol"] for r in rows]
    finally:
        conn.close()


def kegg_enrichment(gene_list, top=15):
    """简化的KEGG通路富集分析（超几何检验）"""
    conn = _get_conn()
    try:
        # 获取背景基因集大小
        bg_row = conn.execute(
            "SELECT COUNT(DISTINCT gene) as total FROM tcmnp_kegg"
        ).fetchone()
        bg_total = bg_row["total"] if bg_row else 20000
        
        input_set = set(gene_list)
        input_size = len(input_set)
        
        if input_size < 2:
            return []
        
        # 获取每个通路的基因
        rows = conn.execute(
            "SELECT term, gene FROM tcmnp_kegg"
        ).fetchall()
        
        from collections import defaultdict
        pathway_genes = defaultdict(set)
        for r in rows:
            pathway_genes[r["term"]].add(r["gene"])
        
        # 超几何检验
        from math import comb
        results = []
        for term, genes in pathway_genes.items():
            overlap = input_set & genes
            if len(overlap) < 2:
                continue
            
            k = len(overlap)  # 命中数
            M = bg_total       # 背景基因总数
            K = len(genes)     # 通路基因数
            n = input_size     # 输入基因数
            
            # 超几何p值（上尾）
            p_value = 0.0
            max_i = min(n, K)
            for i in range(k, max_i + 1):
                try:
                    p_value += comb(K, i) * comb(M - K, n - i) / comb(M, n)
                except (ValueError, ZeroDivisionError, OverflowError):
                    pass
            
            results.append({
                "term": term,
                "overlap": k,
                "pathway_size": K,
                "input_size": input_size,
                "p_value": p_value,
                "genes": sorted(list(overlap))
            })
        
        # 按p值排序，取top
        results.sort(key=lambda x: x["p_value"])
        return results[:top]
    finally:
        conn.close()

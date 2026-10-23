# -*- coding: utf-8 -*-
"""
常见中药活性成分 SMILES 字典
用于从成分名称映射到SMILES，供TCM-Flow计算使用
"""

# 中文名 -> SMILES 映射
# 数据来源：PubChem / 文献
COMPOUND_SMILES = {
    # ===== 黄酮类（论文中使用的7种）=====
    "柚皮素": "O=C1CC(c2ccc(O)cc2)Oc2cc(O)ccc12",
    "木犀草素": "O=c1c2cc(O)cc(O)c2oc2ccc(O)c(O)c12",
    "芹菜素": "O=c1c2cc(O)cc(O)c2oc2ccc(O)cc12",
    "黄芩素": "O=c1c2cc(O)c(O)c(O)c2oc2cccc(O)c12",
    "山柰酚": "O=c1c2cc(O)cc(O)c2oc2c(O)cc(O)c12",
    "槲皮素": "O=c1c2cc(O)cc(O)c2oc2c(O)c(O)cc(O)c12",
    "染料木素": "O=c1cc(O)c2cc(O)ccc2o1",
    
    # ===== 其他常见黄酮类 =====
    "黄芩苷": "O=c1c2cc(O)c(O[C@@H]3O[C@H](C(=O)O)[C@@H](O)[C@H]3O)c(O)c2oc2cccc(O)c12",
    "汉黄芩素": "COc1c(O)cc2c(=O)c3cccc(O)c3oc2c1O",
    "汉黄芩苷": "COc1c(O)cc2c(=O)c3cccc(O[C@@H]4O[C@H](C(=O)O)[C@@H](O)[C@H]4O)c3oc2c1O",
    "千层纸素A": "COc1c(O)cc2c(=O)c3c(O)cccc3oc2c1O",
    "野黄芩苷": "O=c1c2cc(O)c(O[C@@H]3O[C@H](C(=O)O)[C@@H](O)[C@H]3O)c(O)c2oc2c(O)ccc12",
    "橙皮苷": "CC1=C(C(=O)C2=C3C(=O)CC(C4=CC=CC=C43)(C)OC2=C1C)C",
    "芦丁": "",
    "甘草素": "O=C1CC(c2ccc(O)cc2)Oc2ccc(O)cc12",
    "异甘草素": "O=C(/C=C/c1ccc(O)cc1)c1ccc(O)cc1O",
    
    # ===== 生物碱类 =====
    "小檗碱": "COc1ccc2c(c1)Cc1c3ccc(OC)c(OC)c3cc[n+]1C2",
    "黄连素": "COc1ccc2c(c1)Cc1c3ccc(OC)c(OC)c3cc[n+]1C2",
    "巴马汀": "COc1ccc2c(c1)Cc1cc3ccc(OC)c(OC)c3cc[n+]1C2",
    "黄连碱": "",
    "川芎嗪": "Cc1nc(C)c(C)nc1C",
    "苦参碱": "C1CC2=CC(=O)C3=C(C2CCN1C3=O)O",
    "麻黄碱": "CN[C@@H](C)[C@H](O)c1ccccc1",
    
    # ===== 三萜类 =====
    "人参皂苷Rg1": "",
    "人参皂苷Rb1": "",
    "熊果酸": "CC1(C)CCC2(C)C(CCC3(C)C2CCC2(C)C3CCC2(O)C(=O)O)C1",
    "齐墩果酸": "CC1(C)CCC2(C)C(CCC3(C)C2CCC2(C)C3CCC2C(=C)C(=O)O)C1",
    
    # ===== 酚酸类 =====
    "丹参酮IIA": "CC1=C(C(=O)C2=C3C(=O)CC(C4=CC=CC=C43)(C)OC2=C1C)C",
    "丹酚酸B": "",
    "阿魏酸": "C=CC(=O)c1ccc(O)cc1O",
    "绿原酸": "O=C(O)C(CC(=O)O)OC(=O)/C=C/c1ccc(O)c(O)c1",
    
    # ===== 其他常见成分 =====
    "姜黄素": "O=C(/C=C/c1ccc(O)c(O)c1)/C=C/c1ccc(O)c(O)c1",
    "大黄素": "CC1=C2C(=O)c3ccc(O)c3C(=O)C2=C(O)c2c(O)ccc(O)c2C1=O",
    "葛根素": "O=c1c2cc(O)c(OC[C@H]3O[C@@H](CO)[C@@H](O)[C@H]3O)c(O)c2oc2ccc(O)c12",
    "穿心莲内酯": "",
    "青蒿素": "",
    "紫杉醇": "",
    
    # ===== 英文名映射 =====
    "Naringenin": "O=C1CC(c2ccc(O)cc2)Oc2cc(O)ccc12",
    "Luteolin": "O=c1c2cc(O)cc(O)c2oc2ccc(O)c(O)c12",
    "Apigenin": "O=c1c2cc(O)cc(O)c2oc2ccc(O)cc12",
    "Baicalein": "O=c1c2cc(O)c(O)c(O)c2oc2cccc(O)c12",
    "Kaempferol": "O=c1c2cc(O)cc(O)c2oc2c(O)cc(O)c12",
    "Quercetin": "O=c1c2cc(O)cc(O)c2oc2c(O)c(O)cc(O)c12",
    "Genistein": "O=c1cc(O)c2cc(O)ccc2o1",
    "Wogonin": "COc1c(O)cc2c(=O)c3cccc(O)c3oc2c1O",
    "Baicalin": "O=c1c2cc(O)c(O[C@@H]3O[C@H](C(=O)O)[C@@H](O)[C@H]3O)c(O)c2oc2cccc(O)c12",
    "Wogonoside": "COc1c(O)cc2c(=O)c3cccc(O[C@@H]4O[C@H](C(=O)O)[C@@H](O)[C@H]4O)c3oc2c1O",
    "Oroxylin A": "COc1c(O)cc2c(=O)c3c(O)cccc3oc2c1O",
    "Oroxylin-A": "COc1c(O)cc2c(=O)c3c(O)cccc3oc2c1O",
    "Scutellarin": "O=c1c2cc(O)c(O[C@@H]3O[C@H](C(=O)O)[C@@H](O)[C@H]3O)c(O)c2oc2c(O)ccc12",
    "Liquiritigenin": "O=C1CC(c2ccc(O)cc2)Oc2ccc(O)cc12",
    "Isoliquiritigenin": "O=C(/C=C/c1ccc(O)cc1)c1ccc(O)cc1O",
    "Berberine": "COc1ccc2c(c1)Cc1c3ccc(OC)c(OC)c3cc[n+]1C2",
    "Palmatine": "COc1ccc2c(c1)Cc1cc3ccc(OC)c(OC)c3cc[n+]1C2",
    "Tetramethylpyrazine": "Cc1nc(C)c(C)nc1C",
    "Ligustrazine": "Cc1nc(C)c(C)nc1C",
    "Matrine": "C1CC2=CC(=O)C3=C(C2CCN1C3=O)O",
    "Ephedrine": "CN[C@@H](C)[C@H](O)c1ccccc1",
    "Curcumin": "O=C(/C=C/c1ccc(O)c(O)c1)/C=C/c1ccc(O)c(O)c1",
    "Emodin": "CC1=C2C(=O)c3ccc(O)c3C(=O)C2=C(O)c2c(O)ccc(O)c2C1=O",
    "Ursolic acid": "CC1(C)CCC2(C)C(CCC3(C)C2CCC2(C)C3CCC2(O)C(=O)O)C1",
    "Oleanolic acid": "CC1(C)CCC2(C)C(CCC3(C)C2CCC2(C)C3CCC2C(=C)C(=O)O)C1",
    "Ferulic acid": "C=CC(=O)c1ccc(O)cc1O",
    "Chlorogenic acid": "O=C(O)C(CC(=O)O)OC(=O)/C=C/c1ccc(O)c(O)c1",
    "beta-sitosterol": "CC1CCC2C(C1)CCC3C2CCC4=CC(CC34)C5CCC(C)CC5(C)CO",
    "beta Sitosterol": "CC1CCC2C(C1)CCC3C2CCC4=CC(CC34)C5CCC(C)CC5(C)CO",
    "Stigmasterol": "CC1CCC2C(C1)CCC3C2CCC4=CC(CC34)C5CCC(C)C(=CC)C5(C)CO",
    "Palmitic acid": "CCCCCCCCCCCCCCCC(=O)O",
    "Hexadecanoic acid": "CCCCCCCCCCCCCCCC(=O)O",
    "Tanshinone IIA": "CC1=C(C(=O)C2=C3C(=O)CC(C4=CC=CC=C43)(C)OC2=C1C)C",
    "Tanshinone IIA sulfonate": "CC1=C(C(=O)C2=C3C(=O)CC(C4=CC=CC=C43)(C)OC2=C1C)C",
    "Puerarin": "O=c1c2cc(O)c(OC[C@H]3O[C@@H](CO)[C@@H](O)[C@H]3O)c(O)c2oc2ccc(O)c12",
    "Hesperidin": "",
    "Rutin": "",
    "Ginsenoside Rg1": "",
    "Ginsenoside Rb1": "",
    "Salvianolic acid B": "",
    "Andrographolide": "",
    "Artemisinin": "",
    "Paclitaxel": "",
    "Caffeic acid": "C=CC(=O)c1ccc(O)c(O)c1",
    "p-Coumaric acid": "C=CC(=O)c1ccc(O)cc1",
    "Gallic acid": "O=C(O)c1cc(O)c(O)c(O)c1",
}


def get_smiles(compound_name):
    """根据成分名称获取SMILES
    
    优先匹配中文名，再匹配英文名，最后返回None
    大小写不敏感
    """
    if not compound_name:
        return None
    
    name = compound_name.strip()
    name_lower = name.lower()
    
    # 直接匹配（大小写不敏感）
    for cn, smi in COMPOUND_SMILES.items():
        if not smi:
            continue
        if cn.lower() == name_lower:
            return smi
    
    # 模糊匹配（包含关系，大小写不敏感）
    for cn, smi in COMPOUND_SMILES.items():
        if not smi:
            continue
        cn_lower = cn.lower()
        if cn_lower in name_lower or name_lower in cn_lower:
            return smi
    
    return None

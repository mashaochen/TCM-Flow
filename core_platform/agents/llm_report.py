# -*- coding: utf-8 -*-
"""
LLM 报告 agent —— 调用本地 Ollama 生成分析报告 report.txt

把全流程的数值结果（对接、xtb 性质、MD 指标）+ 蛋白真实信息，喂给本地大模型，
生成一段中文分析，回答四类问题：
  1. 分子性质（HOMO-LUMO 能隙、偶极矩、电荷等）预示了什么？
  2. 与蛋白的对接亲和力代表什么（结合强弱判断）？
  3. 这个蛋白属于什么类别 / 关联什么疾病与通路？
  4. 分子动力学指标（RMSD / Rg / RMSF）说明了什么？

抗幻觉设计（重点）：
  - 蛋白信息从 PDB 头部真实解析（HEADER/TITLE/COMPND/KEYWDS），不靠模型猜；
  - 数值解读给出教科书级参考区间，让模型基于具体数值引用，而非自由发挥；
  - 疾病/通路关联强制标注为"推测"，禁止编造具体通路编号或未经证实的疾病名；
  - 温度调低（0.3），减少随机发挥。
"""
import os
import json
import urllib.request

try:
    from config import (OLLAMA_URL, OLLAMA_MODEL, OLLAMA_TIMEOUT, OLLAMA_TEMPERATURE)
except Exception:
    OLLAMA_URL = os.environ.get("OLLAMA_URL", "http://localhost:11434")
    OLLAMA_MODEL = os.environ.get("OLLAMA_MODEL", "qwen2.5:7b-instruct")
    OLLAMA_TIMEOUT = int(os.environ.get("OLLAMA_TIMEOUT", "600"))
    OLLAMA_TEMPERATURE = 0.3

TIMEOUT = OLLAMA_TIMEOUT


# ---------------------------------------------------------------------------
# 1) 蛋白身份接地：从 PDB 头部解析真实信息
# ---------------------------------------------------------------------------
def parse_pdb_header(pdb_text):
    """从 PDB 文本头部解析蛋白真实身份信息，返回 dict。"""
    info = {
        "header": "",       # HEADER 分类，如 "HYDROLASE"
        "title": "",        # TITLE 标题
        "compound": "",     # COMPND 分子名（去分号合并）
        "keywords": "",     # KEYWDS 关键词
        "ec": "",           # EC 编号
        "source": "",       # SOURCE 来源生物
        "n_residues": 0,    # 残基数
    }
    if not pdb_text:
        return info

    comp_lines = []
    src_lines = []
    for line in pdb_text.splitlines():
        rec = line[:6].strip()
        if rec == "HEADER" and not info["header"]:
            info["header"] = line[10:50].strip()
        elif rec == "TITLE" and not info["title"]:
            info["title"] = line[10:].strip()
        elif rec == "COMPND":
            txt = line[10:].strip()
            if txt.startswith("EC:"):
                info["ec"] = txt.split(":", 1)[1].strip().rstrip(";")
            else:
                # 只保留 MOLECULE 行
                if "MOLECULE" in txt:
                    comp_lines.append(txt.split(":", 1)[1].strip().rstrip(";"))
        elif rec == "KEYWDS":
            info["keywords"] += line[10:].strip()
        elif rec == "SOURCE":
            txt = line[10:].strip()
            if "ORGANISM_SCIENTIFIC" in txt:
                src_lines.append(txt.split(":", 1)[1].strip().rstrip(";"))
    info["compound"] = "; ".join(comp_lines)
    info["keywords"] = info["keywords"].strip()
    info["source"] = "; ".join(src_lines)
    info["n_residues"] = sum(
        1 for l in pdb_text.splitlines()
        if l.startswith("ATOM") and l[22:26].strip() == "CA"
    )
    return info


def fetch_uniprot_summary(pdb_id, timeout=20):
    """
    尝试通过 RCSB → UniProt 拿蛋白的功能/疾病注释（可选增强接地）。
    失败时返回空 dict，不阻塞主流程。
    """
    out = {"function": "", "diseases": []}
    if not pdb_id:
        return out
    try:
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
        req = urllib.request.Request(url, headers={"User-Agent": "TCMFlow/1.0"})
        data = json.loads(urllib.request.urlopen(req, timeout=timeout).read())
        uniprot_id = next(iter(data.get(pdb_id.lower(), {}).get("UniProt", {})), None)
        if not uniprot_id:
            return out
        url2 = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        req2 = urllib.request.Request(url2, headers={"User-Agent": "TCMFlow/1.0"})
        up = json.loads(urllib.request.urlopen(req2, timeout=timeout).read())
        for c in up.get("comments", []):
            if c.get("commentType") == "FUNCTION":
                out["function"] = " ".join(
                    t.get("value", "") for t in c.get("texts", [])
                )[:500]
        for c in up.get("comments", []):
            if c.get("commentType") == "DISEASE":
                for d in c.get("disease", {}).values():
                    if isinstance(d, dict) and d.get("diseaseId"):
                        out["diseases"].append(d["diseaseId"][:80])
        out["diseases"] = out["diseases"][:5]
    except Exception:
        pass
    return out


# ---------------------------------------------------------------------------
# 2) 数值解读参考（教科书级，让模型引用而非乱编）
# ---------------------------------------------------------------------------
_REFERENCE = """【数值解读参考区间（请据此引用，不要超出）】
- HOMO-LUMO 能隙：<2 eV 化学活泼（软分子，易发生电荷转移/易氧化）；2~4 eV 中等；
  4 eV 以上较稳定（硬分子）。多数口服药物在 3~5 eV，能隙越小通常反应活性越高、越不稳定。
- 偶极矩：>3 Debye 极性大（利于水溶性与氢键形成）；1.5~3 Debye 中等；<1.5 Debye 非极性（脂溶性好、易跨膜）。
- 对接亲和力（AutoDock Vina，越负越强）：<-9 kcal/mol 强结合；-7~-9 较强；-5~-7 中等；>-5 弱结合。
  这是经验打分函数预测值，用于相对排序，不是实验结合自由能。
- RMSD（均方根偏差）：均值 <0.3 nm 通常认为体系在模拟中结构稳定；数值越小越稳定。
- Rg（回转半径）：描述蛋白整体紧凑度，波动越小结构越稳定，显著变化提示构象变化/部分去折叠。
- RMSF（均方根涨落）：描述每个残基的柔性；高 RMSF 区域通常是柔性 loop 或活性位点，低 RMSF 是刚性核心。"""


# ---------------------------------------------------------------------------
# 3) 组装 prompt 并调用 Ollama
# ---------------------------------------------------------------------------
def build_prompt(stages, smiles, ligand_name, pdb_id, protein_info, uniprot_info):
    """把各阶段结果组装成结构化 prompt。stages 为 pipeline 的 result["stages"]。"""
    dk = stages.get("docking", {})
    rd = stages.get("redocking", {})
    xt = stages.get("xtb_opt", {})
    pr = stages.get("properties", {})
    md = stages.get("md", {}).get("metrics", {})
    an = stages.get("analysis", {})
    li = dk.get("ligand_info", {})

    def fmt(v, unit=""):
        return f"{v:.3f} {unit}".strip() if isinstance(v, (int, float)) else "N/A"

    lines = []
    lines.append("你是中医药分子计算的科研分析助手。请基于下面给出的【真实数据】写一份中文分析报告。")
    lines.append("")
    lines.append("== 任务对象 ==")
    lines.append(f"配体分子：{ligand_name or smiles}（SMILES: {smiles}）")
    lines.append(f"分子式：{li.get('formula', 'N/A')}；分子量：{li.get('mw', 'N/A')}；LogP：{li.get('logp', 'N/A')}")
    lines.append(f"靶蛋白 PDB：{pdb_id or '用户上传'}")
    lines.append("")
    lines.append("== 蛋白真实身份（来自 PDB 头部，可采信）==")
    lines.append(f"结构分类(HEADER)：{protein_info.get('header') or '未知'}")
    lines.append(f"标题(TITLE)：{protein_info.get('title') or '未知'}")
    lines.append(f"分子名(COMPND)：{protein_info.get('compound') or '未知'}")
    lines.append(f"关键词(KEYWDS)：{protein_info.get('keywords') or '未知'}")
    lines.append(f"EC 编号：{protein_info.get('ec') or '未知'}")
    lines.append(f"来源生物：{protein_info.get('source') or '未知'}")
    if uniprot_info.get("function"):
        lines.append(f"UniProt 功能注释：{uniprot_info['function']}")
    if uniprot_info.get("diseases"):
        lines.append(f"UniProt 关联疾病：{'; '.join(uniprot_info['diseases'])}")
    lines.append("")
    lines.append("== 计算数值结果（真实计算，可采信）==")
    lines.append(f"初始对接最佳亲和力：{fmt(dk.get('best_affinity'), 'kcal/mol')}")
    lines.append(f"xtb 优化后再对接最佳亲和力：{fmt(rd.get('best_affinity'), 'kcal/mol')}")
    lines.append(f"xtb 优化后总能量：{fmt(xt.get('opt_energy_eh'), 'Hartree')}")
    lines.append(f"HOMO-LUMO 能隙：{fmt(pr.get('homo_lumo_gap_ev'), 'eV')}")
    lines.append(f"偶极矩：{fmt(pr.get('dipole_debye'), 'Debye')}")
    lines.append(f"MD 温度：{md.get('temperature_k', 'N/A')} K；步长 {md.get('dt_fs', 'N/A')} fs")
    lines.append(f"MD 步数 NVT/NPT/生产：{md.get('nvt_steps')}/{md.get('npt_steps')}/{md.get('prod_steps')}")
    lines.append(f"RMSD：均值 {fmt(an.get('rmsd_mean_nm'), 'nm')}，最大 {fmt(an.get('rmsd_max_nm'), 'nm')}")
    lines.append(f"Rg（回转半径）：均值 {fmt(an.get('rg_mean_nm'), 'nm')}")
    lines.append(f"RMSF：均值 {fmt(an.get('rmsf_mean_nm'), 'nm')}，最大 {fmt(an.get('rmsf_max_nm'), 'nm')}")
    lines.append("")
    lines.append(_REFERENCE)
    lines.append("")
    lines.append("== 输出要求 ==")
    lines.append("请按以下四个小节写，每节 80~200 字，用平实的科研语言：")
    lines.append("【1 分子性质解读】这个分子的能隙/偶极矩数值预示了怎样的化学稳定性、反应活性与成药性质？")
    lines.append("【2 对接结果解读】亲和力数值代表结合强弱如何？优化前后变化说明什么？")
    lines.append("【3 蛋白-疾病-通路关联】基于上述蛋白真实身份，它属于什么类别、与哪些生理功能/疾病相关？"
                 "请明确区分『已知事实』与『推测』；若信息不足，直接说明需要进一步查证，不要编造具体通路编号或疾病名。"
                 "【严禁重复】同一句话不要反复说，每节只讲一次核心观点。")
    lines.append("【4 分子动力学解读】RMSD/Rg/RMSF 的数值说明了体系稳定性与残基柔性如何？")
    lines.append("最后加一句免责声明：本报告由半经验/经验方法计算并由 AI 辅助解读，结论需实验验证。")
    lines.append("")
    lines.append("【重要纪律】绝对禁止同一句话反复重复。如果发现自己在重复，请立即结束该小节。")
    return "\n".join(lines)


def deduplicate_report(text):
    """
    后处理：检测并移除 LLM 输出中的重复循环。
    小模型偶尔会陷入同一句话反复输出的死循环，这里做兜底清理。
    """
    if not text:
        return text

    # 策略1：按句子切分，检测连续重复的句子（允许轻微变化）
    import re
    # 按句号、问号、感叹号、换行切分成"句子片段"
    sentences = re.split(r'(?<=[。！？\n])', text)
    sentences = [s.strip() for s in sentences if s.strip()]

    if len(sentences) < 3:
        return text

    # 检测重复模式：如果某一句的"核心内容"在之前已经出现过多次，就截断
    seen = {}
    result = []
    for s in sentences:
        # 归一化：去标点、去空白、去英文标点
        key = re.sub(r'[，。、；：""''（）\s\.,;:!?\(\)"\']+', '', s)
        # 如果句子太短，跳过检测
        if len(key) < 8:
            result.append(s)
            continue
        # 检查是否已经出现过 2 次以上
        count = seen.get(key, 0)
        if count >= 2:
            # 已经重复 2 次以上，跳过
            continue
        # 同时检查子串重复——如果当前句子是之前某句的子串且长度>20，也算重复
        is_sub_dup = False
        for prev_key in seen:
            if len(prev_key) > 15 and key in prev_key and prev_key.count(key) >= 2:
                is_sub_dup = True
                break
        if is_sub_dup:
            continue
        seen[key] = count + 1
        result.append(s)

    cleaned = ''.join(result)

    # 策略2：检测整段级别的重复（相同的长文本块反复出现）
    # 按行切分，找重复段落
    lines = cleaned.split('\n')
    deduped_lines = []
    line_signatures = {}
    for line in lines:
        sig = re.sub(r'\s+', '', line)
        if len(sig) > 10:
            if sig in line_signatures:
                continue  # 整行重复，跳过
            line_signatures[sig] = True
        deduped_lines.append(line)
    cleaned = '\n'.join(deduped_lines)

    return cleaned


def call_ollama(prompt, model=None, temperature=None, max_tokens=2200):
    """调用 Ollama /api/generate，返回生成的文本；失败抛出异常。"""
    model = model or OLLAMA_MODEL
    if temperature is None:
        temperature = OLLAMA_TEMPERATURE
    payload = {
        "model": model,
        "prompt": prompt,
        "stream": False,
        "options": {
            "temperature": temperature,
            "num_predict": max_tokens,
        },
    }
    req = urllib.request.Request(
        f"{OLLAMA_URL}/api/generate",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    resp = urllib.request.urlopen(req, timeout=TIMEOUT)
    data = json.loads(resp.read().decode("utf-8"))
    return data.get("response", "").strip()


def check_ollama():
    """检查 Ollama 服务是否可用，返回 (ok, 信息)。"""
    try:
        req = urllib.request.Request(f"{OLLAMA_URL}/api/tags", headers={"User-Agent": "TCMFlow/1.0"})
        data = json.loads(urllib.request.urlopen(req, timeout=5).read())
        models = [m["name"] for m in data.get("models", [])]
        return True, f"已连接，可用模型：{models}"
    except Exception as e:
        return False, str(e)


def generate_report(stages, smiles, ligand_name, pdb_id, protein_pdb_text,
                    model=None):
    """
    主入口：生成分析报告文本。返回 (report_text, meta_dict)。
    meta 里含蛋白身份、uniprot 信息，供调用方一并归档。
    """
    protein_info = parse_pdb_header(protein_pdb_text)
    uniprot_info = fetch_uniprot_summary(pdb_id) if pdb_id else {}
    prompt = build_prompt(stages, smiles, ligand_name, pdb_id,
                          protein_info, uniprot_info)
    text = call_ollama(prompt, model=model)
    # 后处理：清理小模型可能出现的重复循环
    text = deduplicate_report(text)
    return text, {"protein_info": protein_info, "uniprot": uniprot_info,
                  "model": model or OLLAMA_MODEL}

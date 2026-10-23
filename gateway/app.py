# -*- coding: utf-8 -*-
"""
TCM-Gateway · 中医药入口平台
基于TCM-Flow计算引擎，提供中医药数据库与网络药理学衔接功能
TCM-Flow核心代码完全不修改
"""
import os
from flask import (Flask, render_template, request, jsonify, session)

import config
import tcm_data
import tcmflow_client

app = Flask(__name__, static_folder="static", template_folder="templates")
app.secret_key = config.SECRET_KEY
app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024


# ================ 页面 ================

@app.route("/")
def index():
    return render_template("index.html",
                           platform_title=config.PLATFORM_TITLE,
                           platform_subtitle=config.PLATFORM_SUBTITLE)


# ================ 中药数据 API ================

@app.route("/api/herbs")
def api_herbs():
    """中药列表，支持搜索"""
    keyword = request.args.get("q", "").strip()
    limit = min(int(request.args.get("limit", 50)), 200)
    herbs = tcm_data.list_herbs(keyword=keyword or None, limit=limit)
    return jsonify({"success": True, "herbs": herbs})


@app.route("/api/herbs/<name>")
def api_herb_detail(name):
    """中药详情"""
    herb = tcm_data.get_herb_detail(name)
    if not herb:
        return jsonify({"success": False, "error": "未找到该中药"}), 404
    return jsonify({"success": True, "herb": herb})


@app.route("/api/herbs/<name>/compounds")
def api_herb_compounds(name):
    """中药活性成分列表（含SMILES）"""
    compounds = tcm_data.get_herb_compounds(name)
    return jsonify({"success": True, "compounds": compounds})


@app.route("/api/compounds/<name>/targets")
def api_compound_targets(name):
    """成分靶点列表（需要指定中药）"""
    herb = request.args.get("herb", "").strip()
    if not herb:
        return jsonify({"success": False, "error": "请指定中药"}), 400
    targets = tcm_data.get_compound_targets(herb, name)
    return jsonify({"success": True, "targets": targets})


@app.route("/api/compounds/<name>/smiles")
def api_compound_smiles(name):
    """获取成分SMILES"""
    herb = request.args.get("herb", "").strip()
    smiles = tcm_data.get_compound_smiles(name, None)
    if smiles:
        return jsonify({"success": True, "smiles": smiles})
    else:
        return jsonify({"success": False, "error": "未找到该成分的SMILES"}), 404


@app.route("/api/kegg/enrichment")
def api_kegg_enrichment():
    """KEGG通路富集分析"""
    genes_str = request.args.get("genes", "").strip()
    if not genes_str:
        return jsonify({"success": False, "error": "请提供基因列表"}), 400
    genes = [g.strip() for g in genes_str.split(",") if g.strip()]
    top = min(int(request.args.get("top", 15)), 30)
    results = tcm_data.kegg_enrichment(genes, top=top)
    return jsonify({"success": True, "pathways": results})


# ================ TCM-Flow 计算 API ================

@app.route("/api/compute/run", methods=["POST"])
def api_compute_run():
    """提交全流程计算任务"""
    data = request.get_json() or {}
    smiles = (data.get("smiles") or "").strip()
    ligand_name = (data.get("ligand_name") or "").strip()
    pdb_id = (data.get("pdb_id") or "").strip().upper()
    pdb_text = data.get("pdb_text") or ""
    dock_params = data.get("dock_params") or {}
    manual_box = data.get("manual_box") or None
    herb_name = (data.get("herb_name") or "").strip()
    
    if not smiles:
        return jsonify({"success": False, "error": "缺少配体SMILES"}), 400
    
    if not pdb_id and not pdb_text:
        return jsonify({"success": False, "error": "请提供PDB ID或PDB文本"}), 400
    
    # 调用TCM-Flow
    client = tcmflow_client.get_client()
    success, result = client.run_full_pipeline(
        smiles=smiles,
        ligand_name=ligand_name,
        pdb_id=pdb_id or None,
        pdb_text=pdb_text or None,
        dock_params=dock_params,
        manual_box=manual_box
    )
    
    if success:
        session["last_herb"] = herb_name
        session["last_compound"] = ligand_name
        return jsonify({"success": True, "job_id": result})
    else:
        return jsonify({"success": False, "error": result}), 500


@app.route("/api/compute/<job_id>")
def api_compute_status(job_id):
    """查询计算任务状态"""
    client = tcmflow_client.get_client()
    success, result = client.get_job_status(job_id)
    if success:
        return jsonify({"success": True, "job": result})
    else:
        return jsonify({"success": False, "error": result}), 500


@app.route("/api/health")
def api_health():
    """健康检查"""
    tcmflow_ok, tcmflow_info = tcmflow_client.get_client().health_check()
    db_ok = os.path.exists(config.TCM_DB_PATH)
    return jsonify({
        "status": "ok",
        "platform": config.PLATFORM_TITLE,
        "database": {"ok": db_ok, "path": config.TCM_DB_PATH},
        "tcmflow": {"ok": tcmflow_ok, "info": tcmflow_info}
    })


if __name__ == "__main__":
    print("=" * 60)
    print(f"  {config.PLATFORM_TITLE}")
    print(f"  {config.PLATFORM_SUBTITLE}")
    print(f"  本机访问：http://localhost:{config.PORT}")
    print(f"  计算引擎：TCM-Flow (http://localhost:5589)")
    print("=" * 60)
    app.run(host="0.0.0.0", port=config.PORT, debug=False, threaded=True)

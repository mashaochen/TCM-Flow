# -*- coding: utf-8 -*-
"""
TCM-Flow · 全流程分子计算平台 — Flask 后端

独立应用（与 DockLab 分离）：邀请码注册（含有效期）、管理员后台、
全流程计算任务（对接→xtb优化→再对接→MD→性质→分析→报告）、文字报告。
用户界面只出文字报告，原始数据归档在本机，由管理员私下分发。
"""
import os
import json
import uuid
import threading
import datetime
import re as _re

from flask import (Flask, render_template, request, jsonify, session,
                   redirect, url_for, send_file, abort, make_response)

import config
import db
import pipeline

app = Flask(__name__, static_folder="static", template_folder="templates")
app.secret_key = config.SECRET_KEY
app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024
app.config["SEND_FILE_MAX_AGE_DEFAULT"] = 0

db.init_db()
config.ensure_dirs()

_job_lock = threading.Lock()
_running = 0
_pending = []


def current_user():
    u = session.get("username")
    return db.get_user(u) if u else None


def login_required(fn):
    from functools import wraps
    @wraps(fn)
    def w(*a, **kw):
        if not current_user():
            if request.path.startswith("/api/"):
                return jsonify({"success": False, "error": "未登录"}), 401
            return redirect(url_for("login"))
        return fn(*a, **kw)
    return w


def admin_required(fn):
    from functools import wraps
    @wraps(fn)
    def w(*a, **kw):
        u = current_user()
        if not u:
            if request.path.startswith("/api/"):
                return jsonify({"success": False, "error": "未登录"}), 401
            return redirect(url_for("login"))
        if u.get("role") != "admin":
            if request.path.startswith("/api/"):
                return jsonify({"success": False, "error": "无管理员权限"}), 403
            # 普通用户访问后台页面：重定向回首页（避免英文 Forbidden 页面）
            return redirect(url_for("index"))
        return fn(*a, **kw)
    return w


def re_username(s):
    return bool(_re.match(r"^[A-Za-z0-9_]{3,20}$", s))


# ---------------- 页面 ----------------
@app.route("/")
@login_required
def index():
    return render_template("gateway.html",
                           platform_title=config.PLATFORM_TITLE,
                           platform_subtitle=config.PLATFORM_SUBTITLE,
                           username=session.get("username", ""))


@app.route("/compute")
@login_required
def compute_page():
    """经典TCM-Flow计算页面（原首页，保留作为高级功能入口）"""
    return render_template("index.html", profile=config.RESEARCHER_PROFILE,
                           vina_version=config.VINA_VERSION)


@app.route("/about")
def about():
    return render_template("about.html", profile=config.RESEARCHER_PROFILE,
                           title=config.PLATFORM_TITLE, subtitle=config.PLATFORM_SUBTITLE)


@app.route("/login")
def login():
    if current_user():
        return redirect(url_for("index"))
    return render_template("login.html", profile=config.RESEARCHER_PROFILE)


@app.route("/admin")
@admin_required
def admin_page():
    return render_template("admin.html", profile=config.RESEARCHER_PROFILE)


# ---------------- 认证 ----------------
@app.route("/api/auth/register", methods=["POST"])
def api_register():
    data = request.get_json() or {}
    username = (data.get("username") or "").strip()
    password = data.get("password") or ""
    invite_code = (data.get("invite_code") or "").strip().upper()
    if not username or not password:
        return jsonify({"success": False, "error": "用户名和密码不能为空"})
    if len(password) < 6:
        return jsonify({"success": False, "error": "密码至少 6 位"})
    if not re_username(username):
        return jsonify({"success": False, "error": "用户名仅支持字母、数字、下划线（3-20 位）"})
    if db.username_exists(username):
        return jsonify({"success": False, "error": "用户名已存在"})
    ok, err = db.validate_invite_code(invite_code)
    if not ok:
        return jsonify({"success": False, "error": err})
    db.create_user(username, password, invite_code)
    db.mark_invite_code_used(invite_code, username)
    return jsonify({"success": True, "message": "注册成功，请登录"})


@app.route("/api/auth/login", methods=["POST"])
def api_login():
    data = request.get_json() or {}
    username = (data.get("username") or "").strip()
    password = data.get("password") or ""
    user = db.verify_user(username, password)
    if not user:
        return jsonify({"success": False, "error": "用户名或密码错误"})
    session["username"] = user["username"]
    session["role"] = user["role"]
    return jsonify({"success": True, "role": user["role"], "username": user["username"]})


@app.route("/api/auth/logout", methods=["POST"])
def api_logout():
    session.clear()
    return jsonify({"success": True})


@app.route("/api/auth/me")
def api_me():
    u = current_user()
    if not u:
        return jsonify({"logged_in": False})
    return jsonify({"logged_in": True, "username": u["username"], "role": u["role"]})


# ---------------- 管理员 ----------------
@app.route("/api/admin/stats")
@admin_required
def api_admin_stats():
    return jsonify({"success": True, "stats": db.stats()})


@app.route("/api/admin/invites", methods=["GET"])
@admin_required
def api_admin_invites():
    return jsonify({"success": True, "codes": db.list_invite_codes()})


@app.route("/api/admin/invites", methods=["POST"])
@admin_required
def api_admin_create_invites():
    data = request.get_json() or {}
    try:
        count = max(1, min(int(data.get("count", 1)), 50))
    except (ValueError, TypeError):
        count = 1
    note = (data.get("note") or "").strip()
    vd = data.get("validity_days")
    validity_days = None
    if vd not in (None, "", "null"):
        try:
            validity_days = int(vd) if int(vd) in (7, 30, 365) else None
        except (ValueError, TypeError):
            validity_days = None
    codes = db.create_invite_codes(count=count, note=note,
                                   created_by=session.get("username", "admin"),
                                   validity_days=validity_days)
    return jsonify({"success": True, "codes": codes})


@app.route("/api/admin/invites/revoke", methods=["POST"])
@admin_required
def api_admin_revoke_invite():
    data = request.get_json() or {}
    db.revoke_invite_code((data.get("code") or "").strip().upper())
    return jsonify({"success": True})


@app.route("/api/admin/users")
@admin_required
def api_admin_users():
    return jsonify({"success": True, "users": db.list_users()})


@app.route("/api/admin/jobs")
@admin_required
def api_admin_jobs():
    return jsonify({"success": True, "jobs": db.list_jobs(limit=200)})


# ---------------- 全流程任务 ----------------
def _download_pdb(pdb_id):
    import requests
    url = config.PDB_DOWNLOAD_URL.format(pdb_id=pdb_id)
    r = requests.get(url, timeout=60)
    if r.status_code != 200:
        raise RuntimeError(f"RCSB 返回 {r.status_code}")
    text = r.text
    if "ATOM" not in text:
        raise RuntimeError("下载内容不是有效的 PDB")
    return text


def _job_worker(job_id, smiles, ligand_name, pdb_id, protein_pdb_text,
                dock_params, manual_box, output_dir):
    try:
        db.update_job(job_id, status="running",
                      started_at=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        last_progress = {"stage": 0, "name": "准备中", "msg": "任务已启动"}
        def cb(idx, name, msg):
            nonlocal last_progress
            last_progress = {"stage": idx, "name": name, "msg": msg}
            db.update_job(job_id, result_json=json.dumps(
                {"progress": last_progress},
                ensure_ascii=False))
        result = pipeline.run_full_pipeline(
            smiles, protein_pdb_text, pdb_id=pdb_id, ligand_name=ligand_name,
            dock_params=dock_params, manual_box=manual_box,
            output_dir=output_dir, progress_cb=cb)
        if result["success"]:
            db.update_job(
                job_id, status="done",
                completed_at=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                best_affinity=result["stages"].get("redocking", {}).get("best_affinity"),
                result_json=json.dumps({
                    "progress": last_progress,
                    "report_text": result.get("report_text", ""),
                    "stages": {k: v for k, v in result["stages"].items() if k != "report"},
                }, ensure_ascii=False))
        else:
            db.update_job(job_id, status="error", error=result.get("error", "未知错误"),
                          completed_at=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    except Exception as e:
        import traceback
        traceback.print_exc()
        db.update_job(job_id, status="error", error=str(e),
                      completed_at=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    finally:
        global _running
        with _job_lock:
            _running -= 1
            while _pending and _running < config.MAX_CONCURRENT_JOBS:
                args = _pending.pop(0)
                _running += 1
                threading.Thread(target=_job_worker, args=args, daemon=True).start()


@app.route("/api/full/run", methods=["POST"])
@login_required
def api_full_run():
    data = request.get_json() or {}
    smiles = (data.get("smiles") or "").strip()
    ligand_name = (data.get("ligand_name") or "").strip()
    pdb_id = (data.get("pdb_id") or "").strip().upper()
    receptor_source = data.get("receptor_source", "pdb_id")
    pdb_text = data.get("pdb_text") or ""
    dock_params = data.get("dock_params") or {}
    manual_box = data.get("manual_box") or None

    if not smiles:
        return jsonify({"success": False, "error": "请输入配体 SMILES"})

    if receptor_source == "pdb_id":
        if not pdb_id:
            return jsonify({"success": False, "error": "请输入蛋白 PDB ID"})
        try:
            protein_pdb_text = _download_pdb(pdb_id)
        except Exception as e:
            return jsonify({"success": False, "error": f"PDB 下载失败: {e}"})
    else:
        if not pdb_text or len(pdb_text) < 100:
            return jsonify({"success": False, "error": "请上传或粘贴蛋白 PDB 内容"})
        protein_pdb_text = pdb_text

    try:
        from docking import clean_protein_pdb
        clean_protein_pdb(protein_pdb_text)
    except Exception as e:
        return jsonify({"success": False, "error": str(e)})

    params = dict(config.DEFAULT_DOCK_PARAMS)
    for k in ("exhaustiveness", "num_modes", "energy_range", "cpu", "seed", "scoring"):
        if k in dock_params:
            try:
                if k in ("exhaustiveness", "num_modes", "cpu", "seed"):
                    params[k] = int(dock_params[k])
                else:
                    params[k] = float(dock_params[k]) if k == "energy_range" else dock_params[k]
            except (ValueError, TypeError):
                pass
    params["exhaustiveness"] = max(1, min(params["exhaustiveness"], 64))
    params["num_modes"] = max(1, min(params["num_modes"], 20))
    if params["scoring"] not in ("vina", "vinardo", "ad4"):
        params["scoring"] = "vina"

    job_id = db.next_job_id()
    date_str = datetime.datetime.now().strftime("%Y-%m-%d")
    output_dir = os.path.join(config.RESULTS_DIR, date_str, session["username"], job_id)
    os.makedirs(output_dir, exist_ok=True)

    db.create_job(job_id, session["username"], smiles, ligand_name,
                  pdb_id if receptor_source == "pdb_id" else "upload",
                  receptor_source, json.dumps(params, ensure_ascii=False), output_dir)

    args = (job_id, smiles, ligand_name,
            pdb_id if receptor_source == "pdb_id" else "upload",
            protein_pdb_text, params, manual_box, output_dir)

    global _running
    with _job_lock:
        if _running < config.MAX_CONCURRENT_JOBS:
            _running += 1
            threading.Thread(target=_job_worker, args=args, daemon=True).start()
        else:
            _pending.append(args)
            db.update_job(job_id, status="queued")

    return jsonify({"success": True, "job_id": job_id})


@app.route("/api/full/<job_id>")
@login_required
def api_full_status(job_id):
    job = db.get_job(job_id)
    if not job:
        return jsonify({"success": False, "error": "任务不存在"}), 404
    if job["username"] != session["username"] and session.get("role") != "admin":
        return jsonify({"success": False, "error": "无权访问"}), 403
    resp = {
        "id": job["id"], "status": job["status"], "smiles": job["smiles"],
        "ligand_name": job["ligand_name"], "pdb_id": job["pdb_id"],
        "created_at": job["created_at"], "best_affinity": job["best_affinity"],
        "error": pipeline._friendly_error_str(job["error"]),
    }
    if job["result_json"]:
        resp["result"] = json.loads(job["result_json"])
    return jsonify({"success": True, "job": resp})


@app.route("/api/jobs")
@login_required
def api_my_jobs():
    jobs = db.list_jobs(username=session["username"], limit=50)
    return jsonify({"success": True, "jobs": [
        {k: j[k] for k in ("id", "status", "smiles", "ligand_name", "pdb_id",
                           "created_at", "best_affinity")} for j in jobs]})


# ---------------- 管理员数据包 ----------------
_KEY_FILES = ("report.md", "report.txt", "summary.txt",
              # 对接
              "complex.pdb", "complex_opt.pdb", "ligand.pdbqt", "ligand_opt.pdbqt",
              "docked_ligand.pdbqt", "docked_ligand_opt.pdbqt",
              "receptor.pdbqt", "receptor_clean.pdb",
              # xtb
              "xtbopt.xyz",
              # 分子动力学（轨迹 + 拓扑 + 最终结构 + 能量/温度日志，供用户画图/可视化）
              "trajectory.dcd", "system.pdb", "md_final.pdb", "md.log",
              # 轨迹分析时序数据（RMSD/Rg/RMSF，可直接画图）
              "analysis_timeseries.csv")


@app.route("/api/admin/packages")
@admin_required
def api_admin_packages():
    jobs = db.list_jobs(limit=1000)
    tree = {}
    for j in jobs:
        date = (j.get("created_at") or "")[:10] or "未知日期"
        user = j.get("username") or "?"
        files = []
        d = j.get("output_dir")
        if d and os.path.isdir(d):
            for root, _dirs, fnames in os.walk(d):
                for fn in fnames:
                    if fn in _KEY_FILES:
                        files.append(os.path.relpath(os.path.join(root, fn), d))
        entry = {
            "job_id": j["id"], "smiles": j["smiles"], "ligand_name": j["ligand_name"],
            "pdb_id": j["pdb_id"], "status": j["status"],
            "best_affinity": j["best_affinity"], "created_at": j["created_at"],
            "files": sorted(files),
        }
        tree.setdefault(date, {}).setdefault(user, []).append(entry)
    ordered = [{"date": d, "users": [
        {"username": u, "jobs": sorted(js, key=lambda x: x["created_at"] or "", reverse=True)}
        for u, js in sorted(users.items())]} for d, users in sorted(tree.items(), reverse=True)]
    return jsonify({"success": True, "packages": ordered})


@app.route("/api/admin/download/<job_id>/<path:filename>")
@admin_required
def api_admin_download(job_id, filename):
    job = db.get_job(job_id)
    if not job or not job.get("output_dir"):
        return jsonify({"success": False, "error": "任务不存在"}), 404
    base = os.path.realpath(job["output_dir"])
    path = os.path.realpath(os.path.join(base, filename))
    if not (path == base or path.startswith(base + os.sep)) or not os.path.isfile(path):
        return jsonify({"success": False, "error": "文件不存在"}), 404
    return send_file(path, as_attachment=True,
                     download_name=os.path.basename(path))


# ---------------- 中医药数据库 API（只读，来自中医方药智能解析平台） ----------------
@app.route("/api/tcm/herbs")
def tcm_herbs():
    """中药列表，支持搜索"""
    import tcm_data
    keyword = request.args.get("q", "").strip()
    limit = min(int(request.args.get("limit", 50)), 200)
    herbs = tcm_data.list_herbs(keyword=keyword or None, limit=limit)
    return jsonify({"success": True, "herbs": herbs})


@app.route("/api/tcm/herbs/<name>/compounds")
def tcm_herb_compounds(name):
    """中药活性成分列表（含SMILES）"""
    import tcm_data
    compounds = tcm_data.get_herb_compounds(name)
    return jsonify({"success": True, "compounds": compounds})


@app.route("/api/tcm/compounds/<name>/targets")
def tcm_compound_targets(name):
    """成分靶点列表"""
    import tcm_data
    herb = request.args.get("herb", "").strip()
    if not herb:
        return jsonify({"success": False, "error": "请指定中药"}), 400
    targets = tcm_data.get_compound_targets(herb, name)
    return jsonify({"success": True, "targets": targets})


@app.route("/api/tcm/kegg/enrichment")
def tcm_kegg_enrichment():
    """KEGG通路富集分析"""
    import tcm_data
    genes_str = request.args.get("genes", "").strip()
    if not genes_str:
        return jsonify({"success": False, "error": "请提供基因列表"}), 400
    genes = [g.strip() for g in genes_str.split(",") if g.strip()]
    top = min(int(request.args.get("top", 15)), 30)
    results = tcm_data.kegg_enrichment(genes, top=top)
    return jsonify({"success": True, "pathways": results})


@app.route("/api/health")
def health():
    ollama_ok = False
    ollama_msg = "未检测"
    try:
        from agents import llm_report
        ollama_ok, ollama_msg = llm_report.check_ollama()
    except Exception as e:
        ollama_msg = str(e)
    return jsonify({"status": "ok", "platform": config.PLATFORM_TITLE,
                    "vina_exe": os.path.exists(config.VINA_EXE),
                    "xtb_exe": os.path.exists(config.XTB_EXE),
                    "ollama": {"ok": ollama_ok, "info": ollama_msg}})


if __name__ == "__main__":
    print("=" * 62)
    print(f"  {config.PLATFORM_TITLE}")
    print(f"  {config.PLATFORM_SUBTITLE}")
    print(f"  管理员账号：{config.ADMIN_USERNAME}")
    print(f"  本机访问：http://localhost:{config.PORT}")
    print("=" * 62)
    app.run(host="0.0.0.0", port=config.PORT, debug=False, threaded=True)

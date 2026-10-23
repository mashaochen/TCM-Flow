# -*- coding: utf-8 -*-
"""
案例二：柚皮素 + 1CRN 七阶段全流程（通过 API 提交）
"""
import urllib.request
import json
import http.cookiejar
import time
import sys
import os

BASE_URL = "http://localhost:5589"
USERNAME = "admin"
PASSWORD = "admin123456"

# 柚皮素 SMILES
SMILES = "CC(CC1=CC(=C(C=C1)O)O)C(=O)C2=CC=C(C=C2)O"
LIGAND_NAME = "柚皮素 (Naringenin) - 论文案例二"
PDB_ID = "1CRN"


def login(opener):
    data = json.dumps({"username": USERNAME, "password": PASSWORD}).encode('utf-8')
    req = urllib.request.Request(f"{BASE_URL}/api/auth/login",
        data=data, headers={'Content-Type': 'application/json'})
    with opener.open(req) as resp:
        return json.loads(resp.read().decode('utf-8'))


def submit_job(opener):
    data = json.dumps({
        "smiles": SMILES,
        "ligand_name": LIGAND_NAME,
        "receptor_source": "pdb_id",
        "pdb_id": PDB_ID,
        "dock_params": {}
    }).encode('utf-8')
    req = urllib.request.Request(f"{BASE_URL}/api/full/run",
        data=data, headers={'Content-Type': 'application/json'})
    with opener.open(req) as resp:
        return json.loads(resp.read().decode('utf-8'))


def check_status(opener, job_id):
    req = urllib.request.Request(f"{BASE_URL}/api/full/{job_id}")
    with opener.open(req) as resp:
        return json.loads(resp.read().decode('utf-8'))


def main():
    print("=" * 60)
    print("  案例二：柚皮素 + 1CRN 七阶段全流程")
    print("=" * 60)

    cj = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))

    print("\n[1/3] 登录...")
    r = login(opener)
    print(f"  成功: {r.get('success')}, 角色: {r.get('role')}")

    print("\n[2/3] 提交任务...")
    r = submit_job(opener)
    job_id = r.get("job_id")
    print(f"  任务ID: {job_id}")

    print("\n[3/3] 监控进度...")
    last_stage = None
    start = time.time()

    while True:
        r = check_status(opener, job_id)
        job = r.get("job", {})
        status = job.get("status", "unknown")
        result = job.get("result", {})
        prog = result.get("progress", {})
        stage_name = prog.get("name", "等待中")
        stage_num = prog.get("stage", "?")

        elapsed = int(time.time() - start)
        elapsed_str = f"{elapsed//60}分{elapsed%60}秒"

        if stage_name != last_stage:
            print(f"\n  [{elapsed_str}] → 阶段 {stage_num}/7: {stage_name}")
            if prog.get("msg"):
                print(f"         {prog['msg']}")
            last_stage = stage_name

        if status == "done":
            affinity = job.get("best_affinity", "N/A")
            print(f"\n🎉 任务完成！")
            print(f"  最佳亲和力: {affinity} kcal/mol")
            print(f"  总用时: {elapsed_str}")
            print(f"\n  查看结果: {BASE_URL}")
            break
        elif status == "error":
            print(f"\n❌ 任务失败: {job.get('error')}")
            break

        time.sleep(10)


if __name__ == "__main__":
    main()

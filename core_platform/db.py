# -*- coding: utf-8 -*-
"""
DockLab 数据库模块
SQLite 存储：用户、邀请码、对接任务
- 密码使用 PBKDF2 加盐哈希（werkzeug）
- 邀请码由管理员生成，单次使用
"""
import os
import sqlite3
import secrets
import datetime
from werkzeug.security import generate_password_hash, check_password_hash

from config import DB_PATH, ensure_dirs, INVITE_CODE_LENGTH, INVITE_CODE_CHARSET

_SCHEMA = """
CREATE TABLE IF NOT EXISTS users (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    username TEXT UNIQUE NOT NULL,
    password_hash TEXT NOT NULL,
    role TEXT NOT NULL DEFAULT 'user',       -- 'admin' | 'user'
    invite_code TEXT,
    created_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS invite_codes (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    code TEXT UNIQUE NOT NULL,
    note TEXT DEFAULT '',
    status TEXT NOT NULL DEFAULT 'unused',   -- 'unused' | 'used' | 'revoked' | 'expired'
    created_by TEXT DEFAULT 'admin',
    created_at TEXT NOT NULL,
    expires_at TEXT,                         -- 有效期截止时间（NULL = 永久）
    used_by TEXT,
    used_at TEXT
);

CREATE TABLE IF NOT EXISTS jobs (
    id TEXT PRIMARY KEY,                       -- 任务 ID（UUID）
    username TEXT NOT NULL,
    smiles TEXT NOT NULL,
    ligand_name TEXT DEFAULT '',
    pdb_id TEXT DEFAULT '',                    -- 蛋白 PDB ID（或 'upload'）
    receptor_source TEXT DEFAULT '',           -- 'pdb_id' | 'upload'
    status TEXT NOT NULL DEFAULT 'queued',     -- queued|running|done|error
    created_at TEXT NOT NULL,
    started_at TEXT,
    completed_at TEXT,
    best_affinity REAL,
    num_atoms_ligand INTEGER,
    num_residues INTEGER,
    params_json TEXT,                          -- 对接参数（JSON）
    result_json TEXT,                          -- 结果摘要（JSON）
    output_dir TEXT,                           -- 任务输出目录
    error TEXT
);
"""


def _now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def get_conn():
    ensure_dirs()
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    return conn


def init_db():
    """初始化数据库 + 默认管理员账号"""
    conn = get_conn()
    conn.executescript(_SCHEMA)
    # 迁移：旧库补充 expires_at 字段
    cols = [r["name"] for r in conn.execute("PRAGMA table_info(invite_codes)").fetchall()]
    if "expires_at" not in cols:
        conn.execute("ALTER TABLE invite_codes ADD COLUMN expires_at TEXT")
        conn.commit()
    from config import ADMIN_USERNAME, ADMIN_PASSWORD
    row = conn.execute("SELECT id FROM users WHERE username=?", (ADMIN_USERNAME,)).fetchone()
    if not row:
        conn.execute(
            "INSERT INTO users(username, password_hash, role, invite_code, created_at) "
            "VALUES (?,?,?,?,?)",
            (ADMIN_USERNAME, generate_password_hash(ADMIN_PASSWORD),
             "admin", "ADMIN", _now()),
        )
        conn.commit()
        print(f"[TCM-Flow] 已创建管理员账号：{ADMIN_USERNAME}（请尽快修改默认密码）")
    conn.commit()
    conn.close()


# ---------------- 邀请码 ----------------

def generate_invite_code(n=INVITE_CODE_LENGTH):
    return "".join(secrets.choice(INVITE_CODE_CHARSET) for _ in range(n))


def _expires_at(validity_days):
    """根据有效期天数返回截止时间字符串；None=永久。"""
    if validity_days is None:
        return None
    return (datetime.datetime.now() + datetime.timedelta(days=validity_days)) \
        .strftime("%Y-%m-%d %H:%M:%S")


def create_invite_codes(count=1, note="", created_by="admin", validity_days=None):
    """批量生成邀请码（validity_days: None=永久 / 7 / 30 / 365）。返回列表。"""
    conn = get_conn()
    codes = []
    for _ in range(count):
        code = generate_invite_code()
        conn.execute(
            "INSERT INTO invite_codes(code, note, status, created_by, created_at, expires_at) "
            "VALUES (?,?,?,?,?,?)",
            (code, note, "unused", created_by, _now(), _expires_at(validity_days)),
        )
        codes.append(code)
    conn.commit()
    conn.close()
    return codes


def validate_invite_code(code):
    """校验邀请码是否可用。返回 (ok, 错误信息)"""
    conn = get_conn()
    row = conn.execute(
        "SELECT * FROM invite_codes WHERE code=? AND status='unused'",
        (code.strip().upper(),),
    ).fetchone()
    conn.close()
    if not row:
        return False, "邀请码无效或已被使用"
    # 检查有效期
    exp = row["expires_at"]
    if exp:
        now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if exp < now:
            _mark_expired(code.strip().upper())
            return False, "邀请码已过期"
    return True, ""


def _mark_expired(code):
    conn = get_conn()
    conn.execute("UPDATE invite_codes SET status='expired' WHERE code=?", (code,))
    conn.commit()
    conn.close()


def list_invite_codes():
    conn = get_conn()
    rows = conn.execute("SELECT * FROM invite_codes ORDER BY id DESC").fetchall()
    conn.close()
    return [dict(r) for r in rows]


def revoke_invite_code(code):
    conn = get_conn()
    conn.execute("UPDATE invite_codes SET status='revoked' WHERE code=?", (code,))
    conn.commit()
    conn.close()


def mark_invite_code_used(code, username):
    conn = get_conn()
    conn.execute(
        "UPDATE invite_codes SET status='used', used_by=?, used_at=? WHERE code=?",
        (username, _now(), code),
    )
    conn.commit()
    conn.close()


# ---------------- 用户 ----------------

def create_user(username, password, invite_code):
    conn = get_conn()
    conn.execute(
        "INSERT INTO users(username, password_hash, role, invite_code, created_at) "
        "VALUES (?,?,?,?,?)",
        (username, generate_password_hash(password), "user", invite_code, _now()),
    )
    conn.commit()
    conn.close()


def get_user(username):
    conn = get_conn()
    row = conn.execute("SELECT * FROM users WHERE username=?", (username,)).fetchone()
    conn.close()
    return dict(row) if row else None


def verify_user(username, password):
    user = get_user(username)
    if not user:
        return None
    if check_password_hash(user["password_hash"], password):
        return user
    return None


def username_exists(username):
    return get_user(username) is not None


def list_users():
    conn = get_conn()
    rows = conn.execute("SELECT id, username, role, invite_code, created_at "
                        "FROM users ORDER BY id ASC").fetchall()
    conn.close()
    return [dict(r) for r in rows]


def change_password(username, new_password):
    conn = get_conn()
    conn.execute("UPDATE users SET password_hash=? WHERE username=?",
                 (generate_password_hash(new_password), username))
    conn.commit()
    conn.close()


# ---------------- 任务 ----------------

def next_job_id():
    """生成下一个数字任务 ID（从 00001 递增，5 位补零；忽略旧的 uuid 任务）。"""
    conn = get_conn()
    rows = conn.execute("SELECT id FROM jobs").fetchall()
    conn.close()
    max_n = 0
    for r in rows:
        if r["id"].isdigit():
            max_n = max(max_n, int(r["id"]))
    return f"{max_n + 1:05d}"


def create_job(job_id, username, smiles, ligand_name, pdb_id, receptor_source,
               params_json, output_dir):
    conn = get_conn()
    conn.execute(
        "INSERT INTO jobs(id, username, smiles, ligand_name, pdb_id, receptor_source, "
        "status, created_at, params_json, output_dir) VALUES (?,?,?,?,?,?,?,?,?,?)",
        (job_id, username, smiles, ligand_name, pdb_id, receptor_source,
         "queued", _now(), params_json, output_dir),
    )
    conn.commit()
    conn.close()


def update_job(job_id, **fields):
    if not fields:
        return
    conn = get_conn()
    sets = ", ".join(f"{k}=?" for k in fields)
    vals = list(fields.values()) + [job_id]
    conn.execute(f"UPDATE jobs SET {sets} WHERE id=?", vals)
    conn.commit()
    conn.close()


def get_job(job_id):
    conn = get_conn()
    row = conn.execute("SELECT * FROM jobs WHERE id=?", (job_id,)).fetchone()
    conn.close()
    return dict(row) if row else None


def list_jobs(username=None, limit=100):
    conn = get_conn()
    if username:
        rows = conn.execute("SELECT * FROM jobs WHERE username=? ORDER BY created_at DESC LIMIT ?",
                            (username, limit)).fetchall()
    else:
        rows = conn.execute("SELECT * FROM jobs ORDER BY created_at DESC LIMIT ?",
                            (limit,)).fetchall()
    conn.close()
    return [dict(r) for r in rows]


def count_active_jobs():
    conn = get_conn()
    n = conn.execute("SELECT COUNT(*) FROM jobs WHERE status IN ('queued','running')").fetchone()[0]
    conn.close()
    return n


def stats():
    conn = get_conn()
    n_users = conn.execute("SELECT COUNT(*) FROM users WHERE role='user'").fetchone()[0]
    n_jobs = conn.execute("SELECT COUNT(*) FROM jobs").fetchone()[0]
    n_done = conn.execute("SELECT COUNT(*) FROM jobs WHERE status='done'").fetchone()[0]
    n_codes = conn.execute("SELECT COUNT(*) FROM invite_codes").fetchone()[0]
    n_unused = conn.execute("SELECT COUNT(*) FROM invite_codes WHERE status='unused'").fetchone()[0]
    conn.close()
    return {
        "n_users": n_users, "n_jobs": n_jobs, "n_done": n_done,
        "n_codes": n_codes, "n_unused": n_unused,
    }

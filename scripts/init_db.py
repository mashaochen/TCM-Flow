# -*- coding: utf-8 -*-
"""
初始化 TCM-Flow 平台数据库
运行后会在 data 目录下创建空的 tcmflow.db（用户和任务数据库）
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "core_platform"))

import db

print("正在初始化 TCM-Flow 数据库...")
print(f"数据库路径: {db.DB_PATH}")

# 确保表存在
db.init_db()

# 创建管理员账号
try:
    db.create_user("admin", "admin123456", role="admin")
    print("✓ 管理员账号已创建: admin / admin123456")
except Exception as e:
    print(f"  管理员账号已存在: {e}")

# 创建邀请码
try:
    code = db.create_invite_code(plan="pro", duration_days=365, created_by="admin")
    print(f"✓ Pro 邀请码已生成: {code}")
except Exception as e:
    print(f"  邀请码生成跳过: {e}")

print("\n数据库初始化完成！")
print(f"启动平台: cd core_platform && python app.py")

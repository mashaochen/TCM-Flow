# -*- coding: utf-8 -*-
"""
TCM-Gateway · 中医药入口平台
配置文件
"""
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# ===== 平台信息 =====
PLATFORM_TITLE = "中医方药智能计算平台"
PLATFORM_SUBTITLE = "基于TCM-Flow的中医药分子计算入口"
SECRET_KEY = "tcm-gateway-secret-key-change-in-production"
PORT = 5000

# ===== TCM-Flow 配置（计算引擎，完全不修改）=====
TCMFLOW_BASE_URL = "http://localhost:5589"
TCMFLOW_ADMIN_USER = "admin"
TCMFLOW_ADMIN_PASS = "admin123456"

# ===== 中药数据库（只读，来自中医方药智能解析平台）=====
TCM_DB_PATH = os.environ.get("TCM_GATEWAY_DB_PATH", os.path.join(BASE_DIR, "..", "tcm-platform", "backend", "tcm.db"))

# ===== 结果缓存 =====
CACHE_DIR = os.path.join(BASE_DIR, "cache")
os.makedirs(CACHE_DIR, exist_ok=True)

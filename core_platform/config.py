# -*- coding: utf-8 -*-
"""
TCM-Flow 配置模块
全流程分子计算平台（对接 → xtb 优化 → 再对接 → OpenMM MD → 性质 → 分析 → 报告）
"""
import os
import sys

# 项目根目录
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# 关键可执行文件路径（支持环境变量覆盖，便于第三方环境复现）
VINA_EXE = os.environ.get("VINA_EXE", os.path.join(BASE_DIR, "bin", "vina.exe"))
XTB_EXE = os.environ.get("XTB_EXE", os.path.join(BASE_DIR, "bin", "xtb-6.7.1", "bin", "xtb.exe"))

# meeko 受体/配体准备 CLI
def _scripts_dir():
    d = os.path.join(os.path.dirname(sys.executable), "Scripts")
    if os.path.isdir(d):
        return d
    return os.path.join(os.path.expanduser("~"),
                        ".workbuddy", "binaries", "python", "envs", "default", "Scripts")

SCRIPTS_DIR = _scripts_dir()
MK_PREPARE_RECEPTOR = os.path.join(SCRIPTS_DIR, "mk_prepare_receptor.exe")

# 数据目录（支持环境变量覆盖，便于部署到不同环境）
DATA_DIR = os.environ.get("TCMFLOW_DATA_DIR", os.path.join(BASE_DIR, "data"))
DB_PATH = os.environ.get("TCMFLOW_DB_PATH", os.path.join(DATA_DIR, "tcmflow.db"))
RESULTS_DIR = os.environ.get("TCMFLOW_RESULTS_DIR", os.path.join(DATA_DIR, "results"))

# 中医药数据库（只读，来自中医方药智能解析平台）
TCM_DB_PATH = os.environ.get("TCM_DB_PATH",
    os.path.join(BASE_DIR, "..", "..", "tcm-platform", "backend", "tcm.db"))

# 平台信息
PLATFORM_NAME = "TCM-Flow"
PLATFORM_TITLE = "TCM-Flow · 全流程分子计算平台"
PLATFORM_SUBTITLE = "对接 · xtb 优化 · 分子动力学 · 性质 · 报告 一体化流水线"
VINA_VERSION = "1.2.7"
AFFILIATION = "南京中医药大学 · 人工智能与信息技术学院"

# 会话/安全
SECRET_KEY = os.environ.get("TCMFLOW_SECRET_KEY",
    "tcmflow-li-xiong-group-2026-secret-key-please-change")

# 管理员账号
ADMIN_USERNAME = os.environ.get("TCMFLOW_ADMIN_USER", "admin")
ADMIN_PASSWORD = os.environ.get("TCMFLOW_ADMIN_PASSWORD", "admin123456")

# 服务端口
PORT = int(os.environ.get("PORT", 5589))

# 对接默认参数
DEFAULT_DOCK_PARAMS = {
    "exhaustiveness": 8,
    "num_modes": 9,
    "energy_range": 3,
    "cpu": 0,
    "seed": 42,
    "scoring": "vina",
    "box_padding": 8.0,
}

# 最大并发任务数
MAX_CONCURRENT_JOBS = 1   # 全流程任务独占资源，一次一个

# 单任务超时（秒）
JOB_TIMEOUT = 3600

# 邀请码
INVITE_CODE_LENGTH = 8
INVITE_CODE_CHARSET = "ABCDEFGHJKLMNPQRSTUVWXYZ23456789"

# PDB 下载源
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# 测试模式开关：True=4核/0.5GB（本机验证），False=24核/1GB（生产）
TEST_MODE = os.environ.get("TCMFLOW_TEST_MODE", "1") == "1"

RES = {
    "test": {"ncores": 4, "mem_mb": 500},
    "prod": {"ncores": 24, "mem_mb": 1024},
}

XTB_GFN = 2
XTB_SOLVENT = None

# 本地 LLM（Ollama）报告配置
OLLAMA_URL = os.environ.get("OLLAMA_URL", "http://localhost:11434")
OLLAMA_MODEL = os.environ.get("OLLAMA_MODEL", "qwen2.5:1.5b-instruct")
OLLAMA_TIMEOUT = int(os.environ.get("OLLAMA_TIMEOUT", "600"))
OLLAMA_TEMPERATURE = float(os.environ.get("OLLAMA_TEMPERATURE", "0.3"))
# 是否启用 LLM 报告（False 时 pipeline 直接走模板兜底）
LLM_REPORT_ENABLED = os.environ.get("TCMFLOW_LLM_ENABLED", "1") == "1"

MD_PARAMS = {
    "test": {
        "temperature": 300, "dt_fs": 2.0,
        "nvt_steps": 500, "npt_steps": 500, "prod_steps": 1000,
        "report_interval": 100,
        "padding_nm": 1.0, "ionic_strength": 0.15,
    },
    "prod": {
        "temperature": 300, "dt_fs": 2.0,
        "nvt_steps": 50000, "npt_steps": 50000, "prod_steps": 500000,
        "report_interval": 1000,
        "padding_nm": 1.0, "ionic_strength": 0.15,
    },
}


def get_res():
    return RES["test" if TEST_MODE else "prod"]


def get_md_params():
    return MD_PARAMS["test" if TEST_MODE else "prod"]


# 李雄个人简介（平台"关于"页/主页展示）
RESEARCHER_PROFILE = {
    "name": "李雄",
    "title": "博士",
    "affiliation": "南京中医药大学 人工智能与信息技术学院",
    "degree": "东南大学 物理学博士（2018）",
    "email": "xli@njucm.edu.cn",
    "publications": "已发表 23 篇 SCI 第一作者论文，持 2 项发明专利，4项软件著作权。",
    "research_areas": [
        "中医四诊（望闻问切）客观化、数字化与现代化",
        "计算物理与量子化学（GFN2-xTB、DFT）",
        "网络医学与 DDI/PPI 网络、靶点数据库",
        "分子对接与结构生物学辅助的药物发现",
        "新型物理传感器与中医客观化结合",
    ],
    "summary": (
        "李雄博士长期从事计算物理、量子化学与中医药现代化的交叉研究，"
        "将分子模拟、网络医学与临床数据相结合，推动中医四诊客观化与中药活性成分的机制阐释。"
        "本平台（TCM-Flow）提供蛋白质-配体全流程分子计算服务：分子对接、半经验优化、"
        "分子动力学模拟、性质提取与轨迹分析，服务于中医药活性成分的靶点发现与机制研究。"
    ),
}


def ensure_dirs():
    for d in (DATA_DIR, RESULTS_DIR):
        os.makedirs(d, exist_ok=True)

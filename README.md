# TCM-Flow · 中医药全流程分子计算平台

> **T**raditional **C**hinese **M**edicine **Flow** —— 面向中医药研究的中药活性成分全流程分子计算平台

TCM-Flow 将**中药活性成分检索 → 半经验量子化学优化 → 分子对接 → 分子动力学模拟 → 性质分析 → 大语言模型智能报告**整合为七阶段自动化管线，研究者仅需输入中药名称与靶点蛋白，即可获得完整的计算分析报告。

> 📖 **详细安装指南请见 [INSTALL.md](INSTALL.md)**

---

## ✨ 特色功能

- 🌿 **中医药专属**：内置 3,800+ 味中药、7,200+ 种活性成分数据库，支持中药→成分→靶点→通路一站式检索
- 🚀 **七阶段自动化管线**：对接 → xtb 优化 → 再对接 → OpenMM MD → 性质 → 分析 → AI 报告，一键提交全自动完成
- 🛡️ **复合物 MD 稳定方案**：HBonds 约束 + 强摩擦 Langevin (γ=10/ps) + 0.5 fs 步长 + 5 阶段逐步升温 + 自动重试
- 🧠 **本地 LLM 智能报告**：基于 Ollama 的本地大模型，三层反幻觉约束（事实接地 + 数值锚定 + 推测标注），数据不出本地
- 🔐 **邀请制三级用户系统**：免费版 / Pro 版 / 企业版，管理员后台统一管控
- 📦 **开箱即用**：内置 Vina 可执行文件，Python 依赖 pip 一键安装

---

## 🏗️ 系统架构

```
┌─────────────────────────────────────────────────────────┐
│                    中医药入口网关 (可选)                  │
│              中药搜索 · 成分筛选 · 靶点通路              │
└────────────────────────────┬────────────────────────────┘
                             │ API 调用
                             ▼
┌─────────────────────────────────────────────────────────┐
│                 TCM-Flow 核心计算平台                    │
│  ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐ │
│  │ 对接  │ │ xtb  │ │再对接│ │  MD  │ │性质  │ │分析  │ │
│  └──┬───┘ └──┬───┘ └──┬───┘ └──┬───┘ └──┬───┘ └──┬───┘ │
│     │        │        │        │        │        │      │
│  ┌──▼────────▼────────▼────────▼────────▼────────▼───┐  │
│  │              七阶段管线调度引擎                    │  │
│  └──────────────────────┬───────────────────────────┘  │
│                         │                              │
│  ┌──────────────────────▼───────────────────────────┐  │
│  │           LLM 智能报告生成 (Ollama)               │  │
│  └──────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────┘
```

---

## 🚀 快速开始（5 分钟跑通案例一）

### 前置条件

- Python 3.9 或 3.10
- Windows 10/11 x64
- 内存 ≥ 4 GB

### 三步跑通虚拟筛选

```bash
# 1. 克隆仓库
git clone https://github.com/mashaochen/TCM-Flow.git
cd TCM-Flow

# 2. 安装依赖
pip install -r requirements.txt

# 3. 运行案例一（7种黄酮对 2J0L 虚拟筛选）
cd examples/case1_virtual_screening
python run_case1.py
```

> Vina 已内置在 `core_platform/bin/vina.exe`，无需额外安装。
> 如需完整七流程计算，还需安装 xtb 和 OpenMM，详见 [INSTALL.md](INSTALL.md)。

### 启动 Web 平台

```bash
cd core_platform
python app.py
```

访问 http://localhost:5589 ，使用以下账号登录：
- 管理员：`admin` / `admin123456`

---

## 📦 依赖清单

### Python 包（必需）

| 包名 | 版本要求 | 用途 |
|------|---------|------|
| flask | ≥ 2.0 | Web 服务框架 |
| requests | ≥ 2.25 | HTTP 请求 |
| rdkit | ≥ 2021.09 | 分子处理（SMILES→3D、性质计算等） |
| meeko | ≥ 0.4.0 | 配体/受体 PDBQT 准备 |
| openmm | ≥ 7.7 | 分子动力学模拟 |
| pdbfixer | ≥ 1.8 | PDB 修复 |
| numpy | ≥ 1.20 | 数值计算 |

### 外部计算引擎

| 引擎 | 版本 | 用途 | 是否内置 | 大小 |
|------|------|------|---------|------|
| AutoDock Vina | 1.2.7 | 分子对接 | ✅ 内置 | 1.2 MB |
| xtb | 6.7.1 | 量子化学优化 | ❌ 需下载 | ~123 MB |
| OpenMM | ≥ 7.7 | 分子动力学 | ❌ pip/conda 安装 | - |
| Ollama | 最新 | AI 报告 | ❌ 可选 | ~40 MB + 模型 |

### 可选依赖

| 依赖 | 用途 | 不安装的影响 |
|------|------|-------------|
| mdtraj | 轨迹分析加速 | 自动降级为纯 NumPy 实现 |
| Ollama + qwen2.5 模型 | AI 报告生成 | 自动降级为模板报告 |
| 中药数据库 tcm.db | 中医药网关页面 | 网关功能不可用，核心计算不受影响（[Zenodo 下载](https://doi.org/10.5281/zenodo.23600814)） |

> **完整安装步骤见 [INSTALL.md](INSTALL.md)**

---

## 🧪 论文案例复现

本仓库包含论文中两个案例的完整复现代码：

### 案例一：7 种黄酮对 FAK (PDB: 2J0L) 虚拟筛选

```bash
cd examples/case1_virtual_screening
python run_case1.py
```

- **所需时间**：约 2-5 分钟
- **依赖**：仅需 rdkit + meeko + Vina（已内置）
- **输出**：7 个成分的对接打分排序表

### 案例二：柚皮素 + Crambin (PDB: 1CRN) 七阶段全流程

```bash
# 先启动平台
cd core_platform
python app.py

# 另开终端，提交任务
cd examples/case2_full_pipeline
python run_case2_api.py
```

- **所需时间**：测试档约 5-10 分钟，生产档约 1-2 小时
- **依赖**：全部计算引擎（Vina + xtb + OpenMM）
- **输出**：七阶段完整结果（对接报告、xtb 能量、MD 轨迹、性质数据、AI 报告）

> 每个案例目录下有 `README.md` 详细说明运行方法和预期结果。

---

## 📁 目录结构

```
TCM-Flow/
├── README.md                     # 本文件
├── INSTALL.md                    # 详细安装指南
├── QUICKSTART.md                 # 快速入门（Windows）
├── LICENSE                       # MIT 开源协议
├── requirements.txt              # Python 依赖
├── .gitignore                    # Git 忽略规则
│
├── core_platform/                # 核心计算平台（主程序）
│   ├── app.py                    # Flask 主应用
│   ├── config.py                 # 配置文件（支持环境变量覆盖）
│   ├── db.py                     # 数据库操作
│   ├── docking.py                # 分子对接模块
│   ├── pipeline.py               # 七阶段管线调度
│   ├── tcm_data.py               # 中药数据访问
│   ├── compound_smiles.py        # 化合物 SMILES 库
│   ├── ollama_server.py          # Ollama LLM 服务
│   ├── bin/                      # 外部计算引擎
│   │   ├── vina.exe              # AutoDock Vina（已内置）
│   │   └── README.md             # 其他引擎安装说明
│   ├── agents/                   # 计算代理
│   │   ├── llm_report.py         # LLM 报告生成
│   │   ├── xtb_opt.py            # xtb 几何优化
│   │   ├── openmm_md.py          # OpenMM 分子动力学
│   │   ├── analysis.py           # 轨迹分析
│   │   └── ligand_param.py       # 配体参数化
│   ├── templates/                # HTML 模板
│   │   ├── gateway.html          # 中医药首页
│   │   ├── index.html            # 高级计算页
│   │   ├── login.html            # 登录/注册
│   │   ├── admin.html            # 管理后台
│   │   └── about.html            # 关于页
│   └── static/                   # 静态资源（CSS + JS）
│
├── gateway/                      # 独立中医药入口网关（可选部署）
│   ├── app.py
│   ├── tcmflow_client.py
│   ├── templates/
│   └── static/
│
├── examples/                     # 论文案例复现代码
│   ├── case1_virtual_screening/  # 案例一：7种黄酮虚拟筛选
│   │   ├── run_case1.py
│   │   ├── ligands.smi
│   │   └── README.md
│   └── case2_full_pipeline/      # 案例二：七阶段全流程
│       ├── run_case2_api.py
│       └── README.md
│
├── data/                         # 数据目录
│   └── README.md                 # 数据库获取说明
│
└── scripts/                      # 工具脚本
    └── init_db.py                # 初始化数据库
```

---

## 🔧 环境变量配置

所有关键路径均支持环境变量覆盖，便于在不同环境部署：

| 环境变量 | 说明 | 默认值 |
|---------|------|--------|
| `TCMFLOW_DATA_DIR` | 数据目录 | `core_platform/data` |
| `TCMFLOW_DB_PATH` | 平台数据库路径 | `$DATA_DIR/tcmflow.db` |
| `TCM_DB_PATH` | 中药数据库路径 | - |
| `TCMFLOW_TEST_MODE` | 测试模式（`1`=测试档，`0`=生产档） | `1` |
| `VINA_EXE` | Vina 可执行文件路径 | `core_platform/bin/vina.exe` |
| `XTB_EXE` | xtb 可执行文件路径 | `core_platform/bin/xtb-6.7.1/bin/xtb.exe` |
| `OLLAMA_URL` | Ollama 服务地址 | `http://localhost:11434` |
| `OLLAMA_MODEL` | LLM 模型名 | `qwen2.5:1.5b-instruct` |
| `PORT` | 服务端口 | `5589` |

---

## 📊 两档配置

平台提供两档 MD 模拟配置，通过环境变量 `TCMFLOW_TEST_MODE` 切换：

| 参数 | 测试档（默认） | 生产档 |
|------|---------------|--------|
| NVT 平衡 | 500 步 | 50,000 步 |
| NPT 平衡 | 500 步 | 50,000 步 |
| 生产 MD | 1,000 步 / 0.5 ps | 500,000 步 / 250 ps |
| 轨迹保存间隔 | 100 步 | 1,000 步 |
| 适用场景 | 演示、教学、快速验证 | 正式科研计算 |

> **切换到生产档**：`set TCMFLOW_TEST_MODE=0` 后重新启动平台

---

## 📝 论文引用

如果您在研究中使用了 TCM-Flow，请引用我们的论文：

```
[待补充] 作者. 标题. 期刊, 年份.
```

数据集已在 Zenodo 公开：[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23600814.svg)](https://doi.org/10.5281/zenodo.23600814)

---

## 📄 许可证

MIT License - 详见 [LICENSE](LICENSE)

---

## 🤝 联系方式

如有问题或建议，欢迎提交 Issue 或联系作者。

- 南京中医药大学 人工智能与信息技术学院
- 邮箱：xli@njucm.edu.cn

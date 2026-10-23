# TCM-Flow · 中医药全流程分子计算平台

> **T**raditional **C**hinese **M**edicine **Flow** —— 面向中医药研究的中药活性成分全流程分子计算平台

TCM-Flow 将**中药活性成分检索 → 半经验量子化学优化 → 分子对接 → 分子动力学模拟 → 性质分析 → 大语言模型智能报告**整合为七阶段自动化管线，研究者仅需输入中药名称与靶点蛋白，即可获得完整的计算分析报告。

---

## ✨ 特色功能

- 🌿 **中医药专属**：内置 3,800+ 味中药、7,200+ 种活性成分数据库，支持中药→成分→靶点→通路一站式检索
- 🚀 **七阶段自动化管线**：对接 → xtb 优化 → 再对接 → OpenMM MD → 性质 → 分析 → AI 报告，一键提交全自动完成
- 🛡️ **复合物 MD 稳定方案**：HBonds 约束 + 强摩擦 Langevin (γ=10/ps) + 0.5 fs 步长 + 5 阶段逐步升温 + 自动重试
- 🧠 **本地 LLM 智能报告**：基于 Ollama 的本地大模型，三层反幻觉约束（事实接地 + 数值锚定 + 推测标注），数据不出本地
- 🔐 **邀请制三级用户系统**：免费版 / Pro 版 / 企业版，管理员后台统一管控
- 📦 **零配置部署**：内置 xtb、Vina 等计算引擎二进制，用户无需安装计算软件

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

## 🚀 快速开始

### 环境要求

- Python 3.8 ~ 3.11
- Windows / Linux / macOS
- 内存 ≥ 4 GB（生产档 MD 建议 ≥ 16 GB）

### 1. 克隆仓库

```bash
git clone https://github.com/yourname/TCM-Flow.git
cd TCM-Flow
```

### 2. 安装依赖

```bash
pip install -r requirements.txt
```

### 3. 初始化数据库

```bash
python scripts/init_db.py
```

### 4. 下载中药数据库

完整数据库（约 100 MB）请从论文补充材料提供的链接下载，
解压后设置环境变量：

```bash
# Windows
set TCM_DB_PATH=C:\path\to\tcm.db

# Linux / macOS
export TCM_DB_PATH=/path/to/tcm.db
```

### 5. 安装计算引擎

平台需要以下计算引擎支持：

| 引擎 | 用途 | 下载地址 |
|------|------|---------|
| AutoDock Vina | 分子对接 | https://vina.scripps.edu/ |
| xtb | 量子化学优化 | https://github.com/grimme-lab/xtb |
| OpenMM | 分子动力学 | `pip install openmm pdbfixer` |
| Ollama（可选） | AI 报告 | https://ollama.ai/ |

将 Vina 和 xtb 的可执行文件路径添加到系统 PATH，
或放入 `core_platform/bin/` 目录。

### 6. 启动平台

```bash
cd core_platform
python app.py
```

访问 http://localhost:5589 ，使用以下账号登录：
- 管理员：`admin` / `admin123456`

---

## 📁 目录结构

```
TCM-Flow/
├── README.md                     # 本文件
├── requirements.txt              # Python 依赖
├── .gitignore                    # Git 忽略规则
│
├── core_platform/                # 核心计算平台（主程序）
│   ├── app.py                    # Flask 主应用
│   ├── config.py                 # 配置文件
│   ├── db.py                     # 数据库操作
│   ├── docking.py                # 分子对接模块
│   ├── pipeline.py               # 七阶段管线调度
│   ├── tcm_data.py               # 中药数据访问
│   ├── agents/                   # 计算代理
│   │   ├── llm_report.py         # LLM 报告生成
│   │   ├── xtb_opt.py            # xtb 几何优化
│   │   ├── openmm_md.py          # OpenMM 分子动力学
│   │   ├── analysis.py           # 轨迹分析
│   │   └── ligand_param.py       # 配体参数化
│   ├── templates/                # HTML 模板
│   └── static/                   # 静态资源
│
├── gateway/                      # 独立中医药入口网关
│   ├── app.py                    # Flask 应用
│   ├── tcmflow_client.py         # TCM-Flow API 客户端
│   ├── templates/                # 前端页面
│   └── static/                   # 静态资源
│
├── examples/                     # 示例
│   ├── case1_virtual_screening/  # 案例一：7种黄酮虚拟筛选
│   └── case2_full_pipeline/      # 案例二：柚皮素七阶段全流程
│
├── data/                         # 数据目录
│   └── README.md                 # 数据获取说明
│
└── scripts/                      # 工具脚本
    └── init_db.py                # 初始化数据库
```

---

## 🧪 使用示例

### 案例一：虚拟筛选

```bash
cd examples/case1_virtual_screening
python run_case1.py
```

### 案例二：七阶段全流程

```bash
# 先启动平台
cd core_platform
python app.py

# 另开终端，提交任务
cd examples/case2_full_pipeline
python run_case2_api.py
```

更多示例请查看 [`examples/`](examples/) 目录。

---

## 📊 两档配置

平台提供两档 MD 模拟配置，可在 `config.py` 中通过 `TEST_MODE` 切换：

| 参数 | 测试档（默认） | 生产档 |
|------|---------------|--------|
| NVT 平衡 | 500 步 | 50,000 步 |
| NPT 平衡 | 500 步 | 50,000 步 |
| 生产 MD | 1,000 步 / 0.5 ps | 500,000 步 / 250 ps |
| 轨迹保存间隔 | 100 步 | 1,000 步 |
| 适用场景 | 演示、教学、快速验证 | 正式科研计算 |

---

## 📝 论文引用

如果您在研究中使用了 TCM-Flow，请引用我们的论文：

```
[待补充] 作者. 标题. 期刊, 年份.
```

---

## 📄 许可证

MIT License

---

## 🤝 联系方式

如有问题或建议，欢迎提交 Issue 或联系作者。

- 南京中医药大学 人工智能与信息技术学院
- 邮箱：xli@njucm.edu.cn

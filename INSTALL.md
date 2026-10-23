# 安装指南

本文档详细说明如何从零开始安装 TCM-Flow 的所有依赖，确保可复现。

---

## 系统要求

| 项目 | 最低要求 | 推荐配置 |
|------|---------|---------|
| 操作系统 | Windows 10 x64 | Windows 10/11 x64 |
| Python | 3.9 | 3.10 |
| 内存 | 4 GB | 8 GB+ |
| 磁盘空间 | 2 GB | 5 GB+ |
| CPU | 4 核 | 8 核+ |

> **注意**：本平台目前仅支持 Windows。Linux/macOS 用户可通过修改可执行文件路径适配。

---

## 第一步：安装 Python

### 方式 A：Miniconda（推荐，OpenMM 安装更方便）

1. 下载 Miniconda：https://docs.conda.io/en/latest/miniconda.html
2. 安装时勾选 "Add Miniconda3 to my PATH environment variable"
3. 打开命令行，创建环境：
   ```bash
   conda create -n tcmflow python=3.10 -y
   conda activate tcmflow
   ```

### 方式 B：原生 Python

1. 下载 Python 3.10：https://www.python.org/downloads/release/python-31011/
2. 安装时勾选 "Add Python to PATH"

---

## 第二步：安装 Python 依赖

```bash
# 进入项目目录
cd TCM-Flow

# 安装核心依赖
pip install -r requirements.txt
```

### 依赖说明

| 包名 | 用途 | 必需 |
|------|------|------|
| flask | Web 服务框架 | ✅ |
| requests | HTTP 请求 | ✅ |
| rdkit | 分子处理（SMILES→3D、性质计算等） | ✅ |
| meeko | 配体/受体 PDBQT 准备 | ✅ |
| numpy | 数值计算 | ✅（rdkit 依赖） |
| openmm | 分子动力学模拟 | ✅ |
| pdbfixer | PDB 修复 | ✅ |
| mdtraj | 轨迹分析加速 | ❌（可选，自动降级） |
| ollama | LLM 调用 | ❌（可选，API 模式不需要） |

### Conda 环境安装 OpenMM（推荐）

```bash
conda install -c conda-forge openmm pdbfixer rdkit -y
pip install flask requests meeko
```

---

## 第三步：安装外部计算引擎

### 3.1 AutoDock Vina（已内置 ✅）

Vina 可执行文件已内置在 `core_platform/bin/vina.exe`，无需额外安装。

版本：AutoDock Vina v1.2.7（Apache 2.0 协议）

如需使用自己的版本：
```bash
set VINA_EXE=C:\path\to\vina.exe
```

### 3.2 xtb（必需）

xtb 用于第二阶段的量子化学几何优化。

**下载安装：**
1. 访问：https://github.com/grimme-lab/xtb/releases
2. 下载 `xtb-6.7.1-windows-x86_64.zip`
3. 解压到 `core_platform/bin/xtb-6.7.1/`
4. 确保 `core_platform/bin/xtb-6.7.1/bin/xtb.exe` 存在

**或设置环境变量：**
```bash
set XTB_EXE=C:\path\to\xtb.exe
```

**验证安装：**
```bash
xtb --version
```

### 3.3 OpenMM + PDBFixer（必需）

用于第四阶段分子动力学模拟。

**Conda 安装（推荐）：**
```bash
conda install -c conda-forge openmm pdbfixer -y
```

**Pip 安装（仅 Linux）：**
```bash
# Windows 上 pip 安装 openmm 可能失败，建议用 conda
pip install openmm pdbfixer
```

**验证安装：**
```bash
python -c "import openmm; print('OpenMM version:', openmm.__version__)"
python -c "import pdbfixer; print('PDBFixer OK')"
```

### 3.4 Ollama（可选）

用于第七阶段 AI 报告生成。不安装也能跑，报告会自动降级为模板模式。

1. 下载安装：https://ollama.com/
2. 拉取模型：
   ```bash
   ollama pull qwen2.5:1.5b-instruct
   ```
3. 验证：
   ```bash
   ollama list
   ```

---

## 第四步：准备数据

### 4.1 中药数据库（可选，用于中医药网关页面）

完整数据库（约 100 MB）包含 3,819 味中药、7,284 种活性成分。

**获取方式：**
- 从论文补充材料下载
- 或从 GitHub Release 下载
- 下载后放入 `data/tcm.db`

**或设置环境变量：**
```bash
set TCM_DB_PATH=C:\path\to\tcm.db
```

### 4.2 平台数据库（自动创建）

`tcmflow.db`（用户、任务、邀请码等）会在首次启动时自动创建，无需手动准备。

默认位置：`data/tcmflow.db`

---

## 第五步：启动平台

### Windows 一键启动

双击 `start_windows.bat`

或手动启动：
```bash
cd core_platform
python app.py
```

### 访问地址

- 主页面：http://localhost:5589
- 管理员账号：`admin` / `admin123456`
  > ⚠️ 首次登录后请立即修改默认密码！

### 环境变量配置（可选）

| 变量名 | 说明 | 默认值 |
|--------|------|--------|
| `TCMFLOW_DATA_DIR` | 数据目录 | `core_platform/data` |
| `TCMFLOW_DB_PATH` | 平台数据库路径 | `$DATA_DIR/tcmflow.db` |
| `TCM_DB_PATH` | 中药数据库路径 | `../tcm-platform/backend/tcm.db` |
| `TCMFLOW_TEST_MODE` | 测试模式（1=测试档，0=生产档） | `1` |
| `VINA_EXE` | Vina 可执行文件路径 | `bin/vina.exe` |
| `XTB_EXE` | xtb 可执行文件路径 | `bin/xtb-6.7.1/bin/xtb.exe` |
| `OLLAMA_URL` | Ollama 服务地址 | `http://localhost:11434` |
| `OLLAMA_MODEL` | LLM 模型名 | `qwen2.5:1.5b-instruct` |
| `PORT` | 服务端口 | `5589` |

---

## 常见问题

### Q: `ImportError: DLL load failed while importing rdkit`
A: 建议使用 conda 安装 rdkit：`conda install -c conda-forge rdkit`

### Q: xtb 运行报错缺少 DLL
A: 确保 xtb 目录下的 `bin/libiomp5md.dll` 存在，或将其加入系统 PATH

### Q: OpenMM 找不到平台
A: Windows 上建议用 conda 安装 openmm，pip 安装可能缺少 CUDA/OpenCL 平台

### Q: 对接报错 "vina.exe not found"
A: 检查 `core_platform/bin/vina.exe` 是否存在，或设置 `VINA_EXE` 环境变量

### Q: 中药数据加载失败
A: 检查 `TCM_DB_PATH` 是否指向正确的 tcm.db 文件路径

### Q: 测试模式和生产模式的区别
A: 测试模式 MD 只跑 1000 步（约几秒），用于验证流程；生产模式跑 500,000 步（约几十分钟到几小时），用于正式计算。

---

## 验证安装是否成功

运行以下命令，依次检查每个阶段：

```bash
cd core_platform

# 检查对接模块
python -c "import docking; print('✅ Docking OK')"

# 检查 xtb 模块
python -c "from agents.xtb_opt import optimize_geometry; print('✅ xTB OK')"

# 检查 OpenMM 模块
python -c "from agents.openmm_md import run_md; print('✅ OpenMM OK')"

# 检查分析模块
python -c "from agents.analysis import analyze_trajectory; print('✅ Analysis OK')"

# 启动平台
python app.py
```

如果全部打印 ✅，说明安装成功！

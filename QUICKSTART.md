# 快速入门

## Windows 用户一键启动

### 第一步：安装 Python
- 下载 Python 3.9 或 3.10：https://www.python.org/downloads/
- 安装时勾选 **Add Python to PATH**

### 第二步：安装依赖
双击运行 `install_windows.bat`

### 第三步：准备数据
1. 下载中药数据库 `tcm.db`（从论文补充材料或 GitHub Release）
2. 放入 `data/` 目录，或设置环境变量 `TCMFLOW_DB_PATH`

### 第四步：安装计算引擎

| 引擎 | 用途 | 安装方法 |
|------|------|---------|
| AutoDock Vina | 分子对接 | 下载 vina.exe，加入系统 PATH 或放入 core_platform/bin/ |
| xtb | 量子化学优化 | 下载 xtb，加入系统 PATH 或放入 core_platform/bin/ |
| OpenMM | 分子动力学 | `pip install openmm pdbfixer` |
| Ollama（可选） | AI 报告 | 安装 Ollama，拉取 qwen2.5:1.5b 模型 |

### 第五步：启动平台
双击运行 `start_windows.bat`

浏览器访问 http://localhost:5589
- 管理员账号：`admin` / `admin123456`

## 论文案例复现

### 案例一：7种黄酮对 FAK 虚拟筛选
```bash
cd examples/case1_virtual_screening
python run_case1.py
```

### 案例二：柚皮素七阶段全流程
```bash
# 先启动平台（双击 start_windows.bat）
cd examples/case2_full_pipeline
python run_case2_api.py
```

## 常见问题

**Q: 提示找不到 vina / xtb**
A: 确保 vina.exe 和 xtb.exe 所在目录已加入系统 PATH。

**Q: OpenMM 安装失败**
A: 建议用 conda 安装：`conda install -c conda-forge openmm pdbfixer`

**Q: 数据库加载失败**
A: 检查 tcm.db 是否在 data/ 目录下，或设置环境变量 TCMFLOW_DB_PATH。

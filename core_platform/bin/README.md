[Uploading README.md…]()
# 二进制工具目录

本目录存放外部计算引擎的可执行文件。

## 已包含

- ina.exe — AutoDock Vina v1.2.7（Apache 2.0 协议，已内置）

## 需要自行下载

### xtb（必需，用于量子化学优化）

- 版本：xtb 6.7.1
- 下载地址：https://github.com/grimme-lab/xtb/releases
- 安装方法：
  1. 下载 xtb-6.7.1-windows-x86_64.zip
  2. 解压到 in/xtb-6.7.1/ 目录
  3. 确保 in/xtb-6.7.1/bin/xtb.exe 存在
- 或设置环境变量：
  `
  set XTB_EXE=C:\path\to\xtb.exe
  `

### Ollama（可选，用于 AI 报告）

- 下载地址：https://ollama.com/
- 安装后拉取模型：
  `
  ollama pull qwen2.5:1.5b-instruct
  `
- 不安装也能跑，报告自动降级为模板模式。

@echo off
chcp 65001 >nul
echo ====================================
echo   TCM-Flow 一键安装依赖
echo ====================================
echo.

echo [1/3] 检查 Python...
python --version
if errorlevel 1 (
    echo 错误: 未找到 Python，请先安装 Python 3.8-3.11
    echo 下载地址: https://www.python.org/downloads/
    pause
    exit /b 1
)

echo.
echo [2/3] 安装 Python 依赖...
pip install -r requirements.txt
if errorlevel 1 (
    echo 警告: 部分依赖安装失败
    echo 可以尝试: pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple
)

echo.
echo [3/3] 初始化数据库...
python scripts\init_db.py

echo.
echo ====================================
echo   安装完成！
echo ====================================
echo.
echo 下一步:
echo   1. 下载中药数据库 tcm.db 放入 data/ 目录
echo   2. 安装计算引擎（Vina, xtb, OpenMM）
echo   3. 双击 start_windows.bat 启动平台
echo.
pause
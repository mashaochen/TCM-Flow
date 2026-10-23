@echo off
chcp 65001 >nul
echo ====================================
echo   TCM-Flow 中医药分子计算平台
echo ====================================
echo.

echo [1/1] 启动 TCM-Flow 核心计算平台...
cd /d "%~dp0core_platform"
python app.py

pause
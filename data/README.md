# 数据目录说明

## 中药活性成分数据库

完整的中药活性成分数据库（含 3,819 味中药、7,284 种活性成分、靶点、通路信息）
体积约 100 MB，超出 GitHub 文件大小限制。

### 获取方式

**方式一：从百度网盘 / 阿里云盘下载**
> 请在论文补充材料中提供下载链接

**方式二：使用脚本从公开数据源构建**
```bash
python ../scripts/build_tcm_db.py
```
该脚本从 TCMSP、BATMAN-TCM 等公开数据库抓取数据并整合。

**方式三：使用示例数据（仅含少量中药）**
仓库中包含 `sample_tcm.db`（约 1MB），包含 10 味常用中药的
活性成分数据，可用于快速测试平台功能。

## 数据格式

数据库为 SQLite 格式，包含以下主要表：

| 表名 | 说明 |
|------|------|
| herbs | 中药基本信息（名称、拼音、性味归经、功效） |
| compounds | 活性成分信息（名称、SMILES、OB、DL） |
| herb_compound | 中药-成分关联 |
| targets | 靶点信息（基因名、蛋白名、UniProt ID） |
| compound_target | 成分-靶点关联 |
| kegg_pathways | KEGG 通路信息 |
| target_pathway | 靶点-通路关联 |

## 环境变量

通过环境变量指定数据库路径：

```bash
# Windows
set TCM_DB_PATH=C:\path\to\tcm.db

# Linux / macOS
export TCM_DB_PATH=/path/to/tcm.db
```

# tlsfinder

**高效空间TLS分块分析与可视化工具包**

## 主要功能

- 自动遍历所有样本，筛选高TLS_score区域
- 支持聚类分块（DBSCAN）空间定位和最大聚类识别
- 批量导出高分区域聚类标签、坐标、表达矩阵等
- 支持分块空间可视化（HE底图叠加/分块编号/轮廓填充等）
- 可灵活扩展分块参数、分位数阈值、聚类方法、可视化色板

---

## 安装

```bash
pip install tlsfinder
# 或直接源码导入
```

---

## 快速上手

### 1. 数据准备

```python
import tlsfinder
import scanpy as sc
# 读取 AnnData 对象
adata = sc.read_h5ad('your_data.h5ad')
```

### 2. 分块参数设置（每个样本可自定义）

```python
dbscan_params = {
    'GSM5924041_ffpe_c_51': {'eps': 100, 'min_samples': 5},
    'GSM5924046_frozen_b_1': {'eps': 250, 'min_samples': 5},
    'GSM5924047_frozen_b_7': {'eps': 350, 'min_samples': 3}
}
```

### 3. 批量空间TLS高分区域分块分析

```python
from tlsfinder.gr import tls_dbscan_blocks

block_dict = tls_dbscan_blocks(
    adata,
    score_key='TLS_score',
    cluster_key='spatial_cluster',
    sample_key='sample',
    dbscan_params=dbscan_params,
    quantile=0.88
)
```

### 4. 高分区域可视化

```python
from tlsfinder.pl import plot_tls_blocks_on_HE, plot_tls_blocks_labeled

# 可视化分块轮廓+HE底图
plot_tls_blocks_on_HE(block_dict)

# 可视化分块编号着色+HE底图
plot_tls_blocks_labeled(block_dict)
```

---

## 核心API

### `tls_dbscan_blocks(...)`
- 输入：AnnData对象、分数/聚类/样本标签、分块参数等
- 输出：分块分析结果字典（坐标、分块轮廓、mask、图片等）
- 自动去除 AnnData.raw 以节省内存，适合大数据批量处理

### `plot_tls_blocks_on_HE(block_dict, ...)`
- 输入：`tls_dbscan_blocks` 的结果
- 功能：批量可视化高TLS区域轮廓/HE底图/高分点

### `plot_tls_blocks_labeled(block_dict, ...)`
- 输入：`tls_dbscan_blocks` 的结果
- 功能：分块区域自动编号着色+HE底图

---

## 进阶功能

- 导出高分区域表达矩阵
- 支持自定义分块参数、分位数阈值
- 支持多样本批量分析与可视化
- 兼容 cellcharter/squidpy 聚类结果

---


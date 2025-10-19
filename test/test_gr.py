import pytest

# 假设你已在 src/celltry/gr/__init__.py 导出了 enrichment_analysis
from celltry.gr import enrichment_analysis, tls_dbscan_blocks

@pytest.fixture
def adata_tiny():
    """生成一个最小的AnnData对象用于测试。"""
    import anndata
    import numpy as np
    import pandas as pd
    X = np.random.rand(4, 3)
    obs = pd.DataFrame({
        "sample": ["A", "A", "B", "B"],
        "spatial_cluster": [0, 1, 0, 1],
        "cell_type": ["T", "B", "T", "B"]
    })
    var = pd.DataFrame(index=["gene1", "gene2", "gene3"])
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    # 必需的空间属性
    adata.obsm["spatial"] = np.random.rand(4, 2)
    return adata

def test_enrichment_analysis_runs(adata_tiny):
    """测试 enrichment_analysis 能否顺利运行（仅不抛错误）"""
    enrichment_analysis(adata_tiny, group_key="spatial_cluster", label_key="cell_type")

def test_tls_dbscan_blocks_runs(adata_tiny):
    """测试 tls_dbscan_blocks 能否输出结果字典"""
    result = tls_dbscan_blocks(adata_tiny)
    assert isinstance(result, dict)
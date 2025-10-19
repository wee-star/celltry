import scanpy as sc
import numpy as np
from scipy import sparse

def grouped_scale(adata, group_key='sample'):
    """
    对 AnnData 按 group_key 分组分别标准化表达矩阵（sc.pp.scale），再合并回去。
    """
    scaled_X = []
    obs_order = []
    for sample in adata.obs[group_key].cat.categories:
        mask = (adata.obs[group_key] == sample).values
        adata_sub = adata[mask].copy()
        sc.pp.scale(adata_sub)
        scaled_X.append(
            adata_sub.X if isinstance(adata_sub.X, np.ndarray) else adata_sub.X.toarray()
        )
        obs_order.extend(adata_sub.obs_names)
    adata = adata[obs_order, :].copy()
    adata.X = np.vstack(scaled_X)
    adata.raw = adata.copy()
    for sample in adata.obs[group_key].cat.categories:
        mask = (adata.obs[group_key] == sample).values
        adata.X[mask, :] = sc.pp.scale(adata[mask], copy=True).X
    return adata
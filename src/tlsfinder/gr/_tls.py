import os
import numpy as np
import pandas as pd

def export_tls_high_score_regions(
    adata,
    score_key='TLS_score',
    cluster_key='spatial_cluster',
    sample_key='sample',
    out_dir='data',
    quantile=0.85,
    expr=True
):
    """
    自动遍历所有样本，导出每个样本高TLS_score区域的聚类、坐标、表达信息。
    """
    os.makedirs(out_dir, exist_ok=True)
    samples = adata.obs[sample_key].cat.categories.tolist() if pd.api.types.is_categorical_dtype(adata.obs[sample_key]) else adata.obs[sample_key].unique().tolist()
    for lib_id in samples:
        print(f'处理样本: {lib_id}')
        adata_sub = adata[adata.obs[sample_key] == lib_id].copy()
        # 自动设定TLS_score阈值（如85分位）
        score_thr = np.percentile(adata_sub.obs[score_key], int(quantile*100))
        print(f"{lib_id} {score_key}阈值: {score_thr:.3f}")
        high_tls_spots_idx = adata_sub.obs.index[adata_sub.obs[score_key] > score_thr]
        print(f'高{score_key} spot数量: {len(high_tls_spots_idx)}')
        if len(high_tls_spots_idx) == 0:
            continue
        adata_tls = adata_sub[high_tls_spots_idx, :].copy()
        selected_data = adata_tls.obs.loc[:, [score_key, cluster_key, sample_key]]
        spatial_coords = adata_tls.obsm['spatial']
        selected_data['spatial_1'] = spatial_coords[:, 0]
        selected_data['spatial_2'] = spatial_coords[:, 1]
        if expr:
            if hasattr(adata_tls.X, "toarray"):
                expr_df = pd.DataFrame(adata_tls.X.toarray(), index=adata_tls.obs_names, columns=adata_tls.var_names)
            else:
                expr_df = pd.DataFrame(adata_tls.X, index=adata_tls.obs_names, columns=adata_tls.var_names)
            output_df = pd.concat([selected_data, expr_df], axis=1)
        else:
            output_df = selected_data
        out_path = os.path.join(out_dir, f"{lib_id}_TLS_high_spots_with_cluster.csv")
        output_df.to_csv(out_path, index=True)
        print(f"已导出: {out_path}")
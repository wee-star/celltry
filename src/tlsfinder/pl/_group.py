import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

def plot_spatial_clusters_by_slice(
    adata,
    cluster_key='spatial_cluster',
    library_key='sample',
    size=1.5,
    legend_loc='right margin'
):
    # 自动保证聚类标签为 categorical
    if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
        adata.obs[cluster_key] = pd.Categorical(adata.obs[cluster_key])

    # 自动保证颜色条存在且长度正确
    color_key = f'{cluster_key}_colors'
    if color_key not in adata.uns or len(adata.uns[color_key]) != adata.obs[cluster_key].nunique():
        n_clusters = adata.obs[cluster_key].nunique()
        adata.uns[color_key] = [matplotlib.colors.to_hex(plt.cm.tab20(i / n_clusters)) for i in range(n_clusters)]

    # 分样本画图
    if pd.api.types.is_categorical_dtype(adata.obs[library_key]):
        lib_ids = adata.obs[library_key].cat.categories.tolist()
    else:
        lib_ids = adata.obs[library_key].unique().tolist()

    for lib in lib_ids:
        sq.pl.spatial_scatter(
            adata,
            color=cluster_key,
            library_key=library_key,
            library_id=[lib],
            img=None,
            size=size,
            legend_loc=legend_loc,
            title=f"{lib}"
        )

def plot_enrichment(
    adata,
    group_key='spatial_cluster',
    label_key='cell_type',
    figsize=(4,4),
    fontsize=6,
    dot_scale=1
):
    """
    领域聚类的富集分析+点图可视化，一步完成。
    """
    import cellcharter as cc
    cc.gr.enrichment(adata, group_key=group_key, label_key=label_key)
    cc.pl.enrichment(
        adata,
        group_key=group_key,
        label_key=label_key,
        figsize=figsize,
        fontsize=fontsize,
        dot_scale=dot_scale
    )
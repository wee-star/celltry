import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import squidpy as sq

def plot_spatial_clusters_by_slice(
    adata,
    cluster_key='spatial_cluster',
    library_key='sample',
    size=1.5,
    legend_loc='right margin'
):
    """
    按切片（library）分别画出空间聚类结果。
    """
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

def enrichment_analysis(
    adata,
    group_key='spatial_cluster',
    label_key='cell_type',
    figsize=(4,4),
    fontsize=6,
    dot_scale=1
):
    """
    进行领域聚类的富集分析，并画富集点图。
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
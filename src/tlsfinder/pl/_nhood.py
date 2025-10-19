import pandas as pd
import squidpy as sq

def plot_spatial_score_or_neighbors(
    adata,
    color='TLS_score',
    sample_key='sample',
    cmap='Reds',
    size=1.5,
    alpha=0.8,
    connectivity_key=None,
    edges_width=0.3,
    legend_loc=None,
    per_sample=True,
    show_img=False,   # 新增参数，控制是否只画底图
):
    """
    通用空间属性/邻域可视化。
    - color: 要展示的属性（如 TLS_score、或 sample、或 cluster）
    - connectivity_key: 若为 None，则只画点，不连线；否则画邻接（如'spatial_connectivities'）。
    - per_sample: 若True则分样本逐个展示；若False则整体画一张。
    - show_img: 若True只画底图，不画点和连线
    """
    if per_sample:
        # 自动获取所有样本名
        if pd.api.types.is_categorical_dtype(adata.obs[sample_key]):
            samples = adata.obs[sample_key].cat.categories.tolist()
        else:
            samples = adata.obs[sample_key].unique().tolist()
        for lib_id in samples:
            adata_sub = adata[adata.obs[sample_key] == lib_id].copy()
            # 保证uns['spatial']只保留当前样本的空间信息
            if 'spatial' in adata.uns and lib_id in adata.uns['spatial']:
                adata_sub.uns['spatial'] = {lib_id: adata.uns['spatial'][lib_id]}
            else:
                print(f"Warning: 空间信息缺失 {lib_id}")
                continue
            sq.pl.spatial_scatter(
                adata_sub,
                color=None if show_img else color,
                cmap=cmap,
                size=size if not show_img else None,
                alpha=alpha if not show_img else None,
                connectivity_key=None if show_img else connectivity_key,
                edges_width=edges_width if not show_img else None,
                legend_loc=legend_loc,
                title=str(lib_id),
                img=lib_id if show_img else None,
            )
    else:
        sq.pl.spatial_scatter(
            adata,
            color=None if show_img else color,
            cmap=cmap,
            size=size if not show_img else None,
            alpha=alpha if not show_img else None,
            connectivity_key=None if show_img else connectivity_key,
            edges_width=edges_width if not show_img else None,
            legend_loc=legend_loc,
            library_key=sample_key,
            img=None if not show_img else list(adata.uns['spatial'].keys())[0],
        )
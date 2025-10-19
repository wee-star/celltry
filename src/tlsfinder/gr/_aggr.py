def enrichment_analysis(
    adata,
    group_key='spatial_cluster',
    label_key='cell_type'
):
    """
    进行领域聚类的富集分析（仅计算，不可视化）。
    """
    import cellcharter as cc
    cc.gr.enrichment(adata, group_key=group_key, label_key=label_key)
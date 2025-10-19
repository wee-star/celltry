import squidpy as sq

def compute_spatial_neighbors(adata, library_key='sample', coord_type='generic', delaunay=True):
    """
    计算空间邻接矩阵（空间图）。
    """
    sq.gr.spatial_neighbors(adata, library_key=library_key, coord_type=coord_type, delaunay=delaunay)
    return adata

def remove_long_links(adata):
    """
    移除空间邻域的长边（如cellcharter方法）。
    """
    import cellcharter as cc
    cc.gr.remove_long_links(adata)
    return adata
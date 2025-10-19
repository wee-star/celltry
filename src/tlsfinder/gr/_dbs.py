import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from scipy.ndimage import binary_dilation, gaussian_filter
from skimage.measure import find_contours

def tls_dbscan_blocks(
    adata,
    score_key='TLS_score',
    cluster_key='spatial_cluster',
    sample_key='sample',
    dbscan_params=None,
    quantile=0.88,
):
    """
    对所有样本高TLS区域的最大聚类做DBSCAN分块，并返回各样本的spot坐标、contour、mask等。
    切片后主动 adata_sub.raw = None 节省内存，不影响后续空间分析和可视化。
    """
    samples = adata.obs[sample_key].cat.categories.tolist() if pd.api.types.is_categorical_dtype(adata.obs[sample_key]) else adata.obs[sample_key].unique().tolist()
    all_pts, all_contours, all_imgs, all_extents, all_lib_ids, all_mask_bin, all_adata_tls_dbscan = [], [], [], [], [], [], []
    for lib_id in samples:
        adata_sub = adata[adata.obs[sample_key] == lib_id].copy()
        adata_sub.raw = None  # 关键：立刻去掉raw，大幅降低内存
        img_HE = adata.uns['spatial'][lib_id]['images']['hires']
        scalef = adata.uns['spatial'][lib_id]['scalefactors']['tissue_hires_scalef']
        h, w = img_HE.shape[:2]
        extent = (0, w, h, 0)
        tls_scores = adata_sub.obs[score_key]
        tls_thr = np.percentile(tls_scores, int(quantile*100))
        high_tls_spots_idx = adata_sub.obs.index[tls_scores >= tls_thr]
        if len(high_tls_spots_idx) == 0:
            continue
        adata_tls_spots = adata_sub[high_tls_spots_idx, :]
        coords = adata_tls_spots.obsm['spatial']
        center_coord = coords.mean(axis=0)
        dists = np.linalg.norm(coords - center_coord, axis=1)
        dist_thr = np.percentile(dists, 90)
        core_spots_mask = dists <= dist_thr
        core_spots_idx = adata_tls_spots.obs.index[core_spots_mask]
        if len(core_spots_idx) == 0:
            continue
        core_tls_scores = adata_sub.obs.loc[core_spots_idx, score_key]
        core_spots_idx_final = core_spots_idx[core_tls_scores >= tls_thr]
        if len(core_spots_idx_final) == 0:
            continue
        adata_tls = adata_sub[core_spots_idx_final, :].copy()
        adata_tls.raw = None  # 再次节省内存
        cluster_counts = adata_tls.obs[cluster_key].value_counts()
        if cluster_counts.empty:
            continue
        max_cluster = cluster_counts.idxmax()
        adata_tls_max_cluster = adata_tls[adata_tls.obs[cluster_key] == max_cluster].copy()
        adata_tls_max_cluster.raw = None  # 再次节省内存
        if dbscan_params is not None and lib_id in dbscan_params:
            eps = dbscan_params[lib_id]['eps']
            min_samples = dbscan_params[lib_id]['min_samples']
        else:
            eps, min_samples = 200, 5
        coords_all = adata_tls_max_cluster.obsm['spatial']
        db = DBSCAN(eps=eps, min_samples=min_samples)
        db_labels = db.fit_predict(coords_all)
        adata_tls_max_cluster.obs['dbscan_label'] = db_labels
        adata_tls_dbscan = adata_tls_max_cluster[adata_tls_max_cluster.obs['dbscan_label'] != -1].copy()
        adata_tls_dbscan.raw = None  # 清理raw
        all_adata_tls_dbscan.append(adata_tls_dbscan)
        pts = adata_tls_dbscan.obsm['spatial'] * scalef
        if len(pts) < 3:
            continue
        grid_size = 800
        Xg = np.linspace(0, w, grid_size)
        Yg = np.linspace(0, h, grid_size)
        mask_img = np.zeros((grid_size, grid_size), dtype=np.uint8)
        x_idx = np.clip(np.searchsorted(Xg, pts[:,0]), 0, grid_size-1)
        y_idx = np.clip(np.searchsorted(Yg, pts[:,1]), 0, grid_size-1)
        mask_img[y_idx, x_idx] = 1
        mask_img = binary_dilation(mask_img, iterations=3)
        mask_blur = gaussian_filter(mask_img.astype(float), sigma=6)
        mask_bin = (mask_blur > mask_blur.max() * 0.20).astype(np.uint8)
        contours = find_contours(mask_bin, 0.5)
        all_pts.append(pts)
        all_contours.append(contours)
        all_imgs.append(img_HE)
        all_extents.append(extent)
        all_lib_ids.append(lib_id)
        all_mask_bin.append(mask_bin)
    return dict(
        pts=all_pts,
        contours=all_contours,
        imgs=all_imgs,
        extents=all_extents,
        lib_ids=all_lib_ids,
        mask_bin=all_mask_bin,
        adata_tls_dbscan=all_adata_tls_dbscan
    )
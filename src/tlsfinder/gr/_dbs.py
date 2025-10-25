import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from skimage.measure import find_contours
from scipy.ndimage import binary_dilation, gaussian_filter

def tls_dbscan_blocks(
    adata,
    score_key='TLS_score',
    cluster_key='spatial_cluster',
    sample_key='sample',
    dbscan_params=None,
    quantile=0.88,
    include_all=False,   # 新增：是否为被跳过的样本也返回占位数据
    verbose=False
):
    """
    对所有样本高TLS区域的最大聚类做DBSCAN分块，并返回各样本的spot坐标、contour、mask等。
    include_all=True 时，即便样本在某一步被过滤（无高分点/无核心点/DBSCAN无点等），
    也会为该样本添加占位项（空 pts、空 contours、全0 mask），以保证绘图时每个样本都有对应子图。
    """
    samples = adata.obs[sample_key].cat.categories.tolist() if pd.api.types.is_categorical_dtype(adata.obs[sample_key]) else adata.obs[sample_key].unique().tolist()
    all_pts, all_contours, all_imgs, all_extents, all_lib_ids, all_mask_bin, all_adata_tls_dbscan = [], [], [], [], [], [], []
    for lib_id in samples:
        if verbose: print("处理样本:", lib_id)
        adata_sub = adata[adata.obs[sample_key] == lib_id].copy()
        adata_sub.raw = None  # 关键：立刻去掉raw，大幅降低内存
        img_HE = adata.uns['spatial'][lib_id]['images']['hires']
        scalef = adata.uns['spatial'][lib_id]['scalefactors']['tissue_hires_scalef']
        h, w = img_HE.shape[:2]
        extent = (0, w, h, 0)

        # 计算高分阈值
        tls_scores = adata_sub.obs[score_key]
        tls_thr = np.percentile(tls_scores, int(quantile*100))
        high_tls_spots_idx = adata_sub.obs.index[tls_scores >= tls_thr]
        if verbose: print(f"  TLS_thr={tls_thr:.4f} high_tls_spots={len(high_tls_spots_idx)}")
        if len(high_tls_spots_idx) == 0:
            if include_all:
                all_pts.append(np.empty((0,2)))
                all_contours.append([])
                all_imgs.append(img_HE)
                all_extents.append(extent)
                all_lib_ids.append(lib_id)
                all_mask_bin.append(np.zeros((800,800), dtype=np.uint8))
                all_adata_tls_dbscan.append(None)
            else:
                if verbose: print("  -> skip: no high-tls spots")
            continue

        adata_tls_spots = adata_sub[high_tls_spots_idx, :]
        coords = adata_tls_spots.obsm['spatial']
        center_coord = coords.mean(axis=0)
        dists = np.linalg.norm(coords - center_coord, axis=1)
        dist_thr = np.percentile(dists, 90)
        core_spots_mask = dists <= dist_thr
        core_spots_idx = adata_tls_spots.obs.index[core_spots_mask]
        if verbose: print(f"  core_spots(before score filter)={len(core_spots_idx)} dist_thr={dist_thr:.2f}")
        if len(core_spots_idx) == 0:
            if include_all:
                all_pts.append(np.empty((0,2)))
                all_contours.append([])
                all_imgs.append(img_HE)
                all_extents.append(extent)
                all_lib_ids.append(lib_id)
                all_mask_bin.append(np.zeros((800,800), dtype=np.uint8))
                all_adata_tls_dbscan.append(None)
            else:
                if verbose: print("  -> skip: no core spots")
            continue

        core_tls_scores = adata_sub.obs.loc[core_spots_idx, score_key]
        core_spots_idx_final = core_spots_idx[core_tls_scores >= tls_thr]
        if verbose: print(f"  core_spots(after score filter)={len(core_spots_idx_final)}")
        if len(core_spots_idx_final) == 0:
            if include_all:
                all_pts.append(np.empty((0,2)))
                all_contours.append([])
                all_imgs.append(img_HE)
                all_extents.append(extent)
                all_lib_ids.append(lib_id)
                all_mask_bin.append(np.zeros((800,800), dtype=np.uint8))
                all_adata_tls_dbscan.append(None)
            else:
                if verbose: print("  -> skip: no core spots passing score threshold")
            continue

        adata_tls = adata_sub[core_spots_idx_final, :].copy()
        adata_tls.raw = None  # 再次节省内存

        if cluster_key not in adata_tls.obs.columns:
            if include_all:
                all_pts.append(np.empty((0,2)))
                all_contours.append([])
                all_imgs.append(img_HE)
                all_extents.append(extent)
                all_lib_ids.append(lib_id)
                all_mask_bin.append(np.zeros((800,800), dtype=np.uint8))
                all_adata_tls_dbscan.append(None)
            else:
                if verbose: print("  -> skip: cluster_key missing")
            continue

        cluster_counts = adata_tls.obs[cluster_key].value_counts()
        if verbose: print(f"  cluster_counts:\n{cluster_counts}")
        if cluster_counts.empty:
            if include_all:
                all_pts.append(np.empty((0,2)))
                all_contours.append([])
                all_imgs.append(img_HE)
                all_extents.append(extent)
                all_lib_ids.append(lib_id)
                all_mask_bin.append(np.zeros((800,800), dtype=np.uint8))
                all_adata_tls_dbscan.append(None)
            else:
                if verbose: print("  -> skip: no clusters in core TLS spots")
            continue

        max_cluster = cluster_counts.idxmax()
        adata_tls_max_cluster = adata_tls[adata_tls.obs[cluster_key] == max_cluster].copy()
        adata_tls_max_cluster.raw = None  # 再次节省内存

        # DBSCAN 参数
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
        if verbose: print(f"  dbscan: raw_points={len(coords_all)} kept_after_dbscan={adata_tls_dbscan.n_obs if hasattr(adata_tls_dbscan,'n_obs') else len(adata_tls_dbscan.obs)}")

        pts = adata_tls_dbscan.obsm['spatial'] * scalef
        # 若点太少，根据 include_all 决定是否返回占位
        if len(pts) < 3:
            if include_all:
                all_pts.append(pts)  # 可能是空或少量点
                all_contours.append([])
                all_imgs.append(img_HE)
                all_extents.append(extent)
                all_lib_ids.append(lib_id)
                all_mask_bin.append(np.zeros((800,800), dtype=np.uint8))
                all_adata_tls_dbscan.append(adata_tls_dbscan)
            else:
                if verbose: print("  -> skip: pts < 3 after dbscan")
            continue

        # 生成 mask/contour（同你原代码）
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
        all_adata_tls_dbscan.append(adata_tls_dbscan)

    return dict(
        pts=all_pts,
        contours=all_contours,
        imgs=all_imgs,
        extents=all_extents,
        lib_ids=all_lib_ids,
        mask_bin=all_mask_bin,
        adata_tls_dbscan=all_adata_tls_dbscan
    )
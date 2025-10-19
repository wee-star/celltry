import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np

def plot_tls_blocks_on_HE(
    block_dict,
    fill_color=(225/255, 80/255, 160/255, 0.5),
    grid_size=800
):
    """
    可视化高分区域聚类分块与HE底图轮廓。
    block_dict: gr/_nhood.py tls_dbscan_blocks返回的结果
    """
    all_pts = block_dict['pts']
    all_contours = block_dict['contours']
    all_imgs = block_dict['imgs']
    all_extents = block_dict['extents']
    all_lib_ids = block_dict['lib_ids']
    num = len(all_pts)
    fig, axes = plt.subplots(1, num, figsize=(5*num, 5))
    if num == 1:
        axes = [axes]
    for i in range(num):
        ax = axes[i]
        pts = all_pts[i]
        contours = all_contours[i]
        img_HE = all_imgs[i]
        extent = all_extents[i]
        lib_id = all_lib_ids[i]
        h, w = img_HE.shape[:2]
        Xg = np.linspace(0, w, grid_size)
        Yg = np.linspace(0, h, grid_size)
        ax.imshow(img_HE, extent=extent, origin='upper')
        for contour in contours:
            contour = np.flip(contour, axis=1)
            contour_x = np.interp(contour[:,0], np.arange(grid_size), Xg)
            contour_y = np.interp(contour[:,1], np.arange(grid_size), Yg)
            poly = Polygon(np.c_[contour_x, contour_y], closed=True, facecolor=fill_color, edgecolor='none')
            ax.add_patch(poly)
        ax.scatter(pts[:,0], pts[:,1], c='magenta', s=8, alpha=0.7, edgecolors='w', linewidths=0.2)
        ax.set_title(f"{lib_id} TLS region", fontsize=14)
        ax.axis('off')
    plt.tight_layout()
    plt.show()

def plot_tls_blocks_labeled(
    block_dict,
    grid_size=800,
    cmap_name='tab10'
):
    """
    分块区域着色并编号可视化。
    """
    from skimage.measure import label, regionprops
    from matplotlib import cm
    all_mask_bin = block_dict['mask_bin']
    all_imgs = block_dict['imgs']
    all_extents = block_dict['extents']
    all_pts = block_dict['pts']
    all_lib_ids = block_dict['lib_ids']
    num = len(all_mask_bin)
    fig, axes = plt.subplots(1, num, figsize=(5*num, 5))
    if num == 1:
        axes = [axes]
    colormap = cm.get_cmap(cmap_name, 10)
    for i in range(num):
        ax = axes[i]
        img_HE = all_imgs[i]
        extent = all_extents[i]
        mask_bin = all_mask_bin[i]
        pts = all_pts[i]
        lib_id = all_lib_ids[i]
        h, w = img_HE.shape[:2]
        Xg = np.linspace(0, w, grid_size)
        Yg = np.linspace(0, h, grid_size)
        ax.imshow(img_HE, extent=extent, origin='upper')
        label_img = label(mask_bin, connectivity=2)
        regions = regionprops(label_img)
        for region in regions:
            region_mask = (label_img == region.label).astype(float)
            from skimage.measure import find_contours
            region_contours = find_contours(region_mask, 0.5)
            color = colormap((region.label-1)%10)
            for contour in region_contours:
                contour = np.flip(contour, axis=1)
                contour_x = np.interp(contour[:,0], np.arange(grid_size), Xg)
                contour_y = np.interp(contour[:,1], np.arange(grid_size), Yg)
                poly = Polygon(np.c_[contour_x, contour_y], closed=True, facecolor=color, edgecolor='black', linewidth=1, alpha=0.5)
                ax.add_patch(poly)
            cy, cx = region.centroid
            cx_disp = np.interp(cx, np.arange(grid_size), Xg)
            cy_disp = np.interp(cy, np.arange(grid_size), Yg)
            ax.text(cx_disp, cy_disp, str(region.label), color='black', fontsize=16, ha='center', va='center', weight='bold')
        ax.scatter(pts[:,0], pts[:,1], c='magenta', s=8, alpha=0.7, edgecolors='w', linewidths=0.2)
        ax.set_title(f"{lib_id} TLS blocks", fontsize=14)
        ax.axis('off')
    plt.tight_layout()
    plt.show()
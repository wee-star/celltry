import scanpy as sc
import pandas as pd
import numpy as np
import os
import anndata
from PIL import Image
import glob

def merge_spatial_data(datapath, output_h5ad):
    """
    批量读取原始空间转录组数据，预处理并合并为一个AnnData h5ad文件。
    每个样本的 obs['sample'] 为类别型变量，保证可视化区分颜色。
    """
    adatas = []
    all_h5_files = glob.glob(os.path.join(datapath, "*_filtered_feature_bc_matrix.h5"))
    prefixes = [os.path.basename(f).replace("_filtered_feature_bc_matrix.h5", "") for f in all_h5_files]

    for prefix in prefixes:
        expr_h5 = os.path.join(datapath, f"{prefix}_filtered_feature_bc_matrix.h5")
        pos_csv = os.path.join(datapath, f"{prefix}_tissue_positions_list.csv")
        hires_img = os.path.join(datapath, f"{prefix}_tissue_hires_image.png")
        lowres_img = os.path.join(datapath, f"{prefix}_tissue_lowres_image.png")
        scalefactors_json = os.path.join(datapath, f"{prefix}_scalefactors_json.json")
        anno_csv = os.path.join(datapath, f"{prefix}_TLS_annotation.csv")

        files_needed = [expr_h5, pos_csv, hires_img, lowres_img, scalefactors_json, anno_csv]
        if not all(os.path.exists(f) for f in files_needed):
            print(f"缺少文件，跳过样本: {prefix}")
            continue

        adata = sc.read_10x_h5(expr_h5)
        adata.var_names_make_unique()
        adata.var = pd.DataFrame(index=adata.var_names)
        adata.obs_names = [f"{bc}_{prefix}" for bc in adata.obs_names]
        if adata.raw is not None:
            adata.raw._obs_names = adata.obs_names

        positions = pd.read_csv(pos_csv, header=None)
        positions.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
        positions['barcode'] = positions['barcode'].astype(str) + f"_{prefix}"
        adata = adata[adata.obs_names.isin(positions['barcode'])].copy()
        adata.obs = adata.obs.join(positions.set_index('barcode'), how='left')

        anno = pd.read_csv(anno_csv)
        anno['Barcode'] = anno['Barcode'].astype(str) + f"_{prefix}"
        anno = anno.set_index('Barcode')
        adata.obs['cell_type'] = adata.obs.index.map(anno['TLS_2_cat']).fillna("NA")

        for col in ['i-niche', 'tile', 'area', 'dataset', 'stage']:
            adata.obs[col] = "NA"
        adata.obs['sample'] = prefix
        adata.obs = adata.obs[['cell_type', 'i-niche', 'tile', 'area', 'dataset', 'stage', 'sample',
                               'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']]

        spatial_coords = positions.set_index('barcode').loc[adata.obs_names, ['pxl_row_in_fullres', 'pxl_col_in_fullres']].astype(float).values
        adata.obsm['spatial'] = spatial_coords
        adata.obsm['blanks'] = np.zeros((adata.n_obs, 2))

        try:
            scalefactors = pd.read_json(scalefactors_json, typ='series').to_dict()
        except Exception:
            scalefactors = {}
        hires_img_arr = np.array(Image.open(hires_img))
        lowres_img_arr = np.array(Image.open(lowres_img))
        adata.uns['spatial'] = {
            prefix: {
                'images': {
                    'hires': hires_img_arr,
                    'lowres': lowres_img_arr
                },
                'scalefactors': scalefactors
            }
        }
        adata.uns['spatial_cluster_colors'] = []
        adatas.append(adata)

    if len(adatas) == 0:
        raise ValueError("未找到任何可合并的样本！")
    adata_merged = anndata.concat(
        adatas,
        axis=0,
        join='outer',
        label=None,
        merge='unique',
    )
    adata_merged.var = pd.DataFrame(index=adata_merged.var_names)

    # 关键：sample设为类别型，保证可视化不同颜色
    adata_merged.obs['sample'] = pd.Categorical(adata_merged.obs['sample'])

    adata_merged.uns['spatial'] = {}
    for prefix in prefixes:
        for ad in adatas:
            if ad.obs['sample'].iloc[0] == prefix:  # 已用iloc[0]避免FutureWarning
                adata_merged.uns['spatial'][prefix] = ad.uns['spatial'][prefix]
    adata_merged.uns['spatial_cluster_colors'] = []

    adata_merged.write(output_h5ad)
    print(f"Saved merged AnnData to {output_h5ad}")
    return adata_merged

# 用法示例
if __name__ == "__main__":
    datapath = "raw_data"
    output_h5ad = "data/merged_spatial_data_formatted.h5ad"
    adata_merged = merge_spatial_data(datapath, output_h5ad)
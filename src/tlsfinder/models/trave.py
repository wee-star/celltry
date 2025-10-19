import numpy as np
import pandas as pd
from scipy import sparse
import os
from cellcharter.tl import TRVAE

def filter_nan_spots(adata):
    """
    过滤掉表达矩阵中含NaN的spot
    """
    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    nan_spots = np.isnan(X).any(axis=1)
    adata = adata[~nan_spots, :].copy()
    print("过滤后NaN数量:", np.isnan(adata.X.toarray()).sum() if sparse.issparse(adata.X) else np.isnan(adata.X).sum())
    return adata

def ensure_float32(adata):
    """
    保证表达矩阵是float32类型
    """
    if adata.X.dtype != np.float32:
        adata.X = adata.X.astype(np.float32)
    return adata

def align_genes_and_get_latent(adata, model, condition_key='dataset'):
    """
    按模型基因顺序对齐表达矩阵，推理并写入adata.obsm['X_trVAE']
    """
    # 获取基因顺序
    if hasattr(model, "genes"):
        reference_genes = list(model.genes)
    elif hasattr(model, "var_names"):
        reference_genes = list(model.var_names)
    elif hasattr(model, "adata") and hasattr(model.adata, "var_names"):
        reference_genes = list(model.adata.var_names)
    else:
        raise ValueError("找不到模型训练时的基因顺序，请手动指定 reference_genes！")

    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    df_aligned = df.reindex(columns=reference_genes, fill_value=0)
    X_aligned = df_aligned.values.astype(np.float32)
    adata.obsm['X_trVAE'] = model.get_latent(X_aligned, adata.obs[condition_key])
    return adata

def train_trvae(
    adata,
    model_path,
    condition_key='dataset',
    hidden_layer_sizes=[128, 128],
    latent_dim=10,
    dr_rate=0.05,
    use_mmd=True,
    recon_loss='mse',
    use_bn=False,
    use_ln=True,
    n_epochs=250,
    lr=1e-3,
    map_location='cpu',
    verbose=True
):
    """
    加载已保存的TRVAE模型或（如无）新训练并保存。
    返回模型对象。
    """
    if os.path.exists(model_path):
        if verbose:
            print(f"发现已有模型：{model_path}，直接加载。")
        model = TRVAE.load(model_path, adata, map_location=map_location)
    else:
        if verbose:
            print("未发现已有模型，开始训练新模型。")
        model = TRVAE(
            adata,
            condition_key=condition_key,
            hidden_layer_sizes=hidden_layer_sizes,
            latent_dim=latent_dim,
            dr_rate=dr_rate,
            use_mmd=use_mmd,
            recon_loss=recon_loss,
            use_bn=use_bn,
            use_ln=use_ln
        )
        model.train(n_epochs=n_epochs, lr=lr)
        model.save(model_path)
        if verbose:
            print(f"模型已保存到：{model_path}")
    return model
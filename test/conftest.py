import pytest
import numpy as np
import pandas as pd
import anndata

@pytest.fixture
def adata_tiny():
    X = np.random.rand(10, 5)
    obs = pd.DataFrame({'sample': ['s1']*5 + ['s2']*5, 'spatial_cluster': [0,1]*5, 'cell_type': ['T', 'B']*5})
    var = pd.DataFrame(index=[f"g{i}" for i in range(5)])
    return anndata.AnnData(X=X, obs=obs, var=var)
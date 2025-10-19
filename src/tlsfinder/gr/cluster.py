class Cluster:
    def __init__(self, n_clusters, random_state=12345, trainer_params=None):
        import cellcharter as cc
        self.model = cc.tl.Cluster(
            n_clusters=n_clusters, 
            random_state=random_state,
            trainer_params=trainer_params
        )
    def fit(self, adata, use_rep='X_trVAE'):
        self.model.fit(adata, use_rep=use_rep)
    def predict(self, adata, use_rep='X_trVAE'):
        return self.model.predict(adata, use_rep=use_rep)
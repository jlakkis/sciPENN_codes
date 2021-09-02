from scanpy import read_10x_h5, read
from pandas import read_csv, DataFrame

def read_malt_data(path):
    adata_malt_gene = read_10x_h5(path + "/malt_10k_protein_v3_filtered_feature_bc_matrix.h5")
    adata_malt = read(path + "/filtered_feature_bc_matrix/matrix.mtx").T
    
    malt_features =  read_csv(path + "/filtered_feature_bc_matrix/features.tsv", sep="\t", header=None)

    adata_malt.var["feature_type"] = list(malt_features[2])
    adata_malt.obs_names = adata_malt_gene.obs_names
    adata_malt.var['protein_names'] = list(malt_features[0])
    adata_malt.var_names = list(malt_features[0])
    
    adata_malt_protein = adata_malt[:,adata_malt.var['feature_type'] == 'Antibody Capture']

    adata_malt_gene.var_names_make_unique()

    adata_malt_gene.obs['sample'] = [1]*8412
    
    gene_matrix = DataFrame(adata_malt_gene.X.toarray(), index = adata_malt_gene.obs.index, columns = adata_malt_gene.var.index)
    protein_matrix = DataFrame(adata_malt_protein.X.toarray(), index = adata_malt_protein.obs.index, columns = adata_malt_protein.var.index)
    
    return gene_matrix.T, protein_matrix.T
import scanpy as sc
from scanpy import read

from gzip import open as gzip_open
from shutil import copyfileobj
import os
from pandas import read_csv, DataFrame
from numpy import count_nonzero
import numpy as np

def read_data(train = [], test = [], dir_path = 'weights', read_unpaired = False, subset_hvg = True, cell_normalize = True, log_normalize = True, feature_normalize = True, n_protein_HVG = None, n_gene_HVG = None, min_cells = 1, min_genes = 1):
    patients = train + test
    if len(patients) == 0:
        patients = ['RPM211A', 'RPM211B', 'RPM232A', 'RPM232B', 'RPM215A', 'RPM215B', 'RPM218A', 'RPM218B']
        
    dir_paths = [os.path.join(dir_path, patient) for patient in patients]
    i = 0
    
    for dir_path in dir_paths:
        patient = os.path.split(dir_path)[1]
        gzipped_data = dir_path + "/matrix.mtx.gz"

        if not os.path.isfile(gzipped_data[:-3]):
            with gzip_open(gzipped_data, 'rb') as f_in:
                with open(gzipped_data[:-3], 'wb') as f_out:
                    copyfileobj(f_in, f_out)

        adata = read(gzipped_data[:-3]).T
        
        genes = read_csv(dir_path + "/features.tsv", sep='\t', header = None)
        cells = read_csv(dir_path + "/barcodes.tsv", sep='\t', header = None).iloc[:,0].values

        adata.var.index, adata.obs.index = genes.iloc[:,0].values, cells
        adata.var['Common Name'] = genes.iloc[:,1].values
        adata.obs['patient'] = [patient] * adata.shape[0]
        adata.var['expression_type'] = genes.iloc[:,2].values
        adata.X = adata.X.toarray()
        adata.layers['raw'] = adata.X.copy()
        
        adata_gene = adata[:, adata.var['expression_type'] != 'Antibody Capture'].copy()
        adata_protein = adata[:, adata.var['expression_type'] == 'Antibody Capture'].copy()

        index = adata.var['expression_type'] != 'Antibody Capture'
        index_p = adata.var['expression_type'] == 'Antibody Capture'
        
        cell_filter = count_nonzero(adata_gene.X, axis = 1) >= min_cells        
        adata_gene, adata_protein = adata_gene[cell_filter].copy(), adata_protein[cell_filter].copy()
        
        if cell_normalize:
            sc.pp.normalize_total(adata_protein)
        if log_normalize:
            sc.pp.log1p(adata_protein)
        
        if i == 0:
            adatat_gene = adata_gene.copy()
            adatat_protein = adata_protein.copy()
            
        else:
            adata_gene.var = adatat_gene.var
            adata_protein.var = adatat_protein.var
            adatat_gene = adatat_gene.concatenate(adata_gene, batch_key = None, index_unique = None)
            adatat_protein = adatat_protein.concatenate(adata_protein, batch_key = None, index_unique = None)    
            
        i += 1
        
    gene_list = adatat_gene.var['Common Name'].values

    gene_dict = {}
    for gene in gene_list:
        if gene not in gene_dict.keys():
            gene_dict[gene] = 0

        gene_dict[gene] += 1

    gene_list = [gene_dict[gene] < 2 for gene in gene_list]
    adatat_gene = adatat_gene[:, gene_list]
    
    adatat_gene.var['Scientific Name'] = adatat_gene.var.index
    adatat_gene.var.index = adatat_gene.var['Common Name'].values
    
    adatat_gene.obs_names_make_unique()
    adatat_protein.obs_names_make_unique()
        
    if read_unpaired:
        adatat_gene.obs = adatat_gene.obs[['patient']]
        scrna_seq = sc.read("../../Data/monocytes_mingyao/scrna_seq/raw_cnt_mat.csv").T
        rnaseq_metadata = read_csv("../../Data/monocytes_mingyao/scrna_seq/monocyte_integrated_20_metadata.txt", sep = "\t")
        scrna_seq.obs = rnaseq_metadata[['orig.ident']]
        scrna_seq.obs.columns = ['patient']
    
        cell_filter = count_nonzero(scrna_seq.X, axis = 1) >= min_cells 
        scrna_seq = scrna_seq[cell_filter].copy()
        if cell_normalize:
            sc.pp.normalize_total(scrna_seq)
        if log_normalize:
            sc.pp.log1p(scrna_seq)

        adatat_gene = adatat_gene.concatenate(scrna_seq, batch_key = None)
    
    features = adatat_gene.X.sum(axis = 0) >= min_genes
    adatat_gene = adatat_gene[:, features]
    
    if cell_normalize:
        sc.pp.normalize_total(adatat_gene)
    if log_normalize:
        sc.pp.log1p(adatat_gene)
    
    sc.pp.filter_genes(adatat_protein, min_counts = 1)
    
    if n_gene_HVG is not None:
        sc.pp.highly_variable_genes(adatat_gene, min_mean = 0.0125, max_mean = 3, min_disp = 0.5, 
                                      n_bins = 20, subset = False, batch_key = 'patient', n_top_genes = n_gene_HVG)
    if n_protein_HVG is not None:
        train_proteins = adatat_protein[[x in train for x in adatat_protein.obs['patient']]]
        sc.pp.highly_variable_genes(train_proteins, subset = False, batch_key = 'patient', n_top_genes = n_protein_HVG)
        adatat_protein.var['highly_variable'] = train_proteins.var['highly_variable']
    
    if subset_hvg:
        adatat_gene = adatat_gene[:, adatat_gene.var['highly_variable'].values].copy()
        adatat_protein = adatat_protein[:, adatat_protein.var['highly_variable'].values].copy()
    
    patients_rna = np.unique(adatat_gene.obs['patient'])
    patients_protein = np.unique(adatat_protein.obs['patient'])
    
    if feature_normalize:
        for patient in patients_rna:
            indices = [x == patient for x in adatat_gene.obs['patient']]
            sub_adata = adatat_gene[indices]

            sc.pp.scale(sub_adata)
            adatat_gene[indices] = sub_adata.X

        if feature_normalize:
            for patient in patients_protein:
                indices = [x == patient for x in adatat_protein.obs['patient']]
                sub_adata = adatat_protein[indices]

                sc.pp.scale(sub_adata)
                adatat_protein[indices] = sub_adata.X 
    
    return adatat_gene, adatat_protein

def read_preprocess_R(dir_path):
    adata_gene, adata_protein = read_data(cell_normalize = False, log_normalize = False, feature_normalize = False, 
                   dir_path = dir_path, subset_hvg = False)
    
    tmp = adata_protein.copy()
    sc.pp.normalize_total(tmp)
    sc.pp.log1p(tmp)

    sums = tmp.X.sum(axis = 0)
    samples = ((tmp.X > 0.0001).sum(axis = 0))

    expression = sums/samples
    adata_protein = adata_protein[:, expression > 0.8].copy()
    
    train_patientset = ['RPM211A', 'RPM211B', 'RPM232A', 'RPM232B']
    test_patientset = ['RPM215A', 'RPM215B', 'RPM218A', 'RPM218B']

    train_patients = [x in train_patientset for x in adata_gene.obs['patient']]
    test_patients = [x in test_patientset for x in adata_gene.obs['patient']]

    adata_gene, adata_gene_test = adata_gene[train_patients], adata_gene[test_patients]
    adata_protein, adata_protein_test = adata_protein[train_patients], adata_protein[test_patients]
    
    train_metadata, test_metadata = adata_gene.obs, adata_gene_test.obs
    
    adata_gene = DataFrame(adata_gene.X, index = adata_gene.obs.index, columns = adata_gene.var.index)
    adata_gene_test = DataFrame(adata_gene_test.X, index = adata_gene_test.obs.index, columns = adata_gene_test.var.index)
    adata_protein = DataFrame(adata_protein.X, index = adata_protein.obs.index, columns = adata_protein.var.index)
    adata_protein_test = DataFrame(adata_protein_test.X, index = adata_protein_test.obs.index, columns = adata_protein_test.var.index)
    
    return adata_gene.T, adata_protein.T, adata_gene_test.T, adata_protein_test.T, train_metadata, test_metadata
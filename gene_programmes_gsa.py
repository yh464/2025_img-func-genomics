# Gene set enrichment tests for the cNMF programmes
import pandas as pd
import numpy as np
rng = np.random.default_rng(19260817)

cnmf_k = 16; cnmf_dt = '0_15'
cnmf_weights = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/wang_2025/wang_2025_neocortex/' +
    f'wang_2025_neocortex.spectra.k_{cnmf_k}.dt_{cnmf_dt}.consensus.txt',index_col = 0).T
cnmf_weights.columns = [f'F{i+1}' for i in range(cnmf_k)]
replication_k = 14; replication_dt = '0_1'
replication_weights = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/polioudakis_2019/neocx_wang_2025_genes/'+
    f'neocx_wang_2025_genes.spectra.k_{replication_k}.dt_{replication_dt}.consensus.txt',
    index_col = 0
).T
replication_weights.columns = [f'F{i+1}' for i in range(replication_k)]

gene_sets = pd.read_table('./braun_2023_s4.txt')
ref = pd.read_table('/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/params/genes_ref.txt', index_col = 'LABEL')
neuronal_like_ipc_genes = [gene for gene in gene_sets.index[:138] if gene in ref.index]
rg_like_ipc_genes = [gene for gene in gene_sets.index[138:] if gene in ref.index]
neuronal_like_ipc_genes = ref.loc[neuronal_like_ipc_genes, 'GENE'].values
rg_like_ipc_genes = ref.loc[rg_like_ipc_genes, 'GENE'].values

def gene_set_enrichment(weights, gene_set, n_perm = 10000):
    n_overlap = weights.index.intersection(gene_set).size
    mean_weight = weights.loc[weights.index.intersection(gene_set),:].mean(axis = 0).values
    perm_weights = np.stack([
        weights.sample(n = n_overlap, replace = False, random_state = rng).mean(axis = 0).values
        for _ in range(n_perm)
    ], axis = 0)
    p_values = (mean_weight > perm_weights).mean(axis = 0)
    out = pd.DataFrame({'mean_weight': mean_weight, 'p_value': p_values}, index = weights.columns)
    print(out)
    return out

cnmf_gsea_neuronal_like_ipc = gene_set_enrichment(cnmf_weights, neuronal_like_ipc_genes)
cnmf_gsea_rg_like_ipc = gene_set_enrichment(cnmf_weights, rg_like_ipc_genes)
replication_gsea_neuronal_like_ipc = gene_set_enrichment(replication_weights, neuronal_like_ipc_genes)
replication_gsea_rg_like_ipc = gene_set_enrichment(replication_weights, rg_like_ipc_genes)
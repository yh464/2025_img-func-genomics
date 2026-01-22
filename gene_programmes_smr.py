import pandas as pd

gene_sets = pd.read_table('./braun_2023_s4.txt', index_col = 0)
ref = pd.read_table('/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/params/genes_ref.txt', index_col = 'LABEL')
neuronal_like_ipc_genes = [gene for gene in gene_sets.index[:138] if gene in ref.index]
rg_like_ipc_genes = [gene for gene in gene_sets.index[138:] if gene in ref.index]
neuronal_like_ipc_genes = ref.loc[neuronal_like_ipc_genes, 'GENE'].values
rg_like_ipc_genes = ref.loc[rg_like_ipc_genes, 'GENE'].values

cnmf_k = 16; cnmf_dt = '0_15'
cnmf_weights = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/wang_2025/wang_2025_neocortex/' +
    f'wang_2025_neocortex.spectra.k_{cnmf_k}.dt_{cnmf_dt}.consensus.txt',index_col = 0).T
cnmf_weights.columns = [f'F{i+1}' for i in range(cnmf_k)]

smr_files = [
    '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/annot/smr/structural_factors/cortical_expansion.smr',
    '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/annot/smr/disorders/adhd2022.smr',
    '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/annot/smr/disorders/adhd2025.smr'
]

def extract_smr(smr_file, gene_set):
    smr = pd.read_table(smr_file)
    gene_set = pd.Index(gene_set)
    genes = gene_set.intersection(smr['probe'])
    out = smr.loc[smr['probe'].isin(genes),:]
    print(out)
    return out

for smr_file in smr_files:
    print(smr_file)
    print('Neuronal-like IPC genes:')
    extract_smr(smr_file, neuronal_like_ipc_genes)
    print('RG-like IPC genes:')
    extract_smr(smr_file, rg_like_ipc_genes)
    for programme in ['F9','F2','F1','F8']:
        print(f'Top 200 of {programme}:')
        top_200_genes = cnmf_weights.nlargest(200, programme).index
        extract_smr(smr_file, top_200_genes)
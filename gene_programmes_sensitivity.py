import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.cluster.hierarchy import linkage, dendrogram
output_dir = '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/test'

# compare different parameters
cnmf_k = 16; cnmf_dt = '0_15'; cnmf_factor_order = [f'F{i}' for i in [12,15,5,9,14,2,1,8,4,6,7,11,13,3,10,16]]
# cnmf_k = 10; cnmf_dt = '0_1'; factor_order = [f'F{i}' for i in [9,10,4,5,6,7,8,1,2,3]]
cnmf_weights = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/wang_2025/wang_2025_neocortex/' +
    f'wang_2025_neocortex.spectra.k_{cnmf_k}.dt_{cnmf_dt}.consensus.txt',index_col = 0).T
cnmf_weights.columns = [f'F{i+1}' for i in range(cnmf_k)]

sensitivity_k = [13, 14, 15, 16, 17, 18, 19, 20]
sensitivity_weights = [
    pd.read_table(
        f'/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/wang_2025/wang_2025_neocortex/' +
        f'wang_2025_neocortex.spectra.k_{k}.dt_0_15.consensus.txt', index_col = 0).T.rename(columns = {i: f'F{i}' for i in range(1, k+1)})
    for k in sensitivity_k
]

for test_weights,k in zip(sensitivity_weights, sensitivity_k):
    linkage_matrix = linkage(test_weights.T, method = 'average', metric = 'correlation')
    dendrogram_res = dendrogram(linkage_matrix, no_plot = True)
    corr = pd.concat([test_weights.corrwith(cnmf_weights[f'F{i+1}']).to_frame(name = f'F{i+1}') for i in range(cnmf_k)], axis = 1)
    corr = corr.loc[:, cnmf_factor_order].iloc[dendrogram_res['leaves'], :]
    out_prefix = f'{output_dir}/cnmf_corr_k{cnmf_k}_to_k{k}'
    fig = plt.figure(figsize = (10,10))
    sns.heatmap(corr, cmap = 'vlag', center = 0, yticklabels = True, xticklabels = True, square = True)
    fig.savefig(f'{out_prefix}.pdf', bbox_inches = 'tight')
    plt.close()

# rgipc_k = 12; rgipc_dt = '0_1'; rgipc_factor_order = [f'F{i}' for i in [11,10,8,4,3,5,1,2,6,12,7,9]]
rgipc_k = 9; rgipc_dt = '0_1'; rgipc_factor_order = [f'F{i}' for i in [6,1,2,5,3,4,8,7,9]]
rgipc_weights = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/wang_2025/wang_2025_rgipc/' +
    f'wang_2025_rgipc.spectra.k_{rgipc_k}.dt_{rgipc_dt}.consensus.txt', index_col = 0).T
rgipc_weights.columns = [f'F{i+1}' for i in range(rgipc_k)]
n_overlapping_genes = rgipc_weights.index.intersection(cnmf_weights.index).size
print(f'Number of overlapping genes between RG-IPC and main dataset: {n_overlapping_genes}')

ipcenn_k = 11; ipcenn_dt = '0_1'; ipcenn_factor_order = [f'F{i}' for i in [11,4,9,10,6,8,1,2,3,5,7]]
ipcenn_weights = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/wang_2025/wang_2025_ipcenn/' +
    f'wang_2025_ipcenn.spectra.k_{ipcenn_k}.dt_{ipcenn_dt}.consensus.txt', index_col = 0).T
ipcenn_weights.columns = [f'F{i+1}' for i in range(ipcenn_k)]
n_overlapping_genes = ipcenn_weights.index.intersection(cnmf_weights.index).size
print(f'Number of overlapping genes between IPC-EN and main dataset: {n_overlapping_genes}')

sensitivity_k = [9, 10, 16, 36]; sensitivity_dt = ['0_1','0_1','0_15', '0_2']
sensitivity_weights = [
    pd.read_table(
        f'/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/wang_2025/wang_2025_neocortex/' +
        f'wang_2025_neocortex.spectra.k_{k}.dt_{dt}.consensus.txt', index_col = 0).T.rename(columns = {i: f'F{i}' for i in range(1, k+1)})
    for k, dt in zip(sensitivity_k, sensitivity_dt)
]

for test_weights,k in zip(sensitivity_weights, sensitivity_k):
    linkage_matrix = linkage(test_weights.T, method = 'average', metric = 'correlation')
    dendrogram_res = dendrogram(linkage_matrix, no_plot = True)
    corr = pd.concat([test_weights.corrwith(rgipc_weights[f'F{i+1}']).to_frame(name = f'F{i+1}') for i in range(rgipc_k)], axis = 1)
    corr = corr.loc[:, rgipc_factor_order].iloc[dendrogram_res['leaves'], :]
    out_prefix = f'{output_dir}/rgipc_corr_k{rgipc_k}_to_k{k}'
    fig = plt.figure(figsize = (8,8))
    sns.heatmap(corr, cmap = 'vlag', center = 0, yticklabels = True, xticklabels = True, square = True)
    fig.savefig(f'{out_prefix}.pdf', bbox_inches = 'tight')
    plt.close()

    corr = pd.concat([test_weights.corrwith(ipcenn_weights[f'F{i+1}']).to_frame(name = f'F{i+1}') for i in range(ipcenn_k)], axis = 1)
    corr = corr.loc[:, ipcenn_factor_order].iloc[dendrogram_res['leaves'], :]
    out_prefix = f'{output_dir}/ipcenn_corr_k{ipcenn_k}_to_k{k}'
    fig = plt.figure(figsize = (8,8))
    sns.heatmap(corr, cmap = 'vlag', center = 0, yticklabels = True, xticklabels = True, square = True)
    fig.savefig(f'{out_prefix}.pdf', bbox_inches = 'tight')
    plt.close()

replication_k = 10; replication_dt = '0_1'
replication_weights = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/polioudakis_2019/neocx_wang_2025_genes/'+
    f'neocx_wang_2025_genes.spectra.k_{replication_k}.dt_{replication_dt}.consensus.txt',
    index_col = 0
).T
n_overlapping_genes = replication_weights.index.intersection(cnmf_weights.index).size
print(f'Number of overlapping genes between Polioudakis et al. 2019 and main dataset: {n_overlapping_genes}')

linkage_matrix = linkage(replication_weights.T, method = 'average', metric = 'correlation')
dendrogram_res = dendrogram(linkage_matrix, no_plot = True)
corr = pd.concat([replication_weights.corrwith(cnmf_weights[f'F{i+1}']).to_frame(name = f'F{i+1}') for i in range(cnmf_k)], axis = 1)
corr = corr.loc[:, cnmf_factor_order].iloc[dendrogram_res['leaves'], :]
out_prefix = f'{output_dir}/cnmf_corr_k{cnmf_k}_to_polioudakis_k{replication_k}'
fig = plt.figure(figsize = (10,10))
sns.heatmap(corr, cmap = 'vlag', center = 0, yticklabels = True, xticklabels = True, square = True)
fig.savefig(f'{out_prefix}.pdf', bbox_inches = 'tight')
plt.close()
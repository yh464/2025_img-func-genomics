import pandas as pd
from scipy.stats import ttest_ind
import numpy as np
import scanpy as sc
import os

if not os.path.isfile('/rds/project/rds-Nl99R8pHODQ/multiomics/scdrs/wang_2025/transient_state_score.txt') or \
    not os.path.isfile('/rds/project/rds-Nl99R8pHODQ/multiomics/scdrs/wang_2025/transient_states.txt'):
    adata = sc.read_h5ad('/rds/project/rds-Nl99R8pHODQ/multiomics/scdrs/wang_2025/wang_2025_neocortex.h5ad')
    ts_score = pd.DataFrame(index = adata.obs_names, columns = ['vRG','oRG','tRG','IPC','EN','RG','S','G2M','M','G1S','MG1'])
    phase_markers = pd.read_table('/rds/project/rds-Nl99R8pHODQ/multiomics/scdrs/wang_2025/macosco_2015_cell_cycle_genes.txt')
    for phase in ['S','G2/M','M','G1/S','M/G1']:
        markers = [gene for gene in phase_markers[phase] if gene in adata.var_names]
        ts_score[phase.replace('/','')] = adata[:,markers].X.mean(axis = 1)

    cell_type_markers = pd.read_table('/rds/project/rds-Nl99R8pHODQ/multiomics/scdrs/wang_2025/marker_genes.txt')
    for abbr, orig in zip(
        ['vRG','oRG','tRG','IPC','EN'],
        ['RG-vRG', 'RG-oRG', 'RG-tRG', 'IPC-EN', 'EN-Newborn']
        ):
        markers = [gene for gene in cell_type_markers.loc[cell_type_markers['Cell type'] == orig, 'ENSG'] if gene in adata.var_names]
        ts_score[abbr] = adata[:,markers].X.mean(axis = 1)
    markers = [gene for gene in cell_type_markers.loc[cell_type_markers['Cell type'].str.startswith('RG'), 'ENSG'] if gene in adata.var_names]
    markers = list(set(markers))
    ts_score['RG'] = adata[:,markers].X.mean(axis = 1)
    ts_score.to_csv('/rds/project/rds-Nl99R8pHODQ/multiomics/scdrs/wang_2025/transient_state_score.txt', sep = '\t')

    ts_orig = pd.DataFrame(index = adata.obs_names, columns = ['S','G2M','M'])
    ts_score['phase_max'] = ts_score[['S','G2M','M','MG1','G1S']].idxmax(axis = 1)
    for phase in ['S','G2M','M']:
        ts_orig[phase] = ts_score['phase_max'] == phase
    ts_orig = ts_orig.loc[ts_orig.any(axis = 1),:]
    ts_orig.to_csv('/rds/project/rds-Nl99R8pHODQ/multiomics/scdrs/wang_2025/transient_states.txt', sep = '\t')

ts_orig = pd.read_table('/rds/project/rds-Nl99R8pHODQ/multiomics/scdrs/wang_2025/transient_states.txt', index_col = 0)
ts_score =  pd.read_table('/rds/project/rds-Nl99R8pHODQ/multiomics/scdrs/wang_2025/transient_state_score.txt', index_col = 0)

def threshold_ts(thr, ts = ts_orig):
    threshold = np.log(100) + thr
    ts = ts.copy()
    ts.loc[:,['RG','vRG','oRG','tRG','IPC','EN']] = ts_score.loc[:,['RG','vRG','oRG','tRG','IPC','EN']] >= threshold
    for col in ['RG','vRG','oRG','tRG']:
        ts[f'{col}EN'] = ts[col] & ts.EN
        ts[f'{col}IPC'] = ts[col] & ts.IPC & ~ts.EN
        ts[f'{col}_'] = ts[col] & ~ts.EN
    ts['IPCEN'] = ts.IPC & ts.EN & ~ts.RG
    ts['IPC_'] = ts.IPC & ~ts.RG & ~ts.EN
    ts['differentiating'] = ts.RGEN | ts.IPCEN
    ts['amplifying'] = ts.RG_ | ts.IPC_
    return ts

def scdrs_assoc_test(df, group_filter = True):
    df = df.loc[group_filter, :]
    ctrl_columns = [col for col in df.columns if col.startswith('ctrl_norm_score')]
    n_ctrl = len(ctrl_columns)
    score_q95 = np.quantile(df['norm_score'], 0.95)
    v_ctrl_score_q95 = np.quantile(df[ctrl_columns], 0.95, axis = 0)
    mc_p = ((v_ctrl_score_q95 >= score_q95).sum() + 1) / (n_ctrl + 1)
    mc_z = (score_q95 - v_ctrl_score_q95.mean()) / v_ctrl_score_q95.std(ddof = 1)
    return mc_p, mc_z


def test(input_file):
    df = pd.read_table(input_file, index_col = 0)
    
    for thr in [0, 0.1, 0.2, 0.4, 0.5, 0.7, 1]:
        ts = threshold_ts(thr, ts_orig)
        tests = pd.DataFrame(
            index = [
                'all_precursors','all_S','all_G2M','all_M',
                'RG','RG_S','RG_G2M','RG_M',
                'vRG','vRG_S','vRG_G2M','vRG_M',
                'oRG','oRG_S','oRG_G2M','oRG_M',
                'tRG','tRG_S','tRG_G2M','tRG_M',
                'IPC','IPC_S','IPC_G2M','IPC_M'
                ],
            columns = ['diff','t','p','df', 'n_amp', 'n_diff', 
                       'amplifying_mcp','amplifying_mcz', 'differentiating_mcp','differentiating_mcz']
            )
        
        type_filters = [
            (ts.amplifying, ts.differentiating),
        #    (ts.RG, ts.RGEN | ts.RGIPC),
            (ts.RG_, ts.RGEN),
            (ts.vRG_, ts.vRGEN),
            (ts.oRG_, ts.oRGEN),
            (ts.tRG_, ts.tRGEN),
            (ts.IPC_, ts.IPCEN)
            ]
        
        phase_filters = [
            True, ts.S, ts.G2M, ts.M
            ]
        
        idx = 0
        for tf in type_filters:
            for pf in phase_filters:
                res = ttest_ind(
                    df.loc[ts.index[tf[0] & pf], 'norm_score'],
                    df.loc[ts.index[tf[1] & pf], 'norm_score']
                    )
                tests.iloc[idx, 0] = df.loc[ts.index[tf[0] & pf], 'norm_score'].mean() - df.loc[ts.index[tf[1] & pf], 'norm_score'].mean()
                tests.iloc[idx, 1:4] = [res.statistic, res.pvalue, res.df]
                tests.iloc[idx, 4] = (tf[0]&pf).sum()
                tests.iloc[idx, 5] = (tf[1]&pf).sum()
                amp_mcp, amp_mcz = scdrs_assoc_test(df.loc[ts.index[tf[0] & pf], :])
                diff_mcp, diff_mcz = scdrs_assoc_test(df.loc[ts.index[tf[1] & pf], :])
                tests.iloc[idx, 6:10] = [amp_mcp, amp_mcz, diff_mcp, diff_mcz]
                idx += 1
        print(tests)

input_files = ['/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/sc/scdrs/disorders/adhd2025/wang_2025/adhd2025.wang_2025_neocortex.scdrs.score.txt',
               '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/sc/scdrs/structural_factors/cortical_expansion/wang_2025/cortical_expansion.wang_2025_neocortex.scdrs.score.txt']

for file in input_files: test(file)



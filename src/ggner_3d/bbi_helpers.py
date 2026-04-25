import bbi
import numpy as np
import matplotlib.pyplot as plt

def stackup_bw_dict(bw_dict, df, nbins=150, summary='mean'):
    bw_dict_rep_avg = {}
    for k, bw_list in bw_dict.items():
        bw_dict_rep_avg[k] = []
        for bw_path in bw_list:
            current_stackup = bbi.stackup(
                bw_path,
                df.chrom,
                df.flank_start,
                df.flank_end,
                bins=nbins,
                summary=summary
            )
            bw_dict_rep_avg[k].append(current_stackup)
    bw_dict_rep_avg = {k: np.nanmean(np.array(v), axis=0) for k, v in bw_dict_rep_avg.items()}
    return bw_dict_rep_avg

def stackup_bw_dict_q(bw_list, df, nbins=150, q_col='Q', mean_per_q=True, summary='mean'):
    _rep_avg = []
    for bw_path in bw_list:
        current_stackup = bbi.stackup(
            bw_path,
            df.chrom,
            df.flank_start,
            df.flank_end,
            bins=nbins,
            summary=summary
        )
        _rep_avg.append(current_stackup)
    _rep_avg = np.nanmean(np.array(_rep_avg), axis=0)
    _rep_avg_q = {}
    if mean_per_q:
        for q in sorted(df[q_col].unique()):
            _rep_avg_q[q] = np.nanmean(_rep_avg[df[q_col] == q], axis=0)
    else:
        for q in sorted(df[q_col].unique()):
            _rep_avg_q[q] = _rep_avg[df[q_col] == q]
    return _rep_avg_q
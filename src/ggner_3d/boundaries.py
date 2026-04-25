import numpy as np
import pandas as pd
from typing import Tuple, Dict
from functools import partial
from coolpuppy import coolpup
import pandas as pd
import numpy as np
from coolpuppy.lib.numutils import get_domain_score
from more_itertools import collapse
import cooler
from typing import Callable, Tuple


def make_tads(insul_df, minlen=100_000, maxlen=1_500_000):
    df = insul_df.sort_values(["chrom", "start", "end"]).reset_index(drop=True).copy()

    next_df = (
        df.groupby("chrom", sort=False)[["start", "end"]]
        .shift(-1)
        .rename(columns={"start": "start2", "end": "end2"})
    )

    tads = pd.concat([df, next_df], axis=1)
    tads = tads[tads["start2"].notna()].copy()

    lengths = tads["start2"] - tads["start"]
    tads = tads[(lengths >= minlen) & (lengths <= maxlen)].copy()

    tads["start"] = (tads["start"] + tads["end"]) // 2
    tads["end"] = (tads["start2"] + tads["end2"]) // 2

    return (
        tads[["chrom", "start", "end"]]
        .astype({"start": "int64", "end": "int64"})
        .reset_index(drop=True)
    )

def accumulate_values_safe(dict1, dict2, key):
    # if the incoming pup doesn't have the key, just keep what we have
    if key not in dict2:
        return dict1

    v2 = dict2[key]
    if key in dict1:
        dict1[key] = list(collapse([dict1[key], v2]))
    else:
        dict1[key] = list(collapse([[v2]]))
    return dict1

def extra_sum_func(dict1, dict2):
    return accumulate_values_safe(dict1, dict2, 'domain_score')

def add_domain_score(snippet, flank=1):
    snippet['domain_score'] = [  # list-like on purpose
        f"{get_domain_score(snippet['data'], flank=flank)}_"
        f"{snippet['chrom1']}_{snippet['start1']}_{snippet['end1']}_"
        f"{snippet['chrom2']}_{snippet['start2']}_{snippet['end2']}"
    ]
    return snippet

def postprocess_str(scores):
    df = {
        'chrom': [],
        'start': [],
        'end': [],
        'domain_score': [],
    }
    for score in scores:
        domain_score, chrom1, start1, end1, _, _, _ = score.split('_')
        df['chrom'].append(chrom1)
        df['start'].append(int(start1))
        df['end'].append(int(end1))
        df['domain_score'].append(float(domain_score))
    return pd.DataFrame(df)
    
def assign_domain_score(
        tads_df: pd.DataFrame,
        clr: cooler.Cooler,
        expected: pd.DataFrame,
        view_df: pd.DataFrame,
        resolution: int = 10_000,
        nproc: int = 6,
        clr_weight_name: str = 'sweight',
        ignore_diags=0,
        rescale_flank: float = 1,
        rescale_size: int = 99,
        score_flank: int = 1,
):
    cc = coolpup.CoordCreator(tads_df, resolution=resolution, features_format='bed', local=True, rescale_flank=rescale_flank)
    pu = coolpup.PileUpper(clr, cc, expected=expected, view_df=view_df, ignore_diags=ignore_diags, rescale_size=rescale_size, rescale=True, nproc=nproc, clr_weight_name=clr_weight_name)
    postprocess_func = partial(add_domain_score, flank=score_flank)
    pup = pu.pileupsWithControl(postprocess_func=postprocess_func, extra_sum_funcs={'domain_score': extra_sum_func})
    return postprocess_str(pup['domain_score'][0])

from typing import Iterable, Optional, Tuple, Dict, Any, Literal, Union
Method = Literal["otsu", "li"]
def multi_upper_threshold(
    samples: Iterable[np.ndarray],
    *,
    method: Method = "otsu",
    bins: int = 256,
    value_range: Optional[Tuple[float, float]] = None,
    weights: Optional[np.ndarray] = None,
    # Li-specific
    max_iter: int = 200,
    tol: float = 1e-6,
    shift_if_needed: bool = True,
    eps: float = 1e-12,
    return_debug: bool = False,
) -> Union[float, Tuple[float, Dict[str, Any]]]:
    """
    Find a single global *upper* threshold from N continuous distributions.

    Supports:
      - method="otsu": maximize the *average* Otsu between-class variance objective.
      - method="li":   solve a *global* Li minimum cross-entropy fixed-point by averaging
                      the per-distribution Li updates.

    Parameters
    ----------
    samples:
        Iterable of 1D arrays; each array is one distribution's samples.
    method:
        "otsu" or "li".
    bins:
        Histogram bins for Otsu (continuous -> discretized). Ignored by Li.
    value_range:
        (min, max) for Otsu histograms. If None, inferred from all samples.
    weights:
        Optional shape (N,) weights for distributions (normalized internally).
    max_iter, tol:
        Iteration controls for Li.
    shift_if_needed:
        Li uses logs of class means; if data may be <= 0, we shift all samples by a constant
        to make them positive, then shift threshold back.
    eps:
        Small constant to avoid division-by-zero / log(0).
    return_debug:
        If True, returns (threshold, debug_dict).

    Returns
    -------
    threshold : float
        Global threshold.
    debug : dict (optional)
    """
    xs = [np.asarray(x, dtype=float).ravel() for x in samples]
    if len(xs) == 0:
        raise ValueError("samples must contain at least one array.")
    if any(x.size == 0 for x in xs):
        raise ValueError("All sample arrays must be non-empty.")
    if any(~np.isfinite(x).any() for x in xs):
        # Not strictly required, but avoids silent issues
        xs = [x[np.isfinite(x)] for x in xs]
        if any(x.size == 0 for x in xs):
            raise ValueError("After removing non-finite values, at least one array became empty.")

    N = len(xs)
    if weights is None:
        w_dist = np.ones(N, dtype=float) / N
    else:
        w_dist = np.asarray(weights, dtype=float).ravel()
        if w_dist.shape[0] != N:
            raise ValueError(f"weights must have length {N}")
        s = w_dist.sum()
        if s <= 0:
            raise ValueError("weights must sum to a positive value")
        w_dist = w_dist / s

    method = method.lower().strip()
    if method not in ("otsu", "li"):
        raise ValueError("method must be 'otsu' or 'li'")

    if method == "otsu":
        # Determine common histogram range
        if value_range is None:
            lo = min(np.min(x) for x in xs)
            hi = max(np.max(x) for x in xs)
        else:
            lo, hi = map(float, value_range)

        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            raise ValueError(f"Invalid value_range inferred/provided: ({lo}, {hi})")

        edges = np.linspace(lo, hi, bins + 1)
        centers = (edges[:-1] + edges[1:]) / 2.0

        # Thresholds are at edges[1:] up to edges[-2] boundaries -> bins-1 candidates
        J_total = np.zeros(bins - 1, dtype=float)

        for wd, x in zip(w_dist, xs):
            hist, _ = np.histogram(x, bins=edges)
            p = hist.astype(float)
            p_sum = p.sum()
            if p_sum <= 0:
                continue
            p /= p_sum

            w0 = np.cumsum(p)          # class0 probability up to bin k
            m = np.cumsum(p * centers) # cumulative mean numerator
            muT = m[-1]                # total mean

            w0_k = w0[:-1]
            w1_k = 1.0 - w0_k

            valid = (w0_k > eps) & (w1_k > eps)

            mu0 = np.zeros_like(w0_k)
            mu1 = np.zeros_like(w0_k)
            mu0[valid] = m[:-1][valid] / w0_k[valid]
            mu1[valid] = (muT - m[:-1][valid]) / w1_k[valid]

            J_i = np.zeros_like(w0_k)
            J_i[valid] = w0_k[valid] * w1_k[valid] * (mu0[valid] - mu1[valid]) ** 2

            J_total += wd * J_i

        k_star = int(np.nanargmax(J_total))
        T_star = float(edges[k_star + 1])  # boundary between bins

        if not return_debug:
            return T_star
        dbg = {
            "method": "otsu",
            "threshold": T_star,
            "k_star": k_star,
            "objective_avg": J_total,
            "objective_avg_max": float(J_total[k_star]),
            "bin_edges": edges,
            "bin_centers": centers,
            "range": (float(lo), float(hi)),
            "bins": bins,
            "dist_weights": w_dist,
        }
        return T_star, dbg

    # -------------------
    # method == "li"
    # -------------------
    # Li update per distribution:
    #   t_i = (mu0 - mu1) / (log(mu0) - log(mu1))
    # where mu0 is mean of values <= T and mu1 mean of values > T.
    #
    # We solve for a *global* T by averaging per-distribution updates:
    #   T_{new} = sum_i w_i * t_i(T)
    #
    # Note: Li requires mu0, mu1 > 0 due to logs. We'll shift if needed.
    all_vals = np.concatenate(xs)
    shift = 0.0

    if shift_if_needed:
        min_val = float(np.min(all_vals))
        if min_val <= 0:
            # shift so that min becomes eps
            shift = (-min_val) + eps
            xs_li = [x + shift for x in xs]
        else:
            xs_li = xs
    else:
        xs_li = xs

    # Initialize global threshold (weighted mean works well)
    T = float(np.sum([wd * np.mean(x) for wd, x in zip(w_dist, xs_li)]))

    invalid_counts = 0
    for it in range(max_iter):
        t_updates = []
        w_used = []
        per_dist_stats = []

        for wd, x in zip(w_dist, xs_li):
            # Split by current threshold
            m0 = x <= T
            m1 = ~m0
            if not m0.any() or not m1.any():
                # can't compute update if one side empty
                invalid_counts += 1
                continue

            mu0 = float(np.mean(x[m0]))
            mu1 = float(np.mean(x[m1]))

            # Ensure positivity for logs (guard against zeros due to degenerate data)
            mu0p = max(mu0, eps)
            mu1p = max(mu1, eps)
            denom = (np.log(mu0p) - np.log(mu1p))

            if abs(denom) < eps:
                invalid_counts += 1
                continue

            t_i = (mu0p - mu1p) / denom
            if not np.isfinite(t_i):
                invalid_counts += 1
                continue

            t_updates.append(float(t_i))
            w_used.append(float(wd))
            per_dist_stats.append((mu0, mu1))

        if len(t_updates) == 0:
            raise RuntimeError(
                "Li failed: could not compute any valid per-distribution updates. "
                "Try shift_if_needed=True (default) or check data degeneracy."
            )

        w_used = np.asarray(w_used, dtype=float)
        w_used /= w_used.sum()

        T_new = float(np.sum(w_used * np.asarray(t_updates, dtype=float)))

        if abs(T_new - T) <= tol * max(1.0, abs(T)):
            T = T_new
            break
        T = T_new

    # shift threshold back to original scale
    T_star = float(T - shift)

    if not return_debug:
        return T_star

    dbg = {
        "method": "li",
        "threshold": T_star,
        "threshold_internal_shifted": float(T),
        "shift_added": float(shift),
        "iterations": it + 1,
        "max_iter": max_iter,
        "tol": tol,
        "invalid_update_count": int(invalid_counts),
        "dist_weights": w_dist,
    }
    return T_star, dbg

def assign_e1_to_tad_bins(
    ev: pd.DataFrame,
    tads: pd.DataFrame,
    ev_res: int = 100_000,
    bin_res: int = 10_000,
    value_col: str = "E1",
    chrom_col: str = "chrom",
    start_col: str = "start",
    end_col: str = "end",
    tad_id_col: str = "tad_id",
    agg: Callable = np.mean,   # e.g. np.mean, np.median, np.nanmean, np.nanmedian
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Tile each TAD into fixed-width bins (default 10kb) and assign each bin an
    overlap-weighted average of an EV track (default 100kb). Then aggregate
    per-TAD using `agg` across bin-level values.

    Notes
    -----
    - Bin-level assignment is overlap-weighted across up to two EV bins.
    - Bin-level NaNs happen when there is no EV value covering the bin.
    - Per-TAD aggregation uses `agg` on the array of bin values for that TAD.
      If you want NaNs ignored, pass np.nanmean / np.nanmedian (or your own).
    - Assumes 0-based half-open intervals: [start, end).

    Returns
    -------
    bins_df : DataFrame
        Columns: tad_id, chrom, start, end, E1_bin
    tad_df : DataFrame
        Original tads (with tad_id) + per-TAD aggregated value column: E1_agg
    """
    # --- checks ---
    for df, name in [(ev, "ev"), (tads, "tads")]:
        for c in [chrom_col, start_col, end_col]:
            if c not in df.columns:
                raise ValueError(f"{name} is missing required column '{c}'")
    if value_col not in ev.columns:
        raise ValueError(f"ev is missing required column '{value_col}'")

    ev = ev.copy()
    tads = tads.copy().reset_index(drop=True)

    if tad_id_col not in tads.columns:
        tads[tad_id_col] = np.arange(len(tads), dtype=np.int64)

    # --- 1) build bins inside each TAD ---
    starts = tads[start_col].to_numpy(np.int64)
    ends   = tads[end_col].to_numpy(np.int64)
    lens   = ends - starts
    if (lens <= 0).any():
        bad = np.where(lens <= 0)[0][:10]
        raise ValueError(f"Found non-positive TAD lengths at indices {bad.tolist()}")

    n_bins = (lens + bin_res - 1) // bin_res  # ceil

    rep_idx = np.repeat(tads.index.to_numpy(np.int64), n_bins)
    bins_df = tads.loc[rep_idx, [tad_id_col, chrom_col, start_col, end_col]].copy()

    k = np.concatenate([np.arange(n, dtype=np.int64) for n in n_bins])
    bs = bins_df[start_col].to_numpy(np.int64) + k * bin_res
    be = np.minimum(bs + bin_res, bins_df[end_col].to_numpy(np.int64))

    bins_df[start_col] = bs
    bins_df[end_col]   = be
    bins_df = bins_df[[tad_id_col, chrom_col, start_col, end_col]]

    # --- 2) overlap-weighted value from EV ---
    ev["ev_id"] = (ev[start_col].to_numpy(np.int64) // ev_res).astype(np.int64)
    ev_s = ev.set_index([chrom_col, "ev_id"])[value_col]

    left_id  = (bins_df[start_col].to_numpy(np.int64) // ev_res).astype(np.int64)
    right_id = ((bins_df[end_col].to_numpy(np.int64) - 1) // ev_res).astype(np.int64)

    mi_left  = pd.MultiIndex.from_arrays([bins_df[chrom_col].to_numpy(), left_id])
    mi_right = pd.MultiIndex.from_arrays([bins_df[chrom_col].to_numpy(), right_id])

    left_val  = ev_s.reindex(mi_left).to_numpy(dtype=float)
    right_val = ev_s.reindex(mi_right).to_numpy(dtype=float)

    bin_start = bins_df[start_col].to_numpy(np.int64)
    bin_end   = bins_df[end_col].to_numpy(np.int64)

    left_end_boundary    = (left_id + 1) * ev_res
    right_start_boundary = right_id * ev_res

    w_left  = (np.minimum(bin_end, left_end_boundary) - bin_start).astype(float)
    w_right = (bin_end - np.maximum(bin_start, right_start_boundary)).astype(float)

    same = (left_id == right_id)
    w_right[same] = 0.0
    w_left[same]  = (bin_end[same] - bin_start[same]).astype(float)

    num = np.zeros(len(bins_df), dtype=float)
    den = np.zeros(len(bins_df), dtype=float)

    m = ~np.isnan(left_val)
    num[m] += left_val[m] * w_left[m]
    den[m] += w_left[m]

    m = ~np.isnan(right_val)
    num[m] += right_val[m] * w_right[m]
    den[m] += w_right[m]

    bins_df[f"{value_col}_bin"] = np.where(den > 0, num / den, np.nan)

    # --- 3) per-TAD aggregation using agg ---
    def _apply_agg(x: pd.Series):
        arr = x.to_numpy(dtype=float)
        return float(agg(arr))

    tad_agg = (
        bins_df.groupby(tad_id_col)[f"{value_col}_bin"]
        .apply(_apply_agg)
        .reset_index()
        .rename(columns={f"{value_col}_bin": f"{value_col}_agg"})
    )

    tad_df = tads.merge(tad_agg, on=tad_id_col, how="left")
    return bins_df, tad_df

from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


def _expand_indices(idxs: Sequence[int], accepted_range: int) -> Tuple[List[int], List[int]]:
    expanded: List[int] = []
    expanded_original: List[int] = []
    for idx in idxs:
        rng = list(range(int(idx) - accepted_range, int(idx) + accepted_range + 1))
        expanded.extend(rng)
        expanded_original.extend([int(idx)] * len(rng))
    return expanded, expanded_original


def preserved_boundaries_only(
    comparisons: Sequence[Tuple[int, int]],
    insulation_table_: Sequence[pd.DataFrame],
    sample_names_str: Sequence[str],
    window: Any,
    accepted_range: int = 1,
    comparisons_str: Optional[Sequence[str]] = None,
    verbose: bool = True,
) -> Tuple[List[Dict[str, List[Tuple[int, int]]]], List[pd.DataFrame]]:
    """
    Only computes PRESERVED boundaries and splits them into Increase vs Decrease.

    Returns
    -------
    preserved_results : list of dict
        Per comparison:
          {'Increase': [(idx1, idx2), ...], 'Decrease': [(idx1, idx2), ...]}

    preserved_dfs : list of DataFrames
        Per comparison tidy table with columns:
          sample1, sample2, dir, idx1, idx2
    """
    preserved_results: List[Dict[str, List[Tuple[int, int]]]] = []
    preserved_dfs: List[pd.DataFrame] = []

    for comp_idx, (i, j) in enumerate(comparisons):
        ins_df1 = insulation_table_[i].copy()
        ins_df2 = insulation_table_[j].copy()

        b1 = ins_df1.loc[ins_df1[f"is_boundary_{window}"] == True].index.to_numpy().astype(int)
        b2 = ins_df2.loc[ins_df2[f"is_boundary_{window}"] == True].index.to_numpy().astype(int)

        b1_exp, b1_exp_orig = _expand_indices(b1, accepted_range)

        preserved_pairs: List[Tuple[int, int]] = []
        for idx_2 in b2:
            idx_2 = int(idx_2)
            if idx_2 in b1_exp:
                k = b1_exp.index(idx_2)  # keeps your original "first match" behavior
                idx_1 = int(b1_exp_orig[k])
                preserved_pairs.append((idx_1, idx_2))

        res = {"Increase": [], "Decrease": []}

        for idx_1, idx_2 in preserved_pairs:
            bs1 = ins_df1.iloc[idx_1][f"boundary_strength_{window}"]
            bs2 = ins_df2.iloc[idx_2][f"boundary_strength_{window}"]
            if bs2 > bs1:
                res["Increase"].append((idx_1, idx_2))
            else:
                res["Decrease"].append((idx_1, idx_2))

        preserved_results.append(res)

        # dataframe for this comparison
        rows = []
        for idx1, idx2 in res["Increase"]:
            rows.append(
                dict(
                    sample1=sample_names_str[i],
                    sample2=sample_names_str[j],
                    dir="Increase",
                    idx1=int(idx1),
                    idx2=int(idx2),
                )
            )
        for idx1, idx2 in res["Decrease"]:
            rows.append(
                dict(
                    sample1=sample_names_str[i],
                    sample2=sample_names_str[j],
                    dir="Decrease",
                    idx1=int(idx1),
                    idx2=int(idx2),
                )
            )

        df_now = pd.DataFrame(rows, columns=["sample1", "sample2", "dir", "idx1", "idx2"])
        preserved_dfs.append(df_now)

    if verbose:
        for comp_idx, (i, j) in enumerate(comparisons):
            tag = (
                comparisons_str[comp_idx]
                if comparisons_str is not None
                else f"{sample_names_str[i]} vs {sample_names_str[j]}"
            )
            inc_n = len(preserved_results[comp_idx]["Increase"])
            dec_n = len(preserved_results[comp_idx]["Decrease"])
            print(f"{tag}: Preserved / Increase ::: {inc_n}")
            print(f"{tag}: Preserved / Decrease ::: {dec_n}")

    return preserved_results, preserved_dfs

def get_insulation_fixed(insul_df, coord_df, insul_column, flank=100_000, resolution=10_000):
    
    insul_df = insul_df.copy()
    coord_df = coord_df.copy()

    insul_df['mid'] = (insul_df['start'] + insul_df['end']) // 2
    coord_df['mid'] = (coord_df['start'] + coord_df['end']) // 2

    bins_each_side = flank // resolution
    window_size = 2 * bins_each_side + 1

    results = []

    for chrom in coord_df['chrom'].unique():
        insul_chrom = insul_df[insul_df['chrom'] == chrom].reset_index(drop=True)
        coord_chrom = coord_df[coord_df['chrom'] == chrom]

        if insul_chrom.empty or coord_chrom.empty:
            continue

        insul_mids = insul_chrom['mid'].values
        insul_scores = insul_chrom[insul_column].values

        for _, row in coord_chrom.iterrows():
            mid = row['mid']

            # find closest insulation bin index
            idx = np.searchsorted(insul_mids, mid)

            start = idx - bins_each_side
            end = idx + bins_each_side + 1

            if start < 0 or end > len(insul_scores):
                continue

            results.append(insul_scores[start:end])

    return np.array(results)


def genomewide_full_tad2tad_from_pixel_field(
    clr: cooler.Cooler,
    tads_df: pd.DataFrame,
    *,
    field: str = "scaled",      # pixel field to aggregate (e.g., "scaled" or "oe")
    agg: str = "mean",          # "mean" or "sum" within each TAD×TAD block
    ignore_diagonal: int = 0,   # drop pixels with dist_bins < ignore_diagonal
    symmetrize: bool = True,    # mirror upper triangle to lower
    return_per_chrom: bool = True,  # if True, also return per-chrom matrices
):
    """
    Build FULL cis TAD×TAD matrices from a cooler pixel field, and optionally a single
    genome-wide block-diagonal matrix.

    This is appropriate when your cooler already contains the quantity you want to aggregate
    (e.g., an O/E field 'oe' or a scaled balanced field 'scaled').

    Notes
    -----
    - This is cis-only by construction (we fetch per chromosome).
    - 'sum' means sum of the *field values* inside each TAD block.
      For scaled balanced fields this is often meaningful; for OE, 'mean' is usually preferred.

    Returns
    -------
    G : np.ndarray
        Genome-wide block-diagonal TAD×TAD matrix (total_TADs × total_TADs)
    meta : dict
        chrom_order, offsets, sizes, tads_by_chrom
    per_chrom : dict (optional)
        per_chrom[chrom] = {"matrix": M, "tads": tads_chr}
    """

    if agg not in {"mean", "sum"}:
        raise ValueError("agg must be 'mean' or 'sum'")
    if ignore_diagonal < 0:
        raise ValueError("ignore_diagonal must be >= 0")

    bins = clr.bins()[:]  # indexed by global bin_id
    chroms = list(bins["chrom"].unique())

    tads = tads_df[["chrom", "start", "end"]].copy()
    tads = tads.sort_values(["chrom", "start", "end"]).reset_index(drop=True)

    def tad_map_for_chrom(bins_chr: pd.DataFrame, tads_chr: pd.DataFrame):
        """Map each bin midpoint to a TAD id (0..n_tads-1) or -1 if outside."""
        tads_chr = tads_chr.sort_values(["start", "end"]).reset_index(drop=True)
        mid = ((bins_chr["start"].to_numpy() + bins_chr["end"].to_numpy()) // 2).astype(np.int64)
        starts = tads_chr["start"].to_numpy(np.int64)
        ends = tads_chr["end"].to_numpy(np.int64)

        idx = np.searchsorted(starts, mid, side="right") - 1
        ok = (idx >= 0) & (mid < ends[idx])
        out = np.full(len(mid), -1, dtype=np.int32)
        out[ok] = idx[ok]
        return out, tads_chr

    per_chrom = {}
    sizes = {}
    tads_by_chrom = {}

    for chrom in chroms:
        tads_chr = tads[tads["chrom"] == chrom]
        if tads_chr.empty:
            continue

        bins_chr = bins[bins["chrom"] == chrom]
        if bins_chr.empty:
            continue

        tad_map, tads_chr = tad_map_for_chrom(bins_chr, tads_chr)
        n_tads = len(tads_chr)
        if n_tads == 0:
            continue

        # cis pixels for this chromosome
        pix = clr.matrix(field=field, as_pixels=True, join=False).fetch(chrom)
        if pix.empty:
            continue

        # upper triangle + dist (as column, stays aligned)
        pix = pix[pix["bin2_id"] >= pix["bin1_id"]].copy()
        pix["dist_bins"] = (pix["bin2_id"] - pix["bin1_id"]).astype(np.int64)

        if ignore_diagonal > 0:
            pix = pix[pix["dist_bins"] >= ignore_diagonal].copy()
            if pix.empty:
                continue

        # map global bin_id -> local index within chrom to index tad_map
        bin_ids_chr = bins_chr.index.to_numpy()
        local_index = pd.Series(np.arange(len(bin_ids_chr)), index=bin_ids_chr)

        b1_local = local_index.loc[pix["bin1_id"]].to_numpy()
        b2_local = local_index.loc[pix["bin2_id"]].to_numpy()

        pix["tad1"] = tad_map[b1_local]
        pix["tad2"] = tad_map[b2_local]
        pix = pix[(pix["tad1"] >= 0) & (pix["tad2"] >= 0)].copy()
        if pix.empty:
            continue

        v = pix[field].to_numpy()
        good = np.isfinite(v)
        if not np.any(good):
            continue

        tad1 = pix["tad1"].to_numpy()[good]
        tad2 = pix["tad2"].to_numpy()[good]
        v = v[good]

        # aggregate pixels -> FULL per-chrom TAD×TAD matrix
        g = pd.DataFrame({"tad1": tad1, "tad2": tad2, "v": v})
        if agg == "mean":
            gb = g.groupby(["tad1", "tad2"], sort=False)["v"].mean().reset_index()
        else:
            gb = g.groupby(["tad1", "tad2"], sort=False)["v"].sum().reset_index()

        M = np.full((n_tads, n_tads), np.nan, dtype=float)
        M[gb["tad1"].to_numpy(), gb["tad2"].to_numpy()] = gb["v"].to_numpy()

        if symmetrize:
            M = np.where(np.isnan(M), M.T, M)

        per_chrom[chrom] = {"matrix": M, "tads": tads_chr}
        sizes[chrom] = n_tads
        tads_by_chrom[chrom] = tads_chr

    # genome-wide block diagonal
    chrom_order = [c for c in chroms if c in per_chrom]
    total = int(sum(sizes[c] for c in chrom_order))
    G = np.full((total, total), np.nan, dtype=float)

    offsets = {}
    cur = 0
    for c in chrom_order:
        offsets[c] = cur
        n = sizes[c]
        G[cur:cur + n, cur:cur + n] = per_chrom[c]["matrix"]
        cur += n

    meta = {
        "field": field,
        "agg": agg,
        "ignore_diagonal": ignore_diagonal,
        "chrom_order": chrom_order,
        "offsets": offsets,
        "sizes": sizes,
        "tads_by_chrom": tads_by_chrom,
    }

    if return_per_chrom:
        return G, meta, per_chrom
    return G, meta


def genomewide_full_tad2tad_from_pixel_field_full(
    clr: cooler.Cooler,
    tads_df: pd.DataFrame,
    *,
    field: str = "scaled",      # pixel field to aggregate (e.g., "scaled" or "oe")
    agg: str = "mean",          # "mean" or "sum" within each TAD×TAD block
    ignore_diagonal: int = 0,   # only applies to cis: drop pixels with dist_bins < ignore_diagonal
    symmetrize: bool = True,    # mirror upper triangle to lower
    return_per_chrom: bool = True,  # return diagonal cis blocks as well
):
    """
    Build a FULL genome-wide TAD×TAD matrix from a cooler pixel field.

    Compared with the cis-only version, this includes:
      - cis TAD×TAD blocks on the diagonal
      - trans TAD×TAD blocks off the diagonal

    Notes
    -----
    - We iterate over chromosome pairs, so this avoids pulling the entire pixel table
      into memory at once.
    - `ignore_diagonal` is only meaningful for cis blocks.
    - For trans blocks, pixels are aggregated exactly the same way, just without a
      genomic-distance filter.
    - If `symmetrize=True`, only the upper triangle / upper chrom-pair blocks are
      computed and then mirrored.

    Returns
    -------
    G : np.ndarray
        Genome-wide TAD×TAD matrix (total_TADs × total_TADs), including trans blocks.
    meta : dict
        chrom_order, offsets, sizes, tads_by_chrom
    per_chrom : dict (optional)
        per_chrom[chrom] = {"matrix": cis_block, "tads": tads_chr}
    """

    if agg not in {"mean", "sum"}:
        raise ValueError("agg must be 'mean' or 'sum'")
    if ignore_diagonal < 0:
        raise ValueError("ignore_diagonal must be >= 0")

    bins = clr.bins()[:].copy()   # global bin ids are in the index
    chroms = list(bins["chrom"].unique())

    tads = tads_df[["chrom", "start", "end"]].copy()
    tads = tads.sort_values(["chrom", "start", "end"]).reset_index(drop=True)

    def tad_map_for_chrom(bins_chr: pd.DataFrame, tads_chr: pd.DataFrame):
        """Map each bin midpoint to a chromosome-local TAD id (0..n_tads-1), or -1."""
        tads_chr = tads_chr.sort_values(["start", "end"]).reset_index(drop=True)

        mid = ((bins_chr["start"].to_numpy() + bins_chr["end"].to_numpy()) // 2).astype(np.int64)
        starts = tads_chr["start"].to_numpy(np.int64)
        ends = tads_chr["end"].to_numpy(np.int64)

        idx = np.searchsorted(starts, mid, side="right") - 1
        ok = (idx >= 0) & (mid < ends[idx])

        out = np.full(len(mid), -1, dtype=np.int32)
        out[ok] = idx[ok]
        return out, tads_chr

    # ------------------------------------------------------------------
    # First pass: define genome-wide TAD indexing
    # ------------------------------------------------------------------
    chrom_order = []
    offsets = {}
    sizes = {}
    tads_by_chrom = {}

    cur = 0
    for chrom in chroms:
        tads_chr = tads[tads["chrom"] == chrom]
        bins_chr = bins[bins["chrom"] == chrom]

        if tads_chr.empty or bins_chr.empty:
            continue

        tads_chr = tads_chr.sort_values(["start", "end"]).reset_index(drop=True)
        if len(tads_chr) == 0:
            continue

        chrom_order.append(chrom)
        offsets[chrom] = cur
        sizes[chrom] = len(tads_chr)
        tads_by_chrom[chrom] = tads_chr
        cur += len(tads_chr)

    total = cur
    G = np.full((total, total), np.nan, dtype=float)

    if total == 0:
        meta = {
            "field": field,
            "agg": agg,
            "ignore_diagonal": ignore_diagonal,
            "chrom_order": chrom_order,
            "offsets": offsets,
            "sizes": sizes,
            "tads_by_chrom": tads_by_chrom,
            "global": True,
        }
        if return_per_chrom:
            return G, meta, {}
        return G, meta

    # ------------------------------------------------------------------
    # Second pass: map every global bin_id -> global TAD id
    # ------------------------------------------------------------------
    max_bin_id = int(bins.index.max())
    bin_to_global_tad = np.full(max_bin_id + 1, -1, dtype=np.int64)

    for chrom in chrom_order:
        bins_chr = bins[bins["chrom"] == chrom]
        tad_map_local, _ = tad_map_for_chrom(bins_chr, tads_by_chrom[chrom])

        bin_ids_chr = bins_chr.index.to_numpy(np.int64)
        valid = tad_map_local >= 0
        bin_to_global_tad[bin_ids_chr[valid]] = offsets[chrom] + tad_map_local[valid]

    # ------------------------------------------------------------------
    # Aggregate pixels chromosome-pair by chromosome-pair
    # ------------------------------------------------------------------
    pixel_selector = clr.matrix(field=field, as_pixels=True, join=False)

    for i, chrom1 in enumerate(chrom_order):
        for j in range(i, len(chrom_order)):
            chrom2 = chrom_order[j]

            # cis block
            if chrom1 == chrom2:
                pix = pixel_selector.fetch(chrom1)
                if pix.empty:
                    continue

                pix = pix[pix["bin2_id"] >= pix["bin1_id"]].copy()
                if pix.empty:
                    continue

                if ignore_diagonal > 0:
                    dist = (pix["bin2_id"] - pix["bin1_id"]).to_numpy(np.int64)
                    pix = pix[dist >= ignore_diagonal].copy()
                    if pix.empty:
                        continue

            # trans block
            else:
                pix = pixel_selector.fetch(chrom1, chrom2)
                if pix.empty:
                    continue
                pix = pix.copy()

            b1 = pix["bin1_id"].to_numpy(np.int64)
            b2 = pix["bin2_id"].to_numpy(np.int64)

            tad1 = bin_to_global_tad[b1]
            tad2 = bin_to_global_tad[b2]
            v = pix[field].to_numpy()

            good = (tad1 >= 0) & (tad2 >= 0) & np.isfinite(v)
            if not np.any(good):
                continue

            g = pd.DataFrame({
                "tad1": tad1[good],
                "tad2": tad2[good],
                "v": v[good],
            })

            if agg == "mean":
                gb = g.groupby(["tad1", "tad2"], sort=False)["v"].mean().reset_index()
            else:
                gb = g.groupby(["tad1", "tad2"], sort=False)["v"].sum().reset_index()

            ii = gb["tad1"].to_numpy(np.int64)
            jj = gb["tad2"].to_numpy(np.int64)
            G[ii, jj] = gb["v"].to_numpy(float)

    if symmetrize:
        G = np.where(np.isnan(G), G.T, G)

    meta = {
        "field": field,
        "agg": agg,
        "ignore_diagonal": ignore_diagonal,
        "chrom_order": chrom_order,
        "offsets": offsets,
        "sizes": sizes,
        "tads_by_chrom": tads_by_chrom,
        "global": True,
    }

    if return_per_chrom:
        per_chrom = {
            chrom: {
                "matrix": G[
                    offsets[chrom]:offsets[chrom] + sizes[chrom],
                    offsets[chrom]:offsets[chrom] + sizes[chrom],
                ],
                "tads": tads_by_chrom[chrom],
            }
            for chrom in chrom_order
        }
        return G, meta, per_chrom

    return G, meta

def oe_cool_dump(clr, view_df=None, weight='weight', nproc=1, dump_path=None, aggregate_trans=False):

    fast_expected_df = obs_over_exp_cooler.expected_full_fast(
    clr,
    view_df=view_df,
    smooth_cis=False,
    aggregate_trans=aggregate_trans,
    expected_column_name="expected",
    nproc=nproc,
    clr_weight_name=weight,
    )
    print("Expected matrix calculated, now calculating obs/exp...")
    
    results = []
    for oe_chunk in obs_over_exp_cooler.obs_over_exp_generator(
            clr,
            fast_expected_df,
            view_df=view_df,
            expected_column_name="expected",
            oe_column_name='oe',
            chunksize=1_000_000,
            clr_weight_name=weight,
        ):
        results.append(oe_chunk)
    res_df = pd.concat(results, ignore_index=True)

    bins_oe = clr.bins()[:].copy()
    _bad_mask = bins_oe["weight"].isna()
    bins_oe["weight"] = 1.
    bins_oe.loc[_bad_mask,"weight"] = np.nan

    if dump_path is not None:
        print(f"Dumping cooler to {dump_path}...")
        cooler.create_cooler(
            cool_uri = dump_path,
            bins = bins_oe,
            pixels = res_df,
            columns=["oe"],
            dtypes={"oe":np.float64},
            )

    return res_df

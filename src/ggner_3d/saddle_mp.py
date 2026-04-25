import numpy as np
import pandas as pd
import cooler

from cooltools.api import snipping
from cooltools.lib import common

try:
    from tqdm.auto import tqdm
except Exception:
    tqdm = None


def saddle_by_interval_category_cis(
    clr,
    expected_df,
    intervals_df,
    category_col="ds_diff_Q",
    categories=None,
    expected_value_col="balanced.avg",
    balance="weight",
    midpoint_dist_range=None,      # (min_bp, max_bp) on interval midpoint distance
    expected_region_col=None,      # e.g. chromosome-arm name if expected is arm-based
    view_df=None,                  # optional matching view for expected_df
    weight_mode="pixel",           # "pixel" or "pair"
    min_pixels_per_block=1,
    ignore_diags=0,                # optional diagonal mask in bins
    use_tqdm=False,
    tqdm_leave=False,
):
    """
    Faster saddle-like O/E aggregation for interval/domain annotations.

    Strategy:
      - pre-bin intervals once
      - group intervals by (expected_region, chrom)
      - load one sparse balanced matrix per region/chrom
      - load one expected selector per region
      - slice both in memory for all interval pairs
      - accumulate into category x category sum/count matrices

    Returns
    -------
    dict with keys:
      sum, count, oe, pair_count, categories
    """
    required_interval_cols = {"chrom", "start", "end", category_col}
    missing = required_interval_cols - set(intervals_df.columns)
    if missing:
        raise ValueError(f"intervals_df is missing required columns: {sorted(missing)}")

    if weight_mode not in {"pixel", "pair"}:
        raise ValueError("weight_mode must be 'pixel' or 'pair'")

    intervals = intervals_df.copy()
    intervals["chrom"] = intervals["chrom"].astype(str)

    if expected_region_col is None:
        intervals["_expected_region"] = intervals["chrom"]
        expected_region_col = "_expected_region"
    elif expected_region_col not in intervals.columns:
        raise ValueError(f"{expected_region_col!r} not found in intervals_df")

    intervals = intervals.loc[intervals[category_col].notna()].copy()

    if categories is None:
        categories = sorted(pd.unique(intervals[category_col]))
    categories = list(categories)
    cat_to_idx = {cat: i for i, cat in enumerate(categories)}
    n_cat = len(categories)

    if midpoint_dist_range is None:
        min_mid_bp, max_mid_bp = 0, np.inf
    else:
        min_mid_bp, max_mid_bp = midpoint_dist_range
        min_mid_bp = 0 if min_mid_bp is None else min_mid_bp
        max_mid_bp = np.inf if max_mid_bp is None else max_mid_bp

    # Build / validate view_df for expected selector
    if view_df is None:
        view_df = common.make_cooler_view(clr).copy()
    else:
        view_df = view_df.copy()

    if "name" not in view_df.columns:
        # whole-chrom expected case
        if view_df["chrom"].duplicated().any():
            raise ValueError("view_df must have a unique 'name' column when regions repeat per chromosome")
        view_df["name"] = view_df["chrom"].astype(str)

    view_df["chrom"] = view_df["chrom"].astype(str)
    view_df = view_df[["chrom", "start", "end", "name"]].copy()

    # Precompute expected selectors like coolpuppy does
    exp_snipper = snipping.ExpectedSnipper(
        clr,
        expected_df,
        view_df=view_df,
        expected_value_col=expected_value_col,
    )
    expected_selectors = {
        region_name: exp_snipper.select(region_name, region_name)
        for region_name in view_df["name"]
    }

    # Bin table once
    bins = clr.bins()[:][["chrom", "start", "end"]].copy()
    bins["chrom"] = bins["chrom"].astype(str)

    chrom_bin_cache = {}
    for chrom, b in bins.groupby("chrom", sort=False, observed=False):
        chrom_bin_cache[chrom] = {
            "starts": b["start"].to_numpy(),
            "ends": b["end"].to_numpy(),
        }

    def _interval_to_chr_local_bin_span(chrom, start, end):
        if chrom not in chrom_bin_cache:
            return None
        starts = chrom_bin_cache[chrom]["starts"]
        ends = chrom_bin_cache[chrom]["ends"]

        # first bin with end > start
        lo = np.searchsorted(ends, start, side="right")
        # first bin with start >= end
        hi = np.searchsorted(starts, end, side="left")

        if lo >= hi:
            return None
        return lo, hi

    spans = []
    mids = []
    for row in intervals.itertuples(index=False):
        chrom = row.chrom
        start = int(row.start)
        end = int(row.end)
        spans.append(_interval_to_chr_local_bin_span(chrom, start, end))
        mids.append((start + end) / 2.0)

    intervals["_bin_span_chr"] = spans
    intervals["_mid"] = mids
    intervals = intervals.loc[intervals["_bin_span_chr"].notna()].copy()

    # Precompute region metadata
    view_df_idx = view_df.set_index("name", drop=False)

    sum_mat = np.zeros((n_cat, n_cat), dtype=float)
    count_mat = np.zeros((n_cat, n_cat), dtype=float)
    pair_count_mat = np.zeros((n_cat, n_cat), dtype=int)

    outer_groups = list(intervals.groupby([expected_region_col, "chrom"], sort=False, observed=False))
    if use_tqdm and tqdm is not None:
        outer_groups = tqdm(outer_groups, desc="region/chrom", leave=tqdm_leave)

    for (region_name, chrom), chrom_df in outer_groups:
        if region_name not in expected_selectors:
            continue
        if region_name not in view_df_idx.index:
            continue

        region_row = view_df_idx.loc[region_name]
        if region_row["chrom"] != chrom:
            # expected region and interval chromosome must match for cis
            continue

        region_tuple = (region_row["chrom"], int(region_row["start"]), int(region_row["end"]))

        # region extents in chromosome-local bin coordinates
        lo_global, hi_global = clr.extent(region_tuple)
        chrom_offset = clr.offset(chrom)
        region_lo_chr = lo_global - chrom_offset
        region_hi_chr = hi_global - chrom_offset
        region_nbins = region_hi_chr - region_lo_chr
        if region_nbins <= 0:
            continue

        # coolpuppy-style: one sparse fetch for the whole region
        bigdata = clr.matrix(balance=balance, sparse=True).fetch(region_tuple, region_tuple).tocsr()
        exp_selector = expected_selectors[region_name]

        # bad bins from balancing mask
        if balance:
            bad_bins = np.isnan(clr.bins()[balance].fetch(region_tuple).to_numpy())
        else:
            bad_bins = np.zeros(region_nbins, dtype=bool)

        chrom_df = chrom_df.sort_values("_mid").reset_index(drop=True)

        mids = chrom_df["_mid"].to_numpy()
        cats = chrom_df[category_col].map(cat_to_idx).to_numpy()

        lo_chr = np.array([x[0] for x in chrom_df["_bin_span_chr"]], dtype=int)
        hi_chr = np.array([x[1] for x in chrom_df["_bin_span_chr"]], dtype=int)

        # region-local coordinates for slicing bigdata / expected selector
        lo_reg = lo_chr - region_lo_chr
        hi_reg = hi_chr - region_lo_chr

        # keep only intervals fully inside this expected region
        inside = (lo_reg >= 0) & (hi_reg <= region_nbins)
        if not np.any(inside):
            continue

        mids = mids[inside]
        cats = cats[inside]
        lo_chr = lo_chr[inside]
        hi_chr = hi_chr[inside]
        lo_reg = lo_reg[inside]
        hi_reg = hi_reg[inside]

        n = len(mids)
        if n == 0:
            continue

        anchor_iter = range(n)
        if use_tqdm and tqdm is not None:
            anchor_iter = tqdm(anchor_iter, desc=f"{region_name}", leave=tqdm_leave)

        for a in anchor_iter:
            ia = cats[a]
            if np.isnan(ia):
                continue
            ia = int(ia)

            # use searchsorted instead of scanning all later intervals
            left_mid = mids[a] + min_mid_bp
            right_mid = mids[a] + max_mid_bp

            b0 = np.searchsorted(mids, left_mid, side="left")
            b1 = n if np.isinf(max_mid_bp) else np.searchsorted(mids, right_mid, side="right")
            b0 = max(a, b0)

            a_lo_reg, a_hi_reg = lo_reg[a], hi_reg[a]
            a_lo_chr, a_hi_chr = lo_chr[a], hi_chr[a]

            for b in range(b0, b1):
                ib = cats[b]
                if np.isnan(ib):
                    continue
                ib = int(ib)

                b_lo_reg, b_hi_reg = lo_reg[b], hi_reg[b]
                b_lo_chr, b_hi_chr = lo_chr[b], hi_chr[b]

                obs = bigdata[a_lo_reg:a_hi_reg, b_lo_reg:b_hi_reg].toarray().astype(float)
                if obs.size == 0:
                    continue

                # apply balancing NaN mask
                bad1 = bad_bins[a_lo_reg:a_hi_reg]
                bad2 = bad_bins[b_lo_reg:b_hi_reg]
                if bad1.any():
                    obs[bad1, :] = np.nan
                if bad2.any():
                    obs[:, bad2] = np.nan

                exp = np.asarray(exp_selector[a_lo_reg:a_hi_reg, b_lo_reg:b_hi_reg], dtype=float)

                # optional diagonal masking
                if ignore_diags > 0:
                    d = np.abs(
                        np.arange(a_lo_chr, a_hi_chr)[:, None]
                        - np.arange(b_lo_chr, b_hi_chr)[None, :]
                    )
                    obs[d < ignore_diags] = np.nan

                with np.errstate(divide="ignore", invalid="ignore"):
                    oe = obs / exp

                finite = np.isfinite(oe)
                n_valid = int(finite.sum())
                if n_valid < min_pixels_per_block:
                    continue

                if weight_mode == "pixel":
                    contrib_sum = float(np.nansum(oe))
                    contrib_count = float(n_valid)
                else:
                    contrib_sum = float(np.nanmean(oe))
                    contrib_count = 1.0

                sum_mat[ia, ib] += contrib_sum
                count_mat[ia, ib] += contrib_count
                pair_count_mat[ia, ib] += 1

                if ia != ib:
                    sum_mat[ib, ia] += contrib_sum
                    count_mat[ib, ia] += contrib_count
                    pair_count_mat[ib, ia] += 1

    with np.errstate(divide="ignore", invalid="ignore"):
        oe_mat = sum_mat / count_mat

    return {
        "sum": sum_mat,
        "count": count_mat,
        "oe": oe_mat,
        "pair_count": pair_count_mat,
        "categories": categories,
    }
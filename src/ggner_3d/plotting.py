# plotting.py
import seaborn as sns
import matplotlib as mpl  # better than plt for global params
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.colors import Normalize, LogNorm
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from cytoolz import merge
from functools import partial
import matplotlib.colors as mcolors

import cooltools
from cooltools.api.saddle import saddle_strength

_RC = {
    # -----------------
    # Output / saving
    # -----------------
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,

    # -----------------
    # Lines / spines
    # -----------------
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.2,

    # -----------------
    # Ticks
    # -----------------
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.minor.size": 2,
    "ytick.minor.size": 2,

    # -----------------
    # Fonts (paper)
    # -----------------
    "font.family": "sans-serif",
    "font.size": 9.0,          # base
    "axes.titlesize": 12.0,
    "axes.labelsize": 12.0,
    "xtick.labelsize": 9.0,
    "ytick.labelsize": 9.0,
    "legend.fontsize": 8.0,
    "legend.title_fontsize": 8.0,

    # -----------------
    # Padding / spacing
    # -----------------
    "axes.labelpad": 4.0,      # distance between axis labels and axes
    "axes.titlepad": 6.0,      # distance between title and plot
    "legend.borderaxespad": 0.3,
    "legend.handlelength": 2.2,
    "legend.handletextpad": 0.4,
    "legend.labelspacing": 0.3,

    # -----------------
    # Math / text rendering
    # -----------------
    "mathtext.default": "it",
}

def update_rcparams():
    sns.set_theme(context="paper", style="ticks", rc=_RC)
    mpl.rcParams.update(_RC)

def despine(ax):
    sns.despine(ax=ax, top=True, right=True, offset=2, trim=True)

COLORS = ['#465775', '#A63446', '#F5B841', '#9DBBAE']


def saddleplot(
    track,
    saddledata,
    n_bins,
    vrange=None,
    qrange=(0.0, 1.0),
    cmap="coolwarm",
    scale="log",
    vmin=0.5,
    vmax=2,
    color=None,
    title=None,
    xlabel=None,
    ylabel=None,
    clabel=None,
    fig=None,
    fig_kws=None,
    heatmap_kws=None,
    margin_kws=None,
    cbar_kws=None,
    subplot_spec=None,
):
    """
    Seaborn-based saddle plot (heatmap + marginals).
    Expects `saddledata` shape (n_bins+2, n_bins+2) OR (n_bins, n_bins).
    """
    # --- digitize track into bins (same idea as cooltools saddle workflow) ---
    track_value_col = track.columns[3]

    # NOTE: requires cooltools in your environment (same as your original)
    digitized_track, binedges = cooltools.digitize(
        track, n_bins, vrange=vrange, qrange=qrange
    )

    # mean of original track values per digitized bin
    groupmean = (
        track[track_value_col]
        .groupby(digitized_track[digitized_track.columns[3]])
        .mean()
    )

    # if qrange is used, treat axes as quantiles 0..1 with uniform bins
    if qrange is not None:
        lo, hi = qrange
        binedges = np.linspace(lo, hi, n_bins + 1)
    else:
        lo, hi = float(binedges.min()), float(binedges.max())

    # --- handle flanking outlier bins in saddledata (+2) ---
    C = np.asarray(saddledata)
    if C.shape[0] == n_bins + 2 and C.shape[1] == n_bins + 2:
        C = C[1:-1, 1:-1]
        # groupmean likely includes outlier groups too; try to drop if present
        if len(groupmean) == n_bins + 2:
            groupmean = groupmean.iloc[1:-1]

    # Ensure groupmean aligns to n_bins
    if len(groupmean) != n_bins:
        # best-effort: reindex to 0..n_bins-1 if possible
        try:
            groupmean = groupmean.reindex(range(n_bins))
        except Exception:
            groupmean = np.asarray(groupmean)[:n_bins]
    groupmean = np.asarray(groupmean)

    # --- layout (3x3) ---
    if subplot_spec is not None:
        GS = partial(GridSpecFromSubplotSpec, subplot_spec=subplot_spec)
    else:
        GS = GridSpec

    gs = GS(
        nrows=3,
        ncols=3,
        width_ratios=[0.2, 1, 0.1],
        height_ratios=[0.2, 1, 0.1],
        wspace=0.05,
        hspace=0.05,
    )

    if fig is None:
        fig_kws_default = dict(figsize=(5, 5))
        fig_kws = merge(fig_kws_default, fig_kws if fig_kws is not None else {})
        fig = plt.figure(**fig_kws)

    grid = {}

    # --- axes ---
    ax_top = fig.add_subplot(gs[1])
    ax_left = fig.add_subplot(gs[3])
    ax_hm = fig.add_subplot(gs[4])
    ax_cbar = fig.add_subplot(gs[5])

    grid["ax_margin_x"] = ax_top
    grid["ax_margin_y"] = ax_left
    grid["ax_heatmap"] = ax_hm
    grid["ax_cbar"] = ax_cbar

    # --- normalization ---
    if scale == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif scale == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("Only 'linear' and 'log' are supported for scale.")

    # --- seaborn heatmap ---
    heatmap_kws_default = dict(
        cmap=cmap,
        square=True,
        xticklabels=False,
        yticklabels=False,
        cbar=True,
        cbar_ax=ax_cbar,
    )
    heatmap_kws = merge(heatmap_kws_default, heatmap_kws or {})

    # seaborn passes norm via "norm" to matplotlib
    hm = sns.heatmap(C, ax=ax_hm, norm=norm, **heatmap_kws)

    # colorbar label / ticks
    if clabel:
        hm.collections[0].colorbar.set_label(clabel)

    if scale == "linear" and vmin is not None and vmax is not None:
        decimal = 10
        nsegments = 5
        cd_ticks = np.trunc(np.linspace(vmin, vmax, nsegments) * decimal) / decimal
        hm.collections[0].colorbar.set_ticks(cd_ticks)

    # --- marginals ---
    margin_kws_default = dict(edgecolor="k", linewidth=1)
    margin_kws = merge(margin_kws_default, margin_kws or {})
    facecolor = color if color is not None else sns.color_palette()[0]

    # Use bin centers for bar positions (cleaner than edges for seaborn layout)
    centers = (binedges[:-1] + binedges[1:]) / 2
    width = (binedges[1] - binedges[0]) if len(binedges) > 1 else 1.0

    # top bars (x)
    ax_top.bar(
        centers, groupmean, width=width, align="center",
        color=facecolor, **margin_kws
    )
    ax_top.set_xlim(lo, hi)
    ax_top.set_xticks([])
    ax_top.set_yticks([])
    for s in ["top", "right", "left"]:
        ax_top.spines[s].set_visible(False)

    # left bars (y) — horizontal, and flip x-axis like your original
    ax_left.barh(
        centers, groupmean, height=width, align="center",
        color=facecolor, **margin_kws
    )
    ax_left.invert_xaxis()
    ax_left.set_ylim(hi, lo)  # match your inverted y
    ax_left.set_xticks([])
    ax_left.set_yticks([])
    for s in ["top", "bottom", "left"]:
        ax_left.spines[s].set_visible(False)

    # --- labels / title ---
    if title is not None:
        ax_top.set_title(title)
    if xlabel is not None:
        ax_hm.set_xlabel(xlabel)
    if ylabel is not None:
        ax_left.set_ylabel(ylabel)

    return grid

def plot_saddle_strength_profiles(
    samples,
    *,
    n_groups=None,
    labels=None,
    ax=None,
    fig_kws=None,
    line_kws=None,
    title="saddle strength profile",
    xlabel="extent",
    ylabel="(AA + BB) / (AB + BA)",
    baseline=1.0,
    colors=COLORS,
):
    """
    Plot multiple saddle strength profiles on one figure using seaborn.

    Parameters
    ----------
    samples : sequence
        A sequence of (int_s, int_c) tuples, or dict-like label -> (int_s, int_c).
        - int_s, int_c are the saddle 'signal' and 'count' matrices expected by cooltools.api.saddle.saddle_strength.
    n_groups : int, optional
        Number of groups (N_GROUPS). If None, inferred from int_s shape as (n-2).
    labels : sequence of str, optional
        Labels for each sample (only used if `samples` is a list/tuple).
    ax : matplotlib Axes, optional
        Plot onto an existing axes.
    fig_kws : dict, optional
        Passed to plt.figure() when ax is None.
    line_kws : dict, optional
        Passed to seaborn.lineplot() (e.g., linewidth, alpha, linestyle).
    baseline : float
        Horizontal reference line (typically 1.0 for ratios).

    Returns
    -------
    ax : matplotlib Axes
    df : pandas.DataFrame
        Tidy data used for plotting.
    """
    import pandas as pd

    # Normalize input to (labels, tuples)
    if hasattr(samples, "items"):  # dict-like
        items = list(samples.items())
        sample_labels = [k for k, _ in items]
        sample_pairs = [v for _, v in items]
    else:
        sample_pairs = list(samples)
        if labels is None:
            sample_labels = [f"sample_{i+1}" for i in range(len(sample_pairs))]
        else:
            if len(labels) != len(sample_pairs):
                raise ValueError("labels must be the same length as samples")
            sample_labels = list(labels)

    # Infer N_GROUPS if needed
    if n_groups is None:
        int_s0, _ = sample_pairs[0]
        n = np.asarray(int_s0).shape[0]
        if n < 3:
            raise ValueError("int_s looks too small to be (N_GROUPS+2, N_GROUPS+2)")
        n_groups = n - 2

    x = np.arange(n_groups + 2)

    # Build tidy dataframe
    rows = []
    for lab, (int_s, int_c) in zip(sample_labels, sample_pairs):
        y = saddle_strength(int_s, int_c)
        if len(y) != len(x):
            raise ValueError(
                f"{lab}: saddle_strength length {len(y)} != expected {len(x)} "
                f"(n_groups={n_groups})"
            )
        rows.extend({"extent": int(ext), "value": float(val), "sample": lab} for ext, val in zip(x, y))

    df = pd.DataFrame(rows)

    # Plot
    if ax is None:
        fig_kws = fig_kws or {}
        plt.figure(**fig_kws)
        ax = plt.gca()

    line_kws_default = dict(drawstyle="steps-pre")  # matches your plt.step(..., where='pre')
    line_kws = {**line_kws_default, **(line_kws or {})}

    sns.lineplot(
        data=df,
        x="extent",
        y="value",
        hue="sample",
        ax=ax,
        **line_kws,
        palette=colors,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if baseline is not None:
        ax.axhline(baseline, color="grey", linestyle="--", linewidth=1)

    ax.set_xlim(0, len(x) - 1)
    ax.set_ylim(bottom=0.75)

    if ax.legend_ is not None:
        ax.legend_.set_title("Sample")

    return ax, df

def hic_upper_triangle(
    mat,
    ax=None,
    cmap="RdBu_r",
    vmin=None,
    vmax=None,
    center=None,
    robust=False,
    robust_pct=(2, 98),
    symmetric=False,
    colorbar=True,
    cbar_kws=None,
    hide_axes=True,
    triangle_only=True,
):
    """
    Hi-C style 45° rotated (diamond) heatmap showing ONLY the upper triangle.

    Parameters
    ----------
    mat : (N, N) array
        Square matrix.
    center : float or None
        If not None, center the colormap at this value (seaborn-like center=0).
        Uses matplotlib.colors.TwoSlopeNorm.
    robust : bool
        If True and vmin/vmax not provided, use percentiles robust_pct.
    robust_pct : tuple(float, float)
        Percentiles used when robust=True (e.g., (2, 98)).
    symmetric : bool
        If True and center is not None and vmin/vmax not provided, use symmetric
        limits around center based on max abs deviation.
    triangle_only : bool
        If True, shows a single triangle (like common Hi-C panels) instead of full diamond.
    """
    mat = np.asarray(mat)
    n, m = mat.shape
    if n != m:
        raise ValueError("mat must be square (N x N).")

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))
    else:
        fig = ax.figure

    # Mask lower triangle (keep i <= j)
    mask = np.tril(np.ones_like(mat, dtype=bool), k=-1)
    Z = np.ma.array(mat, mask=mask)

    # Choose scale limits if not provided
    data = mat[~mask]  # only upper-triangle values
    data = data[np.isfinite(data)]

    if vmin is None or vmax is None:
        if robust and data.size:
            lo, hi = np.percentile(data, robust_pct)
        else:
            lo, hi = (data.min(), data.max()) if data.size else (0.0, 1.0)

        if vmin is None:
            vmin = lo
        if vmax is None:
            vmax = hi

    # Build normalization
    norm = None
    if center is not None:
        if symmetric and (vmin is None or vmax is None):
            # (kept for completeness; vmin/vmax already set above)
            v = np.nanmax(np.abs(data - center))
            vmin, vmax = center - v, center + v
        elif symmetric:
            v = max(abs(vmin - center), abs(vmax - center))
            vmin, vmax = center - v, center + v

        norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)

    # Pixel-edge grid then rotate coords by 45°: x=(i+j), y=(j-i)
    e = np.arange(n + 1)
    I, J = np.meshgrid(e, e, indexing="ij")
    X = I + J
    Y = J - I

    pc = ax.pcolormesh(
        X, Y, Z,
        shading="auto",
        cmap=cmap,
        norm=norm,
        vmin=None if norm is not None else vmin,
        vmax=None if norm is not None else vmax,
    )

    ax.set_aspect("equal")
    ax.set_xlim(0, 2 * n)

    # Show either a single triangle (common) or full diamond
    if triangle_only:
        ax.set_ylim(0, n)
    else:
        ax.set_ylim(-n, n)

    if hide_axes:
        ax.axis("off")

    if colorbar:
        cbar_kws = cbar_kws or {}
        fig.colorbar(pc, ax=ax, **cbar_kws)

    return ax

def plot_flanks_start_end_stackup(
    df,
    bigwigs,
    flank_bp=500_000,
    bins=200,
    agg_across_regions="mean",   # "mean" or "median"
    summary_per_bin="mean",      # 'mean', 'min', 'max', 'cov', 'std', or 'sum'
    shade="sem",                 # "sem", "std", or None (across regions after replicate-mean)
    hue=None,                    # e.g. "Q"
    hue_order=None,
    palette="tab10",             # None | colormap name | list of colors | dict {group: color}
    smooth_window=None,          # int (bins), moving average
    smooth_sigma=None,           # float (bins), gaussian sigma
    delta_keys=None,             # if bigwigs is dict: ("treated","control") => treated - control
    boundary_mode="concat",      # "concat" OR "aggregate"
    znorm=None,                  # None | "region" | "bin" | "global"
    gap_bp=250_000,              # extra gap between 5' and 3' segments in concat mode (bp)
    vline_color="k",             # vertical marker color
    vline_alpha=0.8,             # vertical marker alpha
    color_alpha=None,            # force alpha onto any palette/provided color (lines + fill base rgb)
    line_alpha=1.0,              # line alpha
    lw=1.0,                      # line width
    shade_alpha=None,            # fill alpha (defaults: 0.20 no-hue, 0.15 with hue)
    ax=None,                     # plot into existing axes
    figsize=(10, 4),
    do_tight_layout=True,
    site_type="boundary",
):
    """
    In concat mode:
      - 5' and 3' are plotted as two separate line segments (no connecting line).
      - gap_bp controls the horizontal spacing between segments.
      - duplicate 5' coordinates are deduplicated by (chrom, start) before aggregation.
      - duplicate 3' coordinates are deduplicated by (chrom, end) before aggregation.
    """

    # ---------------- helpers ----------------
    def _force_alpha(color, a):
        if color is None or a is None:
            return color
        r, g, b, _ = to_rgba(color)
        return (r, g, b, float(a))

    def _rgb(color):
        r, g, b, _ = to_rgba(color)
        return (r, g, b)

    def _smooth_1d(y, window=None, sigma=None):
        y = np.asarray(y, dtype=float)
        if (window is None or window <= 1) and (sigma is None or sigma <= 0):
            return y

        nanmask = np.isnan(y)
        v = np.where(nanmask, 0.0, y)
        w = np.where(nanmask, 0.0, 1.0)

        if sigma is not None and sigma > 0:
            rad = int(np.ceil(3 * sigma))
            xk = np.arange(-rad, rad + 1)
            k = np.exp(-(xk**2) / (2 * sigma**2))
            k = k / k.sum()
        else:
            window = int(max(2, window))
            k = np.ones(window, dtype=float) / window

        v_s = np.convolve(v, k, mode="same")
        w_s = np.convolve(w, k, mode="same")
        return v_s / np.where(w_s == 0, np.nan, w_s)

    def _get_color_map(keys, palette):
        if palette is None:
            return {k: None for k in keys}
        if isinstance(palette, dict):
            return {k: palette.get(k, None) for k in keys}
        if isinstance(palette, (list, tuple)):
            cols = list(palette)
            return {k: cols[i % len(cols)] if cols else None for i, k in enumerate(keys)}
        cmap = plt.get_cmap(str(palette))
        n = max(len(keys), 1)
        return {k: cmap(i / max(n - 1, 1)) for i, k in enumerate(keys)}

    def _pick_single_color_from_palette(palette):
        """If only one color is needed (no hue), use index 0."""
        if palette is None:
            return None
        if isinstance(palette, dict):
            return next(iter(palette.values()), None)
        if isinstance(palette, (list, tuple)):
            return palette[0] if len(palette) else None
        cmap = plt.get_cmap(str(palette))
        return cmap(0.0)

    def _aggregate(mat):
        n = mat.shape[0]
        prof = np.nanmedian(mat, axis=0) if agg_across_regions == "median" else np.nanmean(mat, axis=0)

        if shade is None:
            return prof, None, None, n

        spread = np.nanstd(mat, axis=0)
        if shade == "sem":
            band = spread / np.sqrt(max(n, 1))
        elif shade == "std":
            band = spread
        else:
            raise ValueError("shade must be 'sem', 'std', or None")

        return prof, prof - band, prof + band, n

    def _znorm_mat(mat, mode):
        mat = np.asarray(mat, dtype=float)
        if mode is None:
            return mat
        mode = str(mode).lower()

        if mode in ("region", "row"):
            mean = np.nanmean(mat, axis=1, keepdims=True)
            std  = np.nanstd(mat, axis=1, keepdims=True)
            z = (mat - mean) / np.where(std < 1e-12, np.nan, std)
            z = np.where(np.broadcast_to(std < 1e-12, z.shape), 0.0, z)
            return z

        if mode in ("bin", "col", "column"):
            mean = np.nanmean(mat, axis=0, keepdims=True)
            std  = np.nanstd(mat, axis=0, keepdims=True)
            z = (mat - mean) / np.where(std < 1e-12, np.nan, std)
            z = np.where(np.broadcast_to(std < 1e-12, z.shape), 0.0, z)
            return z

        if mode == "global":
            mean = np.nanmean(mat)
            std  = np.nanstd(mat)
            if std < 1e-12:
                return np.zeros_like(mat)
            return (mat - mean) / std

        raise ValueError("znorm must be None, 'region', 'bin', or 'global'.")

    def _apply_smoothing_concat(arr):
        if arr is None:
            return None
        arr = np.asarray(arr, dtype=float).copy()
        if smooth_sigma is not None and smooth_sigma > 0:
            arr[:bins] = _smooth_1d(arr[:bins], sigma=smooth_sigma)
            arr[bins:] = _smooth_1d(arr[bins:], sigma=smooth_sigma)
        elif smooth_window is not None and smooth_window > 1:
            arr[:bins] = _smooth_1d(arr[:bins], window=smooth_window)
            arr[bins:] = _smooth_1d(arr[bins:], window=smooth_window)
        return arr

    def _apply_smoothing_single(arr):
        if arr is None:
            return None
        arr = np.asarray(arr, dtype=float).copy()
        if smooth_sigma is not None and smooth_sigma > 0:
            return _smooth_1d(arr, sigma=smooth_sigma)
        if smooth_window is not None and smooth_window > 1:
            return _smooth_1d(arr, window=smooth_window)
        return arr

    def _dedup_rows_by_coord(mat, chrom_arr, pos_arr):
        """
        Deduplicate rows by exact genomic coordinate while preserving first-seen order.
        Used only in concat mode:
          - 5' side: (chrom, start)
          - 3' side: (chrom, end)
        """
        keep_idx = []
        seen = set()
        for i, (c, p) in enumerate(zip(chrom_arr, pos_arr)):
            key = (str(c), int(p))
            if key not in seen:
                seen.add(key)
                keep_idx.append(i)

        keep_idx = np.asarray(keep_idx, dtype=int)
        return mat[keep_idx], keep_idx

    # ---------------- columns / windows ----------------
    df = df.rename(columns={c: c.strip() for c in df.columns})
    for col in ("chrom", "start", "end"):
        if col not in df.columns:
            raise ValueError(f"df missing required column: '{col}'")

    chroms = df["chrom"].astype(str).to_numpy()
    starts = df["start"].astype(int).to_numpy()
    ends   = df["end"].astype(int).to_numpy()

    s0 = np.maximum(starts - flank_bp, 0)
    e0 = starts + flank_bp
    s1 = np.maximum(ends - flank_bp, 0)
    e1 = ends + flank_bp

    try:
        import bbi
    except ImportError as e:
        raise ImportError("Need `bbi` installed for bbi.stackup.") from e

    def _stackup_repmean(bw_list, summary='mean'):
        sm, em = [], []
        for bw in bw_list:
            sm.append(bbi.stackup(bw, chroms, s0, e0, bins=bins, oob=np.nan, summary=summary))
            em.append(bbi.stackup(bw, chroms, s1, e1, bins=bins, oob=np.nan, summary=summary))
        sm = np.nanmean(np.stack(sm, axis=0), axis=0)
        em = np.nanmean(np.stack(em, axis=0), axis=0)
        return sm, em

    # ---------------- start/end matrices ----------------
    plot_label = "replicate-mean signal"
    if isinstance(bigwigs, dict):
        if len(bigwigs) != 2:
            raise ValueError("bigwigs dict must have exactly 2 keys to compute delta.")
        keys = list(bigwigs.keys())

        if delta_keys is None:
            pos_key, neg_key = keys[1], keys[0]
        else:
            if not (isinstance(delta_keys, (tuple, list)) and len(delta_keys) == 2):
                raise ValueError("delta_keys must be a 2-tuple like ('treated','control').")
            pos_key, neg_key = delta_keys
            if pos_key not in bigwigs or neg_key not in bigwigs:
                raise ValueError(f"delta_keys must match dict keys. Got {delta_keys}, keys={keys}.")

        s_pos, e_pos = _stackup_repmean(bigwigs[pos_key], summary=summary_per_bin)
        s_neg, e_neg = _stackup_repmean(bigwigs[neg_key], summary=summary_per_bin)
        start_mat = s_pos - s_neg
        end_mat   = e_pos - e_neg
        plot_label = f"delta: {pos_key} - {neg_key}"
    else:
        start_mat, end_mat = _stackup_repmean(list(bigwigs), summary=summary_per_bin)

    if znorm is not None:
        start_mat = _znorm_mat(start_mat, znorm)
        end_mat   = _znorm_mat(end_mat, znorm)

    if boundary_mode not in ("concat", "aggregate"):
        raise ValueError("boundary_mode must be 'concat' or 'aggregate'.")

    # x axis
    x_rel = np.linspace(-flank_bp, flank_bp, bins, endpoint=False) + (flank_bp * 2 / bins) / 2
    if boundary_mode == "aggregate":
        boundary_mat = np.nanmean(np.stack([start_mat, end_mat], axis=0), axis=0)
        x = x_rel / 1e6
    else:
        boundary_mat = None
        x = np.concatenate([x_rel, x_rel + 2 * flank_bp + gap_bp]) / 1e6  # length 2*bins

    # axes setup
    created_ax = ax is None
    if created_ax:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    if shade_alpha is None:
        shade_alpha = 0.15 if hue is not None else 0.20

    def _plot_profile(label, prof, lo, hi, requested_color=None):
        # smoothing
        if boundary_mode == "concat":
            prof = _apply_smoothing_concat(prof)
            lo   = _apply_smoothing_concat(lo)
            hi   = _apply_smoothing_concat(hi)
        else:
            prof = _apply_smoothing_single(prof)
            lo   = _apply_smoothing_single(lo)
            hi   = _apply_smoothing_single(hi)

        base_color = _force_alpha(requested_color, color_alpha) if requested_color is not None else None

        if boundary_mode == "concat":
            x1, x2 = x[:bins], x[bins:]
            y1, y2 = prof[:bins], prof[bins:]
            lo1 = None if lo is None else lo[:bins]
            hi1 = None if hi is None else hi[:bins]
            lo2 = None if lo is None else lo[bins:]
            hi2 = None if hi is None else hi[bins:]

            (line1,) = ax.plot(x1, y1, color=base_color, label=label, lw=lw, alpha=line_alpha)
            actual = line1.get_color()
            actual = _force_alpha(actual, color_alpha) if color_alpha is not None else actual

            ax.plot(x2, y2, color=actual, label="_nolegend_", lw=lw, alpha=line_alpha)

            if lo is not None and hi is not None:
                fill = (*_rgb(actual), float(shade_alpha))
                ax.fill_between(x1, lo1, hi1, color=fill)
                ax.fill_between(x2, lo2, hi2, color=fill)

        else:
            (line,) = ax.plot(x, prof, color=base_color, label=label, lw=lw, alpha=line_alpha)
            actual = line.get_color()
            actual = _force_alpha(actual, color_alpha) if color_alpha is not None else actual

            if lo is not None and hi is not None:
                fill = (*_rgb(actual), float(shade_alpha))
                ax.fill_between(x, lo, hi, color=fill)

    # ---------------- plot profiles ----------------
    concat_n_start = None
    concat_n_end = None

    if hue is None:
        single_color = _pick_single_color_from_palette(palette)

        if boundary_mode == "aggregate":
            p, lo, hi, n = _aggregate(boundary_mat)
            _plot_profile(None, p, lo, hi, requested_color=single_color)
            total_n = n

        else:
            start_mat_dedup, keep_start = _dedup_rows_by_coord(start_mat, chroms, starts)
            end_mat_dedup, keep_end     = _dedup_rows_by_coord(end_mat, chroms, ends)

            ps, los, his, n_start = _aggregate(start_mat_dedup)
            pe, loe, hie, n_end   = _aggregate(end_mat_dedup)

            prof = np.concatenate([ps, pe])
            lo   = None if los is None else np.concatenate([los, loe])
            hi   = None if his is None else np.concatenate([his, hie])

            _plot_profile(None, prof, lo, hi, requested_color=single_color)

            concat_n_start = n_start
            concat_n_end   = n_end
            total_n = None

    else:
        if hue not in df.columns:
            raise ValueError(f"hue='{hue}' not found in df.")
        g = df[hue].to_numpy().astype(str)

        if hue_order is None:
            try:
                uniq = np.unique(df[hue].astype(float))
                group_keys = [str(int(u)) if float(u).is_integer() else str(u) for u in np.sort(uniq)]
            except Exception:
                group_keys = sorted(np.unique(g))
        else:
            group_keys = [str(k) for k in hue_order]

        colors = _get_color_map(group_keys, palette)
        total_n = start_mat.shape[0]

        for k in group_keys:
            mask = (g == k)
            if not np.any(mask):
                continue

            c = colors.get(k, None)
            c = _force_alpha(c, color_alpha) if c is not None else None

            if boundary_mode == "aggregate":
                mat = np.nanmean(np.stack([start_mat[mask], end_mat[mask]], axis=0), axis=0)
                p, lo, hi, n_g = _aggregate(mat)
                _plot_profile(f"{hue}={k} (n={n_g})", p, lo, hi, requested_color=c)

            else:
                start_g_dedup, _ = _dedup_rows_by_coord(start_mat[mask], chroms[mask], starts[mask])
                end_g_dedup, _   = _dedup_rows_by_coord(end_mat[mask], chroms[mask], ends[mask])

                ps, los, his, n5 = _aggregate(start_g_dedup)
                pe, loe, hie, n3 = _aggregate(end_g_dedup)

                prof = np.concatenate([ps, pe])
                lo   = None if los is None else np.concatenate([los, loe])
                hi   = None if his is None else np.concatenate([his, hie])

                _plot_profile(f"{hue}={k} (5′ n={n5}, 3′ n={n3})", prof, lo, hi, requested_color=c)

        ax.legend(frameon=False, fontsize=9)

    # ---------------- vertical markers ----------------
    vcol = _force_alpha(vline_color, None)

    if boundary_mode == "aggregate":
        ax.axvline(0, ls="--", lw=1, color=vcol, alpha=vline_alpha)
        tick_pos = np.array([-flank_bp, 0, flank_bp]) / 1e6
        tick_lab = [f"-{flank_bp/1e3:.0f} kb", site_type, f"+{flank_bp/1e3:.0f} kb"]
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_lab)
        ax.set_title(f"Aggregated 5′ & 3′ {site_type} ±{flank_bp:,} bp (total n={total_n}, bins={bins})")

    else:
        # positions in Mb
        end_5p = flank_bp / 1e6
        start_3p = (flank_bp + gap_bp) / 1e6
        anchor_3p = (2 * flank_bp + gap_bp) / 1e6

        ax.axvline(0, ls="--", lw=1, color=vcol, alpha=vline_alpha)           # 5' anchor
        ax.axvline(end_5p, ls=":", lw=1, color=vcol, alpha=vline_alpha)       # end of 5' segment
        ax.axvline(start_3p, ls=":", lw=1, color=vcol, alpha=vline_alpha)     # start of 3' segment
        ax.axvline(anchor_3p, ls="--", lw=1, color=vcol, alpha=vline_alpha)   # 3' anchor

        tick_pos = np.array([
            -flank_bp, 0, flank_bp,
            flank_bp + gap_bp, 2 * flank_bp + gap_bp, 3 * flank_bp + gap_bp
        ]) / 1e6
        tick_lab = [
            f"-{flank_bp/1e3:.0f} kb", f"5′ {site_type}", f"+{flank_bp/1e3:.0f} kb",
            f"-{flank_bp/1e3:.0f} kb", f"3′ {site_type}", f"+{flank_bp/1e3:.0f} kb",
        ]
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_lab)

        if hue is None and concat_n_start is not None and concat_n_end is not None:
            ax.set_title(
                f"5′ ±{flank_bp:,} bp | gap {gap_bp:,} bp | 3′ ±{flank_bp:,} bp "
                f"(unique 5′ n={concat_n_start}, unique 3′ n={concat_n_end}, bins={bins})"
            )
        else:
            ax.set_title(
                f"5′ ±{flank_bp:,} bp | gap {gap_bp:,} bp | 3′ ±{flank_bp:,} bp "
                f"(concat mode dedups 5′/3′ coordinates separately, bins={bins})"
            )

    ax.set_xlabel(f"Position around {site_type}")
    ax.set_ylabel(plot_label + (f" (z={znorm})" if znorm else ""))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if created_ax and do_tight_layout:
        fig.tight_layout()

    return fig, ax, {
        "start_mat": start_mat,
        "end_mat": end_mat,
        "boundary_mat": boundary_mat,
        "label": plot_label,
        "boundary_mode": boundary_mode,
        "znorm": znorm,
        "gap_bp": int(gap_bp) if boundary_mode == "concat" else 0,
        "concat_dedup": bool(boundary_mode == "concat"),
        "concat_unique_start_n": int(concat_n_start) if concat_n_start is not None else None,
        "concat_unique_end_n": int(concat_n_end) if concat_n_end is not None else None,
    }

def lineplot_stackup_dict(
    stackup_dict,
    ax=None,
    *,
    x=None,
    agg="mean",          # "mean" or "median"
    err="sem",           # "sem", "sd", "ci95", or None
    err_axis=0,          # axis to compute error over (0 = rows/regions if v is 2D)
    alpha=0.2,           # fill_between alpha
    labels=None,         # optional mapping: key -> label
    plot_kwargs=None,    # forwarded to ax.plot
    fill_kwargs=None,    # forwarded to ax.fill_between
    structure_type="TAD",
    colors=None,         # mapping key->color OR list/tuple of colors OR single color
    legend_title=None,
    y_title=None,
    despine_fn=None,     # e.g., seaborn.despine or your plotting.despine
):
    """
    Plot 1D profiles from a dict of stackups.
    Values can be:
      - 1D: shape (nbins,)
      - 2D: shape (n, nbins) -> aggregates over err_axis to produce center + error band

    Shaded band defaults to the same color as the line unless overridden via fill_kwargs.
    Returns: ax
    """
    if ax is None:
        _, ax = plt.subplots()

    plot_kwargs = {} if plot_kwargs is None else dict(plot_kwargs)
    fill_kwargs = {} if fill_kwargs is None else dict(fill_kwargs)

    items = list(stackup_dict.items())

    for i, (k, v) in enumerate(items):
        arr = np.asarray(v)

        # pick color for this series
        if colors is None:
            current_color = plot_kwargs.get("color", None)  # let mpl cycle if None
        elif isinstance(colors, dict):
            current_color = colors.get(k, plot_kwargs.get("color", None))
        elif isinstance(colors, (list, tuple)):
            current_color = colors[i % len(colors)]
        else:
            current_color = colors  # single color string

        # x
        nbins = arr.shape[-1]
        xx = np.arange(nbins) if x is None else np.asarray(x)
        if xx.shape[0] != nbins:
            raise ValueError(f"x has length {xx.shape[0]} but expected {nbins} for key={k}")

        # center + band
        if arr.ndim == 1:
            center = arr
            band = None
        elif arr.ndim == 2:
            if agg == "mean":
                center = np.nanmean(arr, axis=err_axis)
            elif agg == "median":
                center = np.nanmedian(arr, axis=err_axis)
            else:
                raise ValueError("agg must be 'mean' or 'median'")

            band = None
            if err is not None:
                arr2 = np.moveaxis(arr, err_axis, 0) if err_axis != 0 else arr

                if err == "sd":
                    band = np.nanstd(arr2, axis=0, ddof=1)
                elif err in ("sem", "ci95"):
                    n = np.sum(~np.isnan(arr2), axis=0)
                    sd = np.nanstd(arr2, axis=0, ddof=1)
                    sem = sd / np.sqrt(np.maximum(n, 1))
                    band = sem if err == "sem" else 1.96 * sem
                else:
                    raise ValueError("err must be one of: None, 'sem', 'sd', 'ci95'")
        else:
            raise ValueError(f"Value for key={k} must be 1D or 2D, got shape {arr.shape}")

        # label
        if isinstance(labels, dict):
            label = labels.get(k, k)
        elif labels is None:
            label = k
        else:
            label = labels

        # plot line (capture actual line color if we let mpl cycle)
        (line,) = ax.plot(xx, center, label=label, color=current_color, **plot_kwargs)
        line_color = line.get_color()

        # band defaults to line color unless fill_kwargs overrides
        if band is not None:
            fb_kwargs = dict(fill_kwargs)
            fb_kwargs.setdefault("color", line_color)   # default = line color
            fb_kwargs.setdefault("linewidth", 0)
            ax.fill_between(xx, center - band, center + band, alpha=alpha, **fb_kwargs)

    # optional despine
    if despine_fn is not None:
        despine_fn(ax=ax)

    # guide lines + xticks
    ax.axvline(25,  color="grey", linestyle="--", linewidth=1.0, alpha=0.7)
    ax.axvline(125, color="grey", linestyle="--", linewidth=1.0, alpha=0.7)

    ax.set_xticks([25, 75, 125])
    if structure_type == "TAD":
        ax.set_xticklabels(["5' Boundary", "TAD", "3' Boundary"])
    elif structure_type == "Loop":
        ax.set_xticklabels(["5' Anchor", "Loop", "3' Anchor"])

    if y_title is not None:
        ax.set_ylabel(y_title)

    if legend_title is not None:
        ax.legend(title=legend_title)
    else:
        ax.legend()

    return ax
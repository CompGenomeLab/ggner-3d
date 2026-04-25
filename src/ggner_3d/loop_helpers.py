import numpy as np
from ggner_3d import plotting
from ggner_3d.dots_common import compare_loops_with_strength
from matplotlib import pyplot as plt

def format_loops_columns(df, how='sample_a'):
    df = df.rename(columns={
        f'chrom1_{how}': 'chrom1',
        f'start1_{how}': 'start1',
        f'end1_{how}': 'end1',
        f'chrom2_{how}': 'chrom2',
        f'start2_{how}': 'start2',
        f'end2_{how}': 'end2',
    })
    return df

def plot_anchor_loop_counts(
    s1,
    s2,
    *,
    clr_,
    dots_dfs_,
    clustering_radius,
    compute_strength=False,
    how="sample_a",
    labels=None,
    colors=None,
    x_positions=(0.0, 0.4),
    figsize=(3.5 * 1.2, 2.5 * 1.2),
    bar_width=0.35,
    ax=None,
    format_loops_columns=format_loops_columns,
    compare_loops_with_strength=compare_loops_with_strength,
    plotting=plotting
):
    """
    Compare s1 vs s2 and plot stacked bars for counts of Anchors and Loops.

    Returns:
        res (dict): output from compare_loops_with_strength (with formatted common_loops)
        anchor_counts (list[int]): [common, s1-specific, s2-specific]
        loop_counts (list[int]):   [common, s1-specific, s2-specific]
        fig, ax: matplotlib objects
    """
    keys = list(clr_.keys())
    if s1 not in keys or s2 not in keys:
        missing = [s for s in (s1, s2) if s not in keys]
        raise KeyError(f"Sample(s) not found in clr_.keys(): {missing}")

    s1_i = keys.index(s1)
    s2_i = keys.index(s2)

    res = compare_loops_with_strength(
        df_a=dots_dfs_[s1_i],
        df_b=dots_dfs_[s2_i],
        compute_strength=compute_strength,
        clustering_radius=clustering_radius,
    )

    # format common loops (same as your snippet)
    res["common_loops"] = format_loops_columns(res["common_loops"], how=how)

    # counts (note: anchors use *_anchors keys; loops use loops keys)
    loop_counts = [
        len(res["common_loops"]),
        len(res["sample_a_specific"]),
        len(res["sample_b_specific"]),
    ]
    anchor_counts = [
        len(res["common_anchors"]),
        len(res["sample_a_specific_anchors"]),
        len(res["sample_b_specific_anchors"]),
    ]

    print("Anchors counts:", anchor_counts)
    print("Loops counts:", loop_counts)

    if labels is None:
        labels = ["Common", f"{s1} - Specific", f"{s2} - Specific"]
    if colors is None:
        colors = getattr(plotting, "COLORS", None)

    x = np.array(x_positions, dtype=float)  # 0=Anchors, 1=Loops (just spaced however you like)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    bottom = np.zeros_like(x, dtype=float)
    for i, lab in enumerate(labels):
        heights = [anchor_counts[i], loop_counts[i]]
        c = colors[i] if colors is not None else None
        ax.bar(x, heights, bottom=bottom, label=lab, color=c, width=bar_width)
        bottom += heights

    ax.set_xticks(x)
    ax.set_xticklabels(["Anchors", "Loops"])
    ax.set_ylabel("Count")
    ax.legend(frameon=False)
    plotting.despine(ax)
    ax.margins(x=0.2)
    plt.tight_layout()

    return res, anchor_counts, loop_counts, fig, ax
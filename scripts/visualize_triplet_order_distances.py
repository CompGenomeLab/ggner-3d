#!/usr/bin/env python3
"""Visualize XPC/CTCF/PolII/promoter order summaries.

The input is the descriptive order table produced at:
misc/xpc_ctcf_pol2_triplets_promoter_orders_with_distances.tsv

The table may contain multiple coordinate systems, including active-TSS
transcription-normalized order and CTCF motif-normalized order.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = {
    "Promoter": "#222222",
    "PolII": "#DBB82A",
    "POLR2A": "#DBB82A",
    "XPC": "#D95F02",
    "CTCF": "#1B9E77",
    "RAD21": "#6A3D9A",
}


def base_landmark(label: str) -> str:
    """Map labels such as CTCF(+) to their base landmark."""
    return label.split("(", 1)[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot triplet order frequencies and median spacing for one coordinate system."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("misc/xpc_ctcf_pol2_triplets_promoter_orders_with_distances.tsv"),
        help="Readable order-distance TSV.",
    )
    parser.add_argument(
        "--cutoff",
        type=int,
        default=5000,
        help="pairwise_cutoff_bp value to plot.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Number of top-ranked order patterns to show.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=None,
        help="Output prefix without extension. Defaults to figs/triplet_orders_cutoff<CUTOFF>_top<TOP>.",
    )
    parser.add_argument(
        "--max-position-bp",
        type=float,
        default=None,
        help=(
            "Clip the median-position schematic x-axis to +/- this many bp. "
            "Points outside the range are drawn at the edge with triangle markers."
        ),
    )
    parser.add_argument(
        "--order-coordinate-system",
        default=None,
        help=(
            "If the input contains multiple order systems, choose one. Examples: "
            "tx_active_tss, tx_active_tss_ctcf_motif_label, ctcf_motif_normalized."
        ),
    )
    parser.add_argument(
        "--title-prefix",
        default="XPC/CTCF/PolII triplet orders",
        help="Title prefix for the figure.",
    )
    return parser.parse_args()


def ordered_group_positions(row: pd.Series, distance_stat: str = "median") -> list[tuple[str, float]]:
    """Return order groups and signed positions relative to the reference landmark.

    The reference landmark is fixed at 0. Adjacent edge distances are
    reconstructed from the edge median or mean columns.
    """
    groups = str(row["order"]).split(" -> ")
    gaps = []
    for i in range(1, len(groups)):
        val = row.get(f"edge{i}_{distance_stat}_bp", np.nan)
        gaps.append(float(val) if pd.notna(val) else 0.0)

    reference_landmark = str(row.get("reference_landmark", "Promoter"))
    reference_idx = None
    for idx, group in enumerate(groups):
        members = [base_landmark(member) for member in group.split("=")]
        if reference_landmark in members:
            reference_idx = idx
            break
    if reference_idx is None:
        raise ValueError(f"Order has no {reference_landmark} group: {row['order']}")

    positions = [np.nan] * len(groups)
    positions[reference_idx] = 0.0

    for idx in range(reference_idx + 1, len(groups)):
        positions[idx] = positions[idx - 1] + gaps[idx - 1]
    for idx in range(reference_idx - 1, -1, -1):
        positions[idx] = positions[idx + 1] - gaps[idx]

    return list(zip(groups, positions))


def landmarks_in_orders(orders: pd.Series) -> list[str]:
    landmarks: list[str] = []
    for order in orders.astype(str):
        for group in order.split(" -> "):
            for member in group.split("="):
                landmark = base_landmark(member)
                if landmark not in landmarks:
                    landmarks.append(landmark)
    return landmarks


def plot_orders(
    df: pd.DataFrame,
    cutoff: int,
    top: int,
    output_prefix: Path,
    max_position_bp: float | None = None,
    order_coordinate_system: str | None = None,
    title_prefix: str = "XPC/CTCF/PolII triplet orders",
) -> None:
    subset = df[df["pairwise_cutoff_bp"].eq(cutoff)].copy()
    if subset.empty:
        available = ", ".join(map(str, sorted(df["pairwise_cutoff_bp"].unique())))
        raise SystemExit(f"No rows for cutoff={cutoff}. Available cutoffs: {available}")

    if "order_coordinate_system" in subset.columns:
        available_systems = list(dict.fromkeys(subset["order_coordinate_system"].astype(str).tolist()))
        if order_coordinate_system is None:
            order_coordinate_system = available_systems[0]
        if order_coordinate_system not in available_systems:
            raise SystemExit(
                f"No rows for order_coordinate_system={order_coordinate_system!r} at cutoff={cutoff}. "
                f"Available systems: {', '.join(available_systems)}"
            )
        subset = subset[subset["order_coordinate_system"].eq(order_coordinate_system)].copy()

    subset = subset.sort_values(["rank"]).head(top).copy()
    subset = subset.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(subset))
    ylabels = [
        f"{int(row.rank)}. {row.order}  ({int(row.n)}/{int(row.denominator_triplets)}, {row.pct:.2f}%)"
        for row in subset.itertuples(index=False)
    ]

    fig_height = max(7, 0.45 * len(subset) + 1.8)
    fig, (ax_bar, ax_pos) = plt.subplots(
        1,
        2,
        figsize=(19, fig_height),
        gridspec_kw={"width_ratios": [1.15, 1.85], "wspace": 0.08},
        sharey=True,
    )

    ax_bar.barh(y, subset["pct"], color="#5B6C83", height=0.72)
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(ylabels, fontsize=8.5)
    ax_bar.set_xlabel("Percent of triplets at this cutoff")
    ax_bar.set_title("Order frequency")
    ax_bar.grid(axis="x", color="#D7DCE2", linewidth=0.8)
    ax_bar.set_axisbelow(True)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.spines["left"].set_visible(False)

    all_positions: list[float] = []
    clipped_any = False
    for row_idx, (_, row) in enumerate(subset.iterrows()):
        groups = ordered_group_positions(row, "median")
        group_positions = [pos for _, pos in groups]
        all_positions.extend(group_positions)

        if max_position_bp is not None:
            plot_positions = [float(np.clip(pos, -max_position_bp, max_position_bp)) for pos in group_positions]
            clipped_flags = [pos < -max_position_bp or pos > max_position_bp for pos in group_positions]
            clipped_any = clipped_any or any(clipped_flags)
        else:
            plot_positions = group_positions
            clipped_flags = [False] * len(group_positions)

        ax_pos.plot(plot_positions, [row_idx] * len(plot_positions), color="#B8C0C9", linewidth=1.6)

        for (group, pos), plot_pos, clipped in zip(groups, plot_positions, clipped_flags):
            members = group.split("=")
            offsets = np.linspace(-0.10, 0.10, len(members)) if len(members) > 1 else [0.0]
            for member, offset in zip(members, offsets):
                marker = ">" if pos > plot_pos else "<" if pos < plot_pos else "o"
                ax_pos.scatter(
                    plot_pos,
                    row_idx + offset,
                    s=58,
                    color=COLORS.get(base_landmark(member), "#666666"),
                    edgecolor="white",
                    linewidth=0.7,
                    marker=marker if clipped else "o",
                    zorder=3,
                )

    if max_position_bp is not None:
        max_abs = float(max_position_bp)
        xpad = max(50.0, max_abs * 0.03)
    else:
        max_abs = max(abs(float(x)) for x in all_positions if np.isfinite(x))
        xpad = max(250.0, max_abs * 0.08)
    ax_pos.set_xlim(-max_abs - xpad, max_abs + xpad)
    ax_pos.axvline(0, color="#222222", linewidth=1.0)
    reference_landmark = str(subset["reference_landmark"].iloc[0]) if "reference_landmark" in subset.columns else "Promoter"
    if reference_landmark == "Promoter":
        xlabel = "Median position relative to promoter/TSS (bp; transcription-normalized)"
    elif reference_landmark == "CTCF":
        xlabel = "Median position relative to CTCF summit (bp; CTCF motif-forward direction)"
    else:
        xlabel = f"Median position relative to {reference_landmark} (bp)"
    if max_position_bp is not None:
        xlabel += f"\nclipped to +/- {max_position_bp:,.0f} bp; triangles indicate off-scale medians"
    ax_pos.set_xlabel(xlabel)
    ax_pos.set_title("Median spacing schematic")
    ax_pos.grid(axis="x", color="#D7DCE2", linewidth=0.8)
    ax_pos.set_axisbelow(True)
    ax_pos.spines["top"].set_visible(False)
    ax_pos.spines["right"].set_visible(False)
    ax_pos.spines["left"].set_visible(False)

    present_landmarks = landmarks_in_orders(subset["order"])
    legend_order = [label for label in COLORS if label in present_landmarks]
    legend_order.extend(label for label in present_landmarks if label not in COLORS)
    legend_handles = [
        ax_pos.scatter(
            [],
            [],
            s=70,
            color=COLORS.get(label, "#666666"),
            edgecolor="white",
            linewidth=0.7,
            label=label,
        )
        for label in legend_order
    ]
    ax_pos.legend(
        handles=legend_handles,
        loc="lower right",
        frameon=True,
        framealpha=0.95,
        title="Landmark",
    )

    denominator = int(subset["denominator_triplets"].iloc[0])
    system_label = ""
    if "order_coordinate_system" in subset.columns:
        system_label = f" | {subset['order_coordinate_system'].iloc[0]}"
    fig.suptitle(
        f"{title_prefix} | all pairwise <= {cutoff:,} bp | n={denominator:,}{system_label}",
        fontsize=14,
        y=0.985,
    )
    if clipped_any:
        ax_pos.text(
            0.995,
            1.015,
            "Off-scale points shown as triangles",
            transform=ax_pos.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.5,
            color="#4B5563",
        )
    fig.subplots_adjust(left=0.34, right=0.985, top=0.90, bottom=0.08, wspace=0.08)

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(output_prefix.with_suffix(".svg"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output_prefix = args.output_prefix
    if output_prefix is None:
        output_prefix = Path("figs") / f"triplet_orders_cutoff{args.cutoff}_top{args.top}"

    df = pd.read_csv(args.input, sep="\t")
    plot_orders(
        df,
        args.cutoff,
        args.top,
        output_prefix,
        args.max_position_bp,
        args.order_coordinate_system,
        args.title_prefix,
    )
    print(f"Wrote {output_prefix.with_suffix('.png')}")
    print(f"Wrote {output_prefix.with_suffix('.svg')}")


if __name__ == "__main__":
    main()

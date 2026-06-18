#!/usr/bin/env python3
"""Analyze XPC/CTCF/PolII triplet order relative to ATAC-active TSSs."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import bioframe as bf


PEAK_COLS = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "signalValue",
    "pValue",
    "qValue",
    "peak",
]
TSS_COLS = ["chrom", "start", "end", "gene", "score", "strand"]
MOTIF_COLS = [
    "chrom",
    "start",
    "end",
    "length",
    "strand",
    "motif",
    "motif_score",
    "pvalue",
    "qvalue",
    "sequence",
]
KEEP_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use ATAC-active TSSs for XPC/CTCF/PolII triplet promoter-order analysis."
    )
    parser.add_argument("--xpc", type=Path, default=Path("misc/XPC_all.narrowPeak"))
    parser.add_argument("--ctcf", type=Path, default=Path("misc/CTCF_all.narrowPeak"))
    parser.add_argument(
        "--pol2",
        type=Path,
        default=Path("misc/GSM7162996_NCRNAPolIIChIP_peaks.ucsc.narrowPeak"),
    )
    parser.add_argument("--tss", type=Path, default=Path("misc/tss.bed"))
    parser.add_argument("--atac", type=Path, default=Path("misc/WT_noUV.mRp.clN_peaks.narrowPeak"))
    parser.add_argument("--ctcf-motifs", type=Path, default=Path("misc/CTCF_hg38.bed"))
    parser.add_argument(
        "--ctcf-motif-max-dist",
        type=int,
        default=100,
        help="Maximum distance from CTCF summit to matched overlapping motif center for high-confidence orientation.",
    )
    parser.add_argument(
        "--require-ctcf-motif",
        action="store_true",
        help="Retain only triplets whose nearest CTCF peak has a high-confidence motif orientation.",
    )
    parser.add_argument(
        "--promoter-window",
        type=int,
        default=1000,
        help="TSS +/- this many bp must overlap an ATAC peak to be called active.",
    )
    parser.add_argument(
        "--pol2-tss-window",
        type=int,
        default=None,
        help=(
            "Optional promoter-proximal Pol II filter. If set, retain only triplets "
            "where the Pol II summit is within this many bp of the nearest active TSS "
            "in transcription-normalized coordinates."
        ),
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("misc/xpc_ctcf_pol2_active_tss"),
        help="Prefix for output files.",
    )
    return parser.parse_args()


def read_peaks(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=PEAK_COLS)
    df = df[df["chrom"].isin(KEEP_CHROMS)].copy().reset_index(drop=True)
    df[f"{label}_id"] = df.index
    peak = df["peak"].astype(int)
    df[f"{label}_summit"] = np.where(
        peak >= 0,
        df["start"].astype(int) + peak,
        ((df["start"].astype(int) + df["end"].astype(int)) // 2),
    )
    return df


def read_tss(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=TSS_COLS)
    df = df[df["chrom"].isin(KEEP_CHROMS)].copy().reset_index(drop=True)
    df["tss_id"] = df.index
    df["tss"] = df["start"].astype(int)
    return df


def read_ctcf_motifs(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=MOTIF_COLS)
    df = df[df["chrom"].isin(KEEP_CHROMS)].copy().reset_index(drop=True)
    df["motif_id"] = df.index
    df["motif_center"] = ((df["start"].astype(int) + df["end"].astype(int)) // 2)
    return df


def assign_ctcf_motifs(ctcf: pd.DataFrame, motifs: pd.DataFrame, max_dist: int) -> pd.DataFrame:
    """Assign each CTCF peak the closest overlapping motif to its summit."""
    overlaps = bf.overlap(
        ctcf[["chrom", "start", "end", "ctcf_id", "ctcf_summit"]],
        motifs[["chrom", "start", "end", "motif_id", "motif_center", "strand", "motif_score"]],
        return_index=False,
        suffixes=("_ctcf", "_motif"),
    )
    overlaps = overlaps.dropna(subset=["motif_id_motif"]).copy()
    overlaps["motif_id_motif"] = overlaps["motif_id_motif"].astype(int)
    overlaps["motif_center_motif"] = overlaps["motif_center_motif"].astype(int)
    overlaps["motif_score_motif"] = overlaps["motif_score_motif"].astype(float)
    overlaps["ctcf_motif_dist_bp"] = (
        overlaps["motif_center_motif"] - overlaps["ctcf_summit_ctcf"].astype(int)
    ).abs()

    best = (
        overlaps.sort_values(
            ["ctcf_id_ctcf", "ctcf_motif_dist_bp", "motif_score_motif"],
            ascending=[True, True, False],
        )
        .drop_duplicates("ctcf_id_ctcf")
        .copy()
    )

    out = ctcf.copy()
    out["ctcf_motif_strand"] = pd.NA
    out["ctcf_motif_center"] = np.nan
    out["ctcf_motif_score"] = np.nan
    out["ctcf_motif_dist_bp"] = np.nan

    ids = best["ctcf_id_ctcf"].to_numpy(dtype=int)
    out.loc[ids, "ctcf_motif_strand"] = best["strand_motif"].to_numpy()
    out.loc[ids, "ctcf_motif_center"] = best["motif_center_motif"].to_numpy()
    out.loc[ids, "ctcf_motif_score"] = best["motif_score_motif"].to_numpy()
    out.loc[ids, "ctcf_motif_dist_bp"] = best["ctcf_motif_dist_bp"].to_numpy()
    out["ctcf_has_overlapping_motif"] = out["ctcf_motif_strand"].notna()
    out["ctcf_high_conf_motif"] = out["ctcf_motif_dist_bp"].le(max_dist)
    return out


def nearest(query: pd.DataFrame, target: pd.DataFrame, qpos: str, tpos: str) -> tuple[np.ndarray, np.ndarray]:
    query = query.reset_index(drop=True)
    nearest_idx = np.full(len(query), -1, dtype=np.int64)
    nearest_dist = np.full(len(query), np.nan)
    for chrom, qsub in query.groupby("chrom", sort=False):
        tsub = target[target["chrom"].eq(chrom)].sort_values(tpos)
        if tsub.empty:
            continue
        target_pos = tsub[tpos].to_numpy(dtype=np.int64)
        target_idx = tsub.index.to_numpy(dtype=np.int64)
        query_pos = qsub[qpos].to_numpy(dtype=np.int64)
        insert = np.searchsorted(target_pos, query_pos, side="left")
        left = np.clip(insert - 1, 0, len(target_pos) - 1)
        right = np.clip(insert, 0, len(target_pos) - 1)
        left_dist = np.abs(query_pos - target_pos[left])
        right_dist = np.abs(query_pos - target_pos[right])
        best = np.where(right_dist < left_dist, right, left)
        rows = qsub.index.to_numpy(dtype=np.int64)
        nearest_idx[rows] = target_idx[best]
        nearest_dist[rows] = np.abs(query_pos - target_pos[best])
    return nearest_idx, nearest_dist


def interval_overlap_any(query: pd.DataFrame, target: pd.DataFrame) -> np.ndarray:
    overlaps = np.zeros(len(query), dtype=bool)
    query = query.reset_index(drop=True)
    for chrom, qsub in query.groupby("chrom", sort=False):
        tsub = target[target["chrom"].eq(chrom)].sort_values("start")
        if tsub.empty:
            continue
        starts = tsub["start"].to_numpy(dtype=np.int64)
        ends = tsub["end"].to_numpy(dtype=np.int64)
        order = np.argsort(starts)
        starts = starts[order]
        ends = ends[order]
        max_end = np.maximum.accumulate(ends)
        for idx, start, end in zip(
            qsub.index.to_numpy(dtype=np.int64),
            qsub["start"].to_numpy(dtype=np.int64),
            qsub["end"].to_numpy(dtype=np.int64),
        ):
            last_possible = np.searchsorted(starts, end, side="left") - 1
            if last_possible >= 0 and max_end[last_possible] > start:
                overlaps[idx] = True
    return overlaps


def active_tss_from_atac(tss: pd.DataFrame, atac: pd.DataFrame, window: int) -> pd.DataFrame:
    promoters = tss.copy()
    promoters["start"] = (promoters["tss"] - window).clip(lower=0)
    promoters["end"] = promoters["tss"] + window + 1
    promoters["active_by_atac_pm_window"] = interval_overlap_any(promoters, atac[["chrom", "start", "end"]])
    return tss.loc[promoters["active_by_atac_pm_window"].to_numpy()].copy().reset_index(drop=True)


def order_label(coords: list[int], labels: list[str]) -> str:
    items = sorted(zip(coords, labels), key=lambda item: (item[0], item[1]))
    groups: list[list[object]] = []
    for coord, label in items:
        if groups and groups[-1][0] == coord:
            groups[-1][1].append(label)
        else:
            groups.append([coord, [label]])
    return " -> ".join("=".join(group[1]) for group in groups)


def side_class(vals: list[int]) -> str:
    neg = sum(val < 0 for val in vals)
    pos = sum(val > 0 for val in vals)
    zero = sum(val == 0 for val in vals)
    if neg == 3:
        return "all three upstream of promoter"
    if pos == 3:
        return "all three downstream of promoter"
    if zero == 3:
        return "all three at promoter"
    if neg > 0 and pos > 0:
        return "promoter between factors"
    if zero > 0 and neg > 0 and pos == 0:
        return "some at promoter, rest upstream"
    if zero > 0 and pos > 0 and neg == 0:
        return "some at promoter, rest downstream"
    if zero > 0 and neg > 0 and pos > 0:
        return "promoter/tied within mixed cluster"
    return "other"


def make_triplets(xpc: pd.DataFrame, ctcf: pd.DataFrame, pol2: pd.DataFrame, tss: pd.DataFrame) -> pd.DataFrame:
    ctcf_idx, xpc_ctcf_dist = nearest(xpc, ctcf, "xpc_summit", "ctcf_summit")
    pol2_idx, xpc_pol2_dist = nearest(xpc, pol2, "xpc_summit", "pol2_summit")
    ctcf_match = ctcf.loc[ctcf_idx].reset_index(drop=True)
    pol2_match = pol2.loc[pol2_idx].reset_index(drop=True)

    triplets = xpc[["chrom", "xpc_id", "xpc_summit"]].copy().reset_index(drop=True)
    triplets["ctcf_id"] = ctcf_match["ctcf_id"].to_numpy()
    triplets["ctcf_summit"] = ctcf_match["ctcf_summit"].to_numpy()
    triplets["ctcf_motif_strand"] = ctcf_match["ctcf_motif_strand"].to_numpy()
    triplets["ctcf_motif_center"] = ctcf_match["ctcf_motif_center"].to_numpy()
    triplets["ctcf_motif_score"] = ctcf_match["ctcf_motif_score"].to_numpy()
    triplets["ctcf_motif_dist_bp"] = ctcf_match["ctcf_motif_dist_bp"].to_numpy()
    triplets["ctcf_has_overlapping_motif"] = ctcf_match["ctcf_has_overlapping_motif"].to_numpy()
    triplets["ctcf_high_conf_motif"] = ctcf_match["ctcf_high_conf_motif"].to_numpy()
    triplets["pol2_id"] = pol2_match["pol2_id"].to_numpy()
    triplets["pol2_summit"] = pol2_match["pol2_summit"].to_numpy()
    triplets["xpc_to_ctcf_bp"] = xpc_ctcf_dist.astype(int)
    triplets["xpc_to_pol2_bp"] = xpc_pol2_dist.astype(int)
    triplets["ctcf_to_pol2_bp"] = np.abs(
        triplets["ctcf_summit"].to_numpy(dtype=np.int64) - triplets["pol2_summit"].to_numpy(dtype=np.int64)
    )
    triplets["span_start"] = np.minimum.reduce(
        [
            triplets["xpc_summit"].to_numpy(dtype=np.int64),
            triplets["ctcf_summit"].to_numpy(dtype=np.int64),
            triplets["pol2_summit"].to_numpy(dtype=np.int64),
        ]
    )
    triplets["span_end"] = np.maximum.reduce(
        [
            triplets["xpc_summit"].to_numpy(dtype=np.int64),
            triplets["ctcf_summit"].to_numpy(dtype=np.int64),
            triplets["pol2_summit"].to_numpy(dtype=np.int64),
        ]
    )
    triplets["triplet_span_bp"] = triplets["span_end"] - triplets["span_start"]
    triplets["triplet_midpoint"] = ((triplets["span_start"] + triplets["span_end"]) // 2).astype(int)

    tss_idx, tss_dist = nearest(triplets, tss, "triplet_midpoint", "tss")
    tss_match = tss.loc[tss_idx].reset_index(drop=True)
    triplets["nearest_active_tss_idx"] = tss_idx
    triplets["nearest_active_tss"] = tss_match["tss"].to_numpy()
    triplets["nearest_active_gene"] = tss_match["gene"].to_numpy()
    triplets["nearest_active_tss_strand"] = tss_match["strand"].to_numpy()
    triplets["triplet_midpoint_to_active_tss_abs_bp"] = tss_dist.astype(int)

    sign = np.where(triplets["nearest_active_tss_strand"].eq("-").to_numpy(), -1, 1)
    for label in ["xpc", "ctcf", "pol2"]:
        triplets[f"{label}_signed_from_active_promoter_tx_bp"] = sign * (
            triplets[f"{label}_summit"].to_numpy(dtype=np.int64)
            - triplets["nearest_active_tss"].to_numpy(dtype=np.int64)
        )

    motif_sign = np.where(triplets["ctcf_motif_strand"].eq("-").to_numpy(), -1, 1)
    for label in ["xpc", "pol2"]:
        triplets[f"{label}_signed_from_ctcf_motif_bp"] = motif_sign * (
            triplets[f"{label}_summit"].to_numpy(dtype=np.int64)
            - triplets["ctcf_summit"].to_numpy(dtype=np.int64)
        )
    triplets["promoter_signed_from_ctcf_motif_bp"] = motif_sign * (
        triplets["nearest_active_tss"].to_numpy(dtype=np.int64)
        - triplets["ctcf_summit"].to_numpy(dtype=np.int64)
    )
    triplets["ctcf_signed_from_ctcf_motif_bp"] = 0

    triplets["cluster_order_tx_active_tss"] = [
        order_label(
            [
                row.xpc_signed_from_active_promoter_tx_bp,
                row.ctcf_signed_from_active_promoter_tx_bp,
                row.pol2_signed_from_active_promoter_tx_bp,
                0,
            ],
            ["XPC", "CTCF", "PolII", "Promoter"],
        )
        for row in triplets.itertuples(index=False)
    ]
    triplets["cluster_order_tx_active_tss_ctcf_motif_label"] = [
        order_label(
            [
                row.xpc_signed_from_active_promoter_tx_bp,
                row.ctcf_signed_from_active_promoter_tx_bp,
                row.pol2_signed_from_active_promoter_tx_bp,
                0,
            ],
            ["XPC", f"CTCF({row.ctcf_motif_strand})", "PolII", "Promoter"],
        )
        for row in triplets.itertuples(index=False)
    ]
    triplets["cluster_order_ctcf_motif_normalized"] = [
        order_label(
            [
                row.xpc_signed_from_ctcf_motif_bp,
                row.ctcf_signed_from_ctcf_motif_bp,
                row.pol2_signed_from_ctcf_motif_bp,
                row.promoter_signed_from_ctcf_motif_bp,
            ],
            ["XPC", "CTCF", "PolII", "Promoter"],
        )
        for row in triplets.itertuples(index=False)
    ]
    triplets["factor_order_tx_no_promoter_active_tss"] = [
        order_label(
            [
                row.xpc_signed_from_active_promoter_tx_bp,
                row.ctcf_signed_from_active_promoter_tx_bp,
                row.pol2_signed_from_active_promoter_tx_bp,
            ],
            ["XPC", "CTCF", "PolII"],
        )
        for row in triplets.itertuples(index=False)
    ]
    triplets["active_promoter_side_class"] = [
        side_class(
            [
                row.xpc_signed_from_active_promoter_tx_bp,
                row.ctcf_signed_from_active_promoter_tx_bp,
                row.pol2_signed_from_active_promoter_tx_bp,
            ]
        )
        for row in triplets.itertuples(index=False)
    ]
    return triplets


def summarize_triplets(triplets: pd.DataFrame, cutoffs: list[int]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary_rows = []
    order_rows = []
    side_rows = []
    for cutoff in cutoffs:
        mask = (
            (triplets["xpc_to_ctcf_bp"] <= cutoff)
            & (triplets["xpc_to_pol2_bp"] <= cutoff)
            & (triplets["ctcf_to_pol2_bp"] <= cutoff)
        )
        subset = triplets[mask].copy()
        n = len(subset)
        if n == 0:
            continue
        abs_mid = subset["triplet_midpoint_to_active_tss_abs_bp"].to_numpy(dtype=float)
        factor_abs_min = np.minimum.reduce(
            [
                np.abs(subset["xpc_signed_from_active_promoter_tx_bp"].to_numpy(dtype=float)),
                np.abs(subset["ctcf_signed_from_active_promoter_tx_bp"].to_numpy(dtype=float)),
                np.abs(subset["pol2_signed_from_active_promoter_tx_bp"].to_numpy(dtype=float)),
            ]
        )
        factor_abs_max = np.maximum.reduce(
            [
                np.abs(subset["xpc_signed_from_active_promoter_tx_bp"].to_numpy(dtype=float)),
                np.abs(subset["ctcf_signed_from_active_promoter_tx_bp"].to_numpy(dtype=float)),
                np.abs(subset["pol2_signed_from_active_promoter_tx_bp"].to_numpy(dtype=float)),
            ]
        )
        summary_rows.append(
            {
                "pairwise_cutoff_bp": cutoff,
                "triplets_all_three_pairwise_within_cutoff": n,
                "triplets_with_ctcf_high_conf_motif": int(subset["ctcf_high_conf_motif"].sum()),
                "ctcf_high_conf_motif_pct": 100 * np.mean(subset["ctcf_high_conf_motif"]),
                "ctcf_motif_plus_pct": 100 * np.mean(subset["ctcf_motif_strand"].eq("+")),
                "ctcf_motif_minus_pct": 100 * np.mean(subset["ctcf_motif_strand"].eq("-")),
                "median_triplet_span_bp": float(np.median(subset["triplet_span_bp"])),
                "median_triplet_midpoint_to_active_promoter_bp": float(np.median(abs_mid)),
                "pct_midpoint_within_1kb_active_promoter": 100 * np.mean(abs_mid <= 1000),
                "pct_midpoint_within_2kb_active_promoter": 100 * np.mean(abs_mid <= 2000),
                "pct_midpoint_within_5kb_active_promoter": 100 * np.mean(abs_mid <= 5000),
                "pct_midpoint_within_10kb_active_promoter": 100 * np.mean(abs_mid <= 10000),
                "pct_any_factor_within_1kb_active_promoter": 100 * np.mean(factor_abs_min <= 1000),
                "pct_all_factors_within_1kb_active_promoter": 100 * np.mean(factor_abs_max <= 1000),
                "pct_active_promoter_between_factors": 100
                * np.mean(subset["active_promoter_side_class"].eq("promoter between factors")),
                "pct_all_upstream_of_active_promoter": 100
                * np.mean(subset["active_promoter_side_class"].eq("all three upstream of promoter")),
                "pct_all_downstream_of_active_promoter": 100
                * np.mean(subset["active_promoter_side_class"].eq("all three downstream of promoter")),
            }
        )

        order_specs = [
            (
                "tx_active_tss",
                "Promoter",
                subset,
                "cluster_order_tx_active_tss",
            ),
            (
                "tx_active_tss_ctcf_motif_label",
                "Promoter",
                subset[subset["ctcf_high_conf_motif"]],
                "cluster_order_tx_active_tss_ctcf_motif_label",
            ),
            (
                "ctcf_motif_normalized",
                "CTCF",
                subset[subset["ctcf_high_conf_motif"]],
                "cluster_order_ctcf_motif_normalized",
            ),
        ]
        for coordinate_system, reference_landmark, order_subset, order_col in order_specs:
            denom = len(order_subset)
            if denom == 0:
                continue
            orders = order_subset[order_col].value_counts().rename_axis("order").reset_index(name="n")
            orders["pct"] = orders["n"] / denom * 100
            orders["rank"] = np.arange(1, len(orders) + 1)
            orders.insert(0, "pairwise_cutoff_bp", cutoff)
            orders.insert(1, "denominator_triplets", denom)
            orders.insert(2, "order_coordinate_system", coordinate_system)
            orders.insert(3, "reference_landmark", reference_landmark)
            orders = orders[
                [
                    "pairwise_cutoff_bp",
                    "denominator_triplets",
                    "order_coordinate_system",
                    "reference_landmark",
                    "rank",
                    "order",
                    "n",
                    "pct",
                ]
            ]
            order_rows.append(orders)

        sides = (
            subset["active_promoter_side_class"]
            .value_counts()
            .rename_axis("active_promoter_side_class")
            .reset_index(name="n")
        )
        sides["pct"] = sides["n"] / n * 100
        sides.insert(0, "pairwise_cutoff_bp", cutoff)
        sides.insert(1, "denominator_triplets", n)
        side_rows.append(sides)

    return (
        pd.DataFrame(summary_rows),
        pd.concat(order_rows, ignore_index=True),
        pd.concat(side_rows, ignore_index=True),
    )


def descriptive_orders_with_distances(triplets: pd.DataFrame, cutoffs: list[int]) -> pd.DataFrame:
    records = []
    for cutoff in cutoffs:
        mask = (
            (triplets["xpc_to_ctcf_bp"] <= cutoff)
            & (triplets["xpc_to_pol2_bp"] <= cutoff)
            & (triplets["ctcf_to_pol2_bp"] <= cutoff)
        )
        subset = triplets[mask].copy()
        order_specs = [
            ("tx_active_tss", "Promoter", subset),
            ("tx_active_tss_ctcf_motif_label", "Promoter", subset[subset["ctcf_high_conf_motif"]]),
            ("ctcf_motif_normalized", "CTCF", subset[subset["ctcf_high_conf_motif"]]),
        ]
        for coordinate_system, reference_landmark, order_subset in order_specs:
            denominator = len(order_subset)
            if denominator == 0:
                continue
            for row in order_subset.itertuples(index=False):
                if coordinate_system == "ctcf_motif_normalized":
                    items = [
                        ("Promoter", int(row.promoter_signed_from_ctcf_motif_bp)),
                        ("PolII", int(row.pol2_signed_from_ctcf_motif_bp)),
                        ("XPC", int(row.xpc_signed_from_ctcf_motif_bp)),
                        ("CTCF", 0),
                    ]
                elif coordinate_system == "tx_active_tss_ctcf_motif_label":
                    items = [
                        ("Promoter", 0),
                        ("PolII", int(row.pol2_signed_from_active_promoter_tx_bp)),
                        ("XPC", int(row.xpc_signed_from_active_promoter_tx_bp)),
                        (f"CTCF({row.ctcf_motif_strand})", int(row.ctcf_signed_from_active_promoter_tx_bp)),
                    ]
                else:
                    items = [
                        ("Promoter", 0),
                        ("PolII", int(row.pol2_signed_from_active_promoter_tx_bp)),
                        ("XPC", int(row.xpc_signed_from_active_promoter_tx_bp)),
                        ("CTCF", int(row.ctcf_signed_from_active_promoter_tx_bp)),
                    ]

                items = sorted(items, key=lambda item: (item[1], item[0]))
                groups = []
                for label, coord in items:
                    if groups and groups[-1]["coord"] == coord:
                        groups[-1]["labels"].append(label)
                    else:
                        groups.append({"coord": coord, "labels": [label]})
                for group in groups:
                    group["label"] = "=".join(group["labels"])

                rec = {
                    "pairwise_cutoff_bp": cutoff,
                    "order_coordinate_system": coordinate_system,
                    "reference_landmark": reference_landmark,
                    "denominator_triplets": denominator,
                    "order": " -> ".join(group["label"] for group in groups),
                    "edge0_label": groups[0]["label"],
                    "cluster_width_bp": groups[-1]["coord"] - groups[0]["coord"],
                    "reference_to_nearest_other_landmark_bp": min(
                        abs(coord) for label, coord in items if reference_landmark not in label.split("=")
                    ),
                }
                for edge_idx in range(1, 4):
                    if edge_idx < len(groups):
                        rec[f"edge{edge_idx}_label"] = f"{groups[edge_idx - 1]['label']} -> {groups[edge_idx]['label']}"
                        rec[f"edge{edge_idx}_gap_bp"] = groups[edge_idx]["coord"] - groups[edge_idx - 1]["coord"]
                    else:
                        rec[f"edge{edge_idx}_label"] = pd.NA
                        rec[f"edge{edge_idx}_gap_bp"] = np.nan
                records.append(rec)

    edges = pd.DataFrame(records)
    grouped = edges.groupby(["pairwise_cutoff_bp", "order_coordinate_system", "reference_landmark", "order"], sort=False).agg(
        denominator_triplets=("denominator_triplets", "first"),
        edge0_label=("edge0_label", "first"),
        edge1_label=("edge1_label", "first"),
        edge2_label=("edge2_label", "first"),
        edge3_label=("edge3_label", "first"),
        n=("order", "size"),
        edge1_mean_bp=("edge1_gap_bp", "mean"),
        edge1_median_bp=("edge1_gap_bp", "median"),
        edge2_mean_bp=("edge2_gap_bp", "mean"),
        edge2_median_bp=("edge2_gap_bp", "median"),
        edge3_mean_bp=("edge3_gap_bp", "mean"),
        edge3_median_bp=("edge3_gap_bp", "median"),
        cluster_width_mean_bp=("cluster_width_bp", "mean"),
        cluster_width_median_bp=("cluster_width_bp", "median"),
        reference_to_nearest_other_landmark_mean_bp=("reference_to_nearest_other_landmark_bp", "mean"),
        reference_to_nearest_other_landmark_median_bp=("reference_to_nearest_other_landmark_bp", "median"),
    )
    grouped = grouped.reset_index()
    grouped["pct"] = grouped["n"] / grouped["denominator_triplets"] * 100
    grouped = grouped.sort_values(
        ["pairwise_cutoff_bp", "order_coordinate_system", "n", "order"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)
    grouped["rank"] = grouped.groupby(["pairwise_cutoff_bp", "order_coordinate_system"]).cumcount() + 1
    for col in [col for col in grouped.columns if col.endswith("_bp")]:
        grouped[col] = grouped[col].round(2)
    grouped["pct"] = grouped["pct"].round(2)
    grouped["descriptive_order_with_mean_median_gaps"] = grouped.apply(format_descriptive_order, axis=1)
    return grouped[
        [
            "pairwise_cutoff_bp",
            "order_coordinate_system",
            "reference_landmark",
            "rank",
            "order",
            "n",
            "denominator_triplets",
            "pct",
            "descriptive_order_with_mean_median_gaps",
            "edge1_label",
            "edge1_mean_bp",
            "edge1_median_bp",
            "edge2_label",
            "edge2_mean_bp",
            "edge2_median_bp",
            "edge3_label",
            "edge3_mean_bp",
            "edge3_median_bp",
            "cluster_width_mean_bp",
            "cluster_width_median_bp",
            "reference_to_nearest_other_landmark_mean_bp",
            "reference_to_nearest_other_landmark_median_bp",
        ]
    ]


def format_bp(value: object) -> str:
    if pd.isna(value):
        return "NA"
    return f"{int(round(float(value))):,} bp"


def format_descriptive_order(row: pd.Series) -> str:
    parts = [str(row["edge0_label"])]
    for edge_idx in range(1, 4):
        edge_label = row.get(f"edge{edge_idx}_label")
        if pd.isna(edge_label) or edge_label == "":
            continue
        target = str(edge_label).split(" -> ")[1]
        parts.append(
            f"--[mean {format_bp(row.get(f'edge{edge_idx}_mean_bp'))}; "
            f"median {format_bp(row.get(f'edge{edge_idx}_median_bp'))}]--> {target}"
        )
    return " ".join(parts)


def write_active_tss_bed(active_tss: pd.DataFrame, path: Path) -> None:
    bed = active_tss[["chrom", "start", "end", "gene", "score", "strand"]].copy()
    bed.to_csv(path, sep="\t", header=False, index=False)


def main() -> None:
    args = parse_args()
    cutoffs = [100, 250, 500, 1000, 2000, 5000]

    xpc = read_peaks(args.xpc, "xpc")
    ctcf = read_peaks(args.ctcf, "ctcf")
    motifs = read_ctcf_motifs(args.ctcf_motifs)
    ctcf = assign_ctcf_motifs(ctcf, motifs, args.ctcf_motif_max_dist)
    pol2 = read_peaks(args.pol2, "pol2")
    tss = read_tss(args.tss)
    atac = read_peaks(args.atac, "atac")

    active_tss = active_tss_from_atac(tss, atac, args.promoter_window)
    write_active_tss_bed(active_tss, args.out_prefix.with_name(f"{args.out_prefix.name}_active_tss_pm{args.promoter_window}.bed"))

    triplets = make_triplets(xpc, ctcf, pol2, active_tss)
    all_triplets_n = len(triplets)
    triplets_before_motif_filter = len(triplets)
    if args.require_ctcf_motif:
        triplets = triplets[triplets["ctcf_high_conf_motif"]].copy().reset_index(drop=True)
    triplets_after_motif_filter = len(triplets)
    triplets_before_pol2_filter = len(triplets)
    if args.pol2_tss_window is not None:
        triplets = triplets[
            triplets["pol2_signed_from_active_promoter_tx_bp"].abs() <= args.pol2_tss_window
        ].copy().reset_index(drop=True)
    triplets_after_pol2_filter = len(triplets)

    summary, orders, sides = summarize_triplets(triplets, cutoffs)
    distance_orders = descriptive_orders_with_distances(triplets, cutoffs)

    triplets.to_csv(args.out_prefix.with_name(f"{args.out_prefix.name}_triplets_all_xpc_anchored.tsv"), sep="\t", index=False)
    summary.to_csv(args.out_prefix.with_name(f"{args.out_prefix.name}_triplets_promoter_summary.tsv"), sep="\t", index=False)
    orders.to_csv(args.out_prefix.with_name(f"{args.out_prefix.name}_triplets_promoter_orders.tsv"), sep="\t", index=False)
    sides.to_csv(args.out_prefix.with_name(f"{args.out_prefix.name}_triplets_promoter_side_classes.tsv"), sep="\t", index=False)
    distance_orders.to_csv(
        args.out_prefix.with_name(f"{args.out_prefix.name}_triplets_promoter_orders_with_distances.tsv"),
        sep="\t",
        index=False,
    )

    active_summary = pd.DataFrame(
        [
            {
                "tss_total_chr1_22_x": len(tss),
                "atac_peaks_chr1_22_x": len(atac),
                "active_tss_definition": f"TSS +/- {args.promoter_window} bp overlaps ATAC peak",
                "active_tss": len(active_tss),
                "active_tss_pct": 100 * len(active_tss) / len(tss),
                "ctcf_peaks_chr1_22_x": len(ctcf),
                "ctcf_peaks_with_overlapping_motif": int(ctcf["ctcf_has_overlapping_motif"].sum()),
                "ctcf_peaks_with_overlapping_motif_pct": 100 * np.mean(ctcf["ctcf_has_overlapping_motif"]),
                "ctcf_peaks_with_high_conf_motif": int(ctcf["ctcf_high_conf_motif"].sum()),
                "ctcf_peaks_with_high_conf_motif_pct": 100 * np.mean(ctcf["ctcf_high_conf_motif"]),
                "ctcf_motif_max_dist_bp": args.ctcf_motif_max_dist,
                "require_ctcf_motif": args.require_ctcf_motif,
                "xpc_anchored_triplets_before_ctcf_motif_filter": triplets_before_motif_filter,
                "xpc_anchored_triplets_after_ctcf_motif_filter": triplets_after_motif_filter,
                "xpc_anchored_triplets_kept_after_ctcf_motif_filter_pct": 100
                * triplets_after_motif_filter
                / triplets_before_motif_filter,
                "pol2_tss_filter_bp": args.pol2_tss_window if args.pol2_tss_window is not None else "none",
                "xpc_anchored_triplets_before_pol2_tss_filter": triplets_before_pol2_filter,
                "xpc_anchored_triplets_after_pol2_tss_filter": triplets_after_pol2_filter,
                "xpc_anchored_triplets_kept_after_pol2_tss_filter_pct": 100
                * triplets_after_pol2_filter
                / triplets_before_pol2_filter,
                "xpc_anchored_triplets_final_kept_pct": 100 * len(triplets) / all_triplets_n,
            }
        ]
    )
    active_summary.to_csv(args.out_prefix.with_name(f"{args.out_prefix.name}_active_tss_summary.tsv"), sep="\t", index=False)

    print("Active TSS definition:", active_summary.loc[0, "active_tss_definition"])
    print(f"TSS chr1-22,X: {len(tss):,}")
    print(f"ATAC peaks chr1-22,X: {len(atac):,}")
    print(f"Active TSS: {len(active_tss):,} ({100 * len(active_tss) / len(tss):.1f}%)")
    print(
        f"CTCF peaks with high-confidence motif <= {args.ctcf_motif_max_dist:,} bp from summit: "
        f"{ctcf['ctcf_high_conf_motif'].sum():,} / {len(ctcf):,} "
        f"({100 * np.mean(ctcf['ctcf_high_conf_motif']):.1f}%)"
    )
    if args.require_ctcf_motif:
        print(
            f"CTCF motif filter: kept {triplets_after_motif_filter:,} / {triplets_before_motif_filter:,} "
            f"XPC-anchored triplets ({100 * triplets_after_motif_filter / triplets_before_motif_filter:.1f}%)"
        )
    if args.pol2_tss_window is not None:
        print(
            f"Pol II promoter-proximal filter: abs(PolII - active TSS, tx-oriented) <= "
            f"{args.pol2_tss_window:,} bp"
        )
        print(
            f"XPC-anchored triplets kept after Pol II filter: {triplets_after_pol2_filter:,} / "
            f"{triplets_before_pol2_filter:,} "
            f"({100 * triplets_after_pol2_filter / triplets_before_pol2_filter:.1f}%)"
        )
    print("\nTriplet summary using nearest active TSS/promoter:")
    print(summary.to_string(index=False))
    print("\nTop orders for pairwise cutoff <= 5000 bp:")
    print(distance_orders[distance_orders["pairwise_cutoff_bp"].eq(5000)].head(10).to_string(index=False))


if __name__ == "__main__":
    main()

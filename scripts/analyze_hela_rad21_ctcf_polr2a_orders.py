#!/usr/bin/env python3
"""Analyze HeLa RAD21/CTCF/POLR2A order relative to POLR2A-active TSSs."""

from __future__ import annotations

import argparse
from pathlib import Path

import bioframe as bf
import numpy as np
import pandas as pd


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
        description=(
            "RAD21-anchored HeLa CTCF/RAD21/POLR2A order analysis. "
            "Active TSSs are inferred by proximity to POLR2A peaks."
        )
    )
    parser.add_argument("--rad21", type=Path, default=Path("data/hela_order/ENCFF239FBO.bed"))
    parser.add_argument("--ctcf", type=Path, default=Path("data/hela_order/ENCFF502CZS.bed"))
    parser.add_argument("--polr2a", type=Path, default=Path("data/hela_order/ENCFF246QVY.bed"))
    parser.add_argument("--tss", type=Path, default=Path("misc/tss.bed"))
    parser.add_argument("--ctcf-motifs", type=Path, default=Path("misc/CTCF_hg38.bed"))
    parser.add_argument(
        "--ctcf-motif-max-dist",
        type=int,
        default=100,
        help="Maximum CTCF summit-to-motif-center distance for high-confidence orientation.",
    )
    parser.add_argument(
        "--require-ctcf-motif",
        action="store_true",
        help="Retain only triplets whose nearest CTCF peak has a high-confidence motif orientation.",
    )
    parser.add_argument(
        "--polr2a-tss-window",
        type=int,
        default=1000,
        help=(
            "Call TSSs active when their nearest POLR2A summit is within this many bp; "
            "also require the triplet POLR2A summit to be within this distance of the assigned active TSS."
        ),
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("data/hela_order/hela_rad21_ctcf_polr2a_polr2aprox1kb"),
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


def active_tss_from_polr2a(tss: pd.DataFrame, polr2a: pd.DataFrame, window: int) -> pd.DataFrame:
    polr2a_idx, polr2a_dist = nearest(tss, polr2a, "tss", "polr2a_summit")
    active = tss.copy()
    active["nearest_polr2a_id"] = polr2a_idx
    active["nearest_polr2a_abs_bp"] = polr2a_dist.astype(int)
    active = active[active["nearest_polr2a_abs_bp"].le(window)].copy().reset_index(drop=True)
    return active


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


def make_triplets(
    rad21: pd.DataFrame,
    ctcf: pd.DataFrame,
    polr2a: pd.DataFrame,
    active_tss: pd.DataFrame,
) -> pd.DataFrame:
    ctcf_idx, rad21_ctcf_dist = nearest(rad21, ctcf, "rad21_summit", "ctcf_summit")
    polr2a_idx, rad21_polr2a_dist = nearest(rad21, polr2a, "rad21_summit", "polr2a_summit")
    ctcf_match = ctcf.loc[ctcf_idx].reset_index(drop=True)
    polr2a_match = polr2a.loc[polr2a_idx].reset_index(drop=True)

    triplets = rad21[["chrom", "rad21_id", "rad21_summit"]].copy().reset_index(drop=True)
    triplets["ctcf_id"] = ctcf_match["ctcf_id"].to_numpy()
    triplets["ctcf_summit"] = ctcf_match["ctcf_summit"].to_numpy()
    triplets["ctcf_motif_strand"] = ctcf_match["ctcf_motif_strand"].to_numpy()
    triplets["ctcf_motif_center"] = ctcf_match["ctcf_motif_center"].to_numpy()
    triplets["ctcf_motif_score"] = ctcf_match["ctcf_motif_score"].to_numpy()
    triplets["ctcf_motif_dist_bp"] = ctcf_match["ctcf_motif_dist_bp"].to_numpy()
    triplets["ctcf_has_overlapping_motif"] = ctcf_match["ctcf_has_overlapping_motif"].to_numpy()
    triplets["ctcf_high_conf_motif"] = ctcf_match["ctcf_high_conf_motif"].to_numpy()
    triplets["polr2a_id"] = polr2a_match["polr2a_id"].to_numpy()
    triplets["polr2a_summit"] = polr2a_match["polr2a_summit"].to_numpy()
    triplets["rad21_to_ctcf_bp"] = rad21_ctcf_dist.astype(int)
    triplets["rad21_to_polr2a_bp"] = rad21_polr2a_dist.astype(int)
    triplets["ctcf_to_polr2a_bp"] = np.abs(
        triplets["ctcf_summit"].to_numpy(dtype=np.int64)
        - triplets["polr2a_summit"].to_numpy(dtype=np.int64)
    )
    triplets["span_start"] = np.minimum.reduce(
        [
            triplets["rad21_summit"].to_numpy(dtype=np.int64),
            triplets["ctcf_summit"].to_numpy(dtype=np.int64),
            triplets["polr2a_summit"].to_numpy(dtype=np.int64),
        ]
    )
    triplets["span_end"] = np.maximum.reduce(
        [
            triplets["rad21_summit"].to_numpy(dtype=np.int64),
            triplets["ctcf_summit"].to_numpy(dtype=np.int64),
            triplets["polr2a_summit"].to_numpy(dtype=np.int64),
        ]
    )
    triplets["triplet_span_bp"] = triplets["span_end"] - triplets["span_start"]
    triplets["triplet_midpoint"] = ((triplets["span_start"] + triplets["span_end"]) // 2).astype(int)

    tss_idx, tss_dist = nearest(triplets, active_tss, "triplet_midpoint", "tss")
    tss_match = active_tss.loc[tss_idx].reset_index(drop=True)
    triplets["nearest_active_tss_idx"] = tss_idx
    triplets["nearest_active_tss"] = tss_match["tss"].to_numpy()
    triplets["nearest_active_gene"] = tss_match["gene"].to_numpy()
    triplets["nearest_active_tss_strand"] = tss_match["strand"].to_numpy()
    triplets["triplet_midpoint_to_active_tss_abs_bp"] = tss_dist.astype(int)

    sign = np.where(triplets["nearest_active_tss_strand"].eq("-").to_numpy(), -1, 1)
    for label in ["rad21", "ctcf", "polr2a"]:
        triplets[f"{label}_signed_from_active_promoter_tx_bp"] = sign * (
            triplets[f"{label}_summit"].to_numpy(dtype=np.int64)
            - triplets["nearest_active_tss"].to_numpy(dtype=np.int64)
        )

    motif_sign = np.where(triplets["ctcf_motif_strand"].eq("-").to_numpy(), -1, 1)
    for label in ["rad21", "polr2a"]:
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
                row.rad21_signed_from_active_promoter_tx_bp,
                row.ctcf_signed_from_active_promoter_tx_bp,
                row.polr2a_signed_from_active_promoter_tx_bp,
                0,
            ],
            ["RAD21", "CTCF", "POLR2A", "Promoter"],
        )
        for row in triplets.itertuples(index=False)
    ]
    triplets["cluster_order_tx_active_tss_ctcf_motif_label"] = [
        order_label(
            [
                row.rad21_signed_from_active_promoter_tx_bp,
                row.ctcf_signed_from_active_promoter_tx_bp,
                row.polr2a_signed_from_active_promoter_tx_bp,
                0,
            ],
            ["RAD21", f"CTCF({row.ctcf_motif_strand})", "POLR2A", "Promoter"],
        )
        for row in triplets.itertuples(index=False)
    ]
    triplets["cluster_order_ctcf_motif_normalized"] = [
        order_label(
            [
                row.rad21_signed_from_ctcf_motif_bp,
                row.ctcf_signed_from_ctcf_motif_bp,
                row.polr2a_signed_from_ctcf_motif_bp,
                row.promoter_signed_from_ctcf_motif_bp,
            ],
            ["RAD21", "CTCF", "POLR2A", "Promoter"],
        )
        for row in triplets.itertuples(index=False)
    ]
    triplets["active_promoter_side_class"] = [
        side_class(
            [
                row.rad21_signed_from_active_promoter_tx_bp,
                row.ctcf_signed_from_active_promoter_tx_bp,
                row.polr2a_signed_from_active_promoter_tx_bp,
            ]
        )
        for row in triplets.itertuples(index=False)
    ]
    return triplets


def cutoff_mask(triplets: pd.DataFrame, cutoff: int) -> pd.Series:
    return (
        (triplets["rad21_to_ctcf_bp"] <= cutoff)
        & (triplets["rad21_to_polr2a_bp"] <= cutoff)
        & (triplets["ctcf_to_polr2a_bp"] <= cutoff)
    )


def pairwise_colocalization_summary(triplets: pd.DataFrame, cutoffs: list[int]) -> pd.DataFrame:
    rows = []
    denom = len(triplets)
    for cutoff in cutoffs:
        rad21_ctcf = triplets["rad21_to_ctcf_bp"].le(cutoff)
        rad21_polr2a = triplets["rad21_to_polr2a_bp"].le(cutoff)
        ctcf_polr2a = triplets["ctcf_to_polr2a_bp"].le(cutoff)
        rows.append(
            {
                "cutoff_bp": cutoff,
                "denominator_rad21_anchored_triplets": denom,
                "rad21_near_ctcf": int(rad21_ctcf.sum()),
                "rad21_near_ctcf_pct": 100 * rad21_ctcf.mean(),
                "rad21_near_polr2a": int(rad21_polr2a.sum()),
                "rad21_near_polr2a_pct": 100 * rad21_polr2a.mean(),
                "rad21_near_both_ctcf_and_polr2a": int((rad21_ctcf & rad21_polr2a).sum()),
                "rad21_near_both_ctcf_and_polr2a_pct": 100 * (rad21_ctcf & rad21_polr2a).mean(),
                "all_three_pairwise": int((rad21_ctcf & rad21_polr2a & ctcf_polr2a).sum()),
                "all_three_pairwise_pct": 100 * (rad21_ctcf & rad21_polr2a & ctcf_polr2a).mean(),
            }
        )
    return pd.DataFrame(rows)


def summarize_triplets(triplets: pd.DataFrame, cutoffs: list[int]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary_rows = []
    order_rows = []
    side_rows = []
    for cutoff in cutoffs:
        subset = triplets[cutoff_mask(triplets, cutoff)].copy()
        n = len(subset)
        if n == 0:
            continue
        abs_mid = subset["triplet_midpoint_to_active_tss_abs_bp"].to_numpy(dtype=float)
        factor_abs_min = np.minimum.reduce(
            [
                np.abs(subset["rad21_signed_from_active_promoter_tx_bp"].to_numpy(dtype=float)),
                np.abs(subset["ctcf_signed_from_active_promoter_tx_bp"].to_numpy(dtype=float)),
                np.abs(subset["polr2a_signed_from_active_promoter_tx_bp"].to_numpy(dtype=float)),
            ]
        )
        factor_abs_max = np.maximum.reduce(
            [
                np.abs(subset["rad21_signed_from_active_promoter_tx_bp"].to_numpy(dtype=float)),
                np.abs(subset["ctcf_signed_from_active_promoter_tx_bp"].to_numpy(dtype=float)),
                np.abs(subset["polr2a_signed_from_active_promoter_tx_bp"].to_numpy(dtype=float)),
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
                "median_rad21_to_ctcf_bp": float(np.median(subset["rad21_to_ctcf_bp"])),
                "median_rad21_to_polr2a_bp": float(np.median(subset["rad21_to_polr2a_bp"])),
                "median_ctcf_to_polr2a_bp": float(np.median(subset["ctcf_to_polr2a_bp"])),
                "median_abs_rad21_to_active_promoter_bp": float(
                    np.median(np.abs(subset["rad21_signed_from_active_promoter_tx_bp"]))
                ),
                "median_abs_ctcf_to_active_promoter_bp": float(
                    np.median(np.abs(subset["ctcf_signed_from_active_promoter_tx_bp"]))
                ),
                "median_abs_polr2a_to_active_promoter_bp": float(
                    np.median(np.abs(subset["polr2a_signed_from_active_promoter_tx_bp"]))
                ),
                "pct_midpoint_within_1kb_active_promoter": 100 * np.mean(abs_mid <= 1000),
                "pct_midpoint_within_2kb_active_promoter": 100 * np.mean(abs_mid <= 2000),
                "pct_midpoint_within_5kb_active_promoter": 100 * np.mean(abs_mid <= 5000),
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
            ("tx_active_tss", "Promoter", subset, "cluster_order_tx_active_tss"),
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
            order_rows.append(
                orders[
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
            )

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
        subset = triplets[cutoff_mask(triplets, cutoff)].copy()
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
                        ("POLR2A", int(row.polr2a_signed_from_ctcf_motif_bp)),
                        ("RAD21", int(row.rad21_signed_from_ctcf_motif_bp)),
                        ("CTCF", 0),
                    ]
                elif coordinate_system == "tx_active_tss_ctcf_motif_label":
                    items = [
                        ("Promoter", 0),
                        ("POLR2A", int(row.polr2a_signed_from_active_promoter_tx_bp)),
                        ("RAD21", int(row.rad21_signed_from_active_promoter_tx_bp)),
                        (f"CTCF({row.ctcf_motif_strand})", int(row.ctcf_signed_from_active_promoter_tx_bp)),
                    ]
                else:
                    items = [
                        ("Promoter", 0),
                        ("POLR2A", int(row.polr2a_signed_from_active_promoter_tx_bp)),
                        ("RAD21", int(row.rad21_signed_from_active_promoter_tx_bp)),
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
    grouped = edges.groupby(
        ["pairwise_cutoff_bp", "order_coordinate_system", "reference_landmark", "order"],
        sort=False,
    ).agg(
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


def write_outputs(
    out_prefix: Path,
    active_tss: pd.DataFrame,
    triplets: pd.DataFrame,
    summary: pd.DataFrame,
    orders: pd.DataFrame,
    sides: pd.DataFrame,
    distance_orders: pd.DataFrame,
    pairwise: pd.DataFrame,
    active_summary: pd.DataFrame,
    polr2a_window: int,
) -> None:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_active_tss_bed(active_tss, out_prefix.with_name(f"{out_prefix.name}_active_tss_polr2a_pm{polr2a_window}.bed"))
    triplets.to_csv(out_prefix.with_name(f"{out_prefix.name}_triplets_all_rad21_anchored.tsv"), sep="\t", index=False)
    summary.to_csv(out_prefix.with_name(f"{out_prefix.name}_triplets_promoter_summary.tsv"), sep="\t", index=False)
    orders.to_csv(out_prefix.with_name(f"{out_prefix.name}_triplets_promoter_orders.tsv"), sep="\t", index=False)
    sides.to_csv(out_prefix.with_name(f"{out_prefix.name}_triplets_promoter_side_classes.tsv"), sep="\t", index=False)
    distance_orders.to_csv(
        out_prefix.with_name(f"{out_prefix.name}_triplets_promoter_orders_with_distances.tsv"),
        sep="\t",
        index=False,
    )
    pairwise.to_csv(out_prefix.with_name(f"{out_prefix.name}_pairwise_colocalization_counts.tsv"), sep="\t", index=False)
    active_summary.to_csv(out_prefix.with_name(f"{out_prefix.name}_active_tss_summary.tsv"), sep="\t", index=False)


def main() -> None:
    args = parse_args()
    cutoffs = [100, 250, 500, 1000, 2000, 5000]

    rad21 = read_peaks(args.rad21, "rad21")
    ctcf = read_peaks(args.ctcf, "ctcf")
    polr2a = read_peaks(args.polr2a, "polr2a")
    tss = read_tss(args.tss)
    motifs = read_ctcf_motifs(args.ctcf_motifs)
    ctcf = assign_ctcf_motifs(ctcf, motifs, args.ctcf_motif_max_dist)
    active_tss = active_tss_from_polr2a(tss, polr2a, args.polr2a_tss_window)

    triplets = make_triplets(rad21, ctcf, polr2a, active_tss)
    all_triplets_n = len(triplets)
    triplets_before_motif_filter = len(triplets)
    if args.require_ctcf_motif:
        triplets = triplets[triplets["ctcf_high_conf_motif"]].copy().reset_index(drop=True)
    triplets_after_motif_filter = len(triplets)
    triplets_before_polr2a_filter = len(triplets)
    triplets = triplets[
        triplets["polr2a_signed_from_active_promoter_tx_bp"].abs() <= args.polr2a_tss_window
    ].copy().reset_index(drop=True)
    triplets_after_polr2a_filter = len(triplets)

    summary, orders, sides = summarize_triplets(triplets, cutoffs)
    distance_orders = descriptive_orders_with_distances(triplets, cutoffs)
    pairwise = pairwise_colocalization_summary(triplets, cutoffs)

    active_summary = pd.DataFrame(
        [
            {
                "tss_total_chr1_22_x": len(tss),
                "rad21_peaks_chr1_22_x": len(rad21),
                "ctcf_peaks_chr1_22_x": len(ctcf),
                "polr2a_peaks_chr1_22_x": len(polr2a),
                "active_tss_definition": f"nearest POLR2A summit within {args.polr2a_tss_window} bp",
                "active_tss": len(active_tss),
                "active_tss_pct": 100 * len(active_tss) / len(tss),
                "ctcf_peaks_with_overlapping_motif": int(ctcf["ctcf_has_overlapping_motif"].sum()),
                "ctcf_peaks_with_overlapping_motif_pct": 100 * np.mean(ctcf["ctcf_has_overlapping_motif"]),
                "ctcf_peaks_with_high_conf_motif": int(ctcf["ctcf_high_conf_motif"].sum()),
                "ctcf_peaks_with_high_conf_motif_pct": 100 * np.mean(ctcf["ctcf_high_conf_motif"]),
                "ctcf_motif_max_dist_bp": args.ctcf_motif_max_dist,
                "require_ctcf_motif": args.require_ctcf_motif,
                "rad21_anchored_triplets_before_ctcf_motif_filter": triplets_before_motif_filter,
                "rad21_anchored_triplets_after_ctcf_motif_filter": triplets_after_motif_filter,
                "rad21_anchored_triplets_kept_after_ctcf_motif_filter_pct": 100
                * triplets_after_motif_filter
                / triplets_before_motif_filter,
                "polr2a_tss_filter_bp": args.polr2a_tss_window,
                "rad21_anchored_triplets_before_polr2a_tss_filter": triplets_before_polr2a_filter,
                "rad21_anchored_triplets_after_polr2a_tss_filter": triplets_after_polr2a_filter,
                "rad21_anchored_triplets_kept_after_polr2a_tss_filter_pct": 100
                * triplets_after_polr2a_filter
                / triplets_before_polr2a_filter,
                "rad21_anchored_triplets_final_kept_pct": 100 * len(triplets) / all_triplets_n,
            }
        ]
    )

    write_outputs(
        args.out_prefix,
        active_tss,
        triplets,
        summary,
        orders,
        sides,
        distance_orders,
        pairwise,
        active_summary,
        args.polr2a_tss_window,
    )

    print("Active TSS definition:", active_summary.loc[0, "active_tss_definition"])
    print(f"TSS chr1-22,X: {len(tss):,}")
    print(f"RAD21 peaks chr1-22,X: {len(rad21):,}")
    print(f"CTCF peaks chr1-22,X: {len(ctcf):,}")
    print(f"POLR2A peaks chr1-22,X: {len(polr2a):,}")
    print(f"Active TSS: {len(active_tss):,} ({100 * len(active_tss) / len(tss):.1f}%)")
    print(
        f"CTCF peaks with high-confidence motif <= {args.ctcf_motif_max_dist:,} bp from summit: "
        f"{ctcf['ctcf_high_conf_motif'].sum():,} / {len(ctcf):,} "
        f"({100 * np.mean(ctcf['ctcf_high_conf_motif']):.1f}%)"
    )
    if args.require_ctcf_motif:
        print(
            f"CTCF motif filter: kept {triplets_after_motif_filter:,} / {triplets_before_motif_filter:,} "
            f"RAD21-anchored triplets ({100 * triplets_after_motif_filter / triplets_before_motif_filter:.1f}%)"
        )
    print(
        f"POLR2A promoter-proximal filter: abs(POLR2A - active TSS, tx-oriented) <= "
        f"{args.polr2a_tss_window:,} bp"
    )
    print(
        f"RAD21-anchored triplets kept after POLR2A filter: {triplets_after_polr2a_filter:,} / "
        f"{triplets_before_polr2a_filter:,} "
        f"({100 * triplets_after_polr2a_filter / triplets_before_polr2a_filter:.1f}%)"
    )
    print("\nPairwise colocalization among final RAD21-anchored triplets:")
    print(pairwise.to_string(index=False))
    print("\nTriplet summary using nearest POLR2A-active TSS/promoter:")
    print(summary.to_string(index=False))
    print("\nTop orders for pairwise cutoff <= 1000 bp:")
    print(distance_orders[distance_orders["pairwise_cutoff_bp"].eq(1000)].head(10).to_string(index=False))


if __name__ == "__main__":
    main()

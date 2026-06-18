#!/usr/bin/env python3
"""TS / NTS metagene lineplots around protein-coding TSSs for ds_bed timepoints.

For each timepoint BED (stranded BED6 nascent reads) build strand-resolved
metagene coverage profiles around TSSs, oriented by transcription direction:

    NTS (non-template / sense)  = reads on the SAME strand as the gene
    TS  (template / antisense)  = reads on the OPPOSITE strand

Window is K kb upstream and L kb downstream of the TSS (default 5 / 20 kb).
Signal is full read coverage (every bp a read spans), library-size normalized
to CPM so timepoints are comparable.

Run:
    uv run scripts/plot_ds_tss_ts_nts_metagene.py
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Reuse repo plotting style helpers.
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
from ggner_3d.plotting import update_rcparams, despine, COLORS, _smooth_1d  # noqa: E402

MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
MAIN_CHROM_SET = set(MAIN_CHROMS)

# GTF chrom -> BED chrom: Ensembl uses "1"/"X", BED uses "chr1"/"chrX".
NTS_COLOR = COLORS[1]  # sense
TS_COLOR = COLORS[0]   # antisense


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed-dir", type=Path, default=Path("data/ds_bed"))
    p.add_argument("--gtf", type=Path, default=Path("misc/Homo_sapiens.GRCh38.115.chr.gtf"))
    p.add_argument("--upstream", type=int, default=5000, help="K bp upstream of TSS.")
    p.add_argument("--downstream", type=int, default=20000, help="L bp downstream of TSS.")
    p.add_argument("--bin", type=int, default=100, help="Bin size (bp). Must divide K+L.")
    p.add_argument("--biotype", default="protein_coding", help="gene_biotype to keep ('all' for any).")
    p.add_argument("--out-dir", type=Path, default=Path("figs/ds_tss_metagene"))
    p.add_argument("--swap-strands", action="store_true", help="Flip NTS/NTS (reverse-stranded library).")
    p.add_argument("--smooth-sigma", type=float, default=1.0, help="Gaussian smoothing sigma in bins (0=off).")
    p.add_argument("--chroms", nargs="*", default=None, help="Restrict to these chroms (e.g. chr1) for quick runs.")
    p.add_argument("--signal", choices=["coverage", "5prime"], default="coverage",
                   help="'coverage' = full read span; '5prime' = single bp at read 5' end (damage-seq lesion site).")
    p.add_argument("--lesion-offset", type=int, default=0,
                   help="Shift 5' position upstream (strand-relative) by N bp to hit the lesion dinucleotide.")
    p.add_argument("--pol2-bw", type=Path, default=Path("misc/GSM7162996_NCRNAPolIIChIP.primary.ucsc.bw"),
                   help="PolII bigwig used to rank gene expression for filtering.")
    p.add_argument("--pol2-top-frac", type=float, default=None,
                   help="Keep only the top fraction of genes by promoter PolII signal (e.g. 0.5). None = no filter.")
    p.add_argument("--pol2-window", type=int, default=1000,
                   help="Promoter window centered on TSS (bp) for --pol2-mode promoter.")
    p.add_argument("--pol2-mode", choices=["promoter", "body", "pausing_index"], default="promoter",
                   help="What pre-UV PolII metric to rank genes by: promoter (paused), "
                        "body (elongating, the pausing-index denominator), or pausing_index.")
    p.add_argument("--pol2-min", type=float, default=None,
                   help="Keep only genes with PolII (in --pol2-mode units) >= this value. "
                        "Selects high-transcription / high-TCR genes to widen the 0h->3h gap.")
    p.add_argument("--pol2-max", type=float, default=None,
                   help="Keep only genes with PolII <= this value (e.g. exclude the top "
                        "body-ChIP turnover quantile).")
    p.add_argument("--pol2-quantiles", type=int, default=None,
                   help="Stratify genes into N PolII quantiles and emit dose-response figures "
                        "(supersedes --pol2-top-frac).")
    p.add_argument("--body-start", type=int, default=2000,
                   help="Downstream start (bp from TSS) of the gene-body window used for the "
                        "dose-response metric (avoids the TSS spike).")
    p.add_argument("--min-gene-length", type=int, default=0,
                   help="Drop genes whose genomic span is shorter than this (bp). Removes the "
                        "gene-length confound when the downstream window overruns short genes.")
    p.add_argument("--sim-norm", action="store_true",
                   help="If a *<sim-tag>* bed matches a timepoint, also emit sim-normalized "
                        "figures: log2(observed/sim) per strand and sim-corrected FC = FC_obs - FC_sim.")
    p.add_argument("--sim-tag", default="sim",
                   help="Filename substring marking simulated (sequence-only) damage beds.")
    return p.parse_args()


def load_tss(gtf: Path, biotype: str, keep_chroms: set[str]) -> pd.DataFrame:
    """One TSS per gene from GTF gene features. Returns chrom/tss/strand/gene_id."""
    rows = []
    want_biotype = None if biotype == "all" else biotype
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom = "chr" + f[0]
            if chrom not in keep_chroms:
                continue
            attrs = f[8]
            if want_biotype is not None and f'gene_biotype "{want_biotype}"' not in attrs:
                continue
            start = int(f[3]) - 1  # GTF is 1-based inclusive -> 0-based
            end = int(f[4])
            strand = f[6]
            tss = start if strand == "+" else end - 1
            gid = ""
            key = 'gene_id "'
            i = attrs.find(key)
            if i != -1:
                j = attrs.find('"', i + len(key))
                gid = attrs[i + len(key):j]
            rows.append((chrom, tss, strand, gid, start, end, end - start))
    tss = pd.DataFrame(rows, columns=["chrom", "tss", "strand", "gene_id",
                                      "gene_start", "gene_end", "length"])
    print(f"  loaded {len(tss):,} TSSs ({biotype}) across {tss.chrom.nunique()} chroms")
    return tss


def make_signal(starts: np.ndarray, ends: np.ndarray, strand: str, size: int,
                mode: str, offset: int) -> np.ndarray:
    """Per-bp signal array (length `size`) for one strand.

    coverage : full read span via diff/cumsum over half-open [start, end).
    5prime   : single count at the read 5' end (start for +, end-1 for -),
               shifted `offset` bp upstream (strand-relative) toward the lesion.
    """
    if mode == "coverage":
        diff = np.zeros(size + 1, dtype=np.int32)
        np.add.at(diff, np.clip(starts, 0, size), 1)
        np.add.at(diff, np.clip(ends, 0, size), -1)
        return np.cumsum(diff[:-1])
    # 5prime
    if strand == "+":
        pos = starts - offset
    else:
        pos = (ends - 1) + offset
    pos = np.clip(pos, 0, size - 1)
    arr = np.zeros(size, dtype=np.int32)
    np.add.at(arr, pos, 1)
    return arr


def _interval_signal(bw_path: Path, cs: dict, chroms: np.ndarray,
                     g0: np.ndarray, g1: np.ndarray) -> np.ndarray:
    """Mean bw signal over each [g0, g1). NaN for off-chrom or non-positive intervals."""
    import bbi
    out = np.full(len(chroms), np.nan, dtype=float)
    s0 = np.maximum(g0.astype(np.int64), 0)
    s1 = g1.astype(np.int64)
    ok = np.array([c in cs for c in chroms]) & (s1 > s0)
    if ok.any():
        ends = np.array([min(e, cs[c]) for e, c in zip(s1[ok], chroms[ok])])
        starts = np.minimum(s0[ok], ends - 1)
        out[ok] = bbi.stackup(str(bw_path), list(chroms[ok]), list(starts), list(ends),
                              bins=1, summary="mean")[:, 0]
    return out


def pol2_scores(tss: pd.DataFrame, bw_path: Path, mode: str, prom_window: int,
                body_start: int = 500, pause_lo: int = 50, pause_hi: int = 300) -> np.ndarray:
    """Per-gene pre-UV PolII score. NaN for genes off the bw chroms.

    promoter      : mean over TSS +/- prom_window/2 (promoter-proximal / paused).
    body          : mean over the elongation region, +body_start..gene-end (transcription
                    oriented), i.e. the pausing-index denominator.
    pausing_index : log2( pause-peak density / body density ); high = more paused.
    """
    import bbi
    cs = bbi.open(str(bw_path)).chromsizes
    chrom = tss.chrom.to_numpy()
    tssv = tss.tss.to_numpy()
    gstart = tss.gene_start.to_numpy()
    gend = tss.gene_end.to_numpy()
    plus = (tss.strand.to_numpy() == "+")

    def body_iv():
        # downstream of TSS by body_start, to the far gene end (transcription oriented)
        g0 = np.where(plus, tssv + body_start, gstart)
        g1 = np.where(plus, gend, tssv - body_start + 1)
        return g0, g1

    def pause_iv():
        g0 = np.where(plus, tssv - pause_lo, tssv - pause_hi)
        g1 = np.where(plus, tssv + pause_hi + 1, tssv + pause_lo + 1)
        return g0, g1

    if mode == "promoter":
        pad = prom_window // 2
        return _interval_signal(bw_path, cs, chrom, tssv - pad, tssv + pad)
    if mode == "body":
        return _interval_signal(bw_path, cs, chrom, *body_iv())
    if mode == "pausing_index":
        eps = 1e-3
        pause = _interval_signal(bw_path, cs, chrom, *pause_iv())
        body = _interval_signal(bw_path, cs, chrom, *body_iv())
        return np.log2((pause + eps) / (body + eps))
    raise ValueError(f"unknown pol2 mode: {mode}")


def extract_window(cov: np.ndarray, g0: int, g1: int) -> np.ndarray:
    """Coverage over genome [g0, g1) as float, NaN-padded where out of bounds."""
    n = g1 - g0
    out = np.full(n, np.nan, dtype=np.float32)
    lo = max(g0, 0)
    hi = min(g1, cov.shape[0])
    if hi > lo:
        out[lo - g0:hi - g0] = cov[lo:hi]
    return out


def stackup_for_timepoint(bed: Path, tss: pd.DataFrame, K: int, L: int, chroms: list[str],
                          swap: bool, mode: str, offset: int, binsize: int):
    """Return binned (NTS, TS) stackups of shape (n_genes, (K+L)//binsize) in CPM.

    Bins each gene's window on the fly so the full per-bp matrix is never allocated
    (keeps memory ~n_genes*nbins instead of n_genes*(K+L))."""
    print(f"  reading {bed.name} ...")
    reads = pd.read_csv(
        bed, sep="\t", header=None, usecols=[0, 1, 2, 5],
        names=["chrom", "start", "end", "strand"],
        dtype={"chrom": "category", "start": np.int64, "end": np.int64, "strand": "category"},
    )
    total = len(reads)
    cpm = 1e6 / total
    print(f"    {total:,} reads, CPM scale = {cpm:.4g}")

    win = K + L
    nbins = win // binsize
    n_genes = len(tss)
    NTS = np.full((n_genes, nbins), np.nan, dtype=np.float32)
    TS = np.full((n_genes, nbins), np.nan, dtype=np.float32)
    idx = np.arange(n_genes)

    def to_bins(row):
        return np.nanmean(row.reshape(nbins, binsize), axis=1)

    for chrom in chroms:
        gmask = tss.chrom.values == chrom
        rmask = reads.chrom.values == chrom
        if not gmask.any() or not rmask.any():
            continue
        rsub = reads[rmask]
        g_idx = idx[gmask]
        g_tss = tss.tss.values[gmask]
        g_strand = tss.strand.values[gmask]

        size = int(max(rsub.end.max(), g_tss.max() + L)) + 2
        cov = {}
        for s in ("+", "-"):
            ss = rsub[rsub.strand.values == s]
            cov[s] = make_signal(ss.start.to_numpy(), ss.end.to_numpy(), s, size, mode, offset)

        for row, t, st in zip(g_idx, g_tss, g_strand):
            sense = cov[st]                      # same strand as gene
            anti = cov["-" if st == "+" else "+"]  # opposite
            nts_cov, ts_cov = (anti, sense) if swap else (sense, anti)
            if st == "+":
                g0, g1, rev = t - K, t + L, False
            else:
                g0, g1, rev = t - L, t + K, True
            nrow = extract_window(nts_cov, g0, g1)
            trow = extract_window(ts_cov, g0, g1)
            if rev:
                nrow = nrow[::-1]
                trow = trow[::-1]
            NTS[row] = to_bins(nrow)
            TS[row] = to_bins(trow)
        del cov
        print(f"    {chrom}: {gmask.sum()} genes, {rmask.sum():,} reads")

    return NTS * cpm, TS * cpm


def bin_rows(arr: np.ndarray, binsize: int) -> np.ndarray:
    """Reduce (n, win) -> (n, win/binsize) by mean (NaN-aware)."""
    n, win = arr.shape
    nbins = win // binsize
    return np.nanmean(arr.reshape(n, nbins, binsize), axis=2)


def mean_band(arr: np.ndarray, sigma: float):
    """NaN-aware mean +/- SEM across genes, optionally smoothed (bins)."""
    center = np.nanmean(arr, axis=0)
    nval = np.sum(~np.isnan(arr), axis=0)
    sd = np.nanstd(arr, axis=0, ddof=1)
    sem = sd / np.sqrt(np.maximum(nval, 1))
    if sigma and sigma > 0:
        center = _smooth_1d(center, sigma=sigma)
        sem = _smooth_1d(sem, sigma=sigma)
    return center, sem


def foldchange_stackup(nts: np.ndarray, ts: np.ndarray, eps: float = 1e-3) -> np.ndarray:
    """Per-gene log2((NTS+eps)/(TS+eps)). >0 = sense-biased."""
    return np.log2((nts + eps) / (ts + eps))


def plot_foldchange(ax, x, fc, sigma, title, color=COLORS[2], label=None,
                    ylabel="log2(NTS / TS)"):
    c, e = mean_band(fc, sigma)
    ax.plot(x, c, color=color, lw=1.6, label=label)
    ax.fill_between(x, c - e, c + e, color=color, alpha=0.2, lw=0)
    ax.axhline(0, color="grey", lw=0.8, alpha=0.6)
    ax.axvline(0, color="grey", ls="--", lw=1.0, alpha=0.7)
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    despine(ax)


def plot_timepoint(ax, x, nts, ts, sigma, title):
    for arr, color, label in ((nts, NTS_COLOR, "NTS (sense)"), (ts, TS_COLOR, "TS (antisense)")):
        c, e = mean_band(arr, sigma)
        ax.plot(x, c, color=color, label=label, lw=1.6)
        ax.fill_between(x, c - e, c + e, color=color, alpha=0.2, lw=0)
    ax.axvline(0, color="grey", ls="--", lw=1.0, alpha=0.7)
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel("Coverage (CPM)")
    ax.set_title(title)
    ax.legend(frameon=False)
    despine(ax)


def main() -> None:
    args = parse_args()
    assert (args.upstream + args.downstream) % args.bin == 0, "K+L must be divisible by --bin"
    K, L, binsize = args.upstream, args.downstream, args.bin
    chroms = args.chroms if args.chroms else MAIN_CHROMS
    keep = set(chroms)

    update_rcparams()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("Parsing TSS from GTF ...")
    tss = load_tss(args.gtf, args.biotype, keep).reset_index(drop=True)

    if args.min_gene_length > 0:
        before = len(tss)
        tss = tss[tss.length >= args.min_gene_length].reset_index(drop=True)
        print(f"  gene-length filter >= {args.min_gene_length:,} bp: kept {len(tss):,}/{before:,}")

    if args.pol2_min is not None or args.pol2_max is not None:
        print(f"PolII filter (mode={args.pol2_mode}) ...")
        tss["pol2"] = pol2_scores(tss, args.pol2_bw, args.pol2_mode, args.pol2_window)
        before = len(tss)
        keep_m = tss.pol2.notna()
        if args.pol2_min is not None:
            keep_m &= tss.pol2 >= args.pol2_min
        if args.pol2_max is not None:
            keep_m &= tss.pol2 <= args.pol2_max
        tss = tss[keep_m].reset_index(drop=True)
        print(f"  kept {len(tss):,}/{before:,} genes "
              f"(PolII in [{args.pol2_min}, {args.pol2_max}])")

    stratify = args.pol2_quantiles is not None
    if stratify:
        nq = args.pol2_quantiles
        print(f"Stratifying genes into {nq} PolII quantiles "
              f"({args.pol2_bw.name}, mode={args.pol2_mode}) ...")
        tss["pol2"] = pol2_scores(tss, args.pol2_bw, args.pol2_mode, args.pol2_window)
        tss = tss[tss.pol2.notna()].reset_index(drop=True)
        tss["Q"] = pd.qcut(tss.pol2.rank(method="first"), nq, labels=range(1, nq + 1)).astype(int)
        for q in range(1, nq + 1):
            sub = tss[tss.Q == q]
            print(f"  Q{q}: {len(sub):,} genes, mean PolII({args.pol2_mode}) {sub.pol2.mean():.3g}")
    elif args.pol2_top_frac is not None:
        print(f"Ranking genes by PolII ({args.pol2_bw.name}, mode={args.pol2_mode}) ...")
        tss["pol2"] = pol2_scores(tss, args.pol2_bw, args.pol2_mode, args.pol2_window)
        valid = tss.pol2.notna()
        thr = tss.loc[valid, "pol2"].quantile(1 - args.pol2_top_frac)
        tss = tss[valid & (tss.pol2 >= thr)].reset_index(drop=True)
        print(f"  kept {len(tss):,} genes (top {args.pol2_top_frac:.0%}, PolII >= {thr:.3g})")

    beds = sorted(args.bed_dir.glob("*.bed"))
    if not beds:
        raise SystemExit(f"No .bed files in {args.bed_dir}")

    nbins = (K + L) // binsize
    x = np.arange(nbins) * binsize - K + binsize / 2  # bin centers, bp from TSS

    def tp_key(name: str) -> str:
        m = re.search(r"(\d+)\s*h", name)
        return f"{m.group(1)}h" if m else name.split(".")[0]

    def is_sim(name: str) -> bool:
        return args.sim_tag.lower() in name.lower()

    # Compute binned (NTS, TS) stackups for every bed once.
    raw: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    label_of: dict[str, str] = {}
    for bed in beds:
        lbl = bed.name.split(".")[0]
        label_of[bed.name] = lbl
        print(f"\n=== {lbl} (signal={args.signal}) ===")
        nts_b, ts_b = stackup_for_timepoint(
            bed, tss, K, L, chroms, args.swap_strands, args.signal, args.lesion_offset, binsize)
        raw[bed.name] = (nts_b, ts_b)

    real_beds = [b for b in beds if not is_sim(b.name)]
    sim_by_key = {tp_key(b.name): raw[b.name] for b in beds if is_sim(b.name)}
    do_sim = args.sim_norm and bool(sim_by_key)
    if args.sim_norm and not sim_by_key:
        print(f"  [sim-norm] no beds matching '*{args.sim_tag}*'; skipping normalization")

    profiles = {}      # label -> (nts, ts, fc)  (observed)
    norm_profiles = {}  # label -> fc_norm        (sim-corrected)
    for bed in real_beds:
        lbl = label_of[bed.name]
        nts, ts = raw[bed.name]
        fc = foldchange_stackup(nts, ts)
        profiles[lbl] = (nts, ts, fc)

        # observed TS/NTS coverage
        fig, ax = plt.subplots(figsize=(5.0, 3.4))
        plot_timepoint(ax, x, nts, ts, args.smooth_sigma, f"{lbl}: TS / NTS around TSS")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(args.out_dir / f"{lbl}_ts_nts_tss.{ext}", dpi=300)
        plt.close(fig)

        # observed fold change
        fig, ax = plt.subplots(figsize=(5.0, 3.4))
        plot_foldchange(ax, x, fc, args.smooth_sigma, f"{lbl}: NTS/TS fold change")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(args.out_dir / f"{lbl}_ts_nts_foldchange.{ext}", dpi=300)
        plt.close(fig)

        nc, ne = mean_band(nts, args.smooth_sigma)
        tc, te = mean_band(ts, args.smooth_sigma)
        fcc, fce = mean_band(fc, args.smooth_sigma)
        np.savez(args.out_dir / f"{lbl}_stackup.npz", NTS=nts, TS=ts, log2FC=fc, x=x)
        cols = {"pos_bp": x, "NTS_mean": nc, "NTS_sem": ne, "TS_mean": tc, "TS_sem": te,
                "log2FC_mean": fcc, "log2FC_sem": fce}

        # sim-normalized variants
        key = tp_key(bed.name)
        if do_sim and key in sim_by_key:
            snts, sts = sim_by_key[key]
            nts_norm = foldchange_stackup(nts, snts)  # log2(obs/sim), NTS
            ts_norm = foldchange_stackup(ts, sts)     # log2(obs/sim), TS
            fc_norm = fc - foldchange_stackup(snts, sts)  # FC_obs - FC_sim
            norm_profiles[lbl] = fc_norm

            fig, ax = plt.subplots(figsize=(5.0, 3.4))
            for arr, color, lab in ((nts_norm, NTS_COLOR, "NTS (sense)"),
                                    (ts_norm, TS_COLOR, "TS (antisense)")):
                c, e = mean_band(arr, args.smooth_sigma)
                ax.plot(x, c, color=color, lw=1.6, label=lab)
                ax.fill_between(x, c - e, c + e, color=color, alpha=0.2, lw=0)
            ax.axhline(0, color="grey", lw=0.8, alpha=0.6)
            ax.axvline(0, color="grey", ls="--", lw=1.0, alpha=0.7)
            ax.set_xlabel("Distance from TSS (bp)")
            ax.set_ylabel("log2(observed / sim)")
            ax.set_title(f"{lbl}: damage enrichment over expected")
            ax.legend(frameon=False)
            despine(ax)
            fig.tight_layout()
            for ext in ("png", "pdf"):
                fig.savefig(args.out_dir / f"{lbl}_ts_nts_simnorm.{ext}", dpi=300)
            plt.close(fig)

            fig, ax = plt.subplots(figsize=(5.0, 3.4))
            plot_foldchange(ax, x, fc_norm, args.smooth_sigma,
                            f"{lbl}: sim-corrected NTS/TS fold change",
                            ylabel=r"$\Delta$log2(NTS/TS) vs sim")
            fig.tight_layout()
            for ext in ("png", "pdf"):
                fig.savefig(args.out_dir / f"{lbl}_ts_nts_foldchange_simnorm.{ext}", dpi=300)
            plt.close(fig)

            nnc, nne = mean_band(nts_norm, args.smooth_sigma)
            tnc, tne = mean_band(ts_norm, args.smooth_sigma)
            fnc, fne = mean_band(fc_norm, args.smooth_sigma)
            cols.update({"NTS_simnorm_mean": nnc, "NTS_simnorm_sem": nne,
                         "TS_simnorm_mean": tnc, "TS_simnorm_sem": tne,
                         "log2FC_simnorm_mean": fnc, "log2FC_simnorm_sem": fne})

        pd.DataFrame(cols).to_csv(args.out_dir / f"{lbl}_mean_profile.tsv", sep="\t", index=False)
        print(f"  wrote {lbl} figures + matrices to {args.out_dir}")

    tp_colors = COLORS[:max(len(profiles), 1)]

    # combined: observed TS/NTS coverage, NTS solid / TS dashed
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for (lbl, (nts, ts, _fc)), col in zip(profiles.items(), tp_colors):
        nc, _ = mean_band(nts, args.smooth_sigma)
        tc, _ = mean_band(ts, args.smooth_sigma)
        ax.plot(x, nc, color=col, lw=1.6, label=f"{lbl} NTS")
        ax.plot(x, tc, color=col, lw=1.3, ls="--", label=f"{lbl} TS")
    ax.axvline(0, color="grey", ls="--", lw=1.0, alpha=0.7)
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel("Coverage (CPM)")
    ax.set_title("TS / NTS around TSS: timepoints")
    ax.legend(frameon=False, fontsize=7)
    despine(ax)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(args.out_dir / f"combined_ts_nts_tss.{ext}", dpi=300)
    plt.close(fig)

    # combined observed fold change
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for (lbl, (_nts, _ts, fc)), col in zip(profiles.items(), tp_colors):
        c, e = mean_band(fc, args.smooth_sigma)
        ax.plot(x, c, color=col, lw=1.6, label=lbl)
        ax.fill_between(x, c - e, c + e, color=col, alpha=0.15, lw=0)
    ax.axhline(0, color="grey", lw=0.8, alpha=0.6)
    ax.axvline(0, color="grey", ls="--", lw=1.0, alpha=0.7)
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel("log2(NTS / TS)")
    ax.set_title("NTS/TS fold change around TSS: timepoints")
    ax.legend(frameon=False, fontsize=7)
    despine(ax)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(args.out_dir / f"combined_ts_nts_foldchange.{ext}", dpi=300)
    plt.close(fig)

    # combined sim-corrected fold change
    if norm_profiles:
        fig, ax = plt.subplots(figsize=(5.4, 3.6))
        for (lbl, fc_norm), col in zip(norm_profiles.items(), tp_colors):
            c, e = mean_band(fc_norm, args.smooth_sigma)
            ax.plot(x, c, color=col, lw=1.6, label=lbl)
            ax.fill_between(x, c - e, c + e, color=col, alpha=0.15, lw=0)
        ax.axhline(0, color="grey", lw=0.8, alpha=0.6)
        ax.axvline(0, color="grey", ls="--", lw=1.0, alpha=0.7)
        ax.set_xlabel("Distance from TSS (bp)")
        ax.set_ylabel(r"$\Delta$log2(NTS/TS) vs sim")
        ax.set_title("Sim-corrected NTS/TS fold change: timepoints")
        ax.legend(frameon=False, fontsize=7)
        despine(ax)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(args.out_dir / f"combined_ts_nts_foldchange_simnorm.{ext}", dpi=300)
        plt.close(fig)

    # ---- PolII-quantile dose-response ----
    if stratify:
        nq = args.pol2_quantiles
        Qs = list(range(1, nq + 1))
        qmask = {q: (tss.Q.values == q) for q in Qs}
        meanpol2 = {q: float(tss.pol2.values[qmask[q]].mean()) for q in Qs}
        region = (x >= args.body_start) & (x <= L)  # downstream body, skips TSS spike
        qcmap = plt.cm.viridis(np.linspace(0, 0.9, nq))
        ylab = r"$\Delta$log2(NTS/TS) vs sim" if do_sim else "log2(NTS / TS)"
        records = []
        pergene_body = {}   # label -> per-gene body sim-corrected log2(NTS/TS)
        fc_stacks = {}      # label -> sim-corrected log2(NTS/TS) stackup (n_genes x nbins)

        for bed in real_beds:
            lbl = label_of[bed.name]
            nts, ts = raw[bed.name]
            key = tp_key(bed.name)
            if do_sim and key in sim_by_key:
                snts, sts = sim_by_key[key]
                fc = foldchange_stackup(nts, ts) - foldchange_stackup(snts, sts)
            else:
                fc = foldchange_stackup(nts, ts)
            fc_stacks[lbl] = fc

            # per-quantile FC profile (low->high PolII)
            fig, ax = plt.subplots(figsize=(5.4, 3.6))
            for q, col in zip(Qs, qcmap):
                c, _ = mean_band(fc[qmask[q]], args.smooth_sigma)
                ax.plot(x, c, color=col, lw=1.4, label=f"Q{q}")
            ax.axhline(0, color="grey", lw=0.8, alpha=0.6)
            ax.axvline(0, color="grey", ls="--", lw=1.0, alpha=0.7)
            ax.set_xlabel("Distance from TSS (bp)")
            ax.set_ylabel(ylab)
            ax.set_title(f"{lbl}: NTS/TS by PolII quantile")
            ax.legend(frameon=False, fontsize=7, ncol=2, title="PolII (low->high)")
            despine(ax)
            fig.tight_layout()
            for ext in ("png", "pdf"):
                fig.savefig(args.out_dir / f"{lbl}_simnorm_fc_by_pol2q.{ext}", dpi=300)
            plt.close(fig)

            # dose-response metric: per-gene body mean, then mean +/- SEM per quantile
            body = np.nanmean(fc[:, region], axis=1)
            pergene_body[lbl] = body
            for q in Qs:
                v = body[qmask[q]]
                v = v[~np.isnan(v)]
                records.append({"timepoint": lbl, "Q": q, "mean_pol2": meanpol2[q],
                                "n": v.size, "body_FC": float(v.mean()),
                                "body_FC_sem": float(v.std(ddof=1) / np.sqrt(max(v.size, 1)))})

        dose = pd.DataFrame.from_records(records)
        dose.to_csv(args.out_dir / "dose_response_pol2q.tsv", sep="\t", index=False)

        # per-gene body table (one column per timepoint)
        pg = tss[["gene_id", "chrom", "strand", "length", "pol2", "Q"]].copy()
        for lbl, arr in pergene_body.items():
            pg[f"bodyFC_{lbl}"] = arr
        pg.to_csv(args.out_dir / "per_gene_body_fc.tsv", sep="\t", index=False)

        # ---- TC-NER test: paired 3h vs 0h within transcription strata ----
        tcner_test(args, x, L, Qs, qmask, meanpol2, dose, pergene_body, fc_stacks,
                   tss, region, ylab)
        print(f"  wrote dose-response + TC-NER test ({nq} PolII quantiles) to {args.out_dir}")

    # ---- repair from timepoints + strand asymmetry (0h-normalized) ----
    repair_asymmetry(args, x, L, raw, real_beds, label_of, tss)

    print(f"\nDone. Outputs in {args.out_dir}")


def _hour(label: str) -> float:
    m = re.search(r"(\d+)\s*h", label)
    return float(m.group(1)) if m else float("inf")


def repair_asymmetry(args, x, L, raw, real_beds, label_of, tss, eps=1e-3):
    """Repair computed from the two real timepoints, 0h as the initial-damage baseline.

    Per strand, per gene/bin: relative change r = log2(late/early) (CPM-relative).
    Strand asymmetry of repair A = r_NTS - r_TS = Dlog2(NTS/TS)_{late-early}; the
    per-library CPM factor cancels in A, so A>0 = TS lost more signal = TC-NER."""
    beds = sorted(real_beds, key=lambda b: _hour(label_of[b.name]))
    if len(beds) < 2:
        print("  [repair] need >=2 timepoints; skipping")
        return
    b0, b1 = beds[0], beds[-1]
    l0, l1 = label_of[b0.name], label_of[b1.name]
    nts0, ts0 = raw[b0.name]
    nts3, ts3 = raw[b1.name]

    if "pol2" in tss.columns:
        g = tss.pol2.values >= np.nanmedian(tss.pol2.values)
        gtag = "active genes (PolII>=median)"
    else:
        g = np.ones(len(tss), dtype=bool)
        gtag = "all genes"

    rN = np.log2((nts3 + eps) / (nts0 + eps))   # per-strand relative change
    rT = np.log2((ts3 + eps) / (ts0 + eps))
    A = rN - rT                                  # CPM-invariant repair asymmetry

    # 1) per-strand relative change (CPM-relative): TS drops more than NTS
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for arr, color, name in ((rN, NTS_COLOR, "NTS (sense)"), (rT, TS_COLOR, "TS (antisense)")):
        c, e = mean_band(arr[g], args.smooth_sigma)
        ax.plot(x, c, color=color, lw=1.8, label=name)
        ax.fill_between(x, c - e, c + e, color=color, alpha=0.18, lw=0)
    ax.axhline(0, color="grey", lw=0.8, alpha=0.6)
    ax.axvline(0, color="grey", ls="--", lw=1.0, alpha=0.7)
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel(f"log2({l1}/{l0}) per strand  (CPM-relative)")
    ax.set_title(f"Relative damage change, {gtag}: TS lost more")
    ax.legend(frameon=False)
    despine(ax)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(args.out_dir / f"repair_perstrand_relchange.{ext}", dpi=300)
    plt.close(fig)

    # 1b) absolute per-strand repair = 0h - 3h (CPM removed; positive = damage lost)
    dN = nts0 - nts3
    dT = ts0 - ts3
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for arr, color, name in ((dN, NTS_COLOR, "NTS (sense)"), (dT, TS_COLOR, "TS (antisense)")):
        c, e = mean_band(arr[g], args.smooth_sigma)
        ax.plot(x, c, color=color, lw=1.8, label=name)
        ax.fill_between(x, c - e, c + e, color=color, alpha=0.18, lw=0)
    ax.axhline(0, color="grey", lw=0.8, alpha=0.6)
    ax.axvline(0, color="grey", ls="--", lw=1.0, alpha=0.7)
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel(f"{l0} - {l1} coverage (CPM removed)")
    ax.set_title(f"Absolute repair (0h-3h), {gtag}: TS removed more")
    ax.legend(frameon=False)
    despine(ax)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(args.out_dir / f"repair_perstrand_absdiff.{ext}", dpi=300)
    plt.close(fig)

    # 2) repair strand asymmetry A (CPM-invariant) + paired body test
    from scipy.stats import wilcoxon
    region = (x >= args.body_start) & (x <= L)
    body = np.nanmean(A[g][:, region], axis=1)
    body = body[~np.isnan(body)]
    p = wilcoxon(body)[1] if body.size > 10 else np.nan
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    c, e = mean_band(A[g], args.smooth_sigma)
    ax.plot(x, c, color=COLORS[2], lw=1.8)
    ax.fill_between(x, c - e, c + e, color=COLORS[2], alpha=0.2, lw=0)
    ax.axhline(0, color="grey", lw=0.8, alpha=0.6)
    ax.axvline(0, color="grey", ls="--", lw=1.0, alpha=0.7)
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel(r"repair asymmetry $\Delta$log2(NTS/TS)")
    ax.set_title(f"TS repaired faster than NTS ({l1} vs {l0}; body p={p:.1e}, n={body.size})")
    despine(ax)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(args.out_dir / f"repair_strand_asymmetry.{ext}", dpi=300)
    plt.close(fig)
    print(f"  repair asymmetry ({l1} vs {l0}, {gtag}): body median={np.median(body):+.4f} p={p:.1e}")


def _stars(p: float) -> str:
    return ("***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns")


def tcner_test(args, x, L, Qs, qmask, meanpol2, dose, pergene_body, fc_stacks,
               tss, region, ylab):
    """Show TC-NER active by the late timepoint: TS preferentially repaired vs NTS,
    increasing with transcription. Paired (late vs early)."""
    from scipy.stats import wilcoxon

    tps = sorted(pergene_body, key=_hour)
    if len(tps) < 2:
        print("  [tcner] need >=2 timepoints; skipping test")
        return
    t0, t1 = tps[0], tps[-1]   # earliest vs latest
    b0, b1 = pergene_body[t0], pergene_body[t1]

    # per-quantile paired Wilcoxon (late vs early)
    qrows = []
    for q in Qs:
        m = qmask[q]
        d = (b1 - b0)[m]
        d = d[~np.isnan(d)]
        stat, p = wilcoxon(d) if d.size > 10 else (np.nan, np.nan)
        qrows.append({"Q": q, "mean_pol2": meanpol2[q], "n": int(d.size),
                      "median_delta": float(np.median(d)) if d.size else np.nan,
                      "wilcoxon_p": float(p)})
    qstats = pd.DataFrame(qrows)
    qstats.to_csv(args.out_dir / "tcner_paired_by_quantile.tsv", sep="\t", index=False)

    # annotate dose-response with per-quantile significance (late vs early)
    fig, ax = plt.subplots(figsize=(5.2, 3.7))
    for lbl, col in zip(dose.timepoint.unique(), COLORS):
        d = dose[dose.timepoint == lbl]
        ax.errorbar(d.Q, d.body_FC, yerr=d.body_FC_sem, color=col, marker="o",
                    ms=4, lw=1.5, capsize=2, label=lbl)
    ax.axhline(0, color="grey", lw=0.8, alpha=0.6)
    ytop = dose.body_FC.max() + dose.body_FC_sem.max()
    for _, r in qstats.iterrows():
        ax.text(r.Q, ytop * 1.04, _stars(r.wilcoxon_p), ha="center", va="bottom", fontsize=8)
    ax.set_ylim(top=ytop * 1.15)
    ax.set_xlabel("PolII quantile (low -> high)")
    ax.set_ylabel(f"Gene-body {ylab}\n(+{args.body_start}..+{L} bp)")
    ax.set_title(f"TC-NER dose-response ({t1} vs {t0}, paired)")
    ax.set_xticks(Qs)
    ax.legend(frameon=False)
    despine(ax)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(args.out_dir / f"dose_response_pol2q.{ext}", dpi=300)
    plt.close(fig)

    # transcription strata. Body-ChIP occupancy turns over at the very top quantile
    # (occupancy != elongation), so "active" = expressed (>= median), not the top tail.
    med, p33 = np.nanpercentile(tss.pol2.values, [50, 33])
    active = (tss.pol2.values >= med)   # expressed / transcribed (TC-NER substrates)
    low = (tss.pol2.values <= p33)      # low transcription (weak TCR expected)
    groups = [("active (PolII>=median)", active), ("low (bottom tercile)", low)]
    pvals = {}
    for gname, gmask in groups:
        d0, d1 = pergene_body[t0][gmask], pergene_body[t1][gmask]
        ok = ~(np.isnan(d0) | np.isnan(d1))
        _, p = wilcoxon((d1 - d0)[ok]) if ok.sum() > 10 else (np.nan, np.nan)
        pvals[gname] = (p, int(ok.sum()))

    # ---- HEADLINE: aggregate sim-corrected log2(NTS/TS) metagene, active genes,
    #      early vs late. Strand bias rises over the gene body by the late point. ----
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for lbl, col in ((t0, COLORS[0]), (t1, COLORS[1])):
        c, e = mean_band(fc_stacks[lbl][active], args.smooth_sigma)
        ax.plot(x, c, color=col, lw=1.8, label=lbl)
        ax.fill_between(x, c - e, c + e, color=col, alpha=0.18, lw=0)
    ax.axhline(0, color="grey", lw=0.8, alpha=0.6)
    ax.axvline(0, color="grey", ls="--", lw=1.0, alpha=0.7)
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel(ylab)
    pA = pvals["active (PolII>=median)"]
    ax.set_title(f"Active genes: NTS/TS bias rises by {t1} (paired p={pA[0]:.1e}, n={pA[1]})")
    ax.legend(frameon=False, title="timepoint")
    despine(ax)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(args.out_dir / f"tcner_active_fc_metagene.{ext}", dpi=300)
    plt.close(fig)

    # ---- supplementary boxplot: per-gene body bias, active vs low, early vs late ----
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    positions, data, colors, xticklab = [], [], [], []
    pos = 0
    for gname, gmask in groups:
        for lbl, col in ((t0, COLORS[0]), (t1, COLORS[1])):
            v = pergene_body[lbl][gmask]; v = v[~np.isnan(v)]
            positions.append(pos); data.append(v); colors.append(col)
            xticklab.append(lbl); pos += 1
        pos += 1
    bp = ax.boxplot(data, positions=positions, widths=0.7, showfliers=False,
                    patch_artist=True, medianprops=dict(color="black"))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.6)
    ax.axhline(0, color="grey", lw=0.8, alpha=0.6)
    ax.set_xticks(positions); ax.set_xticklabels(xticklab, fontsize=8)
    ax.set_ylabel(f"Gene-body {ylab}")
    ytxt = max(np.nanpercentile(np.concatenate(data), 97), 0.0)
    for (gname, _), gx in zip(groups, [1.5, 4.5]):
        p, n = pvals[gname]
        ax.text(gx, ytxt * 1.18, f"{gname.split()[0]}\n{_stars(p)} (n={n})",
                ha="center", va="bottom", fontsize=8)
    ax.set_ylim(top=ytxt * 1.5 if ytxt > 0 else None)
    ax.set_title(f"Per-gene body bias {t0}->{t1}")
    despine(ax)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(args.out_dir / f"tcner_active_vs_low_box.{ext}", dpi=300)
    plt.close(fig)

    print(f"  TC-NER paired ({t1} vs {t0}): "
          + ", ".join(f"{g.split()[0]} p={pvals[g][0]:.1e}" for g in pvals))


if __name__ == "__main__":
    main()

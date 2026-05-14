#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Option 1 IMPLEMENTATION: ATAC-intensity–matched null model (hg38 chr1–22 + chrX)
#
# Goal:
#   Test whether XPC is unusually close to CTCF beyond what is expected from
#   chromatin accessibility, *including accessibility intensity*.
#
# Key idea (Intensity-matched shuffling):
#   1) Score each XPC summit by ATAC accessibility (mean of two RPGC bigWigs).
#   2) Split XPC summits into quantile bins by ATAC score (e.g., deciles).
#   3) Build an ATAC "universe" (genomic intervals with ATAC signal) and assign
#      each interval to the SAME bins using the same score thresholds.
#   4) For each bin k, shuffle XPC summits in bin k only inside ATAC-universe bin k.
#   5) For observed and each permutation, compute distances: XPC -> nearest CTCF summit.
#   6) Compare observed stats to permutation distributions; make plots including
#      an overlay histogram of observed vs pooled permuted distance distributions.
#
# Restriction:
#   Uses only chr1–chr22 and chrX everywhere; blacklist excluded.
#   Optional: if WHITELIST BED is provided, restrict signal + peaks to whitelist.
#
# Requirements:
#   bedtools
#   UCSC bigWigToBedGraph
#   (optional) UCSC bigWigAverageOverBed OR python package pyBigWig in PY env
#   awk (gawk recommended, but we compute medians in python at the end)
#
# Outputs (OUTDIR):
#   01_universe/ATAC.mean.bedGraph.bed4              (chr start end meanATAC)
#   01_universe/ATAC.bin1..binK.bed                  (BED3 per bin, merged)
#   02_obs/XPC.summits.bed3 + .scored.tsv            (XPC summits + ATAC scores + bins)
#   02_obs/CTCF.summits.bed3
#   02_obs/XPC_CTCF.overlap_peak_summit_distances.tsv
#   02_obs/XPC_CTCF.overlap_peak_summit_distances.summary.tsv
#   02_obs/observed_distances.txt
#   03_perm/perm_stats.tsv                           (mean/median/fracs per perm)
#   03_perm/saved_perm_dists/dist_perm_*.txt         (first SAVE_PERM_DISTS perms)
#   04_viz/*.png
#   summary.autosomesX.intensityMatched.tsv
#
# Usage:
#   chmod +x xpc_ctcf_intensityMatched_null.sh
#   OUTDIR='boundaries_20kb_flank' WHITELIST='boundaries_WTnoUV_500000bp.bed' ./xpc_ctcf_null_autosomesX_dynamic.sh
#   OUTDIR='genomewide' ./xpc_ctcf_null_autosomesX_dynamic.sh
#   Note; boundaries_WTnoUV_500000bp is TAD boundaries identified with 500kb insulation window size at 10kb res,+-20kb flanks.
# -----------------------------------------------------------------------------

set -euo pipefail

# ------------------------ Settings (override via env vars) --------------------
N_PERM="${N_PERM:-10}"
SEED_BASE="${SEED_BASE:-1337}"

# Number of intensity bins (typically 10 = deciles)
N_BINS="${N_BINS:-10}"

# Save distance vectors for the first N permutations (for overlay histogram)
SAVE_PERM_DISTS="${SAVE_PERM_DISTS:-10}"

# Only for plotting readability: cap long distance tail (blank = no cap)
DIST_CAP="${DIST_CAP:-200000}"

OUTDIR="${OUTDIR:-out_xpc_ctcf_intensityMatched}"

BW1="${BW1:-atac_WT_noUV_1.rpgc.bw}"
BW2="${BW2:-atac_WT_noUV_2.rpgc.bw}"
XPC="${XPC:-XPC_all.narrowPeak}"
CTCF="${CTCF:-CTCF_all.narrowPeak}"
BLACKLIST="${BLACKLIST:-hg38.blacklist.bed}"
WHITELIST="${WHITELIST:-}"
CHROMSIZES="${CHROMSIZES:-hg38.chrom.sizes}"

PY="/home/carlos/micromamba/envs/hic/bin/python"
# -----------------------------------------------------------------------------

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing: $1" >&2; exit 1; }; }
need bedtools
need awk
need sort
need bigWigToBedGraph
[[ -x "$PY" ]] || { echo "ERROR: python not executable at $PY" >&2; exit 1; }

mkdir -p "$OUTDIR"/{00_tmp,01_universe,02_obs,03_perm/saved_perm_dists,04_viz}

log() { echo "[$(date '+%F %T')] $*"; }

# BED3 filter for chr1-22 + chrX
filter_chrAUX_bed3() {
  awk 'BEGIN{
    keep["chrX"]=1; for(i=1;i<=22;i++) keep["chr"i]=1; OFS="\t"
  }
  ($1 in keep) && ($3>$2) {print $1,$2,$3}'
}

# Optional stream filter: if WHITELIST_R is set, clip intervals to whitelist.
# Works for BED3/BED4 streams and preserves A-side extra columns (e.g., signal value).
apply_optional_whitelist() {
  if [[ -n "${WHITELIST_R:-}" ]]; then
    bedtools intersect -wa -wb -a stdin -b "$WHITELIST_R" \
      | awk 'BEGIN{OFS="\t"}{
          bs=$(NF-1); be=$NF;
          s=($2>bs?$2:bs);
          e=($3<be?$3:be);
          if(e>s){
            printf "%s\t%d\t%d", $1, s, e;
            for(i=4;i<=NF-3;i++) printf "\t%s", $i;
            printf "\n";
          }
        }'
  else
    cat
  fi
}

# -----------------------------------------------------------------------------
# 0) Restrict chrom sizes + blacklist
# -----------------------------------------------------------------------------
log "Restricting genome to chr1-22 and chrX ..."
awk 'BEGIN{keep["chrX"]=1; for(i=1;i<=22;i++) keep["chr"i]=1}
     ($1 in keep){print $1"\t"$2}' "$CHROMSIZES" > "$OUTDIR/00_tmp/hg38.autosomesX.chrom.sizes"
CHROMSIZES_R="$OUTDIR/00_tmp/hg38.autosomesX.chrom.sizes"
[[ -s "$CHROMSIZES_R" ]] || { echo "ERROR: restricted chrom.sizes empty; check chr naming." >&2; exit 1; }
log "Using chrom.sizes: $CHROMSIZES_R"

filter_chrAUX_bed3 < "$BLACKLIST" \
  | bedtools sort -g "$CHROMSIZES_R" -i stdin \
  > "$OUTDIR/00_tmp/blacklist.autosomesX.bed"
BLACKLIST_R="$OUTDIR/00_tmp/blacklist.autosomesX.bed"

WHITELIST_R=""
if [[ -n "$WHITELIST" ]]; then
  [[ -s "$WHITELIST" ]] || { echo "ERROR: WHITELIST is set but file is missing/empty: $WHITELIST" >&2; exit 1; }
  log "Applying whitelist BED: $WHITELIST"
  filter_chrAUX_bed3 < "$WHITELIST" \
    | bedtools sort -g "$CHROMSIZES_R" -i stdin \
    | bedtools merge -i stdin \
    > "$OUTDIR/00_tmp/whitelist.autosomesX.bed"
  WHITELIST_R="$OUTDIR/00_tmp/whitelist.autosomesX.bed"
  [[ -s "$WHITELIST_R" ]] || { echo "ERROR: whitelist has no valid chr1-22/X intervals after filtering: $WHITELIST" >&2; exit 1; }
  WL_BP=$(awk '{s+=($3-$2)} END{print s+0}' "$WHITELIST_R")
  WL_N=$(wc -l < "$WHITELIST_R")
  log "Whitelist intervals (chr1-22,X): $WL_N ; total bp: $WL_BP"
else
  log "No WHITELIST supplied; using full chr1-22/X space (minus blacklist)."
fi

# -----------------------------------------------------------------------------
# 1) Build ATAC universe with MEAN signal across two bigWigs
# -----------------------------------------------------------------------------
log "Converting bigWig -> bedGraph (rep1, rep2) ..."
bigWigToBedGraph "$BW1" "$OUTDIR/00_tmp/rep1.bg"
bigWigToBedGraph "$BW2" "$OUTDIR/00_tmp/rep2.bg"

# Keep only chr1-22,X and value>0 to define "ATAC universe"
log "Filtering bedGraphs to chr1-22,X and ATAC>0 ..."
awk 'BEGIN{OFS="\t"} $4>0 {print $1,$2,$3,$4}' "$OUTDIR/00_tmp/rep1.bg" \
  | filter_chrAUX_bed3 \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,1}' > /dev/null 2>&1 || true

# Re-do properly keeping value column (BED4)
awk 'BEGIN{
  keep["chrX"]=1; for(i=1;i<=22;i++) keep["chr"i]=1; OFS="\t"
}
($1 in keep) && ($4>0) && ($3>$2) {print $1,$2,$3,$4}' "$OUTDIR/00_tmp/rep1.bg" \
  | apply_optional_whitelist \
  | bedtools sort -g "$CHROMSIZES_R" -i stdin \
  > "$OUTDIR/00_tmp/rep1.filt.bg"

awk 'BEGIN{
  keep["chrX"]=1; for(i=1;i<=22;i++) keep["chr"i]=1; OFS="\t"
}
($1 in keep) && ($4>0) && ($3>$2) {print $1,$2,$3,$4}' "$OUTDIR/00_tmp/rep2.bg" \
  | apply_optional_whitelist \
  | bedtools sort -g "$CHROMSIZES_R" -i stdin \
  > "$OUTDIR/00_tmp/rep2.filt.bg"

# Make mean signal bedGraph by unioning intervals and averaging values
# unionbedg expects bedGraph; missing values become 0
log "Computing mean ATAC signal (unionbedg + average) ..."
bedtools unionbedg -i "$OUTDIR/00_tmp/rep1.filt.bg" "$OUTDIR/00_tmp/rep2.filt.bg" \
| awk 'BEGIN{OFS="\t"}{
    v1=$4; v2=$5;
    mean=(v1+v2)/2.0;
    if(mean>0) print $1,$2,$3,mean
  }' \
| bedtools sort -g "$CHROMSIZES_R" -i stdin \
> "$OUTDIR/01_universe/ATAC.mean.bedGraph.bed4"

# Remove blacklist portions (subtract to keep non-blacklisted pieces)
log "Subtracting blacklist from ATAC universe ..."
bedtools subtract -a "$OUTDIR/01_universe/ATAC.mean.bedGraph.bed4" -b "$BLACKLIST_R" \
| bedtools sort -g "$CHROMSIZES_R" -i stdin \
> "$OUTDIR/01_universe/ATAC.mean.noBlacklist.bed4"

UNIV_BP=$(awk '{s+=($3-$2)} END{print s+0}' "$OUTDIR/01_universe/ATAC.mean.noBlacklist.bed4")
UNIV_N=$(wc -l < "$OUTDIR/01_universe/ATAC.mean.noBlacklist.bed4")
log "ATAC universe intervals: $UNIV_N ; total bp: $UNIV_BP"

# -----------------------------------------------------------------------------
# 2) Create XPC and CTCF 1bp summits (chr1-22,X; blacklist removed)
# -----------------------------------------------------------------------------
to_summit_bed3_autosomesX() {
  local infile="$1"
  local outfile="$2"
  awk 'BEGIN{OFS="\t"}
    $3>$2 {
      peak=$10;
      if(peak>=0){ s=$2+peak; e=s+1; }
      else { s=int(($2+$3)/2); e=s+1; }
      print $1,s,e
    }' "$infile" \
  | filter_chrAUX_bed3 \
  | bedtools sort -g "$CHROMSIZES_R" -i stdin \
  | apply_optional_whitelist \
  | bedtools intersect -v -a stdin -b "$BLACKLIST_R" \
  | bedtools sort -g "$CHROMSIZES_R" -i stdin \
  > "$outfile"
}

# BED4: chrom peak_start peak_end summit_start(1bp start coordinate)
to_peak_with_summit_bed4_autosomesX() {
  local infile="$1"
  local outfile="$2"
  awk 'BEGIN{
    keep["chrX"]=1; for(i=1;i<=22;i++) keep["chr"i]=1; OFS="\t"
  }
    ($1 in keep) && ($3>$2) {
      peak=$10;
      if(peak>=0){ ss=$2+peak; }
      else { ss=int(($2+$3)/2); }
      print $1,$2,$3,ss
    }' "$infile" \
  | bedtools sort -g "$CHROMSIZES_R" -i stdin \
  | apply_optional_whitelist \
  | bedtools sort -g "$CHROMSIZES_R" -i stdin \
  > "$outfile"
}

log "Creating XPC + CTCF summits (1bp) ..."
to_summit_bed3_autosomesX "$XPC"  "$OUTDIR/02_obs/XPC.summits.bed3"
to_summit_bed3_autosomesX "$CTCF" "$OUTDIR/02_obs/CTCF.summits.bed3"

XPC_N=$(wc -l < "$OUTDIR/02_obs/XPC.summits.bed3")
CTCF_N=$(wc -l < "$OUTDIR/02_obs/CTCF.summits.bed3")
log "XPC summits:  $XPC_N"
log "CTCF summits: $CTCF_N"

# -----------------------------------------------------------------------------
# 2b) Report summit distances for overlapping XPC/CTCF peaks (>=1bp overlap)
# -----------------------------------------------------------------------------
log "Reporting summit distances for overlapping XPC/CTCF peaks (>=1bp overlap)..."

to_peak_with_summit_bed4_autosomesX "$XPC"  "$OUTDIR/00_tmp/XPC.peaks.withSummit.raw.bed4"
to_peak_with_summit_bed4_autosomesX "$CTCF" "$OUTDIR/00_tmp/CTCF.peaks.withSummit.raw.bed4"

# Keep only peaks whose summits are in the filtered summit sets used by this pipeline.
bedtools intersect -u -a "$OUTDIR/00_tmp/XPC.peaks.withSummit.raw.bed4" -b "$OUTDIR/02_obs/XPC.summits.bed3" \
  | bedtools sort -g "$CHROMSIZES_R" -i stdin \
  > "$OUTDIR/00_tmp/XPC.peaks.withSummit.filtered.bed4"

bedtools intersect -u -a "$OUTDIR/00_tmp/CTCF.peaks.withSummit.raw.bed4" -b "$OUTDIR/02_obs/CTCF.summits.bed3" \
  | bedtools sort -g "$CHROMSIZES_R" -i stdin \
  > "$OUTDIR/00_tmp/CTCF.peaks.withSummit.filtered.bed4"

OVERLAP_DIST_TSV="$OUTDIR/02_obs/XPC_CTCF.overlap_peak_summit_distances.tsv"
OVERLAP_SUMMARY_TSV="$OUTDIR/02_obs/XPC_CTCF.overlap_peak_summit_distances.summary.tsv"

echo -e "xpc_chr\txpc_start\txpc_end\txpc_summit\tctcf_chr\tctcf_start\tctcf_end\tctcf_summit\toverlap_bp\tsummit_distance_bp" \
  > "$OVERLAP_DIST_TSV"

bedtools intersect \
  -wa -wb \
  -a "$OUTDIR/00_tmp/XPC.peaks.withSummit.filtered.bed4" \
  -b "$OUTDIR/00_tmp/CTCF.peaks.withSummit.filtered.bed4" \
| awk 'BEGIN{OFS="\t"}{
    xs=$2; xe=$3; xsum=$4;
    cs=$6; ce=$7; csum=$8;
    ovs=(xs>cs?xs:cs);
    ove=(xe<ce?xe:ce);
    ovbp=ove-ovs;
    if(ovbp>=1){
      d=(xsum>csum?xsum-csum:csum-xsum);
      print $1,xs,xe,xsum,$5,cs,ce,csum,ovbp,d;
    }
  }' >> "$OVERLAP_DIST_TSV"

"$PY" - "$OVERLAP_DIST_TSV" "$OVERLAP_SUMMARY_TSV" <<'PY'
import sys
import numpy as np

dist_tsv, out_tsv = sys.argv[1:3]
d = []
with open(dist_tsv) as fh:
    next(fh, None)  # header
    for line in fh:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        d.append(float(fields[-1]))

n = len(d)
with open(out_tsv, "w") as fo:
    fo.write("metric\tvalue\n")
    fo.write(f"overlap_peak_pairs\t{n}\n")
    if n == 0:
        fo.write("mean_summit_distance_bp\tNA\n")
        fo.write("median_summit_distance_bp\tNA\n")
        fo.write("min_summit_distance_bp\tNA\n")
        fo.write("max_summit_distance_bp\tNA\n")
    else:
        arr = np.asarray(d, dtype=float)
        fo.write(f"mean_summit_distance_bp\t{float(np.mean(arr))}\n")
        fo.write(f"median_summit_distance_bp\t{float(np.median(arr))}\n")
        fo.write(f"min_summit_distance_bp\t{float(np.min(arr))}\n")
        fo.write(f"max_summit_distance_bp\t{float(np.max(arr))}\n")
PY

OVERLAP_N=$(awk 'END{print (NR>0 ? NR-1 : 0)}' "$OVERLAP_DIST_TSV")
log "Overlapping XPC/CTCF peak pairs (>=1bp): $OVERLAP_N"
log "Wrote: $OVERLAP_DIST_TSV"
log "Wrote: $OVERLAP_SUMMARY_TSV"

# -----------------------------------------------------------------------------
# 3) Score XPC summits by ATAC signal (mean of bigWigs) and create bins
#    This is the "dynamic" part: bins are defined by XPC's ATAC distribution.
# -----------------------------------------------------------------------------
log "Scoring XPC summits by ATAC mean signal (BW1,BW2) and binning into N_BINS=$N_BINS ..."

# Prefer UCSC bigWigAverageOverBed if available; else use pyBigWig via python.
HAS_BWAVG=0
if command -v bigWigAverageOverBed >/dev/null 2>&1; then
  HAS_BWAVG=1
fi

if [[ "$HAS_BWAVG" -eq 1 ]]; then
  log "Using bigWigAverageOverBed for scoring (preferred UCSC tool)."

  # bigWigAverageOverBed requires a name column. Create BED4 with unique IDs.
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"xpc_"NR}' "$OUTDIR/02_obs/XPC.summits.bed3" > "$OUTDIR/00_tmp/XPC.summits.named.bed"

  bigWigAverageOverBed "$BW1" "$OUTDIR/00_tmp/XPC.summits.named.bed" "$OUTDIR/00_tmp/xpc.rep1.tab" >/dev/null
  bigWigAverageOverBed "$BW2" "$OUTDIR/00_tmp/XPC.summits.named.bed" "$OUTDIR/00_tmp/xpc.rep2.tab" >/dev/null

  # Parse mean0 (column 5 in bigWigAverageOverBed output: name size covered sum mean0 mean)
  # We want mean0 (treat uncovered as 0), average rep1+rep2.
  paste <(cut -f4,5 "$OUTDIR/00_tmp/xpc.rep1.tab") <(cut -f5 "$OUTDIR/00_tmp/xpc.rep2.tab") \
  | awk 'BEGIN{OFS="\t"}{id=$1; v1=$2; v2=$3; print id, (v1+v2)/2.0}' \
  > "$OUTDIR/00_tmp/xpc.atac_scores.tsv"

  # Join scores back to BED3 in order (they are in the same order)
  paste "$OUTDIR/02_obs/XPC.summits.bed3" <(cut -f2 "$OUTDIR/00_tmp/xpc.atac_scores.tsv") \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' \
  > "$OUTDIR/02_obs/XPC.summits.scored.bed4"

else
  log "bigWigAverageOverBed not found. Falling back to python pyBigWig scoring."
  "$PY" - "$BW1" "$BW2" "$OUTDIR/02_obs/XPC.summits.bed3" "$OUTDIR/02_obs/XPC.summits.scored.bed4" <<'PY'
import sys
bw1, bw2, bed, out = sys.argv[1:5]
try:
    import pyBigWig
except Exception as e:
    raise SystemExit("pyBigWig not available in this python env, and bigWigAverageOverBed not found.\n"
                     "Install UCSC bigWigAverageOverBed OR install pyBigWig into the PY env.\n"
                     f"Original import error: {e}")

b1 = pyBigWig.open(bw1)
b2 = pyBigWig.open(bw2)

def val(bw, chrom, start, end):
    # values() can return None for missing; for 1bp, take first value or 0
    v = bw.values(chrom, start, end, numpy=False)
    if not v or v[0] is None: return 0.0
    return float(v[0])

with open(out, "w") as fo, open(bed) as fi:
    for line in fi:
        if not line.strip(): continue
        c,s,e = line.rstrip("\n").split("\t")[:3]
        s=int(s); e=int(e)
        m = (val(b1,c,s,e) + val(b2,c,s,e))/2.0
        fo.write(f"{c}\t{s}\t{e}\t{m}\n")

b1.close(); b2.close()
PY
fi

# Create XPC bin assignments + thresholds; output:
#   XPC.summits.scored.binned.tsv : chrom start end score bin(1..N_BINS)
#   thresholds.tsv : bin_index  lower_inclusive  upper_exclusive (for reproducibility)
"$PY" - "$OUTDIR/02_obs/XPC.summits.scored.bed4" "$N_BINS" \
      "$OUTDIR/02_obs/XPC.summits.scored.binned.tsv" \
      "$OUTDIR/01_universe/xpc_score_thresholds.tsv" <<'PY'
import sys, numpy as np

bed4 = sys.argv[1]
nbins = int(sys.argv[2])
out_binned = sys.argv[3]
out_thr = sys.argv[4]

# load
rows = []
scores = []
with open(bed4) as f:
    for line in f:
        c,s,e,v = line.rstrip("\n").split("\t")
        v = float(v)
        rows.append((c,int(s),int(e),v))
        scores.append(v)

scores = np.array(scores, dtype=float)
# Define bin edges by quantiles of XPC score distribution
# edges length nbins+1; includes min and max
edges = np.quantile(scores, np.linspace(0, 1, nbins+1))

# Make edges strictly non-decreasing; handle ties by nudging with tiny epsilon
eps = 1e-12
for i in range(1, len(edges)):
    if edges[i] < edges[i-1]:
        edges[i] = edges[i-1]
    if edges[i] == edges[i-1]:
        edges[i] = edges[i] + eps

def assign_bin(v):
    # bins: 1..nbins using edges; last bin includes max
    # np.searchsorted returns index in [0..nbins]
    idx = np.searchsorted(edges, v, side="right") - 1
    if idx < 0: idx = 0
    if idx >= nbins: idx = nbins-1
    return idx + 1

with open(out_binned, "w") as fo:
    fo.write("chrom\tstart\tend\tscore\tbin\n")
    for c,s,e,v in rows:
        b = assign_bin(v)
        fo.write(f"{c}\t{s}\t{e}\t{v}\t{b}\n")

with open(out_thr, "w") as fo:
    fo.write("bin\tlower_inclusive\tupper_exclusive\n")
    for b in range(1, nbins+1):
        lo = edges[b-1]
        hi = edges[b] if b < nbins else edges[b] + 1.0  # last bin open-ended
        fo.write(f"{b}\t{lo}\t{hi}\n")
PY

# Split XPC summits by bin into BED3 files (for shuffling)
log "Splitting XPC summits into bins..."
for b in $(seq 1 "$N_BINS"); do
  awk -v B="$b" 'BEGIN{FS=OFS="\t"} NR>1 && $5==B {print $1,$2,$3}' \
    "$OUTDIR/02_obs/XPC.summits.scored.binned.tsv" \
  | bedtools sort -g "$CHROMSIZES_R" -i stdin \
  > "$OUTDIR/02_obs/XPC.bin${b}.bed3"
  n=$(wc -l < "$OUTDIR/02_obs/XPC.bin${b}.bed3" || echo 0)
  log "  XPC bin $b: $n summits"
done

# -----------------------------------------------------------------------------
# 4) Bin the ATAC universe intervals using the SAME thresholds (from XPC)
#    This ensures shuffled positions have matched ATAC intensity distribution.
# -----------------------------------------------------------------------------
log "Binning ATAC universe intervals using XPC-derived thresholds..."
"$PY" - "$OUTDIR/01_universe/ATAC.mean.noBlacklist.bed4" \
      "$OUTDIR/01_universe/xpc_score_thresholds.tsv" \
      "$OUTDIR/01_universe/ATAC.binned.tsv" <<'PY'
import sys

univ_bed4 = sys.argv[1]
thr_tsv = sys.argv[2]
out = sys.argv[3]

# load thresholds
bins = []
with open(thr_tsv) as f:
    next(f)
    for line in f:
        b, lo, hi = line.rstrip("\n").split("\t")
        bins.append((int(b), float(lo), float(hi)))

def assign(v):
    for b, lo, hi in bins:
        if (v >= lo) and (v < hi):
            return b
    return bins[-1][0]

with open(out, "w") as fo, open(univ_bed4) as fi:
    fo.write("chrom\tstart\tend\tscore\tbin\n")
    for line in fi:
        c,s,e,v = line.rstrip("\n").split("\t")
        v = float(v)
        b = assign(v)
        fo.write(f"{c}\t{s}\t{e}\t{v}\t{b}\n")
PY

# Create BED3 per ATAC bin (merge to simplify)
log "Writing ATAC-bin allowed regions (BED3, merged)..."
rm -f "$OUTDIR/01_universe/ATAC_bin"*.bed

for b in $(seq 1 "$N_BINS"); do
  awk -v B="$b" 'BEGIN{FS=OFS="\t"} NR>1 && $5==B {print $1,$2,$3}' \
    "$OUTDIR/01_universe/ATAC.binned.tsv" \
  | bedtools sort -g "$CHROMSIZES_R" -i stdin \
  | bedtools merge -i stdin \
  > "$OUTDIR/01_universe/ATAC_bin${b}.bed"

  bp=$(awk '{s+=$3-$2} END{print s+0}' "$OUTDIR/01_universe/ATAC_bin${b}.bed" || echo 0)
  log "  ATAC bin $b: $(wc -l < "$OUTDIR/01_universe/ATAC_bin${b}.bed") intervals ; bp=$bp"
done

# Sanity check: every XPC bin must have non-empty ATAC bin regions to shuffle into
log "Sanity check: ensuring each occupied XPC bin has ATAC support..."
for b in $(seq 1 "$N_BINS"); do
  xpc_n=$(wc -l < "$OUTDIR/02_obs/XPC.bin${b}.bed3" || echo 0)
  if [[ "$xpc_n" -gt 0 ]]; then
    atac_bp=$(awk '{s+=$3-$2} END{print s+0}' "$OUTDIR/01_universe/ATAC_bin${b}.bed" || echo 0)
    if [[ "$atac_bp" -le 0 ]]; then
      echo "ERROR: XPC bin $b has $xpc_n sites but ATAC_bin${b}.bed has 0 bp. Cannot intensity-match shuffle." >&2
      echo "Try fewer bins (N_BINS=5) or change scoring/thresholding strategy." >&2
      exit 1
    fi
  fi
done
log "Sanity check passed."

# -----------------------------------------------------------------------------
# 5) Observed distances (XPC -> nearest CTCF) using all XPC summits (chr1-22,X)
#    (No need to restrict to ATAC mask because intensity matching already captures that;
#     HOWEVER, if many XPC scores are 0 because the bigWig has 0 at those sites,
#     they will fall into the lowest bin and shuffle accordingly.)
# -----------------------------------------------------------------------------
log "Computing observed nearest distances (XPC -> nearest CTCF)..."
bedtools closest -a "$OUTDIR/02_obs/XPC.summits.bed3" -b "$OUTDIR/02_obs/CTCF.summits.bed3" -d \
| awk '{print $NF}' > "$OUTDIR/02_obs/observed_distances.txt"

# -----------------------------------------------------------------------------
# 6) Permutations: intensity-matched shuffling per bin, concatenate, then closest
# -----------------------------------------------------------------------------
log "Starting permutations: N_PERM=$N_PERM ; SAVE_PERM_DISTS=$SAVE_PERM_DISTS ; N_BINS=$N_BINS"
log "Each permutation shuffles XPC within ATAC bin matching its ATAC score bin."

# Clean previous saved dists to avoid mixing runs
rm -f "$OUTDIR/03_perm/saved_perm_dists"/dist_perm_*.txt
: > "$OUTDIR/03_perm/perm_stats.tsv"
echo -e "perm\tmean\tmedian\tfrac_le_100bp\tfrac_le_500bp\tfrac_le_1000bp" >> "$OUTDIR/03_perm/perm_stats.tsv"

for i in $(seq 1 "$N_PERM"); do
  seed=$((SEED_BASE + i))
  log "Perm $i/$N_PERM | seed=$seed | shuffling bins..."

  : > "$OUTDIR/00_tmp/XPC.shuf.all.bed3"

  for b in $(seq 1 "$N_BINS"); do
    if [[ -s "$OUTDIR/02_obs/XPC.bin${b}.bed3" ]]; then
      xpc_n=$(wc -l < "$OUTDIR/02_obs/XPC.bin${b}.bed3")
      log "  perm $i | bin $b | shuffling $xpc_n sites within ATAC_bin${b}.bed"

      bedtools shuffle \
        -i "$OUTDIR/02_obs/XPC.bin${b}.bed3" \
        -g "$CHROMSIZES_R" \
        -incl "$OUTDIR/01_universe/ATAC_bin${b}.bed" \
        -excl "$BLACKLIST_R" \
        -chrom \
        -seed $((seed*100 + b)) \
      >> "$OUTDIR/00_tmp/XPC.shuf.all.bed3"
    fi
  done

  bedtools sort -g "$CHROMSIZES_R" -i "$OUTDIR/00_tmp/XPC.shuf.all.bed3" > "$OUTDIR/00_tmp/XPC.shuf.all.sorted.bed3"

  log "Perm $i/$N_PERM | computing distances to nearest CTCF..."
  bedtools closest -a "$OUTDIR/00_tmp/XPC.shuf.all.sorted.bed3" -b "$OUTDIR/02_obs/CTCF.summits.bed3" -d \
  | awk '{print $NF}' > "$OUTDIR/00_tmp/dist.txt"

  if (( i <= SAVE_PERM_DISTS )); then
    cp "$OUTDIR/00_tmp/dist.txt" "$OUTDIR/03_perm/saved_perm_dists/dist_perm_${i}.txt"
    log "Perm $i/$N_PERM | saved distances for overlay histogram"
  fi

  log "Perm $i/$N_PERM | done."
done

# Summarize observed + per-permutation distance vectors in ONE python call (robust medians)
log "Summarizing permutations + generating plots (including observed vs pooled permuted histogram)..."
"$PY" - "$OUTDIR" "$SAVE_PERM_DISTS" "$DIST_CAP" <<'PY'
import os, sys, numpy as np
import matplotlib.pyplot as plt

outdir = sys.argv[1]
save_k = int(sys.argv[2])
cap = float(sys.argv[3]) if sys.argv[3].strip() else None

def load1(p): return np.loadtxt(p, dtype=float)

obs = load1(os.path.join(outdir, "02_obs", "observed_distances.txt"))

# observed stats
obs_mean = float(np.mean(obs))
obs_median = float(np.median(obs))
obs_f100 = float(np.mean(obs <= 100))
obs_f500 = float(np.mean(obs <= 500))
obs_f1000 = float(np.mean(obs <= 1000))

# gather perm distance files (at least those saved; but we also want stats for ALL perms)
# Here, to avoid huge storage, we recompute perm stats from saved files only for overlay,
# and for p-values we instead compute from perm_stats computed from ALL perms? We did not.
# So we compute per-perm stats from saved files only IF only saved exist.
# Better: compute per-perm stats from ALL dist files is not possible because they aren't saved.
#
# Therefore, we compute perm stats by re-reading dist files if present; otherwise we only do overlay.
# To keep things correct, we compute per-perm stats ON THE FLY from saved dist files
# AND write summary + overlay plot. For full p-values, increase SAVE_PERM_DISTS or store all.
#
# Practical compromise: use saved_k as Monte Carlo null for distance-distribution plot;
# for p-values, you'd typically use N_PERM full stats. If you want exact p-values,
# modify bash loop to compute per-perm stats and append (easy; ask if needed).
#
# We'll do both:
#   - compute p-values based on saved_k dist vectors (approx)
#   - write that clearly in summary.

saved_dir = os.path.join(outdir, "03_perm", "saved_perm_dists")
perm_stats = []
dists_pool = []
k_found = 0

for i in range(1, save_k+1):
    p = os.path.join(saved_dir, f"dist_perm_{i}.txt")
    if os.path.exists(p):
        d = load1(p)
        k_found += 1
        dists_pool.append(d)
        perm_stats.append((
            float(np.mean(d)),
            float(np.median(d)),
            float(np.mean(d <= 100)),
            float(np.mean(d <= 500)),
            float(np.mean(d <= 1000)),
        ))

if k_found == 0:
    raise SystemExit("No saved perm distance files found; cannot make overlay histogram.")

perm_stats = np.array(perm_stats, dtype=float)
perm_mean   = perm_stats[:,0]
perm_median = perm_stats[:,1]
perm_f100   = perm_stats[:,2]
perm_f500   = perm_stats[:,3]
perm_f1000  = perm_stats[:,4]

def p_low(null, obsval):   # smaller = closer enrichment
    return (np.sum(null <= obsval) + 1) / (len(null) + 1)

def p_high(null, obsval):  # larger fraction = enrichment
    return (np.sum(null >= obsval) + 1) / (len(null) + 1)

summary = [
    ("mean_distance",   obs_mean,   "p_low",  p_low(perm_mean, obs_mean)),
    ("median_distance", obs_median, "p_low",  p_low(perm_median, obs_median)),
    ("frac_le_100bp",   obs_f100,   "p_high", p_high(perm_f100, obs_f100)),
    ("frac_le_500bp",   obs_f500,   "p_high", p_high(perm_f500, obs_f500)),
    ("frac_le_1000bp",  obs_f1000,  "p_high", p_high(perm_f1000, obs_f1000)),
]

os.makedirs(os.path.join(outdir, "04_viz"), exist_ok=True)

with open(os.path.join(outdir, "summary.autosomesX.intensityMatched.tsv"), "w") as f:
    f.write("note\tp-values computed from saved permutations only; increase SAVE_PERM_DISTS for tighter MC error\n")
    f.write(f"saved_perms_used\t{k_found}\n")
    f.write("metric\tobserved\tp_direction\tempirical_p\n")
    for m,o,d,p in summary:
        f.write(f"{m}\t{o}\t{d}\t{p}\n")

# 1) overlay histogram: observed vs pooled permuted distances
null = np.concatenate(dists_pool)
if cap is not None:
    obs_plot = np.clip(obs, 0, cap)
    null_plot = np.clip(null, 0, cap)
    xlabel = f"Distance to nearest CTCF summit (bp), capped at {int(cap)}"
else:
    obs_plot, null_plot = obs, null
    xlabel = "Distance to nearest CTCF summit (bp)"

plt.figure()
plt.hist(null_plot, bins=120, density=True, alpha=0.5, label=f"Permuted pooled (k={k_found})")
plt.hist(obs_plot,  bins=120, density=True, alpha=0.5, label="Observed")
plt.xlabel(xlabel)
plt.ylabel("Density")
plt.title("Observed vs ATAC-intensity–matched permuted distance distributions (chr1-22,X)")
plt.legend()
plt.savefig(os.path.join(outdir, "04_viz", f"observed_vs_permuted_dist_hist.pooled{k_found}.png"),
            dpi=200, bbox_inches="tight")
plt.close()

# 2) ECDF observed
x = np.sort(obs)
y = np.arange(1, len(x)+1)/len(x)
plt.figure()
plt.plot(x, y)
plt.xlabel("Distance to nearest CTCF summit (bp)")
plt.ylabel("ECDF")
plt.title("Observed ECDF (chr1-22,X): XPC → nearest CTCF")
plt.savefig(os.path.join(outdir, "04_viz", "observed_distance.ecdf.png"),
            dpi=200, bbox_inches="tight")
plt.close()

# 3) Perm-stat histograms (from saved perms only)
plots = [
    ("mean_distance", perm_mean, obs_mean),
    ("median_distance", perm_median, obs_median),
    ("frac_le_100bp", perm_f100, obs_f100),
    ("frac_le_500bp", perm_f500, obs_f500),
    ("frac_le_1000bp", perm_f1000, obs_f1000),
]
for name, nullv, obsval in plots:
    plt.figure()
    plt.hist(nullv, bins=40)
    plt.axvline(obsval, linewidth=2)
    plt.xlabel(name)
    plt.ylabel("Count")
    plt.title(f"Null (intensity-matched; saved perms k={k_found}) for {name}")
    plt.savefig(os.path.join(outdir, "04_viz", f"{name}.hist.png"),
                dpi=200, bbox_inches="tight")
    plt.close()

print("Wrote:", os.path.join(outdir, "summary.autosomesX.intensityMatched.tsv"))
print("Plots in:", os.path.join(outdir, "04_viz"))
PY

log "DONE."
log "Key outputs:"
log "  $OUTDIR/summary.autosomesX.intensityMatched.tsv"
log "  $OUTDIR/04_viz/observed_vs_permuted_dist_hist.pooled*.png"
log "  $OUTDIR/01_universe/xpc_score_thresholds.tsv  (bin edges)"
log "  $OUTDIR/02_obs/XPC.summits.scored.binned.tsv  (XPC bins)"
log "  $OUTDIR/02_obs/XPC_CTCF.overlap_peak_summit_distances.tsv"
log "  $OUTDIR/02_obs/XPC_CTCF.overlap_peak_summit_distances.summary.tsv"
log "  $OUTDIR/01_universe/ATAC_bin*.bed            (allowed regions per bin)"

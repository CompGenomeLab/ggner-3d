# This script processes chromosome-level contact information from multiple `.mcool`
# files using GENOVA. For each sample, a balanced contact map is loaded at 1 Mb
# resolution, and a chromosome-by-chromosome contact matrix is generated. The
# observed contact frequencies are then normalized against an expected contact
# matrix to calculate log2 observed/expected values.
#
# The resulting matrices summarize whether contacts between chromosome pairs occur
# more or less frequently than expected. Positive values indicate contacts enriched
# relative to expectation, while negative values indicate contacts depleted relative
# to expectation. The processed matrices are saved for each sample in both RDS
# format, which preserves the R object structure, and TSV format, which can be
# inspected or used by other tools.
#
# Here, expected values are calculated using the "sums" mode in GENOVA. This uses
# a chi-squared-style null model, where chromosome-pair expectations are conditional
# on the total contact sums for individual chromosomes. In other words, given the
# total number of contacts associated with each chromosome, the expected matrix
# estimates how many contacts would be expected between each chromosome pair under
# this chromosome-level null model.

# --- packages ---
library(GENOVA)

# --- inputs ---
RESOLUTION <- 1000000

MCOOLS_PATH_DICT <- c(
  # --- WT ---
  WTnoUV = "/home/carlos/oldies/ner_collab/mcools_merged_biorep/WTnoUV.mcool",
  WT1h   = "/home/carlos/oldies/ner_collab/mcools_merged/WT_1H_1.mcool",
  WT3h   = "/home/carlos/oldies/ner_collab/mcools_merged_biorep/WT3h.mcool",
  WT6h   = "/home/carlos/oldies/ner_collab/mcools_merged/WT_6H_6.mcool",

  # --- XPA ---
  XPAnoUV = "/home/carlos/oldies/ner_collab/mcools_merged_biorep/XPAnoUV.mcool",
  XPA1h   = "/home/carlos/oldies/ner_collab/mcools_merged/XPA_1H_1.mcool",
  XPA3h   = "/home/carlos/oldies/ner_collab/mcools_merged_biorep/XPA3h.mcool",
  XPA6h   = "/home/carlos/oldies/ner_collab/mcools_merged/XPA_6H_1.mcool",

  # --- XPC ---
  XPCnoUV = "/home/carlos/oldies/ner_collab/mcools_merged_biorep/XPCnoUV.mcool",
  XPC3h   = "/home/carlos/oldies/ner_collab/mcools_merged_biorep/XPC3h.mcool"
)

# --- helper: compute log2(O/E) in the same "proportions of total" framework ---
chrommat_log2_obsexp <- function(cm, expected_mode = c("bins","sums","trans","cis","regress"),
                                 finite_to_na = TRUE) {
  expected_mode <- match.arg(expected_mode)

  obs <- cm$obs
  exp <- cm$exp

  if (length(dim(obs)) == 3L) {
    obs <- obs[,,1, drop = TRUE]
    exp <- exp[,,1, drop = TRUE]
  }

  tmp <- obs
  if (expected_mode == "trans") diag(tmp) <- 0
  obs_prop <- obs / sum(tmp, na.rm = TRUE)

  out <- log2(obs_prop / exp)

  if (finite_to_na) out[!is.finite(out)] <- NA_real_
  out
}

# --- run ---
EXPECTED_MODE <- "sums"
OUTDIR <- paste0("/home/carlos/Clone/ggner-3d/data/genova_trans/log2_obsexp_", EXPECTED_MODE)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

clrs_ <- list()
cm_list <- list()
log2oe_list <- list()

for (s_name in names(MCOOLS_PATH_DICT)) {

  clr <- load_contacts(
    signal_path = MCOOLS_PATH_DICT[[s_name]],
    sample_name = s_name,
    balancing   = TRUE,
    resolution  = RESOLUTION
  )

  cm <- chromosome_matrix(clr, expected = EXPECTED_MODE)

  log2oe <- chrommat_log2_obsexp(cm, expected_mode = EXPECTED_MODE)

  clrs_[[s_name]]     <- clr
  cm_list[[s_name]]   <- cm
  log2oe_list[[s_name]] <- log2oe

  saveRDS(
    log2oe,
    file = file.path(OUTDIR, paste0(s_name, "_log2_obsexp.rds"))
  )

  write.table(
    log2oe,
    file = file.path(OUTDIR, paste0(s_name, "_log2_obsexp.tsv")),
    sep = "\t", quote = FALSE, col.names = NA
  )
}

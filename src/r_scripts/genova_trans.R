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
# chromosome_matrix() returns obs as counts and exp as proportions-of-total expected :contentReference[oaicite:0]{index=0}
# For expected="trans", GENOVA zeros the diagonal before computing totals/expected :contentReference[oaicite:1]{index=1}
chrommat_log2_obsexp <- function(cm, expected_mode = c("bins","sums","trans","cis","regress"),
                                 finite_to_na = TRUE) {
  expected_mode <- match.arg(expected_mode)

  obs <- cm$obs
  exp <- cm$exp

  # if you passed a single contacts object, these are still 3D arrays with 3rd dim = 1 :contentReference[oaicite:2]{index=2}
  if (length(dim(obs)) == 3L) {
    obs <- obs[,,1, drop = TRUE]
    exp <- exp[,,1, drop = TRUE]
  }

  # compute obs as proportion-of-total (with the same diagonal handling used for "trans" expected) :contentReference[oaicite:3]{index=3}
  tmp <- obs
  if (expected_mode == "trans") diag(tmp) <- 0
  obs_prop <- obs / sum(tmp, na.rm = TRUE)

  out <- log2(obs_prop / exp)

  if (finite_to_na) out[!is.finite(out)] <- NA_real_
  out
}

# --- run ---
EXPECTED_MODE <- "sums"   # you used expected="trans" in chromosome_matrix()
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

  # store in memory (optional)
  clrs_[[s_name]]     <- clr
  cm_list[[s_name]]   <- cm
  log2oe_list[[s_name]] <- log2oe

  # save as RDS (keeps dimnames clean)
  saveRDS(
    log2oe,
    file = file.path(OUTDIR, paste0(s_name, "_log2_obsexp.rds"))
  )

  # optional: also save as TSV matrix
  write.table(
    log2oe,
    file = file.path(OUTDIR, paste0(s_name, "_log2_obsexp.tsv")),
    sep = "\t", quote = FALSE, col.names = NA
  )
}

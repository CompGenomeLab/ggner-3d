#!/usr/bin/env Rscript

# -----------------------------
# Parallel Grid Search + Best Scenario Data Export
# -----------------------------

# Load the argparse library
suppressPackageStartupMessages(library(argparse))

parse_numeric_seq_arg <- function(raw_values, arg_name) {
  if (is.null(raw_values) || length(raw_values) == 0L) {
    stop(sprintf("%s must contain at least one value.", arg_name), call. = FALSE)
  }

  raw_text <- paste(raw_values, collapse = ",")
  raw_text <- gsub("^\\s*c\\((.*)\\)\\s*$", "\\1", raw_text, perl = TRUE)

  tokens <- unlist(strsplit(raw_text, "[,\\s]+", perl = TRUE), use.names = FALSE)
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) == 0L) {
    stop(sprintf("%s must contain at least one value.", arg_name), call. = FALSE)
  }

  parsed_values <- suppressWarnings(as.numeric(tokens))
  bad_values <- is.na(parsed_values) | !is.finite(parsed_values)
  if (any(bad_values)) {
    stop(
      sprintf(
        "%s contains non-numeric values: %s",
        arg_name,
        paste(unique(tokens[bad_values]), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  parsed_values
}

# Create a parser object
parser <- ArgumentParser(description = "Process input, summary, and SMD files with covariate options and export best scenario data.")

# Add command-line arguments
parser$add_argument(
  "--input_file",
  required = TRUE,
  help = "Path to the input CSV file."
)

parser$add_argument(
  "--summary_file",
  required = TRUE,
  help = "Path to the summary CSV file (grid search results)."
)

parser$add_argument(
  "--smd_file",
  required = TRUE,
  help = "Path to the SMD CSV file."
)

parser$add_argument(
  "--output_data",
  required = TRUE,
  help = "Path to save the rows/features for the best scenario (CSV)."
)

parser$add_argument(
  "--covariates",
  nargs = "+",
  default = c("n_ctcf_peaks", "looped_bin"),
  help = "List of covariate names (space-separated)."
)

parser$add_argument(
  "--critical_exact_match_cols",
  nargs = "*",
  default = c("n_ctcf_peaks", "looped_bin"),
  help = "Columns that require exact matching. Pass the flag with no values to disable exact matching."
)

parser$add_argument(
  "--drop_failed_strata",
  action = "store_true",
  default = FALSE,
  help = "If set, any stratum without a viable match invalidates the scenario. Otherwise, those strata are skipped."
)

parser$add_argument(
  "--quiet_grid_search",
  action = "store_true",
  default = FALSE,
  help = "If set, suppress per-scenario grid-search logs."
)

parser$add_argument(
  "--workers",
  type = "integer",
  default = 0L,
  help = "Parallel workers for grid search. Use 0 for auto (detected cores - 1)."
)

parser$add_argument(
  "--focal_keep_frac_seq",
  nargs = "+",
  default = c("1.00", "0.90", "0.80", "0.70", "0.60", "0.50", "0.40", "0.30", "0.20"),
  help = "Focal keep fractions to test. Accepts space-separated values, comma-separated values, or quoted c(...)."
)

# Parse the arguments
args <- parser$parse_args()

# Assign to variables
input_file <- args$input_file
summary_file <- args$summary_file
smd_file <- args$smd_file
output_data_file <- args$output_data
covariate_names <- args$covariates
critical_exact_match_cols <- args$critical_exact_match_cols
drop_failed_strata <- args$drop_failed_strata
verbose_grid_search <- !args$quiet_grid_search
requested_workers <- as.integer(args$workers)
if (is.na(requested_workers) || requested_workers < 0L) {
  requested_workers <- 0L
}
focal_keep_frac_seq_raw <- parse_numeric_seq_arg(args$focal_keep_frac_seq, "--focal_keep_frac_seq")
if (any(focal_keep_frac_seq_raw <= 0)) {
  stop("--focal_keep_frac_seq values must be > 0.", call. = FALSE)
}
if (any(focal_keep_frac_seq_raw > 1)) {
  cat("Warning: --focal_keep_frac_seq values > 1 were capped at 1.\n")
}
focal_keep_frac_seq <- pmin(1, focal_keep_frac_seq_raw)

detected_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}
auto_workers <- max(1L, detected_cores - 1L)

cat("Running with the following arguments:\n")
cat("Input file:", input_file, "\n")
cat("Summary file:", summary_file, "\n")
cat("SMD file:", smd_file, "\n")
cat("Output Data file:", output_data_file, "\n")
cat("Covariates:", paste(covariate_names, collapse = ", "), "\n")
exact_match_cols_label <- if (length(critical_exact_match_cols) == 0L) "(none)" else paste(critical_exact_match_cols, collapse = ", ")
cat("Exact-match columns:", exact_match_cols_label, "\n")
cat("Drop failed strata: ", drop_failed_strata, "\n\n")
cat("Verbose grid search:", verbose_grid_search, "\n\n")
cat("Requested workers (0=auto):", requested_workers, "\n")
cat("Detected CPU cores:", detected_cores, "\n")
cat("Focal keep fraction seq:", paste(focal_keep_frac_seq, collapse = ", "), "\n\n")

# define quantile grid
q_lo_seq = seq(0.60, 0.75, by = 0.05)
q_hi_seq = seq(0.75, 0.90, by = 0.05)
lo_hi_diff_min = 0.10

set.seed(12345)
options(dplyr.summarise.inform = FALSE)

required_cran <- c("dplyr", "readr", "purrr", "tibble", "data.table", "ks", "tidyr", "forcats", "stringr")
required_bioc <- c("GenomicRanges", "nullranges")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(GenomicRanges)
  library(nullranges)
})

# -----------------------------
# Helpers
# -----------------------------
safe_div <- function(num, den, na_on_zero = TRUE, eps = NULL) {
  if (!is.null(eps)) {
    return(num / (den + eps))
  }
  if (na_on_zero) {
    return(ifelse(den == 0 | is.na(den), NA_real_, num / den))
  }
  ifelse(den == 0, 0, num / den)
}

format_duration <- function(seconds) {
  if (is.na(seconds) || !is.finite(seconds)) {
    return("NA")
  }
  total_sec <- as.integer(round(max(0, seconds)))
  hrs <- total_sec %/% 3600L
  mins <- (total_sec %% 3600L) %/% 60L
  secs <- total_sec %% 60L
  sprintf("%02d:%02d:%02d", hrs, mins, secs)
}

scale_covariate <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }
  sd_val <- stats::sd(x, na.rm = TRUE)
  if (is.na(sd_val) || sd_val == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

compute_smd_table <- function(data, covariates, treat_col = "treatment") {
  treat_vec <- as.numeric(data[[treat_col]])
  stopifnot(all(treat_vec %in% c(0, 1)))
  purrr::map_dfr(covariates, function(var_name) {
    values <- data[[var_name]]
    treated_vals <- values[treat_vec == 1]
    control_vals <- values[treat_vec == 0]
    mean_t <- mean(treated_vals, na.rm = TRUE)
    mean_c <- mean(control_vals, na.rm = TRUE)
    var_t <- stats::var(treated_vals, na.rm = TRUE)
    var_c <- stats::var(control_vals, na.rm = TRUE)
    pooled_sd <- sqrt((var_t + var_c) / 2)
    smd <- dplyr::case_when(
      is.na(pooled_sd) ~ NA_real_,
      pooled_sd == 0 ~ NA_real_,
      TRUE ~ (mean_t - mean_c) / pooled_sd
    )
    tibble::tibble(
      covariate = var_name,
      mean_treated = mean_t,
      mean_control = mean_c,
      smd = smd
    )
  })
}

# -----------------------------
# Config
# -----------------------------
matching_methods <- c("rejection", "stratified", "nearest") # order of preference

quantile_grid <- expand.grid(
  q_lo = q_lo_seq,
  q_hi = q_hi_seq,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    q_lo = as.numeric(q_lo),
    q_hi = as.numeric(q_hi)
  ) %>%
  filter(q_lo < q_hi - lo_hi_diff_min)

if (verbose_grid_search) {
  cat(
    "[grid] candidates=", nrow(quantile_grid),
    " q_lo_range=[", min(quantile_grid$q_lo), ", ", max(quantile_grid$q_lo), "]",
    " q_hi_range=[", min(quantile_grid$q_hi), ", ", max(quantile_grid$q_hi), "]",
    " min_diff>", lo_hi_diff_min, "\n",
    sep = ""
  )
}

# -----------------------------
# Data load + preprocessing
# -----------------------------
df_raw <- readr::read_csv(
  input_file,
  col_types = cols(
    chrom = col_character(),
    start = col_double(),
    end = col_double(),
    log2_insulation_score = col_double(),
    xpc_signal = col_double(),
    ctcf_signal = col_double(),
    atac_signal = col_double(),
    n_ctcf_motif = col_double(),
    looped = col_logical(),
    n_ctcf_peaks = col_double(),
    is_boundary = col_logical(),
    h3k4me1_signal = col_double(),
    h3k4me3_signal = col_double(),
    h3k27ac_signal = col_double(),
    pc1 = col_double(),
    damageseq = col_double(),
    pol2 = col_double()
  )
)

analysis_df <- df_raw %>%
  mutate(
    n_ctcf_peaks = as.integer(n_ctcf_peaks),
    ins_eff = -log2_insulation_score,
    looped_bin = as.integer(looped),
    is_boundary_bin = as.integer(is_boundary)
  ) %>%
  filter(
    is.finite(n_ctcf_peaks),
    is.finite(ins_eff),
    is.finite(ins_eff),
    is.finite(xpc_signal),
    is.finite(ctcf_signal),
    n_ctcf_peaks > 0
  ) %>%
  filter(!is.na(looped_bin), !is.na(is_boundary_bin))

if (nrow(analysis_df) == 0) {
  cat("No loci remain after filtering for matching.\n")
  quit(save = "no", status = 0)
}

# -----------------------------
# Matching preparation
# -----------------------------
prepare_matching_pool <- function(base_df, q_lo, q_hi, covariates, exact_cols) {
  q_vals <- stats::quantile(base_df$xpc_signal, probs = c(q_lo, q_hi), na.rm = TRUE, names = FALSE)
  dataset <- base_df %>%
    mutate(
      xpc_group = dplyr::case_when(
        xpc_signal >= q_vals[2] ~ "high",
        xpc_signal <= q_vals[1] ~ "low",
        TRUE ~ NA_character_
      ),
      treatment = dplyr::if_else(xpc_group == "high", 1L, 0L, missing = NA_integer_)
    ) %>%
    filter(!is.na(xpc_group))
  if (nrow(dataset) == 0L || length(unique(dataset$treatment)) < 2L) {
    return(NULL)
  }
  dataset <- dataset %>%
    mutate(
      row_id = dplyr::row_number(),
      across(all_of(covariates), as.numeric),
      across(all_of(covariates), scale_covariate, .names = "scaled_{.col}")
    )
  if (length(exact_cols) > 0) {
    exact_df <- dataset %>%
      select(all_of(exact_cols)) %>%
      mutate(across(everything(), as.factor))
    critical_key <- interaction(exact_df, drop = TRUE)
  } else {
    critical_key <- factor(rep("all", nrow(dataset)))
  }
  dataset <- dataset %>%
    mutate(critical_key = critical_key)
  scaled_covariate_names <- paste0("scaled_", covariates)
  covar_formula <- if (length(scaled_covariate_names) > 0) {
    stats::as.formula(paste("~", paste(scaled_covariate_names, collapse = " + ")))
  } else {
    ~1
  }
  gr <- GenomicRanges::makeGRangesFromDataFrame(
    dataset %>%
      mutate(
        start_bp = as.integer(start),
        end_bp = as.integer(end),
        start_1based = start_bp + 1L
      ) %>%
      select(-start, -end) %>%
      select(chrom, start_1based, end_bp, everything()),
    seqnames.field = "chrom",
    start.field = "start_1based",
    end.field = "end_bp",
    keep.extra.columns = TRUE
  )
  mcols(gr)$critical_key <- dataset$critical_key
  list(
    data = dataset,
    gr = gr,
    q_vals = q_vals,
    covar_formula = covar_formula
  )
}

match_strata <- function(gr, covar_formula, methods, drop_failed_strata) {
  focal_ranges <- gr[gr$treatment == 1L]
  pool_ranges  <- gr[gr$treatment == 0L]
  if (length(focal_ranges) == 0L || length(pool_ranges) == 0L) {
    return(NULL)
  }
  strata_vals <- sort(unique(as.character(gr$critical_key)))
  matched_focal_list <- list()
  matched_control_list <- list()
  pair_offset <- 0L
  method_track <- list()

  attempt_match <- function(focal_stratum, pool_stratum, method, replace_flag) {
    tryCatch(
      {
        match_obj <- nullranges::matchRanges(
          focal = focal_stratum,
          pool = pool_stratum,
          covar = covar_formula,
          method = method,
          replace = replace_flag
        )
        list(
          focal_df = tibble::as_tibble(nullranges::focal(match_obj)),
          control_df = tibble::as_tibble(nullranges::matched(match_obj)),
          method = method,
          replace = replace_flag
        )
      },
      error = function(e) NULL
    )
  }

  for (stratum in strata_vals) {
    focal_stratum <- focal_ranges[focal_ranges$critical_key == stratum]
    pool_stratum  <- pool_ranges[pool_ranges$critical_key == stratum]
    if (length(focal_stratum) == 0L) {
      next
    }
    if (length(pool_stratum) == 0L) {
      if (drop_failed_strata) {
        return(NULL)
      } else {
        message("Skipping stratum ", stratum, " due to zero available controls.")
        next
      }
    }
    match_success <- NULL
    for (replace_flag in c(FALSE, TRUE)) {
      if (!is.null(match_success)) break
      for (method in methods) {
        match_attempt <- attempt_match(focal_stratum, pool_stratum, method, replace_flag)
        if (is.null(match_attempt)) {
          next
        }
        focal_df <- match_attempt$focal_df
        control_df <- match_attempt$control_df
        if (nrow(focal_df) == 0 || nrow(control_df) == 0) next
        if (nrow(focal_df) != nrow(control_df)) next
        match_success <- match_attempt
        break
      }
    }
    if (is.null(match_success)) {
      if (drop_failed_strata) {
        return(NULL)
      } else {
        message("Skipping stratum ", stratum, " after all matching methods failed.")
        next
      }
    }
    focal_df <- match_success$focal_df
    control_df <- match_success$control_df
    if (nrow(focal_df) == 0 || nrow(control_df) == 0) {
      if (drop_failed_strata) return(NULL)
      message("Skipping stratum ", stratum, " because matched data were empty.")
      next
    }
    if (nrow(focal_df) != nrow(control_df)) {
      if (drop_failed_strata) return(NULL)
      message("Skipping stratum ", stratum, " due to unequal matched counts.")
      next
    }
    pair_ids <- seq_len(nrow(focal_df)) + pair_offset
    focal_df <- focal_df %>% mutate(pair_id = pair_ids, critical_key = stratum)
    control_df <- control_df %>% mutate(pair_id = pair_ids, critical_key = stratum)
    matched_focal_list[[length(matched_focal_list) + 1L]]   <- focal_df
    matched_control_list[[length(matched_control_list) + 1L]] <- control_df
    method_track[[length(method_track) + 1L]] <- tibble::tibble(
      critical_key = stratum,
      method = match_success$method,
      replace = match_success$replace
    )
    pair_offset <- pair_offset + nrow(focal_df)
  }
  if (pair_offset == 0L) return(NULL)

  matched_focal   <- dplyr::bind_rows(matched_focal_list)
  matched_control <- dplyr::bind_rows(matched_control_list)
  method_df <- if (length(method_track) > 0) dplyr::bind_rows(method_track) else {
    tibble::tibble(critical_key = character(), method = character(), replace = logical())
  }
  list(
    matched_focal = matched_focal,
    matched_control = matched_control,
    method_per_stratum = method_df,
    total_pairs_raw = pair_offset
  )
}

run_matching_for_scenario <- function(prep, methods) {
  match_strata(prep$gr, prep$covar_formula, methods, drop_failed_strata)
}

summarize_scenario <- function(prep, match_res, q_lo, q_hi, covariates) {
  matched_focal <- match_res$matched_focal
  matched_control <- match_res$matched_control
  method_per_stratum <- match_res$method_per_stratum
  method_label <- paste(unique(method_per_stratum$method), collapse = ";")
  method_label <- ifelse(nchar(method_label) > 0, method_label, NA_character_)
  replace_flag <- any(method_per_stratum$replace)
  total_controls <- nrow(matched_control)
  unique_controls <- dplyr::n_distinct(matched_control$row_id)

  pairs_df <- matched_focal %>%
    transmute(
      pair_id = pair_id,
      treated_row_id = row_id,
      treated_xpc = xpc_signal,
      critical_key = critical_key
    ) %>%
    inner_join(
      matched_control %>%
        transmute(
          pair_id = pair_id,
          control_row_id = row_id,
          control_xpc = xpc_signal,
          critical_key = critical_key
        ),
      by = c("pair_id","critical_key")
    )

  if (nrow(pairs_df) == 0L) return(NULL)

  pairs_no_repl <- pairs_df %>%
    arrange(control_row_id, pair_id) %>%
    distinct(control_row_id, .keep_all = TRUE)

  matched_pairs <- nrow(pairs_no_repl)
  if (matched_pairs == 0L) return(NULL)

  # treated / control sets used for SMD
  pair_ids_tbl <- tibble::tibble(pair_id = pairs_no_repl$pair_id)

  treated_no_repl <- matched_focal %>%
    semi_join(pair_ids_tbl, by = "pair_id") %>%
    mutate(treatment = 1L)

  control_no_repl <- matched_control %>%
    semi_join(pair_ids_tbl, by = "pair_id") %>%
    mutate(treatment = 0L)

  matched_for_smd <- dplyr::bind_rows(treated_no_repl, control_no_repl)

  smd_tbl <- compute_smd_table(matched_for_smd, covariates, treat_col = "treatment")
  max_abs_smd  <- if (nrow(smd_tbl) == 0 || all(is.na(smd_tbl$smd))) NA_real_ else max(abs(smd_tbl$smd), na.rm = TRUE)
  mean_abs_smd <- if (nrow(smd_tbl) == 0 || all(is.na(smd_tbl$smd))) NA_real_ else mean(abs(smd_tbl$smd), na.rm = TRUE)

  xpc_treated_mean <- mean(pairs_no_repl$treated_xpc, na.rm = TRUE)
  xpc_control_mean <- mean(pairs_no_repl$control_xpc, na.rm = TRUE)
  xpc_diff <- xpc_treated_mean - xpc_control_mean

  duplicates <- total_controls - unique_controls
  score <- if (is.na(max_abs_smd)) NA_real_ else xpc_diff - max_abs_smd

  used_strata <- unique(pairs_no_repl$critical_key)

  list(
    q_lo = q_lo,
    q_hi = q_hi,
    xpc_lo_value = prep$q_vals[[1]],
    xpc_hi_value = prep$q_vals[[2]],
    method = method_label,
    replace = replace_flag,
    matched_pairs = matched_pairs,
    total_pairs_raw = match_res$total_pairs_raw,
    unique_controls = unique_controls,
    total_controls = total_controls,
    duplicates = duplicates,
    xpc_treated_mean = xpc_treated_mean,
    xpc_control_mean = xpc_control_mean,
    xpc_diff = xpc_diff,
    smd_table = smd_tbl,
    max_abs_smd = max_abs_smd,
    mean_abs_smd = mean_abs_smd,
    score = score,
    focal_keep_frac = prep$focal_keep_frac,
    focal_kept_n = prep$focal_kept_n,
    focal_total_n = prep$focal_total_n,
    # IDs for extraction
    treated_ids = treated_no_repl$row_id,
    matched_control_ids = control_no_repl$row_id,
    used_strata = used_strata
  )
}

shrink_prepared_focal <- function(prep, focal_keep_frac) {
  focal_keep_frac <- as.numeric(focal_keep_frac)
  if (is.na(focal_keep_frac) || focal_keep_frac <= 0) {
    return(NULL)
  }
  focal_keep_frac <- min(1, focal_keep_frac)

  focal_df <- prep$data %>%
    filter(treatment == 1L) %>%
    arrange(desc(xpc_signal), row_id)

  focal_total_n <- nrow(focal_df)
  if (focal_total_n == 0L) {
    return(NULL)
  }

  focal_kept_n <- as.integer(max(1L, ceiling(focal_total_n * focal_keep_frac)))
  focal_kept_n <- min(focal_total_n, focal_kept_n)
  keep_ids <- focal_df$row_id[seq_len(focal_kept_n)]

  data_subset <- prep$data %>%
    filter(treatment == 0L | row_id %in% keep_ids)

  gr_subset <- prep$gr[prep$gr$treatment == 0L | prep$gr$row_id %in% keep_ids]

  if (nrow(data_subset) == 0L || length(gr_subset) == 0L) {
    return(NULL)
  }

  treated_n <- sum(data_subset$treatment == 1L, na.rm = TRUE)
  pool_n <- sum(data_subset$treatment == 0L, na.rm = TRUE)
  if (treated_n == 0L || pool_n == 0L) {
    return(NULL)
  }

  prep_subset <- prep
  prep_subset$data <- data_subset
  prep_subset$gr <- gr_subset
  prep_subset$focal_keep_frac <- focal_kept_n / focal_total_n
  prep_subset$focal_kept_n <- focal_kept_n
  prep_subset$focal_total_n <- focal_total_n
  prep_subset
}

is_better_match_quality <- function(candidate_metrics, current_best_metrics) {
  if (is.null(current_best_metrics)) {
    return(TRUE)
  }
  cand_mean <- dplyr::coalesce(candidate_metrics$mean_abs_smd, Inf)
  best_mean <- dplyr::coalesce(current_best_metrics$mean_abs_smd, Inf)
  if (cand_mean < best_mean) {
    return(TRUE)
  }
  if (cand_mean > best_mean) {
    return(FALSE)
  }
  cand_max <- dplyr::coalesce(candidate_metrics$max_abs_smd, Inf)
  best_max <- dplyr::coalesce(current_best_metrics$max_abs_smd, Inf)
  if (cand_max < best_max) {
    return(TRUE)
  }
  if (cand_max > best_max) {
    return(FALSE)
  }
  cand_xpc <- dplyr::coalesce(candidate_metrics$xpc_diff, -Inf)
  best_xpc <- dplyr::coalesce(current_best_metrics$xpc_diff, -Inf)
  if (cand_xpc > best_xpc) {
    return(TRUE)
  }
  if (cand_xpc < best_xpc) {
    return(FALSE)
  }
  cand_frac <- dplyr::coalesce(candidate_metrics$focal_keep_frac, 0)
  best_frac <- dplyr::coalesce(current_best_metrics$focal_keep_frac, 0)
  cand_frac > best_frac
}

run_matching_with_focal_shrink <- function(prep, q_lo, q_hi, covariates, methods, focal_keep_fracs) {
  best_metrics <- NULL
  best_meta <- NULL
  attempts <- list()

  for (focal_keep_frac in focal_keep_fracs) {
    prep_subset <- shrink_prepared_focal(prep, focal_keep_frac)
    if (is.null(prep_subset)) {
      attempts[[length(attempts) + 1L]] <- tibble::tibble(
        focal_keep_frac = focal_keep_frac,
        focal_kept_n = NA_integer_,
        focal_total_n = NA_integer_,
        status = "skip_prep"
      )
      next
    }

    match_res <- run_matching_for_scenario(prep_subset, methods)
    if (is.null(match_res)) {
      attempts[[length(attempts) + 1L]] <- tibble::tibble(
        focal_keep_frac = prep_subset$focal_keep_frac,
        focal_kept_n = prep_subset$focal_kept_n,
        focal_total_n = prep_subset$focal_total_n,
        status = "skip_match"
      )
      next
    }

    metrics <- summarize_scenario(prep_subset, match_res, q_lo, q_hi, covariates)
    if (is.null(metrics)) {
      attempts[[length(attempts) + 1L]] <- tibble::tibble(
        focal_keep_frac = prep_subset$focal_keep_frac,
        focal_kept_n = prep_subset$focal_kept_n,
        focal_total_n = prep_subset$focal_total_n,
        status = "skip_metrics"
      )
      next
    }

    attempts[[length(attempts) + 1L]] <- tibble::tibble(
      focal_keep_frac = metrics$focal_keep_frac,
      focal_kept_n = metrics$focal_kept_n,
      focal_total_n = metrics$focal_total_n,
      status = "success",
      mean_abs_smd = metrics$mean_abs_smd,
      max_abs_smd = metrics$max_abs_smd,
      xpc_diff = metrics$xpc_diff
    )

    if (is_better_match_quality(metrics, best_metrics)) {
      best_metrics <- metrics
      best_meta <- list(prep = prep_subset, match_res = match_res)
    }
  }

  attempt_tbl <- if (length(attempts) > 0L) {
    dplyr::bind_rows(attempts)
  } else {
    tibble::tibble()
  }

  if (is.null(best_metrics)) {
    return(list(
      status = "no_success",
      attempts = attempt_tbl,
      best_metrics = NULL,
      best_meta = NULL
    ))
  }

  list(
    status = "success",
    attempts = attempt_tbl,
    best_metrics = best_metrics,
    best_meta = best_meta
  )
}

# -----------------------------
# Grid search over scenarios
# -----------------------------
scenario_results <- list()
scenario_meta <- list()
grid_total <- nrow(quantile_grid)
grid_success <- 0L
grid_skip_prep <- 0L
grid_skip_match <- 0L
grid_skip_metrics <- 0L
grid_error <- 0L
grid_processed <- 0L
grid_best_score <- -Inf
grid_best_idx <- NA_integer_
grid_start_time <- Sys.time()
grid_workers <- if (requested_workers == 0L) auto_workers else requested_workers
grid_workers <- max(1L, grid_workers)
if (grid_total > 0L) {
  grid_workers <- min(grid_workers, grid_total)
}

if (verbose_grid_search) {
  cat(
    "[grid] execution_mode=", ifelse(grid_workers > 1L, "parallel", "sequential"),
    " workers=", grid_workers,
    " requested=", requested_workers,
    " auto=", auto_workers, "\n",
    sep = ""
  )
}

evaluate_grid_scenario <- function(idx) {
  tryCatch(
    {
      q_lo <- quantile_grid$q_lo[idx]
      q_hi <- quantile_grid$q_hi[idx]

      prep <- prepare_matching_pool(analysis_df, q_lo, q_hi, covariate_names, critical_exact_match_cols)
      if (is.null(prep)) {
        return(list(
          scenario_index = idx,
          q_lo = q_lo,
          q_hi = q_hi,
          status = "skip_prep"
        ))
      }

      prep_stats <- list(
        xpc_lo = prep$q_vals[[1]],
        xpc_hi = prep$q_vals[[2]],
        rows = nrow(prep$data),
        treated = sum(prep$data$treatment == 1L, na.rm = TRUE),
        pool = sum(prep$data$treatment == 0L, na.rm = TRUE),
        strata = dplyr::n_distinct(prep$data$critical_key)
      )

      shrink_eval <- run_matching_with_focal_shrink(
        prep = prep,
        q_lo = q_lo,
        q_hi = q_hi,
        covariates = covariate_names,
        methods = matching_methods,
        focal_keep_fracs = focal_keep_frac_seq
      )
      if (shrink_eval$status != "success") {
        return(list(
          scenario_index = idx,
          q_lo = q_lo,
          q_hi = q_hi,
          status = "skip_match",
          prep_stats = prep_stats,
          shrink_attempts = shrink_eval$attempts
        ))
      }

      metrics <- shrink_eval$best_metrics
      metrics$scenario_index <- idx
      list(
        scenario_index = idx,
        q_lo = q_lo,
        q_hi = q_hi,
        status = "success",
        prep_stats = prep_stats,
        shrink_attempts = shrink_eval$attempts,
        metrics = metrics,
        meta = shrink_eval$best_meta
      )
    },
    error = function(e) {
      list(
        scenario_index = idx,
        q_lo = quantile_grid$q_lo[idx],
        q_hi = quantile_grid$q_hi[idx],
        status = "error",
        error_message = conditionMessage(e)
      )
    }
  )
}

process_grid_result <- function(res) {
  idx <- res$scenario_index
  on.exit(
    {
      grid_processed <<- grid_processed + 1L
      if (verbose_grid_search) {
        elapsed_sec <- as.numeric(difftime(Sys.time(), grid_start_time, units = "secs"))
        remaining <- max(0L, grid_total - grid_processed)
        eta_sec <- if (grid_processed > 0L && elapsed_sec > 0) {
          (elapsed_sec / grid_processed) * remaining
        } else {
          NA_real_
        }
        cat(
          "[grid] progress=", grid_processed, "/", grid_total,
          " elapsed=", format_duration(elapsed_sec),
          " eta=", format_duration(eta_sec),
          " status=", res$status,
          " idx=", idx, "\n",
          sep = ""
        )
      }
    },
    add = TRUE
  )
  if (verbose_grid_search) {
    cat("[grid][", idx, "/", grid_total, "] q_lo=", res$q_lo, " q_hi=", res$q_hi, "\n", sep = "")
    if (!is.null(res$prep_stats)) {
      cat(
        "[grid][", idx, "/", grid_total, "] ",
        "xpc_lo=", res$prep_stats$xpc_lo,
        " xpc_hi=", res$prep_stats$xpc_hi,
        " rows=", res$prep_stats$rows,
        " treated=", res$prep_stats$treated,
        " pool=", res$prep_stats$pool,
        " strata=", res$prep_stats$strata, "\n",
        sep = ""
      )
    }
    if (!is.null(res$shrink_attempts) && nrow(res$shrink_attempts) > 0L) {
      shrink_lines <- res$shrink_attempts %>%
        mutate(
          frac_label = ifelse(is.na(focal_keep_frac), "NA", sprintf("%.2f", focal_keep_frac)),
          kept_label = ifelse(is.na(focal_kept_n), "NA", as.character(focal_kept_n)),
          total_label = ifelse(is.na(focal_total_n), "NA", as.character(focal_total_n)),
          status_label = status
        ) %>%
        transmute(label = paste0(frac_label, "(", kept_label, "/", total_label, "):", status_label)) %>%
        pull(label)
      cat(
        "[grid][", idx, "/", grid_total, "] focal_attempts=",
        paste(shrink_lines, collapse = " | "), "\n",
        sep = ""
      )
    }
  }

  if (res$status == "skip_prep") {
    grid_skip_prep <<- grid_skip_prep + 1L
    if (verbose_grid_search) {
      cat("[grid][", idx, "/", grid_total, "] skip=prep(no valid treated/pool split)\n", sep = "")
    }
    return(invisible(NULL))
  }

  if (res$status == "skip_match") {
    grid_skip_match <<- grid_skip_match + 1L
    if (verbose_grid_search) {
      cat("[grid][", idx, "/", grid_total, "] skip=matching(no valid matched strata)\n", sep = "")
    }
    return(invisible(NULL))
  }

  if (res$status == "error") {
    grid_error <<- grid_error + 1L
    if (verbose_grid_search) {
      cat("[grid][", idx, "/", grid_total, "] skip=error(", res$error_message, ")\n", sep = "")
    }
    return(invisible(NULL))
  }

  metrics <- res$metrics
  grid_success <<- grid_success + 1L
  if (verbose_grid_search) {
    cat(
      "[grid][", idx, "/", grid_total, "] success ",
      "pairs=", metrics$matched_pairs,
      " focal_keep=", sprintf("%.2f", metrics$focal_keep_frac),
      " focal_n=", metrics$focal_kept_n, "/", metrics$focal_total_n,
      " method=", metrics$method,
      " replace=", metrics$replace,
      " xpc_diff=", metrics$xpc_diff,
      " max_abs_smd=", metrics$max_abs_smd,
      " mean_abs_smd=", metrics$mean_abs_smd,
      " score=", metrics$score, "\n",
      sep = ""
    )
  }
  if (!is.na(metrics$score) && (is.na(grid_best_idx) || metrics$score > grid_best_score)) {
    grid_best_score <<- metrics$score
    grid_best_idx <<- idx
    if (verbose_grid_search) {
      cat("[grid] best_update idx=", idx, " score=", metrics$score, "\n", sep = "")
    }
  }
  scenario_results[[length(scenario_results) + 1L]] <<- metrics
  scenario_meta[[length(scenario_meta) + 1L]] <<- res$meta
  invisible(NULL)
}

grid_indices <- seq_len(grid_total)
if (length(grid_indices) == 0L) {
  # no-op
} else if (grid_workers <= 1L) {
  if (verbose_grid_search) {
    cat("[grid] processing sequentially across ", length(grid_indices), " scenarios\n", sep = "")
  }
  for (idx in grid_indices) {
    process_grid_result(evaluate_grid_scenario(idx))
  }
} else {
  cl <- parallel::makePSOCKcluster(grid_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  parallel::clusterSetRNGStream(cl, iseed = 12345)
  invisible(parallel::clusterEvalQ(
    cl,
    suppressPackageStartupMessages({
      library(dplyr)
      library(readr)
      library(purrr)
      library(tibble)
      library(tidyr)
      library(forcats)
      library(stringr)
      library(GenomicRanges)
      library(nullranges)
    })
  ))
  parallel::clusterExport(
    cl,
    varlist = c(
      "quantile_grid",
      "analysis_df",
      "covariate_names",
      "critical_exact_match_cols",
      "matching_methods",
      "drop_failed_strata",
      "scale_covariate",
      "compute_smd_table",
      "prepare_matching_pool",
      "match_strata",
      "run_matching_for_scenario",
      "shrink_prepared_focal",
      "is_better_match_quality",
      "run_matching_with_focal_shrink",
      "summarize_scenario",
      "evaluate_grid_scenario",
      "focal_keep_frac_seq"
    ),
    envir = environment()
  )

  chunk_size <- grid_workers
  n_chunks <- as.integer(ceiling(length(grid_indices) / chunk_size))
  for (chunk_id in seq_len(n_chunks)) {
    start_pos <- (chunk_id - 1L) * chunk_size + 1L
    end_pos <- min(length(grid_indices), chunk_id * chunk_size)
    chunk_indices <- grid_indices[start_pos:end_pos]
    if (verbose_grid_search) {
      cat(
        "[grid] dispatch chunk=", chunk_id, "/", n_chunks,
        " size=", length(chunk_indices),
        " scenarios=", min(chunk_indices), "-", max(chunk_indices), "\n",
        sep = ""
      )
    }
    chunk_results <- parallel::parLapply(cl, chunk_indices, evaluate_grid_scenario)
    for (res in chunk_results) {
      process_grid_result(res)
    }
    if (verbose_grid_search) {
      elapsed_sec <- as.numeric(difftime(Sys.time(), grid_start_time, units = "secs"))
      remaining <- max(0L, grid_total - grid_processed)
      eta_sec <- if (grid_processed > 0L && elapsed_sec > 0) {
        (elapsed_sec / grid_processed) * remaining
      } else {
        NA_real_
      }
      cat(
        "[grid] complete chunk=", chunk_id, "/", n_chunks,
        " processed=", end_pos, "/", length(grid_indices),
        " elapsed=", format_duration(elapsed_sec),
        " eta=", format_duration(eta_sec), "\n",
        sep = ""
      )
    }
  }
  parallel::stopCluster(cl)
}

if (verbose_grid_search) {
  grid_elapsed_sec <- as.numeric(difftime(Sys.time(), grid_start_time, units = "secs"))
  cat(
    "[grid] complete total=", grid_total,
    " success=", grid_success,
    " skip_prep=", grid_skip_prep,
    " skip_matching=", grid_skip_match,
    " skip_metrics=", grid_skip_metrics,
    " errors=", grid_error,
    " elapsed=", format_duration(grid_elapsed_sec), "\n",
    sep = ""
  )
  if (!is.na(grid_best_idx)) {
    cat("[grid] best_seen idx=", grid_best_idx, " score=", grid_best_score, "\n", sep = "")
  }
}

if (length(scenario_results) == 0) {
  cat("No matching configuration succeeded for the supplied grid.\n")
  quit(save = "no", status = 0)
}

# -----------------------------
# Build summary tables and write CSVs
# -----------------------------
summary_tbl <- purrr::map_dfr(scenario_results, function(res) {
  tibble::tibble(
    q_lo = res$q_lo,
    q_hi = res$q_hi,
    xpc_lo_value = res$xpc_lo_value,
    xpc_hi_value = res$xpc_hi_value,
    focal_keep_frac = res$focal_keep_frac,
    focal_kept_n = res$focal_kept_n,
    focal_total_n = res$focal_total_n,
    method = res$method,
    replace = res$replace,
    matched_pairs = res$matched_pairs,
    total_pairs_raw = res$total_pairs_raw,
    unique_controls = res$unique_controls,
    total_controls = res$total_controls,
    duplicates = res$duplicates,
    xpc_treated_mean = res$xpc_treated_mean,
    xpc_control_mean = res$xpc_control_mean,
    xpc_diff = res$xpc_diff,
    max_abs_smd = res$max_abs_smd,
    mean_abs_smd = res$mean_abs_smd,
    score = res$score,
    smd_table = list(res$smd_table),
    scenario_index = res$scenario_index
  )
})

summary_tbl <- summary_tbl %>%
  mutate(
    score_order = dplyr::coalesce(score, -Inf),
    max_smd_order = dplyr::coalesce(max_abs_smd, Inf)
  ) %>%
  arrange(desc(score_order), desc(xpc_diff), max_smd_order) %>%
  mutate(scenario_rank = row_number()) %>%
  select(-score_order, -max_smd_order)

summary_flat <- summary_tbl %>%
  select(-smd_table)

smd_long <- purrr::map2_dfr(
  seq_len(nrow(summary_tbl)),
  summary_tbl$smd_table,
  function(idx, smd_tbl) {
    if (is.null(smd_tbl) || nrow(smd_tbl) == 0) {
      return(tibble::tibble(
        scenario_rank = summary_tbl$scenario_rank[[idx]],
        q_lo = summary_tbl$q_lo[[idx]],
        q_hi = summary_tbl$q_hi[[idx]],
        method = summary_tbl$method[[idx]],
        replace = summary_tbl$replace[[idx]],
        covariate = character(0),
        mean_treated = numeric(0),
        mean_control = numeric(0),
        smd = numeric(0)
      ))
    }
    smd_tbl %>%
      mutate(
        scenario_rank = summary_tbl$scenario_rank[[idx]],
        q_lo = summary_tbl$q_lo[[idx]],
        q_hi = summary_tbl$q_hi[[idx]],
        method = summary_tbl$method[[idx]],
        replace = summary_tbl$replace[[idx]]
      ) %>%
      select(
        scenario_rank, q_lo, q_hi, method, replace,
        covariate, mean_treated, mean_control, smd
      )
  }
)

readr::write_csv(summary_flat, summary_file)
readr::write_csv(smd_long, smd_file)

cat("Grid search summary written to:", summary_file, "\n")
cat("Covariate SMD details written to:", smd_file, "\n")

# -----------------------------
# Select best scenario and Export Data
# -----------------------------
best_row <- summary_tbl %>%
  mutate(
    mean_abs_smd_order = dplyr::coalesce(mean_abs_smd, Inf),
    max_abs_smd_order = dplyr::coalesce(max_abs_smd, Inf),
    xpc_diff_order = -dplyr::coalesce(xpc_diff, -Inf)
  ) %>%
  arrange(mean_abs_smd_order, max_abs_smd_order, xpc_diff_order) %>%
  dplyr::slice_head(n = 1)

best_idx <- best_row$scenario_index[[1]]
best_rank <- best_row$scenario_rank[[1]]

cat("\nBest Scenario Selected (Rank", best_rank, "):\n")
cat("  Index:", best_idx, "\n")
cat("  Focal keep:", best_row$focal_kept_n, "/", best_row$focal_total_n, " (", best_row$focal_keep_frac, ")\n", sep = "")
cat("  Mean |SMD|:", best_row$mean_abs_smd, "\n")
cat("  XPC Diff:", best_row$xpc_diff, "\n")

# Retrieve metadata for the best scenario
best_metrics <- scenario_results[[ which(vapply(scenario_results, function(x) x$scenario_index, integer(1)) == best_idx) ]]
best_meta    <- scenario_meta[[ which(vapply(scenario_results, function(x) x$scenario_index, integer(1)) == best_idx) ]]

prep_best    <- best_meta$prep
match_best   <- best_meta$match_res
data_best    <- prep_best$data

used_strata <- best_metrics$used_strata
treated_ids <- best_metrics$treated_ids
matched_control_ids <- best_metrics$matched_control_ids

# Extract raw data rows for various groups
# 1. Focal (Treated)
# 2. Pool (All available controls in used strata)
# 3. Matched Controls (Subset of pool)
# 4. Unmatched Controls (Subset of pool)

# Get Pool (Controls in used strata)
pool_in_used <- data_best %>%
  filter(treatment == 0L, critical_key %in% used_strata)

# Identify Unmatched
unmatched_controls <- pool_in_used %>%
  filter(!(row_id %in% matched_control_ids))

# Prepare focal data (with pair_ids)
focal_pairs_info <- match_best$matched_focal %>%
  select(row_id, pair_id) %>%
  distinct()

df_focal <- data_best %>%
  filter(row_id %in% treated_ids) %>%
  left_join(focal_pairs_info, by = "row_id") %>%
  mutate(group_label = "Focal", treatment = 1L)

# Prepare matched control data (with pair_ids)
control_pairs_info <- match_best$matched_control %>%
  select(row_id, pair_id) %>%
  distinct()

df_matched <- data_best %>%
  filter(row_id %in% matched_control_ids) %>%
  left_join(control_pairs_info, by = "row_id") %>%
  mutate(group_label = "Matched Control", treatment = 0L)

# Prepare unmatched control data (no pairs)
df_unmatched <- unmatched_controls %>%
  mutate(pair_id = NA_integer_, group_label = "Unmatched Control", treatment = 0L)

# Combine all into one dataframe
final_export_df <- bind_rows(df_focal, df_matched, df_unmatched) %>%
  select(
    chrom, start, end,
    row_id, pair_id, group_label, treatment,
    xpc_signal, ins_eff, n_ctcf_peaks, looped_bin, is_boundary_bin,
    everything()
  ) %>%
  # drop internal scaling columns if desired, or keep them.
  # We will keep everything for completeness, but sorting helps.
  arrange(group_label, pair_id)

# Write to file
readr::write_csv(final_export_df, output_data_file)

cat("\nBest scenario data written to:", output_data_file, "\n")
cat("This file contains the following groups in the 'group_label' column:\n")
cat(" - Focal (High XPC)\n")
cat(" - Matched Control\n")
cat(" - Unmatched Control\n")
cat("Note: 'Pool Control' is the combination of 'Matched Control' and 'Unmatched Control'.\n")

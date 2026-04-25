import os
import cooler
import pandas as pd
import bioframe
import numpy as np
from cooltools.sandbox import obs_over_exp_cooler

MCOOLS_PATH_DICT = {
    # --- WT ---
    "WTnoUV": "/home/carlos/oldies/ner_collab/mcools_merged_biorep/WTnoUV.mcool",
    "WT1h":   "/home/carlos/oldies/ner_collab/mcools_merged/WT_1H_1.mcool",
    "WT3h":   "/home/carlos/oldies/ner_collab/mcools_merged_biorep/WT3h.mcool",
    "WT6h":   "/home/carlos/oldies/ner_collab/mcools_merged/WT_6H_6.mcool",


    # --- XPA ---
    "XPAnoUV": "/home/carlos/oldies/ner_collab/mcools_merged_biorep/XPAnoUV.mcool",
    "XPA1h":   "/home/carlos/oldies/ner_collab/mcools_merged/XPA_1H_1.mcool",
    "XPA3h":   "/home/carlos/oldies/ner_collab/mcools_merged_biorep/XPA3h.mcool",
    "XPA6h":   "/home/carlos/oldies/ner_collab/mcools_merged/XPA_6H_1.mcool",

    # --- XPC ---
    "XPCnoUV": "/home/carlos/oldies/ner_collab/mcools_merged_biorep/XPCnoUV.mcool",
    "XPC3h":   "/home/carlos/oldies/ner_collab/mcools_merged_biorep/XPC3h.mcool",

    # # --- Pol2 ---
    # "POL2_control": "/home/carlos/Clone/ggner-3d/data/degron_pol2/GSM5394172_control_mergeRep_500bp.mcool",
    # "POL2_degron": "/home/carlos/Clone/ggner-3d/data/degron_pol2/GSM5394173_degron_mergeRep_500bp.mcool",

    # # --- Wapl ---
    # "WAPL_control": "/home/carlos/Clone/ggner-3d/data/degron_wapl/GSE178982_WAPL-UT_pool.mcool",
    # "WAPL_degron": "/home/carlos/Clone/ggner-3d/data/degron_wapl/GSE178982_WAPL-AID_pool.mcool",
}

def CLR_CONNS(resolution, genotype=None, timepoint=None, contains=None, biomerged=True):
    """
    Return dict of Cooler objects filtered by genotype, timepoint, or custom substrings.

    Parameters
    ----------
    resolution : int
        Resolution (bp) for selecting the cooler matrix.
    genotype : str or list or None
        e.g. "WT", ["WT","XPA"]
    timepoint : str or list or None
        e.g. "noUV", "1h", ["1h","3h"]
    contains : str or list or None
        Arbitrary substring filters. Example: contains="noUV"

    Returns
    -------
    dict
        e.g. {"WT1h": Cooler(...), ...}
    """
    def to_list(x):
        if x is None:
            return None
        return x if isinstance(x, (list, tuple)) else [x]

    genotype = to_list(genotype)
    timepoint = to_list(timepoint)
    contains = to_list(contains)

    filtered = {}

    for key, path in (MCOOLS_PATH_DICT if biomerged else MCOOL_PATH_DICT_BIOREPS).items():

        # Filter by genotype
        if genotype and not any(g in key for g in genotype):
            continue

        # Filter by timepoint
        if timepoint and not any(tp in key for tp in timepoint):
            continue

        # Arbitrary substring filters
        if contains and not all(sub in key for sub in contains):
            continue

        filtered[key] = cooler.Cooler(f"{path}::/resolutions/{resolution}")

    return filtered

MCOOL_PATH_DICT_BIOREPS = {
    # --- WT ---
    "WT_noUV_r1": "/home/carlos/oldies/ner_collab/mcools_merged/WTnoUV1.mcool",
    "WT_noUV_r2": "/home/carlos/oldies/ner_collab/mcools_merged/WTnoUV2.mcool",

    "WT_1h_r1": "/home/carlos/oldies/ner_collab/mcools_merged/WT_1H_1.mcool",

    "WT_3h_r1": "/home/carlos/oldies/ner_collab/mcools_merged/WT3h1.mcool",
    # "WT3h2": "/home/carlos/oldies/ner_collab/mcools_merged/WT3h2.mcool",
    "WT_3h_r2": "/home/carlos/oldies/ner_collab/mcools_merged/WT_3H_3.mcool",

    "WT_6h_r1": "/home/carlos/oldies/ner_collab/mcools_merged/WT_6H_6.mcool",

    # --- XPA ---
    "XPA_noUV_r1": "/home/carlos/oldies/ner_collab/mcools_merged/XPAnoUV1.mcool",
    "XPA_noUV_r2": "/home/carlos/oldies/ner_collab/mcools_merged/XPAnoUV2.mcool",

    "XPA_1h_r1": "/home/carlos/oldies/ner_collab/mcools_merged/XPA_1H_1.mcool",

    "XPA_3h_r1": "/home/carlos/oldies/ner_collab/mcools_merged/XPA3h1.mcool",
    "XPA_3h_r2": "/home/carlos/oldies/ner_collab/mcools_merged/XPA3h2.mcool",

    "XPA_6h_r1": "/home/carlos/oldies/ner_collab/mcools_merged/XPA_6H_1.mcool",

    # --- XPC ---
    "XPC_noUV_r1": "/home/carlos/oldies/ner_collab/mcools_merged/XPCnoUV1.mcool",
    "XPC_noUV_r2": "/home/carlos/oldies/ner_collab/mcools_merged/XPCnoUV2.mcool",

    "XPC_3h_r1": "/home/carlos/oldies/ner_collab/mcools_merged/XPC3h1.mcool",
    "XPC_3h_r2": "/home/carlos/oldies/ner_collab/mcools_merged/XPC3h2.mcool",
}

def assert_all_cools_has_weights(clr_dict, weight_name):
    for k, v in clr_dict.items():
        bins_table = v.bins()[:]
        if weight_name not in bins_table.columns:
            raise ValueError(f"Cooler {k} does not have weight column '{weight_name}'")
    print(f"All coolers have the required '{weight_name}' column; {[k for k in clr_dict.keys()]}")

def prepare_view_df(arm=True):
    # pick one clr
    clr = cooler.Cooler(list(MCOOLS_PATH_DICT.values())[0] + "::/resolutions/10000")
    if arm:
        hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
        hg38_cens = bioframe.fetch_centromeres('hg38')
        hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)
        view_df = hg38_arms.set_index("chrom").loc[clr.chromnames].reset_index()
    else:
        view_df = pd.DataFrame({'chrom': clr.chromnames,
                        'start': 0,
                        'end': clr.chromsizes.values,
                        'name': clr.chromnames}
                      )
    return view_df

def get_eig_data(sample_name, resolution):
    # return compartment eigenvalues and eigenvectors dataframes
    compartment_data_dir = '/home/carlos/Clone/ggner-3d/data/compartment'
    return pd.read_csv(
        f'{compartment_data_dir}/compartment_eigvals_{sample_name}_{resolution}.tsv', sep='\t'
    ), pd.read_csv(f'{compartment_data_dir}/compartment_eigvecs_{sample_name}_{resolution}.tsv', sep='\t')

EXPECTED_DATA_DIR = '/home/carlos/Clone/ggner-3d/data/expected'
def get_expected_cis_df(condition, resolution_kb, label='arm'):
    # label is either 'arm' or 'chrom'
    # assert file exists
    filename = f"expected_cis_{label}_{condition}_res{resolution_kb}kb.tsv"
    path = f"{EXPECTED_DATA_DIR}/{filename}"
    if not os.path.exists(path):
        raise FileNotFoundError(f"Expected cis file not found: {path}")
    return pd.read_csv(path, sep='\t')

def get_expected_trans_df(condition, resolution_kb):
    # assert file exists
    filename = f"expected_trans_{condition}_res{resolution_kb}kb.tsv"
    path = f"{EXPECTED_DATA_DIR}/{filename}"
    if not os.path.exists(path):
        raise FileNotFoundError(f"Expected trans file not found: {path}")
    return pd.read_csv(path, sep='\t')

DOTS_DATA_DIR = '/home/carlos/Clone/ggner-3d/data/dots'
def get_dots_df(condition, resolution):
    # assert file exists
    filename = f"dots_{condition}_{resolution}.csv"
    path = f"{DOTS_DATA_DIR}/{filename}"
    if not os.path.exists(path):
        raise FileNotFoundError(f"Dots file not found: {path}")
    return pd.read_csv(path)

def chrom_to_bin_indices(df, bins_table, bin_column = 'idx1'):
    df = df.copy()
    bins_table['indices'] = bins_table.index
    return pd.merge(df, bins_table, left_on=bin_column, right_on='indices')

def scaled_balanced_cool_dump_streaming(
    clr,
    *,
    target_sum=1e6,
    dump_path=None,
    field_name="scaled",
    cis_only=None,
    include_bad_bins_as_nan=True,
    weight_name="weight",
    divisive_weights=None,
    chunksize=1_000_000,
    return_pixels=False,
    count_dtype=np.int32,
):
    """
    Stream a scaled balanced value column from an input cooler.

    For each pixel:
      - balanced = count * w1 * w2            (multiplicative weights)
      - balanced = count / (w1 * w2)          (divisive weights)

    Scaling behavior depends on `cis_only`:

      - cis_only is None or False:
          Use all cis + trans pixels to compute one global scale factor.
          Output all pixels.

      - cis_only is True or "full":
          Use all cis pixels to compute one global scale factor.
          Output cis pixels only.

      - cis_only == "chrom":
          Use cis pixels within each chromosome to compute one scale factor
          per chromosome.
          Output cis pixels only.

    Parameters
    ----------
    clr : cooler.Cooler
        Input cooler.
    target_sum : float, optional
        Desired sum after scaling.
        For `cis_only="chrom"`, this is the desired sum per chromosome.
    dump_path : str or None, optional
        If provided, write a new cooler to this URI/path.
    field_name : str, optional
        Name of the output value column to write.
    cis_only : {None, False, True, "full", "chrom"}, optional
        Controls normalization scope and output selection.

        - None / False: normalize on all pixels, output all pixels
        - True / "full": normalize on all cis pixels, output cis only
        - "chrom": normalize per chromosome on cis pixels, output cis only
    include_bad_bins_as_nan : bool, optional
        Preserve NaN weights in output bins table when weight_name exists.
    weight_name : str, optional
        Bin weight column name.
    divisive_weights : bool or None, optional
        If None, mimic cooler.matrix() defaults:
        KR / VC / SQRT_VC are treated as divisive, others multiplicative.
    chunksize : int, optional
        Number of pixels per chunk.
    return_pixels : bool, optional
        If True, return (pixels_iter(), scale) instead of (None, scale).
        `scale` is a float for global modes and a dict for `cis_only="chrom"`.
    count_dtype : numpy dtype, optional
        Dtype to use for the required output 'count' column.

    Returns
    -------
    (generator or None, float or dict)
        (pixels iterator, scale) if return_pixels=True, else (None, scale).

        - float for global scaling modes
        - dict {chrom_name: scale} for `cis_only="chrom"`
    """
    import numpy as np
    import pandas as pd
    import cooler

    bins = clr.bins()[:]

    # Normalize/validate mode while preserving backward compatibility.
    if cis_only is True:
        cis_mode = "full"
    elif cis_only in (None, False):
        cis_mode = None
    elif cis_only in {"full", "chrom"}:
        cis_mode = cis_only
    else:
        raise ValueError(
            "cis_only must be one of: None, False, True, 'full', 'chrom'"
        )

    # Load weights.
    if weight_name in bins.columns:
        w = bins[weight_name].to_numpy(dtype=np.float64, copy=False)
    else:
        w = np.ones(len(bins), dtype=np.float64)

    # Match cooler.matrix() behavior by default.
    if divisive_weights is None:
        divisive = weight_name in {"KR", "VC", "SQRT_VC"}
    else:
        divisive = bool(divisive_weights)

    # Fast cis mask via chromosome codes.
    chrom_codes = pd.Categorical(bins["chrom"]).codes.astype(np.int32, copy=False)
    chrom_names = pd.Index(pd.Categorical(bins["chrom"]).categories).tolist()
    n_chroms = len(chrom_names)

    # Number of nonzero pixels.
    nnz = int(clr.info.get("nnz", 0)) or len(clr.pixels(join=False))
    pixsel = clr.pixels(join=False)

    def _balanced_from_chunk(df):
        b1 = df["bin1_id"].to_numpy(dtype=np.int64, copy=False)
        b2 = df["bin2_id"].to_numpy(dtype=np.int64, copy=False)
        cnt = df["count"].to_numpy(dtype=np.float64, copy=False)

        pair_w = w[b1] * w[b2]

        if divisive:
            bal = np.empty_like(cnt, dtype=np.float64)
            ok = np.isfinite(pair_w) & (pair_w != 0.0)
            np.divide(cnt, pair_w, out=bal, where=ok)
            bal[~ok] = np.nan
        else:
            bal = cnt * pair_w

        cis_mask = chrom_codes[b1] == chrom_codes[b2]
        return b1, b2, bal, cis_mask

    # Pass 1: compute scale(s).
    if cis_mode == "chrom":
        totals = np.zeros(n_chroms, dtype=np.float64)
        seen = np.zeros(n_chroms, dtype=bool)

        for lo in range(0, nnz, chunksize):
            hi = min(lo + chunksize, nnz)
            df = pixsel[lo:hi]
            b1, _, bal, cis_mask = _balanced_from_chunk(df)

            mask = cis_mask & np.isfinite(bal)
            if not np.any(mask):
                continue

            codes = chrom_codes[b1[mask]]
            totals += np.bincount(codes, weights=bal[mask], minlength=n_chroms)
            seen[codes] = True

        bad = seen & (~np.isfinite(totals) | (totals <= 0))
        if np.any(bad):
            bad_chroms = [chrom_names[i] for i in np.flatnonzero(bad)]
            raise ValueError(
                f"Per-chromosome balanced sums are invalid for: {bad_chroms}"
            )

        scale_arr = np.full(n_chroms, np.nan, dtype=np.float64)
        scale_arr[seen] = float(target_sum) / totals[seen]
        scale = {chrom_names[i]: float(scale_arr[i]) for i in np.flatnonzero(seen)}

    else:
        total = 0.0

        for lo in range(0, nnz, chunksize):
            hi = min(lo + chunksize, nnz)
            df = pixsel[lo:hi]
            _, _, bal, cis_mask = _balanced_from_chunk(df)

            if cis_mode == "full":
                total += float(np.nansum(bal[cis_mask]))
            else:
                total += float(np.nansum(bal))

        if not np.isfinite(total) or total <= 0:
            raise ValueError(f"Total balanced sum is {total}. Check weights / selection.")

        scale = float(target_sum) / total

    # Pass 2: stream scaled pixels.
    def pixels_iter():
        for lo in range(0, nnz, chunksize):
            hi = min(lo + chunksize, nnz)
            df = pixsel[lo:hi]
            b1, b2, bal, cis_mask = _balanced_from_chunk(df)

            if cis_mode is None:
                out_b1 = b1
                out_b2 = b2
                out_val = np.asarray(bal * scale, dtype=np.float64)

            elif cis_mode == "full":
                out_b1 = b1[cis_mask]
                out_b2 = b2[cis_mask]
                out_val = np.asarray(bal[cis_mask] * scale, dtype=np.float64)

            else:  # cis_mode == "chrom"
                out_b1 = b1[cis_mask]
                out_b2 = b2[cis_mask]
                out_codes = chrom_codes[out_b1]
                out_val = np.asarray(bal[cis_mask] * scale_arr[out_codes], dtype=np.float64)

            yield pd.DataFrame(
                {
                    "bin1_id": out_b1,
                    "bin2_id": out_b2,
                    "count": np.zeros(len(out_b1), dtype=count_dtype),
                    field_name: out_val,
                }
            )

    if dump_path is not None:
        bins_out = bins.copy()

        if include_bad_bins_as_nan and weight_name in bins_out.columns:
            bad = bins_out[weight_name].isna()
            bins_out[weight_name] = 1.0
            bins_out.loc[bad, weight_name] = np.nan

        cooler.create_cooler(
            cool_uri=dump_path,
            bins=bins_out,
            pixels=pixels_iter(),
            columns=[field_name],
            dtypes={"count": count_dtype, field_name: np.float64},
            ordered=True,
            symmetric_upper=(clr.storage_mode == "symmetric-upper"),
        )

    if return_pixels:
        return pixels_iter(), scale

    return None, scale
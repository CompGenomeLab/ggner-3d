import numpy as np
import bioframe
import cooltools
from cooltools.lib.common import align_track_with_cooler

import multiprocessing as mp
from itertools import combinations

import numpy as np
import pandas as pd
import bioframe
import cooler
from cooltools.api import saddle as _saddle_mod


def _saddle_worker(args):
    (
        clr_uri,
        expected,
        view_df,
        supports,
        digitized_tracks,
        contact_type,
        n_bins,
        clr_weight_name,
        expected_value_col,
        view_name_col,
        min_diag,
        max_diag,
        verbose,
    ) = args

    clr = cooler.Cooler(clr_uri)

    if contact_type == "cis":
        getmatrix = _saddle_mod._make_cis_obsexp_fetcher(
            clr,
            expected,
            view_df,
            clr_weight_name=clr_weight_name,
            expected_value_col=expected_value_col,
            view_name_col=view_name_col,
        )
    elif contact_type == "trans":
        getmatrix = _saddle_mod._make_trans_obsexp_fetcher(
            clr,
            expected,
            view_df,
            clr_weight_name=clr_weight_name,
            expected_value_col=expected_value_col,
            view_name_col=view_name_col,
        )
    else:
        raise ValueError("Allowed values for contact_type are 'cis' or 'trans'.")

    S = np.zeros((n_bins + 2, n_bins + 2), dtype=float)
    C = np.zeros((n_bins + 2, n_bins + 2), dtype=float)

    for reg1, reg2 in supports:
        _saddle_mod._accumulate(
            S,
            C,
            getmatrix,
            digitized_tracks,
            reg1,
            reg2,
            min_diag=min_diag,
            max_diag=max_diag,
            verbose=verbose,
        )

    return S, C


def saddle_mp(
    clr,
    expected,
    track,
    contact_type,
    n_bins,
    vrange=None,
    qrange=None,
    view_df=None,
    clr_weight_name="weight",
    expected_value_col="balanced.avg",
    view_name_col="name",
    min_diag=3,
    max_diag=-1,
    trim_outliers=False,
    verbose=False,
    drop_track_na=False,
    nproc=1,
):
    """
    Multiprocessing version of cooltools.api.saddle.saddle with the same args
    plus `nproc`.

    Returns
    -------
    interaction_sum : 2D array
    interaction_count : 2D array
    """
    if not isinstance(nproc, int) or nproc < 1:
        raise ValueError("nproc must be a positive integer")

    # Keep exact upstream behavior when nproc == 1
    if nproc == 1:
        return _saddle_mod.saddle(
            clr=clr,
            expected=expected,
            track=track,
            contact_type=contact_type,
            n_bins=n_bins,
            vrange=vrange,
            qrange=qrange,
            view_df=view_df,
            clr_weight_name=clr_weight_name,
            expected_value_col=expected_value_col,
            view_name_col=view_name_col,
            min_diag=min_diag,
            max_diag=max_diag,
            trim_outliers=trim_outliers,
            verbose=verbose,
            drop_track_na=drop_track_na,
        )

    # ---- copied closely from upstream saddle() setup ----
    if isinstance(n_bins, int):
        track = _saddle_mod.align_track_with_cooler(
            track,
            clr,
            view_df=view_df,
            clr_weight_name=clr_weight_name,
            mask_clr_bad_bins=True,
            drop_track_na=drop_track_na,
        )
        digitized_track, _ = _saddle_mod.digitize(
            track.iloc[:, :4],
            n_bins,
            vrange=vrange,
            qrange=qrange,
            digitized_suffix=".d",
        )
        digitized_col = digitized_track.columns[3]
    elif n_bins is None:
        digitized_track = track
        digitized_col = digitized_track.columns[3]
        _saddle_mod.is_track(track.astype({digitized_col: "float"}), raise_errors=True)
        if not isinstance(digitized_track.dtypes.iloc[3], pd.CategoricalDtype):
            raise ValueError(
                "when n_bins=None, saddle assumes the track has been "
                "pre-digitized and the value column is a pandas categorical. "
                "See cooltools.api.saddle.digitize()."
            )
        cats = digitized_track[digitized_col].dtype.categories.values
        n_bins = len(cats[cats > -1]) - 2
    else:
        raise ValueError("n_bins must be provided as int or None")

    if view_df is None:
        view_df = _saddle_mod.view_from_track(digitized_track)
    else:
        try:
            _saddle_mod.is_compatible_viewframe(
                view_df,
                clr,
                check_sorting=True,
                raise_errors=True,
            )
        except Exception as e:
            raise ValueError("view_df is not a valid viewframe or incompatible") from e

    try:
        _saddle_mod.is_valid_expected(
            expected,
            contact_type,
            view_df,
            verify_cooler=clr,
            expected_value_cols=[expected_value_col],
            raise_errors=True,
        )
    except Exception as e:
        raise ValueError("provided expected is not compatible") from e

    if clr_weight_name:
        try:
            _saddle_mod.is_cooler_balanced(clr, clr_weight_name, raise_errors=True)
        except Exception as e:
            raise ValueError(
                f"provided cooler is not balanced or {clr_weight_name} is missing"
            ) from e

    digitized_tracks = {}
    for _, reg in view_df.iterrows():
        digitized_reg = bioframe.select(digitized_track, reg)
        digitized_tracks[reg[view_name_col]] = digitized_reg[digitized_col]

    if contact_type == "cis":
        supports = list(zip(view_df[view_name_col], view_df[view_name_col]))
    elif contact_type == "trans":
        chrom_by_name = view_df.set_index(view_name_col)["chrom"]
        supports = [
            pair
            for pair in combinations(view_df[view_name_col], 2)
            if chrom_by_name[pair[0]] != chrom_by_name[pair[1]]
        ]
    else:
        raise ValueError("Allowed values for contact_type are 'cis' or 'trans'.")

    interaction_sum = np.zeros((n_bins + 2, n_bins + 2), dtype=float)
    interaction_count = np.zeros((n_bins + 2, n_bins + 2), dtype=float)

    if len(supports) == 0:
        interaction_sum += interaction_sum.T
        interaction_count += interaction_count.T
        if trim_outliers:
            interaction_sum = interaction_sum[1:-1, 1:-1]
            interaction_count = interaction_count[1:-1, 1:-1]
        return interaction_sum, interaction_count

    n_jobs = min(nproc, len(supports))
    support_chunks = [supports[i::n_jobs] for i in range(n_jobs)]

    worker_args = [
        (
            clr.uri,
            expected,
            view_df,
            chunk,
            digitized_tracks,
            contact_type,
            n_bins,
            clr_weight_name,
            expected_value_col,
            view_name_col,
            min_diag,
            max_diag,
            verbose,
        )
        for chunk in support_chunks
    ]

    # h5py recommends independent open-per-process readers; use spawn.
    ctx = mp.get_context("spawn")
    with ctx.Pool(n_jobs) as pool:
        parts = pool.map(_saddle_worker, worker_args)

    for S_part, C_part in parts:
        interaction_sum += S_part
        interaction_count += C_part

    interaction_sum += interaction_sum.T
    interaction_count += interaction_count.T

    if trim_outliers:
        interaction_sum = interaction_sum[1:-1, 1:-1]
        interaction_count = interaction_count[1:-1, 1:-1]

    return interaction_sum, interaction_count
    
def pq_saddle(
    clr,
    expected_pq,
    track,
    view_df_arm,
    clr_weight_name="sweight",
    expected_value_col="balanced.avg.smoothed.agg",
    n_bins=25,
    qrange=(0.02, 0.98),
    view_name_col="name",
    drop_track_na=False,
    symmetrize=True,
    trim_outliers=False,
):
    binsize = clr.binsize

    # cooltools wants a 4-col bedGraph-like track: chrom,start,end,value
    track4 = track.iloc[:, :4].copy()
    if track4.columns[3] != "value":
        track4 = track4.rename(columns={track4.columns[3]: "value"})

    # 1) align track to cooler bins
    track_aligned = align_track_with_cooler(
        track4,
        clr,
        view_df=view_df_arm,
        clr_weight_name=clr_weight_name,
        mask_clr_bad_bins=True,
        drop_track_na=drop_track_na,
    )

    # 2) digitize exactly like cooltools
    digitized_track, binedges = cooltools.digitize(
        track_aligned[["chrom", "start", "end", "value"]],
        n_bins=n_bins,
        qrange=qrange,
    )
    dcol = digitized_track.columns[3]

    # 3) store digitized values per arm
    arm_info = view_df_arm.set_index(view_name_col)[["chrom", "start", "end"]]
    digitized = {}
    for _, reg in view_df_arm.iterrows():
        name = reg[view_name_col]
        # keep integer labels: -1, 0, 1..n, n+1
        digitized[name] = bioframe.select(digitized_track, reg)[dcol].astype(int).to_numpy()

    # 4) expected lookup
    exp_lookup = {
        (r1, r2): df.set_index("dist")[expected_value_col]
        for (r1, r2), df in expected_pq.groupby(["region1", "region2"])
    }

    # n regular bins + 2 outlier bins (0 and n+1), matching cooltools
    nb = n_bins + 2
    S = np.zeros((nb, nb), dtype=float)
    C = np.zeros((nb, nb), dtype=float)

    for chrom in view_df_arm["chrom"].unique():
        p = f"{chrom}_p"
        q = f"{chrom}_q"

        if p not in digitized or q not in digitized:
            continue

        if (p, q) in exp_lookup:
            exp_series = exp_lookup[(p, q)]
            row_arm, col_arm = p, q
        elif (q, p) in exp_lookup:
            exp_series = exp_lookup[(q, p)]
            row_arm, col_arm = q, p
        else:
            continue

        reg1 = tuple(arm_info.loc[row_arm])
        reg2 = tuple(arm_info.loc[col_arm])

        obs = clr.matrix(balance=clr_weight_name).fetch(reg1, reg2)

        start1_bins = arm_info.loc[row_arm, "start"] // binsize
        start2_bins = arm_info.loc[col_arm, "start"] // binsize

        rr = np.arange(obs.shape[0])[:, None]
        cc = np.arange(obs.shape[1])[None, :]
        dist_mat = (start2_bins - start1_bins) + cc - rr

        exp_mat = exp_series.reindex(dist_mat.ravel()).to_numpy().reshape(obs.shape)
        oe = obs / exp_mat

        rcat = digitized[row_arm]
        ccat = digitized[col_arm]

        # Drop only missing-data category (-1), keep outliers 0 and n+1
        keep_r = rcat >= 0
        keep_c = ccat >= 0
        oe = oe[np.ix_(keep_r, keep_c)]
        rcat = rcat[keep_r]
        ccat = ccat[keep_c]

        pair_codes = rcat[:, None] * nb + ccat[None, :]
        good = np.isfinite(oe)

        S += np.bincount(
            pair_codes[good].ravel(),
            weights=oe[good].ravel(),
            minlength=nb * nb,
        ).reshape(nb, nb)

        C += np.bincount(
            pair_codes[good].ravel(),
            minlength=nb * nb,
        ).reshape(nb, nb)

    # match cooltools: symmetrize first
    if symmetrize:
        S = S + S.T
        C = C + C.T

    # match cooltools: trim outlier bins last
    if trim_outliers:
        S = S[1:-1, 1:-1]
        C = C[1:-1, 1:-1]

    M = S / C
    M[C == 0] = np.nan

    return {
        "saddle": M,
        "sum": S,
        "count": C,
        "binedges": binedges,
        "digitized_track": digitized_track,
        "track_aligned": track_aligned,
    }
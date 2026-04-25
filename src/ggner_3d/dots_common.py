import numpy as np
import pandas as pd
import cooler
import cooltools

def compare_loops_with_strength(
    df_a,
    df_b,
    clr_a=None,
    clr_b=None,
    expected_a=None,
    expected_b=None,
    view_df_a=None,
    view_df_b=None,
    resolution=None,
    clustering_radius=20_000,
    flank=100_000,
    clr_weight_name="weight",
    nproc=1,
    sample_a_name="sample_a",
    sample_b_name="sample_b",
    deduplicate=True,
    compute_strength=True,
):
    req = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    miss_a = [c for c in req if c not in df_a.columns]
    miss_b = [c for c in req if c not in df_b.columns]
    if miss_a or miss_b:
        raise ValueError(f"Missing required columns. {sample_a_name}: {miss_a}; {sample_b_name}: {miss_b}")

    if compute_strength:
        if clr_a is None or clr_b is None:
            raise ValueError("clr_a and clr_b are required when compute_strength=True")
        if expected_a is None:
            raise ValueError("expected_a is required when compute_strength=True")
        if view_df_a is None:
            raise ValueError("view_df_a is required when compute_strength=True")
        if expected_b is None:
            expected_b = expected_a
        if view_df_b is None:
            view_df_b = view_df_a

    def _to_cooler(clr_like):
        if isinstance(clr_like, cooler.Cooler):
            return clr_like
        uri = str(clr_like)
        if "::" in uri:
            return cooler.Cooler(uri)
        if resolution is None:
            raise ValueError("resolution is required when clr path is .mcool without ::resolutions/<res>")
        return cooler.Cooler(f"{uri}::resolutions/{resolution}")

    def _prepare(df, id_col):
        x = df[req].copy()
        if deduplicate:
            x = x.drop_duplicates()
        x = x.reset_index(drop=True)
        x[id_col] = np.arange(len(x), dtype=np.int64)
        x["_bin1"] = (x["start1"] // clustering_radius).astype(np.int64)
        x["_bin2"] = (x["start2"] // clustering_radius).astype(np.int64)
        return x

    # --- NEW: anchor extraction + matching helpers ---------------------------------
    def _prepare_anchors(loops_df, id_col):
        # raw anchors (duplicates allowed so we can compute n_loops)
        a1 = loops_df[["chrom1", "start1", "end1"]].copy()
        a1.columns = ["chrom", "start", "end"]
        a2 = loops_df[["chrom2", "start2", "end2"]].copy()
        a2.columns = ["chrom", "start", "end"]
        raw = pd.concat([a1, a2], ignore_index=True)

        if deduplicate:
            # anchor degree (# loops touching anchor)
            anchors = (
                raw.groupby(["chrom", "start", "end"], as_index=False)
                   .size()
                   .rename(columns={"size": "n_loops"})
            )
        else:
            anchors = raw.copy()
            anchors["n_loops"] = 1

        anchors = anchors.reset_index(drop=True)
        anchors[id_col] = np.arange(len(anchors), dtype=np.int64)

        # use midpoint for binning/matching
        anchors["_pos"] = ((anchors["start"].astype(np.int64) + anchors["end"].astype(np.int64)) // 2).astype(np.int64)
        anchors["_abin"] = (anchors["_pos"] // clustering_radius).astype(np.int64)
        return anchors

    def _match_anchors(anc_a, anc_b, id_a, id_b):
        a_small = anc_a[["chrom", "_abin", id_a, "_pos"]]
        b_small = anc_b[["chrom", "_abin", id_b, "_pos"]]
        keys = ["chrom", "_abin"]

        candidates = []
        for d in (-1, 0, 1):
            b_shift = b_small.copy()
            b_shift["_abin"] += d
            m = a_small.merge(b_shift, on=keys, how="inner", suffixes=("_a", "_b"))
            if m.empty:
                continue
            dpos = (m["_pos_a"] - m["_pos_b"]).abs()
            keep = dpos <= clustering_radius
            if not keep.any():
                continue
            mm = m.loc[keep, [f"{id_a}", f"{id_b}"]].copy()
            mm["delta_pos"] = dpos[keep].to_numpy()
            mm["dist_bp"] = mm["delta_pos"]
            candidates.append(mm)

        if not candidates:
            return (
                pd.DataFrame(columns=[id_a, id_b, "delta_pos", "dist_bp"]),
                anc_a.copy(),
                anc_b.copy(),
            )

        cands = pd.concat(candidates, ignore_index=True).drop_duplicates([id_a, id_b])
        cands = cands.sort_values(["dist_bp", "delta_pos", id_a, id_b])

        used_a, used_b, keep_rows = set(), set(), []
        for _ida, _idb, dpos, dist in cands[[id_a, id_b, "delta_pos", "dist_bp"]].itertuples(index=False, name=None):
            if _ida in used_a or _idb in used_b:
                continue
            used_a.add(_ida)
            used_b.add(_idb)
            keep_rows.append((_ida, _idb, dpos, dist))

        matched = pd.DataFrame(keep_rows, columns=[id_a, id_b, "delta_pos", "dist_bp"])

        matched_a = set(matched[id_a].tolist())
        matched_b = set(matched[id_b].tolist())

        a_spec = anc_a.loc[~anc_a[id_a].isin(matched_a), ["chrom", "start", "end", "n_loops"]].reset_index(drop=True)
        b_spec = anc_b.loc[~anc_b[id_b].isin(matched_b), ["chrom", "start", "end", "n_loops"]].reset_index(drop=True)

        common = (
            matched
            .merge(anc_a[[id_a, "chrom", "start", "end", "n_loops"]], on=id_a, how="left")
            .merge(
                anc_b[[id_b, "chrom", "start", "end", "n_loops"]],
                on=id_b,
                how="left",
                suffixes=(f"_{sample_a_name}", f"_{sample_b_name}"),
            )
        )
        # rename the left-side columns to carry sample suffix too
        common = common.rename(
            columns={
                "chrom": f"chrom_{sample_a_name}",
                "start": f"start_{sample_a_name}",
                "end": f"end_{sample_a_name}",
                "n_loops": f"n_loops_{sample_a_name}",
            }
        )
        return common.reset_index(drop=True), a_spec, b_spec
    # -----------------------------------------------------------------------------


    # Same loop-strength definition as scripts/calculate_loop_strength.py
    def _quantify_loops(mtx):
        sq_size = mtx.shape[0]
        midpoint = int(np.floor(sq_size / 2))
        mid_9pixels_mean = np.nanmean(mtx[midpoint - 1:midpoint + 2, midpoint - 1:midpoint + 2])
        neighboring_size = int(np.ceil(0.3 * sq_size) // 2 * 2 + 1)
        upper_left_mean = np.nanmean(mtx[:neighboring_size, :neighboring_size])
        upper_right_mean = np.nanmean(mtx[:neighboring_size, (sq_size - neighboring_size):])
        lower_right_mean = np.nanmean(mtx[(sq_size - neighboring_size):, (sq_size - neighboring_size):])
        bg = np.nanmean(np.array([upper_left_mean, upper_right_mean, lower_right_mean]))
        return mid_9pixels_mean / bg

    def _strength_from_stack(stack, n_loops):
        if n_loops == 0:
            return np.array([], dtype=float)

        stack = np.asarray(stack)
        if stack.ndim != 3:
            raise ValueError(f"Unexpected pileup output shape: {stack.shape}")

        if stack.shape[0] == n_loops:        # cooltools >= 0.6
            loop_view = stack
        elif stack.shape[2] == n_loops:      # older cooltools
            loop_view = np.moveaxis(stack, 2, 0)
        else:
            raise ValueError(f"Cannot infer loop axis from shape {stack.shape}, n_loops={n_loops}")

        return np.array([_quantify_loops(loop_view[i, :, :]) for i in range(n_loops)], dtype=float)

    def _build_strength_runner(clr, expected_df, view_df):
        def _run(loops_df):
            loops_df = loops_df.reset_index(drop=True)
            if loops_df.empty:
                return pd.Series([], dtype=float)
            stack = cooltools.pileup(
                clr,
                loops_df[req],
                view_df=view_df,
                expected_df=expected_df,
                flank=flank,
                clr_weight_name=clr_weight_name,
                nproc=nproc,
            )
            return pd.Series(_strength_from_stack(stack, len(loops_df)), index=loops_df.index)
        return _run

    a = _prepare(df_a, "_id_a")
    b = _prepare(df_b, "_id_b")

    keys = ["chrom1", "chrom2", "_bin1", "_bin2"]
    a_small = a[["chrom1", "chrom2", "_bin1", "_bin2", "_id_a", "start1", "start2"]]
    b_small = b[["chrom1", "chrom2", "_bin1", "_bin2", "_id_b", "start1", "start2"]]

    candidates = []
    for d1 in (-1, 0, 1):
        for d2 in (-1, 0, 1):
            b_shift = b_small.copy()
            b_shift["_bin1"] += d1
            b_shift["_bin2"] += d2
            m = a_small.merge(b_shift, on=keys, how="inner", suffixes=("_a", "_b"))
            if m.empty:
                continue
            d_start1 = (m["start1_a"] - m["start1_b"]).abs()
            d_start2 = (m["start2_a"] - m["start2_b"]).abs()
            keep = (d_start1 <= clustering_radius) & (d_start2 <= clustering_radius)
            if not keep.any():
                continue
            m = m.loc[keep, ["_id_a", "_id_b"]].copy()
            m["delta_start1"] = d_start1[keep].to_numpy()
            m["delta_start2"] = d_start2[keep].to_numpy()
            m["dist_bp"] = np.maximum(m["delta_start1"], m["delta_start2"])
            candidates.append(m)

    if candidates:
        cands = pd.concat(candidates, ignore_index=True).drop_duplicates(["_id_a", "_id_b"])
        cands = cands.sort_values(["dist_bp", "delta_start1", "delta_start2", "_id_a", "_id_b"])
        used_a, used_b, keep_rows = set(), set(), []
        for _id_a, _id_b, ds1, ds2, dist in cands[
            ["_id_a", "_id_b", "delta_start1", "delta_start2", "dist_bp"]
        ].itertuples(index=False, name=None):
            if _id_a in used_a or _id_b in used_b:
                continue
            used_a.add(_id_a)
            used_b.add(_id_b)
            keep_rows.append((_id_a, _id_b, ds1, ds2, dist))
        matched = pd.DataFrame(keep_rows, columns=["_id_a", "_id_b", "delta_start1", "delta_start2", "dist_bp"])
    else:
        matched = pd.DataFrame(columns=["_id_a", "_id_b", "delta_start1", "delta_start2", "dist_bp"])

    matched_a = set(matched["_id_a"].tolist())
    matched_b = set(matched["_id_b"].tolist())

    a_specific = a.loc[~a["_id_a"].isin(matched_a), req].reset_index(drop=True)
    b_specific = b.loc[~b["_id_b"].isin(matched_b), req].reset_index(drop=True)

    common = (
        matched
        .merge(a[["_id_a"] + req], on="_id_a", how="left")
        .merge(b[["_id_b"] + req], on="_id_b", how="left", suffixes=(f"_{sample_a_name}", f"_{sample_b_name}"))
        .reset_index(drop=True)
    )

    # --- NEW: compute anchor common/specific sets ---------------------------------
    anchors_a = _prepare_anchors(a, "_id_anchor_a")
    anchors_b = _prepare_anchors(b, "_id_anchor_b")
    common_anchors, a_specific_anchors, b_specific_anchors = _match_anchors(
        anchors_a, anchors_b, "_id_anchor_a", "_id_anchor_b"
    )
    # -----------------------------------------------------------------------------

    if compute_strength:
        clrA = _to_cooler(clr_a)
        clrB = _to_cooler(clr_b)
        run_a = _build_strength_runner(clrA, expected_a, view_df_a)
        run_b = _build_strength_runner(clrB, expected_b, view_df_b)

        a_specific[f"{sample_a_name}_loop_strength"] = run_a(a_specific).values
        b_specific[f"{sample_b_name}_loop_strength"] = run_b(b_specific).values

        if not common.empty:
            loops_a_common = common[
                [f"chrom1_{sample_a_name}", f"start1_{sample_a_name}", f"end1_{sample_a_name}",
                 f"chrom2_{sample_a_name}", f"start2_{sample_a_name}", f"end2_{sample_a_name}"]
            ].copy()
            loops_a_common.columns = req

            loops_b_common = common[
                [f"chrom1_{sample_b_name}", f"start1_{sample_b_name}", f"end1_{sample_b_name}",
                 f"chrom2_{sample_b_name}", f"start2_{sample_b_name}", f"end2_{sample_b_name}"]
            ].copy()
            loops_b_common.columns = req

            common[f"{sample_a_name}_loop_strength"] = run_a(loops_a_common).values
            common[f"{sample_b_name}_loop_strength"] = run_b(loops_b_common).values
            common["log2_strength_fc"] = np.log2(
                (common[f"{sample_a_name}_loop_strength"] + 1e-12) /
                (common[f"{sample_b_name}_loop_strength"] + 1e-12)
            )

    return {
        f"{sample_a_name}_specific": a_specific,
        f"{sample_b_name}_specific": b_specific,
        "common_loops": common,
        "common_anchors": common_anchors,
        f"{sample_a_name}_specific_anchors": a_specific_anchors,
        f"{sample_b_name}_specific_anchors": b_specific_anchors,
    }

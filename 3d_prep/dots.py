import cooltools
import pandas as pd
from ggner_3d.cm import CLR_CONNS, prepare_view_df, get_expected_cis_df, assert_all_cools_has_weights

RESOLUTION = 10000
WEIGHT = 'sweight'
NPROC = 6
DATA_DIR = "/home/carlos/Clone/ggner-3d/data/dots"

if __name__ == "__main__":
    # open connections to cooler files
    clrs_ = CLR_CONNS(RESOLUTION)
    assert_all_cools_has_weights(clrs_, WEIGHT)

    # prepare view dataframe
    view_df = prepare_view_df(arm=True)

    # compute expected cis
    expected_ = [get_expected_cis_df(condition=k, resolution_kb=RESOLUTION//1000, label='arm') for k in clrs_.keys()]

    # compute and save dots
    for (k, v), exp in zip(clrs_.items(), expected_):

        print(f"Processing {k}...")

        dots_df = cooltools.dots(
            v,
            expected=exp,
            view_df=view_df,
            max_loci_separation=10_000_000,
            clr_weight_name=WEIGHT,
            nproc=NPROC,
        )

        dots_df.to_csv(
            f"{DATA_DIR}/dots_{k}_{RESOLUTION}.csv",
            index=False
        )
        print(
            f"Dots for {k} at {RESOLUTION}bp saved."
        )
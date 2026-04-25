import cooler
import cooltools

from ggner_3d.cm import MCOOL_PATH_DICT_BIOREPS, assert_all_cools_has_weights, prepare_view_df

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

NPROC = 6
RESOLUTION = 10_000
WEIGHT_NAME = "sweight"
DATA_DIR = "/home/carlos/Clone/ggner-3d/data/expected_reps"
ARM = False
DOTRANS = False

SAMPLES_OF_INTEREST = [
    'WT_noUV_r1',
    'WT_noUV_r2',
    'WT_3h_r1',
    'WT_3h_r2',
    'XPC_noUV_r1',
    'XPC_noUV_r2',
    'XPC_3h_r1',
    'XPC_3h_r2',
    'XPA_noUV_r1',
    'XPA_noUV_r2',
    'XPA_3h_r1',
    'XPA_3h_r2',
]


if __name__ == "__main__":

    # open connections to cooler files
    CLR_ = {k: cooler.Cooler(f'{v}::resolutions/{RESOLUTION}') for k, v in MCOOL_PATH_DICT_BIOREPS.items() if k in SAMPLES_OF_INTEREST}

    # prepare view df
    view_df_arm = prepare_view_df(arm=ARM)

    # check that all coolers have the required weight column
    assert_all_cools_has_weights(CLR_, WEIGHT_NAME)

    # compute and save expected cis
    for k, v in CLR_.items():

        print(f"Processing {k}...")

        expected_df_cis = cooltools.expected_cis(
            clr=v,
            nproc=NPROC,
            view_df = view_df_arm,
            clr_weight_name=WEIGHT_NAME,
        )
        if ARM:
            label = 'arm'
        else:
            label = 'chrom'
        save_path = f"{DATA_DIR}/expected_cis_{label}_{k}_res{RESOLUTION//1000}kb.tsv"
        expected_df_cis.to_csv(save_path, sep="\t", index=False)
        print(f"Saved to {save_path}")

        if not DOTRANS:
            continue

        expected_df_trans = cooltools.expected_trans(
            clr=v,
            nproc=NPROC,
            view_df=None,
            clr_weight_name=WEIGHT_NAME,
        )
        save_path = f"{DATA_DIR}/expected_trans_{k}_res{RESOLUTION//1000}kb.tsv"
        expected_df_trans.to_csv(save_path, sep="\t", index=False)
        print(f"Saved to {save_path}")


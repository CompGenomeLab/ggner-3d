import cooler
import cooltools

from ggner_3d.cm import MCOOL_PATH_DICT_BIOREPS, assert_all_cools_has_weights

NPROC = 8
RESOLUTION = 10000
WEIGHT_NAME = "sweight"
DATA_DIR = "/home/carlos/Clone/ggner-3d/data/insulation_reps_2"
WINDOWS_FACTOR = [50] # 10, 25, 50, 100
INSULATION_WINDOWS = [
    RESOLUTION * w for w in WINDOWS_FACTOR
]

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

    # check that all coolers have the required weight column
    assert_all_cools_has_weights(CLR_, WEIGHT_NAME)

    # compute and save insulation
    for k, v in CLR_.items():

        print(f"Processing {k}...")

        insulation_df = cooltools.insulation(
            clr = v,
            window_bp = INSULATION_WINDOWS,
            clr_weight_name='sweight',
            nproc=NPROC
        )
        insulation_df.to_csv(
            f"{DATA_DIR}/insulation_{k}_{RESOLUTION}bp.csv",
            index=False
        )
        print(
            f"Insulation for {k} at {RESOLUTION}bp saved."
        )


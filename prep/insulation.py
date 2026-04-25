import cooltools

from ggner_3d.cm import CLR_CONNS, assert_all_cools_has_weights

NPROC = 6
RESOLUTION = 10000
WEIGHT_NAME = "sweight"
DATA_DIR = "/home/carlos/Clone/ggner-3d/data/insulation"
INSULATION_WINDOWS = [
    RESOLUTION * w for w in [
        10, 25, 50, 100
    ]
]
if __name__ == "__main__":
    
    # open connections to cooler files
    CLR_ = CLR_CONNS(RESOLUTION)

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


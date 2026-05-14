# GG-NER and 3D Genome Analysis

This repository contains analysis code for the paper "XPC organises the 3D genome for nucleotide excision repair".

All analyses were performed on an Arch Linux system. To facilitate reproducibility of the analyses and manuscript figures, the corresponding code has been deposited primarily as Python-based Jupyter notebooks, with cell outputs preserved, including the generated figures. Additional shell and R scripts used in the analysis are also provided.

# Package Management

Python packages were managed using [uv](https://docs.astral.sh/uv/). The required dependencies can be installed using the following command:

```bash
uv sync --dev
```

R packages were managed using [renv](https://rstudio.github.io/renv/). The package versions used in this analysis are recorded in `r_env/renv.lock`.

To install the required R packages and activate the R environment, navigate to the `r_env/` directory and start R:

```bash
cd r_env
R
```

Then restore the package environment from the lockfile:

```r
renv::restore()
```

After restoration, the `renv` environment will be loaded automatically whenever R is started from the `r_env/` directory.


# Repository structure
- 'src/ggner_3d/': Contains modules for 
- `3d_prep/`: Contains scripts for identifying architectural features on contact maps, including compartment identification, insulation based TAD calling, and loop calling. 
- `notebooks/`: Contains Jupyter notebooks for each analysis step, with outputs and figures included.
    - `data_processing/`: Notebooks for processing raw data, including quality control and normalization steps.

# Data availability
Please see the data availability statement in the manuscript for details on how to access the datasets used in this study.

# Contact
Please email [vogulcan@sabanciuniv.edu](vogulcan@sabanciuniv.edu) or raise an issue in the github repository with any questions about installation or usage.
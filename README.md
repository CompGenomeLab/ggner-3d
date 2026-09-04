# NER and 3D Genome Analysis

This repository contains analysis code for the paper "RNAPII and NER stall loop extrusion at UV lesions, shaping the 3D genome during repair".

All analyses were performed on an Arch Linux system. To facilitate reproducibility of the analyses and manuscript figures, the corresponding code has been deposited primarily as Python-based Jupyter notebooks, with cell outputs preserved, including the generated figures. Additional shell and R scripts used in the analysis are also provided.

Please see the other repositories related to this study for additional research code:
- [1D Simulation Engine](https://github.com/CompGenomeLab/ner_rnapii_3d)
- [PHACE-MP](https://github.com/CompGenomeLab/phace-mp)

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
- `src/ggner_3d/`: Contains helper modules for the analysis, including; data loading, processing, and plotting functions.
- `3d_prep/`: Contains scripts for identifying architectural features on contact maps, including; expected frequency calculation, compartment identification, insulation based TAD calling, and loop calling. 
- `notebooks/`: Contains Jupyter notebooks for each analysis step, with outputs and figures included.
- `r_env/`: Contains the R environment and package lockfile for reproducibility.
- `scripts/`: Contains additional shell and R scripts used in the analysis.
  - `scripts/atacseq_nfcore/`: Runs the nf-core/atacseq pipeline and generates RPGC-normalized ATAC-seq signal tracks.
  - `scripts/ctcf_chipseq_nfcore/`: Runs the nf-core/chipseq pipeline for CTCF ChIP-seq and generates normalized signal and input-ratio tracks.
  - `scripts/rnaseq_nfcore/`: Runs the nf-core/rnaseq pipeline and generates CPM-normalized RNA-seq signal tracks.
  - `scripts/xpc_chipseq_nfcore/`: Runs the nf-core/chipseq pipeline for XPC ChIP-seq and generates RPGC-normalized signal tracks.

# Todo:
- Prior to publication, we will reorganize the repository, remove temporary and redundant files, and provide a more detailed overview of its structure, including the purpose of each script and notebook.

# Data availability
Please see the data availability statement in the manuscript for details on how to access the datasets used in this study.

# Contact
Please email [vogulcan@sabanciuniv.edu](vogulcan@sabanciuniv.edu) or raise an issue in the github repository with any questions about installation or usage.

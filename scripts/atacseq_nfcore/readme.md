The scripts in this folder are used to run the nf-core/atacseq pipeline and to generate signal tracks from the resulting alignment files.

1. **run_nf.sh**: Runs the nf-core/atacseq pipeline for the ATAC-seq samples.
2. **atac_rpgc.sh**: Generates RPGC-normalized bigWig files from the pipeline alignment files for downstream analysis in the manuscript, with optional Tn5 shifting and configurable bin size, smoothing length, and minimum mapping quality.

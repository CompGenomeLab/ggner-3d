import cooltools
import bioframe
import os
import subprocess
import pandas as pd
import multiprocess as mp

from ggner_3d.cm import CLR_CONNS, prepare_view_df, assert_all_cools_has_weights

NPROC = 6
RESOLUTION = 100000
WEIGHT_NAME = "sweight"
DATA_DIR = "/home/carlos/Clone/ggner-3d/data/compartment"

if __name__ == "__main__":
    
    # open connections to cooler files
    CLR_ = CLR_CONNS(RESOLUTION)

    # prepare view df
    view_df = prepare_view_df(arm=False)

    # check that all coolers have the required weight column
    assert_all_cools_has_weights(CLR_, WEIGHT_NAME)

    # download fasta if not present
    print("Checking for hg38 fasta...")
    fasta_path = '/home/carlos/Clone/ggner-3d/data/hg38.fa'
    if not os.path.isfile(fasta_path):
        subprocess.call(f'wget --progress=bar:force:noscroll https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -O {fasta_path}.gz', shell=True)
        subprocess.call(f'gunzip {fasta_path}.gz', shell=True)

    # prepare GC content
    print("Preparing GC content...")
    gc_path = os.path.join(DATA_DIR, f'hg38_gc_cov_{RESOLUTION//1000}kb.tsv')
    if not os.path.isfile(gc_path):
        bins = list(CLR_.values())[0].bins()[:]
        hg38_genome = bioframe.load_fasta(fasta_path)
        gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], hg38_genome)
        gc_cov.to_csv(gc_path, index=False, sep='\t')
    else:
        gc_cov = pd.read_csv(gc_path, sep='\t')
    
    # compute and save expected cis
    for k, v in CLR_.items():

        print(f"Processing {k}...")
        pool = mp.Pool(NPROC)
        
        eigvals, eigvec_table = cooltools.eigs_cis(
            clr=v,
            phasing_track=gc_cov,
            view_df=view_df,
            n_eigs=3,
            clr_weight_name=WEIGHT_NAME,
            sort_metric='MAD_explained',
            map=pool.map
        )

        eigvals.to_csv(
            os.path.join(DATA_DIR, f'compartment_eigvals_{k}_{RESOLUTION}.tsv'),
            sep='\t',
            index=False
        )

        eigvec_table.to_csv(
            os.path.join(DATA_DIR, f'compartment_eigvecs_{k}_{RESOLUTION}.tsv'),
            sep='\t',
            index=False
        )
        print(f"Saved eigvals and eigvecs for {k}.")
        pool.close()

    print("All done.")

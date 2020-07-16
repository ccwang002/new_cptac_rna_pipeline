export PATH="/opt/conda/bin:$PATH"

snakemake \
    --profile=/storage1/fs1/dinglab/Active/Projects/CPTAC3/Analysis/bulk_rna_pipeline_development/new_cptac_rna_pipeline/compute1_lsf \
    --configfile=snakemake_config.json \
    -s /storage1/fs1/dinglab/Active/Projects/CPTAC3/Analysis/bulk_rna_pipeline_development/new_cptac_rna_pipeline/star.smk \
    -j 8 \
    -d $PWD \
    --restart-times 2 \
    --resources io_heavy=4 \
    -- \
    star_align_all_samples

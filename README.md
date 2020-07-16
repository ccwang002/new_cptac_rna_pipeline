# CPTAC Bulk RNA Expression pipeline
The pipeline aligns bulk RNA-seq reads and produces aligned BAMs, readcount, and FPKM.

Based on [GTEx RNA pipeline](https://github.com/broadinstitute/gtex-pipeline) compatible to GDC reference.


## Release
- Still developing

## Output
Each batch of execution will produce a TSV `analysis_summary.dat` containing all results.

Possible result types are:
- genomic_bam
- transcriptomic_bam
- chimeric_sam
- fpkm_tsv
- splic_junction_tab



## Run the pipeline

### Setup
The easiest way is to start a conda environment with all the dependencies:
```
conda create -n cptac_gtex_rna python=3.8 \
    snakemake-minimal=5.20.1 \
    pandas=1.0.5 \
    star=2.6.1d \
    samtools=1.10 htslib=1.10 \
    subread=2.0.1 \
    picard=2.21.4 \
    rsem=1.3.1
```


### Create a new batch
1. Copy `example_batch/` to the desired location to store the output
2. Modify the `snakemake_config.json` to ensure all file paths exist
3. Define the file map and the list of samples to run the pipeline (same format as the example)

```
# Create the result summary of the alignment outputs and readcount TSVs
snakemake --configfile=snakemake_config.json -s ../pipeline_workflow/Snakefile \
    --cores 54 -p \
    --resouces io_heavy=5 -- \
    make_analysis_summary

# Only the alignment
snakemake ... star_align_all_samples

# All readcount and FPKMs
snakemake ... all_fpkms
```


## Processing description

### Genome alignment
STAR v2.6.1d with GDC's genome reference [GRCh38.d1.vd1][GDC Reference Files] and [GENCODE v22][gencode-gtf].

[GDC Reference Files]: https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
[gencode-gtf]: https://api.gdc.cancer.gov/data/25aa497c-e615-4cb7-8751-71f744f9691f
[GDC's formula]: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#upper-quartile-fpkm





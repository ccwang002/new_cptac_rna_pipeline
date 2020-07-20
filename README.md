# CPTAC Bulk RNA Expression pipeline (experimental)
The pipeline aligns bulk RNA-seq reads using STAR and quantifies gene level expression using RNASeQC and isoform level expression using RSEM. It's based on [GTEx RNA pipeline](https://github.com/broadinstitute/gtex-pipeline), but the annotation is compatible to GDC's reference.

You probably want to use the original pipeline. This repo was created for me to understand the details of the pipeline.


## Release
- Still developing


## Run the pipeline

### Setup
The following command won't work due to package conflicts. See `Dockerfile`:

```
conda create -n cptac_gtex_rna python=3.8 \
    snakemake-minimal=5.20.1 \
    pandas=1.0.5 \
    star=2.6.1d \
    samtools=1.10 htslib=1.10 \
    picard=2.21.4 \
    rsem=1.3.1
```


### Create a new batch
1. Copy `example_batch/` to the desired location to store the output
2. Modify the `snakemake_config.json` to ensure all file paths exist
3. When using LSF, modifiy `lsf.yaml` to ensure all variables and paths are valid
4. When not using LSF, ensure `$TMPDIR` is set to a sufficiently large space (> 300GB)
5. Define the file map and the sample list to run the pipeline

```
# When using LSF,
bash run.sh

# When not using LSF,
snakemake --configfile=snakemake_config.json \
    -s /path/to/repo/rnaseq.smk \
    -j 54 -p \
    --resouces io_heavy=5 -- \
    star_align_all_samples expression_all_samples

# Generate the alignment manifest
snakemake ... make_alignment_output_manifest
```

Clean up the outputs:

```
# Further compress the outputs and logs
rm -rf .snakemake/

cd rsem/
parallel -j4 --bar 'tar czf {}.tar.gz {}' ::: *.rsem.stat
rm -rf *.rsem.stat

cd ../logs
rm -rf cluster rnaseqc
gzip -9 */*.log
```


## Annotation
STAR v2.6.1d with GDC's genome reference [GRCh38.d1.vd1][GDC Reference Files] and [GENCODE v22][gencode-gtf].

[GDC Reference Files]: https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
[gencode-gtf]: https://api.gdc.cancer.gov/data/25aa497c-e615-4cb7-8751-71f744f9691f

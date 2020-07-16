bsub \
    -n 4 -q general -G compute-dinglab \
    -M 8GB -R 'select[mem>8000] rusage[mem=8000] span[hosts=1]' \
    -N -u 'liang-bo.wang@wustl.edu' \
    -a 'docker(lbwang/cptac_rna_expression)' \
    bash run_snakemake_main_job.sh

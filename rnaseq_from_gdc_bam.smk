import csv
from pathlib import Path

from find_input import BAMFileMapping as FileMapping

WORKFLOW_ROOT = config['workflow_root']
GENE_LEVEL_GTF_PTH = config['gene_level_gtf']
STAR_INDEX_FOLDER = config['star_index']  # Path to the STAR index
RSEM_REF_PREFIX = config['rsem_ref']  # Path to the RSEM index prefix


filemap = FileMapping(
    # The mapping of sample name to other information.
    config['file_map'],
    config['catalog'],
    # List of sample names to run the pipeline
    config['sample_list'],
)


def find_sample_genomic_bam(wildcards):
    return str(filemap.bam_map[wildcards.sample]['genomic'])


def find_sample_transcriptome_bam(wildcards):
    return str(filemap.bam_map[wildcards.sample]['transcriptome'])


rule samtools_index_bam:
    """Index a sorted BAM by samtools."""
    output: '{name}.bam.bai'
    input: '{name}.bam'
    resources:
        io_heavy=1
    shell: 'samtools index {input} {output}'


rule picard_mark_dup:
    output:
        bam='star/{sample}/Aligned.sortedByCoord.out.md.bam',
        metrics='star/{sample}/Aligned.sortedByCoord.out.marked_dup_metrics.txt'
    input: bam=find_sample_genomic_bam
    threads: 8
    resources:
        io_heavy=1,
        mem_mb=lambda wildcards, attempt: 32000 + 16000 * (attempt - 1),
        tmp_mb=32000
    params:
        jar_pth="/usr/local/lib/picard-tools/picard.jar",
        java_mem_mb=lambda wildcards, resources: resources.mem_mb - 500
    log: 'logs/picard_mark_dup/{sample}.log'
    shell:
        "java -jar -Xmx{params.java_mem_mb}m {params.jar_pth} "
        "MarkDuplicates "
        "I={input.bam} O={output.bam} "
        "PROGRAM_RECORD_ID=null "
        "MAX_RECORDS_IN_RAM=500000 "
        "SORTING_COLLECTION_SIZE_RATIO=0.25 "
        "TMP_DIR=$(mktemp -d) "
        "M={output.metrics} "
        "ASSUME_SORT_ORDER=coordinate "
        "TAGGING_POLICY=DontTag "
        "OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 "
        "2>{log} 1>&2"


def expand_to_all_samples(patterns):
    return {
        name: expand(pattern, sample=filemap.samples)
        for name, pattern in patterns.items()
    }


rule rsem_calc_expression:
    output:
        genes=temporary("rsem/{sample}.rsem.genes.results"),
        isoforms=temporary("rsem/{sample}.rsem.isoforms.results")
    input:
        bam=find_sample_transcriptome_bam
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: 36000 + 8000 * (attempt - 1),
        tmp_mb=32000
    params:
        rsem_ref_prefix=RSEM_REF_PREFIX
    log: 'logs/rsem/{sample}.log'
    group: "rsem"
    shell:
        "rsem-calculate-expression "
        "--num-threads {threads} "
        "--temporary-folder $(mktemp -d) "
        "--no-bam-output "
        "--fragment-length-max 1000 "
        "--paired-end "
        "--estimate-rspd "
        "--forward-prob 0.0 "  # Stranded
        "--bam {input.bam} "
        "{params.rsem_ref_prefix} "
        "rsem/{wildcards.sample}.rsem "
        "2>{log} 1>&2"


rule gzip_rsem_outputs:
    output: "rsem/{name}.results.gz"
    input: "rsem/{name}.results"
    group: "rsem"
    shell: "gzip -9n -c {input} > {output}"


rule rnaseqc:
    output:
        metrics="rnaseqc/{sample}.metrics.tsv",
        exon_reads=temporary("rnaseqc/{sample}.exon_reads.gct"),
        gene_fragments=temporary("rnaseqc/{sample}.gene_fragments.gct"),
        gene_reads=temporary("rnaseqc/{sample}.gene_reads.gct"),
        gene_tpm=temporary("rnaseqc/{sample}.gene_tpm.gct"),
    input:
        bam=rules.picard_mark_dup.output.bam,
        bai=rules.picard_mark_dup.output.bam + '.bai'
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 8000 + 8000 * (attempt - 1),
    params:
        gene_level_gtf=GENE_LEVEL_GTF_PTH
    log: 'logs/rnaseqc/{sample}.log'
    group: "rnaseqc"
    shell:
        "rnaseqc "
        "-s {wildcards.sample} "
        "--stranded rf "
        "-vv "
        "-- "
        "{params.gene_level_gtf} "
        "{input.bam} "
        "rnaseqc "
        "2>{log} 1>&2"


rule gzip_rnaseqc_gcts:
    output: "rnaseqc/{name}.gct.gz"
    input: "rnaseqc/{name}.gct"
    group: "rnaseqc"
    shell: "gzip -9n -c {input} > {output}"


rule rsem_all_samples:
    input: **expand_to_all_samples({ \
        'rsem_genes': rules.rsem_calc_expression.output.genes + '.gz', \
        'rsem_isoforms': rules.rsem_calc_expression.output.isoforms + '.gz', \
    })


rule rnaseqc_all_samples:
    input: **expand_to_all_samples({ \
        'rnaseqc_exon_reads': rules.rnaseqc.output.exon_reads + '.gz', \
        'rnaseqc_gene_fragments': rules.rnaseqc.output.gene_fragments + '.gz', \
        'rnaseqc_gene_reads': rules.rnaseqc.output.gene_reads + '.gz', \
        'rnaseqc_gene_tpm': rules.rnaseqc.output.gene_tpm + '.gz', \
    })

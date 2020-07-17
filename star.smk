from textwrap import dedent
from pathlib import Path

from find_input import FileMapping

WORKFLOW_ROOT = config['workflow_root']
GENE_GTF_PTH = config['gdc_gtf']
STAR_INDEX_FOLDER = config['star_index']  # Path to the STAR index


filemap = FileMapping(
    # The mapping of sample name to other information.
    config['file_map'],
    config['catalog'],
    # List of sample names to run the pipeline
    config['sample_list'],
)


def find_sample_fqs(wildcards):
    """Find the FASTQ file paths of a given sample."""
    fqs = filemap.fastq_map[wildcards.sample]
    return {
        'r1_fq': str(fqs['R1']),
        'r2_fq': str(fqs['R2']),
    }


def create_rg_line(wildcards, input):
    fq_name = Path(input.r1_fq).name
    m = re.search(r'_R[12]_\d+\.fastq\.gz$', fq_name)
    rg = fq_name[:m.start()]
    return f"ID:{rg} SM:{wildcards.sample}"


rule star_align:
    """STAR align one sample."""
    output:
        unsorted_bam=temporary('star/{sample}/Aligned.out.bam'),
        # sorted_bam='star/{sample}/Aligned.sortedByCoord.out.bam',
        chimeric_sam=temporary('star/{sample}/Chimeric.out.sam'),
        chimeric_junction='star/{sample}/Chimeric.out.junction',
        quant_tx_bam='star/{sample}/Aligned.toTranscriptome.out.bam',
        quant_gene_count_tab=temporary('star/{sample}/ReadsPerGene.out.tab'),
        sj_count_tab=temporary('star/{sample}/SJ.out.tab'),
    input: unpack(find_sample_fqs)
    params:
        star_ix=STAR_INDEX_FOLDER,
        out_folder='star/{sample}/',
        outSAMattrRGline=create_rg_line
    log: 'logs/star/{sample}.log'
    threads: 8
    resources:
        io_heavy=1,
        mem_mb=lambda wildcards, attempt: 32000 + 16000 * (attempt - 1)
    shell:
        "STAR "
        "--readFilesIn {input.r1_fq} {input.r2_fq} "
        # Most parameters follow GDC
        "--alignIntronMax 1000000 "
        "--alignIntronMin 20 "
        "--alignMatesGapMax 1000000 "
        "--alignSJDBoverhangMin 1 "
        "--alignSJoverhangMin 8 "
        "--alignSoftClipAtReferenceEnds Yes "

        # Follow arriba's recommendation regarding chimera parameters
        # Ref: https://arriba.readthedocs.io/en/latest/workflow/
        "--chimJunctionOverhangMin 10 "
        "--chimMainSegmentMultNmax 1 "
        "--chimOutType Junctions SeparateSAMold WithinBAM SoftClip "
        "--chimOutJunctionFormat 1 "
        "--chimSegmentMin 10 "
        "--chimScoreMin 1"
        "--chimScoreDropMax 30 "
        "--chimScoreJunctionNonGTAG 0 "
        "--chimScoreSeparation 1 "
        "--alignSJstitchMismatchNmax 5 -1 5 5 "
        "--chimSegmentReadGapMax 3 "

        "--genomeDir {params.star_ix} "
        "--genomeLoad NoSharedMemory "
        "--limitBAMsortRAM 0 "
        "--limitSjdbInsertNsj 1200000 "
        "--outFileNamePrefix {params.out_folder} "
        "--outFilterIntronMotifs None "
        "--outFilterMatchNminOverLread 0.33 "
        "--outFilterMismatchNmax 999 "
        "--outFilterMismatchNoverLmax 0.1 "
        "--outFilterMultimapNmax 20 "
        "--outFilterScoreMinOverLread 0.33 "
        "--outFilterType BySJout "
        "--outSAMattributes NH HI AS nM NM ch "
        "--outSAMattrRGline {params.outSAMattrRGline} "
        "--outSAMstrandField intronMotif "
        "--outSAMtype BAM Unsorted "
        "--outSAMunmapped Within "
        "--quantMode TranscriptomeSAM GeneCounts "
        "--readFilesCommand zcat "
        "--runThreadN {threads} "
        "--twopassMode Basic "
        "> {log}"


rule samtools_index_bam:
    """Index a sorted BAM by samtools."""
    output: '{name}.bam.bai'
    input: '{name}.bam'
    resources:
        io_heavy=1
    shell: 'samtools index {input} {output}'


rule samtools_sort_star_bam:
    output: 'star/{sample}/Aligned.sortedByCoord.out.bam'
    input: rules.star_align.output.unsorted_bam
    threads: 8
    resources:
        io_heavy=1,
        mem_mb=lambda wildcards, attempt: 16000 + 8000 * (attempt - 1),
        tmp_mb = 50000
    shell:
        "samtools sort "
        "--threads {threads} -m 1500M "
        "-T $(mktemp -d) "
        "-o {output} {input}"


rule samtools_sort_star_chimeric_bam:
    output: 'star/{sample}/Chimeric.out.sorted.bam'
    input: rules.star_align.output.chimeric_sam
    threads: 4
    resources:
        io_heavy=1,
        mem_mb=lambda wildcards, attempt: 4000 + 8000 * (attempt - 1),
        tmp_mb = 8000
    shell:
        "samtools sort "
        "--threads {threads} -m 1500M "
        "-T $(mktemp -d) "
        "-o {output} {input}"


rule gzip_star_quant_tab:
    output: 'star/{sample}/ReadsPerGene.out.tab.gz'
    input: rules.star_align.output.quant_gene_count_tab
    shell: "gzip -9n -c {input} > {output}"


rule gzip_star_sj_tab:
    output: 'star/{sample}/SJ.pass1.out.tab.gz'
    input: rules.star_align.output.sj_count_tab
    shell: "gzip -9n -c {input} > {output}"


rule gzip_star_chimeric_junction:
    output: 'star/{sample}/Chimeric.out.junction.gz'
    input: rules.star_align.output.chimeric_junction
    shell: "gzip -9n -c {input} > {output}"


def expand_to_all_samples(patterns):
    return {
        name: expand(pattern, sample=filemap.samples)
        for name, pattern in patterns.items()
    }


rule star_align_all_samples:
    """Align all RNA-seq samples."""
    input:
        **expand_to_all_samples({ \
            "sorted_bams": rules.samtools_sort_star_bam.output[0], \
            "sorted_bam_bais": rules.samtools_sort_star_bam.output[0] + '.bai', \
            "chimeric_bams": rules.samtools_sort_star_chimeric_bam.output[0], \
            "chimeric_bam_bais": rules.samtools_sort_star_chimeric_bam.output[0] + '.bai', \
            "chimeric_junction_gzs": rules.gzip_star_chimeric_junction.output[0], \
            "quant_tx_bams": rules.star_align.output.quant_tx_bam, \
            "quant_gene_count_tab_gzs": rules.gzip_star_quant_tab.output[0], \
            "sj_count_tab_gzs": rules.gzip_star_sj_tab.output[0] \
        })

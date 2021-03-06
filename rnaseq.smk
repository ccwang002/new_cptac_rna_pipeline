import csv
from pathlib import Path

from find_input import FASTQFileMapping as FileMapping

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
        chimeric_junction=temporary('star/{sample}/Chimeric.out.junction'),
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
        mem_mb=lambda wildcards, attempt: 40000 + 8000 * (attempt - 1),
        tmp_mb=32000
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
        "--outTmpDir $(mktemp -d)/_STARtmp "
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
        mem_mb=lambda wildcards, attempt: 32000 + 8000 * (attempt - 1),
        tmp_mb=50000
    shell:
        "samtools sort "
        "--threads {threads} "
        # it uses much more memory than what's specified below
        "-m 1400M "
        "-T $(mktemp -d) "
        "-o {output} {input}"


rule samtools_sort_star_chimeric_bam:
    output: 'star/{sample}/Chimeric.out.sorted.bam'
    input: rules.star_align.output.chimeric_sam
    threads: 4
    resources:
        io_heavy=1,
        mem_mb=lambda wildcards, attempt: 4000 + 8000 * (attempt - 1),
        tmp_mb=8000
    shell:
        "samtools sort "
        "--threads {threads} -m 1400M "
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


rule picard_mark_dup:
    output:
        bam='star/{sample}/Aligned.sortedByCoord.out.md.bam',
        metrics='star/{sample}/Aligned.sortedByCoord.out.marked_dup_metrics.txt'
    input: bam=rules.samtools_sort_star_bam.output
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


rule star_align_all_samples:
    """Align all RNA-seq samples."""
    input:
        **expand_to_all_samples({ \
            "sorted_bams": rules.samtools_sort_star_bam.output[0], \
            "sorted_bam_bais": rules.samtools_sort_star_bam.output[0] + '.bai', \
            # "mark_dup_bams": rules.picard_mark_dup.output.bam, \
            # "mark_dup_bam_bais": rules.picard_mark_dup.output.bam + '.bai', \
            # "mark_dup_metric_txts": rules.picard_mark_dup.output.metrics, \
            "chimeric_bams": rules.samtools_sort_star_chimeric_bam.output[0], \
            "chimeric_bam_bais": rules.samtools_sort_star_chimeric_bam.output[0] + '.bai', \
            "chimeric_junction_gzs": rules.gzip_star_chimeric_junction.output[0], \
            "quant_tx_bams": rules.star_align.output.quant_tx_bam, \
            "quant_gene_count_tab_gzs": rules.gzip_star_quant_tab.output[0], \
            "sj_count_tab_gzs": rules.gzip_star_sj_tab.output[0] \
        })


rule make_alignment_output_manifest:
    output:
        manifest='alignment_summary.dat'
    input: rules.star_align_all_samples.input
    run:
        result_file_tpls = {
            ('genomic_bam', 'BAM'): rules.samtools_sort_star_bam.output[0],
            ('transcriptomic_bam', 'BAM'): rules.star_align.output.quant_tx_bam,
            ('chimeric_bam', 'BAM'): rules.samtools_sort_star_chimeric_bam.output[0],
            ('chimeric_junction', 'TSV'): rules.gzip_star_chimeric_junction.output[0],
            ('STAR_gene_count_tab', 'TSV'): rules.gzip_star_quant_tab.output[0],
            ('splic_junction_tab', 'TSV'): rules.gzip_star_sj_tab.output[0],
        }
        with open(output.manifest, 'w') as f:
            writer = csv.writer(f, dialect='excel-tab', lineterminator='\n')
            # Write column header
            cols = [
                '# case', 'disease', 'short_sample_type', 'aliquot',
                'result_type', 'result_path', 'result_format',
                'sample_name_fastq_R1', 'sample_uuid_fastq_R1',
                'sample_name_fastq_R2', 'sample_uuid_fastq_R2',
            ]
            writer.writerow(cols)

            for sample in filemap.samples:
                # Get sample info
                fq_pths = [filemap.fastq_map[sample][strand] for strand in ['R1', 'R2']]

                fq_catalogs = [filemap.pth_to_catalog[pth] for pth in fq_pths]
                fq_uuids = [c['UUID'] for c in fq_catalogs]
                fq_names = [c['# sample_name'] for c in fq_catalogs]
                case = fq_catalogs[0]['case']
                disease = fq_catalogs[0]['disease']
                sample_type = fq_catalogs[0]['short_sample_type']

                for (result_type, result_format), pth_fmt in result_file_tpls.items():
                    # Get absolute file path
                    abs_result_path = Path(pth_fmt.format(sample=sample)).resolve(strict=True)
                    writer.writerow([
                        case, disease, sample_type, sample,
                        result_type, str(abs_result_path), result_format,
                        fq_names[0], fq_uuids[0],
                        fq_names[1], fq_uuids[1],
                    ])


rule rsem_calc_expression:
    output:
        genes=temporary("rsem/{sample}.rsem.genes.results"),
        isoforms=temporary("rsem/{sample}.rsem.isoforms.results")
    input:
        bam=rules.star_align.output.quant_tx_bam
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: 36000 + 8000 * (attempt - 1),
        tmp_mb=32000
    params:
        rsem_ref_prefix=RSEM_REF_PREFIX
    log: 'logs/rsem/{sample}.log'
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
    shell: "gzip -9n -c {input} > {output}"


rule expression_all_samples:
    input: **expand_to_all_samples({ \
        'rsem_genes': rules.rsem_calc_expression.output.genes + '.gz', \
        'rsem_isoforms': rules.rsem_calc_expression.output.isoforms + '.gz', \
        'rnaseqc_exon_reads': rules.rnaseqc.output.exon_reads + '.gz', \
        'rnaseqc_gene_fragments': rules.rnaseqc.output.gene_fragments + '.gz', \
        'rnaseqc_gene_reads': rules.rnaseqc.output.gene_reads + '.gz', \
        'rnaseqc_gene_tpm': rules.rnaseqc.output.gene_tpm + '.gz', \
    })

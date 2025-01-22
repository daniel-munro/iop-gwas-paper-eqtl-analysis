def fastq_paths_1(wildcards):
    files = fastqs.loc[fastqs["sample_id"] == wildcards.sample_id]
    return ["data/fastq/{}1.fq.gz".format(x) for x in files["prefix"]]


def fastq_paths_2(wildcards):
    files = fastqs.loc[fastqs["sample_id"] == wildcards.sample_id]
    return ["data/fastq/{}2.fq.gz".format(x) for x in files["prefix"]]


rule star_index:
    input:
        fasta = "../data/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa",
        gtf = "../data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.gtf"
    output:
        # Among others:
        "data/star_index/SAindex"
    threads: 4
    shell:
        """
        STAR --runMode genomeGenerate \
        --genomeDir data/star_index \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 149 \
        --runThreadN 4
        """

rule star_align_1st_pass:
    # Get novel junctions to use for 2nd pass.
    input:
        fastq1 = fastq_paths_1,
        fastq2 = fastq_paths_2,
        index = "data/star_index/SAindex"
    params:
        fastq_list_1 = lambda wildcards, input: ",".join(input.fastq1),
        fastq_list_2 = lambda wildcards, input: ",".join(input.fastq2),
        prefix = "data/star_sj/{sample_id}."
    output:
        "data/star_sj/{sample_id}.SJ.out.tab"
    resources:
        mem_mb = 120000
    threads: 16
    shell:
        """
        STAR --runMode alignReads \
        --runThreadN 16 \
        --genomeDir data/star_index \
        --readFilesIn {params.fastq_list_1} {params.fastq_list_2} \
        --readFilesCommand zcat \
        --outSAMtype None \
        --outFileNamePrefix {params.prefix}
        """


def read_groups(wildcards, input):
    fastqs = [x.replace("_1.fq.gz", "") for x in input.fastq1]
    rgs = expand("ID:{fq} SM:{sam}", fq=fastqs, sam=wildcards.sample_id)
    return " , ".join(rgs)

rule star_align_2nd_pass:
    input:
        fastq1 = fastq_paths_1,
        fastq2 = fastq_paths_2,
        vcf = "data/genotype/individual/{sample_id}.vcf.gz",
        sj = "data/star_sj/{sample_id}.SJ.out.tab",
        index = "data/star_index/SAindex"
    params:
        fastq_list_1 = lambda wildcards, input: ",".join(input.fastq1),
        fastq_list_2 = lambda wildcards, input: ",".join(input.fastq2),
        read_groups = read_groups,
        prefix = "data/star_out/{sample_id}.",
        extra = expand("data/star_out/{{sample_id}}{ext}",
                       # ext=[".SJ.out.tab", ".Log*", "._STAR*"])
                       ext=[".SJ.out.tab", "._STAR*"])
    output:
        # phASER requires sorted by coord, indexed. RSEM requires transcriptome-aligned.
        coord = "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam",
        trans = "data/star_out/{sample_id}.Aligned.toTranscriptome.out.bam",
        chim = "data/star_out/{sample_id}.Chimeric.out.junction"
    resources:
        mem_mb = 120000
    threads: 16
    shell:
        """
        STAR --runMode alignReads \
        --runThreadN 16 \
        --genomeDir data/star_index \
        --readFilesIn {params.fastq_list_1} {params.fastq_list_2} \
        --readFilesCommand zcat \
        --varVCFfile <(zcat {input.vcf}) \
        --waspOutputMode SAMtag \
        --sjdbFileChrStartEnd {input.sj} \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions WithinBAM SoftClip \
        --chimMainSegmentMultNmax 1 \
        --quantMode TranscriptomeSAM \
        --outSAMattrRGline {params.read_groups} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outFileNamePrefix {params.prefix}
        rm -r {params.extra}
        """

# rule index_bam:
#     input:
#         "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam"
#     output:
#         "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam.bai"
#     group: "align"
#     shell:
#         "samtools index {input}"

rule mark_duplicates:
    input:
        "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam"
    output:
        # bam = "data/markdup_out/{sample_id}.Aligned.sortedByCoord.out.md.bam",
        bam = "data/markdup_out/{sample_id}.bam",
        metrics = "data/markdup_out/{sample_id}.marked_dup_metrics.txt"
    shell:
        """
        picard MarkDuplicates \
        I={input} \
        O={output.bam} \
        M={output.metrics} \
        ASSUME_SORT_ORDER=coordinate \
        PROGRAM_RECORD_ID=null \
        TMP_DIR=$PBSTMPDIR \
        MAX_RECORDS_IN_RAM=2000000
        """

rule retrieve_fastq:
    output:
        "data/fastq/{batch}/{file}"
    shell:
        "rsync -av pejlab:~/data/fastq/{wildcards.batch}/{wildcards.file} data/fastq/{wildcards.batch}/"

rule retrieve_vcf:
    output:
        "data/vcfs_for_star/{vcf}"
    shell:
        "rsync -av pejlab:~/{output} data/vcfs_for_star/"

rule retrieve_sj:
    output:
        "data/star_out_sj/{sj}"
    shell:
        "rsync -av pejlab:~/{output} data/star_out_sj/"

rule bam_to_fastq:
    input:
        lambda wildcards: input_bams[wildcards.sample_id]
    output:
        "data/fastq_for_star/{sample_id}.fastq.gz"
    shell:
        """
        mkfifo read_pipe_{wildcards.sample_id}
        gzip -1 -c < read_pipe_{wildcards.sample_id} > {output} &
        java -jar -Xmx8g tools/picard.jar SamToFastq \
        I={input} FASTQ=read_pipe_{wildcards.sample_id} \
        INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=SILENT
        rm read_pipe_{wildcards.sample_id}
        """

rule star_load_genome:
    output:
        "star_genome_loaded"
    shell:
        """
        STAR --genomeDir data/star_index --genomeLoad LoadAndExit
        touch {output}
        """


def star_unload_genome():
    if os.path.exists("star_genome_loaded"):
        shell("STAR --genomeDir data/star_index --genomeLoad Remove")
        os.remove("star_genome_loaded")


onsuccess:
    star_unload_genome()
onerror:
    star_unload_genome()

rule star_align_1st_pass:
    # Get novel junctions to use for 2nd pass.
    input:
        fastq = fastq_paths,
        index = "data/star_index/SAindex",
        genome = "star_genome_loaded"
    params:
        fastq_list = lambda wildcards, input: ",".join(input.fastq)
        # unsorted = "data/star_out_sj/{sample_id}.Aligned.out.bam",
        # trans = "data/star_out_sj/{sample_id}.Aligned.toTranscriptome.out.bam",
        # chim = "data/star_out_sj/{sample_id}.Chimeric.out.junction"
    output:
        # 2nd pass deletes this file, so output 1st pass in different dir:
        "data/star_out_sj/{sample_id}.SJ.out.tab"
    # threads: 1000       # One STAR job should run alone due to memory.
    shell:
        """
        STAR --runMode alignReads \
        --runThreadN 8 \
        --genomeDir data/star_index \
        --genomeLoad LoadAndKeep \
        --readFilesIn {params.fastq_list} \
        --readFilesCommand zcat \
        --outSAMtype None \
        --outFileNamePrefix data/star_out_sj/{wildcards.sample_id}.
        """

rule star_align_2nd_pass:
    input:
        fastq = fastq_paths,
        vcf = vcf_path,
        # sj = expand("data/star_out_sj/{sample_id}.SJ.out.tab", sample_id=sample_ids),
        sj = "data/star_out_sj/{sample_id}.SJ.out.tab",
        index = "data/star_index/SAindex"
    params:
        fastq_list = lambda wildcards, input: ",".join(input.fastq),
        unsorted = "data/star_out/{sample_id}.Aligned.out.bam",
        sjnew = "data/star_out/{sample_id}.SJ.out.tab"
    output:
        coord = "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam",
        coordi = "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam.bai",
        trans = "data/star_out/{sample_id}.Aligned.toTranscriptome.out.bam",
        chim = "data/star_out/{sample_id}.Chimeric.out.junction"
    threads: 1000       # One STAR job should run alone due to memory.
    shell:
        """
        STAR --runMode alignReads \
        --runThreadN 8 \
        --genomeDir data/star_index \
        --readFilesIn {params.fastq_list} \
        --readFilesCommand zcat \
        --varVCFfile <(zcat {input.vcf}) \
        --waspOutputMode SAMtag \
        --sjdbFileChrStartEnd {input.sj} \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions WithinBAM SoftClip \
        --chimMainSegmentMultNmax 1 \
        --quantMode TranscriptomeSAM \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outFileNamePrefix data/star_out/{wildcards.sample_id}.
        tools/samtools sort --threads 8 -o {output.coord} {params.unsorted}
        tools/samtools index {output.coord}
        rm {params.unsorted} {params.sjnew}
        """

rule star_align_2nd_pass:
    input:
        fastq = fastq_paths,
        vcf = vcf_path,
        # sj = expand("data/star_out_sj/{sample_id}.SJ.out.tab", sample_id=sample_ids),
        sj = "data/star_out_sj/{sample_id}.SJ.out.tab",
        index = "data/star_index/SAindex"
    output:
        coord = "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam",
        trans = "data/star_out/{sample_id}.Aligned.toTranscriptome.out.bam",
        chim = "data/star_out/{sample_id}.Chimeric.out.junction"
    shell:
        """
        rsync -av garibaldi:~/data/star_out/{wildcards.sample_id}.* data/star_out/
        """

# Need to get bamsync binary.  Source code included in gtex pipeline.
rule bamsync:
    # sync BAMs (optional; copy QC flags and read group IDs)
    input:
        lambda wildcards: input_bams[wildcards.sample_id],
        "/secure/dmunro/star_out/{sample_id}.Aligned.sortedByCoord.out.bam"
    output:
        "/secure/dmunro/star_out/{sample_id}.Aligned.sortedByCoord.out.patched.bam"
    shell:
docker run --rm -v /secure/dmunro:/secure/dmunro -t broadinstitute/gtex_rnaseq:V8 \
    /bin/bash -c "/src/run_bamsync.sh \
        $input_bam \
        /secure/dmunro/star_out/${sample_id}.Aligned.sortedByCoord.out.bam \
        /secure/dmunro/star_out/${sample_id}"

rule rnaseqc:
    input:
        bam = "data/markdup_out/{sample_id}.Aligned.sortedByCoord.out.md.bam",
        gtf = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    params:
        interm1 = "data/rnaseqc_out/{sample_id}.exon_reads.gct",
        interm2 = "data/rnaseqc_out/{sample_id}.gene_reads.gct",
        interm3 = "data/rnaseqc_out/{sample_id}.gene_tpm.gct"
    output:
        "data/rnaseqc_out/{sample_id}.exon_reads.gct.gz",
        "data/rnaseqc_out/{sample_id}.gene_reads.gct.gz",
        "data/rnaseqc_out/{sample_id}.gene_tpm.gct.gz",
        "data/rnaseqc_out/{sample_id}.metrics.tsv"
    shell:
        """
        tools/rnaseqc {input.gtf} {input.bam} data/rnaseqc_out \
        -s {wildcards.sample_id} --unpaired
        # Seems identical to gene_reads file:
        rm data/rnaseqc_out/{wildcards.sample_id}.gene_fragments.gct
        gzip {params.interm1} {params.interm2} {params.interm3}
        """

rule combine_gcts:
    input:
        expand("data/rnaseqc_out/{sample_id}.{{quant_type}}.gct.gz",
               sample_id=sample_ids)
    output:
        "data/{quant_type,\w+}.gct.gz"
    shell:
        "python3 src/combine_GCTs.py {input} {wildcards.quant_type} -o data/"


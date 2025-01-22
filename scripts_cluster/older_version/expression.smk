rule rsem:
    input:
        ref = "../data/rsem_reference/rsem_reference.transcripts.fa",
        bam = "data/star_out/{sample_id}.Aligned.toTranscriptome.out.bam"
    output:
        "data/rsem_out/{sample_id}.genes.results.gz"
    params:
        prefix = "data/rsem_out/{sample_id}"
    conda:
        "../envs/rsem.yaml"
    threads: 16
    shell:
        """
        rsem-calculate-expression \
        --paired-end \
        --num-threads 16 \
        --quiet \
        --estimate-rspd \
        --no-bam-output \
        --alignments {input.bam} \
        ../data/rsem_reference/rsem_reference \
        data/rsem_out/{wildcards.sample_id}
        gzip {params.prefix}.genes.results
        rm {params.prefix}.isoforms.results
        rm -r {params.prefix}.stat
        """

# rule gzip_rsem_output:
#     input:
#         "data/rsem_out/{sample_id}.genes.results"
#     output:
#         "data/rsem_out/{sample_id}.genes.results.gz"
#     group:
#         "rsem"
#     shell:
#         "gzip {input}"

rule combine_RSEM:
    input:
        expand("data/rsem_out/{sample_id}.genes.results.gz", sample_id=sample_ids)
    output:
        "data/expression/rsem_{field}.gct.gz"
    shell:
        "python3 ../src/combine_RSEM.py {input} {wildcards.field} {output}"

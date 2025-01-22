localrules:
    phaser_expr_matrix,

rule phaser:
    input:
        # bam = "data/markdup_out/{rat_id}_{region}.Aligned.sortedByCoord.out.md.bam",
        # bai = "data/markdup_out/{rat_id}_{region}.Aligned.sortedByCoord.out.md.bam.bai",
        bam = "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam",
        bai = "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam.bai",
        # bam = "data/markdup_out/{sample_id}.bam",
        # bai = "data/markdup_out/{sample_id}.bam.bai",
        vcf = "data/genotype/eyes.vcf.gz",
        vcfi = "data/genotype/eyes.vcf.gz.tbi"
    output:
        "data/phaser_out/{sample_id}.vcf.gz",
        "data/phaser_out/{sample_id}.vcf.gz.tbi",
        "data/phaser_out/{sample_id}.haplotypic_counts.txt"
    conda:
        "../envs/phaser.yaml"
    # group: "phaser"
    shell:
        # "python tools/phaser/wrapper.py id {input.bam} {input.vcf} {input.gene_models} output --paired-end"
        """
        python ../tools/phaser/phaser/phaser.py \
        --temp_dir $PBSTMPDIR \
        --bam {input.bam} \
        --vcf {input.vcf} \
        --sample {wildcards.sample_id} \
        --baseq 10 \
        --mapq 255 \
        --isize 1e6 \
        --paired_end 0 \
        --o data/phaser_out/{wildcards.sample_id} \
        --include_indels 0 \
        --gw_phase_vcf 1 \
        --threads 8
        """

rule phaser_gene_ae:
    input:
        counts = "data/phaser_out/{sample_id}.haplotypic_counts.txt",
        gene_models = "../data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.bed"
    output:
        "data/phaser_gene_ae_out/{sample_id}.gene_ae.txt"
    conda:
        "../envs/phaser.yaml"
    # group: "phaser"
    shell:
        # --min_haplo_maf 0.05 \
        """
        python ../tools/phaser/phaser_gene_ae/phaser_gene_ae.py \
        --haplotypic_counts {input.counts} \
        --features {input.gene_models} \
        --o {output}
        """

rule phaser_expr_matrix:
    """Note: Due to phaser bug involving bed column order, internal tabix command doesn't work, but index isn't necessary. So I ignore the error and delete the bad index files."""
    input:
        gene_ae = expand("data/phaser_gene_ae_out/{sample_id}.gene_ae.txt",
                         sample_id=sample_ids),
        gene_models = "../data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.bed"
    output:
        bed = "data/phaser_pop_out/expr_matrix.bed.gz",
        # bedi = "data/phaser_pop_out/expr_matrix.bed.gz.tbi",
        bedgw = "data/phaser_pop_out/expr_matrix.gw_phased.bed.gz",
        # bedgwi = "data/phaser_pop_out/expr_matrix.gw_phased.bed.gz.tbi"
    conda:
        "../envs/phaser.yaml"
    shell:
        # tabix {output.bed}
        # tabix {output.bedgw}
        """
        python ../tools/phaser/phaser_pop/phaser_expr_matrix.py \
        --gene_ae_dir data/phaser_gene_ae_out \
        --features {input.gene_models} \
        --o data/phaser_pop_out/expr_matrix || true
        rm {output.bed}.tbi {output.bedgw}.tbi
        """

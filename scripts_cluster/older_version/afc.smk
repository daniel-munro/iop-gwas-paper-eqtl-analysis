localrules:
    phaser_sample_map,
    phaser_gene_var_pairs,

def afc_expr_bed(wildcards):
    return {
        "basic": "data/expression/rsem_expected_count.bed.gz",
        "basic3": "data/expression/rsem_TPM.bed.gz",
        "basic4": "data/expression/rsem_expected_count.bed.gz",
        "main": "data/expression/rsem_expected_count.bed.gz",
        "main3": "data/expression/rsem_TPM.bed.gz",
        "main4": "data/expression/rsem_expected_count.bed.gz",
        "main5": "data/expression/rsem_expected_count.bed.gz",
        "qtl2":  "data/expression/rsem_expected_count.bed.gz",
    }[wildcards.method]

def afc_covar_file(wildcards):
    return {
        "basic": [],
        "basic3": [],
        "basic4": [],
        "main": "data/tensorqtl/main.combined_covariates.txt",
        "main3": "data/tensorqtl/main3.combined_covariates.txt",
        "main4": "data/tensorqtl/main4.combined_covariates.txt",
        "main4": "data/tensorqtl/main4.combined_covariates.txt",
        "main5": "data/tensorqtl/main4.PEER_covariates.txt",
        "qtl2": [],
    }[wildcards.method]

def afc_qtl_file(wildcards):
    return {
        "basic": "data/tensorqtl/basic.cis_qtl.txt.gz",
        "basic3": "data/tensorqtl/basic3.cis_qtl.txt.gz",
        "basic4": "data/tensorqtl/basic4.cis_qtl.txt.gz",
        "main": "data/tensorqtl/main.cis_qtl.txt.gz",
        "main3": "data/tensorqtl/main3.cis_qtl.txt.gz",
        "main4": "data/tensorqtl/main4.cis_qtl.txt.gz",
        # "main5": "data/tensorqtl/main5.cis_qtl.txt.gz",
        "main5": "data/tensorqtl/main5.cis_qtl_signif.txt.gz",
        "qtl2": "data/qtl2/eyes.gene_var_pval.tsv.gz",
    }[wildcards.method]

def aFC_qtl_arg(wildcards, input):
    if wildcards.method == "qtl2":
        qtl_call = "<(python3 ../src/prepare_qtl_for_afc.py {qtl} 1 2)"
    else:
        # qtl_call = "<(python3 ../src/prepare_qtl_for_afc.py {qtl})"
        ## Now I'm getting aFC for all significant eQTLs:
        qtl_call = "<(python3 ../src/prepare_qtl_for_afc.py {qtl} 1 2)"
    return qtl_call.format(qtl=input.qtl)


def aFC_covar_arg(wildcards, input):
    if wildcards.method in {"basic", "basic2", "basic3", "basic4", "qtl2"}:
        return ""
    else:
        return "--cov <(sed 's/-N//g' {covar})".format(covar=input.covar)


rule aFC_from_eQTL_model:
    input:
        vcf = "data/genotype/eyes.vcf.gz",
        vcfi = "data/genotype/eyes.vcf.gz.tbi",
        bed = afc_expr_bed,
        bedi = lambda w: afc_expr_bed(w) + ".tbi",
        qtl = afc_qtl_file,
        covar = afc_covar_file,
    output:
        "data/afc/{method}.aFC.txt"
    params:
        qtl_arg = aFC_qtl_arg,
        covar_arg = aFC_covar_arg
    conda:
        "../envs/afc.yaml"
    resources:
        walltime = 12
    shell:
        """
        python3 ../tools/aFC/aFC.py \
        --vcf {input.vcf} \
        --pheno {input.bed} \
        --qtl {params.qtl_arg} \
        {params.covar_arg} \
        --log_xform 1 \
        --output {output}
        """

# rule aFC_from_eQTL_model_main:
#     input:
#         vcf = "data/genotype/eyes.vcf.gz",
#         vcfi = "data/genotype/eyes.vcf.gz.tbi",
#         bed = "data/expression/rsem_expected_count.bed.gz",
#         bedi = "data/expression/rsem_expected_count.bed.gz.tbi",
#         qtl = "data/tensorqtl/main.cis_qtl.txt.gz",
#         covar = "data/tensorqtl/main.combined_covariates.txt",
#     output:
#         "data/afc/main.aFC.txt"
#     conda:
#         "../envs/afc.yaml"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python3 ../tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 ../src/prepare_qtl_for_afc.py {input.qtl}) \
#         --cov <(sed 's/-N//g' {input.covar}) \
#         --log_xform 1 \
#         --output {output}
#         """

# rule aFC_from_eQTL_model_main3:
#     input:
#         vcf = "data/genotype/eyes.vcf.gz",
#         vcfi = "data/genotype/eyes.vcf.gz.tbi",
#         bed = "data/expression/rsem_TPM.bed.gz",
#         bedi = "data/expression/rsem_TPM.bed.gz.tbi",
#         qtl = "data/tensorqtl/main3.cis_qtl.txt.gz",
#         covar = "data/tensorqtl/main3.combined_covariates.txt",
#     output:
#         "data/afc/main3.aFC.txt"
#     conda:
#         "../envs/afc.yaml"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python3 ../tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 ../src/prepare_qtl_for_afc.py {input.qtl}) \
#         --cov <(sed 's/-N//g' {input.covar}) \
#         --log_xform 1 \
#         --output {output}
#         """

# rule aFC_from_eQTL_model_basic:
#     input:
#         vcf = "data/genotype/eyes.vcf.gz",
#         vcfi = "data/genotype/eyes.vcf.gz.tbi",
#         bed = "data/expression/rsem_expected_count.bed.gz",
#         bedi = "data/expression/rsem_expected_count.bed.gz.tbi",
#         qtl = "data/tensorqtl/basic.cis_qtl.txt.gz",
#     output:
#         "data/afc/basic.aFC.txt"
#     conda:
#         "../envs/afc.yaml"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python3 ../tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 ../src/prepare_qtl_for_afc.py {input.qtl}) \
#         --log_xform 1 \
#         --output {output}
#         """

# rule aFC_from_eQTL_model_basic3:
#     input:
#         vcf = "data/genotype/eyes.vcf.gz",
#         vcfi = "data/genotype/eyes.vcf.gz.tbi",
#         bed = "data/expression/rsem_TPM.bed.gz",
#         bedi = "data/expression/rsem_TPM.bed.gz.tbi",
#         qtl = "data/tensorqtl/basic3.cis_qtl.txt.gz",
#     output:
#         "data/afc/basic3.aFC.txt"
#     conda:
#         "../envs/afc.yaml"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python3 ../tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 ../src/prepare_qtl_for_afc.py {input.qtl}) \
#         --log_xform 1 \
#         --output {output}
#         """

# rule aFC_from_eQTL_model_qtl2:
#     input:
#         vcf = "data/genotype/eyes.vcf.gz",
#         vcfi = "data/genotype/eyes.vcf.gz",
#         bed = "data/expression/rsem_expected_count.bed.gz",
#         bedi = "data/expression/rsem_expected_count.bed.gz.tbi",
#         qtl = "data/qtl2/eyes.gene_var_pval.tsv.gz",
#     output:
#         "data/afc/qtl2.aFC.txt"
#     conda:
#         "../envs/afc.yaml"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python3 ../tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 ../src/prepare_qtl_for_afc.py {input.qtl} 1 2) \
#         --log_xform 1 \
#         --output {output}
#         """

###################
## ASE-based aFC ##
###################

rule phaser_sample_map:
    output:
        "data/afc/sample_map.txt"
    run:
        with open(output[0], "w") as out:
            out.write("vcf_sample\tbed_sample\n")
            for rat_id in rat_ids:
                bam_base = "{}.Aligned.sortedByCoord.out".format(rat_id)
                out.write("{}\t{}\n".format(rat_id, bam_base))
                # out.write("{}\t{}\n".format(rat_id, rat_id))

rule phaser_gene_var_pairs:
    input:
        "data/afc/{version}.aFC.txt"
    output:
        "data/afc/{version}.pairs.txt"
    run:
        inlines = open(input[0], "r").read().splitlines()
        gene_var = [x.split("\t")[:2] for x in inlines[1:]]
        with open(output[0], "w") as out:
            out.write("gene_id\tvar_id\tvar_contig\tvar_pos\tvar_ref\tvar_alt\n")
            for x in gene_var:
                chrom, pos = tuple(x[1].split(":"))
                chrom = chrom.replace("chr", "")
                out.write("{}\t{}\t{}\t{}\t\t\n".format(x[0], x[1], chrom, pos))

rule phaser_cis_var:
    input:
        bed = "data/phaser_pop_out/expr_matrix.gw_phased.bed.gz",
        # bedi = "data/expression/expr_matrix.gw_phased.bed.gz.tbi",
        vcf = "data/genotype/eyes.vcf.gz",
        gene_var_pairs = "data/afc/{method}.pairs.txt",
        sample_map = "data/afc/sample_map.txt"
    output:
        "data/afc/{method}.ASE_aFC.txt"
    conda:
        "../envs/phaser.yaml"
    resources:
        walltime = 12
    threads: 16
    shell:
        """
        python ../tools/phaser/phaser_pop/phaser_cis_var.py \
        --bed {input.bed} \
        --vcf {input.vcf} \
        --pair {input.gene_var_pairs} \
        --map {input.sample_map} \
        --t 16 \
        --o {output}
        """

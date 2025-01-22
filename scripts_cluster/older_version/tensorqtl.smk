localrules:
    vcf_to_plink,
    sample_participant_lookup,
    prepare_expression,
    # expression_tsv_to_bed,
    expression_gct_to_bed,
    # other_covariates,
    combine_covariates,
    empty_covariates,
    tensorqtl_all_signif,

rule vcf_to_plink:
    input:
        "data/genotype/eyes.vcf.gz"
    output:
        multiext("data/genotype/eyes", ".bed", ".bim", ".fam")
    params:
        prefix = "data/genotype/eyes"
    shell:
        """
        plink2 --make-bed \
        --vcf {input} \
        --maf 0.0001 \
        --max-maf 0.9999 \
        --out {params.prefix}
        """

rule sample_participant_lookup:
    input: "data/rat_ids.txt"
    output: "data/expression/samples_participants.txt"
    run:
        with open(output[0], "w") as out:
            out.write("sample_id\tparticipant_id\n")
            for rat_id in rat_ids:
                out.write("{0}\t{1}\n".format(rat_id, rat_id))

rule prepare_expression:
    input:
        tpm_gct = "data/expression/rsem_TPM.gct.gz",
        counts_gct = "data/expression/rsem_expected_count.gct.gz",
        gtf = "../data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf",
        samples = "data/expression/samples_participants.txt",
        chrlist = "data/genotype/eyes.chrlist.txt",
        src = "../src/rnaseqnorm.py"
    output:
        "data/tensorqtl/eyes.expression.bed.gz",
        "data/tensorqtl/eyes.expression.bed.gz.tbi"
    params:
        prefix = "data/tensorqtl/eyes"
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python3 ../src/eqtl_prepare_expression.py \
        {input.tpm_gct} {input.counts_gct} {input.gtf} \
        {input.samples} {input.chrlist} {params.prefix} \
        --tpm_threshold 0.1 \
        --count_threshold 6 \
        --sample_frac_threshold 0.2 \
        --normalization_method tmm
        """

rule expression_gct_to_bed:
    input:
        gct = "data/expression/rsem_{field}.gct.gz",
        gtf = "../data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    output:
        bed = "data/expression/rsem_{field}.bed.gz",
        bedi = "data/expression/rsem_{field}.bed.gz.tbi"
    params:
        interm = "data/expression/rsem_{field}.bed",
    conda:
        "../envs/bioinfo.yaml"
    shell:
        """
        python3 ../src/gct_to_bed.py {input.gct} {input.gtf} {params.interm}
        bgzip {params.interm}
        tabix {output.bed}
        """

def expr_bed(wildcards):
    return {
        "basic": "data/expression/rsem_expected_count.bed.gz",
        "basic3": "data/expression/rsem_TPM.bed.gz",
        "basic4": "data/tensorqtl/eyes.expression.bed.gz",
        "main": "data/expression/rsem_expected_count.bed.gz",
        "main3": "data/expression/rsem_TPM.bed.gz",
        "main4": "data/tensorqtl/eyes.expression.bed.gz",
        "main5": "data/tensorqtl/eyes.expression.bed.gz",
    }[wildcards.method]

def covar_file(wildcards):
    return {
        "basic": "data/tensorqtl/basic.covar_empty.txt",
        "basic3": "data/tensorqtl/basic3.covar_empty.txt",
        "basic4": "data/tensorqtl/basic4.covar_empty.txt",
        "main": "data/tensorqtl/main.combined_covariates.txt",
        "main3": "data/tensorqtl/main3.combined_covariates.txt",
        "main4": "data/tensorqtl/main4.combined_covariates.txt",
        "main5": "data/tensorqtl/main4.PEER_covariates.txt",
    }[wildcards.method]

rule run_peer:
    input:
        # "data/tensorqtl/eyes.expression.bed.gz"
        # "data/expression/rsem_expected_count.bed.gz"
        expr_bed
    output:
        "data/tensorqtl/{method}.PEER_residuals.txt",
        "data/tensorqtl/{method}.PEER_alpha.txt",
        "data/tensorqtl/{method}.PEER_covariates.txt"
    params:
        prefix = "data/tensorqtl/{method}",
        num_peer = 15
    conda:
        "../envs/peer.yaml"
    resources:
        walltime = 6
    shell:
        "Rscript ../src/run_PEER.R {input} {params.prefix} {params.num_peer}"

rule similarity_to_founders:
    input:
        pop_vcf = "data/genotype/eyes.vcf.gz",
        pop_vcfi = "data/genotype/eyes.vcf.gz.tbi",
        founder_vcf = "../data/genotype/founders.vcf.gz",
        founder_vcfi = "../data/genotype/founders.vcf.gz.tbi"
    output:
        "data/tensorqtl/sim_to_founders.txt"
    conda:
        "../envs/biopython.yaml"
    shell:
        "python3 ../src/sim_to_each_founder.py {input.pop_vcf} {input.founder_vcf} {output}"

rule combine_covariates:
    input:
        peer = "data/tensorqtl/{method}.PEER_covariates.txt",
        founder_sim = "data/tensorqtl/sim_to_founders.txt",
        # other = "data/tensorqtl/other_covariates.tsv"
    output:
        "data/tensorqtl/{method}.combined_covariates.txt"
    params:
        prefix = "data/tensorqtl/{method}"
    shell:
        """
        python3 ../src/combine_covariates.py {input.peer} {params.prefix} \
        --add_covariates {input.founder_sim}
        """

rule empty_covariates:
    """tensorQTL currently requires a covariates file."""
    input:
        # "data/tensorqtl/main.expression.bed.gz"
        # "data/expression/rsem_expected_count.bed.gz"
        expr_bed
    output:
        "data/tensorqtl/{method}.covar_empty.txt"
    shell:
        "zcat {input} | head -1 | sed 's/gene_id/ID/' | cut -f4- > {output} || true"

rule tensorqtl_perm:
    input:
        geno = multiext("data/genotype/eyes", ".bed", ".bim", ".fam"),
        bed = expr_bed,
        bedi = lambda w: expr_bed(w) + ".tbi",
        covar = covar_file
    output:
        "data/tensorqtl/{method}.cis_qtl.txt.gz"
    params:
        geno_prefix = "data/genotype/eyes",
        window = 500000
    conda:
        "../envs/tensorqtl.yaml"
    resources:
        walltime = 12,
        # gpu = ",nodes=1:ppn=8:gtx1080ti -q gpu"
    shell:
        """
        # module load cuda
        # pip install -e ../tools/tensorqtl
        python -m tensorqtl --mode cis \
            {params.geno_prefix} \
            {input.bed} \
            {wildcards.method} \
            --covariates {input.covar} \
            --window {params.window} \
            --output_dir data/tensorqtl
        """

rule tensorqtl_nominal:
    input:
        geno = multiext("data/genotype/eyes", ".bed", ".bim", ".fam"),
        bed = expr_bed,
        bedi = lambda w: expr_bed(w) + ".tbi",
        covar = covar_file
    output:
        expand("data/tensorqtl/{{method}}/{{method}}.cis_qtl_pairs.{chrn}.parquet",
               chrn=range(1, 21))
    params:
        geno_prefix = "data/genotype/eyes",
        outdir = "data/tensorqtl/{method}",
        window = 500000
    conda:
        "../envs/tensorqtl.yaml"
    resources:
        walltime = 12,
        # gpu = ",nodes=1:ppn=8:gtx1080ti -q gpu"
    shell:
        """
        # module load cuda
        # pip install -e ../tools/tensorqtl
        mkdir -p {params.outdir}
        python -m tensorqtl --mode cis_nominal \
            {params.geno_prefix} \
            {input.bed} \
            {wildcards.method} \
            --covariates {input.covar} \
            --window {params.window} \
            --output_dir {params.outdir}
        """

rule tensorqtl_all_signif:
    input:
        perm = "data/tensorqtl/{code}.cis_qtl.txt.gz",
        nom = expand("data/tensorqtl/{{code}}/{{code}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    output:
        "data/tensorqtl/{code}.cis_qtl_signif.txt.gz"
    params:
        nom_prefix = "data/tensorqtl/{code}/{code}"
    conda:
        "../envs/tensorqtl.yaml"
    shell:
        "python3 ../src/tensorqtl_all_sig_eqtls.py {input.perm} {params.nom_prefix} {output}"


# rule run_tensorqtl_perm_main:
#     input:
#         geno = multiext("data/genotype/eyes", ".bed", ".bim", ".fam"),
#         bed = "data/expression/rsem_expected_count.bed.gz",
#         bedi = "data/expression/rsem_expected_count.bed.gz.tbi",
#         covar = "data/tensorqtl/main.combined_covariates.txt",
#     output:
#         "data/tensorqtl/main.cis_qtl.txt.gz"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python -m tensorqtl --mode cis \
#         data/genotype/eyes \
#         {input.bed} \
#         main \
#         --covariates {input.covar} \
#         --output_dir data/tensorqtl
#         """

# rule run_tensorqtl_perm_main3:
#     input:
#         geno = multiext("data/genotype/eyes", ".bed", ".bim", ".fam"),
#         bed = "data/expression/rsem_TPM.bed.gz",
#         bedi = "data/expression/rsem_TPM.bed.gz.tbi",
#         covar = "data/tensorqtl/main3.combined_covariates.txt",
#     output:
#         "data/tensorqtl/main3.cis_qtl.txt.gz"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python -m tensorqtl --mode cis \
#         data/genotype/eyes \
#         {input.bed} \
#         main3 \
#         --covariates {input.covar} \
#         --output_dir data/tensorqtl
#         """

# rule run_tensorqtl_perm_basic:
#     input:
#         geno = multiext("data/genotype/eyes", ".bed", ".bim", ".fam"),
#         bed = "data/expression/rsem_expected_count.bed.gz",
#         bedi = "data/expression/rsem_expected_count.bed.gz.tbi",
#         covar = "data/tensorqtl/basic.covar_empty.txt",
#     output:
#         "data/tensorqtl/basic.cis_qtl.txt.gz"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python -m tensorqtl --mode cis \
#         data/genotype/eyes \
#         {input.bed} \
#         basic \
#         --covariates {input.covar} \
#         --output_dir data/tensorqtl
#         """

# rule run_tensorqtl_perm_basic3:
#     input:
#         geno = multiext("data/genotype/eyes", ".bed", ".bim", ".fam"),
#         bed = "data/expression/rsem_TPM.bed.gz",
#         bedi = "data/expression/rsem_TPM.bed.gz.tbi",
#         covar = "data/tensorqtl/basic3.covar_empty.txt",
#     output:
#         "data/tensorqtl/basic3.cis_qtl.txt.gz"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python -m tensorqtl --mode cis \
#         data/genotype/eyes \
#         {input.bed} \
#         basic3 \
#         --covariates {input.covar} \
#         --output_dir data/tensorqtl
#         """

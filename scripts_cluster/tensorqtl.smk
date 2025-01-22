localrules:
    vcf_to_plink,


rule vcf_to_plink:
    """Get SNPs that are not monomorphic in a given set of samples."""
    input:
        vcf = "data/genotype/eyes.vcf.gz",
        samples = "data/rat_ids.txt",
    output:
        multiext("data/genotype/geno", ".bed", ".bim", ".fam")
    params:
        prefix = "data/genotype/geno"
    shell:
        """
        plink2 --make-bed \
        --vcf {input.vcf} \
        --keep {input.samples} \
        --maf 0.01 \
        --mac 2 \
        --max-alleles 2 \
        --out {params.prefix}
        """


rule tensorqtl_perm:
    """Map cis-eQTLs, determining significance using permutations.
    Outputs the top association per gene.
    """
    input:
        geno = multiext("data/genotype/geno", ".bed", ".bim", ".fam"),
        bed = "data/Eye.expr.iqn.filtered.bed.gz",
        bedi = "data/Eye.expr.iqn.filtered.bed.gz.tbi",
        covar = "data/covar.txt",
    output:
        "data/Eye.cis_qtl.txt.gz"
    params:
        geno_prefix = "data/genotype/geno",
        window = window,
    resources:
        walltime = 12
    shell:
        """
        python3 src/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --window {params.window} \
            --mode cis
        """


rule tensorqtl_independent:
    """Use stepwise regression to identify multiple conditionally independent cis-eQTLs per gene."""
    input:
        geno = multiext("data/genotype/geno", ".bed", ".bim", ".fam"),
        bed = "data/Eye.expr.iqn.filtered.bed.gz",
        bedi = "data/Eye.expr.iqn.filtered.bed.gz.tbi",
        covar = "data/covar.txt",
        cis = "data/Eye.cis_qtl.txt.gz"
    output:
        "data/Eye.cis_independent_qtl.txt.gz"
    params:
        geno_prefix = "data/genotype/geno",
        window = window,
    resources:
        walltime = 12
    shell:
        """
        python3 src/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --cis_output {input.cis} \
            --window {params.window} \
            --mode cis_independent
        """


rule tensorqtl_nominal:
    """Get summary statistics for all tested cis-window SNPs per gene."""
    input:
        geno = multiext("data/genotype/geno", ".bed", ".bim", ".fam"),
        bed = "data/Eye.expr.iqn.filtered.bed.gz",
        bedi = "data/Eye.expr.iqn.filtered.bed.gz.tbi",
        covar = "data/covar.txt",
    output:
        expand("data/nominal/Eye.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    params:
        geno_prefix = "data/genotype/geno",
        outdir = "data/nominal",
        window = window,
    resources:
        walltime = 12,
    shell:
        """
        mkdir -p {params.outdir}
        python3 -m tensorqtl \
            {params.geno_prefix} \
            {input.bed} \
            Eye \
            --covariates {input.covar} \
            --window {params.window} \
            --output_dir {params.outdir} \
            --mode cis_nominal
        """


rule tensorqtl_all_signif:
    """Extract all significant cis SNP-gene pairs."""
    input:
        perm = "data/Eye.cis_qtl.txt.gz",
        nom = expand("data/nominal/Eye.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21)),
    output:
        "data/Eye.cis_qtl_signif.txt.gz"
    params:
        nom_prefix = "data/nominal/Eye"
    shell:
        "python3 src/tensorqtl_all_signif.py {input.perm} {params.nom_prefix} {output}"

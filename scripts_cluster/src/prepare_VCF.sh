# I am extracting the Eye samples from the RatGTEx combined VCF, which has
# corrected header info and multiallelic variants removed.

bcftools view \
    ~/ratgtex/geno/ratgtex.vcf.gz \
    --samples-file data/rat_ids.txt \
    --exclude-uncalled \
    -Oz -o data/genotype/eyes.vcf.gz

# ── Rules: ggtyper profiling and genotyping ──────────────────────────────────


rule profile_samples:
    """Build insert-size profiles for all BAM files."""
    input:
        bam_list=f"{OUTDIR}/bam_list.txt"
    output:
        profiles_list=temp(f"{OUTDIR}/profiles/sampleProfiles.txt")
    params:
        outdir=f"{OUTDIR}/profiles",
        threads=THREADS
    shell:
        """
        mkdir -p {params.outdir}
        {GGTYPER} profile-samples \
          {input.bam_list} \
          {output.profiles_list} \
          {params.outdir}/ \
          -T {params.threads} -f
        """


rule profile_variants:
    """Build variant profiles using library parameters from all sample profiles."""
    input:
        variants=f"{OUTDIR}/merged/variants.json",
        sample_profiles=f"{OUTDIR}/profiles/sampleProfiles.txt"
    output:
        profiles_list=temp(f"{OUTDIR}/profiles/variantProfiles.txt")
    params:
        outdir=f"{OUTDIR}/profiles",
        threads=THREADS
    shell:
        """
        {GGTYPER} profile-variants \
          {input.variants} \
          {output.profiles_list} \
          {params.outdir}/ \
          {input.sample_profiles} \
          -T {params.threads}
        """


rule genotype:
    """Genotype all variants across all samples."""
    input:
        variant_profiles=f"{OUTDIR}/profiles/variantProfiles.txt",
        sample_profiles=f"{OUTDIR}/profiles/sampleProfiles.txt"
    output:
        results=f"{OUTDIR}/out_genotype_results.tsv.gz"
    params:
        prefix=f"{OUTDIR}/out",
        threads=THREADS
    shell:
        """
        {GGTYPER} genotype \
          {input.variant_profiles} \
          {input.sample_profiles} \
          {params.prefix} \
          -T {params.threads} -e
        gzip {params.prefix}_genotype_results.tsv
        """


rule annotate_vcf_with_ggtyper:
    """Add GGTyper genotype metrics to the merged VCF as GGT_* INFO/FORMAT fields."""
    input:
        vcf=f"{OUTDIR}/merged/population_sv.vcf",
        results=f"{OUTDIR}/out_genotype_results.tsv.gz"
    output:
        vcf=f"{OUTDIR}/out_ggtyper_annotated.vcf.gz",
        tbi=f"{OUTDIR}/out_ggtyper_annotated.vcf.gz.tbi"
    params:
        script=f"{SCRIPTS}/annotate_vcf_with_ggtyper.py",
        tmp_vcf=f"{OUTDIR}/out_ggtyper_annotated.tmp.vcf"
    conda:
        "../envs/genotyping.yaml"
    shell:
        """
        python3 {params.script} \
          {input.vcf} \
          {input.results} \
          {params.tmp_vcf}
        bcftools view {params.tmp_vcf} -O z -o {output.vcf}
        bcftools index -t {output.vcf}
        rm {params.tmp_vcf}
        """

# ── Rules: ggtyper profiling and genotyping ──────────────────────────────────


rule profile_samples:
    """Build insert-size profiles for all BAM files."""
    input:
        bam_list=f"{OUTDIR}/bam_list.txt"
    output:
        profiles_list=f"{OUTDIR}/profiles/sampleProfiles.txt"
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
          -T {params.threads}
        """


rule profile_variants:
    """Build variant profiles using library parameters from all sample profiles."""
    input:
        variants=f"{OUTDIR}/merged/variants.json",
        sample_profiles=f"{OUTDIR}/profiles/sampleProfiles.txt"
    output:
        profiles_list=f"{OUTDIR}/profiles/variantProfiles.txt"
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

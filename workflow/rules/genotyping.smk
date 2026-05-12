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


rule plot_ggtyper_certainty:
    """Plot GGTyper certainty distributions to help choose post-genotyping thresholds."""
    input:
        results=f"{OUTDIR}/out_genotype_results.tsv.gz"
    output:
        all=f"{OUTDIR}/out_genotype_results.certainty_all.svg",
        variant=f"{OUTDIR}/out_genotype_results.certainty_variant.svg",
        variant_quality_pass=f"{OUTDIR}/out_genotype_results.certainty_variant-quality-pass.svg",
        summary=f"{OUTDIR}/out_genotype_results.certainty_summary.tsv",
        thresholds=f"{OUTDIR}/out_genotype_results.certainty_summary.thresholds.tsv"
    params:
        script=f"{SCRIPTS}/plot_ggtyper_certainty.py",
        prefix=f"{OUTDIR}/out_genotype_results"
    conda:
        "../envs/genotyping.yaml"
    shell:
        """
        python3 {params.script} {input.results} --output-prefix {params.prefix} --mode all
        python3 {params.script} {input.results} --output-prefix {params.prefix} --mode variant
        python3 {params.script} {input.results} --output-prefix {params.prefix} --mode variant-quality-pass
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


rule filter_ggtyper_vcf:
    """Filter the annotated VCF using configurable GGTyper FORMAT thresholds."""
    input:
        vcf=f"{OUTDIR}/out_ggtyper_annotated.vcf.gz"
    output:
        vcf=f"{OUTDIR}/out_ggtyper_filtered.vcf.gz",
        tbi=f"{OUTDIR}/out_ggtyper_filtered.vcf.gz.tbi"
    params:
        script=f"{SCRIPTS}/filter_ggtyper_vcf.py",
        tmp_vcf=f"{OUTDIR}/out_ggtyper_filtered.tmp.vcf",
        require_quality_pass="--require-quality-pass" if GGTYPER_FILTER.get("require_quality_pass", True) else "--no-require-quality-pass",
        min_certainty=GGTYPER_FILTER.get("min_certainty", 0.8),
        min_genotype_quality=GGTYPER_FILTER.get("min_genotype_quality", 20),
        min_total_reads=GGTYPER_FILTER.get("min_total_reads", 10),
        min_avg_mapq=GGTYPER_FILTER.get("min_avg_mapq", 20),
        min_passing_samples=GGTYPER_FILTER.get("min_passing_samples", 1)
    conda:
        "../envs/genotyping.yaml"
    shell:
        """
        python3 {params.script} \
          {input.vcf} \
          {params.tmp_vcf} \
          {params.require_quality_pass} \
          --min-certainty {params.min_certainty} \
          --min-genotype-quality {params.min_genotype_quality} \
          --min-total-reads {params.min_total_reads} \
          --min-avg-mapq {params.min_avg_mapq} \
          --min-passing-samples {params.min_passing_samples}
        bcftools view {params.tmp_vcf} -O z -o {output.vcf}
        bcftools index -t {output.vcf}
        rm {params.tmp_vcf}
        """

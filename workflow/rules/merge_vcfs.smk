# ── Rules: per-sample VCF extraction and population-level merge ─────────────


rule extract_caller_vcf:
    """Extract calls from the configured caller for one sample."""
    input:
        vcf=lambda wc: VCF_FOR[wc.sample]
    output:
        vcf=f"{OUTDIR}/per_sample/{{sample}}_caller.vcf.gz",
        tbi=f"{OUTDIR}/per_sample/{{sample}}_caller.vcf.gz.tbi"
    params:
        caller=CALLER
    conda:
        "../envs/genotyping.yaml"
    shell:
        """
        bcftools view -i 'INFO/svdb_origin="{params.caller}"' \
          {input.vcf} -O z -o {output.vcf}
        bcftools index -t {output.vcf}
        """


rule merge_population_vcf:
    """Merge per-sample VCFs into one non-redundant population variant set."""
    input:
        vcfs=expand(f"{OUTDIR}/per_sample/{{sample}}_caller.vcf.gz", sample=SAMPLES)
    output:
        vcf=f"{OUTDIR}/merged/population_sv.vcf"
    params:
        vcf_list=lambda wc, input: " ".join(input.vcfs)
    conda:
        "../envs/genotyping.yaml"
    shell:
        """
        mkdir -p {OUTDIR}/merged
        svdb --merge --pass_only --vcf {params.vcf_list} > {output.vcf}
        """

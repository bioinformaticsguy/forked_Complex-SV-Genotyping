# ── Rules: per-sample VCF extraction and population-level merge ─────────────


rule extract_caller_vcf:
    """Extract calls from the configured caller for one sample."""
    input:
        vcf=lambda wc: VCF_FOR[wc.sample]
    output:
        vcf=temp(f"{OUTDIR}/per_sample/{{sample}}_caller.vcf.gz"),
        tbi=temp(f"{OUTDIR}/per_sample/{{sample}}_caller.vcf.gz.tbi")
    params:
        caller=CALLER
    conda:
        "../envs/genotyping.yaml"
    shell:
        """
        if bcftools view -h {input.vcf} | grep -q '##INFO=<ID=svdb_origin'; then
            bcftools view -i 'INFO/svdb_origin="{params.caller}"' \
              {input.vcf} -O z -o {output.vcf}
        else
            bcftools view {input.vcf} -O z -o {output.vcf}
        fi
        bcftools index -t {output.vcf}
        """


rule merge_population_vcf:
    """Merge per-sample VCFs into one non-redundant population variant set."""
    input:
        vcfs=expand(f"{OUTDIR}/per_sample/{{sample}}_caller.vcf.gz", sample=SAMPLES)
    output:
        vcf=temp(f"{OUTDIR}/merged/population_sv.vcf")
    params:
        vcf_list=lambda wc, input: " ".join(input.vcfs)
    conda:
        "../envs/genotyping.yaml"
    shell:
        """
        mkdir -p {OUTDIR}/merged
        svdb --merge --pass_only --vcf {params.vcf_list} > {output.vcf}
        """

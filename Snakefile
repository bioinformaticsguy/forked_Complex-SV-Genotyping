import pandas as pd

configfile: "configs/config.yaml"

samples_df = pd.read_csv(config["samples_sheet"], sep="\t")
SAMPLES    = samples_df["sample_id"].tolist()
BAM_FOR    = dict(zip(samples_df["sample_id"], samples_df["bam_path"]))
VCF_FOR    = dict(zip(samples_df["sample_id"], samples_df["vcf_path"]))

OUTDIR  = config["output_dir"]
GGTYPER = config.get("ggtyper", "./ggtyper")
THREADS = config.get("threads", 4)
CALLER  = config.get("caller", "manta")


rule all:
    input:
        f"{OUTDIR}/genotype_results.tsv"


# ── Step 1: Extract calls from the desired caller per sample ────────────────
rule extract_caller_vcf:
    input:
        vcf=lambda wc: VCF_FOR[wc.sample]
    output:
        vcf=f"{OUTDIR}/per_sample/{{sample}}_caller.vcf.gz",
        tbi=f"{OUTDIR}/per_sample/{{sample}}_caller.vcf.gz.tbi"
    params:
        caller=CALLER
    shell:
        """
        bcftools view -i 'INFO/svdb_origin="{params.caller}"' \
          {input.vcf} -O z -o {output.vcf}
        bcftools index -t {output.vcf}
        """


# ── Step 2: Merge per-sample VCFs into one population variant set ───────────
rule merge_population_vcf:
    input:
        vcfs=expand(f"{OUTDIR}/per_sample/{{sample}}_caller.vcf.gz", sample=SAMPLES)
    output:
        vcf=f"{OUTDIR}/merged/population_sv.vcf"
    params:
        vcf_list=lambda wc, input: " ".join(input.vcfs)
    shell:
        """
        mkdir -p {OUTDIR}/merged
        svdb --merge --pass_only --vcf {params.vcf_list} > {output.vcf}
        """


# ── Step 3: Convert merged VCF to GGTyper JSON ──────────────────────────────
rule convert_to_json:
    input:
        vcf=f"{OUTDIR}/merged/population_sv.vcf",
        script="vcf_to_ggtyper.py"
    output:
        json=f"{OUTDIR}/merged/variants_raw.json"
    shell:
        "python3 {input.script} {input.vcf} {output.json} --all-callers"


# ── Step 4: Remove variants on non-standard contigs ─────────────────────────
rule filter_json:
    input:
        json=f"{OUTDIR}/merged/variants_raw.json",
        script="filter_variants_json.py"
    output:
        json=f"{OUTDIR}/merged/variants.json"
    shell:
        "python3 {input.script} {input.json} {output.json}"


# ── Step 5: Write BAM list for profile-samples ──────────────────────────────
rule create_bam_list:
    output:
        txt=f"{OUTDIR}/bam_list.txt"
    run:
        with open(output.txt, "w") as f:
            for s in SAMPLES:
                f.write(BAM_FOR[s] + "\n")


# ── Step 6: Profile all samples ─────────────────────────────────────────────
rule profile_samples:
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


# ── Step 7: Profile variants (uses library params from all sample profiles) ──
rule profile_variants:
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


# ── Step 8: Genotype all samples against all variants ───────────────────────
rule genotype:
    input:
        variant_profiles=f"{OUTDIR}/profiles/variantProfiles.txt",
        sample_profiles=f"{OUTDIR}/profiles/sampleProfiles.txt"
    output:
        results=f"{OUTDIR}/genotype_results.tsv"
    params:
        prefix=f"{OUTDIR}/genotype_results",
        threads=THREADS
    shell:
        """
        {GGTYPER} genotype \
          {input.variant_profiles} \
          {input.sample_profiles} \
          {params.prefix} \
          -T {params.threads} -e
        """

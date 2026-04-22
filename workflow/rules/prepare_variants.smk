# ── Rules: convert VCF to GGTyper JSON and prepare BAM list ─────────────────


rule convert_to_json:
    """Convert merged population VCF to GGTyper junction JSON format."""
    input:
        vcf=f"{OUTDIR}/merged/population_sv.vcf"
    output:
        json=f"{OUTDIR}/merged/variants_raw.json"
    params:
        script=f"{SCRIPTS}/vcf_to_ggtyper.py"
    conda:
        "../envs/genotyping.yaml"
    shell:
        "python3 {params.script} {input.vcf} {output.json} --all-callers"


rule filter_json:
    """Remove variants on non-standard contigs (chrUn_*, *_random, *_alt)."""
    input:
        json=f"{OUTDIR}/merged/variants_raw.json"
    output:
        json=f"{OUTDIR}/merged/variants.json"
    params:
        script=f"{SCRIPTS}/filter_variants_json.py"
    conda:
        "../envs/genotyping.yaml"
    shell:
        "python3 {params.script} {input.json} {output.json}"


rule create_bam_list:
    """Write one BAM path per line for ggtyper profile-samples."""
    output:
        txt=f"{OUTDIR}/bam_list.txt"
    run:
        with open(output.txt, "w") as f:
            for s in SAMPLES:
                f.write(BAM_FOR[s] + "\n")

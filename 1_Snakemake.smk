#!/usr/bin/env runsnakemake
configfile: "config.yaml"
BUILD = config["BUILD"]
SAMPLE = config["SAMPLE"]
FASTA = config["FASTA"]
CHROMS = config["CHROMS"]
METHOD = config["METHOD"]
NAME = "%s_%s_%s" % (BUILD, SAMPLE, METHOD)
FASTQS = config["FASTQS"]
SNP = config["SNP"]
THREADS = 24
MM2_PRESET = config["MM2_PRESET"]
OUTDIR = config["OUTDIR"]

rule all:
    input:
        OUTDIR + "/%s.primary.fasta" % BUILD,
        OUTDIR + "/%s.tandemRepeats.bed" % BUILD,
        OUTDIR + "/%s.map-hifi.mmi" % BUILD,
        OUTDIR + "/%s.mm2.bam" % NAME,
        OUTDIR + "/%s.haplotag.bam" % NAME,
        OUTDIR + "/%s.haplotag.flagstat" % NAME,
        OUTDIR + "/%s.sniffles2.vcf.gz" % NAME,

rule make_fasta:
    input:
        fa = FASTA
    output:
        fa = OUTDIR + "/%s.primary.fasta" % BUILD
    shell:
        """
        samtools faidx {input.fa} {CHROMS} > {output.fa}
        samtools faidx {output.fa}
        """

rule find_tandem_repeats:
    input:
        fa = rules.make_fasta.output.fa
    output:
        bed = OUTDIR + "/tandemRepeats/%s.{chrom}.bed" % BUILD
    log:
        OUTDIR + "/tandemRepeats/%s.{chrom}.log" % BUILD
    shell:
        """
        findTandemRepeats --merge --chrom {wildcards.chrom} {input.fa} {output.bed} &> {log}
        """

rule merge_randem_repeats:
    input:
        beds = expand(rules.find_tandem_repeats.output.bed, chrom=CHROMS)
    output:
        bed = OUTDIR + "/%s.tandemRepeats.bed" % BUILD
    shell:
        """
        cat {input.beds} | sort -k1,1 -k2,2n > {output.bed}
        bgzip -c {output.bed} > {output.bed}.gz
        tabix -p bed {output.bed}.gz
        """

rule build_index:
    input:
        fa = rules.make_fasta.output.fa
    output:
        mmi = OUTDIR + "/%s.%s.mmi" % (BUILD, MM2_PRESET)
    log:
        OUTDIR + "/%s.%s.log" % (BUILD, MM2_PRESET)
    threads:
        THREADS
    shell:
        """
        minimap2 -t {threads} -x {MM2_PRESET} -d {output.mmi} {input.fa} &> {log} 
        """

rule minimap2:
    input:
        fqs = FASTQS,
        mmi = rules.build_index.output.mmi
    output:
        bam = OUTDIR + "/%s.mm2.bam" % NAME
    log:
        OUTDIR + "/%s.mm2.log" % NAME
    params:
        rg = "@RG\\tID:%s\\tLB:%s\\tSM:%s" % (SAMPLE, SAMPLE, SAMPLE)
    threads:
        THREADS
    shell:
        """(
        minimap2 -ax {MM2_PRESET} --MD -t {threads} -R '{params.rg}' --secondary=no {input.mmi} {input.fqs} \
            | samtools view -@ {threads} -u -F 4 - \
            | samtools sort -@ {threads} -T {output.bam}_TMP -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule haplotag:
    input:
        fa = FASTA,
        vcf = SNP,
        bam = rules.minimap2.output.bam
    output:
        bam = OUTDIR + "/%s.haplotag.bam" % NAME
    log:
        OUTDIR + "/%s.haplotag.log" % NAME
    threads:
        4
    shell:
        """
        whatshap haplotag --output-threads {threads} -o {output.bam} -r {input.fa} {input.vcf} {input.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule sniffles2:
    input:
        bam = rules.haplotag.output.bam,
        fa = FASTA,
        bed = rules.merge_randem_repeats.output.bed
    output:
        tmp1 = temp(OUTDIR + "/%s.sniffles2.vcf" % NAME),
        vcf1 = OUTDIR + "/%s.sniffles2.vcf.gz" % NAME,
        tmp2 = temp(OUTDIR + "/%s.sniffles2.filtered.vcf" % NAME),
        vcf2 = OUTDIR + "/%s.sniffles2.filtered.vcf.gz" % NAME
    log:
        OUTDIR + "/%s.sniffles2.log" % NAME
    threads:
        THREADS
    shell:
        """(
        set +u; source activate sniffles2
        sniffles -t {threads} --phase --output-rnames --minsvlen 20 --sample-id {SAMPLE} \
            --tandem-repeats {input.bed} --reference {input.fa} -i {input.bam} -v {output.tmp1}
        bgzip -c {output.tmp1} > {output.vcf1}
        tabix -p vcf -f {output.vcf1}
        cat {output.tmp1} | grep '#' > {output.tmp2}
        cat {output.tmp1} | grep -v '#' | grep -v 'INV' | grep -v 'DUP' | grep -v 'BND' >> {output.tmp2}
        bgzip -c {output.tmp2} > {output.vcf2}
        tabix -p vcf -f {output.vcf2} ) &> {log}
        """

rule flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    threads:
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """

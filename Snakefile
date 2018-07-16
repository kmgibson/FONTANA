#! /usr/bin/env python
from future import standard_library
standard_library.install_aliases()

import snakemake
from snakemake.utils import report


configfile: "config.yaml"

''' Rules the do not need cluster job '''
localrules: all, complete_sample, structure


wildcard_constraints:
    pdir=".*",
    sampid="\w+"
    
''' Constants '''
SAMPLES = [f.strip('\n') for f in open(config['samples'],"rU")]
#NPROC = int(check_output('nproc'))


def dict_args(d):
    """ Convert a dictionary to command line arguments"""
    ret = ''
    for k,v in d.items():
        if len(k) == 1:
            ret += '-{} {} '.format(k,v)
        else:
            ret += '--{} {} '.format(k,v)
    return ret


print(SAMPLES)

rule all:
    input:
        expand("samples/{s}/completed.txt", s=SAMPLES)


rule complete_sample:
    output:
        touch("samples/{sampid}/completed.txt")
    input:
        "haps.vcf",
        "bam_list.txt",
        "samples/{sampid}/02_align/bam_summary.txt",
        "samples/{sampid}/02_align/bt2.sorted.bam.bai",
        "samples/{sampid}/02_align/bam_status.txt"


rule structure:
    output:
        "samples/{sampid}/00_raw/"
    shell:
        "mkdir -p {output}"

''' add check to see if bam files are there. '''

rule bamtofastq:
    ''' Covert IonTorrent .BAM file back to FASTQ '''
    input:
        "samples/{sampid}/00_raw/original.bam"
    output:
        "samples/{sampid}/00_raw/original.fq"
    shell:
        "bedtools bamtofastq -i {input} -fq {output}"

''' add check to see if fastqs are good '''

rule flexbar:
    input:
        "samples/{sampid}/00_raw/original.fq"
    output:
        "samples/{sampid}/01_trim/clean.fastq.gz"
    params:
        outdir = "samples/{sampid}/01_trim"
    threads: snakemake.utils.available_cpu_count()
    run:
        flexbar_args = dict_args(config['flexbar'])
        adapters = config['adapters']['SE']
        shell("module load flexbar/3.0.3 && "
            "flexbar "
            "--threads {threads} "
            "{flexbar_args} "
            "--adapters {adapters} "
            "-reads {input} "
            "--target {params.outdir}/clean "
        )


rule align_hg19:
    ''' Align clean FASTQ to HG19 '''
    input:
        "samples/{sampid}/01_trim/clean.fastq.gz",
        expand(config['bt2idx'] + ".{i}.bt2", i = range(1, 5)),
        expand(config['bt2idx'] + ".rev.{i}.bt2", i = range(1, 3))
    output:
        "samples/{sampid}/02_align/bt2.unsorted.bam",
        "samples/{sampid}/02_align/bt2.summary.txt"
    threads: snakemake.utils.available_cpu_count()
    run:
        bowtie2_args = dict_args(config['bowtie2'])
        shell(
			"(bowtie2 "
			"-p {threads} "
			"{bowtie2_args} "
			"--rg-id $(basename $(dirname $(dirname {input[0]}))) "
			"--rg SM:$(basename $(dirname $(dirname {input[0]}))) "
			"-x {config[bt2idx]} "
			"-U {input[0]} | "
			"samtools view -b > {output[0]}"
            ") 3>&1 1>&2 2>&3 | tee {output[1]} "
        )


rule sortbam:
    ''' sort and index unsorted.bam file '''
    input:
        "samples/{sampid}/02_align/bt2.unsorted.bam"
    output:
        "samples/{sampid}/02_align/bt2.sorted.bam",
        "samples/{sampid}/02_align/bt2.sorted.bam.bai"
    shell:
        "samtools sort {input} -o {output[0]} && "
        "samtools index {output[0]}"


rule checkbam:
    input:
        "samples/{sampid}/02_align/bt2.sorted.bam"
    output:
        "samples/{sampid}/02_align/bam_summary.txt"
    shell:
        "picard ValidateSamFile "
        "I={input} "
        "MODE=SUMMARY > {output}"


rule checkbamoutput:
    input:
        "samples/{sampid}/02_align/bam_summary.txt"
    output:
        "samples/{sampid}/02_align/bam_status.txt"
    params:
        outdir = "samples/{sampid}/02_align"
    shell:
        "bash checkbamoutput.sh {input} {output} {params.outdir} "

rule bamlist:
    output:
        "bam_list.txt"
    shell: 
        "ls samples/*/02_align/bt2.sorted.bam >> {output} "

rule freebayes_variant_call:
    input:
        "bam_list.txt"
    output:
        "vars.vcf"
    run:
        freebayes_args = dict_args(config['freebayes']['variant'])
        shell(
            "freebayes "
            "--fasta-reference {config[fasta_ref]} "
            "{freebayes_args} "
            "--targets {config[mh_bedfile]} "
            "-L {input} | "
            'vcffilter -f "QUAL > {config[vcffilter_qual]}" > {output} '
        )


rule index_vcf:
    input:
        "vars.vcf"
    output:
        "vars.vcf.gz",
        "vars.vcf.gz.tbi"
    shell:
        "bgzip {input} && "
        "tabix {output[0]}"


rule freebayes_hap:
    input:
        "bam_list.txt",
        "vars.vcf.gz"
    output:
        "haps.vcf"
    run:
        freebayes_args = dict_args(config['freebayes']['haps'])
        shell(
            "freebayes "
            "--fasta-reference {config[fasta_ref]} "
            "--targets {config[mh_bedfile]} "
            "--haplotype-basis-alleles {input[1]} "
            "{freebayes_args} "
            "-L {input[0]} > {output} "
        )









# snakemake --forceall --dag | dot -Tpdf > dag3.pdf


#
# for f in Upload/*/*.bam; do name=$(echo $(basename $f) | cut -d"." -f1); ln -s ../../../$f samples/${name}/00_raw/original.bam; done
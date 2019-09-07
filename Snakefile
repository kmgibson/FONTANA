#! /usr/bin/env python
# -*- coding: utf-8 -*-
from future import standard_library
standard_library.install_aliases()

import snakemake
from snakemake.utils import report


configfile: "config.yaml"

''' Rules the do not need cluster job '''
localrules: all, complete_sample


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

print("Samples include: ", SAMPLES)

rule all:
    input:
        expand("samples/{s}/completed.txt", s=SAMPLES)


rule complete_sample:
    output:
        touch("samples/{sampid}/completed.txt")
    input:
        "haps.vcf",
        "bam_list.txt",
        "samples/{sampid}/bam_summary.txt",
        "samples/{sampid}/bt2.sorted.bam.bai",
        "samples/{sampid}/bam_status.txt"


''' HG19 References '''
include: "Snakefile.references"

''' Processing correct data '''
''' This decides whether there is Illumina data or IonTorrent data'''
snfile = "Snakefile.process%s" % config['platform']
print("Using the snakemake file: ", snfile)
include: snfile

''' Haplotype calling '''
include: "Snakefile.haplotype"



# snakemake --forceall --dag | dot -Tpdf > dag3.pdf
#
# for f in Upload/*/*.bam; do name=$(echo $(basename $f) | cut -d"." -f1); ln -s ../../../$f samples/${name}/00_raw/original.bam; done
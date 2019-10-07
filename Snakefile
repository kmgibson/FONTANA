#! /usr/bin/env python
# -*- coding: utf-8 -*-
from future import standard_library
standard_library.install_aliases()

import snakemake
from snakemake.utils import report
from datetime import date
import os.path
import random
import string

configfile: "config.yaml"

''' Rules the do not need cluster job '''
localrules: all, complete_sample


wildcard_constraints:
    pdir=".*",
    sampid="\w+",
    fontana="fontana_.*"
    
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

def randomStringDigits(stringLength=6):
    """Generate a random string of letters and digits """
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))

''' Create a directory '''
rand_string = randomStringDigits(6)
fontana_dir = 'fontana_%s' % date.today()
alt_fontana_dir = '%s_%s' % (fontana_dir, rand_string)
for samp in SAMPLES:
    dirpath = 'samples/%s/%s' % (samp, fontana_dir)
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    elif os.path.exists(dirpath):
         os.makedirs(alt_fontana_dir)
         fontana_dir = alt_fontana_dir
    

print("Samples include: ", SAMPLES)
print("Directory for this run: ", fontana_dir)


rule all:
    input:
        expand("samples/{s}/{f}/completed.txt", s=SAMPLES, f=fontana_dir)


rule complete_sample:
    output:
        touch("samples/{sampid}/{fontana}/completed.txt")
    input:
        "haps.vcf",
        "bam_list.txt",
        "samples/{sampid}/{fontana}/bam_summary.txt",
        "samples/{sampid}/{fontana}/bt2.sorted.bam.bai",
        "samples/{sampid}/{fontana}/bam_status.txt"


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
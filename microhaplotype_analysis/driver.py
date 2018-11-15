#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *
from collections import defaultdict, OrderedDict
import itertools

from Bio import SeqIO

from mhanalysis import *



# Using FASTA file
ref_file = 'references.fasta'

mh_dict = OrderedDict()
regions = defaultdict(list)

ref_iter = (s for s in SeqIO.parse(ref_file, 'fasta'))
for r in ref_iter:
    f = r.id.split('_')
    chrom = f[0]
    mhid = '_'.join(f[1:3])
    start, end = int(f[-3]), int(f[-2])
    mh = Microhaplotype(mhid, chrom, start, end)
    mh_dict[mhid] = mh
    mh.ref = str(r.seq).upper()
    regions[chrom].append((start, end, mhid))


# Parse VCF
vcf_file = 'haps.vcf'
vcfgen = (l.strip('\n').split('\t') for l in open(vcf_file, 'rU') if not l.startswith('##'))
header = next(vcfgen)

sample_ids = header[9:]
sample_genotypes = {sid:defaultdict(list) for sid in sample_ids}

vars = [VCFLine(x) for x in vcfgen]

for varnum, var in enumerate(vars):
    print('[Variant %d]' % (varnum+1))
    print('\t%s %s:%d' % ('Location:'.ljust(15), var.chrom, var.pos))
    print('\t%s %d' % ('Alt. alleles:'.ljust(15), len(var.alt)))
    print('\t%s ?' % ('Indels:'.ljust(15)))    
    found_in = []
    for reg_s, reg_e, mhid in regions[var.chrom]:
        if reg_s <= var.pos < reg_e:
            found_in.append(mhid)
    assert len(found_in) == 1
    print('\t%s %s' % ('Microhaplotype:'.ljust(15), found_in[0]))
    cur_mh = mh_dict[found_in[0]]
    assert cur_mh.overlaps(var)
    assert var.ref == cur_mh.get_ref(var.pos, len(var.ref))
    cur_mh.mutate(var)
    # Sample genotypes
    for sid, gt in zip(sample_ids, var.genotypes):
        sample_genotypes[sid][cur_mh.mhid].append(gt['GT'])

for mhid, cur_mh in mh_dict.items():
    print(cur_mh)

for sid in sample_ids:
    print('[Sample %s]' % sid)
    for mhid, cur_mh in mh_dict.items():
        mhgt = sample_genotypes[sid][mhid]
        _poss = list(itertools.product(*[_.split('/') for _ in mhgt]))
        _poss = sorted(set([''.join(t) for t in _poss]))
        if len(_poss) == 1:
            mh_allele_1 = mh_allele_2 = _poss[0]
            print('\tMH %s: %s (homozygous)' % (mhid, ';'.join(mhgt)))
            print('\t\tAllele %s: %s' % (mh_allele_1, cur_mh.alleles[mh_allele_1]))
            print('\t\tAllele %s: %s' % (mh_allele_2, cur_mh.alleles[mh_allele_2]))
        elif len(_poss) == 2:
            mh_allele_1, mh_allele_2 = _poss
            print('\tMH %s: %s (heterozygous)' % (mhid, ';'.join(mhgt)))
            print('\t\tAllele %s: %s' % (mh_allele_1, cur_mh.alleles[mh_allele_1]))
            print('\t\tAllele %s: %s' % (mh_allele_2, cur_mh.alleles[mh_allele_2]))
        elif len(_poss) > 2:
            print('\tMH %s: %s (ambigous)' % (mhid, ';'.join(mhgt)))
            for mh_allele in _poss:
                print('\t\tAllele %s: %s' % (mh_allele, cur_mh.alleles[mh_allele]))




# Using BED file
"""
reg_file = 'chr1_mh.bed'

mh_dict = {}
regions = defaultdict(list)

reg_iter = (l.strip('\n').split('\t') for l in open(reg_file, 'rU'))
_ = next(reg_iter) # get rid of dummy line
for l in reg_iter:
    mh = Microhaplotype(l[3], l[0], int(l[1]), int(l[2]))
    mh_dict[l[3]] = mh
    regions[l[0]].append((int(l[1]), int(l[2]), l[3]))
"""

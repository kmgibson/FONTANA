#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import argparse
from collections import defaultdict, OrderedDict, Counter

import pysam


def args_params(args):
    """ Returns a dictionary from argparse namespace
        Excludes "func" argument
    """
    d = {k:v for k,v in vars(args).items() if v is not None}
    if 'func' in d: d.pop('func')
    return d


def trim_read_region(read, spos, epos):
    """ Get read sequence aligned to reference region 
    """
    apairs = read.get_aligned_pairs()
    for i in range(len(apairs)):
        qstart, rp = apairs[i]
        if rp is not None and rp >= spos: break
    
    for i in range(len(apairs)-1,-1,-1):
        qend, rp = apairs[i]
        if rp is not None and rp <= epos: break
    return read.query_sequence[qstart:qend]


def filter_spanning_reads(bam_file=None, region_file=None, bam_out=None):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    out_samfile = pysam.AlignmentFile(bam_out, "wb", template=samfile)
    
    regions = defaultdict(list)
    reg_iter = (l.strip('\n').split('\t') for l in open(region_file, 'rU'))
    # _ = next(reg_iter) # get rid of dummy line
    for i,l in enumerate(reg_iter):
        chrom, spos, epos = l[0], int(l[1]), int(l[2])
        name = l[3] if len(l) > 3 else 'region%02d' % (i+1)
        print('%s - %s:%d-%d' % (name, chrom, spos, epos))
        total_reads = good_reads = 0
        seqs = []
        for aseg in samfile.fetch(chrom, spos, epos):
            total_reads += 1
            if aseg.reference_end is not None:
                if aseg.reference_start <= spos and aseg.reference_end >= epos:
                    good_reads += 1
                    seqs.append(trim_read_region(aseg, spos, epos))
                    out_samfile.write(aseg)
        print('%d overlapping reads, %d spanning reads' % (total_reads, good_reads))
        c = Counter(seqs).most_common()
        print('%d unique sequences' % (len(c)))
        print('Top 5:')
        print('\n'.join('%d\t%s' % (t[1],t[0]) for t in c[:5]))
    
    samfile.close()
    out_samfile.close()


def console():
    parser = argparse.ArgumentParser(
        description='Filter reads that span full region.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('bam_file',
                        help="bam file")
    parser.add_argument('region_file',
                        help="Regions (bed file)")
    parser.add_argument('bam_out', default='filtered.bam',
                        help="Filtered bam file")
    args = parser.parse_args()
    filter_spanning_reads(**args_params(args))


if __name__ == '__main__':
    console()

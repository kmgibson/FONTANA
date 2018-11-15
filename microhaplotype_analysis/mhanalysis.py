# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

from collections import OrderedDict,defaultdict

from Bio import SeqIO

def numstr(x):
    """ Convert argument to numeric type.

    Attempts to convert argument to an integer. If this fails, attempts to
    convert to float. If both fail, return as string.

    Args:
        x: Argument to convert.

    Returns:
        The argument as int, float, or string.

    """
    try:
        return int(x)
    except ValueError:
        try:
            return float(x)
        except ValueError:
            return str(x)

def infostr(x):
    """ Convert argument to dictionary.

    Converts a VCF info field to a dictionary. The attribute string
    should be a semicolon-separated list of tag-value pairs. See
    https://samtools.github.io/hts-specs/VCFv4.2.pdf.

    Args:
        x: Argument to convert.

    Returns:
        A dictionary containing the tag-value pairs.

    """
    if isinstance(x, dict):
        return x
    ret = OrderedDict()
    for kv in x.split(';'):
        k, v = kv.split('=')
        if ',' in v:
            ret[k] = list(map(numstr, v.split(',')))
        else:
            ret[k] = numstr(v)
    return ret



class VCFLine(object):
    # COLS is a list of tuples for parsing VCFLine. 
    # The tuple elements are: (field_name, function, default) where `field_name` is the
    # name of the field, `function` is a function to convert a string into the correct 
    # type, and `default` is a default value for the field.
    COLS = [
        ('chrom', str, '.'),
        ('pos', int, -1),
        ('id', str, '.'),
        ('ref', str, '.'),
        ('alt', lambda x:x.split(','), '.'),
        ('qual', numstr, -1),
        ('filter', lambda x:x.split(';'), '.'),
        ('info', infostr, '.'),
    ]
    
    def __init__(self, val=None):
        if val is None:
            # Create an empty VCFLine object
            for (n, t, d) in self.cols:
                setattr(self, n, d)
            self.genotypes = None
        elif isinstance(val, list):
            # Create VCFLine object from list of strings
            for (n, t, d), v in zip(self.COLS, val):
                setattr(self, n, t(v))
            
            # Parse genotype fields            
            if len(val) <= len(self.COLS):
                self.genotypes = None
            else:
                self.genotypes = []
                gt_format = val[8].split(':')
                for gt_str in val[9:]:
                    d = OrderedDict()
                    for k,v in zip(gt_format, gt_str.split(':')):
                        if ',' in v:
                            d[k] = list(map(numstr, v.split(',')))
                        else:
                            d[k] = numstr(v)
                    self.genotypes.append(d)
        elif isinstance(val, VCFLine):
            # WARNING: this is not a deep copy        
            # Create VCFLine object from another VCFLine      
            for (n, t, d) in self.COLS:
                setattr(self, n, getattr(val, n))
            self.genotypes = val.genotypes
    
    def num_genotypes(self):
        return len(self.genotypes)
    
    def _fmt(self):
        """ Internal format function for first 8 columns """
        ret = []
        for n,t,d in self.COLS:
            attval = getattr(self, n)
            if isinstance(attval, list):
                ret.append(','.join(str(_) for _ in attval))
            elif isinstance(attval, dict):
                dstr = []
                for k,v in attval.items():
                    if isinstance(v, list):
                        dstr.append('%s=%s' % (k, ','.join(str(_) for _ in v)))
                    else:
                        dstr.append('%s=%s' % (k,v))
                ret.append(';'.join(dstr))
            else:
                ret.append(str(attval))
        return ret

    def _fmt_genotypes(self):
        """ Internal format function for genotype columns """
        if self.genotypes is None:
            return []
        
        ret = []
        gt_format = self.genotypes[0].keys()
        ret.append(':'.join(gt_format))
        for gt in self.genotypes:
            gstr = []
            for k in gt_format:
                if isinstance(gt[k], list):
                    gstr.append(','.join(str(_) for _ in gt[k]))
                else:
                    gstr.append(str(gt[k]))
            ret.append(':'.join(gstr))
        return ret
    
    def __str__(self):
        return '\t'.join(self._fmt() + self._fmt_genotypes())
    
class Microhaplotype(object):
    def __init__(self, mhid, chrom, start, end):
        self.mhid = mhid
        self.chrom = chrom
        self.start = start if isinstance(start, int) else int(start)
        self.end = end if isinstance(end, int) else int(end)
        self.ref = None
        self.alleles = None
        self.varcount = 0
    
    def overlaps(self, var):
        if isinstance(var, VCFLine):
            return self.chrom == var.chrom and (self.start <= var.pos < self.end)
        else:
            raise TypeError("Expected VCFLine")
    
    def get_ref(self, pos, length):
        relpos = (pos - self.start - 1)
        return self.ref[relpos:(relpos+length)]
    
    def mutate(self, var):    
        """ Mutate the reference allele to alternate allele(s)"""
        if not self.overlaps(var):
            return
        
        self.varcount += 1
        relpos = (var.pos - self.start -1)
        assert self.ref[relpos:(relpos+len(var.ref))] == var.ref
                                
        if self.alleles is None:
            self.alleles = {}
            
            # Get the left and right sides that are not part of allele
            lside = self.ref[:relpos]
            rside = self.ref[(relpos+len(var.ref)): ]
            
            assert (lside + var.ref + rside) == self.ref
            self.alleles['0'] = lside + var.ref + rside
            
            for i, aa in enumerate(var.alt):
                self.alleles['%d' % (i+1)] = lside + aa + rside
        else:
            prev_alleles = self.alleles
            self.alleles = {}
            for akey in sorted(prev_alleles):
                lside = prev_alleles[akey][:relpos]
                rside = prev_alleles[akey][(relpos+len(var.ref)): ]
                self.alleles['%s0' % akey] = lside + var.ref + rside
                for i, aa in enumerate(var.alt):
                    self.alleles['%s%d' % (akey, i+1)] = lside + aa + rside
    
    def __str__(self):
        ret = []
        ret.append('[Microhaplotype %s] %s:%d-%d; variants=%d' % (self.mhid, self.chrom, self.start, self.end, self.varcount))
        if self.alleles is not None:
            for akey in sorted(self.alleles.keys()):
                ret.append('\t%s: %s' % (akey, self.alleles[akey]))
        return '\n'.join(ret)

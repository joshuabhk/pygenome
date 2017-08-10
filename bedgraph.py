# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 16:06:21 2017

@author: Bong-Hyun Kim

Read bedgraph file and produce various summaries.
"""

from numpy import zeros
#from pybedtools import bedtool
import sys
import matplotlib.pyplot as pp
from gzip import GzipFile

class BedgraphCounter:
    def __init__( self, bdgfn = '', sizefn = 'mm9.chrom.sizes' ) :
        self.sizefn = sizefn
        
        self.lengths = { l.split()[0]:int(l.split()[1]) for l in open( self.sizefn ) if l.split() }
        
        self.counts = {k:zeros(v) for k,v in self.lengths.items()}
        self.bdgfn = bdgfn 
        
        if self.bdgfn :
            if bdgfn.endswith('gz'):
                fp = GzipFile(bdgfn)
            else :
                fp = open(bdgfn)
        
            self.read_bedgraph(fp)
            
    def read_bedgraph( self, fp ) :
        for l in fp :
            if l.startswith(b'track'):
                continue
            
            line = l.split()
            chrom = line[0].decode()
            start = int(line[1])
            end = int(line[2])
            val = float(line[3])
            
            if chrom in self.counts :
                self.counts[chrom][start:end] += val
    
    def normalize( self, factor=1000000000.0 ) :
        total = sum([v.sum() for v in self.counts.values()])
        mf = factor/total
        for k,v in self.counts.items() :
            self.counts[k] = v*mf
    
    def sum( self ) :
        return sum([v.sum() for v in self.counts.values()])
    
    def get_value( self, chrom, start, end, add_chr=False ) :
        if chrom not in self.counts and add_chr and not chrom.startswith('chr') :
            chrom = 'chr'+chrom
        
        return self.counts[chrom][start:end].sum()
    
    def bed2file( self, bed, ofp=sys.stdout ):
        for rec in bed :
            value = self.get_value(rec[0], rec[1], rec[2])
            if len(rec) < 4 :
                print( rec[0], rec[1], rec[2], value, sep='\t', file=ofp )
            else :
                print( rec[0], rec[1], rec[2], value, '\t'.join(rec[3:]), sep='\t', file=ofp)
    
    def bed2average( self, bed, strand_col=None, plot=True ):
        # 'strand_col == None' means no strandedness counted.
        # 'strand_col > 3' means strandedness counted. -> not checked.
        # 'strandedness counted' means 
        #    that the regions are symmetrically defined from the center of the region
        # All the region should have the same length.
        profile = None
        for rec in bed :
            if profile == None :
                length = rec[2]-rec[1]
                profile = zeros(length)
            
            if strand_col != None :
                strand = rec[strand_col]
            else :
                strand = '.'
            
            if strand == '+' or strand == '.':
                temp = self.counts[rec[0]][rec[1]:rec[2]]
            elif strand == '-' :
                temp = self.counts[rec[0]][rec[1]:rec[2]][::-1]
            else :
                assert False
                
            if len(temp) != len(profile):
                print(rec, 'skipped')
                continue
                
            profile += temp
        
        if plot :
            pp.plot( profile )
            pp.show()
        
        return profile
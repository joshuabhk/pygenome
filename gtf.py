#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 16:06:21 2017

@author: Bong-Hyun Kim

My own stripped down GTF file reader.
Focused on easy and flexible parsing on GENCODE GTF file.
"""

class Gene :
    def __init__(self, gene_id="", name="", line="") :
        self.name = name
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None
        self.gene_id = id
        self.info = []
        self.type = None
        
        if line :
            self.add_info( line )
            
            l = line.split('\t')
            if l[2] == 'gene' :
                info = l[-1].split(';')
                temp_id = [i.split()[1].replace('"','') for i in info if i.startswith('gene_id') ]
                if temp_id :
                    self.gene_id = temp_id[0]

                temp_name = [i.split()[1].replace('"','') for i in info if i.startswith(' gene_name') ]
                if temp_name :
                    self.name = temp_name[0]
                
                temp_type = [i.split()[1].replace('"','') for i in info if i.startswith(' transcript_type') ]
                if temp_type :
                    self.type = temp_type[0]
                    
                self.chrom = l[0]
                self.start = int(l[3])
                self.end = int(l[4])
                self.strand = l[6]
                
                if self.strand == '-' :
                    self.start, self.end = self.end, self.start
                
    def add_info( self, infoline ):
        self.info.append(infoline)
        
    def get_promoter_region( self, *args, up=1000, down=1000 ) :
        if self.strand == '-':
            return self.chrom, self.start -down, self.start + 1 + up, '\t'.join(args), self.strand
        else :
            return self.chrom, self.start -up, self.start + 1 + down, '\t'.join(args), self.strand
        
        
class GTF :
    def __init__( self, fn ) :
        self.fn = fn
        self.header = []
        
        self.genes = {}
        self.n2g = {}
        
        if fn :
            self._parse() 
        
    def _parse( self ) :
        fp = open(self.fn)
        for l in fp :
            #header line
            if l.startswith( '#' ) :
                self.header.append(l)
                continue
            
            self._add_info( l )
    
    def _parse_gene_id( self, infoline ):
        info = infoline.split('\t')[-1].split(';')
        temp_id = [i.split()[1].replace('"','') for i in info if i.startswith('gene_id') ]
        return temp_id[0]
    
    def _add_info( self, infoline ):
        gene_id = self._parse_gene_id( infoline )
        
        if gene_id in self.genes :
            gene = self.genes[gene_id]
            gene.add_info(infoline)
        else :
            gene = Gene( line = infoline )
            self.genes[gene_id] = gene
            
    def name2gene(self, name):
        #return one gene 
        if not self.n2g :
            self.n2g = {v.name:v for v in self.genes.values()}
        
        if name in self.n2g :
                return self.n2g[name]    
       
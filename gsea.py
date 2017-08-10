#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 16:09:59 2017

@author: kimb8

Reimplementation of BROAD Gene Set Enrichment Analysis
in python with ES (Enrichment Score) 
and AAES (Absolute Area Enrichment Score)

"""
import numpy as np
from numpy import zeros
from scipy.stats import gumbel_r, norm
import matplotlib.pyplot as pp
from seaborn import distplot #, jointplot


class GSEA :
    def __init__(self, 
                 input_list, 
                 nrepeats=100, 
                 pplot=True, 
                 splot=True, 
                 weights=1,
                 weight_type=1) :
        
        self.input_list = input_list
        self.weights = weights
        self.weight_type = weight_type
        
        sum1, sum2, maxs, mins = self.get_random_scores(input_list, 
                                                        n=nrepeats)
        v = self.get_v( input_list )
        s1 = self.get_area(v)
        s2 = self.get_integral(v)
        mx = self.get_max(v)
        mn = self.get_min(v)
        
        if pplot :
            self.fit_distribution(sum1, 
                              distribution=gumbel_r, 
                              plot=pplot)
            pp.show()
            self.fit_distribution(sum2, 
                              distribution=norm, 
                              plot=pplot)
            pp.show()
            self.fit_distribution(maxs, 
                              distribution=gumbel_r, 
                              plot=pplot)
            pp.show()
            self.fit_distribution(-mins, 
                              distribution=gumbel_r, 
                              plot=pplot)
            pp.show()
        
        if splot : self.plot()# input_list )
            
        self.get_p_value( s1, sum1, distribution=gumbel_r )
        self.get_p_value( s2, sum2, distribution=norm )
        self.get_p_value( mx, maxs, distribution=gumbel_r )
        self.get_p_value( -mn, -mins, distribution=gumbel_r )
    
    def plot(self, weight_not_apply=True):
        '''Prepare plot to understand the results
        '''
        if weight_not_apply :
            if self.weight_type == 1 : #early weights
                a = self.input_list * self.weights
                x = np.cumsum(a - a.sum()/len(a) )
                pp.plot(list(x))
                pp.show()
            elif self.weight_type == 2: #late weights
                a = self.input_list
                x = np.cumsum(a - a.sum()/len(a) )
                x = x*self.weights
                pp.plot(list(x))
                pp.show()
        else :
            a = self.input_list
            x = np.cumsum(a - a.sum()/len(a) )
            pp.plot(list(x))
            pp.show()
    
    
    def get_v( self, alist ):
        #helper function to reduce calculation
        
        #if self.weights != None :
        if self.weight_type == 1 :
            alist = alist*self.weights
            return np.cumsum( alist - alist.sum()/len(alist) )
        elif self.weight_type == 2 :
            return np.cumsum( alist - alist.sum()/len(alist) )*self.weights
        
    def get_area( self, v ) :
        #v = np.cumsum( alist - alist.sum()/len(expr2) )
        return v.abs().sum()

    def get_integral( self, v ):
        #v = np.cumsum( alist - alist.sum()/len(expr2) )
        return v.sum()
    
    def get_max( self, v ):
        #v = np.cumsum( alist - alist.sum()/len(expr2) )
        return v.max()
    
    def get_min( self, v ):
        #v = np.cumsum( alist - alist.sum()/len(expr2) )
        return v.min()

    def get_random_scores( self, alist, n=5000 ):
        N = len(alist)

        sum1 = zeros(n)
        sum2 = zeros(n)
        maxs = zeros(n)
        mins = zeros(n)

        for i in range(n) :
            a = alist.sample(N)
            v = self.get_v(a)

            sum1[i] = self.get_area(v)
            sum2[i] = self.get_integral(v)

            maxs[i] = self.get_max(v)
            mins[i] = self.get_min(v)

        #distplot(sum1, axlabel="absolute sum area"); pp.show()
        #distplot(sum2, axlabel="sum area"); pp.show()
        #distplot(maxs, axlabel="highs"); pp.show()
        #distplot(mins, axlabel='lows'); pp.show()

        return sum1, sum2, maxs, mins
    
    def fit_distribution( self, 
                         values,
                         distribution=gumbel_r, 
                         plot=True) :
        
        x=np.linspace(min(values), max(values), 500)
        fit = distribution.fit( values )
        y=distribution(*fit).pdf(x)
        distplot(values, kde=False, norm_hist=True)
        if plot : pp.plot(x,y)
        
        
    def get_p_value( self, 
                    score, 
                    values, 
                    distribution=gumbel_r ) :
        
        fit = distribution.fit( values )
        pvalue = 1-distribution.cdf( [score], *fit )[0]
        
        #correction of sum2 with normal distribution
        #test both tails
        if distribution == norm :
            pvalue = min( pvalue, 1-pvalue )
            
        #print(score, pvalue)
        return pvalue
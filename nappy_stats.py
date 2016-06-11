# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 09:24:39 2016

@author: Brad Kossmann
"""
import math
import numpy as np
    
    
"""
Get the avarage and standard dev for each DNA parameter from per-frame data.
params is the output of get_param_data
Returns dictionary ave_std[bp] = [param][ave, std]
"""
def get_ave_std(params):
    
    ave_std={}
    for key, mat in params.iteritems():
        mat = zip( *mat)
        ave_std[key] = [ [np.average(i), np.std(i)] for i in mat ]
        
    return ave_std
    
    
"""
Finds differences in averagea and spread from each paramter between two 
systems.
ave_std1 and 2 can be the outputs of get_ave_std (for each param set)
Returns dictionaries ave_diff[bp] = [param][diff(ave1, ave2)], spread_diff""
"""
def comp_ave_std(ave_std1, ave_std2): 
    
    ave_diff={}
    spread_diff={}
    for key, params in ave_std1.iteritems():
        ave_diff[key] = \
        [ params[i][0] - ave_std2[key][i][0] for i in range(len(params))]
        
        spread_diff[key] = \
        [ params[i][1] - ave_std2[key][i][1] for i in range(len(params))]
        
    return ave_diff, spread_diff
    
    
"""
do_bin is needed for get_hist.
It takes in a data point and a nested list of [bin_min, bin_index] and returns
the bin index of the bin the data point belongs in.
While this is a general method, it's really designed only to work inside
get_hist.
"""
def do_bin(point, windows):
    
    mid = int(math.floor(len(windows) / 2))
    
    if  len(windows) == 1:  
        return windows[0][1]
        
    elif point < windows[mid][0]:
        return do_bin(point, windows[:mid])
        
    elif point > windows[mid][0]:
        return do_bin(point, windows[mid:])
        
    elif point == windows[mid][0]:
        return windows[mid][1]
        
        
"""
Prepare and return histograms for each parameter based on per-frame data.
params is a params list, as from get_param_data, num_bins is dictionary of 
bin_numbers as param_idx:bin_num
If subtract is false, simple binning returns:
1) per-bp dictionary of nested [param][hist] list
2) per-parameter dictionary of lists of window centers
If subtract is true, requires second set of params *with same bp dict indices*
returns:
1) histograms for first set of parameters
2) histograms for second set of parameters
3) set of second - first histograms
4) per-parameter ""
"""
def get_hist(params, num_bins, params2=None, subtract=False):
    
    #transpose each per-bp array for convenience
    for bp, ar in params.iteritems():
        params[bp] = zip( *ar)
        
        if subtract == True:
            params2[bp] = zip( *params2[bp])
    
    #Need to find global (across all bps) min-max for each param      
    bins={}
    for bp, ar in params.iteritems():
        
        for i in range(len(ar)):
            
            if i not in bins:
                bins[i] = [min(ar[i]), max(ar[i])]
                
            else:
                if min(ar[i]) < bins[i][0]:
                    bins[i][0] = min(ar[i])
                
                if max(ar[i]) > bins[i][1]:
                    bins[i][1] = max(ar[i])
                    
    if subtract == True: #Not sure if this can be handled more efficiently
    
        for bp, ar in params2.iteritems():
            
            for i in range(len(ar)):
                
                if i not in bins:
                    return 'Incompatible parameter sets' 
                    
                else:
                    if min(ar[i]) < bins[i][0]:
                        bins[i][0] = min(ar[i])
                    if max(ar[i]) > bins[i][1]:
                        bins[i][1] = max(ar[i])
    
    #Figure out bin widths
    widths=[0]*len(bins)
    for key, pair in bins.iteritems():
        width = (pair[1] - pair[0] ) / num_bins[key]   
        widths[key] = width
        bins[key] = [ [ (pair[0] + i 
                         * width), i] for i in range(num_bins[key]) ]
    #Bins is now a per-param dictionary of [bin_minimum, bin_index] lists
        
    histograms={} #Will store the final histograms
    for bp, ar in params.iteritems(): #For each bp
        histograms[bp]=[] #Each bp gets a new array to store the histograms
        
        for i in range(len(ar)): #For each parameter
            #Not going to use the counter here for convenience
            count={ j:0 for j in range(num_bins[i]) } #Initialize counter
            histograms[bp].append( [0]*(num_bins[i]+1) ) #Initialize histogram list
            hist = [ do_bin(j, bins[i]) for j in ar[i] ] #Make histogram
            
            for j in hist: #Count up histogram populations
                count[j]+=1
                
            for bin_idx, pop in count.iteritems():
                histograms[bp][i][bin_idx] = pop #Move populations to histogram
                
    if subtract == True:
        
        histograms2={}
        for bp, ar in params2.iteritems():
            histograms2[bp]=[]
            
            for i in range(len(ar)):
                count={ j:0 for j in range(num_bins[i]) }
                histograms2[bp].append( [0]*(num_bins[i]+1) )
                hist = [ do_bin(j, bins[i]) for j in ar[i] ]
                
                for j in hist:
                    count[j]+=1
                    
                for bin_idx, pop in count.iteritems():
                    histograms2[bp][i][bin_idx] = pop
                
    #Convert bins to dictionary of window centers lists
    for key, ar in bins.iteritems(): 
        bins[key] = [ j[0] + widths[key] for j in ar  ]
        
    if subtract == True:
        hist_sub={}
        
        for key, ar in histograms.iteritems():
            hist_sub[key] = \
            [ np.subtract(histograms2[key][i], ar[i]).tolist() 
            for i in range(len(ar)) ]
        
    if subtract == False:
        return histograms, bins
        
    else:
        return histograms, histograms2, hist_sub, bins
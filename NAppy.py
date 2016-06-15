# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 12:51:52 2016

@author: Brad Kossmann
"""
import time
import nappy_stats as napstat
import nappy_io as napio
import nappy_plot as napplot

idx, PU1_params = napio.get_param_data('../PU1/PU1.nastruct.dat')
DNA_params = napio.get_param_data('../PU1/DNA.nastruct.dat', indices=False)

DNA_ave_std = napstat.get_ave_std(DNA_params)
PU1_ave_std = napstat.get_ave_std(PU1_params)

ave_diff, spread_diff = napstat.comp_ave_std(PU1_ave_std, DNA_ave_std)

bin_nums = {}
for key, ind in idx.iteritems():
    bin_nums[ind] = 500 #Going w/100 bins per parameter for testing purposes
    
start= time.time()
'''
To use get_hist separately for each param set
DNA_hists, DNA_bins = napstat.get_hist(DNA_params, bin_nums)
PU1_hists, PU1_bins = napstat.get_hist(PU1_params, bin_nums)
'''

PU1_hists, DNA_hists, sub_hist, bins = napstat.get_hist(PU1_params, bin_nums, 
                                                        DNA_params, 
                                                        subtract=True)

biv_ker = napstat.biv_ker_dens(DNA_params, 'Shear', 'Stagger', idx, 100j)
napplot.plot_dens(biv_ker, marker=0.5)
napplot.hist_plot(DNA_hists, 'Minor', idx, normalize='True')
#print sub_hist
print time.time() - start 
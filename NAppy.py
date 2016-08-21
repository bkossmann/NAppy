# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 12:51:52 2016

@author: Brad Kossmann
"""
import time
import nappy_stats as napstat
import nappy_io as napio
import nappy_plot as napplot

#idx, PU1_params = napio.get_param_data_amber('../PU1/PU1.nastruct.dat')
#DNA1_params = napio.get_param_data_amber('../PU1/DNA.nastruct.dat', indices=False)

idx, PU1_params = napio.get_param_data_curves('../PU1/PU_curves/1/1_', 15,)

DNA1_params = napio.get_param_data_curves('../PU1/DNA_curves/1/1_', 15,
                                               indices=False)

#print idx
#idx, PU1_params = napio.get_param_data('../PU1/PU_1.dat')
#DNA1_params = napio.get_param_data('../PU1/DNA_1.dat', indices=False)
'''
PU2_params = napio.get_param_data('../PU1/PU_2.dat', indices=False)
DNA2_params = napio.get_param_data('../PU1/DNA_2.dat', indices=False)

PU3_params = napio.get_param_data('../PU1/PU_3.dat', indices=False)
DNA3_params = napio.get_param_data('../PU1/DNA_3.dat', indices=False)

idx, PU4_params = napio.get_param_data('../PU1/PU_4.dat')
DNA4_params = napio.get_param_data('../PU1/DNA_4.dat', indices=False)

DNA1_ave_std = napstat.get_ave_std(DNA1_params)
PU1_ave_std = napstat.get_ave_std(PU1_params)
ave_diff1, spread_diff1 = napstat.comp_ave_std(PU1_ave_std, DNA1_ave_std)

DNA2_ave_std = napstat.get_ave_std(DNA2_params)
PU2_ave_std = napstat.get_ave_std(PU2_params)
ave_diff2, spread_diff2 = napstat.comp_ave_std(PU2_ave_std, DNA2_ave_std)

DNA3_ave_std = napstat.get_ave_std(DNA3_params)
PU3_ave_std = napstat.get_ave_std(PU3_params)
ave_diff3, spread_diff3 = napstat.comp_ave_std(PU3_ave_std, DNA3_ave_std)

DNA4_ave_std = napstat.get_ave_std(DNA4_params)
PU4_ave_std = napstat.get_ave_std(PU4_params)
ave_diff4, spread_diff4 = napstat.comp_ave_std(PU4_ave_std, DNA4_ave_std)

napio.output_params_csv(ave_diff1, idx, ['Minor'], 'avediff1.csv','bp,Minor')
napio.output_params_csv(ave_diff2, idx, ['Minor'], 'avediff2.csv','bp,Minor')
napio.output_params_csv(ave_diff3, idx, ['Minor'], 'avediff3.csv','bp,Minor')
napio.output_params_csv(ave_diff4, idx, ['Minor'], 'avediff4.csv','bp,Minor')

napio.output_params_csv(DNA1_ave_std, idx, ['Minor'], 'DNA1_minor.csv',
                        'bp,Minor', list_idx=[0])
napio.output_params_csv(DNA2_ave_std, idx, ['Minor'], 'DNA2_minor.csv',
                        'bp,Minor', list_idx=[0])
napio.output_params_csv(DNA3_ave_std, idx, ['Minor'], 'DNA3_minor.csv',
                        'bp,Minor', list_idx=[0])
napio.output_params_csv(DNA4_ave_std, idx, ['Minor'], 'DNA4_minor.csv',
                        'bp,Minor', list_idx=[0])

napio.output_params_csv(PU1_ave_std, idx, ['Minor'], 'PU1_minor.csv',
                        'bp,Minor', list_idx=[0])
napio.output_params_csv(PU2_ave_std, idx, ['Minor'], 'PU2_minor.csv',
                        'bp,Minor', list_idx=[0])
napio.output_params_csv(PU3_ave_std, idx, ['Minor'], 'PU3_minor.csv',
                        'bp,Minor', list_idx=[0])
napio.output_params_csv(PU4_ave_std, idx, ['Minor'], 'PU4_minor.csv',
                        'bp,Minor', list_idx=[0])
                        '''
'''
out1 = open('avediff1.csv', 'w')
out1.write('bp,diff\n')
for key, index in ave_diff1.iteritems():
    out1.write(key.split()[0]+','+str(ave_diff1[key][idx['Minor']])+'\n')
out1.close()

out2 = open('avediff2.csv', 'w')
out2.write('bp,diff\n')
for key, index in ave_diff2.iteritems():
    out2.write(key.split()[0]+','+str(ave_diff2[key][idx['Minor']])+'\n')
out2.close()

out3 = open('avediff3.csv', 'w')
out3.write('bp,diff\n')
for key, index in ave_diff3.iteritems():
    out3.write(key.split()[0]+','+str(ave_diff3[key][idx['Minor']])+'\n')
out3.close()

out4 = open('avediff4.csv', 'w')
out4.write('bp,diff\n')
for key, index in ave_diff4.iteritems():
    out4.write(key.split()[0]+','+str(ave_diff4[key][idx['Minor']])+'\n')
out4.close()
'''

bin_nums = {}
for key, ind in idx.iteritems():
    bin_nums[ind] = 500 #Going w/100 bins per parameter for testing purposes
    
start= time.time()

PU1_hists, DNA1_hists, sub1_hist, bins = napstat.get_hist(PU1_params, bin_nums,
                                                        DNA1_params,
                                                        subtract=True)
                            
'''
PU2_hists, DNA2_hists, sub2_hist, bins = napstat.get_hist(PU2_params, bin_nums, 
                                                        DNA2_params, 
                                                        subtract=True)
                                                        
PU3_hists, DNA3_hists, sub3_hist, bins = napstat.get_hist(PU3_params, bin_nums, 
                                                        DNA3_params, 
                                                        subtract=True)
                                                        
PU4_hists, DNA4_hists, sub4_hist, bins = napstat.get_hist(PU4_params, bin_nums, 
                                                        DNA4_params, 
                                                        subtract=True)
                                                        
#biv_ker = napstat.biv_ker_dens(DNA_params, 'Shear', 'Stagger', idx, 100j)
#napplot.plot_dens(biv_ker, marker=0.5)

'''                                                        '''
napplot.hist_plot(DNA1_hists, 'Minor', idx, 'DNA1.png', normalize='True')
napplot.hist_plot(PU1_hists, 'Minor', idx, 'PU1.png', normalize='True')
napplot.hist_plot(DNA1_hists, 'Minor', idx, 'DNA1.png', normalize='True')
napplot.hist_plot(PU2_hists, 'Minor', idx, 'PU2.png', normalize='True')
napplot.hist_plot(DNA2_hists, 'Minor', idx, 'DNA2.png', normalize='True')
napplot.hist_plot(PU3_hists, 'Minor', idx, 'PU3.png', normalize='True')
napplot.hist_plot(DNA3_hists, 'Minor', idx, 'DNA3.png', normalize='True')
napplot.hist_plot(PU4_hists, 'Minor', idx, 'PU4.png', normalize='True')
napplot.hist_plot(DNA4_hists, 'Minor', idx, 'DNA4.png', normalize='True')
napplot.hist_plot(sub1_hist, 'Minor', idx, 'sub1.png', normalize='True')
napplot.hist_plot(sub2_hist, 'Minor', idx, 'sub2.png', normalize='True')
napplot.hist_plot(sub3_hist, 'Minor', idx, 'sub3.png', normalize='True')
napplot.hist_plot(sub4_hist, 'Minor', idx, 'sub4.png', normalize='True')
#print sub_hist
print time.time() - start 
'''
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 11:59:28 2016

@author: Brad Kossmann
"""
import numpy as np
import matplotlib.pyplot as plt


"""
dens_plot takes in the output from napstat.biv_ker_dens (directly) and plots
the kernel density profile
If a per-bp dictionary is passed, each bp will be plotted separately. 
By default, only one array is passed in.
"""
def plot_dens(ker_dens, marker=0.1):
    
    if not isinstance(ker_dens, dict):
        z = ker_dens[0]
        x = ker_dens[1]
        y = ker_dens[2]
        xmin = ker_dens[3]
        xmax = ker_dens[4]
        ymin = ker_dens[5]
        ymax = ker_dens[6]
        fig, ax = plt.subplots()
        ax.imshow(np.rot90(z), cmap=plt.cm.rainbow,
                  extent=[xmin, xmax, ymin, ymax])
        ax.plot(x, y, 'k.', markersize=marker)
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        plt.show()
        
    else:
        for key, args in ker_dens.iteritems():
            z = args[0]
            x = args[1]
            y = args[2]
            xmin = args[3]
            xmax = args[4]
            ymin = args[5]
            ymax = args[6]
            fig, ax = plt.subplots()
            ax.imshow(np.rot90(z), cmap=plt.cm.rainbow,
                      extent=[xmin, xmax, ymin, ymax])
            ax.plot(x, y, 'k.', markersize=marker)
            ax.set_xlim([xmin, xmax])
            ax.set_ylim([ymin, ymax])
            plt.show()
            
            
"""
hist_plot takes in a histogram dictionary (e.g. from get_hist) and plots the 
parm from parm_idx for base pair bp
normalize automatically normalizes the histogram between 0 and 1 if set to 
'True'
cap removes any bin that is consistently 0 across all rows (bp usually)
Future implementations will include support for x-axis labeling of bins instead
of bin numbers
"""
def hist_plot(hist_dict, parm, parm_idx, name, normalize='False', cap='True'):
    
    hist_sort = {}
    for key, vals, in hist_dict.iteritems(): #Order bases by sequence
        hist_sort[int(key.split()[0])] = vals
    
    hist=[]
    for key, vals in hist_sort.iteritems():
        hist.append(vals[parm_idx[parm]])

    if normalize == 'True':
        import nappy_stats as napstat
        histogram = napstat.norm_hist(hist)
        
    if cap == 'True':      
        
        low_bound = 'None'
        high_bound = 'None'
        histogram = zip(*histogram)
        for i in range(len(histogram)): #Lowest bin index with all 0's
            if low_bound != 'None':
                pass
            elif float(0) not in histogram[i]:
                low_bound = i -1
        
        histogram.reverse()
        for i in range(len(histogram)): #Highest bin index with all 0's
            if high_bound != 'None':
                pass
            elif float(0) not in histogram[i]:
                high_bound = len(histogram) - i
        
        histogram.reverse()
        histogram = zip(*histogram)
        for i in range(len(histogram)):
            histogram[i] = histogram[i][low_bound:high_bound]
            
    #plt.imshow(histogram)
    plt.imsave(name, np.array(histogram))
    plt.close()
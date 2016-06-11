# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 09:23:32 2016

@author: Brad Kossmann
"""


'''
Import the data from a raw ptraj output.
Returns idx as dictionary of param_name:index, params as dictionary of
[bp][frame][parameter]
'''
def get_param_data(data, indices=True):
    
    idx={}
    params={}
    with open(data, 'r') as filein:
    
        for i, x in enumerate(next(filein).strip('\n').strip('#').split()):
            idx[x]=i
            
        for line in filein:
            lines=line.split()
            
            if len(lines) > 0:
                
                if lines[1]+'-'+lines[2] not in params:
                    params[lines[1]+'-'+lines[2]] = \
                    [ [ float(i) for i in lines[3:] ] ]
                    
                else:
                    params[lines[1]+'-'+lines[2]].append \
                    ([ float(i) for i in lines[3:] ])
                
    filein.close()
    
    if indices:
        return idx, params
        
    else:
        return params
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
def get_param_data_amber(data, indices=True):
    
    print "Importing cpptraj data from", data, "\n"
    idx={}
    params={}
    with open(data, 'r') as filein:
    
        for i, x in enumerate(next(filein).strip('\n').strip('#').split()):
            idx[x]=i - 3 #For convenience. "Frame, base1, base2 aren't useful.
            
        for line in filein:
            lines=line.split()
            
            if len(lines) > 0:
                
                if lines[1]+' '+lines[2] not in params:
                    params[lines[1]+' '+lines[2]] = \
                    [[float(i) for i in lines[3:]]]
                    
                else:
                    params[lines[1]+' '+lines[2]].append \
                    ([float(i) for i in lines[3:]])
                
    filein.close()
    
    if indices:
        return idx, params
        
    else:
        return params

'''
Import *.ser data from curves+/canal output.
Prefix is the file path prefix + filename prefix.
All standard output parameters are included, unless a truncated parameters list
is passed into parameters. bp is is the number of base pairs being passed in.
Returns same data as get_param_data_amber.
Warning: any NAs in the data get converted to 0 for simplicity. The user must
know what curves parameters are not calculated for all base pairs and analyze
accordingly. Tbend currently not supported.
'''
def get_param_data_curves(prefix, bp, parameters=[], indices=True):
    
    print "Importing curves+/canal data for set", prefix, "\n"
    idx={}
    
    bp1 = range(1,bp+1)
    bp2 = range(2*bp,bp,-1)
        
    temp = {}
    for x in bp1:
        temp[x] = []
        
    parameters=['alphaC','alphaW','ampC','ax-bend','betaC','betaW','buckle',
                'chiC','chiW','deltaC','deltaW','epsilC','epsilW','gammaC',
                'gammaW','h-ris','h-twi','inclin','majd','majw','mind','minw',
                'opening','phaseC','phaseW','propel','rise','roll','shear',
                'shift','slide','stagger','stretch','tilt','tip',
                'twist','xdisp','ydisp','zetaC','zetaW'] 
                
    for p in range(len(parameters)):
        idx[parameters[p]] = p
        
        for key, index in temp.iteritems():
            temp[key].append([])
            
        filein = open(prefix+parameters[p]+'.ser', 'r')
        for line in filein:
            lines=line.split()
            for i in range(1,len(lines)):
                val = lines[i]
                if val != 'NA':
                    temp[i][len(temp[i])-1].append(float(val))
                else:
                    temp[i][len(temp[i])-1].append(float(0))
        filein.close()
        
    params={}
    
    for key, index in temp.iteritems():
        params[str(bp1[key - 1])+' '+str(bp2[key - 1])] = zip(*index)
    
    if indices:
        return idx, params
        
    else:
        return params
        
'''
params is parameters list, idx is index list (both perhaps from data import),
filename is desired output file name, header is a comma-separated list of 
column headings, list_idx is a list of indices to write in case writing from,
e.g. ave_std where only one component maybe desired for each field.
Warning: only prints out first base in basepair for bp for now
'''
def output_params_csv(params, idx, data_out, filename, header, list_idx=False):
    outfile = open(filename, 'w')
    outfile.write(header+'\n')
    
    for key, index in params.iteritems():
        outfile.write(key.split()[0])
        
        for dat in data_out:
            
            if not list_idx:
                outfile.write(','+str(params[key][idx[dat]]))
            else:
                for idc in list_idx:
                    outfile.write(','+str(params[key][idx[dat]][idc]))
                
        outfile.write('\n')
    outfile.close()
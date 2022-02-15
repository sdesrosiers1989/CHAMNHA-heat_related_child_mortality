#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Dcomposition of infant mortality model components

Das Gupta 1993 method

Created on Mon Apr 20 11:12:29 2020

@author: earsch
"""

#%%set wd and import packages
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec

import iris
import iris.plot as iplt
import iris.quickplot as qplt
import iris.coord_categorisation

import numpy as np
import numpy.ma as ma

import pandas as pd

import math

import cartopy.crs as ccrs

import copy

import glob

from iris.experimental.equalise_cubes import equalise_attributes

proj = ccrs.PlateCarree(central_longitude = 38)

#Import my functions
import sys
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/Tanga/Plot_functions')
import tanzania1 as tp

def add_code_from_filename(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[97:-3]
    
    code = file.partition('_')[2].partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(code, long_name = 'code'))

def ens_mean(cube_list):  
    cor_cubes = iris.cube.CubeList()
    #add ensemble as dimension so can use cube statistics
    for item in cube_list:
        gcm = item.coord('gcm').points[0]
        try:
            item.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'ens'))
        except:
            print('Has ens coord.')
        try:
            item.remove_coord('month')
        except:
            print('No month coord')
        try:
            item.remove_coord('time')
        except:
            print('No time coord')
        try:
            item.remove_coord('year')
        except:
            print('No month coord')
        cor_cubes.append(item)
    ens_cube = cor_cubes.merge_cube()
    ens_cube = ens_cube.collapsed('ens', iris.analysis.MEAN)
    return ens_cube

def find_dif(his_list, histnat_list):
    #his - histnat
    dif_list = []
    difper_list = []
    
    for cube in his_list:
        mod = tp.gcm(cube)
        comp_cube =  [x for x in histnat_list if tp.gcm(x) == mod][0] # checked, len always 1, as it should be
        print(mod, tp.gcm(comp_cube))
        
        dif_cube = copy.deepcopy(cube)
        dif_cube.data = dif_cube.data - comp_cube.data
        dif_list.append(dif_cube)
        
        difper_cube = copy.deepcopy(cube)
        difper_cube.data = (difper_cube.data - comp_cube.data) / difper_cube.data
        difper_list.append(difper_cube)
    
    return dif_list, difper_list

#%%
    
def calc_comp(cube_list, return_list =  False):
    '''order should be:
        'e2p1m1', 'e2p2m1', 'e2p1m2',
        'e1p1m1', 'e1p2m1', 'e1p1m2'
        'e2p2m2', 'e1p2m2'
    '''
    
    p1m1e1_list = cube_list[3]
    p1m2e2_list = cube_list[2]
    p1m1e2_list = cube_list[0]
    p1m2e1_list = cube_list[5]
    p2m1e1_list = cube_list[4]
    p2m2e2_list = cube_list[6]
    p2m1e2_list = cube_list[1]
    p2m2e1_list = cube_list[7]
    
    p_list = []
    m_list = []
    e_list = []
    
    for i in np.arange(len(p1m1e1_list)):
        
        p1m1e1 = p1m1e1_list[i]
        
        #find same gcm
        gcm = tp.gcm(p1m1e1)
        
        p1m2e2 = [x for x in p1m2e2_list if tp.gcm(x) == gcm][0]
        p1m1e2 = [x for x in p1m1e2_list if tp.gcm(x) == gcm][0]
        p1m2e1 = [x for x in p1m2e1_list if tp.gcm(x) == gcm][0]
        p2m1e1 = [x for x in p2m1e1_list if tp.gcm(x) == gcm][0]
        p2m2e2 = [x for x in p2m2e2_list if tp.gcm(x) == gcm][0]
        p2m1e2 = [x for x in p2m1e2_list if tp.gcm(x) == gcm][0]
        p2m2e1 = [x for x in p2m2e1_list if tp.gcm(x) == gcm][0]
    
        peq1 = (p1m1e1 + p1m2e2) / 3
        peq2 = (p1m1e2 + p1m2e1) / 6
        peq3 = (p2m1e1 + p2m2e2) / 3
        peq4 = (p2m1e2 + p2m2e1) / 6
        
        p_effect = (peq1 + peq2) - (peq3 + peq4)
        
        meq1 = (p1m1e1 + p2m1e2) / 3
        meq2 = (p1m1e2 + p2m1e1) / 6
        meq3 = (p1m2e1 + p2m2e2) / 3
        meq4 = (p1m2e2 + p2m2e1) / 6
        
        m_effect = (meq1 + meq2) - (meq3 + meq4)
        
        eeq1 = (p1m1e1 + p2m2e1) / 3
        eeq2 = (p1m2e1 + p2m1e1) / 6
        eeq3 = (p1m1e2 + p2m2e2) / 3
        eeq4 = (p1m2e2 + p2m1e2) / 6
        
        e_effect = (eeq1 + eeq2) - (eeq3 + eeq4)
        
        p_list.append(p_effect)
        m_list.append(m_effect)
        e_list.append(e_effect)
        
    p_ens = ens_mean(p_list)
    m_ens = ens_mean(m_list)
    e_ens = ens_mean(e_list)
        
    if return_list == True:
        return p_ens, m_ens, e_ens, p_list, m_list, e_list
    else:
        return p_ens, m_ens, e_ens
    #return p_ens, m_ens, e_ens, p_list, m_list, e_list



#%% Import outputs of ifnant mortality health burden model
    # compare hist-nat T1 and T2
    
codes = ['e2p1m1', 'e2p2m1', 'e2p1m2',
        'e1p1m1', 'e1p2m1', 'e1p1m2',
        'e2p2m2', 'e1p2m2']

path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/decomp_attribution2/2010/'
damip_mods = []
filenames = glob.glob(path + '*.nc')
cube_list = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file, callback = add_code_from_filename)
    damip_mods.append(tp.gcm(x))
    cube_list.append(x)


#'e2p2m2'
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/'
filenames = glob.glob(path + '*hist-nat*20052014*.nc')
for file in filenames:
    x = iris.load_cube(file)
    x.add_aux_coord(iris.coords.AuxCoord('e2p2m2', long_name = 'code'))
    if tp.gcm(x) in damip_mods:
        cube_list.append(x)


code_cubes = []
for code in codes:
    new_list = [x for x in cube_list if x.coord('code').points == code]
    code_cubes.append(new_list)
    
        
#%% What is the cause of the differences between hist-nat T1 and T2?
    
#1995 - 2004: pop and mortality the same (onyl due to climate)
period1 = [x for x in code_cubes if x[0].coord('code').points == 'e1p1m1'][0]
period2 = [x for x in code_cubes if x[0].coord('code').points == 'e2p2m2'][0]

dif, dif_per = find_dif(period1, period2)

ens_dif = ens_mean(dif)
ens_perdif = ens_mean(dif_per)

period1_ens = ens_mean(period1)
period2_ens = ens_mean(period2)


#%% Cause of difference between Histnat T1 and T2?
    
#codes = ['e2p1m1', 'e2p2m1', 'e2p1m2',
#        'e1p1m1', 'e1p2m1', 'e1p1m2',
#        'e2p2m2', 'e1p2m2']

p_effect, m_effect, e_effect, plist, mlist, elist = calc_comp(code_cubes, return_list = True)
    
#%%
print(np.nanmean(period1_ens.data))
print(np.nanmean(period2_ens.data))
print(np.nanmean(ens_dif.data))


np.nanmean(p_effect.data)
np.nanmean(m_effect.data)
np.nanmean(e_effect.data)

np.nanmean(p_effect.data) + np.nanmean(m_effect.data) +  np.nanmean(e_effect.data) 


#%%

i, j, k = 40,40, 0 

print(period1[k][i,j].data)
print(period2[k][i,j].data)
print(dif[k][i,j].data)


print(plist[k][i,j].data)
print(mlist[k][i,j].data)
print(elist[k][i,j].data)

plist[k][i,j].data + mlist[k][i,j].data +  elist[k][i,j].data

#%% What percentage of the change is due to climate, pop and mortality?

total_change = np.nanmean(ens_dif.data)

p_per = np.nanmean(p_effect.data) / total_change * 100
m_per = np.nanmean(m_effect.data) / total_change * 100
e_per = np.nanmean(e_effect.data) / total_change * 100



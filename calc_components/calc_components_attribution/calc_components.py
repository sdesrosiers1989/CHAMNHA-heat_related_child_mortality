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

def ens_mean(cube_list):  
    cor_cubes = iris.cube.CubeList()
    #add ensemble as dimension so can use cube statistics
    for i in np.arange(len(cube_list)):
        item = cube_list[i]
        gcm = i
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
    
def calc_comp(p_list, m_list, e_list, return_list =  False):
    '''order should be:
    period1, period2
    list of models within each period only for e (pop and mor same each model)
    '''
    
    p1 = p_list[0].data
    p2 = p_list[1].data
    m1 = m_list[0].data
    m2 = m_list[1].data 
    e1_mlist = e_list[0]
    e2_mlist = e_list[1]
    
    p_outlist = iris.cube.CubeList()
    m_outlist = iris.cube.CubeList()
    e_outlist = iris.cube.CubeList()
    
    for i in np.arange(len(e1_mlist)):
        
        e1 = e1_mlist[i].data
        
        #find same gcm
        gcm = tp.gcm(e1_mlist[i])
        e2 = [x for x in e2_mlist if tp.gcm(x) == gcm][0].data
        
        
        peq1 = (m1 * e1 + m2 * e2)  / 3
        peq2 = (m1 * e2 + m2 * e1) / 6
        
        p_effect = copy.deepcopy(p_list[0])
        p_effect.data = (p1 - p2)*(peq1 + peq2)
        
        meq1 = (p1 * e1 + p2 * e2)  / 3
        meq2 = (p1 * e2 + p2 * e1) / 6
                
        m_effect = copy.deepcopy(m_list[0])
        m_effect.data = (m1 - m2)*(meq1 + meq2)
        
        eeq1 = (m1 * p1 + m2 * p2)  / 3
        eeq2 = (p1 * m2 + p2 * m1) / 6
                
        e_effect = copy.deepcopy(e1_mlist[i])
        e_effect.data = (e1 - e2)*(eeq1 + eeq2)
        
        p_outlist.append(p_effect)
        m_outlist.append(m_effect)
        e_outlist.append(e_effect)
        
    p_ens = ens_mean(p_outlist)
    m_ens = ens_mean(m_outlist)
    e_ens = ens_mean(e_outlist)
        
    if return_list == True:
        return p_ens, m_ens, e_ens, p_outlist, m_outlist, e_outlist
    else:
        return p_ens, m_ens, e_ens
    #return p_ens, m_ens, e_ens, p_list, m_list, e_list

#%% import e model data
    
path = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_061/historical/'

#hist-nat
damip_mods = []
filenames = glob.glob(path + '*hist-nat*.nc')
damip_list = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    damip_mods.append(tp.gcm(x))
    damip_list.append(x)


#historical
filenames = glob.glob(path + '*historical*.nc')
hist_list = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    if tp.gcm(x) in damip_mods:
        hist_list.append(x)
  
#%% Import pop and mor data 

#population -from world pop
# 0 - 4 year olds
pop_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000_regrid.nc')
pop_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2010_regrid.nc')
pop_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2019_regrid.nc')

#daily mortality
dmor_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2000_regrid.nc')
dmor_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2010_regrid.nc')
dmor_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2019_regrid.nc')


pop_list = [pop_2000, pop_2010, pop_2019]
mor_list = [dmor_2000, dmor_2010, dmor_2019]

#%% actual mortality

path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/'

filenames = glob.glob(path + '*hist-nat*')
damip_mor = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    damip_mor.append(x)

#%% calc pop ratio as sued in health burden model
#pop data

pop_ratio_list = []
for i in np.arange(len(pop_list)):

    pop_ratio = pop_list[i]/pop_2000 #ratio future pop to baseline pop   
    #0 in denominator causing problems
    pop_ratio.data = np.where(pop_2000.data == 0, 1, pop_ratio.data)
    pop_ratio.data = ma.masked_array(pop_ratio.data, mask = pop_2000.data.mask)  
    pop_ratio_list.append(pop_ratio)
    
#%% get mdoels into year lists
    
years = np.unique([x.coord('year').points[0] for x in hist_list]) #damip and his years same

damip_years =  iris.cube.CubeList()
his_years = iris.cube.CubeList()

for y in years:
    d_list = [x for x in damip_list if x.coord('year').points[0] == y]
    d_mean_list = iris.cube.CubeList()
    for cube in d_list:
        d_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))

    h_list = [x for x in hist_list if x.coord('year').points[0] == y]
    h_mean_list = iris.cube.CubeList()
    for cube in h_list:
        h_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))
    
    damip_years.append(d_mean_list)
    his_years.append(h_mean_list)
    
#%% mortality years

years = np.unique([x.coord('year').points[0] for x in damip_mor]) #damip and his years same

damip_mor_years =  iris.cube.CubeList()

for y in years:
    d_list = [x for x in damip_mor if x.coord('year').points[0] == y]

    
    damip_mor_years.append(d_list)


#%% Cause of difference between Histnat T1 and T2?
    
p_inlist = [pop_ratio_list[0], pop_ratio_list[1]]
m_inlist = [dmor_2000, dmor_2010]
e_inlist = [damip_years[0], damip_years[1]]

p_effect, m_effect, e_effect, plist, mlist, elist = calc_comp(p_inlist, m_inlist, e_inlist, return_list = True)
    
#%% Check against actual dif

#1995 - 2004: pop and mortality the same (onyl due to climate)
period1 = damip_mor_years[0]
period2 = damip_mor_years[1]

dif, dif_per = find_dif(period1, period2)

ens_dif = ens_mean(dif)
ens_perdif = ens_mean(dif_per)

period1_ens = ens_mean(period1)
period2_ens = ens_mean(period2)

#%%

print(np.nanmean(period1_ens.data))
print(np.nanmean(period2_ens.data))
print(np.nanmean(ens_dif.data))

np.nanmean(p_effect.data)
np.nanmean(m_effect.data)
np.nanmean(e_effect.data)

np.nanmean(p_effect.data) +  np.nanmean(m_effect.data) +   np.nanmean(e_effect.data) 


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



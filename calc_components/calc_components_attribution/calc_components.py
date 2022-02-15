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
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib


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

fs = 16
matplotlib.rcParams.update({'font.size': fs})

#Import my functions
import sys
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/Tanga/Plot_functions')
import tanzania1 as tp

def ens_mean(cube_list):  
    cor_cubes = iris.cube.CubeList()
    #add ensemble as dimension so can use cube statistics
    for i in np.arange(len(cube_list)):
        item = cube_list[i]
        gcm = tp.gcm(item)
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
        
        p_effect = copy.deepcopy(e1_mlist[i])
        p_effect.data = (p1 - p2)*(peq1 + peq2)
        
        meq1 = (p1 * e1 + p2 * e2)  / 3
        meq2 = (p1 * e2 + p2 * e1) / 6
                
        m_effect = copy.deepcopy(e1_mlist[i])
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



def calc_comp_num(p_list, m_list, e_list):
    '''order should be:
    period1, period2
    list of models within each period only for e (pop and mor same each model)
    '''
    
    p1 = p_list[0]
    p2 = p_list[1]
    m1 = m_list[0]
    m2 = m_list[1]
    e1 = e_list[0]
    e2 = e_list[1]
    
        
    peq1 = (m1 * e1 + m2 * e2)  / 3
    peq2 = (m1 * e2 + m2 * e1) / 6
    
    p_effect = (p1 - p2)*(peq1 + peq2)
    
    meq1 = (p1 * e1 + p2 * e2)  / 3
    meq2 = (p1 * e2 + p2 * e1) / 6
            
    m_effect = (m1 - m2)*(meq1 + meq2)
    
    eeq1 = (m1 * p1 + m2 * p2)  / 3
    eeq2 = (p1 * m2 + p2 * m1) / 6
            
    e_effect = (e1 - e2)*(eeq1 + eeq2)
            
    return p_effect, m_effect, e_effect
    #return p_ens, m_ens, e_ens, p_list, m_list, e_list

#%%

p_list = [5, 10]
m_list = [8, 6]
e_list = [2, 4]

p, m, e = calc_comp_num(p_list, m_list, e_list)

print(p, m, e)


p_out_list = []
m_out_list = []
e_out_list = []

for i in np.arange(5, 20):
    p_list = [5, i]
    p, m, e = calc_comp_num(p_list, m_list, e_list)
    p_out_list.append(p * -1)
    m_out_list.append(m * -1)
    e_out_list.append(e * -1)


plt.plot(np.arange(5, 20), p_out_list, label = 'p effect')
plt.plot(np.arange(5, 20), m_out_list, label = 'm effect')
plt.plot(np.arange(5, 20), e_out_list, label = 'e effect')
plt.legend()
plt.xlabel('Value P in 2nd time period (m, e change held constant)')
plt.ylabel('Effect size')

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
        
path = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_1/historical/'

#hist-nat
filenames = glob.glob(path + '*hist-nat*.nc')
damip_list1 = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    damip_list1.append(x)


#historical
filenames = glob.glob(path + '*historical*.nc')
hist_list1 = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    if tp.gcm(x) in damip_mods:
        hist_list1.append(x)
  
#%% Import pop and mor data 

#population -from world pop
# 0 - 4 year olds
pop_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000_regrid.nc')
pop_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2010_regrid.nc')
pop_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2019_regrid.nc')

pop_2015 = pop_2010 + ((pop_2019 - pop_2010) / 2)


#daily mortality
dmor_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2000_regrid.nc')
dmor_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2010_regrid.nc')
dmor_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2019_regrid.nc')

dmor_2015 = dmor_2010 + ((dmor_2019 - dmor_2010) / 2)


pop_list = [pop_2000, pop_2010, pop_2015, pop_2019]
mor_list = [dmor_2000, dmor_2010, dmor_2015, dmor_2019]

#%% actual mortality

path = '/nfs/a321/earsch/CHAMNHA/output/annual_mortality/coeff_061/thres_hismodel/'

filenames = glob.glob(path + '*hist-nat*')
damip_mor = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    damip_mor.append(x)
    
filenames = glob.glob(path + '*historical*')
his_mor = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    if tp.gcm(x) in damip_mods:
       his_mor.append(x)
       
path = '/nfs/a321/earsch/CHAMNHA/output/annual_mortality/coeff_1/thres_hismodel/'

filenames = glob.glob(path + '*hist-nat*')
damip_mor1 = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    damip_mor1.append(x)
    
filenames = glob.glob(path + '*historical*')
his_mor1 = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    if tp.gcm(x) in damip_mods:
       his_mor1.append(x)
       
#x = [tp.gcm(y) for y in his_mor1]

#%% restrict to afric only
countries = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/africa_countires_only.nc')


country_lookup = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/pop/country_lookup.csv')
country_lookup['Name'][country_lookup['Name'] == 'Republic of the Congo'] = 'Congo'
country_lookup['Name'][country_lookup['Name'] == 'eSwatini'] = 'Swaziland'


not_in_afr = ['Albania', 'Armenia', 'Azerbaijan', 'Cyprus', 'France', 'Greece', 
              'Iran', 'Iraq', 'Israel', 'Italy', 'Jordan', 'Kuwait', 'Lebanon',
              'Malta', 'Northern Cyprus', 'Oman', 'Palestine', 'Portugal', 'Qatar',
              'Saudi Arabia', 'Spain', 'Syria', 'Turkey', 'Turkmenistan', 'United Arab Emirates', 'Yemen']

c_dict = dict(zip(country_lookup['Value'], country_lookup['Name']))

#%% restrict africa only

#obs - used to get mask
cru_tas = iris.load('/nfs/a321/earsch/Tanga/Data/CRU/tmp/*.nc',
                    iris.Constraint(cube_func = lambda cube: cube.var_name == 'tmp'))

cru_tas = cru_tas[0][0]
cru_tas = cru_tas.regrid(damip_list[0], iris.analysis.Linear())

#%%

for cube in damip_list:
    for i in np.arange(cube.shape[0]):
        cube.data[i][np.isnan(countries.data)] = np.nan
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)
    
for cube in hist_list:
    for i in np.arange(cube.shape[0]):
        cube.data[i][np.isnan(countries.data)] = np.nan
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)

for cube in damip_list1:
    for i in np.arange(cube.shape[0]):
        cube.data[i][np.isnan(countries.data)] = np.nan
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)
        
for cube in hist_list1:
    for i in np.arange(cube.shape[0]):
        cube.data[i][np.isnan(countries.data)] = np.nan
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)


for cube in damip_mor:
    for i in np.arange(cube.shape[0]):
        cube.data[i][np.isnan(countries.data)] = np.nan
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)
    
for cube in his_mor:
    for i in np.arange(cube.shape[0]):
        cube.data[i][np.isnan(countries.data)] = np.nan
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)

for cube in damip_mor1:
    for i in np.arange(cube.shape[0]):
        cube.data[i][np.isnan(countries.data)] = np.nan
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)
    
for cube in his_mor1:
    for i in np.arange(cube.shape[0]):
        cube.data[i][np.isnan(countries.data)] = np.nan
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)
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
        #will compare 2011 - 2020 to 1995 - 2004

def split_years(in_list):
    hyears = [np.arange(1995, 2005),
              np.arange(2011, 2021)]

    his_years = iris.cube.CubeList()
    
    mods = np.unique([tp.gcm(x) for x in in_list])
    
    for year_group in hyears:
        if 1995 in year_group:
            h_list = [x for x in in_list if x.coord('year').points[0] == 1995]
            
        else:
            #2011 - 2020 are split over 2 cube lists
            h_list0 = [x for x in in_list if x.coord('year').points[0] == 2005]
            h_list1 = [x for x in in_list if x.coord('year').points[0] == 2015]
            h_list = iris.cube.CubeList()
            
            for mod in mods:
                h_list0m = [y for y in h_list0 if tp.gcm(y) == mod][0]
                h_list1m = [y for y in h_list1 if tp.gcm(y) == mod]
                
                if len(h_list1m) > 0:
                    #future mods only don't have 2015 - 2020 historical period, as that's actually climate change
                    h_list1m = h_list1m[0]
                    x0 = h_list0m.extract(iris.Constraint(year = lambda cell: cell >= 2011))
                    x1 = h_list1m.extract(iris.Constraint(year = lambda cell: cell <= 2020))
                    
                    x = iris.cube.CubeList([x0,x1])
                    x = x.concatenate_cube()
                    
                    h_list.append(x)
                else:
                    x0 = h_list0m.extract(iris.Constraint(year = lambda cell: cell >= 2011))
                    h_list.append(x0)
            
        h_mean_list = iris.cube.CubeList()
        for cube in h_list:
            h_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))
        
        his_years.append(h_mean_list)
    return his_years


damip_years = split_years(damip_list)
his_years = split_years(hist_list)
       
#coeff 1
damip_years1 =  split_years(damip_list1)
his_years1 = split_years(hist_list1)


    
#%% mortality years

def split_years_mor(in_list):
    hyears = [np.arange(1995, 2005),
              np.arange(2011, 2021)]
    his_mor_years = iris.cube.CubeList()
    
    mods = np.unique([tp.gcm(x) for x in in_list])
    
    for year_group in hyears:
        if 1995 in year_group:
            h_list = [x for x in in_list if x.coord('year').points[0] == 1995]
            
        else:
            #2011 - 2020 are split over 2 cube lists
            h_list0 = [x for x in in_list if x.coord('year').points[0] == 2005]
            h_list1 = [x for x in in_list if x.coord('year').points[0] == 2015]
            h_list = iris.cube.CubeList()
            
            for mod in mods:
                h_list0m = [y for y in h_list0 if tp.gcm(y) == mod][0]
                h_list1m = [y for y in h_list1 if tp.gcm(y) == mod]
                
                if len(h_list1m) > 0:
                    #future mods only don't have 2015 - 2020 historical period, as that's actually climate change
                    h_list1m = h_list1m[0]
                    x0 = h_list0m.extract(iris.Constraint(year = lambda cell: cell >= 2011))
                    x1 = h_list1m.extract(iris.Constraint(year = lambda cell: cell <= 2020))
                    
                    x = iris.cube.CubeList([x0,x1])
                    x = x.concatenate_cube()
                    
                    h_list.append(x)
                else:
                    x0 = h_list0m.extract(iris.Constraint(year = lambda cell: cell >= 2011))
                    h_list.append(x0)
            
        his_mor_years.append(h_list)
    return his_mor_years


damip_mor_years =  split_years_mor(damip_mor)
his_mor_years = split_years_mor(his_mor)

damip_mor_years1 =  split_years_mor(damip_mor1)
his_mor_years1 = split_years_mor(his_mor1)



#%% Cause of difference between Histnat T1 and T2?
#T1= 1995 - 2004
#T2 = 2011 - 2020    


p_inlist = [pop_ratio_list[0], pop_ratio_list[2]]
m_inlist = [dmor_2000, dmor_2015]
e_inlist = [damip_years[0], damip_years[1]]
e_his_inlist = [his_years[0], his_years[1]]

p_effect, m_effect, e_effect, plist, mlist, elist = calc_comp(p_inlist, m_inlist, e_inlist, return_list = True)
#his T2 vs T1
p_effect_h, m_effect_h, e_effect_h, plist_h, mlist_h, elist_h = calc_comp(p_inlist, m_inlist, e_his_inlist, return_list = True)
 

#coeff 1
   
p_inlist = [pop_ratio_list[0], pop_ratio_list[2]]
m_inlist = [dmor_2000, dmor_2015]
e_inlist = [damip_years1[0], damip_years1[1]]
e_his_inlist = [his_years1[0], his_years1[1]]

p_effect1, m_effect1, e_effect1, plist1, mlist1, elist1 = calc_comp(p_inlist, m_inlist, e_inlist, return_list = True)
#his T2 vs T1
p_effect_h1, m_effect_h1, e_effect_h1, plist_h1, mlist_h1, elist_h1 = calc_comp(p_inlist, m_inlist, e_his_inlist, return_list = True)
  

  
#%% Check against actual dif

#1995 - 2004: pop and mortality the same (onyl due to climate)
period1 = damip_mor_years[0]
period2 = damip_mor_years[1]
period1 = [x.collapsed(['time'], iris.analysis.MEAN) for x in period1]
period2 = [x.collapsed(['time'], iris.analysis.MEAN) for x in period2]


dif, dif_per = find_dif(period1, period2)

ens_dif = ens_mean(dif)
ens_perdif = ens_mean(dif_per)

period1_ens = ens_mean(period1)
period2_ens = ens_mean(period2)


#historical difs

period1_h = his_mor_years[0]
period2_h = his_mor_years[1]
period1_h = [x.collapsed(['time'], iris.analysis.MEAN) for x in period1_h]
period2_h = [x.collapsed(['time'], iris.analysis.MEAN) for x in period2_h]

dif_h, dif_per_h = find_dif(period1_h, period2_h)

ens_dif_h = ens_mean(dif_h)
ens_perdif_h = ens_mean(dif_per_h)

period1_ens_h = ens_mean(period1_h)
period2_ens_h = ens_mean(period2_h)


# coeff 1


#1995 - 2004: pop and mortality the same (onyl due to climate)
period11 = damip_mor_years1[0]
period21 = damip_mor_years1[1]
period11 = [x.collapsed(['time'], iris.analysis.MEAN) for x in period11]
period21 = [x.collapsed(['time'], iris.analysis.MEAN) for x in period21]

dif1, dif_per1 = find_dif(period11, period21)

ens_dif1 = ens_mean(dif1)
ens_perdif1 = ens_mean(dif_per1)

period1_ens1 = ens_mean(period11)
period2_ens1 = ens_mean(period21)


#historical difs

period1_h1 = his_mor_years1[0]
period2_h1 = his_mor_years1[1]
period1_h1 = [x.collapsed(['time'], iris.analysis.MEAN) for x in period1_h1]
period2_h1 = [x.collapsed(['time'], iris.analysis.MEAN) for x in period2_h1]

dif_h1, dif_per_h1 = find_dif(period1_h1, period2_h1)

ens_dif_h1 = ens_mean(dif_h1)
ens_perdif_h1 = ens_mean(dif_per_h1)

period1_ens_h1 = ens_mean(period1_h1)
period2_ens_h1 = ens_mean(period2_h1)



#%% Check decomp

print(np.nanmean(period1_ens.data))
print(np.nanmean(period2_ens.data))
print(np.nanmean(ens_dif.data))

np.nanmean(p_effect.data)
np.nanmean(m_effect.data)
np.nanmean(e_effect.data)

np.nanmean(p_effect.data) +  np.nanmean(m_effect.data) +   np.nanmean(e_effect.data) 
 

#historical

print(np.nanmean(ens_dif_h.data))
np.nanmean(p_effect_h.data) +  np.nanmean(m_effect_h.data) +   np.nanmean(e_effect_h.data) 

#%%check individual historical models - order mixed up


np.nanmean(p_effect_h.data)
np.nanmean(m_effect_h.data)
np.nanmean(e_effect_h.data)

print(np.nanmean(dif_h[2].data))
np.nanmean(plist_h[0].data) +  np.nanmean(mlist_h[0].data) +   np.nanmean(elist_h3[0].data) 


#%% What percentage of the change is due to climate, pop and mortality?

def cube_to_frame(p_effect_list, m_effect_list, e_effect_list, total_dif, period, coeff):
    
    df = pd.DataFrame(columns = ['model', 'period', 'coeff', 'p', 'm', 'e', 'p_per', 'm_per', 'e_per', 'total'])

    for i in np.arange(len(p_effect_list)):
        p = p_effect_list[i]
        gcm = tp.gcm(p)
        
        m = [x for x in m_effect_list if x.coord('ens').points[0] == gcm][0]
        e = [x for x in e_effect_list if x.coord('ens').points[0] == gcm][0]
        d = [x for x in total_dif if x.coord('ens').points[0] == gcm][0]
        
        p_val =  np.nansum(p.data)
        m_val =  np.nansum(m.data)
        e_val =  np.nansum(e.data)
        d_val =  np.nansum(d.data)
        
        p_per = p_val / d_val * 100
        m_per = m_val / d_val * 100
        e_per = e_val / d_val * 100
    
            
        y = pd.DataFrame(data = {'model': gcm,
                                 'period': period,
                                 'coeff': coeff,
                                 'p': p_val,
                                 'm': m_val,
                                 'e': e_val,
                                 'p_per': p_per,
                                 'm_per': m_per,
                                 'e_per': e_per,
                                 'total': p_per + m_per + e_per}, index = [0])

        
        df = df.append(y)
    
    return df


damip_p2 = cube_to_frame(plist, mlist, elist, dif, 'period2', '0.61')
his_p2 = cube_to_frame(plist_h, mlist_h, elist_h, dif_h, 'period2', '0.61')

damip_p21 = cube_to_frame(plist1, mlist1, elist1, dif1, 'period2', '1.0')
his_p21 = cube_to_frame(plist_h1, mlist_h1, elist_h1, dif_h1, 'period2', '1.0')


damip = damip_p2.append(damip_p21)
his = his_p2.append(his_p21)


damip_agg = damip.groupby(['period', 'coeff'], as_index = True).mean()
his_agg = his.groupby(['period', 'coeff'], as_index = True).mean()

d_per = damip_agg.drop(columns = ['p' , 'm', 'e', 'total'])
h_per = his_agg.drop(columns = ['p' , 'm', 'e', 'total'])

d_abs = damip_agg.drop(columns = ['p_per' , 'm_per', 'e_per', 'total'])
h_abs = his_agg.drop(columns = ['p_per' , 'm_per', 'e_per', 'total'])

#otherwise looks like population growth makes mortlaity go down
d_abs = d_abs * -1
h_abs = h_abs * -1


#d_abs.to_csv('/nfs/see-fs-02_users/earsch/Documents/Leeds/CHAMNHA/output_csvs/d_abs.csv')
#h_abs.to_csv('/nfs/see-fs-02_users/earsch/Documents/Leeds/CHAMNHA/output_csvs/h_abs.csv')

#%% https://stackoverflow.com/questions/22787209/how-to-have-clusters-of-stacked-bars-with-python-pandas

def plot_clustered_stacked(dfall, cols, labels=None, title="",  H="/", width = 0.5, z = 1, **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe"""

    n_df = len(dfall)
    n_col = len(dfall[0].columns) 
    n_ind = len(dfall[0].index)
    axe = plt.subplot(1,2, z)

    for df in dfall : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ax=axe,
                      legend=False,
                      grid=False,
                      color = cols,
                      **kwargs)  # make bar plots

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + width / float(n_df + width) * i / float(n_col))
                rect.set_hatch(H * int(i / n_col)) #edited part     
                rect.set_width(width / float(n_df + width))

    axe.set_xticks(np.arange(n_ind))
    axe.set_xticklabels(df.index, rotation = 0)
    #axe.set_title(title)

    # Add invisible data to add another legend
    #n=[]        
    #for i in range(n_df):
    #    n.append(axe.bar(0, 0, color="gray", hatch=H * i))

    #l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0.5])
    #if labels is not None:
     #   l2 = plt.legend(n, labels, loc=[1.01, 0.1]) 
    #axe.add_artist(l1)
    return axe


#%% Plot absoltue mortality
    
h_abs_low = h_abs.iloc[h_abs.index.get_level_values('coeff') == '0.61']
h_abs_high = h_abs.iloc[h_abs.index.get_level_values('coeff') == '1.0']
d_abs_low = d_abs.iloc[d_abs.index.get_level_values('coeff') == '0.61']
d_abs_high = d_abs.iloc[d_abs.index.get_level_values('coeff') == '1.0']

fig = plt.figure(figsize=(12,6))

ax = plot_clustered_stacked([h_abs_low, d_abs_low],labels = ["his", "damip"],
                            cols = ['#1b9e77', '#7570b3', '#d95f02'])
ax2 = plot_clustered_stacked([h_abs_high, d_abs_high],labels = ["his", "damip"],
                            cols = ['#1b9e77', '#7570b3', '#d95f02'], z = 2)                                    
                                    

ax_list = [ax, ax2]
for a in ax_list:
    
    a.set_ylim([-6000, 13000])
    a.set_xlabel('Period')
    a.set_xticklabels(['2005 - 2014', '2015 - 2020'])
    a.axhline(0, c = 'k', linewidth = 0.5)

ax2.set_yticklabels([''] * 6)

ax.set_ylabel('Contribution to change in heat mortality')


plt.draw()
x, y = 0.4,1.05

ax.annotate('Coeff: 0.61', (x, y), xycoords = 'axes fraction',
             rotation = 0, ha = 'left')
ax2.annotate('Coeff: 1.00', (x, y), xycoords = 'axes fraction',
             rotation = 0, ha = 'left')


mor_patch = Patch(facecolor='#7570b3')
pop_patch = Patch(facecolor='#1b9e77')
e_patch = Patch(facecolor='#d95f02')
his_patch = Patch(facecolor='grey')
damip_patch = Patch(facecolor='grey', hatch = '//')

custom_lines = [pop_patch, his_patch, mor_patch, damip_patch, e_patch]


ax.legend(custom_lines, ['Population', 'Historical', 'Mortality', 
                         'Hist-nat', 'Climate'], 
           bbox_to_anchor=(1.9, -0.2),ncol=3,frameon = False, handletextpad = 0.5)


#have to run plt.savefig at same time as create figure to get it to save correctly

#plt.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/decomp_abs_total_bothcoeffs.png',
#            bbox_inches = 'tight', pad_inches = 0.3)

#%% Plot per mortality


ax = plot_clustered_stacked([h_per, d_per],labels = ["his", "damip"], cols = ['#1b9e77', '#7570b3', '#d95f02'])

plt.xlabel('Period')

ax.set_ylabel('% contribution to change in heat mortality')
ax.set_xticklabels(['2005 - 2014', '2015 - 2020'])

ax.axhline(0, c = 'k', linewidth = 0.5)



mor_patch = Patch(facecolor='#7570b3')
pop_patch = Patch(facecolor='#1b9e77')
e_patch = Patch(facecolor='#d95f02')
his_patch = Patch(facecolor='grey')
damip_patch = Patch(facecolor='grey', hatch = '//')

custom_lines = [pop_patch, his_patch, mor_patch, damip_patch, e_patch]


ax.legend(custom_lines, ['Population', 'Historical', 'Mortality', 
                         'Hist-nat', 'Climate'], 
           bbox_to_anchor=(0.9, -0.2),ncol=3,frameon = False, handletextpad = 0.5)


#plt.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/decomp_per.png',
#            bbox_inches = 'tight', pad_inches = 0.3)
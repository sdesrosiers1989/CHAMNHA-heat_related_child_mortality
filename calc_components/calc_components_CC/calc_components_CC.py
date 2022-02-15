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

def find_dif(his_list, fut_list):
    #fut - his
    dif_list = []
    difper_list = []
    
    for cube in fut_list:
        mod = tp.gcm(cube)
        comp_cube =  [x for x in his_list if tp.gcm(x) == mod][0] # checked, len always 1, as it should be
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
        
        p_effect = copy.deepcopy(e2_mlist[i])
        p_effect.data = (p1 - p2)*(peq1 + peq2)
        
        meq1 = (p1 * e1 + p2 * e2)  / 3
        meq2 = (p1 * e2 + p2 * e1) / 6
                
        m_effect = copy.deepcopy(e2_mlist[i])
        m_effect.data = (m1 - m2)*(meq1 + meq2)
        
        eeq1 = (m1 * p1 + m2 * p2)  / 3
        eeq2 = (p1 * m2 + p2 * m1) / 6
                
        e_effect = copy.deepcopy(e2_mlist[i])
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



#%% import e model data
    
path = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_061/historical/'

#historical
filenames = glob.glob(path + '*historical*.nc')
hist_list = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    hist_list.append(x)
        
path = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_1/historical/'

#historical
filenames = glob.glob(path + '*historical*.nc')

hist_list1 = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    hist_list1.append(x)
    
#future
path = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_061/future/'

#historical
filenames = glob.glob(path + '*PBC_MBC.nc')
fut_list = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    fut_list.append(x)
        
path = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_1/future/'

#historical
filenames = glob.glob(path + '*PBC_MBC.nc')
fut_list1 = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    fut_list1.append(x)
    
 #%% actual mortality

path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/'

filenames = glob.glob(path + '*historical*')

his = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    his.append(x)
    
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_1/thres_hismodel/'

filenames = glob.glob(path + '*historical*')

his1 = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    his1.append(x)
    
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/future/'

filenames = glob.glob(path + '*PBC_MBC.nc')
fut = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    fut.append(x)


path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_1/thres_hismodel/future/'

filenames = glob.glob(path + '*PBC_MBC.nc')
fut1 = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    fut1.append(x)
  
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

#future more and pop

years = np.arange(2025, 2050, 10)
pop_path = '/nfs/a321/earsch/CHAMNHA/input_data/pop/future/processed/'
mor_path = '/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/processed/'
mor_list_fut = []

for y in years:
    p_name = pop_path + 'ssp2_' + str(y) + '_04population_mf_BIASCORR.nc'
    m_name = mor_path + 'ref_' + str(y) + '_04_totalmor_mf_BIASCORR.nc'
    
    pop_list.append(iris.load_cube(p_name))
    
    
    mor_list_fut.append(iris.load_cube(m_name))


mor_list_fut = [x/365 for x in mor_list_fut]


#%% calc pop ratio as sued in health burden model
#pop data

pop_ratio_list = []
for i in np.arange(len(pop_list)):

    pop_ratio = pop_list[i]/pop_2000 #ratio future pop to baseline pop   
    #0 in denominator causing problems
    pop_ratio.data = np.where(pop_2000.data == 0, 1, pop_ratio.data)
    pop_ratio.data = ma.masked_array(pop_ratio.data, mask = pop_2000.data.mask)  
    pop_ratio_list.append(pop_ratio)
    
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
cru_tas = cru_tas.regrid(fut_list[0], iris.analysis.Linear())
    

#%%

for cube in fut_list:
    for i in np.arange(cube.shape[0]):
        cube.data.mask[i][np.isnan(countries.data)] = True
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)
    
for cube in hist_list:
    for i in np.arange(cube.shape[0]):
        cube.data.mask[i][np.isnan(countries.data)] = True
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)

for cube in fut_list1:
    for i in np.arange(cube.shape[0]):
        cube.data.mask[i][np.isnan(countries.data)] = True
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)
        
for cube in hist_list1:
    for i in np.arange(cube.shape[0]):
        cube.data.mask[i][np.isnan(countries.data)] = True
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)

for cube in his:
    cube.data.mask[np.isnan(countries.data)] = True
    cube.data = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube.data)
    
for cube in his1:
    cube.data.mask[np.isnan(countries.data)] = True
    cube.data = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube.data)

for cube in fut:
    cube.data.mask[np.isnan(countries.data)] = True
    cube.data = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube.data)
    
for cube in fut1:
    cube.data.mask[np.isnan(countries.data)] = True
    cube.data = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube.data)

#%% get mdoels into year lists
    
hyears = np.unique([x.coord('year').points[0] for x in hist_list]) 

his_years = iris.cube.CubeList()

for y in hyears:
    h_list = [x for x in hist_list if x.coord('year').points[0] == y]
    h_mean_list = iris.cube.CubeList()
    for cube in h_list:
        h_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))
    
    his_years.append(h_mean_list)
    
 
    
#coeff 1
his_years1 = iris.cube.CubeList()

for y in hyears:
 

    h_list = [x for x in hist_list1 if x.coord('year').points[0] == y]
    h_mean_list = iris.cube.CubeList()
    for cube in h_list:
        h_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))
    

    his_years1.append(h_mean_list)
    
    
    
# future
    
years = np.unique([x.coord('year').points[0] for x in fut_list]) 

fut_years = iris.cube.CubeList()

for y in years:
    f_list = [x for x in fut_list if x.coord('year').points[0] == y]
    f_mean_list = iris.cube.CubeList()
    for cube in f_list:
        f_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))
    
    fut_years.append(f_mean_list)
    
 
    
#coeff 1
fut_years1 = iris.cube.CubeList()

for y in years:
    f_list = [x for x in fut_list1 if x.coord('year').points[0] == y]
    f_mean_list = iris.cube.CubeList()
    for cube in f_list:
        f_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))
    
    fut_years1.append(f_mean_list)

#%% actual mortality years
    
years = np.unique([x.coord('year').points[0] for x in his]) 

his_mor_years = iris.cube.CubeList()

for y in years:
    
    h_list = [x for x in his if x.coord('year').points[0] == y]
    his_mor_years.append(h_list)

his_mor_years1 = iris.cube.CubeList()

for y in years:

    
    h_list = [x for x in his1 if x.coord('year').points[0] == y]
    his_mor_years1.append(h_list)


#fut
years = np.unique([x.coord('year').points[0] for x in fut]) 

fut_mor_years = iris.cube.CubeList()

for y in years:
    
    f_list = [x for x in fut if x.coord('year').points[0] == y]
    fut_mor_years.append(f_list)

fut_mor_years1 = iris.cube.CubeList()

for y in years:

    
    f_list = [x for x in fut1 if x.coord('year').points[0] == y]
    fut_mor_years1.append(f_list)


#%% Cause of difference between Hist  and future?
# comp_period = 2005 - 2014
    
i_h = 1 # index corresponding to 2005 - 2014

#compare 2020, 2030, 2040

#2020
p_inlist = [pop_ratio_list[i_h], pop_ratio_list[3]]
m_inlist = [dmor_2010, mor_list_fut[0]]
e_inlist = [his_years[i_h], fut_years[0]]

p_effect, m_effect, e_effect, plist, mlist, elist = calc_comp(p_inlist, m_inlist, e_inlist, return_list = True)

#2030
p_inlist = [pop_ratio_list[i_h], pop_ratio_list[4]]
m_inlist = [dmor_2010, mor_list_fut[1]]
e_inlist = [his_years[i_h], fut_years[1]]

p_effect3, m_effect3, e_effect3, plist3, mlist3, elist3 = calc_comp(p_inlist, m_inlist, e_inlist, return_list = True)
  
#2040

p_inlist = [pop_ratio_list[i_h], pop_ratio_list[5]]
m_inlist = [dmor_2010, mor_list_fut[2]]
e_inlist = [his_years[i_h], fut_years[2]]

p_effect4, m_effect4, e_effect4, plist4, mlist4, elist4 = calc_comp(p_inlist, m_inlist, e_inlist, return_list = True)
  

#coeff 1
   
#2020
p_inlist = [pop_ratio_list[i_h], pop_ratio_list[3]]
m_inlist = [dmor_2010, mor_list_fut[0]]
e_inlist = [his_years1[i_h], fut_years1[0]]

p_effect1, m_effect1, e_effect1, plist1, mlist1, elist1 = calc_comp(p_inlist, m_inlist, e_inlist, return_list = True)

#2030
p_inlist = [pop_ratio_list[i_h], pop_ratio_list[4]]
m_inlist = [dmor_2010, mor_list_fut[1]]
e_inlist = [his_years1[i_h], fut_years1[1]]

p_effect31, m_effect31, e_effect31, plist31, mlist31, elist31 = calc_comp(p_inlist, m_inlist, e_inlist, return_list = True)
  
#2040

p_inlist = [pop_ratio_list[i_h], pop_ratio_list[5]]
m_inlist = [dmor_2010, mor_list_fut[2]]
e_inlist = [his_years1[i_h], fut_years1[2]]

p_effect41, m_effect41, e_effect41, plist41, mlist41, elist41 = calc_comp(p_inlist, m_inlist, e_inlist, return_list = True)
  
#%% Check values

#1995 - 2004: his
base_period = his_mor_years[i_h]
fperiod = fut_mor_years[0]
fperiod3 = fut_mor_years[1]
fperiod4 = fut_mor_years[2]

dif, dif_per = find_dif(base_period, fperiod)
dif3, dif_per3 = find_dif(base_period, fperiod3)
dif4, dif_per4 = find_dif(base_period, fperiod4)

ens_dif = ens_mean(dif)
ens_perdif = ens_mean(dif_per)
ens_dif3 = ens_mean(dif3)
ens_perdif3 = ens_mean(dif_per3)
ens_dif4 = ens_mean(dif4)
ens_perdif4 = ens_mean(dif_per4)

bperiod_ens = ens_mean(base_period)
fperiod_ens = ens_mean(fperiod)
fperiod3_ens = ens_mean(fperiod3)
fperiod4_ens = ens_mean(fperiod4)





# coeff 1

base_period1 = his_mor_years1[i_h]
fperiod1 = fut_mor_years1[0]
fperiod31 = fut_mor_years1[1]
fperiod41 = fut_mor_years1[2]

dif1, dif_per1 = find_dif(base_period1, fperiod1)
dif31, dif_per31 = find_dif(base_period1, fperiod31)
dif41, dif_per41 = find_dif(base_period1, fperiod41)

ens_dif1 = ens_mean(dif1)
ens_perdif1 = ens_mean(dif_per1)
ens_dif31 = ens_mean(dif31)
ens_perdif31 = ens_mean(dif_per31)
ens_dif41 = ens_mean(dif41)
ens_perdif41 = ens_mean(dif_per41)

bperiod_ens1 = ens_mean(base_period1)
fperiod_ens1 = ens_mean(fperiod1)
fperiod3_ens1 = ens_mean(fperiod31)
fperiod4_ens1 = ens_mean(fperiod41)

#%% Check decomp

print(np.nanmean(bperiod_ens.data))
print(np.nanmean(fperiod_ens.data))
print(np.nanmean(ens_dif.data))

np.nanmean(p_effect.data)
np.nanmean(m_effect.data)
np.nanmean(e_effect.data)

np.nanmean(p_effect.data) +  np.nanmean(m_effect.data) +   np.nanmean(e_effect.data) 

#coeff 1
print(np.nanmean(bperiod_ens1.data))
print(np.nanmean(fperiod_ens1.data))
print(np.nanmean(ens_dif1.data))

np.nanmean(p_effect1.data)
np.nanmean(m_effect1.data)
np.nanmean(e_effect1.data)

np.nanmean(p_effect1.data) +  np.nanmean(m_effect1.data) +   np.nanmean(e_effect1.data) 


#%% What percentage of the change is due to climate, pop and mortality?

def cube_to_frame(p_effect_list, m_effect_list, e_effect_list, total_dif, period, coeff):
    
    df = pd.DataFrame(columns = ['model', 'scen', 'period', 'coeff', 'p', 'm', 'e', 'p_per', 'm_per', 'e_per', 'total'])

    scens = np.unique([x.coord('sim').points[0] for x in p_effect_list])
    
    for s in scens:
        p_effect_list_s = [x for x in p_effect_list if x.coord('sim').points[0] == s]
        m_effect_list_s = [x for x in m_effect_list if x.coord('sim').points[0] == s]
        e_effect_list_s = [x for x in e_effect_list if x.coord('sim').points[0] == s]
        total_dif_s = [x for x in total_dif if x.coord('sim').points[0] == s]

        for i in np.arange(len(p_effect_list_s)):
            p = p_effect_list_s[i]
            gcm = tp.gcm(p)
            
            m = [x for x in m_effect_list_s if x.coord('ens').points[0] == gcm][0]
            e = [x for x in e_effect_list_s if x.coord('ens').points[0] == gcm][0]
            d = [x for x in total_dif_s if x.coord('ens').points[0] == gcm][0]
            
                        
            p_val =  np.nansum(p.data)
            m_val =  np.nansum(m.data)
            e_val =  np.nansum(e.data)
            d_val =  np.nansum(d.data)
            
            p_per = p_val / d_val * 100
            m_per = m_val / d_val * 100
            e_per = e_val / d_val * 100
        
                
            y = pd.DataFrame(data = {'model': gcm,
                                     'scen': s,
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


df_p = cube_to_frame(plist, mlist, elist, dif, '2020s', '0.61')
df_p3 = cube_to_frame(plist3, mlist3, elist3, dif3, '2030s', '0.61')
df_p4 = cube_to_frame(plist4, mlist4, elist4, dif4, '2040s', '0.61')

df_p1 = cube_to_frame(plist1, mlist1, elist1, dif1, '2020s', '1.0')
df_p31 = cube_to_frame(plist31, mlist31, elist31, dif31, '2030s', '1.0')
df_p41 = cube_to_frame(plist41, mlist41, elist41, dif41, '2040s', '1.0')


df_p = df_p.append(df_p1)
df_p3 = df_p3.append(df_p31)
df_p4 = df_p4.append(df_p41)

df = df_p.append(df_p3)
df = df.append(df_p4)


df_agg = df.groupby(['scen', 'period', 'coeff'], as_index = True).mean()


df_per = df_agg.drop(columns = ['p' , 'm', 'e', 'total'])

df_abs = df_agg.drop(columns = ['p_per' , 'm_per', 'e_per', 'total'])

#otherwise looks like population growth makes mortlaity go down
df_abs = df_abs * -1

print(df_agg.iloc[df_agg.index.get_level_values('coeff') == '0.61'])

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
    
    a.set_ylim([-6000, 15000])
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
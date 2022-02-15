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
from matplotlib.gridspec import GridSpec

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
    
    for i in np.arange(len(e2_mlist)):
        
        e2 = e2_mlist[i].data
        
        #find same gcm
        gcm = tp.gcm(e2_mlist[i])
        e1 = [x for x in e1_mlist if tp.gcm(x) == gcm][0].data
        
                
        
        
        
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
 
#future
path = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_1/future/'


filenames = glob.glob(path + '*PBC2_MBC.nc')
ssp119_list = iris.cube.CubeList()
ssp245_list = iris.cube.CubeList()
ssp585_list = iris.cube.CubeList()

for file in filenames:
    x = iris.load_cube(file)
    if 'ssp119' in file:
        ssp119_list.append(x)
    elif 'ssp245' in file:
        ssp245_list.append(x)
    elif 'ssp585' in file:
        ssp585_list.append(x)
        

mods = np.unique([tp.gcm(x) for x in ssp119_list])
    
    
    
path = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_1/historical/'

#historical
filenames = glob.glob(path + '*historical*.nc')
hist_list = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    if tp.gcm(x) in mods:
        hist_list.append(x)
        

ssp245_list = [x for x in ssp245_list if tp.gcm(x) in mods]    
ssp585_list = [x for x in ssp585_list if tp.gcm(x) in mods]    

        


    
 #%% actual mortality

path = '/nfs/a321/earsch/CHAMNHA/output/annual_mortality/coeff_1/thres_hismodel/'

filenames = glob.glob(path + '*historical*')

his = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    if tp.gcm(x) in mods:
        his.append(x)
    
    
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_1/thres_hismodel/future/'

filenames = glob.glob(path + '*PBC2_MBC.nc')
fut119_list = iris.cube.CubeList()
fut245_list = iris.cube.CubeList()
fut585_list = iris.cube.CubeList()

for file in filenames:
    x = iris.load_cube(file)
    if tp.gcm(x) in mods:

        if 'ssp119' in file:
            fut119_list.append(x)
        elif 'ssp245' in file:
            fut245_list.append(x)
        elif 'ssp585' in file:
            fut585_list.append(x)



  
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
    p_name = pop_path + 'ssp2_' + str(y) + '_04population_mf_BIASCORR2.nc'
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
cru_tas = cru_tas.regrid(ssp119_list[0], iris.analysis.Linear())
    

#%%

fut_list = [ssp119_list, ssp245_list, ssp585_list]

for cube_list in fut_list:
    for cube in cube_list:
        for i in np.arange(cube.shape[0]):
            cube.data.mask[i][np.isnan(countries.data)] = True
            cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)
    
for cube in hist_list:
    for i in np.arange(cube.shape[0]):
        cube.data.mask[i][np.isnan(countries.data)] = True
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)


for cube in his:
    for i in np.arange(cube.shape[0]):
        cube.data.mask[i][np.isnan(countries.data)] = True
        cube.data[i] = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube[i,:,:].data)
    
futmor_list = [fut119_list, fut245_list, fut585_list]

for cube_list in futmor_list:
    for cube in cube_list:
        cube.data.mask[np.isnan(countries.data)] = True
        cube.data = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube.data)


#%% get mdoels into year lists
    
hyears = [np.arange(1995, 2005)]

his_years = iris.cube.CubeList()

mods = np.unique([tp.gcm(x) for x in hist_list])

for year_group in hyears:
    if 1995 in year_group:
        h_list = [x for x in hist_list if x.coord('year').points[0] == 1995]
        
        
    h_mean_list = iris.cube.CubeList()
    for cube in h_list:
        h_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))
    
    his_years.append(h_mean_list)

 
        
    
    
# future
    
years = np.unique([x.coord('year').points[0] for x in ssp119_list]) 

ssp119_years = iris.cube.CubeList()

for y in years:
    f_list = [x for x in ssp119_list if x.coord('year').points[0] == y]
    f_mean_list = iris.cube.CubeList()
    for cube in f_list:
        f_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))
    
    ssp119_years.append(f_mean_list)
    

#ssp245
years = np.unique([x.coord('year').points[0] for x in ssp245_list]) 

ssp245_years = iris.cube.CubeList()

for y in years:
    f_list = [x for x in ssp245_list if x.coord('year').points[0] == y]
    f_mean_list = iris.cube.CubeList()
    for cube in f_list:
        f_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))
    
    ssp245_years.append(f_mean_list)
  
    
#585
years = np.unique([x.coord('year').points[0] for x in ssp585_list]) 

ssp585_years = iris.cube.CubeList()

for y in years:
    f_list = [x for x in ssp585_list if x.coord('year').points[0] == y]
    f_mean_list = iris.cube.CubeList()
    for cube in f_list:
        f_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))
    
    ssp585_years.append(f_mean_list)

#%% actual mortality years (need 1995 -2004 and 2011 - 2020 for historical)
    
years = np.unique([x.coord('year').points[0] for x in his]) 

hyears = [np.arange(1995, 2005)]
his_mor_years = iris.cube.CubeList()

mods = np.unique([tp.gcm(x) for x in his])

for year_group in hyears:
    if 1995 in year_group:
        h_list = [x for x in his if x.coord('year').points[0] == 1995]
        
    
        
    his_mor_years.append(h_list)



#fut
years = np.unique([x.coord('year').points[0] for x in fut119_list]) 

ssp119_mor_years = iris.cube.CubeList()

for y in years:
    
    f_list = [x for x in fut119_list if x.coord('year').points[0] == y]
    ssp119_mor_years.append(f_list)

years = np.unique([x.coord('year').points[0] for x in fut245_list]) 

ssp245_mor_years = iris.cube.CubeList()

for y in years:
    
    f_list = [x for x in fut245_list if x.coord('year').points[0] == y]
    ssp245_mor_years.append(f_list)
    
years = np.unique([x.coord('year').points[0] for x in fut585_list]) 

ssp585_mor_years = iris.cube.CubeList()

for y in years:
    
    f_list = [x for x in fut585_list if x.coord('year').points[0] == y]
    ssp585_mor_years.append(f_list)




#%% Cause of difference between Hist  and future?
# comparing to index 0 = 1995 - 2004
    
i_h = 0 # index corresponding to 1995 - 2004

#compare 2020, 2030, 2040

scen_list = [ssp119_years, ssp245_years, ssp585_years]

#2020
p_outlist = []
m_outlist = []
e_outlist = []

for i in np.arange(len(scen_list)):
    p_inlist = [pop_ratio_list[i_h], pop_ratio_list[3]]
    m_inlist = [dmor_2000, mor_list_fut[0]]
    e_inlist = [his_years[i_h], scen_list[i][0]]

    p_effect, m_effect, e_effect = calc_comp(p_inlist, m_inlist, e_inlist, return_list = False)
    
    p_outlist.append(p_effect)
    m_outlist.append(m_effect)
    e_outlist.append(e_effect)




#2030

p_outlist3 = []
m_outlist3 = []
e_outlist3 = []

for i in np.arange(len(scen_list)):
    p_inlist = [pop_ratio_list[i_h], pop_ratio_list[4]]
    m_inlist = [dmor_2000, mor_list_fut[1]]
    e_inlist = [his_years[i_h], scen_list[i][1]]

    p_effect, m_effect, e_effect = calc_comp(p_inlist, m_inlist, e_inlist, return_list = False)
    
    p_outlist3.append(p_effect)
    m_outlist3.append(m_effect)
    e_outlist3.append(e_effect)

#2040

p_outlist4 = []
m_outlist4 = []
e_outlist4 = []

for i in np.arange(len(scen_list)):
    p_inlist = [pop_ratio_list[i_h], pop_ratio_list[5]]
    m_inlist = [dmor_2000, mor_list_fut[2]]
    e_inlist = [his_years[i_h], scen_list[i][2]]

    p_effect, m_effect, e_effect = calc_comp(p_inlist, m_inlist, e_inlist, return_list = False)
    
    p_outlist4.append(p_effect)
    m_outlist4.append(m_effect)
    e_outlist4.append(e_effect)


#%% Check values

#ssp119
#1995 - 2004: his
base_period = his_mor_years[i_h]
base_period = [x.collapsed(['time'], iris.analysis.MEAN) for x in base_period]
fperiod = ssp119_mor_years[0]
fperiod3 = ssp119_mor_years[1]
fperiod4 = ssp119_mor_years[2]

#245
fperiod = ssp245_mor_years[0]
fperiod3 = ssp245_mor_years[1]
fperiod4 = ssp245_mor_years[2]

#585
fperiod = ssp585_mor_years[0]
fperiod3 = ssp585_mor_years[1]
fperiod4 = ssp585_mor_years[2]

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


#
base_period = his_mor_years[i_h]
base_period = [x.collapsed(['time'], iris.analysis.MEAN) for x in base_period]
scen_list = [ssp119_mor_years, ssp245_mor_years, ssp585_mor_years]
ens_dif_list = []
ens_dif3_list = []
ens_dif4_list = []
for i in np.arange(len(scen_list)):
    
    fperiod = scen_list[i][0]
    fperiod3 = scen_list[i][1]
    fperiod4 = scen_list[i][2]

    
    dif, dif_per = find_dif(base_period, fperiod)
    dif3, dif_per3 = find_dif(base_period, fperiod3)
    dif4, dif_per4 = find_dif(base_period, fperiod4)
    
    ens_dif_list.append(ens_mean(dif))
    ens_dif3_list.append(ens_mean(dif3))
    ens_dif4_list.append(ens_mean(dif4))


#%% Check decomp - run scenario specific ens from above first

print(np.nanmean(bperiod_ens.data))
print(np.nanmean(fperiod_ens.data))
print(np.nanmean(ens_dif.data))

np.nanmean(p_outlist[0].data)
np.nanmean(m_outlist[0].data)
np.nanmean(e_outlist[0].data)

#ssp119
print(np.nanmean(ens_dif.data))
np.nanmean(p_outlist[0].data) +  np.nanmean(m_outlist[0].data) +   np.nanmean(e_outlist[0].data) 
print(np.nanmean(ens_dif3.data))
np.nanmean(p_outlist3[0].data) +  np.nanmean(m_outlist3[0].data) +   np.nanmean(e_outlist3[0].data) 

#245
print(np.nanmean(ens_dif.data))
np.nanmean(p_outlist[1].data) +  np.nanmean(m_outlist[1].data) +   np.nanmean(e_outlist[1].data) 

#585
print(np.nanmean(ens_dif3.data))
np.nanmean(p_outlist3[2].data) +  np.nanmean(m_outlist3[2].data) +   np.nanmean(e_outlist3[2].data) 



#%% What percentage of the change is due to climate, pop and mortality?

def cube_to_frame(p_effect_list, m_effect_list, e_effect_list, total_dif, period):
    
    df = pd.DataFrame(columns = ['scen', 'period', 'p', 'm', 'e', 'p_per', 'm_per', 'e_per', 'total'])

    scens = np.unique([x.coord('sim').points[0] for x in p_effect_list])
    
    for s in scens:
        p = [x for x in p_effect_list if x.coord('sim').points[0] == s][0]
        m = [x for x in m_effect_list if x.coord('sim').points[0] == s][0]
        e = [x for x in e_effect_list if x.coord('sim').points[0] == s][0]
        d = [x for x in total_dif if x.coord('sim').points[0] == s][0]
                       
        p_val =  np.nansum(p.data)
        m_val =  np.nansum(m.data)
        e_val =  np.nansum(e.data)
        d_val =  np.nansum(d.data)
        
        p_per = p_val / d_val * 100
        m_per = m_val / d_val * 100
        e_per = e_val / d_val * 100
    
            
        y = pd.DataFrame(data = {'scen': s,
                                 'period': period,
                                 'p': p_val,
                                 'm': m_val,
                                 'e': e_val,
                                 'p_per': p_per,
                                 'm_per': m_per,
                                 'e_per': e_per,
                                 'total': p_per + m_per + e_per}, index = [0])

        
        df = df.append(y)
    
    return df


df_p = cube_to_frame(p_outlist, m_outlist, e_outlist, ens_dif_list, '2020s')
df_p3 = cube_to_frame(p_outlist3, m_outlist3, e_outlist3, ens_dif3_list,'2030s')
df_p4 = cube_to_frame(p_outlist4, m_outlist4, e_outlist4, ens_dif4_list, '2040s')


df_p = df_p.append(df_p3)
df_p = df_p.append(df_p4)
#df_p.set_index('period', inplace=  True)

df_per = df_p.drop(columns = ['p' , 'm', 'e', 'total'])

df_abs = df_p.drop(columns = ['p_per' , 'm_per', 'e_per', 'total'])

#otherwise looks like population growth makes mortlaity go down
df_abs['p'] = df_abs['p'] * -1
df_abs['m'] = df_abs['m'] * -1
df_abs['e'] = df_abs['e'] * -1

df_abs['scen'] = pd.Categorical(df_abs['scen'], ['ssp119', 'ssp245', 'ssp585'])
df_abs = df_abs[df_abs['period'] != '2020s']

#load historical

d_abs = pd.read_csv('/nfs/see-fs-02_users/earsch/Documents/Leeds/CHAMNHA/output_csvs/d_abs_2.csv')
h_abs = pd.read_csv('/nfs/see-fs-02_users/earsch/Documents/Leeds/CHAMNHA/output_csvs/h_abs_2.csv')

#combine historical and future df for plotting

d_abs = d_abs[d_abs['coeff'] == 1.0]
d_abs['scen'] = 'hist-nat'
h_abs = h_abs[h_abs['coeff'] == 1.0]
h_abs['scen'] = 'historical'

his_df = d_abs.append(h_abs)
his_df = his_df.drop(columns = ['coeff'])


# calculate net
df_abs_net = df_abs[df_abs['period'] != '2030s']
df_abs_net['net'] = df_abs_net['p'] + df_abs_net['m'] + df_abs_net['e']

his_df_net = copy.deepcopy(his_df)
his_df_net['net'] = his_df_net['p'] + his_df_net['m'] + his_df_net['e']


#%% his + future

import matplotlib
fs = 14
matplotlib.rcParams.update({'font.size': fs})

####set up figure and axes
fig = plt.figure(figsize = (9,6))
gs = GridSpec(ncols = 4, nrows = 1, figure = fig) #3 rows, 3 cols

ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1:])

ax_list = [ax1, ax2]
for ax in ax_list:
    ax.set_ylim(-15000, 25000)

###labels and lists
#scenarios and years to cycle through
scen_his = ['hist-nat', 'historical']
scen_fut = ['ssp119', 'ssp245', 'ssp585']
years_his = ['period2']
years_fut = ['2030s', '2040s']

#labels and colours
scen_labs_his = ['Hist-nat', 'Historical']
scen_labs_fut = ['SSP119', 'SSP245', 'SSP585']    
p_labs_his = ['P2']
p_labs_fut = ['30s', '40s']

lw = 2.0

#sort scens ind ataframe so plot in correct order
df_abs['scen'] = pd.Categorical(df_abs['scen'], scen_fut)
his_df['scen'] = pd.Categorical(his_df['scen'], scen_his)


#### set up bar locations for axes 1
barWidth1 = 0.12
barWidth2 = barWidth1 / 2.0
#ax1.set_xlim(-0.1, 0.7)

# Set position of bar on X axis (position for each rcm bar within each gcm category)
n = len(scen_his)
positions = []
r1 = np.arange(n)
r1 = r1 * 0.2 # decrease space between groups
positions.append(r1)

 
#### set up bar locations for axes 2
# Set position of bar on X axis (position for each rcm bar within each gcm category)

n = len(scen_fut)
positions2 = []
r1 = np.arange(n)
r1 = r1 * 0.2 # decrease space between groups
positions2.append(r1)

for i in np.arange(1,len(years_fut)):
    next_pos = [x + barWidth2 for x in positions2[i-1]]
    positions2.append(next_pos)


## Plot ax1 - his
dat = his_df  
dat = his_df.sort_values('scen')

#label = scen_labs[i]
pos = positions[0]

ax1.bar(pos, dat['p'], width = barWidth1, color = '#1b9e77', label = 'Population',
       edgecolor = 'k')
ax1.bar(pos, dat['m'], width = barWidth1, color = '#7570b3', label = 'Mortality',
       edgecolor = 'k')
ax1.bar(pos, dat['e'], width = barWidth1, color = '#d95f02', label = 'Climate',
       bottom = dat['p'], edgecolor = 'k')



## plot ax2 - future
for i in np.arange(len(years_fut)):
    
    pos = positions2[i]
    
    #scen = scens[i]
    dat = df_abs[df_abs['period'] == years_fut[i]]   
    dat = dat.sort_values('scen')
    
    #label = scen_labs[i]
    
    ax2.bar(pos, dat['p'], width = barWidth2, color = '#1b9e77', label = 'Population',
           edgecolor = 'k')
    ax2.bar(pos, dat['m'], width = barWidth2, color = '#7570b3', label = 'Mortality',
           edgecolor = 'k')
    ax2.bar(pos, dat['e'], width = barWidth2, color = '#d95f02', label = 'Climate',
           bottom = dat['p'], edgecolor = 'k')
    
#remove y ticks
ax2.set_yticklabels([''] * 9)
    

#set x labels

x1_locs = [0.0, 0.2]
ax1.set_xticks(x1_locs)
ax1.set_xticklabels(['  2011-\n2020', '  2011-\n2020'], fontsize = 12)

x_locs = [0.0, 0.07, 
          0.2, 0.27,
          0.4, 0.47]

ax2.set_xticks(np.sort(x_locs))
xlabs = ['  2030-\n2039', '  2040-\n2049']
ax2.set_xticklabels(np.tile(xlabs,3), fontsize = 12)


#set horizontal line
ax1.axhline(0, c = 'k', linewidth = 0.5)
ax2.axhline(0, c = 'k', linewidth = 0.5)

#add scenario labels
y = -0.19
ax1.annotate(scen_labs_his[0], (0.02, y), xycoords = 'axes fraction',
                 rotation = 0, ha = 'left', fontsize = 12)
ax1.annotate(scen_labs_his[1], (0.56, y), xycoords = 'axes fraction',
                 rotation = 0, ha = 'left', fontsize = 12)
ax2.annotate(scen_labs_fut[0], (0.10, y), xycoords = 'axes fraction',
                 rotation = 0, ha = 'left', fontsize = 12)
ax2.annotate(scen_labs_fut[1], (0.45, y), xycoords = 'axes fraction',
                 rotation = 0, ha = 'left', fontsize = 12)
ax2.annotate(scen_labs_fut[2], (0.8, y), xycoords = 'axes fraction',
                 rotation = 0, ha = 'left', fontsize = 12)

#lines
width = .065
trans = ax.get_xaxis_transform()
y = -0.13
x = -0.25
ax1.plot([x, x + width], [y, y], color = 'k', transform = trans, clip_on = False)
x = -0.155
ax1.plot([x, x + width], [y, y], color = 'k', transform = trans, clip_on = False)


width = .065
trans = ax2.get_xaxis_transform()
x = 0
ax2.plot([x, x + width], [y, y], color = 'k', transform = trans, clip_on = False)
x = 0.20
ax2.plot([x, x + width], [y, y], color = 'k', transform = trans, clip_on = False)
x = 0.4
ax2.plot([x, x + width], [y, y], color = 'k', transform = trans, clip_on = False)


##
ax1.set_ylabel('Contribution to change in heat mortality')  

#
ax1.annotate('A)', (0.01, 1.02), xycoords = 'axes fraction',
             rotation = 0, ha = 'left', weight = 'bold')
ax2.annotate('B)', (0.01, 1.02), xycoords = 'axes fraction',
             rotation = 0, ha = 'left', weight = 'bold')


plt.draw()

#create legend

mor_patch = Patch(facecolor='#7570b3', edgecolor = 'k')
pop_patch = Patch(facecolor='#1b9e77', edgecolor = 'k')
e_patch = Patch(facecolor='#d95f02', edgecolor = 'k')
         

custom_lines = [pop_patch, mor_patch, e_patch]


ax.legend(custom_lines, ['Population growth', 'Mortality rate declines', 'Climate change'], 
           bbox_to_anchor=(1.1, -0.25),ncol=3,frameon = False, handletextpad = 0.5)


#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/cc_components_withhis_2.png',
#            bbox_inches = 'tight', pad_inches = 0.3)

#%% Net


import matplotlib
fs = 14
matplotlib.rcParams.update({'font.size': fs})

####set up figure and axes
fig = plt.figure(figsize = (9,6))
gs = GridSpec(ncols = 3, nrows = 1, figure = fig) #3 rows, 3 cols

ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1:])

ax_list = [ax1, ax2]
for ax in ax_list:
    ax.set_ylim(-25000, 40000)

###labels and lists
#scenarios and years to cycle through
scen_his = ['hist-nat', 'historical']
scen_fut = ['ssp119', 'ssp245', 'ssp585']
years_his = ['period2', 'net']
years_fut = ['2040s', 'net']

#labels and colours
scen_labs_his = ['Hist-nat', 'Historical']
scen_labs_fut = ['SSP119', 'SSP245', 'SSP585']    
p_labs_his = ['P2', 'net']
p_labs_fut = ['40s', 'net']

lw = 2.0

#sort scens ind ataframe so plot in correct order
df_abs['scen'] = pd.Categorical(df_abs['scen'], scen_fut)
his_df['scen'] = pd.Categorical(his_df['scen'], scen_his)


#### set up bar locations for axes 1
barWidth1 = 0.12
barWidth2 = barWidth1 
#ax1.set_xlim(-0.1, 0.7)

# Set position of bar on X axis (position for each rcm bar within each gcm category)
n = len(scen_his)
positions = []
r1 = np.arange(n)
r1 = r1 * 0.3 # decrease space between groups
positions.append(r1)

for i in np.arange(1,len(years_his)):
    next_pos = [x + barWidth1 for x in positions[i-1]]
    positions.append(next_pos)
 
#### set up bar locations for axes 2
# Set position of bar on X axis (position for each rcm bar within each gcm category)

n = len(scen_fut)
positions2 = []
r1 = np.arange(n)
r1 = r1 * 0.4 # decrease space between groups
positions2.append(r1)

for i in np.arange(1,len(years_fut)):
    next_pos = [x + barWidth2 for x in positions2[i-1]]
    positions2.append(next_pos)


## Plot ax1 - his
dat = his_df_net
dat = dat.sort_values('scen')

#label = scen_labs[i]
pos = positions[0]

ax1.bar(pos, dat['p'], width = barWidth1, color = '#1b9e77', label = 'Population',
       edgecolor = 'k')
ax1.bar(pos, dat['m'], width = barWidth1, color = '#7570b3', label = 'Mortality',
       edgecolor = 'k')
ax1.bar(pos, dat['e'], width = barWidth1, color = '#d95f02', label = 'Climate',
       bottom = dat['p'], edgecolor = 'k')
#net
pos = positions[1]
ax1.bar(pos, dat['net'], width = barWidth1, color = 'grey', label = 'Net change',
       edgecolor = 'k')


## plot ax2 - future
pos = positions2[0]
    
#scen = scens[i]
dat = df_abs_net   
dat = dat.sort_values('scen')

#label = scen_labs[i]

ax2.bar(pos, dat['p'], width = barWidth2, color = '#1b9e77', label = 'Population',
       edgecolor = 'k')
ax2.bar(pos, dat['m'], width = barWidth2, color = '#7570b3', label = 'Mortality',
       edgecolor = 'k')
ax2.bar(pos, dat['e'], width = barWidth2, color = '#d95f02', label = 'Climate',
       bottom = dat['p'], edgecolor = 'k')

#net
pos = positions2[1]
ax2.bar(pos, dat['net'], width = barWidth1, color = 'grey', label = 'Net change',
       edgecolor = 'k')
    
#remove y ticks
ax2.set_yticklabels([''] * 9)
    

#set x labels

x1_locs = [0.06, 0.36]
ax1.set_xticks(x1_locs)
ax1.set_xticklabels(['  2011-\n2020', '  2011-\n2020'], fontsize = 12)

x_locs = [0.06, 
          0.47, 
          0.87]

ax2.set_xticks(np.sort(x_locs))
xlabs = ['  2040-\n2049']
ax2.set_xticklabels(np.tile(xlabs,3), fontsize = 12)


#set horizontal line
ax1.axhline(0, c = 'k', linewidth = 0.5)
ax2.axhline(0, c = 'k', linewidth = 0.5)

#add scenario labels
y = -0.19
ax1.annotate(scen_labs_his[0], (0.11, y), xycoords = 'axes fraction',
                 rotation = 0, ha = 'left', fontsize = 12)
ax1.annotate(scen_labs_his[1], (0.58, y), xycoords = 'axes fraction',
                 rotation = 0, ha = 'left', fontsize = 12)
ax2.annotate(scen_labs_fut[0], (0.09, y), xycoords = 'axes fraction',
                 rotation = 0, ha = 'left', fontsize = 12)
ax2.annotate(scen_labs_fut[1], (0.45, y), xycoords = 'axes fraction',
                 rotation = 0, ha = 'left', fontsize = 12)
ax2.annotate(scen_labs_fut[2], (0.79, y), xycoords = 'axes fraction',
                 rotation = 0, ha = 'left', fontsize = 12)

#lines
width = .15
trans = ax1.get_xaxis_transform()
y = -0.13
x = 0.0
ax1.plot([x, x + width], [y, y], color = 'k', transform = trans, clip_on = False)

x = 0.29
ax1.plot([x, x + width], [y, y], color = 'k', transform = trans, clip_on = False)


width = .15
trans = ax2.get_xaxis_transform()
x = 0
ax2.plot([x, x + width], [y, y], color = 'k', transform = trans, clip_on = False)
x = 0.41
ax2.plot([x, x + width], [y, y], color = 'k', transform = trans, clip_on = False)
x = 0.81
ax2.plot([x, x + width], [y, y], color = 'k', transform = trans, clip_on = False)


##
ax1.set_ylabel('Contribution to change in heat mortality \n (number of deaths per year)')  

#
ax1.annotate('A)', (0.01, 1.02), xycoords = 'axes fraction',
             rotation = 0, ha = 'left', weight = 'bold')
ax2.annotate('B)', (0.01, 1.02), xycoords = 'axes fraction',
             rotation = 0, ha = 'left', weight = 'bold')


plt.draw()

#create legend

mor_patch = Patch(facecolor='#7570b3', edgecolor = 'k')
pop_patch = Patch(facecolor='#1b9e77', edgecolor = 'k')
e_patch = Patch(facecolor='#d95f02', edgecolor = 'k')        
net_patch = Patch(facecolor='grey', edgecolor = 'k')

custom_lines = [pop_patch, mor_patch, e_patch, net_patch]


ax.legend(custom_lines, ['Population growth', 'All-cause mortality rate declines', 'Climate change',
                         'Net change'], 
           bbox_to_anchor=(0.95, -0.25),ncol=2,frameon = False, handletextpad = 0.5)


#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/cc_components_withhis_4_net_coeff1.png',
#            bbox_inches = 'tight', pad_inches = 0.3)
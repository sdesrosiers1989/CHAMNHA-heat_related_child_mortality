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
    
path = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_061/historical/'

#historical
filenames = glob.glob(path + '*historical*.nc')
hist_list = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    hist_list.append(x)
        
path = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_1/historical/'


    
#future
path = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_061/future/'

#historical
filenames = glob.glob(path + '*PBC_MBC.nc')
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
        


    
 #%% actual mortality

path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/'

filenames = glob.glob(path + '*historical*')

his = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    his.append(x)
    
    
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/future/'

filenames = glob.glob(path + '*PBC_MBC.nc')
fut119_list = iris.cube.CubeList()
fut245_list = iris.cube.CubeList()
fut585_list = iris.cube.CubeList()

for file in filenames:
    x = iris.load_cube(file)
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
    cube.data.mask[np.isnan(countries.data)] = True
    cube.data = np.ma.masked_where(np.ma.getmask(cru_tas.data), cube.data)
    
futmor_list = [fut119_list, fut245_list, fut585_list]

for cube_list in futmor_list:
    for cube in cube_list:
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

#%% actual mortality years
    
years = np.unique([x.coord('year').points[0] for x in his]) 

his_mor_years = iris.cube.CubeList()

for y in years:
    
    h_list = [x for x in his if x.coord('year').points[0] == y]
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
# comp_period = 2005 - 2014
    
i_h = 1 # index corresponding to 2005 - 2014

#compare 2020, 2030, 2040

scen_list = [ssp119_years, ssp245_years, ssp585_years]

#2020
p_outlist = []
m_outlist = []
e_outlist = []

for i in np.arange(len(scen_list)):
    p_inlist = [pop_ratio_list[i_h], pop_ratio_list[3]]
    m_inlist = [dmor_2010, mor_list_fut[0]]
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
    m_inlist = [dmor_2010, mor_list_fut[1]]
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
    m_inlist = [dmor_2010, mor_list_fut[2]]
    e_inlist = [his_years[i_h], scen_list[i][2]]

    p_effect, m_effect, e_effect = calc_comp(p_inlist, m_inlist, e_inlist, return_list = False)
    
    p_outlist4.append(p_effect)
    m_outlist4.append(m_effect)
    e_outlist4.append(e_effect)


#%% Check values

#ssp119
#1995 - 2004: his
base_period = his_mor_years[i_h]
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


#%% Check decomp

print(np.nanmean(bperiod_ens.data))
print(np.nanmean(fperiod_ens.data))
print(np.nanmean(ens_dif.data))

np.nanmean(p_outlist[0].data)
np.nanmean(m_outlist[0].data)
np.nanmean(e_outlist[0].data)

np.nanmean(p_outlist[0].data) +  np.nanmean(m_outlist[0].data) +   np.nanmean(e_outlist[0].data) 

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

#load historical

d_abs = pd.read_csv('/nfs/see-fs-02_users/earsch/Documents/Leeds/CHAMNHA/output_csvs/d_abs.csv')
h_abs = pd.read_csv('/nfs/see-fs-02_users/earsch/Documents/Leeds/CHAMNHA/output_csvs/h_abs.csv')

#%% https://stackoverflow.com/questions/22787209/how-to-have-clusters-of-stacked-bars-with-python-pandas

def plot_clustered_stacked(dfall, cols, labels=None, title="",  H=['', "/", '.'], width = 0.25, z = 1, **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe"""

    n_df = len(dfall) # number of dfs
    n_col = len(dfall[0].columns) #columns
    n_ind = len(dfall[0].index) #index - groups 
    axe = plt.subplot(1,1, z)
    # x = subplot number

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
                rect.set_hatch(H[int(i / n_col)]) #edited part     
                rect.set_width(width / float(n_df + width))

    axe.set_xticks(np.arange(n_ind))
    axe.set_xticklabels(df.index, rotation = 0)

    return axe


#%% Plot absoltue mortality
    
df_119 = df_abs[df_abs['scen'] == 'ssp119']
df_245 = df_abs[df_abs['scen'] == 'ssp245']
df_585 = df_abs[df_abs['scen'] == 'ssp585']


fig = plt.figure(figsize=(12,6))

ax = plot_clustered_stacked([df_119, df_245, df_585],labels = ["ssp119", "ssp245", 'ssp585'],
                            cols = ['#1b9e77', '#7570b3', '#d95f02'])
                                   
                                    
#%%
ax_list = [ax]
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

#%% Plot mortality

import matplotlib
fs = 14
matplotlib.rcParams.update({'font.size': fs})

scens = ['ssp119', 'ssp245', 'ssp585']
years = ['2020s', '2030s', '2040s']


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(1,1,1)


barWidth = 0.1

#ax1.set_xlim(-0.1, 0.7)


# Set position of bar on X axis (position for each rcm bar within each gcm category)
n = len(years)
positions = []
r1 = np.arange(n)
#r1 = r1 * 0.5 # decrease space between groups
positions.append(r1)

for i in np.arange(1,len(scens)):
    next_pos = [x + barWidth for x in positions[i-1]]
    positions.append(next_pos)
 
#set colors
hatch_list = ['', '//', '..']
labels = ["his1", "his2", 'histnat1', 'histnat3']
scen_labs = ['SSP119', 'SSP245', 'SSP585']    
p_labs = ['2020 - 2029', '2030 - 2039', '2040 - 2049']

lw = 2.0
for i in np.arange(len(scens)):
    
    pos = positions[i]
    
    scen = scens[i]
    dat = df_abs[df_abs['scen'] == scen]   
    
    label = scen_labs[i]
    
    ax.bar(pos, dat['p'], width = barWidth, color = '#1b9e77', label = 'Population',
           edgecolor = 'k', hatch = hatch_list[i])
    ax.bar(pos, dat['m'], width = barWidth, color = '#7570b3', label = 'Mortality',
           edgecolor = 'k', hatch = hatch_list[i])
    ax.bar(pos, dat['e'], width = barWidth, color = '#d95f02', label = 'Climate',
           bottom = dat['p'], edgecolor = 'k', hatch = hatch_list[i])
    

ax.set_xticks(positions[1]) #positions 5 is the middle rcm
ax.set_xticklabels(p_labs)

ax.axhline(0, c = 'k', linewidth = 0.5)



ax.set_ylabel('Contribution to change in heat mortality')  

plt.draw()

#create legend

mor_patch = Patch(facecolor='#7570b3', edgecolor = 'k')
pop_patch = Patch(facecolor='#1b9e77', edgecolor = 'k')
e_patch = Patch(facecolor='#d95f02', edgecolor = 'k')
                
ssp119_patch = Patch(facecolor='white', edgecolor = 'k')
ssp245_patch = Patch(facecolor='white', hatch = '//', edgecolor = 'k')
ssp585_patch = Patch(facecolor='white', hatch = '..', edgecolor = 'k')


custom_lines = [pop_patch, ssp119_patch, mor_patch, ssp245_patch, e_patch, ssp585_patch]


ax.legend(custom_lines, ['Population', 'SSP119', 'Mortality', 
                         'SSP245', 'Climate', 'SSP585'], 
           bbox_to_anchor=(1.2, -0.2),ncol=3,frameon = False, handletextpad = 0.5)


#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/totalrate_mort_coeff061_range_largefont.png',
#            bbox_inches = 'tight', pad_inches = 0.3)

#%% group by scenario

import matplotlib
fs = 14
matplotlib.rcParams.update({'font.size': fs})

scens = ['ssp119', 'ssp245', 'ssp585']
years = ['2020s', '2030s', '2040s']


fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(1,1,1)


barWidth = 0.2

#ax1.set_xlim(-0.1, 0.7)


# Set position of bar on X axis (position for each rcm bar within each gcm category)
n = len(scens)
positions = []
r1 = np.arange(n)
#r1 = r1 * 0.5 # decrease space between groups
positions.append(r1)

for i in np.arange(1,len(years)):
    next_pos = [x + barWidth for x in positions[i-1]]
    positions.append(next_pos)
 
#set colors
hatch_list = ['', '//', '..']
labels = ["his1", "his2", 'histnat1', 'histnat3']
scen_labs = ['SSP119', 'SSP245', 'SSP585']    
p_labs = ['20s', '30s', '40s']

lw = 2.0
for i in np.arange(len(years)):
    
    pos = positions[i]
    
    scen = scens[i]
    dat = df_abs[df_abs['period'] == years[i]]   
    
    label = scen_labs[i]
    
    ax.bar(pos, dat['p'], width = barWidth, color = '#1b9e77', label = 'Population',
           edgecolor = 'k', hatch = hatch_list[i])
    ax.bar(pos, dat['m'], width = barWidth, color = '#7570b3', label = 'Mortality',
           edgecolor = 'k', hatch = hatch_list[i])
    ax.bar(pos, dat['e'], width = barWidth, color = '#d95f02', label = 'Climate',
           bottom = dat['p'], edgecolor = 'k', hatch = hatch_list[i])
    

x_locs = np.append(positions[0], positions[1])
x_locs = np.append(x_locs, positions[2])


ax.set_xticks(np.sort(x_locs))
ax.set_xticklabels(np.tile(p_labs,3))


#ax.set_xticks(positions[1]) #positions 5 is the middle rcm
#ax.set_xticklabels(scen_labs)

ax.axhline(0, c = 'k', linewidth = 0.5)


ax.set_ylabel('Contribution to change in heat mortality')  

plt.draw()

#create legend

mor_patch = Patch(facecolor='#7570b3', edgecolor = 'k')
pop_patch = Patch(facecolor='#1b9e77', edgecolor = 'k')
e_patch = Patch(facecolor='#d95f02', edgecolor = 'k')
                
#ssp119_patch = Patch(facecolor='white', edgecolor = 'k')
#ssp245_patch = Patch(facecolor='white', hatch = '//', edgecolor = 'k')
#ssp585_patch = Patch(facecolor='white', hatch = '..', edgecolor = 'k')


custom_lines = [pop_patch, ssp119_patch, mor_patch, ssp245_patch, e_patch, ssp585_patch]


ax.legend(custom_lines, ['Population', '2020s', 'Mortality', 
                         '2030s', 'Climate', '2040s'], 
           bbox_to_anchor=(1.2, -0.2),ncol=3,frameon = False, handletextpad = 0.5)


#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/totalrate_mort_coeff061_range_largefont.png',
#            bbox_inches = 'tight', pad_inches = 0.3)
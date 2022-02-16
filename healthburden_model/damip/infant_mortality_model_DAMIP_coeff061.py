#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Infant mortality model (due to temp)

CORDEX and CP4A

Created on Mon Apr 20 11:12:29 2020

@author: earsch
"""

#%%set wd and import packages

import iris
import iris.quickplot as qplt
import iris.coord_categorisation
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units

import numpy as np
import numpy.ma as ma

import math


import copy

import glob


#Import my functions
import sys
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/Tanga/Plot_functions')
import tanzania1 as tp

        

#%% Import CMIP6 bias corrected data
    
#### temp
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/bias_corr/CHAMNHA/tas/his/'
filenames = glob.glob(path + '*.nc')
tas = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    tas.append(x)
            
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/bias_corr/CHAMNHA/tas/end/'
tas_damip = iris.cube.CubeList()
filenames = glob.glob(path + '*hist-nat.nc')
for file in filenames:
    x = iris.load_cube(file)
    tas_damip.append(x)

tas_245 = iris.cube.CubeList()
filenames = glob.glob(path + '*ssp245.nc')
for file in filenames:
    x = iris.load_cube(file)
    tas_245.append(x)
    
#obs - used to get threshold
cru_tas = iris.load('/nfs/a321/earsch/Tanga/Data/CRU/tmp/*.nc',
                    iris.Constraint(cube_func = lambda cube: cube.var_name == 'tmp'))

cru_tas = cru_tas.concatenate_cube()
iris.coord_categorisation.add_year(cru_tas, 'time')
cru_tas = cru_tas.extract(iris.Constraint(year= lambda cell: 1995 <= cell <= 2010))
cru_tas = cru_tas.regrid(tas[0][0], iris.analysis.Linear())

#%% Limit historical and ssp245 models to same models in DAMIP

damip_mods = [tp.gcm(x) for x in tas_damip]

tas = [x for x in tas if tp.gcm(x) in damip_mods]
tas_245 = [x for x in tas_245 if tp.gcm(x) in damip_mods]


#%% get historical period and concatendate (historical + ssp245) 1995 - 2020

#for historical period  (1995 - 2005)
    #concatenate historical and ssp245

def find_mod(his, fut):
    
    his_mod = tp.gcm(his) 
    
    mod_list=  iris.cube.CubeList()
    for i in np.arange(len(fut)):
        fut_mod = tp.gcm(fut[i]) 
        
        if fut_mod == his_mod:
            mod_list.append(fut[i])
    
    mod_list.append(his)
    return mod_list

tas_concat = []

for i in np.arange(len(tas)):
    fut_cube = tas_245[i]
    mod_list = find_mod(fut_cube, tas)
    for cube in mod_list:
        cube.remove_coord('sim')
    equalise_attributes(mod_list)
    unify_time_units(mod_list)
    
    for cube in mod_list:
        cube.coord('month').attributes = ''
        cube.coord('year').attributes = ''
        cube.coord('time').attributes = ''
    
    tas_hisear = mod_list.concatenate_cube()
    tas_hisear.add_aux_coord(iris.coords.AuxCoord('historical-245', long_name = 'sim'))
    tas_concat.append(tas_hisear)
    
#%% Select years of interest 1995 - 2020

tas_damip_years = []
for cube in tas_damip:
    x = cube.extract(iris.Constraint(year = lambda cell: 1995 <= cell <= 2020))
    tas_damip_years.append(x)
    
tas_his_years = []
for cube in tas_concat:
    x = cube.extract(iris.Constraint(year = lambda cell: 1995 <= cell <= 2020))
    tas_his_years.append(x)

#%% Import input data 

#population -from world pop
# 0 - 4 year olds
pop_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000_regrid.nc')
pop_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2010_regrid.nc')
pop_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2019_regrid.nc')

#daily mortality
dmor_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2000_regrid.nc')
dmor_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2010_regrid.nc')
dmor_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2019_regrid.nc')

#%% Create coefficient data

def place_holder(base, dat):
    new_dat = copy.deepcopy(base)
    new_dat.data = np.where(~np.isnan(new_dat.data), dat, new_dat.data)
    new_dat.data = ma.masked_array(new_dat.data, mask = base.data.mask)
    return new_dat

#Coeff
per_c = 0.61
c = math.log((per_c/100) + 1)
coeff = place_holder(tas[0][0], c) 

#Threshold
#thres = cru_tas.collapsed('time', iris.analysis.PERCENTILE, percent = 75)
#thres.units = 'celsius'

#%% calculate temp diff with threshold

def calc_tdif(tavg, thres):
    
    tdif = iris.analysis.maths.subtract(tavg, thres.data)
    tdif.data = np.ma.where(tdif.data < 0, 0, tdif.data)
    return tdif

tdif_hislist = []
thres_list = []
for cube in tas_his_years:
    thres_years = cube.extract(iris.Constraint(year = lambda cell: 1995 <= cell <= 2010))
    thres = thres_years.collapsed('time', iris.analysis.PERCENTILE, percent = 75)
    tdif = calc_tdif(cube, thres)
    thres_list.append(thres)
    tdif_hislist.append(tdif)

tdif_damiplist = []
for cube in tas_damip_years:
    thres = [x for x in thres_list if tp.gcm(x) == tp.gcm(cube)]
    if len(thres) > 1:
        print('too long')
    tdif = calc_tdif(cube, thres[0])
    tdif_damiplist.append(tdif)
 


#%% calculate attributable death per decade


def ann_death_per_decade(temp, dec_start, dec_end, pop_ratio, davg_mort):
    b_dec = temp.extract(iris.Constraint(year= lambda cell: dec_start <= cell < dec_end))
    dims =  b_dec.shape
    
    #cycle through each year in decade
    #for each gridcell, extract daily data
    #from daily tdif (difference between threshold and tavg), calcualte daily att deaths
    #sum daily att deaths to get total annual deaths from heat per year
    #calculate mean annual deaths per year (for the decade of interest)
    
    year_output = iris.cube.CubeList()
    years = np.unique(b_dec.coord('year').points)
    
    for y in years:
        b_year = b_dec.extract(iris.Constraint(year= lambda cell: cell == y))
        output = copy.deepcopy(b_year)
        print(y)
        
        
        for j in np.arange(dims[1]):
            for k in np.arange(dims[2]):
                b_day = b_year[:, j, k]
                   
                if np.isnan(b_day[0].data) == False: # if not masked
                    
                    
                    p_b = pop_ratio[j,k].data #ratio of future to baseline pop for gridcell
                    c = coeff[j,k].data # coeff for gridcell
                    m_loc = davg_mort[j, k].data
                    
                    #sensitive to placement of brackets
                    daily_att_deaths = p_b * (m_loc / np.exp(c * b_day.data)) * (np.exp(c * b_day.data)- 1)
                
                    
                    output.data[:,j,k] = daily_att_deaths
        #sum total heat deaths
        out_sum = output.collapsed('time', iris.analysis.SUM)
        year_output.append(out_sum)
        
    #calculate annual heat deaths
    year_merge =  year_output.merge_cube()
    ann_avg_heat_death = year_merge.collapsed('time', iris.analysis.MEAN)
    
    return year_merge, ann_avg_heat_death



#%% Run model



def apply_model(dec_list, pop_list, mor_list, cube_list, path, path_indyears, code = ''):
     
    for i in np.arange(0, len(dec_list)):
        dstart = dec_list[i] # dec start and dec end used for subsetting climate data
        dec_end = dstart + 10  
        period = str(dstart) + str(dec_end - 1)
        
        print(period)
        
        #pop data
        pop_ratio = pop_list[i]/pop_2000 #ratio future pop to baseline pop   
        #0 in denominator causing problems
        pop_ratio.data = np.where(pop_2000.data == 0, 1, pop_ratio.data)
        pop_ratio.data = ma.masked_array(pop_ratio.data, mask = pop_2000.data.mask)
    
        #mor data
        mor = mor_list[i]
    
    
        for j in np.arange(0, len(cube_list)):
            cube = cube_list[j]
            
            print(j, tp.gcm(cube))
            
            ahd_indyear, ahd_mean = ann_death_per_decade(cube, dstart, dec_end, pop_ratio, mor)
            
            sim =  ahd_mean.coord('sim').points[0]
            #save data
                    
            save_name = sim  + '_' + tp.gcm(ahd_mean) + '_' + period + '_' + code
    
            iris.save(ahd_indyear, path_indyears + save_name + '.nc')
            iris.save(ahd_mean, path + save_name + '.nc')
    
            print(save_name, 'saved')

    
#%% HIST NAT compare to HISTNAT
#save path   
            
path0 = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/decomp_attribution2/'
path_indyears0 = '/nfs/a321/earsch/CHAMNHA/output/annual_mortality/coeff_061/thres_hismodel/decomp_attribution2/'
cube_list = tdif_damiplist

periods = ['2010', '2020']
code_list = ['e2p1m1', 'e2p2m1', 'e2p1m2',
             'e1p1m1', 'e1p2m1', 'e1p1m2',
             'e1p2m2']

pop_lists = [pop_2010, pop_2019]
dmor_lists = [dmor_2010, dmor_2019]
d_lists = [2005, 2015]

for p in np.arange(1, len(periods)):
    period = periods[p]
    path = path0 + period + '/'
    path_indyears = path_indyears0 + period + '/'

    for i in np.arange(0,2):
        code = code_list[i]
        print(i, code)
        
        if 'p2' in code:
            pop_list = [pop_lists[p]]
            #print('Using p2', np.nanmean(pop_list[0].data), np.nanmean(pop_list[1].data))      
        if 'p1' in code:
            pop_list = [pop_2000]
            #print('Using p1', np.nanmean(pop_list[0].data), np.nanmean(pop_list[1].data))
                   
        if 'm2' in code:
            mor_list = [dmor_lists[p]]
            #print('Using m2', np.nanmean(mor_list[0].data), np.nanmean(mor_list[1].data))
        if 'm1' in code:
            mor_list = [dmor_2000]
            #print('Using m1', np.nanmean(mor_list[0].data), np.nanmean(mor_list[1].data))
            
        if 'e2' in code:
            dec_start = [d_lists[p]]
            #print('Using e2', cube_list[0].coord('sim').points[0])
        if 'e1' in code:
            dec_start = [1995]
            #print('Using e1', cube_list[0].coord('sim').points[0])
             
            
        apply_model(dec_start, pop_list, mor_list, cube_list, path, path_indyears, code)



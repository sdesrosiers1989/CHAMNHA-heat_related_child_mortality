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
#scen = '585'
#scen = '245'
scen = '119'

#importing just to geg et hist-nat mdoels - already ran his for them
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/bias_corr/CHAMNHA/tas/end/'
tas = iris.cube.CubeList()
filenames = glob.glob(path + '*' + scen + '.nc')
for file in filenames:
    x = iris.load_cube(file)
    tas.append(x)

    
#obs - used to get threshold
path = '/nfs/a321/earsch/CHAMNHA/output/thres_his/'
filenames = glob.glob(path + 'historical' + '*' + '1995_2010.nc')
thres_list = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    thres_list.append(x)
    
#%% Select years of interest 1995 - 2020

tas_years = []
for cube in tas:
    x = cube.extract(iris.Constraint(year = lambda cell: 2020 <= cell <= 2100)) 
    tas_years.append(x)

#%% Import input data 

#population -from world pop
# 0 - 4 year olds
pop_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000_regrid.nc')
    #pop used for last decade of historical

pop_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2010_regrid.nc')

#daily mortality
    #mor used for last decade of historical
dmor_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2010_regrid.nc')

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



#%% calculate temp diff with threshold

def calc_tdif(tavg, thres):
    
    tdif = iris.analysis.maths.subtract(tavg, thres.data)
    tdif.data = np.ma.where(tdif.data < 0, 0, tdif.data)
    return tdif

    

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

#save path    
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/future/'
path_indyears = '/nfs/a321/earsch/CHAMNHA/output/annual_mortality/coeff_061/thres_hismodel/future/'

dec_start = np.arange(2020, 2100, 10)

for i in np.arange(2, len(dec_start)):
    dstart = dec_start[i] # dec start and dec end used for subsetting climate data
    dec_end = dstart + 10  #goes to 2015, but extracted as < not <=, so will be same time period as period 2 of damip historical mods
    period = str(dstart) + str(dec_end - 1)
    
    print(period)
    
    #pop data
    pop_ratio = pop_2010/pop_2000 #ratio future pop to baseline pop   
    #0 in denominator causing problems
    pop_ratio.data = np.where(pop_2000.data == 0, 1, pop_ratio.data)
    pop_ratio.data = ma.masked_array(pop_ratio.data, mask = pop_2000.data.mask)

    #mor data
    mor = dmor_2010


    for j in np.arange(0, len(tas_years)):
        gcm = tp.gcm(cube)
        
        thres = [x for x in thres_list if tp.gcm(x) == gcm]
        if len(thres) > 1:
            print('thres too long')
        cube = calc_tdif(tas_years[j], thres[0])
        
        print(j, gcm)
        
        ahd_indyear, ahd_mean = ann_death_per_decade(cube, dstart, dec_end, pop_ratio, mor)
        
        sim = ahd_mean.coord('sim').points[0]
        #save data
                
        save_name = sim  + '_' + tp.gcm(ahd_mean) + '_' + period + '_' + 'P2010' + '_' + 'M2010'

        iris.save(ahd_indyear, path_indyears + save_name + '.nc')
        iris.save(ahd_mean, path + save_name + '.nc')

        print(save_name, 'saved')

    

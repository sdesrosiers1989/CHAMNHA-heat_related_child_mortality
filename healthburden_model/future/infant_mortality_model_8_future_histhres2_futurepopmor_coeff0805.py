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
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/CHAMNHA/functions')
import functions as f
        

#%% Import CMIP6 bias corrected data
    
#### temp
scen = '585'
scen = '245'
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


years = np.arange(2025, 2060, 10)
pop_path = '/nfs/a321/earsch/CHAMNHA/input_data/pop/future/processed/'
mor_path = '/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/processed/'
pop_list = []
mor_list = []

for y in years:
    p_name = pop_path + 'ssp2_' + str(y) + '_04population_mf_BIASCORR2.nc'
    m_name = mor_path + 'ref_' + str(y) + '_04_totalmor_mf_BIASCORR.nc'
    
    pop_list.append(iris.load_cube(p_name))
    mor_list.append(iris.load_cube(m_name))

#%% turn into daily mortlaity

dmor_list = [x/365 for x in mor_list]



#%% Create coefficient data


#Coeff
per_c = 0.805
c = math.log((per_c/100) + 1)
coeff = f.place_holder(tas[0][0], c) 



#%% Run model

#save path    
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_0805/thres_hismodel/future/'
path_indyears = '/nfs/a321/earsch/CHAMNHA/output/annual_mortality/coeff_0805/thres_hismodel/future/'
path_e = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_0805/future/'


dec_start = np.arange(2020, 2060, 10)

for i in np.arange(3, len(dec_start)):
    dstart = dec_start[i] # dec start and dec end used for subsetting climate data
    dec_end = dstart + 10  #goes to 2015, but extracted as < not <=, so will be same time period as period 2 of damip historical mods
    period = str(dstart) + str(dec_end - 1)
    
    print(period)
    
    #pop data
    p = pop_list[i]
    pop_ratio = p/pop_2000 #ratio future pop to baseline pop   
    #0 in denominator causing problems
    pop_ratio.data = np.where(pop_2000.data == 0, 1, pop_ratio.data)
    pop_ratio.data = ma.masked_array(pop_ratio.data, mask = pop_2000.data.mask)

    #mor data
    mor = dmor_list[i]


    for j in np.arange(0, len(tas_years)):
        gcm = tp.gcm(tas_years[j])
        
        thres = [x for x in thres_list if tp.gcm(x) == gcm]
        if len(thres) > 1:
            print('thres too long')
        cube = f.calc_tdif(tas_years[j], thres[0])
        
        print(j, gcm)
        
        ahd_indyear, ahd_mean, e = f.ann_death_per_decade(cube, dstart, dec_end, pop_ratio, mor, coeff)
        
        sim = ahd_mean.coord('sim').points[0]
        #save data
                
        save_name = sim  + '_' + tp.gcm(ahd_mean) + '_' + period + '_' + 'PBC2' + '_' + 'MBC'

        iris.save(ahd_indyear, path_indyears + save_name + '.nc')
        iris.save(ahd_mean, path + save_name + '.nc')
        iris.save(e, path_e + save_name + '.nc')

        print(save_name, 'saved')

    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Infant mortality model (due to temp)

CMIP6 models

Created on Mon Apr 20 11:12:29 2020

@author: earsch
"""

#%%set wd and import packages

import iris
import iris.coord_categorisation

import numpy as np
import numpy.ma as ma

import math

import glob


#Import my functions
import sys
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/Tanga/Plot_functions')
import tanzania1 as tp
import sys
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/CHAMNHA/functions')
import functions as f

        

#%% Import CMIP6 bias corrected data
    
#### temp
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/bias_corr/CHAMNHA/tas/his/'
filenames = glob.glob(path + '*.nc')
tas = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    tas.append(x)
            
    #importing just to geg et hist-nat mdoels - already ran his for them
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/bias_corr/CHAMNHA/tas/end/'
tas_damip = iris.cube.CubeList()
filenames = glob.glob(path + '*hist-nat.nc')
for file in filenames:
    x = iris.load_cube(file)
    tas_damip.append(x)


#%% Limit historical to models NOT in DAMIP (already ran) 

damip_mods = [tp.gcm(x) for x in tas_damip]

tas = [x for x in tas if tp.gcm(x) not in damip_mods]


    
#%% Select years of interest 1995 - 2020

tas_his_years = []
for cube in tas:
    x = cube.extract(iris.Constraint(year = lambda cell: 1995 <= cell <= 2014)) # last year
    tas_his_years.append(x)

#%% Import input data 

#population -from world pop
# 0 - 4 year olds
pop_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000_regrid.nc')
pop_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2010_regrid.nc')
#pop_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2019_regrid.nc')

#daily mortality
dmor_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2000_regrid.nc')
dmor_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2010_regrid.nc')
#dmor_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2019_regrid.nc')

#%% Create coefficient data


#Coeff
#per_c = 0.61
per_c= 1.0
c = math.log((per_c/100) + 1)
coeff = f.place_holder(tas[0][0], c) 


#%% calculate temp diff with threshold


tdif_hislist = []
thres_list = []
for cube in tas_his_years:
    thres_years = cube.extract(iris.Constraint(year = lambda cell: 1995 <= cell <= 2010))
    thres = thres_years.collapsed('time', iris.analysis.PERCENTILE, percent = 75)
    tdif = f.calc_tdif(cube, thres)
    thres_list.append(thres)
    tdif_hislist.append(tdif)

#%% save thres

path = '/nfs/a321/earsch/CHAMNHA/output/thres_his/'
for cube in thres_list:
    save_name = 'historical'  + '_' + tp.gcm(cube) + '_' + '1995_2010'
    iris.save(cube, path + save_name + '.nc')    



#%% Run model

#save path    
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_1/thres_hismodel/'
path_indyears = '/nfs/a321/earsch/CHAMNHA/output/annual_mortality/coeff_1/thres_hismodel/'
path_e = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_1/historical/'


dec_start = [1995, 2005]
pop_list = [pop_2000, pop_2010]
mor_list = [dmor_2000, dmor_2010]

for i in np.arange(0, len(dec_start)):
    dstart = dec_start[i] # dec start and dec end used for subsetting climate data
    dec_end = dstart + 10  #goes to 2015, but extracted as < not <=, so will be same time period as period 2 of damip historical mods
    period = str(dstart) + str(dec_end - 1)
    
    print(period)
    
    #pop data
    pop_ratio = pop_list[i]/pop_2000 #ratio future pop to baseline pop   
    #0 in denominator causing problems
    pop_ratio.data = np.where(pop_2000.data == 0, 1, pop_ratio.data)
    pop_ratio.data = ma.masked_array(pop_ratio.data, mask = pop_2000.data.mask)

    #mor data
    mor = mor_list[i]


    for j in np.arange(0, len(tdif_hislist)):
        cube = tdif_hislist[j]
        
        print(j, tp.gcm(cube))
        
        ahd_indyear, ahd_mean, e = f.ann_death_per_decade(cube, dstart, dec_end, pop_ratio, mor, coeff)
        
        sim =  ahd_mean.coord('sim').points[0]
        #save data
                
        save_name = sim  + '_' + tp.gcm(ahd_mean) + '_' + period

        iris.save(ahd_indyear, path_indyears + save_name + '.nc')
        iris.save(ahd_mean, path + save_name + '.nc')
        iris.save(e, path_e + save_name + '.nc')

        print(save_name, 'saved')

    

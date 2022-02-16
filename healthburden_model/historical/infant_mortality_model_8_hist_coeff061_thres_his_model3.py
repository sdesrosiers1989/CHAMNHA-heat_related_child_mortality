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
    if 'FGOALS' in file:
        x = x.extract(iris.Constraint(year= lambda cell: cell <= 2014))
    
    tas.append(x)
            
    #importing just to geg et hist-nat mdoels - already ran his for them
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/bias_corr/CHAMNHA/tas/end/'
tas_damip = iris.cube.CubeList()
filenames = glob.glob(path + '*hist-nat.nc')
for file in filenames:
    x = iris.load_cube(file)
    tas_damip.append(x)


    
#obs - used to get threshold
cru_tas = iris.load('/nfs/a321/earsch/Tanga/Data/CRU/tmp/*.nc',
                    iris.Constraint(cube_func = lambda cube: cube.var_name == 'tmp'))

cru_tas = cru_tas.concatenate_cube()
iris.coord_categorisation.add_year(cru_tas, 'time')
cru_tas = cru_tas.extract(iris.Constraint(year= lambda cell: 1995 <= cell <= 2010))
cru_tas = cru_tas.regrid(tas[0][0], iris.analysis.Linear())

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
#dmor_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2000_regrid.nc')


countries = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/africa_countryname.nc')
popfrac_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/popfrac_2019_regrid.nc')
countries_regrid = countries.regrid(popfrac_2019, iris.analysis.Nearest())
countries_regrid.coord('latitude').guess_bounds()
countries_regrid.coord('longitude').guess_bounds()


popfrac_2019.coord('latitude').guess_bounds()
popfrac_2019.coord('longitude').guess_bounds()

p2000_regrid = pop_2000.regrid(popfrac_2019, iris.analysis.AreaWeighted())
p2010_regrid = pop_2010.regrid(popfrac_2019, iris.analysis.AreaWeighted())

p2010_regrid.data = np.ma.masked_array(p2010_regrid.data, countries_regrid.data.mask)
p2000_regrid.data = np.ma.masked_array(p2000_regrid.data, countries_regrid.data.mask)

m2000_regrid = dmor_2000.regrid(popfrac_2019, iris.analysis.AreaWeighted())
m2010_regrid = dmor_2010.regrid(popfrac_2019, iris.analysis.AreaWeighted())

m2010_regrid.data = np.ma.masked_array(m2010_regrid.data, countries_regrid.data.mask)
m2000_regrid.data = np.ma.masked_array(m2000_regrid.data, countries_regrid.data.mask)




p2000_regrid = p2000_regrid.regrid(tas[0][0], iris.analysis.AreaWeighted())
p2010_regrid = p2010_regrid.regrid(tas[0][0], iris.analysis.AreaWeighted())
m2000_regrid = m2000_regrid.regrid(tas[0][0], iris.analysis.AreaWeighted())
m2010_regrid = m2010_regrid.regrid(tas[0][0], iris.analysis.AreaWeighted())
#%% Create coefficient data


#Coeff
per_c = 0.61
#per_c= 1.0
c = math.log((per_c/100) + 1)
coeff = f.place_holder(tas[0][0], c) 

#Threshold
thres = cru_tas.collapsed('time', iris.analysis.PERCENTILE, percent = 75)
thres.units = 'celsius'

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

#path = '/nfs/a321/earsch/CHAMNHA/output/thres_his/'
#for cube in thres_list:
#    save_name = 'historical'  + '_' + tp.gcm(cube) + '_' + '1995_2010'
#    iris.save(cube, path + save_name + '.nc')    



#%% Run model

#save path    
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/'
path_indyears = '/nfs/a321/earsch/CHAMNHA/output/annual_mortality/coeff_061/thres_hismodel/'
path_e = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_061/historical/'


dec_start = [1995, 2005]
pop_list = [p2000_regrid, p2010_regrid]
mor_list = [m2000_regrid, m2010_regrid]

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

        iris.save(ahd_indyear, path_indyears + save_name + '_2.nc')
        iris.save(ahd_mean, path + save_name + '_2.nc')
        #iris.save(e, path_e + save_name + '_2.nc')

        print(save_name, 'saved')

    

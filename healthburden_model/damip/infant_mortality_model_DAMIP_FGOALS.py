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
import iris.coord_categorisation
from iris.util import equalise_attributes
from iris.util import unify_time_units

import numpy as np
import numpy.ma as ma

import math

import glob


#Import my functions
import sys
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/Tanga/Plot_functions')
import tanzania1 as tp
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/CHAMNHA/functions')
import functions as f
        

#%% Import CMIP6 bias corrected data
    
#### temp
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/bias_corr/CHAMNHA/tas/his/'
filenames = glob.glob(path + '*FGOALS*.nc')
tas = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    x = x.extract(iris.Constraint(year= lambda cell: cell <= 2014))
    tas.append(x)
            
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/bias_corr/CHAMNHA/tas/end/'
tas_damip = iris.cube.CubeList()
filenames = glob.glob(path + '*FGOALS**hist-nat.nc')
for file in filenames:
    x = iris.load_cube(file)
    tas_damip.append(x)

tas_245 = iris.cube.CubeList()
filenames = glob.glob(path + '*FGOALS*ssp245.nc')
for file in filenames:
    x = iris.load_cube(file)
    x = x.extract(iris.Constraint(year= lambda cell: cell >= 2015))
    tas_245.append(x)
    

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



#Coeff
per_c = 0.61
per_c = 0.805
#per_c = 1.0
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

tdif_damiplist = []
for cube in tas_damip_years:
    thres = [x for x in thres_list if tp.gcm(x) == tp.gcm(cube)]
    if len(thres) > 1:
        print('too long')
    tdif = f.calc_tdif(cube, thres[0])
    tdif_damiplist.append(tdif)
 
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


dec_start = [1995, 2005, 2015]
pop_list = [pop_2000, pop_2010, pop_2019]
mor_list = [dmor_2000, dmor_2010, dmor_2019]

for i in np.arange(2, len(dec_start)):
    dstart = dec_start[i] # dec start and dec end used for subsetting climate data
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


    for j in np.arange(0, len(tdif_hislist)):
        cube = tdif_hislist[j]
        
        print(j, tp.gcm(cube))
        
        ahd_indyear, ahd_mean, e  = f.ann_death_per_decade(cube, dstart, dec_end, pop_ratio, mor, coeff)
        
        sim =  ahd_mean.coord('sim').points[0]
        #save data
                
        save_name = sim  + '_' + tp.gcm(ahd_mean) + '_' + period

        iris.save(ahd_indyear, path_indyears + save_name + '.nc')
        iris.save(ahd_mean, path + save_name + '.nc')
        iris.save(e, path_e + save_name + '.nc')

        print(save_name, 'saved')

    

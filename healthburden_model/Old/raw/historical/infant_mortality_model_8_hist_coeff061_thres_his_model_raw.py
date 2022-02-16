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

def add_att_from_filename(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[52:-3]
    
    gcm = file.partition('_')[0]
    sim = file.partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

#%% Import CMIP6 bias corrected data
    
#### temp
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/Processed/tas/'
filenames = glob.glob(path + '*.nc')

tas_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'tas')


tas = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file,  tas_constraint, callback = add_att_from_filename)
    tas.append(x)

tas_damip = [x for x in tas if x.coord('sim').points[0] == 'hist-nat']
tas_his = [x for x in tas if x.coord('sim').points[0] == 'historical']

    
#obs - used to get threshold
cru_tas = iris.load('/nfs/a321/earsch/Tanga/Data/CRU/tmp/*.nc',
                    iris.Constraint(cube_func = lambda cube: cube.var_name == 'tmp'))

cru_tas = cru_tas.concatenate_cube()
iris.coord_categorisation.add_year(cru_tas, 'time')
cru_tas = cru_tas.extract(iris.Constraint(year= lambda cell: 1995 <= cell <= 2010))
cru_tas = cru_tas.regrid(tas[0][0], iris.analysis.Linear())

#%%
damip_mods = [tp.gcm(x) for x in tas_damip]

tas_his = [x for x in tas_his if tp.gcm(x) in damip_mods]
    
#%% Select years of interest 1995 - 2020

tas_his_years = []
for cube in tas_his:
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

#%% Create coefficient data

def place_holder(base, dat):
    new_dat = copy.deepcopy(base)
    new_dat.data = np.where(~np.isnan(new_dat.data), dat, new_dat.data)
    new_dat.data = ma.masked_array(new_dat.data, mask = base.data.mask)
    return new_dat

#Coeff
per_c = 0.61
#per_c= 1.0
c = math.log((per_c/100) + 1)
coeff = place_holder(tas[0][0], c) 



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

#%% save thres

path = '/nfs/a321/earsch/CHAMNHA/output/thres_his/raw/'
for cube in thres_list:
    save_name = 'historical'  + '_' + tp.gcm(cube) + '_' + '1995_2010'
    iris.save(cube, path + save_name + '.nc')    

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
    
    e_list = iris.cube.CubeList()
    
    for y in years:
        b_year = b_dec.extract(iris.Constraint(year= lambda cell: cell == y))
        output = copy.deepcopy(b_year)
        e_output = copy.deepcopy(b_year)
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
                    e = (1 / np.exp(c * b_day.data)) * (np.exp(c * b_day.data)- 1)
                    
                    output.data[:,j,k] = daily_att_deaths
                    e_output.data[:,j,k] = e
        #sum total heat deaths
        out_sum = output.collapsed('time', iris.analysis.SUM)
        year_output.append(out_sum)
        e_mean = e_output.collapsed('time', iris.analysis.MEAN)
        e_list.append(e_mean)
        
    #calculate annual heat deaths
    year_merge =  year_output.merge_cube()
    ann_avg_heat_death = year_merge.collapsed('time', iris.analysis.MEAN)
    e_merge = e_list.merge_cube()
    
    return year_merge, ann_avg_heat_death, e_merge



#%% Run model

#save path    
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/raw/'
path_indyears = '/nfs/a321/earsch/CHAMNHA/output/annual_mortality/coeff_061/thres_hismodel/raw/'
path_e = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_061/raw/historical/'


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
        
        ahd_indyear, ahd_mean, e = ann_death_per_decade(cube, dstart, dec_end, pop_ratio, mor)
        
        sim =  ahd_mean.coord('sim').points[0]
        #save data
                
        save_name = sim  + '_' + tp.gcm(ahd_mean) + '_' + period

        iris.save(ahd_indyear, path_indyears + save_name + '.nc')
        iris.save(ahd_mean, path + save_name + '.nc')
        iris.save(e, path_e + save_name + '.nc')

        print(save_name, 'saved')

    

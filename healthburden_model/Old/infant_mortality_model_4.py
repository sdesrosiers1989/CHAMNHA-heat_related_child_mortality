#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Infant mortality model (due to temp)

CORDEX and CP4A

Created on Mon Apr 20 11:12:29 2020

@author: earsch
"""

#%%set wd and import packages
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec

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

#Import my functions
import sys
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/Tanga/Plot_functions')
import tanzania1 as tp

def add_att_from_filename_tas(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[67:-3]
    gcm = file.partition('_')[0]
    mod = file.partition('_')[2].partition('_')[0]
    sim = file.partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

def add_att_from_filename_tashis(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[70:-3]
    gcm = file.partition('_')[0]
    mod = file.partition('_')[2].partition('_')[0]
    sim = file.partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

def add_att_from_filename_tas45(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[59:-3]
    gcm = file.partition('_')[0]
    mod = file.partition('_')[2].partition('_')[0]
    sim = file.partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))
    
def add_att_from_filename_tasmax(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[68:-3]
    gcm = file.partition('_')[0]
    mod = file.partition('_')[2].partition('_')[0]
    sim = file.partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))


def add_att_from_filename_tasmax45(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[65:-3]
    gcm = file.partition('_')[0]
    mod = file.partition('_')[2].partition('_')[0]
    sim = file.partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

def add_att_from_filename_pr(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[63:-3]
    gcm = file.partition('_')[0]
    mod = file.partition('_')[2].partition('_')[0]
    sim = file.partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

def add_att_from_filename_prhist(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[64:-3]
    gcm = file.partition('_')[0]
    mod = file.partition('_')[2].partition('_')[0]
    sim = file.partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

def add_att_from_filename_pr45(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[54:-3]
    gcm = file.partition('_')[0]
    mod = file.partition('_')[2].partition('_')[0]
    sim = file.partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def add_time(cube_list):
    for cube in cube_list:
        iris.coord_categorisation.add_year(cube, 'time')
        iris.coord_categorisation.add_month(cube, 'time')
        iris.coord_categorisation.add_day_of_month(cube, 'time')
        
def extract_box(cube_list, box_con):
    output = iris.cube.CubeList()
    for cube in cube_list:
        x = cube.extract(box_con)
        output.append(x)
    return output

#%% Import CORDEX data
    
#### temp
path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/CHAMNHA/tas/his/'
filenames = glob.glob(path + '*.nc')
tas = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    tas.append(x)
    
path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/CHAMNHA/tas/ear/'
filenames = glob.glob(path + '*.nc')
for file in filenames:
    x = iris.load_cube(file)
    tas.append(x)
        
path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/CHAMNHA/tas/mid/'
tas_fut = iris.cube.CubeList()
filenames = glob.glob(path + '*.nc')
for file in filenames:
    x = iris.load_cube(file)
    tas_fut.append(x)
        
path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/CHAMNHA/tas/end/'
filenames = glob.glob(path + '*.nc')
for file in filenames:
    x = iris.load_cube(file)
    tas_fut.append(x)

cube_orog = iris.load('/nfs/a277/IMPALA/data/4km/ANCILS/orog_combined_POSTac144_ancil_4km.nc')
orog = cube_orog[9] 

ulon = orog[0,0,:,:].coord('grid_longitude')
transform = ulon.coord_system.as_cartopy_projection()

#obs
cru_tas = iris.load('/nfs/a321/earsch/Tanga/Data/CRU/tmp/*.nc',
                    iris.Constraint(cube_func = lambda cube: cube.var_name == 'tmp'))

cru_tas = cru_tas.concatenate_cube()
iris.coord_categorisation.add_year(cru_tas, 'time')
cru_tas = cru_tas.extract(iris.Constraint(year= lambda cell: 1971 <= cell <= 2000))
cru_tas = cru_tas.regrid(tas[0][0], iris.analysis.Linear())

#%% get historical period and concatendate (historical + historicalear) 1995 - 2005

#for historical period  (1995 - 2005)
    #concatenate historical and historical ear

def find_mod(his, fut):
    
    his_mod = tp.gcm(his) + tp.model(his)
    
    mod_list=  iris.cube.CubeList()
    for i in np.arange(len(fut)):
        fut_mod = tp.gcm(fut[i]) + tp.model(fut[i])
        
        if fut_mod == his_mod:
            mod_list.append(fut[i])
    
    mod_list.append(his)
    return mod_list

tas_his = [x for x in tas if x.coord('sim').points == 'historical']
tas_ear = [x for x in tas if x.coord('sim').points == 'historicalear']

tas_concat = []

for i in np.arange(len(tas_ear)):
    ear_cube = tas_ear[i]
    mod_list = find_mod(ear_cube, tas_his)
    for cube in mod_list:
        cube.remove_coord('sim')
    tas_hisear = mod_list.concatenate_cube()
    tas_hisear.add_aux_coord(iris.coords.AuxCoord('historical', long_name = 'sim'))
    tas_concat.append(tas_hisear)

#%% Import input data

#population -from world pop
# 0 - 4 year olds
pop_f = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/global_f_0_2000_1km.nc')
pop_m = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/global_m_0_2000_1km.nc')
pop_2000 = pop_f + pop_m


#%% Placeholder data

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
thres = place_holder(tas[0][0], 30.68)
#thres = cru_tas.collapsed('time', iris.analysis.PERCENTILE, percent = 75)
thres.units = 'celsius'

#baseline pop
bpop = place_holder(tas[0][0], 2519000)

#pop mid 2000
mpop = place_holder(tas[0][0], 2551915)

#import daily avg mortality palceholder
davg_mort = np.loadtxt(open('/nfs/a321/earsch/CHAMNHA/daily_avg_mortality_placeholder.csv', 'rb'),
                       delimiter = ',', skiprows = 1)


#%% regrid input data to tas

cs = tas[0][0].coord_system(iris.coord_systems.CoordSystem)
pop_2000.coord('longitude').coord_system = cs
pop_2000.coord('latitude').coord_system = cs

pop_2000 = pop_2000.regrid(tas[0][0], iris.analysis.Linear())
pop_2000 = pop_2000 #original unit is people per pixel (1km2) - 
                    #as its a a ratio only needs to be same unit as future

#%% calculate temp diff with threshold

def calc_tdif(tavg, thres):
    
    tdif = iris.analysis.maths.subtract(tavg, thres.data)
    tdif.data = np.where(tdif.data < 0, 0, tdif.data)
    tdif.data = ma.masked_array(tdif.data, mask = tavg.data.mask)
    return tdif


tdif_hislist = []
for cube in tas_concat[0:1]:
    tdif = calc_tdif(cube, thres)
    tdif_hislist.append(tdif)

tdif_futlist = []
for cube in tas_fut[0:1]:
    tdif = calc_tdif(cube, thres)
    tdif_futlist.append(tdif)
#%% calculate attributable death per decade



def ann_death_per_decade(base, dec_start, dec_end, pop_ratio):
    b_dec = base.extract(iris.Constraint(year= lambda cell: dec_start <= cell < dec_end))
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
        
        for j in np.arange(dims[1]):
            for k in np.arange(dims[2]):
                b_day = b_year[:, j, k]
    
                if np.isnan(b_day[0].data) == False: # if not masked
                    
                    ndays = b_day.shape[0] # number of days in temperature data
                                             #used to get rigth number of days from daily avg mortality info
                    p_b = pop_ratio[j,k].data #ratio of future to baseline pop for gridcell
                    c = coeff[j,k].data # coeff for gridcell
                    
                    #sensitive to placement of brackets
                    daily_att_deaths = p_b * (davg_mort[0:ndays] / np.exp(c * b_day.data)) * (np.exp(c * b_day.data)- 1)
                
                    
                    output.data[:,j,k] = daily_att_deaths
        #sum total heat deaths
        out_sum = output.collapsed('time', iris.analysis.SUM)
        year_output.append(out_sum)
        
    #calculate annual heat deaths
    year_merge =  year_output.merge_cube()
    ann_avg_heat_death = year_merge.collapsed('time', iris.analysis.MEAN)
    
    return ann_avg_heat_death



#%%

#saev path    
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_30/'
end =  'his'
#historical period (1995 - 2005)
pop_ratio = pop_2000 / pop_2000 #ratio future pop to baseline pop

dec_start = 1995
dec_end = 2005    
period = str(dec_start) + str(dec_end)


for cube in tdif_hislist:
    ahd = ann_death_per_decade(cube, dec_start, dec_end, pop_ratio)
    
    sim =  ahd.coord('sim').points[0]
    #save data
    spath = path + end + '/'
    
    save_name = tp.gcm(ahd) + '_' + tp.model(ahd) + '_' + sim  + '_' + period
    save_path = spath + save_name + '.nc'
    iris.save(ahd, save_path)
    print(save_name, 'saved')

    
#%%
#furue period
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_30/'
end =  'fut'

#historical period (1995 - 2005)
pop_ratio = pop_2000 / pop_2000 #ratio future pop to baseline pop

dec_start = 2025
dec_end = 2035    
period = str(dec_start) + str(dec_end)

for cube in tdif_futlist:
    ahd = ann_death_per_decade(cube, dec_start, dec_end, pop_ratio)
    
    sim =  ahd.coord('sim').points[0]
    #save data
    spath = path + end + '/'
    
    save_name = tp.gcm(ahd) + '_' + tp.model(ahd) + '_' + sim + '_' + period
    
    save_path = spath + save_name + '.nc'
    iris.save(ahd, save_path)
    print(save_name, 'saved')
    

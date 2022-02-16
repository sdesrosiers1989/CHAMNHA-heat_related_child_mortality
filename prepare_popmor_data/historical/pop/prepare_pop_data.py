#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process WorldPop data for use

-already converted to netcdf
-downloaded 0 (0mth - 12mth) and 1(12mth - 4 years inclusive) m and f
-Combine for all years and subset to africa only

Created on Tue Aug 17 11:44:49 2021

@author: earsch
"""
#%% import libraris

import iris
import iris.quickplot as qplt

import glob

import numpy as np

#%%

path = '/nfs/a321/earsch/CHAMNHA/input_data/pop/'

filenames = glob.glob(path + '*2000*.nc')
clist_2000 = iris.cube.CubeList()
for file in filenames:
    clist_2000.append(iris.load_cube(file))

filenames = glob.glob(path + '*2010*.nc')
clist_2010 = iris.cube.CubeList()
for file in filenames:
    clist_2010.append(iris.load_cube(file))
    
filenames = glob.glob(path + '*2019*.nc')
clist_2019 = iris.cube.CubeList()
for file in filenames:
    clist_2019.append(iris.load_cube(file))
    
#%% Subset to africa
    
min_lat = -40.0
max_lat = 40.0
max_lon = 60.0
min_lon = -25.0

def safrica_lat(input):
    return min_lat  <= input <= max_lat 

def safrica_long(input):
    return min_lon  <= input <= max_lon 

afr_con = iris.Constraint(latitude = safrica_lat, longitude = safrica_long)

africa_list = []
list_of_lists = [clist_2000, clist_2010, clist_2019]
for clist in list_of_lists:
    new_list = iris.cube.CubeList()
    for cube in clist:
        x = cube.extract(afr_con)
        #x.units = ('1')
        new_list.append(x)
    africa_list.append(new_list)
    
#%% Combine data together for each year
    
    
comb_2000 = africa_list[0][0]
for i in np.arange(1, len(africa_list[0])):
    comb_2000 = iris.analysis.maths.add(comb_2000, africa_list[0][i])

comb_2010 = africa_list[1][0]
for i in np.arange(1, len(africa_list[1])):
    comb_2010 = iris.analysis.maths.add(comb_2010, africa_list[1][i])
    
comb_2019 = africa_list[2][0]
for i in np.arange(1, len(africa_list[2])):
    comb_2019 = iris.analysis.maths.add(comb_2019, africa_list[2][i])
    
#%% save output
    
iris.save(comb_2000, '/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000.nc')
iris.save(comb_2010, '/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2010.nc')
iris.save(comb_2019, '/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2019.nc')



#%% Import saved pop and mort data, regrid to tas data and mulitple by size of gridcells

p2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000.nc')
p2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2010.nc')
p2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2019.nc')

mor_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2000.nc')
mor_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2010.nc')
mor_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2019.nc')

dat_list = [p2000, p2010, p2019, mor_2000, mor_2010, mor_2019]

#import model data for regridding
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/bias_corr/CHAMNHA/tas/his/'
filenames = glob.glob(path + '*.nc')
tas = iris.load_cube(filenames[0])



#%%
cs = tas[0][0].coord_system(iris.coord_systems.CoordSystem)

regrid_list = []
for cube in dat_list:

    cube.coord('longitude').coord_system = cs
    cube.coord('latitude').coord_system = cs

    new_cube = cube.regrid(tas[0], iris.analysis.Linear())
    #original unit is people per pixel (1km2) - 
    #as its a a ratio only needs to be same unit as future

    regrid_list.append(new_cube)

#%%
iris.save(regrid_list[0], '/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000_regrid.nc')
iris.save(regrid_list[1], '/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2010_regrid.nc')
iris.save(regrid_list[2], '/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2019_regrid.nc')

iris.save(regrid_list[3], '/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2000_regrid.nc')
iris.save(regrid_list[4], '/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2000_regrid.nc')
iris.save(regrid_list[5], '/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2000_regrid.nc')

x = regrid_list[3] / 365 # daily_mor
iris.save(x, '/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2000_regrid.nc')
x = regrid_list[4] / 365 # daily_mor
iris.save(x, '/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2000_regrid.nc')
x = regrid_list[5] / 365 # daily_mor
iris.save(x, '/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2000_regrid.nc')


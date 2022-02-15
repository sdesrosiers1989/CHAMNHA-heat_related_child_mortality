#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Data Preparation  - CMIP6 tas
Date: 28/01/2019
Author: Sarah C
Purpose: Concantenate separeate model files
1. Import model data 
    CMIP6 tas - DAMIP hist-nat, historical, SCENARIOMIP ssp119 ssp585
2. For each model, concatenate years into one cube
3. Regrid onto coarsest 
4. Save data
            
'''

#%% Import packages

import matplotlib.pyplot as plt 

import iris
import iris.plot as iplt
import iris.quickplot as qplt
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units
import iris.coord_categorisation
from cf_units import Unit

import cartopy.crs as ccrs

import glob

import numpy as np

import itertools

import sys
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/Tanga/Plot_functions')
import tanzania1 as tp

proj = ccrs.PlateCarree(central_longitude = 38)

#%% Import model data  on native grids
# https://github.com/ESMValGroup/ESMValCore/issues/576

#CMIP6 data (tas = near surface air temperature, which is between 1.5 to 10m)
tas_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'tas')
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/ScenarioMIP/'
filenames = glob.glob(path+ 'ssp585/tas*day*.nc')
filenames_ssp119 = glob.glob(path+ 'ssp119/tas*day*.nc')
filenames_ssp245 = glob.glob(path+ 'ssp245/tas*day*.nc')

filenames.extend(filenames_ssp119)
filenames.extend(filenames_ssp245)

filenames = [x for x in filenames if 'CESM' in x]


filenames = sorted(filenames) # need to sort before splitting otherwise will end up with duplicate groups

#split filenames into groups based on gcm_rcm combo - will need to load by individual sim as years
#will be loaded separately, and will need to concatenate only the same simulation together

def find_gcmsim(filename):
    file = filename[61:]
    gcm = file.partition('_')[0]
    sim = file.partition('_')[2].partition('_')[0]
    gcm_sim = gcm + sim 
    return gcm_sim

new_files = [list(k) for g, k in itertools.groupby(filenames,key = lambda x:find_gcmsim(x))]

#load file by groups
cmip6_cubes = iris.cube.CubeList()
for file_list in new_files:
    cubes = iris.load(file_list, tas_constraint)
    cmip6_cubes.append(cubes)

#%%
#Get coord system to assign to all cordex

mpi_hr_hist = iris.load('/nfs/a321/earsch/Tanga/Data/CMIP6/ScenarioMIP/ssp585/tas_day_MPI-ESM1-2-HR_ssp585_*.nc', tas_constraint)
cs = mpi_hr_hist[0].coord_system(iris.coord_systems.CoordSystem)

def add_time_dim(cube):
    tcoord0_auxtime = cube.coord('time')
    tcoord0_dimtime = iris.coords.DimCoord(tcoord0_auxtime.points, 
                                               standard_name = tcoord0_auxtime.standard_name, 
                                               units = tcoord0_auxtime.units, 
                                               var_name = tcoord0_auxtime.var_name, 
                                               attributes = tcoord0_auxtime.attributes)
    cube.remove_coord('time')
    cube.add_dim_coord(tcoord0_dimtime, 0)
    return cube

def add_lat_lon_dim(cube, base):
    latcoord = base.coord('latitude')
    loncoord = base.coord('longitude')
    
    cube.remove_coord('latitude')
    cube.remove_coord('longitude')
    cube.add_dim_coord(latcoord, 1)
    cube.add_dim_coord(loncoord, 2)
    
    latitude = cube.coord('latitude')
    latitude.bounds = latitude.bounds.astype(np.float64)
    latitude.bounds=np.round(latitude.bounds,4)
    longitude = cube.coord('longitude')
    longitude.bounds = longitude.bounds.astype(np.float64)
    
    return cube 

def md(cc, cc2): #cc = cube coord
    print(cc.standard_name, cc2.standard_name)
    print(cc.long_name, cc2.long_name)
    print(cc.var_name, cc2.var_name)
    print(cc.units, cc2.units)
    print(cc.attributes, cc2.attributes)
    print(cc.shape, cc2.shape)
    print(cc.points.dtype, cc2.points.dtype)
    print(cc.bounds_dtype, cc2.bounds_dtype)
    
#md(item[0].coord('latitude'), item[1].coord('latitude'))
#md(item[0].coord('longitude'), item[1].coord('longitude'))


    #%%

#concatenate and ensure all have same coord system
cmip6_cat = iris.cube.CubeList()
for i in np.arange(len(cmip6_cubes)):
    item =  cmip6_cubes[i]
    equalise_attributes(item)
    unify_time_units(item)
    for j in item:
        try:
            j.remove_coord('height')
        except:
            print('no height')
    if i == 0:
        item[0] = add_time_dim(item[0])
        
        x = item[0].data
        time_points = item[0].coord('time').points
        for k in np.arange(1, len(item)):
            x = np.concatenate((x, item[k].data), axis = 0)
            time_points = np.append(time_points, item[k].coord('time').points)
        
        x = np.dstack((item[0].data, item[1].data))
        
        for i in item:
            base = item[1]
            i = i.regrid(base, iris.analysis.AreaWeighted())
            
            
            latitude = i.coord('latitude')
            i.coord('latitude').bounds = latitude.bounds.astype(np.float64)
            i.coord('latitude').bounds=np.round(latitude.bounds,4)
            longitude = i.coord('longitude')
            i.coord('longitude').bounds = longitude.bounds.astype(np.float64)
            i.coord('longitude').bounds=np.round(longitude.bounds,4)
            
    x = item.concatenate_cube()
    x.coord('longitude').coord_system = cs
    x.coord('latitude').coord_system = cs
    cmip6_cat.append(x)

    
# get data of interest only, 1971 - 1999, and 2071 - 2099
    #drop 2099 from HadGEM2 CCLM4 rcp8.5 as not full year
goodyears = []
for cube in cmip6_cat:
    try:
        iris.coord_categorisation.add_year(cube, 'time')
    except:
        print('Has year.')
    mod = cube.attributes['source_id']
    sim = cube.attributes['experiment_id']
    if sim == 'historical':
        x = cube.extract(iris.Constraint(year = lambda cell: 1960 <= cell ))
    else:
        x = cube

    goodyears.append(x)
           
for cube in goodyears:
    #print (cube.coord('time').units)
    atts = cube.attributes
    print(atts['experiment_id'], atts['source_id'], cube.shape)
    

#%% Regrid first so can do area-weighted, then extract

min_lat = -40.0
max_lat = 40.0

p25_sub = iris.Constraint(latitude = lambda cell: min_lat < cell < max_lat)
 
base = mpi_hr_hist[0][0]
#base.coord('longitude').guess_bounds()
#base.coord('latitude').guess_bounds()
base = base.extract(p25_sub)

base_subset = base.intersection(longitude=(-25, 60)) 
#because ciruclar coords have a split in Africa - can't extract using constraint


for cube in goodyears:
    try:
        cube.coord('longitude').guess_bounds()
        cube.coord('latitude').guess_bounds()
    except:
        print('Has bounds.')
        pass

regridded = iris.cube.CubeList()

for cube in goodyears:
    x = cube.regrid(base_subset, iris.analysis.AreaWeighted())
    regridded.append(x)



#%% save cmip6 datadata

for cube in regridded:
    try:
        atts = cube.attributes
        save_name = 'tas' + '_' + atts['source_id'] + '_' + atts['experiment_id']
        save_path = '/nfs/a321/earsch/Tanga/Data/CMIP6/Processed/tas/' + save_name + '.nc'
        iris.save(cube, save_path)
    except:
        print('Will remove some metadata when saving...')
        pass

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BIAS CORRECTION
Purpose: Bias correct rainfall and temperature using linear-scaling
Steps:
    1. Import Data
    2. Calculate corretion factors
    3. Apply correction factors

Created on Fri April  26 12:08 2019
@author: earsch
"""

# set wd and import packages
import iris
import iris.quickplot as qplt
import iris.coord_categorisation

import numpy as np

import cartopy.crs as ccrs

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
    

def add_att_from_filename_masks(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[49:-3]
    mod_type = file.partition('_')[0]
    if mod_type == 'corgrid':
        #gcm = 'um'
        mod = 'p25'
    elif mod_type == 'trmmsize':
        #gcm = 'none'
        mod = 'trmm'
    elif mod_type == 'cp4size':
        #gcm = 'um'
        mod = 'cp4'
    else:
        #gcm = file.partition('_')[0]
        mod = file.partition('_')[2]
    
    #cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))

def find_mod(his, fut):
    
    his_mod = tp.gcm(his) 
    
    mod_list=  iris.cube.CubeList()
    for i in np.arange(len(fut)):
        fut_mod = tp.gcm(fut[i]) 
        
        if fut_mod == his_mod:
            mod_list.append(fut[i])
    
    mod_list.append(his)
        
    return mod_list

#%% Define bias corr functions

def regrid_obs(obs, ls, base):
    #regrid model and mask to observations
    
    cs = base.coord_system(iris.coord_systems.CoordSystem)
    obs.coord('latitude').guess_bounds()
    obs.coord('longitude').guess_bounds()
    
    obs.coord('longitude').coord_system = cs
    obs.coord('latitude').coord_system = cs
        #x = cube.extract(iris.Constraint(grid_latitude = lambda cell: cell <= 12.0))
    obs_regrid = obs.regrid(base, iris.analysis.AreaWeighted())
        
    ls.coord('longitude').coord_system = cs
    ls.coord('latitude').coord_system = cs
    ls_regrid = ls.regrid(base, iris.analysis.Linear())
    
    return obs_regrid, ls_regrid

def extract_years(obs, mod_list, obs_years, myears ):
    #extract to time periods of interest
    # dd year and month coords
    # convert to celsius
    
    for cube in mod_list:
        cube.convert_units('celsius')
            
    # add year and month to obs
    try:
        iris.coord_categorisation.add_year(obs, 'time')
    except:
        print('Has year.')
                   
    try:
        iris.coord_categorisation.add_month(obs, 'time', name = 'month')
    except:
        print('Has month.')
        
    # add year and month to mods            
    for cube in mod_list:
        try:
            iris.coord_categorisation.add_year(cube, 'time')
        except:
            print('Has year.')
            

    for cube in mod_list:
        try:
            iris.coord_categorisation.add_month(cube, 'time')
        except:
            print('Has month.')
    
    #extract years of interst from cru
    cru_years = obs.extract(iris.Constraint(year = lambda cell: obs_years[0] <= cell <= obs_years[1]))
    
    tas_future = iris.cube.CubeList()
    tas_calc_corr = iris.cube.CubeList()
    tas_his = iris.cube.CubeList()
    for cube in mod_list:
        sim = cube.coord('sim').points[0]
        if sim == 'historical':
            x = cube.extract(iris.Constraint(year = lambda cell: obs_years[0] <= cell <= obs_years[1]))
            #x = masking(cube, ls_cor, small_mask = True)
            tas_calc_corr.append(x)
            x = cube.extract(iris.Constraint(year = lambda cell: myears[0] <= cell ))
            #x = masking(cube, ls_cor, small_mask = True)
            tas_his.append(x)
        else:
            
            x = cube.extract(iris.Constraint(year = lambda cell: myears[0] <= cell <= myears[1]))
            #x = masking(cube, ls_cor, small_mask = True)
            tas_future.append(x)
    
    return cru_years, tas_calc_corr, tas_his, tas_future

def get_avg_month(cru_years, tas_calc_corr):
    
    #print('getting cru month')
    cru_month = cru_years.aggregated_by(['month'], iris.analysis.MEAN)
    
    #print('getting mod months', len(tas_his))
    tas_months = iris.cube.CubeList();
    for cube in tas_calc_corr:
        x = cube.aggregated_by(['month'], iris.analysis.MEAN)
        tas_months.append(x)
    
    return cru_month, tas_months

def calc_corr_fact(obs_month, mod_month):
 
    p_cor = copy.deepcopy(mod_month)
    obs = copy.deepcopy(obs_month)
    
    #If historical model data is zero will lead to correction factor being NA and introducting
    # NAS. If replace future corrected data with 0 where 0 in historical means can't get more rainfall
    # in very dry areas.
    p_cor.data = obs.data - p_cor.data
    
    return p_cor

def apply_corr(p_cor, mod_list, path, end):
    #apply correction factor
    
    month_dict = {'Jan': 0, 'Feb': 1, 'Mar':2, 'Apr':3, 'May':4, 'Jun': 5, 'Jul':6, 
               'Aug':7, 'Sep':8, 'Oct': 9, 'Nov': 10, 'Dec': 11}
    
    for cube in mod_list:
        sim = cube.coord('sim').points[0]

        dims = cube.shape
        corrected_cube = copy.deepcopy(cube)
        #cycle through days and apply monthly correction factor to each day
        for t in np.arange(0, dims[0]):
            month = cube[t].coord('month').points[0] #find month
            month_num = month_dict[month] #get number associated with month
            corr_month = p_cor[month_num] #get correction factor associated with month
            if corr_month.coord('month').points[0] == cube[t].coord('month').points[0]:
                new_day = copy.deepcopy(cube[t])
                new_day.data = new_day.data + corr_month.data
                corrected_cube.data[t] = new_day.data
            else:
                print('Month mismatch')
        
        #save data
        spath = path + end + '/'
        
        save_name = 'tas' + '_' + tp.gcm(cube) + '_' + sim
        save_path = spath + save_name + '.nc'
        iris.save(corrected_cube, save_path)
        print('Saving ', save_name)
        

def bias_corr(obs_regrid, mod_list, ls_regrid, obs_years, myears, path):
    
    print('Extracting years...')
    #extract years and time period of interest
    #convert to celsius
    cru_years, tas_calc_corr, tas_his, tas_future = extract_years(obs_regrid, mod_list, obs_years, myears)

    print('Averaging months...')
    #get avg months (for creation of correction factors)
    cru_month, tas_months = get_avg_month(cru_years, tas_calc_corr)
    
    print('Calculating correction factors..')
    p_cor = calc_corr_fact(cru_month, tas_months[0])
    
    print('Applying correction factors..')
    apply_corr(p_cor, tas_his, path, end = 'his')
    apply_corr(p_cor, tas_future, path, end = 'end')

#%%
    
    ''' Temperature '''
    
    
'''
Tmodel(day) = Tmodel(day) + MEANmonth(Pobs) - MEANmonth(Pmodel)

For each model and each month calculate a correction factor        
'''

## define variables

# years of interest
obs_years = [1970, 1989] # use these years for bias-correction (create correction factor)
myears = [1970, 2100] #start year for historical data, end year future data

#%%
tas_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'tas')

##load obs data
# Get CRU
cru_tas = iris.load('/nfs/a321/earsch/Tanga/Data/CRU/tmp/*.nc',
                    iris.Constraint(cube_func = lambda cube: cube.var_name == 'tmp'))

cru_tas = cru_tas.concatenate_cube()


#masks
#mask
ls = iris.load_cube('/nfs/a277/IMPALA/data/4km/ANCILS/landseamask_ancil_4km_regrid.nc')
ls = ls[0,0,:,:] #remove first two time, surface coords that are unused

ls.coord('longitude').points = ls.coord('longitude').points - 360

#tas future
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/Processed/tas/'
filenames = glob.glob(path + '*.nc')
#filenames = [ x for x in filenames if 'CESM' not in x]
tas = iris.cube.CubeList()

for file in filenames:
    x = iris.load_cube(file,  tas_constraint, callback = add_att_from_filename)
    tas.append(x)

tas_his = [x for x in tas if x.coord('sim').points[0] == 'historical']
tas = [x for x in tas if x.coord('sim').points[0] == 'ssp585']



    
#%% Run bias corr function for each historical model seperately

#save path
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/bias_corr/CHAMNHA/tas/'

#regrid obs to cmip6
obs_regrid, regrid_ls = regrid_obs(cru_tas, ls, tas[0][0])

for i in np.arange(22, 24):
    his_cube = tas_his[i]
    mod_list = find_mod(his_cube, tas)
    print(i, tp.gcm(his_cube))
    
    bias_corr(obs_regrid, mod_list, regrid_ls, obs_years, myears, path)
    


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

import matplotlib.pyplot as plt 

import iris
import iris.plot as iplt
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
sys.path.append('/nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/Tanga/Onset_functions')
from onset_functions import masking

proj = ccrs.PlateCarree(central_longitude = 38)

'''Set up plotting variables '''
col = 'GnBu'
col_dif = 'RdBu'
col_dif_o = 'RdBu_r' # for onset, where makes sense ti have reversed colorbar
cmap = 'viridis' #cmap for plotting onsets/cessations/durations

levs = np.arange(0,400,60) # Rainfall (monthly)
levs_dif = np.arange(-195,220,30) #Monthly rainfall difference
levs_onsets = np.arange(0, 390, 30)
levs_dif_rcp = np.arange(-90,110,20)
levs_dif_rcp_per = np.arange(-90,110,20)
levs_dif_rcp_o = np.arange(-45,55,10) #rcp differences, onset
#fourier
levs_dif_f = [-1.5, -0.5, 0.5, 1.5]
levs_f = [-0.5, 0.5, 1.5]


xticks = [10,25,40] # if need to change ticks for plotting Safrica

#add model and sim data to files
def add_att_from_filename_um(cube, field, filename):
    file = filename[62:]
    mod = filename[62:65]
    #split filename into sections by '_', find second partition and split again
    sim = file.partition('_')[2].partition('_')[2].partition('_')[0]
    
    cube.add_aux_coord(iris.coords.AuxCoord('um', long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

def add_att_from_filename_rain(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[45:-3]
    mod_type = file.partition('_')[0]
    if mod_type == 'trmm':
        gcm = 'none'
        mod = mod_type
        sim = 'historical'
    elif mod_type == 'chirps':
        gcm = 'none'
        mod = mod_type
        sim = 'historical'
    elif mod_type == 'pr':
        gcm = file.partition('_')[2].partition('_')[0]
        mod = file.partition('_')[2].partition('_')[2].partition('_')[0]
        sim = file.partition('_')[2].partition('_')[2].partition('_')[2]
    else:
        gcm = 'um'
        mod = mod_type
        sim = file.partition('_')[2].partition('_')[2].partition('_')[2]
        if sim[0] == 'h':
            sim = 'historical'
        if sim[0] == 'r':
            sim = 'rcp85'
    
    #print(gcm, mod, sim)

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

def add_att_from_filename_tas(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[67:-3]
    mod_type = file.partition('_')[0]
    if (mod_type == 'cp4') or (mod_type == 'p25'):
        gcm = 'um'
        mod = mod_type
        sim = file.partition('_')[2].partition('_')[2]
        if sim == 'hist':
            sim = 'historical'
    else:
        gcm = mod_type
        mod = file.partition('_')[2].partition('_')[0]
        sim = file.partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))
    
def add_att_from_filename_tas45(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[65:-3]
    mod_type = file.partition('_')[0]
    if (mod_type == 'cp4') or (mod_type == 'p25'):
        gcm = 'um'
        mod = mod_type
        sim = file.partition('_')[2].partition('_')[2]
        if sim == 'hist':
            sim = 'historical'
    else:
        gcm = mod_type
        mod = file.partition('_')[2].partition('_')[0]
        sim = file.partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

def add_att_from_filename_tasmax(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[60:-3]
    mod = file.partition('_')[0]
    sim = file.partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord('um', long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

def add_att_from_filename_cor(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[62:-13]
    gcm = file.partition('_')[0]
    mod = file.partition('_')[2].partition('_')[0]
    sim = file.partition('_')[2].partition('_')[2]
    
    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
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

def add_att_from_filename_tashis(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[53:-3]
    mod_type = file.partition('_')[0]
    if (mod_type == 'cp4') or (mod_type == 'p25'):
        gcm = 'um'
        mod = mod_type
        sim = file.partition('_')[2].partition('_')[2]
    else:
        gcm = mod_type
        mod = file.partition('_')[2].partition('_')[0]
        sim = file.partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
    cube.add_aux_coord(iris.coords.AuxCoord(mod, long_name = 'model'))
    cube.add_aux_coord(iris.coords.AuxCoord(sim, long_name = 'sim'))

def find_mod(his, fut):
    
    his_mod = tp.gcm(his) + tp.model(his)
    
    mod_list=  iris.cube.CubeList()
    for i in np.arange(len(fut)):
        fut_mod = tp.gcm(fut[i]) + tp.model(fut[i])
        
        if fut_mod == his_mod:
            mod_list.append(fut[i])
    
    mod_list.append(his)
        
    return mod_list

#%% Define bias corr functions

def regrid(obs, mod_list, area):
    #regrid obs to model
    
    def safrica_lat(input):
        return area[0] <= input <= area[1]

    def safrica_long(input):
        return area[2] <= input <= area[3]

    saf_con2 = iris.Constraint(latitude = safrica_lat, longitude = safrica_long)
    
    safrica_mods = iris.cube.CubeList()
    for cube in mod_list:
        try:
            cube.coord('latitude').guess_bounds()
            cube.coord('longitude').guess_bounds()
        except:
            print('Has bounds')
        x = cube.extract(saf_con2)
        safrica_mods.append(x)
        
    #regrid obs
    cs = safrica_mods[0].coord_system(iris.coord_systems.CoordSystem)
    
    obs.coord('longitude').coord_system = cs
    obs.coord('latitude').coord_system = cs  
    obs.coord('latitude').guess_bounds()
    obs.coord('longitude').guess_bounds()
    safrica_cru = obs.regrid(safrica_mods[0], iris.analysis.AreaWeighted())
    
    return safrica_cru, safrica_mods

def add_year_coords(obs, mod_list, obs_years):
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
    
    return cru_years, mod_list

def get_avg_month(cru_years, tas_his):
    
    cru_month = cru_years.aggregated_by(['month'], iris.analysis.MEAN)
    
    tas_months = iris.cube.CubeList();
    for cube in tas_his:
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

def apply_corr(p_cor, mod_list, path):
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
        
        save_name = 'tas' + '_' + tp.gcm(cube) + '_' + tp.model(cube) + '_' + sim
        save_path = path + save_name + '.nc'
        iris.save(corrected_cube, save_path)
        print('Saving ', save_name)
        


def bias_corr(obs, mod_list, area,  obs_years, path):
    print('Regridding...')
    #regrid mod, mask to model
    safrica_cru, safrica_mods = regrid(obs, mod_list, area)
    
    print('Add year coords...')
    #add year, month coords
    #convert to celsius
    cru_years, mod_list = add_year_coords(safrica_cru, safrica_mods, obs_years)

    print('Averaging months...')
    #get avg months (for creation of correction factors)
    tas_his = [x for x in mod_list if x.coord('sim').points[0] == 'historical']
    tas_future = [x for x in mod_list if x.coord('sim').points[0] == 'rcp85']
    cru_month, tas_months = get_avg_month(cru_years, tas_his)
    
    print('Calculating correction factors..')
    p_cor = calc_corr_fact(cru_month, tas_months[0])
    
    print('Applying correction factors..')
    apply_corr(p_cor, tas_his, path)
    apply_corr(p_cor, tas_future, path)

#%%
    
    ''' Temperature '''
    
    
'''
Tmodel(day) = Tmodel(day) + MEANmonth(Pobs) - MEANmonth(Pmodel)

For each model and each month calculate a correction factor        
'''

## define variables
# define area
min_lat = -35.0
max_lat = 17.0
min_lon = -20.0
max_lon = 53.0
area = [min_lat, max_lat, min_lon, max_lon]

obs_years = [1997, 2006]

#%%

##load obs data
# Get CRU
cru_tmin = iris.load('/nfs/a321/earsch/Tanga/Data/CRU/tmn/*.nc',
                    iris.Constraint(cube_func = lambda cube: cube.var_name == 'tmn'))
cru_tmax = iris.load('/nfs/a321/earsch/Tanga/Data/CRU/tmx/*.nc',
                    iris.Constraint(cube_func = lambda cube: cube.var_name == 'tmx'))

cru_tmin = cru_tmin.concatenate_cube()
cru_tmax = cru_tmax.concatenate_cube()


#tasmax
path = '/nfs/a321/earsch/Tanga/Data/CP4_Processed/tasmax/'
filenames = glob.glob(path + '*.nc')
tasmax = iris.cube.CubeList()
tasmax_his = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file,  callback = add_att_from_filename_tasmax)
    if 'historical' in file:
        tasmax_his.append(x)
    else:
        tasmax.append(x)

  
#tasmin
path = '/nfs/a321/earsch/Tanga/Data/CP4_Processed/tasmin/'
filenames = glob.glob(path + '*.nc')
tasmin = iris.cube.CubeList()
tasmin_his = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file,   callback = add_att_from_filename_tasmax)
    if 'historical' in file:
        tasmin_his.append(x)
    else:
        tasmin.append(x)
    

   
    
#%% Run bias corr function for each historical model seperately

var = 'tasmax'
#var = 'tasmin'
if var == 'tasmax':
    tas_his = tasmax_his
    tas = tasmax
    cru_tas = cru_tmax
    #save path
    path = '/nfs/a321/earsch/Tanga/Data/CP4_Processed/tasmax/bias_corr/'
elif var == 'tasmin':
    tas_his = tasmin_his
    tas = tasmin
    cru_tas = cru_tmax
    #save path
    path = '/nfs/a321/earsch/Tanga/Data/CP4_Processed/tasmin/bias_corr/'
    

for i in np.arange(0, len(tas_his)):
    print(i)
    his_cube = tas_his[i]
    mod_list = find_mod(his_cube, tas)
    
    #cru doesn't have ocean, so bias corr will automatically remove oceans
    bias_corr(cru_tas, mod_list, area,  obs_years, path)
    


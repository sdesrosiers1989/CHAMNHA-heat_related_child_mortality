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
    
def add_att_from_filename_tashis_1971(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[70:-3]
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
    cube.add_aux_coord(iris.coords.AuxCoord('historicalear', long_name = 'sim')) #otherwise will extract his years
    
def add_att_from_filename_tas45(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[59:-3]
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
    file = filename[59:-3]
    mod_type = file.partition('_')[0]
    if (mod_type == 'cp4') or (mod_type == 'p25'):
        gcm = 'um'
        mod = mod_type
        sim = file.partition('_')[2]
    else:
        gcm = mod_type
        mod = file.partition('_')[2].partition('_')[0]
        sim = file.partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'gcm'))
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

def regrid(obs, mod_list, area, ls):
    #regrid model and mask to observations
    
    def safrica_lat(input):
        return area[0] <= input <= area[1]

    def safrica_long(input):
        return area[2] <= input <= area[3]

    saf_con2 = iris.Constraint(latitude = safrica_lat, longitude = safrica_long)
    
    safrica_cru = obs.extract(saf_con2)
    
    cs = safrica_cru.coord_system(iris.coord_systems.CoordSystem)
    safrica_cru.coord('latitude').guess_bounds()
    safrica_cru.coord('longitude').guess_bounds()
    
    safrica_mods = iris.cube.CubeList()
    for cube in mod_list:
        try:
            cube.coord('grid_latitude').guess_bounds()
            cube.coord('grid_longitude').guess_bounds()
        except:
            print('Has bounds')
        cube.coord('grid_longitude').coord_system = cs
        cube.coord('grid_latitude').coord_system = cs
        #x = cube.extract(iris.Constraint(grid_latitude = lambda cell: cell <= 12.0))
        x = cube.regrid(safrica_cru, iris.analysis.AreaWeighted())
        safrica_mods.append(x)
        
    ls.coord('longitude').coord_system = cs
    ls.coord('latitude').coord_system = cs
    ls_cor = ls.regrid(safrica_mods[0], iris.analysis.Linear())
    
    return safrica_cru, safrica_mods, ls_cor

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
    tas_mid = iris.cube.CubeList()
    tas_his = iris.cube.CubeList()
    for cube in mod_list:
        sim = cube.coord('sim').points[0]
        if sim == 'historical':
            x = cube.extract(iris.Constraint(year = lambda cell: obs_years[0] <= cell <= obs_years[1]))
            #x = masking(cube, ls_cor, small_mask = True)
            tas_his.append(x)
        else:
            x = cube.extract(iris.Constraint(year = lambda cell: myears[0] <= cell <= myears[1]))
            #x = masking(cube, ls_cor, small_mask = True)
            tas_mid.append(x)
            x = cube.extract(iris.Constraint(year = lambda cell: myears[2] <= cell <= myears[3]))
            #x = masking(cube, ls_cor, small_mask = True)
            tas_future.append(x)
    
    return cru_years, tas_his, tas_mid, tas_future

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
        
        save_name = 'tas' + '_' + tp.gcm(cube) + '_' + tp.model(cube) + '_' + sim
        save_path = spath + save_name + '.nc'
        iris.save(corrected_cube, save_path)
        print('Saving ', save_name)
        


def bias_corr(obs, mod_list, area, ls, obs_years, myears, path):
    print('Regridding...')
    #regrid mod, mask to model
    safrica_cru, safrica_mods, ls_cor = regrid(obs, mod_list, area, ls)
    
    print('Extracting years...')
    #extract years and time period of interest
    #convert to celsius
    cru_years, tas_his, tas_mid, tas_future = extract_years(safrica_cru, safrica_mods, obs_years, myears)

    print('Averaging months...')
    #get avg months (for creation of correction factors)
    cru_month, tas_months = get_avg_month(cru_years, tas_his)
    
    print('Calculating correction factors..')
    p_cor = calc_corr_fact(cru_month, tas_months[0])
    
    print('Applying correction factors..')
    #apply_corr(p_cor, tas_his, path, end = 'his')
    apply_corr(p_cor, tas_mid, path, end = 'ear')
    #apply_corr(p_cor, tas_future, path, end = 'end')

#%%
    
    ''' Temperature '''
    
    
'''
Tmodel(day) = Tmodel(day) + MEANmonth(Pobs) - MEANmonth(Pmodel)

For each model and each month calculate a correction factor        
'''

## define variables
# define area
min_lat = -35.0
max_lat = 30.0
min_lon = -20.0
max_lon = 53.0
area = [min_lat, max_lat, min_lon, max_lon]

# years of interest
obs_years = [1971, 2000]
#obs_yers = [1997, 2006]
myears = [2001, 2005, 2061, 2062] #start end year mid, start end year end-century

#%% !!!Just to get 2001 - 2005 (actually from historical sim)

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


#tas future (1971 - 2006)
path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Processed/tas_hist_1971onwards/'
filenames = glob.glob(path + '*historical.nc')
filenames= [x for x in filenames if 'p25' not in x]
tas = iris.cube.CubeList()
tas_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'tas')
tas2_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'd03236')
for file in filenames:
    print(file)
    if ('p25' in file) or ('cp4' in file):
        print('Skip um')
    else:
        x = iris.load_cube(file,  tas_constraint, callback = add_att_from_filename_tashis_1971)
    tas.append(x)
    


#tas his
path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Processed/tas/'
filenames = glob.glob(path + '*historical.nc')
filenames= [x for x in filenames if 'p25' not in x]
tas_his = iris.cube.CubeList()
tas_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'tas')
tas2_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'd03236')
for file in filenames:
    print(file)
    if ('p25' in file) or ('cp4' in file):
        print('Skip um')
    else:
        x = iris.load_cube(file,  tas_constraint, callback = add_att_from_filename_tashis)
    tas_his.append(x)
    
    
#%% Run bias corr function for each historical model seperately to get 2000 - 2005 only
#historical 2001 - 2005 data put in mid position

#save path
path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/CHAMNHA/tas/'

for i in np.arange(0, len(tas_his)):
    print(i)
    his_cube = tas_his[i]
    mod_list = find_mod(his_cube, tas)
    
    bias_corr(cru_tas, mod_list, area, ls, obs_years, myears, path)
   
    
    
    
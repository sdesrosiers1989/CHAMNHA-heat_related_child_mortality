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


#%%
    
    ''' Temperature '''
    
    
'''
Tmodel(day) = Tmodel(day) + MEANmonth(Pobs) - MEANmonth(Pmodel)

For each model and each month calculate a correction factor        
'''

#load temperature data

#tas
path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Processed/tas_rcp85allyears/'
filenames = glob.glob(path + '*.nc')
tas = iris.cube.CubeList()
tas_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'tas')
tas2_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'd03236')
for file in filenames:
    if ('p25' in file) or ('cp4' in file):
        print('skipping um') # already done
    else:
        x = iris.load_cube(file,  tas_constraint, callback = add_att_from_filename_tas)
    tas.append(x)
    
#add rcp45
#path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Processed/tas_rcp45/'
#filenames = glob.glob(path + '*.nc')
#for file in filenames:
#    x = iris.load_cube(file,  tas_constraint, callback = add_att_from_filename_tas45)
#    tas.append(x)

#tasmax
path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Processed/tasmax/'
filenames = glob.glob(path + '*.nc')
tasmax = iris.cube.CubeList()
tasmax_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'tasmax')
tasmax2_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'd03236')
for file in filenames:
    if ('p25' in file) or ('cp4' in file):
        x = iris.load_cube(file,  tasmax2_constraint, callback = add_att_from_filename_tasmax)
    else:
        x = iris.load_cube(file,  tasmax_constraint, callback = add_att_from_filename_tasmax)
    tasmax.append(x)
    
#tasmin
path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Processed/tasmin/'
filenames = glob.glob(path + '*.nc')
tasmin = iris.cube.CubeList()
tasmin_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'tasmin')
tasmin2_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'd03236')
for file in filenames:
    if ('p25' in file) or ('cp4' in file):
        x = iris.load_cube(file,  tasmin2_constraint, callback = add_att_from_filename_tasmax)
    else:
        x = iris.load_cube(file,  tasmin_constraint, callback = add_att_from_filename_tasmax)
    tasmin.append(x)
    
tas_mod = [tas, tasmax, tasmin]
tas_mod = [tas]

#masks
#mask
ls = iris.load_cube('/nfs/a277/IMPALA/data/4km/ANCILS/landseamask_ancil_4km_regrid.nc')
ls = ls[0,0,:,:] #remove first two time, surface coords that are unused

ls.coord('longitude').points = ls.coord('longitude').points - 360

# Get CRU
cru_tas = iris.load('/nfs/a321/earsch/Tanga/Data/CRU/tmp/*.nc',
                    iris.Constraint(cube_func = lambda cube: cube.var_name == 'tmp'))
cru_tmin = iris.load('/nfs/a321/earsch/Tanga/Data/CRU/tmn/*.nc',
                    iris.Constraint(cube_func = lambda cube: cube.var_name == 'tmn'))
cru_tmax = iris.load('/nfs/a321/earsch/Tanga/Data/CRU/tmx/*.nc',
                    iris.Constraint(cube_func = lambda cube: cube.var_name == 'tmx'))

cru_tas = cru_tas.concatenate_cube()
cru_tmin = cru_tmin.concatenate_cube()
cru_tmax = cru_tmax.concatenate_cube()

cru = [cru_tas]##, cru_tmin, cru_tmax]

#%%
#Get sub-Saharan Africa
min_lat = -35.0
max_lat = 17.0
min_lon = -20.0
max_lon = 53.0

def safrica_lat(input):
    return min_lat <= input <= max_lat

def safrica_long(input):
    return min_lon <= input <= max_lon 

saf_con2 = iris.Constraint(latitude = safrica_lat, longitude = safrica_long)

safrica_cru = iris.cube.CubeList()
for cube in cru:
    x = cube.extract(saf_con2)
    safrica_cru.append(x)

cs = safrica_cru[0].coord_system(iris.coord_systems.CoordSystem)
safrica_cru[0].coord('latitude').guess_bounds()
safrica_cru[0].coord('longitude').guess_bounds()

safrica_tas = iris.cube.CubeList()
for cube_list in tas_mod:
    new_list = iris.cube.CubeList()
    for cube in cube_list:
        try:
            cube.coord('latitude').guess_bounds()
            cube.coord('longitude').guess_bounds()
        except:
            print('Has bounds')
        cube.coord('grid_longitude').coord_system = cs
        cube.coord('grid_latitude').coord_system = cs
        #x = cube.extract(iris.Constraint(grid_latitude = lambda cell: cell <= 12.0))
        x = cube.regrid(safrica_cru[0], iris.analysis.AreaWeighted())
        new_list.append(x)
    safrica_tas.append(new_list)

#regrid ls
ls.coord('longitude').coord_system = cs
ls.coord('latitude').coord_system = cs
ls_cor = ls.regrid(safrica_tas[0][0], iris.analysis.Linear())
    
#After regridding, make smaller

for cube_list in safrica_tas:
    for cube in cube_list:
        cube.convert_units('celsius')
            

for cube in safrica_cru:
    try:
        iris.coord_categorisation.add_year(cube, 'time')
    except:
        print('Has year.')
               
for cube in safrica_cru:
    try:
        iris.coord_categorisation.add_month(cube, 'time', name = 'month')
    except:
        print('Has month.')
        
for cube_list in safrica_tas:
    for cube in cube_list:
        try:
            iris.coord_categorisation.add_year(cube, 'time')
        except:
            print('Has year.')
            
for cube_list in safrica_tas:
    for cube in cube_list:
        try:
            iris.coord_categorisation.add_month(cube, 'time')
        except:
            print('Has month.')

cru_um = iris.cube.CubeList()
cru_cor = iris.cube.CubeList()
for cube in safrica_cru:
    um = cube.extract(iris.Constraint(year = lambda cell: 1997 <= cell <= 2006))
    cor = cube.extract(iris.Constraint(year = lambda cell: 1971 <= cell <= 2000))
    cru_um.append(um)
    cru_cor.append(cor)
cru_years = [cru_um, cru_cor]

tas_future = iris.cube.CubeList()
tas_mid = iris.cube.CubeList()
tas_his = iris.cube.CubeList()
for cube_list in safrica_tas:
    new_list_his = iris.cube.CubeList()
    new_list_mid = iris.cube.CubeList()
    new_list_future = iris.cube.CubeList()
    for cube in cube_list:
        sim = cube.coord('sim').points[0]
        if sim == 'historical':
            x = cube.extract(iris.Constraint(year = lambda cell: 1971 <= cell <= 2000))
            #x = masking(cube, ls_cor, small_mask = True)
            new_list_his.append(x)
        else:
            x = cube.extract(iris.Constraint(year = lambda cell: 2031 <= cell <= 2060))
            #x = masking(cube, ls_cor, small_mask = True)
            new_list_mid.append(x)
            x = cube.extract(iris.Constraint(year = lambda cell: 2071 <= cell <= 2100))
            #x = masking(cube, ls_cor, small_mask = True)
            new_list_future.append(x)
    tas_his.append(new_list_his)
    tas_mid.append(new_list_mid)
    tas_future.append(new_list_future)

cru_month = iris.cube.CubeList()
for cube_list in cru_years:
    new_list = iris.cube.CubeList()
    for cube in cube_list:
        x = cube.aggregated_by(['month'], iris.analysis.MEAN)
        new_list.append(x)
    cru_month.append(new_list) 

tas_months = iris.cube.CubeList();
for cube_list in tas_his:
    new_list = iris.cube.CubeList()
    for cube in cube_list:
        x = cube.aggregated_by(['month'], iris.analysis.MEAN)
        new_list.append(x)
    tas_months.append(new_list)

#%%

p_cor = iris.cube.CubeList()

for cube in tas_months[0]:
    if tp.gcm(cube) == 'um':
        obs_orig = cru_month[0][0]
        if np.min(obs_orig.coord('year').bounds) != 1997:
            print('Wrong cru file.')
        print(tp.model(cube), 'using 1997 chirps')
    else:
        obs_orig = cru_month[1][0]
        if np.min(obs_orig.coord('year').bounds) != 1971:
            print('Wrong cru file.')
    
    model = copy.deepcopy(cube)
    obs = copy.deepcopy(obs_orig)
    
    #If historical model data is zero will lead to correction factor being NA and introducting
    # NAS. If replace future corrected data with 0 where 0 in historical means can't get more rainfall
    # in very dry areas.
    model.data = obs.data - model.data
    p_cor.append(model)
    
#apply correction factor
    
month_dict = {'Jan': 0, 'Feb': 1, 'Mar':2, 'Apr':3, 'May':4, 'Jun': 5, 'Jul':6, 
           'Aug':7, 'Sep':8, 'Oct': 9, 'Nov': 10, 'Dec': 11}

tas_his_corrected = iris.cube.CubeList()

for cube in tas_his[0]:
    #Get correction cube
    mod = tp.gcm(cube) + tp.model(cube)
    for corr in p_cor:
        cor_mod = tp.gcm(corr) + tp.model(corr)
        if cor_mod == mod:
            cor_cube = corr
            print('Correction cube found:', mod, cor_mod)
            break
    #Apply correction factor to daily data
    if cor_mod == mod:
        dims = cube.shape
        corrected_cube = copy.deepcopy(cube)
        #cycle through days and apply monthly correction factor to each day
        for t in np.arange(0, dims[0]):
            month = cube[t].coord('month').points[0] #find month
            month_num = month_dict[month] #get number associated with month
            corr_month = cor_cube[month_num] #get correction factor associated with month
            if corr_month.coord('month').points[0] == cube[t].coord('month').points[0]:
                new_day = copy.deepcopy(cube[t])
                new_day.data = new_day.data + corr_month.data
                corrected_cube.data[t] = new_day.data
            else:
                print('Month mismatch')
        tas_his_corrected.append(corrected_cube)
    else:
        print('Model mismatch:', cor_mod, mod)
       
def apply_cor(tas_future_list, path):
    tas_fut_corrected = iris.cube.CubeList()
    for cube in tas_future[0]:
        #Get correction cube
        mod = tp.gcm(cube) + tp.model(cube)
        for corr in p_cor:
            cor_mod = tp.gcm(corr) + tp.model(corr)
            if cor_mod == mod:
                cor_cube = corr
                print('Correction cube found:', mod, cor_mod)
                break
        if cor_mod == mod:
            #Apply correction factor to daily data
            dims = cube.shape
            corrected_cube = copy.deepcopy(cube)
            #cycle through days and apply monthly correction factor to each day
            for t in np.arange(0, dims[0]):
                month = cube[t].coord('month').points[0] #find month
                month_num = month_dict[month] #get number associated with month
                corr_month = cor_cube[month_num] #get correction factor associated with month
                if corr_month.coord('month').points[0] == cube[t].coord('month').points[0]:
                    new_day = copy.deepcopy(cube[t])
                    new_day.data = new_day.data + corr_month.data
                    corrected_cube.data[t] = new_day.data
                else:
                    print('Month mismatch')
                    
            #save data

            save_name = 'tas' + '_' + tp.gcm(cube) + '_' + tp.model(cube) + '_' + 'rcp85'
            save_path = path + save_name + '.nc'
            iris.save(cube, save_path)
            print('Saving ', save_name)
            tas_fut_corrected.append(corrected_cube)
    return tas_fut_corrected

e_path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/tas/end/'
m_path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/tas/mid/'

tas_mid = apply_cor(tas_mid, m_path)

tas_end_rcp45 = [x for x in tas_future if x.coord('sim').points[0] == 'rcp45']
tas_end = apply_cor(tas_end_rcp45, e_path)
#%%
# check
#new = corrected_cube.aggregated_by(['month'], iris.analysis.MEAN)
#old = cube.aggregated_by(['month'], iris.analysis.MEAN)

#his tas already saved
#for cube in tas_his_corrected:
#    try:
#        save_name = 'tas' + '_' + tp.gcm(cube) + '_' + tp.model(cube) + '_' + 'historical'
#        save_path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/tas/' + save_name + '.nc'
#        iris.save(cube, save_path)
#    except:
#        print('Will remove some metadata when saving...')
#        pass
#    
for cube in tas_fut_corrected:
    try:
        save_name = 'tas' + '_' + tp.gcm(cube) + '_' + tp.model(cube) + '_' + 'rcp85'
        save_path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/tas/' + save_name + '.nc'
        iris.save(cube, save_path)
    except:
        print('Will remove some metadata when saving...')
        pass

#%% Check tas correction
#% Plot correction factors
        
cube_list = p_cor[24]

#cube_list = dif_cor
rows = 3
cols = 4
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)


fig = plt.figure(figsize=(9,9))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, 12):
    ax = ax_list[i]
    cube = cube_list[i]
    
    im_dif = iplt.contourf(cube, axes = ax, cmap = 'RdBu', levels = np.arange(-5, 6, 1), extend = 'both')
    
    #title = cube.coord('model').points[0] 
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True)
       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = 'Jan', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = 'Feb', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = 'March', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[3].annotate(s = 'April', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[4].annotate(s = 'May', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[5].annotate(s = 'June', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[6].annotate(s = 'July', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[7].annotate(s = 'Aug', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[8].annotate(s = 'Sept', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[9].annotate(s = 'Oct', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[10].annotate(s = 'Nov', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[11].annotate(s = 'Dec', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')



cbax_dif = tp.get_cbax(fig, ax = ax_list[8], last_ax = [ax_list[-1]], orientation = 'horizontal',
                       dif = 0.06, h_w = 0.02)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'horizontal')
#cbar_dif.set_ticks([-40,  -30,  -20, -10, 0, 10, 20, 30, 40])
cbar_dif.set_label('Correction Factor')

#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Tanga Project/Figures/CORDEX/Hist/Monthly_rainfall/historical_trmmdif_bias_correted_ens_levs2.png',
#            bbox_inches = 'tight', pad_inches = 0.3)

#%% Check - monthly avgs should match CRU tas
        
avg_tas = iris.cube.CubeList()
for cube in tas_his_corrected:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas.append(x)
  
cru_cor = cru_month[1][0]
cru_um = cru_month[0][0]

#split cordex into list for each month
cor = [x for x in avg_tas if tp.gcm(x) != 'um']
cor_months = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in cor:
        x = cube[i]
        new_month.append(x)
    cor_months.append(new_month)

dif_cor = iris.cube.CubeList()
for i in np.arange(0, len(cor_months)):
    month = cor_months[i]
    month_dif = iris.cube.CubeList()
    for cube in month:
        dif = copy.deepcopy(cube)
        dif.data = cube.data - cru_cor[i].data
        month_dif.append(dif)
    dif_cor.append(month_dif)

for cube_list in dif_cor:
    for cube in cube_list:
        print(np.min(cube.data), np.max(cube.data))

#check cpr and p25
um = [x for x in avg_tas if tp.gcm(x) == 'um']
um_months = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in um:
        x = cube[i]
        new_month.append(x)
    um_months.append(new_month)

dif_um = iris.cube.CubeList()
for i in np.arange(0, len(um_months)):
    month = um_months[i]
    month_dif = iris.cube.CubeList()
    for cube in month:
        dif = copy.deepcopy(cube)
        dif.data = cube.data - cru_um[i].data
        month_dif.append(dif)
    dif_um.append(month_dif)

for cube_list in dif_um:
    for cube in cube_list:
        print(np.min(cube.data), np.max(cube.data))
        
#check no negative numbers in bias-corrected data
for cube in dif_cor:
    for x in cube:
        print(np.min(x.data), np.max(x.data))

#%% Create plot
cube_list = iris.cube.CubeList()
for i in np.arange(0, 12):
    for cube in dif_um[i]:
        try:
            cube.remove_coord('time')
        except:
            print('no time')
        cube.cell_methods = ''
    ens = ens_mean(dif_cor[i])
    cube_list.append(ens)

#cube_list = dif_cor
rows = 3
cols = 4
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)


fig = plt.figure(figsize=(9,9))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, len(cube_list)):
    ax = ax_list[i]
    cube = cube_list[i]
    
    im_dif = iplt.contourf(cube, axes = ax, cmap = 'RdBu', levels = np.arange(-0.01, 0.01,0.001), extend = 'both')
    
    #title = cube.coord('model').points[0] 
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True)
       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = 'Jan', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = 'Feb', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = 'March', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[3].annotate(s = 'April', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[4].annotate(s = 'May', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[5].annotate(s = 'June', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[6].annotate(s = 'July', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[7].annotate(s = 'Aug', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[8].annotate(s = 'Sept', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[9].annotate(s = 'Oct', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[10].annotate(s = 'Nov', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[11].annotate(s = 'Dec', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')



cbax_dif = tp.get_cbax(fig, ax = ax_list[8], last_ax = [ax_list[-1]], orientation = 'horizontal',
                       dif = 0.06, h_w = 0.02)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'horizontal')
#cbar_dif.set_ticks([-40,  -30,  -20, -10, 0, 10, 20, 30, 40])
cbar_dif.set_label('Model - CRU')

#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Tanga Project/Figures/CORDEX/Hist/Monthly_rainfall/historical_trmmdif_bias_correted_ens_levs2.png',
#            bbox_inches = 'tight', pad_inches = 0.3)

#%% Compare monthly change in temp

tas_his_corrected = sorted(tas_his_corrected, key = lambda o: (o.coord('gcm').points[0],
                                                               o.coord('model').points[0]))
tas_fut_corrected = sorted(tas_fut_corrected, key = lambda o: (o.coord('gcm').points[0],
                                                               o.coord('model').points[0]))
tas_his[0] = sorted(tas_his[0], key = lambda o: (o.coord('gcm').points[0],
                                                               o.coord('model').points[0]))
tas_future[0] = sorted(tas_future[0], key = lambda o: (o.coord('gcm').points[0],
                                                               o.coord('model').points[0]))

cube_list = [tas_his_corrected, tas_fut_corrected, tas_his[0], tas_future[0]]
for item in cube_list:
    for cube in item:
        cube.convert_units('celsius')

avg_tas_his_bc = iris.cube.CubeList()
for cube in tas_his_corrected:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas_his_bc.append(x)
  
avg_tas_fut_bc = iris.cube.CubeList()
for cube in tas_fut_corrected:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas_fut_bc.append(x)   
    
#split models into list for each month
his_months_bc = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in avg_tas_his_bc:
        x = cube[i]
        new_month.append(x)
    his_months_bc.append(new_month)

fut_months_bc = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in avg_tas_fut_bc:
        x = cube[i]
        new_month.append(x)
    fut_months_bc.append(new_month)

dif_bc = iris.cube.CubeList()
for i in np.arange(0, len(his_months_bc)):
    pc_month = his_months_bc[i]
    fc_month = fut_months_bc[i]
    month_dif = iris.cube.CubeList()
    for j in np.arange(0, len(pc_month)):
        pc_mod = pc_month[j]
        fc_mod = fc_month[j]
        pc_model = tp.gcm(pc_mod) + tp.model(pc_mod)
        fc_model = tp.gcm(fc_mod) + tp.model(fc_mod)
        if pc_model == fc_model:
            month_dif.append(fc_mod - pc_mod)
        else:
            print('Model mismatch:', i, pc_model, fc_model)
    dif_bc.append(month_dif)

#%%
# non bias-corrected
avg_tas_his = iris.cube.CubeList()
for cube in tas_his[0]:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas_his.append(x)
  
avg_tas_fut = iris.cube.CubeList()
for cube in tas_future[0]:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas_fut.append(x)   
    
#split models into list for each month
his_months = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in avg_tas_his:
        x = cube[i]
        new_month.append(x)
    his_months.append(new_month)

fut_months = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in avg_tas_fut:
        x = cube[i]
        new_month.append(x)
    fut_months.append(new_month)

dif_orig = iris.cube.CubeList()
for i in np.arange(0, len(his_months)):
    pc_month = his_months[i]
    fc_month = fut_months[i]
    month_dif = iris.cube.CubeList()
    if pc_month[0].coord('month').points[0] == fc_month[0].coord('month').points[0]:
        for j in np.arange(0, len(pc_month)):
            pc_mod = pc_month[j]
            fc_mod = fc_month[j]
            pc_model = tp.gcm(pc_mod) + tp.model(pc_mod)
            fc_model = tp.gcm(fc_mod) + tp.model(fc_mod)
            if pc_model == fc_model:
                month_dif.append(fc_mod - pc_mod)
            else:
                print('Model mismatch:', pc_model, fc_model)
        dif_orig.append(month_dif)
    else:
        print('Month mismatch', pc_month[0].coord('month').points[0], fc_month[0].coord('month').points[0])


#%%
cube_list = iris.cube.CubeList()
for i in np.arange(0, 12):
    for cube in dif_orig[i]:
        try:
            cube.remove_coord('time')
        except:
            print('no time')
        cube.cell_methods = ''
    ens = ens_mean(dif_orig[i])
    cube_list.append(ens)


rows = 3
cols = 4
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)


fig = plt.figure(figsize=(9,9))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, len(cube_list)):
    ax = ax_list[i]
    cube = cube_list[i]
    
    im_dif = iplt.contourf(cube, axes = ax, cmap = 'Reds', levels = np.arange(0,8,1), extend = 'max')
    
    #title = cube.coord('model').points[0] 
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True)
       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = 'Jan', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = 'Feb', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = 'March', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[3].annotate(s = 'April', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[4].annotate(s = 'May', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[5].annotate(s = 'June', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[6].annotate(s = 'July', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[7].annotate(s = 'Aug', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[8].annotate(s = 'Sept', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[9].annotate(s = 'Oct', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[10].annotate(s = 'Nov', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[11].annotate(s = 'Dec', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')



cbax_dif = tp.get_cbax(fig, ax = ax_list[8], last_ax = [ax_list[-1]], orientation = 'horizontal',
                       dif = 0.06, h_w = 0.02)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'horizontal')
#cbar_dif.set_ticks([-80, -60, -40, -20,  0,  20, 40, 60, 80])
cbar_dif.set_label('Future - Present')

#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Tanga Project/Figures/CORDEX/BiasCorr/TAS_RCP85_MonthlyComparison_orig_cordexens.png',
#            bbox_inches = 'tight', pad_inches = 0.3)

#%% check cp4 and p25

cube_list = iris.cube.CubeList()
for i in np.arange(0, 12):
    for cube in dif_orig[i]:
        x = [x for x in dif_orig[i] if tp.model(x) == 'p25'][0]
    cube_list.append(x)


rows = 3
cols = 4
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)


fig = plt.figure(figsize=(9,9))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, len(cube_list)):
    ax = ax_list[i]
    cube = cube_list[i]
    
    im_dif = iplt.contourf(cube, axes = ax, cmap = 'Reds', levels = np.arange(0,8,1), extend = 'max')
    
    #title = cube.coord('model').points[0] 
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True)
       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = 'Jan', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = 'Feb', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = 'March', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[3].annotate(s = 'April', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[4].annotate(s = 'May', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[5].annotate(s = 'June', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[6].annotate(s = 'July', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[7].annotate(s = 'Aug', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[8].annotate(s = 'Sept', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[9].annotate(s = 'Oct', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[10].annotate(s = 'Nov', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[11].annotate(s = 'Dec', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')



cbax_dif = tp.get_cbax(fig, ax = ax_list[8], last_ax = [ax_list[-1]], orientation = 'horizontal',
                       dif = 0.06, h_w = 0.02)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'horizontal')
#cbar_dif.set_ticks([-80, -60, -40, -20,  0,  20, 40, 60, 80])
cbar_dif.set_label('Future - Present')

#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Tanga Project/Figures/CORDEX/BiasCorr/TAS_RCP85_MonthlyComparison_orig_cordexens.png',
#            bbox_inches = 'tight', pad_inches = 0.3)


#%%
    
''' TMAX '''
    
    
'''
Tmodel(day) = Tmodel(day) + MEANmonth(Pobs) - MEANmonth(Pmodel)

For each model and each month calculate a correction factor        
'''
p_cor = iris.cube.CubeList()

for cube in tas_months[1]:
    print(tas_months[1][0].var_name)
    if tp.gcm(cube) == 'um':
        obs_orig = cru_month[0][2]
        if np.min(obs_orig.coord('year').bounds) != 1997:
            print('Wrong cru file.')
        print(tp.model(cube), 'using 1997 chirps')
    else:
        obs_orig = cru_month[1][2]
        print(obs_orig.var_name)
        if np.min(obs_orig.coord('year').bounds) != 1971:
            print('Wrong cru file.')
    
    model = copy.deepcopy(cube)
    obs = copy.deepcopy(obs_orig)
    
    #If historical model data is zero will lead to correction factor being NA and introducting
    # NAS. If replace future corrected data with 0 where 0 in historical means can't get more rainfall
    # in very dry areas.
    model.data = obs.data - model.data
    p_cor.append(model)
    
#apply correction factor
    
month_dict = {'Jan': 0, 'Feb': 1, 'Mar':2, 'Apr':3, 'May':4, 'Jun': 5, 'Jul':6, 
           'Aug':7, 'Sep':8, 'Oct': 9, 'Nov': 10, 'Dec': 11}

tas_his_corrected = iris.cube.CubeList()

for cube in tas_his[1]:
    #Get correction cube
    mod = tp.gcm(cube) + tp.model(cube)
    for corr in p_cor:
        cor_mod = tp.gcm(corr) + tp.model(corr)
        if cor_mod == mod:
            cor_cube = corr
            print('Correction cube found:', mod, cor_mod)
            break
    #Apply correction factor to daily data
    if cor_mod == mod:
        dims = cube.shape
        corrected_cube = copy.deepcopy(cube)
        #cycle through days and apply monthly correction factor to each day
        for t in np.arange(0, dims[0]):
            month = cube[t].coord('month').points[0] #find month
            month_num = month_dict[month] #get number associated with month
            corr_month = cor_cube[month_num] #get correction factor associated with month
            if corr_month.coord('month').points[0] == cube[t].coord('month').points[0]:
                new_day = copy.deepcopy(cube[t])
                new_day.data = new_day.data + corr_month.data
                corrected_cube.data[t] = new_day.data
            else:
                print('Month mismatch')
        tas_his_corrected.append(corrected_cube)
    else:
        print('Model mismatch:', cor_mod, mod)
       
        
tas_fut_corrected = iris.cube.CubeList()
for cube in tas_future[1]:
    #Get correction cube
    mod = tp.gcm(cube) + tp.model(cube)
    for corr in p_cor:
        cor_mod = tp.gcm(corr) + tp.model(corr)
        if cor_mod == mod:
            cor_cube = corr
            print('Correction cube found:', mod, cor_mod)
            break
    if cor_mod == mod:
        #Apply correction factor to daily data
        dims = cube.shape
        corrected_cube = copy.deepcopy(cube)
        #cycle through days and apply monthly correction factor to each day
        for t in np.arange(0, dims[0]):
            month = cube[t].coord('month').points[0] #find month
            month_num = month_dict[month] #get number associated with month
            corr_month = cor_cube[month_num] #get correction factor associated with month
            if corr_month.coord('month').points[0] == cube[t].coord('month').points[0]:
                new_day = copy.deepcopy(cube[t])
                new_day.data = new_day.data + corr_month.data
                corrected_cube.data[t] = new_day.data
            else:
                print('Month mismatch')
        tas_fut_corrected.append(corrected_cube)

#%%
# check
#new = corrected_cube.aggregated_by(['month'], iris.analysis.MEAN)
#old = cube.aggregated_by(['month'], iris.analysis.MEAN)

for cube in tas_his_corrected:
    try:
        save_name = 'tasmax' + '_' + tp.gcm(cube) + '_' + tp.model(cube) + '_' + 'historical'
        save_path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/tasmax/' + save_name + '.nc'
        iris.save(cube, save_path)
    except:
        print('Will remove some metadata when saving...')
        pass
    
for cube in tas_fut_corrected:
    try:
        save_name = 'tasmax' + '_' + tp.gcm(cube) + '_' + tp.model(cube) + '_' + 'rcp85'
        save_path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/tasmax/' + save_name + '.nc'
        iris.save(cube, save_path)
    except:
        print('Will remove some metadata when saving...')
        pass

#%% Check tas correction
#% Plot correction factors
        
cube_list = p_cor[24]

#cube_list = dif_cor
rows = 3
cols = 4
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)


fig = plt.figure(figsize=(9,9))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, 12):
    ax = ax_list[i]
    cube = cube_list[i]
    
    im_dif = iplt.contourf(cube, axes = ax, cmap = 'RdBu', levels = np.arange(-5, 6, 1), extend = 'both')
    
    #title = cube.coord('model').points[0] 
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True)
       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = 'Jan', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = 'Feb', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = 'March', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[3].annotate(s = 'April', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[4].annotate(s = 'May', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[5].annotate(s = 'June', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[6].annotate(s = 'July', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[7].annotate(s = 'Aug', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[8].annotate(s = 'Sept', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[9].annotate(s = 'Oct', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[10].annotate(s = 'Nov', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[11].annotate(s = 'Dec', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')



cbax_dif = tp.get_cbax(fig, ax = ax_list[8], last_ax = [ax_list[-1]], orientation = 'horizontal',
                       dif = 0.06, h_w = 0.02)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'horizontal')
#cbar_dif.set_ticks([-40,  -30,  -20, -10, 0, 10, 20, 30, 40])
cbar_dif.set_label('Correction Factor')

#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Tanga Project/Figures/CORDEX/Hist/Monthly_rainfall/historical_trmmdif_bias_correted_ens_levs2.png',
#            bbox_inches = 'tight', pad_inches = 0.3)

#%% Check - monthly avgs should match CRU tas
        
avg_tas = iris.cube.CubeList()
for cube in tas_his_corrected:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas.append(x)
  
cru_cor = cru_month[1][2]
cru_um = cru_month[0][2]

#split cordex into list for each month
cor = [x for x in avg_tas if tp.gcm(x) != 'um']
cor_months = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in cor:
        x = cube[i]
        new_month.append(x)
    cor_months.append(new_month)

dif_cor = iris.cube.CubeList()
for i in np.arange(0, len(cor_months)):
    month = cor_months[i]
    month_dif = iris.cube.CubeList()
    for cube in month:
        dif = copy.deepcopy(cube)
        dif.data = cube.data - cru_cor[i].data
        month_dif.append(dif)
    dif_cor.append(month_dif)

for cube_list in dif_cor:
    for cube in cube_list:
        print(np.min(cube.data), np.max(cube.data))

#check cpr and p25
um = [x for x in avg_tas if tp.gcm(x) == 'um']
um_months = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in um:
        x = cube[i]
        new_month.append(x)
    um_months.append(new_month)

dif_um = iris.cube.CubeList()
for i in np.arange(0, len(um_months)):
    month = um_months[i]
    month_dif = iris.cube.CubeList()
    for cube in month:
        dif = copy.deepcopy(cube)
        dif.data = cube.data - cru_um[i].data
        month_dif.append(dif)
    dif_um.append(month_dif)

for cube_list in dif_um:
    for cube in cube_list:
        print(np.min(cube.data), np.max(cube.data))
        
#check no negative numbers in bias-corrected data
for cube in dif_cor:
    for x in cube:
        print(np.min(x.data), np.max(x.data))

#%% Create plot
cube_list = iris.cube.CubeList()
for i in np.arange(0, 12):
    for cube in dif_cor[i]:
        try:
            cube.remove_coord('time')
        except:
            print('no time')
        cube.cell_methods = ''
    ens = ens_mean(dif_cor[i])
    cube_list.append(ens)

#cube_list = dif_cor
rows = 3
cols = 4
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)


fig = plt.figure(figsize=(9,9))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, len(cube_list)):
    ax = ax_list[i]
    cube = cube_list[i]
    
    im_dif = iplt.contourf(cube, axes = ax, cmap = 'RdBu', levels = np.arange(-0.01, 0.01,0.001), extend = 'both')
    
    #title = cube.coord('model').points[0] 
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True)
       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = 'Jan', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = 'Feb', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = 'March', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[3].annotate(s = 'April', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[4].annotate(s = 'May', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[5].annotate(s = 'June', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[6].annotate(s = 'July', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[7].annotate(s = 'Aug', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[8].annotate(s = 'Sept', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[9].annotate(s = 'Oct', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[10].annotate(s = 'Nov', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[11].annotate(s = 'Dec', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')



cbax_dif = tp.get_cbax(fig, ax = ax_list[8], last_ax = [ax_list[-1]], orientation = 'horizontal',
                       dif = 0.06, h_w = 0.02)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'horizontal')
#cbar_dif.set_ticks([-40,  -30,  -20, -10, 0, 10, 20, 30, 40])
cbar_dif.set_label('Model - CRU')

#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Tanga Project/Figures/CORDEX/Hist/Monthly_rainfall/historical_trmmdif_bias_correted_ens_levs2.png',
#            bbox_inches = 'tight', pad_inches = 0.3)

#%% Compare monthly change in temp

tas_his_corrected = sorted(tas_his_corrected, key = lambda o: (o.coord('gcm').points[0],
                                                               o.coord('model').points[0]))
tas_fut_corrected = sorted(tas_fut_corrected, key = lambda o: (o.coord('gcm').points[0],
                                                               o.coord('model').points[0]))
tas_his[1] = sorted(tas_his[1], key = lambda o: (o.coord('gcm').points[0],
                                                               o.coord('model').points[0]))
tas_future[1] = sorted(tas_future[1], key = lambda o: (o.coord('gcm').points[0],
                                                               o.coord('model').points[0]))

cube_list = [tas_his_corrected, tas_fut_corrected, tas_his[1], tas_future[1]]
for item in cube_list:
    for cube in item:
        cube.convert_units('celsius')

avg_tas_his_bc = iris.cube.CubeList()
for cube in tas_his_corrected:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas_his_bc.append(x)
  
avg_tas_fut_bc = iris.cube.CubeList()
for cube in tas_fut_corrected:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas_fut_bc.append(x)   
    
#split models into list for each month
his_months_bc = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in avg_tas_his_bc:
        x = cube[i]
        new_month.append(x)
    his_months_bc.append(new_month)

fut_months_bc = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in avg_tas_fut_bc:
        x = cube[i]
        new_month.append(x)
    fut_months_bc.append(new_month)

dif_bc = iris.cube.CubeList()
for i in np.arange(0, len(his_months_bc)):
    pc_month = his_months_bc[i]
    fc_month = fut_months_bc[i]
    month_dif = iris.cube.CubeList()
    for j in np.arange(0, len(pc_month)):
        pc_mod = pc_month[j]
        fc_mod = fc_month[j]
        pc_model = tp.gcm(pc_mod) + tp.model(pc_mod)
        fc_model = tp.gcm(fc_mod) + tp.model(fc_mod)
        if pc_model == fc_model:
            month_dif.append(fc_mod - pc_mod)
        else:
            print('Model mismatch:', i, pc_model, fc_model)
    dif_bc.append(month_dif)

#%%
# non bias-corrected
avg_tas_his = iris.cube.CubeList()
for cube in tas_his[1]:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas_his.append(x)
  
avg_tas_fut = iris.cube.CubeList()
for cube in tas_future[1]:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas_fut.append(x)   
    
#split models into list for each month
his_months = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in avg_tas_his:
        x = cube[i]
        new_month.append(x)
    his_months.append(new_month)

fut_months = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in avg_tas_fut:
        x = cube[i]
        new_month.append(x)
    fut_months.append(new_month)

dif_orig = iris.cube.CubeList()
for i in np.arange(0, len(his_months)):
    pc_month = his_months[i]
    fc_month = fut_months[i]
    month_dif = iris.cube.CubeList()
    if pc_month[0].coord('month').points[0] == fc_month[0].coord('month').points[0]:
        for j in np.arange(0, len(pc_month)):
            pc_mod = pc_month[j]
            fc_mod = fc_month[j]
            pc_model = tp.gcm(pc_mod) + tp.model(pc_mod)
            fc_model = tp.gcm(fc_mod) + tp.model(fc_mod)
            if pc_model == fc_model:
                month_dif.append(fc_mod - pc_mod)
            else:
                print('Model mismatch:', pc_model, fc_model)
        dif_orig.append(month_dif)
    else:
        print('Month mismatch', pc_month[0].coord('month').points[0], fc_month[0].coord('month').points[0])


#%%
cube_list = iris.cube.CubeList()
for i in np.arange(0, 12):
    for cube in dif_bc[i]:
        try:
            cube.remove_coord('time')
        except:
            print('no time')
        cube.cell_methods = ''
    ens = ens_mean(dif_bc[i])
    cube_list.append(ens)


rows = 3
cols = 4
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)


fig = plt.figure(figsize=(9,9))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, len(cube_list)):
    ax = ax_list[i]
    cube = cube_list[i]
    
    im_dif = iplt.contourf(cube, axes = ax, cmap = 'Reds', levels = np.arange(0,8,1), extend = 'max')
    
    #title = cube.coord('model').points[0] 
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True)
       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = 'Jan', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = 'Feb', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = 'March', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[3].annotate(s = 'April', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[4].annotate(s = 'May', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[5].annotate(s = 'June', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[6].annotate(s = 'July', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[7].annotate(s = 'Aug', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[8].annotate(s = 'Sept', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[9].annotate(s = 'Oct', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[10].annotate(s = 'Nov', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[11].annotate(s = 'Dec', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')



cbax_dif = tp.get_cbax(fig, ax = ax_list[8], last_ax = [ax_list[-1]], orientation = 'horizontal',
                       dif = 0.06, h_w = 0.02)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'horizontal')
#cbar_dif.set_ticks([-80, -60, -40, -20,  0,  20, 40, 60, 80])
cbar_dif.set_label('Future - Present')

fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Tanga Project/Figures/CORDEX/BiasCorr/TMAX_RCP85_MonthlyComparison_corrected_cordexens.png',
            bbox_inches = 'tight', pad_inches = 0.3)

#%% check cp4 and p25

cube_list = iris.cube.CubeList()
for i in np.arange(0, 12):
    for cube in dif_orig[i]:
        x = [x for x in dif_orig[i] if tp.model(x) == 'cp4'][0]
    cube_list.append(x)


rows = 3
cols = 4
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)


fig = plt.figure(figsize=(9,9))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, len(cube_list)):
    ax = ax_list[i]
    cube = cube_list[i]
    
    im_dif = iplt.contourf(cube, axes = ax, cmap = 'Reds', levels = np.arange(0,8,1), extend = 'max')
    
    #title = cube.coord('model').points[0] 
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True)
       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = 'Jan', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = 'Feb', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = 'March', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[3].annotate(s = 'April', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[4].annotate(s = 'May', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[5].annotate(s = 'June', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[6].annotate(s = 'July', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[7].annotate(s = 'Aug', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[8].annotate(s = 'Sept', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[9].annotate(s = 'Oct', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[10].annotate(s = 'Nov', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[11].annotate(s = 'Dec', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')



cbax_dif = tp.get_cbax(fig, ax = ax_list[8], last_ax = [ax_list[-1]], orientation = 'horizontal',
                       dif = 0.06, h_w = 0.02)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'horizontal')
#cbar_dif.set_ticks([-80, -60, -40, -20,  0,  20, 40, 60, 80])
cbar_dif.set_label('Future - Present')

#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Tanga Project/Figures/CORDEX/BiasCorr/TMAX_RCP85_MonthlyComparison_orig_cordexens.png',
#            bbox_inches = 'tight', pad_inches = 0.3)

#%%
    
''' TMIN '''
    
    
'''
Tmodel(day) = Tmodel(day) + MEANmonth(Pobs) - MEANmonth(Pmodel)

For each model and each month calculate a correction factor        
'''
p_cor = iris.cube.CubeList()

for cube in tas_months[2]:
    print(cube.var_name)
    if tp.gcm(cube) == 'um':
        obs_orig = cru_month[0][1]
        if np.min(obs_orig.coord('year').bounds) != 1997:
            print('Wrong cru file.')
        print(tp.model(cube), 'using 1997 cru')
    else:
        obs_orig = cru_month[1][1]
        print(obs_orig.var_name)
        if np.min(obs_orig.coord('year').bounds) != 1971:
            print('Wrong cru file.')
    
    model = copy.deepcopy(cube)
    obs = copy.deepcopy(obs_orig)
    
    #If historical model data is zero will lead to correction factor being NA and introducting
    # NAS. If replace future corrected data with 0 where 0 in historical means can't get more rainfall
    # in very dry areas.
    model.data = obs.data - model.data
    p_cor.append(model)
    
#apply correction factor
    
month_dict = {'Jan': 0, 'Feb': 1, 'Mar':2, 'Apr':3, 'May':4, 'Jun': 5, 'Jul':6, 
           'Aug':7, 'Sep':8, 'Oct': 9, 'Nov': 10, 'Dec': 11}

tas_his_corrected = iris.cube.CubeList()

for cube in tas_his[2]:
    #Get correction cube
    mod = tp.gcm(cube) + tp.model(cube)
    for corr in p_cor:
        cor_mod = tp.gcm(corr) + tp.model(corr)
        if cor_mod == mod:
            cor_cube = corr
            print('Correction cube found:', mod, cor_mod)
            break
    #Apply correction factor to daily data
    if cor_mod == mod:
        dims = cube.shape
        corrected_cube = copy.deepcopy(cube)
        #cycle through days and apply monthly correction factor to each day
        for t in np.arange(0, dims[0]):
            month = cube[t].coord('month').points[0] #find month
            month_num = month_dict[month] #get number associated with month
            corr_month = cor_cube[month_num] #get correction factor associated with month
            if corr_month.coord('month').points[0] == cube[t].coord('month').points[0]:
                new_day = copy.deepcopy(cube[t])
                new_day.data = new_day.data + corr_month.data
                corrected_cube.data[t] = new_day.data
            else:
                print('Month mismatch')
        tas_his_corrected.append(corrected_cube)
    else:
        print('Model mismatch:', cor_mod, mod)
       
        
tas_fut_corrected = iris.cube.CubeList()
for cube in tas_future[2]:
    #Get correction cube
    mod = tp.gcm(cube) + tp.model(cube)
    for corr in p_cor:
        cor_mod = tp.gcm(corr) + tp.model(corr)
        if cor_mod == mod:
            cor_cube = corr
            print('Correction cube found:', mod, cor_mod)
            break
    if cor_mod == mod:
        #Apply correction factor to daily data
        dims = cube.shape
        corrected_cube = copy.deepcopy(cube)
        #cycle through days and apply monthly correction factor to each day
        for t in np.arange(0, dims[0]):
            month = cube[t].coord('month').points[0] #find month
            month_num = month_dict[month] #get number associated with month
            corr_month = cor_cube[month_num] #get correction factor associated with month
            if corr_month.coord('month').points[0] == cube[t].coord('month').points[0]:
                new_day = copy.deepcopy(cube[t])
                new_day.data = new_day.data + corr_month.data
                corrected_cube.data[t] = new_day.data
            else:
                print('Month mismatch')
        tas_fut_corrected.append(corrected_cube)

#%%
# check
#new = corrected_cube.aggregated_by(['month'], iris.analysis.MEAN)
#old = cube.aggregated_by(['month'], iris.analysis.MEAN)

for cube in tas_his_corrected:
    try:
        save_name = 'tasmin' + '_' + tp.gcm(cube) + '_' + tp.model(cube) + '_' + 'historical'
        save_path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/tasmin/' + save_name + '.nc'
        iris.save(cube, save_path)
    except:
        print('Will remove some metadata when saving...')
        pass
    
for cube in tas_fut_corrected:
    try:
        save_name = 'tasmin' + '_' + tp.gcm(cube) + '_' + tp.model(cube) + '_' + 'rcp85'
        save_path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/tasmin/' + save_name + '.nc'
        iris.save(cube, save_path)
    except:
        print('Will remove some metadata when saving...')
        pass

#%% Check tas correction
#% Plot correction factors
        
cube_list = p_cor[24]

#cube_list = dif_cor
rows = 3
cols = 4
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)


fig = plt.figure(figsize=(9,9))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, 12):
    ax = ax_list[i]
    cube = cube_list[i]
    
    im_dif = iplt.contourf(cube, axes = ax, cmap = 'RdBu', levels = np.arange(-5, 6, 1), extend = 'both')
    
    #title = cube.coord('model').points[0] 
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True)
       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = 'Jan', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = 'Feb', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = 'March', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[3].annotate(s = 'April', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[4].annotate(s = 'May', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[5].annotate(s = 'June', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[6].annotate(s = 'July', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[7].annotate(s = 'Aug', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[8].annotate(s = 'Sept', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[9].annotate(s = 'Oct', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[10].annotate(s = 'Nov', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[11].annotate(s = 'Dec', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')



cbax_dif = tp.get_cbax(fig, ax = ax_list[8], last_ax = [ax_list[-1]], orientation = 'horizontal',
                       dif = 0.06, h_w = 0.02)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'horizontal')
#cbar_dif.set_ticks([-40,  -30,  -20, -10, 0, 10, 20, 30, 40])
cbar_dif.set_label('Correction Factor')

#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Tanga Project/Figures/CORDEX/Hist/Monthly_rainfall/historical_trmmdif_bias_correted_ens_levs2.png',
#            bbox_inches = 'tight', pad_inches = 0.3)

#%% Check - monthly avgs should match CRU tas
        
avg_tas = iris.cube.CubeList()
for cube in tas_his_corrected:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas.append(x)
  
cru_cor = cru_month[1][1]
cru_um = cru_month[0][1]

#split cordex into list for each month
cor = [x for x in avg_tas if tp.gcm(x) != 'um']
cor_months = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in cor:
        x = cube[i]
        new_month.append(x)
    cor_months.append(new_month)

dif_cor = iris.cube.CubeList()
for i in np.arange(0, len(cor_months)):
    month = cor_months[i]
    month_dif = iris.cube.CubeList()
    for cube in month:
        dif = copy.deepcopy(cube)
        dif.data = cube.data - cru_cor[i].data
        month_dif.append(dif)
    dif_cor.append(month_dif)

for cube_list in dif_cor:
    for cube in cube_list:
        print(np.min(cube.data), np.max(cube.data))

#check cpr and p25
um = [x for x in avg_tas if tp.gcm(x) == 'um']
um_months = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in um:
        x = cube[i]
        new_month.append(x)
    um_months.append(new_month)

dif_um = iris.cube.CubeList()
for i in np.arange(0, len(um_months)):
    month = um_months[i]
    month_dif = iris.cube.CubeList()
    for cube in month:
        dif = copy.deepcopy(cube)
        dif.data = cube.data - cru_um[i].data
        month_dif.append(dif)
    dif_um.append(month_dif)

for cube_list in dif_um:
    for cube in cube_list:
        print(np.min(cube.data), np.max(cube.data))
        
#check no negative numbers in bias-corrected data
for cube in dif_cor:
    for x in cube:
        print(np.min(x.data), np.max(x.data))

#%% Create plot
cube_list = iris.cube.CubeList()
for i in np.arange(0, 12):
    for cube in dif_cor[i]:
        try:
            cube.remove_coord('time')
        except:
            print('no time')
        cube.cell_methods = ''
    ens = ens_mean(dif_cor[i])
    cube_list.append(ens)

#cube_list = dif_cor
rows = 3
cols = 4
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)


fig = plt.figure(figsize=(9,9))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, len(cube_list)):
    ax = ax_list[i]
    cube = cube_list[i]
    
    im_dif = iplt.contourf(cube, axes = ax, cmap = 'RdBu', levels = np.arange(-0.01, 0.01,0.001), extend = 'both')
    
    #title = cube.coord('model').points[0] 
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True)
       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = 'Jan', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = 'Feb', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = 'March', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[3].annotate(s = 'April', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[4].annotate(s = 'May', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[5].annotate(s = 'June', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[6].annotate(s = 'July', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[7].annotate(s = 'Aug', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[8].annotate(s = 'Sept', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[9].annotate(s = 'Oct', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[10].annotate(s = 'Nov', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[11].annotate(s = 'Dec', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')



cbax_dif = tp.get_cbax(fig, ax = ax_list[8], last_ax = [ax_list[-1]], orientation = 'horizontal',
                       dif = 0.06, h_w = 0.02)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'horizontal')
#cbar_dif.set_ticks([-40,  -30,  -20, -10, 0, 10, 20, 30, 40])
cbar_dif.set_label('Model - CRU')

#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Tanga Project/Figures/CORDEX/Hist/Monthly_rainfall/historical_trmmdif_bias_correted_ens_levs2.png',
#            bbox_inches = 'tight', pad_inches = 0.3)

#%% Compare monthly change in temp

tas_his_corrected = sorted(tas_his_corrected, key = lambda o: (o.coord('gcm').points[0],
                                                               o.coord('model').points[0]))
tas_fut_corrected = sorted(tas_fut_corrected, key = lambda o: (o.coord('gcm').points[0],
                                                               o.coord('model').points[0]))
tas_his[2] = sorted(tas_his[2], key = lambda o: (o.coord('gcm').points[0],
                                                               o.coord('model').points[0]))
tas_future[2] = sorted(tas_future[2], key = lambda o: (o.coord('gcm').points[0],
                                                               o.coord('model').points[0]))

cube_list = [tas_his_corrected, tas_fut_corrected, tas_his[2], tas_future[2]]
for item in cube_list:
    for cube in item:
        cube.convert_units('celsius')

avg_tas_his_bc = iris.cube.CubeList()
for cube in tas_his_corrected:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas_his_bc.append(x)
  
avg_tas_fut_bc = iris.cube.CubeList()
for cube in tas_fut_corrected:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas_fut_bc.append(x)   
    
#split models into list for each month
his_months_bc = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in avg_tas_his_bc:
        x = cube[i]
        new_month.append(x)
    his_months_bc.append(new_month)

fut_months_bc = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in avg_tas_fut_bc:
        x = cube[i]
        new_month.append(x)
    fut_months_bc.append(new_month)

dif_bc = iris.cube.CubeList()
for i in np.arange(0, len(his_months_bc)):
    pc_month = his_months_bc[i]
    fc_month = fut_months_bc[i]
    month_dif = iris.cube.CubeList()
    for j in np.arange(0, len(pc_month)):
        pc_mod = pc_month[j]
        fc_mod = fc_month[j]
        pc_model = tp.gcm(pc_mod) + tp.model(pc_mod)
        fc_model = tp.gcm(fc_mod) + tp.model(fc_mod)
        if pc_model == fc_model:
            month_dif.append(fc_mod - pc_mod)
        else:
            print('Model mismatch:', i, pc_model, fc_model)
    dif_bc.append(month_dif)

#%%
# non bias-corrected
avg_tas_his = iris.cube.CubeList()
for cube in tas_his[2]:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas_his.append(x)
  
avg_tas_fut = iris.cube.CubeList()
for cube in tas_future[2]:
    x = cube.aggregated_by('month', iris.analysis.MEAN)
    avg_tas_fut.append(x)   
    
#split models into list for each month
his_months = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in avg_tas_his:
        x = cube[i]
        new_month.append(x)
    his_months.append(new_month)

fut_months = iris.cube.CubeList()
for i in np.arange(0, 12):
    new_month = iris.cube.CubeList()
    for cube in avg_tas_fut:
        x = cube[i]
        new_month.append(x)
    fut_months.append(new_month)

dif_orig = iris.cube.CubeList()
for i in np.arange(0, len(his_months)):
    pc_month = his_months[i]
    fc_month = fut_months[i]
    month_dif = iris.cube.CubeList()
    if pc_month[0].coord('month').points[0] == fc_month[0].coord('month').points[0]:
        for j in np.arange(0, len(pc_month)):
            pc_mod = pc_month[j]
            fc_mod = fc_month[j]
            pc_model = tp.gcm(pc_mod) + tp.model(pc_mod)
            fc_model = tp.gcm(fc_mod) + tp.model(fc_mod)
            if pc_model == fc_model:
                month_dif.append(fc_mod - pc_mod)
            else:
                print('Model mismatch:', pc_model, fc_model)
        dif_orig.append(month_dif)
    else:
        print('Month mismatch', pc_month[0].coord('month').points[0], fc_month[0].coord('month').points[0])


#%%
cube_list = iris.cube.CubeList()
for i in np.arange(0, 12):
    for cube in dif_bc[i]:
        try:
            cube.remove_coord('time')
        except:
            print('no time')
        cube.cell_methods = ''
    ens = ens_mean(dif_bc[i])
    cube_list.append(ens)


rows = 3
cols = 4
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)


fig = plt.figure(figsize=(9,9))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, len(cube_list)):
    ax = ax_list[i]
    cube = cube_list[i]
    
    im_dif = iplt.contourf(cube, axes = ax, cmap = 'Reds', levels = np.arange(0,8,1), extend = 'max')
    
    #title = cube.coord('model').points[0] 
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True)
       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = 'Jan', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = 'Feb', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = 'March', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[3].annotate(s = 'April', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[4].annotate(s = 'May', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[5].annotate(s = 'June', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[6].annotate(s = 'July', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[7].annotate(s = 'Aug', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[8].annotate(s = 'Sept', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[9].annotate(s = 'Oct', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[10].annotate(s = 'Nov', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[11].annotate(s = 'Dec', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')



cbax_dif = tp.get_cbax(fig, ax = ax_list[8], last_ax = [ax_list[-1]], orientation = 'horizontal',
                       dif = 0.06, h_w = 0.02)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'horizontal')
#cbar_dif.set_ticks([-80, -60, -40, -20,  0,  20, 40, 60, 80])
cbar_dif.set_label('Future - Present')

fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Tanga Project/Figures/CORDEX/BiasCorr/TMIN_RCP85_MonthlyComparison_corrected_cordexens.png',
            bbox_inches = 'tight', pad_inches = 0.3)

#%% check cp4 and p25

cube_list = iris.cube.CubeList()
for i in np.arange(0, 12):
    for cube in dif_orig[i]:
        x = [x for x in dif_orig[i] if tp.model(x) == 'cp4'][0]
    cube_list.append(x)


rows = 3
cols = 4
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)


fig = plt.figure(figsize=(9,9))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, len(cube_list)):
    ax = ax_list[i]
    cube = cube_list[i]
    
    im_dif = iplt.contourf(cube, axes = ax, cmap = 'Reds', levels = np.arange(0,8,1), extend = 'max')
    
    #title = cube.coord('model').points[0] 
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True)
       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = 'Jan', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = 'Feb', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = 'March', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[3].annotate(s = 'April', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[4].annotate(s = 'May', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[5].annotate(s = 'June', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[6].annotate(s = 'July', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[7].annotate(s = 'Aug', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[8].annotate(s = 'Sept', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[9].annotate(s = 'Oct', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[10].annotate(s = 'Nov', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[11].annotate(s = 'Dec', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')



cbax_dif = tp.get_cbax(fig, ax = ax_list[8], last_ax = [ax_list[-1]], orientation = 'horizontal',
                       dif = 0.06, h_w = 0.02)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'horizontal')
#cbar_dif.set_ticks([-80, -60, -40, -20,  0,  20, 40, 60, 80])
cbar_dif.set_label('Future - Present')

#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Tanga Project/Figures/CORDEX/BiasCorr/TMIN_RCP85_MonthlyComparison_orig_cordexens.png',
#            bbox_inches = 'tight', pad_inches = 0.3)
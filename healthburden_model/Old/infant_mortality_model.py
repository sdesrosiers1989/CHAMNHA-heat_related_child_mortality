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
        
path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/CHAMNHA/tas/mid/'
filenames = glob.glob(path + '*.nc')
for file in filenames:
    x = iris.load_cube(file)
    tas.append(x)
        
path = '/nfs/a321/earsch/Tanga/Data/CORDEX/Bias_corr/CHAMNHA/tas/end/'
filenames = glob.glob(path + '*.nc')
for file in filenames:
    x = iris.load_cube(file)
    tas.append(x)

cube_orog = iris.load('/nfs/a277/IMPALA/data/4km/ANCILS/orog_combined_POSTac144_ancil_4km.nc')
orog = cube_orog[9] 

ulon = orog[0,0,:,:].coord('grid_longitude')
transform = ulon.coord_system.as_cartopy_projection()

#%% Placeholder data

def place_holder(base, dat):
    new_dat = copy.deepcopy(base)
    new_dat.data = np.where(~np.isnan(new_dat.data), dat, new_dat.data)
    new_dat.data = ma.masked_array(new_dat.data, mask = base.data.mask)
    return new_dat

#Coeff
coeff = place_holder(tas[0][0], 0.013346) 

#Threshold
thres = place_holder(tas[0][0], 16.5)

#baseline pop
bpop = place_holder(tas[0][0], 2519000)

#pop mid 2000
mpop = place_holder(tas[0][0], 2551915)

#import daily avg mortality palceholder
davg_mort = np.loadtxt(open('/nfs/a321/earsch/CHAMNHA/daily_avg_mortality_placeholder.csv', 'rb'),
                       delimiter = ',', skiprows = 1)

#%% calculate temp diff with threshold

def calc_tdif(tavg, thres):
    
    tdif = tavg - thres
    tdif.data = np.where(tdif.data < 0, 0, tdif.data)
    tdif.data = ma.masked_array(tdif.data, mask = tavg.data.mask)
    return tdif


tdif = calc_tdif(tas[0], thres)
#%% calculate attributable death per decade

#for earch decade
pop_ratio = mpop / bpop #ratio future pop to baseline pop

dec_start = 1990
dec_end = 1999

def ann_death_per_decade(base, dec_start, dec_end, pop_ratio):
    b_dec = base.extract(iris.Constraint(year= lambda cell: dec_start <= cell <= dec_end))
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

ann_avg_heat_death = ann_death_per_decade(tdif, dec_start, dec_end, pop_ratio)


#%% import station data
    
stat_df = pd.DataFrame({'Name': ['Korogocho', 'Viwandani', 'Nouna'],
                        'Latitude': [-1.246779, -1.304098, 12.2],
                        'Longitude': [36.89, 36.87, -4.02],
                        'Thres': [ 20, 20, 20]})
    
    

#%% Confirm locations on map

func = tp.plot_africa
  
rows = 1
cols = 1
nfigs = rows * cols

xticks = [-10, 10, 30, 50]
yticks = [-40, -20, 0, 20, 40]

fig = plt.figure(figsize=(9, 9))


ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))
    
    
im_dif = iplt.contourf(orog[0,0,:,:], axes = ax_list[0], cmap = 'BrBG', levels = np.arange(0, 2800, 200),
                           extend = 'max')

    
func(ax_list[0], high = False, xticks = xticks, yticks = yticks, Tanga = False)

ax_list[0].scatter(stat_df.Longitude,
                   stat_df.Latitude,
                   zorder = 5,
                   transform = transform,
                   facecolor = 'black')

plt.draw()

cbax_dif = tp.get_cbax(fig, ax = ax_list[0],  orientation = 'vertical',
                       dif = 0.02, h_w = 0.015)
cbar_dif = plt.colorbar(im_dif, cax = cbax_dif, orientation = 'vertical')
cbar_dif.set_label('m ASL')

#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/CHAMNHA_Orography_Stations.png',
#            bbox_inches = 'tight', pad_inches = 0.3)
    
#%% Extract CORDEX data to box and calculate daily ensemble mean
    

for cube in tas:    
    try:
        iris.coord_categorisation.add_day_of_month(cube, 'time', name = 'day')
    except:
        print('has day')
    try:
        iris.coord_categorisation.add_month_number(cube, 'time', name = 'month_num')
    except:
        print('has month num')



tas_list = sorted(tas, key = lambda o: (o.coord('gcm').points[0], o.coord('model').points[0]))

#remove MOHC models as 360 day caelndar
#tas_listx = [x for x in tas_list if 'MOHC' not in tp.gcm(x)]
#tasmax_listx = [x for x in tasmax_list if 'MOHC' not in tp.gcm(x)]
#tasmin_listx = [x for x in tasmin_list if 'MOHC' not in tp.gcm(x)]
#pr_listx = [x for x in pr_list if 'MOHC' not in tp.gcm(x)]#

#%%
    
def model_name(cube):
    name = tp.gcm(cube) + tp.model(cube) 
    return name

def cube_to_frame(tas_list, stat_name, stat_pts, thres):
    df = pd.DataFrame(columns = ['Date', 'Year', 'Month', 'Day', 'Region',
                                 'Scen', 'Model', 'Tmean', 'Threshold', 'Temp_diff'])
    
    hist_mod = np.unique([model_name(x) for x in tas_list if x.coord('sim').points[0] == 'historical'])
    
    for mod in hist_mod:
        
        mod_tas = [x for x in tas_list if model_name(x) == mod]
        
        for cube in mod_tas:
            
            x = cube.interpolate(stat_pts, iris.analysis.Linear())
            
            sim = cube.coord('sim').points[0]
            mod =  tp.gcm(cube) + "_" + tp.model(cube)
            
            y = x.coord('year').points
            m = x.coord('month').points
            d = x.coord('day').points
            
            time = []
            for i in np.arange(len(y)):
                t = str(d[i]) + '-' + str(m[i]) + '-' + str(y[i]) 
                time.append(t)
            
            #calc temp dif (min 0)
            td = x.data - thres
            td = np.where(td < 0, 0, td)
            
            
            mod_df = pd.DataFrame({'Date': time,
                                   'Year': y,
                                   'Month': x.coord('month_num').points,
                                   'Day': d,
                                   'Region': stat_name, 
                                   'Scen': sim,
                                   'Model': mod,
                                   'Tmean': x.data,
                                   'Threshold': thres,
                                   'Temp_diff': td})
            df = df.append(mod_df)
            
    
    

    return df

#%%

save_path = '/nfs/a321/earsch/CHAMNHA/station_output/'

for i in np.arange(0, len(stat_df)):
    stat_name = stat_df['Name'][i].replace(' ', '_')
    
    save_name = stat_name +  '.csv'

    
    print(i, stat_name)
    
    stat_pts = [('latitude', stat_df['Latitude'][i]), ('longitude', stat_df['Longitude'][i])]
    thres = stat_df['Thres'][i]
    
    df = cube_to_frame(tas_list, stat_name, stat_pts, thres)    
    df.to_csv(save_path + save_name, index = False)

    #get daily ensemble mean
    df_agg = df.groupby(['Date', 'Year', 'Month', 'Day', 'Scen', 'Region'], as_index = False).mean()
    df_agg['Temp_diff'] = df_agg['Tmean'] - thres
    df_agg['Threshold'] = thres
    df_agg['Temp_diff'][df_agg['Temp_diff'] < 0] = 0
    
    
    df_agg.to_csv(save_path + 'ENSMEAN_' + save_name, index = False)



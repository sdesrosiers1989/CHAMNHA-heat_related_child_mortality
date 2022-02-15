#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Dcomposition of infant mortality model components

Das Gupta 1993 method

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

def add_code_from_filename(cube, field, filename):
    #split filename into sections by '_', find second partition and split again
    file = filename[97:-3]
    
    code = file.partition('_')[2].partition('_')[2].partition('_')[2]

    cube.add_aux_coord(iris.coords.AuxCoord(code, long_name = 'code'))

def ens_mean(cube_list):  
    cor_cubes = iris.cube.CubeList()
    #add ensemble as dimension so can use cube statistics
    for item in cube_list:
        gcm = item.coord('gcm').points[0]
        try:
            item.add_aux_coord(iris.coords.AuxCoord(gcm, long_name = 'ens'))
        except:
            print('Has ens coord.')
        try:
            item.remove_coord('month')
        except:
            print('No month coord')
        try:
            item.remove_coord('time')
        except:
            print('No time coord')
        try:
            item.remove_coord('year')
        except:
            print('No month coord')
        cor_cubes.append(item)
    ens_cube = cor_cubes.merge_cube()
    ens_cube = ens_cube.collapsed('ens', iris.analysis.MEAN)
    return ens_cube

def find_dif(his_list, histnat_list):
    #his - histnat
    dif_list = []
    difper_list = []
    
    for i in np.arange(len(his_list)):
        cube_list = his_list[i]
        damip_cube_list = histnat_list[i]
        
        new_list = iris.cube.CubeList()
        new_listper = iris.cube.CubeList()
        for cube in cube_list:
            mod = tp.gcm(cube)
            comp_cube =  [x for x in damip_cube_list if tp.gcm(x) == mod][0] # checked, len always 1, as it should be
            print(mod, tp.gcm(comp_cube))
            
            dif_cube = copy.deepcopy(cube)
            dif_cube.data = dif_cube.data - comp_cube.data
            new_list.append(dif_cube)
            
            difper_cube = copy.deepcopy(cube)
            difper_cube.data = (difper_cube.data - comp_cube.data) / difper_cube.data
            new_listper.append(difper_cube)

        dif_list.append(new_list)
        difper_list.append(new_listper)
    
    return dif_list, difper_list

#%%
    
def calc_comp(cube_list):
    '''order should be:
        'e2p1m1', 'e2p2m1', 'e2p1m2',
        'e1p1m1', 'e1p2m1', 'e1p1m2'
        'e2p2m2', 'e1p2m2'
    '''
    
    p1m1e1 = cube_list[3]
    p1m2e2 = cube_list[2]
    p1m1e2 = cube_list[0]
    p1m2e1 = cube_list[5]
    p2m1e1 = cube_list[4]
    p2m2e2 = cube_list[6]
    p2m1e2 = cube_list[1]
    p2m2e1 = cube_list[7]
    
    peq1 = (p1m1e1 + p1m2e2) / 3
    peq2 = (p1m1e2 + p1m2e1) / 6
    peq3 = (p2m1e1 + p2m2e2) / 3
    peq4 = (p2m1e2 + p2m2e1) / 6
    
    p_effect = (peq1 + peq2) - (peq3 + peq4)
    
    meq1 = (p1m1e1 + p2m1e2) / 3
    meq2 = (p1m1e2 + p2m1e1) / 6
    meq3 = (p1m2e1 + p2m2e2) / 3
    meq4 = (p1m2e2 + p2m2e1) / 6
    
    m_effect = (meq1 + meq2) - (meq3 + meq4)
    
    eeq1 = (p1m1e1 + p2m2e1) / 3
    eeq2 = (p1m2e1 + p2m1e1) / 6
    eeq3 = (p1m1e2 + p2m2e2) / 3
    eeq4 = (p1m2e2 + p2m1e2) / 6
    
    e_effect = (eeq1 + eeq2) - (eeq3 + eeq4)
    
    return p_effect, m_effect, e_effect



#%% Import outputs of ifnant mortality health burden model
    
codes = ['e2p1m1', 'e2p2m1', 'e2p1m2',
        'e1p1m1', 'e1p2m1', 'e1p1m2',
        'e2p2m2', 'e1p2m2']

path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/decomp_attribution/'

filenames = glob.glob(path + '*.nc')
cube_list = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file, callback = add_code_from_filename)
    cube_list.append(x)


#'e2p2m2', 'e1p2m2'
    #hist-nat = e1
    #hist = e2
path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/'
filenames = glob.glob(path + '*hist-nat*.nc')
damip_mods = []
for file in filenames:
    x = iris.load_cube(file)
    x.add_aux_coord(iris.coords.AuxCoord('e1p2m2', long_name = 'code'))
    damip_mods.append(tp.gcm(x))
    cube_list.append(x)


filenames = glob.glob(path + '*historical*.nc')
for file in filenames:
    x = iris.load_cube(file)
    x.add_aux_coord(iris.coords.AuxCoord('e2p2m2', long_name = 'code'))
    if tp.gcm(x) in damip_mods:
        cube_list.append(x)


code_cubes = []
for code in codes:
    new_list = [x for x in cube_list if x.coord('code').points == code]
    code_cubes.append(new_list)
    
    
#%% Split into years
    
years = np.unique([x.coord('year').points[0] for x in cube_list]) 

#code_years = []
#for cube_list in code_cubes:
#    new_list = iris.cube.CubeList()
#    for y in years:
#        newy_list = [x for x in cube_list if x.coord('year').points[0] == y]
#    
#        new_list.append(newy_list)
#    code_years.append(new_list)
    
code_years = iris.cube.CubeList()
for y in years:
    new_list = iris.cube.CubeList()
    for cube_list in code_cubes:
        newy_list = [x for x in cube_list if x.coord('year').points[0] == y]

        new_list.append(newy_list)
    code_years.append(new_list)
    
#%% What is the cause of the differences between hist and hist-nat?
    
#1995 - 2004: pop and mortality the same (onyl due to climate)
period1 = code_years[0]
hist_period1 = period1[6]
histnat_period1 = period1[7]


#%% Find difference
    

    
    
#%% Get ens mean

ens_hist = []
ens_damip = []
ens_damip_ht = []

for i in np.arange(len(damip_years)):
    
    x = ens_mean(damip_years[i])
    ens_damip.append(x)
    
    x = ens_mean(damip_ht_years[i])
    ens_damip_ht.append(x)
    
    x = ens_mean(his_years[i])
    ens_hist.append(x)
    

ens_dif = [ens_mean(x) for x in dif_list]
ens_dif_ht = [ens_mean(x) for x in dif_ht_list]

#coeff 1
ens_hist1 = []
ens_damip1 = []
ens_damip1_ht = []

for i in np.arange(len(damip_years1)):
    
    x = ens_mean(damip_years1[i])
    ens_damip1.append(x)
    
    x = ens_mean(damip_ht_years1[i])
    ens_damip1_ht.append(x)
    
    x = ens_mean(his_years1[i])
    ens_hist1.append(x)
    

ens_dif1 = [ens_mean(x) for x in dif_list1]
ens_dif1_ht = [ens_mean(x) for x in dif_ht_list1]


#%%
    
    #coeff 0.61
cube_list = [ens_hist[0], ens_hist[1], ens_hist[2], 
             ens_damip_ht[0], ens_damip_ht[1], ens_damip_ht[2],
             ens_dif_ht[0], ens_dif_ht[1], ens_dif_ht[2]]


# coeff 1
#cube_list = [ens_hist1[0], ens_hist1[1], ens_hist1[2], 
#             ens_damip1[0], ens_damip1[1], ens_damip1[2],
#             ens_dif1[0], ens_dif1[1], ens_dif1[2]]


rows = 3
cols = 3
nfigs = rows * cols

no_grids = np.arange(1,nfigs - cols)[np.mod(np.arange(1, nfigs - cols),cols)!= 0]
y_only = np.arange(0, nfigs - cols, cols)

levels = np.arange(0, 55, 5)
levs_dif = np.arange(-5, 6, 1)
#levels_dif = np.arange(-200, 250, 50)
xticks = np.arange(-10, 60, 20)
yticks = [40, 20, 0, -20, -40]

fig = plt.figure(figsize=(9,10))

ax_list = []
for i in np.arange(1, nfigs+1):
    ax_list.append(plt.subplot(rows, cols, i, projection = proj))

for i in np.arange(0, len(cube_list)):
    ax = ax_list[i]
    cube = cube_list[i]
    
    if i < 6:
        im = iplt.contourf(cube, axes = ax, cmap = 'YlOrRd', levels = levels, extend = 'max')
    else:
        im_dif = iplt.contourf(cube, axes = ax, cmap = 'RdBu_r', levels = levs_dif, extend = 'both')
    
    if i in no_grids:
        tp.plot_africa(ax_list[i], high = False, no_x = True, no_y = True, Tanga = False)
    elif i in y_only:
        tp.plot_africa(ax_list[i], high = False, no_x = True, yticks = yticks, Tanga = False)
    elif i == no_grids[-1] + 1:
        tp.plot_africa(ax_list[i], high = False, xticks = xticks, yticks = yticks, Tanga = False)
    else:
        tp.plot_africa(ax_list[i], high = False, no_y = True, xticks = xticks, Tanga = False)

       
    #ax.set_title(title)
    #tp.plot_africa(ax, high = False, xticks = xticks)

xpos = 0.5
ypos = 1.05
ax_list[0].annotate(s = '1995 - 2004', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[1].annotate(s = '2005 - 2014', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')
ax_list[2].annotate(s = '2015 - 2020', xy = (xpos, ypos), xycoords = 'axes fraction',
             ha = 'center')


xpos = -0.35
ypos = 0.5
ax_list[0].annotate(s = 'Historical', xy = (xpos, ypos), xycoords = 'axes fraction',
             va = 'center', rotation = 90)
ax_list[3].annotate(s = 'Hist-nat', xy = (xpos, ypos), xycoords = 'axes fraction',
             va = 'center', rotation = 90)
ax_list[6].annotate(s = 'Historical minus his-nat', xy = (xpos, ypos), xycoords = 'axes fraction',
             va = 'center', rotation = 90)


plt.draw()

dif = 0.05

cbax = tp.get_cbax(fig, ax = ax_list[3], last_ax = [ax_list[5]], orientation = 'horizontal', 
                   dif = dif/2, h_w = 0.02)
cbar = plt.colorbar(im, cax = cbax, orientation = 'horizontal')
cbar.set_label('Annual avg. heat related deaths for 0 - 4 year olds ')

cbax = tp.get_cbax(fig, ax = ax_list[6], last_ax = [ax_list[-1]], orientation = 'horizontal', 
                   dif = dif, h_w = 0.02)
cbar = plt.colorbar(im_dif, cax = cbax, orientation = 'horizontal')
cbar.set_label('Difference in annual avg. heat related deaths for 0 - 4 year olds ')

#
#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/attributionresulst_ensmean_difens_coeff061_threshismodel.png',
#         bbox_inches = 'tight', pad_inches = 0.3)



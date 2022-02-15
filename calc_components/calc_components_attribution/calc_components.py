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
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

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

def ens_mean(cube_list):  
    cor_cubes = iris.cube.CubeList()
    #add ensemble as dimension so can use cube statistics
    for i in np.arange(len(cube_list)):
        item = cube_list[i]
        gcm = tp.gcm(item)
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
    
    for cube in his_list:
        mod = tp.gcm(cube)
        comp_cube =  [x for x in histnat_list if tp.gcm(x) == mod][0] # checked, len always 1, as it should be
        print(mod, tp.gcm(comp_cube))
        
        dif_cube = copy.deepcopy(cube)
        dif_cube.data = dif_cube.data - comp_cube.data
        dif_list.append(dif_cube)
        
        difper_cube = copy.deepcopy(cube)
        difper_cube.data = (difper_cube.data - comp_cube.data) / difper_cube.data
        difper_list.append(difper_cube)
    
    return dif_list, difper_list

#%%
    
def calc_comp(p_list, m_list, e_list, return_list =  False):
    '''order should be:
    period1, period2
    list of models within each period only for e (pop and mor same each model)
    '''
    
    p1 = p_list[0].data
    p2 = p_list[1].data
    m1 = m_list[0].data
    m2 = m_list[1].data 
    e1_mlist = e_list[0]
    e2_mlist = e_list[1]
    
    p_outlist = iris.cube.CubeList()
    m_outlist = iris.cube.CubeList()
    e_outlist = iris.cube.CubeList()
    
    for i in np.arange(len(e1_mlist)):
        
        e1 = e1_mlist[i].data
        
        #find same gcm
        gcm = tp.gcm(e1_mlist[i])
        e2 = [x for x in e2_mlist if tp.gcm(x) == gcm][0].data
        
        
        peq1 = (m1 * e1 + m2 * e2)  / 3
        peq2 = (m1 * e2 + m2 * e1) / 6
        
        p_effect = copy.deepcopy(e1_mlist[i])
        p_effect.data = (p1 - p2)*(peq1 + peq2)
        
        meq1 = (p1 * e1 + p2 * e2)  / 3
        meq2 = (p1 * e2 + p2 * e1) / 6
                
        m_effect = copy.deepcopy(e1_mlist[i])
        m_effect.data = (m1 - m2)*(meq1 + meq2)
        
        eeq1 = (m1 * p1 + m2 * p2)  / 3
        eeq2 = (p1 * m2 + p2 * m1) / 6
                
        e_effect = copy.deepcopy(e1_mlist[i])
        e_effect.data = (e1 - e2)*(eeq1 + eeq2)
        
        p_outlist.append(p_effect)
        m_outlist.append(m_effect)
        e_outlist.append(e_effect)
        
    p_ens = ens_mean(p_outlist)
    m_ens = ens_mean(m_outlist)
    e_ens = ens_mean(e_outlist)
        
    if return_list == True:
        return p_ens, m_ens, e_ens, p_outlist, m_outlist, e_outlist
    else:
        return p_ens, m_ens, e_ens
    #return p_ens, m_ens, e_ens, p_list, m_list, e_list

#%% import e model data
    
path = '/nfs/a321/earsch/CHAMNHA/output/e/coeff_061/historical/'

#hist-nat
damip_mods = []
filenames = glob.glob(path + '*hist-nat*.nc')
damip_list = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    damip_mods.append(tp.gcm(x))
    damip_list.append(x)


#historical
filenames = glob.glob(path + '*historical*.nc')
hist_list = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    if tp.gcm(x) in damip_mods:
        hist_list.append(x)
  
#%% Import pop and mor data 

#population -from world pop
# 0 - 4 year olds
pop_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000_regrid.nc')
pop_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2010_regrid.nc')
pop_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2019_regrid.nc')

#daily mortality
dmor_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2000_regrid.nc')
dmor_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2010_regrid.nc')
dmor_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/daily_mor_mf_01_2019_regrid.nc')


pop_list = [pop_2000, pop_2010, pop_2019]
mor_list = [dmor_2000, dmor_2010, dmor_2019]

#%% actual mortality

path = '/nfs/a321/earsch/CHAMNHA/output/annual_avg_mortality/coeff_061/thres_hismodel/'

filenames = glob.glob(path + '*hist-nat*')
damip_mor = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    damip_mor.append(x)
    
filenames = glob.glob(path + '*historical*')
his_mor = iris.cube.CubeList()
for file in filenames:
    x = iris.load_cube(file)
    if tp.gcm(x) in damip_mods:
       his_mor.append(x)

#%% calc pop ratio as sued in health burden model
#pop data

pop_ratio_list = []
for i in np.arange(len(pop_list)):

    pop_ratio = pop_list[i]/pop_2000 #ratio future pop to baseline pop   
    #0 in denominator causing problems
    pop_ratio.data = np.where(pop_2000.data == 0, 1, pop_ratio.data)
    pop_ratio.data = ma.masked_array(pop_ratio.data, mask = pop_2000.data.mask)  
    pop_ratio_list.append(pop_ratio)
    
#%% get mdoels into year lists
    
years = np.unique([x.coord('year').points[0] for x in hist_list]) #damip and his years same

damip_years =  iris.cube.CubeList()
his_years = iris.cube.CubeList()

for y in years:
    d_list = [x for x in damip_list if x.coord('year').points[0] == y]
    d_mean_list = iris.cube.CubeList()
    for cube in d_list:
        d_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))

    h_list = [x for x in hist_list if x.coord('year').points[0] == y]
    h_mean_list = iris.cube.CubeList()
    for cube in h_list:
        h_mean_list.append(cube.collapsed('time', iris.analysis.MEAN))
    
    damip_years.append(d_mean_list)
    his_years.append(h_mean_list)
    
#%% mortality years

years = np.unique([x.coord('year').points[0] for x in damip_mor]) #damip and his years same

damip_mor_years =  iris.cube.CubeList()
his_mor_years = iris.cube.CubeList()

for y in years:
    d_list = [x for x in damip_mor if x.coord('year').points[0] == y]
    damip_mor_years.append(d_list)
    
    h_list = [x for x in his_mor if x.coord('year').points[0] == y]
    his_mor_years.append(h_list)


#%% Cause of difference between Histnat T1 and T2?
    
p_inlist = [pop_ratio_list[0], pop_ratio_list[1]]
m_inlist = [dmor_2000, dmor_2010]
e_inlist = [damip_years[0], damip_years[1]]
e_his_inlist = [his_years[0], his_years[1]]

p_effect, m_effect, e_effect, plist, mlist, elist = calc_comp(p_inlist, m_inlist, e_inlist, return_list = True)
#his T2 vs T1
p_effect_h, m_effect_h, e_effect_h, plist_h, mlist_h, elist_h = calc_comp(p_inlist, m_inlist, e_his_inlist, return_list = True)
  
#hist nat T3 vs T1
p_inlist = [pop_ratio_list[0], pop_ratio_list[2]]
m_inlist = [dmor_2000, dmor_2019]
e_inlist = [damip_years[0], damip_years[2]]

p_effect3, m_effect3, e_effect3, plist3, mlist3, elist3 = calc_comp(p_inlist, m_inlist, e_inlist, return_list = True)
  
#hist nat T3 vs T1
p_inlist = [pop_ratio_list[0], pop_ratio_list[2]]
m_inlist = [dmor_2000, dmor_2019]
e_inlist = [damip_years[0], damip_years[2]]
e_his_inlist = [his_years[0], his_years[2]]

p_effect3, m_effect3, e_effect3, plist3, mlist3, elist3 = calc_comp(p_inlist, m_inlist, e_inlist, return_list = True)
p_effect_h3, m_effect_h3, e_effect_h3, plist_h3, mlist_h3, elist_h3 = calc_comp(p_inlist, m_inlist, e_his_inlist, return_list = True)

  
#%% Check against actual dif

#1995 - 2004: pop and mortality the same (onyl due to climate)
period1 = damip_mor_years[0]
period2 = damip_mor_years[1]
period3 = damip_mor_years[2]

dif, dif_per = find_dif(period1, period2)

ens_dif = ens_mean(dif)
ens_perdif = ens_mean(dif_per)

period1_ens = ens_mean(period1)
period2_ens = ens_mean(period2)
period3_ens = ens_mean(period3)


#hist ant T3 s T1 dif
dif3, dif_per3 = find_dif(period1, period3)

ens_dif3 = ens_mean(dif3)
ens_perdif3 = ens_mean(dif_per3)

#historical difs

period1_h = his_mor_years[0]
period2_h = his_mor_years[1]
period3_h = his_mor_years[2]

dif_h, dif_per_h = find_dif(period1_h, period2_h)

ens_dif_h = ens_mean(dif_h)
ens_perdif_h = ens_mean(dif_per_h)

period1_ens_h = ens_mean(period1_h)
period2_ens_h = ens_mean(period2_h)
period3_ens_h = ens_mean(period3_h)


#hist ant T3 s T1 dif
dif_h3, dif_per_h3 = find_dif(period1_h, period3_h)

ens_dif_h3 = ens_mean(dif_h3)
ens_perdif_h3 = ens_mean(dif_per_h3)


#%% Check decomp

print(np.nanmean(period1_ens.data))
print(np.nanmean(period2_ens.data))
print(np.nanmean(ens_dif.data))

np.nanmean(p_effect.data)
np.nanmean(m_effect.data)
np.nanmean(e_effect.data)

np.nanmean(p_effect.data) +  np.nanmean(m_effect.data) +   np.nanmean(e_effect.data) 


np.nanmean(p_effect3.data)
np.nanmean(m_effect3.data)
np.nanmean(e_effect3.data)
print(np.nanmean(ens_dif3.data))
np.nanmean(p_effect3.data) +  np.nanmean(m_effect3.data) +   np.nanmean(e_effect3.data) 

#historical
print(np.nanmean(ens_dif_h3.data))
np.nanmean(p_effect_h3.data) +  np.nanmean(m_effect_h3.data) +   np.nanmean(e_effect_h3.data) 

print(np.nanmean(ens_dif_h.data))
np.nanmean(p_effect_h.data) +  np.nanmean(m_effect_h.data) +   np.nanmean(e_effect_h.data) 

#%%check individual historical models - order mixed up

print(np.nanmean(dif_h3[2].data))
np.nanmean(plist_h3[0].data) +  np.nanmean(mlist_h3[0].data) +   np.nanmean(elist_h3[0].data) 


#%% What percentage of the change is due to climate, pop and mortality?

def cube_to_frame(p_effect_list, m_effect_list, e_effect_list, total_dif, period):
    
    df = pd.DataFrame(columns = ['model', 'period', 'p', 'm', 'e', 'p_per', 'm_per', 'e_per', 'total'])

    for i in np.arange(len(p_effect_list)):
        p = p_effect_list[i]
        gcm = tp.gcm(p)
        
        m = [x for x in m_effect_list if x.coord('ens').points[0] == gcm][0]
        e = [x for x in e_effect_list if x.coord('ens').points[0] == gcm][0]
        d = [x for x in total_dif if x.coord('ens').points[0] == gcm][0]
        
        p_val =  np.nanmean(p.data)
        m_val =  np.nanmean(m.data)
        e_val =  np.nanmean(e.data)
        d_val =  np.nanmean(d.data)
        
        p_per = p_val / d_val * 100
        m_per = m_val / d_val * 100
        e_per = e_val / d_val * 100
    
            
        y = pd.DataFrame(data = {'model': gcm,
                                 'period': period,
                                 'p': p_val,
                                 'm': m_val,
                                 'e': e_val,
                                 'p_per': p_per,
                                 'm_per': m_per,
                                 'e_per': e_per,
                                 'total': p_per + m_per + e_per}, index = [0])

        
        df = df.append(y)
    
    return df


damip_p2 = cube_to_frame(plist, mlist, elist, dif, 'period2')
damip_p3 = cube_to_frame(plist3, mlist3, elist3, dif3, 'period3')
his_p2 = cube_to_frame(plist_h, mlist_h, elist_h, dif_h, 'period2')
his_p3 = cube_to_frame(plist_h3, mlist_h3, elist_h3, dif_h3, 'period3')

damip = damip_p2.append(damip_p3)
his = his_p2.append(his_p3)


damip_agg = damip.groupby(['period'], as_index = True).mean()
his_agg = his.groupby(['period'], as_index = True).mean()

d_per = damip_agg.drop(columns = ['p' , 'm', 'e', 'total'])
h_per = his_agg.drop(columns = ['p' , 'm', 'e', 'total'])

d_abs = damip_agg.drop(columns = ['p_per' , 'm_per', 'e_per', 'total'])
h_abs = his_agg.drop(columns = ['p_per' , 'm_per', 'e_per', 'total'])

#otherwise looks like population growth makes mortlaity go down
d_abs = d_abs * -1
h_abs = h_abs * -1

#%% https://stackoverflow.com/questions/22787209/how-to-have-clusters-of-stacked-bars-with-python-pandas

def plot_clustered_stacked(dfall, cols, labels=None, title="",  H="/", width = 0.5,  **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe"""

    n_df = len(dfall)
    n_col = len(dfall[0].columns) 
    n_ind = len(dfall[0].index)
    axe = plt.subplot(111)

    for df in dfall : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ax=axe,
                      legend=False,
                      grid=False,
                      color = cols,
                      **kwargs)  # make bar plots

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + width / float(n_df + width) * i / float(n_col))
                rect.set_hatch(H * int(i / n_col)) #edited part     
                rect.set_width(width / float(n_df + width))

    axe.set_xticks(np.arange(n_ind))
    axe.set_xticklabels(df.index, rotation = 0)
    #axe.set_title(title)

    # Add invisible data to add another legend
    #n=[]        
    #for i in range(n_df):
    #    n.append(axe.bar(0, 0, color="gray", hatch=H * i))

    #l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0.5])
    #if labels is not None:
     #   l2 = plt.legend(n, labels, loc=[1.01, 0.1]) 
    #axe.add_artist(l1)
    return axe

# create fake dataframes
df1 = pd.DataFrame(np.random.rand(4, 5),
                   index=["A", "B", "C", "D"],
                   columns=["I", "J", "K", "L", "M"])
df2 = pd.DataFrame(np.random.rand(4, 5),
                   index=["A", "B", "C", "D"],
                   columns=["I", "J", "K", "L", "M"])
df3 = pd.DataFrame(np.random.rand(4, 5),
                   index=["A", "B", "C", "D"], 
                   columns=["I", "J", "K", "L", "M"])

# Then, just call :
plot_clustered_stacked([df1, df2],labels = ["df1", "df2"], cols = ['black', 'green', 'blue'])

#%% plot

labels = ['Period 2', 'Period 3', '6', '7']
mor_means = [his_agg['m_per'], damip_agg['m_per'], his_agg3['m_per'], damip_agg3['m_per']]
pop_means = [his_agg['p_per'], damip_agg['p_per'], his_agg3['p_per'], damip_agg3['p_per']]
e_means = [his_agg['e_per'], damip_agg['e_per'], his_agg3['e_per'], damip_agg3['e_per']]
width  = 0.35

plt.bar(labels, mor_means, width, label = 'Mortality')
plt.bar(labels, pop_means, width, label = 'Population')
plt.bar(labels, e_means, width, label = 'Climate')

plt.legend()

#%%

plot_clustered_stacked([h_per, d_per],["his", "damip"])

#%%
ax = plot_clustered_stacked([h_abs, d_abs],labels = ["his", "damip"], cols = ['#1b9e77', '#7570b3', '#d95f02'])

plt.xlabel('Period')

ax.set_ylabel('Mean annual heat mortality per gridcell')
ax.set_xticklabels(['2005 - 2014', '2015 - 2020'])

ax.axhline(0, c = 'k', linewidth = 0.5)



mor_patch = Patch(facecolor='#7570b3')
pop_patch = Patch(facecolor='#1b9e77')
e_patch = Patch(facecolor='#d95f02')
his_patch = Patch(facecolor='grey')
damip_patch = Patch(facecolor='grey', hatch = '//')

custom_lines = [pop_patch, his_patch, mor_patch, damip_patch, e_patch]


ax.legend(custom_lines, ['Population', 'Historical', 'Mortality', 
                         'Hist-nat', 'Climate'], 
           bbox_to_anchor=(0.9, -0.2),ncol=3,frameon = False, handletextpad = 0.5)




#%%


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(1,1,1)

barWidth = 0.1

# Set position of bar on X axis (position for each rcm bar within each gcm category)
positions = []
r1 = np.arange(2)
positions.append(r1)

periods = ['Period2', 'Period3']
scen_dat = [h_abs, d_abs]

for i in np.arange(1,len(periods)):
    next_pos = [x + barWidth for x in positions[i-1]]
    positions.append(next_pos)
 
#set colors
cols = ['green', 'blue', 'red']
var_name = ['m', 'p', 'e']
var_label = ['Mortality', 'Population', 'Climate']
scen_label = ['hst', 'hist-nat']

lw = 2.0


# Make the plot
for i in np.arange(len(scen_dat)):
    pos = positions[i]
    dat = scen_dat[i]
    label = scen_label
    plt.bar(pos, [5, 10], width = barWidth, color = cols[i], label = label)
    plt.bar(pos, [5, 10], width = barWidth, color = cols[2], label = label)
    
    
    #%%
#ax.set_ylim([-70, 70])
#ax.axhline(y = 0, color = 'k', linewidth = 0.5)
#ax.annotate(s = var_name[j], xy = (0.01, 1.02), xycoords = 'axes fraction')
# Add xticks on the middle of the group bars
ax.set_xticks(positions[5]) #positions 5 is the middle rcm
ax.set_xticklabels(np.repeat('', 11))

#add lines between groups
y = [25, 0]
if j == 2: 
    y = [0.025, 0]
elif j == 1:
    y = [15, 0]
    
adjust = 0.05
for i in np.arange(0, len(positions[0])):
    #position[0] is position of the first gcm
    #position[-1] has the position of the last gcm
    x1 = positions[0][i]
    x = [x1 - adjust, x1 - adjust]     
    ax.plot(x, y, color = 'black', linewidth = 0.5)


if j == 2:
    ax.set_xlabel('RCM')
    ax.set_xticklabels(rcm_short)
ax.tick_params(bottom = 'False')

ax.set_ylabel(var_label[j])
    
    # add CP4 and P25
    ax.scatter(7.3, df_cp4[var_name[j].lower()][(df_cp4['scen'] == 'hist')].values, c = 'black',
                            label = 'CP4A')
    ax.scatter(7.3, df_p25[var_name[j].lower()][(df_p25['scen'] == 'hist')].values, c = 'black',
                            marker = 's', label = 'P25')
    #add trmm and obs
    
    ax.axhline(df15_agg_obs[var_name[j].lower()][(df15_agg_obs['mod'] == 'obs_trmm')].values,
                            linestyle = 'dotted', color = '#e52929', label = 'TRMM', linewidth = lw)
                            
    ax.axhline(df38_agg_obs[var_name[j].lower()].values,
                             color = '#e52929', label = 'CHIRPS', linewidth = lw)
    

plt.draw()
#plt.annotate('Change in rainfall (mm/month)', (0.04, 0.65), xycoords = 'figure fraction',
#             rotation = 90)

# Create legend & Show graphic
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25),ncol=4, frameon = False)

#fig.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Barchar_params_30_rcms_obs.png',
#            bbox_inches = 'tight', pad_inches = 0.3)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process SSP future data for use

-in csv format
-link to coutnry raster(dictionary)
-distribution using 2019 pop distribution

Created on Tue Aug 17 11:44:49 2021

@author: earsch
"""
#%% import libraris

import iris
import iris.quickplot as qplt

import matplotlib.pyplot as plt

import glob

import numpy as np

import pandas as pd

import copy
#%%

countries = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/africa_countryname.nc')
popfrac_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/popfrac_2019_regrid.nc')
pop2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2019_regrid.nc')


#change lookup to match country names in GBD data
country_lookup = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/pop/country_lookup.csv')
country_lookup['Name'][country_lookup['Name'] == 'Republic of the Congo'] = 'Congo'

mor_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2019_regrid.nc')

mor_ref = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_reference_both_04.csv')
mor_better = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_better_both_04.csv')
mor_worse = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_worse_both_04.csv')


#%% Create diciotnary of country codes and names

c_dict = dict(zip(country_lookup['Value'], country_lookup['Name']))


#%% regrid countries to poprfac

countries_regrid = countries.regrid(popfrac_2019, iris.analysis.Nearest())


#%% erplace country code with total under 5 population

save_path = '/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/processed/'

def distr_dat(in_table, scenario_name):


    years = np.unique(in_table['Year'])
    years = years[years >= 2020]
    
    output = []
    
    for y in years:
    
        year_dat = in_table[in_table['Year'] == y]
        cnames = np.unique(year_dat['Location'])
        
        mor_year = copy.deepcopy(countries_regrid)
        
        # get all country codes from map
        dat = mor_year.data
        vals = np.unique(dat)
        vals = vals[vals.mask == False]
        
        for val in vals:
            c_name = c_dict[val]
            print(c_name)
            
            if c_name in cnames:
        
                replace_val = year_dat[year_dat['Location'] == c_name]['Value'].values
                print(replace_val)
                
                mor_year.data = np.where(mor_year.data == val, replace_val, mor_year.data)
                
                mor_year_frac = copy.deepcopy(mor_year)
                mor_year_frac.data = mor_year_frac.data*popfrac_2019.data
                
        save_name = scenario_name + '_' + str(y) + '_04_totalmor_mf.nc'
        iris.save(mor_year_frac, save_path + save_name)
        
        output.append(mor_year_frac)
    return output

mor_ref_output = distr_dat(mor_ref, 'ref')
mor_better_output = distr_dat(mor_better, 'better')
mor_worse_output = distr_dat(mor_worse, 'worse')


#%%  check
        
pop_levs = [0, 3000, 20000, 30000, 50000, 100000, 250000, 500000] 
mor_levs = np.array([0, 1, 4, 6, 10, 15, 30]) * 365


qplt.contourf(mor_2019, levels = mor_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(mor_ref_output[0], levels = mor_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(mor_better_output[0], levels = mor_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(mor_worse_output[0], levels = mor_levs, extend = 'max', cmap = 'YlOrRd')

#%%
print(np.nanmean(mor_2019.data))
mmor = [np.nanmean(x.data) for x in mor_ref_output]
mmorb = [np.nanmean(x.data) for x in mor_better_output]
mmorw = [np.nanmean(x.data) for x in mor_worse_output]

mmor[0]
mmorb[0]
mmorw[0]
#mmor3 = [np.nanmean(x.data) for x in pop_output_ssp3]


print(np.nansum(mor_2019.data))
tmor = [np.nansum(x.data) for x in mor_ref_output]
tmorb = [np.nansum(x.data) for x in mor_better_output]
tmorw = [np.nansum(x.data) for x in mor_worse_output]

tmor[0]
tmorb[0]
tmorw[0]

#%%

years = np.unique(mor_ref['Year'])
years = years[years >= 2020]

plt.plot(years, tmor, label = 'Reference')
plt.plot(years, tmorb, label = 'Better')
plt.plot(years, tmorw, label = 'Worse')

plt.scatter(2019, np.nansum(mor_2019.data), label = 'actual 2019')
plt.xlabel('Year')
plt.ylabel('Total under 5 mortality')
plt.legend()

#plt.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/future_mortality.png',
#         bbox_inches = 'tight', pad_inches = 0.3)
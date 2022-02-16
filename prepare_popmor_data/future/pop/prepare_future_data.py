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


country_lookup = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/pop/country_lookup.csv')
ssp2 = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/pop/future/ssp2_allafricacountries_04_bothsexes.csv', 
                   skiprows = 8)

ssp2['Population'] = ssp2['Population'] * 1000
#%% Create diciotnary of country codes and names

c_dict = dict(zip(country_lookup['Value'], country_lookup['Name']))


#%% regrid countries to poprfac

countries_regrid = countries.regrid(popfrac_2019, iris.analysis.Nearest())


#%% erplace country code with total under 5 population

save_path = '/nfs/a321/earsch/CHAMNHA/input_data/pop/future/processed/'

years = np.unique(ssp2['Year'])
years = years[years >= 2020]

pop_output = []

for y in years:

    year_dat = ssp2[ssp2['Year'] == y]
    cnames = np.unique(year_dat['Area'])
    
    pop_year = copy.deepcopy(countries_regrid)
    
    # get all country codes from map
    dat = pop_year.data
    vals = np.unique(dat)
    vals = vals[vals.mask == False]
    
    for val in vals:
        c_name = c_dict[val]
        print(c_name)
        
        if c_name in cnames:
    
            replace_val = year_dat[year_dat['Area'] == c_name]['Population'].values
            print(replace_val)
            
            pop_year.data = np.where(pop_year.data == val, replace_val, pop_year.data)
            
            pop_year_frac = copy.deepcopy(pop_year)
            pop_year_frac.data = pop_year_frac.data*popfrac_2019.data
            
    save_name = 'ssp2_' + str(y) + '_04population_mf.nc'
    #iris.save(pop_year_frac, save_path + save_name)
    
    pop_output.append(pop_year_frac)
    

#%%  check
        
pop_levs = [0, 3000, 20000, 30000, 50000, 100000, 250000, 500000] 

qplt.contourf(pop2019, levels = pop_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(pop_output[4], levels = pop_levs, extend = 'max', cmap = 'YlOrRd')

#%%
print(np.nanmean(pop2019.data))
[print(np.nanmean(x.data)) for x in pop_output]

print(np.nansum(pop2019.data))
tpop = [np.nansum(x.data) for x in pop_output]

#%%

plt.plot(years, tpop, label = 'SSP2')
plt.scatter(2019, np.nansum(pop2019.data), label = 'actual 2019')
plt.xlabel('Year')
plt.ylabel('Total African population (under 5)')
plt.legend()
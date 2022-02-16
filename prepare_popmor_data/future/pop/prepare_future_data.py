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

ssp3 = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/pop/future/ssp3_africa_04_bothsexes.csv', 
                   skiprows = 8)

ssp3['Population'] = ssp3['Population'] * 1000
#%% Create diciotnary of country codes and names

c_dict = dict(zip(country_lookup['Value'], country_lookup['Name']))


#%% regrid countries to poprfac

countries_regrid = countries.regrid(popfrac_2019, iris.analysis.Nearest())


#%% erplace country code with total under 5 population

save_path = '/nfs/a321/earsch/CHAMNHA/input_data/pop/future/processed/'

def distr_pop(pop_table, scenario_name):


    years = np.unique(pop_table['Year'])
    years = years[years >= 2020]
    
    pop_output = []
    
    for y in years:
    
        year_dat = pop_table[pop_table['Year'] == y]
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
                
        save_name = scenario_name + '_' + str(y) + '_04population_mf.nc'
        iris.save(pop_year_frac, save_path + save_name)
        
        pop_output.append(pop_year_frac)
    return pop_output

pop_output = distr_pop(ssp2, 'ssp2')
pop_output_ssp3 = distr_pop(ssp3, 'ssp3')

#%%  check
        
pop_levs = [0, 3000, 20000, 30000, 50000, 100000, 250000, 500000] 

qplt.contourf(pop2019, levels = pop_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(pop_output_ssp3[0], levels = pop_levs, extend = 'max', cmap = 'YlOrRd')

#%%
print(np.nanmean(pop2019.data))
mpop = [np.nanmean(x.data) for x in pop_output]
mpop3 = [np.nanmean(x.data) for x in pop_output_ssp3]


print(np.nansum(pop2019.data))
tpop = [np.nansum(x.data) for x in pop_output]
tpop3 = [np.nansum(x.data) for x in pop_output_ssp3]

#%%

years = np.unique(ssp2['Year'])
years = years[years >= 2020]

plt.plot(years, tpop, label = 'SSP2')
plt.plot(years, tpop3, label = 'SSP3')
plt.scatter(2019, np.nansum(pop2019.data), label = 'actual 2019')
plt.xlabel('Year')
plt.ylabel('Total African population (under 5)')
plt.legend()
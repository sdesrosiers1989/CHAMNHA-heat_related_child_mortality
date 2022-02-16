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
pop2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2010_regrid.nc')
pop2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000_regrid.nc')


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
countries_regrid.coord('latitude').guess_bounds()
countries_regrid.coord('longitude').guess_bounds()


popfrac_2019.coord('latitude').guess_bounds()
popfrac_2019.coord('longitude').guess_bounds()
p2000_regrid = pop2000.regrid(popfrac_2019, iris.analysis.AreaWeighted())
p2010_regrid = pop2010.regrid(popfrac_2019, iris.analysis.AreaWeighted())


print(np.nansum(pop2000.data), np.nansum(p2000_regrid.data))
print(np.nansum(pop2010.data), np.nansum(p2010_regrid.data))


#%% Bias correct -finding dif between 2000, 2010 - avg diff


def country_total(pop, country_raster):
            
    # get all country codes from map
    country = copy.deepcopy(country_raster)
    dat = country.data
    vals = np.unique(dat)
    vals = vals[vals.mask == False]
    
    df = pd.DataFrame(columns = ['Area', 'pop'])
    
    for val in vals:
        val_raster = copy.deepcopy(country)
        val_raster.data = np.where(val_raster.data == val, 1, np.nan)
        
        country_pop = pop * val_raster
        
        cpop_val = np.nansum(country_pop.data)
        cname = c_dict[val]
        
        y = pd.DataFrame({'Area': cname,
                          'pop': cpop_val},
                            index = [0])
        
        df = df.append(y)
        
    return df
    
df_2000 = country_total(p2000_regrid, countries_regrid)
df_2010 = country_total(p2010_regrid, countries_regrid)

df_ssp_2000 = ssp2[ssp2['Year'] == 2000]

both2000 = pd.merge(df_ssp_2000, df_2000, on = 'Area')
both2000['dif'] = both2000['pop'] - both2000['Population']

df_ssp_2010 = ssp2[ssp2['Year'] == 2010]

both2010 = pd.merge(df_ssp_2010, df_2010, on = 'Area')
both2010['dif'] = both2010['pop'] - both2010['Population']

both_years = pd.merge(both2000, both2010, on = ['Area'])
both_years['avg_dif'] = (both_years['dif_x'] + both_years['dif_y']) / 2

#%% erplace country code with total under 5 population

save_path = '/nfs/a321/earsch/CHAMNHA/input_data/pop/future/processed/'

def distr_pop(pop_table, scenario_name, bias_corr):


    years = np.unique(pop_table['Year'])
    years = years[years >= 2000]
    
    pop_output = []
    
    for y in years:
    
        year_dat = pop_table[pop_table['Year'] == y]
        
        #adjust
        year_dat_merge = pd.merge(year_dat, bias_corr, on = 'Area')
        year_dat_merge['Population_corr'] = year_dat_merge['Population'] - year_dat_merge['avg_dif']
        
        
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
        
                replace_val = year_dat_merge[year_dat_merge['Area'] == c_name]['Population_corr'].values
                print(replace_val)
                
                pop_year.data = np.where(pop_year.data == val, replace_val, pop_year.data)
                
                pop_year_frac = copy.deepcopy(pop_year)
                pop_year_frac.data = pop_year_frac.data*popfrac_2019.data
                
        save_name = scenario_name + '_' + str(y) + '_04population_mf_BIASCORR.nc'
        iris.save(pop_year_frac, save_path + save_name)
        
        pop_output.append(pop_year_frac)
    return pop_output

pop_output = distr_pop(ssp2, 'ssp2', both_years)



#%%  check
        
pop_levs = [0, 3000, 20000, 30000, 50000, 100000, 250000, 500000] 

qplt.contourf(pop2019, levels = pop_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(pop_output[4], levels = pop_levs, extend = 'max', cmap = 'YlOrRd')


qplt.contourf(pop2000, levels = pop_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(p2000_regrid, levels = pop_levs, extend = 'max', cmap = 'YlOrRd')

#%%
print(np.nanmean(pop2019.data))
mpop = [np.nanmean(x.data) for x in pop_output]
#mpop3 = [np.nanmean(x.data) for x in pop_output_ssp3]


print(np.nansum(pop2019.data))
tpop = [np.nansum(x.data) for x in pop_output]
#tpop3 = [np.nansum(x.data) for x in pop_output_ssp3]

#%%

years = np.unique(ssp2['Year'])
years = years[years >= 2000]

plt.plot(years, tpop, label = 'SSP2')
#plt.plot(years, tpop3, label = 'SSP3')
plt.scatter(2010, np.nansum(p2010_regrid.data), label = 'actual 2010')
plt.scatter(2000, np.nansum(p2000_regrid.data), label = 'actual 2000')
plt.xlabel('Year')
plt.ylabel('Total African population (under 5)')
plt.legend()
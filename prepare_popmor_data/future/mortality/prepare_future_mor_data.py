#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process future mortality data (GBD)

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
mor_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2010_regrid.nc')
mor_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2000_regrid.nc')

mor_ref = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_reference_both_04.csv')
mor_better = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_better_both_04.csv')
mor_worse = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_worse_both_04.csv')


#%% Create diciotnary of country codes and names

c_dict = dict(zip(country_lookup['Value'], country_lookup['Name']))

#countires in countries not in mor_ref
names = country_lookup['Name']
ref_names = mor_ref['Location']

[x for x in names if x not in ref_names]


#%% regrid countries to poprfac

countries_regrid = countries.regrid(popfrac_2019, iris.analysis.Nearest())
countries_regrid.coord('latitude').guess_bounds()
countries_regrid.coord('longitude').guess_bounds()


popfrac_2019.coord('latitude').guess_bounds()
popfrac_2019.coord('longitude').guess_bounds()
m2000_regrid = mor_2000.regrid(popfrac_2019, iris.analysis.AreaWeighted())
m2010_regrid = mor_2010.regrid(popfrac_2019, iris.analysis.AreaWeighted())
m2019_regrid = mor_2019.regrid(popfrac_2019, iris.analysis.AreaWeighted())

#%%

m2019_regrid.data = np.ma.masked_array(m2019_regrid.data, countries_regrid.data.mask)
m2010_regrid.data = np.ma.masked_array(m2010_regrid.data, countries_regrid.data.mask)
m2000_regrid.data = np.ma.masked_array(m2000_regrid.data, countries_regrid.data.mask)


#%% Bias correct -finding dif between 2000, 2010 - avg diff


def country_total(mor, country_raster):
            
    # get all country codes from map
    country = copy.deepcopy(country_raster)
    dat = country.data
    vals = np.unique(dat)
    vals = vals[vals.mask == False]
    
    df = pd.DataFrame(columns = ['Location', 'mor'])
    
    for val in vals:
        val_raster = copy.deepcopy(country)
        val_raster.data = np.where(val_raster.data == val, 1, np.nan)
        
        country_mor = mor * val_raster
        
        cmor_val = np.nansum(country_mor.data)
        cname = c_dict[val]
        
        y = pd.DataFrame({'Location': cname,
                          'mor': cmor_val},
                            index = [0])
        
        df = df.append(y)
        
    return df
    
df_2000 = country_total(m2000_regrid, countries_regrid)
df_2010 = country_total(m2010_regrid, countries_regrid)
df_2019 = country_total(m2019_regrid, countries_regrid)


df_ref_2000 = mor_ref[mor_ref['Year'] == 2000]

both2000 = pd.merge(df_ref_2000, df_2000, on = 'Location')
both2000['dif'] = both2000['mor'] - both2000['Value']

df_ref_2010 = mor_ref[mor_ref['Year'] == 2010]

both2010 = pd.merge(df_ref_2010, df_2010, on = 'Location')
both2010['dif'] = both2010['mor'] - both2010['Value']

both_years = pd.merge(both2000, both2010, on = ['Location'])
both_years['avg_dif'] = (both_years['dif_x'] + both_years['dif_y']) / 2

cols = both_years.columns
keep = ['Location', 'Year_x', 'Value_x', 'Value_y', 'Year_y', 'mor_x', 'mor_y', 'avg_dif']
cols = [x for x in cols if x not in keep]
both_years.drop(columns = cols, inplace =  True)

#%% erplace country code with total under 5 population

save_path = '/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/processed/'

def distr_dat(in_table, scenario_name, bias_corr):


    years = np.unique(in_table['Year'])
    years = years[years >= 2000]
    
    output = []
    
    for y in years:
    
        year_dat = in_table[in_table['Year'] == y]
        #adjust
        year_dat_merge = pd.merge(year_dat, bias_corr, on = 'Location')
        year_dat_merge['Mor_corr'] = year_dat_merge['Value'] + year_dat_merge['avg_dif']
        
        
        
        cnames = np.unique(year_dat['Location'])
        
        mor_year = copy.deepcopy(countries_regrid)
        
        # get all country codes from map
        dat = mor_year.data
        vals = np.unique(dat)
        vals = vals[vals.mask == False]
        
        for val in vals:
            c_name = c_dict[val]
            #print(c_name)
            
            if c_name in cnames:
        
                replace_val = year_dat_merge[year_dat_merge['Location'] == c_name]['Mor_corr'].values
                #print(replace_val)
                if replace_val < 0:
                    print( y, c_name, replace_val)
                
                mor_year.data = np.where(mor_year.data == val, replace_val, mor_year.data)
                
                mor_year_frac = copy.deepcopy(mor_year)
                mor_year_frac.data = mor_year_frac.data*popfrac_2019.data
                
        save_name = scenario_name + '_' + str(y) + '_04_totalmor_mf_BIASCORR.nc'
        iris.save(mor_year_frac, save_path + save_name)
        
        output.append(mor_year_frac)
    return output

mor_ref_output = distr_dat(mor_ref, 'ref', both_years)
#pop_output_ssp3 = distr_pop(ssp3, 'ssp3')

#%%  check
        
pop_levs = [0, 3000, 20000, 30000, 50000, 100000, 250000, 500000] 
mor_levs = np.array([0, 1, 4, 6, 10, 15, 30]) * 365


qplt.contourf(mor_2019, levels = mor_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(mor_ref_output[19], levels = mor_levs, extend = 'max', cmap = 'YlOrRd')

#%%
print(np.nanmean(mor_2019.data))
mmor = [np.nanmean(x.data) for x in mor_ref_output]
#mmor3 = [np.nanmean(x.data) for x in pop_output_ssp3]


print(np.nansum(mor_2019.data))
tmor = [np.nansum(x.data) for x in mor_ref_output]
#tmor3 = [np.nansum(x.data) for x in pop_output_ssp3]

#%%

years = np.unique(mor_ref['Year'])
years = years[years >= 2000]

plt.plot(years, tmor, label = 'Ref')
plt.scatter(2019, np.nansum(mor_2019.data), label = 'actual 2019')
plt.xlabel('Year')
plt.ylabel('Total African mortality (under 5)')
plt.legend()
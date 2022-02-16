#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process SSP future data for use

-in csv format
-link to coutnry raster(dictionary)
-distribution using 2019 pop distribution

Biascorrect all countries based on specific country correction factor
-need to be regridded and correctly masked prior to correction

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


gpw_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/gpw_both_densityperkm_0_4_2010_regrid.nc')
gpw_grid_area = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/gpw_gridarea.nc')
gpw_2010 = gpw_2010 * gpw_grid_area / 1000000

#change llok up anmes to match ssp
country_lookup = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/pop/country_lookup.csv')
country_lookup['Name'][country_lookup['Name'] == 'Republic of the Congo'] = 'Congo'
country_lookup['Name'][country_lookup['Name'] == 'eSwatini'] = 'Swaziland'
country_lookup['Name'][country_lookup['Name'] == 'Tanzania'] = 'United Republic of Tanzania'
country_lookup['Name'][country_lookup['Name'] == 'The Gambia'] = 'Gambia'
country_lookup['Name'][country_lookup['Name'] == 'Libya'] = 'Libyan Arab Jamahiriya'


ssp2 = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/pop/future/ssp2_allcountries_04_bothsexes.csv', 
                   skiprows = 8)

ssp2['Population'] = ssp2['Population'] * 1000




#regrid to CMIP6
scen = '119'
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/bias_corr/CHAMNHA/tas/end/'
filenames = glob.glob(path + '*' + scen + '.nc')
tas = iris.load_cube(filenames[0])[0]

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
p2019_regrid = pop2019.regrid(popfrac_2019, iris.analysis.AreaWeighted())


print(np.nansum(pop2000.data), np.nansum(p2000_regrid.data))
print(np.nansum(pop2010.data), np.nansum(p2010_regrid.data))

#%%apply countries mask to pop data

p2019_regrid.data = np.ma.masked_array(p2019_regrid.data, countries_regrid.data.mask)
p2010_regrid.data = np.ma.masked_array(p2010_regrid.data, countries_regrid.data.mask)
p2000_regrid.data = np.ma.masked_array(p2000_regrid.data, countries_regrid.data.mask)

gpw_2010.data = np.ma.masked_array(gpw_2010.data, countries_regrid.data.mask)

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
df_2019 = country_total(p2019_regrid, countries_regrid)

df_g2010 = country_total(gpw_2010, countries_regrid)

df_ssp_2000 = ssp2[ssp2['Year'] == 2000]

both2000 = pd.merge(df_ssp_2000, df_2000, on = 'Area')
both2000['dif'] = both2000['pop'] - both2000['Population']

df_ssp_2010 = ssp2[ssp2['Year'] == 2010]

both2010 = pd.merge(df_ssp_2010, df_2010, on = 'Area')
both2010['dif'] = both2010['pop'] - both2010['Population']

both_years = pd.merge(both2000, both2010, on = ['Area'])
both_years['avg_dif'] = (both_years['dif_x'] + both_years['dif_y']) / 2

x = country_total(popfrac_2019, countries_regrid)

#%%

c = df_ssp_2000['Area']
c2 = df_2000['Area'].values

[x for x in c if x not in c2]



#%%

path = '/nfs/a321/earsch/CHAMNHA/input_data/pop/future/processed/'
file_names = glob.glob(path + 'ssp2*BIASCORR2.nc')
pop_output = iris.load(file_names)

file_names = glob.glob(path + 'ssp2*mf.nc')
pop_output_raw = iris.load(file_names)

#%% restrict to afric only

not_in_afr = ['Albania', 'Armenia', 'Azerbaijan', 'Cyprus', 'France', 'Greece', 
              'Iran', 'Iraq', 'Israel', 'Italy', 'Jordan', 'Kuwait', 'Lebanon',
              'Malta', 'Northern Cyprus', 'Oman', 'Palestine', 'Portugal', 'Qatar',
              'Saudi Arabia', 'Spain', 'Syria', 'Turkey', 'Turkmenistan', 'United Arab Emirates', 'Yemen']



#%%

tpop = []
tpop_raw =[]

cregrid = countries_regrid.regrid(pop_output[0], iris.analysis.Nearest())
cregrid_raw = countries_regrid.regrid(pop_output_raw[0], iris.analysis.Nearest())

for i in np.arange(len(pop_output)):
    df = country_total(pop_output[i], cregrid)
    df = df[~df['Area'].isin(not_in_afr)]
    
    x = np.nansum(df['pop'])
    tpop.append(x)
    
    df = country_total(pop_output_raw[i], cregrid_raw)
    df = df[~df['Area'].isin(not_in_afr)]
    
    x = np.nansum(df['pop'])
    tpop_raw.append(x)


#actual 200 and 2010
df = country_total(p2000_regrid, countries_regrid)
df = df[~df['Area'].isin(not_in_afr)]
x2000 = np.nansum(df['pop'])

df = country_total(p2010_regrid, countries_regrid)
df = df[~df['Area'].isin(not_in_afr)]
x2010 = np.nansum(df['pop'])

df = country_total(gpw_2010, countries_regrid)
df = df[~df['Area'].isin(not_in_afr)]
gx2010 = np.nansum(df['pop'])

UN = [131107 * 1000, 165704.243 * 1000]
#from WPP2019_POP_F07_1_Population_By_AGE_BOTH_SEXES

#%%

fig = plt.figure(figsize=(9,9))

years = np.unique(ssp2['Year'])
years = years[years >= 2000]

plt.plot(years, tpop, label = 'SSP2 Adjusted')
plt.plot(years, tpop_raw, label = 'SSP2')

#plt.plot(years, tpop3, label = 'SSP3')
plt.scatter([2000, 2010], [x2000, x2010], label = 'WorldPop')
plt.scatter([2010], [gx2010], label = 'GPW', c = 'black')
plt.scatter([2000, 2010], UN, label = 'UN', c = 'red')


plt.xlabel('Year')
plt.ylabel('Total African population (under 5)')

plt.legend(bbox_to_anchor=(0.85, -0.1),ncol=5,frameon = False, handletextpad = 0.5)

#plt.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Inputdata_SSP2_corrected_tasgrid_afonly_morestimates.png',
#         bbox_inches = 'tight', pad_inches = 0.3)

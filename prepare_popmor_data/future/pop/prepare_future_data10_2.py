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
pop2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2019_regrid.nc')
pop2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2010_regrid.nc')
pop2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000_regrid.nc')

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

countries_regrid = countries.regrid(pop2019, iris.analysis.Nearest())



print(np.nansum(pop2000.data))
print(np.nansum(pop2010.data))

#%%apply countries mask to pop data

p2019_regrid = copy.deepcopy(pop2019)
p2010_regrid = copy.deepcopy(pop2010)
p2000_regrid = copy.deepcopy(pop2000)

p2019_regrid.data = np.ma.masked_array(p2019_regrid.data, countries_regrid.data.mask)
p2010_regrid.data = np.ma.masked_array(p2010_regrid.data, countries_regrid.data.mask)
p2000_regrid.data = np.ma.masked_array(p2000_regrid.data, countries_regrid.data.mask)


#%% Find country totals


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
    


#%% recreate popfraoc

country_totalpop = country_total(p2019_regrid, countries_regrid)
cnames = country_totalpop['Area'].values
   
popfrac_2019 = copy.deepcopy(p2019_regrid)
dims = popfrac_2019.shape    

# get all country codes from map
dat = popfrac_2019.data
vals = np.unique(dat)
vals = vals[vals.mask == False]

        
for i in np.arange(dims[0]):
    for j in np.arange(dims[1]):
        
        val = countries_regrid[i,j].data
        if np.ma.is_masked(val) == False:
            c_name = c_dict[int(val)]
            if c_name in cnames:
    
    
                p = country_totalpop[country_totalpop['Area'] == c_name]['pop'].values
                replace_val = p2019_regrid[i,j].data / p
                popfrac_2019.data[i,j] = replace_val


#%% erplace country code with total under 5 population

save_path = '/nfs/a321/earsch/CHAMNHA/input_data/pop/future/processed/'


def distr_pop(pop_table, scenario_name):


    years = np.unique(pop_table['Year'])
    years = years[years >= 2000]
    
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
            
        pop_year_frac_regrid = pop_year_frac.regrid(tas, iris.analysis.AreaWeighted())
        pop_output.append(pop_year_frac_regrid)
    return pop_output

pop_output = distr_pop(ssp2, 'ssp2')

#%% Modify pop2019
 #calc frac chagne from 2000 for each ssp2 year

base = pop_output[0]

fracchange =  iris.cube.CubeList()
for cube in pop_output:
    x = cube/base
    fracchange.append(x)

#apply frac change to correctly gridded input pop2090 data
actual_pop2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000_regrid.nc')
actual_pop2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2019_regrid.nc')


fracchange_regrid =  iris.cube.CubeList()
for cube in fracchange:
    x = cube.regrid(actual_pop2019, iris.analysis.AreaWeighted())
    fracchange_regrid.append(x)

#calc avg change (do it this way as 0s in egypt messing up avg change otherwise)
total_pop = [np.nansum(x.data) for x in pop_output]
base_pop = np.nansum(base.data) 

avg_change = total_pop / base_pop
    

act_mask = actual_pop2000.data.mask
new_mask = fracchange_regrid[0].data.mask

act_mask = np.where(act_mask == True, 1, 0)
new_mask = np.where(new_mask == True, 1, 0)
mask_dif = act_mask - new_mask

new_dat = iris.cube.CubeList()
for i in np.arange(len(fracchange_regrid)):
    cube = fracchange_regrid[i]
    
    #change all gridcells, even masked ones
    new_cube = actual_pop2000 * avg_change[i]
    #change specifically ones have data for
    x = cube * actual_pop2000

    #where mask si the same, use fracchange * pop2000, otherwise use
    # pop2000 * avg_change
    new_cube.data = np.ma.where(mask_dif == 0, x.data, new_cube.data)    
    new_dat.append(new_cube)
    

years = np.unique(ssp2['Year'])
years = years[years >= 2000]

save_path = '/nfs/a321/earsch/CHAMNHA/input_data/pop/future/processed/'


for i in np.arange(len(new_dat)):
    y = years[i]
    save_name = 'ssp2' + '_' + str(y) + '_04population_mf_BIASCORR4.nc'
    #iris.save(new_dat[i], save_path + save_name)



#%% restrict to afric only

not_in_afr = ['Albania', 'Armenia', 'Azerbaijan', 'Cyprus', 'France', 'Greece', 
              'Iran', 'Iraq', 'Israel', 'Italy', 'Jordan', 'Kuwait', 'Lebanon',
              'Malta', 'Northern Cyprus', 'Oman', 'Palestine', 'Portugal', 'Qatar',
              'Saudi Arabia', 'Spain', 'Syria', 'Turkey', 'Turkmenistan', 'United Arab Emirates', 'Yemen']



#%%

tpop = []
tpop_raw = []

cregrid = countries_regrid.regrid(pop_output[0], iris.analysis.Nearest())

for i in np.arange(len(new_dat)):
    df = country_total(new_dat[i], cregrid)
    df = df[~df['Area'].isin(not_in_afr)]
    
    x = np.nansum(df['pop'])
    tpop.append(x)
    
    df = country_total(pop_output[i], cregrid)
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

#%%

fig = plt.figure(figsize=(9,9))

years = np.unique(ssp2['Year'])
years = years[years >= 2000]

plt.plot(years, tpop, label = 'SSP2 Corrected')
plt.plot(years, tpop_raw, label = 'SSP2')

#plt.plot(years, tpop3, label = 'SSP3')
plt.scatter([2000, 2010], [x2000, x2010], label = 'Actual')
plt.xlabel('Year')
plt.ylabel('Total African population (under 5)')

plt.legend(bbox_to_anchor=(0.95, -0.1),ncol=3,frameon = False, handletextpad = 0.5)

#plt.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Inputdata_SSP2_corrected_tasgrid_afonly.png',
#         bbox_inches = 'tight', pad_inches = 0.3)

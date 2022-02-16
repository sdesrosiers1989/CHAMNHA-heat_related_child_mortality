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

#not regridded to CANeSM - units per km2
countries = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/africa_countryname.nc')
pop2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2000.nc')
pop2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2010.nc')
pop2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/pop/processed/afr_01_mf_2019.nc')

gs = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/africa_cname_gridarea.nc')
gs.convert_units('km2')

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



#%% regrid data

min_lat = np.nanmin(pop2019.coord('latitude').points)
max_lat = np.nanmax(pop2019.coord('latitude').points)
min_lon = np.nanmin(pop2019.coord('longitude').points)
max_lon = np.nanmax(pop2019.coord('longitude').points)

def safrica_lat(input):
    return min_lat  <= input <= max_lat 

def safrica_long(input):
    return min_lon  <= input <= max_lon 

af_con = iris.Constraint(latitude = safrica_lat, longitude = safrica_long)

countries = countries.extract(af_con)
gs = gs.extract(af_con)

countries.coord('latitude').guess_bounds()
countries.coord('longitude').guess_bounds()
gs.coord('latitude').guess_bounds()
gs.coord('longitude').guess_bounds()

cs = countries.coord_system(iris.coord_systems.CoordSystem)

pop_regrid = []
pop_list = [pop2000, pop2010, pop2019]
for cube in pop_list:
    try:
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
    except:
        print('Has bounds')
    cube.coord('longitude').coord_system = cs
    cube.coord('latitude').coord_system = cs
    x = cube.regrid(countries, iris.analysis.AreaWeighted())
    x = x * gs
    pop_regrid.append(x)



#%%apply countries mask to pop data - cuts bits of some countries

p2019_regrid = copy.deepcopy(pop_regrid[2])
p2010_regrid = copy.deepcopy(pop_regrid[1])
p2000_regrid = copy.deepcopy(pop_regrid[0])


p2019_regrid.data = np.ma.masked_array(p2019_regrid.data, countries.data.mask)
p2010_regrid.data = np.ma.masked_array(p2010_regrid.data, countries.data.mask)
p2000_regrid.data = np.ma.masked_array(p2000_regrid.data, countries.data.mask)



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
    
df_2000 = country_total(p2000_regrid, countries)
df_2010 = country_total(p2010_regrid, countries)
df_2019 = country_total(p2019_regrid, countries)

df_ssp_2000 = ssp2[ssp2['Year'] == 2000]

both2000 = pd.merge(df_ssp_2000, df_2000, on = 'Area')
both2000['dif'] = both2000['pop'] / both2000['Population']

df_ssp_2010 = ssp2[ssp2['Year'] == 2010]

both2010 = pd.merge(df_ssp_2010, df_2010, on = 'Area')
both2010['dif'] = both2010['pop'] / both2010['Population']

both_years = pd.merge(both2000, both2010, on = ['Area'])
both_years['avg_dif'] = (both_years['dif_x'] + both_years['dif_y']) / 2


#%%

c = df_ssp_2000['Area']
c2 = df_2000['Area'].values

[x for x in c if x not in c2]

#%% recreate popfraoc

country_totalpop = country_total(p2019_regrid, countries)
cnames = country_totalpop['Area'].values
   
popfrac_2019 = copy.deepcopy(p2019_regrid)
dims = popfrac_2019.shape    

# get all country codes from map
dat = popfrac_2019.data
vals = np.unique(dat)
vals = vals[vals.mask == False]

        
for i in np.arange(dims[0]):
    for j in np.arange(dims[1]):
        
        val = countries[i,j].data
        if np.ma.is_masked(val) == False:
            c_name = c_dict[int(val)]
            if c_name in cnames:
    
    
                p = country_totalpop[country_totalpop['Area'] == c_name]['pop'].values
                replace_val = p2019_regrid[i,j].data / p
                popfrac_2019.data[i,j] = replace_val


#%% erplace country code with total under 5 population

save_path = '/nfs/a321/earsch/CHAMNHA/input_data/pop/future/processed/'

def distr_pop_bc(pop_table, scenario_name, bias_corr):


    years = np.unique(pop_table['Year'])
    years = years[years >= 2000]
    
    pop_output = []
    
    for y in years:
    
        year_dat = pop_table[pop_table['Year'] == y]
        
        #adjust
        year_dat_merge = pd.merge(year_dat, bias_corr, on = 'Area')
        year_dat_merge['Population_corr'] = year_dat_merge['Population'] * year_dat_merge['avg_dif']
        
        
        cnames = np.unique(year_dat['Area'])
        
        pop_year = copy.deepcopy(countries)
        
        # get all country codes from map
        dat = pop_year.data
        vals = np.unique(dat)
        vals = vals[vals.mask == False]
        
        for val in vals:
            c_name = c_dict[val]
            #print(c_name)
            
            if c_name in cnames:
        
                replace_val = year_dat_merge[year_dat_merge['Area'] == c_name]['Population_corr'].values
                if replace_val < 0:
                    print(replace_val)
                
                pop_year.data = np.where(pop_year.data == val, replace_val, pop_year.data)
                
                pop_year_frac = copy.deepcopy(pop_year)
                pop_year_frac.data = pop_year_frac.data*popfrac_2019.data
                
        #save_name = scenario_name + '_' + str(y) + '_04population_mf_BIASCORR3.nc'
        #pop_year_frac_regrid = pop_year_frac.regrid(tas, iris.analysis.AreaWeighted())
        
        #iris.save(pop_year_frac_regrid, save_path + save_name)
        
        pop_output.append(pop_year_frac)
    return pop_output

def distr_pop(pop_table, scenario_name):


    years = np.unique(pop_table['Year'])
    years = years[years >= 2000]
    
    pop_output = []
    
    for y in years:
    
        year_dat = pop_table[pop_table['Year'] == y]
        cnames = np.unique(year_dat['Area'])
        
        pop_year = copy.deepcopy(countries)
        
        # get all country codes from map
        dat = pop_year.data
        vals = np.unique(dat)
        vals = vals[vals.mask == False]
        
        for val in vals:
            c_name = c_dict[val]
            #print(c_name)
            
            if c_name in cnames:
        
                replace_val = year_dat[year_dat['Area'] == c_name]['Population'].values
                #print(replace_val)
                
                pop_year.data = np.where(pop_year.data == val, replace_val, pop_year.data)
                
                pop_year_frac = copy.deepcopy(pop_year)
                pop_year_frac.data = pop_year_frac.data*popfrac_2019.data
            
        #pop_year_frac_regrid = pop_year_frac.regrid(tas, iris.analysis.AreaWeighted())
        pop_output.append(pop_year_frac)
    return pop_output


pop_output = distr_pop_bc(ssp2, 'ssp2', both_years)
pop_output_raw = distr_pop(ssp2, 'ssp2')

#%%

check = country_total(pop_output[9], countries)

#%%  check
        
pop_levs = [0, 3000, 20000, 30000, 50000, 100000, 250000, 500000] 

qplt.contourf(pop2019, levels = pop_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(p2019_regrid, levels = pop_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(pop_output[4], levels = pop_levs, extend = 'max', cmap = 'YlOrRd')


qplt.contourf(pop2000, levels = pop_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(p2000_regrid, levels = pop_levs, extend = 'max', cmap = 'YlOrRd')

#%% restrict to afric only
countries = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/africa_countires_only.nc')


country_lookup = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/pop/country_lookup.csv')
country_lookup['Name'][country_lookup['Name'] == 'Republic of the Congo'] = 'Congo'
country_lookup['Name'][country_lookup['Name'] == 'eSwatini'] = 'Swaziland'


not_in_afr = ['Albania', 'Armenia', 'Azerbaijan', 'Cyprus', 'France', 'Greece', 
              'Iran', 'Iraq', 'Israel', 'Italy', 'Jordan', 'Kuwait', 'Lebanon',
              'Malta', 'Northern Cyprus', 'Oman', 'Palestine', 'Portugal', 'Qatar',
              'Saudi Arabia', 'Spain', 'Syria', 'Turkey', 'Turkmenistan', 'United Arab Emirates', 'Yemen']

c_dict = dict(zip(country_lookup['Value'], country_lookup['Name']))

#obs - used to get threshold
cru_tas = iris.load('/nfs/a321/earsch/Tanga/Data/CRU/tmp/*.nc',
                    iris.Constraint(cube_func = lambda cube: cube.var_name == 'tmp'))

cru_tas = cru_tas[0][0]
cru_tas = cru_tas.regrid(pop2019_regrid, iris.analysis.Linear())

#%%

pop_2019.data.mask[np.isnan(countries.data)] = True
pop_2019.data = np.ma.masked_where(np.ma.getmask(cru_tas.data), pop2019.data)

pop_2045.data.mask[np.isnan(countries.data)] = True
pop_2045.data = np.ma.masked_where(np.ma.getmask(cru_tas.data), pop2045.data)


mor_2019.data.mask[np.isnan(countries.data)] = True
mor_2019.data = np.ma.masked_where(np.ma.getmask(cru_tas.data), mor_2019.data)

mor_2045.data.mask[np.isnan(countries.data)] = True
mor_2045.data = np.ma.masked_where(np.ma.getmask(cru_tas.data), mor_2045.data)

rmor_2019.data.mask[np.isnan(countries.data)] = True
rmor_2019.data = np.ma.masked_where(np.ma.getmask(cru_tas.data), rmor_2019.data)

rmor_2045.data.mask[np.isnan(countries.data)] = True
rmor_2045.data = np.ma.masked_where(np.ma.getmask(cru_tas.data), rmor_2045.data)
    

#%%
print(np.nanmean(pop2019.data))
print(np.nanmean(p2019_regrid.data))
mpop = [np.nanmean(x.data) for x in pop_output]
#mpop3 = [np.nanmean(x.data) for x in pop_output_ssp3]


print(np.nansum(pop2019.data))
print(np.nansum(p2019_regrid.data))
print(np.nansum(df_2019['pop']))
tpop = [np.nansum(x.data) for x in pop_output]
tpop_raw = [np.nansum(x.data) for x in pop_output_raw]

#tpop3 = [np.nansum(x.data) for x in pop_output_ssp3]

#%%

years = np.unique(ssp2['Year'])
years = years[years >= 2000]

plt.plot(years, tpop, label = 'SSP2 Corrected')
plt.plot(years, tpop_raw, label = 'SSP2')

#plt.plot(years, tpop3, label = 'SSP3')
plt.scatter([2000, 2010], [np.nansum(p2000_regrid.data), np.nansum(p2010_regrid.data)], label = 'Actual')
plt.xlabel('Year')
plt.ylabel('Total African population (under 5)')
plt.legend()

#plt.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Inputdata_SSP2_corrected_tasgrid.png',
#         bbox_inches = 'tight', pad_inches = 0.3)

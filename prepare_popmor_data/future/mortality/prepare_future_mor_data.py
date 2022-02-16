#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process future mortality data (GBD)

-in csv format
-link to coutnry raster(dictionary)
-distribution using 2019 pop distribution

-extrapolate to 2050

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
country_lookup['Name'][country_lookup['Name'] == 'eSwatini'] = 'Swaziland'
#country_lookup['Name'][country_lookup['Name'] == 'Northern Cyprus'] = 'Cyprus'
    #already cyprus val, which in other dates is cyprus + nrorthern cyprus

mor_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2019_regrid.nc')
mor_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2010_regrid.nc')
mor_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2000_regrid.nc')

mor_ref = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_reference_both_04.csv')
mor_better = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_better_both_04.csv')
mor_worse = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_worse_both_04.csv')

mor_ref_ex = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_reference_both_04_extracountries.csv')
mor_better_ex = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_better_both_04_extracountries.csv')
mor_worse_ex = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_worse_both_04_extracountries.csv')

#%%

n = np.unique(mor_ref['Location'].values)
mor_ref_ex = mor_ref_ex[~mor_ref_ex['Location'].isin(n)]

mor_ref = mor_ref.append(mor_ref_ex)

#%% Create diciotnary of country codes and names

c_dict = dict(zip(country_lookup['Value'], country_lookup['Name']))

#countires in countries not in mor_ref
names = country_lookup['Name']
ref_names = np.unique(mor_ref['Location'])

[x for x in names if x not in ref_names]

#[x for x in ref_names if x not in ref_names]

#mor_ref = mor_ref.drop(mor_ref[mor_ref['Location']  not in names])


#%% regrid countries to poprfac

countries_regrid = countries.regrid(popfrac_2019, iris.analysis.Nearest())
countries_regrid.coord('latitude').guess_bounds()
countries_regrid.coord('longitude').guess_bounds()


popfrac_2019.coord('latitude').guess_bounds()
popfrac_2019.coord('longitude').guess_bounds()

pop2019 = pop2019.regrid(popfrac_2019, iris.analysis.AreaWeighted())

m2000_regrid = mor_2000.regrid(popfrac_2019, iris.analysis.AreaWeighted())
m2010_regrid = mor_2010.regrid(popfrac_2019, iris.analysis.AreaWeighted())
m2019_regrid = mor_2019.regrid(popfrac_2019, iris.analysis.AreaWeighted())


#%%

popfrac_2019.data = np.ma.masked_array(popfrac_2019.data, countries_regrid.data.mask)
pop2019.data = np.ma.masked_array(pop2019.data, countries_regrid.data.mask)


m2019_regrid.data = np.ma.masked_array(m2019_regrid.data, countries_regrid.data.mask)
m2010_regrid.data = np.ma.masked_array(m2010_regrid.data, countries_regrid.data.mask)
m2000_regrid.data = np.ma.masked_array(m2000_regrid.data, countries_regrid.data.mask)



#%% extrapolate to 2050

def extend(in_table, interp_sy = 2020, end_year = 2050, country_names = names):
    
    for n in country_names:
        c_df = in_table[in_table['Location'] == n]
        
        #calc line
        
        c_df_projperiod = c_df[c_df['Year'] >= interp_sy]
        #m, b = np.polyfit(df['rp'].values, df['rl'].values, 1)
        z = np.polyfit(c_df_projperiod['Year'].values, c_df_projperiod['Value'].values, 2)
        p = np.poly1d(z)
        
        years = np.arange(2041, end_year +1)
        proj = p(years)
                
        
        df = pd.DataFrame({'Location': n,
                           'Year': years,
                           'Value': proj })
        in_table = in_table.append(df)
        
    return in_table


mor_ref_ex = extend(mor_ref, country_names = [x for x in names if x != 'Northern Cyprus'])
    

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
both2000['frac'] =  both2000['mor'] / both2000['Value'] 

df_ref_2010 = mor_ref[mor_ref['Year'] == 2010]

both2010 = pd.merge(df_ref_2010, df_2010, on = 'Location')
both2010['frac'] = both2010['mor']/ both2010['Value'] 

both_years = pd.merge(both2000, both2010, on = ['Location'])
both_years['avg_frac'] = (both_years['frac_x'] + both_years['frac_y']) / 2

cols = both_years.columns
keep = ['Location', 'Year_x', 'Value_x', 'Value_y', 'Year_y', 'mor_x', 'mor_y', 'frac_x', 'frac_y', 'avg_frac']
cols = [x for x in cols if x not in keep]
both_years.drop(columns = cols, inplace =  True)

np.nansum(both_years['Value_x'])
np.nansum(df_2000['mor'])
np.nansum(both_years['Value_x'] * both_years['frac_x'])
np.nansum(both_years['Value_x'] * both_years['avg_frac'])
both_years['check_2000'] = both_years['Value_x'] * both_years['avg_frac']

np.nansum(both_years['Value_y'])
np.nansum(df_2010['mor'])
np.nansum(both_years['Value_y'] * both_years['frac_y'])
np.nansum(both_years['Value_y'] * both_years['avg_frac'])




#%% recreate popfraoc

country_totalpop = country_total(pop2019, countries_regrid)
cnames = country_totalpop['Location'].values
   
popfrac_2019 = copy.deepcopy(pop2019)
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
    
    
                p = country_totalpop[country_totalpop['Location'] == c_name]['mor'].values
                replace_val = pop2019[i,j].data / p
                popfrac_2019.data[i,j] = replace_val




#%% erplace country code with total under 5 population

save_path = '/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/processed/'

def distr_dat_bc(in_table, scenario_name, bias_corr):


    years = np.unique(in_table['Year'])
    years = years[years >= 2000]
    
    output = []
    
    for y in years:
    
        year_dat = in_table[in_table['Year'] == y]
        #adjust
        year_dat_merge = pd.merge(year_dat, bias_corr, on = 'Location')
        year_dat_merge['Mor_corr'] = year_dat_merge['Value'] * year_dat_merge['avg_frac']
        
        
        
        cnames = np.unique(year_dat['Location'])
        
        mor_year = copy.deepcopy(countries_regrid)
        
        # get all country codes from map
        dat = mor_year.data
        vals = np.unique(dat)
        vals = vals[vals.mask == False]
        
        for val in vals:
            c_name = c_dict[val]
            
            
            if c_name in cnames:
        
                replace_val = year_dat_merge[year_dat_merge['Location'] == c_name]['Mor_corr'].values
                
                if replace_val < 0:
                    print( y, c_name, replace_val)
                    #replace_val = 0
                
                if y == 2000:
                    print(c_name, replace_val)
                
                mor_year.data = np.where(mor_year.data == val, replace_val, mor_year.data)
                
                mor_year_frac = copy.deepcopy(mor_year)
                mor_year_frac.data = mor_year_frac.data*popfrac_2019.data
            else:
                print(c_name)
                
        save_name = scenario_name + '_' + str(int(y)) + '_04_totalmor_mf_BIASCORR.nc'
        iris.save(mor_year_frac, save_path + save_name)
        
        output.append(mor_year_frac)
    return output

def distr_dat(in_table, scenario_name):


    years = np.unique(in_table['Year'])
    years = years[years >= 2000]
    
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
            #print(c_name)
            
            if c_name in cnames:
        
                replace_val = year_dat[year_dat['Location'] == c_name]['Value'].values
                #print(replace_val)
                
                mor_year.data = np.where(mor_year.data == val, replace_val, mor_year.data)
                
                mor_year_frac = copy.deepcopy(mor_year)
                mor_year_frac.data = mor_year_frac.data*popfrac_2019.data
                
        
        output.append(mor_year_frac)
    return output

mor_ref_output = distr_dat_bc(mor_ref_ex, 'ref', both_years)
mor_ref_output_raw = distr_dat(mor_ref_ex, 'ref')

#pop_output_ssp3 = distr_pop(ssp3, 'ssp3')



#%%  check
        
pop_levs = [0, 3000, 20000, 30000, 50000, 100000, 250000, 500000] 
mor_levs = np.array([0, 1, 4, 6, 10, 15, 30]) * 365


qplt.contourf(mor_2019, levels = mor_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(mor_ref_output[19], levels = mor_levs, extend = 'max', cmap = 'YlOrRd')


x = country_total(mor_ref_output[0], countries_regrid)
x_check = pd.merge(x, df_2000, on = 'Location')

#%%
print(np.nanmean(mor_2019.data))
mmor = [np.nanmean(x.data) for x in mor_ref_output]
#mmor3 = [np.nanmean(x.data) for x in pop_output_ssp3]


print(np.nansum(mor_2019.data))
print(np.nansum(m2019_regrid.data))
print(np.nansum(df_2019['mor']))
tmor = [np.nansum(x.data) for x in mor_ref_output]
tmor_raw = [np.nansum(x.data) for x in mor_ref_output_raw]
#tmor3 = [np.nansum(x.data) for x in pop_output_ssp3]

#%%

actual_x = [2000, 2010, 2019]
actual_y = [np.nansum(m2000_regrid.data),
            np.nansum(m2010_regrid.data),
            np.nansum(m2019_regrid.data)]

years = np.arange(2000, 2051)

plt.plot(years, tmor, label = 'GBD Ref Corrected')
plt.plot(years, tmor_raw, label = 'GBD Ref')

plt.scatter(actual_x, actual_y, label = 'Actual')
plt.xlabel('Year')
plt.ylabel('Total African mortality (under 5)')
plt.legend()

plt.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Inputdata_GBDRef_corrected.png',
         bbox_inches = 'tight', pad_inches = 0.3)
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


#regrid to CMIP6
scen = '119'
path = '/nfs/a321/earsch/Tanga/Data/CMIP6/bias_corr/CHAMNHA/tas/end/'
filenames = glob.glob(path + '*' + scen + '.nc')
tas = iris.load_cube(filenames[0])[0]
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

countries_regrid = countries.regrid(mor_2019, iris.analysis.Nearest())





#%%

pop2019.data = np.ma.masked_array(pop2019.data, countries_regrid.data.mask)

m2019_regrid = copy.deepcopy(mor_2019)
m2010_regrid = copy.deepcopy(mor_2010)
m2000_regrid = copy.deepcopy(mor_2000)

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
        
        if np.nanmin(proj) <= 0:
            new_min = np.min(proj[proj> 0])
            proj[proj <= 0] = new_min
                
        
        df = pd.DataFrame({'Location': n,
                           'Year': years,
                           'Value': proj })
        in_table = in_table.append(df)
        
    return in_table


mor_ref_ex = extend(mor_ref, country_names = [x for x in names if x != 'Northern Cyprus'], end_year = 2060)
    

x = mor_ref_ex[mor_ref_ex['Year'] == 2055]
y = mor_ref_ex[mor_ref_ex['Location'] == 'Equatorial Guinea']
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
        mor_year_frac_regrid = mor_year_frac.regrid(tas, iris.analysis.AreaWeighted())
        
        #if int(y) > 2050:
        #    save_name = scenario_name + '_' + str(int(y)) + '_04_totalmor_mf_BIASCORR.nc'      
        #    iris.save(mor_year_frac_regrid, save_path + save_name)
        
        output.append(mor_year_frac_regrid)
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

        mor_year_frac_regrid = mor_year_frac.regrid(tas, iris.analysis.AreaWeighted())
    
        output.append(mor_year_frac_regrid)
        
    return output

mor_ref_output = distr_dat_bc(mor_ref_ex, 'ref', both_years)
mor_ref_output_raw = distr_dat(mor_ref_ex, 'ref')

#pop_output_ssp3 = distr_pop(ssp3, 'ssp3')

#%% Modify pop2019
 #calc frac chagne from 2000 for each ssp2 year

base =  mor_ref_output[0]

fracchange=  iris.cube.CubeList()
for cube in mor_ref_output:
    x = cube/base
    fracchange.append(x)

#apply frac change to correctly gridded input pop2090 data
actual_mor2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2019_regrid.nc')
actual_mor2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2010_regrid.nc')
actual_mor2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2000_regrid.nc')


fracchange_regrid =  iris.cube.CubeList()
for cube in fracchange:
    x = cube.regrid(actual_mor2019, iris.analysis.AreaWeighted())
    fracchange_regrid.append(x)

#calc avg change (do it this way as 0s in egypt messing up avg change otherwise)
total_mor = [np.nansum(x.data) for x in mor_ref_output]
base_mor = np.nansum(base.data) 

avg_change = total_mor / base_mor
    

act_mask = actual_mor2000.data.mask
new_mask = fracchange_regrid[0].data.mask

act_mask = np.where(act_mask == True, 1, 0)
new_mask = np.where(new_mask == True, 1, 0)
mask_dif = act_mask - new_mask

new_dat = iris.cube.CubeList()
for i in np.arange(len(fracchange_regrid)):
    cube = fracchange_regrid[i]
    
    #change all gridcells, even masked ones
    new_cube = actual_mor2000 * avg_change[i]
    #change specifically ones have data for
    x = cube * actual_mor2000

    #where mask si the same, use fracchange * pop2000, otherwise use
    # pop2000 * avg_change
    new_cube.data = np.ma.where(mask_dif == 0, x.data, new_cube.data)    
    new_dat.append(new_cube)
    

years = np.unique(mor_ref_ex['Year'])
years = years[years >= 2000]

save_path = '/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/processed/'


for i in np.arange(len(new_dat)):
    y = years[i]
    save_name = 'ref' + '_' + str(y)[0:4] + '_04_totalmor_mf_BIASCORR2.nc' 
    iris.save(new_dat[i], save_path + save_name)


#%%  check
        
pop_levs = [0, 3000, 20000, 30000, 50000, 100000, 250000, 500000] 
mor_levs = np.array([0, 1, 4, 6, 10, 15, 30]) * 365


qplt.contourf(mor_2019, levels = mor_levs, extend = 'max', cmap = 'YlOrRd')
qplt.contourf(mor_ref_output[19], levels = mor_levs, extend = 'max', cmap = 'YlOrRd')


x = country_total(mor_ref_output[0], countries_regrid)
x_check = pd.merge(x, df_2000, on = 'Location')

#%%

years = np.unique(mor_ref_ex['Year'])
years = years[years >= 2000]

mor_ref_output[35]= mor_ref_output[35].regrid(countries_regrid, iris.analysis.AreaWeighted())

check = country_total(mor_ref_output[35], countries_regrid)


#%% restrict to afric only

not_in_afr = ['Albania', 'Armenia', 'Azerbaijan', 'Cyprus', 'France', 'Greece', 
              'Iran', 'Iraq', 'Israel', 'Italy', 'Jordan', 'Kuwait', 'Lebanon',
              'Malta', 'Northern Cyprus', 'Oman', 'Palestine', 'Portugal', 'Qatar',
              'Saudi Arabia', 'Spain', 'Syria', 'Turkey', 'Turkmenistan', 'United Arab Emirates', 'Yemen']



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

tmor = []
tmor_raw =[]

cregrid = countries_regrid.regrid(mor_ref_output[0], iris.analysis.Nearest())

for i in np.arange(len(mor_ref_output)):
    print(i)
    df = country_total(mor_ref_output[i], cregrid)
    df = df[~df['Location'].isin(not_in_afr)]
    
    x = np.nansum(df['mor'])
    tmor.append(x)
    
    df = country_total(mor_ref_output_raw[i], cregrid)
    df = df[~df['Location'].isin(not_in_afr)]
    
    x = np.nansum(df['mor'])
    tmor_raw.append(x)


#actual 200 and 2010
df = country_total(m2000_regrid, countries_regrid)
df = df[~df['Location'].isin(not_in_afr)]
x2000 = np.nansum(df['mor'])

df = country_total(m2010_regrid, countries_regrid)
df = df[~df['Location'].isin(not_in_afr)]
x2010 = np.nansum(df['mor'])

df = country_total(m2019_regrid, countries_regrid)
df = df[~df['Location'].isin(not_in_afr)]
x2019 = np.nansum(df['mor'])

#%%

fig = plt.figure(figsize=(9,9))

actual_x = [2000, 2010, 2019]
actual_y = [x2000,
            x2010,
            x2019]

years = np.arange(2000, 2061)

plt.plot(years, tmor, label = 'GBD Ref Corrected')
plt.plot(years, tmor_raw, label = 'GBD Ref')

plt.scatter(actual_x, actual_y, label = 'Actual')
plt.xlabel('Year')
plt.ylabel('Total African mortality (under 5)')
plt.legend(bbox_to_anchor=(0.95, -0.1),ncol=3,frameon = False, handletextpad = 0.5)


#plt.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Inputdata_GBDRef_corrected_regrid_afonly2060.png',
#         bbox_inches = 'tight', pad_inches = 0.3)

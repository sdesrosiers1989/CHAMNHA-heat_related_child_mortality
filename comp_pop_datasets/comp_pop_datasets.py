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
#pop data

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


#mortality data
mor_2019 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2019_regrid.nc')
mor_2010 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2010_regrid.nc')
mor_2000 = iris.load_cube('/nfs/a321/earsch/CHAMNHA/input_data/mortality/processed/total_mor_mf_01_2000_regrid.nc')

mor_ref = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_reference_both_04.csv')
mor_ref_ex = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/vis_totaldeaths_reference_both_04_extracountries.csv')


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




#%% regrid countries to poprfac

countries_regrid = countries.regrid(popfrac_2019, iris.analysis.Nearest())
countries_regrid.coord('latitude').guess_bounds()
countries_regrid.coord('longitude').guess_bounds()


popfrac_2019.coord('latitude').guess_bounds()
popfrac_2019.coord('longitude').guess_bounds()
p2000_regrid = pop2000.regrid(popfrac_2019, iris.analysis.AreaWeighted())
p2010_regrid = pop2010.regrid(popfrac_2019, iris.analysis.AreaWeighted())
p2019_regrid = pop2019.regrid(popfrac_2019, iris.analysis.AreaWeighted())

m2000_regrid = mor_2000.regrid(popfrac_2019, iris.analysis.AreaWeighted())
m2010_regrid = mor_2010.regrid(popfrac_2019, iris.analysis.AreaWeighted())
m2019_regrid = mor_2019.regrid(popfrac_2019, iris.analysis.AreaWeighted())

print(np.nansum(pop2000.data), np.nansum(p2000_regrid.data))
print(np.nansum(pop2010.data), np.nansum(p2010_regrid.data))

#%%apply countries mask to pop data

p2019_regrid.data = np.ma.masked_array(p2019_regrid.data, countries_regrid.data.mask)
p2010_regrid.data = np.ma.masked_array(p2010_regrid.data, countries_regrid.data.mask)
p2000_regrid.data = np.ma.masked_array(p2000_regrid.data, countries_regrid.data.mask)

gpw_2010.data = np.ma.masked_array(gpw_2010.data, countries_regrid.data.mask)


m2019_regrid.data = np.ma.masked_array(m2019_regrid.data, countries_regrid.data.mask)
m2010_regrid.data = np.ma.masked_array(m2010_regrid.data, countries_regrid.data.mask)
m2000_regrid.data = np.ma.masked_array(m2000_regrid.data, countries_regrid.data.mask)

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

#change lookup to match country names in GBD data
country_lookup = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/pop/country_lookup.csv')
country_lookup['Name'][country_lookup['Name'] == 'Republic of the Congo'] = 'Congo'
country_lookup['Name'][country_lookup['Name'] == 'eSwatini'] = 'Swaziland'

#countires in countries not in mor_ref
names = country_lookup['Name']
ref_names = np.unique(mor_ref['Location'])

[x for x in names if x not in ref_names]

mor_ref_ex = extend(mor_ref, country_names = [x for x in names if x != 'Northern Cyprus'], end_year = 2060)
    

x = mor_ref_ex[mor_ref_ex['Year'] == 2055]
y = mor_ref_ex[mor_ref_ex['Location'] == 'Equatorial Guinea']

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

path = '/nfs/a321/earsch/CHAMNHA/input_data/mortality/future/processed/'
file_names = glob.glob(path + '*_04_totalmor_mf_BIASCORR.nc')
mor_output = iris.load(file_names)

file_names = glob.glob(path + 'ref*_04_totalmor_mf.nc')
file_names = [x for x in file_names if '.0' not in x]
mor_output_raw = iris.load(file_names)

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


tmor = []
tmor_raw =[]

cregrid = countries_regrid.regrid(mor_output[0], iris.analysis.Nearest())
cregrid_raw = countries_regrid.regrid(mor_output_raw[0], iris.analysis.Nearest())

for i in np.arange(len(mor_output)):
    df = country_total(mor_output[i], cregrid)
    df = df[~df['Area'].isin(not_in_afr)]
    
    x = np.nansum(df['pop'])
    tmor.append(x)
    
    df = country_total(mor_output_raw[i], cregrid_raw)
    df = df[~df['Area'].isin(not_in_afr)]
    
    x = np.nansum(df['pop'])
    tmor_raw.append(x)

#actual 200 and 2010
df = country_total(p2000_regrid, countries_regrid)
df = df[~df['Area'].isin(not_in_afr)]
x2000 = np.nansum(df['pop'])

df = country_total(p2010_regrid, countries_regrid)
df = df[~df['Area'].isin(not_in_afr)]
x2010 = np.nansum(df['pop'])

df = country_total(p2019_regrid, countries_regrid)
df = df[~df['Area'].isin(not_in_afr)]
x2019 = np.nansum(df['pop'])

df = country_total(gpw_2010, countries_regrid)
df = df[~df['Area'].isin(not_in_afr)]
gx2010 = np.nansum(df['pop'])

#mor 200 and 2010
df = country_total(m2000_regrid, countries_regrid)
df = df[~df['Area'].isin(not_in_afr)]
mx2000 = np.nansum(df['pop'])

df = country_total(m2010_regrid, countries_regrid)
df = df[~df['Area'].isin(not_in_afr)]
mx2010 = np.nansum(df['pop'])

df = country_total(m2019_regrid, countries_regrid)
df = df[~df['Area'].isin(not_in_afr)]
mx2019 = np.nansum(df['pop'])

df = country_total(gpw_2010, countries_regrid)
df = df[~df['Area'].isin(not_in_afr)]
gx2010 = np.nansum(df['pop'])

#%%
un_csv = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES_SMALL.csv')

country_lookup = pd.read_csv('/nfs/a321/earsch/CHAMNHA/input_data/pop/country_lookup.csv')
country_lookup['Name'][country_lookup['Name'] == 'The Gambia'] = 'Gambia'
country_lookup['Name'][country_lookup['Name'] == 'Republic of the Congo'] = 'Congo'
country_lookup['Name'][country_lookup['Name'] == 'eSwatini'] = 'Eswatini'
country_lookup['Name'][country_lookup['Name'] == 'Tanzania'] = 'United Republic of Tanzania'

un_csv = un_csv[un_csv['Area'].isin(country_lookup['Name'])]
un_csv =un_csv[~un_csv['Area'].isin(not_in_afr)]

un_csv['new_pop'] = un_csv['Pop'].str.replace(' ', '')
un_csv['new_pop'] = un_csv['new_pop'].astype(int)


x = np.unique(un_csv['Area'])

un_2000 = un_csv[un_csv['Year'] == 2000]['new_pop']
un_2010 = un_csv[un_csv['Year'] == 2010]['new_pop']


UN = [np.nansum(un_2000) * 1000, np.nansum(un_2010) * 1000]
#from WPP2019_POP_F07_1_Population_By_AGE_BOTH_SEXES

#%%

fig = plt.figure(figsize = (12,14))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)


ax_list = [ax1, ax2]
for ax in ax_list:
    ax.set_xlim([1995, 2105])

#pop
years = np.unique(ssp2['Year'])
years = years[years >= 2000]

ax1.plot(years, tpop_raw, label = 'SSP2', c= '#d95f02')

ax1.plot(years, tpop, label = 'SSP2 Adjusted', c= '#7570b3')

#plt.plot(years, tpop3, label = 'SSP3')
ax1.scatter([2000, 2010], [x2000, x2010], label = 'WorldPop')
ax1.scatter([2010], [gx2010], label = 'GPWv4', c = 'black')
ax1.scatter([2000, 2010], UN, label = 'United Nations', c = 'red')


plt.xlabel('Year')
ax1.set_ylabel('Total under 5 African population')

#pmor
actual_x = [2000, 2010]
actual_y = np.array([mx2000,
                     mx2010])

actual_y = actual_y / 1000000

years = np.arange(2000, 2061)
ax2.plot(years[0:41], np.array(tmor_raw[0:41])/ 1000000, label = 'GBD Ref', c= '#d95f02')

ax2.plot(years[0:41], np.array(tmor[0:41]) / 1000000, label = 'GBD Ref Adjusted', c= '#7570b3')
ax2.plot(years[41:], np.array(tmor_raw[41:])/ 1000000, label = 'Extrapolation', linestyle = '--',
          c= '#d95f02')
ax2.plot(years[41:], np.array(tmor[41:]) / 1000000, label = 'Extrapolation', linestyle = '--',
          c= '#7570b3')


ax1.annotate('A)', (-0.03, 1.01), xycoords = 'axes fraction',
             rotation = 0, ha = 'left', weight = 'bold')
ax2.annotate('B)', (-0.03, 1.01), xycoords = 'axes fraction',
             rotation = 0, ha = 'left', weight = 'bold')

ax2.scatter(actual_x, actual_y, label = 'UNICEF')
ax2.set_xlabel('Year')
ax2.set_ylabel('Total under 5 African mortality')

ax2.annotate(s = '1e6', xy = (0.02, 1.01), xycoords = 'axes fraction',ha = 'center')

ax1.legend(bbox_to_anchor=(1.0, -0.05),ncol=5,frameon = False, handletextpad = 0.5)
ax2.legend(bbox_to_anchor=(0.8, -0.1),ncol=3,frameon = False, handletextpad = 0.5)

#plt.savefig('/nfs/see-fs-02_users/earsch/Documents/Leeds/Inputdata_POP_MOR.png',
#         bbox_inches = 'tight', pad_inches = 0.3)

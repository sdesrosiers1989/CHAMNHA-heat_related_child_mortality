#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Functions for infant mortality model (due to temp)

Created on Mon Apr 20 11:12:29 2020

@author: earsch
"""

#%%set wd and import packages

import iris
import iris.quickplot as qplt
import iris.coord_categorisation
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units

import numpy as np
import numpy.ma as ma

import math


import copy

import glob

#%%

def place_holder(base, dat):
    new_dat = copy.deepcopy(base)
    new_dat.data = np.where(~np.isnan(new_dat.data), dat, new_dat.data)
    new_dat.data = ma.masked_array(new_dat.data, mask = base.data.mask)
    return new_dat

def calc_tdif(tavg, thres):
    
    tdif = iris.analysis.maths.subtract(tavg, thres.data)
    tdif.data = np.ma.where(tdif.data < 0, 0, tdif.data)
    return tdif

def ann_death_per_decade(temp, dec_start, dec_end, pop_ratio, davg_mort, coeff):
    b_dec = temp.extract(iris.Constraint(year= lambda cell: dec_start <= cell < dec_end))
    dims =  b_dec.shape
    
    #cycle through each year in decade
    #for each gridcell, extract daily data
    #from daily tdif (difference between threshold and tavg), calcualte daily att deaths
    #sum daily att deaths to get total annual deaths from heat per year
    #calculate mean annual deaths per year (for the decade of interest)
    
    year_output = iris.cube.CubeList()
    years = np.unique(b_dec.coord('year').points)
    
    e_list = iris.cube.CubeList()
    
    for y in years:
        b_year = b_dec.extract(iris.Constraint(year= lambda cell: cell == y))
        output = copy.deepcopy(b_year)
        e_output = copy.deepcopy(b_year)
        print(y)
        
        
        for j in np.arange(dims[1]):
            for k in np.arange(dims[2]):
                b_day = b_year[:, j, k]
                   
                if np.isnan(b_day[0].data) == False: # if not masked
                    
                    
                    p_b = pop_ratio[j,k].data #ratio of future to baseline pop for gridcell
                    c = coeff[j,k].data # coeff for gridcell
                    m_loc = davg_mort[j, k].data
                    
                    #sensitive to placement of brackets
                    daily_att_deaths = p_b * (m_loc / np.exp(c * b_day.data)) * (np.exp(c * b_day.data)- 1)
                    e = (1 / np.exp(c * b_day.data)) * (np.exp(c * b_day.data)- 1)
                    
                    output.data[:,j,k] = daily_att_deaths
                    e_output.data[:,j,k] = e
        #sum total heat deaths
        out_sum = output.collapsed('time', iris.analysis.SUM)
        year_output.append(out_sum)
        e_mean = e_output.collapsed('time', iris.analysis.SUM)
        e_list.append(e_mean)
        
    #calculate annual heat deaths
    year_merge =  year_output.merge_cube()
    ann_avg_heat_death = year_merge.collapsed('time', iris.analysis.MEAN)
    e_merge = e_list.merge_cube()
    
    return year_merge, ann_avg_heat_death, e_merge



    

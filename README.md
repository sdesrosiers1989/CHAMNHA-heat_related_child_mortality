# CHAMNHA Projection: Heat related child mortality
Scripts for calculating heat-relatated child mortality (HRCM) due to climate change, and changes in population size and total (all-cause mortality), and for decomposing changes in HRCM due to changes in climate, population and all-cause mortality. <br /> 

Children are defined here as under 5, and the scripts and input data are based on that. However, the scripts will work with any age group.

## Inputs
- daily temperature data (CMIP6)
- population data
  - historical: WorldPop under-5s 1km gridded dataset
  - future: SSP2 under-5s country totals 
- total (all-cause) mortality (not as a rate)
  - historical: UNICEF under-5s all-cause mortality, country totals
  - future: GBD Reference scenario under-5s, country totals

## Process

1. Prepare temperature (tas) input data
   - Assuming that CMIP6 data already downloaded from ESGF or CEDA, script will concatenate model files into one file per model and scenario, and ensure all on same grid and cover same time frame.
   - Scripts (in CMIP6_temp_prepare):
     - CMIP6_temp_prepare.py - for SSPs and historical, not including CESM models
     - CMIP6_temp_prepare_DAMIP.py - for DAMIP scenario
     - CMIP6_temp_prepare_CESM.py - CESM models not concatentating properly - fix them in this script
3. Bias correct input tas data (OPTIONAL)
   - Bias-correct temperature data using linear-scaling and monthly correction factors, and the CRU temperature dataset as the reference observations
   - Scripts (in bias_corr)
     - tas - sripts are available for CORDEX, CMIP6 and the CP4A and P25 UM models, however the rest of these scripts are based on using the CMIP6 models (though any netcdf file in same format should work)
4. Preprae input pop/mortality data
   - Health burden model requires gridded total population, and gridded total annual mortality.
   - Compare different historical population and mortality input datasets (comp_pop_datasets) (OPTIONAL)
   - Scripts (in prepare_popmor_data)
     - historical
       - WorldPop data already available gridded - download individual files and converted to netCDF using NCO
       - Using prepare_pop_data.py concatenate into one file and regrid to CMIP6 grid
       - Note that the historical mortality data was available as country totals, but was distribtued to individual gridcells in ArcGIS using the 2019 WorldPop gridded data. If ArcGIS is not available, future scripts could be adapted for the same purpose.
     - future (prepare_popmor_data/future)
       - Population script (pop/prepare_future_data.py):
         - Adjust SSP2 country totals so that historical data in-line with historical WorldPop data (linear-scaling)
         - Distribute adjusted country totals to individual gridcells, based on 2019 population distribution
       - Mortality script (mortality/prepare_future_mor_data.py):
         - Extrapolate from 2040 to 2050 using a 2nd order polynomial
         - Adjust country totals so GBD Reference Scenario historical values in line with UNICEF historical values
         - Distributed adjusted country totals to individual grid cells, based on 2019 population distribution
5. Run health burden model
   - Health burden model based on a linear threshold model, whereby temperature impacts on mortality increases linearly above a certain threshold. Model described in Hajat (2014). In the implementation here, the threshold used is the 75th percentile based on the historical scenario for each individual climate model. Two coefficients are used based on previous studies in Africa (Azongo 2012, Egondi 2012). Note that the relationship between mortality and temperature above the threshold in these papers is presented as percentage. The scripts convert that into the coefficient required by the health burden model. 
   - Adjusting coefficients: historical model has a separate script for each coefficient, however other scenarios coefficient can be adjusted within the #Coeff section of the script (change file paths accordingly)
   - The scripts take in climate data, as well as the population and mortality data. Population and mortality data available per decade, temperature data daily. The health burden model requires daily all-cause mortality, however as this data was not available annual mortality is simply divided by 365.
   - Outputs are the annual heat related mortality (individual years), and as a mean per decade, and * *e*, which is the fraction of total mortality due to heat, and is required in the decomposition scripts
   - Scripts (healthburden_model):
     - historical: 1995 - 2014, historical scenario from CMIP6, and the appropriate population and mortality data. Separate scripts for coeff = 0.61 and 1.0
     - damip: 1995 - 2020, with 1995 - 2014 being the hist-nat scenario, and 2015 - 2020 being ssp245. Separate file for FGOALs model.
     - future: 2020 - 2050, using future population and mortality
7. Decompose into components
   - Decompose components of change in HRCM into change due to climate change (temperature increases), population (population growth) and all-cause mortality (declining)
   - Based on the method of Das Gupta (1993)
   - Note that you need to be quite careful with the inputs to this script, and check that all the components add up to the total change. 
   - Scripts (in calc_components folder): <br />
     - calc_components.py - for attribution experiment. Run this first to get csv files of results, which will be used as inputs in calc_components_cc.py for creating plots
     - calc_components_CC.py - for climate change experiments, but also takes csv file input from calc_components.py to create plots
## Branches
- master: use this one for all analysis
- pop2019:
  - used for checking impact of using 2019 population for 2020 - 2050
- popmor_dif_biascorr: <br />
  - Created when checking the impact of different bias-correction methods for future population and mortality on the HRCM and decomposition results. Impact of different methods minimal, so stayed used original, simpler method.
- mor_rate:
  - Created when checking if future mortality should be based on mortality rate rather than total. Total worked better.

## References

Azongo, D. K., Awine, T., Wak, G., Binka, F. N. & Oduro, A. R. A time series analysis of weather variability and all-cause mortality in the Kasena-Nankana districts of Northern Ghana, 1995-2010. Glob. Health Action 5, 14–22 (2012).  <br />

Das Gupta, P. Decomposition of Rates: A User’s Manual. Stand. Decompos. Rates A User’s Man. P23-186 (1993). <br />

Egondi, T. et al. Time-series analysis of weather and mortality patterns in Nairobi’s informal settlements. Glob. Health Action 5, 23–32 (2012).  <br />

Hajat, S., Vardoulakis, S., Heaviside, C. & Eggen, B. Climate change effects on human health: Projections of temperature-related mortality for the UK during the 2020s, 2050s and 2080s. J. Epidemiol. Community Health 68, 641–648 (2014).  <br />

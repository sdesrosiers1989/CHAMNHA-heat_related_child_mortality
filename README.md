# CHAMNHA-heat_related_child_mortality
Scripts for calculating heat-relatated child mortality (HRCM) due to climate change. Also includes mortality and population change. 
Children are defined here as under 5, and the scripts and input data are based on that. However, the scripts will work with any age group.

## Inputs
- daily temperature data (CMIP6)
- population data
- total (all-cause) mortality (not as a rate)

## Process

1. Prepare temperature (tas) input data
   - CMIP6 data already downloaded from ESGF or CEDA
   - Need to concatenate into one file per model, scenario
   - Scripts (in CMIP6_temp_prepare):
     - CMIP6_temp_prepare.py - for SSPs and historical, not including CESM models
     - CMIP6_temp_prepare_DAMIP.py - for DAMIP scenario
     - CMIP6_temp_prepare_CESM.py - CESM models not concatentating properly - fix them in this script
3. Bias correct input tas data (potentially could be skipped)
   - biascorr/tas/tas_biascorr_CMIP6 (note that scripts are available for tmax, tmin, as well as the CORDEX and CP4A/P25 UM models, however the scripts here all are based on daily tas, and CMIP6 data
   -  
4. Preprae input pop/mortality data
   - point 1
   - Compare different population and mortality input datasets (comp_pop_datasets) - Optional
5. Run health burden model
6. Decompose into components
   - Decompose components of change in HRCM into change due to climate change (temperature increases), population (population growth) and all-cause mortality (declining)
   - Based on the method of Das Gupta (1993)
   - Note that you need to be quite careful with the inputs to this script, and check that all the components add up to the total change. 
   - Scripts (in calc_components folder): <br />
     - calc_components.py - for attribution experiment. Run this first to get csv files of results, which will be used as inputs in calc_components_cc.py for creating plots
     - calc_components_CC.py - for climate change experiments, but also takes csv file input from calc_components.py to create plots
## Branches
- master: use this one for all analysis
- popmor_dif_biascorr: <br />
  - Created when checking the impact of different bias-correction methods for future population and mortality on the HRCM and decomposition results. Impact of different methods minimal, so stayed used original, simpler method.
- mor_rate:
  - Created when checking if future mortality should be based on mortality rate rather than total. Total worked better.

## References

Das Gupta, P. Decomposition of Rates: A User’s Manual. Stand. Decompos. Rates A User’s Man. P23-186 (1993). <br />

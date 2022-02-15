# CHAMNHA-heat_related_child_mortality
Scripts for calculating heat-relatated child mortality due to climate change. Also includes mortality and population change. 
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
4. Next step

## References

# rgsm_integrated

This project contains data and code required to reproduce analysis and create figures in a manuscript descibing an integrated population model for rio grande silvery minnow that has been submitted to Ecological Applications.

Below are brief descriptions of the code and data in this project

## R and Stan code required to reproduce results and figures in manuscript
code_together.r – R code that uses data to fit models based on Stan code and produce Figures 2 – 9 and S1 – S4 in the manuscript
int3.stan – Stan code for fitting integrated model with constant mortality rates
int3re.stan – Stan code for fitting integrated model with time varying mortality rates

##CSV files containing information required to reproduce results and figures in manuscript
abq_gage_08330000.csv – daily flow data for USGS gage 08330000 in long format
aQ.csv – daily flow data for USGS gage 08330000 from March 1st to September 30th organized by year
ASIRQuery4Combined.csv – RGSM monitoring dataset - subset for analysis, updated version available online
Fish Rescue 2009-2020.csv – Fish rescue data
fws2br_rmconverter.csv - file for converting between river mile system used by USFWS and USBOR
mesohab_sum.csv – summary of habitat measurements from Braun et al., 2015
oos_catch.csv – catch data for October 2019 and October 2020.
oos_releases.csv – information on augmentation relevant to out of sample analysis
oos_rivereyes.csv – information on river drying during 2019 and 2020
releases 8_26_2019.csv – information on augmentation from 2002 to 2018
rescue_surv.csv – information on survival of rescued fish from Archdeacon et al., 2020
rivereyes_v2.csv - information on river drying from 2002 to 2018
SanAcacia_gage_0835490.csv 0 – daily flow data for USGS gage 08354900 in long format
sQ.csv – daily flow data for USGS gage 08354900from March 1st to September 30th organized by year
sum_ee1.csv - Expert elicitation results from expert 1
sum_ee2.csv - Expert elicitation results from expert 2
sum_ee3.csv - Expert elicitation results from expert 3
sum_ee4.csv - Expert elicitation results from expert 4
sum_ee5.csv - Expert elicitation results from expert 5

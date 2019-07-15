%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Memo for MATLAB version of the model   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1) The function to run the model can be found in the file run_sig_iso_dtchange.m and can be executed via:

run_sig_iso_dtchange(cc_rain_flx_in, rainratio_in, dep_in, dt_in, oxonly_in, folder) 

and example call would be: run_sig_iso_dtchange(6.0e-5, 1.5, 0.24, 1d8, true, './1207_test')


2) To run the lysocline experiments (i.e., Section 3.1 in the manuscript) execuet the functions in the file lysocline_exp.m


3) Subroutines of the main code can be found in caco3_main.m


4) CaCO3 therdomynamic subroutines & functions are in caco3_therm.m 


5) Test functions used during model development can be found in caco3_test.m 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Boundary/initial conditions and  model options

Boundary/initial conditions can be changed in caco3_main.caco3_set_boundary_cond(), such as bottom water concentrations/fluxes

Global properties such as depth of simulated sediment column, grid numbers, densities, threshold values, rate constants and model options (e.g. OM degradation method, type of mixing to be used, enable signal tracking?)
are defined under properties of the class caco3_main. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plotting of the results is cuurently possible with the following python scripts (folder of model results and number of CaCO3 classes has to be specified in these python scripts):

1) matlab_caco3_profiles_sum_multi_v3.py	:	plots the geochemical profiles, e.g. as in Fig. 3 of the masnuscript

2) matlab_caco3_lys_oxanox_sum_v3.py		:	plots the lysocline results for oxic-only and oxic-anoxic OM degradation model, e.g. as in Figs. 5+6 of the masnuscript

3) matlab_caco3_signals.py			:	plots the time change of the proxy signals, e.g. as in Fig. 8 of the masnuscript



Memo for MATLAB version of the model

CaCO3 therdomynamic subroutines & functions in caco3_therm.m 
Main code caco3_main.m
Test functions in caco3_test.m 

Need BLAS libraries 

To run the model without signal tracking, e.g.:
caco3_test.chk_through(2.0, 35.0, 4.0)


To run the model with signal tracking, e.g.:
caco3_test.chk_through_signal(2.0, 35.0, 4.0, 5.0)

Change the boundary conditions in caco3_test.caco3_set_boundary_cond()

Change main properties (e.g. sediment characteristics, constants, preprocessor variables as in defines.h) in properties of class caco3_main.m

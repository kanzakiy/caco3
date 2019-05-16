Memo for MATLAB version of the model

CaCO3 therdomynamic subroutines & functions in caco3_therm.m 
Main code caco3_main.m
Test functions in caco3_test.m 

Need BLAS libraries - follow instructions on https://www.mathworks.com/help/matlab/matlab_external/calling-lapack-and-blas-functions-from-mex-files.html 
to build and copy the MEX file matrixDivide.c to your local matlab working-directory:

copyfile(fullfile(matlabroot,'extern','examples','refbook','matrixDivide.c'),'.')
fileattrib('matrixDivide.c','+w')
mex -v matrixDivide.c -lmwlapack

then can be used as:
X = matrixDivide(A,B)


To run the model without signal tracking, e.g.:
caco3_test.chk_through(2.0, 35.0, 4.0)


To run the model with signal tracking, e.g.:
caco3_test.chk_through_signal(2.0, 35.0, 4.0, 5.0)

Change the boundary conditions in caco3_test.caco3_set_boundary_cond()

Change main properties (e.g. sediment characteristics, constants, preprocessor variables as in defines.h) in properties of class caco3_main.m

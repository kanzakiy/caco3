Signal tracking diagenesis by Python2.7 

This is another way to run the model with Python
    by creating a python module from Fortran code through f2py* and just importing and using the module
* You need numpy 

1. Create caco3mod.dll
    (a) Compile caco3_therm.f90 and caco3_test_mod_v5_6.f90 as in readme file for Fortran (see readme.txt)
    (b) Type 'python -m numpy.f2py -c -m caco3mod caco3_python.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas'  
        (or 'python -m numpy.f2py -c -m caco3mod caco3_python.f90 caco3_test_mod_v5_6.o caco3_therm.o umf4_f77wrapper.o
            -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -llapack -lopenblas' 
            if you choose to use sparse matrix solver)
2. Import and use the module
    (a) Change the variables in caco3_fortran.py (rain flux of caco3, om/caco3 rain ratio, file name etc.)
        (you can change other variables in define.h before compiling fortran codes) 
    (b) Type 'python caco3_fortran.py'  
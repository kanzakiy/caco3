Memo

CaCO3 therdomynamic subroutines & functions in caco3_therm.f90 
Main code caco3_test_mod_v5_5.f90
caco3_test_mod_v5_6.f90 works as well consisting of lots of subroutines 

Need BLAS (& UMFPACK if the number of CaCO3 species is large) libraries 

Following http://www.hnagata.net/archives/212
And http://www5.hp-ez.com/hp/calculations/page17
http://qiita.com/AnchorBlues/items/69c1744de818b5e045ab
   OpenBLAS (+lapack)
1. download: http://www.openblas.net/ -> TAR
2. tar zxvf OpenBLAS-0.2.20.tar.gz
3. cd OpenBLAS-0.2.20
4. make BINARY=64 CC="gcc -m64" FC="gfortran -m64"
5. su
6. make PREFIX=/usr/local install

If you get the error 'libopenblas.so.0: cannot open shared object file: No such file or directory' then type
1) sudo apt-get install libopenblas-base
2) export LD_LIBRARY_PATH=/usr/lib/openblas-base/

See https://github.com/PetterS/SuiteSparse/tree/master/UMFPACK for UMFPACK

Note that UMFPACK is usually not necessary. 

Simulation can be run by following steps:
(a) specify directory where results are stored in caco3_test_mod_v5_5.f90 (lines 235 and 2650) 
(b) comment/comment out macros in defines.h 
(c) compile by typing: 
    [if you do not have UMFPACK]
        a) gfortran -c caco3_therm.f90
        b) gfortran -c -cpp -I/path/to/working/directory caco3_test_mod_v5_6.f90
        c) gfortran caco3_fortran.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas -g -fcheck=all
    [if you have UMFPACK] 
        a) gfortran -c caco3_therm.f90
        b) gfortran -c -cpp -I/path/to/working/directory caco3_test_mod_v5_6.f90
        c) gfortran caco3_fortran.f90 caco3_test_mod_v5_6.o caco3_therm.o umf4_f77wrapper.o -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lopenblas -g -fcheck=all
(d) run by './a.exe cc caco3_rainflux_value rr om/caco3_rain_ratio_value dep water_depth_value dt time_step_value fl simulation_name'
    where caco3_rainflux_value [umol cm2 yr-1], om/caco3_rain_ratio_value, water_depth_value [km], time_step_value [yr] and simulation_name are your inputs. 
    The water_depth_value [km] represents water depth when simulation does not track proxy signals ('sense' macro is defined in defines.h), and 
    maximum water depth when proxy signal is tracked (i.e., when 'sense' is not defined in defines.h). 
    
In step (b) above, you need to define macros in defines.h to simulate what you want. 
E.g.,
(1) With defining 'sense' with/without 'oxonly', you can run the model to predict lysoclines and caco3 burial fluxes.  
(2) With defining 'biotest', you can test bioturbation effect in 5 kyr signal tracking. 
(3) With defining 'track2', you can test different signal tracking method.
(4) With defining 'size', you can track two different types of caco3 species. 
(5) You need switch on 'sparse' to use sparse matrix solver (you need UMFPACK in this case). 

Default run can be made with switching off all options (commenting out all options in defines.h), 
which track 2 isotope signals with 4 caco3 species for 50 kyr isotope shift event assuming Fickian mixing for bioturbation. 
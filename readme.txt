Memo

CaCO3 therdomynamic subroutines & functions in caco3_therm.f90 
Main code caco3_test_mod_v5_5.f90

Need BLAS (& UMFPACK if the number of CaCO3 species is large) libraries 

Following http://www.hnagata.net/archives/212
And http://www5.hp-ez.com/hp/calculations/page17
http://qiita.com/AnchorBlues/items/69c1744de818b5e045ab
[] OpenBLAS (+lapack)
1. download: http://www.openblas.net/ -> TAR
2. tar zxvf OpenBLAS-0.2.20.tar.gz
3. cd OpenBLAS-0.2.20
4. make BINARY=64 CC="gcc -m64" FC="gfortran -m64"
5. su
6. make PREFIX=/usr/local install

See https://github.com/PetterS/SuiteSparse/tree/master/UMFPACK for UMFPACK
program chk_poro
! just check the calculation of transition matrix using subroutines in caco3_test_mod_v5_6.f90

! (1) gfortran -c caco3_therm.f90
! (2) gfortran -c -cpp -I/path/to/working/directory caco3_test_mod_v5_6.f90 
! (3) gfortran chk_poro.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas -g -fcheck=all
! (4) ./a.out

! you can check porosity as function of depth  
! please copy and paste results to whatever file to be compared with results with MATLAB version 
use globalvariables
implicit none 

beta = 1.00000000005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
call makegrid(beta,nz,ztot,dz,z)

call getporosity() ! assume porosity profile 

! showing porosity profile as function of depth
! porosity              ---> poro
! volume solid fraction ---> sporo = 1 - poro
! z                     ---> depth 
do iz=1,nz
    print*,z(iz),poro(iz),sporo(iz)
enddo

endprogram 
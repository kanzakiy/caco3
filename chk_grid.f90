program chk_grid
! just check the calculation of sediment grid cells caco3_test_mod_v5_6.f90

! to run type 
! (1) gfortran -c caco3_therm.f90
! (2) gfortran -c -cpp -I/path/to/working/directory_test_mod_v5_6.f90 
! (3) gfortran chk_grid.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas -g -fcheck=all
! (4) ./a.out

! you can check grid No. vs. dz and z (cm) 
! please copy and paste results to whatever file to be compared with results with MATLAB version 

implicit none 
integer(kind=4),parameter :: nz=100
real(kind=8) beta,dz(nz),z(nz),ztot 
integer(kind=4) iz

ztot=500d0
beta = 1.00000000005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
call makegrid(beta,nz,ztot,dz,z)

do iz=1,nz
    print*,iz,dz(iz),z(iz)
enddo

endprogram 
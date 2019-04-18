program chk_caco3_therm_sbrtns
! just check the calculation of aqueous co2 species (and its derivatives) using subroutines in caco3_therm.f90

! (1) gfortran -c caco3_therm.f90
! (2) gfortran -c -cpp -I/path/to/working/directory caco3_test_mod_v5_6.f90 
! (3) gfortran chk_therm_sb.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas -g -fcheck=all
! (4) ./a.out

! you can check dic,alk,co2,hco3,co3,ph,dco3_ddic,dco3_dalk   
! please copy and paste results to whatever file to be compared with results with MATLAB version 

implicit none 
integer(kind=4),parameter :: nz = 100
real(kind=8) dep,temp,sal,dic(nz),alk(nz),ph(nz),co2(nz),hco3(nz),co3(nz),dco3_dalk(nz),dco3_ddic(nz)
integer(kind=4) iz,info

temp = 2d0 ! temperature in Celsius 
sal = 35d0 ! salinity wt o/oo 
dep = 4.0d0 ! depth in km
dic = 2211d0*1d-6/1d3  ! 2211 uM converted to mol/cm3; all the same value at different grids
alk = 2285d0*1d-6/1d3
! calling subroutine to calculate all aqueous co2 species and pH
call calcspecies(dic,alk,temp,sal,dep,ph,co2,hco3,co3,nz,info)
if (info/=0) then ! if ph cannot be calculated in the subroutine, info=1 is returned 
    print*, 'error in calcspecies'
    stop
endif 
! calling subroutine to calculate derivatives of co3 conc. wrt dic and alk (defined as dco3_ddic and dco3_dalk, respectively)
call calcdevs(dic,alk,temp,sal,dep,nz,info,dco3_dalk,dco3_ddic)
if (info/=0) then ! if ph cannot be calculated in the subroutine, info=1 is returned 
    print*, 'error in calcdevs'
    stop
endif 
! printing results on screen; if you want check with MATLAB version by copy and paste
do iz=1,nz
    print*,dic(iz),alk(iz),co2(iz),hco3(iz),co3(iz),ph(iz),dco3_ddic(iz),dco3_dalk(iz)
    ! note that the above concentrations (alk,dic,co2,hco3,co3) are all in mol/cm3
    ! ph is actually H+ concentration in mol/L 
    ! dco3_dalk and dco3_ddic should be dimensionless
enddo

endprogram 
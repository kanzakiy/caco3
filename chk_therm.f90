program chk_caco3_therm
! just check the calcite and aragonite solubility

! to run type 'gfortran -c caco3_therm.f90; gfortran chk_therm.f90 caco3_therm.f90;./a.exe' 
! or 'gfortran -c caco3_therm.f90; gfortran chk_therm.f90 caco3_therm.f90;./a.out'

! you can check depth vs. equillibrium consts (mol^2 kg^-2) for calcite and aragonite on screen
! please copy and paste results to whatever file to be compared with results with MATLAB version 

implicit none 
real(kind=8) dep,temp,sal,keqcc,keqag,calceqcc,calceqag
integer(kind=4) i,n

temp = 2d0
sal = 35d0
print*,'dep, keqcc, keqag'
n = 100
do i=1,n
    dep=6d0*i/real(n,kind=8)
    keqcc=calceqcc(temp,sal,dep)
    keqag=calceqag(temp,sal,dep)
    print*,dep,keqcc,keqag
enddo

endprogram 
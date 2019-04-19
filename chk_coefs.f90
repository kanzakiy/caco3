program chk_coefs
! just check the calculation of kinetic and thermodynamic coefficients using subroutines in caco3_test_mod_v5_6.f90

! (1) gfortran -c caco3_therm.f90
! (2) gfortran -c -cpp -I/path/to/working/directory caco3_test_mod_v5_6.f90 
! (3) gfortran -cpp -I/path/to/working/directory chk_coefs.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas -g -fcheck=all
! (4) ./a.out

! <<<< Maybe first only consider Fickian mixing so switch off all macros in defines.h (you can switch on 'test') >>>>>

! you can check rate and thermodynamic coefficients  
! please copy and paste results to whatever file to be compared with results with MATLAB version 

use globalvariables
implicit none 
#include <defines.h>

#ifdef allnobio 
nobio = .true.
#elif defined allturbo2 
turbo2 = .true.
#elif defined alllabs 
labs = .true.
#endif 

#ifdef oxonly
anoxic = .false. 
#endif

beta = 1.00000000005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
call makegrid(beta,nz,ztot,dz,z)

call getporosity() ! assume porosity profile 

call coefs(temp,sal,dep)  ! determining diffusion coeffs, rate constants for om decomposition and cc dissolution
                           ! , dissociation constants of co2 and hco3, and co3 conc. at calcite saturation 
                            ! Note that porosity is needed for diffusion coefficient to reflect tortuosity

! showing thermodynamic consts which are not functions of depth
print*,'keq1    ,keq2   ,keqcc  ,co3sat '
print*,keq1,keq2,keqcc,co3sat ! dissociation const. of h2co3, dissociation const. of hco3, solubility of calcite, and co3 conc. at calcite saturation 
print*,''
print*,''
print*,''

print*,'z   ,dif_dic    ,dif_alk    ,dif_o2 ,kom'
! showing coeffs against depth on screen 
do iz=1,nz 
    print*,z(iz),dif_dic(iz),dif_alk(iz),dif_o2(iz),kom(iz) 
    ! depth [cm], dic, alk and o2 diffusion coefficients [cm2 yr-1], om decompposition rate const. [ yr-1] 
enddo
print*,''
print*,''
print*,''

do iz=1,nz  ! dissolution rate constant are species specific so recording on txt files 
    if (iz==1) print*,'z', (isp ,isp=1,nspcc)
    print*,z(iz),(kcc(iz,isp),isp=1,nspcc)
enddo

endprogram 
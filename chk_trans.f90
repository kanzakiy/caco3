program chk_trans
! just check the calculation of transition matrix using subroutines in caco3_test_mod_v5_6.f90

! (1) gfortran -c caco3_therm.f90
! (2) gfortran -c -cpp -I/path/to/working/directory caco3_test_mod_v5_6.f90 
! (3) gfortran -cpp -I/path/to/working/directory chk_trans.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas -g -fcheck=all
! (4) ./a.out

! you can check transition matrix for individual species 
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

! call getporosity() ! assume porosity profile 
call getporosity(  &
     poro,porof,sporo,sporof,sporoi & ! output
     ,z,nz  & ! input
     )

! call make_transmx()
call make_transmx(  &
    trans,izrec,izrec2,izml,nonlocal  & ! output 
    ,labs,nspcc,turbo2,nobio,dz,sporo,nz,z  & ! input
    )

! recording transition matrices for individual species (isp=1 --> om, 2-->clay,3~2+nspcc --> caco3 species)
! in default, transition matrices are the same for all the species 
! file will be created your working directory of a name 'chk_trans_sp-xxx.txt' where xxx denotes species number 
! sizes of matrices are all (nz, nz )
do isp=1,nspcc+2
    write(dumchr(1),'(i3.3)') isp 
    open(unit=file_tmp,file='./chk_trans_sp-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace') 
    do iz=1,nz
        write(file_tmp,*)(trans(iz,iiz,isp),iiz=1,nz)
    enddo
    close(file_tmp)
enddo

endprogram 
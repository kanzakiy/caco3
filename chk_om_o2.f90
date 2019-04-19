program chk_om_o2
! just check the calculation of conc. and flx of om & o2 using subroutines in caco3_test_mod_v5_6.f90

! (1) gfortran -c caco3_therm.f90
! (2) gfortran -c -cpp -I/path/to/working/directory caco3_test_mod_v5_6.f90 
! (3) gfortran -cpp -I/path/to/working/directory chk_om_o2.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas -g -fcheck=all
! (4) ./a.out

! <<<< Maybe first only consider Fickian mixing so switch off all macros in defines.h (you can switch on 'test') >>>>>

! you can check calculation of concs. and flxes of om & o2 
! please copy and paste results to whatever file to be compared with results with MATLAB version 
! NOTE: This tests o2 calculation, but o2 calculation should be boring without om, so the code here does iterations at individual time steps as in the whole code.  

#include <defines.h>
use globalvariables
implicit none 
integer(kind=4) interval  ! choose value between 1 to nz 

interval =10 ! choose a value between 1 to nz; om depth profile is shown with this interval; e.g., if inteval = nz, om conc. at all depths are shown
! e.g., if interval = 5, om conc. at 5 depths are shown   

write(dumchr(2),'(i0)') interval
dumchr(1)="(A,"//trim(adjustl(dumchr(2)))//"E11.3"//")"

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

call make_transmx()

call coefs(temp,sal,dep)  ! need to specify diffusion coefficient as well as om decomposition rate const. etc.

om = 1d-8  ! assume an arbitrary low conc. 
o2 = o2i*1d-6/1d3 ! o2 conc. in uM converted to mol/cm3

omx = om
o2x = o2

time = 0d0 ! model time [yr]
it = 1 ! integration count 
nt = 10 ! total integration 
dt = 100d0 ! time step [yr]

rho = 2.5d0 ! assume here 

!!!  addition to chk_om.f90 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
zox = 10d0  ! initial assumption on oxygen penetaration depth [cm]
!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do it=1,nt
    print'(A,i0,A,E11.3,A,E11.3,A)','(it,dt,time)  (',it,',',dt,',',time,')'
    
    itr = 0  ! iteration number for om and o2 calcuation 
    error = 1d4 ! error in ieration for zox 
    minerr= 1d4  ! recording minimum relative difference in zox from previously considered zox 
    
    do while (error > tol)
    
        call omcalc() ! om conc. calculation 
        ! calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations )
        call calcflxom()
        
        ! print*,'~~~~ conc ~~~~'
        ! print dumchr(1), 'z  :',(z(iz),iz=1,nz,nz/interval)
        ! print dumchr(1), 'om :',(omx(iz)*mom/rho(iz)*100d0,iz=1,nz,nz/interval)
        ! print*,'++++ flx ++++'
        ! print'(7A11)', 'tflx','adv','dif','omrxn','ccrxn','rain','res'
        ! print'(A,7E11.3)', 'om :', omtflx, omadv,  omdif, omdec,0d0,omrain, omres
        
        ! print*,'izox',izox
        ! sb omcalc calculates izox, which is the deepest grid where o2 >=0. 
        
        if (izox == nz) then ! fully oxic; lower boundary condition ---> no diffusive out flow  
            call o2calc_ox()  ! o2 calculation when o2 penetration depth (zox) is the same as bottom depth. 
            call calcflxo2_ox() !  fluxes relevant to o2 (at the same time checking the satisfaction of difference equations) 
        else  !! if oxygen is depleted within calculation domain, lower boundary changes to zero concs.
            call o2calc_sbox() ! o2 calculation when o2 is depleted within the calculation domain.
            call calcflxo2_sbox() ! fluxes relevant to oxygen 
        endif
        
        ! print*,'~~~~ conc ~~~~'
        ! print dumchr(1), 'z  :',(z(iz),iz=1,nz,nz/interval)
        ! print dumchr(1), 'o2 :',(o2x(iz)*1d3,iz=1,nz,nz/interval)  ! o2 in mol/L
        ! print*,'++++ flx ++++'
        ! print'(7A11)', 'tflx','adv','dif','omrxn','ccrxn','rain','res'
        ! print'(A,7E11.3)', 'o2 :',o2tflx,0d0, o2dif,o2dec, 0d0,0d0,o2res

        ! update of zox 
        zoxx = 0d0
        do iz=1,nz
            if (o2x(iz)<=0d0) exit
        enddo

        if (iz==nz+1) then ! oxygen never gets less than 0 
            zoxx = ztot ! zox is the bottom depth 
        else if (iz==1) then ! calculating zox interpolating at z=0 with SWI conc. and at z=z(iz) with conc. o2x(iz)
            zoxx = (z(iz)*o2i*1d-6/1d3 + 0d0*abs(o2x(iz)))/(o2i*1d-6/1d3+abs(o2x(iz)))
        else     ! calculating zox interpolating at z=z(iz-1) with o2x(iz-1) and at z=z(iz) with conc. o2x(iz)
            zoxx = (z(iz)*o2x(iz-1) + z(iz-1)*abs(o2x(iz)))/(o2x(iz-1)+abs(o2x(iz)))
        endif
        
        ! error evaluation as relative difference of zox
        error = abs((zox -zoxx)/zox)   
        
        ! print*, 'itr,zox, zoxx, error',itr,zox, zoxx, error
        ! print*,'~~~~~~~~~~~////~~~~~~~~~~~~~'
        
        if (zox==zoxx) exit 
         
        zox = 0.5d0*(zox + zoxx)  ! new zox 
        
        ! if iteration reaches 100, error in zox is tested assuming individual grid depths as zox and find where error gets minimized 
        if (itr>=100 .and. itr <= nz+99) then 
            zox = z(itr-99) ! zox value in next test 
            if (minerr >=error ) then ! if this time error is less than last adopt as optimum 
                if (itr/=100) then 
                    izox_minerr = itr -100
                    minerr = error 
                endif 
            endif
        elseif (itr ==nz+100) then ! check last test z(nz)
            if (minerr >=error ) then 
                izox_minerr = itr -100
                minerr = error 
            endif
            zox = z(izox_minerr)  ! determine next test which should be most optimum 
        elseif (itr ==nz+101) then  ! results should be optimum and thus exit 
            exit
        endif 

        if (itr >nz+101) then 
            stop
        endif

        itr = itr + 1
    enddo
    
    ! showing results on screen
    print*,'~~~~ conc ~~~~'
    print dumchr(1), 'z  :',(z(iz),iz=1,nz,nz/interval)
    print dumchr(1), 'om :',(omx(iz)*mom/rho(iz)*100d0,iz=1,nz,nz/interval)
    print dumchr(1), 'o2 :',(o2x(iz)*1d3,iz=1,nz,nz/interval)  ! o2 in mol/L
    print*,'++++ flx ++++'
    print'(7A11)', 'tflx','adv','dif','omrxn','ccrxn','rain','res'
    print'(A,7E11.3)', 'om :', omtflx, omadv,  omdif, omdec,0d0,omrain, omres
    print'(A,7E11.3)', 'o2 :',o2tflx,0d0, o2dif,o2dec, 0d0,0d0,o2res
    
    print*,''
    print*,''
    print*,''
    
    om = omx
    o2 = o2x
    time = time +dt
    
enddo

endprogram 
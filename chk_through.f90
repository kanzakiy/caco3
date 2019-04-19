program chk_through
! just check the calculation of conc. and flx of om, o2, cc & clay using subroutines in caco3_test_mod_v5_6.f90
! burial rate is updated at a give time step until it converges well 

! (1) gfortran -c caco3_therm.f90
! (2) gfortran -c -cpp -I/path/to/working/directory caco3_test_mod_v5_6.f90 
! (3) gfortran -cpp -I/path/to/working/directory chk_through.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas -g -fcheck=all
! (4) ./a.out

! <<<< Maybe first only consider Fickian mixing so switch off all macros in defines.h (you can switch on 'test') >>>>>

! you can check calculation of concs. and flxes of om, o2, cc & clay, with burial modified at individual time steps 
! please copy and paste results to whatever file to be compared with results with MATLAB version 
! NOTE: Burial is not modified according to reactions and non-local mixing so that solid conc. may go above 100 wt %.  

#include <defines.h>
use globalvariables
implicit none 
integer(kind=4) interval  ! choose value between 1 to nz 

interval =10 ! choose a value between 1 to nz; om depth profile is shown with this interval; e.g., if inteval = nz, om conc. at all depths are shown
! e.g., if interval = 5, om conc. at 5 depths are shown   

workdir = './' ! working directory

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

dep = 4.0d0 ! km water depth; note that temperature and salinitiy has initially assumed values in globalvariables.mod  

beta = 1.00000000005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
call makegrid(beta,nz,ztot,dz,z)

call getporosity() ! assume porosity profile 

!!!!!!!!!!!!! flx assignement and initial guess for burial rate !!!!!!!!!!!!!!!!!!!!!!
call flxstat()  ! assume fluxes of om, cc and clay, required to calculate burial velocity
! print*,om2cc,ccflxi,detflx,omflx,sum(ccflx)
! molar volume (cm3 mol-1) needed for burial rate calculation 
mvom = mom/rhoom  ! om
mvsed = msed/rhosed ! clay 
mvcc = mcc/rhocc ! caco3
call burial_pre() ! initial guess of burial profile, requiring porosity profile  
call dep2age() ! depth -age conversion 
call calcupwindscheme() ! determine factors for upwind scheme to represent burial advection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call make_transmx()

call coefs(temp,sal,dep)  ! need to specify diffusion coefficient as well as om decomposition rate const. etc.

om = 1d-8  ! assume an arbitrary low conc. 
o2 = o2i*1d-6/1d3 ! o2 conc. in uM converted to mol/cm3
cc = 1d-8   ! assume an arbitrary low conc. 
dic = dici*1d-6/1d3 ! mol/cm3; factor is added to change uM to mol/cm3 
alk = alki*1d-6/1d3 ! mol/cm3
pt = 1d-8  ! assume an arbitrary low conc.

call calcspecies(dic,alk,temp,sal,dep,pro,co2,hco3,co3,nz,infosbr)   ! calculation of individual co2 species at initial conditions

omx = om
o2x = o2
ccx = cc
dicx = dic
alkx = alk 
! this may not be necessary as these individual species assume equilibrium 
co2x = co2
hco3x = hco3
co3x = co3
ptx = pt

time = 0d0 ! model time [yr]
it = 1 ! integration count 
nt = 10 ! total integration 
dt = 100d0 ! time step [yr]

rho = 2.5d0 ! assume here 

oxco2 = 0d0  ! oxic degradation of om; here assumed 0 at all time and depth  
anco2 = 0d0  ! anoxic degradation of om; here assumed 0 at all time and depth  
!!!  addition to chk_om.f90 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
zox = 10d0  ! initial assumption on oxygen penetaration depth [cm]
!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do it=1,nt
    
    itr_w = 0  ! # of iteration made to converge w 
    err_w_min = 1d4 ! minimum relative difference of w compared to the previous w 
    err_f = 0d0  !! relative different of total vol. fraction of solids wrt the previous value 
    
300 continue ! point of restart when burial velocity does not converge
    
    print'(A,i0,A,E11.3,A,E11.3,A)','(it,dt,time)  (',it,',',dt,',',time,')'
        
    dw = 0d0 ! change in burial rate caused by reaction and non-local mixing 
    
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
    
    
    !~~  OM & O2 calculation END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! calculation of oxic and anoxic degradation of om (oxco2 and anco2, respectively)
    do iz = 1,nz 
        if (o2x(iz) > o2th) then
            oxco2(iz) = (1d0-poro(iz))*kom(iz)*omx(iz)  ! aerobic respiration 
        else 
            ! o2x(iz) = o2th
            if (anoxic) then 
                anco2(iz) = (1d0-poro(iz))*kom(iz)*omx(iz)  ! anaerobic respiration 
            endif
        endif
    enddo

    do iz=1,nz
        dw(iz) = dw(iz) -(1d0-poro(iz))*mvom*kom(iz)*omx(iz)  !! burial rate change need reflect volume change caused by chemical reactions 
        ! as well as non-local mixing 
        if (turbo2(1).or.labs(1)) then 
            do iiz = 1, nz
                if (trans(iiz,iz,1)==0d0) cycle
                dw(iz) = dw(iz) - mvom*(-trans(iiz,iz,1)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*omx(iiz))
            enddo
        else 
            if (nonlocal(1)) then 
                do iiz = 1, nz
                    if (trans(iiz,iz,1)==0d0) cycle
                    dw(iz) = dw(iz) - mvom*(-trans(iiz,iz,1)/dz(iz)*omx(iiz))
                enddo
            endif
        endif
    enddo

    do iz=1,nz
        if (omx(iz)<omx_th) omx(iz)=omx_th  !! truncated at minimum value 
    enddo
    
    ! calculation of caco3 system 
    call calccaco3sys()
    if (flg_500) then
        print*,'error after calccaco3sys'
        stop
    endif 
    
    call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,infosbr)
    if (infosbr==1) then 
        print*,'error after calcspecies after calccaco3sys'
        stop
    endif 
    
    ! calculation of fluxes relevant to caco3 and co2 system
    call calcflxcaco3sys()
    ! end of caco3 calculations 
    
    ! clay calculation 
    call claycalc()
    call calcflxclay()
    ! end of clay calculation 
    
    ! checking for total volume of solids, density and burial velocity 

    call getsldprop() ! get solid property, rho (density) and frt (total vol.frac)

    err_f = maxval(abs(frt - 1d0))  ! new error in total vol. fraction (must be 1 in theory) 
    !! ========= calculation of burial velocity =============================

    wx = w  ! recording previous burial velocity 

    call burialcalc() ! get new burial velocity

    !  determine calculation scheme for advection 
    call calcupwindscheme()

    ! error and iteration evaluation 
    itr_w = itr_w + 1  ! counting iteration for w 
    err_w = maxval(abs((w-wx)/wx))  ! relative difference of w 
    if (err_w<err_w_min) then 
        err_w_min= err_w  ! recording minimum relative difference of  w 
        wxx = wx  ! recording w which minimizes deviation of total sld fraction from 1 
    endif 
    if (itr_w>100) then   ! if iteration gets too many (100), force to end with optimum w where error is minimum
        if (itr_w==101) then 
            w = wxx   
            go to 300
        elseif (itr_w==102) then 
            w = wxx
            print*, 'not converging w',time, err_w, err_w_min
            pause
            go to 400
        endif 
    endif
    if (err_w > tol) go to 300

400 continue
    
    ! showing results on screen
    write(dumchr(2),'(i0)') interval
    dumchr(1)="(A,"//trim(adjustl(dumchr(2)))//"E11.3"//")"
    print*,'~~~~ conc ~~~~'
    print dumchr(1), 'z  :',(z(iz),iz=1,nz,nz/interval)
    print dumchr(1), 'om :',(omx(iz)*mom/rho(iz)*100d0,iz=1,nz,nz/interval)
    print dumchr(1), 'o2 :',(o2x(iz)*1d3,iz=1,nz,nz/interval)  ! o2 in mol/L
    print dumchr(1), 'cc :',(sum(ccx(iz,:)*mcc)/rho(iz)*100d0,iz=1,nz,nz/interval)
    print dumchr(1), 'dic:',(dicx(iz)*1d3,iz=1,nz,nz/interval)
    print dumchr(1), 'alk:',(alkx(iz)*1d3,iz=1,nz,nz/interval)
    print dumchr(1), 'sed:',(ptx(iz)*msed/rho(iz)*100d0,iz=1,nz,nz/interval)
    print*, '   ..... multiple cc species ..... '
    write(dumchr(2),'(i0)') interval
    dumchr(1)="(i0.3,':',"//trim(adjustl(dumchr(2)))//"E11.3"//")"
    do isp=1,nspcc 
        print dumchr(1),isp,(ccx(iz,isp)*mcc/rho(iz)*100d0,iz=1,nz,nz/interval)
    enddo
    print*,'++++ flx ++++'
    print'(7A11)', 'tflx','adv','dif','omrxn','ccrxn','rain','res'
    print'(A,7E11.3)', 'om :', omtflx, omadv,  omdif, omdec,0d0,omrain, omres
    print'(A,7E11.3)', 'o2 :',o2tflx,0d0, o2dif,o2dec, 0d0,0d0,o2res
    print'(A,7E11.3)', 'cc :',sum(cctflx),  sum(ccadv), sum(ccdif),0d0,sum(ccdis), sum(ccrain), sum(ccres) 
    print'(A,7E11.3)', 'dic:',dictflx, 0d0,dicdif, dicdec,  dicdis, 0d0,dicres 
    print'(A,7E11.3)', 'alk:',alktflx, 0d0, alkdif, alkdec, alkdis, 0d0, alkres 
    print'(A,7E11.3)', 'sed:',pttflx, ptadv,ptdif,  0d0, 0d0, ptrain, ptres
    print*, '   ..... multiple cc species ..... '
    do isp=1,nspcc 
        print'(i0.3,":",7E11.3)',isp,cctflx(isp), ccadv(isp), ccdif(isp),0d0,ccdis(isp), ccrain(isp), ccres(isp) 
    enddo
    write(dumchr(2),'(i0)') interval
    dumchr(1)="(A,"//trim(adjustl(dumchr(2)))//"E11.3"//")"
    print*,'==== burial etc ===='
    print dumchr(1), 'z  :',(z(iz),iz=1,nz,nz/interval)
    print dumchr(1), 'w  :',(w(iz),iz=1,nz,nz/interval)
    print dumchr(1), 'rho:',(rho(iz),iz=1,nz,nz/interval)
    print dumchr(1), 'frc:',(frt(iz),iz=1,nz,nz/interval)
    
    print*,''
    print*,''
    print*,''
    
    om = omx
    o2 = o2x
    cc = ccx
    dic = dicx
    alk = alkx
    pt = ptx
    time = time +dt
    
enddo

endprogram 
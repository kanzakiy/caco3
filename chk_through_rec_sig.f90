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
! NOTE: Now this code is almost the same as the whole code. The difference is only that this code does not track any signals   

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

!********************************************************************************************************************************  ADDED-START
dep = 5.0d0
ccflxi = 12d-6
om2cc = 0.7d0
open(unit=file_sigmly,file=trim(adjustl(workdir))//'sigmly.txt',action='write',status='unknown')! recording signals etc at just below mixed layer 
open(unit=file_sigmlyd,file=trim(adjustl(workdir))//'sigmlyd.txt',action='write',status='unknown') ! recording signals etc at depths of 2x mixed layer thickness 
open(unit=file_sigbtm,file=trim(adjustl(workdir))//'sigbtm.txt',action='write',status='unknown')! ! recording signals etc at bottom of sediment 
! <<<<<<<<<<<<<<<<<<<<<  NEW: 05/13/2019 <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<
open(unit=file_bound,file=trim(adjustl(workdir))//'bound.txt',action='write',status='unknown')! recording boundary conditions changes 
! <<<<<<<<<<<<<<<<<<<<<  NEW: 05/13/2019  <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<

!********************************************************************************************************************************  ADDED-END

beta = 1.00000000005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
call makegrid(beta,nz,ztot,dz,z)

! call getporosity() ! assume porosity profile 
call getporosity(  &
     poro,porof,sporo,sporof,sporoi & ! output
     ,z,nz  & ! input
     )

!!!!!!!!!!!!! flx assignement and initial guess for burial rate !!!!!!!!!!!!!!!!!!!!!!
! call flxstat()  ! assume fluxes of om, cc and clay, required to calculate burial velocity
call flxstat(  &
    omflx,detflx,ccflx  & ! output
    ,om2cc,ccflxi,mcc,nspcc  & ! input 
    )
! print*,om2cc,ccflxi,detflx,omflx,sum(ccflx)
! molar volume (cm3 mol-1) needed for burial rate calculation 
mvom = mom/rhoom  ! om
mvsed = msed/rhosed ! clay 
mvcc = mcc/rhocc ! caco3
! call burial_pre() ! initial guess of burial profile, requiring porosity profile  
call burial_pre(  &
    w,wi  & ! output
    ,detflx,ccflx,nspcc,nz  & ! input 
    )  
! call dep2age() ! depth -age conversion 
call dep2age(  &
    age &  ! output 
    ,dz,w,nz  &  ! input
   )
! call calcupwindscheme() ! determine factors for upwind scheme to represent burial advection
call calcupwindscheme(  &
    up,dwn,cnr,adf & ! output 
    ,w,nz   & ! input &
    )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!********************************************************************************************************************************  ADDED-START


!!! ~~~~~~~~~~~~~~ set recording time 
call recordtime()

depi = 4d0  ! depth before event 
depf = dep   ! max depth to be changed to  

flxfini = 0.5d0  !  total caco3 rain flux for fine species assumed before event 
flxfinf = 0.9d0 !  maximum changed value 

! ///////////// isotopes  ////////////////
d13c_ocni = 2d0  ! initial ocean d13c value 
d13c_ocnf = -1d0 ! ocean d13c value with maximum change  
d18o_ocni = 1d0 ! initial ocean d18o value 
d18o_ocnf = -1d0 ! ocean d18o value with maximum change 

call sig2sp_pre() ! end-member signal assignment 




!********************************************************************************************************************************  ADDED-END






! call make_transmx()
call make_transmx(  &
    trans,izrec,izrec2,izml,nonlocal  & ! output 
    ,labs,nspcc,turbo2,nobio,dz,sporo,nz,z  & ! input
    )

! call coefs(temp,sal,dep)  ! need to specify diffusion coefficient as well as om decomposition rate const. etc.
call coefs(  &
    dif_dic,dif_alk,dif_o2,kom,kcc,co3sat & ! output 
    ,temp,sal,dep,nz,nspcc,poro,cai  & !  input 
    )

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
ptx = pt
! this may not be necessary as these individual species assume equilibrium 
co2x = co2
hco3x = hco3
co3x = co3

time = 0d0 ! model time [yr]
it = 1 ! integration count 
nt = 10 ! total integration 
dt = 1d2 ! time step [yr]

rho = 2.5d0 ! assume here 

oxco2 = 0d0  ! oxic degradation of om; here initially assumed 0   
anco2 = 0d0  ! anoxic degradation of om; here initially assumed 0   
!!!  addition to chk_om.f90 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
zox = 10d0  ! initial assumption on oxygen penetaration depth [cm]
!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! do it=1,nt
do   ! <<<--------------------------------------------------------------------------------------------------------------------------  Here change to unknwon iteration 


!********************************************************************************************************************************  ADDED-START

    !! ///////// isotopes & fluxes settings ////////////// 
#ifndef sense    
    call timestep(time,800,5000,1000,dt) ! determine time step dt by calling timestep(time,nt_spn,nt_trs,nt_aft,dt) where nt_xx denotes total iteration number  
    call signal_flx(time)
    call bdcnd(time,dep)
#endif 

    ! isotope signals represented by caco3 rain fluxes 
    d18o_flx = sum(d18o_sp(:)*ccflx(:))/ccflxi
    d13c_flx = sum(d13c_sp(:)*ccflx(:))/ccflxi

#ifndef track2
    if (abs(d13c_flx - d13c_ocn)>tol .or. abs(d18o_flx - d18o_ocn)>tol) then ! check comparability with input signals 
        print*,'error in assignment of proxy'
        write(file_err,*)'error in assignment of proxy',d18o_ocn,d13c_ocn,d18o_flx,d13c_flx
        stop
    endif 
#endif
 
! <<<<<<<<<<<<<<<<<<<<<  NEW: 05/13/2019  <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<
#ifndef size 
    !  recording fluxes of two types of caco3 separately 
    write(file_bound,*) time, d13c_ocn, d18o_ocn, (ccflx(isp),isp=1,nspcc),temp, dep, sal,dici,alki, o2i
#else 
    !  do not record separately 
    write(file_bound,*) time, d13c_ocn, d18o_ocn, sum(ccflx(1:4)),sum(ccflx(5:8)),(ccflx(isp),isp=1,nspcc),temp, dep, sal,dici,alki, o2i
#endif  
! <<<<<<<<<<<<<<<<<<<<<  NEW: 05/13/2019  <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<

    !! === temperature & pressure and associated boundary changes ====
    ! if temperature is changed during signal change event this affect diffusion coeff etc. 
    ! call coefs(temp,sal,dep)
    call coefs(  &
        dif_dic,dif_alk,dif_o2,kom,kcc,co3sat & ! output 
        ,temp,sal,dep,nz,nspcc,poro,cai  & !  input 
        )
    !! /////////////////////
!********************************************************************************************************************************  ADDED-END

    
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
    
        ! call omcalc() ! om conc. calculation 
        call omcalc( &
            omx,izox  & ! output 
            ,kom   &  ! in&output
            ,oxic,anoxic,o2x,om,nz,sporo,sporoi,sporof &! input 
            ,w,wi,dt,up,dwn,cnr,adf,trans,nspcc,labs,turbo2,nonlocal,omflx,poro,dz &! input 
            ) 
        ! calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations )
        ! call calcflxom()
        call calcflxom(  &
            omadv,omdec,omdif,omrain,omres,omtflx  & ! output 
            ,sporo,om,omx,dt,w,dz,z,nz,turbo2,labs,nonlocal,poro,up,dwn,cnr,adf,rho,mom,trans,kom,sporof,sporoi,wi,nspcc,omflx  & ! input 
            )
        
        ! print*,'~~~~ conc ~~~~'
        ! print dumchr(1), 'z  :',(z(iz),iz=1,nz,nz/interval)
        ! print dumchr(1), 'om :',(omx(iz)*mom/rho(iz)*100d0,iz=1,nz,nz/interval)
        ! print*,'++++ flx ++++'
        ! print'(7A11)', 'tflx','adv','dif','omrxn','ccrxn','rain','res'
        ! print'(A,7E11.3)', 'om :', omtflx, omadv,  omdif, omdec,0d0,omrain, omres
        
        ! print*,'izox',izox
        ! sb omcalc calculates izox, which is the deepest grid where o2 >=0. 
        
        if (izox == nz) then ! fully oxic; lower boundary condition ---> no diffusive out flow  
            ! call o2calc_ox()  ! o2 calculation when o2 penetration depth (zox) is the same as bottom depth. 
            call o2calc_ox(  &
                o2x  & ! output
                ,izox,nz,poro,o2,kom,omx,sporo,dif_o2,dz,dt & ! input
                )
            ! call calcflxo2_ox() !  fluxes relevant to o2 (at the same time checking the satisfaction of difference equations) 
            call calcflxo2_ox( &
                o2dec,o2dif,o2tflx,o2res  & ! output 
                ,nz,sporo,kom,omx,dz,poro,dif_o2,dt,o2,o2x  & ! input
                )
        else  !! if oxygen is depleted within calculation domain, lower boundary changes to zero concs.
            ! call o2calc_sbox() ! o2 calculation when o2 is depleted within the calculation domain.
            call o2calc_sbox(  &
                o2x  & ! output
                ,izox,nz,poro,o2,kom,omx,sporo,dif_o2,dz,dt & ! input
                )
            ! call calcflxo2_sbox() ! fluxes relevant to oxygen 
            call calcflxo2_sbox( &
                o2dec,o2dif,o2tflx,o2res  & ! output 
                ,nz,sporo,kom,omx,dz,poro,dif_o2,dt,o2,o2x,izox  & ! input
                )
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
    ! call calccaco3sys()
    call calccaco3sys(  &
        ccx,dicx,alkx,rcc,dt  & ! in&output
        ,nspcc,dic,alk,dep,sal,temp,labs,turbo2,nonlocal,sporo,sporoi,sporof,poro,dif_alk,dif_dic & ! input
        ,w,up,dwn,cnr,adf,dz,trans,cc,oxco2,anco2,co3sat,kcc,ccflx,ncc,omega,nz  & ! input
        )
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
    ! call calcflxcaco3sys()
    call calcflxcaco3sys(  &
         cctflx,ccflx,ccdis,ccdif,ccadv,ccrain,ccres,alktflx,alkdis,alkdif,alkdec,alkres & ! output
         ,dictflx,dicdis,dicdif,dicres,dicdec   & ! output
         ,dw & ! inoutput
         ,nspcc,ccx,cc,dt,dz,rcc,adf,up,dwn,cnr,w,dif_alk,dif_dic,dic,dicx,alk,alkx,oxco2,anco2,trans    & ! input
         ,turbo2,labs,nonlocal,sporof,it,nz,poro,sporo        & ! input
         )
    ! end of caco3 calculations 
    
    ! clay calculation 
    ! call claycalc()
    call claycalc(  &   
        ptx                  &  ! output
        ,nz,sporo,pt,dt,w,dz,detflx,adf,up,dwn,cnr,trans  &  ! input
        ,nspcc,labs,turbo2,nonlocal,poro,sporof     &  !  intput
        )
    ! call calcflxclay()
    call calcflxclay( &
        pttflx,ptdif,ptadv,ptres,ptrain  & ! output
        ,dw          &  ! in&output
        ,nz,sporo,ptx,pt,dt,dz,detflx,w,adf,up,dwn,cnr,sporof,trans,nspcc,turbo2,labs,nonlocal,poro           &  !  input
        )
    ! end of clay calculation 
    
    ! checking for total volume of solids, density and burial velocity 

    ! call getsldprop() ! get solid property, rho (density) and frt (total vol.frac)
    call getsldprop(  &
        rho,frt,       &  ! output
        nz,omx,ptx,ccx,nspcc,w,up,dwn,cnr,adf,z      & ! input
        )

    err_f = maxval(abs(frt - 1d0))  ! new error in total vol. fraction (must be 1 in theory) 
    !! ========= calculation of burial velocity =============================

    wx = w  ! recording previous burial velocity 

    ! call burialcalc() ! get new burial velocity
    call burialcalc(  &
        w,wi         & !  output
        ,detflx,ccflx,nspcc,omflx,dw,dz,poro,nz    & ! input
        )

    !  determine calculation scheme for advection 
    ! call calcupwindscheme()
    call calcupwindscheme(  &
        up,dwn,cnr,adf & ! output 
        ,w,nz   & ! input &
        )

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
            ! pause  ! <------------------------------------------------------------------------------------------------------------ comment out 
            go to 400
        endif 
    endif
    if (err_w > tol) go to 300

400 continue

!********************************************************************************************************************************  ADDED-START

    call dep2age(  &
        age &  ! output 
        ,dz,w,nz  &  ! input
       )

    ! ---------------------
    !/////// ISOTOPES /////
    !  calculating bulk isotopic composition
    do iz=1,nz 
        d18o_blk(iz) = sum(d18o_sp(:)*ccx(iz,:))/sum(ccx(iz,:))
        d13c_blk(iz) = sum(d13c_sp(:)*ccx(iz,:))/sum(ccx(iz,:))
#ifdef size
        d18o_blkf(iz) = sum(d18o_sp(1:4)*ccx(iz,1:4))/sum(ccx(iz,1:4))
        d13c_blkf(iz) = sum(d13c_sp(1:4)*ccx(iz,1:4))/sum(ccx(iz,1:4))
        d18o_blkc(iz) = sum(d18o_sp(5:8)*ccx(iz,5:8))/sum(ccx(iz,5:8))
        d13c_blkc(iz) = sum(d13c_sp(5:8)*ccx(iz,5:8))/sum(ccx(iz,5:8))
#endif 
    enddo
    
    ! recording 
    
    if (time>=rectime(cntrec)) then 
        call recordprofile(cntrec )
        
        cntrec = cntrec + 1
        if (cntrec == nrec+1) exit
    endif 
    
    
    call sigrec()  ! recording signals at 3 different depths (btm of mixed layer, 2xdepths of btm of mixed layer and btm depth of calculation domain)
    
    
!********************************************************************************************************************************  ADDED-END

    print'(A,E11.3)', 'error in frt:', maxval(abs(frt - 1d0))
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
    
    ! call recordprofile(it)  !! <----------------------------------------------------------------------------------- comment out 
    it = 1+it            ! <---------------------------------------------------------------------------------------------------- added 
enddo


!********************************************************************************************************************************  ADDED-START
close(file_sigmly)! recording signals etc at just below mixed layer 
close(file_sigmlyd) ! recording signals etc at depths of 2x mixed layer thickness 
close(file_sigbtm)! ! recording signals etc at bottom of sediment  
! <<<<<<<<<<<<<<<<<<<<<  NEW: 05/13/2019 <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<
close(file_bound)! recording boundary conditions changes 
! <<<<<<<<<<<<<<<<<<<<<  NEW: 05/13/2019  <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<< <<<<<<<<<<<<<<<<<<<<<
!********************************************************************************************************************************  ADDED-END
endprogram 
program caco3 
! trying a simple diagenesis
! irregular grid
! attempt to allow negative burial from v3_2 (which works until non-local mixing is considered) 
! allow multiple carbonate species developed from v3_3 
! allow different transition matrix and reactivity for different species, developed from v5_2 
! enabling separate signal tracking developed from v5_3
! enabling large number of species with sparse matrix solver
! to compile, type 'gfortran -c caco3_therm.f90 -g -fcheck=all;
! gfortran -cpp -Dtest  caco3_test_mod_v5_5.f90 caco3_therm.o umf4_f77wrapper.o -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lopenblas -g -fcheck=all'
! allow fraction to advection scheme e.g., 0.5 upwind + 0.5 central scheme
implicit none

integer(kind=4),parameter :: nz = 100
#ifdef sense
integer(kind=4),parameter :: nspcc = 1
#elif defined track2
integer(kind=4),parameter :: nspcc = 22
#elif defined size
integer(kind=4),parameter :: nspcc = 8
#else
integer(kind=4),parameter :: nspcc = 4
#endif
real(kind=8) cc(nz,nspcc),ccx(nz,nspcc)  ! mol cm-3 sld
real(kind=8) om(nz),omx(nz)  ! mol cm-3 sld
real(kind=8) pt(nz), ptx(nz) ! # of specific particles/# of total particles ??
real(kind=8) ccflx(nspcc), d13c_sp(nspcc),d18o_sp(nspcc), d13c_blk(nz), d18o_blk(nz)
real(kind=8) d13c_blkf(nz), d18o_blkf(nz),  d13c_blkc(nz), d18o_blkc(nz)
real(kind=8) d13c_ocni,d13c_ocnf,d18o_ocni,d18o_ocnf,d13c_ocn, d18o_ocn, flxfrc(nspcc),flxfrc2(nspcc)
real(kind=8) d13c_flx, d18o_flx
real(kind=8) :: ccflxi = 10d-6 ! mol (CaCO3) cm-2 yr-1  ! Emerson and Archer (1990) 
real(kind=8) :: omflx = 12d-6 ! mol cm-2 yr-1       ! Emerson and Archer (1990)
real(kind=8) :: detflx = 180d-6 ! g cm-2 yr-1  ! MUDS input http://forecast.uchicago.edu/Projects/muds.html
real(kind=8) :: alki =  2285d0 ! uM
real(kind=8) :: dici = 2211d0 ! uM
real(kind=8) :: o2i = 165d0 ! uM
! real(kind=8) :: komi = 2d0  ! /yr  ! MUDS 
real(kind=8) :: komi = 0.5d0  ! /yr  ! arbitrary 
! real(kind=8) :: komi = 0.1d0  ! /yr  ! Canfield 1994
! real(kind=8) :: komi = 0.06d0  ! /yr  ! ?? Emerson 1985? who adopted relatively slow decomposition rate 
! real(kind=8) :: kcci = 10.0d0*365.25d0  ! /yr  0.15 to 30 d-1 Emerson and Archer (1990) 0.1 to 10 d-1 in Archer 1991
real(kind=8) :: kcci = 1d0*365.25d0  ! /yr  0.15 to 30 d-1 Emerson and Archer (1990) 0.1 to 10 d-1 in Archer 1991
! real(kind=8) :: kcci = 0d0*365.25d0  ! /yr  0.15 to 30 d-1 Emerson and Archer (1990) 0.1 to 10 d-1 in Archer 1991
real(kind=8) :: poroi = 0.8d0  
real(kind=8) :: keqcc = 4.4d-7   ! mol2 kg-2 (Mucci 1983 cited by Emerson and Archer 1990)
real(kind=8) :: ncc = 4.5d0   ! (Archer et al. 1989)
real(kind=8) :: temp = 2d0  ! C 
real(kind=8) :: sal = 35d0  ! wt o/oo 
real(kind=8) :: cai = 10.3d-3 ! mol kg-1
real(kind=8) :: fact = 1d-3 ! w/r factor to facilitate calculation 
! real(kind=8) :: rhosed = 2.09d0 ! g/cm3 sediment particle density assuming opal
real(kind=8) :: rhosed = 2.6d0 ! g/cm3 sediment particle density assming kaolinite 
real(kind=8) :: rhoom = 1.2d0 ! g/cm3 organic particle density 
real(kind=8) :: rhocc = 2.71d0 ! g/cm3 organic particle density 
real(kind=8) :: mom = 30d0 ! g/mol OM assuming CH2O
! real(kind=8) :: msed = 87.11d0 ! g/mol arbitrary sediment g/mol assuming opal (SiO2•n(H2O) )
real(kind=8) :: msed = 258.16d0 ! g/mol arbitrary sediment g/mol assuming kaolinite ( 	Al2Si2O5(OH)4 )
real(kind=8) :: mcc = 100d0 ! g/mol CaCO3 
real(kind=8) :: ox2om = 1.3d0 ! o2/om ratio for om decomposition (Emerson 1985; Archer 1991)
real(kind=8) :: om2cc = 0.666d0  ! rain ratio of organic matter to calcite
! real(kind=8) :: om2cc = 0.5d0  ! rain ratio of organic matter to calcite
real(kind=8) :: o2th = 0d0 ! threshold oxygen level below which not to calculate 
real(kind=8) :: dev = 1d-6 ! deviation addumed 
real(kind=8) :: zml_ref = 12d0 ! mixed layer depth
real(kind=8) dif_alk0, dif_dic0, dif_o20, zml(nspcc+2) , zrec, zrec2, chgf, flxfin, flxfini, flxfinf
real(kind=8) pore_max, exp_pore, calgg  ! parameters to determine porosity in Archer
real(kind=8) mvom, mvsed, mvcc  ! molar volumes (cm3 mol-1) mv_i = m_i/rho_i
real(kind=8) keq1, keq2, calceq1,calceq2, co3sat, keqag, dep, zox, zoxx, depi,depf, corrf, df, err_f, err_fx
real(kind=8) dic(nz), dicx(nz), alk(nz), alkx(nz), o2(nz), o2x(nz) !  mol cm-3 porewater 
real(kind=8) dbio(nz),dif_dic(nz), dif_alk(nz), dif_o2(nz), ff(nz)
real(kind=8) co2(nz), hco3(nz), co3(nz), pro(nz)
real(kind=8) co2x(nz), hco3x(nz), co3x(nz), prox(nz),co3i
real(kind=8) poro(nz), rho(nz), frt(nz), sporo(nz), sporoi, porof, sporof
real(kind=8) rcc(nz,nspcc), drcc_dcc(nz,nspcc), drcc_ddic(nz,nspcc), drcc_dalk(nz,nspcc), drcc_dco3(nz,nspcc) 
real(kind=8) dco3_ddic(nz), dco3_dalk(nz), ddum(nz), dpro_dalk(nz), dpro_ddic(nz)
real(kind=8) kcc(nz,nspcc) 
real(kind=8) kom(nz), oxco2(nz),anco2(nz) 
real(kind=8) w(nz) , wi, dw(nz), wx(nz), err_w, wxx(nz),err_f_min, dfrt_df, d2frt_df2, dfrt_dfx, err_w_min
real(kind=8) z(nz), dz(nz), eta(nz), beta, dage(nz), age(nz) 
real(kind=8) :: ztot = 500d0 ! cm 
integer(kind=4) :: nsp = 3  ! independent chemical variables 
integer(kind=4) :: nmx 
real(kind=8),allocatable :: amx(:,:),ymx(:),emx(:)
integer(kind=4),allocatable :: ipiv(:),dumx(:,:)
integer(kind=4) infobls, infosbr
real(kind=8) error, error2, minerr
real(kind=8) :: tol = 1d-6
integer(kind=4) iz, row, col, itr  , it, iiz, itr_w, cntsp, izrec, izrec2, itr_f
integer(kind=4) :: nt = 1000000
integer(kind=4),parameter :: nrec = 15
integer(kind=4) :: cntrec, itrec, izox_minerr =0
real(kind=4) up(nz), dwn(nz), cnr(nz) ! advection calc. schemes; up or down wind, or central schemes if 1
real(kind=8) adf(nz)  ! factor to make sure mass conversion 
! integer(kind=4),parameter :: nt = 1
real(kind=8) :: time, dt = 1d2
real(kind=8) :: rectime(nrec), time_max, dumreal, time_spn, time_trs, time_aft
! real(kind=8) :: dtcc=1d5, dtom = 1d5
real(kind=8) :: omadv, omdec, omdif, omrain, omres, omtflx 
real(kind=8) :: o2tflx, o2dec, o2dif, o2res 
real(kind=8) :: cctflx(nspcc), ccdis(nspcc), ccdif(nspcc), ccadv(nspcc), ccres(nspcc),ccrain(nspcc) 
real(kind=8) :: dictflx, dicdis, dicdif, dicdec, dicres 
real(kind=8) :: alktflx, alkdis, alkdif, alkdec, alkres 
real(kind=8) :: pttflx, ptdif, ptadv, ptres, ptrain 
real(kind=8) :: trans(nz,nz,nspcc+2), transdbio(nz,nz), translabs(nz,nz),transturbo2(nz,nz), translabs_tmp(nz,nz)
character*512 workdir, filechr
character*10 dumchr(3)
character*25 arg, chr(3,4)
integer(kind=4) dumint(8), idp, izox, narg, ia, izml, isp, ilabs, nlabs
integer(kind=4) :: file_tmp=100,file_ccflx=101,file_omflx=102,file_o2flx=103,file_dicflx=104,file_alkflx=105,file_ptflx=106  
integer(kind=4) :: file_err=107,file_ccflxes(nspcc), file_bound=108, file_totfrac=109, file_sigmly=110,file_sigmlyd=111
integer(kind=4) :: file_sigbtm=112 
external dgesv
logical :: oxic = .true.  ! oxic only model of OM degradation by Emerson (1985) 
logical :: anoxic = .true.  ! oxic-anoxic model of OM degradation by Archer (1991) 
logical :: nobio(nspcc+2) = .false.  ! no biogenic reworking assumed 
logical :: turbo2(nspcc+2) = .false.  ! random mixing 
logical :: labs(nspcc+2) = .false.  ! mixing info from LABS 
logical :: nonlocal(nspcc+2)  ! ON if assuming non-local mixing (i.e., if labs or turbo2 is ON)
real(kind=8) :: calceqcc, calceqag
real(kind=8), parameter :: pi=4.0d0*atan(1.0d0) 
! when using sparse matrix solver 
integer(kind=4) n, nnz
integer(kind=4), allocatable :: ai(:), ap(:) 
real(kind=8), allocatable :: ax(:), bx(:) 
real(kind=8) :: control(20)
integer(kind=4) i,j,cnt,cnt2
real(kind=8) info(90)
integer(kind=8) numeric
integer(kind=4) status
integer(kind=8) symbolic
integer(kind=4) sys
real(kind=8), allocatable :: kai(:)

call date_and_time(dumchr(1),dumchr(2),dumchr(3),dumint)

!!! get variables !!!
narg = iargc()
do ia = 1, narg,2
    call getarg(ia,arg)
    select case(trim(arg))
        case('cc','CC','Cc')
            call getarg(ia+1,arg)
            read(arg,*)ccflxi 
        case('rr','RR','Rr')
            call getarg(ia+1,arg)
            read(arg,*)om2cc
        case('dep','DEP','Dep')
            call getarg(ia+1,arg)
            read(arg,*)dep
        case('dt','DT','Dt')
            call getarg(ia+1,arg)
            read(arg,*)dt
        case('fl','FL','Fl')
            call getarg(ia+1,arg)
            read(arg,*)filechr
    end select
enddo
print'(3A,3E11.3)','ccflxi','om2cc','dep:',ccflxi,om2cc, dep
do ia = 1,3
    if (ia==1) dumreal=ccflxi
    if (ia==2) dumreal=om2cc
    if (ia==3) dumreal=dep
    if (dumreal/=0d0) then 
        write(chr(ia,1),'(i0)') floor(log10(dumreal))
        write(chr(ia,2),'(i0)') int(dumreal/(10d0**(floor(log10(dumreal)))))
        write(chr(ia,3),'(i0)') int((dumreal/(10d0**(floor(log10(dumreal)))) &
            - int(dumreal/(10d0**(floor(log10(dumreal))))))*10d0)
    else 
        write(chr(ia,1),'(i0)') 0
        write(chr(ia,2),'(i0)') 0
        write(chr(ia,3),'(i0)') 0
    endif
    
    chr(ia,4) = trim(adjustl(chr(ia,2)))//'_'//trim(adjustl(chr(ia,3)))//'E'//trim(adjustl(chr(ia,1)))
    
enddo 
print'(6A)','ccflx','om2cc','dep:',(chr(ia,4),ia=1,3)
! pause

!!!!!!!!!!!!!!
#ifndef nonrec
!! FILES !!!!!!!!!
workdir = 'C:/Users/YK/Desktop/Sed_res/'
workdir = trim(adjustl(workdir))//'test-translabs/profiles/'
workdir = trim(adjustl(workdir))//'multi/'
#ifdef test 
workdir = trim(adjustl(workdir))//'test/'
#endif
if (.not. anoxic) then 
    workdir = trim(adjustl(workdir))//'ox'
else 
    workdir = trim(adjustl(workdir))//'oxanox'
endif
if (any(labs)) workdir = trim(adjustl(workdir))//'_labs'
if (any(turbo2)) workdir = trim(adjustl(workdir))//'_turbo2'
if (any(nobio)) workdir = trim(adjustl(workdir))//'_nobio'
workdir = trim(adjustl(workdir))//'/'
workdir = trim(adjustl(workdir))//'cc-'//trim(adjustl(chr(1,4)))//'_rr-'//trim(adjustl(chr(2,4)))  &
#ifdef sense
    //'_dep-'//trim(adjustl(chr(3,4)))
#else
    //'_'//trim(adjustl(filechr))
#endif
! workdir = trim(adjustl(workdir))//'-'//trim(adjustl(dumchr(1)))  ! adding date
call system ('mkdir -p '//trim(adjustl(workdir)))
workdir = trim(adjustl(workdir))//'/'

open(unit=file_ptflx,file=trim(adjustl(workdir))//'ptflx.txt',action='write',status='unknown')
open(unit=file_ccflx,file=trim(adjustl(workdir))//'ccflx.txt',action='write',status='unknown')
open(unit=file_omflx,file=trim(adjustl(workdir))//'omflx.txt',action='write',status='unknown')
open(unit=file_o2flx,file=trim(adjustl(workdir))//'o2flx.txt',action='write',status='unknown')
open(unit=file_dicflx,file=trim(adjustl(workdir))//'dicflx.txt',action='write',status='unknown')
open(unit=file_alkflx,file=trim(adjustl(workdir))//'alkflx.txt',action='write',status='unknown')
open(unit=file_err,file=trim(adjustl(workdir))//'errlog.txt',action='write',status='unknown')
open(unit=file_bound,file=trim(adjustl(workdir))//'bound.txt',action='write',status='unknown')
open(unit=file_totfrac,file=trim(adjustl(workdir))//'frac.txt',action='write',status='unknown')
open(unit=file_sigmly,file=trim(adjustl(workdir))//'sigmly.txt',action='write',status='unknown')
open(unit=file_sigmlyd,file=trim(adjustl(workdir))//'sigmlyd.txt',action='write',status='unknown')
open(unit=file_sigbtm,file=trim(adjustl(workdir))//'sigbtm.txt',action='write',status='unknown')
do isp = 1,nspcc
    file_ccflxes(isp)=800+(isp-1)
    write(dumchr(1),'(i3.3)') isp
    open(unit=file_ccflxes(isp),file=trim(adjustl(workdir))//'ccflx-sp_'//trim(adjustl(dumchr(1))) &
        //'.txt',action='write',status='unknown')
enddo
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!  MAKING GRID !!!!!!!!!!!!!!!!!
beta = 1.00000000005d0

do iz = 1, nz 
    ! z(iz) = (iz-1)*ztot/(nz-1)
    z(iz) = iz*ztot/nz
    eta(iz) = log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0))
    if (iz==1) then
        dz(iz) = ztot*log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0))
    endif
    if (iz/=1) then 
        dz(iz) = ztot*log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0)) - sum(dz(:iz-1))
    endif
enddo

do iz=1,nz
    if (iz==1) z(iz)=dz(iz)*0.5d0
    if (iz/=1) z(iz) = z(iz-1)+dz(iz-1)*0.5d0 + 0.5d0*dz(iz)
enddo

!~~~~~~~~~~~~~ saving grid for LABS ~~~~~~~~~~~~~~~~~~~~~~
#ifdef recgrid
open(unit=file_tmp, file='C:/cygwin64/home/YK/LABS/1dgrid.txt',action='write',status='unknown')
do iz = 1, nz
    write(file_tmp,*) dz(iz)
enddo
close(file_tmp)
#endif
! stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! FUNDAMENTAL PARAMETERS !!!!!!!!!!!!!

!! Archer's (1991) parameterization 
! ccflx = 18d-6 ! mol cm-2 yr-1
! omflx = 12d-6 ! mol cm-2 yr-1

! om2cc = 0.666d0
! ccflx = 60d-6
! dep = 6d0   ! sea floor depth in km

omflx = om2cc*ccflxi
detflx = (1d0/9d0)*ccflxi*mcc ! 90% of mass flux becomes inorganic C; g cm-2 yr-1
ccflx = ccflxi/nspcc

mvom = mom/rhoom
mvsed = msed/rhosed
mvcc = mcc/rhocc
    

poro = poroi  ! constant porosity 
! ----------- Archer's parameterization 
calgg = 0.0d0
pore_max =  1d0 - ( 0.483d0 + 0.45d0 * calgg) / 2.5d0
exp_pore = 0.25d0*calgg + 3.d0 *(1d0-calgg)
poro = EXP(-z/exp_pore) * (1.d0-pore_max) + pore_max 
porof = pore_max
porof = poro(nz)
sporof = 1d0-porof
! ------------------
! print*,poro
! stop
! cf., poro = poro_0*exp(-z(iz)/poro_scale)  ! Hydrate modeling parameterization where poro_0 = 0.69 & poro_scale = 2000 (m)
sporoi = 1d0-poroi
sporo = 1d0 - poro

w = 0.003d0 ! cm/yr 

wi = (detflx/msed*mvsed + sum(ccflx)*mvcc +omflx*mvom)/(1d0-poroi)
wi = (detflx/msed*mvsed + sum(ccflx)*mvcc            )/(1d0-poroi)
! wi = (detflx/msed*msed + ccflx*mcc + omflx*mom)/(1d0-poroi)/2.5d0

! print*,wi
! pause

500 continue

w= wi

!! depth -age 
dage = dz/w  
age = 0d0
do iz=1,nz
    if (iz==1) age(iz)=dage(iz)*0.5d0
    if (iz/=1) age(iz) = age(iz-1)+dage(iz-1)*0.5d0 + 0.5d0*dage(iz)
enddo

! ------------ determine calculation scheme for advection 
up = 0
dwn=0
cnr =0
adf = 1d0
do iz=1,nz 
    if (iz==1) then 
        if (w(iz)>=0d0 .and. w(iz+1)>=0d0) then
            up(iz) = 1
        elseif (w(iz)<=0d0 .and. w(iz+1)<=0d0) then
            dwn(iz) = 1
        else 
            if (.not.(w(iz)*w(iz+1) <=0d0)) then 
                print*,'error'
                stop
            endif
            cnr(iz) = 1
        endif
    elseif (iz==nz) then 
        if (w(iz)>=0d0 .and. w(iz-1)>=0d0) then
            up(iz) = 1
        elseif (w(iz)<=0d0 .and. w(iz-1)<=0d0) then
            dwn(iz) = 1
        else 
            if (.not.(w(iz)*w(iz-1) <=0d0)) then 
                print*,'error'
                stop
            endif
            cnr(iz) = 1
        endif
    else 
        if (w(iz) >=0d0) then 
            if (w(iz+1)>=0d0 .and. w(iz-1)>=0d0) then
                up(iz) = 1
            else
                cnr(iz) = 1
            endif
        else  
            if (w(iz+1)<=0d0 .and. w(iz-1)<=0d0) then
                dwn(iz) = 1
            else
                cnr(iz) = 1
            endif
        endif
    endif
enddo        
! up = 0
! dwn=0
! cnr =1

if (sum(up(:)+dwn(:)+cnr(:))/=nz) then
    print*,error,sum(up),sum(dwn),sum(cnr)
    stop
endif
! ---------------------

zox = 10d0 ! priori assumed oxic zone 

!!! ~~~~~~~~~~~~~~ set recording time 
!! simulation to steady state 
! time_max = ztot / wi ! yr
! do itrec=1,nrec 
    ! rectime(itrec)=itrec*time_max/real(nrec)*10d0
! enddo
!!!!
!!! +++ tracing experiment 
time_spn = ztot / wi *50d0 ! yr
time_trs = 50d3 !  5 kyr transition
! time_aft = time_spn 
time_aft = time_trs*3d0 
#ifdef sense
time_trs = 0d0
time_aft = 0d0 
#endif 
do itrec=1,nrec/3
    rectime(itrec)=itrec*time_spn/real(nrec/3)
enddo
do itrec=nrec/3+1,nrec/3*2
    rectime(itrec)=rectime(nrec/3)+(itrec-nrec/3)*time_trs/real(nrec/3)
enddo
do itrec=nrec/3*2+1,nrec
    rectime(itrec)=rectime(nrec/3*2)+(itrec-nrec/3*2)*time_aft/real(nrec/3)
enddo
!!! ++++
cntrec = 1
#ifndef nonrec
open(unit=file_tmp,file=trim(adjustl(workdir))//'rectime.txt',action='write',status='unknown')
do itrec=1,nrec 
    write(file_tmp,*) rectime(itrec)
enddo
close(file_tmp)
#endif

depi = 4d0
depf = 5d0

flxfini = 0.5d0
flxfinf = 0.9d0

! ///////////// isotopes  ////////////////
d13c_ocni = 2d0
d13c_ocnf = -1d0
d18o_ocni = 1d0
d18o_ocnf = -1d0

! do isp = 1,nspcc 
    ! d13c_sp(isp) = d13c_ocni + (d13c_ocnf - d13c_ocni)*(isp-1)/real(nspcc-1)
    ! d18o_sp(isp) = d18o_ocni + (d18o_ocnf - d18o_ocni)*(isp-1)/real(nspcc-1)
! enddo
#ifndef sense
d13c_sp(1)=d13c_ocni
d18o_sp(1)=d18o_ocni
d13c_sp(2)=d13c_ocni
d18o_sp(2)=d18o_ocnf
d13c_sp(3)=d13c_ocnf
d18o_sp(3)=d18o_ocni
d13c_sp(4)=d13c_ocnf
d18o_sp(4)=d18o_ocnf
#else
d18o_sp=0d0
d13c_sp=0d0
#endif 

#ifdef size 
d13c_sp(5)=d13c_ocni
d18o_sp(5)=d18o_ocni
d13c_sp(6)=d13c_ocni
d18o_sp(6)=d18o_ocnf
d13c_sp(7)=d13c_ocnf
d18o_sp(7)=d18o_ocni
d13c_sp(8)=d13c_ocnf
d18o_sp(8)=d18o_ocnf
#endif 
!! //////////////////////////////////
!!!~~~~~~~~~~~~~~~~~

!!!! TRANSITION MATRIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~~~~~~~~~~~~ loading transition matrix from LABS ~~~~~~~~~~~~~~~~~~~~~~~~
if (any(labs)) then
    translabs = 0d0
    nlabs = 7394
    do ilabs=1,nlabs
        translabs_tmp=0d0
        write(dumchr(1),'(i0)') ilabs*2000
        open(unit=file_tmp, file='C:/Users/YK/Desktop/biot-res/trans-test-1kyr-frq-20190315/mix/' &
            //'transmtx-'//trim(adjustl(dumchr(1)))//'.txt',action='read',status='unknown')
        ! open(unit=file_tmp, file='C:/Users/YK/Desktop/biot-res/'// &
            ! 'trans-test-100yr-20190108' &
            ! //'/data/'//&
            ! 'transmtx-219000' &
            ! //'.txt',action='read',status='unknown')
        do iz = 1, nz
            read(file_tmp,*) translabs_tmp(iz,:)
        enddo
        close(file_tmp)
        translabs = translabs + translabs_tmp
    enddo
    translabs = translabs/real(nlabs)
endif

if (.true.) then 
! if (.false.) then 
    ! translabs = translabs *365.25d0/10d0*1d0/2.3d0
    translabs = translabs *365.25d0/10d0*1d0/10d0
endif

if (.false.) then 
! if (.true.) then 
    do iz=1,nz
        ! translabs(iz,iz) = 0d0
        do iiz=1,nz
            if (iiz==iz) cycle
            ! translabs(iz,iiz) = translabs(iz,iiz)*(1d0-poro(iz))/dz(iz)*dz(iiz)/(1d0-poro(iiz))
            translabs(iiz,iz) = translabs(iiz,iz)*(1d0-poro(iz))/dz(iz)*dz(iiz)/(1d0-poro(iiz))
        enddo
    enddo
endif 

if (.false.) then 
! if (.true.) then 
    do isp=1,nspcc+2
        do iiz = 1, nz
            ! trans(iiz,:) = trans(iiz,:)*dz(:)
            trans(iiz,:,isp) = trans(iiz,:,isp)/(1d0-poro(:))
        enddo
    enddo
endif
! stop
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
zml=zml_ref

zrec = 1.1d0*maxval(zml)
zrec2 = 2.0d0*maxval(zml)

#ifdef size 
zml(2+1)=20d0 
zml(2+2)=20d0 
zml(2+3)=20d0 
zml(2+4)=20d0 
zrec = 1.1d0*minval(zml)
zrec2 = 1.1d0*maxval(zml)
#endif 

do iz=1,nz
    if (z(iz)<=zrec) izrec = iz
    if (z(iz)<=zrec2) izrec2 = iz
enddo

nonlocal = .false.
do isp=1,nspcc+2
    if (turbo2(isp) .or. labs(isp)) nonlocal(isp)=.true.
    
    dbio=0d0
    do iz = 1, nz
        if (z(iz) <=zml(isp)) then
            dbio(iz) =  0.15d0
            izml = iz
        else
            dbio(iz) =  0d0
        endif
    enddo

    ! dbio = 0.15d0
    ! if (nobio) dbio = 0d0

    transdbio = 0d0
    do iz = 1, izml
        if (iz==1) then 
            transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        elseif (iz==izml) then 
            transdbio(iz,iz) = 0.5d0*(sporo(Iz)*dbio(iz)+sporo(Iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
        else 
            transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))  &
                + 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        endif
    enddo

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    ! transition matrix for random mixing 
    transturbo2 = 0d0
    transturbo2(:izml,:izml) = 0.0015d0/1d0
    do iz=1,izml
       transturbo2(iz,iz)=-0.0015d0*(izml-1)/1d0
    enddo

    if (turbo2(isp)) translabs = transturbo2

    trans(:,:,isp) = transdbio(:,:)

    if (nonlocal(isp)) trans(:,:,isp) = translabs(:,:)
    if (nobio(isp)) trans(:,:,isp) = 0d0
enddo

if (all(.not.nonlocal)) then 
    do isp=1,nspcc+2-1
        if (any(trans(:,:,isp+1)/=trans(:,:,isp))) nonlocal=.true.
    enddo
endif 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!~~~~~~~~diffusion & reaction~~~~~~~~~~~~~~~~~~~~~~
! ff = poro*poro*poro  ! Archer 
ff = poro*poro       ! Huelse et al., 2018

dif_dic0 = (151.69d0 + 7.93d0*temp) ! cm2/yr at 2 oC (Huelse et al. 2018)
dif_alk0 = (151.69d0 + 7.93d0*temp) ! cm2/yr  (Huelse et al. 2018)
dif_o20 =  (348.62d0 + 14.09d0*temp) 

dif_dic = dif_dic0*ff
dif_alk = dif_alk0*ff
dif_o2 = dif_o20*ff

kom = komi 
kcc = kcci

#ifdef size 
kcc(:,1) = kcci*10d0
kcc(:,2) = kcci*10d0
kcc(:,3) = kcci*10d0
kcc(:,4) = kcci*10d0
#endif 

keq1 = calceq1(temp,sal,dep)
keq2 = calceq2(temp,sal,dep)

keqcc = calceqcc(temp,sal,dep)
co3sat = keqcc/cai

!!   INITIAL CONDITIONS !!!!!!!!!!!!!!!!!!!    
cc = 1d-8   
! cc = ccflx/100d0/(1d0-poroi)/w 
dic = dici*1d-6/1d3 ! mol/cm3 
alk = alki*1d-6/1d3 ! mol/cm3

call calcspecies(dic,alk,temp,sal,dep,pro,co2,hco3,co3,nz,infosbr)  
    
pt = 1d-8
! do iz = 1, nz
    ! if (z(iz) < 2d0) pt(iz) = 1d0
! enddo
om = 1d-8
o2 = o2i*1d-6/1d3 ! mol/cm3

! ~~~ passing to temporary variables ~~~~~~~~~~~
ccx = cc
dicx = dic
alkx = alk 

co2x = co2
hco3x = hco3
co3x = co3

co3i=co3(1)

ptx = pt

omx = om
o2x = o2

do isp=1,nspcc 
    rcc(:,isp) = kcc(:,isp)*ccx(:,isp)*abs(1d0-co3x(:)*1d3/co3sat)**ncc*merge(1d0,0d0,(1d0-co3x(:)*1d3/co3sat)>0d0)
enddo

oxco2 = 0d0
anco2 = 0d0

! ~~~ saving initial conditions 
#ifndef nonrec
write(dumchr(1),'(i3.3)') 0

open(unit=file_tmp,file=trim(adjustl(workdir))//'ptx-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace') 
do iz = 1,nz
    write(file_tmp,*) z(iz),age(iz),pt(iz)*msed/2.5d0*100,0d0,1d0  ,wi
enddo
close(file_tmp)

open(unit=file_tmp,file=trim(adjustl(workdir))//'ccx-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace') 
do iz = 1,nz
    write(file_tmp,*) z(iz),age(iz),sum(cc(iz,:))*100d0/2.5d0*100d0, dic(iz)*1d3, alk(iz)*1d3, co3(iz)*1d3-co3sat &
        , sum(rcc(iz,:)),-log10(pro(iz)) 
enddo
close(file_tmp)

open(unit=file_tmp,file=trim(adjustl(workdir))//'o2x-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace') 
do iz = 1,nz
    write(file_tmp,*) z(iz),age(iz),o2x(iz)*1d3, oxco2(iz), anco2(iz)
enddo
close(file_tmp)

open(unit=file_tmp,file=trim(adjustl(workdir))//'omx-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace') 
do iz = 1,nz
    write(file_tmp,*) z(iz),age(iz),om(iz)*mom/2.5d0*100d0
enddo
close(file_tmp)
    
open(unit=file_tmp,file=trim(adjustl(workdir))//'ccx_sp-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace') 
do iz = 1,nz
    write(file_tmp,*) z(iz),age(iz),(cc(iz,isp)*mcc/2.5d0*100d0,isp=1,nspcc) 
enddo
close(file_tmp)

open(unit=file_tmp,file=trim(adjustl(workdir))//'sig-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace') 
do iz = 1,nz
    write(file_tmp,*) z(iz),age(iz),d13c_ocni,d18o_ocni
enddo
close(file_tmp)

open(unit=file_tmp,file=trim(adjustl(workdir))//'bur-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace') 
do iz = 1,nz
    write(file_tmp,*) z(iz),age(iz),w(iz),up(iz),dwn(iz),cnr(iz),adf(iz)
enddo
close(file_tmp)
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~
!! START OF TIME INTEGLATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

time = 0d0
it = 1
do 

! dt = 0.25d0  ! yr
! dt = 1d0  ! yr
! dt = dtom  ! yr

!! ///////// isotopes & fluxes settings ////////////// 
#ifndef sense    
if (time <= time_spn) then 
    if (time+dt> time_spn) then
        dt=time_trs/5000d0 
    else
        dt = time_spn/ 800d0
    endif
    ! dt=time_trs/100d0 
    ! dt = 100d0
    d13c_ocn = d13c_ocni
    d18o_ocn = d18o_ocni
    ccflx = 0d0
    ccflx(1) = ccflxi
    dep = depi
#ifdef track2
    cntsp=1
    d18o_sp(cntsp)=d18o_ocn
    d13c_sp(cntsp)=d13c_ocn
    ccflx = 0d0
    ccflx(cntsp) = ccflxi
#endif 
#ifdef size 
    ccflx = 0d0
    flxfin = flxfini
    ccflx(1) = flxfin*ccflxi
    ccflx(5) = (1d0-flxfin)*ccflxi
#endif 
elseif (time>time_spn .and. time<=time_spn+time_trs) then 
    dt = time_trs/5000d0
    d13c_ocn = d13c_ocni + (time-time_spn)*(d13c_ocnf-d13c_ocni)/time_trs
    d18o_ocn = d18o_ocni + (time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs 
    if (time-time_spn<=time_trs/2d0) then
        ! d13c_ocn = d13c_ocni + (time-time_spn)*(d13c_ocnf-d13c_ocni)/time_trs*2d0
        d18o_ocn = d18o_ocni + (time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs*2d0
        flxfin = flxfini + (time-time_spn)*(flxfinf-flxfini)/time_trs*2d0
    else 
        ! d13c_ocn = 2d0*d13c_ocnf - d13c_ocni  - (time-time_spn)*(d13c_ocnf-d13c_ocni)/time_trs*2d0
        d18o_ocn = 2d0*d18o_ocnf - d18o_ocni - (time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs*2d0
        flxfin = 2d0*flxfinf - flxfini - (time-time_spn)*(flxfinf-flxfini)/time_trs*2d0
    endif
    
    if (time-time_spn<=time_trs/10d0) then
        d13c_ocn = d13c_ocni + (time-time_spn)*(d13c_ocnf-d13c_ocni)/time_trs*10d0
        ! d18o_ocn = d18o_ocni + (time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs*2d0
        dep = depi + (depf-depi)*(time-time_spn)/time_trs*10d0
    elseif (time-time_spn>time_trs/10d0 .and. time-time_spn<=time_trs/10d0*9d0) then
        d13c_ocn = d13c_ocnf
        dep = depf
    elseif  (time-time_spn>time_trs/10d0*9d0) then 
        d13c_ocn = 10d0*d13c_ocnf - 9d0*d13c_ocni  - (time-time_spn)*(d13c_ocnf-d13c_ocni)/time_trs*10d0
        ! d18o_ocn = 2d0*d18o_ocnf - d18o_ocni - (time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs*2d0
        dep = 10d0*depf-9d0*depi - (depf-depi)*(time-time_spn)/time_trs*10d0
    endif
    
    ! d18o_ocn = (d18o_ocni+d18o_ocnf)*0.5d0 + (d18o_ocni-d18o_ocnf)*0.5*cos((time-time_spn)*2d0*pi/(time_trs/4d0)) 
    
    ! print*,d13c_ocn
    if (.not.(d13c_ocn>=d13c_ocnf .and.d13c_ocn<=d13c_ocni)) then 
        print*,'error in d13c',d13c_ocn
        stop
    endif
    ! do isp=1,nspcc-1
        ! if ((d13c_sp(isp)-d13c_ocn)*(d13c_sp(isp+1)-d13c_ocn)<=0d0) exit 
    ! enddo 
    ! ccflx = 0d0
    ! ccflx(isp) = ccflxi*abs((d13c_sp(isp+1)-d13c_ocn))/(abs(d13c_sp(isp)-d13c_ocn)+abs(d13c_sp(isp+1)-d13c_ocn))
    ! ccflx(isp+1) = ccflxi*abs((d13c_sp(isp)-d13c_ocn))/(abs(d13c_sp(isp)-d13c_ocn)+abs(d13c_sp(isp+1)-d13c_ocn))
    flxfrc(1) = abs(d13c_ocnf-d13c_ocn)/(abs(d13c_ocnf-d13c_ocn)+abs(d13c_ocni-d13c_ocn))
    flxfrc(2) = abs(d13c_ocni-d13c_ocn)/(abs(d13c_ocnf-d13c_ocn)+abs(d13c_ocni-d13c_ocn))
    flxfrc(3) = abs(d18o_ocnf-d18o_ocn)/(abs(d18o_ocnf-d18o_ocn)+abs(d18o_ocni-d18o_ocn))
    flxfrc(4) = abs(d18o_ocni-d18o_ocn)/(abs(d18o_ocnf-d18o_ocn)+abs(d18o_ocni-d18o_ocn))
    do 
        ! case flxfrc2(1)=0
        flxfrc2=0d0
        flxfrc2(2) = flxfrc(1)
        flxfrc2(3) = flxfrc(3)
        flxfrc2(4) = flxfrc(2)-flxfrc2(3)
        if (all(flxfrc2>=0d0))exit 
        flxfrc2=0d0
        ! case flxfrc2(2)=0
        flxfrc2(1) = flxfrc(1)
        flxfrc2(3) = flxfrc(3)-flxfrc2(1)
        flxfrc2(4) = flxfrc(4)
        if (all(flxfrc2>=0d0))exit 
        flxfrc2=0d0
        ! case flxfrc2(3)=0
        flxfrc2(1) = flxfrc(3)
        flxfrc2(2) = flxfrc(1)-flxfrc2(1)
        flxfrc2(4) = flxfrc(2)
        if (all(flxfrc2>=0d0))exit 
        flxfrc2=0d0
        ! case flxfrc2(4)=0
        flxfrc2(2) = flxfrc(4)
        flxfrc2(1) = flxfrc(1)-flxfrc2(2)
        flxfrc2(3) = flxfrc(2)
        if (all(flxfrc2>=0d0))exit 
        print*,'error'
        stop
    enddo 
#ifdef track2    
    if (mod(it,int(5000/(nspcc-2)))==0) then  
        cntsp=cntsp+1
        d18o_sp(cntsp)=d18o_ocn
        d13c_sp(cntsp)=d13c_ocn
        ccflx = 0d0
        ccflx(cntsp) = ccflxi
    endif 
#else 
    do isp=1,nspcc
        ccflx(isp)=flxfrc2(isp)*ccflxi
    enddo
#endif 
    
#ifdef size
    do isp=1,4
        ccflx(isp)=flxfrc2(isp)*ccflxi*flxfin
    enddo
    do isp=5,8
        ccflx(isp)=flxfrc2(isp-4)*ccflxi*(1d0-flxfin)
    enddo
#endif 
    
    ! dep = 4.0d0 + (6d0-4d0)*(time-time_spn)/time_trs
    ! if (time-time_spn<=time_trs/2d0) then
        ! dep = 4.0d0 + (6d0-4d0)*(time-time_spn)/time_trs*2d0
    ! else 
        ! dep = 2*6d0 - 4.0d0 - (6d0-4d0)*(time-time_spn)/time_trs*2d0
    ! endif
    ! dep = 4.0d0
    ! dep = 6.0d0
elseif (time>time_spn+time_trs) then 
    dt = time_aft/1000d0
    dt=time_trs/1000d0 
    d13c_ocn = d13c_ocni 
    d18o_ocn = d18o_ocni
    ccflx = 0d0
    ! ccflx(3) = ccflxi
    ccflx(1) = ccflxi
    ! dep = 6.0d0
#ifdef track2
    if (cntsp+1/=nspcc) then
        print*,'fatal error in counting',cntsp
        stop
    endif 
    d18o_sp(cntsp+1)=d18o_ocn
    d13c_sp(cntsp+1)=d13c_ocn
    ccflx = 0d0
    ccflx(cntsp+1) = ccflxi
#endif 
#ifdef size 
    ccflx = 0d0
    flxfin=flxfini 
    ccflx(1) = flxfin*ccflxi
    ccflx(5) = (1d0-flxfin)*ccflxi
#endif 
    ! d13c_ocn = d13c_ocni 
    ! d18o_ocn = d18o_ocni
    ! ccflx = 0d0
    ! ccflx(1) = ccflxi
    dep = depi
endif
#endif 

d18o_flx = sum(d18o_sp(:)*ccflx(:))/ccflxi
d13c_flx = sum(d13c_sp(:)*ccflx(:))/ccflxi

#ifndef track2
if (abs(d13c_flx - d13c_ocn)>tol .or. abs(d18o_flx - d18o_ocn)>tol) then 
    print*,'error in assignment of proxy'
    write(file_err,*)'error in assignment of proxy',d18o_ocn,d13c_ocn,d18o_flx,d13c_flx
    stop
endif 
#endif
! dt = 50d0  ! when difficult to converge with non-local mixing

!! === temperature & pressure and associated boundary changes ====

dif_dic0 = (151.69d0 + 7.93d0*temp) ! cm2/yr at 2 oC (Huelse et al. 2018)
dif_alk0 = (151.69d0 + 7.93d0*temp) ! cm2/yr  (Huelse et al. 2018)
dif_o20 =  (348.62d0 + 14.09d0*temp) 
dif_dic = dif_dic0*ff
dif_alk = dif_alk0*ff
dif_o2 = dif_o20*ff

keqcc = calceqcc(temp,sal,dep)
co3sat = keqcc/cai
!! /////////////////////

#ifndef nonrec 
if (it==1) then  
    write(file_bound,*) '#time  d13c_ocn  d18o_ocn, fluxes of cc:',(isp,isp=1,nspcc),'temp  dep  sal  dici  alki  o2i'
endif 
#ifndef size 
write(file_bound,*) time, d13c_ocn, d18o_ocn, (ccflx(isp),isp=1,nspcc),temp, dep, sal,dici,alki, o2i
#else 
write(file_bound,*) time, d13c_ocn, d18o_ocn, sum(ccflx(1:4)),sum(ccflx(5:8)),(ccflx(isp),isp=1,nspcc),temp, dep, sal,dici,alki, o2i
#endif 
#endif 

itr_w = 0 
err_w_min = 1d4

itr_f = 0
err_f_min = 1d4
dfrt_df = 0d0
d2frt_df2 = 0d0

err_f = 0d0
err_fx = 0d0

300 continue 

#ifndef nondisp
print*,'it :',it,dt
#endif

dw = 0d0
oxco2 = 0d0
anco2 = 0d0

itr = 0
error = 1d4
minerr= 1d4
! ~~~~~~~~~ OM & O2 iteration wrt zox ~~~~~~~~~~~~~~~
do while (error > tol)

!~~~~~~ OM calculation ~~~~~~~~~~~~~~~~~

if (oxic) then 
    do iz=1,nz
        ! if (z(iz) < zox) 
        if (o2x(iz) > o2th) then
            kom(iz) = komi
            ! dbio(iz) =  0.15d0
            izox = iz
        else
            kom(iz) = 0d0
            if (anoxic) kom(iz) = komi
            ! dbio(iz) =  0d0
        endif
    enddo
    
! NOTE: I first assumed that mixed layer = oxic zone. The above and following commented out parts are thus implemented.
!       After checking Archer's code, he assumed a fixed mixed layer depth of 12 cm despite oxygen conditions. 
!       Accordingly, dio-disturbances are fixed all the time. 
    
    ! transdbio = 0d0
    ! do iz = 1, izox
        ! if (iz==1) then 
            ! transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            ! transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        ! elseif (iz==izox) then 
            ! transdbio(iz,iz) = 0.5d0*(sporo(Iz)*dbio(iz)+sporo(Iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            ! transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
        ! else 
            ! transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))  &
                ! + 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            ! transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            ! transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        ! endif
    ! enddo
    ! trans = transdbio
    ! if (nonlocal) trans = translabs

endif
    
nsp=1
nmx = nz*nsp
if (allocated(amx)) deallocate(amx)
if (allocated(ymx)) deallocate(ymx)
if (allocated(emx)) deallocate(emx)
if (allocated(ipiv)) deallocate(ipiv)
allocate(amx(nmx,nmx),ymx(nmx),emx(nmx),ipiv(nmx))

amx = 0d0
ymx = 0d0

do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then 
        ymx(row) = &
            + sporo(iz)*(-om(iz))/dt &
            - omflx/dz(1)
            ! - (1d0-poro(iz))*((dbio(iz)+dbio(iz+1))*0.5d0*(omx(iz+1)-omx(iz))/(0.5d0*(dz(1)+dz(2)))  &
            ! - 0d0  &  !  no bioturbation loss at the top boundary 
            ! - (dbio(iz)+dbio(iz))*0.5d0*(omx(iz)-0d0)/(0.5d0*(dz(1)+dz(1)))  &  !  assumes loss via bioturbation; this makes it difficult to converge 
            ! )/dz(1)  &
            ! + (1d0-poro(iz))*w(iz)*(omx(iz)-0d0)/dz(1)  &
            ! + kom(iz)*omx(iz)
        amx(row,row) = (&
            + sporo(iz)*(1d0)/dt &
            ! - (1d0-poro(iz))*((dbio(iz)+dbio(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(1)+dz(2)))-0d0)/dz(1)  &
            ! + w(iz)*(1d0-0d0)/dz(1)   &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-sporoi*wi*0d0)/dz(1)   &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz)*w(iz)*1d0)/dz(1)   &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporoi*wi*0d0)/dz(1)   &
            + sporo(iz)*kom(iz)   &
            ) 
        amx(row,row+nsp) =  (&
            ! - (1d0-poro(iz))*((dbio(iz)+dbio(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(1)+dz(2)))-0d0)/dz(1)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz)*w(iz)*0d0)/dz(1)   &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporoi*wi*0d0)/dz(1)   &
            )
    else if (iz == nz) then 
        ymx(row) = 0d0   &
            + sporo(iz)*(-om(iz))/dt 
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(omx(iz)-omx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! + w(iz)*(omx(iz)-omx(iz-1))/dz(iz)&
            ! + kom(iz)*omx(iz)
        amx(row,row) = (&
            + sporo(iz)*(1d0)/dt &
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! + w(iz)*(1d0)/dz(iz)  &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(sporof*w(iz)*1d0-sporo(iz)*w(iz)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporof*w(iz)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  &
            + sporo(iz)*kom(iz)   &
            )
        amx(row,row-nsp) = ( &
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! + w(iz)*(-1d0)/dz(iz)  &
            + adf(iz)*up(iz)*(sporof*w(iz)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporof*w(iz)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            )
    else 
        ymx(row) = 0d0  &
            + sporo(iz)*(0d0-om(iz))/dt 
            ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(omx(iz+1)-omx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(omx(iz)-omx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! + w(iz)*(omx(iz)-omx(iz-1))/dz(iz)&
            ! + kom(iz)*omx(iz)
        amx(row,row) = (&
            + sporo(Iz)*(1d0)/dt &
            ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(-1d0)/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! + w(iz)*(1d0)/dz(iz)  &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz)*w(iz)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  &
            + sporo(iz)*kom(iz)   &
            )
        amx(row,row+nsp) =  (&
            ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(1d0)/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - 0d0)/dz(iz)  & 
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz)*w(iz)*0d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  &
            )
        amx(row,row-nsp) =  (&
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! + w(iz)*(-1d0)/dz(iz)  &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            )
    endif
    ! if (nonlocal(1)) then 
    if (turbo2(1).or.labs(1)) then 
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp 
            if (trans(iiz,iz,1)==0d0) cycle
            ! amx(row,col) = amx(row,col) -trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))
            amx(row,col) = amx(row,col) -trans(iiz,iz,1)/dz(iz)*dz(iiz)*(1d0-poro(iiz))
            ! amx(row,col) = amx(row,col) -trans(iiz,iz,1)/dz(iz)*dz(iiz)
        enddo
    else 
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp 
            if (trans(iiz,iz,1)==0d0) cycle
            amx(row,col) = amx(row,col) -trans(iiz,iz,1)/dz(iz)
        enddo
    endif
enddo

ymx = - ymx

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls) 

omx = ymx

omadv = 0d0
omdec = 0d0
omdif = 0d0
omrain = 0d0
omtflx = 0d0

do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then 
        omtflx = omtflx + sporo(iz)*(omx(iz)-om(iz))/dt*dz(iz) 
        omadv = omadv + adf(iz)*up(iz)*(sporo(iz)*w(iz)*omx(iz)-0d0)/dz(iz)*dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*omx(iz+1)-sporo(iz)*w(iz)*omx(iz))/dz(iz)*dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*omx(iz+1)-0d0)/dz(iz)*dz(iz)  
            ! + ((1d0-poro(iz))*w(iz)-(1d0-poroi)*wi)*om(iz)/dz(iz)*dz(Iz)  
        ! omdif = omdif &
            ! - (1d0-poro(iz))*((dbio(iz)+dbio(iz+1))*0.5d0*(omx(Iz+1)-omx(iz))/(0.5d0*(dz(1)+dz(2)))-0d0)/dz(1)*dz(iz)
        omdec = omdec + sporo(iz)*kom(iz)*omx(iz)*dz(iz)
        omrain = omrain - omflx/dz(1)*dz(iz)
    else if (iz == nz) then 
        omtflx = omtflx + sporo(iz)*(omx(iz)-om(iz))/dt*dz(iz)
        omadv = omadv + adf(iz)*up(iz)*(sporo(iz)*w(iz)*omx(iz)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz) &
            + adf(iz)*dwn(iz)*(sporof*w(iz)*omx(iz)-sporo(iz)*w(iz)*omx(iz))/dz(iz)*dz(iz) &
            + adf(iz)*cnr(iz)*(sporof*w(iz)*omx(iz)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz) 
            ! + ((1d0-poro(iz))*w(iz)-(1d0-poro(iz-1))*w(Iz-1))*omx(iz)/dz(iz)*dz(Iz)
        ! omdif = omdif &
            ! - (1d0-poro(iz))*(0d0  &
            ! -(dbio(iz)+dbio(iz-1))*0.5d0*(omx(Iz)-omx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)*dz(iz)
        omdec = omdec + sporo(iz)*kom(iz)*omx(iz)*dz(iz)
    else 
        omtflx = omtflx + sporo(iz)*(omx(iz)-om(iz))/dt*dz(iz)
        omadv = omadv + adf(iz)*up(iz)*(sporo(iz)*w(iz)*omx(iz)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*omx(iz+1)-sporo(iz)*w(iz)*omx(iz))/dz(iz)*dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*omx(iz+1)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz) 
        ! omdif = omdif &
            ! - (1d0-poro(iz))*((dbio(iz)+dbio(iz+1))*0.5d0*(omx(Iz+1)-omx(iz))/(0.5d0*(dz(iz)+dz(iz+1)))  &
            ! -(dbio(iz)+dbio(iz-1))*0.5d0*(omx(Iz)-omx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)*dz(iz)
        omdec = omdec + sporo(iz)*kom(iz)*omx(iz)*dz(iz)
    endif
    ! if (nonlocal(1)) then 
    if (turbo2(1).or.labs(1)) then 
        do iiz = 1, nz
            if (trans(iiz,iz,1)==0d0) cycle
            ! omdif = omdif -trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))*dz(iz)*omx(iiz)
            omdif = omdif -trans(iiz,iz,1)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*dz(iz)*omx(iiz)
            ! omdif = omdif -trans(iiz,iz,1)/dz(iz)*dz(iiz)*dz(iz)*omx(iiz)
        enddo
    else
        do iiz = 1, nz
            if (trans(iiz,iz,1)==0d0) cycle
            omdif = omdif -trans(iiz,iz,1)/dz(iz)*dz(iz)*omx(iiz)  ! check previous versions 
        enddo
    endif
enddo

omres = omadv + omdec + omdif + omrain + omtflx 

if (any(omx<0d0)) then
    print*,'negative om, stop'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'NEGATIVE_OM.txt',status = 'unknown')
    do iz = 1, nz
        write (file_tmp,*) z(iz),omx(iz)*mom/rho(iz)*100d0,w(iz),up(iz),dwn(iz),cnr(iz),adf(iz)
    enddo
    close(file_tmp)
    stop
endif 
if (any(isnan(omx))) then
    print*,'nan om, stop'
    print*,omx
    stop
endif 

if (it==1) write(file_omflx,*)'time, omtflx, omadv, omdec, omdif, omrain, omres'
write(file_omflx,*)time, omtflx, omadv, omdec, omdif, omrain, omres

! stop

!~~~~~~~~~~~~~~~~~ O2 calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

amx = 0d0
ymx = 0d0

if (izox == nz) then 

do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then 
        ymx(row) = ( &
            + poro(iz)*(0d0-o2(iz))/dt & 
            - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(0d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_o2(iz)*(0d0-o2i*1d-6/1d3)/dz(iz))/dz(iz)  &
            + sporo(iz)*ox2om*kom(iz)*omx(iz)  &
            )
        amx(row,row) = (& 
            + poro(iz)*(1d0)/dt & 
            - ((poro(iz)*dif_o2(iz)+poro(Iz+1)*dif_o2(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_o2(iz)*(1d0)/dz(iz))/dz(iz)&
            )
        amx(row,row+nsp) = (& 
            - ((poro(Iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - 0d0)/dz(iz)&
            )
    else if (iz == nz) then 
        ymx(row) = (0d0 & 
            + poro(iz)*(0d0-o2(iz))/dt &
            - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(Iz-1))*(0d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            + sporo(iz)*ox2om*kom(iz)*omx(iz)  &
            )
        amx(row,row) = ( & 
            + poro(iz)*(1d0)/dt &
            - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(Iz-1)*dif_o2(Iz-1))*(1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            )
        amx(row,row-nsp) = ( & 
            - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(Iz-1))*(-1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            ) 
    else 
        ymx(row) = ( 0d0& 
            + poro(iz)*(0d0-o2(iz))/dt & 
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(0d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(0d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            + sporo(iz)*ox2om*kom(iz)*omx(iz)  &
            )
        amx(row,row) = (& 
            + poro(iz)*(1d0)/dt & 
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(-1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            )
        amx(row,row+nsp) = (& 
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0d0)/dz(iz)  &
            )
        amx(row,row-nsp) = (& 
            - (0d0 &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) &
            )
    endif
enddo

ymx = - ymx

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls) 

o2x = ymx

o2dec = 0d0 
o2dif = 0d0
o2tflx = 0d0

do iz = 1,nz 
    if (iz == 1) then 
        o2dec = o2dec + sporo(iz)*ox2om*kom(iz)*omx(iz)*dz(iz)
        o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz)
        o2dif = o2dif - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(o2x(iz+1)-o2x(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_o2(iz)*(o2x(iz)-o2i*1d-6/1d3)/dz(iz))/dz(iz) *dz(iz)
    else if (iz == nz) then 
        o2dec = o2dec + (1d0-poro(iz))*ox2om*kom(iz)*omx(iz)/poro(iz)*dz(iz)*poro(iz)
        o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz)
        o2dif = o2dif & 
            - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(Iz-1))*(o2x(iz)-o2x(iz-1))/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            *dz(iz)
    else 
        o2dec = o2dec + (1d0-poro(iz))*ox2om*kom(iz)*omx(iz)/poro(iz)*dz(iz)*poro(iz)
        o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz)
        o2dif = o2dif &
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(o2x(iz+1)-o2x(iz))/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(o2x(Iz)-o2x(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            *dz(iz)
    endif
enddo

else 

do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then 
        ymx(row) = ( &
            + poro(iz)*(0d0-o2(iz))/dt & 
            - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(0d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_o2(iz)*(0d0-o2i*1d-6/1d3)/dz(iz))/dz(iz)  &
            + sporo(iz)*ox2om*kom(iz)*omx(iz)  &
            )
        amx(row,row) = (& 
            + poro(iz)*(1d0)/dt & 
            - ((poro(iz)*dif_o2(iz)+poro(Iz+1)*dif_o2(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_o2(iz)*(1d0)/dz(iz))/dz(iz)&
            )
        amx(row,row+nsp) = (& 
            - ((poro(Iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - 0d0)/dz(iz)&
            )
    ! else if (iz == izox) then 
        ! ymx(row) = (0d0 & 
            ! + poro(iz)*(0d0-o2(iz))/dt &
            ! - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(Iz-1))*(0d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            ! + sporo(iz)*ox2om*kom(iz)*omx(iz)  &
            ! )
        ! amx(row,row) = ( & 
            ! + poro(iz)*(1d0)/dt &
            ! - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(Iz-1)*dif_o2(Iz-1))*(1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            ! )
        ! amx(row,row-nsp) = ( & 
            ! - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(Iz-1))*(-1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            ! ) 
    else if (iz>1 .and. iz<= izox) then 
        ymx(row) = ( 0d0& 
            + poro(iz)*(0d0-o2(iz))/dt & 
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(0d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(0d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            + sporo(iz)*ox2om*kom(iz)*omx(iz)  &
            )
        amx(row,row) = (& 
            + poro(iz)*(1d0)/dt & 
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(-1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            )
        amx(row,row+nsp) = (& 
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0d0)/dz(iz)  &
            )
        amx(row,row-nsp) = (& 
            - (0d0 &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) &
            )
    else if (iz> izox) then  
        ymx(row) = ( 0d0& 
            )
        amx(row,row) = (& 
            + 1d0 &
            )
    endif
enddo

ymx = - ymx

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls) 

o2x = ymx

o2dec = 0d0 
o2dif = 0d0
o2tflx = 0d0

do iz = 1,nz 
    if (iz == 1) then 
        o2dec = o2dec + sporo(iz)*ox2om*kom(iz)*omx(iz)*dz(iz)
        o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz)
        o2dif = o2dif - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(o2x(iz+1)-o2x(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_o2(iz)*(o2x(iz)-o2i*1d-6/1d3)/dz(iz))/dz(iz) *dz(iz)
    ! else if (iz == izox) then 
        ! o2dec = o2dec + (1d0-poro(iz))*ox2om*kom(iz)*omx(iz)/poro(iz)*dz(iz)*poro(iz)
        ! o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz)
        ! o2dif = o2dif & 
            ! - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(Iz-1))*(o2x(iz)-o2x(iz-1))/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            ! *dz(iz)
    else if (iz>1 .and. iz<=izox) then 
        o2dec = o2dec + (1d0-poro(iz))*ox2om*kom(iz)*omx(iz)/poro(iz)*dz(iz)*poro(iz)
        o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz)
        o2dif = o2dif &
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(o2x(iz+1)-o2x(iz))/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(o2x(Iz)-o2x(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            *dz(iz)
    endif
enddo

endif

o2res = o2dec + o2dif + o2tflx 

if (it==1) write(file_o2flx,*)'time, o2dec, o2dif, o2tflx, o2res'
write(file_o2flx,*)time,o2dec, o2dif,o2tflx,o2res

zoxx = 0d0
do iz=1,nz
    if (o2x(iz)<=0d0) exit
enddo

if (iz==nz+1) then 
    zoxx = ztot 
else if (iz==1) then 
    zoxx = (z(iz)*o2i*1d-6/1d3 + 0d0*abs(o2x(iz)))/(o2i*1d-6/1d3+abs(o2x(iz)))
else     
    zoxx = (z(iz)*o2x(iz-1) + z(iz-1)*abs(o2x(iz)))/(o2x(iz-1)+abs(o2x(iz)))
endif

error = abs((zox -zoxx)/zox)
#ifdef showiter 
print*, 'zox',itr,zox, zoxx, error
#endif
if (zox==zoxx) exit 
 
zox = 0.5d0*(zox + zoxx)
! zox = zoxx

if (itr>=100 .and. itr <= nz+99) then 
    zox = z(itr-99) ! next test 
    if (minerr >=error ) then ! if this time is good adopt as optimum 
        if (itr/=100) then 
            izox_minerr = itr -100
            minerr = error 
        endif 
    endif
    ! print*,izox_minerr, minerr ,error
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

! print*,omx
! print*,o2x
! stop
! pause

!~~  O2 calculation END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

! if (.not.anoxic) anco2 = 0d0

! if (oxic) then 
    ! do iz=1,nz 
        ! if (o2x(iz)<=o2th) exit 
    ! enddo
    ! if (iz==nz+1) iz=nz
    ! zox=z(iz)
    ! izox = iz
    ! anco2 = 0d0
! endif

do iz=1,nz
    dw(iz) = dw(iz) -(1d0-poro(iz))*mvom*kom(iz)*omx(iz)
    
    ! if (nonlocal(1)) then 
    if (turbo2(1).or.labs(1)) then 
        do iiz = 1, nz
            if (trans(iiz,iz,1)==0d0) cycle
            ! dw(iz) = dw(iz) - mvom*(-trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))*omx(iiz))
            dw(iz) = dw(iz) - mvom*(-trans(iiz,iz,1)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*omx(iiz))
            ! dw(iz) = dw(iz) - mvom*(-trans(iiz,iz,1)/dz(iz)*dz(iiz)*omx(iiz))
        enddo
    ! endif
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
    if (omx(iz)<1d-30) omx(iz)=1d-30  !! minimum value 
enddo

! print*,oxco2
! print*,anco2

! oxco2 = oxco2*0.1d0

! pause

! stop 

!!  ~~~~~~~~~~~~~~~~~~~~~~ CaCO3 solid, ALK and DIC  calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! do it = 1,nt

error = 1d4
itr = 0

! dt = 0.25d0
! dt = dtcc
! dt = 1d100

nsp = 2 + nspcc
nmx = nz*nsp
deallocate(amx,ymx,emx,ipiv)
allocate(amx(nmx,nmx),ymx(nmx),emx(nmx),ipiv(nmx))

if (allocated(dumx))deallocate(dumx)
allocate(dumx(nmx,nmx))

do while (error > tol)
    
amx = 0d0
ymx = 0d0

call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,infosbr)
if (infosbr==1) then 
    dt=dt/10d0
    ! go to 300
#ifdef sense
    go to 500
#else
    stop
#endif 
endif 
call calcdevs(dicx,alkx,temp,sal,dep,nz,infosbr,dco3_dalk,dco3_ddic)
if (infosbr==1) then 
    dt=dt/10d0
    ! go to 300
#ifdef sense
    go to 500
#else
    stop
#endif 
endif 
! print*, dco3_dalk
! stop
do isp=1,nspcc
    rcc(:,isp) = kcc(:,isp)*ccx(:,isp)*abs(1d0-co3x(:)*1d3/co3sat)**ncc*merge(1d0,0d0,(1d0-co3x(:)*1d3/co3sat)>0d0)
    drcc_dcc(:,isp) = kcc(:,isp)*abs(1d0-co3x(:)*1d3/co3sat)**ncc*merge(1d0,0d0,(1d0-co3x(:)*1d3/co3sat)>0d0)
    drcc_dco3(:,isp) = kcc(:,isp)*ccx(:,isp)*ncc*abs(1d0-co3x(:)*1d3/co3sat)**(ncc-1d0)  &
        *merge(1d0,0d0,(1d0-co3x(:)*1d3/co3sat)>0d0)*(-1d3/co3sat)
    drcc_ddic(:,isp) = drcc_dco3(:,isp)*dco3_ddic(:)
    drcc_dalk(:,isp) = drcc_dco3(:,isp)*dco3_dalk(:)
enddo

do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then 
        do isp = 1,nspcc  ! multiple caco3 species 
            ymx(row+isp-1) = &
                + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt &
                - ccflx(isp)/dz(1) &
                ! - ((dbio(iz)+dbio(iz+1))*0.5d0*(ccx(iz+1)-ccx(iz))/(0.5d0*(dz(1)+dz(2)))  &
                ! - 0d0  &  !  no bioturbation loss at the top boundary  
                ! )/dz(1)  &
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-0d0)/dz(1)  &
                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(1)  &
                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-0d0)/dz(1)  &
                + sporo(iz)*rcc(iz,isp)
            amx(row+isp-1,row+isp-1) = (&
                + sporo(iz)*(1d0)/dt &
                ! - ((dbio(iz)+dbio(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(1)+dz(2)))-0d0)/dz(1)  &
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(1)   &
                + adf(iz)*dwn(iz)*(0d0-sporo(iz)*w(iz)*1d0)/dz(1)  &
                + sporo(iz)* drcc_dcc(iz,isp)  &
                )* ccx(iz,isp) 
            amx(row+isp-1,row+isp-1+nsp) =  (&
                ! - ((dbio(iz)+dbio(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(1)+dz(2)))-0d0)/dz(1)  &
                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(1)  &
                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(1)  &
                )*ccx(iz+1,isp)
            amx(row+isp-1,row+nspcc) = (&
                + sporo(iz)*drcc_ddic(iz,isp)  &
                )*dicx(iz)
            amx(row+isp-1,row+nspcc+1) = (&
                + sporo(iz)*drcc_dalk(iz,isp)  &
                )*alkx(iz)
            ! DIC
            amx(row+nspcc,row+isp-1) = (&
                - (1d0-poro(Iz))*drcc_dcc(iz,isp)  &
                )*ccx(iz,isp)*fact
            ! ALK 
            amx(row+nspcc+1,row+isp-1) = (&
                - 2d0* (1d0-poro(Iz))*drcc_dcc(iz,isp)  &
                )*ccx(iz,isp)*fact
        enddo 
        !  DIC 
        ymx(row+nspcc) = ( &
            + poro(iz)*(dicx(iz)-dic(iz))/dt & 
            - ((poro(iz)*dif_dic(iz)+poro(iz+1)*dif_dic(iz+1))*0.5d0*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_dic(iz)*(dicx(iz)-dici*1d-6/1d3)/dz(iz))/dz(iz)  &
            - oxco2(iz) &
            - anco2(iz) &
            - (1d0-poro(Iz))*sum(rcc(iz,:))  &
            )*fact
        amx(row+nspcc,row+nspcc) = (& 
            + poro(iz)*(1d0)/dt & 
            - ((poro(iz)*dif_dic(iz)+poro(Iz+1)*dif_dic(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(Iz)*dif_dic(iz)*(1d0)/dz(iz))/dz(iz)&
            - (1d0-poro(Iz))*sum(drcc_ddic(iz,:))  &
            )*dicx(iz)*fact
        amx(row+nspcc,row+nspcc+nsp) = (& 
            - ((poro(iz)*dif_dic(iz)+poro(Iz+1)*dif_dic(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - 0d0)/dz(iz)&
            )*dicx(iz+1)*fact
        amx(row+nspcc,row+nspcc+1) = ( &
            - (1d0-poro(Iz))*sum(drcc_dalk(iz,:))  &
            )*alkx(iz)*fact
        ! ALK
        ymx(row+nspcc+1) = (& 
            + poro(iz)*(alkx(iz)-alk(iz))/dt & 
            - ((poro(iz)*dif_alk(iz)+poro(Iz+1)*dif_alk(iz+1))*0.5d0*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_alk(iz)*(alkx(iz)-alki*1d-6/1d3)/dz(iz))/dz(iz) &
            - anco2(iz) &
            - 2d0* (1d0-poro(Iz))*sum(rcc(iz,:))  &
            )*fact
        amx(row+nspcc+1,row+nspcc+1) = (& 
            + poro(iz)*(1d0)/dt & 
            - ((poro(iz)*dif_alk(iz)+poro(iz+1)*dif_alk(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_alk(iz)*(1d0)/dz(iz))/dz(iz)  &
            - 2d0* (1d0-poro(Iz))*sum(drcc_dalk(iz,:))  &
            )*alkx(iz)*fact
        amx(row+nspcc+1,row+nspcc+1+nsp) = (& 
            - ((poro(Iz)*dif_alk(iz)+poro(Iz+1)*dif_alk(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - 0d0)/dz(iz)&
            )*alkx(iz+1)*fact
        amx(row+nspcc+1,row+nspcc) = (&
            - 2d0* (1d0-poro(Iz))*sum(drcc_ddic(iz,:))  &
            )*dicx(iz)*fact
    else if (iz == nz) then 
        do isp=1,nspcc
            ymx(row+isp-1) = & 
                + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt &
                ! - (0d0 &
                ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ccx(iz)-ccx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz)  &
                + adf(iz)*cnr(iz)*(sporof*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz)  &
                + adf(iz)*dwn(iz)*(sporof*w(iz)*ccx(iz,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(iz)  &
                + sporo(iz)*rcc(iz,isp)
            amx(row+isp-1,row+isp-1) = (&
                + sporo(iz)*(1d0)/dt &
                ! - (0d0 &
                ! - (dbio(iz)+dbio(iz-1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  &
                + adf(iz)*cnr(iz)*(sporof*w(iz)*1d0-0d0)/dz(iz)  &
                + adf(iz)*dwn(iz)*(sporof*w(iz)*1d0-sporo(iz)*w(iz)*1d0)/dz(iz)  &
                + sporo(iz)*drcc_dcc(iz,isp)   &
                )*ccx(iz,isp)
            amx(row+isp-1,row+isp-1-nsp) = ( &
                ! - (0d0 &
                ! - (dbio(iz)+dbio(iz-1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
                + adf(iz)*up(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
                + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
                )*ccx(iz-1,isp)
            amx(row+isp-1,row+nspcc) = (&
                + sporo(iz)*drcc_ddic(iz,isp) &
                )*dicx(iz)
            amx(row+isp-1,row+nspcc+1) = (&
                + sporo(iz)*drcc_dalk(iz,isp) &
                )*alkx(iz)
            ! ymx(row) = ccx(iz) - ccx(iz-1)
            ! amx(row, row) = 1d0*ccx(iz)
            ! amx(row,row-1)= -1d0*ccx(iz-1)
            
            !DIC 
            amx(row+nspcc,row+isp-1) = (&
                - sporo(Iz)*drcc_dcc(iz,isp)  &
                )*ccx(iz,isp)*fact
            !ALK 
            amx(row+nspcc+1,row+isp-1) = (&
                - 2d0*sporo(Iz)*drcc_dcc(iz,isp)  &
                )*ccx(Iz,isp)*fact
        enddo
        ! DIC
        ymx(row+nspcc) = (& 
            + poro(iz)*(dicx(iz)-dic(iz))/dt &
            - (0d0 - 0.5d0*(poro(iz)*dif_dic(iz)+poro(Iz-1)*dif_dic(Iz-1))*(dicx(iz)-dicx(iz-1))  &
                /(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            - oxco2(iz) &
            - anco2(iz) &
            - sporo(iz)*sum(rcc(iz,:))  &
            )*fact
        amx(row+nspcc,row+nspcc) = ( & 
            + poro(iz)*(1d0)/dt &
            - (0d0 - 0.5d0*(poro(iz)*dif_dic(iz)+poro(Iz-1)*dif_dic(Iz-1))*(1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            - sporo(Iz)*sum(drcc_ddic(iz,:))  &
            )*dicx(iz)*fact
        amx(row+nspcc,row+nspcc-nsp) = ( & 
            - (0d0 - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(Iz-1))*(-1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            ) * dicx(iz-1)*fact
        amx(row+nspcc,row+nspcc+1) = (&
            - sporo(Iz)*sum(drcc_dalk(iz,:))  &
            )*alkx(iz)*fact
        ! ALK 
        ymx(row+nspcc+1) = ( & 
            + poro(iz)*(alkx(iz)-alk(iz))/dt &
            - (0d0 - 0.5d0*(poro(iz)*dif_alk(iz)+poro(Iz-1)*dif_alk(Iz-1))*(alkx(iz)-alkx(iz-1))/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            - anco2(iz) &
            - 2d0*sporo(Iz)*sum(rcc(iz,:))  &
            )*fact
        amx(row+nspcc+1,row+nspcc+1) = ( & 
            + poro(iz)*(1d0)/dt &
            - (0d0 - 0.5d0*(poro(Iz)*dif_alk(iz)+poro(iz-1)*dif_alk(Iz-1))*(1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            - 2d0*sporo(Iz)*sum(drcc_dalk(iz,:))  &
            )*alkx(iz)*fact
        amx(row+nspcc+1,row+nspcc+1-nsp) = ( & 
            - (0d0 - 0.5d0*(poro(iz)*dif_alk(iz)+poro(Iz-1)*dif_alk(Iz-1))*(-1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            ) * alkx(iz-1)*fact
        amx(row+nspcc+1,row+nspcc) = (&
            - 2d0*sporo(Iz)*sum(drcc_ddic(iz,:))  &
            )*dicx(Iz)*fact
    else 
        do isp=1,nspcc
            ymx(row+isp-1) = & 
                + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt &
                ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(ccx(iz+1)-ccx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) &
                ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ccx(iz)-ccx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-sporo(Iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz)  &
                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(Iz)*w(iz)*ccx(iz,isp))/dz(iz)  &
                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(Iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz)  &
                + sporo(iz)*rcc(iz,isp)
            amx(row+isp-1,row+isp-1) = (&
                + sporo(iz)*(1d0)/dt &
                ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(-1d0)/(0.5d0*(dz(iz+1)+dz(iz))) &
                ! - (dbio(iz)+dbio(iz-1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  &
                + adf(iz)*dwn(iz)*(0d0-sporo(Iz)*w(iz)*1d0)/dz(iz)  &
                + sporo(iz)*drcc_dcc(iz,isp)  &
                )*ccx(iz,isp)
            amx(row+isp-1,row+isp-1+nsp) =  (&
                ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(1d0)/(0.5d0*(dz(iz+1)+dz(iz))) &
                ! - 0d0)/dz(iz)  & 
                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(iz)  &
                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(iz)  &
                )*ccx(iz+1,isp)
            amx(row+isp-1,row+isp-1-nsp) =  (&
                ! - (0d0 &
                ! - (dbio(iz)+dbio(iz-1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
                + adf(iz)*up(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
                + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
                )*ccx(iz-1,isp)
            amx(row+isp-1,row+nspcc) = (& 
                + sporo(Iz)*drcc_ddic(iz,isp)  &
                )*dicx(Iz)
            amx(row+isp-1,row+nspcc+1) = (&
                + sporo(Iz)*drcc_dalk(iz,isp) &
                )*alkx(iz)
            ! DIC 
            amx(row+nspcc,row+isp-1) = (&
                - sporo(Iz)*drcc_dcc(iz,isp)  &
                )*ccx(iz,isp)*fact
            ! ALK 
            amx(row+nspcc+1,row+isp-1) = (&
                - 2d0*sporo(Iz)*drcc_dcc(iz,isp)  &
                )*ccx(iz,isp)*fact 
        enddo
        ! DIC 
        ymx(row+nspcc) = ( & 
            + poro(iz)*(dicx(iz)-dic(iz))/dt & 
            - (0.5d0*(poro(iz+1)*dif_dic(iz+1)+poro(Iz)*dif_dic(iz))*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(dicx(Iz)-dicx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            - oxco2(iz) &
            - anco2(iz) &
            - sporo(Iz)*sum(rcc(iz,:))  &
            )*fact
        amx(row+nspcc,row+nspcc) = (& 
            + poro(iz)*(1d0)/dt & 
            - (0.5d0*(poro(iz+1)*dif_dic(iz+1)+poro(iz)*dif_dic(iz))*(-1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            - sporo(iz)*sum(drcc_ddic(iz,:))  &
            )*dicx(iz)*fact
        amx(row+nspcc,row+nspcc+nsp) = (& 
            - (0.5d0*(poro(iz+1)*dif_dic(iz+1)+poro(iz)*dif_dic(iz))*(1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0d0)/dz(iz)  &
            )*dicx(iz+1)*fact
        amx(row+nspcc,row+nspcc-nsp) = (& 
            - (0d0 &
            - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) &
            )*dicx(iz-1)*fact
        amx(row+nspcc,row+nspcc+1) = (&
            - sporo(Iz)*sum(drcc_dalk(iz,:))  &
            )*alkx(iz)*fact
        ! ALK 
        ymx(row+nspcc+1) = (& 
            + poro(iz)*(alkx(iz)-alk(iz))/dt & 
            - (0.5d0*(poro(iz+1)*dif_alk(iz+1)+poro(iz)*dif_alk(iz))*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(Iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(alkx(iz)-alkx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) &
            - anco2(iz) &
            - 2d0*sporo(Iz)*sum(rcc(iz,:))  &
            ) *fact
        amx(row+nspcc+1,row+nspcc+1) = (& 
            + poro(iz)*(1d0)/dt & 
            - (0.5d0*(poro(iz+1)*dif_alk(iz+1)+poro(iz)*dif_alk(iz))*(-1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(Iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            - 2d0*sporo(Iz)*sum(drcc_dalk(iz,:))  &
            )*alkx(iz)*fact
        amx(row+nspcc+1,row+nspcc+1+nsp) = ( & 
            - (0.5d0*(poro(iz+1)*dif_alk(iz+1)+poro(iz)*dif_alk(iz))*(1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0d0)/dz(iz)  &
            )*alkx(iz+1)*fact
        amx(row+nspcc+1,row+nspcc+1-nsp) = (& 
            - (0d0 &
            - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) &
            )*alkx(iz-1)*fact
        amx(row+nspcc+1,row+nspcc) = (&
            - 2d0*sporo(Iz)*sum(drcc_ddic(iz,:))  &
            )*dicx(iz)*fact 
    endif
    do isp=1,nspcc
        ! if (nonlocal(isp+2)) then
        if (turbo2(isp+2).or.labs(isp+2)) then
            do iiz = 1, nz
                col = 1 + (iiz-1)*nsp
                if (trans(iiz,iz,isp+2)==0d0) cycle
                ! amx(row+isp-1,col+isp-1) = amx(row+isp-1,col+isp-1) &
                    ! - trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))*ccx(iiz,isp)
                ! ymx(row+isp-1) = ymx(row+isp-1) &
                    ! - trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))*ccx(iiz,isp)
                amx(row+isp-1,col+isp-1) = amx(row+isp-1,col+isp-1) &
                    - trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*ccx(iiz,isp)
                ymx(row+isp-1) = ymx(row+isp-1) &
                    - trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*ccx(iiz,isp)
                ! amx(row+isp-1,col+isp-1) = amx(row+isp-1,col+isp-1) &
                    ! - trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*ccx(iiz,isp)
                ! ymx(row+isp-1) = ymx(row+isp-1) &
                    ! - trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*ccx(iiz,isp)
            enddo
        else
            do iiz = 1, nz
                col = 1 + (iiz-1)*nsp
                if (trans(iiz,iz,isp+2)==0d0) cycle
                amx(row+isp-1,col+isp-1) = amx(row+isp-1,col+isp-1) -trans(iiz,iz,isp+2)/dz(iz)*ccx(iiz,isp)
                ymx(row+isp-1) = ymx(row+isp-1) - trans(iiz,iz,isp+2)/dz(iz)*ccx(iiz,isp)
            enddo
        endif
    enddo
enddo

ymx = - ymx

#ifndef nonrec
if (any(isnan(ymx))) then 
    print*,'NAN in ymx'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_ymx_pre.txt',status = 'unknown')
    do iz = 1, nmx
        write (file_tmp,*) ymx(iz)
    enddo
    close(file_tmp)
    ! stop
endif
#endif 

! call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls) 

!  slowest way of using sparse matrix solver
n = nmx
where(amx/=0d0)
    dumx=1
elsewhere
    dumx=0
endwhere
! dumx = merge(1,0,amx/=0d0) 
nnz = sum(dumx)

if (allocated(ai)) deallocate(ai)
if (allocated(ap)) deallocate(ap)
if (allocated(ax)) deallocate(ax)
if (allocated(bx)) deallocate(bx)
if (allocated(kai)) deallocate(kai)

allocate(ai(nnz))
allocate(ap(n+1))
allocate(ax(nnz))
allocate(bx(n))
allocate(kai(n))

ai = 0
ap = 0
ax = 0d0
bx = 0d0
kai = 0d0

ap(1)=0
cnt2=0
do i=1,n
    ap(i+1)=ap(i)+sum(dumx(:,i))
    if (ap(i+1)==0) cycle
    cnt=0
    do j=1,n
        if (dumx(j,i)==0) cycle
        cnt=cnt+1
        cnt2=cnt2+1
        ai(cnt2)=j-1
        ax(cnt2)=amx(j,i)
        if (cnt==sum(dumx(:,i)))exit 
    enddo
enddo
if (cnt2/=nnz) then
    print*,'fatal error'
    stop
endif 
bx = ymx
        
! solving matrix with UMFPACK (following is pasted from umfpack_simple.f90)
! Set the default control parameters.
call umf4def( control )
! From the matrix data, create the symbolic factorization information.
call umf4sym ( n, n, ap, ai, ax, symbolic, control, info )
if ( info(1) < 0.0D+00 ) then
    write ( *, * ) ''
    write ( *, *) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, * ) '  UMF4SYM returns INFO(1) = ', info(1)
    stop 1
end if
! From the symbolic factorization information, carry out the numeric factorization.
call umf4num ( ap, ai, ax, symbolic, numeric, control, info )
if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4NUM returns INFO(1) = ', info(1)
    stop 1
end if
!  Free the memory associated with the symbolic factorization.
call umf4fsym ( symbolic )
! Solve the linear system.
sys = 0
call umf4sol ( sys, kai, bx, numeric, control, info )
if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4SOL returns INFO(1) = ', info(1)
    stop 1
end if
! Free the memory associated with the numeric factorization.
call umf4fnum ( numeric )
!  Print the solution.
! write ( *, * ) ''
! write ( *, * ) '  Computed solution'
! write ( *, * ) ''

ymx = kai 

! if (infobls/=0) then 
    ! print*,'nonzero info',infobls
    ! stop
! endif 

! #ifdef nonrec
! if (any(isnan(amx))) then
    ! print*,'NAN in amx'
    ! open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_amx.txt',status = 'unknown')
    ! do iz = 1, nmx
        ! write (file_tmp,*) amx(iz,:)
    ! enddo
    ! close(file_tmp)
    ! stop
! endif

! if (any(isnan(ymx))) then 
    ! print*,'NAN in ymx'
    ! open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_ymx.txt',status = 'unknown')
    ! do iz = 1, nmx
        ! write (file_tmp,*) ymx(iz)
    ! enddo
    ! close(file_tmp)
    ! stop
! endif
! #endif
! stop

do iz = 1, nz 
    row = 1+(iz-1)*nsp
    do isp=1,nspcc
        if (ymx(row+isp-1)>10d0) then 
            ccx(iz,isp) = ccx(iz,isp)*1.5d0
        elseif (ymx(row+isp-1)<-10d0) then 
            ccx(iz,isp) = ccx(iz,isp)*0.5d0
        else
            ccx(iz,isp) = ccx(iz,isp)*exp(ymx(row+isp-1))
        endif
        if (ccx(iz,isp)<1d-30) then
            ccx(iz,isp)=1d-30
            ymx(row+isp-1) = 0d0
        endif
    enddo
    if (ymx(row+nspcc)>10d0) then 
        dicx(iz)=dicx(iz)*1.5d0
    elseif (ymx(row+nspcc)<-10d0) then 
        dicx(iz)=dicx(iz)*0.5d0
    else 
        dicx(iz) = dicx(iz)*exp(ymx(row+nspcc))
    endif
    if (ymx(row+nspcc+1)>10d0) then 
        alkx(Iz) = alkx(iz)*1.5d0
    elseif (ymx(row+nspcc+1)<-10d0) then 
        alkx(iz) = alkx(iz)*0.5d0
    else 
        alkx(iz) = alkx(iz)*exp(ymx(row+nspcc+1))
    endif
    if (dicx(iz)<1d-100) ymx(row+nspcc) = 0d0
    if (alkx(iz)<1d-100) ymx(row+nspcc) = 0d0
    ! print *,exp(ymx(row))
enddo

error = maxval(exp(abs(ymx))) - 1d0
itr = itr + 1
#ifdef showiter
print*,'co2 iteration',itr,error,infobls
#endif
! print*,(ccx(iz),iz=1,nz,20)
! print*,(dicx(iz),iz=1,nz,20)
! print*,(alkx(iz),iz=1,nz,20)

! stop

if (any(ccx<0d0)) then
    print*,'negative ccx, stop'
    print*,ccx
    stop
endif
if (any(isnan(ccx))) then
    print*,'nan om, stop'
    print*,ccx
    stop
endif

if (any(dicx<0d0)) then
    print*,'negative dicx, stop'
    print*,dicx
    stop
endif 
if (any(isnan(dicx))) then
    print*,'nan dic, stop'
    print*,dicx
    stop
endif 

if (any(alkx<0d0)) then
    print*,'negative alk, stop'
    print*,alkx
    stop
endif
if (any(isnan(alkx))) then
    print*,'nan alk, stop'
    print*,alkx
    stop
endif

enddo

! ~~~~  End of calculation iteration for CO2 species ~~~~~~~~~~~~~~~~~~~~

call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,infosbr)
if (infosbr==1) then 
    dt=dt/10d0
    ! go to 300
#ifdef sense
    go to 500
#else
    stop
#endif 
endif 

cctflx =0d0 
ccdis = 0d0 
ccdif = 0d0 
ccadv = 0d0 
ccrain = 0d0
ccres = 0d0 

dictflx = 0d0 
dicdis = 0d0 
dicdif = 0d0 
dicdec = 0d0 
dicres = 0d0

alktflx = 0d0 
alkdis = 0d0 
alkdif = 0d0 
alkdec = 0d0 
alkres = 0d0

do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then 
        do isp=1,nspcc
            cctflx(isp) = cctflx(isp) + (1d0-poro(iz))*(ccx(iz,isp)-cc(iz,isp))/dt *dz(iz)
            ccdis(isp) = ccdis(isp)  + (1d0-poro(Iz))*rcc(iz,isp) *dz(iz)
            ! ccdif = ccdif &
                ! - (1d0-poro(iz))*((dbio(iz)+dbio(iz+1))*0.5d0*(ccx(iz+1)-ccx(iz))/(0.5d0*(dz(1)+dz(2)))  &
                ! - 0d0 )/dz(1)* dz(iz)
            ccrain(isp) = ccrain(isp) - ccflx(isp)/dz(1)*dz(iz)
            ccadv(isp) = ccadv(Isp) + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-0d0)/dz(1) * dz(iz) &
                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(1) * dz(iz)  &
                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-0d0)/dz(1) * dz(iz)
            ! print*, iz,(1d0-poro(iz))*(ccx(iz)-cc(iz))/dt *dz(iz) &
                ! + (1d0-poro(Iz))*rcc(iz) *dz(iz) &
                ! - (1d0-poro(iz))*((dbio(iz)+dbio(iz+1))*0.5d0*(ccx(iz+1)-ccx(iz))/(0.5d0*(dz(1)+dz(2)))  &
                ! - 0d0 )/dz(1)* dz(iz) &
                ! - ccflx/(40d0+60d0)/dz(1)*dz(iz) &
                ! + (1d0-poro(iz))*w(iz)*(ccx(iz)-0d0)/dz(1) * dz(iz)
        enddo
        !  DIC 
        dictflx = dictflx +(dicx(iz)-dic(iz))/dt*dz(iz)*poro(iz) 
        dicdif = dicdif - ((poro(iz)*dif_dic(iz)+poro(iz+1)*dif_dic(iz+1))*0.5d0*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_dic(iz)*(dicx(iz)-dici*1d-6/1d3)/dz(iz))/dz(iz)*dz(iz)
        dicdec = dicdec - oxco2(iz)*dz(iz) - anco2(iz)*dz(iz) 
        dicdis = dicdis - sum(rcc(iz,:))*sporo(iz)*dz(iz) 
        ! ALK
        alktflx = alktflx + (alkx(iz)-alk(iz))/dt*dz(iz)*poro(iz)
        alkdif = alkdif - ((poro(iz)*dif_alk(iz)+poro(iz+1)*dif_alk(iz+1))*0.5d0*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_alk(iz)*(alkx(iz)-alki*1d-6/1d3)/dz(iz))/dz(iz)*dz(iz)
        alkdec = alkdec - anco2(iz)*dz(iz) 
        alkdis = alkdis - 2d0* sporo(Iz)*sum(rcc(iz,:))*dz(iz) 
    else if (iz == nz) then 
        do isp=1,nspcc
            cctflx(isp) = cctflx(isp) + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt *dz(iz)
            ccdis(isp) = ccdis(isp)  + sporo(Iz)*rcc(iz,isp) *dz(iz)
            ! ccdif = ccdif &
                ! - (1d0-poro(iz))*(0d0  &
                ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ccx(iz)-ccx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))) )/dz(iz)* dz(iz)
            ccadv(isp) = ccadv(isp) &
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz) * dz(iz)  &
                + adf(iz)*cnr(iz)*(sporof*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz) * dz(iz)  &
                + adf(iz)*dwn(iz)*(sporof*w(iz)*ccx(iz,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(iz) * dz(iz)  
            ! print*,iz, (1d0-poro(iz))*(ccx(iz)-cc(iz))/dt *dz(iz) &
                ! + (1d0-poro(Iz))*rcc(iz) *dz(iz) &
                ! - (1d0-poro(iz))*(0d0  &
                ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ccx(iz)-ccx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))) )/dz(iz)* dz(iz) &
                ! + (1d0-poro(iz))*w(iz)*(ccx(iz)-ccx(iz-1))/dz(iz) * dz(iz)
        enddo
        ! DIC
        dictflx = dictflx +(dicx(iz)-dic(iz))/dt*dz(iz)*poro(iz) 
        dicdif = dicdif - (0d0 &
            - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(Iz-1))*(dicx(iz)-dicx(iz-1))/(0.5d0*(dz(iz-1)+dz(iz))) &
            )/dz(iz)*dz(iz)
        dicdec = dicdec - oxco2(iz)*dz(iz) - anco2(iz)*dz(iz) 
        dicdis = dicdis - sporo(Iz)*sum(rcc(iz,:))*dz(iz) 
        ! ALK 
        alktflx = alktflx + (alkx(iz)-alk(iz))/dt*dz(iz)*poro(iz)
        alkdif = alkdif - (0d0 &
            - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(Iz-1))*(alkx(iz)-alkx(iz-1))/(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz)*dz(iz)
        alkdec = alkdec - anco2(iz)*dz(iz)
        alkdis = alkdis - 2d0* Sporo(Iz)*sum(rcc(iz,:))*dz(iz)
    else 
        do isp=1,nspcc
            cctflx(isp) = cctflx(isp) + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt *dz(iz)
            ccdis(isp) = ccdis(isp)  + sporo(Iz)*rcc(iz,isp) *dz(iz)
            ! ccdif = ccdif &
                ! - (1d0-poro(iz))*((dbio(iz+1)+dbio(iz))*0.5d0*(ccx(iz+1)-ccx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) &
                ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ccx(iz)-ccx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))) )/dz(iz)* dz(iz)
            ccadv(isp) = ccadv(isp) &
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz) * dz(iz)  &
                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(iz) * dz(iz)  &
                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz) * dz(iz)  
            ! print*,iz,(1d0-poro(iz))*(ccx(iz)-cc(iz))/dt *dz(iz) &
                ! + (1d0-poro(Iz))*rcc(iz) *dz(iz)  &
                ! - (1d0-poro(iz))*((dbio(iz+1)+dbio(iz))*0.5d0*(ccx(iz+1)-ccx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) &
                ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ccx(iz)-ccx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))) )/dz(iz)* dz(iz) &
                ! + (1d0-poro(iz))*w(iz)*(ccx(iz)-ccx(iz-1))/dz(iz) * dz(iz)
        enddo
        ! DIC 
        dictflx = dictflx +(dicx(iz)-dic(iz))/dt*dz(iz)*poro(iz) 
        dicdif = dicdif - (0.5d0*(poro(iz+1)*dif_dic(iz+1)+poro(iz)*dif_dic(iz))*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(Iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(dicx(Iz)-dicx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))) &
            )/dz(iz)*dz(iz)
        dicdec = dicdec - oxco2(iz)*dz(iz) - anco2(iz)*dz(iz) 
        dicdis = dicdis - sporo(Iz)*sum(rcc(iz,:))*dz(iz) 
        ! ALK 
        alktflx = alktflx + (alkx(iz)-alk(iz))/dt*dz(iz)*poro(iz)
        alkdif = alkdif - (0.5d0*(poro(iz+1)*dif_alk(iz+1)+poro(iz)*dif_alk(iz))*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(alkx(iz)-alkx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)*dz(iz)
        alkdec = alkdec - anco2(iz)*dz(iz)
        alkdis = alkdis - 2d0* sporo(iz)*sum(rcc(iz,:))*dz(iz) 
    endif
    do isp=1,nspcc
        ! if (nonlocal(isp+2)) then 
        if (labs(isp+2).or. turbo2(isp+2)) then 
            do iiz = 1, nz
                if (trans(iiz,iz,isp+2)==0d0) cycle
                ! ccdif(isp) = ccdif(isp) -trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))*dz(iz)*ccx(iiz,isp)
                ccdif(isp) = ccdif(isp) -trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*dz(iz)*ccx(iiz,isp)
                ! ccdif(isp) = ccdif(isp) -trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*dz(iz)*ccx(iiz,isp)
            enddo
        else 
            do iiz = 1, nz
                if (trans(iiz,iz,isp+2)==0d0) cycle
                ccdif(isp) = ccdif(isp) -trans(iiz,iz,isp+2)/dz(iz)*dz(iz)*ccx(iiz,isp)
            enddo
        endif
        ! if (nonlocal(isp+2)) then 
        if (labs(isp+2).or. turbo2(isp+2)) then 
            do iiz = 1, nz
                if (trans(iiz,iz,isp+2)==0d0) cycle
                ! dw(iz) = dw(iz) - mvcc*(-trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))*ccx(iiz,isp))
                dw(iz) = dw(iz) - mvcc*(-trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*ccx(iiz,isp))
                ! dw(iz) = dw(iz) - mvcc*(-trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*ccx(iiz,isp))
            enddo
        ! endif
        else 
            if (nonlocal(isp+2)) then 
                do iiz = 1, nz
                    if (trans(iiz,iz,isp+2)==0d0) cycle
                    dw(iz) = dw(iz) -mvcc*(-trans(iiz,iz,isp+2)/dz(iz)*ccx(iiz,isp))
                enddo
            endif 
        endif 
    enddo 
    dw(iz) = dw(iz) -(1d0-poro(iz))*mvcc*sum(rcc(iz,:))
enddo

ccres = cctflx +  ccdis +  ccdif + ccadv + ccrain
dicres = dictflx + dicdis + dicdif + dicdec 
alkres = alktflx + alkdis + alkdif + alkdec 

if (abs(alkres)/max(alktflx,alkdis ,alkdif , alkdec) > tol*10d0) then 
    print*,'not enough accuracy in co2 calc:stop',abs(alkres)/max(alktflx,alkdis ,alkdif , alkdec)
    write(file_err,*)'not enough accuracy in co2 calc:stop',abs(alkres)/max(alktflx,alkdis ,alkdif , alkdec)
    ! stop
endif

#ifndef nonrec
if (it==1) then 
    write(file_ccflx,*) 'time, cctflx, ccdis, ccdif, ccadv, ccrain, ccres' 
    write(file_dicflx,*) 'time, dictflx, dicdis, dicdif, dicdec,  dicres' 
    write(file_alkflx,*) 'time, alktflx, alkdis, alkdif, alkdec, alkres' 
    do isp=1,nspcc
        write(file_ccflxes(isp),*) 'time, cctflx, ccdis, ccdif, ccadv, ccrain, ccres' 
    enddo
endif
write(file_ccflx,*) time,sum(cctflx), sum(ccdis), sum(ccdif), sum(ccadv), sum(ccrain), sum(ccres)
write(file_dicflx,*) time,dictflx, dicdis, dicdif, dicdec,  dicres 
write(file_alkflx,*) time,alktflx, alkdis, alkdif, alkdec, alkres 
do isp=1,nspcc
    write(file_ccflxes(isp),*) time,cctflx(isp), ccdis(isp), ccdif(isp), ccadv(isp), ccrain(isp), ccres(isp)
enddo
#endif
! ~~~~ calculation particle density ~~~~~~~~~~~~~~~~~~
error = 1d4
itr = 0

nsp = 1
nmx = nz*nsp
deallocate(amx,ymx,emx,ipiv)
allocate(amx(nmx,nmx),ymx(nmx),emx(nmx),ipiv(nmx))

! do while (error > tol)
    
amx = 0d0
ymx = 0d0

do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then 
        ymx(row) = &
            + sporo(iz)*(-pt(iz))/dt &
            - detflx/msed/dz(iz)
            ! - ((dbio(iz)+dbio(iz+1))*0.5d0*(ptx(iz+1)-ptx(iz))/(0.5d0*(dz(iz)+dz(iz+1)))  &
            ! - 0d0  &  !  no bioturbation loss at the top boundary 
            ! )/dz(iz)  &
            ! + w(iz)*(ptx(iz)-0d0)/dz(1) 
        amx(row,row) = (&
            + sporo(iz)*(1d0)/dt &
            ! - ((dbio(iz)+dbio(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))-0d0)/dz(iz)  &
            ! - trans(iz,iz)  &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)   &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz)*w(iz)*1d0)/dz(iz)   &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-0d0)/dz(iz)   &
            )            
        amx(row,row+nsp) =  (&
            ! - ((dbio(iz)+dbio(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))-0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz)*w(iz)*0d0)/dz(iz)   &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(iz)   &
            )
        ! do iiz = 1, nz
            ! if (iiz==iz) cycle
            ! col = 1 + (iiz-1)*nsp
            ! amx(row,col) = amx(row,col) -trans(iz,iiz)
        ! enddo
    else if (iz == nz) then 
        ymx(row) = & 
            + sporo(iz)*(-pt(iz))/dt 
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ptx(iz)-ptx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! + w(iz)*(ptx(iz)-ptx(iz-1))/dz(iz)  
        amx(row,row) = (&
            + sporo(iz)*(1d0)/dt &
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! - trans(iz,iz)  &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporof*w(iz)*1d0-0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(sporof*w(iz)*1d0-sporo(iz)*w(iz)*1d0)/dz(iz)  &
            )
        amx(row,row-nsp) = ( &
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            + adf(iz)*up(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            )
        ! do iiz = 1, nz
            ! if (iiz==iz) cycle
            ! col = 1 + (iiz-1)*nsp
            ! amx(row,col) = amx(row,col) -trans(iz,iiz)
        ! enddo
    else 
        ymx(row) = & 
            + sporo(iz)*(-pt(iz))/dt 
            ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(ptx(iz+1)-ptx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ptx(iz)-ptx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! + w(iz)*(ptx(iz)-ptx(iz-1))/dz(iz)  
        amx(row,row) = (&
            + sporo(iz)*(1d0)/dt &
            ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(-1d0)/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! - trans(iz,iz)  &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz)*w(iz)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-0d0)/dz(iz)  &
            )
        amx(row,row+nsp) =  (&
            ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(1d0)/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - 0d0)/dz(iz)  & 
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz)*w(iz)*0d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(iz)  &
            )
        amx(row,row-nsp) =  (&
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            + adf(iz)*up(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            )
    endif
    ! if (nonlocal(2)) then 
    if (labs(2).or.turbo2(2)) then 
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp
            if (trans(iiz,iz,2)==0d0) cycle
            ! amx(row,col) = amx(row,col) -trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))
            amx(row,col) = amx(row,col) -trans(iiz,iz,2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))
            ! amx(row,col) = amx(row,col) -trans(iiz,iz,2)/dz(iz)*dz(iiz)
        enddo
    else 
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp
            if (trans(iiz,iz,2)==0d0) cycle
            amx(row,col) = amx(row,col) -trans(iiz,iz,2)/dz(iz)
        enddo
    endif
enddo

ymx = - ymx

#ifndef nonrec
if (any(isnan(ymx))) then 
    print*,'NAN in ymx:pt'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_ymx_pre_pt.txt',status = 'unknown')
    do iz = 1, nmx
        write (file_tmp,*) ymx(iz)
    enddo
    close(file_tmp)
    stop
endif
#endif 

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls) 

#ifndef nonrec
if (any(isnan(amx))) then
    print*,'NAN in amx:pt'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_amx_pt.txt',status = 'unknown')
    do iz = 1, nmx
        write (file_tmp,*) amx(iz,:)
    enddo
    close(file_tmp)
    stop
endif

if (any(isnan(ymx))) then 
    print*,'NAN in ymx:pt'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_ymx_pt.txt',status = 'unknown')
    do iz = 1, nmx
        write (file_tmp,*) ymx(iz)
    enddo
    close(file_tmp)
    stop
endif
#endif

ptx = ymx

pttflx = 0d0 
ptdif = 0d0 
ptadv = 0d0 
ptres = 0d0
ptrain = 0d0

do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then
        pttflx = pttflx + sporo(iz)*(ptx(iz)-pt(iz))/dt*dz(iz)
        ptrain = ptrain - detflx/msed
        ! ptdif = ptdif &
            ! - ((dbio(iz)+dbio(iz+1))*0.5d0*(ptx(iz+1)-ptx(iz))/(0.5d0*(dz(iz)+dz(iz+1)))  &
            ! - 0d0  &  !  no bioturbation loss at the top boundary 
            ! )/dz(iz)  &
            ! - trans(iz,iz)*ptx(iz) &
            ! *(1d0-poro(iz))*rhosed/msed*dz(iz)
        ptadv = ptadv &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ptx(iz)-0d0)/dz(iz)*dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ptx(iz+1)-sporo(iz)*w(iz)*ptx(iz))/dz(iz)*dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ptx(iz+1)-0d0)/dz(iz)*dz(iz)  
        ! do iiz = 1, nz
            ! if (iiz==iz) cycle
            ! ptdif = ptdif -trans(iz,iiz)*ptx(iiz)&
                ! *(1d0-poro(iz))*rhosed/msed*dz(iz)
        ! enddo
    else if (iz == nz) then 
        pttflx = pttflx + (1d0-poro(iz))*(ptx(iz)-pt(iz))/dt*dz(iz)
        ! ptdif = ptdif &
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ptx(iz)-ptx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! - trans(iz,iz)*ptx(iz)  &
            ! *(1d0-poro(iz))*rhosed/msed*dz(iz)
        ptadv = ptadv &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ptx(iz)-sporo(iz-1)*w(iz-1)*ptx(iz-1))/dz(iz)*dz(iz)  &
            + adf(iz)*cnr(iz)*(sporof*w(iz)*ptx(iz)-sporo(iz-1)*w(iz-1)*ptx(iz-1))/dz(iz)*dz(iz)  &
            + adf(iz)*dwn(iz)*(sporof*w(iz)*ptx(iz)-sporo(iz)*w(iz)*ptx(iz))/dz(iz)*dz(iz)  
        ! do iiz = 1, nz
            ! if (iiz==iz) cycle
            ! ptdif = ptdif -trans(iz,iiz)*ptx(iiz)&
                ! *(1d0-poro(iz))*rhosed/msed*dz(iz)
        ! enddo
    else 
        pttflx = pttflx + (1d0-poro(iz))*(ptx(iz)-pt(iz))/dt*dz(iz)
        ! ptdif = ptdif &
            ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(ptx(iz+1)-ptx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ptx(iz)-ptx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! - trans(iz,iz)*ptx(iz)  &
            ! *(1d0-poro(iz))*rhosed/msed*dz(iz)
        ptadv = ptadv &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ptx(iz)-sporo(iz-1)*w(Iz-1)*ptx(iz-1))/dz(iz)*dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ptx(iz+1)-sporo(iz)*w(Iz)*ptx(iz))/dz(iz)*dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ptx(iz+1)-sporo(iz-1)*w(Iz-1)*ptx(iz-1))/dz(iz)*dz(iz)
    endif
    ! if(nonlocal(2)) then 
    if(turbo2(2).or.labs(2)) then 
        do iiz = 1, nz
            ! if (iiz==iz) cycle
            if (trans(iiz,iz,2)==0d0) cycle
            ! ptdif = ptdif -trans(iiz,iz)*ptx(iiz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))*dz(iz)
            ptdif = ptdif -trans(iiz,iz,2)*ptx(iiz)/dz(iz)*dz(iiz)*dz(iz)
        enddo
    else 
        do iiz = 1, nz
            ! if (iiz==iz) cycle
            if (trans(iiz,iz,2)==0d0) cycle
            ptdif = ptdif -trans(iiz,iz,2)*ptx(iiz)/dz(iz)    &
                *dz(iz)
        enddo
    endif
    ! if(nonlocal(2)) then 
    if(turbo2(2).or.labs(2)) then 
        do iiz = 1, nz
            if (trans(iiz,iz,2)==0d0) cycle
            ! dw(iz) = dw(iz) - mvsed*(-trans(iiz,iz)*ptx(iiz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz)))
            dw(iz) = dw(iz) - mvsed*(-trans(iiz,iz,2)*ptx(iiz)/dz(iz)*dz(iiz)*(1d0-poro(iiz)))
            ! dw(iz) = dw(iz) - mvsed*(-trans(iiz,iz,2)*ptx(iiz)/dz(iz)*dz(iiz))
        enddo
    ! endif
    else 
        if (nonlocal(2)) then 
            do iiz = 1, nz
                if (trans(iiz,iz,2)==0d0) cycle
                dw(iz) = dw(iz) - mvsed*(-trans(iiz,iz,2)*ptx(iiz)/dz(iz))
            enddo
        endif 
    endif 
enddo

ptres = pttflx + ptdif + ptadv + ptrain

#ifndef nonrec
if (it==1) write(file_ptflx,*) 'time, pttflx, ptdif, ptadv, ptrain, ptres'
write(file_ptflx,*) time, pttflx, ptdif, ptadv, ptrain, ptres
#endif
! stop
!! ~~~~~~~~~End of detrital particle calculation 


!!!!! PRINTING RESULTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! write(dumchr(1),'(i7.7)') it

! if (it/=1) then 
    ! if (all(omx==0d0) .or. all(ptx==0d0) .or.all(ccx==0d0)) then 
        ! print*,it,'error:all zeros',all(omx==0d0),all(ptx==0d0) ,all(ccx==0d0)
        ! stop
    ! endif
! endif
err_fx = maxval(abs(frt - 1d0))
do iz=1,nz 
    rho(iz) = omx(iz)*mom + ptx(iz)*msed +  sum(ccx(iz,:))*mcc
    frt(iz) = omx(Iz)*mvom + ptx(iz)*mvsed + sum(ccx(iz,:))*mvcc
enddo 
err_f = maxval(abs(frt - 1d0))
if (err_f < err_fx) err_f_min = err_f
!! calculation of burial velocity =============================

wx = w

wi = (detflx/msed*mvsed + sum(ccflx)*mvcc +omflx*mvom + sum(dw*dz))/(1d0-poroi)
wi = (detflx/msed*mvsed + sum(ccflx)*mvcc +omflx*mvom             )/(1d0-poroi)

do iz=1,nz
    if (iz==1) then 
        w(iz)=((1d0-poroi)*wi + dw(iz)*dz(iz))/(1d0-poro(iz))
    else 
        w(iz)=((1d0-poro(iz-1))*w(iz-1) + dw(iz)*dz(iz))/(1d0-poro(iz))
    endif
    
    ! if (w(iz) < 0d0) w(iz) = wi*1d-3
    
enddo

! wxx = w

! do iz=nz-1,1,-1 
    ! if (w(iz)<0d0) then 
        ! w(iz+1) = w(iz)
    ! endif 
! enddo


error = 1d4
! itr = 0

! do while (error>tol) 

! wxx = w

! up=1
! dwn=0
! cnr=0
! adf=1

! nsp = 1
! nmx = nz*nsp
! deallocate(amx,ymx,emx,ipiv)
! allocate(amx(nmx,nmx),ymx(nmx),emx(nmx),ipiv(nmx))

! amx = 0d0
! ymx = 0d0

! do iz = 1,nz 
    ! row = 1 + (iz-1)*nsp 
    ! if (iz == 1) then 
        ! ymx(row) = &
            ! - (detflx/msed*mvsed + sum(ccflx)*mvcc +omflx*mvom)/dz(iz) 
        ! amx(row,row) = (&
            ! + adf(iz)*up(iz)*(sporo(iz)*1d0-0d0)/dz(iz)   &
            ! + adf(iz)*dwn(iz)*(sporo(iz+1)*0d0-sporo(iz)*1d0)/dz(iz)   &
            ! + adf(iz)*cnr(iz)*(sporo(iz+1)*0d0-0d0)/dz(iz)   &
            ! )            
        ! amx(row,row+nsp) =  (&
            ! + adf(iz)*dwn(iz)*(sporo(iz+1)*1d0-sporo(iz)*0d0)/dz(iz)   &
            ! + adf(iz)*cnr(iz)*(sporo(iz+1)*1d0-0d0)/dz(iz)   &
            ! )
    ! else if (iz == nz) then 
        ! amx(row,row) = (&
            ! + adf(iz)*up(iz)*(sporo(iz)*1d0-0d0)/dz(iz)  &
            ! + adf(iz)*cnr(iz)*(sporof*1d0-0d0)/dz(iz)  &
            ! + adf(iz)*dwn(iz)*(sporof*1d0-sporo(iz)*1d0)/dz(iz)  &
            ! )
        ! amx(row,row-nsp) = ( & 
            ! + adf(iz)*up(iz)*(0d0-sporo(iz-1)*1d0)/dz(iz)  &
            ! + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*1d0)/dz(iz)  &
            ! )
    ! else 
        ! amx(row,row) = (&
            ! + adf(iz)*up(iz)*(sporo(iz)*1d0-0d0)/dz(iz)  &
            ! + adf(iz)*dwn(iz)*(sporo(iz+1)*0d0-sporo(iz)*1d0)/dz(iz)  &
            ! + adf(iz)*cnr(iz)*(sporo(iz+1)*0d0-0d0)/dz(iz)  &
            ! )
        ! amx(row,row+nsp) =  (&
            ! + adf(iz)*dwn(iz)*(sporo(iz+1)*1d0-sporo(iz)*0d0)/dz(iz)  &
            ! + adf(iz)*cnr(iz)*(sporo(iz+1)*1d0-0d0)/dz(iz)  &
            ! )
        ! amx(row,row-nsp) =  (&
            ! + adf(iz)*up(iz)*(0d0-sporo(iz-1)*1d0)/dz(iz)  &
            ! + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*1d0)/dz(iz)  &
            ! )
    ! endif
    ! ymx(row) = ymx(row) - dw(iz)
! enddo

! ymx = -ymx

! #ifndef nonrec
! if (any(isnan(ymx))) then 
    ! print*,'NAN in ymx:w'
    ! open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_ymx_pre_w.txt',status = 'unknown')
    ! do iz = 1, nmx
        ! write (file_tmp,*) ymx(iz)
    ! enddo
    ! close(file_tmp)
    ! stop
! endif
! #endif 

! call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,info) 

! #ifndef nonrec
! if (any(isnan(amx))) then
    ! print*,'NAN in amx:w'
    ! open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_amx_w.txt',status = 'unknown')
    ! do iz = 1, nmx
        ! write (file_tmp,*) amx(iz,:)
    ! enddo
    ! close(file_tmp)
    ! stop
! endif

! if (any(isnan(ymx))) then 
    ! print*,'NAN in ymx:w'
    ! open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_ymx_w.txt',status = 'unknown')
    ! do iz = 1, nmx
        ! write (file_tmp,*) ymx(iz)
    ! enddo
    ! close(file_tmp)
    ! stop
! endif
! #endif

! w = ymx


! if (time>time_spn .and. time <= time_spn + time_trs) w = w*frt

!!!! IF modifying porosity each time !!!!!!!!!!!!!!
! #ifdef sense 
! calgg = sum(ccx(nz,:))*mcc/rho(nz)
! pore_max =  1d0 - ( 0.483d0 + 0.45d0 * calgg) / 2.5d0
! exp_pore = 0.25d0*calgg + 3.d0 *(1d0-calgg)
! poro = EXP(-z/exp_pore) * (1.d0-pore_max) + pore_max 

! ff = poro*poro  ! Archer 

! dif_dic = dif_dic0*ff
! dif_alk = dif_alk0*ff
! dif_o2 = dif_o20*ff

! sporo = 1d0 - poro
! #endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! if (any(w<0d0)) then
    ! print*,'negative burial velocity'
    ! print*, w
    ! stop
! endif 

! ------------ determine calculation scheme for advection 
up = 0
dwn=0
cnr =0
adf=1d0
do iz=1,nz 
    if (iz==1) then 
        if (w(iz)>=0d0 .and. w(iz+1)>=0d0) then
            up(iz) = 1
        elseif (w(iz)<=0d0 .and. w(iz+1)<=0d0) then
            dwn(iz) = 1
        else 
            if (.not.(w(iz)*w(iz+1) <=0d0)) then 
                print*,'error'
                stop
            endif
            cnr(iz) = 1
        endif
    elseif (iz==nz) then 
        if (w(iz)>=0d0 .and. w(iz-1)>=0d0) then
            up(iz) = 1
        elseif (w(iz)<=0d0 .and. w(iz-1)<=0d0) then
            dwn(iz) = 1
        else 
            if (.not.(w(iz)*w(iz-1) <=0d0)) then 
                print*,'error'
                stop
            endif
            cnr(iz) = 1
        endif
    else 
        if (w(iz) >=0d0) then 
            if (w(iz+1)>=0d0 .and. w(iz-1)>=0d0) then
                up(iz) = 1
            else
                cnr(iz) = 1
            endif
        else  
            if (w(iz+1)<=0d0 .and. w(iz-1)<=0d0) then
                dwn(iz) = 1
            else
                cnr(iz) = 1
            endif
        endif
    endif
enddo        

if (sum(up(:)+dwn(:)+cnr(:))/=nz) then
    print*,error,sum(up),sum(dwn),sum(cnr)
    stop
endif

do iz=1,nz-1
    if (cnr(iz)==1 .and. cnr(iz+1)==1) then 
        if (w(iz)>=0d0 .and. w(iz+1) < 0d0) then
            ! cnr(iz+1)=0
            ! dwn(iz+1)=1
            ! cnr(iz+1)=abs(frt(iz)-1d0)/(abs(frt(iz+1)-1d0)+abs(frt(iz)-1d0))
            ! cnr(iz)=abs(frt(iz+1)-1d0)/(abs(frt(iz+1)-1d0)+abs(frt(iz)-1d0))
            
            ! if (itr_f == 0) then 
                ! corrf = 2d0
                ! corrf = abs(w(iz)**corrf)/(abs(w(iz+1)**corrf)+abs(w(iz)**corrf))
            ! elseif (itr_f == 1) then 
                ! if (corrf>0.5d0) then
                    ! df = -0.5d0*corrf
                ! else 
                    ! df = 0.5d0*corrf
                ! endif 
                ! if (corrf+df >=1d0) df = 1d0-corrf
                ! if (corrf+df <=0d0) df = -corrf
                ! corrf = corrf + df
            ! elseif (itr_f == 2) then 
                ! dfrt_df = (err_f -err_fx)/df
                ! if (dfrt_df <0d0) then 
                    ! df = 0.5d0*corrf
                ! else 
                    ! df = -0.5d0*corrf
                ! endif 
                ! if (corrf+df >=1d0) df = 1d0-corrf
                ! if (corrf+df <=0d0) df = -corrf
                ! corrf = corrf + df
            ! elseif (itr_f >=3) then 
                ! dfrt_dfx = dfrt_df
                ! dfrt_df = (err_f -err_fx)/df
                ! d2frt_df2 = (dfrt_df - dfrt_dfx)/df
                ! if (dfrt_df <0d0) then 
                    ! df = 0.5d0*corrf
                ! else 
                    ! df = -0.5d0*corrf
                ! endif 
                ! if (corrf+df >=1d0) df = 1d0-corrf
                ! if (corrf+df <=0d0) df = -corrf
                ! corrf = corrf + df
            ! endif 
            ! itr_f = itr_f + 1
            ! cnr(iz+1) = corrf
            ! cnr(iz) = 1d0-corrf
            
            corrf = 2d0
            cnr(iz+1)=abs(w(iz)**corrf)/(abs(w(iz+1)**corrf)+abs(w(iz)**corrf))
            cnr(iz)=abs(w(iz+1)**corrf)/(abs(w(iz+1)**corrf)+abs(w(iz)**corrf))
            ! cnr(iz+1)=abs(w(iz))/(abs(w(iz+1))+abs(w(iz)))
            ! cnr(iz)=abs(w(iz+1))/(abs(w(iz+1))+abs(w(iz)))
            dwn(iz+1)=1d0-cnr(iz+1)
            up(iz)=1d0-cnr(iz)
        endif 
    endif 
    if (cnr(iz)==1 .and. cnr(iz+1)==1) then 
        if (w(iz)< 0d0 .and. w(iz+1) >= 0d0) then
            cnr(iz+1)=0
            cnr(iz)=0
            up(iz+1)=1
            dwn(iz)=1
            ! adf(iz)=0.5d0
            ! adf(iz+1)=0.5d0
            adf(iz)=abs(w(iz+1))/(abs(w(iz+1))+abs(w(iz)))
            adf(iz+1)=abs(w(iz))/(abs(w(iz+1))+abs(w(iz)))
        endif 
    endif 
enddo       
! up = 0
! dwn=0
! cnr =1
! adf=1d0

! error = maxval(abs((w-wxx)/wxx))
! enddo
itr_w = itr_w + 1
err_w = maxval(abs((w-wx)/wx))
! if (maxval(abs(frt - 1d0))<err_f_min) then 
if (err_w<err_w_min) then 
    ! err_f_min= maxval(abs(frt - 1d0))
    err_w_min= err_w
    wxx = wx  ! recording w which minimizes deviation of total sld fraction from 1 
endif 
if (itr_w>100) then 
    if (itr_w==101) then 
        w = wxx
        go to 300
    elseif (itr_w==102) then 
        w = wxx
        write(file_err,*) 'not converging w',time, err_w, err_w_min
        go to 400
    endif 
endif
if (err_w > tol) go to 300
! if (err_w > 1d-2) go to 300

400 continue

!! depth -age 
dage=0d0
age = 0d0
dage = dz/w  
do iz=1,nz
    if (iz==1) age(iz)=dage(iz)*0.5d0
    if (iz/=1) age(iz) = age(iz-1)+dage(iz-1)*0.5d0 + 0.5d0*dage(iz)
enddo

! ---------------------

if (any(rho<0d0)) then
    print*,'negative density'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'NEGATIVE_RHO.txt',status = 'unknown')
    do iz = 1, nz
        write (file_tmp,*) z(iz),rho(iz),w(iz),up(iz),dwn(iz),cnr(iz),adf(iz)
    enddo
    close(file_tmp)
    stop
endif 

!/////// ISOTOPES /////
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

#ifndef nonrec
if (time>=rectime(cntrec)) then 
    write(dumchr(1),'(i3.3)') cntrec 
    
    open(unit=file_tmp,file=trim(adjustl(workdir))//'ptx-'//trim(adjustl(dumchr(1)))//'.txt' &
        ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),ptx(iz)*msed/rho(iz)*100d0,rho(iz),frt(iz)  ,w(iz)
    enddo
    close(file_tmp)

    open(unit=file_tmp,file=trim(adjustl(workdir))//'ccx-'//trim(adjustl(dumchr(1)))//'.txt' &
        ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),sum(ccx(iz,:))*mcc/rho(iz)*100d0, dicx(iz)*1d3, alkx(iz)*1d3  &
            , co3x(iz)*1d3-co3sat, sum(rcc(iz,:)),-log10(prox(iz)) 
    enddo
    close(file_tmp)

    open(unit=file_tmp,file=trim(adjustl(workdir))//'omx-'//trim(adjustl(dumchr(1)))//'.txt'  &
        ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),omx(iz)*mom/rho(iz)*100d0
    enddo
    close(file_tmp)

    open(unit=file_tmp,file=trim(adjustl(workdir))//'o2x-'//trim(adjustl(dumchr(1)))//'.txt'  &
        ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),o2x(iz)*1d3, oxco2(iz), anco2(iz)
    enddo
    close(file_tmp)
    
    open(unit=file_tmp,file=trim(adjustl(workdir))//'ccx_sp-'  &
        //trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),(ccx(iz,isp)*mcc/rho(iz)*100d0,isp=1,nspcc) 
    enddo
    close(file_tmp)
    
    open(unit=file_tmp,file=trim(adjustl(workdir))//'sig-'  &
        //trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),d13c_blk(iz),d18o_blk(iz)
    enddo
    close(file_tmp)
    
    open(unit=file_tmp,file=trim(adjustl(workdir))//'bur-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),w(iz),up(iz),dwn(iz),cnr(iz),adf(iz)
    enddo
    close(file_tmp)
    
    cntrec = cntrec + 1
    if (cntrec == nrec+1) exit
endif 
#endif

#ifndef nondisp    
print*, 'time   :',time, maxval(abs(frt - 1d0))
print*,'~~~~ conc ~~~~'
print'(A,5E11.3)', 'z  :',(z(iz),iz=1,nz,nz/5)
print'(A,5E11.3)', 'om :',(omx(iz)*mom/rho(iz)*100d0,iz=1,nz,nz/5)
print'(A,5E11.3)', 'o2 :',(o2x(iz)*1d3,iz=1,nz,nz/5)
print'(A,5E11.3)', 'cc :',(sum(ccx(iz,:))*mcc/rho(iz)*100d0,iz=1,nz,nz/5)
print'(A,5E11.3)', 'dic:',(dicx(iz)*1d3,iz=1,nz,nz/5)
print'(A,5E11.3)', 'alk:',(alkx(iz)*1d3,iz=1,nz,nz/5)
print'(A,5E11.3)', 'sed:',(ptx(iz)*msed/rho(iz)*100d0,iz=1,nz,nz/5)
print*, '   ..... multiple cc species ..... '
do isp=1,nspcc 
    print'(i0.3,":",5E11.3)',isp,(ccx(iz,isp)*mcc/rho(iz)*100d0,iz=1,nz,nz/5)
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

print*,'==== burial etc ===='
print'(A,5E11.3)', 'z  :',(z(iz),iz=1,nz,nz/5)
print'(A,5E11.3)', 'w  :',(w(iz),iz=1,nz,nz/5)
print'(A,5E11.3)', 'rho:',(rho(iz),iz=1,nz,nz/5)
print'(A,5E11.3)', 'frc:',(frt(iz),iz=1,nz,nz/5)

print*,''
#endif

!! in theory, o2dec/ox2om + alkdec = dicdec = omdec (in absolute value)
! if (dicdec /= 0d0) then 
if (om2cc /= 0d0) then 
    if ( abs((o2dec/ox2om - alkdec + dicdec)/dicdec) > tol) then 
        print*, abs((o2dec/ox2om + alkdec - dicdec)/dicdec) ,o2dec/ox2om,alkdec,dicdec
        write(file_err,*) trim(adjustl(dumchr(1))), time, dt &
            , abs((o2dec/ox2om + alkdec - dicdec)/dicdec),o2dec/ox2om,alkdec,dicdec
        ! stop
    endif
endif 

#ifndef nonrec 
write(file_totfrac,*) time,maxval(abs(frt - 1d0))
#ifndef size 
! if (.not.any(w<0d0)) then  ! not recording when burial is negative 
if (all(w>=0d0)) then  ! not recording when burial is negative 
    write(file_sigmly,*) time-age(izrec),d13c_blk(izrec),d18o_blk(izrec) &
        ,sum(ccx(izrec,:))*mcc/rho(izrec)*100d0,ptx(izrec)*msed/rho(izrec)*100d0
    write(file_sigmlyd,*) time-age(izrec2),d13c_blk(izrec2),d18o_blk(izrec2) &
        ,sum(ccx(izrec2,:))*mcc/rho(izrec2)*100d0,ptx(izrec2)*msed/rho(izrec2)*100d0
    write(file_sigbtm,*) time-age(nz),d13c_blk(nz),d18o_blk(nz) &
        ,sum(ccx(nz,:))*mcc/rho(nz)*100d0,ptx(nz)*msed/rho(nz)*100d0
endif 
#else 
if (all(w>=0d0)) then  ! not recording when burial is negative 
    write(file_sigmly,*) time-age(izrec),d13c_blk(izrec),d18o_blk(izrec) &
        ,sum(ccx(izrec,:))*mcc/rho(izrec)*100d0,ptx(izrec)*msed/rho(izrec)*100d0  &
        ,d13c_blkf(izrec),d18o_blkf(izrec),sum(ccx(izrec,1:4))*mcc/rho(izrec)*100d0  &
        ,d13c_blkc(izrec),d18o_blkc(izrec),sum(ccx(izrec,5:8))*mcc/rho(izrec)*100d0  
    write(file_sigmlyd,*) time-age(izrec2),d13c_blk(izrec2),d18o_blk(izrec2) &
        ,sum(ccx(izrec2,:))*mcc/rho(izrec2)*100d0,ptx(izrec2)*msed/rho(izrec2)*100d0  &
        ,d13c_blkf(izrec2),d18o_blkf(izrec2),sum(ccx(izrec2,1:4))*mcc/rho(izrec2)*100d0  &
        ,d13c_blkc(izrec2),d18o_blkc(izrec2),sum(ccx(izrec2,5:8))*mcc/rho(izrec2)*100d0  
    write(file_sigbtm,*) time-age(nz),d13c_blk(nz),d18o_blk(nz) &
        ,sum(ccx(nz,:))*mcc/rho(nz)*100d0,ptx(nz)*msed/rho(nz)*100d0 &
        ,d13c_blkf(nz),d18o_blkf(nz),sum(ccx(nz,1:4))*mcc/rho(nz)*100d0  &
        ,d13c_blkc(nz),d18o_blkc(nz),sum(ccx(nz,5:8))*mcc/rho(nz)*100d0  
endif 
#endif
#endif 

! pause

time = time + dt
it = it + 1

o2 = o2x
om = omx

cc = ccx
dic = dicx
alk = alkx

pt = ptx

enddo

#ifndef nonrec
open(unit=file_tmp,file=trim(adjustl(workdir))//'sp-trace.txt',action='write',status='replace') 
do isp  = 1,nspcc
    write(file_tmp,*) isp,d13c_sp(isp),d18o_sp(isp)
enddo
close(file_tmp)
#endif 

#ifndef nonrec
close(file_ptflx)
close(file_ccflx)
close(file_omflx)
close(file_o2flx)
close(file_dicflx)
close(file_alkflx)
close(file_err)
close(file_bound)
close(file_totfrac)
close(file_sigmly)
close(file_sigmlyd)
close(file_sigbtm)
do isp=1,nspcc 
    close(file_ccflxes(isp))
enddo
#endif

workdir = 'C:/Users/YK/Desktop/Sed_res/'
workdir = trim(adjustl(workdir))//'test-translabs/res/'
workdir = trim(adjustl(workdir))//'multi/'
#ifdef test 
workdir = trim(adjustl(workdir))//'test/'
#endif
if (.not. anoxic) then 
    workdir = trim(adjustl(workdir))//'ox'
else 
    workdir = trim(adjustl(workdir))//'oxanox'
endif

if (any(labs)) workdir = trim(adjustl(workdir))//'-labs'
if (any(turbo2)) workdir = trim(adjustl(workdir))//'-turbo2'
if (any(nobio)) workdir = trim(adjustl(workdir))//'-nobio'

workdir = trim(adjustl(workdir))//'/'

call system ('mkdir -p '//trim(adjustl(workdir)))

open(unit=file_tmp,file=trim(adjustl(workdir))//'lys_sense_'// &
    'cc-'//trim(adjustl(chr(1,4)))//'_rr-'//trim(adjustl(chr(2,4)))  &
    //'.txt',action='write',status='unknown',access='append') 
write(file_tmp,*) 1d6*(co3i*1d3-co3sat), sum(ccx(1,:))*mcc/rho(1)*100d0, frt(1)  &
    ,sum(ccx(nz,:))*mcc/rho(nz)*100d0, frt(nz),sum(ccx(izml,:))*mcc/rho(izml)*100d0, frt(izml)
close(file_tmp)

open(unit=file_tmp,file=trim(adjustl(workdir))//'ccbur_sense_'// &
    'cc-'//trim(adjustl(chr(1,4)))//'_rr-'//trim(adjustl(chr(2,4)))  &
    //'.txt',action='write',status='unknown',access='append') 
write(file_tmp,*) 1d6*(co3i*1d3-co3sat), 1d6*sum(ccadv)
close(file_tmp)

end program 
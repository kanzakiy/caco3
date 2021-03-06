
!**************************************************************************************************************************************
subroutine caco3(ccflxi,om2cc,depi2,dti2,filechr,oxonly,biotmode) 
! a signal tracking diagenesis
implicit none 

#include <defines.h>
integer(kind=4)nz      ! grid number
integer(kind=4)nspcc   ! number of CaCO3 species 

#ifndef nzinput
#define nzinput 100
#endif 

#ifndef nspccinput
#ifdef sense
#define nspccinput 1
#elif defined track2
#define nspccinput 42
#elif defined size
#define nspccinput 8
#elif defined isotrack
#define nspccinput 5
#else
#define nspccinput 4
#endif
#endif 

parameter(nz=nzinput)  
parameter(nspcc=nspccinput)

real(kind=8),intent(in)::ccflxi,om2cc,depi2,dti2
character*255,intent(in)::filechr,biotmode
logical,intent(in)::oxonly

integer(kind=4),parameter :: nsig = 2
integer(kind=4) nt_spn,nt_trs,nt_aft
real(kind=8) cc(nz,nspcc),ccx(nz,nspcc)  ! mol cm-3 sld; concentration of caco3, subscript x denotes dummy variable used during iteration 
real(kind=8) om(nz),omx(nz)  ! mol cm-3 sld; om conc. 
real(kind=8) pt(nz), ptx(nz) ! mol cm-3 sld; clay conc.
real(kind=8) ccflx(nspcc), d13c_sp(nspcc),d18o_sp(nspcc) ! flux of caco3, d13c signal of caco3, d18o signal of caco3
real(kind=8) d13c_blk(nz), d18o_blk(nz)  ! d13c signal of bulk caco3, d18o signal of bulk caco3 
real(kind=8) d13c_blkf(nz), d18o_blkf(nz),  d13c_blkc(nz), d18o_blkc(nz) ! subscripts f and c denotes variables of fine and coarse caco3 species, respectively 
real(kind=8) d13c_flx, d18o_flx  ! d13c signal averaged over flux values, d18o counterpart 
real(kind=8) d13c_ocni, d13c_ocnf, d13c_ocn  ! initial value of ocean d13c, final value of ocean d13c, ocean d13c  
real(kind=8) d18o_ocni, d18o_ocnf, d18o_ocn  ! the same as above expect 18o insted of 13c
!!!!!!!!!!!! used only when isotrack is on !!!!!!!!!!!!!!!!!!!!!!!!
real(kind=8) capd47_ocni,capd47_ocnf
real(kind=8) capd47_ocn
real(kind=8) r13c_blk(nz),r18o_blk(nz),r17o_blk(nz),d17o_blk(nz),d14c_age(nz),capd47(nz)
real(kind=8) r45,r46,r47,r45s,r46s,r47s
real(kind=8) :: r18o_pdb = 0.0020672d0 ! Fry (2006)
real(kind=8) :: r17o_pdb = 0.0003859d0 ! Fry (2006) cf., 0.000379 by Hoef (2015) saying after Hayes (1983)
real(kind=8) :: r18o_smow = 0.0020052d0 ! Fry (2006)
real(kind=8) :: r17o_smow = 0.0003799d0 ! Fry (2006), cf., 0.000373 by Hoef (2015) saying after Hayes (1983)
real(kind=8) :: r13c_pdb = 0.011180d0 ! Fry (2006)
real(kind=8) :: d13c_om = -25d0 ! e.g., Ridgwell and Arndt (2014) (probably vs PDB)
real(kind=8) :: d18o_o2 = 23.5d0 ! Kroopnick and Craig (1972) vs SMOW; 18O16O/16O16O
real(kind=8) :: d18o_so4 = 9.5d0 ! Longinelli and Craig (1967) vs SMOW
real(kind=8) :: f13c18o, f12c18o,f13c16o, f12c16o, f13c17o, f12c17o  !  relative to whole species 
real(kind=8) :: f12c17o18o,f12c17o17o,f12c18o18o
integer(kind=4) :: i12c16o=1,i12c18o=2,i13c16o=3,i13c18o=4,i14c=5
integer(kind=4) :: i12c17o=6,i12c17o18o=7,i12c18o18o=8,i13c17o=9,i12c17o17o=10 ! not implemented yet
real(kind=8) krad(nz,nspcc)  ! caco3 decay consts (for 14c alone)
real(kind=8) deccc(nz,nspcc)  ! radio-active decay rate of caco3 
real(kind=8) ddeccc_dcc(nz,nspcc)  ! radio-active decay rate of caco3 
real(kind=8) ccrad(nspcc),alkrad
real(kind=8) :: k14ci = 1d0/8033d0 ! [yr-1], Aloisi et al. 2004
real(kind=8) :: r14ci = 1.2d-12 ! c14/c12 in modern, Aloisi et al. 2004, citing Kutschera 2000
real(kind=8) r2d,d2r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
real(kind=8) flxfrc(nspcc),flxfrc2(nspcc)  ! flux fractions used for assigning flux values to realized isotope input changes 
! real(kind=8) :: ccflxi = 10d-6 ! mol (CaCO3) cm-2 yr-1  ! a reference caco3 flux; Emerson and Archer (1990) 
real(kind=8) :: omflx = 12d-6 ! mol cm-2 yr-1       ! a reference om flux; Emerson and Archer (1990)
real(kind=8) :: detflx = 180d-6 ! g cm-2 yr-1  ! a reference detrital flux; MUDS input http://forecast.uchicago.edu/Projects/muds.html
#ifndef mocsy
real(kind=8) :: alki =  2285d0 ! uM  ! a reference ALK; MUDS
real(kind=8) :: dici = 2211d0 ! uM   ! a reference DIC; MUDS 
#else 
! real(kind=8) :: alki =  2295d0 ! uM  ! a reference ALK from mocsy
! real(kind=8) :: dici = 2154d0 ! uM   ! a reference DIC from mocsy 
real(kind=8) :: alki =  2285d0 ! uM  ! a reference ALK; MUDS
real(kind=8) :: dici = 2211d0 ! uM   ! a reference DIC; MUDS 
#endif 
real(kind=8) :: o2i = 165d0 ! uM     ! a reference O2 ; MUDS
! real(kind=8) :: komi = 2d0  ! /yr  ! a reference om degradation rate const.; MUDS 
! real(kind=8) :: komi = 0.5d0  ! /yr  ! arbitrary 
! real(kind=8) :: komi = 0.1d0  ! /yr  ! Canfield 1994
real(kind=8) :: komi = 0.06d0  ! /yr  ! ?? Emerson 1985? who adopted relatively slow decomposition rate 
! real(kind=8) :: kcci = 10.0d0*365.25d0  ! /yr; a reference caco3 dissolution rate const. 
#ifndef nodissolve
real(kind=8) :: kcci = 1d0*365.25d0  ! /yr  ;cf., 0.15 to 30 d-1 Emerson and Archer (1990) 0.1 to 10 d-1 in Archer 1991
#else
real(kind=8) :: kcci = 0d0*365.25d0  ! /yr 
#endif 
real(kind=8) :: poroi = 0.8d0  ! a reference porosity 
real(kind=8) :: keqcc = 4.4d-7   ! mol2 kg-2 caco3 solutiblity (Mucci 1983 cited by Emerson and Archer 1990)
real(kind=8) :: ncc = 4.5d0   ! (Archer et al. 1989) reaction order for caco3 dissolution 
real(kind=8) :: temp = 2d0  ! C a refernce temperature 
real(kind=8) :: sal = 35d0  ! wt o/oo  salinity 
real(kind=8) :: cai = 10.3d-3 ! mol kg-1 calcium conc. seawater 
real(kind=8) :: fact = 1d-3 ! w/r factor to facilitate calculation 
! real(kind=8) :: rhosed = 2.09d0 ! g/cm3 sediment particle density assuming opal
real(kind=8) :: rhosed = 2.6d0 ! g/cm3 sediment particle density assming kaolinite 
real(kind=8) :: rhoom = 1.2d0 ! g/cm3 organic particle density 
real(kind=8) :: rhocc = 2.71d0 ! g/cm3 organic particle density 
real(kind=8) :: mom = 30d0 ! g/mol OM assuming CH2O
! real(kind=8) :: msed = 87.11d0 ! g/mol arbitrary sediment g/mol assuming opal (SiO2•n(H2O) )
real(kind=8) :: msed = 258.16d0 ! g/mol arbitrary sediment g/mol assuming kaolinite ( 	Al2Si2O5(OH)4 )
real(kind=8) :: mcc(nspcc) = 100d0 ! g/mol CaCO3 
real(kind=8) :: ox2om = 1.3d0 ! o2/om ratio for om decomposition (Emerson 1985; Archer 1991)
! real(kind=8) :: om2cc = 0.666d0  ! rain ratio of organic matter to calcite
! real(kind=8) :: om2cc = 0.5d0  ! rain ratio of organic matter to calcite
real(kind=8) :: o2th = 0d0 ! threshold oxygen level below which not to calculate 
real(kind=8) :: dev = 1d-6 ! deviation addumed 
real(kind=8) :: zml_ref = 12d0 ! a referece mixed layer depth
real(kind=8) :: ccx_th = 1d-300 ! threshold caco3 conc. (mol cm-3) below which calculation is not conducted 
real(kind=8) :: omx_th = 1d-300 ! threshold om    conc. (mol cm-3) below which calculation is not conducted 
real(kind=8) dif_alk0, dif_dic0, dif_o20  ! diffusion coefficient of alk, dik and o2 in seawater 
real(kind=8) zml(nspcc+2) , zrec, zrec2  ! mixed layer depth, sediment depth where proxy signal is read 1 & 2
real(kind=8) chgf  ! variable to check change in total fraction of solid materials
real(kind=8) flxfin, flxfini, flxfinf  !  flux ratio of fine particles: i and f denote initial and final values  
real(kind=8) pore_max, exp_pore, calgg  ! parameters to determine porosity in Archer (1991) 
real(kind=8) mvom, mvsed, mvcc(nspcc)  ! molar volumes (cm3 mol-1) mv_i = m_i/rho_i where i = om, sed and cc for organic matter, clay and caco3, respectively
real(kind=8) keq1, keq2  ! equilibrium const. for h2co3 and hco3 dissociations and functions to calculate them  
real(kind=8) co3sat, keqag  ! co3 conc. at caco3 saturation and solubility product of aragonite  
real(kind=8) zox, zoxx  ! oxygen penetration depth (cm) and its dummy variable 
real(kind=8) dep,depi,depf ! water depth, i and f denote initial and final values 
real(kind=8) corrf, df, err_f, err_fx  !  variables to help total fraction of solid materials to converge 1
real(kind=8) dic(nz), dicx(nz), alk(nz), alkx(nz), o2(nz), o2x(nz) !  mol cm-3 porewater; dic, alk and o2 concs., x denotes dummy variables 
real(kind=8) dif_dic(nz), dif_alk(nz), dif_o2(nz) ! dic, alk and o2 diffusion coeffs inclueing effect of tortuosity
real(kind=8) dbio(nz), ff(nz)  ! biodiffusion coeffs, and formation factor 
real(kind=8) co2(nz), hco3(nz), co3(nz), pro(nz) ! co2, hco3, co3 and h+ concs. 
real(kind=8) co2x(nz), hco3x(nz), co3x(nz), prox(nz),co3i  ! dummy variables and initial co3 conc.  
real(kind=8) ohmega(nz),dohmega_ddic(nz),dohmega_dalk(nz)
real(kind=8) poro(nz), rho(nz), frt(nz)  ! porositiy, bulk density and total volume fraction of solid materials 
real(kind=8) sporo(nz), sporoi, porof, sporof  ! solid volume fraction (1 - poro), i and f denote the top and bottom values of variables  
real(kind=8) rcc(nz,nspcc)  ! dissolution rate of caco3 
real(kind=8) drcc_dcc(nz,nspcc), drcc_ddic(nz,nspcc)! derivatives of caco3 dissolution rate wrt caco3 and dic concs.
real(kind=8) drcc_dalk(nz,nspcc), drcc_dco3(nz,nspcc),drcc_dohmega(nz,nspcc) ! derivatives of caco3 dissolution rate wrt alk and co3 concs.
real(kind=8) dco3_ddic(nz), dco3_dalk(nz)  ! derivatives of co3 conc. wrt dic and alk concs. 
real(kind=8) ddum(nz)  ! dummy variable 
real(kind=8) dpro_dalk(nz), dpro_ddic(nz)  ! derivatives of h+ conc. wrt alk and dic concs. 
real(kind=8) kcc(nz,nspcc)  ! caco3 dissolution rate consts. 
real(kind=8) kom(nz)  ! degradation rate consts. 
real(kind=8) oxco2(nz),anco2(nz)  ! oxic and anoxic om degradation rate 
real(kind=8) w(nz) , wi, dw(nz), wx(nz) ! burial rate, burial rate initial guess, burial rate change, burial rate dummy 
real(kind=8) err_w, wxx(nz) ! err in burial rate, dummy dummy burial rate  
real(kind=8) err_f_min, dfrt_df, d2frt_df2, dfrt_dfx, err_w_min  ! variables to minimize errors in burial rate and total fractions of solid phases 
real(kind=8) z(nz), dz(nz)! depth, individual sediment layer thickness
real(kind=8) eta(nz), beta  ! parameters to make a grid 
real(kind=8) dage(nz), age(nz)  ! individual time span and age of sediment grids  
#ifdef sense
real(kind=8) :: ztot = 50d0 ! cm , total sediment thickness 
#else
real(kind=8) :: ztot = 500d0 ! cm 
#endif
integer(kind=4) :: nsp = 3  ! independent chemical variables, this does not have to be decided here    
integer(kind=4) :: nmx      ! row (and col) number of matrix created to solve linear difference equations 
real(kind=8),allocatable :: amx(:,:),ymx(:),emx(:) ! amx and ymx correspond to A and B in Ax = B, but ymx is also x when Ax = B is solved. emx is array of error 
integer(kind=4),allocatable :: ipiv(:),dumx(:,:) ! matrix used to solve linear system Ax = B 
integer(kind=4) infobls, infosbr  ! variables used to tell errors when calling a subroutine to solve matrix 
real(kind=8) error, error2, minerr  !  errors in iterations and minimum error produced 
real(kind=8) :: tol = 1d-6  ! tolerance of error 
integer(kind=4) iz, row, col, itr  , it, iiz, itr_w, itr_f ! integers for sediment grid, matrix row and col and iteration numbers 
integer(kind=4) cntsp  ! counting caco3 species numbers 
integer(kind=4) izrec, izrec2  ! grid number where signal is read 
integer(kind=4) :: nt = 1000000  ! maximum interation for time (do not have to be defined ) 
integer(kind=4),parameter :: nrec = 15  ! total recording time of sediment profiles 
integer(kind=4) cntrec, itrec  ! integers used for counting recording time 
integer(kind=4) :: izox_minerr =0  ! grid number where error in zox is minimum 
real(kind=8) up(nz), dwn(nz), cnr(nz) ! advection calc. schemes; up or down wind, or central schemes if 1
real(kind=8) adf(nz)  ! factor to make sure mass conversion 
real(kind=8) :: time, dt = 1d2  ! time and time step 
real(kind=8) rectime(nrec), time_max ! recording time and maximum time 
real(kind=8) dumreal  !  dummy variable 
real(kind=8) time_spn, time_trs, time_aft  ! time durations of spin-up, signal transition and after transition  
! fluxes, adv, dec, dis, dif, res, t and rain denote burial, decomposition, dissoution, diffusion, residual, time change and rain fluxes, respectively  
real(kind=8) :: omadv, omdec, omdif, omrain, omres, omtflx ! fluxes of om
real(kind=8) :: o2tflx, o2dec, o2dif, o2res  ! o2 fluxes 
real(kind=8) :: cctflx(nspcc), ccdis(nspcc), ccdif(nspcc), ccadv(nspcc), ccres(nspcc),ccrain(nspcc) ! caco3 fluxes 
real(kind=8) :: dictflx, dicdis, dicdif, dicdec, dicres  ! dic fluxes 
real(kind=8) :: alktflx, alkdis, alkdif, alkdec, alkres  ! alk fluxes 
real(kind=8) :: pttflx, ptdif, ptadv, ptres, ptrain  ! clay fluxes 
real(kind=8) :: trans(nz,nz,nspcc+2)  ! transition matrix 
real(kind=8) :: transdbio(nz,nz), translabs(nz,nz) ! transition matrices created assuming Fickian mixing and LABS simulation
real(kind=8) :: transturbo2(nz,nz), translabs_tmp(nz,nz) ! transition matrices assuming random mixing and LABS simulation 
character*512 workdir  ! work directory and created file names 
! character*255 filechr  ! work directory and created file names 
character*25 dumchr(3)  ! character dummy variables 
character*25 arg, chr(3,4)  ! used for reading variables and dummy variables
integer(kind=4) dumint(8)  ! dummy integer 
integer(kind=4) idp, izox  ! integer for depth and grid number of zox 
integer(kind=4) narg, ia  ! integers for getting input variables 
integer(kind=4) izml, isp, ilabs, nlabs ! grid # of bottom of mixed layer, # of caco3 species, # of labs simulation and total # of labs simulations  
integer(kind=4) :: file_tmp=100,file_ccflx=101,file_omflx=102,file_o2flx=103,file_dicflx=104,file_alkflx=105,file_ptflx=106  !  file #
integer(kind=4) :: file_err=107,file_ccflxes(nspcc), file_bound=108, file_totfrac=109, file_sigmly=110,file_sigmlyd=111  ! file #
integer(kind=4) :: file_sigbtm=112 ! file #
external dgesv  ! subroutine in BALS library 
logical :: oxic = .true.  ! oxic only model of OM degradation by Emerson (1985) 
logical :: anoxic = .true.  ! oxic-anoxic model of OM degradation by Archer (1991) 
logical :: nobio(nspcc+2) = .false.  ! no biogenic reworking assumed 
logical :: turbo2(nspcc+2) = .false.  ! random mixing 
logical :: labs(nspcc+2) = .false.  ! mixing info from LABS 
logical :: nonlocal(nspcc+2)  ! ON if assuming non-local mixing (i.e., if labs or turbo2 is ON)
logical :: flg_500,izox_calc_done
real(kind=8) dt_om_o2,error_o2min ,dti
integer(kind=4) iizox, iizox_errmin, w_save(nz)
integer(kind=4) :: itr_w_max = 20
real(kind=8) :: kom_ox(nz),kom_an(nz),kom_dum(nz,3)
logical :: warmup_done = .false.
logical :: all_oxic

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

print*
print*,'ccflxi    = ',ccflxi
print*,'om2cc     = ',om2cc
print*,'dep       = ',depi2
print*,'dt        = ',dti2
print*,'runname   = ',trim(adjustl(filechr))
print*,'oxonly?   = ',oxonly
print*,'biot      = ',trim(adjustl(biotmode))
print*
! print*,'If above variables are OK, please press [enter] to start'
! read *

select case(trim(adjustl(biotmode))) 
    case('nobio','NOBIO','Nobio')
        nobio = .true.
        print*
        print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        print*,'XXXXXXXXXXXXXXX No bioturbation XXXXXXXXXXXXX'
        print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        print*
    case('turbo2','TURBO2','Turbo2')
        turbo2 = .true.
        print*
        print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        print*,'XXXXXXXXXXX Homogeneous mixing XXXXXXXXXXXXXX'
        print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        print*
    case('labs','LABS','Labs')
        labs = .true. 
        print*
        print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        print*,'XXXXXXXXXXXXXXXX LABS mixing XXXXXXXXXXXXXXXX'
        print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        print*
    case default  
        print*
        print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        print*,'XXXXXXXXXXXXX Fickian mixing XXXXXXXXXXXXXXXX'
        print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        print*
endselect 

if (oxonly) anoxic = .false. 

flg_500 = .false.

! dt = dti2
dti = dti2
dep = depi2

call date_and_time(dumchr(1),dumchr(2),dumchr(3),dumint)  ! get date and time in character 

!!! get variables !!!
! get total caco3 flux, om/caco3 rain ratio and water depth   
! call getinput(ccflxi,om2cc,dti,filechr,dep,anoxic,nobio,labs,turbo2,nspcc)
print'(3A,3E11.3)','ccflxi','om2cc','dep:',ccflxi,om2cc, dep  ! printing read data 

! prepare directory to store result files 
call makeprofdir(  &  ! make profile files and a directory to store them 
    workdir   &
    ,filechr,anoxic,labs,turbo2,nobio,nspcc  &
    ,file_ptflx,file_ccflx,file_omflx,file_o2flx,file_dicflx,file_alkflx,file_err  &
    ,file_bound,file_totfrac,file_sigmly,file_sigmlyd,file_sigbtm,file_ccflxes  &
    ,ccflxi,dep,om2cc,chr &
    )
!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!  MAKING GRID !!!!!!!!!!!!!!!!! 
beta = 1.00000000005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
call makegrid(beta,nz,ztot,dz,z)
! stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! FUNDAMENTAL PARAMETERS !!!!!!!!!!!!!
!! assume steady state flux for detrital material 
call flxstat(  &
    omflx,detflx,ccflx  & ! output
    ,om2cc,ccflxi,mcc,nspcc  & ! input 
    )
    
    
!!!!!!!!! specific boundary conditions if needed

! o2i =  162.35434104144247
! dici = 2298.6964696813252
! alki = 2429.3198745074496
! dep =  5.0000000000000000 
! sal =  34.957461001340832
! temp =  1.3566180092130367
! ccflxi = 7.6549260517794933E-005
! omflx =  5.8994023885595817E-005
! detflx = 3.3969233213725432E-003
! ccflx = ccflxi/nspcc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! calculate molar volume (cm3/mol) = molar mass (g/mol) / density (g/cm3)
#ifdef isotrack
mcc(i12c16o) = 40.078d0+12d0+16d0*3d0
mcc(i12c18o) = 40.078d0+12d0+16d0*2d0+18d0
mcc(i13c16o) = 40.078d0+13d0+16d0*3d0
mcc(i13c18o) = 40.078d0+13d0+16d0*2d0+18d0
mcc(i14c) = 40.078d0+14d0+16d0*3d0
#endif 
mvom = mom/rhoom  ! om
mvsed = msed/rhosed ! clay 
mvcc = mcc/rhocc ! caco3
    
!  assume porosity profile 
call getporosity(  &
     poro,porof,sporo,sporof,sporoi & ! output
     ,z,nz,poroi  & ! input
     )

! below is the point when the calculation is re-started with smaller time step
500 continue

! initial guess of burial profile 
call burial_pre(  &
    w,wi  & ! output
    ,detflx,ccflx,nspcc,nz,poroi,msed,mvsed,mvcc  & ! input 
    )
! depth -age conversion 
call dep2age(  &
    age &  ! output 
    ,dz,w,nz  &  ! input
   )
! determine factors for upwind scheme to represent burial advection
call calcupwindscheme(  &
    up,dwn,cnr,adf & ! output 
    ,w,nz   & ! input &
    )
! ---------------------

zox = 10d0 ! priori assumed oxic zone 

!!! ~~~~~~~~~~~~~~ set recording time 
call recordtime(  &
    rectime,time_spn,time_trs,time_aft,cntrec  &
    ,ztot,wi,file_tmp,workdir,nrec  &
    )
    
! rectime(1) = 1000d0
! rectime(2) = 10000d0
! rectime(3) = 100000d0
! rectime(4) = 1000000d0
! rectime(5) = 10000000d0

! depi = 4d0  ! depth before event 
depi = 3.5d0  ! depth before event 
depf = dep   ! max depth to be changed to  

flxfini = 0.5d0  !  total caco3 rain flux for fine species assumed before event 
flxfinf = 0.9d0 !  maximum changed value 

! ///////////// isotopes  ////////////////
d13c_ocni = 2d0  ! initial ocean d13c value 
d13c_ocnf = -1d0 ! ocean d13c value with maximum change  
d18o_ocni = 1d0 ! initial ocean d18o value 
d18o_ocnf = -1d0 ! ocean d18o value with maximum change 
capd47_ocni = 0.6d0
capd47_ocnf = 0.5d0

call sig2sp_pre(  &  ! end-member signal assignment 
    d13c_sp,d18o_sp  &
    ,d13c_ocni,d13c_ocnf,d18o_ocni,d18o_ocnf,nspcc  &
    ) 
!! //////////////////////////////////
!!!~~~~~~~~~~~~~~~~~

!!!! TRANSITION MATRIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call make_transmx(  &
    trans,izrec,izrec2,izml,nonlocal  & ! output 
    ,labs,nspcc,turbo2,nobio,dz,sporo,nz,z,file_tmp,zml_ref  & ! input
    )

!~~~~~~~~diffusion & reaction~~~~~~~~~~~~~~~~~~~~~~
call coefs(  &
    dif_dic,dif_alk,dif_o2,kom,kcc,co3sat,krad & ! output 
    ,temp,sal,dep,nz,nspcc,poro,cai,komi,kcci,k14ci,i14c  & !  input 
    ,i13c18o  &
    )

!!   INITIAL CONDITIONS !!!!!!!!!!!!!!!!!!! 
cc = 1d-8   ! assume an arbitrary low conc. 
! cc = 0.9d0/(mcc/rhocc)   ! assume 90 wt%  
dic = dici*1d-6/1d3 ! mol/cm3; factor is added to change uM to mol/cm3 
alk = alki*1d-6/1d3 ! mol/cm3

! call subroutine to calculate all aqueous co2 species reflecting initial assumption on dic and alk 
#ifndef mocsy
call calcspecies(dic,alk,temp,sal,dep,pro,co2,hco3,co3,nz,infosbr)  
#else
call co2sys_mocsy(nz,alk*1d6,dic*1d6,temp,dep*1d3,sal  &
                        ,co2,hco3,co3,pro,ohmega,dohmega_ddic,dohmega_dalk) ! using mocsy
co2 = co2/1d6
hco3 = hco3/1d6
co3 = co3/1d6
#endif     
pt = 1d-8  ! assume an arbitrary low conc. 
om = 1d-8  ! assume an arbitrary low conc. 
o2 = o2i*1d-6/1d3 ! mol/cm3  ; factor is added to change uM to mol/cm3 
! o2(2:) = 0d0

! ~~~ passing to temporary variables with subscript x ~~~~~~~~~~~
ccx = cc
dicx = dic
alkx = alk 

co2x = co2
hco3x = hco3
co3x = co3

co3i=co3(1) ! recording seawater conc. of co3 
#ifdef mocsy 
cai = (0.02128d0/40.078d0) * sal/1.80655d0
co3sat = co3i*1d3/ohmega(1)
#endif  

ptx = pt

omx = om
o2x = o2

! calculating initial dissolution rate of caco3 for all caco3 species 
do isp=1,nspcc 
    rcc(:,isp) = kcc(:,isp)*ccx(:,isp)*abs(1d0-co3x(:)*1d3/co3sat)**ncc*merge(1d0,0d0,(1d0-co3x(:)*1d3/co3sat)>0d0)
enddo

oxco2 = 0d0  ! initial oxic degradation 
anco2 = 0d0  ! anoxic counterpart 

! ~~~ saving initial conditions 
#ifndef nonrec
call recordprofile(  &
    0,file_tmp,workdir,nz,z,age,pt,rho,cc,ccx,dic,dicx,alk,alkx,co3,co3x,nspcc,msed,wi,co3sat,rcc  &
    ,pro,o2x,oxco2,anco2,om,mom,mcc,d13c_ocni,d18o_ocni,up,dwn,cnr,adf,ptx,w,frt,prox,omx,d13c_blk,d18o_blk  &
    ,d17o_blk,d14c_age,capd47  &
    )
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~
!! START OF TIME INTEGLATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
warmup_done = .false.
time = 0d0
it = 1
do 

    !! ///////// isotopes & fluxes settings ////////////// 
#ifndef sense    
    nt_spn =400
    ! if (.not.warmup_done) then 
        ! if (it<11) then 
            ! nt_spn = 80000
        ! elseif (it<21) then
            ! nt_spn = 8000
        ! else 
            ! warmup_done = .true.
            ! time = 0d0
            ! it = 1
            ! cycle
        ! endif 
    ! endif 
    nt_trs = 5000
    nt_aft = 1000
    700 continue
    call timestep(time,nt_spn,nt_trs,nt_aft,dt,time_spn,time_trs,time_aft)
    call signal_flx(  &
        d13c_ocn,d18o_ocn,ccflx,d18o_sp,d13c_sp,cntsp  &
        ,time,time_spn,time_trs,time_aft,d13c_ocni,d13c_ocnf,d18o_ocni,d18o_ocnf,nspcc,ccflxi,it,flxfini,flxfinf  &
        ,r14ci,capd47_ocni,capd47_ocnf,capd47_ocn,r13c_pdb,r18o_pdb,r17o_pdb,tol,nt_trs  &
        ) 
    call bdcnd(   &
        time,dep,time_spn,time_trs,time_aft,depi,depf  &
        ) 

    ! isotope signals represented by caco3 rain fluxes 
    d18o_flx = sum(d18o_sp(:)*ccflx(:))/ccflxi
    d13c_flx = sum(d13c_sp(:)*ccflx(:))/ccflxi

#ifndef track2
#ifndef isotrack
    if (abs(d13c_flx - d13c_ocn)>tol .or. abs(d18o_flx - d18o_ocn)>tol) then ! check comparability with input signals 
        print*,'error in assignment of proxy'
        write(file_err,*)'error in assignment of proxy',d18o_ocn,d13c_ocn,d18o_flx,d13c_flx
        stop
    endif 
#endif 
#endif

#else 
    dt = dti
    600 continue 
#endif 

    !! === temperature & pressure and associated boundary changes ====
    ! if temperature is changed during signal change event this affect diffusion coeff etc. 
    call coefs(  &
        dif_dic,dif_alk,dif_o2,kom,kcc,co3sat,krad & ! output 
        ,temp,sal,dep,nz,nspcc,poro,cai,komi,kcci,k14ci,i14c  & !  input 
        ,i13c18o  &
        )
    !! /////////////////////

#ifndef nonrec 
    if (it==1) then    !! recording boundary conditions 
        if (.not. flg_500) then 
            write(file_bound,*) '#time  d13c_ocn  d18o_ocn, D47, fluxes of cc:',(isp,isp=1,nspcc) & 
                ,'temp  dep  sal  dici  alki  o2i'
        endif 
    endif 
#ifndef size 
    !  recording fluxes of two types of caco3 separately 
    if (.not. flg_500) write(file_bound,*) time, d13c_ocn, d18o_ocn, capd47_ocn,(ccflx(isp),isp=1,nspcc)  &
        ,temp, dep, sal,dici,alki, o2i
#else 
    !  do not record separately 
    if (.not. flg_500) write(file_bound,*) time, d13c_ocn, d18o_ocn, capd47_ocn,sum(ccflx(1:4)),sum(ccflx(5:8)) &
        ,(ccflx(isp),isp=1,nspcc),temp, dep, sal,dici,alki, o2i
#endif 
#endif 

    itr_w = 0  ! # of iteration made to converge w 
    err_w_min = 1d4 ! minimum relative difference of w compared to the previous w 

    !  the followings are currently not used !!
    itr_f = 0  ! # of iteration made to converge total vol. fraction of solids 
    err_f_min = 1d4 ! minimum relative difference of total vol. fraction of solids compared to the previous value 
    dfrt_df = 0d0 ! change in total volume fraction divided by change in variable f that is used to tune how advection is calculated where burial changes its signs 
    d2frt_df2 = 0d0  ! change of dfrt_dt by change in variable f 
    !  the above are currently not used !!!!

    err_f = 0d0  !! relative different of total vol. fraction of solids wrt the previous value 
    err_fx = 0d0 !! previous err_f
    
    w_save = w !! saving previous w 

    ! point where iteration for w is conducted 
    300 continue 

#ifndef nondisp
    ! displaying time step 
    print*,'it :',it,dt
#endif

    ! initializing
    dw = 0d0 ! change in burial rate caused by reaction and non-local mixing 

    ! ~~~~~~~~~ OM & O2 iteration wrt zox ~~~~~~~~~~~~~~~

    oxco2 = 0d0 ! oxic degradation of om 
    anco2 = 0d0 ! anoxic degradation of om 
    itr = 0  ! iteration # for om and o2 calcuation 
    error = 1d4 ! error in ieration for om and o2 
    minerr= 1d4  ! recording minimum relative difference in zox from previously considered zox 

    do !  while (error > tol)
        !~~~~~~ OM calculation ~~~~~~~~~~~~~~~~~
        dt_om_o2 = 1d8 
        dt_om_o2 = dt 
        
        ! calculating zox from assumed/previous o2 profiles
        call calc_zox( &
            izox,kom,zox,kom_ox,kom_an   &  ! output 
            ,oxic,anoxic,nz,o2x,o2th,komi,ztot,z,o2i,dz  & ! input
            )
        
        call omcalc( &
            omx  & ! output 
            ,kom   &  ! input
            ,om,nz,sporo,sporoi,sporof &! input 
            ,w,wi,dt_om_o2,up,dwn,cnr,adf,trans,nspcc,labs,turbo2,nonlocal,omflx,poro,dz &! input 
            ) 
        ! calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations )
        call calcflxom(  &
            omadv,omdec,omdif,omrain,omres,omtflx  & ! output 
            ,sporo,om,omx,dt_om_o2,w,dz,z,nz,turbo2,labs,nonlocal,poro,up,dwn,cnr,adf,rho,mom  &
            ,trans,kom,sporof,sporoi,wi,nspcc,omflx  & ! input 
            ,file_tmp,workdir &
            ,flg_500  &
            )
        if (flg_500) then 
            dt = dt/10d0
            w = w_save  
            call calcupwindscheme(  &
                up,dwn,cnr,adf & ! output 
                ,w,nz   & ! input &
                )            
            go to 300
        endif
        
        !~~~~~~~~~~~~~~~~~ O2 calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! if (izox == nz) then ! fully oxic; lower boundary condition ---> no diffusive out flow  
        call o2calc_ox(  &
            o2x  & ! output
            ,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt_om_o2,ox2om,o2i & ! input
            )
        !  fluxes relevant to o2 (at the same time checking the satisfaction of difference equations) 
        call calcflxo2_ox( &
            o2dec,o2dif,o2tflx,o2res  & ! output 
            ,nz,sporo,kom_ox,omx,dz,poro,dif_o2,dt_om_o2,o2,o2x,ox2om,o2i  & ! input
            )
        ! endif 
        ! else  !! if oxygen is depleted within calculation domain, lower boundary changes to zero concs.
            ! izox = nz
            ! call o2calc_ox(  &
                ! o2x  & ! output
                ! ,izox,nz,poro,o2,kom,omx,sporo,dif_o2,dz,dt_om_o2,ox2om,o2i & ! input
                ! )
            ! izox_calc_done = .false.
            ! if (oxic) then 
                ! do iz=1,nz
                    ! if (o2x(iz) > o2th) then
                        ! if (.not. izox_calc_done) izox = iz
                    ! else! unless anoxi degradation is allowed, om cannot degradate below zox
                        ! izox_calc_done = .true.
                    ! endif
                ! enddo
            ! endif
        ! print'(A,5E11.3)', 'o2 :',(o2x(iz)*1d3,iz=1,nz,nz/5)
        ! if (all(o2x>=0d0)) then 
        if (all(o2x>=0d0).and.izox==nz) then 
            iizox_errmin = nz
            ! print *,'all oxic',iizox_errmin
            all_oxic = .true.
        elseif (any(o2x<0d0)) then 
            all_oxic = .false.
            error_o2min = 1d4
            iizox_errmin = izox
            do iizox = 1,nz   
            ! do iizox = izox,nz   
            ! do iizox = max(1,izox-20),min(nz,izox+20)   
                if (iizox<nz) then 
                    call o2calc_sbox(  &
                        o2x  & ! output
                        ,iizox,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt_om_o2,ox2om,o2i & ! input
                        )
                elseif (izox==nz) then  
                    call o2calc_ox(  &
                        o2x  & ! output
                        ,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt_om_o2,ox2om,o2i & ! input
                        )
                endif 
                ! print'(A,I0,10E11.3)', 'o2 :',iizox,(o2x(iz)*1d3,iz=1,nz,nz/10)
                ! fluxes relevant to oxygen 
                ! call calcflxo2_sbox( &
                    ! o2dec,o2dif,o2tflx,o2res  & ! output 
                    ! ,nz,sporo,kom,omx,dz,poro,dif_o2,dt_om_o2,o2,o2x,iizox,ox2om,o2i  & ! input
                    ! )
                if (all(o2x>=0d0)) then 
                    if (abs(o2x(max(iizox-1,1)))<error_o2min) then 
                        error_o2min = abs(o2x(max(iizox-1,1)))
                        iizox_errmin = iizox
                        ! print*,'find smaller difference',iizox_errmin 
                    endif 
                endif 
            enddo 
            ! print*,iizox_errmin 
            if (iizox_errmin < nz) then 
                call o2calc_sbox(  &
                    o2x  & ! output
                    ,iizox_errmin,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt_om_o2,ox2om,o2i & ! input
                    )
                ! print'(A,I0,10E11.3)', 'o2 :',iizox_errmin,(o2x(iz)*1d3,iz=1,nz,nz/10)
                ! fluxes relevant to oxygen 
                call calcflxo2_sbox( &
                    o2dec,o2dif,o2tflx,o2res  & ! output 
                    ,nz,sporo,kom_ox,omx,dz,poro,dif_o2,dt_om_o2,o2,o2x,iizox_errmin,ox2om,o2i  & ! input
                    )
            elseif (iizox_errmin ==nz) then 
                call o2calc_ox(  &
                    o2x  & ! output
                    ,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt_om_o2,ox2om,o2i & ! input
                    )
                !  fluxes relevant to o2 (at the same time checking the satisfaction of difference equations) 
                call calcflxo2_ox( &
                    o2dec,o2dif,o2tflx,o2res  & ! output 
                    ,nz,sporo,kom_ox,omx,dz,poro,dif_o2,dt_om_o2,o2,o2x,ox2om,o2i  & ! input
                    )
            endif 
                
            ! izox_calc_done = .false.
            ! if (oxic) then 
                ! do iz=1,nz
                    ! if (o2x(iz) > o2th) then
                        ! if (.not. izox_calc_done) iizox_errmin = iz
                    ! else! unless anoxi degradation is allowed, om cannot degradate below zox
                        ! izox_calc_done = .true.
                    ! endif
                ! enddo
            ! endif
            
            call calc_zox( &
                iizox_errmin,kom_dum(:,1),zox,kom_dum(:,2),kom_dum(:,3)   &  ! output 
                ,oxic,anoxic,nz,o2x,o2th,komi,ztot,z,o2i,dz  & ! input
                )
            ! print*,iizox_errmin             
        endif

        ! update of zox 
        ! zoxx = 0d0
        ! do iz=1,nz
            ! if (o2x(iz)<=0d0) exit
        ! enddo

        ! if (iz==nz+1) then ! oxygen never gets less than 0 
            ! zoxx = ztot ! zox is the bottom depth 
        ! else if (iz==1) then ! calculating zox interpolating at z=0 with SWI conc. and at z=z(iz) with conc. o2x(iz)
            ! zoxx = (z(iz)*o2i*1d-6/1d3 + 0d0*abs(o2x(iz)))/(o2i*1d-6/1d3+abs(o2x(iz)))
        ! else     ! calculating zox interpolating at z=z(iz-1) with o2x(iz-1) and at z=z(iz) with conc. o2x(iz)
            ! zoxx = (z(iz)*o2x(iz-1) + z(iz-1)*abs(o2x(iz)))/(o2x(iz-1)+abs(o2x(iz)))
        ! endif

        ! error = abs((zox -zoxx)/zox)  ! relative difference 
        error = abs(izox-iizox_errmin)  ! relative difference 
#ifdef showiter 
        print*, 'zox',itr,izox, iizox_errmin
#endif
        ! if (zox==zoxx) exit 
         
        ! zox = 0.5d0*(zox + zoxx)  ! new zox 

        ! if iteration reaches 100 error in zox is tested assuming individual grid depths as zox and find where error gets minimized 
        ! if (itr>=100 .and. itr <= nz+99) then 
            ! zox = z(itr-99) ! zox value in next test 
            ! if (minerr >=error ) then ! if this time error is less than last adopt as optimum 
                ! if (itr/=100) then 
                    ! izox_minerr = itr -100
                    ! minerr = error 
                ! endif 
            ! endif
        ! elseif (itr ==nz+100) then ! check last test z(nz)
            ! if (minerr >=error ) then 
                ! izox_minerr = itr -100
                ! minerr = error 
            ! endif
            ! zox = z(izox_minerr)  ! determine next test which should be most optimum 
        ! elseif (itr ==nz+101) then  ! results should be optimum and thus exit 
            ! exit
        ! endif 

        ! if (itr >nz+101) then 
            ! stop
        ! endif
        if (izox==iizox_errmin) then 
            if (all_oxic) then            
                exit 
            else 
                if (izox < nz) then 
                    call o2calc_sbox(  &
                        o2x  & ! output
                        ,izox,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt_om_o2,ox2om,o2i & ! input
                        )
                    ! print'(A,I0,10E11.3)', 'o2 :',izox,(o2x(iz)*1d3,iz=1,nz,nz/10)
                    ! fluxes relevant to oxygen 
                    call calcflxo2_sbox( &
                        o2dec,o2dif,o2tflx,o2res  & ! output 
                        ,nz,sporo,kom_ox,omx,dz,poro,dif_o2,dt_om_o2,o2,o2x,izox,ox2om,o2i  & ! input
                        )       
                elseif (izox == nz) then  
                    call o2calc_ox(  &
                        o2x  & ! output
                        ,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt_om_o2,ox2om,o2i & ! input
                        )
                    !  fluxes relevant to o2 (at the same time checking the satisfaction of difference equations) 
                    call calcflxo2_ox( &
                        o2dec,o2dif,o2tflx,o2res  & ! output 
                        ,nz,sporo,kom_ox,omx,dz,poro,dif_o2,dt_om_o2,o2,o2x,ox2om,o2i  & ! input
                        )
                endif 
                exit 
            endif  
        endif 
        
        if (error < minerr ) then 
            minerr = error 
        else 
            if (izox < nz .and. iizox_errmin == nz) then 
            
                call o2calc_sbox(  &
                    o2x  & ! output
                    ,izox,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt_om_o2,ox2om,o2i & ! input
                    )
                ! print'(A,I0,10E11.3)', 'o2 :',izox,(o2x(iz)*1d3,iz=1,nz,nz/10)
                ! fluxes relevant to oxygen 
                call calcflxo2_sbox( &
                    o2dec,o2dif,o2tflx,o2res  & ! output 
                    ,nz,sporo,kom_ox,omx,dz,poro,dif_o2,dt_om_o2,o2,o2x,izox,ox2om,o2i  & ! input
                    )
                exit 
                
            endif 
        endif 
        
        itr = itr + 1

    enddo 

    ! print'(A,I0,10E11.3)', 'o2 :',izox,(o2x(iz)*1d3,iz=1,nz,nz/10)
    !~~  OM & O2 calculation END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! do iz = 1,nz 
        ! if (o2x(iz) > o2th) then
            ! oxco2(iz) = (1d0-poro(iz))*kom(iz)*omx(iz)  ! aerobic respiration 
        ! else 
            ! if (anoxic) then 
                ! anco2(iz) = (1d0-poro(iz))*kom(iz)*omx(iz)  ! anaerobic respiration 
            ! endif
        ! endif
    ! enddo
    
    oxco2(:) = (1d0-poro(:))*kom_ox(:)*omx(:)
    anco2(:) = (1d0-poro(:))*kom_an(:)*omx(:)

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

    !!  ~~~~~~~~~~~~~~~~~~~~~~ CaCO3 solid, ALK and DIC  calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call calccaco3sys(  &
        ccx,dicx,alkx,rcc,dt  & ! in&output
        ,nspcc,dic,alk,dep,sal,temp,labs,turbo2,nonlocal,sporo,sporoi,sporof,poro,dif_alk,dif_dic & ! input
        ,w,up,dwn,cnr,adf,dz,trans,cc,oxco2,anco2,co3sat,kcc,ccflx,ncc,ohmega,nz  & ! input
        ! ,dum_sfcsumocn  & ! input for genie geochemistry
        ,tol,poroi,flg_500,fact,file_tmp,alki,dici,ccx_th,workdir  &
        ,krad,deccc  & 
        )
    if (flg_500) then 
        dt = dt/10d0
        ! go to 500
        ccx = cc
        alkx = alk
        dicx = dic
        w = w_save 
        call calcupwindscheme(  &
            up,dwn,cnr,adf & ! output 
            ,w,nz   & ! input &
            )
#ifdef sense
        go to 600
#else 
        go to 700
#endif 
    endif 
    ! ~~~~  End of calculation iteration for CO2 species ~~~~~~~~~~~~~~~~~~~~
    ! update aqueous co2 species 
#ifndef mocsy
    call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,infosbr)
    if (infosbr==1) then 
        dt=dt/10d0
#ifdef sense
        ! go to 500
        ccx = cc
        alkx = alk
        dicx = dic
        w = w_save 
        call calcupwindscheme(  &
            up,dwn,cnr,adf & ! output 
            ,w,nz   & ! input &
            )
        go to 600
#else
        stop
#endif 
    endif 
#else 
    call co2sys_mocsy(nz,alkx*1d6,dicx*1d6,temp,dep*1d3,sal  &
                            ,co2x,hco3x,co3x,prox,ohmega,dohmega_ddic,dohmega_dalk) ! using mocsy
    co2x = co2x/1d6
    hco3x = hco3x/1d6
    co3x = co3x/1d6
#endif 

    ! calculation of fluxes relevant to caco3 and co2 system
    call calcflxcaco3sys(  &
         cctflx,ccflx,ccdis,ccdif,ccadv,ccrain,ccres,alktflx,alkdis,alkdif,alkdec,alkres & ! output
         ,dictflx,dicdis,dicdif,dicres,dicdec   & ! output
         ,dw & ! inoutput
         ,nspcc,ccx,cc,dt,dz,rcc,adf,up,dwn,cnr,w,dif_alk,dif_dic,dic,dicx,alk,alkx,oxco2,anco2,trans    & ! input
         ,turbo2,labs,nonlocal,sporof,it,nz,poro,sporo        & ! input
         ,dici,alki,file_err,mvcc,tol,flg_500  &
         ,ccrad,alkrad,deccc  &  
         )
    ! if (flg_500) go to 500
    
    ! ~~~~ calculation clay  ~~~~~~~~~~~~~~~~~~
    call claycalc(  &   
        ptx                  &  ! output
        ,nz,sporo,pt,dt,w,dz,detflx,adf,up,dwn,cnr,trans  &  ! input
        ,nspcc,labs,turbo2,nonlocal,poro,sporof     &  !  intput
        ,msed,file_tmp,workdir &
        )
    call calcflxclay( &
        pttflx,ptdif,ptadv,ptres,ptrain  & ! output
        ,dw          &  ! in&output
        ,nz,sporo,ptx,pt,dt,dz,detflx,w,adf,up,dwn,cnr,sporof,trans,nspcc,turbo2,labs,nonlocal,poro           &  !  input
        ,msed,mvsed  &
        )
    !! ~~~~~~~~~End of clay calculation 

    err_fx = maxval(abs(frt - 1d0))  ! recording previous error in total vol. fraction of solids 

    ! get solid property, rho (density) and frt (total vol.frac)
    call getsldprop(  &
        rho,frt,       &  ! output
        nz,omx,ptx,ccx,nspcc,w,up,dwn,cnr,adf,z      & ! input
        ,mom,msed,mcc,mvom,mvsed,mvcc,file_tmp,workdir  &
        )

    err_f = maxval(abs(frt - 1d0))  ! new error in total vol. fraction (must be 1 in theory) 
    if (err_f < err_fx) err_f_min = err_f  ! recording minimum error 
#ifdef sense
    if (err_f < tol) exit  ! if total vol. fraction is near enough to 1, steady-state solution is obtained 
#endif 
    !! calculation of burial velocity =============================

    wx = w  ! recording previous burial velocity 

    !  get new burial velocity
    call burialcalc(  &
        w,wi        & !  output
        ,detflx,ccflx,nspcc,omflx,dw,dz,poro,nz    & ! input
        ,msed,mvsed,mvcc,mvom,poroi &
        )
    ! call burialcalc_fdm(  &
        ! w       & !  output
        ! ,detflx,ccflx,nspcc,omflx,dw,dz,poro,nz    & ! input
        ! ,msed,mvsed,mvcc,mvom,poroi,up,dwn,cnr,adf,sporo,sporof,file_tmp,workdir &
        ! )

    ! ------------ determine calculation scheme for advection 
    call calcupwindscheme(  &
        up,dwn,cnr,adf & ! output 
        ,w,nz   & ! input &
        )

    ! error and iteration evaluation 
    itr_w = itr_w + 1  ! counting iteration for w 
    err_w = maxval(abs((w-wx)/wx))  ! relative difference of w 
#ifdef showiter
    print*,'itr_w,err_w',itr_w,err_w 
#endif 
    if (err_w<err_w_min) then 
        err_w_min= err_w  ! recording minimum relative difference of  w 
        wxx = wx  ! recording w which minimizes deviation of total sld fraction from 1 
    endif 
    if (itr_w>itr_w_max) then   ! if iteration gets too many, force to end with optimum w where error is minimum
        if (itr_w==itr_w_max+1) then 
            w = wxx   
            go to 300
        elseif (itr_w==itr_w_max+2) then 
            w = wxx
            write(file_err,*) 'not converging w',time, err_w, err_w_min
            go to 400
        endif 
    endif
    if (err_w > tol) go to 300

    400 continue

#ifndef nonrec
    if (it==1) write(file_omflx,*)'time, omtflx, omadv, omdec, omdif, omrain, omres'
    write(file_omflx,*)time, omtflx, omadv, omdec, omdif, omrain, omres
    if (it==1) write(file_o2flx,*)'time, o2dec, o2dif, o2tflx, o2res'
    write(file_o2flx,*)time,o2dec, o2dif,o2tflx,o2res
    
    if (it==1) then 
        write(file_ccflx,*) 'time, cctflx, ccdis, ccrad, ccdif, ccadv, ccrain, ccres' 
        write(file_dicflx,*) 'time, dictflx, dicdis, dicdif, dicdec,  dicres' 
        write(file_alkflx,*) 'time, alktflx, alkdis, alkrad, alkdif, alkdec, alkres' 
        do isp=1,nspcc
            write(file_ccflxes(isp),*) 'time, cctflx, ccdis, ccrad, ccdif, ccadv, ccrain, ccres' 
        enddo
    endif
    write(file_ccflx,*) time,sum(cctflx), sum(ccdis), sum(ccrad), sum(ccdif), sum(ccadv), sum(ccrain), sum(ccres)
    write(file_dicflx,*) time,dictflx, dicdis, dicdif, dicdec,  dicres 
    write(file_alkflx,*) time,alktflx, alkdis, alkrad, alkdif, alkdec, alkres 
    do isp=1,nspcc
        write(file_ccflxes(isp),*) time,cctflx(isp), ccdis(isp), ccrad(isp),ccdif(isp), ccadv(isp), ccrain(isp), ccres(isp)
    enddo
    
    if (it==1) write(file_ptflx,*) 'time, pttflx, ptdif, ptadv, ptrain, ptres'
    write(file_ptflx,*) time, pttflx, ptdif, ptadv, ptrain, ptres
#endif

    !! depth -age conversion 
    ! call dep2age()
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
!!!!!  direct tracking 
#ifdef isotrack 
        r18o_blk(iz) = sum((/ccx(iz,i12c18o),ccx(iz,i13c18o)/))  &
            /sum((/3d0*ccx(iz,i12c16o),3d0*ccx(iz,i13c16o),2d0*ccx(iz,i12c18o),2d0*ccx(iz,i13c18o)/))
        r13c_blk(iz) = sum((/ccx(iz,i13c16o),ccx(iz,i13c18o)/))  &
            /sum((/ccx(iz,i12c18o),ccx(iz,i12c16o)/))
        r17o_blk(iz) = 0d0
        d18o_blk(iz) = r2d(r18o_blk(iz),r18o_pdb)
        d17o_blk(iz) = r2d(r17o_blk(iz),r17o_pdb)
        d13c_blk(iz) = r2d(r13c_blk(iz),r13c_pdb)
        d14c_age(iz) = -8033d0*log(ccx(iz,i14c)   &
            /sum((/ccx(iz,i12c18o),ccx(iz,i12c16o)/)) &
            /r14ci) ! Stuiver and Polach (1977)
        
        r47 = (ccx(iz,i13c18o))/ccx(iz,i12c16o)
        r47s = r13c_blk(iz)*r18o_blk(iz) 
        
        capd47(iz) = ((r47/r47s-1d0) )*1d3

#ifdef fullclump
        !!! not yet implemented
        r18o_blk(iz) = sum((/ccx(iz,i12c18o),ccx(iz,i13c18o),ccx(iz,i12c17o18o),1d0*ccx(iz,i12c18o18o)/))  &
            /sum((/ccx(iz,i12c18o),1d0*ccx(iz,i12c16o),ccx(iz,i12c17o) &
            ,1d0*ccx(iz,i13c16o),ccx(iz,i13c17o),ccx(iz,i13c18o)/))
        r17o_blk(iz) = sum((/ccx(iz,i12c17o),ccx(iz,i13c17o),ccx(iz,i12c17o18o),1d0*ccx(iz,i12c17o17o)/))  &
            /sum((/ccx(iz,i12c18o),1d0*ccx(iz,i12c16o),ccx(iz,i12c17o) &
            ,1d0*ccx(iz,i13c16o),ccx(iz,i13c17o),ccx(iz,i13c18o)/))
        r13c_blk(iz) = sum((/ccx(iz,i13c16o),ccx(iz,i13c18o),ccx(iz,i13c17o)/))  &
            /sum((/ccx(iz,i12c18o),ccx(iz,i12c16o),ccx(iz,i12c17o) &
            ,ccx(iz,i12c17o18o),ccx(iz,i12c17o17o),ccx(iz,i12c18o18o)/)) 
        d18o_blk(iz) = r2d(r18o_blk(iz),r18o_pdb)
        d17o_blk(iz) = r2d(r17o_blk(iz),r17o_pdb)
        d13c_blk(iz) = r2d(r13c_blk(iz),r13c_pdb)
        d14c_age(iz) = -8033d0*log(ccx(iz,i14c)   &
            /sum((/ccx(iz,i12c18o),ccx(iz,i12c16o),ccx(iz,i12c17o) &
            ,ccx(iz,i12c17o18o),ccx(iz,i12c17o17o),ccx(iz,i12c18o18o)/)) &
            /r14ci) ! Stuiver and Polach (1977)
        r45 = (ccx(iz,i13c16o)+ccx(iz,i12c17o))/ccx(iz,i12c16o)
        r45s = r13c_blk(iz) + 1d0*r17o_blk(iz)
        r46 = (ccx(iz,i13c17o)+ccx(iz,i12c18o)+ccx(iz,i12c17o17o))/ccx(iz,i12c16o)
        r46s = 1d0*r13c_blk(iz)*r17o_blk(iz) +1d0*r18o_blk(iz) + r17o_blk(iz)**2d0
        r47 = (ccx(iz,i13c18o)+ccx(iz,i12c17o18o))/ccx(iz,i12c16o)
        r47s = 1d0*r13c_blk(iz)*r18o_blk(iz) + 1d0*r18o_blk(iz)*r17o_blk(iz)
        capd47(iz) = ((r47/r47s-1d0) - (r46/r46s-1d0) - (r45/r45s-1d0))*1d3
#endif 
#endif 

    enddo

    !!!!! PRINTING RESULTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef nonrec
    if (time>=rectime(cntrec)) then
#ifndef nondisp
        print*,'****** TIME TO RECORD *********'
        print*,time,cntrec
        call recordprofile(  &
            cntrec,file_tmp,workdir,nz,z,age,pt,rho,cc,ccx,dic,dicx,alk,alkx,co3,co3x,nspcc,msed,wi,co3sat,rcc  &
            ,pro,o2x,oxco2,anco2,om,mom,mcc,d13c_ocni,d18o_ocni,up,dwn,cnr,adf,ptx,w,frt,prox,omx,d13c_blk,d18o_blk  &
            ,d17o_blk,d14c_age,capd47  &
            )
        print*,'****** RECORD FINISHED ********'
        print*
        print*
        print*
#endif 
        cntrec = cntrec + 1
        if (cntrec == nrec+1) exit
    endif 
#endif

    ! displaying results 
#ifndef nondisp    
    call resdisplay(  &
        nz,nspcc,it &
        ,z,frt,omx,rho,o2x,dicx,alkx,ptx,w &
        ,ccx  &
        ,cctflx,ccadv,ccdif,ccdis,ccrain,ccres  &
        ,time,omtflx,omadv,omdif,omdec,omrain,omres,o2tflx,o2dif,o2dec,o2res  &
        ,dictflx,dicdif,dicdec,dicdis,dicres,alktflx,alkdif,alkdec,alkdis,alkres  &
        ,pttflx,ptadv,ptdif,ptrain,ptres,mom,mcc,msed  &
        ,alkrad,ccrad  &
        ) 
#endif

    !! in theory, o2dec/ox2om + alkdec = dicdec = omdec (in absolute value)
    if (om2cc /= 0d0) then 
        if ( abs((o2dec/ox2om - alkdec + dicdec)/dicdec) > tol) then 
            print*, 'chk_om_dec',abs((o2dec/ox2om + alkdec - dicdec)/dicdec) ,o2dec/ox2om,alkdec,dicdec
            write(file_err,*) trim(adjustl(dumchr(1))), time, dt &
                , abs((o2dec/ox2om + alkdec - dicdec)/dicdec),o2dec/ox2om,alkdec,dicdec
        endif
    endif 

#ifndef nonrec 
    write(file_totfrac,*) time,maxval(abs(frt - 1d0))
    ! recording signals at 3 different depths (btm of mixed layer, 2xdepths of btm of mixed layer and btm depth of calculation domain)
    call sigrec(  &
        nz,file_sigmly,file_sigmlyd,file_sigbtm,w,time,age,izrec,d13c_blk,d13c_blkc &
        ,d13c_blkf,d18o_blk,d18o_blkc,d18o_blkf,ccx,mcc,rho,ptx,msed,izrec2,nspcc  & 
        ,d14c_age,capd47   &  
        )
#endif 

    ! before going to next time step, update variables 

    time = time + dt
    it = it + 1

    o2 = o2x
    om = omx

    cc = ccx
    dic = dicx
    alk = alkx

    pt = ptx
    
    ! if (time>=100d0) exit 

enddo

#ifndef nonrec
open(unit=file_tmp,file=trim(adjustl(workdir))//'sp-trace.txt',action='write',status='replace') 
do isp  = 1,nspcc
    write(file_tmp,*) isp,d13c_sp(isp),d18o_sp(isp)
enddo
close(file_tmp)
#endif 

#ifndef nonrec
call closefiles(  &
    file_ptflx,file_ccflx,file_omflx,file_o2flx,file_dicflx,file_alkflx,file_err  &
    ,file_bound,file_totfrac,file_sigmly,file_sigmlyd,file_sigbtm,file_ccflxes,nspcc  &
    )
#endif

! recording end results for lysoclines and caco3 burial fluxes
call resrec(  &
    workdir,anoxic,nspcc,labs,turbo2,nobio,co3i,co3sat,mcc,ccx,nz,rho,frt,ccadv,file_tmp,izml,chr,dt,it,time  &
    )

endsubroutine caco3
!**************************************************************************************************************************************




!**************************************************************************************************************************************
subroutine makegrid(beta,nz,ztot,dz,z)  !  making grid, after Hoffmann & Chiang, 2000
implicit none
integer(kind=4),intent(in) :: nz
real(kind=8),intent(in)::beta,ztot
real(kind=8),intent(out)::dz(nz),z(nz)
integer(kind=4) iz

do iz = 1, nz 
    z(iz) = iz*ztot/nz  ! regular grid 
    if (iz==1) then
        dz(iz) = ztot*log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0))
    endif
    if (iz/=1) then 
        dz(iz) = ztot*log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0)) - sum(dz(:iz-1))
    endif
enddo

do iz=1,nz  ! depth is defined at the middle of individual layers 
    if (iz==1) z(iz)=dz(iz)*0.5d0  
    if (iz/=1) z(iz) = z(iz-1)+dz(iz-1)*0.5d0 + 0.5d0*dz(iz)
enddo

!~~~~~~~~~~~~~ saving grid for LABS ~~~~~~~~~~~~~~~~~~~~~~
#ifdef recgrid
open(unit=100, file='C:/cygwin64/home/YK/LABS/1dgrid.txt',action='write',status='unknown')
do iz = 1, nz
    write(100,*) dz(iz)
enddo
close(100)
#endif

endsubroutine makegrid
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine getinput(ccflxi,om2cc,dti,filechr,dep,anoxic,nobio,labs,turbo2,nspcc)
implicit none 
integer(kind=4),intent(in)::nspcc
integer(kind=4) narg,ia
character*25 arg
real(kind=8),intent(out)::ccflxi,om2cc,dti,dep
logical,intent(out)::anoxic,nobio(nspcc+2),labs(nspcc+2),turbo2(nspcc+2)
character*255,intent(out):: filechr
logical ox_tmp 
character*10:: biotchr

narg = iargc()
do ia = 1, narg,2
    call getarg(ia,arg)
    select case(trim(arg))
        case('cc','CC','Cc')
            call getarg(ia+1,arg)
            read(arg,*)ccflxi   ! reading caco3 total rain flux
        case('rr','RR','Rr')
            call getarg(ia+1,arg)
            read(arg,*)om2cc  ! reading om/caco3 rain ratio 
        case('dep','DEP','Dep')
            call getarg(ia+1,arg)
            read(arg,*)dep  ! reading water depth in km 
        case('dt','DT','Dt')
            call getarg(ia+1,arg)
            read(arg,*)dti   ! reading time step used in calculation 
        case('fl','FL','Fl')
            call getarg(ia+1,arg)
            read(arg,*)filechr  ! reading file name that store sediment profiles etc. 
        case('ox','OX','Ox')
            call getarg(ia+1,arg)
            read(arg,*)ox_tmp  ! reading file name that store sediment profiles etc.
            anoxic = .not. ox_tmp
        case('biot','BIOT','Biot')
            call getarg(ia+1,arg)
            read(arg,*)biotchr  ! character to define styles of bioturbation
            select case(trim(biotchr))
                case('labs','LABS','Labs')
                    print*
                    print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    print*,'XXXXXXXXXXXXXXXX LABS mixing XXXXXXXXXXXXXXXX'
                    print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    print*
                    labs = .true.
                case('nobio','NOBIO','Nobio')
                    print*
                    print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    print*,'XXXXXXXXXXXXXXX No bioturbation XXXXXXXXXXXXX'
                    print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    print*
                    nobio=.true.
                case('turbo2','TURBO2','Turbo2')
                    print*
                    print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    print*,'XXXXXXXXXXX Homogeneous mixing XXXXXXXXXXXXXX'
                    print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    print*
                    turbo2 = .true.
                case default 
                    print*
                    print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    print*,'XXXXXXXXXXXXX Fickian mixing XXXXXXXXXXXXXXXX'
                    print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    print*
            endselect
    end select
enddo

! stop

endsubroutine getinput 
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine getinput_v2(ccflxi,om2cc,dti,filechr,dep,oxonly,biotchr)
implicit none 
integer(kind=4) narg,ia
character*25 arg
real(kind=8),intent(out)::ccflxi,om2cc,dti,dep
logical,intent(out)::oxonly
character*255,intent(out):: filechr,biotchr

narg = iargc()
do ia = 1, narg,2
    call getarg(ia,arg)
    select case(trim(arg))
        case('cc','CC','Cc')
            call getarg(ia+1,arg)
            read(arg,*)ccflxi   ! reading caco3 total rain flux
        case('rr','RR','Rr')
            call getarg(ia+1,arg)
            read(arg,*)om2cc  ! reading om/caco3 rain ratio 
        case('dep','DEP','Dep')
            call getarg(ia+1,arg)
            read(arg,*)dep  ! reading water depth in km 
        case('dt','DT','Dt')
            call getarg(ia+1,arg)
            read(arg,*)dti   ! reading time step used in calculation 
        case('fl','FL','Fl')
            call getarg(ia+1,arg)
            read(arg,*)filechr  ! reading file name that store sediment profiles etc. 
        case('ox','OX','Ox')
            call getarg(ia+1,arg)
            read(arg,*)oxonly  ! reading file name that store sediment profiles etc.
        case('biot','BIOT','Biot')
            call getarg(ia+1,arg)
            read(arg,*)biotchr  ! character to define styles of bioturbation
    end select
enddo

! stop

endsubroutine getinput_v2 
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine makeprofdir(  &  ! make profile files and a directory to store them 
    workdir   &
    ,filechr,anoxic,labs,turbo2,nobio,nspcc  &
    ,file_ptflx,file_ccflx,file_omflx,file_o2flx,file_dicflx,file_alkflx,file_err  &
    ,file_bound,file_totfrac,file_sigmly,file_sigmlyd,file_sigbtm,file_ccflxes  &
    ,ccflxi,dep,om2cc,chr &
    )
implicit none 
integer(kind=4),intent(in)::nspcc
integer(kind=4),intent(in)::file_ptflx,file_ccflx,file_omflx,file_o2flx,file_dicflx,file_alkflx,file_err
integer(kind=4),intent(in)::file_bound,file_totfrac,file_sigmly,file_sigmlyd,file_sigbtm
integer(kind=4),intent(inout)::file_ccflxes(nspcc)
real(kind=8),intent(in)::ccflxi,dep,om2cc
character*512,intent(inout)::workdir
character*255,intent(in)::filechr
logical,intent(in)::anoxic
logical,dimension(nspcc+2),intent(in)::labs,turbo2,nobio
real(kind=8) dumreal
integer(kind=4) ia,isp
character*25 dumchr(3),chr_tmp
character*25,intent(out)::chr(3,4)

do ia = 1,3  !  creating file name based on read caco3 rain flux, rain ratio and water depth 
    if (ia==1) dumreal=ccflxi  
    if (ia==2) dumreal=om2cc
    if (ia==3) dumreal=dep
    ! if (dumreal/=0d0) then 
        ! write(chr(ia,1),'(i0)') floor(log10(dumreal))  !! reading order 
        ! write(chr(ia,2),'(i0)') int(dumreal/(10d0**(floor(log10(dumreal))))) !  first digit 
        ! write(chr(ia,3),'(i0)') int((dumreal/(10d0**(floor(log10(dumreal)))) &         
            ! - int(dumreal/(10d0**(floor(log10(dumreal))))))*10d0)   ! second digit 
    ! else   ! if value is 0 set every character 0 
        ! write(chr(ia,1),'(i0)') 0 
        ! write(chr(ia,2),'(i0)') 0
        ! write(chr(ia,3),'(i0)') 0
    ! endif
    write(chr_tmp,'(E7.2e1)') dumreal 
    ! chr(ia,4) = trim(adjustl(chr(ia,2)))//'_'//trim(adjustl(chr(ia,3)))//'E'//trim(adjustl(chr(ia,1)))  
    ! chr_tmp(2:2)='_'
    chr(ia,4) = trim(adjustl(chr_tmp))
    ! x_yEz where x,y and z are read first and second digits and z is the order 
enddo 
print'(6A)','ccflx','om2cc','dep:',(chr(ia,4),ia=1,3)
! pause
!! FILES !!!!!!!!!
workdir = 'C:/Users/YK/Desktop/Sed_res/'
workdir = trim(adjustl(workdir))//'caco3_output/profiles/'
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
#ifndef nonrec
call system ('mkdir -p '//trim(adjustl(workdir)))
workdir = trim(adjustl(workdir))//'/'
open(unit=file_ptflx,file=trim(adjustl(workdir))//'ptflx.txt',action='write',status='unknown') ! recording fluxes of clay
open(unit=file_ccflx,file=trim(adjustl(workdir))//'ccflx.txt',action='write',status='unknown')! recording fluxes of caco3
open(unit=file_omflx,file=trim(adjustl(workdir))//'omflx.txt',action='write',status='unknown')! recording fluxes of om
open(unit=file_o2flx,file=trim(adjustl(workdir))//'o2flx.txt',action='write',status='unknown')! recording fluxes of o2
open(unit=file_dicflx,file=trim(adjustl(workdir))//'dicflx.txt',action='write',status='unknown')! recording fluxes of dic
open(unit=file_alkflx,file=trim(adjustl(workdir))//'alkflx.txt',action='write',status='unknown')! recording fluxes of alk
open(unit=file_err,file=trim(adjustl(workdir))//'errlog.txt',action='write',status='unknown')! recording errors 
open(unit=file_bound,file=trim(adjustl(workdir))//'bound.txt',action='write',status='unknown')! recording boundary conditions changes 
open(unit=file_totfrac,file=trim(adjustl(workdir))//'frac.txt',action='write',status='unknown') ! recording total fractions of solids 
open(unit=file_sigmly,file=trim(adjustl(workdir))//'sigmly.txt',action='write',status='unknown')! recording signals etc at just below mixed layer 
open(unit=file_sigmlyd,file=trim(adjustl(workdir))//'sigmlyd.txt',action='write',status='unknown') ! recording signals etc at depths of 2x mixed layer thickness 
open(unit=file_sigbtm,file=trim(adjustl(workdir))//'sigbtm.txt',action='write',status='unknown')! ! recording signals etc at bottom of sediment  
do isp = 1,nspcc
    file_ccflxes(isp)=800+(isp-1)  ! assigning intergers to files that record fluxes of individual caco3 species 
    write(dumchr(1),'(i3.3)') isp 
    open(unit=file_ccflxes(isp),file=trim(adjustl(workdir))//'ccflx-sp_'//trim(adjustl(dumchr(1))) &
        //'.txt',action='write',status='unknown')
enddo
#endif 
endsubroutine makeprofdir 
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine flxstat(  &
    omflx,detflx,ccflx  & ! output
    ,om2cc,ccflxi,mcc,nspcc  & ! input 
    )
implicit none 
integer(kind=4),intent(in)::nspcc
real(kind=8),intent(in)::om2cc,ccflxi,mcc(nspcc)
real(kind=8),intent(out)::omflx,detflx,ccflx(nspcc)

omflx = om2cc*ccflxi  ! om rain = rain ratio x caco3 rain 
! detflx = (1d0/9d0)*ccflxi*mcc ! 90% of mass flux becomes inorganic C; g cm-2 yr-1
detflx = (1d0/9d0)*ccflxi*100d0 ! 90% of mass flux becomes inorganic C; g cm-2 yr-1
ccflx = ccflxi/nspcc  !  rains of individual caco3 species is equivalently distributed as default 

endsubroutine flxstat
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine getporosity(  &
     poro,porof,sporo,sporof,sporoi & ! output
     ,z,nz,poroi  & ! input
     )
implicit none
integer(kind=4),intent(in)::nz
real(kind=8),dimension(nz),intent(in)::z
real(kind=8),dimension(nz),intent(out)::poro,sporo
real(kind=8),intent(in)::poroi 
real(kind=8),intent(out)::porof,sporoi,sporof
real(kind=8) calgg,pore_max,exp_pore

poro = poroi  ! constant porosity 
! ----------- Archer's parameterization 
calgg = 0.0d0  ! caco3 in g/g (here 0 is assumed )
pore_max =  1d0 - ( 0.483d0 + 0.45d0 * calgg) / 2.5d0  ! porosity at the bottom 
exp_pore = 0.25d0*calgg + 3.d0 *(1d0-calgg)  ! scale depth of e-fold decrease of porosity 
poro = EXP(-z/exp_pore) * (1.d0-pore_max) + pore_max 
porof = pore_max  ! porosity at the depth 
porof = poro(nz)  ! this assumes zero-porosity gradient at the depth; these choices do not affect the calculation 
sporof = 1d0-porof  !  volume fraction of solids at bottom depth  
! ------------------
! cf., poro = poro_0*exp(-z(iz)/poro_scale)  ! Hydrate modeling parameterization where poro_0 = 0.69 & poro_scale = 2000 (m)
sporoi = 1d0-poroi  ! volume fraction of solids at the seawater-sediment interface (SWI)
sporo = 1d0 - poro  !  volume fraction of solids 

endsubroutine getporosity
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine burial_pre(  &
    w,wi  & ! output
    ,detflx,ccflx,nspcc,nz,poroi,msed,mvsed,mvcc  & ! input 
    )
implicit none
integer(kind=4),intent(in)::nspcc,nz
real(kind=8),intent(in)::detflx,ccflx(nspcc),poroi,msed,mvsed,mvcc(nspcc)
real(kind=8),intent(out)::w(nz),wi

! burial rate w from rain fluxes represented by volumes
! initial guess assuming a box representation (this guess is accurate when there is no caco3 dissolution occurring) 
! om is not considered as it gets totally depleted 

wi = (detflx/msed*mvsed + sum(ccflx*mvcc)            )/(1d0-poroi)
w = wi

endsubroutine burial_pre
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine dep2age(  &
    age &  ! output 
    ,dz,w,nz  &  ! input
   )
implicit none
integer(kind=4),intent(in)::nz
real(kind=8),intent(in)::dz(nz),w(nz)
real(kind=8),intent(out)::age(nz)
real(kind=8)::dage(nz)
integer(kind=4) iz

dage = dz/w  ! time spans of individual sediment layers 
age = 0d0
do iz=1,nz  ! assigning ages to depth in the same way to assign depths to individual grids 
    if (iz==1) age(iz)=dage(iz)*0.5d0  
    if (iz/=1) age(iz) = age(iz-1)+dage(iz-1)*0.5d0 + 0.5d0*dage(iz)
enddo

endsubroutine dep2age
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine calcupwindscheme(  &
    up,dwn,cnr,adf & ! output 
    ,w,nz   & ! input &
    )
implicit none
integer(kind=4),intent(in)::nz
real(kind=8),intent(in)::w(nz)
real(kind=8),dimension(nz),intent(out)::up,dwn,cnr,adf
real(kind=8) corrf
integer(kind=4) iz


! ------------ determine variables to realize advection 
!  upwind scheme 
!  up  ---- burial advection at grid i = sporo(i)*w(i)*(some conc. at i) - sporo(i-1)*w(i-1)*(some conc. at i - 1) 
!  dwn ---- burial advection at grid i = sporo(i+1)*w(i+1)*(some conc. at i+1) - sporo(i)*w(i)*(some conc. at i) 
!  cnr ---- burial advection at grid i = sporo(i+1)*w(i+1)*(some conc. at i+1) - sporo(i-1)*w(i-1)*(some conc. at i - 1) 
!  when burial rate is positive, scheme need to choose up, i.e., up = 1.  
!  when burial rate is negative, scheme need to choose dwn, i.e., dwn = 1.  
!  where burial change from positive to negative or vice versa, scheme chooses cnr, i.e., cnr = 1. for the mass balance sake 

up = 0
dwn=0
cnr =0
adf=1d0
do iz=1,nz 
    if (iz==1) then 
        if (w(iz)>=0d0 .and. w(iz+1)>=0d0) then  ! positive burial 
            up(iz) = 1
        elseif (w(iz)<=0d0 .and. w(iz+1)<=0d0) then  ! negative burial 
            dwn(iz) = 1
        else   !  where burial sign changes  
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
    print*,'error',sum(up),sum(dwn),sum(cnr)
    stop
endif

do iz=1,nz-1
    if (cnr(iz)==1 .and. cnr(iz+1)==1) then 
        if (w(iz)>=0d0 .and. w(iz+1) < 0d0) then
            corrf = 5d0  !  This assignment of central advection term helps conversion especially when assuming turbo2 mixing 
            cnr(iz+1)=abs(w(iz)**corrf)/(abs(w(iz+1)**corrf)+abs(w(iz)**corrf))
            cnr(iz)=abs(w(iz+1)**corrf)/(abs(w(iz+1)**corrf)+abs(w(iz)**corrf))
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
            adf(iz)=abs(w(iz+1))/(abs(w(iz+1))+abs(w(iz)))
            adf(iz+1)=abs(w(iz))/(abs(w(iz+1))+abs(w(iz)))
        endif 
    endif 
enddo       

endsubroutine calcupwindscheme
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine recordtime(  &
    rectime,time_spn,time_trs,time_aft,cntrec  &
    ,ztot,wi,file_tmp,workdir,nrec  &
    )
implicit none 
integer(kind=4),intent(in)::nrec,file_tmp
real(kind=8),intent(in)::ztot,wi
real(kind=8),intent(out)::rectime(nrec),time_spn,time_trs,time_aft
integer(kind=4),intent(out)::cntrec
character*255,intent(in)::workdir
integer(kind=4) itrec
 
!!! +++ tracing experiment 
time_spn = ztot / wi *50d0 ! yr  ! spin-up duration, 50 times the shortest residence time possible (assuming no caco3 dissolution) 
time_trs = 50d3 !  duration of signal change event 
time_aft = time_trs*3d0  ! duration of simulation after the event 
#ifdef biotest
time_trs = 5d3   !  smaller duration of event assumed 
time_aft = time_trs*10d0  
#endif 
! distributing recording time in 3 different periods 
do itrec=1,nrec/3
    rectime(itrec)=itrec*time_spn/real(nrec/3)
enddo
do itrec=nrec/3+1,nrec/3*2
    rectime(itrec)=rectime(nrec/3)+(itrec-nrec/3)*time_trs/real(nrec/3)
enddo
do itrec=nrec/3*2+1,nrec
    rectime(itrec)=rectime(nrec/3*2)+(itrec-nrec/3*2)*time_aft/real(nrec/3)
enddo
#ifdef sense
time_trs = 0d0   !  there is no event needed 
time_aft = 0d0 
! time_spn = 10d0*time_spn  ! longer
do itrec=1,nrec
    rectime(itrec)=itrec*time_spn/real(nrec)
enddo
#endif 
!!! ++++
cntrec = 1  ! rec number (increasing with recording done )
#ifndef nonrec
open(unit=file_tmp,file=trim(adjustl(workdir))//'rectime.txt',action='write',status='unknown')
do itrec=1,nrec 
    write(file_tmp,*) rectime(itrec)  ! recording when records are made 
enddo
close(file_tmp)
#endif

endsubroutine recordtime
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine sig2sp_pre(  &  ! end-member signal assignment 
    d13c_sp,d18o_sp  &
    ,d13c_ocni,d13c_ocnf,d18o_ocni,d18o_ocnf,nspcc  &
    )
implicit none 
integer(kind=4),intent(in)::nspcc
real(kind=8),dimension(nspcc),intent(out)::d13c_sp,d18o_sp
real(kind=8),intent(in)::d13c_ocni,d13c_ocnf,d18o_ocni,d18o_ocnf 

#ifndef sense
!  four end-member caco3 species interpolation 
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
! assuming two differently sized caco3 species and giving additional 4 species their end-member isotopic compositions 
d13c_sp(5)=d13c_ocni
d18o_sp(5)=d18o_ocni
d13c_sp(6)=d13c_ocni
d18o_sp(6)=d18o_ocnf
d13c_sp(7)=d13c_ocnf
d18o_sp(7)=d18o_ocni
d13c_sp(8)=d13c_ocnf
d18o_sp(8)=d18o_ocnf
#endif 

endsubroutine sig2sp_pre
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine make_transmx(  &
    trans,izrec,izrec2,izml,nonlocal  & ! output 
    ,labs,nspcc,turbo2,nobio,dz,sporo,nz,z,file_tmp,zml_ref  & ! input
    )
implicit none
integer(kind=4),intent(in)::nspcc,nz,file_tmp
real(kind=8),intent(in)::dz(nz),sporo(nz),z(nz),zml_ref
logical,intent(in)::labs(nspcc+2),turbo2(nspcc+2),nobio(nspcc+2)
real(kind=8),intent(out)::trans(nz,nz,nspcc+2)
logical,intent(out)::nonlocal(nspcc+2)
integer(kind=4),intent(out)::izrec,izrec2,izml
integer(kind=4) nlabs,ilabs,iz,isp
real(kind=8) :: translabs(nz,nz),translabs_tmp(nz,nz),dbio(nz),transdbio(nz,nz),transturbo2(nz,nz)
real(kind=8) :: zml(nspcc+2),zrec,zrec2,probh
character*25 dumchr(3)

!~~~~~~~~~~~~ loading transition matrix from LABS ~~~~~~~~~~~~~~~~~~~~~~~~
if (any(labs)) then
    translabs = 0d0

    open(unit=file_tmp,file='./labs-mtx.txt',action='read',status='unknown')
    do iz=1,nz
        read(file_tmp,*) translabs(iz,:)  ! writing 
    enddo
    close(file_tmp)

endif

if (.true.) then  ! devided by the time duration when transition matrices are created in LABS and weakening by a factor
! if (.false.) then 
    translabs = translabs *365.25d0/10d0*1d0/3d0  
    ! translabs = translabs *365.25d0/10d0*1d0/10d0
endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
zml=zml_ref ! mixed layer depth assumed to be a reference value at first 

zrec = 1.1d0*maxval(zml)  ! depth where recording is made, aimed at slightly below the bottom of mixed layer  
! zrec = maxval(zml)  ! depth where recording is made, aimed at slightly below the bottom of mixed layer  
zrec2 = 2.0d0*maxval(zml)  ! depth where recording is made ver. 2, aimed at 2 time bottom depth of mixed layer 

#ifdef size 
zml(2+1)=20d0   ! fine species have larger mixed layers 
zml(2+2)=20d0   ! note that total number of solid species is 2 + nspcc including om, clay and nspcc of caco3; thus index has '2+'
zml(2+3)=20d0 
zml(2+4)=20d0 
zrec = 1.1d0*minval(zml)  ! first recording is made below minimum depth of mixed layer 
zrec2 = 1.1d0*maxval(zml) ! second recording is made below maximum depth of mixed layer
! zrec = minval(zml)  ! first recording is made below minimum depth of mixed layer 
! zrec2 = maxval(zml) ! second recording is made below maximum depth of mixed layer
#endif 

do iz=1,nz ! determine grid locations where signal recording is made 
    if (z(iz)<=zrec) izrec = iz  
    if (z(iz)<=zrec2) izrec2 = iz
enddo
! izrec = min(nz,izrec+1)
! izrec2 = min(nz,izrec2+1)

nonlocal = .false. ! initial assumption 
do isp=1,nspcc+2
    if (turbo2(isp) .or. labs(isp)) nonlocal(isp)=.true. ! if mixing is made by turbo2 or labs, then nonlocal 
    
    dbio=0d0
    do iz = 1, nz
        if (z(iz) <=zml(isp)) then
            dbio(iz) =  0.15d0   !  within mixed layer 150 cm2/kyr (Emerson, 1985) 
            izml = iz   ! determine grid of bottom of mixed layer 
        else
            dbio(iz) =  0d0 ! no biodiffusion in deeper depths 
        endif
    enddo

    transdbio = 0d0   ! transition matrix to realize Fickian mixing with biodiffusion coefficient dbio which is defined just above 
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
    probh = 0.0010d0
    transturbo2(:izml,:izml) = probh  ! arbitrary assumed probability 
    do iz=1,izml  ! when i = j, transition matrix contains probabilities with which particles are moved from other layers of sediment   
       transturbo2(iz,iz)=-probh*(izml-1)  
    enddo

    if (turbo2(isp)) translabs = transturbo2   ! translabs temporarily used to represents nonlocal mixing 

    trans(:,:,isp) = transdbio(:,:)  !  firstly assume local mixing implemented by dbio 

    if (nonlocal(isp)) trans(:,:,isp) = translabs(:,:)  ! if nonlocal, replaced by either turbo2 mixing or labs mixing 
    if (nobio(isp)) trans(:,:,isp) = 0d0  ! if assuming no bioturbation, transition matrix is set at zero  
enddo
! even when all are local Fickian mixing, mixing treatment must be the same as in case of nonlocal 
! if mixing intensity and depths are different between different species  
if (all(.not.nonlocal)) then  
    do isp=1,nspcc+2-1
        if (any(trans(:,:,isp+1)/=trans(:,:,isp))) nonlocal=.true.
    enddo
endif 

endsubroutine make_transmx
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine coefs(  &
    dif_dic,dif_alk,dif_o2,kom,kcc,co3sat,krad & ! output 
    ,temp,sal,dep,nz,nspcc,poro,cai,komi,kcci,k14ci,i14c  & !  input
    ,i13c18o  &
    )
integer(kind=4),intent(in)::nz,nspcc
real(kind=8),intent(in)::temp,sal,dep,poro(nz),cai,komi,kcci,k14ci
real(kind=8),intent(out)::dif_dic(nz),dif_alk(nz),dif_o2(nz),kom(nz),kcc(nz,nspcc),co3sat,krad(nz,nspcc)
real(kind=8) dif_dic0,dif_alk0,dif_o20,ff(nz),keq1,keq2,keqcc
real(kind=8) calceqcc,calceq1,calceq2

ff = poro*poro       ! representing tortuosity factor 

dif_dic0 = (151.69d0 + 7.93d0*temp) ! cm2/yr at 2 oC (Huelse et al. 2018)
dif_alk0 = (151.69d0 + 7.93d0*temp) ! cm2/yr  (Huelse et al. 2018)
dif_o20 =  (348.62d0 + 14.09d0*temp) 

dif_dic = dif_dic0*ff  ! reflecting tortuosity factor 
dif_alk = dif_alk0*ff
dif_o2 = dif_o20*ff

kom = komi  ! assume reference values for all reaction terms 
kcc = kcci
krad = 0d0
#ifdef isotrack
krad(:,i14c) = k14ci 
#ifdef kie 
kcc(:,i13c18o) = kcci*(1d0-5d-5) 
#endif 
#endif 

#ifdef size 
! assume stronger dissolution for fine species (1-4) 
kcc(:,1) = kcci*10d0
kcc(:,2) = kcci*10d0
kcc(:,3) = kcci*10d0
kcc(:,4) = kcci*10d0
#endif 

keq1 = calceq1(temp,sal,dep) ! carbonic acid dissociation const. function called from caco3_therm.f90 
keq2 = calceq2(temp,sal,dep) ! bicarbonate dissociation const. function called from caco3_therm.f90

keqcc = calceqcc(temp,sal,dep) ! calcite solubility function called from caco3_therm.f90
co3sat = keqcc/cai ! co3 conc. at calcite saturation 

! print*,cai,keqcc,co3sat

endsubroutine coefs
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine recordprofile(  &
    itrec,file_tmp,workdir,nz,z,age,pt,rho,cc,ccx,dic,dicx,alk,alkx,co3,co3x,nspcc,msed,wi,co3sat,rcc  &
    ,pro,o2x,oxco2,anco2,om,mom,mcc,d13c_ocni,d18o_ocni,up,dwn,cnr,adf,ptx,w,frt,prox,omx,d13c_blk,d18o_blk  &
    ,d17o_blk,d14c_age,capd47  &
    )
implicit none 
integer(kind=4),intent(in):: itrec,file_tmp,nz,nspcc
real(kind=8),dimension(nz),intent(in)::z,age,pt,rho,dic,dicx,alk,alkx,co3,co3x,pro,o2x,oxco2,anco2,om,up,dwn,cnr,adf
real(kind=8),dimension(nz),intent(in)::ptx,w,frt,prox,omx,d13c_blk,d18o_blk,d17o_blk,d14c_age,capd47
real(kind=8),dimension(nz,nspcc),intent(in)::cc,ccx,rcc
real(kind=8),intent(in)::msed,wi,co3sat,mom,mcc(nspcc),d13c_ocni,d18o_ocni
character*255,intent(in)::workdir
character*25 dumchr(3)
integer(kind=4) iz,isp

write(dumchr(1),'(i3.3)') itrec  

if (itrec==0) then 
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
        write(file_tmp,*) z(iz),age(iz),(cc(iz,isp)*mcc(isp)/2.5d0*100d0,isp=1,nspcc) 
    enddo
    close(file_tmp)

    open(unit=file_tmp,file=trim(adjustl(workdir))//'sig-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),d13c_ocni,d18o_ocni,0d0,0d0,0d0
    enddo
    close(file_tmp)

    open(unit=file_tmp,file=trim(adjustl(workdir))//'bur-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),w(iz),up(iz),dwn(iz),cnr(iz),adf(iz)
    enddo
    close(file_tmp)
else 
    open(unit=file_tmp,file=trim(adjustl(workdir))//'ptx-'//trim(adjustl(dumchr(1)))//'.txt' &
        ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),ptx(iz)*msed/rho(iz)*100d0,rho(iz),frt(iz)  ,w(iz)
    enddo
    close(file_tmp)

    open(unit=file_tmp,file=trim(adjustl(workdir))//'ccx-'//trim(adjustl(dumchr(1)))//'.txt' &
        ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),sum(ccx(iz,:)*mcc(:))/rho(iz)*100d0, dicx(iz)*1d3, alkx(iz)*1d3  &
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
        write(file_tmp,*) z(iz),age(iz),(ccx(iz,isp)*mcc(isp)/rho(iz)*100d0,isp=1,nspcc) 
    enddo
    close(file_tmp)

    open(unit=file_tmp,file=trim(adjustl(workdir))//'sig-'  &
        //trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),d13c_blk(iz),d18o_blk(iz),d17o_blk(iz),d14c_age(iz),capd47(iz)
    enddo
    close(file_tmp)

    open(unit=file_tmp,file=trim(adjustl(workdir))//'bur-'//trim(adjustl(dumchr(1)))//'.txt' ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),age(iz),w(iz),up(iz),dwn(iz),cnr(iz),adf(iz)
    enddo
    close(file_tmp)
endif 

endsubroutine recordprofile
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine signal_flx(  &
    d13c_ocn,d18o_ocn,ccflx,d18o_sp,d13c_sp,cntsp  &
    ,time,time_spn,time_trs,time_aft,d13c_ocni,d13c_ocnf,d18o_ocni,d18o_ocnf,nspcc,ccflxi,it,flxfini,flxfinf  &
    ,r14ci,capd47_ocni,capd47_ocnf,capd47_ocn,r13c_pdb,r18o_pdb,r17o_pdb,tol,nt_trs   &
    ) 
implicit none 
integer(kind=4),intent(in)::nspcc,it,nt_trs
integer(kind=4),intent(inout)::cntsp
real(kind=8),intent(in)::time,time_spn,time_trs,time_aft,d13c_ocni,d13c_ocnf,d18o_ocni,d18o_ocnf,ccflxi
real(kind=8),intent(in)::flxfini,flxfinf
real(kind=8),intent(out)::d13c_ocn,d18o_ocn,ccflx(nspcc)
real(kind=8),dimension(nspcc),intent(inout)::d18o_sp,d13c_sp
real(kind=8) flxfin
real(kind=8),dimension(nspcc)::flxfrc,flxfrc2 
integer(kind=4) isp
! used only when isotrack is ON
real(kind=8),intent(in)::capd47_ocni,capd47_ocnf,r14ci,r13c_pdb,r18o_pdb,r17o_pdb,tol
real(kind=8),intent(out)::capd47_ocn
real(kind=8) r13c_ocn,r12c_ocn,r14c_ocn,f13c_ocn,f12c_ocn,f14c_ocn
real(kind=8) r18o_ocn,r16o_ocn,r17o_ocn,f18o_ocn,f16o_ocn,f17o_ocn
integer(kind=4),allocatable::ipiv(:)
real(kind=8),allocatable :: amx(:,:),ymx(:),emx(:)
integer(kind=4) infobls,nmx
real(kind=8) d2r,r2d

if (time <= time_spn) then   ! spin-up
    d13c_ocn = d13c_ocni  ! take initial values 
    d18o_ocn = d18o_ocni
    capd47_ocn=capd47_ocni
    ccflx = 0d0  
    ccflx(1) = ccflxi  ! raining only species with initial signal values 
#ifdef track2
    cntsp=1  ! when signal is tracked by method 2, signal is assigned species as time goes so counting the number of species whose signal has been assigned   
    d18o_sp(cntsp)=d18o_ocn
    d13c_sp(cntsp)=d13c_ocn
    ccflx = 0d0
    ccflx(cntsp) = ccflxi
#endif 
#ifdef size 
    ccflx = 0d0
    flxfin = flxfini
    ccflx(1) = (1d0-flxfin)*ccflxi ! fine species with initial signals 
    ccflx(5) = flxfin*ccflxi       ! coarse species with initial signals 
#endif 
elseif (time>time_spn .and. time<=time_spn+time_trs) then ! during event 
    d13c_ocn = d13c_ocni + (time-time_spn)*(d13c_ocnf-d13c_ocni)/time_trs ! single step
    d18o_ocn = d18o_ocni + (time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs  ! single step 
    capd47_ocn = capd47_ocni + (time-time_spn)*(capd47_ocnf-capd47_ocni)/time_trs  ! single step 
    ! creating spike (go from initial to final and come back to initial again taking half time of event duration )
    ! this shape is assumed for d18o and fine (& coarse) caco3 flux changes 
    if (time-time_spn<=time_trs/2d0) then
        d18o_ocn = d18o_ocni + (time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs*2d0
        capd47_ocn = capd47_ocni + (time-time_spn)*(capd47_ocnf-capd47_ocni)/time_trs*2d0      
        flxfin = flxfini + (time-time_spn)*(flxfinf-flxfini)/time_trs*2d0           
    else 
        d18o_ocn = 2d0*d18o_ocnf - d18o_ocni - (time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs*2d0
        capd47_ocn = 2d0*capd47_ocnf - capd47_ocni - (time-time_spn)*(capd47_ocnf-capd47_ocni)/time_trs*2d0
        flxfin = 2d0*flxfinf - flxfini - (time-time_spn)*(flxfinf-flxfini)/time_trs*2d0
    endif
#ifndef biotest    
    ! assuming a d13 excursion with shifts occurring with 1/10 times event duration
    ! the same assumption for depth change 
    if (time-time_spn<=time_trs/10d0) then
        d13c_ocn = d13c_ocni + (time-time_spn)*(d13c_ocnf-d13c_ocni)/time_trs*10d0
    elseif (time-time_spn>time_trs/10d0 .and. time-time_spn<=time_trs/10d0*9d0) then
        d13c_ocn = d13c_ocnf
    elseif  (time-time_spn>time_trs/10d0*9d0) then 
        d13c_ocn = 10d0*d13c_ocnf - 9d0*d13c_ocni  - (time-time_spn)*(d13c_ocnf-d13c_ocni)/time_trs*10d0
    endif
#endif    
    if (.not.(d13c_ocn>=d13c_ocnf .and.d13c_ocn<=d13c_ocni)) then ! check if calculated d13c and d18o are within the assumed ranges  
        print*,'error in d13c',d13c_ocn
        stop
    endif
    ! calculating fractions of flux assigned to individual caco3 species in case of interpolation of 2 signal inputs by 4 species 
    ! NEED to extend to allow tacking of any number of signals 
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
        print*,'error' ! should not come here 
        stop
    enddo 
#ifdef track2    
    ! tracking as time goes 
    ! assignment of new caco3 species is conducted so that nspcc species is used up within nt_trs interations
    ! timing of new species assignment can change at different time period during the event 
    if (time-time_spn<=time_trs/10d0) then
        if (mod(it,int(nt_trs/(nspcc-2)))==0) then   
            cntsp=cntsp+1 ! new species assinged 
            d18o_sp(cntsp)=d18o_ocn
            d13c_sp(cntsp)=d13c_ocn
            ccflx = 0d0
            ccflx(cntsp) = ccflxi
        endif 
    elseif (time-time_spn>time_trs/10d0 .and. time-time_spn<=time_trs/10d0*9d0) then
        if (mod(it,int(nt_trs/(nspcc-2)))==0) then  
            cntsp=cntsp+1
            d18o_sp(cntsp)=d18o_ocn
            d13c_sp(cntsp)=d13c_ocn
            ccflx = 0d0
            ccflx(cntsp) = ccflxi
        endif 
    elseif  (time-time_spn>time_trs/10d0*9d0) then 
        if (mod(it,int(nt_trs/(nspcc-2)))==0) then  
            cntsp=cntsp+1
            d18o_sp(cntsp)=d18o_ocn
            d13c_sp(cntsp)=d13c_ocn
            ccflx = 0d0
            ccflx(cntsp) = ccflxi
        endif 
    endif 
#else 
    ! usually using 4 species interpolation for 2 signals 
    ! NEED to be extended so that any number of signals can be tracked with corresponding minimum numbers of caco3 species 
    do isp=1,nspcc
        ccflx(isp)=flxfrc2(isp)*ccflxi
    enddo
#endif 
    
#ifdef size
    do isp=1,4 ! fine species 
        ccflx(isp)=flxfrc2(isp)*ccflxi*(1d0-flxfin)
    enddo
    do isp=5,8 ! coarse species 
        ccflx(isp)=flxfrc2(isp-4)*ccflxi*flxfin
    enddo
#endif 

elseif (time>time_spn+time_trs) then 
    d13c_ocn = d13c_ocni ! now again initial values 
    d18o_ocn = d18o_ocni
    capd47_ocn=capd47_ocni
    ccflx = 0d0
    ccflx(1) = ccflxi  ! allowing only caco3 flux with initial signal values 
#ifdef track2
    if (cntsp+1/=nspcc) then  ! checking used caco3 species number is enough and necessary 
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
    ccflx(1) = (1d0-flxfin)*ccflxi  ! fine species 
    ccflx(5) = flxfin*ccflxi  ! coarse species 
#endif 
#ifdef biotest 
    d13c_ocn = d13c_ocnf   ! finish with final value 
    d18o_ocn = d18o_ocni
    capd47_ocn=capd47_ocni
    ccflx = 0d0
    ccflx(3) = ccflxi  ! caco3 species with d13c_ocnf and d18o_ocni
#endif 
endif

! only when directly tracking isotopes 
#ifdef isotrack
r13c_ocn = d2r(d13c_ocn,r13c_pdb) 
r12c_ocn = 1d0
r14c_ocn = r14ci 
f12c_ocn = r12c_ocn/(r12c_ocn+r13c_ocn+r14c_ocn)
f13c_ocn = r13c_ocn/(r12c_ocn+r13c_ocn+r14c_ocn)
f14c_ocn = r14c_ocn/(r12c_ocn+r13c_ocn+r14c_ocn)
r18o_ocn = d2r(d18o_ocn,r18o_pdb)
r16o_ocn = 1d0
r17o_ocn = ((17d0-16d0)/(18d0-16d0)*18d0*16d0/(17d0*16d0)*(r18o_ocn/r18o_pdb-1d0)+1d0)*r17o_pdb
#ifndef fullclump
r17o_ocn = 0d0  ! when not considering 17o
#endif 
f16o_ocn = r16o_ocn/(r16o_ocn+r17o_ocn+r18o_ocn)
f17o_ocn = r17o_ocn/(r16o_ocn+r17o_ocn+r18o_ocn)
f18o_ocn = r18o_ocn/(r16o_ocn+r17o_ocn+r18o_ocn)

! 4 species 

nmx = 5 !  
if (allocated(amx)) deallocate(amx)
if (allocated(ymx)) deallocate(ymx)
if (allocated(emx)) deallocate(emx)
if (allocated(ipiv)) deallocate(ipiv)
allocate(amx(nmx,nmx),ymx(nmx),emx(nmx),ipiv(nmx))

amx = 0d0
ymx = 0d0

amx(1,1)=1d0
amx(1,2)=1d0
ymx(1)=f12c_ocn*ccflxi 

amx(2,3)=1d0
amx(2,4)=1d0
ymx(2)=f13c_ocn*ccflxi 

amx(3,1)=3d0
amx(3,2)=2d0
amx(3,3)=3d0
amx(3,4)=2d0
amx(3,5)=3d0*f16o_ocn
ymx(3)=3d0*f16o_ocn*ccflxi

amx(4,4)=1d0
amx(4,1)=-1d0*r13c_ocn*r18o_ocn*(capd47_ocn*1d-3+1d0)

amx(5,5)=1d0
ymx(5)=f14c_ocn*ccflxi

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls) 

! print*,ymx,sum(ymx),infobls

flxfrc2 = 0d0
flxfrc2 = ymx/ccflxi
#ifndef fullclump
if (abs(sum(flxfrc2)-1d0)>tol) then 
    print*,'error in flx',flxfrc2
    stop
endif 
if (any(flxfrc2<0d0)) then 
    print*,'negative flx',flxfrc2
    stop
endif 
#endif
flxfrc = 0d0
flxfrc(:) = flxfrc2(:)!/sum(flxfrc2)

ccflx = 0d0
ccflx = flxfrc*ccflxi
#endif 

endsubroutine signal_flx
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine bdcnd(   &
    time,dep,time_spn,time_trs,time_aft,depi,depf  &
    ) 
implicit none 
real(kind=8),intent(in):: time,time_spn,time_trs,time_aft,depi,depf
real(kind=8),intent(out):: dep

if (time <= time_spn) then   ! spin-up
    dep = depi  ! initial depth 
elseif (time>time_spn .and. time<=time_spn+time_trs) then ! during event 
#ifndef biotest    
    ! assuming a d13 excursion with shifts occurring with 1/10 times event duration
    ! the same assumption for depth change 
    if (time-time_spn<=time_trs/10d0) then
        dep = depi + (depf-depi)*(time-time_spn)/time_trs*10d0
    elseif (time-time_spn>time_trs/10d0 .and. time-time_spn<=time_trs/10d0*9d0) then
        dep = depf
    elseif  (time-time_spn>time_trs/10d0*9d0) then 
        dep = 10d0*depf-9d0*depi - (depf-depi)*(time-time_spn)/time_trs*10d0
    endif
#else     
    dep = depi
#endif
elseif (time>time_spn+time_trs) then 
    dep = depi
endif

endsubroutine bdcnd
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine timestep(time,nt_spn,nt_trs,nt_aft,dt,time_spn,time_trs,time_aft)
implicit none 
real(kind=8),intent(in)::time,time_spn,time_trs,time_aft
integer(kind=4),intent(in)::nt_spn,nt_trs,nt_aft
real(kind=8),intent(out)::dt
  
if (time <= time_spn) then   ! spin-up
    if (time+dt> time_spn) then
        dt=time_trs/real(nt_trs,kind=8) !5000d0   ! when close to 'event', time step needs to get smaller   
    else
        dt = time_spn/real(nt_spn,kind=8)! 800d0 ! otherwise larger time step is better to fasten calculation 
    endif
elseif (time>time_spn .and. time<=time_spn+time_trs) then ! during event 
    dt = time_trs/real(nt_trs,kind=8) !5000d0
elseif (time>time_spn+time_trs) then 
    dt=time_trs/real(nt_aft,kind=8) !1000d0 ! not too large time step
endif

endsubroutine timestep
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine calc_zox( &
    izox,kom,zox,kom_ox,kom_an   &  ! output 
    ,oxic,anoxic,nz,o2x,o2th,komi,ztot,z,o2i,dz  & ! input
    )
implicit none
integer(kind=4),intent(in)::nz
real(kind=8),intent(in)::o2x(nz),o2th,komi,ztot,z(nz),dz(nz),o2i
logical,intent(in)::oxic,anoxic
integer(kind=4),intent(out)::izox
real(kind=8),intent(out)::kom(nz),zox,kom_ox(nz),kom_an(nz)
logical izox_calc_done 
integer(kind=4) iz
real(kind=8)::tol=1d-6

! izox_calc_done = .false.
! if (oxic) then 
    ! do iz=1,nz
        ! if (o2x(iz) > o2th) then
            ! kom(iz) = komi
            ! if (.not. izox_calc_done) izox = iz
        ! else! unless anoxi degradation is allowed, om cannot degradate below zox
            ! kom(iz) = 0d0
            ! if (anoxic) kom(iz) = komi
            ! izox_calc_done = .true.
        ! endif
    ! enddo
! endif

!! calculation of zox 
zox = 0d0
do iz=1,nz
    if (o2x(iz)<=0d0) exit
enddo

if (iz==nz+1) then ! oxygen never gets less than 0 
    zox = ztot ! zox is the bottom depth 
else if (iz==1) then ! calculating zox interpolating at z=0 with SWI conc. and at z=z(iz) with conc. o2x(iz)
    zox = (z(iz)*o2i*1d-6/1d3 + 0d0*abs(o2x(iz)))/(o2i*1d-6/1d3+abs(o2x(iz)))
else if (iz==2) then 
    zox = z(iz-1) - o2x(iz-1)/((o2i*1d-6/1d3 - o2x(iz-1))/(0d0-z(iz-1)))
else     ! calculating zox interpolating at z=z(iz-1) with o2x(iz-1) and at z=z(iz) with conc. o2x(iz)
    ! zox = (z(iz)*o2x(iz-1) + z(iz-1)*abs(o2x(iz)))/(o2x(iz-1)+abs(o2x(iz)))
    zox = z(iz-1) - o2x(iz-1)/((o2x(iz-2) - o2x(iz-1))/(z(iz-2)-z(iz-1)))
endif

! calculation of kom 
kom = 0d0
kom_ox = 0d0
kom_an = 0d0
izox = 0 
if (anoxic) then 
    kom = komi
    do iz=1,nz
        if (z(iz)+0.5d0*dz(iz)<=zox) then 
            kom_ox(iz)=komi
            if (iz> izox ) izox = iz
        elseif (z(iz)+0.5d0*dz(iz)>zox .and. z(iz)-0.5d0*dz(iz)< zox) then 
            kom_ox(iz)=komi* (1d0- ( (z(iz)+0.5d0*dz(iz)) - zox)/dz(iz))
            kom_an(iz)=komi* (( (z(iz)+0.5d0*dz(iz)) - zox)/dz(iz))
            if (iz> izox ) izox = iz
        elseif (z(iz)-0.5d0*dz(iz)>=zox) then 
            kom_an(iz)=komi
        endif 
    enddo 
    if (.not.all(abs(kom_ox+kom_an-kom)/komi<tol)) then 
        print*,'error: calc kom',kom_ox,kom_an,kom
        stop
    endif 
else
    do iz=1,nz
        if (z(iz)+0.5d0*dz(iz)<=zox) then 
            kom_ox(iz)=komi
            if (iz> izox ) izox = iz
        elseif (z(iz)+0.5d0*dz(iz)>zox .and. z(iz)-0.5d0*dz(iz)< zox) then 
            kom_ox(iz)=komi* (1d0- ( (z(iz)+0.5d0*dz(iz)) - zox)/dz(iz))
            if (iz> izox ) izox = iz
        elseif (z(iz)-0.5d0*dz(iz)>=zox) then 
            continue
        endif 
    enddo 
    kom = kom_ox
endif 

endsubroutine calc_zox
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine omcalc( &
    omx  & ! output 
    ,kom   &  ! input
    ,om,nz,sporo,sporoi,sporof &! input 
    ,w,wi,dt,up,dwn,cnr,adf,trans,nspcc,labs,turbo2,nonlocal,omflx,poro,dz &! input 
    ) 
implicit none 
integer(kind=4),intent(in)::nz,nspcc
real(kind=8),dimension(nz),intent(in)::om,sporo,w,up,dwn,cnr,adf,poro,dz
real(kind=8),intent(in)::sporoi,sporof,wi,dt,trans(nz,nz,nspcc+2),omflx
logical,dimension(nspcc+2),intent(in)::labs,turbo2,nonlocal
real(kind=8),intent(out)::omx(nz)
real(kind=8),intent(in)::kom(nz)
integer(kind=4) :: nsp,nmx,iz,row,col,iiz,infobls
integer(kind=4),allocatable::ipiv(:)
real(kind=8),allocatable :: amx(:,:),ymx(:),emx(:)

!  amx and ymx correspond to A and B in Ax = B
!  dgesv subroutine of BLAS returns the solution x in ymx 
!  emx stores errors (not required for dgesv); ipiv is required for dgesv
!  E.g., difference form of governing equation at grid 2 can be expressed as 
!
!               sporo(2)*(omx(2)-om(2))/dt + (sporo(2)*w(2)*omx(2)-sporo(1)*w(1)*omx(1))/dz(2) + sporo(2)*kom(2)*omx(2) = 0 
!
!  In above difference equation, w(2) and w(1) are assumed to be positive for illustration purpose (and dropping adf, up, dwn and cnr terms)
!  and omx(2) and omx(1) are unknowns to be solved. 
!  x contains omx(1), omx(2), ...., omx(nz), i.e., nz unknowns  
!  and thus the above equation fills the matrix A as  
!
!               A(2,1) =  (-sporo(1)*w(1)*1)/dz(2)
!               A(2,2) =  sporo(2)*(1)/dt + (sporo(2)*w(2)*1)/dz(2) + sporo(2)*kom(2)*1
!
!  and the matrix B as 
!
!               - B(2) = sporo(2)*(-om(2))/dt
!
!  Matrices A and B are filled in this way. Note again amx and ymx correspond A and B, respectively. 
    
nsp=1 ! number of species considered here; 1, only om 
nmx = nz*nsp ! # of col (& row) of matrix A to in linear equations Ax = B to be solved, each species has nz (# of grids) unknowns 
if (allocated(amx)) deallocate(amx)
if (allocated(ymx)) deallocate(ymx)
if (allocated(emx)) deallocate(emx)
if (allocated(ipiv)) deallocate(ipiv)
allocate(amx(nmx,nmx),ymx(nmx),emx(nmx),ipiv(nmx))

amx = 0d0
ymx = 0d0

do iz = 1,nz 
    row = 1 + (iz-1)*nsp ! row number is obtained from grid number; here simply gird 1 corresponds to row 1 
    if (iz == 1) then ! need to reflect upper boundary, rain flux; and be careful that iz - 1 does not exit  
        ymx(row) = &
            ! time change term 
            + sporo(iz)*(-om(iz))/dt &
            ! rain flux term 
            - omflx/dz(1)
        amx(row,row) = (&
            ! time change term 
            + sporo(iz)*(1d0)/dt &
            ! advection terms 
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-sporoi*wi*0d0)/dz(1)   &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz)*w(iz)*1d0)/dz(1)   &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporoi*wi*0d0)/dz(1)   &
            !  rxn term 
            + sporo(iz)*kom(iz)   &
            ) 
        ! matrix filling at grid iz but for unknwon at grid iz + 1 (here only advection terms) 
        amx(row,row+nsp) =  (&
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz)*w(iz)*0d0)/dz(1)   &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporoi*wi*0d0)/dz(1)   &
            )
    else if (iz == nz) then ! need to reflect lower boundary; none; but must care that iz + 1 does not exist 
        ymx(row) = 0d0   &
            ! time change term 
            + sporo(iz)*(-om(iz))/dt 
        amx(row,row) = (&
            ! time change term 
            + sporo(iz)*(1d0)/dt &
            ! advection terms 
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(sporof*w(iz)*1d0-sporo(iz)*w(iz)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporof*w(iz)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  &
            ! rxn term 
            + sporo(iz)*kom(iz)   &
            )
        ! filling matrix at grid iz but for unknown at grid iz-1 (only advection terms) 
        amx(row,row-nsp) = ( &
            + adf(iz)*up(iz)*(sporof*w(iz)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporof*w(iz)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            )
    else ! do not have to care about boundaries 
        ymx(row) = 0d0  &
            ! time change term 
            + sporo(iz)*(0d0-om(iz))/dt 
        amx(row,row) = (&
            ! time change term 
            + sporo(Iz)*(1d0)/dt &
            ! advection terms 
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz)*w(iz)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  &
            ! rxn term 
            + sporo(iz)*kom(iz)   &
            )
        ! filling matrix at grid iz but for unknown at grid iz+1 (only advection terms) 
        amx(row,row+nsp) =  (&
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz)*w(iz)*0d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  &
            )
        ! filling matrix at grid iz but for unknown at grid iz-1 (only advection terms) 
        amx(row,row-nsp) =  (&
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            )
    endif
    ! diffusion terms are reflected with transition matrices 
    if (turbo2(1).or.labs(1)) then 
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp 
            if (trans(iiz,iz,1)==0d0) cycle
            amx(row,col) = amx(row,col) -trans(iiz,iz,1)/dz(iz)*dz(iiz)*(1d0-poro(iiz))
        enddo
    else 
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp 
            if (trans(iiz,iz,1)==0d0) cycle
            amx(row,col) = amx(row,col) -trans(iiz,iz,1)/dz(iz)
        enddo
    endif
enddo

ymx = - ymx  ! I have filled matrix B without changing signs; here I change signs at once 

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls) 

omx = ymx ! now passing the solution to unknowns omx 

endsubroutine omcalc
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine calcflxom(  &
    omadv,omdec,omdif,omrain,omres,omtflx  & ! output 
    ,sporo,om,omx,dt,w,dz,z,nz,turbo2,labs,nonlocal,poro,up,dwn,cnr,adf,rho,mom,trans,kom,sporof,sporoi,wi,nspcc,omflx  & ! input 
    ,file_tmp,workdir &
    ,flg_500  &
    )
implicit none 
integer(kind=4),intent(in)::nz,nspcc,file_tmp
real(kind=8),dimension(nz),intent(in)::sporo,om,omx,poro,up,dwn,cnr,adf,rho,kom,w,dz,z
real(kind=8),intent(in)::dt,mom,trans(nz,nz,nspcc+2),sporof,sporoi,wi
real(kind=8),intent(out)::omadv,omdec,omdif,omrain,omflx,omres,omtflx
logical,dimension(nspcc+2),intent(in)::turbo2,labs,nonlocal
character*255,intent(in)::workdir
logical,intent(out)::flg_500
integer(kind=4) :: iz,row,iiz,col,isp,nsp=1

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
        omdec = omdec + sporo(iz)*kom(iz)*omx(iz)*dz(iz)
        omrain = omrain - omflx/dz(1)*dz(iz)
    else if (iz == nz) then 
        omtflx = omtflx + sporo(iz)*(omx(iz)-om(iz))/dt*dz(iz)
        omadv = omadv + adf(iz)*up(iz)*(sporo(iz)*w(iz)*omx(iz)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz) &
            + adf(iz)*dwn(iz)*(sporof*w(iz)*omx(iz)-sporo(iz)*w(iz)*omx(iz))/dz(iz)*dz(iz) &
            + adf(iz)*cnr(iz)*(sporof*w(iz)*omx(iz)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz) 
        omdec = omdec + sporo(iz)*kom(iz)*omx(iz)*dz(iz)
    else 
        omtflx = omtflx + sporo(iz)*(omx(iz)-om(iz))/dt*dz(iz)
        omadv = omadv + adf(iz)*up(iz)*(sporo(iz)*w(iz)*omx(iz)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*omx(iz+1)-sporo(iz)*w(iz)*omx(iz))/dz(iz)*dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*omx(iz+1)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz) 
        omdec = omdec + sporo(iz)*kom(iz)*omx(iz)*dz(iz)
    endif
    if (turbo2(1).or.labs(1)) then 
        do iiz = 1, nz
            if (trans(iiz,iz,1)==0d0) cycle
            omdif = omdif -trans(iiz,iz,1)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*dz(iz)*omx(iiz)
        enddo
    else
        do iiz = 1, nz
            if (trans(iiz,iz,1)==0d0) cycle
            omdif = omdif -trans(iiz,iz,1)/dz(iz)*dz(iz)*omx(iiz)  ! check previous versions 
        enddo
    endif
enddo

omres = omadv + omdec + omdif + omrain + omtflx ! this is residual flux should be zero equations are exactly satisfied 

flg_500 = .false.
if (any(omx<0d0)) then  ! if negative om conc. is detected, need to stop  
    print*,'negative om, stop'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'NEGATIVE_OM.txt',status = 'unknown')
    do iz = 1, nz
        write (file_tmp,*) z(iz),omx(iz)*mom/rho(iz)*100d0,w(iz),up(iz),dwn(iz),cnr(iz),adf(iz)
    enddo
    close(file_tmp)
    ! stop
    flg_500 = .true.
endif 
if (any(isnan(omx))) then  ! if NAN, ... the same ... stop
    print*,'nan om, stop'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'NAN_OM.txt',status = 'unknown')
    do iz = 1, nz
        write (file_tmp,*) z(iz),omx(iz),kom(iz),rho(iz),w(iz),up(iz),dwn(iz),cnr(iz),adf(iz)
    enddo
    close(file_tmp)
    ! print*,omx
    ! stop
    flg_500 = .true.
endif 

endsubroutine calcflxom
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine o2calc_ox(  &
    o2x  & ! output
    ,nz,poro,o2,kom,omx,sporo,dif_o2,dz,dt,ox2om,o2i & ! input
    )
implicit none
integer(kind=4),intent(in)::nz
real(kind=8),dimension(nz),intent(in)::poro,o2,kom,omx,sporo,dif_o2,dz
real(kind=8),intent(in)::dt,ox2om,o2i
real(kind=8),intent(out)::o2x(nz)
integer(kind=4) row,nmx,nsp,iz,infobls
real(kind=8),allocatable :: amx(:,:),ymx(:),emx(:)
integer(kind=4),allocatable::ipiv(:)
    
nsp=1 ! number of species considered here; 1, only om 
nmx = nz*nsp ! # of col (& row) of matrix A to in linear equations Ax = B to be solved, each species has nz (# of grids) unknowns 
if (allocated(amx)) deallocate(amx)
if (allocated(ymx)) deallocate(ymx)
if (allocated(emx)) deallocate(emx)
if (allocated(ipiv)) deallocate(ipiv)
allocate(amx(nmx,nmx),ymx(nmx),emx(nmx),ipiv(nmx))

! reset matrices 
amx = 0d0
ymx = 0d0

do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then ! be careful about upper boundary 
        ymx(row) = ( &
            ! time change term 
            + poro(iz)*(0d0-o2(iz))/dt & 
            ! diffusion term 
            - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(0d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_o2(iz)*(0d0-o2i*1d-6/1d3)/dz(iz))/dz(iz)  &
            ! rxn term 
            + sporo(iz)*ox2om*kom(iz)*omx(iz)  &
            )
        amx(row,row) = (& 
            ! time change term 
            + poro(iz)*(1d0)/dt &
            ! diffusion term 
            - ((poro(iz)*dif_o2(iz)+poro(Iz+1)*dif_o2(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_o2(iz)*(1d0)/dz(iz))/dz(iz)&
            )
        ! filling matrix at grid iz but for unknown at grid iz+1 (only diffusion term) 
        amx(row,row+nsp) = (& 
            - ((poro(Iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - 0d0)/dz(iz)&
            )
    else if (iz == nz) then ! be careful about lower boundary 
        ymx(row) = (0d0 & 
            ! time change term 
            + poro(iz)*(0d0-o2(iz))/dt &
            ! diffusion term 
            - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(Iz-1))*(0d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            ! rxn term 
            + sporo(iz)*ox2om*kom(iz)*omx(iz)  &
            )
        amx(row,row) = ( & 
            ! time change term 
            + poro(iz)*(1d0)/dt &
            ! diffusion term 
            - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(Iz-1)*dif_o2(Iz-1))*(1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            )
        ! filling matrix at grid iz but for unknown at grid iz-1 (only diffusion term) 
        amx(row,row-nsp) = ( & 
            - (0d0 - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(Iz-1))*(-1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            ) 
    else 
        ymx(row) = ( 0d0& 
            ! time change term 
            + poro(iz)*(0d0-o2(iz))/dt & 
            ! diffusion term 
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(0d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(0d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            ! rxn term 
            + sporo(iz)*ox2om*kom(iz)*omx(iz)  &
            )
        amx(row,row) = (& 
            ! time change term 
            + poro(iz)*(1d0)/dt & 
            ! diffusion term 
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(-1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            )
        ! filling matrix at grid iz but for unknown at grid iz+1 (only diffusion term) 
        amx(row,row+nsp) = (& 
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0d0)/dz(iz)  &
            )
        ! filling matrix at grid iz but for unknown at grid iz-1 (only diffusion term) 
        amx(row,row-nsp) = (& 
            - (0d0 &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) &
            )
    endif
enddo

ymx = - ymx  ! sign change; see above for the case of om 

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls) ! solving 

o2x = ymx ! passing solutions to unknowns
 
endsubroutine o2calc_ox
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine calcflxo2_ox( &
    o2dec,o2dif,o2tflx,o2res  & ! output 
    ,nz,sporo,kom,omx,dz,poro,dif_o2,dt,o2,o2x,ox2om,o2i  & ! input
    )
implicit none
integer(kind=4),intent(in)::nz
real(kind=8),dimension(nz),intent(in)::sporo,kom,omx,dz,poro,dif_o2,o2,o2x
real(kind=8),intent(in)::dt,ox2om,o2i
real(kind=8),intent(out)::o2dec,o2dif,o2tflx,o2res
integer(kind=4) iz

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

o2res = o2dec + o2dif + o2tflx  ! residual flux

endsubroutine calcflxo2_ox 
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine o2calc_sbox(  &
    o2x  & ! output
    ,izox,nz,poro,o2,kom,omx,sporo,dif_o2,dz,dt,ox2om,o2i & ! input
    )
implicit none
integer(kind=4),intent(in)::nz,izox
real(kind=8),dimension(nz),intent(in)::poro,o2,kom,omx,sporo,dif_o2,dz
real(kind=8),intent(in):: dt,ox2om,o2i
real(kind=8),intent(out)::o2x(nz)
integer(kind=4) :: row,nmx,nsp,iz,infobls
real(kind=8),allocatable :: amx(:,:),ymx(:),emx(:)
integer(kind=4),allocatable::ipiv(:)
    
nsp=1 ! number of species considered here; 1, only om 
nmx = nz*nsp ! # of col (& row) of matrix A to in linear equations Ax = B to be solved, each species has nz (# of grids) unknowns 
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
        ymx(row) = ( &
            ! time change 
            + poro(iz)*(0d0-o2(iz))/dt & 
            ! diffusion 
            - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(0d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_o2(iz)*(0d0-o2i*1d-6/1d3)/dz(iz))/dz(iz)  &
            ! rxn 
            + sporo(iz)*ox2om*kom(iz)*omx(iz)  &
            )
        amx(row,row) = (& 
            ! time change 
            + poro(iz)*(1d0)/dt & 
            ! diffusion 
            - ((poro(iz)*dif_o2(iz)+poro(Iz+1)*dif_o2(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_o2(iz)*(1d0)/dz(iz))/dz(iz)&
            )
        ! filling matrix at grid iz but for unknown at grid iz+1 (only diffusion term) 
        amx(row,row+nsp) = (& 
            - ((poro(Iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - 0d0)/dz(iz)&
            )
    else if (iz>1 .and. iz<= izox) then 
        ymx(row) = ( 0d0& 
            ! time change 
            + poro(iz)*(0d0-o2(iz))/dt & 
            ! diffusion 
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(0d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(0d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            ! rxn 
            + sporo(iz)*ox2om*kom(iz)*omx(iz)  &
            )
        amx(row,row) = (& 
            ! time change 
            + poro(iz)*(1d0)/dt & 
            ! diffusion
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(-1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            )
        ! filling matrix at grid iz but for unknown at grid iz+1 (only diffusion term) 
        amx(row,row+nsp) = (& 
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0d0)/dz(iz)  &
            )
        ! filling matrix at grid iz but for unknown at grid iz-1 (only diffusion term) 
        amx(row,row-nsp) = (& 
            - (0d0 &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) &
            )
    else if (iz> izox) then  ! at lower than zox; zero conc. is forced 
        ymx(row) = ( 0d0& 
            )
        amx(row,row) = (& 
            + 1d0 &
            )
    endif
enddo

ymx = - ymx  ! change signs 

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls) ! solving 

o2x = ymx ! passing solution to variable 

endsubroutine o2calc_sbox
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine calcflxo2_sbox( &
    o2dec,o2dif,o2tflx,o2res  & ! output 
    ,nz,sporo,kom,omx,dz,poro,dif_o2,dt,o2,o2x,izox,ox2om,o2i  & ! input
    )
implicit none
integer(kind=4),intent(in)::nz,izox
real(kind=8),dimension(nz),intent(in)::sporo,kom,omx,dz,poro,dif_o2,o2,o2x
real(kind=8),intent(in)::dt,ox2om,o2i
real(kind=8),intent(out)::o2dec,o2dif,o2tflx,o2res
integer(kind=4) iz

o2dec = 0d0 
o2dif = 0d0
o2tflx = 0d0

do iz = 1,nz 
    if (iz == 1) then 
        o2dec = o2dec + sporo(iz)*ox2om*kom(iz)*omx(iz)*dz(iz)
        o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz)
        o2dif = o2dif - ((poro(iz)*dif_o2(iz)+poro(iz+1)*dif_o2(iz+1))*0.5d0*(o2x(iz+1)-o2x(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_o2(iz)*(o2x(iz)-o2i*1d-6/1d3)/dz(iz))/dz(iz) *dz(iz)
    else if (iz>1 .and. iz<=izox) then 
        o2dec = o2dec + (1d0-poro(iz))*ox2om*kom(iz)*omx(iz)/poro(iz)*dz(iz)*poro(iz)
        o2tflx = o2tflx + (o2x(iz)-o2(iz))/dt*dz(iz)*poro(iz)
        o2dif = o2dif &
            - (0.5d0*(poro(iz+1)*dif_o2(iz+1)+poro(iz)*dif_o2(iz))*(o2x(iz+1)-o2x(iz))/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_o2(iz)+poro(iz-1)*dif_o2(iz-1))*(o2x(Iz)-o2x(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            *dz(iz)
    endif
enddo

o2res = o2dec + o2dif + o2tflx  ! residual flux 

endsubroutine calcflxo2_sbox 
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine calccaco3sys(  &
    ccx,dicx,alkx,rcc,dt  & ! in&output
    ,nspcc,dic,alk,dep,sal,temp,labs,turbo2,nonlocal,sporo,sporoi,sporof,poro,dif_alk,dif_dic & ! input
    ,w,up,dwn,cnr,adf,dz,trans,cc,oxco2,anco2,co3sat,kcc,ccflx,ncc,ohmega,nz  & ! input
    ! ,dum_sfcsumocn  & ! input for genie geochemistry
    ,tol,poroi,flg_500,fact,file_tmp,alki,dici,ccx_th,workdir  &
    ,krad,deccc  & 
    )
implicit none 
integer(kind=4),intent(in)::nspcc,nz,file_tmp
real(kind=8),dimension(nz),intent(in)::dic,alk,sporo,poro,dif_alk,dif_dic,w,up,dwn,cnr,adf,dz,oxco2,anco2
! real(kind=8),dimension(n_ocn),intent(in)::dum_sfcsumocn
real(kind=8),intent(in)::dep,sal,temp,sporoi,sporof,trans(nz,nz,nspcc+2),cc(nz,nspcc),kcc(nz,nspcc),ccflx(nspcc)
real(kind=8),intent(in)::ncc,tol,poroi,fact,alki,dici,ccx_th
logical,dimension(nspcc+2),intent(in)::labs,turbo2,nonlocal
logical,intent(inout)::flg_500
character*255,intent(in)::workdir
real(kind=8),intent(inout)::dicx(nz),alkx(nz),ccx(nz,nspcc),rcc(nz,nspcc),dt
integer(kind=4)::itr,nsp,nmx,infosbr,iiz,n,nnz,infobls,cnt2,cnt,sys,status,isp,iz,row,col,i,j
integer(kind=4),allocatable :: ipiv(:),ap(:),ai(:)
integer(kind=8) symbolic,numeric
real(kind=8)::loc_error,prox(nz),co2x(nz),hco3x(nz),co3x(nz),dco3_ddic(nz),dco3_dalk(nz),drcc_dcc(nz,nspcc)  
real(kind=8)::drcc_dco3(nz,nspcc),drcc_ddic(nz,nspcc),drcc_dalk(nz,nspcc),info(90),control(20)
real(kind=8)::drcc_dohmega(nz,nspcc),dohmega_dalk(nz),dohmega_ddic(nz),ohmega(nz)
real(kind=8),allocatable :: amx(:,:),ymx(:),emx(:),dumx(:,:),ax(:),kai(:),bx(:)
real(kind=8),intent(in)::co3sat
! only when directly tracking isotopes 
real(kind=8),intent(in)::krad(nz,nspcc)  ! caco3 decay consts (for 14c alone)
real(kind=8),intent(out)::deccc(nz,nspcc)  ! radio-active decay rate of caco3 
real(kind=8)::ddeccc_dcc(nz,nspcc)  ! radio-active decay rate of caco3 
integer(kind=4)::iter_max = 500


! for genie geochemistry
! REAL,DIMENSION(n_carbconst)::dum_carbconst
! REAL,DIMENSION(n_carb)::dum_carb
! REAL,DIMENSION(n_carbalk)::dum_carbalk
! real dum_DIC,dum_ALK,dum_Ca,dum_PO4tot,dum_SiO2tot,dum_Btot,dum_SO4tot,dum_Ftot,dum_H2Stot,dum_NH4tot
! real dev
! real,intent(inout)::co3sat

!       Here the system is non-linear and thus Newton's method is used (e.g., Steefel and Lasaga, 1994).
! 
!       Problem: f(x) = 0  (note that x is array, containing N unknowns, here N = nmx)
!       Expanding around x0 (initial/previous guess), 
!           f(x) =  f(x0) + f'(x0)*(x-x0) + O((x-x0)**2)  
!       where f'(x0) represents a Jacobian matrix.
!       Solution is sought by iteration 
!           x = x0 - f'(x0)^-1*f(x0) or x- x0 = - f'(x0)^-1*f(x0)
!       More practically by following steps. 
!           (1) solving Ax = B where A = f'(x0) and B = -f(x0), which gives delta (i.e., x - x0) (again these are arrays)
!           (2) update solution by x = x0 + delta  
!           (3) replace x as x0
!       Three steps are repeated until relative solution difference (i.e., delta) becomes negligible.
!
!       Now matrices A and B are Jacobian and array that contains f(x), respectively, represented by amx and ymx in this code  
!
!       E.g., if equation at grid iz for caco3 is given by (for simplicity it cuts off several terms)
!           (sporo(iz)*ccx(iz)-sporo(iz)*cc(iz))/dt + (sporo(iz)*w(iz)*ccx(iz)-sporo(iz-1)*w(iz-1)*ccx(iz-1))/dz(iz) + rcc(iz) = 0
!       Then, 
!           B(row) = - left-hand side
!       where row is the row number of caco3 at grid iz, i.e., row = 1+(iz-1)*nsp + isp -1 where isp = 1,..., nspcc,
!       and   
!           A(row,row) = -dB(row)/dccx(iz) = (sporo(iz))/dt + (sporo(iz)*w(iz)*1)/dz(iz) + drcc_dcc(iz)    
!           A(row,row-nsp) = -dB(row)/dccx(iz-1) = (-sporo(iz-1)*w(iz-1)*1)/dz(iz)
!           A(row,row+nspcc-1+1) = -dB(row)/ddic(iz) = drcc_ddic(iz) 
!           A(row,row+nspcc-1+2) = -dB(row)/dalk(iz) = drcc_dalk(iz) 
!       Something like this.
!       Note, however, the present code uses ln (conc.) as x following Steefel and Lasaga (1994). So treatment is correspondingly a bit different.
!       E.g.,   
!           dB(row)/dln(alk(iz)) = dB(row)/dalk(iz)*dalk(iz)/dln(alk(iz)) =  dB(row)/dalk(iz) * alkx(iz) = drcc_dalk(iz)*alkx(iz)
!           ln x = ln x0 + delta, or, x = x0*exp(delta)
!
!       See e.g., Steefel and Lasaga (1994) for more details. 

flg_500 = .false.
loc_error = 1d4
itr = 0

nsp = 2 + nspcc  ! now considered species are dic, alk and nspcc of caco3 
nmx = nz*nsp  ! col (and row) of matrix; the same number of unknowns 

! deallocate(amx,ymx,emx,ipiv)
if (allocated(amx))deallocate(amx)
if (allocated(ymx))deallocate(ymx)
if (allocated(emx))deallocate(emx)
if (allocated(ipiv))deallocate(ipiv)
allocate(amx(nmx,nmx),ymx(nmx),emx(nmx),ipiv(nmx))

if (allocated(dumx))deallocate(dumx)  ! used for sparse matrix solver 
allocate(dumx(nmx,nmx))

do while (loc_error > tol)
    
amx = 0d0
ymx = 0d0

! calling subroutine from caco3_therm.f90 to calculate aqueous co2 species 
#ifndef mocsy
call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,infosbr)
if (infosbr==1) then ! which means error in calculation 
    ! dt=dt/10d0
#ifdef sense
    ! go to 500
    flg_500=.true.
    return
#else
    stop
#endif 
endif 
! calling subroutine from caco3_therm.f90 to calculate derivatives of co3 wrt alk and dic 
call calcdevs(dicx,alkx,temp,sal,dep,nz,infosbr,dco3_dalk,dco3_ddic)
if (infosbr==1) then ! if error in calculation 
    ! dt=dt/10d0
#ifdef sense
    ! go to 500
    flg_500=.true.
    return
#else
    stop
#endif 
endif 
do isp=1,nspcc
    ! calculation of dissolution rate for individual species 
    rcc(:,isp) = kcc(:,isp)*ccx(:,isp)*abs(1d0-co3x(:)*1d3/co3sat)**ncc*merge(1d0,0d0,(1d0-co3x(:)*1d3/co3sat)>0d0)
    ! calculation of derivatives of dissolution rate wrt conc. of caco3 species, dic and alk 
    drcc_dcc(:,isp) = kcc(:,isp)*abs(1d0-co3x(:)*1d3/co3sat)**ncc*merge(1d0,0d0,(1d0-co3x(:)*1d3/co3sat)>0d0)
    drcc_dco3(:,isp) = kcc(:,isp)*ccx(:,isp)*ncc*abs(1d0-co3x(:)*1d3/co3sat)**(ncc-1d0)  &
        *merge(1d0,0d0,(1d0-co3x(:)*1d3/co3sat)>0d0)*(-1d3/co3sat)
    drcc_ddic(:,isp) = drcc_dco3(:,isp)*dco3_ddic(:)
    drcc_dalk(:,isp) = drcc_dco3(:,isp)*dco3_dalk(:)
    deccc(:,isp) = krad(:,isp)*ccx(:,isp)
    ddeccc_dcc(:,isp) = krad(:,isp)
enddo
#else
call co2sys_mocsy(nz,alkx*1d6,dicx*1d6,temp,dep*1d3,sal  &
                        ,co2x,hco3x,co3x,prox,ohmega,dohmega_ddic,dohmega_dalk) ! using mocsy
co2x = co2x/1d6
hco3x = hco3x/1d6
co3x = co3x/1d6
dohmega_ddic = dohmega_ddic*1d6
dohmega_dalk = dohmega_dalk*1d6
do isp=1,nspcc
    ! calculation of dissolution rate for individual species 
    rcc(:,isp) = kcc(:,isp)*ccx(:,isp)*abs(1d0-ohmega(:))**ncc*merge(1d0,0d0,(1d0-ohmega(:))>0d0)
    ! calculation of derivatives of dissolution rate wrt conc. of caco3 species, dic and alk 
    drcc_dcc(:,isp) = kcc(:,isp)*abs(1d0-ohmega(:))**ncc*merge(1d0,0d0,(1d0-ohmega(:))>0d0)
    drcc_dohmega(:,isp) = kcc(:,isp)*ccx(:,isp)*ncc*abs(1d0-ohmega(:))**(ncc-1d0)  &
        *merge(1d0,0d0,(1d0-ohmega(:))>0d0)*(-1d0)
    drcc_ddic(:,isp) = drcc_dohmega(:,isp)*dohmega_ddic(:)
    drcc_dalk(:,isp) = drcc_dohmega(:,isp)*dohmega_dalk(:)
    deccc(:,isp) = krad(:,isp)*ccx(:,isp)
    ddeccc_dcc(:,isp) = krad(:,isp)
enddo
#endif 

do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then ! when upper condition must be taken account; *** comments for matrix filling are given only in this case 
        do isp = 1,nspcc  ! multiple caco3 species 
            ! put f(x) for isp caco3 species 
            ymx(row+isp-1) = (&
                + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt &
                - ccflx(isp)/dz(1) &
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-0d0)/dz(1)  &
                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(1)  &
                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-0d0)/dz(1)  &
                + sporo(iz)*rcc(iz,isp) &
                + sporo(iz)*deccc(iz,isp) &
                )
            ! derivative of f(x) wrt isp caco3 conc. at grid iz in ln 
            amx(row+isp-1,row+isp-1) = (&
                + sporo(iz)*(1d0)/dt &
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(1)   &
                + adf(iz)*dwn(iz)*(0d0-sporo(iz)*w(iz)*1d0)/dz(1)  &
                + sporo(iz)* drcc_dcc(iz,isp)  &
                + sporo(iz)*ddeccc_dcc(iz,isp)  &
                )* ccx(iz,isp) 
            ! derivative of f(x) wrt isp caco3 conc. at grid iz+1 in ln 
            amx(row+isp-1,row+isp-1+nsp) =  (&
                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(1)  &
                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(1)  &
                )*ccx(iz+1,isp)
            ! derivative of f(x) wrt dic conc. at grid iz in ln
            amx(row+isp-1,row+nspcc) = (&
                + sporo(iz)*drcc_ddic(iz,isp)  &
                )*dicx(iz)
            ! derivative of f(x) wrt alk conc. at grid iz in ln
            amx(row+isp-1,row+nspcc+1) = (&
                + sporo(iz)*drcc_dalk(iz,isp)  &
                )*alkx(iz)
            ! DIC
            ! derivative of f(x) for dic at iz wrt isp caco3 conc. at grid iz in ln
            amx(row+nspcc,row+isp-1) = (&
                - (1d0-poro(Iz))*drcc_dcc(iz,isp)  &
                )*ccx(iz,isp)*fact
            ! ALK 
            ! derivative of f(x) for alk at iz wrt isp caco3 conc. at grid iz in ln
            amx(row+nspcc+1,row+isp-1) = (&
                - 2d0* (1d0-poro(Iz))*drcc_dcc(iz,isp)  &
                - sporo(iz)*ddeccc_dcc(iz,isp)    &
                )*ccx(iz,isp)*fact
        enddo 
        !  DIC 
        ! put f(x) for dic at iz  
        ymx(row+nspcc) = ( &
            + poro(iz)*(dicx(iz)-dic(iz))/dt & 
            - ((poro(iz)*dif_dic(iz)+poro(iz+1)*dif_dic(iz+1))*0.5d0*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_dic(iz)*(dicx(iz)-dici*1d-6/1d3)/dz(iz))/dz(iz)  &
            - oxco2(iz) &
            - anco2(iz) &
            - (1d0-poro(Iz))*sum(rcc(iz,:))  &
            )*fact
        ! put derivative of f(x) for dic at iz wrt dic at iz in ln 
        amx(row+nspcc,row+nspcc) = (& 
            + poro(iz)*(1d0)/dt & 
            - ((poro(iz)*dif_dic(iz)+poro(Iz+1)*dif_dic(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(Iz)*dif_dic(iz)*(1d0)/dz(iz))/dz(iz)&
            - (1d0-poro(Iz))*sum(drcc_ddic(iz,:))  &
            )*dicx(iz)*fact
        ! put derivative of f(x) for dic at iz wrt dic at iz+1 in ln 
        amx(row+nspcc,row+nspcc+nsp) = (& 
            - ((poro(iz)*dif_dic(iz)+poro(Iz+1)*dif_dic(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - 0d0)/dz(iz)&
            )*dicx(iz+1)*fact
        ! put derivative of f(x) for dic at iz wrt alk at iz in ln 
        amx(row+nspcc,row+nspcc+1) = ( &
            - (1d0-poro(Iz))*sum(drcc_dalk(iz,:))  &
            )*alkx(iz)*fact
        ! ALK
        ! put f(x) for alk at iz  
        ymx(row+nspcc+1) = (& 
            + poro(iz)*(alkx(iz)-alk(iz))/dt & 
            - ((poro(iz)*dif_alk(iz)+poro(Iz+1)*dif_alk(iz+1))*0.5d0*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_alk(iz)*(alkx(iz)-alki*1d-6/1d3)/dz(iz))/dz(iz) &
            - anco2(iz) &
            - 2d0* (1d0-poro(Iz))*sum(rcc(iz,:))  &
            - sporo(iz)*sum(deccc(iz,:)) &
            )*fact
        ! put derivative of f(x) for alk at iz wrt alk at iz in ln 
        amx(row+nspcc+1,row+nspcc+1) = (& 
            + poro(iz)*(1d0)/dt & 
            - ((poro(iz)*dif_alk(iz)+poro(iz+1)*dif_alk(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_alk(iz)*(1d0)/dz(iz))/dz(iz)  &
            - 2d0* (1d0-poro(Iz))*sum(drcc_dalk(iz,:))  &
            )*alkx(iz)*fact
        ! put derivative of f(x) for alk at iz wrt alk at iz+1 in ln 
        amx(row+nspcc+1,row+nspcc+1+nsp) = (& 
            - ((poro(Iz)*dif_alk(iz)+poro(Iz+1)*dif_alk(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - 0d0)/dz(iz)&
            )*alkx(iz+1)*fact
        ! put derivative of f(x) for alk at iz wrt dic at iz in ln 
        amx(row+nspcc+1,row+nspcc) = (&
            - 2d0* (1d0-poro(Iz))*sum(drcc_ddic(iz,:))  &
            )*dicx(iz)*fact
    else if (iz == nz) then ! need be careful about lower boundary condition; no diffusive flux from the bottom  
        do isp=1,nspcc
            ymx(row+isp-1) = (& 
                + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt &
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz)  &
                + adf(iz)*cnr(iz)*(sporof*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz)  &
                + adf(iz)*dwn(iz)*(sporof*w(iz)*ccx(iz,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(iz)  &
                + sporo(iz)*rcc(iz,isp)  &
                + sporo(iz)*deccc(iz,isp) &
                )
            amx(row+isp-1,row+isp-1) = (&
                + sporo(iz)*(1d0)/dt &
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  &
                + adf(iz)*cnr(iz)*(sporof*w(iz)*1d0-0d0)/dz(iz)  &
                + adf(iz)*dwn(iz)*(sporof*w(iz)*1d0-sporo(iz)*w(iz)*1d0)/dz(iz)  &
                + sporo(iz)*drcc_dcc(iz,isp)   &
                + sporo(iz)*ddeccc_dcc(iz,isp)  &
                )*ccx(iz,isp)
            amx(row+isp-1,row+isp-1-nsp) = ( &
                + adf(iz)*up(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
                + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
                )*ccx(iz-1,isp)
            amx(row+isp-1,row+nspcc) = (&
                + sporo(iz)*drcc_ddic(iz,isp) &
                )*dicx(iz)
            amx(row+isp-1,row+nspcc+1) = (&
                + sporo(iz)*drcc_dalk(iz,isp) &
                )*alkx(iz)
            
            !DIC 
            amx(row+nspcc,row+isp-1) = (&
                - sporo(Iz)*drcc_dcc(iz,isp)  &
                )*ccx(iz,isp)*fact
            !ALK 
            amx(row+nspcc+1,row+isp-1) = (&
                - 2d0*sporo(Iz)*drcc_dcc(iz,isp)  &
                - sporo(iz)*ddeccc_dcc(iz,isp)    &
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
            - sporo(iz)*sum(deccc(iz,:)) &
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
    else !  do not have to be careful abount boundary conditions 
        do isp=1,nspcc
            ymx(row+isp-1) = (& 
                + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt &
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-sporo(Iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz)  &
                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(Iz)*w(iz)*ccx(iz,isp))/dz(iz)  &
                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(Iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz)  &
                + sporo(iz)*rcc(iz,isp)  &
                + sporo(iz)*deccc(iz,isp) &
                )
            amx(row+isp-1,row+isp-1) = (&
                + sporo(iz)*(1d0)/dt &
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  &
                + adf(iz)*dwn(iz)*(0d0-sporo(Iz)*w(iz)*1d0)/dz(iz)  &
                + sporo(iz)*drcc_dcc(iz,isp)  &
                + sporo(iz)*ddeccc_dcc(iz,isp)  &
                )*ccx(iz,isp)
            amx(row+isp-1,row+isp-1+nsp) =  (&
                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(iz)  &
                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(iz)  &
                )*ccx(iz+1,isp)
            amx(row+isp-1,row+isp-1-nsp) =  (&
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
                - sporo(iz)*ddeccc_dcc(iz,isp)    &
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
            - sporo(iz)*sum(deccc(iz,:)) &
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
    ! diffusion terms are filled with transition matrices 
    do isp=1,nspcc
        if (turbo2(isp+2).or.labs(isp+2)) then
            do iiz = 1, nz
                col = 1 + (iiz-1)*nsp
                if (trans(iiz,iz,isp+2)==0d0) cycle
                amx(row+isp-1,col+isp-1) = amx(row+isp-1,col+isp-1) &
                    - trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*ccx(iiz,isp)
                ymx(row+isp-1) = ymx(row+isp-1) &
                    - trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*ccx(iiz,isp)
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

ymx = - ymx  ! because I put f(x) into ymx (=B), minus sign need be added 

#ifndef nonrec
if (any(isnan(ymx))) then 
    print*,'NAN in ymx'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_ymx_pre.txt',status = 'unknown')
    do iz = 1, nmx
        write (file_tmp,*) ymx(iz)
    enddo
    close(file_tmp)
    ! stop
    flg_500 = .true.
    exit 
endif
#endif 

#ifndef sparse 
! using non-sparse solver 
call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,infobls) 
#else 
!  slowest way of using sparse matrix solver
n = nmx
where(amx/=0d0)
    dumx=1
elsewhere
    dumx=0
endwhere
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
#endif 

do iz = 1, nz 
    row = 1+(iz-1)*nsp
    do isp=1,nspcc
        if (ymx(row+isp-1)>10d0) then ! this help conversion 
            ccx(iz,isp) = ccx(iz,isp)*1.5d0
        elseif (ymx(row+isp-1)<-10d0) then ! this help conversion  
            ccx(iz,isp) = ccx(iz,isp)*0.5d0
        else
            ccx(iz,isp) = ccx(iz,isp)*exp(ymx(row+isp-1))
        endif
        if (ccx(iz,isp)<ccx_th) then ! too small trancate value and not be accounted for error 
            ccx(iz,isp)=ccx_th
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
    if (alkx(iz)<1d-100) ymx(row+nspcc+1) = 0d0
enddo

if (flg_500) exit 

loc_error = maxval(exp(abs(ymx))) - 1d0
itr = itr + 1
#ifdef showiter
print*,'co2 iteration',itr,loc_error,infobls
#endif

!  if negative or NAN calculation stops 
if (any(ccx<0d0)) then
    print*,'negative ccx, stop'
    print*,ccx
    ! stop
    flg_500 = .true.
    exit 
endif
if (any(isnan(ccx))) then
    print*,'nan ccx, stop'
    print*,ccx
    ! stop
    flg_500 = .true.
    exit 
endif

if (any(dicx<0d0)) then
    print*,'negative dicx, stop'
    print*,dicx
    ! stop
    flg_500 = .true.
    exit 
endif 
if (any(isnan(dicx))) then
    print*,'nan dic, stop'
    print*,dicx
    ! stop
    flg_500 = .true.
    exit 
endif 

if (any(alkx<0d0)) then
    print*,'negative alk, stop'
    print*,alkx
    ! stop
    flg_500 = .true.
    exit 
endif
if (any(isnan(alkx))) then
    print*,'nan alk, stop'
    print*,alkx
    ! stop
    flg_500 = .true.
    exit 
endif
#ifdef sense
if (itr > iter_max) then 
    print*, '....cannot converge...'
    flg_500 = .true.
    exit
endif 
#endif
enddo

endsubroutine calccaco3sys
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine calcflxcaco3sys(  &
     cctflx,ccflx,ccdis,ccdif,ccadv,ccrain,ccres,alktflx,alkdis,alkdif,alkdec,alkres & ! output
     ,dictflx,dicdis,dicdif,dicres,dicdec   & ! output
     ,dw & ! inoutput
     ,nspcc,ccx,cc,dt,dz,rcc,adf,up,dwn,cnr,w,dif_alk,dif_dic,dic,dicx,alk,alkx,oxco2,anco2,trans    & ! input
     ,turbo2,labs,nonlocal,sporof,it,nz,poro,sporo        & ! input
     ,dici,alki,file_err,mvcc,tol,flg_500  &
     ,ccrad,alkrad,deccc  &  
     )
implicit none 
integer(kind=4),intent(in)::nz,nspcc,it,file_err
real(kind=8),dimension(nz),intent(in)::poro,dz,adf,up,dwn,cnr,w,dif_alk,dif_dic,dic,dicx,alk,alkx,oxco2,anco2
real(kind=8),dimension(nz),intent(in)::sporo
real(kind=8),intent(in)::ccx(nz,nspcc),cc(nz,nspcc),dt,rcc(nz,nspcc),trans(nz,nz,nspcc+2),sporof,dici,alki,mvcc(nspcc),tol
real(kind=8),intent(inout)::dw(nz)
logical,dimension(nspcc+2),intent(in)::turbo2,labs,nonlocal
real(kind=8),dimension(nspcc),intent(out)::cctflx,ccflx,ccdis,ccdif,ccadv,ccrain,ccres
real(kind=8),intent(out)::dictflx,dicdis,dicdif,dicres,alktflx,alkdis,alkdif,alkdec,alkres,dicdec
logical,intent(inout)::flg_500
integer(kind=4)::iz,row,nsp,isp,iiz,col
! when switching on isotrack
real(kind=8),intent(out)::ccrad(nspcc),alkrad
real(kind=8),intent(in)::deccc(nz,nspcc)

nsp = nspcc+2

cctflx =0d0 
ccdis = 0d0 
ccdif = 0d0 
ccadv = 0d0 
ccrain = 0d0
ccrad = 0d0
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
alkrad = 0d0
alkres = 0d0

do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then 
        do isp=1,nspcc
            cctflx(isp) = cctflx(isp) + (1d0-poro(iz))*(ccx(iz,isp)-cc(iz,isp))/dt *dz(iz)
            ccdis(isp) = ccdis(isp)  + (1d0-poro(Iz))*rcc(iz,isp) *dz(iz)
            ccrad(isp) = ccrad(isp) + sporo(iz)*deccc(iz,isp)*dz(iz)
            ccrain(isp) = ccrain(isp) - ccflx(isp)/dz(1)*dz(iz)
            ccadv(isp) = ccadv(Isp) + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-0d0)/dz(1) * dz(iz) &
                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(1) * dz(iz)  &
                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-0d0)/dz(1) * dz(iz)
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
        alkrad = alkrad - sporo(iz)*sum(deccc(iz,:))*dz(iz)
    else if (iz == nz) then 
        do isp=1,nspcc
            cctflx(isp) = cctflx(isp) + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt *dz(iz)
            ccdis(isp) = ccdis(isp)  + sporo(Iz)*rcc(iz,isp) *dz(iz)
            ccrad(isp) = ccrad(isp) + sporo(iz)*deccc(iz,isp)*dz(iz)
            ccadv(isp) = ccadv(isp) &
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz) * dz(iz)  &
                + adf(iz)*cnr(iz)*(sporof*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz) * dz(iz)  &
                + adf(iz)*dwn(iz)*(sporof*w(iz)*ccx(iz,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(iz) * dz(iz)  
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
        alkrad = alkrad - sporo(iz)*sum(deccc(iz,:))*dz(iz)
    else 
        do isp=1,nspcc
            cctflx(isp) = cctflx(isp) + sporo(iz)*(ccx(iz,isp)-cc(iz,isp))/dt *dz(iz)
            ccdis(isp) = ccdis(isp)  + sporo(Iz)*rcc(iz,isp) *dz(iz)
            ccrad(isp) = ccrad(isp) + sporo(iz)*deccc(iz,isp)*dz(iz)
            ccadv(isp) = ccadv(isp) &
                + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ccx(iz,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz) * dz(iz)  &
                + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz)*w(iz)*ccx(iz,isp))/dz(iz) * dz(iz)  &
                + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ccx(iz+1,isp)-sporo(iz-1)*w(iz-1)*ccx(iz-1,isp))/dz(iz) * dz(iz)  
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
        alkrad = alkrad - sporo(iz)*sum(deccc(iz,:))*dz(iz)
    endif
    do isp=1,nspcc
        if (labs(isp+2).or. turbo2(isp+2)) then 
            do iiz = 1, nz
                if (trans(iiz,iz,isp+2)==0d0) cycle
                ccdif(isp) = ccdif(isp) -trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*dz(iz)*ccx(iiz,isp)
            enddo
        else 
            do iiz = 1, nz
                if (trans(iiz,iz,isp+2)==0d0) cycle
                ccdif(isp) = ccdif(isp) -trans(iiz,iz,isp+2)/dz(iz)*dz(iz)*ccx(iiz,isp)
            enddo
        endif
        if (labs(isp+2).or. turbo2(isp+2)) then 
            do iiz = 1, nz
                if (trans(iiz,iz,isp+2)==0d0) cycle
                dw(iz) = dw(iz) - mvcc(isp)*(-trans(iiz,iz,isp+2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*ccx(iiz,isp))
            enddo
        else 
            if (nonlocal(isp+2)) then 
                do iiz = 1, nz
                    if (trans(iiz,iz,isp+2)==0d0) cycle
                    dw(iz) = dw(iz) -mvcc(isp)*(-trans(iiz,iz,isp+2)/dz(iz)*ccx(iiz,isp))
                enddo
            endif 
        endif 
    enddo 
    dw(iz) = dw(iz) -(1d0-poro(iz))*sum(mvcc(:)*rcc(iz,:)) - sporo(iz)*sum(mvcc(:)*deccc(iz,:))
enddo

! residual fluxes 
ccres = cctflx +  ccdis +  ccdif + ccadv + ccrain + ccrad
dicres = dictflx + dicdis + dicdif + dicdec 
alkres = alktflx + alkdis + alkdif + alkdec + alkrad

flg_500 = .false.
! #ifdef sense
! if (abs(alkres)/max(alktflx,alkdis ,alkdif , alkdec) > tol*10d0) then   ! if residual fluxes are relatively large, record just in case  
    ! print*,'not enough accuracy in co2 calc:stop',abs(alkres)/max(alktflx,alkdis ,alkdif , alkdec)
    ! write(file_err,*)it,'not enough accuracy in co2 calc:stop',abs(alkres)/max(alktflx,alkdis ,alkdif , alkdec)  &
        ! ,alkres, alktflx,alkdis , alkdif , alkdec 
    ! flg_500 = .true.
! endif
! #endif 
        
if (abs(alkres)/maxval(abs(ccflx)) > tol*10d0) then   ! if residual fluxes are relatively large, record just in case  
    print*,'not enough accuracy in co2 calc:stop',abs(alkres)/maxval(abs(ccflx))
    write(file_err,*)it,'not enough accuracy in co2 calc:stop',abs(alkres)/maxval(abs(ccflx))  &
        ,alkres, alktflx,alkdis , alkdif , alkdec, alkrad 
    flg_500 = .true.
endif

endsubroutine calcflxcaco3sys 
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine claycalc(  &   
    ptx                  &  ! output
    ,nz,sporo,pt,dt,w,dz,detflx,adf,up,dwn,cnr,trans  &  ! input
    ,nspcc,labs,turbo2,nonlocal,poro,sporof     &  !  intput
    ,msed,file_tmp,workdir &
    )
implicit none
integer(kind=4),intent(in)::nz,nspcc,file_tmp
real(kind=8),dimension(nz),intent(in)::sporo,pt,w,dz,adf,up,dwn,cnr,poro
real(kind=8),intent(in)::dt,detflx,trans(nz,nz,nspcc+2),sporof,msed
logical,dimension(nspcc+2),intent(in)::labs,turbo2,nonlocal
real(kind=8),intent(out)::ptx(nz)
character*255,intent(in)::workdir
integer(kind=4)::nsp,nmx,iz,row,iiz,infobls,col
integer(kind=4),allocatable::ipiv(:)
real(kind=8),allocatable::amx(:,:),ymx(:),emx(:)

nsp = 1 !  only consider clay
nmx = nz*nsp  ! matrix is linear and solved like om and o2, so see comments there for calculation procedures 
if (allocated(amx))deallocate(amx)
if (allocated(ymx))deallocate(ymx)
if (allocated(emx))deallocate(emx)
if (allocated(ipiv))deallocate(ipiv)
! deallocate(amx,ymx,emx,ipiv)
allocate(amx(nmx,nmx),ymx(nmx),emx(nmx),ipiv(nmx))
    
amx = 0d0
ymx = 0d0

do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then 
        ymx(row) = &
            + sporo(iz)*(-pt(iz))/dt &
            - detflx/msed/dz(iz)
        amx(row,row) = (&
            + sporo(iz)*(1d0)/dt &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)   &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz)*w(iz)*1d0)/dz(iz)   &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-0d0)/dz(iz)   &
            )            
        amx(row,row+nsp) =  (&
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz)*w(iz)*0d0)/dz(iz)   &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(iz)   &
            )
    else if (iz == nz) then 
        ymx(row) = & 
            + sporo(iz)*(-pt(iz))/dt 
        amx(row,row) = (&
            + sporo(iz)*(1d0)/dt &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporof*w(iz)*1d0-0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(sporof*w(iz)*1d0-sporo(iz)*w(iz)*1d0)/dz(iz)  &
            )
        amx(row,row-nsp) = ( &
            + adf(iz)*up(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            )
    else 
        ymx(row) = & 
            + sporo(iz)*(-pt(iz))/dt 
        amx(row,row) = (&
            + sporo(iz)*(1d0)/dt &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*0d0-sporo(iz)*w(iz)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*0d0-0d0)/dz(iz)  &
            )
        amx(row,row+nsp) =  (&
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*1d0-sporo(iz)*w(iz)*0d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*1d0-0d0)/dz(iz)  &
            )
        amx(row,row-nsp) =  (&
            + adf(iz)*up(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            )
    endif
    if (labs(2).or.turbo2(2)) then 
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp
            if (trans(iiz,iz,2)==0d0) cycle
            amx(row,col) = amx(row,col) -trans(iiz,iz,2)/dz(iz)*dz(iiz)*(1d0-poro(iiz))
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

endsubroutine claycalc
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine calcflxclay( &
    pttflx,ptdif,ptadv,ptres,ptrain  & ! output
    ,dw          &  ! in&output
    ,nz,sporo,ptx,pt,dt,dz,detflx,w,adf,up,dwn,cnr,sporof,trans,nspcc,turbo2,labs,nonlocal,poro           &  !  input
    ,msed,mvsed  &
    )
implicit none 
integer(kind=4),intent(in)::nz,nspcc
real(kind=8),dimension(nz),intent(in)::sporo,ptx,pt,dz,w,adf,up,dwn,cnr,poro
real(kind=8),intent(in)::dt,detflx,sporof,trans(nz,nz,nspcc+2),msed,mvsed
logical,dimension(nspcc+2),intent(in)::turbo2,labs,nonlocal
real(kind=8),intent(inout)::dw(nz)
real(kind=8),intent(out)::pttflx,ptdif,ptadv,ptres,ptrain
integer(kind=4)::iz,row,nsp=1,col,iiz

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
        ptadv = ptadv &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ptx(iz)-0d0)/dz(iz)*dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ptx(iz+1)-sporo(iz)*w(iz)*ptx(iz))/dz(iz)*dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ptx(iz+1)-0d0)/dz(iz)*dz(iz)  
    else if (iz == nz) then 
        pttflx = pttflx + (1d0-poro(iz))*(ptx(iz)-pt(iz))/dt*dz(iz)
        ptadv = ptadv &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ptx(iz)-sporo(iz-1)*w(iz-1)*ptx(iz-1))/dz(iz)*dz(iz)  &
            + adf(iz)*cnr(iz)*(sporof*w(iz)*ptx(iz)-sporo(iz-1)*w(iz-1)*ptx(iz-1))/dz(iz)*dz(iz)  &
            + adf(iz)*dwn(iz)*(sporof*w(iz)*ptx(iz)-sporo(iz)*w(iz)*ptx(iz))/dz(iz)*dz(iz)  
    else 
        pttflx = pttflx + (1d0-poro(iz))*(ptx(iz)-pt(iz))/dt*dz(iz)
        ptadv = ptadv &
            + adf(iz)*up(iz)*(sporo(iz)*w(iz)*ptx(iz)-sporo(iz-1)*w(Iz-1)*ptx(iz-1))/dz(iz)*dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*w(iz+1)*ptx(iz+1)-sporo(iz)*w(Iz)*ptx(iz))/dz(iz)*dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*w(iz+1)*ptx(iz+1)-sporo(iz-1)*w(Iz-1)*ptx(iz-1))/dz(iz)*dz(iz)
    endif
    if(turbo2(2).or.labs(2)) then 
        do iiz = 1, nz
            if (trans(iiz,iz,2)==0d0) cycle
            ptdif = ptdif -trans(iiz,iz,2)*ptx(iiz)/dz(iz)*dz(iiz)*dz(iz)
        enddo
    else 
        do iiz = 1, nz
            if (trans(iiz,iz,2)==0d0) cycle
            ptdif = ptdif -trans(iiz,iz,2)*ptx(iiz)/dz(iz)    &
                *dz(iz)
        enddo
    endif
    if(turbo2(2).or.labs(2)) then 
        do iiz = 1, nz
            if (trans(iiz,iz,2)==0d0) cycle
            dw(iz) = dw(iz) - mvsed*(-trans(iiz,iz,2)*ptx(iiz)/dz(iz)*dz(iiz)*(1d0-poro(iiz)))
        enddo
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

endsubroutine calcflxclay
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine getsldprop(  &
    rho,frt,       &  ! output
    nz,omx,ptx,ccx,nspcc,w,up,dwn,cnr,adf,z      & ! input
    ,mom,msed,mcc,mvom,mvsed,mvcc,file_tmp,workdir  &
    )
implicit none 
integer(kind=4),intent(in)::nz,nspcc,file_tmp
real(kind=8),dimension(nz),intent(in)::omx,ptx,w,up,dwn,cnr,adf,z
real(kind=8),intent(in)::ccx(nz,nspcc),mom,msed,mcc(nspcc),mvom,mvsed,mvcc(nspcc)
real(kind=8),intent(out)::rho(nz),frt(nz)
character*255,intent(in)::workdir
integer(kind=4)::iz

do iz=1,nz 
    rho(iz) = omx(iz)*mom + ptx(iz)*msed +  sum(ccx(iz,:)*mcc(:))  ! calculating bulk density 
    frt(iz) = omx(Iz)*mvom + ptx(iz)*mvsed + sum(ccx(iz,:)*mvcc(:))  ! calculation of total vol. fraction of solids 
enddo 

! check error for density (rho)
if (any(rho<0d0)) then  ! if negative density stop ....
    print*,'negative density'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'NEGATIVE_RHO.txt',status = 'unknown')
    do iz = 1, nz
        write (file_tmp,*) z(iz),rho(iz),w(iz),up(iz),dwn(iz),cnr(iz),adf(iz)
    enddo
    close(file_tmp)
    stop
endif 

endsubroutine getsldprop
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine burialcalc(  &
    w,wi        & !  output
    ,detflx,ccflx,nspcc,omflx,dw,dz,poro,nz    & ! input
    ,msed,mvsed,mvcc,mvom,poroi &
    )
implicit none
integer(kind=4),intent(in)::nz,nspcc
real(kind=8),intent(in)::detflx,ccflx(nspcc),omflx,dw(nz),dz(nz),poro(nz),msed,mvsed,mvcc(nspcc),mvom,poroi
real(kind=8),intent(out)::wi,w(nz)
integer(kind=4)::iz

! w is up dated by solving 
!           d(1 - poro)*w/dz = dw
! note that dw has recorded volume changes by reactions and non-local mixing (see Eqs. B2 and B6 in ms)
! finite difference form is 
!           if (iz/=1) {(1-poro(iz))*w(iz)-(1-poro(iz-1))*w(iz-1)}/dz(iz) = dw(iz)          
!           if (iz==1) (1-poro(iz))*w(iz) = total volume flux + dw(iz)*dz(iz)          
! which leads to the following calculations 

wi = (detflx/msed*mvsed + sum(ccflx*mvcc) +omflx*mvom)/(1d0-poroi)  ! upper value; (1d0-poroi) is almost meaningless, see below 
do iz=1,nz
    if (iz==1) then 
        w(iz)=((1d0-poroi)*wi + dw(iz)*dz(iz))/(1d0-poro(iz))
    else 
        w(iz)=((1d0-poro(iz-1))*w(iz-1) + dw(iz)*dz(iz))/(1d0-poro(iz))
    endif
enddo

endsubroutine burialcalc
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine burialcalc_fdm(  &
    w       & !  output
    ,detflx,ccflx,nspcc,omflx,dw,dz,poro,nz    & ! input
    ,msed,mvsed,mvcc,mvom,poroi,up,dwn,cnr,adf,sporo,sporof,file_tmp,workdir &
    )
implicit none
integer(kind=4),intent(in)::nz,nspcc,file_tmp
real(kind=8),intent(in)::detflx,ccflx(nspcc),omflx,dw(nz),dz(nz),poro(nz),msed,mvsed,mvcc(nspcc),mvom,poroi
real(kind=8),intent(out)::w(nz)
real(kind=8),dimension(nz),intent(in)::up,dwn,cnr,adf,sporo
real(kind=8),intent(in)::sporof
character*255,intent(in)::workdir
integer(kind=4)::iz

integer(kind=4)::nsp,nmx,row,infobls
integer(kind=4),allocatable::ipiv(:)
real(kind=8),allocatable::amx(:,:),ymx(:),emx(:)

! w is up dated by solving 
!           d(1 - poro)*w/dz = dw
! note that dw has recorded volume changes by reactions and non-local mixing (see Eqs. B2 and B6 in ms)
! finite difference form is 
!           if (iz/=1) {(1-poro(iz))*w(iz)-(1-poro(iz-1))*w(iz-1)}/dz(iz) = dw(iz)          
!           if (iz==1) (1-poro(iz))*w(iz) = total volume flux + dw(iz)*dz(iz)          
! which leads to the following calculations 

nsp = 1 !  only consider burial rate
nmx = nz*nsp  ! matrix is linear and solved like om and o2, so see comments there for calculation procedures 
if (allocated(amx))deallocate(amx)
if (allocated(ymx))deallocate(ymx)
if (allocated(emx))deallocate(emx)
if (allocated(ipiv))deallocate(ipiv)
! deallocate(amx,ymx,emx,ipiv)
allocate(amx(nmx,nmx),ymx(nmx),emx(nmx),ipiv(nmx))
    
amx = 0d0
ymx = 0d0

do iz=1,nz
    row=iz
    if (iz==1) then 
        amx(row,row) = adf(iz)*up(iz)*(sporo(iz)*1d0-0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(0d0-sporo(iz)*1d0)/dz(iz)  
        amx(row,row+1) =   &
            + adf(iz)*cnr(iz)*(sporo(iz+1)*1d0-0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(sporo(iz+1)*1d0-0d0)/dz(iz)  
        ymx(row) = adf(iz)*up(iz)*(-(detflx/msed*mvsed + sum(ccflx*mvcc) +omflx*mvom))/dz(iz)  &  
            + adf(iz)*cnr(iz)*(-(detflx/msed*mvsed + sum(ccflx*mvcc) +omflx*mvom))/dz(iz)  &  
            - dw(iz)
    elseif (iz==nz) then 
        amx(row,row) = adf(iz)*up(iz)*(sporo(iz)*1d0-0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(sporof*1d0-sporo(iz)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporof*1d0-0d0)/dz(iz)  
        amx(row,row-1) = adf(iz)*up(iz)*(0d0-sporo(iz-1)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*1d0)/dz(iz)  
        ymx(row) = - dw(iz)
    else  
        amx(row,row) = adf(iz)*up(iz)*(sporo(iz)*1d0-0d0)/dz(iz)  &
            + adf(iz)*dwn(iz)*(0d0-sporo(iz)*1d0)/dz(iz) 
        amx(row,row-1) = adf(iz)*up(iz)*(0d0-sporo(iz-1)*1d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(0d0-sporo(iz-1)*1d0)/dz(iz)  
        amx(row,row+1) = adf(iz)*dwn(iz)*(sporo(iz+1)*1d0-0d0)/dz(iz)  &
            + adf(iz)*cnr(iz)*(sporo(iz)*1d0-0d0)/dz(iz)  
        ymx(row) = - dw(iz)
    endif 
enddo

ymx=-ymx

#ifndef nonrec
if (any(isnan(ymx))) then 
    print*,'NAN in ymx:w'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_ymx_pre_w.txt',status = 'unknown')
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
    print*,'NAN in amx:w'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_amx_w.txt',status = 'unknown')
    do iz = 1, nmx
        write (file_tmp,*) amx(iz,:)
    enddo
    close(file_tmp)
    stop
endif

if (any(isnan(ymx))) then 
    print*,'NAN in ymx:w'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_ymx_w.txt',status = 'unknown')
    do iz = 1, nmx
        write (file_tmp,*) ymx(iz)
    enddo
    close(file_tmp)
    stop
endif
#endif

w = ymx

endsubroutine burialcalc_fdm
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine resdisplay(  &
    nz,nspcc,it &
    ,z,frt,omx,rho,o2x,dicx,alkx,ptx,w &
    ,ccx  &
    ,cctflx,ccadv,ccdif,ccdis,ccrain,ccres  &
    ,time,omtflx,omadv,omdif,omdec,omrain,omres,o2tflx,o2dif,o2dec,o2res  &
    ,dictflx,dicdif,dicdec,dicdis,dicres,alktflx,alkdif,alkdec,alkdis,alkres  &
    ,pttflx,ptadv,ptdif,ptrain,ptres,mom,mcc,msed  &
    ,alkrad,ccrad  &
    )   
implicit none 
integer(kind=4),intent(in)::nz,nspcc,it
real(kind=8),dimension(nz),intent(in)::z,frt,omx,rho,o2x,dicx,alkx,ptx,w
real(kind=8),dimension(nz,nspcc),intent(in)::ccx
real(kind=8),dimension(nspcc),intent(in)::cctflx,ccadv,ccdif,ccdis,ccrain,ccres
real(kind=8),intent(in)::time,omtflx,omadv,omdif,omdec,omrain,omres,o2tflx,o2dif,o2dec,o2res  
real(kind=8),intent(in)::dictflx,dicdif,dicdec,dicdis,dicres,alktflx,alkdif,alkdec,alkdis,alkres  
real(kind=8),intent(in)::pttflx,ptadv,ptdif,ptrain,ptres,mom,mcc(nspcc),msed
integer(kind=4) iz,isp
!!!
real(kind=8),intent(in)::alkrad,ccrad(nspcc)

print*, 'time   :',time, maxval(abs(frt - 1d0))
print*,'~~~~ conc ~~~~'
print'(A,5E11.3)', 'z  :',(z(iz),iz=1,nz,nz/5)
print'(A,5E11.3)', 'om :',(omx(iz)*mom/rho(iz)*100d0,iz=1,nz,nz/5)
print'(A,5E11.3)', 'o2 :',(o2x(iz)*1d3,iz=1,nz,nz/5)
print'(A,5E11.3)', 'cc :',(sum(ccx(iz,:)*mcc(:))/rho(iz)*100d0,iz=1,nz,nz/5)
print'(A,5E11.3)', 'dic:',(dicx(iz)*1d3,iz=1,nz,nz/5)
print'(A,5E11.3)', 'alk:',(alkx(iz)*1d3,iz=1,nz,nz/5)
print'(A,5E11.3)', 'sed:',(ptx(iz)*msed/rho(iz)*100d0,iz=1,nz,nz/5)
print*, '   ..... multiple cc species ..... '
do isp=1,nspcc 
    print'(i0.3,":",5E11.3)',isp,(ccx(iz,isp)*mcc(isp)/rho(iz)*100d0,iz=1,nz,nz/5)
enddo
print*,'++++ flx ++++'
print'(8A11)', 'tflx','adv','dif','omrxn','ccrxn','rad','rain','res'
print'(A,8E11.3)', 'om :', omtflx, omadv,  omdif, omdec,0d0,0d0,omrain, omres
print'(A,8E11.3)', 'o2 :',o2tflx,0d0, o2dif,o2dec, 0d0,0d0,0d0,o2res

print'(A,8E11.3)', 'cc :',sum(cctflx),  sum(ccadv), sum(ccdif),0d0,sum(ccdis),sum(ccrad), sum(ccrain), sum(ccres) 
print'(A,8E11.3)', 'dic:',dictflx, 0d0,dicdif, dicdec,  dicdis, 0d0,0d0,dicres 
print'(A,8E11.3)', 'alk:',alktflx, 0d0, alkdif, alkdec, alkdis, alkrad,0d0, alkres 
print'(A,8E11.3)', 'sed:',pttflx, ptadv,ptdif,  0d0, 0d0,0d0, ptrain, ptres

print*, '   ..... multiple cc species ..... '
do isp=1,nspcc 
    print'(i0.3,":",8E11.3)',isp,cctflx(isp), ccadv(isp), ccdif(isp),0d0,ccdis(isp),ccrad(isp), ccrain(isp), ccres(isp) 
enddo

print*,'==== burial etc ===='
print'(A,5E11.3)', 'z  :',(z(iz),iz=1,nz,nz/5)
print'(A,5E11.3)', 'w  :',(w(iz),iz=1,nz,nz/5)
print'(A,5E11.3)', 'rho:',(rho(iz),iz=1,nz,nz/5)
print'(A,5E11.3)', 'frc:',(frt(iz),iz=1,nz,nz/5)

print*,''

endsubroutine resdisplay 
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine sigrec(  &
    nz,file_sigmly,file_sigmlyd,file_sigbtm,w,time,age,izrec,d13c_blk,d13c_blkc &
    ,d13c_blkf,d18o_blk,d18o_blkc,d18o_blkf,ccx,mcc,rho,ptx,msed,izrec2,nspcc  &
    ,d14c_age,capd47   &  
    )
implicit none 
integer(kind=4),intent(in)::nz,file_sigmly,file_sigmlyd,file_sigbtm,izrec,izrec2,nspcc
real(kind=8),dimension(nz),intent(in)::w,age,d13c_blk,d13c_blkc,d13c_blkf,d18o_blk,d18o_blkc,d18o_blkf
real(kind=8),dimension(nz),intent(in)::rho,ptx
real(kind=8),dimension(nz,nspcc),intent(in)::ccx
real(kind=8),intent(in)::time,mcc(nspcc),msed 
!!!!!!
real(kind=8),intent(in)::d14c_age(nz),capd47(nz)  

#ifndef size 
if (all(w>=0d0)) then  ! not recording when burial is negative 
    write(file_sigmly,*) time-age(izrec),d13c_blk(izrec),d18o_blk(izrec) &
        ,sum(ccx(izrec,:)*mcc(:))/rho(izrec)*100d0,ptx(izrec)*msed/rho(izrec)*100d0  &
        ,d14c_age(izrec),capd47(izrec)  
    write(file_sigmlyd,*) time-age(izrec2),d13c_blk(izrec2),d18o_blk(izrec2) &
        ,sum(ccx(izrec2,:)*mcc(:))/rho(izrec2)*100d0,ptx(izrec2)*msed/rho(izrec2)*100d0  &
        ,d14c_age(izrec2),capd47(izrec2)  
    write(file_sigbtm,*) time-age(nz),d13c_blk(nz),d18o_blk(nz) &
        ,sum(ccx(nz,:)*mcc(:))/rho(nz)*100d0,ptx(nz)*msed/rho(nz)*100d0  &
        ,d14c_age(nz),capd47(nz)  
endif  
#else 
if (all(w>=0d0)) then  ! not recording when burial is negative 
    write(file_sigmly,*) time-age(izrec),d13c_blk(izrec),d18o_blk(izrec) &
        ,sum(ccx(izrec,:)*mcc(:))/rho(izrec)*100d0,ptx(izrec)*msed/rho(izrec)*100d0  &
        ,d13c_blkf(izrec),d18o_blkf(izrec),sum(ccx(izrec,1:4)*mcc(1:4))/rho(izrec)*100d0  &
        ,d13c_blkc(izrec),d18o_blkc(izrec),sum(ccx(izrec,5:8)*mcc(5:8))/rho(izrec)*100d0  
    write(file_sigmlyd,*) time-age(izrec2),d13c_blk(izrec2),d18o_blk(izrec2) &
        ,sum(ccx(izrec2,:)*mcc(:))/rho(izrec2)*100d0,ptx(izrec2)*msed/rho(izrec2)*100d0  &
        ,d13c_blkf(izrec2),d18o_blkf(izrec2),sum(ccx(izrec2,1:4)*mcc(1:4))/rho(izrec2)*100d0  &
        ,d13c_blkc(izrec2),d18o_blkc(izrec2),sum(ccx(izrec2,5:8)*mcc(5:8))/rho(izrec2)*100d0  
    write(file_sigbtm,*) time-age(nz),d13c_blk(nz),d18o_blk(nz) &
        ,sum(ccx(nz,:)*mcc(:))/rho(nz)*100d0,ptx(nz)*msed/rho(nz)*100d0 &
        ,d13c_blkf(nz),d18o_blkf(nz),sum(ccx(nz,1:4)*mcc(1:4))/rho(nz)*100d0  &
        ,d13c_blkc(nz),d18o_blkc(nz),sum(ccx(nz,5:8)*mcc(5:8))/rho(nz)*100d0  
endif 
#endif

endsubroutine sigrec
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine closefiles(  &
    file_ptflx,file_ccflx,file_omflx,file_o2flx,file_dicflx,file_alkflx,file_err  &
    ,file_bound,file_totfrac,file_sigmly,file_sigmlyd,file_sigbtm,file_ccflxes,nspcc  &
    )
implicit none 
integer(kind=4),intent(in)::file_ptflx,file_ccflx,file_omflx,file_o2flx,file_dicflx,file_alkflx,file_err,nspcc  
integer(kind=4),intent(in)::file_bound,file_totfrac,file_sigmly,file_sigmlyd,file_sigbtm,file_ccflxes(nspcc)
integer(kind=4) isp

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

endsubroutine closefiles 
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine resrec(  &
    workdir,anoxic,nspcc,labs,turbo2,nobio,co3i,co3sat,mcc,ccx,nz,rho,frt,ccadv,file_tmp,izml,chr,dt,it,time  &
    )
implicit none
integer(kind=4),intent(in)::nspcc,file_tmp,nz,izml
real(kind=8),intent(in)::co3i,co3sat,mcc(nspcc)
real(kind=8),dimension(nz,nspcc),intent(in)::ccx
real(kind=8),dimension(nz),intent(in)::rho,frt
real(kind=8),dimension(nspcc),intent(in)::ccadv
character*255,intent(inout)::workdir
logical,intent(in)::anoxic
logical,dimension(nspcc+2),intent(in)::labs,turbo2,nobio
character*25,intent(in)::chr(3,4)
real(kind=8),intent(in)::dt,time
integer(kind=4),intent(in)::it

workdir = 'C:/Users/YK/Desktop/Sed_res/'
workdir = trim(adjustl(workdir))//'caco3_output/res/'
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

open(unit=file_tmp,file=trim(adjustl(workdir))//'lys_sense_'//    &
    'cc-'//trim(adjustl(chr(1,4)))//'_rr-'//trim(adjustl(chr(2,4)))  &
    //'.res',action='write',status='unknown',access='append') 
write(file_tmp,*) 1d6*(co3i*1d3-co3sat), sum(ccx(1,:)*mcc(:))/rho(1)*100d0, frt(1)  &
    ,sum(ccx(nz,:)*mcc(:))/rho(nz)*100d0, frt(nz),sum(ccx(izml,:)*mcc(:))/rho(izml)*100d0, frt(izml)  &
    ,dt,it,time 
close(file_tmp)

open(unit=file_tmp,file=trim(adjustl(workdir))//'ccbur_sense_'// &
    'cc-'//trim(adjustl(chr(1,4)))//'_rr-'//trim(adjustl(chr(2,4)))  &
    //'.res',action='write',status='unknown',access='append') 
write(file_tmp,*) 1d6*(co3i*1d3-co3sat), 1d6*sum(ccadv),dt,it,time 
close(file_tmp)

endsubroutine resrec
!**************************************************************************************************************************************


!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function d2r(delta,rstd)
implicit none
real(kind=8) d2r,delta,rstd
d2r=(delta*1d-3+1d0)*rstd
endfunction d2r
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function r2d(ratio,rstd)
implicit none
real(kind=8) r2d,ratio,rstd
r2d=(ratio/rstd-1d0)*1d3
endfunction r2d
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

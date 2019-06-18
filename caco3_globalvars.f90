

!**************************************************************************************************************************************
module globalvariables
! Module to define variables 
implicit none
#include <defines.h>
integer(kind=4),parameter :: nz = 100  ! grid number 
#ifdef sense
integer(kind=4),parameter :: nspcc = 12  ! number of CaCO3 species 
#elif defined track2
integer(kind=4),parameter :: nspcc = 42
#elif defined size
integer(kind=4),parameter :: nspcc = 8
#else
integer(kind=4),parameter :: nspcc = 4
#endif
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
real(kind=8) sigocn(nsig),sigocni(nsig),sigocnf(nsig),sigflx(nsig),sigblk(nsig,nz),sigblkf(nsig,nz),sigblkc(nsig,nz)
real(kind=8) flxfrc(nspcc),flxfrc2(nspcc)  ! flux fractions used for assigning flux values to realized isotope input changes 
real(kind=8) :: ccflxi = 10d-6 ! mol (CaCO3) cm-2 yr-1  ! a reference caco3 flux; Emerson and Archer (1990) 
real(kind=8) :: omflx = 12d-6 ! mol cm-2 yr-1       ! a reference om flux; Emerson and Archer (1990)
real(kind=8) :: detflx = 180d-6 ! g cm-2 yr-1  ! a reference detrital flux; MUDS input http://forecast.uchicago.edu/Projects/muds.html
#ifndef mocsy
real(kind=8) :: alki =  2285d0 ! uM  ! a reference ALK; MUDS
real(kind=8) :: dici = 2211d0 ! uM   ! a reference DIC; MUDS 
#else 
real(kind=8) :: alki =  2295d0 ! uM  ! a reference ALK from mocsy
real(kind=8) :: dici = 2154d0 ! uM   ! a reference DIC from mocsy 
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
real(kind=8) :: om2cc = 0.666d0  ! rain ratio of organic matter to calcite
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
real(kind=8) dep, depi,depf ! water depth, i and f denote initial and final values 
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
real(kind=8) :: tol = 1d-8  ! tolerance of error 
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
character*255 workdir, filechr  ! work directory and created file names 
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
logical :: flg_500

interface  ! functions in caco3_therm.f90 
 
    function calceq1(tmp,sal,dep)
    implicit none
    real(kind=8) calceq1,tmp,sal,dep
    real(kind=8) coeff(6)
    real(kind=8) tmp_k,pres
    endfunction calceq1
    
    function calceq2(tmp,sal,dep)
    implicit none
    real(kind=8) calceq2,tmp,sal,dep
    real(kind=8) coeff(6)
    real(kind=8) tmp_k,pres
    endfunction calceq2
    
    function calceqcc(tmp,sal,dep)
    implicit none
    real(kind=8) calceqcc,tmp,sal,dep 
    real(kind=8) tmp_k,pres
    endfunction calceqcc
    
    function calceqag(tmp,sal,dep) 
    implicit none
    real(kind=8) calceqag,tmp,sal,pres
    real(kind=8) tmp_k,dep
    endfunction calceqag
    
endinterface
    
endmodule globalvariables
!**************************************************************************************************************************************


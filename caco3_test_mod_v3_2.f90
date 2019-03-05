program caco3 
! trying a simple diagenesis
! irregular grid
! to compile, type 'gfortran -c caco3_therm.f90 -g -fcheck=all; gfortran -cpp -Dsense caco3_test_mod_v3.f90 caco3_therm.o -lopenblas -g -fcheck=all'
implicit none

integer(kind=4),parameter :: nz = 100
real(kind=8) cc(nz),ccx(nz)  ! mol cm-3 sld
real(kind=8) om(nz),omx(nz)  ! mol cm-3 sld
real(kind=8) pt(nz), ptx(nz) ! # of specific particles/# of total particles ??
real(kind=8) :: ccflx = 10d-6 ! mol (CaCO3) cm-2 yr-1  ! Emerson and Archer (1990) 
real(kind=8) :: omflx = 12d-6 ! mol cm-2 yr-1       ! Emerson and Archer (1990)
real(kind=8) :: detflx = 180d-6 ! g cm-2 yr-1  ! MUDS input http://forecast.uchicago.edu/Projects/muds.html
real(kind=8) :: alki =  2285d0 ! uM
real(kind=8) :: dici = 2211d0 ! uM
real(kind=8) :: o2i = 165d0 ! uM
! real(kind=8) :: komi = 2d0  ! /yr  ! MUDS 
real(kind=8) :: komi = 0.1d0  ! /yr  ! Canfield 1994
! real(kind=8) :: komi = 0.01d0  ! /yr  ! ?? Emerson 1985? who adopted relatively slow decomposition rate 
! real(kind=8) :: kcci = 10.0d0*365.25d0  ! /yr  0.15 to 30 d-1 Emerson and Archer (1990) 0.1 to 10 d-1 in Archer 1991
real(kind=8) :: kcci = 1d0*365.25d0  ! /yr  0.15 to 30 d-1 Emerson and Archer (1990) 0.1 to 10 d-1 in Archer 1991
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
real(kind=8) mvom, mvsed, mvcc  ! molar volumes (cm3 mol-1) mv_i = m_i/rho_i
real(kind=8) keq1, keq2, calceq1,calceq2, co3sat, keqag, dep, zox, zoxx
real(kind=8) dic(nz), dicx(nz), alk(nz), alkx(nz), o2(nz), o2x(nz) !  mol cm-3 porewater 
real(kind=8) dbio(nz),dif_dic(nz), dif_alk(nz), dif_o2(nz)
real(kind=8) co2(nz), hco3(nz), co3(nz), pro(nz)
real(kind=8) co2x(nz), hco3x(nz), co3x(nz), prox(nz),co3i
real(kind=8) poro(nz), rho(nz), frt(nz), sporo(nz), sporoi
real(kind=8) rcc(nz), drcc_dcc(nz), drcc_ddic(nz), drcc_dalk(nz), dco3_ddic(nz), dco3_dalk(nz), drcc_dco3(nz), ddum(nz) 
real(kind=8) dpro_dalk(nz), dpro_ddic(nz)
real(kind=8) kcc(nz) 
real(kind=8) kom(nz), oxco2(nz),anco2(nz) 
real(kind=8) w(nz) , wi, dw(nz)
real(kind=8) z(nz), dz(nz), eta(nz), beta 
real(kind=8) :: ztot = 500d0 ! cm 
integer(kind=4) :: nsp = 3  ! independent chemical variables 
integer(kind=4) :: nmx 
real(kind=8),allocatable :: amx(:,:),ymx(:),emx(:)
integer(kind=4),allocatable :: ipiv(:)
integer(kind=4) info
real(kind=8) error, error2
real(kind=8) :: tol = 1d-6
integer(kind=4) iz, row, col, itr  , it, iiz 
integer(kind=4) :: nt = 1000000
integer(kind=4),parameter :: nrec = 10
integer(kind=4) :: cntrec, itrec
integer(kind=4) up(nz), dwn(nz), cnr(nz) ! advection calc. schemes; up or down wind, or central schemes if 1
! integer(kind=4),parameter :: nt = 1
real(kind=8) :: time, dt = 1d2
real(kind=8) :: rectime(nrec), time_max, dumreal
! real(kind=8) :: dtcc=1d5, dtom = 1d5
real(kind=8) :: omadv, omdec, omdif, omrain, omres, omtflx 
real(kind=8) :: o2tflx, o2dec, o2dif, o2res 
real(kind=8) :: cctflx, ccdis, ccdif, ccadv, ccres 
real(kind=8) :: dictflx, dicdis, dicdif, dicdec, ccrain, dicres 
real(kind=8) :: alktflx, alkdis, alkdif, alkdec, alkres 
real(kind=8) :: pttflx, ptdif, ptadv, ptres, ptrain 
real(kind=8) :: trans(nz,nz), transdbio(nz,nz), translabs(nz,nz)
character*256 workdir
character*10 dumchr(3)
character*25 arg, chr(3,4)
integer(kind=4) dumint(8), idp, izox, narg, ia
integer(kind=4) :: file_tmp=100,file_ccflx=101,file_omflx=102,file_o2flx=103,file_dicflx=104,file_alkflx=105,file_ptflx=106 
external dgesv
logical :: nonlocal = .false. ! use transition matrix from LABS for non-local mixing
logical :: oxic = .true.  ! oxic only model of OM degradation by Emerson (1985) 
logical :: anoxic = .false.  ! oxic-anoxic model of OM degradation by Archer (1991) 
real(kind=8) :: calceqcc, calceqag

call date_and_time(dumchr(1),dumchr(2),dumchr(3),dumint)

!!! get variables !!!
narg = iargc()
do ia = 1, narg,2
    call getarg(ia,arg)
    select case(trim(arg))
        case('cc','CC','Cc')
            call getarg(ia+1,arg)
            read(arg,*)ccflx 
        case('rr','RR','Rr')
            call getarg(ia+1,arg)
            read(arg,*)om2cc
        case('dep','DEP','Dep')
            call getarg(ia+1,arg)
            read(arg,*)dep
        case('dt','DT','Dt')
            call getarg(ia+1,arg)
            read(arg,*)dt
    end select
enddo
print'(3A,3E11.3)','ccflx','om2cc','dep:',ccflx,om2cc, dep
do ia = 1,3
    if (ia==1) dumreal=ccflx
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

!! SWITCHES !!!!

nonlocal = .false.   ! local mixing 
! nonlocal = .true.  ! non-local mixing from LABS

!!!!!!!!!!!!!!
#ifndef nonrec
!! FILES !!!!!!!!!
workdir = 'C:/Users/YK/Desktop/Sed_res/'
workdir = trim(adjustl(workdir))//'test-translabs/profiles/'
workdir = trim(adjustl(workdir))//'cc-'//trim(adjustl(chr(1,4)))//'_rr-'//trim(adjustl(chr(2,4)))  &
    //'_dep-'//trim(adjustl(chr(3,4)))
if (.not. anoxic) then 
    workdir = trim(adjustl(workdir))//'_ox'
else 
    workdir = trim(adjustl(workdir))//'_oxanox'
endif
! workdir = trim(adjustl(workdir))//'-'//trim(adjustl(dumchr(1)))  ! adding date
call system ('mkdir -p '//trim(adjustl(workdir)))
workdir = trim(adjustl(workdir))//'/'

open(unit=file_ptflx,file=trim(adjustl(workdir))//'ptflx.txt',action='write',status='unknown')
open(unit=file_ccflx,file=trim(adjustl(workdir))//'ccflx.txt',action='write',status='unknown')
open(unit=file_omflx,file=trim(adjustl(workdir))//'omflx.txt',action='write',status='unknown')
open(unit=file_o2flx,file=trim(adjustl(workdir))//'o2flx.txt',action='write',status='unknown')
open(unit=file_dicflx,file=trim(adjustl(workdir))//'dicflx.txt',action='write',status='unknown')
open(unit=file_alkflx,file=trim(adjustl(workdir))//'alkflx.txt',action='write',status='unknown')
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
open(unit=file_tmp, file='C:/cygwin64/home/YK/LABS/1dgrid.txt',action='write',status='unknown')
do iz = 1, nz
    write(file_tmp,*) dz(iz)
enddo
close(file_tmp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! FUNDAMENTAL PARAMETERS !!!!!!!!!!!!!

!! Archer's (1991) parameterization 
! ccflx = 18d-6 ! mol cm-2 yr-1
! omflx = 12d-6 ! mol cm-2 yr-1

! om2cc = 0.666d0
! ccflx = 60d-6
! dep = 6d0   ! sea floor depth in km

omflx = om2cc*ccflx
detflx = (1d0/9d0)*ccflx*mcc ! 90% of mass flux becomes inorganic C; g cm-2 yr-1

mvom = mom/rhoom
mvsed = msed/rhosed
mvcc = mcc/rhocc
    

poro = poroi
sporoi = 1d0-poroi
sporo = 1d0 - poro

w = 0.003d0 ! cm/yr 

wi = (detflx/msed*mvsed + ccflx*mvcc +omflx*mvom)/(1d0-poroi)
wi = (detflx/msed*mvsed + ccflx*mvcc            )/(1d0-poroi)
! wi = (detflx/msed*msed + ccflx*mcc + omflx*mom)/(1d0-poroi)/2.5d0

! print*,wi
! pause

300 w= wi

zox = 10d0 ! priori assumed oxic zone 

!!! ~~~~~~~~~~~~~~ set recording time 
#ifndef nonrec
open(unit=file_tmp,file=trim(adjustl(workdir))//'rectime.txt',action='write',status='unknown')
time_max = ztot / wi ! yr
do itrec=1,nrec 
    rectime(itrec)=itrec*time_max/real(nrec)*10d0
    write(file_tmp,*) rectime(itrec)
enddo
close(file_tmp)
! pause
cntrec = 1
#endif
!!!~~~~~~~~~~~~~~~~~

!!!! TRANSITION MATRIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~~~~~~~~~~~~ loading transition matrix from LABS ~~~~~~~~~~~~~~~~~~~~~~~~
if (nonlocal) then
    ! open(unit=file_tmp, file='C:/Users/YK/Desktop/biot-res/test_makematrix-20181206/data/'&
        ! //'transmtx-5000.txt',action='read',status='unknown')
    open(unit=file_tmp, file='C:/Users/YK/Desktop/biot-res/'// &
        'trans-test-100yr-20190108' &
        //'/data/'//&
        'transmtx-219000' &
        //'.txt',action='read',status='unknown')
    do iz = 1, nz
        read(file_tmp,*) translabs(iz,:)
    enddo
    close(file_tmp)
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
    do iiz = 1, nz
        trans(iiz,:) = trans(iiz,:)*dz(:)
    enddo
endif

! stop
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dbio(:) = max(0d0,min(0.15d0,-0.15d0*(z(:)-10.d0)/6d0))  ! 0.15 cm2/yr (Emerson, 1985)

do iz = 1, nz
    if (z(iz) <=zox) then
        dbio(iz) =  0.15d0
        izox = iz
    else
        dbio(iz) =  0d0
    endif
enddo

! dbio = 0.15d0
! dbio = 0d0

transdbio = 0d0
do iz = 1, izox
    if (iz==1) then 
        transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
    elseif (iz==izox) then 
        transdbio(iz,iz) = 0.5d0*(sporo(Iz)*dbio(iz)+sporo(Iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
        transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
    else 
        transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))  &
            + 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
        transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
    endif
enddo

trans = transdbio
if (nonlocal) trans = translabs
! trans = 0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!~~~~~~~~diffusion & reaction~~~~~~~~~~~~~~~~~~~~~~
dif_dic = (151.69d0 + 7.93d0*temp)*poroi*poroi ! cm2/yr at 2 oC (Huelse et al. 2018)
dif_alk = (151.69d0 + 7.93d0*temp)*poroi*poroi ! cm2/yr  (Huelse et al. 2018)
dif_o2 =  (348.62d0 + 14.09d0*temp)*poroi*poroi

kom = komi 
kcc = kcci

keq1 = calceq1(temp,sal,dep)
keq2 = calceq2(temp,sal,dep)

keqcc = calceqcc(temp,sal,dep)
co3sat = keqcc/cai

!!   INITIAL CONDITIONS !!!!!!!!!!!!!!!!!!!    
cc = 1d-8   
! cc = ccflx/100d0/(1d0-poroi)/w 
dic = dici*1d-6/1d3 ! mol/cm3 
alk = alki*1d-6/1d3 ! mol/cm3

call calcspecies(dic,alk,temp,sal,dep,pro,co2,hco3,co3,nz,info)  
    
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

rcc = kcc*ccx*abs(1d0-co3x*1d3/co3sat)**ncc*merge(1d0,0d0,(1d0-co3x*1d3/co3sat)>0d0)

oxco2 = 0d0
anco2 = 0d0

! ~~~ saving initial conditions 
#ifndef nonrec
write(dumchr(1),'(i3.3)') 0

open(unit=file_tmp,file=trim(adjustl(workdir))//'ptx-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace') 
do iz = 1,nz
    write(file_tmp,*) z(iz),pt(iz)*msed/2.5d0*100,0d0,1d0  ,wi
enddo
close(file_tmp)

open(unit=file_tmp,file=trim(adjustl(workdir))//'ccx-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace') 
do iz = 1,nz
    write(file_tmp,*) z(iz),cc(iz)*100d0/2.5d0*100d0, dic(iz)*1d3, alk(iz)*1d3, co3(iz)*1d3-co3sat, rcc(iz),-log10(pro(iz)) 
enddo
close(file_tmp)

open(unit=file_tmp,file=trim(adjustl(workdir))//'o2x-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace') 
do iz = 1,nz
    write(file_tmp,*) z(iz),o2x(iz)*1d3, oxco2(iz), anco2(iz)
enddo
close(file_tmp)

open(unit=file_tmp,file=trim(adjustl(workdir))//'omx-'//trim(adjustl(dumchr(1)))//'.txt',action='write',status='replace') 
do iz = 1,nz
    write(file_tmp,*) z(iz),om(iz)*mom/2.5d0*100d0
enddo
close(file_tmp)
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~
!! START OF TIME INTEGLATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

time = 0d0

do it = 1,nt

! dt = 0.25d0  ! yr
! dt = 1d0  ! yr
! dt = dtom  ! yr
#ifndef nondisp
print*,'it :',it,dt
#endif
dw = 0d0
oxco2 = 0d0
anco2 = 0d0

itr = 0
error = 1d4

! ~~~~~~~~~ OM & O2 iteration wrt zox ~~~~~~~~~~~~~~~
do while (error > tol)

!~~~~~~ OM calculation ~~~~~~~~~~~~~~~~~

if (oxic) then 
    do iz=1,nz
        if (z(iz) <=zox) then
            kom(iz) = komi
            dbio(iz) =  0.15d0
            izox = iz
        else
            kom(iz) = 0d0
            if (anoxic) kom(iz) = komi
            dbio(iz) =  0d0
        endif
    enddo
    
    transdbio = 0d0
    do iz = 1, izox
        if (iz==1) then 
            transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        elseif (iz==izox) then 
            transdbio(iz,iz) = 0.5d0*(sporo(Iz)*dbio(iz)+sporo(Iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
        else 
            transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))  &
                + 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        endif
    enddo
    trans = transdbio
    if (nonlocal) trans = translabs

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
            + (sporo(iz)*w(iz)*1d0-sporoi*wi*0d0)/dz(1)   &
            + sporo(iz)*kom(iz)   &
            ) 
        ! amx(row,row+nsp) =  (&
            ! - (1d0-poro(iz))*((dbio(iz)+dbio(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(1)+dz(2)))-0d0)/dz(1)  &
            ! )
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
            + (sporo(iz)*w(iz)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  &
            + sporo(iz)*kom(iz)   &
            )
        amx(row,row-nsp) = ( &
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! + w(iz)*(-1d0)/dz(iz)  &
            + (sporo(iz)*w(iz)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
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
            + (sporo(iz)*w(iz)*1d0-sporo(iz-1)*w(iz-1)*0d0)/dz(iz)  &
            + sporo(iz)*kom(iz)   &
            )
        ! amx(row,row+nsp) =  (&
            ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(1d0)/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - 0d0)/dz(iz)  & 
            ! )
        amx(row,row-nsp) =  (&
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            ! + w(iz)*(-1d0)/dz(iz)  &
            + (sporo(iz)*w(iz)*0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            )
    endif
    if (nonlocal) then 
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp 
            if (trans(iiz,iz)==0d0) cycle
            amx(row,col) = amx(row,col) -trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))
        enddo
    else 
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp 
            if (trans(iiz,iz)==0d0) cycle
            amx(row,col) = amx(row,col) -trans(iiz,iz)/dz(iz)
        enddo
    endif
enddo

ymx = - ymx

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,info) 

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
        omadv = omadv + (sporo(iz)*w(iz)*omx(iz)-0d0)/dz(iz)*dz(iz)  
            ! + ((1d0-poro(iz))*w(iz)-(1d0-poroi)*wi)*om(iz)/dz(iz)*dz(Iz)  
        ! omdif = omdif &
            ! - (1d0-poro(iz))*((dbio(iz)+dbio(iz+1))*0.5d0*(omx(Iz+1)-omx(iz))/(0.5d0*(dz(1)+dz(2)))-0d0)/dz(1)*dz(iz)
        omdec = omdec + sporo(iz)*kom(iz)*omx(iz)*dz(iz)
        omrain = omrain - omflx/dz(1)*dz(iz)
    else if (iz == nz) then 
        omtflx = omtflx + sporo(iz)*(omx(iz)-om(iz))/dt*dz(iz)
        omadv = omadv + (sporo(iz)*w(iz)*omx(iz)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz) 
            ! + ((1d0-poro(iz))*w(iz)-(1d0-poro(iz-1))*w(Iz-1))*omx(iz)/dz(iz)*dz(Iz)
        ! omdif = omdif &
            ! - (1d0-poro(iz))*(0d0  &
            ! -(dbio(iz)+dbio(iz-1))*0.5d0*(omx(Iz)-omx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)*dz(iz)
        omdec = omdec + sporo(iz)*kom(iz)*omx(iz)*dz(iz)
    else 
        omtflx = omtflx + sporo(iz)*(omx(iz)-om(iz))/dt*dz(iz)
        omadv = omadv + (sporo(iz)*w(iz)*omx(iz)-sporo(iz-1)*w(iz-1)*omx(iz-1))/dz(iz)*dz(iz) 
        ! omdif = omdif &
            ! - (1d0-poro(iz))*((dbio(iz)+dbio(iz+1))*0.5d0*(omx(Iz+1)-omx(iz))/(0.5d0*(dz(iz)+dz(iz+1)))  &
            ! -(dbio(iz)+dbio(iz-1))*0.5d0*(omx(Iz)-omx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)*dz(iz)
        omdec = omdec + sporo(iz)*kom(iz)*omx(iz)*dz(iz)
    endif
    if (nonlocal) then 
        do iiz = 1, nz
            if (trans(iiz,iz)==0d0) cycle
            omdif = omdif -trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))*(1d0-poro(iz))*dz(iz)
        enddo
    else
        do iiz = 1, nz
            if (trans(iiz,iz)==0d0) cycle
            omdif = omdif -trans(iiz,iz)/dz(iz)*dz(iz)
        enddo
    endif
enddo

omres = omadv + omdec + omdif + omrain + omtflx 

if (any(omx<0d0)) then
    print*,'negative om, stop'
    print*,omx
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

if (zox == ztot) then 

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

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,info) 

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
    else if (iz>izox) then  
        ymx(row) = ( 0d0& 
            )
        amx(row,row) = (& 
            + 1d0 &
            )
    endif
enddo

ymx = - ymx

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,info) 

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
 
! print*, 'zox',zox, zoxx
if (zox==zoxx) exit 

error = abs((zox -zoxx)/zox)
 
zox = 0.5d0*(zox + zoxx)
! zox = zoxx

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
        o2x(iz) = o2th
        anco2(iz) = (1d0-poro(iz))*kom(iz)*omx(iz)  ! anaerobic respiration 
    endif
enddo

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

nsp = 3
nmx = nz*nsp
deallocate(amx,ymx,emx,ipiv)
allocate(amx(nmx,nmx),ymx(nmx),emx(nmx),ipiv(nmx))

do while (error > tol)
    
amx = 0d0
ymx = 0d0

rcc = kcc*ccx*abs(1d0-co3x*1d3/co3sat)**ncc*merge(1d0,0d0,(1d0-co3x*1d3/co3sat)>0d0)
drcc_dcc = kcc*abs(1d0-co3x*1d3/co3sat)**ncc*merge(1d0,0d0,(1d0-co3x*1d3/co3sat)>0d0)
drcc_dco3 = kcc*ccx*ncc*abs(1d0-co3x*1d3/co3sat)**(ncc-1d0)*merge(1d0,0d0,(1d0-co3x*1d3/co3sat)>0d0)  &
    *(-1d3/co3sat)
ddum = dicx*dev
call calcspecies(dicx+ddum,alkx,temp,sal,dep,prox,co2x,hco3x,dco3_ddic,nz,info)  
if (info==1) then
    dt = dt/10d0
    go to 300
endif 
call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,info)
if (info==1) then 
    dt=dt/10d0
    go to 300
endif 
dco3_ddic = (dco3_ddic-co3x)/ddum
! print*, dco3_ddic
! stop
ddum = alkx*dev
call calcspecies(dicx,alkx+ddum,temp,sal,dep,prox,co2x,hco3x,dco3_dalk,nz,info)  
if (info==1) then 
    dt=dt/10d0
    go to 300
endif 
call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,info)
if (info==1) then 
    dt=dt/10d0
    go to 300
endif 
dco3_dalk = (dco3_dalk-co3x)/ddum
! print*, dco3_dalk
! stop
drcc_ddic = drcc_dco3*dco3_ddic
drcc_dalk = drcc_dco3*dco3_dalk


do iz = 1,nz 
    row = 1 + (iz-1)*nsp 
    if (iz == 1) then 
        ymx(row) = &
            + sporo(iz)*(ccx(iz)-cc(iz))/dt &
            - ccflx/dz(1) &
            ! - ((dbio(iz)+dbio(iz+1))*0.5d0*(ccx(iz+1)-ccx(iz))/(0.5d0*(dz(1)+dz(2)))  &
            ! - 0d0  &  !  no bioturbation loss at the top boundary  
            ! )/dz(1)  &
            + (sporo(iz)*w(iz)*ccx(iz)-0d0)/dz(1)  &
            + sporo(iz)*rcc(iz)
        amx(row,row) = (&
            + sporo(iz)*(1d0)/dt &
            ! - ((dbio(iz)+dbio(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(1)+dz(2)))-0d0)/dz(1)  &
            + (sporo(iz)*w(iz)*1d0-0d0)/dz(1)   &
            + sporo(iz)* drcc_dcc(iz)  &
            )* ccx(iz) 
        ! amx(row,row+nsp) =  (&
            ! - ((dbio(iz)+dbio(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(1)+dz(2)))-0d0)/dz(1)  &
            ! )*ccx(iz+1)
        amx(row,row+1) = (&
            + sporo(iz)*drcc_ddic(iz)  &
            )*dicx(iz)
        amx(row,row+2) = (&
            + sporo(iz)*drcc_dalk(iz)  &
            )*alkx(iz)
        !  DIC 
        ymx(row+1) = ( &
            + poro(iz)*(dicx(iz)-dic(iz))/dt & 
            - ((poro(iz)*dif_dic(iz)+poro(iz+1)*dif_dic(iz+1))*0.5d0*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_dic(iz)*(dicx(iz)-dici*1d-6/1d3)/dz(iz))/dz(iz)  &
            - oxco2(iz) &
            - anco2(iz) &
            - (1d0-poro(Iz))*rcc(iz)  &
            )*fact
        amx(row+1,row+1) = (& 
            + poro(iz)*(1d0)/dt & 
            - ((poro(iz)*dif_dic(iz)+poro(Iz+1)*dif_dic(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(Iz)*dif_dic(iz)*(1d0)/dz(iz))/dz(iz)&
            - (1d0-poro(Iz))*drcc_ddic(iz)  &
            )*dicx(iz)*fact
        amx(row+1,row+1+nsp) = (& 
            - ((poro(iz)*dif_dic(iz)+poro(Iz+1)*dif_dic(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - 0d0)/dz(iz)&
            )*dicx(iz+1)*fact
        amx(row+1,row+2) = ( &
            - (1d0-poro(Iz))*drcc_dalk(iz)  &
            )*alkx(iz)*fact
        amx(row+1,row) = (&
            - (1d0-poro(Iz))*drcc_dcc(iz)  &
            )*ccx(iz)*fact
        ! ALK
        ymx(row+2) = (& 
            + poro(iz)*(alkx(iz)-alk(iz))/dt & 
            - ((poro(iz)*dif_alk(iz)+poro(Iz+1)*dif_alk(iz+1))*0.5d0*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_alk(iz)*(alkx(iz)-alki*1d-6/1d3)/dz(iz))/dz(iz) &
            - anco2(iz) &
            - 2d0* (1d0-poro(Iz))*rcc(iz)  &
            )*fact
        amx(row+2,row+2) = (& 
            + poro(iz)*(1d0)/dt & 
            - ((poro(iz)*dif_alk(iz)+poro(iz+1)*dif_alk(iz+1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_alk(iz)*(1d0)/dz(iz))/dz(iz)  &
            - 2d0* (1d0-poro(Iz))*drcc_dalk(iz)  &
            )*alkx(iz)*fact
        amx(row+2,row+2+nsp) = (& 
            - ((poro(Iz)*dif_alk(iz)+poro(Iz+1)*dif_alk(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1))) &
            - 0d0)/dz(iz)&
            )*alkx(iz+1)*fact
        amx(row+2,row+1) = (&
            - 2d0* (1d0-poro(Iz))*drcc_ddic(iz)  &
            )*dicx(iz)*fact
        amx(row+2,row) = (&
            - 2d0* (1d0-poro(Iz))*drcc_dcc(iz)  &
            )*ccx(iz)*fact
    else if (iz == nz) then 
        ymx(row) = & 
            + sporo(iz)*(ccx(iz)-cc(iz))/dt &
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ccx(iz)-ccx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            + (sporo(iz)*w(iz)*ccx(iz)-sporo(iz-1)*w(iz-1)*ccx(iz-1))/dz(iz)  &
            + sporo(iz)*rcc(iz)
        amx(row,row) = (&
            + sporo(iz)*(1d0)/dt &
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            + (sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  &
            + sporo(iz)*drcc_dcc(iz)   &
            )*ccx(iz)
        amx(row,row-nsp) = ( &
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            + (0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            )*ccx(iz-1)
        amx(row,row+1) = (&
            + sporo(iz)*drcc_ddic(iz) &
            )*dicx(iz)
        amx(row,row+2) = (&
            + sporo(iz)*drcc_dalk(iz) &
            )*alkx(iz)
        ! ymx(row) = ccx(iz) - ccx(iz-1)
        ! amx(row, row) = 1d0*ccx(iz)
        ! amx(row,row-1)= -1d0*ccx(iz-1)
        
        ! DIC
        ymx(row+1) = (& 
            + poro(iz)*(dicx(iz)-dic(iz))/dt &
            - (0d0 - 0.5d0*(poro(iz)*dif_dic(iz)+poro(Iz-1)*dif_dic(Iz-1))*(dicx(iz)-dicx(iz-1))  &
                /(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            - oxco2(iz) &
            - anco2(iz) &
            - sporo(iz)*rcc(iz)  &
            )*fact
        amx(row+1,row+1) = ( & 
            + poro(iz)*(1d0)/dt &
            - (0d0 - 0.5d0*(poro(iz)*dif_dic(iz)+poro(Iz-1)*dif_dic(Iz-1))*(1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            - sporo(Iz)*drcc_ddic(iz)  &
            )*dicx(iz)*fact
        amx(row+1,row+1-nsp) = ( & 
            - (0d0 - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(Iz-1))*(-1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            ) * dicx(iz-1)*fact
        amx(row+1,row+2) = (&
            - sporo(Iz)*drcc_dalk(iz)  &
            )*alkx(iz)*fact
        amx(row+1,row) = (&
            - sporo(Iz)*drcc_dcc(iz)  &
            )*ccx(iz)*fact
        ! ALK 
        ymx(row+2) = ( & 
            + poro(iz)*(alkx(iz)-alk(iz))/dt &
            - (0d0 - 0.5d0*(poro(iz)*dif_alk(iz)+poro(Iz-1)*dif_alk(Iz-1))*(alkx(iz)-alkx(iz-1))/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            - anco2(iz) &
            - 2d0*sporo(Iz)*rcc(iz)  &
            )*fact
        amx(row+2,row+2) = ( & 
            + poro(iz)*(1d0)/dt &
            - (0d0 - 0.5d0*(poro(Iz)*dif_alk(iz)+poro(iz-1)*dif_alk(Iz-1))*(1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            - 2d0*sporo(Iz)*drcc_dalk(iz)  &
            )*alkx(iz)*fact
        amx(row+2,row+2-nsp) = ( & 
            - (0d0 - 0.5d0*(poro(iz)*dif_alk(iz)+poro(Iz-1)*dif_alk(Iz-1))*(-1d0)/(0.5d0*(dz(iz-1)+dz(iz))))/dz(Iz) &
            ) * alkx(iz-1)*fact
        amx(row+2,row+1) = (&
            - 2d0*sporo(Iz)*drcc_ddic(iz)  &
            )*dicx(Iz)*fact
        amx(row+2,row) = (&
            - 2d0*sporo(Iz)*drcc_dcc(iz)  &
            )*ccx(Iz)*fact
    else 
        ymx(row) = & 
            + sporo(iz)*(ccx(iz)-cc(iz))/dt &
            ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(ccx(iz+1)-ccx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ccx(iz)-ccx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            + (sporo(iz)*w(iz)*ccx(iz)-sporo(Iz-1)*w(iz-1)*ccx(iz-1))/dz(iz)  &
            + sporo(iz)*rcc(iz)
        amx(row,row) = (&
            + sporo(iz)*(1d0)/dt &
            ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(-1d0)/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            + (sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  &
            + sporo(iz)*drcc_dcc(iz)  &
            )*ccx(iz)
        ! amx(row,row+nsp) =  (&
            ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(1d0)/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - 0d0)/dz(iz)  & 
            ! )*ccx(iz+1)
        amx(row,row-nsp) =  (&
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            + (0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            )*ccx(iz-1)
        amx(row,row+1) = (& 
            + sporo(Iz)*drcc_ddic(iz)  &
            )*dicx(Iz)
        amx(row,row+2) = (&
            + sporo(Iz)*drcc_dalk(iz) &
            )*alkx(iz)
        ! DIC 
        ymx(row+1) = ( & 
            + poro(iz)*(dicx(iz)-dic(iz))/dt & 
            - (0.5d0*(poro(iz+1)*dif_dic(iz+1)+poro(Iz)*dif_dic(iz))*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(dicx(Iz)-dicx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            - oxco2(iz) &
            - anco2(iz) &
            - sporo(Iz)*rcc(iz)  &
            )*fact
        amx(row+1,row+1) = (& 
            + poro(iz)*(1d0)/dt & 
            - (0.5d0*(poro(iz+1)*dif_dic(iz+1)+poro(iz)*dif_dic(iz))*(-1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            - sporo(iz)*drcc_ddic(iz)  &
            )*dicx(iz)*fact
        amx(row+1,row+1+nsp) = (& 
            - (0.5d0*(poro(iz+1)*dif_dic(iz+1)+poro(iz)*dif_dic(iz))*(1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0d0)/dz(iz)  &
            )*dicx(iz+1)*fact
        amx(row+1,row+1-nsp) = (& 
            - (0d0 &
            - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) &
            )*dicx(iz-1)*fact
        amx(row+1,row+2) = (&
            - sporo(Iz)*drcc_dalk(iz)  &
            )*alkx(iz)*fact
        amx(row+1,row) = (&
            - sporo(Iz)*drcc_dcc(iz)  &
            )*ccx(iz)*fact
        ! ALK 
        ymx(row+2) = (& 
            + poro(iz)*(alkx(iz)-alk(iz))/dt & 
            - (0.5d0*(poro(iz+1)*dif_alk(iz+1)+poro(iz)*dif_alk(iz))*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(Iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(alkx(iz)-alkx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) &
            - anco2(iz) &
            - 2d0*sporo(Iz)*rcc(iz)  &
            ) *fact
        amx(row+2,row+2) = (& 
            + poro(iz)*(1d0)/dt & 
            - (0.5d0*(poro(iz+1)*dif_alk(iz+1)+poro(iz)*dif_alk(iz))*(-1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(Iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  &
            - 2d0*sporo(Iz)*drcc_dalk(iz)  &
            )*alkx(iz)*fact
        amx(row+2,row+2+nsp) = ( & 
            - (0.5d0*(poro(iz+1)*dif_alk(iz+1)+poro(iz)*dif_alk(iz))*(1d0)/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0d0)/dz(iz)  &
            )*alkx(iz+1)*fact
        amx(row+2,row+2-nsp) = (& 
            - (0d0 &
            - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz) &
            )*alkx(iz-1)*fact
        amx(row+2,row+1) = (&
            - 2d0*sporo(Iz)*drcc_ddic(iz)  &
            )*dicx(iz)*fact 
        amx(row+2,row) = (&
            - 2d0*sporo(Iz)*drcc_dcc(iz)  &
            )*ccx(iz)*fact 
    endif
    if (nonlocal) then
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp
            if (trans(iiz,iz)==0d0) cycle
            amx(row,col) = amx(row,col) -trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))*ccx(iiz)
            ymx(row) = ymx(row) - trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))*ccx(iiz)
        enddo
    else
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp
            if (trans(iiz,iz)==0d0) cycle
            amx(row,col) = amx(row,col) -trans(iiz,iz)/dz(iz)*ccx(iiz)
            ymx(row) = ymx(row) - trans(iiz,iz)/dz(iz)*ccx(iiz)
        enddo
    endif
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

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,info) 

if (any(isnan(amx))) then
    print*,'NAN in amx'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_amx.txt',status = 'unknown')
    do iz = 1, nmx
        write (file_tmp,*) amx(iz,:)
    enddo
    close(file_tmp)
    ! stop
endif

if (any(isnan(ymx))) then 
    print*,'NAN in ymx'
    open(unit=file_tmp,file=trim(adjustl(workdir))//'chk_ymx.txt',status = 'unknown')
    do iz = 1, nmx
        write (file_tmp,*) ymx(iz)
    enddo
    close(file_tmp)
    ! stop
endif
#endif
! stop

do iz = 1, nz 
    row = 1+(iz-1)*nsp
    if (ymx(row)>10d0) then 
        ccx(iz) = ccx(iz)*1.5d0
    elseif (ymx(row)<-10d0) then 
        ccx(iz) = ccx(iz)*0.5d0
    else
        ccx(iz) = ccx(iz)*exp(ymx(row))
    endif
    if (ymx(row+1)>10d0) then 
        dicx(iz)=dicx(iz)*1.5d0
    elseif (ymx(row+1)<-10d0) then 
        dicx(iz)=dicx(iz)*0.5d0
    else 
        dicx(iz) = dicx(iz)*exp(ymx(row+1))
    endif
    if (ymx(row+2)>10d0) then 
        alkx(Iz) = alkx(iz)*1.5d0
    elseif (ymx(row+2)<-10d0) then 
        alkx(iz) = alkx(iz)*0.5d0
    else 
        alkx(iz) = alkx(iz)*exp(ymx(row+2))
    endif
    if (ccx(iz)<1d-30) then
        ccx(iz)=1d-30
        ymx(row) = 0d0
    endif
    if (dicx(iz)<1d-100) ymx(row+1) = 0d0
    if (alkx(iz)<1d-100) ymx(row+2) = 0d0
    ! print *,exp(ymx(row))
enddo

error = maxval(exp(abs(ymx))) - 1d0
itr = itr + 1

! print*,'co2 iteration',itr,error,info
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

call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,info)
if (info==1) then 
    dt=dt/10d0
    go to 300
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
        cctflx = cctflx + (1d0-poro(iz))*(ccx(iz)-cc(iz))/dt *dz(iz)
        ccdis = ccdis  + (1d0-poro(Iz))*rcc(iz) *dz(iz)
        ! ccdif = ccdif &
            ! - (1d0-poro(iz))*((dbio(iz)+dbio(iz+1))*0.5d0*(ccx(iz+1)-ccx(iz))/(0.5d0*(dz(1)+dz(2)))  &
            ! - 0d0 )/dz(1)* dz(iz)
        ccrain = ccrain - ccflx/dz(1)*dz(iz)
        ccadv = ccadv + (sporo(iz)*w(iz)*ccx(iz)-0d0)/dz(1) * dz(iz)
        ! print*, iz,(1d0-poro(iz))*(ccx(iz)-cc(iz))/dt *dz(iz) &
            ! + (1d0-poro(Iz))*rcc(iz) *dz(iz) &
            ! - (1d0-poro(iz))*((dbio(iz)+dbio(iz+1))*0.5d0*(ccx(iz+1)-ccx(iz))/(0.5d0*(dz(1)+dz(2)))  &
            ! - 0d0 )/dz(1)* dz(iz) &
            ! - ccflx/(40d0+60d0)/dz(1)*dz(iz) &
            ! + (1d0-poro(iz))*w(iz)*(ccx(iz)-0d0)/dz(1) * dz(iz)
        !  DIC 
        dictflx = dictflx +(dicx(iz)-dic(iz))/dt*dz(iz)*poro(iz) 
        dicdif = dicdif - ((poro(iz)*dif_dic(iz)+poro(iz+1)*dif_dic(iz+1))*0.5d0*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_dic(iz)*(dicx(iz)-dici*1d-6/1d3)/dz(iz))/dz(iz)*dz(iz)
        dicdec = dicdec - oxco2(iz)*dz(iz) - anco2(iz)*dz(iz) 
        dicdis = dicdis - rcc(iz)*sporo(iz)*dz(iz) 
        ! ALK
        alktflx = alktflx + (alkx(iz)-alk(iz))/dt*dz(iz)*poro(iz)
        alkdif = alkdif - ((poro(iz)*dif_alk(iz)+poro(iz+1)*dif_alk(iz+1))*0.5d0*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz)+dz(iz+1))) &
            - poro(iz)*dif_alk(iz)*(alkx(iz)-alki*1d-6/1d3)/dz(iz))/dz(iz)*dz(iz)
        alkdec = alkdec - anco2(iz)*dz(iz) 
        alkdis = alkdis - 2d0* sporo(Iz)*rcc(iz)*dz(iz) 
    else if (iz == nz) then 
        cctflx = cctflx + sporo(iz)*(ccx(iz)-cc(iz))/dt *dz(iz)
        ccdis = ccdis  + sporo(Iz)*rcc(iz) *dz(iz)
        ! ccdif = ccdif &
            ! - (1d0-poro(iz))*(0d0  &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ccx(iz)-ccx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))) )/dz(iz)* dz(iz)
        ccadv = ccadv + (sporo(iz)*w(iz)*ccx(iz)-sporo(iz-1)*w(iz-1)*ccx(iz-1))/dz(iz) * dz(iz)
        ! print*,iz, (1d0-poro(iz))*(ccx(iz)-cc(iz))/dt *dz(iz) &
            ! + (1d0-poro(Iz))*rcc(iz) *dz(iz) &
            ! - (1d0-poro(iz))*(0d0  &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ccx(iz)-ccx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))) )/dz(iz)* dz(iz) &
            ! + (1d0-poro(iz))*w(iz)*(ccx(iz)-ccx(iz-1))/dz(iz) * dz(iz)
        ! DIC
        dictflx = dictflx +(dicx(iz)-dic(iz))/dt*dz(iz)*poro(iz) 
        dicdif = dicdif - (0d0 &
            - 0.5d0*(poro(iz)*dif_dic(iz)+poro(iz-1)*dif_dic(Iz-1))*(dicx(iz)-dicx(iz-1))/(0.5d0*(dz(iz-1)+dz(iz))) &
            )/dz(iz)*dz(iz)
        dicdec = dicdec - oxco2(iz)*dz(iz) - anco2(iz)*dz(iz) 
        dicdis = dicdis - sporo(Iz)*rcc(iz)*dz(iz) 
        ! ALK 
        alktflx = alktflx + (alkx(iz)-alk(iz))/dt*dz(iz)*poro(iz)
        alkdif = alkdif - (0d0 &
            - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(Iz-1))*(alkx(iz)-alkx(iz-1))/(0.5d0*(dz(iz-1)+dz(iz))))/dz(iz)*dz(iz)
        alkdec = alkdec - anco2(iz)*dz(iz)
        alkdis = alkdis - 2d0* Sporo(Iz)*rcc(iz)*dz(iz)
    else 
        cctflx = cctflx + sporo(iz)*(ccx(iz)-cc(iz))/dt *dz(iz)
        ccdis = ccdis  + sporo(Iz)*rcc(iz) *dz(iz)
        ! ccdif = ccdif &
            ! - (1d0-poro(iz))*((dbio(iz+1)+dbio(iz))*0.5d0*(ccx(iz+1)-ccx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ccx(iz)-ccx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))) )/dz(iz)* dz(iz)
        ccadv = ccadv + (sporo(iz)*w(iz)*ccx(iz)-sporo(iz-1)*w(iz-1)*ccx(iz-1))/dz(iz) * dz(iz)
        ! print*,iz,(1d0-poro(iz))*(ccx(iz)-cc(iz))/dt *dz(iz) &
            ! + (1d0-poro(Iz))*rcc(iz) *dz(iz)  &
            ! - (1d0-poro(iz))*((dbio(iz+1)+dbio(iz))*0.5d0*(ccx(iz+1)-ccx(iz))/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(ccx(iz)-ccx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))) )/dz(iz)* dz(iz) &
            ! + (1d0-poro(iz))*w(iz)*(ccx(iz)-ccx(iz-1))/dz(iz) * dz(iz)
        ! DIC 
        dictflx = dictflx +(dicx(iz)-dic(iz))/dt*dz(iz)*poro(iz) 
        dicdif = dicdif - (0.5d0*(poro(iz+1)*dif_dic(iz+1)+poro(iz)*dif_dic(iz))*(dicx(iz+1)-dicx(iz))/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(Iz)*dif_dic(iz)+poro(iz-1)*dif_dic(iz-1))*(dicx(Iz)-dicx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))) &
            )/dz(iz)*dz(iz)
        dicdec = dicdec - oxco2(iz)*dz(iz) - anco2(iz)*dz(iz) 
        dicdis = dicdis - sporo(Iz)*rcc(iz)*dz(iz) 
        ! ALK 
        alktflx = alktflx + (alkx(iz)-alk(iz))/dt*dz(iz)*poro(iz)
        alkdif = alkdif - (0.5d0*(poro(iz+1)*dif_alk(iz+1)+poro(iz)*dif_alk(iz))*(alkx(iz+1)-alkx(iz))/(0.5d0*(dz(iz+1)+dz(Iz))) &
            - 0.5d0*(poro(iz)*dif_alk(iz)+poro(iz-1)*dif_alk(iz-1))*(alkx(iz)-alkx(iz-1))/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)*dz(iz)
        alkdec = alkdec - anco2(iz)*dz(iz)
        alkdis = alkdis - 2d0* sporo(iz)*rcc(iz)*dz(iz) 
    endif
    if (nonlocal) then 
        do iiz = 1, nz
            if (trans(iiz,iz)==0d0) cycle
            ccdif = ccdif -trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))*(1d0-poro(iz))*dz(iz)*ccx(iiz)
        enddo
    else 
        do iiz = 1, nz
            if (trans(iiz,iz)==0d0) cycle
            ccdif = ccdif -trans(iiz,iz)/dz(iz)*dz(iz)*ccx(iiz)
        enddo
    endif
    dw(iz) = dw(iz) -(1d0-poro(iz))*mvcc*rcc(iz)
enddo

ccres = cctflx +  ccdis +  ccdif + ccadv + ccrain
dicres = dictflx + dicdis + dicdif + dicdec 
alkres = alktflx + alkdis + alkdif + alkdec 

if (it==1) then 
    write(file_ccflx,*) 'time, cctflx, ccdis, ccdif, ccadv, ccrain, ccres' 
    write(file_dicflx,*) 'time, dictflx, dicdis, dicdif, dicdec,  dicres' 
    write(file_alkflx,*) 'time, alktflx, alkdis, alkdif, alkdec, alkres' 
endif
write(file_ccflx,*) time,cctflx, ccdis, ccdif, ccadv, ccrain, ccres 
write(file_dicflx,*) time,dictflx, dicdis, dicdif, dicdec,  dicres 
write(file_alkflx,*) time,alktflx, alkdis, alkdif, alkdec, alkres 

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
            + (sporo(iz)*w(iz)*1d0-0d0)/dz(iz)   &
            )            
        ! amx(row,row+nsp) =  (&
            ! - ((dbio(iz)+dbio(iz+1))*0.5d0*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))-0d0)/dz(iz)  &
            ! )
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
            + (sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  &
            )
        amx(row,row-nsp) = ( &
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            + (0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
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
            + (sporo(iz)*w(iz)*1d0-0d0)/dz(iz)  &
            )
        ! amx(row,row+nsp) =  (&
            ! - ((dbio(iz+1)+dbio(iz))*0.5d0*(1d0)/(0.5d0*(dz(iz+1)+dz(iz))) &
            ! - 0d0)/dz(iz)  & 
            ! )
        amx(row,row-nsp) =  (&
            ! - (0d0 &
            ! - (dbio(iz)+dbio(iz-1))*0.5d0*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1))))/dz(iz)  & 
            + (0d0-sporo(iz-1)*w(iz-1)*1d0)/dz(iz)  &
            )
    endif
    if (nonlocal) then 
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp
            if (trans(iiz,iz)==0d0) cycle
            amx(row,col) = amx(row,col) -trans(iiz,iz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))
        enddo
    else 
        do iiz = 1, nz
            col = 1 + (iiz-1)*nsp
            if (trans(iiz,iz)==0d0) cycle
            amx(row,col) = amx(row,col) -trans(iiz,iz)/dz(iz)
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

call dgesv(nmx,int(1),amx,nmx,ipiv,ymx,nmx,info) 

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
            + (sporo(iz)*w(iz)*ptx(iz)-0d0)/dz(iz)*dz(iz)
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
            + (sporo(iz)*w(iz)*ptx(iz)-sporo(iz-1)*w(iz-1)*ptx(iz-1))/dz(iz)*dz(iz)
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
            + (sporo(iz)*w(iz)*ptx(iz)-sporo(iz-1)*w(Iz-1)*ptx(iz-1))/dz(iz)*dz(iz)
    endif
    if(nonlocal) then 
        do iiz = 1, nz
            ! if (iiz==iz) cycle
            if (trans(iiz,iz)==0d0) cycle
            ptdif = ptdif -trans(iiz,iz)*ptx(iiz)/dz(iz)*(1d0-poro(iz))*dz(iiz)/(1d0-poro(iiz))    &
                *(1d0-poro(iz))*dz(iz)
        enddo
    else 
        do iiz = 1, nz
            ! if (iiz==iz) cycle
            if (trans(iiz,iz)==0d0) cycle
            ptdif = ptdif -trans(iiz,iz)*ptx(iiz)/dz(iz)    &
                *dz(iz)
        enddo
    endif
enddo

ptres = pttflx + ptdif + ptadv + ptrain

if (it==1) write(file_ptflx,*) 'time, pttflx, ptdif, ptadv, ptrain, ptres'
write(file_ptflx,*) time, pttflx, ptdif, ptadv, ptrain, ptres

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

rho = omx*mom + ptx*msed +  ccx*mcc
frt = omx*mvom + ptx*mvsed + ccx*mvcc

wi = (detflx/msed*mvsed + ccflx*mvcc +omflx*mvom + sum(dw*dz))/(1d0-poroi)
wi = (detflx/msed*mvsed + ccflx*mvcc +omflx*mvom             )/(1d0-poroi)

do iz=1,nz
    if (iz==1) then 
        w(iz)=((1d0-poroi)*wi + dw(iz)*dz(iz))/(1d0-poro(iz))
    else 
        w(iz)=((1d0-poro(iz-1))*w(iz-1) + dw(iz)*dz(iz))/(1d0-poro(iz))
    endif
    
    ! if (w(iz) < 0d0) w(iz) = wi*1d-3
    
enddo

if (any(w<0d0)) then
    print*,'negative burial velocity'
    print*, w
    stop
endif 

if (any(rho<0d0)) then
    print*,'negative density'
    print*, rho
    stop
endif 

#ifndef nonrec
if (time>=rectime(cntrec)) then 
    write(dumchr(1),'(i3.3)') cntrec 
    
    open(unit=file_tmp,file=trim(adjustl(workdir))//'ptx-'//trim(adjustl(dumchr(1)))//'.txt' &
        ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),ptx(iz)*msed/rho(iz)*100d0,rho(iz),frt(iz)  ,w(iz)
    enddo
    close(file_tmp)

    open(unit=file_tmp,file=trim(adjustl(workdir))//'ccx-'//trim(adjustl(dumchr(1)))//'.txt' &
        ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),ccx(iz)*mcc/rho(iz)*100d0, dicx(iz)*1d3, alkx(iz)*1d3  &
            , co3x(iz)*1d3-co3sat, rcc(iz),-log10(prox(iz)) 
    enddo
    close(file_tmp)

    open(unit=file_tmp,file=trim(adjustl(workdir))//'omx-'//trim(adjustl(dumchr(1)))//'.txt'  &
        ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),omx(iz)*mom/rho(iz)*100d0
    enddo
    close(file_tmp)

    open(unit=file_tmp,file=trim(adjustl(workdir))//'o2x-'//trim(adjustl(dumchr(1)))//'.txt'  &
        ,action='write',status='replace') 
    do iz = 1,nz
        write(file_tmp,*) z(iz),o2x(iz)*1d3, oxco2(iz), anco2(iz)
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
print'(A,5E11.3)', 'cc :',(ccx(iz)*mcc/rho(iz)*100d0,iz=1,nz,nz/5)
print'(A,5E11.3)', 'dic:',(dicx(iz)*1d3,iz=1,nz,nz/5)
print'(A,5E11.3)', 'alk:',(alkx(iz)*1d3,iz=1,nz,nz/5)
print'(A,5E11.3)', 'sed:',(ptx(iz)*msed/rho(iz)*100d0,iz=1,nz,nz/5)
print*,'++++ flx ++++'
print'(7A11)', 'tflx','adv','dif','omrxn','ccrxn','rain','res'
print'(A,7E11.3)', 'om :', omtflx, omadv,  omdif, omdec,0d0,omrain, omres
print'(A,7E11.3)', 'o2 :',o2tflx,0d0, o2dif,o2dec, 0d0,0d0,o2res

print'(A,7E11.3)', 'cc :',cctflx,  ccadv, ccdif,0d0,ccdis, ccrain, ccres 
print'(A,7E11.3)', 'dic:',dictflx, 0d0,dicdif, dicdec,  dicdis, 0d0,dicres 
print'(A,7E11.3)', 'alk:',alktflx, 0d0, alkdif, alkdec, alkdis, 0d0, alkres 
print'(A,7E11.3)', 'sed:',pttflx, ptadv,ptdif,  0d0, 0d0, ptrain, ptres

print*,'==== burial etc ===='
print'(A,5E11.3)', 'z  :',(z(iz),iz=1,nz,nz/5)
print'(A,5E11.3)', 'w  :',(w(iz),iz=1,nz,nz/5)
print'(A,5E11.3)', 'rho:',(rho(iz),iz=1,nz,nz/5)
print'(A,5E11.3)', 'frc:',(frt(iz),iz=1,nz,nz/5)

print*,''
#endif

!! in theory, o2dec/ox2om + alkdec = dicdec = omdec (in absolute value)
if (dicdec /= 0d0) then 
    if ( abs((o2dec/ox2om - alkdec + dicdec)/dicdec) > tol) then 
        print*, abs((o2dec/ox2om + alkdec - dicdec)/dicdec) 
        stop
    endif
endif 

! pause

time = time + dt

o2 = o2x
om = omx

cc = ccx
dic = dicx
alk = alkx

pt = ptx

enddo

#ifndef nonrec
close(file_ptflx)
close(file_ccflx)
close(file_omflx)
close(file_o2flx)
close(file_dicflx)
close(file_alkflx)
#endif

workdir = 'C:/Users/YK/Desktop/Sed_res/'
workdir = trim(adjustl(workdir))//'test-translabs/res/'

call system ('mkdir -p '//trim(adjustl(workdir)))


if (.not. anoxic) then 
    workdir = trim(adjustl(workdir))//'ox-'
else 
    workdir = trim(adjustl(workdir))//'oxanox-'
endif

open(unit=file_tmp,file=trim(adjustl(workdir))//'lys_sense_'// &
    'cc-'//trim(adjustl(chr(1,4)))//'_rr-'//trim(adjustl(chr(2,4)))  &
    //'.txt',action='write',status='unknown',access='append') 
write(file_tmp,*) 1d6*(co3i*1d3-co3sat), ccx(1)*mcc/rho(1)*100d0, frt(1)  &
    ,ccx(nz)*mcc/rho(nz)*100d0, frt(nz)
close(file_tmp)

end program 
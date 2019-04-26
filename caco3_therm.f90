! functions and subroutines to calculate caco3 thermodynamics

function calceq1(tmp,sal,dep) ! Millero et al. (2006, MC); Pressure effect from Zeebe and Wolf-Gladrow (2001) 
implicit none
real(kind=8) calceq1,tmp,sal,dep
real(kind=8) coeff(6)
real(kind=8) tmp_k,pres

tmp_k=tmp+273.15d0
pres = dep*100d0 ! bar

coeff(1)=13.4191D0
coeff(2)=0.0331D0
coeff(3)=-5.33D-5
coeff(4)=-530.1228D0
coeff(5)=-6.103d0
coeff(6)=-2.06950d0

calceq1=-126.34048d0+6320.813d0/tmp_k+19.568224d0*log(tmp_k)
calceq1=calceq1+coeff(1)*sal**0.5d0+coeff(2)*sal+coeff(3)*sal**2d0 &
    +(coeff(4)*sal**0.5d0+coeff(5)*sal)/tmp_k   &
    +coeff(6)*sal**0.5d0*log(tmp_k)
calceq1=10d0**(-calceq1)

calceq1=calceq1*exp((-(-15.82d0-0.0219d0*tmp)*pres+0.5d0*(1.13d-3-0.1475d-3*tmp)*pres*pres)/83.131d0/tmp_k)

endfunction calceq1

function calceq2(tmp,sal,dep) ! Millero et al. (2006, MC); Pressure effect from Zeebe and Wolf-Gladrow (2001) 
implicit none
real(kind=8) calceq2,tmp,sal,dep ! dep in km
real(kind=8) coeff(6)
real(kind=8) tmp_k,pres

tmp_k=tmp+273.15d0
pres = dep*100d0 ! bar

coeff(1)=21.0894D0
coeff(2)=0.1248D0
coeff(3)=-0.0003687D0
coeff(4)=-772.483d0
coeff(5)=-20.051D0
coeff(6)=-3.32254d0

calceq2=-90.18333d0+5143.692d0/tmp_k+14.613358d0*log(tmp_k)
calceq2=calceq2+coeff(1)*sal**0.5d0+coeff(2)*sal+coeff(3)*sal**2d0 &
    +(coeff(4)*sal**0.5d0+coeff(5)*sal)/tmp_k   &
    +coeff(6)*sal**0.5d0*log(tmp_k)
calceq2=10d0**(-calceq2)

calceq2=calceq2*exp((-(-25.50d0+0.1271d0*tmp)*pres+0.5d0*(-3.08d-3+0.0877d-3*tmp)*pres*pres)/83.131d0/tmp_k)

endfunction calceq2

function calceqcc(tmp,sal,dep) ! Mucci (1983) cited by Zeebe and Wolf-Gladrow (2001)
implicit none
real(kind=8) calceqcc,tmp,sal,dep ! depth in km
real(kind=8) tmp_k,pres

tmp_k=tmp+273.15d0
pres = dep*100d0 ! bar

calceqcc = -171.9065d0 - 0.077993d0*tmp_k + 2839.319d0/tmp_k &
    +71.595d0*log10(tmp_k) &
    +(-0.77712d0+0.0028426d0*tmp_k+178.34d0/tmp_k)*sal**0.5d0 &
    -0.07711d0*sal+0.0041249d0*sal**1.5d0 
calceqcc = 10d0**calceqcc

calceqcc=calceqcc*exp((-(-48.76d0+0.5304d0*tmp)*pres+0.5d0*(-11.76d-3+0.3692d-3*tmp)*pres*pres)/83.131d0/tmp_k)

endfunction calceqcc

function calceqag(tmp,sal,dep) ! Mucci (1983) cited by Zeebe and Wolf-Gladrow (2001)
implicit none
real(kind=8) calceqag,tmp,sal,pres
real(kind=8) tmp_k,dep

tmp_k=tmp+273.15d0
pres = dep*100d0 ! bar

calceqag = -171.945d0 - 0.077993d0*tmp_k + 2903.293d0/tmp_k &
    +71.595d0*log10(tmp_k) &
    +(-0.068393d0+0.0017276d0*tmp_k+88.135d0/tmp_k)*sal**0.5d0 &
    -0.10018d0*sal+0.0059415d0*sal**1.5d0 
calceqag = 10d0**calceqag

calceqag=calceqag*exp((-(-46.00d0+0.5304d0*tmp)*pres+0.5d0*(-11.76d-3+0.3692d-3*tmp)*pres*pres)/83.131d0/tmp_k)

endfunction calceqag

subroutine calcspecies(dic,alk,tmp,sal,dep,ph,co2,hco3,co3,nz,info) ! returning the same units as inputs (?)
implicit none
integer(kind=4),intent(in) :: nz
integer(kind=4),intent(out)::info
real(kind=8),intent(in)::dic(nz),alk(nz),tmp,sal,dep
real(kind=8),intent(out)::ph(nz),co2(nz),hco3(nz),co3(nz)
real(kind=8)::a(nz),b(nz),c(nz)
real(kind=8)::calceq1,calceq2,k1,k2

info=0
k1=calceq1(tmp,sal,dep)
k2=calceq2(tmp,sal,dep)
a=1d0
b=(1d0-dic/alk)*k1
c=(1d0-2d0*dic/alk)*k1*k2
ph = (-b+(b*b-4d0*a*c)**0.5d0)/2d0/a
if (any(ph<0d0)) then
    print*,'... unsable to calculate ph'
    ! print*,dic
    ! print*,alk
    ! print*,ph
    ! stop
    info=1
    return
endif
co2 = alk/(k1/ph+2d0*k1*k2/ph/ph)
! hco3=co2*k1/ph
hco3=alk/(1d0+2d0*k2/ph)
! co3=hco3*k2/ph
co3=alk/(ph/k2+2d0)
! ph=-log10(ph)
endsubroutine calcspecies

subroutine calcdevs(dic,alk,tmp,sal,dep,nz,info,dco3_dalk,dco3_ddic) ! returning the same units as inputs (?)
implicit none
integer(kind=4),intent(in) :: nz
integer(kind=4),intent(out)::info
real(kind=8),intent(in)::dic(nz),alk(nz),tmp,sal,dep
real(kind=8),intent(out)::dco3_dalk(nz),dco3_ddic(nz)
real(kind=8)::ph(nz),co2(nz),hco3(nz),co3(nz)
real(kind=8)::a(nz),b(nz),c(nz),db_dalk(nz),dc_dalk(nz)
real(kind=8)::db_ddic(nz),dc_ddic(nz),dph_dalk(nz),dph_ddic(nz)
real(kind=8)::calceq1,calceq2,k1,k2

info=0
k1=calceq1(tmp,sal,dep)
k2=calceq2(tmp,sal,dep)
a=1d0
b=(1d0-dic/alk)*k1
c=(1d0-2d0*dic/alk)*k1*k2
ph = (-b+(b*b-4d0*a*c)**0.5d0)/2d0/a
if (any(ph<0d0)) then
    print*,'... unsable to calculate ph'
    ! print*,dic
    ! print*,alk
    ! print*,ph
    ! stop
    info=1
    return
endif
co2 = alk/(k1/ph+2d0*k1*k2/ph/ph)
! hco3=co2*k1/ph
hco3=alk/(1d0+2d0*k2/ph)
! co3=hco3*k2/ph
co3=alk/(ph/k2+2d0)
! ph=-log10(ph)

db_dalk = k1*(-1d0)*dic*(-1d0)/alk/alk
dc_dalk = k1*k2*(-2d0)*dic*(-1d0)/alk/alk
db_ddic = k1*(-1d0/alk)
dc_ddic = k1*k2*(-2d0/alk)
dph_dalk = -0.5d0*db_dalk + 0.5d0*0.5d0*(b*b-4d0*c)**(-0.5d0)*(2d0*b*db_dalk - 4d0*dc_dalk)
dph_ddic = -0.5d0*db_ddic + 0.5d0*0.5d0*(b*b-4d0*c)**(-0.5d0)*(2d0*b*db_ddic - 4d0*dc_ddic)
dco3_dalk = 1d0/(ph/k2+2d0) + alk*(-1d0)/((ph/k2+2d0)**2d0)*(1d0/k2)*dph_dalk
dco3_ddic = 0d0/(ph/k2+2d0) + alk*(-1d0)/((ph/k2+2d0)**2d0)*(1d0/k2)*dph_ddic

endsubroutine calcdevs
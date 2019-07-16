program caco3_fort
implicit none 
! contains 
real(kind=8)::ccflxi,om2cc,dtinput,dep,ztot
character*255::runname,biotmode,co2chem,runmode
logical::oxonly

call getinput_v2(ccflxi,om2cc,dtinput,runname,dep,oxonly,biotmode)
call caco3(ccflxi,om2cc,dep,dtinput,runname,oxonly,biotmode) 

endprogram 

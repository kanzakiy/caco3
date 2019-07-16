subroutine caco3_python(ccflxi,om2cc,dep,dtinput,runname,oxonly,biotmode)
implicit none 
real(kind=8),intent(in)::ccflxi,om2cc,dtinput,dep
character*255,intent(in)::runname,biotmode
logical,intent(in)::oxonly

call caco3(ccflxi,om2cc,dep,dtinput,runname,oxonly,biotmode) 

endsubroutine

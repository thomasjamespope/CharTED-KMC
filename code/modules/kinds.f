!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% kinds.f
!%
!% module kinds: 
!%  Parameters used to define fortran varaible 
!%  specifications
!%
!% Copyright (C) 2020 Yvelin Giret & Thomas Pope
!% \author Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%% simple variable definitions %%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

module kinds
    use, intrinsic :: ISO_FORTRAN_ENV
    implicit none
    integer, parameter, public :: ch = 100
    integer, parameter, public :: fm = 50
    integer, parameter         :: long_integer = selected_int_kind(16)
    integer, parameter         :: max_integer  = 2_8**31 - 1
    character(ch)              :: delim0 = "=================================================&
                                          &=================================================="
    character(ch)              :: delim1 = "|================================================&
                                          &=================================================|"
    character(ch)              :: delim2 = "|------------------------------------------------&
                                          &-------------------------------------------------|"
    character(ch)              :: delim3 = "|************************************************&
                                          &*************************************************|"
    character(ch)              :: blank1 = "|                                                &
                                          &                                                 |"
    character(ch)              :: rsrtdl = "#################################################&
                                          &##################################################"
    character(fm)              :: fmdlm1  = '(t2,a100)'
    character(fm)              :: fmdlm2  = '(t2,a100,/,t2,a100)'
    character(fm)              :: fmdlm3  = '(t2,a100,/,t2,a100,/,t2,a100)'
    character(fm)              :: fmchr1  = '(t2,"|",t5,a,t100,"|")'
    character(fm)              :: fmchr2  = '(t2,"|",t5,a,2x,a,t100,"|")'
    character(fm)              :: fmchr3  = '(t2,"|",t5,a,2x,a,2x,a,t100,"|")'
    character(fm)              :: fmchfl  = '(t2,"|",t5,a,2x,f20.12,t100,"|")'
    character(fm)              :: fmcflc  = '(t2,"|",t5,a,2x,f20.12,2x,a,t100,"|")'
    character(fm)              :: fmches  = '(t2,"|",t5,a,2x,es20.12,t100,"|")'
    character(fm)              :: fmchin  = '(t2,"|",t5,a,2x,i0,t100,"|")'
    character(fm)              :: fmchrv  = '(t2,"|",t5,a,2x,f5.2,2x,f5.2,2x,f5.2,t100,"|")'
    character(fm)              :: rdchr2  = '(t2,"|",t5,a,t55,":",t60,a,t100,"|")'
    character(fm)              :: erchr2  = '(t2,"|",t5,a,t47,":",t52,a,t100,"|")'
    character(fm)              :: erchfl  = '(t2,"|",t5,a,t47,":",t52,f10.5,t100,"|")'
    character(fm)              :: erches  = '(t2,"|",t5,a,t47,":",t52,es20.12,t100,"|")'
    character(fm)              :: erchin  = '(t2,"|",t5,a,t47,":",t52,i0,t100,"|")'
end module kinds



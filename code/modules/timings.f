!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% module timings.f: 
!%
!% Copyright (C) 2021 Thomas Pope & Yvelin Giret
!% \author Thomas Pope & Yvelin Giret
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    module timings
    use kinds
    use constants
    implicit none
    integer(8)          :: start_run
    double precision            :: start
    double precision, parameter :: mins =   60.0d0
    double precision, parameter :: hrs  = 3600.0d0
    contains 
!=====================================================================!
!   initialize the timing information                                 !
!=====================================================================!
    subroutine timing_init()
    implicit none
    call cpu_time(start)
    call system_clock(start_run)        
    write(screen,fmdlm1) blank1
    write(screen,fmchr1) 'Program started at:'
    write(screen,fmdlm2) delim2, blank1
    call time_now()
    write(screen,fmdlm2) blank1, delim1
    flush(screen)
    endsubroutine timing_init
!=====================================================================!
!   finalize and print the timing information                         !
!=====================================================================!
    subroutine timing_end()
    implicit none
    integer(8) :: finish_run
    double precision   :: finish, time_diff
    integer(8) :: clock_rate, clock_max
    write(screen,fmdlm3) blank1, delim1, blank1
    write(screen,fmchr1)  'PROGRAM ENDED NORMALLY!'
    write(screen,fmchr1)  'Program finished at:'
    call time_now()
    write(screen,fmdlm3) blank1, delim2, blank1
    write(screen,fmchr1)  'Timing information:'
    call cpu_time(finish)
    call system_clock(finish_run,clock_rate,clock_max)
    time_diff = finish - start
    write(screen,10) 'cpu time', int(time_diff / hrs),  &
     &int(int(time_diff / mins) - (int(time_diff / hrs) * mins)), mod(time_diff,mins)
    time_diff = dble(finish_run - start_run) / dble(clock_rate)
    write(screen,10) 'run time', int(time_diff / hrs),  &
     &int(int(time_diff / mins) - (int(time_diff / hrs) * mins)), mod(time_diff,mins)
    write(screen,fmdlm1) blank1
    write(screen,fmdlm1) delim0
    flush(screen)
10  format(t2,"|",t5,a,t14,"=",1x,i5,1x,"[h]",1x,i3,1x,"[m]",1x,f6.2,"[s]",t100,"|")
    endsubroutine timing_end
!=====================================================================!
!   print the time and date                                           !
!=====================================================================!
    subroutine time_now()
    implicit none
    integer, dimension(8) :: values
    character(8)          :: date
    character(10)         :: time
    character(5)          :: zone
    character(9)          :: month
    month = ""
    call date_and_time(date,time,zone,values)
    call date_and_time(ZONE=zone)
    select case(values(2))
      case(1);  month = "January"
      case(2);  month = "February"
      case(3);  month = "March"
      case(4);  month = "April"
      case(5);  month = "May"
      case(6);  month = "June"
      case(7);  month = "July"
      case(8);  month = "August"
      case(9);  month = "September"
      case(10); month = "October"
      case(11); month = "November"
      case(12); month = "December"
    endselect
    write(screen,'(t2,"|",t5,"Date:",2x,i2,1x,A,1x,i4,t100,"|")') values(3), month, values(1)
    write(screen,'(t2,"|",t5,"Time:",2x,i2,"h",i2,"m",i2,"s",t100,"|")') values(5:7)
    write(screen,'(t2,"|",t5,"zone:  UTC",A,t100,"|")') zone
    flush(screen)
    endsubroutine time_now
    endmodule timings

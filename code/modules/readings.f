!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% module readings.f90: ...description...
!%
!% Copyright (C) 2020 Yvelin Giret & Thomas Pope
!% \author Yvelin Giret & Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module readings
    use kinds
    use constants
    implicit none
    character(ch), parameter :: line_error="ERROR::SET:DEFAULT" 
    contains
    
!#####################################################################!
!#### Set an input string to all caps ################################!
!#####################################################################!   
    function capslock(string) result(CAPstrng)
    implicit none
    integer                  :: i,ic
    character(*)             :: string
    character(ch)            :: CAPstrng
    character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    CAPstrng=string
    do i = 1, LEN_TRIM(CAPstrng)
      ic = INDEX(low, CAPstrng(i:i))
      if(ic.gt.0) CAPstrng(i:i) = cap(ic:ic)
    enddo
    endfunction capslock
!#####################################################################!
!#### Search function to find and return the appropriate line ########!
!#####################################################################!   
    function stopat(u,search) result(line)
    implicit none
    integer, intent(in)      :: u
    character(*), intent(in) :: search
    character(ch)            :: line,dum
    integer                  :: ios
    ios = 0
    do while (index(line,search).eq.0 .and. ios.eq.0)
      read(u,'(a100)',iostat=ios) dum
      if(index(dum,"#").eq.0) line=capslock(dum)
    enddo
    if(ios.lt.0) then
      call error_ios_neg(search)
      line=line_error
    elseif(ios.gt.0) then
      call error_ios_pos(search)
      line=line_error
    endif
    endfunction stopat
!#####################################################################!
!#### Error messages for if the search fails #########################!
!#####################################################################!   
    subroutine error_ios_pos(search)
    use constants
    implicit none
    character(*),  intent(in) :: search
    write(screen,fmdlm1) delim3
    write(screen,erchr2) 'WARNING: while trying to read: ', adjustl(trim(search))
    write(screen,fmchr1) '-------- IOSTAT > 0'
    write(screen,fmchr1) '-------- default value is set!'
    write(screen,fmdlm1) delim3
    endsubroutine error_ios_pos
!#####################################################################!   
    subroutine error_ios_neg(search)
    use constants
    implicit none
    character(*),  intent(in) :: search
    write(screen,fmdlm1) delim3
    write(screen,erchr2) 'WARNING: while trying to read: ', adjustl(trim(search))
    write(screen,fmchr1) '-------- IOSTAT < 0'
    write(screen,fmchr1) '-------- default value is set!'
    write(screen,fmdlm1) delim3
    endsubroutine error_ios_neg
!#####################################################################!
!#### Statement Search ###############################################!
!#####################################################################!   
    function read_statements(u,search) result(output)
    implicit none
    integer, intent(in)        :: u
    character(*),  intent(in)  :: search
    logical                    :: output
    ! Local variables:
    character(ch) :: tmp
    integer       :: ios
    rewind(u)
    ios = 0
    do while (index(tmp,search).eq.0 .and. ios.eq.0)
      read(u,'(a100)',iostat=ios) tmp
    enddo
    if(ios.ne.0) then 
      output = .false.
    else
      output = .true.
    endif
    endfunction read_statements 
!#####################################################################!
!#### String Variable Reader #########################################!
!#####################################################################!   
    subroutine read_1_string(u,r,search,output,def_val,val1,val2)
    implicit none
    integer, intent(in)        :: u
    character(*),  intent(in)  :: search, output
    character(*),  intent(in)  :: def_val
    character(ch)              :: dumm
    character(*),  intent(out) :: r
    character(*),  intent(in), optional :: val1, val2
    ! Local variables:
    integer       :: ios
    character(ch) :: tmp, dummy
    logical       :: read_fail
    rewind(u)
    tmp = stopat(u,search) 
    if(trim(adjustl(tmp)).eq.line_error) then
      r = def_val
    else
      read(tmp,*,iostat=ios) dummy, dumm
      if(ios.lt.0) then
        call error_ios_neg(search)
        r = def_val
      elseif(ios.gt.0) then
        call error_ios_pos(search)
        r = def_val
      endif
    endif
    r = adjustl(trim(capslock(dumm)))
    read_fail = .true.
    if(r.eq.def_val) read_fail=.false.
    if(present(val1).and.read_fail) then; if(r.eq.val1) read_fail=.false.; endif
    if(present(val2).and.read_fail) then; if(r.eq.val2) read_fail=.false.; endif
    if(read_fail) then
      write(screen,fmdlm1) delim3
      write(screen,erchr2) 'WARNING: unauthorized value for', adjustl(trim(search))
      write(screen,erchr2) '-------- input value', r
      r = adjustl(trim(def_val))
      write(screen,erchr2) '-------- value replaced by default value', r
      write(screen,fmdlm1) delim3
    endif
    write(screen,rdchr2) output, adjustl(r)
    endsubroutine read_1_string 
!#####################################################################!
!#### Float Variable Reader ##########################################!
!#####################################################################!   
    subroutine read_1_real(u,r,search,output,def_val,frmt)
    implicit none
    integer, intent(in)        :: u
    character(*),  intent(in)  :: search, output
    double precision,      intent(in)  :: def_val
    character(*),  intent(in)  :: frmt
    double precision,      intent(out) :: r
    ! Local variables:
    integer       :: ios
    character(ch) :: tmp, dummy
    character(90) :: char_r
    rewind(u)
    tmp = stopat(u,search) 
    if(trim(adjustl(tmp)).eq.line_error) then
      r = def_val
    else
      read(tmp,*,iostat=ios) dummy, r
      if(ios.lt.0) then
        call error_ios_neg(search)
        r = def_val
      elseif(ios.gt.0) then
        call error_ios_pos(search)
        r = def_val
      endif
    endif
    write(char_r,frmt) r
    write(screen,rdchr2) output, trim(adjustl(char_r))
    endsubroutine read_1_real
!#####################################################################!
!#### Error message for Float Variable Reader ########################!
!#####################################################################!   
    subroutine wrong_value_1_real(r,search,def_val,frmt)
    implicit none
    character(*),  intent(in)    :: search
    double precision,      intent(in)    :: def_val
    character(fm), intent(in)    :: frmt
    double precision,      intent(inout) :: r
    write(screen,fmdlm1) delim3
    write(screen,erchr2) 'WARNING: unauthorized value for', adjustl(trim(search))
    write(screen,frmt) '-------- input value', r
    r = def_val
    write(screen,frmt) '-------- value replaced by default value', r
    write(screen,fmdlm1) delim3
    endsubroutine wrong_value_1_real
!#####################################################################!
!#### Integer Variable Reader ########################################!
!#####################################################################!   
    subroutine read_1_integer(u,r,search,output,def_val)
    implicit none
    integer, intent(in)        :: u
    character(*),  intent(in)  :: search, output
    integer,       intent(in)  :: def_val
    integer,       intent(out) :: r
    ! Local variables:
    integer       :: ios
    character(ch) :: tmp, dummy
    character(90) :: char_r
    rewind(u)
    tmp = stopat(u,search) 
    if(trim(adjustl(tmp)).eq.line_error) then
      r = def_val
    else
      read(tmp,*,iostat=ios) dummy, r
      if(ios.lt.0) then
        call error_ios_neg(search)
        r = def_val
      elseif(ios.gt.0) then
        call error_ios_pos(search)
        r = def_val
      endif
    endif
    write(char_r,'(i0)') r
    write(screen,rdchr2) output, trim(adjustl(char_r))
    endsubroutine read_1_integer
!#####################################################################!
!#### Error message for Integer Variable Reader ######################!
!#####################################################################!   
    subroutine wrong_value_1_int(r,search,def_val)
    implicit none
    character(*),  intent(in)    :: search
    integer,       intent(in)    :: def_val
    integer,       intent(inout) :: r
    write(screen,fmdlm1) delim3
    write(screen,erchr2) 'WARNING: unauthorized value for', adjustl(trim(search))
    write(screen,erchin) '-------- input value', r
    r = def_val
    write(screen,erchin) '-------- value replaced by default value', r
    write(screen,fmdlm1) delim3
    end subroutine wrong_value_1_int
!#####################################################################!
!#### Logical Variable Reader ########################################!
!#####################################################################!   
    subroutine read_1_logical(u,r,search,output,def_val)
    implicit none
    integer, intent(in)        :: u
    character(*),  intent(in)  :: search, output
    logical,       intent(in)  :: def_val
    logical,       intent(out) :: r
    ! Local variables:
    integer       :: ios
    character(ch) :: tmp, dummy
    rewind(u)
    tmp = stopat(u,search) 
    if(trim(adjustl(tmp)).eq.line_error) then
      r = def_val
    else
      read(tmp,*,iostat=ios) dummy, r
      if(ios.lt.0) then
        call error_ios_neg(search)
        r = def_val
      elseif(ios.gt.0) then
        call error_ios_pos(search)
        r = def_val
      endif
    endif
    if(r) then
      write(screen,rdchr2) output, "TRUE"
    else
      write(screen,rdchr2) output, "FALSE"
    endif
    end subroutine read_1_logical
!#####################################################################!
!#### Float Vector Variable Reader ###################################!
!#####################################################################!   
    subroutine read_3_real(u,v,search,output,def)
    implicit none
    integer, intent(in)                 :: u
    character(*),  intent(in)           :: search, output
    double precision, dimension(3), intent(in)  :: def
    double precision, dimension(3), intent(out) :: v
    ! Local variables:
    integer       :: ios
    character(ch) :: tmp, dummy
    character(90) :: char_r
    rewind(u)
    tmp = stopat(u,search) 
    if(trim(adjustl(tmp)).eq.line_error) then
      v = def
    else
      read(tmp,*,iostat=ios) dummy, v
      if(ios.lt.0) then
        call error_ios_neg(search)
        v = def
      elseif(ios.gt.0) then
        call error_ios_pos(search)
        v = def
      endif
    endif
    write(char_r,'(2(f5.2,","),f5.2)') v
    write(screen,rdchr2) output, trim(adjustl(char_r))
    endsubroutine read_3_real
!#####################################################################!
!#### Integer Vector Variable Reader #################################!
!#####################################################################!   
    subroutine read_3_integers(u,x,y,z,search,output,defx,defy,defz)
    implicit none
    integer, intent(in)        :: u
    character(*),  intent(in)  :: search, output
    integer,       intent(in)  :: defx, defy, defz
    integer,       intent(out) :: x, y, z
    ! Local variables:
    integer       :: ios
    character(ch) :: tmp, dummy
    character(90) :: char_r
    rewind(u)
    tmp = stopat(u,search) 
    if(trim(adjustl(tmp)).eq.line_error) then
     x = defx; y = defy; z = defz
    else
      read(tmp,*,iostat=ios) dummy, x, y, z
      if(ios.lt.0) then
        call error_ios_neg(search)
        x = defx; y = defy; z = defz
      elseif(ios.gt.0) then
        call error_ios_pos(search)
        x = defx; y = defy; z = defz
      endif
    endif
    write(char_r,'(3(i0,2x))') x, y, z
    write(screen,rdchr2) output, trim(adjustl(char_r))
    endsubroutine read_3_integers
!#####################################################################!
!#### Float List Variable Reader #####################################!
!#####################################################################!   
    subroutine read_n_real(u,n,l,search,output,def)
    implicit none
    integer, intent(in)                 :: u
    character(*),           intent(in)  :: search, output
    integer,                intent(in)  :: n
    double precision,               intent(in)  :: def
    double precision, dimension(n), intent(out) :: l
    ! Local variables:
    integer       :: ios
    character(ch) :: tmp, dummy
    character(90) :: char_r
    rewind(u)
    tmp = stopat(u,search) 
    if(trim(adjustl(tmp)).eq.line_error) then
      l = def
    else
      read(tmp,*,iostat=ios) dummy, l
      if(ios.lt.0) then
        call error_ios_neg(search)
        l = def
      elseif(ios.gt.0) then
        call error_ios_pos(search)
        l = def
      endif
    endif
    write(char_r,'(10(f5.2,1x))') l
    write(screen,rdchr2) output, trim(adjustl(char_r))
    endsubroutine read_n_real
!#####################################################################!
!#### Crossing Rates Reader ##########################################!
!#####################################################################!
    function read_crossing_rates(u,N_SING,N_TRIP,N_ISC) result(ISC)
    implicit none
    integer, intent(in)        :: u
    integer,       intent(in)        :: N_SING, N_TRIP, N_ISC
    double precision, dimension(N_ISC,N_ISC) :: ISC
    integer                          :: i, j
    character(ch)                    :: text1, text2
    do i=1,N_SING
      do j=1,N_SING
        if(j.ne.i) then
          write(text1,'("S",i0,"->S",i0)') i, j
          write(text2,'("Intersystem Crossing Rate: ",a," [s^-1]")') trim(adjustl(text1))
          call read_1_real(u,ISC(i,j),trim(adjustl(text1)),trim(adjustl(text2)),0.0d0,"(es10.3)")
        endif
      enddo
      do j=1,N_TRIP
        write(text1,'("S",i0,"->T",i0)') i, j
        write(text2,'("Intersystem Crossing Rate: ",a," [s^-1]")') trim(adjustl(text1))
        call read_1_real(u,ISC(i,j+N_SING),trim(adjustl(text1)),trim(adjustl(text2)),0.0d0,"(es10.3)")
      enddo
    enddo
    do i=1,N_TRIP
      do j=1,N_SING
        write(text1,'("T",i0,"->S",i0)') i, j
        write(text2,'("Intersystem Crossing Rate: ",a," [s^-1]")') trim(adjustl(text1))
        call read_1_real(u,ISC(i+N_SING,j),trim(adjustl(text1)),trim(adjustl(text2)),0.0d0,"(es10.3)")
      enddo
      do j=1,N_TRIP
        if(j.ne.i) then
          write(text1,'("T",i0,"->T",i0)') i, j
          write(text2,'("Intersystem Crossing Rate: ",a," [s^-1]")') trim(adjustl(text1))
          call read_1_real(u,ISC(i+N_SING,j+N_SING),trim(adjustl(text1)),trim(adjustl(text2)),0.0d0,"(es10.3)")
        endif
      enddo
    enddo
    endfunction read_crossing_rates

    function read_count_states(u,state) result(N)
    implicit none
    integer, intent(in)             :: u
    character(*), intent(in)        :: state
    logical                         :: found
    character(ch)                   :: search
    integer                         :: N, ios
    N=0; ios=0
    write(search,'(a,"_ENERGY")') state
    found = isitthere(u,search)
    if(.not.found) N = N + 1
    do while(found)
      N = N + 1 
      write(search,'(a,"+",i0,"_ENERGY")') state, N   
      found = isitthere(u,search)
    enddo
    N = N - 1
    endfunction read_count_states
    
    function isitthere(u,search) result(found)
    implicit none
    integer, intent(in)             :: u
    character(ch), intent(in)       :: search
    logical                         :: found
    character(ch)                   :: line, dum
    integer                         :: ios
    found=.false.
    rewind(u)
    line=""; ios=0;
    do while (index(line,trim(adjustl(search))).eq.0.and.ios.eq.0)
      read(u,'(a100)',iostat=ios) dum
      if(index(dum,"#").eq.0) line=capslock(dum)
    enddo
    if(ios.eq.0) found=.true.    
    endfunction isitthere

    function read_E_site(Nx,Ny,Nz,namefile) result(e_site)
        implicit none
        character(18), intent(in)    :: namefile
        integer,       intent(in)    :: Nx, Ny, Nz
        double precision,dimension(Nx,Ny,Nz) :: E_site

        write(screen,fmchr2) 'Reading the file: ', namefile
        call read_real_matrix_18(namefile,Nx,Ny,Nz,E_site)
        write(screen,fmdlm1) delim2
        flush(screen)
    end function read_E_site

    subroutine read_charge(Nx,Ny,Nz,hole_distri,namefile)
        implicit none
        character(26), intent(in) :: namefile
        integer,       intent(in) :: Nx, Ny, Nz
        logical,      intent(out) :: hole_distri(Nx,Ny,Nz)

        write(screen,fmchr2) 'Reading the file: ', namefile
        call read_bool_matrix_26(namefile,Nx,Ny,Nz,hole_distri)
        write(screen,fmdlm1) delim2
        flush(screen)
    end subroutine


    subroutine inquire_inputs(namefile,restart)
       !
       ! Determines whether the simulation is from scratch or from a previous run.
       ! Checks the existence of either input.file or restart.file accordingly.
       ! Checks the existence of stop.file, and create if not in the root folder.
       !
       implicit none
       character(ch), intent(in)  :: namefile
       logical,       intent(out) :: restart
       integer, parameter         :: u=1222

       ! Local variables
       character(ch) :: name1, name2, name3
       logical       :: exist_file

       name1 = "input.file"
       name2 = "CHRG.RESTART"
       name3 = "stop.file"
       exist_file = .false.
       restart = .false.

       write(screen,fmdlm1) blank1
       write(screen,fmchr1)  'Initializing simulation:'
       write(screen,fmdlm2) delim2, blank1

       ! Check for input/restart file:
       if(namefile == name1) then
          write(screen,fmchr1) 'Simulation starts from scratch.'
          inquire(file = name1, exist = exist_file)
          if(exist_file) then
             write(screen,fmchr3) 'File', adjustl(trim(namefile)), 'found in root folder.'
          else
             write(screen,'(2x,a)') '*****************************************************'
             write(screen,'(2x,A,1x,A,1x,A)') 'ERROR: File', adjustl(trim(namefile)), 'not found in root folder.'
             write(screen,'(2x,A)') 'Please provide the required file in the root folder.'
             write(screen,'(2x,A)') 'You should either provide input.file or restart.file!'
             write(screen,'(2x,A)') 'ABORT PROGRAM!'
             write(screen,'(2x,a)') '*****************************************************'
             stop
          end if
       else if(namefile == name2) then
          stop 'ERROR: checkpointing is no longer implemented'
          restart = .true.
          write(screen,fmchr1) 'Simulation starts from a previous run.'
          inquire(file = name2, exist = exist_file)
          if(exist_file) then
             write(screen,fmchr3) 'File', adjustl(trim(namefile)), 'found in root folder.'
          else
             write(screen,'(2x,a)') '*****************************************************'
             write(screen,'(2x,A,1x,A,1x,A)') 'ERROR: File', adjustl(trim(namefile)), 'not found in root folder.'
             write(screen,'(2x,A)') 'Please provide the required file in the root folder.'
             write(screen,'(2x,A)') 'You should either provide input.file or restart.file!'
             write(screen,'(2x,A)') 'ABORT PROGRAM!'
             write(screen,'(2x,a)') '*****************************************************'
             stop
          end if
       else
          write(screen,'(2x,a)') '*****************************************************'
          write(screen,'(2x,A)') 'ERROR: The program expects to read either'
          write(screen,'(2x,A)') 'input.file or restart.file.'
          write(screen,'(2x,A,2x,A)') 'Instead you have asked for:', adjustl(trim(namefile))
          write(screen,'(2x,A)') 'Please change the name of your input file.'
          write(screen,'(2x,A)') 'ABORT PROGRAM!'
          write(screen,'(2x,a)') '*****************************************************'
          stop
       end if

       ! Check for stop file:
       exist_file = .false.
       inquire(file = name3, exist = exist_file)
       if(exist_file) then
          write(screen,fmchr3) 'File', adjustl(trim(name3)), 'found in root folder.'
          open(10, file=name3)
          write(10,102) ' TERMINATE', '.false.'
          close(10)
          write(screen,fmchr1) 'TERMINATE set to .false.'
       else
          write(screen,fmchr3) 'File', adjustl(trim(name3)), 'not found in root folder.'
          open(UNIT=u, file=name3, form='formatted', status='unknown')
              write(u,102) ' TERMINATE', '.false.'
          close(UNIT=u)
          write(screen,fmchr3) 'File', adjustl(trim(name3)), 'created in root folder.'
       end if
       write(screen,fmdlm1) blank1
       flush(screen)
102    format(1x,A,4x,A)
    end subroutine inquire_inputs

    subroutine read_stop(stop_run)
        implicit none
        logical, intent(inout) :: stop_run
        character(ch)          :: tmp
        integer                :: u=100
        open(u, file='stop.file')
        rewind u
        read(u,*) tmp, stop_run
        close(u)
        if(stop_run) then
           write(screen,fmchr1) 'WARNING: User asked to end the simulation! One last iteration will be performed'
           flush(screen)
        endif
    endsubroutine read_stop

    subroutine read_restart_matrices(Nx,Ny,Nz,E_site,hole_distri,f_local, &
                    b_local,name_E_site,name_forward,name_backward,       &
                    name_charge)
        implicit none
        integer,       intent(in)  :: Nx, Ny, Nz
        character(18), intent(in)  :: name_E_site
        character(24), intent(in)  :: name_forward
        character(25), intent(in)  :: name_backward
        character(26), intent(in)  :: name_charge
        integer,       intent(out) :: f_local(Nx,Ny,Nz) 
        integer,       intent(out) :: b_local(Nx,Ny,Nz)
        double precision,      intent(out) :: E_site(Nx,Ny,Nz)
        logical,       intent(out) :: hole_distri(Nx,Ny,Nz)

        write(screen,*) 'Reading the file: ', name_E_site
        call read_real_matrix_18(name_E_site,Nx,Ny,Nz,E_site)

        write(screen,*) 'Reading the file: ', name_forward
        call read_int_matrix_24(name_forward,Nx,Ny,Nz,f_local)

        write(screen,*) 'Reading the file: ', name_backward
        call read_int_matrix_25(name_backward,Nx,Ny,Nz,b_local)

        write(screen,*) 'Reading the file: ', name_charge
        call read_bool_matrix_26(name_charge,Nx,Ny,Nz,hole_distri)

        write(screen,*) trim(adjustl(delim1))
        flush(screen)
    end subroutine read_restart_matrices

    subroutine read_int_matrix_24(namefile,Nx,Ny,Nz,f)
        implicit none
        character(24), intent(in)  :: namefile
        integer,       intent(in)  :: Nx, Ny, Nz
        integer,       intent(out) :: f(Nx,Ny,Nz)

        ! Local variables:
        integer :: ix, iy, iz, u=100

        open(u, file=namefile)
        do iz=1, Nz
           do iy=1, Ny
              do ix=1, Nx
                  read(u,*) f(ix,iy,iz)
              end do
           end do
        end do
        close(u)
    end subroutine read_int_matrix_24

    subroutine read_int_matrix_25(namefile,Nx,Ny,Nz,f)
        implicit none
        character(25), intent(in)  :: namefile
        integer,       intent(in)  :: Nx, Ny, Nz
        integer,       intent(out) :: f(Nx,Ny,Nz)

        ! Local variables:
        integer :: ix, iy, iz, u=100

        open(u, file=namefile)
        do iz=1, Nz
           do iy=1, Ny
              do ix=1, Nx
                  read(u,*) f(ix,iy,iz)
              end do
           end do
        end do
        close(u)
    end subroutine read_int_matrix_25

    subroutine read_real_matrix_18(namefile,Nx,Ny,Nz,f)
        implicit none
        character(18), intent(in)  :: namefile
        integer,       intent(in)  :: Nx, Ny, Nz
        double precision,      intent(out) :: f(Nx,Ny,Nz)

        ! Local variables:
        integer :: ix, iy, iz, tmp1, tmp2, tmp3, u=100

        open(u, file=namefile)
        do iz=1, Nz
           do iy=1, Ny
              do ix=1, Nx
                  read(u,*) tmp1, tmp2, tmp3, f(ix,iy,iz)
              end do
           end do
        end do
        close(u)
        f(:,:,:) = f(:,:,:) * eV_to_Eh
    end subroutine read_real_matrix_18

    subroutine read_bool_matrix_26(namefile,Nx,Ny,Nz,f)
        implicit none
        character(26), intent(in)  :: namefile
        integer,       intent(in)  :: Nx, Ny, Nz
        logical,       intent(out) :: f(Nx,Ny,Nz)

        ! Local variables:
        integer :: ix, iy, iz, tmp1, tmp2, tmp3, u=100

        open(u, file=namefile)
        do iz=1, Nz
           do iy=1, Ny
              do ix=1, Nx
                  read(u,*) tmp1, tmp2, tmp3, f(ix,iy,iz)
              end do
           end do
        end do
        close(u)
    end subroutine read_bool_matrix_26
end module readings

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% module event.f90: ...description...
!%
!% Copyright (C) 2020 Yvelin Giret and Thomas Pope
!% \author Yvelin Giret, Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module event
    use kinds
    use type_charge
    implicit none
!   general indexing arrays    
!   hopping parameters
    integer                                     :: hopboxN      ! dimension of the hopbox
    type(indexes),    allocatable, dimension(:) :: hopbox       ! index array containing the neighbouring hop sites
    double precision, allocatable, dimension(:) :: hopfld       ! change in E-field due to a hop
    double precision, dimension(-1:1,-1:1,-1:1) :: hopdir       ! flag for whether the hop is in the direction of the field
    double precision, dimension(-1:1,-1:1,-1:1) :: erxyz        ! array containing the spacial exponential term of the rate equation 
    double precision                            :: econs        ! the leading constant of the rate equation
!   selection parameters
    double precision, allocatable, dimension(:) :: rate         ! selection array: contains the hopping rates
    type(indexes),    allocatable, dimension(:) :: svo          ! selection array: contains origin of hops
    type(indexes),    allocatable, dimension(:) :: svd          ! selection array: contains destinations of hops
    integer,          allocatable, dimension(:) :: iEd          ! selection array: contains level index of destination
!   grid partition array

    contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%% routine to initialize various constant elements of the events routines %%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine init_events(holes,elecs)
    use constants
    use param
    use my_mpi
    use type_charge
    implicit none
    type(charge)  :: holes, elecs
    type(indexes) :: I
    type(coor)    :: IR
    type(coor)    :: fieldvector
    integer       :: x,y,z,j
    hopboxN = (N_hopbox * 2 + 1)**3 - 1
    max_events = hopboxN * N_chrg * (max(N_HOMO,N_LUMO) + 1) * 2
! allocate the arrays
    allocate(rate(max_events))
    allocate(svo(max_events),svd(max_events),iEd(max_events))
    allocate(hopbox(hopboxN),hopfld(hopboxN))
    if(inc_hole) call holes%rate_init(hopboxN)
    if(inc_elec) call elecs%rate_init(hopboxN)
! set up the constant and the grid arrays
    j             = 0
    fieldvector%x = field_vector(1); fieldvector%y = field_vector(2); fieldvector%z = field_vector(3)
    hopdir = 0
    do x=-N_hopbox,N_hopbox;     I%x = x; IR%x = real(x)
      do y=-N_hopbox,N_hopbox;   I%y = y; IR%y = real(y)
        do z=-N_hopbox,N_hopbox; I%z = z; IR%z = real(z)
          if(x.ne.0 .or. y.ne.0 .or. z.ne.0) then
            erxyz(x,y,z)  = hop_prefactor * exp( - local * I%mag() * lattice)
            j             = j + 1
            hopbox(j)     = I
            hopdir(x,y,z) = (IR.sp.fieldvector)
            hopfld(j)     = - hopdir(x,y,z) * field_step
          else
            erxyz(x,y,z) = 0.0d0
          endif
        enddo
      enddo
    enddo
    if(mrcs_rate) then
      econs = 1.0d0 / (4.0d0 * lambda * kB_in_AU * temperature)
    else
      econs = 1.0d0 / (kB_in_AU * temperature)   
    endif
    endsubroutine init_events

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%% wrapper routine for calculating hopping event %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine calculate_events(holes,elecs,sings,trips,rate_tot)
     use param
     implicit none      
     type(charge)     :: holes, elecs ! Hole and electron related data
     type(charge)     :: sings, trips ! Singlet and Triplet related data
     double precision :: rate_tot
     if(guest_host_flag) then
       call calculate_gh_events_double(holes,elecs,sings,trips,rate_tot)      
     else
      if(inc_excn) then
       call calculate_events_double(holes,elecs,sings,trips,rate_tot)
      elseif(inc_hole) then
       call calculate_events_single(holes,1,rate_tot)
      elseif(inc_elec) then
       call calculate_events_single(elecs,-1,rate_tot)
      else
       stop 'error in events'
      endif
     endif
    endsubroutine calculate_events

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%% routine for calculating hopping rate for a given energy difference %%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function myrate(dE) result(k)
    use param
    implicit none
    double precision, intent(in) :: dE
    double precision             :: k
    if(mrcs_rate) then
     k = exp( - (dE+lambda)**2 * econs)
    else
     if(dE.gt.0.0d0) then
      k = exp( - dE * econs)
     else
      k = 1.0d0
     endif
    endif    
    endfunction myrate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%% routine for calculating hopping event in a sparse one-charge system %%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine calculate_events_single(chrg,sgn,rate_tot)
      use my_mpi
      use pbc
      use coulomb
      use general_functions
      implicit none      
!----- input constant
      type(charge)                        :: chrg                 ! Hole or electron related data
      integer, intent(in)                 :: sgn
!----- calc variables
      type(indexes)                       :: o                    ! origin point
      type(indexes)                       :: s                    ! step index
      type(indexes)                       :: t                    ! target point
      type(indexes)                       :: d                    ! destination point (target adjusted for periodicity)
      integer                             :: me                   ! identity of the hopping charge
      integer                             :: i, j                 ! random indexes
      integer                             :: events, event_keep   ! indexes for events
      integer                             :: i_state              ! index for state levels 
      double precision                    :: delta_E              ! energy difference associated with hop
      double precision                    :: del_E                ! energy difference associated with hop - without multiple levels
      double precision                    :: oEc                  ! Coloumb energy for the origin on a hopping event
      double precision                    :: eta
      integer, dimension(max_events)      :: sai                  ! selection array: general indexes for all events
      double precision, dimension(max_events) :: sum_rate
!----- output variables
      double precision                    :: rate_tot
      events = 0
      do me=1,chrg%N
        o  = chrg%R(me)       ! Origin point
        if(chrg%iE(me).gt.0) then
          events        = events + 1
          rate(events)  = chrg%relax_rate
          svo(events)   = o
          svd(events)   = o 
          iEd(events)   = 0
          sai(events)   = me
        endif
        if(chrg%log_rate(me)) then ! if the last hop has changed things, recalculate the coulomb difference
          oEc = coulomb_energy_single(chrg%N,o,chrg%R,me)
          do j=1,hopboxN
             s = hopbox(j)        ! set up the step index 
             t = o + s            ! find the target point
             d = indx_periodic(t) ! adjust target for periodicity
             if( .not. (chrg.chk.d) ) then      ! if there is nothing at the destination...
!               chrg%coul(j,me) = coulomb_energy_single(chrg%N,d,chrg%R,me) - oEc
               del_E = hopfld(j) + sgn * (coulomb_energy_single(chrg%N,d,chrg%R,me) - oEc)
               do i_state=0,chrg%NE
                 events       = events + 1 ! add an event
                 delta_E      = chrg%E(d%x,d%y,d%z,i_state) - chrg%E(o%x,o%y,o%z,chrg%iE(me)) + del_E
                 rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
                 chrg%rate(i_state,j,me) = rate(events)
                 sai(events)  = me
                 svo(events)  = o
                 svd(events)  = d 
                 iEd(events)  = i_state
               enddo
!             else
!               chrg%rate(:,j,me) = 0.d0
             endif
          enddo
          chrg%log_rate(me) = .false.
        else                       ! otherwise, use the previous coulomb difference
          do j=1,hopboxN
             s = hopbox(j)        ! set up the step index 
             t = o + s            ! find the target point
             d = indx_periodic(t) ! adjust target for periodicity
             if( .not. (chrg.chk.d) ) then      ! if there is nothing at the destination...
               do i_state=0,chrg%NE
                 events       = events + 1 ! add an event
                 rate(events) = chrg%rate(i_state,j,me)
                 sai(events)  = me
                 svo(events)  = o
                 svd(events)  = d 
                 iEd(events)  = i_state
               enddo
             endif
          enddo
        endif
      enddo
!------- nominate an event
      if(events.eq.0) stop 'ERROR: no events!!'
      if(rate(1).gt.1.0d30) then
       write(*,*) "ERROR: You have an infinite rate"
       write(*,*) rate(1)
       write(*,*) svo(1)%x, svo(1)%y, svo(1)%z, chrg%E(svo(1)%x,svo(1)%y,svo(1)%z,0)
       write(*,*) svd(1)%x, svd(1)%y, svd(1)%z, chrg%E(svd(1)%x,svd(1)%y,svd(1)%z,0)
       flush(6)
       stop
      endif
!------- Set up the Partial Sums
      sum_rate(1) = rate(1)
      do i=2,events
       sum_rate(i) = sum_rate(i-1) + rate(i) !Partial sums
      enddo
!------- Normalize the Partial Sums
      rate_tot          = sum_rate(events)
      
!      sum_rate(:events) = sum_rate(:events) / rate_tot
!------- Pick an Event with MC
      call RANDOM_NUMBER(eta)
      event_keep = findeta(sum_rate(:events),eta * rate_tot)
!      do event_keep=1,events; if(sum_rate(event_keep).gt.eta) exit; enddo
!------- update the charge
      chrg%o           = svo(event_keep)
      chrg%d           = svd(event_keep)
      chrg%iE(sai(event_keep)) = iEd(event_keep)
      s                        = chrg%hopid()
      call chrg%update(sai(event_keep),hopdir(s%x,s%y,s%z))   
      chrg%log_rate(sai(event_keep)) = .true.
      do i=1,chrg%N
       if(i.ne.sai(event_keep)) chrg%log_rate(i) = (chrg%d.nearto.chrg%R(i))
      enddo
      endsubroutine calculate_events_single

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%% routine for calculating hopping event in a sparse two-charge system %%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine calculate_events_double(holes,elecs,sings,trips,rate_tot)
      use constants
      use param
      use general_functions
      use my_mpi
      use pbc
      use coulomb
      use checks
      implicit none      
!----- input constant
      type(charge)                            :: holes, elecs       ! Hole and electron related data
      type(charge)                            :: sings, trips       ! Singlet and Triplet related data
!----- calc variables
      type(indexes)                           :: o                  ! origin point
      type(indexes)                           :: s                  ! step index
      type(indexes)                           :: t                  ! target point
      type(indexes)                           :: d                  ! destination point (target adjusted for periodicity)
      integer                                 :: me                 ! identity of the hopping charge
      integer                                 :: i, j               ! random indexes
      integer                                 :: events, event_keep ! indexes for events
      integer                                 :: i_state            ! index for state levels 
      double precision                        :: delta_E            ! energy difference associated with hop
      double precision                        :: del_E              ! energy difference associated with hop - without multiple levels
      double precision                        :: oEc                ! Coloumb energy for the origin on a hopping event
      double precision                        :: eta                ! random number
      integer, dimension(max_events)          :: sai, saj           ! selection array: general indexes for all events
      integer, dimension(max_events)          :: p_i                ! selection array: general indexes for all events
      double precision, dimension(max_events) :: sum_rate           ! cummulative sum of the rates
!      double precision                        :: time(5)
!----- output variables
      double precision                        :: rate_tot           ! sum of all the rates
      events = 0
      do me=1,holes%N
       if(holes%exists(me)) then ! does the hole exist?
        o = holes%R(me)          ! Origin point
        if(holes%iE(me).gt.0) then
         events       = events + 1
         rate(events) = holes%relax_rate
         svo(events)  = o
         svd(events)  = o 
         iEd(events)  = 0
         sai(events)  = me
         p_i(events)  = 0
        endif
        if(holes%log_rate(me)) then ! if the last hop has charged thing, recalculate the coulomb difference
          oEc = coulomb_energy_double(holes%N,elecs%N,o,holes%R,elecs%R,me)
          do j=1,hopboxN ! loop over the hopbox
           s = hopbox(j)        ! set up the step index 
           t = o + s            ! find the target point
           d = indx_periodic(t) ! adjust target for periodicity
           if( .not. (holes.chk.d) .and. .not. (trips.chk.d) ) then  ! if there are no holes at the destination...
            if( .not. (elecs.chk.d) ) then ! if there are no electrons at the destination...
!             if(i_am_root) call CPU_TIME(time(1))
!             holes%coul(j,me) = coulomb_energy_double(holes%N,elecs%N,d,holes%R,elecs%R,me) - oEc
!             holes%coul(j,me) = coulomb_difference_double_sparse(holes%N,elecs%N,o,d,holes%R,elecs%R,me)
             del_E = hopfld(j) + coulomb_energy_double(holes%N,elecs%N,d,holes%R,elecs%R,me) - oEc
!             if(i_am_root) call CPU_TIME(time(2))
             do i_state=0,holes%NE
              events       = events + 1 ! add an event
              delta_E      = holes%E(d%x,d%y,d%z,i_state) - holes%E(o%x,o%y,o%z,holes%iE(me)) + del_E
              rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
              holes%rate(i_state,j,me) = rate(events)
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 0
             enddo
!             if(i_am_root) then 
!              call CPU_TIME(time(3))
!              time_save_test(1) = time(2) - time(1) + time_save_test(1)
!              time_save_test(2) = time(3) - time(2) + time_save_test(2)
!             endif
            else
!             if(i_am_root) call CPU_TIME(time(1))
!             holes%coul(j,me) = coulomb_energy_double(holes%N,elecs%N,d,holes%R,elecs%R,me) - oEc
!             holes%coul(j,me) = coulomb_difference_double_sparse(holes%N,elecs%N,o,d,holes%R,elecs%R,me)
             del_E = hopfld(j) + coulomb_energy_double(holes%N,elecs%N,d,holes%R,elecs%R,me) - oEc
!             if(i_am_root) call CPU_TIME(time(2))
             do i_state=0,holes%NE
              events       = events + 1 ! add an event
              delta_E      = holes%E(d%x,d%y,d%z,i_state) - holes%E(o%x,o%y,o%z,holes%iE(me)) + del_E
              rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
              holes%rate(i_state,j,me) = rate(events)
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 2
              saj(events)  = (elecs.find.d)
             enddo         
!             if(i_am_root) then 
!              call CPU_TIME(time(3))
!              time_save_test(1) = time(2) - time(1) + time_save_test(1)
!              time_save_test(2) = time(3) - time(2) + time_save_test(2)
!             endif
            endif ! if there are no electrons at the destination
           else
            holes%rate(:,j,me) = 0.d0
           endif ! if there are no holes at the destination
          enddo ! loop over the hopbox
          holes%log_rate(me) = .false.
        else
          do j=1,hopboxN ! loop over the hopbox
           s = hopbox(j)        ! set up the step index 
           t = o + s            ! find the target point
           d = indx_periodic(t) ! adjust target for periodicity
           if( .not. (holes.chk.d) .and. .not. (trips.chk.d) ) then  ! if there are no holes at the destination...
            if( .not. (elecs.chk.d) ) then ! if there are no electrons at the destination...
!             if(i_am_root) call CPU_TIME(time(1))
!             del_E = hopfld(j) + holes%coul(j,me)
!             if(i_am_root) call CPU_TIME(time(2))
             do i_state=0,holes%NE
              events       = events + 1 ! add an event
              rate(events) = holes%rate(i_state,j,me)
!              delta_E      = holes%E(d%x,d%y,d%z,i_state) - holes%E(o%x,o%y,o%z,holes%iE(me)) + del_E
!              rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 0
             enddo
!             if(i_am_root) then 
!              call CPU_TIME(time(3))
!              time_save_test(1) = time(2) - time(1) + time_save_test(1)
!              time_save_test(2) = time(3) - time(2) + time_save_test(2)
!             endif
            else
!             if(i_am_root) call CPU_TIME(time(1))
!             del_E = hopfld(j) + holes%coul(j,me)
!             if(i_am_root) call CPU_TIME(time(2))
             do i_state=0,holes%NE
              events       = events + 1 ! add an event
              rate(events) = holes%rate(i_state,j,me)
!              delta_E      = holes%E(d%x,d%y,d%z,i_state) - holes%E(o%x,o%y,o%z,holes%iE(me)) + del_E
!              rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
              
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 2
              saj(events)  = (elecs.find.d)
             enddo         
!             if(i_am_root) then 
!              call CPU_TIME(time(3))
!              time_save_test(1) = time(2) - time(1) + time_save_test(1)
!              time_save_test(2) = time(3) - time(2) + time_save_test(2)
!             endif
            endif ! if there are no electrons at the destination
           endif ! if there are no holes at the destination
          enddo ! loop over the hopbox
        endif
       endif ! does the hole exist?
      enddo ! loop over holes
      
      do me=1,elecs%N
       if(elecs%exists(me)) then  ! does the electron exist?
        o = elecs%R(me)       ! Origin point
        if(elecs%iE(me).gt.0) then
         events       = events + 1
         rate(events) = elecs%relax_rate
         svo(events)  = o
         svd(events)  = o 
         iEd(events)  = 0
         sai(events)  = me
         p_i(events)  = 1
        endif
        if(elecs%log_rate(me)) then ! if the last hop has charged thing, recalculate the coulomb difference
          oEc = coulomb_energy_double(elecs%N,holes%N,o,elecs%R,holes%R,me)        
          do j=1,hopboxN ! loop over the hopbox
           s = hopbox(j)        ! set up the step index 
           t = o + s            ! find the target point
           d = indx_periodic(t) ! adjust target for periodicity
           if( .not. (elecs.chk.d) .and. .not. (trips.chk.d) ) then ! if there are no electrons at the destination...
            if( .not. (holes.chk.d) ) then  ! if there are no holes at the destination...
!             if(i_am_root) call CPU_TIME(time(1))
!             elecs%coul(j,me) = coulomb_energy_double(elecs%N,holes%N,d,elecs%R,holes%R,me) - oEc
!             elecs%coul(j,me) = coulomb_difference_double_sparse(elecs%N,holes%N,o,d,elecs%R,holes%R,me)
             del_E = - hopfld(j) + coulomb_energy_double(elecs%N,holes%N,d,elecs%R,holes%R,me) - oEc
!             if(i_am_root) call CPU_TIME(time(2))
             do i_state=0,elecs%NE
              events       = events + 1 ! add an event
              delta_E      = elecs%E(d%x,d%y,d%z,i_state) - elecs%E(o%x,o%y,o%z,elecs%iE(me)) + del_E
              rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
              elecs%rate(i_state,j,me) = rate(events)
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 1
             enddo
!             if(i_am_root) then 
!              call CPU_TIME(time(3))
!              time_save_test(1) = time(2) - time(1) + time_save_test(1)
!              time_save_test(2) = time(3) - time(2) + time_save_test(2)
!             endif
            else
!             if(i_am_root) call CPU_TIME(time(1))
!             elecs%coul(j,me) = coulomb_energy_double(elecs%N,holes%N,d,elecs%R,holes%R,me) - oEc
!             elecs%coul(j,me) = coulomb_difference_double_sparse(elecs%N,holes%N,o,d,elecs%R,holes%R,me)
             del_E = - hopfld(j) + coulomb_energy_double(elecs%N,holes%N,d,elecs%R,holes%R,me) - oEc
!             if(i_am_root) call CPU_TIME(time(2))
             do i_state=0,elecs%NE
              events       = events + 1 ! add an event
              delta_E      = elecs%E(d%x,d%y,d%z,i_state) - elecs%E(o%x,o%y,o%z,elecs%iE(me)) + del_E
              rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
              elecs%rate(i_state,j,me) = rate(events)
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 3
              saj(events)  = (holes.find.d)
             enddo
!             if(i_am_root) then 
!              call CPU_TIME(time(3))
!              time_save_test(1) = time(2) - time(1) + time_save_test(1)
!              time_save_test(2) = time(3) - time(2) + time_save_test(2)
!             endif
            endif ! if there are no holes at the destination
           else
            elecs%rate(:,j,me) = 0.d0
           endif ! if there are no electrons at the destination
          enddo ! loop over the hopbox
          elecs%log_rate(me) = .false.
        else
          do j=1,hopboxN ! loop over the hopbox
           s = hopbox(j)        ! set up the step index 
           t = o + s            ! find the target point
           d = indx_periodic(t) ! adjust target for periodicity
           if( .not. (elecs.chk.d) .and. .not. (trips.chk.d) ) then ! if there are no electrons at the destination...
            if( .not. (holes.chk.d) ) then  ! if there are no holes at the destination...
!             if(i_am_root) call CPU_TIME(time(1))
!             del_E = - hopfld(j) + elecs%coul(j,me)
!             if(i_am_root) call CPU_TIME(time(2))
             do i_state=0,elecs%NE
              events       = events + 1 ! add an event
              rate(events) = elecs%rate(i_state,j,me)
!              delta_E      = elecs%E(d%x,d%y,d%z,i_state) - elecs%E(o%x,o%y,o%z,elecs%iE(me)) + del_E
!              rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 1
             enddo
!             if(i_am_root) then 
!              call CPU_TIME(time(3))
!              time_save_test(1) = time(2) - time(1) + time_save_test(1)
!              time_save_test(2) = time(3) - time(2) + time_save_test(2)
!             endif
            else
!             if(i_am_root) call CPU_TIME(time(1))
!             del_E = - hopfld(j) + elecs%coul(j,me)
!             if(i_am_root) call CPU_TIME(time(2))
             do i_state=0,elecs%NE
              events       = events + 1 ! add an event
              rate(events) = elecs%rate(i_state,j,me)
!              delta_E      = elecs%E(d%x,d%y,d%z,i_state) - elecs%E(o%x,o%y,o%z,elecs%iE(me)) + del_E
!              rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 3
              saj(events)  = (holes.find.d)
             enddo
!             if(i_am_root) then 
!              call CPU_TIME(time(3))
!              time_save_test(1) = time(2) - time(1) + time_save_test(1)
!              time_save_test(2) = time(3) - time(2) + time_save_test(2)
!             endif
            endif ! if there are no holes at the destination
           endif ! if there are no electrons at the destination
          enddo ! loop over the hopbox
        endif
       endif ! does the electron exist?
      enddo ! loop over electrons
!------- nominate an event
      if(events.eq.0) stop 'ERROR: no events!!'
!------- Set up the Partial Sums
      sum_rate(1)  = rate(1)
      do i=2,events
       sum_rate(i) = sum_rate(i-1) + rate(i) !Partial sums
      enddo
!------- Normalize the Partial Sums
      rate_tot          = sum_rate(events)
!      sum_rate(:events) = sum_rate(:events) / rate_tot
!------- Pick an Event with MC
      call RANDOM_NUMBER(eta)
      event_keep = findeta(sum_rate(:events),eta * rate_tot)
!      do event_keep=1,events; if(sum_rate(event_keep).gt.eta) exit; enddo
!------- update the charge
      if(p_i(event_keep).eq.0) then
       holes%o           = svo(event_keep)
       holes%d           = svd(event_keep)
       holes%iE(sai(event_keep)) = iEd(event_keep)
       s                         = holes%hopid()
       call holes%update(sai(event_keep),hopdir(s%x,s%y,s%z))  
       holes%log_rate(sai(event_keep)) = .true.
       do i=1,holes%N
        if(i.ne.sai(event_keep)) holes%log_rate(i) = (holes%d.nearto.holes%R(i))
       enddo
      elseif(p_i(event_keep).eq.1) then
       elecs%o           = svo(event_keep)
       elecs%d           = svd(event_keep)
       elecs%iE(sai(event_keep)) = iEd(event_keep)
       s                         = elecs%hopid()
       call elecs%update(sai(event_keep),hopdir(s%x,s%y,s%z))   
       elecs%log_rate(sai(event_keep)) = .true.
       do i=1,elecs%N
        if(i.ne.sai(event_keep)) elecs%log_rate(i) = (elecs%d.nearto.elecs%R(i))
       enddo
      elseif(p_i(event_keep).eq.2) then
       call RANDOM_NUMBER(eta)
       holes%o           = svo(event_keep)
       holes%d           = svd(event_keep)
       holes%iE(sai(event_keep)) = iEd(event_keep)
       s                         = holes%hopid()
       holes%flux           = holes%flux + hopdir(s%x,s%y,s%z) 
       holes%R(sai(event_keep))  = holes%d
       call holes%update_flux(holes%o)
       holes%distri(holes%o%x,holes%o%y,holes%o%z) = .false.
       holes%distri(holes%d%x,holes%d%y,holes%d%z) = .true.
       holes%exists(sai(event_keep)) = .false.
       elecs%exists(saj(event_keep)) = .false.
       if(eta .ge. .75d0) then
        sings%N                       = sings%N + 1
        sings%exists(sings%N)         = .true.
        sings%R(sings%N)              = holes%d
        sings%iE(sings%N)             = sings%Eindex(holes%iE(sai(event_keep)),elecs%iE(saj(event_keep)))
       else
        trips%N                       = trips%N + 1
        trips%exists(trips%N)         = .true.
        trips%R(trips%N)              = holes%d
        trips%iE(trips%N)             = trips%Eindex(holes%iE(sai(event_keep)),elecs%iE(saj(event_keep)))
       endif
       holes%log_rate(sai(event_keep)) = .false.
       do i=1,holes%N
        if(i.ne.sai(event_keep)) holes%log_rate(i) = (holes%d.nearto.holes%R(i))
       enddo
       elecs%log_rate(saj(event_keep)) = .false.
       do i=1,elecs%N
        if(i.ne.saj(event_keep)) elecs%log_rate(i) = (holes%d.nearto.elecs%R(i))
       enddo
      elseif(p_i(event_keep).eq.3) then
       call RANDOM_NUMBER(eta)
       elecs%o           = svo(event_keep)
       elecs%d           = svd(event_keep)
       elecs%iE(sai(event_keep)) = iEd(event_keep)
       s                         = elecs%hopid()
       elecs%flux           = elecs%flux + hopdir(s%x,s%y,s%z) 
       elecs%R(sai(event_keep))  = elecs%d
       call elecs%update_flux(elecs%o)
       elecs%distri(elecs%o%x,elecs%o%y,elecs%o%z) = .false.
       elecs%distri(elecs%d%x,elecs%d%y,elecs%d%z) = .true.       
       elecs%exists(sai(event_keep)) = .false.
       holes%exists(saj(event_keep)) = .false.
       if(eta .ge. .75d0) then
        sings%N                       = sings%N + 1
        sings%exists(sings%N)         = .true.
        sings%R(sings%N)              = elecs%d
        sings%iE(sings%N)             = sings%Eindex(holes%iE(saj(event_keep)),elecs%iE(sai(event_keep)))
       else
        trips%N                       = trips%N + 1
        trips%exists(trips%N)         = .true.
        trips%R(trips%N)              = elecs%d
        trips%iE(trips%N)             = trips%Eindex(holes%iE(saj(event_keep)),elecs%iE(sai(event_keep)))
       endif
       elecs%log_rate(sai(event_keep)) = .false.
       do i=1,elecs%N
        if(i.ne.sai(event_keep)) elecs%log_rate(i) = (elecs%d.nearto.elecs%R(i))
       enddo
       holes%log_rate(saj(event_keep)) = .false.
       do i=1,holes%N
        if(i.ne.saj(event_keep)) holes%log_rate(i) = (elecs%d.nearto.holes%R(i))
       enddo
      endif
      endsubroutine calculate_events_double


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%% routine for calculating hopping event in a sparse two-charge guest-host system %%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine calculate_gh_events_double(holes,elecs,sings,trips,rate_tot)
      use constants
      use param
      use general_functions
      use my_mpi
      use pbc
      use coulomb
      use checks
      implicit none      
!----- input constant
      type(charge)                            :: holes, elecs       ! Hole and electron related data
      type(charge)                            :: sings, trips       ! Singlet and Triplet related data
!----- calc variables
      type(indexes)                           :: o                  ! origin point
      type(indexes)                           :: s                  ! step index
      type(indexes)                           :: t                  ! target point
      type(indexes)                           :: d                  ! destination point (target adjusted for periodicity)
      integer                                 :: me                 ! identity of the hopping charge
      integer                                 :: i, j               ! random indexes
      integer                                 :: guest_o, guest_d   ! guest indexes
      integer                                 :: events, event_keep ! indexes for events
      integer                                 :: i_state            ! index for state levels 
      double precision                        :: E_origin
      double precision                        :: delta_E            ! energy difference associated with hop
      double precision                        :: del_E              ! energy difference associated with hop - without multiple levels
      double precision                        :: oEc                ! Coloumb energy for the origin on a hopping event
      double precision                        :: eta                ! random number
      integer,          dimension(max_events) :: sai, saj           ! selection array: general indexes for all events
      integer,          dimension(max_events) :: p_i                ! selection array: general indexes for all events
      double precision, dimension(max_events) :: sum_rate           ! cummulative sum of the rates
!----- output variables
      double precision                        :: rate_tot           ! sum of all the rates
      events = 0
      do me=1,holes%N
       if(holes%exists(me)) then ! does the hole exist?
        o       = holes%R(me)              ! Origin point
        guest_o = holes%guest(o%x,o%y,o%z)
        if(holes%iE(me).gt.0) then
         events       = events + 1
         rate(events) = holes%relax_rate
         svo(events)  = o
         svd(events)  = o 
         iEd(events)  = 0
         sai(events)  = me
         p_i(events)  = 0
        endif
        if(holes%log_rate(me)) then ! if the last hop has charged things, recalculate the coulomb difference
          if(guest_o.eq.0) then
           E_origin = holes%E(o%x,o%y,o%z,holes%iE(me))
          else
           E_origin = holes%guest_E(holes%iE(me),guest_o)
          endif
          oEc = coulomb_energy_double(holes%N,elecs%N,o,holes%R,elecs%R,me)
          do j=1,hopboxN ! loop over the hopbox
           s = hopbox(j)        ! set up the step index 
           t = o + s            ! find the target point
           d = indx_periodic(t) ! adjust target for periodicity
           if( .not. (holes.chk.d) .and. .not. (trips.chk.d) .and. .not. (sings.chk.d) ) then  ! if there are no holes or excitons at the destination...
            guest_d = holes%guest(d%x,d%y,d%z)
            if( .not. (elecs.chk.d) ) then ! if there are no electrons at the destination...
             del_E   = hopfld(j) + coulomb_energy_double(holes%N,elecs%N,d,holes%R,elecs%R,me) - oEc
             if(guest_d.eq.0) then
              do i_state=0,holes%NE
               events       = events + 1 ! add an event
               delta_E      = holes%E(d%x,d%y,d%z,i_state) - E_origin + del_E
               rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
               holes%rate(i_state,j,me) = rate(events)
               ! Load event into the selection arrays
               sai(events)  = me
               svo(events)  = o
               svd(events)  = d 
               iEd(events)  = i_state
               p_i(events)  = 0
              enddo
             else
              do i_state=0,holes%guest_NE
               events       = events + 1 ! add an event
               delta_E      = holes%guest_E(i_state,guest_d) - E_origin + del_E
               rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
               holes%rate(i_state,j,me) = rate(events)
               ! Load event into the selection arrays
               sai(events)  = me
               svo(events)  = o
               svd(events)  = d 
               iEd(events)  = i_state
               p_i(events)  = 0
              enddo
             endif             
            elseif(guest_d.ne.0) then
             del_E = hopfld(j) + coulomb_energy_double(holes%N,elecs%N,d,holes%R,elecs%R,me) - oEc
             do i_state=0,holes%guest_NE
              events       = events + 1 ! add an event
              delta_E      = holes%guest_E(i_state,guest_d) - E_origin + del_E
              rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
              holes%rate(i_state,j,me) = rate(events)
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 2
              saj(events)  = (elecs.find.d)
             enddo
            endif ! if there are no electrons at the destination
           else
            holes%rate(:,j,me) = 0.d0
           endif ! if there are no holes at the destination
          enddo  ! loop over the hopbox
          holes%log_rate(me) = .false.
        else
          do j=1,hopboxN ! loop over the hopbox
           s = hopbox(j)        ! set up the step index 
           t = o + s            ! find the target point
           d = indx_periodic(t) ! adjust target for periodicity
           if( .not. (holes.chk.d) .and. .not. (sings.chk.d) .and. .not. (trips.chk.d) ) then  ! if there are no holes at the destination...
            guest_d = holes%guest(d%x,d%y,d%z)
            if( .not. (elecs.chk.d) ) then ! if there are no electrons at the destination...
             do i_state=0,holes%NE
              events       = events + 1 ! add an event
              rate(events) = holes%rate(i_state,j,me)
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 0
             enddo
            elseif(guest_d.ne.0) then
             do i_state=0,holes%NE
              events       = events + 1 ! add an event
              rate(events) = holes%rate(i_state,j,me)
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 2
              saj(events)  = (elecs.find.d)
             enddo         
            endif ! if there are no electrons at the destination
           endif ! if there are no holes at the destination
          enddo ! loop over the hopbox
        endif
       endif ! does the hole exist?
      enddo ! loop over holes
      
      do me=1,elecs%N
       if(elecs%exists(me)) then  ! does the electron exist?
        o = elecs%R(me)       ! Origin point
        guest_o = elecs%guest(o%x,o%y,o%z)
        if(elecs%iE(me).gt.0) then
         events       = events + 1
         rate(events) = elecs%relax_rate
         svo(events)  = o
         svd(events)  = o 
         iEd(events)  = 0
         sai(events)  = me
         p_i(events)  = 1
        endif
        if(elecs%log_rate(me)) then ! if the last hop has charged things, recalculate the coulomb difference
          if(guest_o.eq.0) then
           E_origin = elecs%E(o%x,o%y,o%z,elecs%iE(me))
          else
           E_origin = elecs%guest_E(elecs%iE(me),guest_o)
          endif
          oEc = coulomb_energy_double(elecs%N,holes%N,o,elecs%R,holes%R,me)        
          do j=1,hopboxN ! loop over the hopbox
           s = hopbox(j)        ! set up the step index 
           t = o + s            ! find the target point
           d = indx_periodic(t) ! adjust target for periodicity
           if( .not. (elecs.chk.d) .and. .not. (trips.chk.d) .and. .not. (sings.chk.d) ) then ! if there are no electrons at the destination...
            guest_d = elecs%guest(d%x,d%y,d%z)
            if( .not. (holes.chk.d) ) then  ! if there are no holes at the destination...
             del_E = - hopfld(j) + coulomb_energy_double(elecs%N,holes%N,d,elecs%R,holes%R,me) - oEc
             if(guest_d.eq.0) then
              do i_state=0,elecs%NE
               events       = events + 1 ! add an event
               delta_E      = elecs%E(d%x,d%y,d%z,i_state) - E_origin + del_E
               rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
               elecs%rate(i_state,j,me) = rate(events)
               ! Load event into the selection arrays
               sai(events)  = me
               svo(events)  = o
               svd(events)  = d 
               iEd(events)  = i_state
               p_i(events)  = 1
              enddo
             else
              do i_state=0,elecs%guest_NE
               events       = events + 1 ! add an event
               delta_E      = elecs%guest_E(i_state,guest_d) - E_origin + del_E
               rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
               elecs%rate(i_state,j,me) = rate(events)
               ! Load event into the selection arrays
               sai(events)  = me
               svo(events)  = o
               svd(events)  = d 
               iEd(events)  = i_state
               p_i(events)  = 1
              enddo
             endif
            elseif(guest_d.ne.0) then
             del_E = - hopfld(j) + coulomb_energy_double(elecs%N,holes%N,d,elecs%R,holes%R,me) - oEc
             do i_state=0,elecs%guest_NE
              events       = events + 1 ! add an event
              delta_E      = elecs%guest_E(i_state,guest_d) - E_origin + del_E
              rate(events) = erxyz(s%x,s%y,s%z) * myrate(delta_E)       ! Rate
              elecs%rate(i_state,j,me) = rate(events)
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 3
              saj(events)  = (holes.find.d)
             enddo
            endif ! if there are no holes at the destination
           else
            elecs%rate(:,j,me) = 0.d0
           endif ! if there are no electrons at the destination
          enddo ! loop over the hopbox
          elecs%log_rate(me) = .false.
        else
          do j=1,hopboxN ! loop over the hopbox
           s = hopbox(j)        ! set up the step index 
           t = o + s            ! find the target point
           d = indx_periodic(t) ! adjust target for periodicity
           if( .not. (elecs.chk.d) .and. .not. (sings.chk.d) .and. .not. (trips.chk.d) ) then ! if there are no electrons at the destination...
            guest_d = holes%guest(d%x,d%y,d%z)
            if( .not. (holes.chk.d) ) then  ! if there are no holes at the destination...
             do i_state=0,elecs%NE
              events       = events + 1 ! add an event
              rate(events) = elecs%rate(i_state,j,me)
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 1
             enddo
            elseif(guest_d.ne.0) then
             do i_state=0,elecs%NE
              events       = events + 1 ! add an event
              rate(events) = elecs%rate(i_state,j,me)
              ! Load event into the selection arrays
              sai(events)  = me
              svo(events)  = o
              svd(events)  = d 
              iEd(events)  = i_state
              p_i(events)  = 3
              saj(events)  = (holes.find.d)
             enddo
            endif ! if there are no holes at the destination
           endif ! if there are no electrons at the destination
          enddo ! loop over the hopbox
        endif
       endif ! does the electron exist?
      enddo ! loop over electrons
!------- nominate an event
      if(events.eq.0) stop 'ERROR: no events!!'
!------- Set up the Partial Sums
      sum_rate(1)  = rate(1)
      do i=2,events
       sum_rate(i) = sum_rate(i-1) + rate(i) !Partial sums
      enddo
!------- Normalize the Partial Sums
      rate_tot          = sum_rate(events)
!------- Pick an Event with MC
      call RANDOM_NUMBER(eta)
      event_keep = findeta(sum_rate(:events),eta * rate_tot)
!------- update the charge
      if(p_i(event_keep).eq.0) then
       holes%o           = svo(event_keep)
       holes%d           = svd(event_keep)
       holes%iE(sai(event_keep)) = iEd(event_keep)
       s                         = holes%hopid()
       call holes%update(sai(event_keep),hopdir(s%x,s%y,s%z))  
       holes%log_rate(sai(event_keep)) = .true.
       do i=1,holes%N
        if(i.ne.sai(event_keep)) holes%log_rate(i) = (holes%d.nearto.holes%R(i))
       enddo
      elseif(p_i(event_keep).eq.1) then
       elecs%o           = svo(event_keep)
       elecs%d           = svd(event_keep)
       elecs%iE(sai(event_keep)) = iEd(event_keep)
       s                         = elecs%hopid()
       call elecs%update(sai(event_keep),hopdir(s%x,s%y,s%z))   
       elecs%log_rate(sai(event_keep)) = .true.
       do i=1,elecs%N
        if(i.ne.sai(event_keep)) elecs%log_rate(i) = (elecs%d.nearto.elecs%R(i))
       enddo
      elseif(p_i(event_keep).eq.2) then
       call RANDOM_NUMBER(eta)
       holes%o           = svo(event_keep)
       holes%d           = svd(event_keep)
       holes%iE(sai(event_keep)) = iEd(event_keep)
       s                         = holes%hopid()
       holes%flux           = holes%flux + hopdir(s%x,s%y,s%z) 
       holes%R(sai(event_keep))  = holes%d
       call holes%update_flux(holes%o)
       holes%distri(holes%o%x,holes%o%y,holes%o%z) = .false.
       holes%distri(holes%d%x,holes%d%y,holes%d%z) = .true.
       holes%exists(sai(event_keep)) = .false.
       elecs%exists(saj(event_keep)) = .false.
       if(eta .ge. .75d0) then
        sings%N                       = sings%N + 1
        sings%exists(sings%N)         = .true.
        sings%R(sings%N)              = holes%d
        sings%iE(sings%N)             = sings%Eindex(holes%iE(sai(event_keep)),elecs%iE(saj(event_keep)))
       else
        trips%N                       = trips%N + 1
        trips%exists(trips%N)         = .true.
        trips%R(trips%N)              = holes%d
        trips%iE(trips%N)             = trips%Eindex(holes%iE(sai(event_keep)),elecs%iE(saj(event_keep)))
       endif
       holes%log_rate(sai(event_keep)) = .false.
       do i=1,holes%N
        if(i.ne.sai(event_keep)) holes%log_rate(i) = (holes%d.nearto.holes%R(i))
       enddo
       elecs%log_rate(saj(event_keep)) = .false.
       do i=1,elecs%N
        if(i.ne.saj(event_keep)) elecs%log_rate(i) = (holes%d.nearto.elecs%R(i))
       enddo
      elseif(p_i(event_keep).eq.3) then
       call RANDOM_NUMBER(eta)
       elecs%o           = svo(event_keep)
       elecs%d           = svd(event_keep)
       elecs%iE(sai(event_keep)) = iEd(event_keep)
       s                         = elecs%hopid()
       elecs%flux           = elecs%flux + hopdir(s%x,s%y,s%z) 
       elecs%R(sai(event_keep))  = elecs%d
       call elecs%update_flux(elecs%o)
       elecs%distri(elecs%o%x,elecs%o%y,elecs%o%z) = .false.
       elecs%distri(elecs%d%x,elecs%d%y,elecs%d%z) = .true.       
       elecs%exists(sai(event_keep)) = .false.
       holes%exists(saj(event_keep)) = .false.
       if(eta .ge. .75d0) then
        sings%N                       = sings%N + 1
        sings%exists(sings%N)         = .true.
        sings%R(sings%N)              = elecs%d
        sings%iE(sings%N)             = sings%Eindex(holes%iE(saj(event_keep)),elecs%iE(sai(event_keep)))
       else
        trips%N                       = trips%N + 1
        trips%exists(trips%N)         = .true.
        trips%R(trips%N)              = elecs%d
        trips%iE(trips%N)             = trips%Eindex(holes%iE(saj(event_keep)),elecs%iE(sai(event_keep)))
       endif
       elecs%log_rate(sai(event_keep)) = .false.
       do i=1,elecs%N
        if(i.ne.sai(event_keep)) elecs%log_rate(i) = (elecs%d.nearto.elecs%R(i))
       enddo
       holes%log_rate(saj(event_keep)) = .false.
       do i=1,holes%N
        if(i.ne.saj(event_keep)) holes%log_rate(i) = (elecs%d.nearto.holes%R(i))
       enddo
      endif
      endsubroutine calculate_gh_events_double

endmodule event






































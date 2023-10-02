!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% program CharTED-KMC.f90: ...description...
!%
!% Copyright (C) 2020 Yvelin Giret & Thomas Pope
!% \authors Yvelin Giret and Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    program kmc_main
    use type_charge
    use my_mpi
    use constants
    use param
    use general_functions
    use timings
    use printings
    use initialization
    use pbc
    use coulomb
    use random
    use checks
    use event
    implicit none
    type(charge)                                    :: holes, elecs
    type(charge)                                    :: sings, trips
    double precision                                :: rate_tot, giga_rho
    logical                                         :: stop_run, simulation
    integer, allocatable, dimension(:)              :: seed
    integer                                         :: nseed
    integer                                         :: i
    double precision, allocatable, dimension(:,:,:) :: noise_map
    double precision                                :: t3,t4

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%% Initialisation %%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call mpi_initialize()
    simulation = .true.
    if(i_am_root) then   ! head process
      call print_logo()  ! Print the logo
      call timing_init() ! Timing initialisation
      time_save_test=0.d0
    endif
    call initialize_input(holes,elecs) ! Read input file
    if(inc_hole.and..not.inc_elec.and.i_am_root) then
     giga_rho = dble(N_hole) / dble((Nx*Ny*Nz))
     write(screen,fmchfl) "Expected Mem Usage [GB/cpu]: ", &
     &(64.d0+N_HOMO*8.d0+(8.d0*(N_HOMO+2)+(24.d0+2080.d0*(N_HOMO+1))*giga_rho)*Nx*Ny*Nz)/(1024.d0**3)
    endif
    if(.not.inc_hole) holes%N = 0
    if(.not.inc_elec) elecs%N = 0
    if(inc_excn) then
     if(guest_host_flag) then
      if(inc_hole) call holes%guest_aloc(guest_density,GUEST_N_HOMO,GUEST_E_HOMO,sigma) 
      if(inc_elec) then
       if(inc_hole) then 
        call elecs%guest_aloc(guest_density,GUEST_N_LUMO,GUEST_E_LUMO,sigma,holes%guest) 
       else
        call elecs%guest_aloc(guest_density,GUEST_N_LUMO,GUEST_E_LUMO,sigma) 
       endif
      endif
      call sings%Eindex_init(Nx,Ny,Nz,GUEST_N_HOMO,holes%guest_E0,GUEST_N_LUMO,elecs%guest_E0,max(holes%N,elecs%N))   
      call trips%Eindex_init(Nx,Ny,Nz,GUEST_N_HOMO,holes%guest_E0,GUEST_N_LUMO,elecs%guest_E0,max(holes%N,elecs%N))   
     else
      call sings%Eindex_init(Nx,Ny,Nz,N_HOMO,holes%E0,N_LUMO,elecs%E0,max(holes%N,elecs%N))   
      call trips%Eindex_init(Nx,Ny,Nz,N_HOMO,holes%E0,N_LUMO,elecs%E0,max(holes%N,elecs%N))   
     endif
    endif
    if(debugging) then !   stop the random numbers from being random (only use for debugging)
      call random_seed(size = nseed); allocate(seed(nseed)) ! allocate seed array
      seed(:) = 1; call random_seed(put = seed)             ! set ALL seeds to 1
    endif
    if(.not.restart) then ! head process
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%%%% Disorder init.              %%%%%%%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(i_am_root) then
       call CPU_TIME(t3)
       write(screen,fmdlm1) blank1
       write(screen,fmchr1) 'Disorder:'
       write(screen,fmdlm2) delim2, blank1
      endif
      allocate(noise_map(Nx,Ny,Nz))
      noise_map = 0.0d0
      if(NN_corr.gt.0) then
       if(dimer_flag) then
        noise_map = noise_map + correlated_dimers(NN_corr,Nx,Ny,Nz,dipole,lattice,epsilon_r,dimer_dipole,dimer_rho)
       else
        noise_map = noise_map + correlated(NN_corr,Nx,Ny,Nz,dipole,lattice,epsilon_r,order_parameter) 
       endif
      endif
      if(inc_hole) then
       do i=0,holes%NE
        if(i_am_root) write(screen,fmchin) "HOMO Level ", i
        if(sigma.gt.0.d0) then
          holes%E(:,:,:,i) = holes%E(:,:,:,i) + noise_map + uncorrelated(Nx,Ny,Nz,sigma)
        else
          holes%E(:,:,:,i) = holes%E(:,:,:,i) + noise_map
        endif
       enddo
      endif
      if(inc_elec) then
       do i=0,elecs%NE
        if(i_am_root) write(screen,fmchin) "LUMO Level ", i
        if(sigma.gt.0.d0) then
          elecs%E(:,:,:,i) = elecs%E(:,:,:,i) + noise_map + uncorrelated(Nx,Ny,Nz,sigma)
        else
          elecs%E(:,:,:,i) = elecs%E(:,:,:,i) + noise_map
        endif
       enddo
      endif   
      deallocate(noise_map) 
      if(i_am_root) then
       call CPU_TIME(t4)
       write(screen,fmdlm1) blank1
       write(screen,fmcflc) "Distribution(s) generated in", t4-t3, "seconds"
       write(screen,fmdlm3) blank1, delim1, blank1
      endif
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%% Hole distribution initialisation %%%%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(i_am_root) call CPU_TIME(t3)
      if(inc_hole) then
        if(i_am_root) then
          write(screen,fmchr1) 'Hole distribution initialisation:'
          write(screen,fmdlm2) delim2, blank1
          flush(screen)
        endif
        call charge_init(holes,E_fermi,kB_in_AU*temperature,"HOMO")
      endif
      if(inc_hole.and.inc_elec.and.i_am_root) write(screen,fmdlm1) delim2
      if(inc_elec) then
        if(i_am_root) then
          write(screen,fmchr1) 'Electron distribution initialisation:'
          write(screen,fmdlm2) delim2, blank1
          flush(screen)
        endif
        if(inc_hole) then
          call charge_init(elecs,E_fermi,kB_in_AU*temperature,"LUMO",holes%distri)
        else
          call charge_init(elecs,E_fermi,kB_in_AU*temperature,"LUMO")
        endif
      endif
      if(i_am_root) then
       call CPU_TIME(t4)
       write(screen,fmdlm1) blank1
       write(screen,fmcflc) "Distribution(s) generated in", t4-t3, "seconds"
       write(screen,fmdlm2) blank1, delim1
      endif
    endif
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%% Set up the electric field vector %%%%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    if(i_am_root) then
      write(screen,fmdlm1) blank1
      write(screen,fmchrv) 'Electric field along v=', field_vector
      write(screen,fmdlm1) delim2, blank1
      write(screen,fmches) 'Applied E  [V/nm]             :', field /Vnm_to_AU
      write(screen,fmchfl) 'Field step [eV]               :', field_step * Eh_to_eV
      write(screen,fmdlm1) blank1
    endif

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%%%%%%%%%%%%%% Open outputs %%%%%%%%%%%%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call open_file(111, "flux",me_ip)
    call open_file(112, "current",me_ip)
    call open_file(114, "mobility_a",me_ip)
    call open_file(115, "mobility_b",me_ip)
    if (inc_excn) call open_file(117, "emissions",me_ip)
    call open_file(118, "trajectory",me_ip,".xyz")
    if(.not. restart) then
      time_tot = 0.0d0 ; iteration = 0
    end if
    
    ! initialize KMC loop
    call init_coloumb_matrix()
    call init_events(holes,elecs)

    ! kmc loop:
    if(i_am_root) then
     write(screen,fmdlm1) blank1
     write(screen,fmchr1) 'KMC loop...'
     write(screen,fmdlm2) delim2, blank1
     flush(screen)
     call initialize_print_timer()
    endif
    call initialize_print_converters(holes)
    do while(simulation)
      call read_stop(stop_run)                          ! check if stopfile is active
      iteration = iteration + 1
      if(inc_hole) call holes%check_mat(iteration,me_ip,"HOLE",0)      ! check that the number of holes is conserved
      if(inc_elec) call elecs%check_mat(iteration,me_ip,"ELEC",0)      ! check that the number of electrons is conserved
      call calculate_events(holes,elecs,sings,trips,rate_tot) ! select an event for the iteration
      time_tot = time_tot - log(myrandom()) / rate_tot      ! Add the randomized hopping time to the time total
      if((holes%N+elecs%N).eq.0 .or.  iteration.eq.Nsteps) stop_run = .true.
      if(mod(iteration,modulo_print).eq.0  .or. stop_run) then
       if(inc_hole.and..not.inc_elec) then
        call print_step_output(time_tot,iteration,holes)
       elseif(inc_elec.and..not.inc_hole) then
        call print_step_output(time_tot,iteration,elecs)
       else
!        if(i_am_root) write(*,*) time_save_test(1:2); flush(6)
        call print_step_output(time_tot,iteration,holes,elecs,sings,trips)
       endif
      endif
      if(mod(iteration,modulo_print_xyz).eq.0 .or. iteration.eq.Nsteps .or. stop_run) then
       if(inc_excn) then
        call print_trajectory2(118,holes,elecs,sings,trips,iteration,time_tot)
        call print_exciton_checks(sings,"singlet-analysis.dat")
        call print_exciton_checks(trips,"triplet-analysis.dat")
       else
        if(inc_hole) then
         call print_trajectory1(118,holes,iteration,time_tot)
        elseif(inc_elec) then
         call print_trajectory1(118,elecs,iteration,time_tot)
        endif
       endif
      endif
        ! If we've hit the end of the calculation, or if we've been told to, stop
      if((iteration.eq.Nsteps) .or. stop_run) simulation = .false.
    enddo
    call mpi_barrier(we_comm, ierr)
    if(inc_excn) then
     call print_exciton_checks(sings,"singlet-analysis.dat")
     call print_exciton_checks(trips,"triplet-analysis.dat")
    endif
    close(111); close(112); close(114); close(115); if(inc_excn) close(117); close(118)
    if(i_am_root.and..not.inc_excn) call print_mobility_analysis()
    if(inc_excn) call print_exciton_analysis(sings,trips)    
    if(inc_hole) then
     deallocate(holes%distri) ; deallocate(holes%E)
    endif
    if(inc_elec) then
     deallocate(elecs%distri) ; deallocate(elecs%E)
    endif
    ! Timing:
    if(i_am_root) then
      close(211)
      close(212)
      call timing_end()
    endif
    call mpi_finalize(ierr)
    end program kmc_main
    


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% module mpi.f: set of routines to handle all 
!%            message passing and job allocations
!%            in the parrallel environment
!%
!% Copyright (C) 2020 Thomas Pope
!% \author Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    module my_mpi
    use mpi
    use kinds
    implicit none
    integer :: ierr         ! Error flag used in message passing
    integer :: we_comm      ! comm world
    integer :: we_np        ! number of overall processes 
    integer :: me_ip        ! global index of the current process
    logical :: i_am_root    ! logical flag, only true on the root node
    integer,dimension(mpi_status_size) :: istat ! status parameter
    integer :: we_real      ! size of the real numbers on disk
    contains
!---------------------------------------------------------------------!
! mpi_initialize:                                                     !
!      initialize the mpi enviroment and tell the root node it's root !
!---------------------------------------------------------------------!
    subroutine mpi_initialize()
!   routine to initialize the mpi enviroment
    implicit none
    call mpi_init(ierr)
    we_comm = mpi_comm_world
    call mpi_comm_size(we_comm,we_np,ierr)
    call mpi_comm_rank(we_comm,me_ip,ierr)
    if(me_ip.eq.0) then
      i_am_root = .true.
    else
      i_am_root = .false.
    endif
    endsubroutine mpi_initialize
!---------------------------------------------------------------------!
! mpi_bcast_parameters:                                               !
!      broadcast the input parameters to all nodes                    !
!---------------------------------------------------------------------!
    subroutine mpi_bcast_parameters()
    use param
    implicit none
    we_real=mpi_double_precision
    !logic flags
    call mpi_bcast(inc_hole,          1,mpi_logical,0,we_comm,ierr)
    call mpi_bcast(inc_elec,          1,mpi_logical,0,we_comm,ierr)
    call mpi_bcast(inc_excn,          1,mpi_logical,0,we_comm,ierr)
    call mpi_bcast(debugging,         1,mpi_logical,0,we_comm,ierr)
    call mpi_bcast(restart,           1,mpi_logical,0,we_comm,ierr)
    call mpi_bcast(energy_init,       1,mpi_logical,0,we_comm,ierr)
    call mpi_bcast(distri_init,       1,mpi_logical,0,we_comm,ierr)
    call mpi_bcast(mill_rate,         1,mpi_logical,0,we_comm,ierr)
    call mpi_bcast(mrcs_rate,         1,mpi_logical,0,we_comm,ierr)
    call mpi_bcast(dimer_flag,        1,mpi_logical,0,we_comm,ierr)
    call mpi_bcast(guest_host_flag,   1,mpi_logical,0,we_comm,ierr)
    call mpi_bcast(test_the_new_thing,1,mpi_logical,0,we_comm,ierr)
    !printing flags
    call mpi_bcast(modulo_print,      1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(modulo_print_xyz,  1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(Nsteps,            1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(N_equilibration,   1,mpi_integer,0,we_comm,ierr)
    !array dimensions
    call mpi_bcast(Nx,                1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(Ny,                1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(Nz,                1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(N_hopbox,          1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(N_hole,            1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(N_elec,            1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(N_HOMO,            1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(N_LUMO,            1,mpi_integer,0,we_comm,ierr)
    if(guest_host_flag) then
     call mpi_bcast(GUEST_N_HOMO,     1,mpi_integer,0,we_comm,ierr)
     call mpi_bcast(GUEST_N_LUMO,     1,mpi_integer,0,we_comm,ierr)
    endif
    call mpi_bcast(N_chrg,            1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(NN_corr,           1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(N_SING,            1,mpi_integer,0,we_comm,ierr)
    call mpi_bcast(N_TRIP,            1,mpi_integer,0,we_comm,ierr)
    if(.not.i_am_root) then
     allocate(E_HOMO(0:N_HOMO),E_LUMO(0:N_LUMO))
     if(guest_host_flag) allocate(GUEST_E_HOMO(0:GUEST_N_HOMO),GUEST_E_LUMO(0:GUEST_N_LUMO))
    endif
    !system parameters
    call mpi_bcast(field,             1,we_real,0,we_comm,ierr)
    call mpi_bcast(field_step,        1,we_real,0,we_comm,ierr)
    call mpi_bcast(field_vector,      3,we_real,0,we_comm,ierr)
    call mpi_bcast(coulomb_factor,    1,we_real,0,we_comm,ierr)
    call mpi_bcast(lattice,           1,we_real,0,we_comm,ierr)
    call mpi_bcast(local,             1,we_real,0,we_comm,ierr)
    call mpi_bcast(hop_prefactor,     1,we_real,0,we_comm,ierr)
    call mpi_bcast(alpha,             1,we_real,0,we_comm,ierr)
    call mpi_bcast(lambda,            1,we_real,0,we_comm,ierr)
    call mpi_bcast(temperature,       1,we_real,0,we_comm,ierr)
    call mpi_bcast(sigma,             1,we_real,0,we_comm,ierr)
    call mpi_bcast(dipole,            1,we_real,0,we_comm,ierr)
    call mpi_bcast(epsilon_r,         1,we_real,0,we_comm,ierr)
    call mpi_bcast(order_parameter,   1,we_real,0,we_comm,ierr)
    call mpi_bcast(dimer_dipole,      1,we_real,0,we_comm,ierr)
    call mpi_bcast(dimer_rho,         1,we_real,0,we_comm,ierr)
    call mpi_bcast(E_fermi,           1,we_real,0,we_comm,ierr)
    call mpi_bcast(E_HOMO,     N_HOMO+1,we_real,0,we_comm,ierr)
    call mpi_bcast(E_LUMO,     N_LUMO+1,we_real,0,we_comm,ierr)
    call mpi_bcast(density,           1,we_real,0,we_comm,ierr)
    call mpi_bcast(time_tot,          1,we_real,0,we_comm,ierr)
    call mpi_bcast(density_real,      1,we_real,0,we_comm,ierr)
    if(guest_host_flag) then
     call mpi_bcast(GUEST_E_HOMO, GUEST_N_HOMO+1,we_real,0,we_comm,ierr)
     call mpi_bcast(GUEST_E_LUMO, GUEST_N_LUMO+1,we_real,0,we_comm,ierr)
     call mpi_bcast(guest_density,             1,we_real,0,we_comm,ierr)    
    endif
    if(inc_hole) call mpi_bcast(hole_hop,1,we_real,0,we_comm,ierr)
    if(inc_elec) call mpi_bcast(elec_hop,1,we_real,0,we_comm,ierr)
    if(inc_hole) call mpi_bcast(hole_rel,1,we_real,0,we_comm,ierr)
    if(inc_elec) call mpi_bcast(elec_rel,1,we_real,0,we_comm,ierr)
    endsubroutine mpi_bcast_parameters
!---------------------------------------------------------------------!
    endmodule my_mpi

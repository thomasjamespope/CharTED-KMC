!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% module charge_distri_xyz.f90: ...description...
!%
!% Copyright (C) 2020 Yvelin Giret & Thomas Pope
!% \author Yvelin Giret & Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module initialization
     use type_charge
     use constants
     use param
     use readings
     implicit none
     

     contains

     subroutine initialize_input(holes,elecs)
     use printings
     use my_mpi
     use checks
     implicit none
     type(charge)   :: holes, elecs!, sings, trips
     character(ch)  :: fname
     logical        :: i_am_insane
     if(i_am_root) then   ! head process
       call GETARG(1,fname)
       call inquire_inputs(fname,restart)
     endif
     call mpi_bcast(restart,1,mpi_logical,0,we_comm,ierr)
     if(i_am_root) then
       call read_input(fname,i_am_insane)
       if(i_am_insane) then
        write(screen,fmdlm1) delim3
        write(screen,fmchr1) "You have requested insanity... Good luck!"
        write(screen,fmdlm1) delim3
       else
        call sanity_checks()
       endif
       call convert_AU()           ! Convert to AU units
       call check_site_number(Nx,Ny,Nz,correlated_type)
     endif
     call mpi_bcast_parameters()
     max_events = 26 * N_chrg * (max(N_HOMO,N_LUMO) + 1) * 2
     if(inc_excn) then
      call allocate_print_arrays(2)
     else
      call allocate_print_arrays(1)
     endif
     ! allocate charge arrays on all nodes
     if(inc_hole) call holes%aloc(Nx,Ny,Nz,N_HOMO,E_HOMO,hole_hop,hole_rel,N_hole)
     if(inc_elec) call elecs%aloc(Nx,Ny,Nz,N_LUMO,E_LUMO,elec_hop,hole_rel,N_elec)
    endsubroutine initialize_input



     subroutine read_input(fname,i_am_insane)
     use my_mpi
     implicit none
     character(ch)      :: fname
     ! Local variable:
     character(ch)      :: E_init, C_init, dum_str1, dum_str2
     character(90)      :: char_r
     logical            :: i_am_insane
     integer, parameter :: u=100
     integer            :: i_state
     write(screen,fmdlm1) blank1
     write(screen,fmchr2)  'Reading input file:', TRIM(fname)
     write(screen,fmdlm2) delim2, blank1
!#####################################################################!
!#### Set Defaults ###################################################!
!#####################################################################!
     energy_init     = .false.
     distri_init     = .false.
     debugging       = .false.
     inc_elec        = .false.
     inc_hole        = .false.
     dimer_flag      = .false.
     dimer_rho       = 0.d0
     dimer_dipole    = 0.d0
     N_elec          = 0
     N_hole          = 0
     N_HOMO          = 0
     N_LUMO          = 0
     N_SING          = 0
     N_TRIP          = 0
     N_hopbox        = 1
     MPI_N_calc      = we_np
     i_am_insane     = .false.
     mill_rate       = .false.
     mrcs_rate       = .false.
     guest_host_flag = .false.
     GUEST_N_HOMO    = 0
     GUEST_N_LUMO    = 0
     test_the_new_thing = .false.

!#####################################################################!
!#### Open the Input File   ##########################################!
!#####################################################################! 
     open(u,file=fname)

!#####################################################################!
!#### Read statements ##########################################!
!#####################################################################! 
     i_am_insane        = read_statements(u,"I am aware that this is insane")
     test_the_new_thing = read_statements(u,"test the new thing")
     guest_host_flag    = read_statements(u,"GUEST-HOST")
    
!#####################################################################!
!#### Read string variables ##########################################!
!#####################################################################! 
     write(screen,fmdlm1) blank1
     write(screen,fmchr1) "String Variables:"
     write(screen,fmdlm1) blank1
     
!     call read_1_string(u,E_init,         'ENERGY_INIT',    'Method for initializing site energies' ,      'CALC','READ')
!     call read_1_string(u,C_init,         'CHARGE_INIT',    'Method for initializing charge distribution' ,'CALC','READ')
     
     call read_1_string(u,run_type,       'RUN_TYPE',       'Calculation type' , 'HOLE',  'ELEC',   'FULL')
     call read_1_string(u,rate_type,      'CHARGE_RATE',    'Rate constant type' , 'MILLER-ABRAHAMS', 'MARCUS')
     call read_1_string(u,correlated_type,'CORRELATED_TYPE','Correlated disorder scheme' , 'RANDOM','ORDERED','MIXED')
     call read_1_integer(u,print_verb,    'VERBOSITY',      'Amount printed to output' , 1) 
     
!     if(E_init == 'READ') energy_init = .true.
!     if(C_init == 'READ') distri_init = .true.
     if(run_type == 'FULL') then 
      inc_elec = .true.; inc_hole = .true.; inc_excn=.true.
     elseif(run_type == 'HOLE') then 
      inc_hole = .true.
     elseif(run_type == 'ELEC') then
      inc_elec = .true.
     endif
     if(rate_type == 'MARCUS') then
       mrcs_rate = .true.
     elseif(rate_type == 'MILLER-ABRAHAMS') then
       mill_rate = .true.
     endif
     
!#####################################################################!
!#### Read remaining simulation: format variables ####################!
!#####################################################################! 
     write(screen,fmdlm1) blank1
     write(screen,fmchr1) "Format Variables:"
     write(screen,fmdlm1) blank1

     call read_1_logical(u,debugging,        'DEBUG',            'Inhibit RNG for debugging?',                 .false.)
     call read_1_integer(u,Nsteps,           'NUMBER_STEPS',     'Number of steps',                            1000000)
     if(Nsteps.lt.-1)          call wrong_value_1_int(Nsteps,          'NUMBER_STEPS',    1000000)
     call read_1_integer(u,modulo_print,     'N_PRINT_STEP',     'Number of steps between standard prints',    1000)
     if(modulo_print.le.0)     call wrong_value_1_int(modulo_print,    'N_PRINT_STEP',    1000)
     call read_1_integer(u,modulo_print_xyz, 'N_PRINT_STEP_XYZ', 'Number of steps between large-array prints', 10000)
     if(modulo_print_xyz.le.0) call wrong_value_1_int(modulo_print_xyz,'N_PRINT_STEP_XYZ',10000)

     call read_1_integer(u,N_equilibration, 'N_EQUILIBRATION', 'Sample Size for most-recent results calculation',  2 * modulo_print)

!#####################################################################!
!#### Read simulation: convergence: variables ########################!
!#####################################################################! 
     write(screen,fmdlm1) blank1
     write(screen,fmchr1) "Convergence Variables:"
     write(screen,fmdlm1) blank1

     call read_3_integers(u,Nx,Ny,Nz,'DIMENSIONS','Number of voxels (X,Y,Z)',50,50,50)
     if(Nx.le.0) call wrong_value_1_int(Nx,'Nx',50)
     if(Ny.le.0) call wrong_value_1_int(Ny,'Ny',50)
     if(Nz.le.0) call wrong_value_1_int(Nz,'Nz',50)
     write(char_r,'(i0)') Nx*Ny*Nz
     write(screen,rdchr2) 'Total number of voxels', adjustl(char_r)
     call read_1_integer(u,NN_corr,'CORRELATED_NEIGHBORS','Number neighbors for correlated disorder',7)
     if(NN_corr.lt.0) call wrong_value_1_int(NN_corr,'CORRELATED_NEIGHBORS',7)
     !call read_1_integer(u,N_hopbox,'HOPBOX','Number of hop-viable neighbors',1)
     N_hopbox = 1
     
!#####################################################################!
!#### Read Physical - system variables ###############################!
!#####################################################################! 
     write(screen,fmdlm1) blank1
     write(screen,fmchr1) "Physical Variables:"
     write(screen,fmdlm1) blank1

     call read_1_real(u,lattice,'LATTICE','Lattice parameter [nm]',1.0d0,"(f10.5)")
     if(lattice.le.0) call wrong_value_1_real(lattice,'LATTICE',1.0d0,erchfl)
     call read_1_real(u,sigma,'UNCORRELATED_DISORDER','Uncorrelated disorder strength [eV]', 0.05d0,"(f10.5)")
     if(sigma.lt.0) call wrong_value_1_real(sigma,'UNCORRELATED_DISORDER',0.05d0,erchfl)
     call read_1_real(u,dipole,'DIPOLE_MOMENT','Average dipole moment [D]', 0.0d0,"(f10.5)")
     if(dipole.lt.0) call wrong_value_1_real(dipole,'DIPOLE_MOMENT',0.0d0,erchfl)
     if(trim(correlated_type).eq.'MIXED') then
       call read_1_real(u,order_parameter,'ORDER_PARAMETER','Order parameter',0.0d0,"(f10.5)")
     elseif(trim(correlated_type).eq.'RANDOM') then
       order_parameter = 0.0d0
     elseif(trim(correlated_type).eq.'ORDERED') then
       order_parameter = 1.0d0
     endif
     if(order_parameter.lt.0.0d0.or.order_parameter.gt.1.0d0) &
      & call wrong_value_1_real(order_parameter,'ORDER_PARAMETER',0.0d0,erchfl)

     call read_1_real(u,dimer_rho,'DIMER_DENSITY','Dimer density [/site]', 0.0d0,"(f10.5)")
     if(dimer_rho.gt.0.0d0) dimer_flag=.true.
     call read_1_real(u,dimer_dipole,'DIMER_DIPOLE','Dipole of the dimers [D]', 0.0d0,"(f10.5)")
     if(mrcs_rate) then 
       call read_1_real(u,lambda,'REORGANIZATION','Reorganization energy [eV]', 0.1d0,"(f10.5)")
       if(lambda.le.0) call wrong_value_1_real(lambda,'REORGANIZATION',0.1d0,erchfl)
     endif
     call read_1_real(u,coul_cutoff,'COUL_CUTOFF','Cutoff radius for the coulomb memory algo.', 15.0d0,"(f10.5)")
     
     call read_1_real(u,local,'LOCAL','Charge Localization Factor', 2.0d0,"(f10.5)")
     call read_1_real(u,E_fermi,'FERMI_ENERGY','Fermi Energy [eV]',0.0d0,"(f10.5)")
     call read_1_real(u,field,'ELECTRIC_FIELD','Electric field [V/nm]',0.05d0,"(f10.5)")
!     call read_3_real(u,field_vector,'ELECTRIC_FIELD_VECTOR','Electric field vector',(/0.0d0,1.0d0,0.0d0/))
     field_vector=(/0.0d0,1.0d0,0.0d0/)
     call read_1_real(u,temperature,'TEMPERATURE','Temperature [K]',298.0d0,"(f10.5)")
     if(temperature.le.0) call wrong_value_1_real(temperature,'TEMPERATURE',298.0d0,erchfl)
     call read_1_real(u,epsilon_r,'PERMITTIVITY','Relative permittivity [e_0]', 3.5d0,"(f10.5)")
     if(epsilon_r.le.0) call wrong_value_1_real(epsilon_r,'PERMITTIVITY',3.5d0,erchfl)
     
!#####################################################################!
!#### Read Physical - hole variables #################################!
!#####################################################################!

     if(inc_hole) then
       write(screen,fmdlm1) blank1
       write(screen,fmchr1) "Hole-specific Variables:"
       write(screen,fmdlm1) blank1
       N_HOMO=read_count_states(u,"HOMO")
       allocate(E_HOMO(0:N_HOMO))
       call read_1_real(u,E_HOMO(0),'HOMO_ENERGY','HOMO Energy [eV]',0.0d0,"(f10.5)")
       do i_state = 1, N_HOMO
         write(dum_str1,'("HOMO+",i0,"_ENERGY")') i_state
         write(dum_str2,'("HOMO+",i0," Energy [eV]")') i_state
         call read_1_real(u,E_HOMO(i_state),trim(adjustl(dum_str1)),trim(adjustl(dum_str2)),0.0d0,"(f10.5)")
       enddo
       E_HOMO = - E_HOMO * eV_to_Eh
       call read_1_real(u,density,'HOLE_DENSITY','Hole density [/site]',0.001d0,"(f10.5)")
       if(density.le.0) call wrong_value_1_real(density,'HOLE_DENSITY',0.001d0,erchfl)
       N_hole = nint(density*dble(Nx*Ny*Nz))
       if(N_hole .LT. 1) then
         write(screen,fmdlm1) delim3
         write(screen,fmchr1) 'WARNING:'
         write(screen,fmchr1) 'Simu. cell too small to contain holes.'
         write(screen,fmchr1) 'ACTION: increase the number of voxels or the density.'
         write(screen,fmchr1) 'Number of holes fixed to 1.'
         N_hole = 1
         write(char_r,'(f10.5)') dble(N_hole) / dble(Nx*Ny*Nz) 
         write(screen,rdchr2) 'Hole density [/site]', adjustl(char_r)
         write(screen,fmdlm1) delim3
       endif
       write(char_r,'(i0)') N_hole 
       write(screen,rdchr2) 'Total number of holes', adjustl(char_r)
       call read_1_real(u,hole_hop,'HOLE_HOP','Hole hopping attempt rate [s^-1]',280.0d12,"(es10.3)")
       call read_1_real(u,hole_rel,'HOLE_RELAX','Hole relaxation rate [/s]',280.0d12,"(es10.3)")
     else
       allocate(E_HOMO(0:N_HOMO))
     endif
     
!#####################################################################!
!#### Read Physical - electron variables #############################!
!#####################################################################! 

     if(inc_elec) then
       write(screen,fmdlm1) blank1
       write(screen,fmchr1) "Electron-specific Variables:"
       write(screen,fmdlm1) blank1
       N_LUMO=read_count_states(u,"LUMO")
       allocate(E_LUMO(0:N_LUMO))
       call read_1_real(u,E_LUMO(0),'LUMO_ENERGY','LUMO Energy [eV]',0.0d0,"(f10.5)")
       do i_state = 1, N_LUMO
         write(dum_str1,'("LUMO+",i0,"_ENERGY")') i_state
         write(dum_str2,'("LUMO+",i0," Energy [eV]")') i_state
         call read_1_real(u,E_LUMO(i_state),trim(adjustl(dum_str1)),trim(adjustl(dum_str2)),0.0d0,"(f10.5)")
       enddo
       E_LUMO = E_LUMO * eV_to_Eh
       call read_1_real(u,density,'ELEC_DENSITY','Electron density [/site]',0.001d0,"(f10.5)")
       if(density.le.0) call wrong_value_1_real(density,'ELEC_DENSITY',0.001d0,erchfl)
       N_elec = nint(density*dble(Nx*Ny*Nz))
       if(N_elec .LT. 1) then
         write(screen,fmdlm1) delim3
         write(screen,fmchr1) 'WARNING:'
         write(screen,fmchr1) 'Simu. cell too small to contain electrons.'
         write(screen,fmchr1) 'ACTION: increase the number of voxels or the density.'
         write(screen,fmchr1) 'Number of electrons fixed to 1.'
         N_elec = 1
         write(char_r,'(f10.5)') dble(N_elec) / dble(Nx*Ny*Nz) 
         write(screen,rdchr2) 'Electron density [/site]', adjustl(char_r)
         write(screen,fmdlm1) delim3
       endif
       write(char_r,'(i0)') N_elec 
       write(screen,rdchr2) 'Total number of electrons', adjustl(char_r)
       call read_1_real(u,elec_hop,'ELEC_HOP','Electron hopping attempt rate [s^-1]',280.0d12,"(es10.3)")
       call read_1_real(u,hole_rel,'ELEC_RELAX','Electron relaxation rate [/s]',280.0d12,"(es10.3)")
     else
       allocate(E_LUMO(0:N_LUMO))
     endif
     N_chrg = N_hole + N_elec
     
!#####################################################################!
!#### Read Guest-Host variables ######################################!
!#####################################################################!      
         
     if(guest_host_flag) then
      write(screen,fmdlm1) blank1
      write(screen,fmchr1) "A Guest Host System has been requested"
      write(screen,fmdlm1) blank1
      call read_1_real(u,guest_density,'GUEST_DENSITY','Guest Density',0.0d0,"(f10.5)")
      if(inc_hole) then
       GUEST_N_HOMO=read_count_states(u,"HGUEST")
       allocate(GUEST_E_HOMO(0:GUEST_N_HOMO))
       call read_1_real(u,GUEST_E_HOMO(0),'HGUEST_ENERGY','HOMO Energy in the Guest [eV]',0.0d0,"(f10.5)")
       do i_state = 1, GUEST_N_HOMO
         write(dum_str1,'("HGUEST+",i0,"_ENERGY")') i_state
         write(dum_str2,'("HOMO+",i0," Energy in the Guest [eV]")') i_state
         call read_1_real(u,GUEST_E_HOMO(i_state),trim(adjustl(dum_str1)),trim(adjustl(dum_str2)),0.0d0,"(f10.5)")
       enddo
       GUEST_E_HOMO = - GUEST_E_HOMO * eV_to_Eh
      else
       allocate(GUEST_E_HOMO(0:N_HOMO))      
      endif
      if(inc_elec) then
       GUEST_N_LUMO=read_count_states(u,"LGUEST")
       allocate(GUEST_E_LUMO(0:GUEST_N_LUMO))
       call read_1_real(u,GUEST_E_LUMO(0),'LGUEST_ENERGY','LUMO Energy in the Guest [eV]',0.0d0,"(f10.5)")
       do i_state = 1, GUEST_N_LUMO
         write(dum_str1,'("LGUEST+",i0,"_ENERGY")') i_state
         write(dum_str2,'("LUMO+",i0," Energy in the Guest [eV]")') i_state
         call read_1_real(u,GUEST_E_LUMO(i_state),trim(adjustl(dum_str1)),trim(adjustl(dum_str2)),0.0d0,"(f10.5)")
       enddo
       GUEST_E_LUMO = GUEST_E_LUMO * eV_to_Eh
      else
       allocate(GUEST_E_LUMO(0:N_LUMO))      
      endif
     endif 

!#####################################################################!
!#### Finished Reading ###############################################!
!#####################################################################! 
     write(screen,fmdlm1) blank1
     flush(screen)
     end subroutine read_input
!#####################################################################!
     subroutine sanity_checks()
     use my_mpi
     implicit none
     logical :: sane, nuts
     sane = .true.
     nuts = .false.
     write(screen,fmdlm1) blank1
     write(screen,fmchr1) "Performing sanity checks on the input parameters"
     write(screen,fmdlm2) delim2,blank1
     if(print_verb.ne.1.and.print_verb.ne.2) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Danger! You have selected an invalid verbosity"
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "This can be either 1 or 2"
       write(screen,fmchr1) "VERBOSITY has been set to 1"
       write(screen,fmdlm1) blank1
       print_verb = 1
       sane = .false.       
     endif
     if(NN_corr.gt.0.and.dipole.eq.0.0d0) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Danger! Your Dipole is zero, but the number of correlated neighbours is greater than zero"
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "We don't want to be messing around with dipole vectors in this case"
       write(screen,fmchr1) "CORRELATED_NEIGHBORS has been set to zero"
       write(screen,fmdlm1) blank1
       NN_corr = 0
       sane = .false.
     endif
     if(order_parameter.gt.0.d0.and.dimer_flag) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Danger! Your order paramter is non-zero, but you have requested dimers"
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "We don't allow both routines at one"
       write(screen,fmchr1) "ORDER_PARAMETER has been set to zero"
       write(screen,fmdlm1) blank1
       order_parameter = 0.d0
       sane = .false.
     endif   
     if(NN_corr.eq.0.and.dimer_flag) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Danger! Your number of correlated neighbours is zero, but you've asked for dimers"
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "The dimers are pointless is they can't talk to their neighbours"
       write(screen,fmchr1) "CORRELATED_NEIGHBORS has been set to 7"
       write(screen,fmdlm1) blank1
       NN_corr = 7
       sane = .false.
     endif
     if(NN_corr.eq.0.and.dipole.ne.0.0d0) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Danger! Your number of correlated neighbours is zero, but you've asked for dipoles"
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "The dipoles are pointless is they can't talk to their neighbours"
       write(screen,fmchr1) "CORRELATED_NEIGHBORS has been set to 7"
       write(screen,fmdlm1) blank1
       NN_corr = 7
       sane = .false.
     endif
     if(NN_corr.gt.min(Nx,Ny,Nz)) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Danger! You are asking for more correlated neighbours than there are grid points"
       write(screen,fmdlm1) delim3
       NN_corr = min(Nx,Ny,Nz) - 1
       write(screen,fmchin) "CORRELATED_NEIGHBORS has been changed to   : ", NN_corr
       write(screen,fmdlm1) blank1
       sane = .false.
     endif
     if(dimer_rho.gt.0.5d0) then     
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Danger! You are asking for too many dimers per site"
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "The maximum number of dimers per site is 1/2"
       dimer_rho = 0.5d0
       write(screen,fmchin) "DIMER_DENSITY has been changed to   : ", dimer_rho
       write(screen,fmdlm1) blank1
       sane = .false.
     endif
     if(N_hopbox.lt.1) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Danger! Your hop box is too small"
       write(screen,fmdlm1) delim3
       N_hopbox = 1
       write(screen,fmchin) "HOPBOX has been changed to   : ", N_hopbox
       write(screen,fmdlm1) blank1
       sane = .false.
     endif
     if(N_hopbox.gt.min(Nx,Ny,Nz)) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Danger! Your hop box is too large"
       write(screen,fmdlm1) delim3
       N_hopbox = min(Nx,Ny,Nz)
       write(screen,fmchin) "HOPBOX has been changed to   : ", N_hopbox
       write(screen,fmdlm1) blank1
       sane = .false.
     endif
     if(mod(Nsteps,modulo_print).ne.0) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Danger! Print step size is indivisible with the number of simulation steps"         
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "NUMBER_STEPS must be a multiple of (N_PRINT_STEP)"
       Nsteps = (int(Nsteps / modulo_print) + 1) * modulo_print
       write(screen,fmchin) "NUMBER_STEPS has been changed to   : ", Nsteps
       sane = .false.
     endif
     if(mod(N_equilibration,modulo_print).ne.0) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Danger! Erroneous Sample Size for most-recent mobility calculation"         
       write(screen,fmdlm1) delim3
       if(N_equilibration.lt.modulo_print * 2) then
         N_equilibration = 2 * modulo_print
         write(screen,fmchr1) "The following inequality MUST be satisfied:"
         write(screen,fmchr1) "(N_equilibration) > 2 x (N_PRINT_STEP)"
       else
         N_equilibration = int(N_equilibration / modulo_print) * modulo_print
         write(screen,fmchr1) "(N_equilibration) must be a multiple of (N_PRINT_STEP)"
       endif
       write(screen,fmchin) "(N_equilibration) has been changed to   : ", N_equilibration
       write(screen,fmdlm1) blank1
       sane=.false.     
     endif
     if(.not.mill_rate.and..not.mrcs_rate) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Danger! You have not declared a valid charge rate type"         
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "(CHARGE_RATE) has been set to Miller-Abrahams"
       mill_rate = .true.
       write(screen,fmdlm1) blank1
       sane=.false.            
     endif
     if(field * lattice .gt. lambda .and. mrcs_rate) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Fatal! The parameters you've chosen exceed the limits of Marcus Rate Theory"
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "The following inequality MUST be satisfied:"
       write(screen,fmchr1) "(ELECTRIC_FIELD) X (LATTICE) < (REORGANIZATION)"
       write(screen,fmdlm1) blank1
       write(screen,fmchfl) 'ELECTRIC_FIELD                : ', field / Vnm_to_AU
       write(screen,fmchfl) 'LATTICE                       : ', lattice / nm_to_a0
       write(screen,fmchfl) 'REORGANIZATION                : ', lambda / eV_to_Eh
       write(screen,fmdlm1) blank1
       nuts=.true.
     endif
     if(guest_host_flag) then
      if(dipole.ne.0.0d0) then
        write(screen,fmdlm1) delim3
        write(screen,fmchr1) "Fatal! You have requested dipoles in a guest-host system"
        write(screen,fmdlm1) delim3
        write(screen,fmchr1) "This is not currently supported. You'll have to chose one"
        write(screen,fmdlm1) blank1
        nuts=.true.
      endif
      if(nint(Nx * Ny * Nz * guest_density).lt.1) then
        write(screen,fmdlm1) delim3
        write(screen,fmchr1) "Fatal! You have requested an unachievably low guest density"
        write(screen,fmdlm1) delim3
        write(screen,fmchr1) "Just run the code without the guest-host parameters in this case"
        write(screen,fmdlm1) blank1
        nuts=.true.
      endif
     endif
     if(Nsteps.eq.0) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Fatal! Number of simulation steps set to zero!"
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "The indefinite run feature has been disabled for compatibility"
       write(screen,fmchr1) "If an indefinite run is needed, please set it to by arbitrarily large"
       write(screen,fmdlm1) blank1
       nuts=.true.
     endif
     if(Nx*Ny*Nz.ge.max_integer) then
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "Fatal! Grid Dimensions are too large!"
       write(screen,fmdlm1) delim3
       write(screen,fmchr1) "The grid you've choosen is so big it is causing integer overflow problems"
       write(screen,fmchr1) "Please make sure that Nx * Ny * Nz < 2^(31) - 1"
       write(screen,fmdlm1) blank1
       nuts=.true.
     endif
     if(nuts) then
       write(screen,fmchr1) "One or more Fatal Error states have been detected!"
       write(screen,fmchr1) "To proceed, please fix the input file"
       write(screen,fmdlm2) blank1, delim0
       stop
     endif
     if(sane) then
       write(screen,fmchr1) "Congratulations, you are sane!!"
     else
       write(screen,fmdlm1) blank1
       write(screen,fmchr1) "Congratulations, you have been made sane!!"
     endif
     write(screen,fmdlm1) blank1
     end subroutine sanity_checks

    subroutine converter(val,con,nm,u1,u2)
    implicit none
    double precision, intent(inout)  :: val
    double precision                 :: val_keep
    double precision, intent(in)     :: con
    character(*), intent(in) :: nm
    character(*), intent(in) :: u1,u2
    val_keep = val;   val = val * con
    write(screen,10) trim(nm), val_keep, trim(u1), val, trim(u2)
10  format(t2,"|",t5,a,t35,":",2x,f10.8,1x,"[",a,"]",t59,"===>",2x,f20.15,1x,"[",a,"]",t100,"|")    
    end subroutine converter

    subroutine converter_b(val,con,nm,u1,u2)
    implicit none
    double precision, intent(inout)  :: val
    double precision                 :: val_keep
    double precision, intent(in)     :: con
    character(*), intent(in) :: nm
    character(*), intent(in) :: u1,u2
    val_keep = val;   val = val * con
    write(screen,10) trim(nm), val_keep, trim(u1), val, trim(u2)
10  format(t2,"|",t5,a,t35,":",2x,e10.4,1x,"[",a,"]",t59,"===>",2x,f20.15,1x,"[",a,"]",t100,"|")    
    end subroutine converter_b


    subroutine convert_AU()
     implicit none
     write(screen,fmdlm1) blank1
     write(screen,fmchr1) 'Conversion to AU units:'
     write(screen,fmdlm2) delim2, blank1
       N_equilibration = N_equilibration / modulo_print
       call converter(lattice,nm_to_a0,"lattice","nm","a0")
       call converter(sigma,eV_to_Eh,"sigma","eV","Eh")
       call converter(dipole,D_to_ea0,"dipole","D","e.a0")
       call converter(dimer_dipole,D_to_ea0,"dimer dipole","D","e.a0")
       call converter(lambda,eV_to_Eh,"lambda","eV","Eh")
       call converter(E_fermi,eV_to_Eh,"E_Fermi","eV","Eh")
       if(inc_hole) call converter_b(hole_hop,Hz_to_AU,"hole_hop","Hz","Eh/hbar")
       if(inc_elec) call converter_b(elec_hop,Hz_to_AU,"elec_hop","Hz","Eh/hbar")
       if(inc_hole) call converter_b(hole_rel,Hz_to_AU,"hole_rel","Hz","Eh/hbar")
       if(inc_elec) call converter_b(elec_rel,Hz_to_AU,"elec_rel","Hz","Eh/hbar")
       call converter(field,Vnm_to_AU,"E_field","V/nm","Eh/(e.a0)")
       field_step   = field * lattice
       density_real = density
       call converter(density_real,lattice**(-3),"density","/site","/a0^3")
       ! elec. coupl. extinction coeff.:
       alpha = 1.0d0/lattice
       write(screen,10) 'Alpha parameter [nm^-1]       :', alpha * nm_to_a0, "/nm", alpha, "/a0"
       ! Coulomb interaction between 2 adjacent sites
       coulomb_factor = 1.0d0/(epsilon_r * lattice)
       write(screen,10) 'Coulomb factor', coulomb_factor * Eh_to_eV, "eV", coulomb_factor, "Eh"
       ! Marcus rate pre-factor:
!       hop_prefactor = (2.0d0 * pi * (coupling**2)) / &
!                          sqrt(4.0d0 * pi * lambda * kB_in_AU * temperature)
       ! Charge localization factor:       
       local         = local * 2.0d0 / nm_to_a0 
       write(screen,10) 'Charge Localization', local * nm_to_a0 / 2.0d0, "/nm", local / 2.0d0, "/a0"
       hop_prefactor = 1.0d0
       if(mrcs_rate) then
!         hop_prefactor = sqrt(pi / (lambda * kB_in_AU * temperature))
         write(screen,11) 'Marcus prefactor', hop_prefactor, "Unitless"
       else
!         hop_prefactor = sqrt(pi / (lambda * kB_in_AU * temperature))
         write(screen,11) 'Miller-Abrahams prefactor', hop_prefactor, "Unitless"
       endif
       if(inc_hole.and..not.inc_elec) then
         hop_prefactor = hop_prefactor * hole_hop
       else
         hop_prefactor = hop_prefactor * elec_hop
       endif
       coul_cutoff = coul_cutoff / lattice
       write(screen,fmdlm2) blank1, blank1

       write(screen,fmchr1) 'Converters used in the calculation:'
       write(screen,fmdlm2) delim2, blank1
       write(screen,10) "Energy", 1.0, "eV",    eV_to_Eh, "Eh"
       write(screen,10) "Length", 1.0, "nm",    nm_to_a0, "a0"
       write(screen,10) "Dipole", 1.0, "D",     D_to_ea0, "e.a0"
       write(screen,12) "Time",   1.0, "s",     1.d0 / hEh_to_s, "hbar/Eh"
       write(screen,12) "Field" , 1.0, "V/nm",   Vnm_to_AU, "Eh/(e.a0)"
       write(screen,fmdlm1) blank1
       write(screen,12) "Current" , 1.0, "AU",  AU_to_A,  "A/m^2"
       write(screen,12) "Mobilty" , 1.0, "AU",  mb_to_SI, "cm^2/(Vs)"
       write(screen,fmdlm2) blank1, blank1

       write(screen,fmchr1) 'Constants used in the calculation:'
       write(screen,fmdlm2) delim2, blank1
       write(screen,fmchfl) 'pi                            : ', pi
       write(screen,fmches) 'e  [C]                        : ', e_to_C
       write(screen,fmches) 'kB [Eh/K]                     : ', kB_in_AU
!       write(screen,fmchfl) 'epsilon_0 [e^2 / (a0.Eh)]     : ', e0_in_AU
       write(screen,fmdlm2) blank1, delim1




       flush(screen)
10  format(t2,"|",t5,a,t35,":",2x,f10.8,1x,"[",a,"]",t59,"===>",2x,f20.15,1x,"[",a,"]",t100,"|") 
11  format(t2,"|",t5,a,t35,":",2x,e10.4,1x,"[",a,"]",t100,"|") 
12  format(t2,"|",t5,a,t35,":",2x,f10.8,1x,"[",a,"]",t59,"===>",2x,e20.11,1x,"[",a,"]",t100,"|") 

    end subroutine convert_AU


       subroutine charge_init(chrg,E_fermi,kBT,state,mask)
       use general_functions
       use my_mpi, only: i_am_root
       implicit none
       type(charge)  :: chrg 
       double precision, intent(in)                          :: E_fermi, kBT
       character(4), intent(in)                              :: state
       logical, dimension(chrg%Nx,chrg%Ny,chrg%Nz), optional :: mask
       ! Local variables:
       type(indexes), dimension(chrg%Nx*chrg%Ny*chrg%Nz)     :: grid
       integer                                               :: placed_charges
       integer                                               :: ix,iy,iz,j,Ntot
       double precision                                      :: eta
       double precision, dimension(chrg%Nx*chrg%Ny*chrg%Nz)  :: p, s
       Ntot = chrg%Nx*chrg%Ny*chrg%Nz
       if(i_am_root) then
         write(screen,fmchin) 'Number of particles to place  :', chrg%N
         write(screen,fmchin) 'Total number of sites         :', Ntot
         write(screen,fmchfl) 'Fermi energy [eV]             :', E_fermi * Eh_to_eV
         write(screen,'(t2,"|",t5,a4," energy [eV]              :  ",f20.12,t100,"|")') state, chrg%E0(0) * Eh_to_eV
         write(screen,fmdlm1) blank1
         flush(screen)
       endif
       chrg%distri    = .false.
       placed_charges = 0
       j              = 0
       if (PRESENT(mask)) then 
        do iz=1, chrg%Nz
         do iy=1, chrg%Ny
          do ix=1, chrg%Nx
           if(mask(ix,iy,iz)) then
            j         = j + 1
            p(j)      = 0.d0
           else
            j         = j + 1
            p(j)      = 1.0d0 / ( 1.0d0 + exp( (chrg%E(ix,iy,iz,0) - E_fermi) / kBT ) )
            grid(j)%x = ix
            grid(j)%y = iy
            grid(j)%z = iz
           endif
          enddo
         enddo
        enddo
       else
        do iz=1, chrg%Nz
         do iy=1, chrg%Ny
          do ix=1, chrg%Nx
           j         = j + 1
           p(j)      = 1.0d0 / ( 1.0d0 + exp( (chrg%E(ix,iy,iz,0) - E_fermi) / kBT ) )
           grid(j)%x = ix
           grid(j)%y = iy
           grid(j)%z = iz
          enddo
         enddo
        enddo
       endif
       s = cumulative_sum(p)
       do while(placed_charges.lt.chrg%N)
        call random_number(eta)
        j                      = findeta(s,eta)
        p(j)                   = 0.0d0
        placed_charges         = placed_charges + 1
        chrg%R(placed_charges) = grid(j)
        s                      = cumulative_sum(p)
       enddo
       do j=1,chrg%N
        chrg%distri(chrg%R(j)%x,chrg%R(j)%y,chrg%R(j)%z) = .true.       
       enddo
       endsubroutine charge_init

end module initialization

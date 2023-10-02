!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Param.f:
!%
!% Contains the global variables
!% 
!% Copyright (C) 2020 Thomas Pope
!% \author Yvelin Giret & Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%% constants        %%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module constants
    use kinds
    implicit none
!------ general parameters
    integer,          parameter, public :: screen   = 6                   ! Output Unit
    double precision, parameter, public :: pi       = 3.1415926535897932384626433d0 ! 
    double precision, parameter, public :: petit    = 1.0d-8              ! small number
!------ atomic unit parameters 
!------- constants
    double precision, parameter, public :: kB_in_AU =  3.16681156349d-6   ! [Eh.K^-1]
!------- converters into AU
    double precision, parameter, public :: eV_to_Eh  =  0.036749322176d0     ! [Eh/eV]
    double precision, parameter, public :: nm_to_a0  = 18.897261246365d0     ! [a0/nm]
    double precision, parameter, public :: D_to_ea0  =  0.393430307364d0     ! [e.a0/D] 1 Debye in e.a0 
    double precision, parameter, public :: Vnm_to_AU =  1.94469038115d-3     ! [Eh(e.a0)^-1 / Vm^-1]
    double precision, parameter, public :: Hz_to_AU  =  2.4188843265857d-17  ! [Eh.hbar^-1 / Hz]
!------- converters into SI
    double precision, parameter, public :: Eh_to_eV = 27.211386245988d0     ! [eV/Eh] Hartree energy in eV
    double precision, parameter, public :: a0_to_nm =  0.052917721090d0     ! [nm/a0] Borh radius in nm
    double precision, parameter, public :: ea0_to_D =  2.541746228703d0     ! [D/e.a0] Dipole
    double precision, parameter, public :: hEh_to_s =  2.4188843265857d-17  ! [s / hbar(Eh)^-1] Time
    double precision, parameter, public :: AU_to_Hz =  4.13413733351d16     ! [Hz / Eh.hbar^-1] Inverse Time
    double precision, parameter, public :: e_to_C   =  1.602176634d-19      ! [C / e] Charge
    double precision, parameter, public :: AU_to_A  =  2.3653370d18         ! [(A/m^2) / (e.Eh.hbar^-1/a0^2)] Current
    double precision, parameter, public :: mb_to_SI =  4.25438215733d-2     ! (ea0^2/hbar) to cm2/Vs -- Mobility 
endmodule constants
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%% input parameters %%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module param
    use kinds
    implicit none  
    integer                                     :: MPI_N_calc                     ! Number of calculations to be run on the same system
    integer                                     :: Nx, Ny, Nz                     ! Number of grid points in xyz directions
    integer                                     :: N_hopbox                       ! Number of nieghbours considered hop-viable
    integer                                     :: Nsteps, iteration              ! Number of steps in the MC loop and index parameter
    integer                                     :: NN_corr                        ! Number of neighbours for correlated disorder
    integer                                     :: modulo_print, modulo_print_xyz ! Number of steps between printing outputs (data and 3d grids)
    integer                                     :: N_equilibration                ! Sample size used for results equilibration
    integer                                     :: N_SING, N_TRIP                 ! Number of singlet/triplets per site
    integer                                     :: print_verb                     ! Verbosity
    double precision                            :: lattice                        ! Lattice parameter
    double precision                            :: coul_cutoff                    ! Coulomb cut-off term
    double precision                            :: local                          ! Charge Localization parameter
    double precision                            :: sigma, dipole                  ! Noise width of uncorrelated/correlated disorder
    double precision                            :: lambda                         ! Reorganization energy for Marcus rate equation
    double precision                            :: E_fermi                        ! Fermi Energy
    integer                                     :: N_HOMO,  N_LUMO                ! Number of HOMO and LUMO levels considered
    double precision, allocatable, dimension(:) :: E_HOMO,  E_LUMO                ! HOMO and LUMO energies
    integer                                     :: GUEST_N_HOMO,  GUEST_N_LUMO    ! Number of HOMO and LUMO levels considered in the guest
    double precision, allocatable, dimension(:) :: GUEST_E_HOMO,  GUEST_E_LUMO    ! HOMO and LUMO energies for the guest
    double precision                            :: guest_density                  ! Density of the Guest Molecules
    double precision                            :: hole_hop, elec_hop             ! Hopping attempt rates for electrons and holes
    double precision                            :: hole_rel, elec_rel             ! Relaxation rates for electrons and holes
    double precision                            :: temperature                    ! Temperature
    double precision                            :: field                          ! Electric field strength
    double precision                            :: field_vector(3)                ! field vector - used only for the direction
    double precision                            :: density                        ! Hole density per site
    double precision                            :: epsilon_r                      ! Relative permittivity
    double precision                            :: order_parameter                ! Weight of anti-ferromagnetic ordering 
    double precision                            :: dimer_dipole                   ! Value of the dipole for the dimer
    double precision                            :: dimer_rho                      ! Density of dimers in the system
    logical                                     :: restart                        ! Flag to check if restarting an old calculation
    logical                                     :: debugging                      ! Inhibit RNG for debugging purposes
    logical                                     :: energy_init                    ! Flag to read or to calculate site energies
    logical                                     :: distri_init                    ! Flag to read or to calculate charge distribution
    logical                                     :: dimer_flag                     ! Flag to generate dimers for the dipole routine
    character(ch)                               :: run_type                       ! Type of calculation (hole, hole + electron, etc)
    character(ch)                               :: spatial                        ! Type of spatial distribution of grid points (hole, hole + electron, etc)
    character(ch)                               :: rate_type                      ! Type of rate calculation
    character(ch)                               :: pbc_type                       ! Type of periodic boundary conditions
    character(ch)                               :: coulomb_type                   ! Type of coulomb calculation
    character(ch)                               :: correlated_type                ! Type of correlation
    double precision                            :: time_save_test(4)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%% input restarters %%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    double precision                            :: time_tot                       ! Total time steps calculated so far
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%% input dependants %%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    integer                                     :: N_hole                         ! Number of holes in the system
    integer                                     :: N_elec                         ! Number of electrons in the system
    integer                                     :: N_exct                         ! Maximum available exicton in the system
    integer                                     :: N_chrg                         ! Total number of charge particles in the system
    integer                                     :: N_particle                     ! Total number of charge particles and possible excitons in the system
    integer                                     :: max_events                     ! Number of possible events per step
    double precision                            :: alpha                          ! Inverse of the Lattice constant
    double precision                            :: coulomb_factor                 ! Coulomb factor for calculating coulomb energy
    double precision                            :: hop_prefactor                  ! Marcus Factor for calculating jump rates
    double precision                            :: density_real                   ! Density in SI units
    double precision                            :: field_step                     ! Field stepsize per lattice point
    logical                                     :: inc_elec                       ! Flag signifying the inclusion of electrons
    logical                                     :: inc_hole                       ! Flag signifying the inclusion of holes
    logical                                     :: inc_excn                       ! Flag signifying the inclusion of excitations
    logical                                     :: mrcs_rate                      ! Flag signifying the use of marcus rates
    logical                                     :: mill_rate                      ! Flag signifying the use of miller-abrahams rates
    logical                                     :: guest_host_flag                ! Flag signifying the use of the guest-host algorithm
    logical                                     :: test_the_new_thing
end module param
              

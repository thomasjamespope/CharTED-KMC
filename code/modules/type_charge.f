!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% type_charge.f
!%
!% module type_charge: 
!%  Definition of the charge derived data type
!%  and its associated procedures
!%
!% Copyright (C) 2020 Yvelin Giret & Thomas Pope
!% \author Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module type_charge
    use kinds
    use type_coor
    use type_indexes
!********************************************************************************************!
! charge data type:                                                                          !
!   contains lots of things, all described below                                             !
!********************************************************************************************!
    type charge
    ! general system information that is good to keep available
      integer                                           :: Nx,Ny,Nz    ! Grid dimensions
      integer                                           :: N           ! Number of charged particles
      integer                                           :: NE          ! Number of States per site
      double precision, allocatable, dimension(:)       :: E0          ! State energy without disorder
      double precision                                  :: omega       ! Hopping attempt rate
      double precision, allocatable, dimension(:,:,:,:) :: E           ! Site energies
    ! Information on the location of the charged particle 
      logical,          allocatable, dimension(:)       :: exists      ! Existence status of each charged particle
      type(indexes),    allocatable, dimension(:)       :: R           ! Location of each charged particle
      logical,          allocatable, dimension(:,:,:)   :: distri      ! Particle distribution
      integer,          allocatable, dimension(:)       :: iE          ! State level index of each charged particle
    ! Information on the rate calculation
      logical,          allocatable, dimension(:)       :: log_rate    ! Flag indicating if the rates needs to be calculated
      double precision, allocatable, dimension(:,:,:)   :: rate        ! rate component memory 
    ! Exciton information, only used by singlets and triplets
      integer,          allocatable, dimension(:)       :: iH, iL      ! Component level indices for each excited state
      double precision                                  :: relax_rate  ! Relaxation rate
      integer                                           :: relax_count
      integer,          allocatable, dimension(:,:)     :: Eindex      ! excitation index array
    ! Guest-Host information, only used for guest-host systems
      integer,          allocatable, dimension(:,:,:)   :: guest       ! Guest distribution
      integer                                           :: Nguest      ! Number of guest molecules
      integer                                           :: guest_NE
      double precision, allocatable, dimension(:)       :: guest_E0    ! Site energies without disorder for the guest
      double precision, allocatable, dimension(:,:)     :: guest_E     ! Site energies for the guest
    ! Hopping information
      type(indexes)                                     :: o,d         ! Storage for hopping event origin and destination
      double precision                                  :: flux        ! charge flux accumulator
    ! Information used in generating output
      integer,          allocatable, dimension(:,:,:)   :: site_flux   ! Total site flux
      integer                                           :: Nd          ! Exciton decay counter (if applicable)
    contains
      procedure :: rate_init        => charge_rate_init
      procedure :: aloc             => charge_allocate
      procedure :: guest_aloc       => guest_allocate
      procedure :: update_flux      => site_flux_update
      procedure :: update           => update_all
      procedure :: check_mat        => check_charge_matrix
      procedure :: hopid            => hops_indexdirection
      procedure :: Eindex_init
      procedure :: charge_indexes_find
      procedure :: charge_index_check
      generic   :: operator(.find.) => charge_indexes_find
      generic   :: operator(.chk.)  => charge_index_check
    endtype charge
    contains
!--------------------------------------------------------------------------------------------!
! charge_rate_init (remampped to the "%rate_init" procedure for charge data types):          !
!     allocates the rate memory arrays                                                       !
! usage, for charge type varbiable a                                                         !
!     call a%rate_init                                                                       !
!--------------------------------------------------------------------------------------------!
    subroutine charge_rate_init(chrg,N)
    implicit none
    class(charge)       :: chrg
    integer, intent(in) :: N
    if(chrg%guest_NE.eq.-1) then
     allocate(chrg%rate(0:chrg%NE,N,chrg%N))
    else
     allocate(chrg%rate(0:max(chrg%NE,chrg%guest_NE),N,chrg%N))
    endif
    allocate(chrg%log_rate(chrg%N))
    chrg%log_rate = .true.
    endsubroutine charge_rate_init
!--------------------------------------------------------------------------------------------!
! charge_allocate (remampped to the "%aloc" procedure for charge data types):                !
!     allocates the distribution arrays in the charge data type                              !
! usage, for charge type varbiable a                                                         !
!     call a%aloc                                                                            !
!--------------------------------------------------------------------------------------------!
    subroutine charge_allocate(chrg,Nx,Ny,Nz,NE,E0,hop_rate,rel_rate,N)
    implicit none
    class(charge)        :: chrg
    integer, intent(in)  :: Nx, Ny, Nz, NE, N
    double precision, intent(in) :: E0(0:NE), hop_rate, rel_rate
    integer              :: i
    chrg%relax_count = 0
    chrg%Nx = Nx; chrg%Ny = Ny; chrg%Nz = Nz; chrg%N = N
    allocate(chrg%distri(chrg%Nx,chrg%Ny,chrg%Nz))
    allocate(chrg%site_flux(chrg%Nx,chrg%Ny,chrg%Nz))
    chrg%site_flux = 0
    chrg%N = N 
    if(chrg%N.gt.0) then
     allocate(chrg%R(chrg%N),chrg%iE(chrg%N),chrg%exists(chrg%N))
     chrg%exists = .true.
    endif
    chrg%NE = NE
    chrg%iE = 0
    allocate(chrg%E0(0:chrg%NE))
    allocate(chrg%E(chrg%Nx,chrg%Ny,chrg%Nz,0:chrg%NE))
    chrg%E0 = E0
    do i=0,chrg%NE
      chrg%E(:,:,:,i) = chrg%E0(i)
    enddo
    chrg%omega      = hop_rate
    chrg%relax_rate = rel_rate
    chrg%flux       = 0.d0
    chrg%guest_NE         = -1
    endsubroutine charge_allocate
!--------------------------------------------------------------------------------------------!
! guest_allocate (remampped to the "%guest_aloc" procedure for charge data types):           !
!     allocates the guest arrays in the charge data type and distributes them on the grid    !
! usage, for charge type varbiable a                                                         !
!     call a%guest_aloc                                                                      !
!--------------------------------------------------------------------------------------------!
    subroutine guest_allocate(chrg,guest_density,NE,E0,sigma,in_grid)
    use general_functions
    implicit none
    class(charge)                                                     :: chrg
    integer,                                     intent(in)           :: NE
    double precision,                            intent(in)           :: guest_density, E0(0:NE), sigma
    integer, dimension(chrg%Nx,chrg%Ny,chrg%Nz), intent(in), optional :: in_grid
    logical, dimension(chrg%Nx,chrg%Ny,chrg%Nz)                       :: guest_grid
    integer                                                           :: ix, iy, iz, ie
    allocate(chrg%guest(chrg%Nx,chrg%Ny,chrg%Nz))
    chrg%Nguest = 0
    chrg%guest  = 0
    if(present(in_grid)) then
     chrg%guest  = in_grid
     chrg%Nguest = maxval(chrg%guest)
    else
     guest_grid = random_sites(chrg%Nx,chrg%Ny,chrg%Nz,guest_density)
     do iz = 1, chrg%Nz; do iy = 1, chrg%Ny; do ix = 1, chrg%Nx
      if(guest_grid(ix,iy,iz)) then
       chrg%Nguest          = chrg%Nguest + 1
       chrg%guest(ix,iy,iz) = chrg%Nguest
      endif
     enddo; enddo; enddo
    endif
    chrg%guest_NE = NE
    allocate(chrg%guest_E0(0:chrg%guest_NE))
    chrg%guest_E0 = E0
    allocate(chrg%guest_E(0:chrg%guest_NE,chrg%Nguest))
    do ie=0,chrg%guest_NE
     chrg%guest_E(ie,:) = E0(ie) + random_noise(chrg%Nguest,sigma)
    enddo
    endsubroutine guest_allocate
!--------------------------------------------------------------------------------------------!
! Eindex_init:                                                                               !
!     allocates and sets up the index array for excitations                                  !
! usage, for charge type varbiable a                                                         !
!     call a%Eindex_init                                                                     !
!--------------------------------------------------------------------------------------------!
    subroutine Eindex_init(chrg,Nx,Ny,Nz,NH,EH,NL,EL,maxN) 
    implicit none
    class(charge)                                :: chrg
    integer                                      :: Nx, Ny, Nz, NH, NL, maxN
    double precision, dimension(0:NH)            :: EH
    double precision, dimension(0:NL)            :: EL
    double precision, dimension(0:(NL+1)*(NH+1)) :: Ediff
    integer                                      :: i,j
    double precision                             :: d1
    integer                                      :: d2,d3
    allocate(chrg%Eindex(0:NH,0:NL))
    allocate(chrg%iH(0:(NL+1)*(NH+1)))
    allocate(chrg%iL(0:(NL+1)*(NH+1)))
    chrg%NE = -1
    do i=0,NH
     do j=0,NL
      chrg%NE          = chrg%NE + 1
      Ediff(chrg%NE)   = EH(i) + EL(j)
      chrg%iH(chrg%NE) = i
      chrg%iL(chrg%NE) = j
     enddo
    enddo
    do i=1,chrg%NE
     d1 = Ediff(i)
     d2 = chrg%iH(i)
     d3 = chrg%iL(i)
     do j=i-1,0,-1
      if(Ediff(j).lt.d1) goto 10
      Ediff(j+1)   = Ediff(j)
      chrg%iH(j+1) = chrg%iH(j)
      chrg%iL(j+1) = chrg%iL(j)
     enddo
     j = -1
10   Ediff(j+1)   = d1
     chrg%iH(j+1) = d2
     chrg%iL(j+1) = d3
    enddo
    do i=0,chrg%NE
     chrg%Eindex(chrg%iH(i),chrg%iL(i)) = i
    enddo
    chrg%Nx = Nx; chrg%Ny = Ny; chrg%Nz = Nz
    allocate(chrg%distri(chrg%Nx,chrg%Ny,chrg%Nz))
    allocate(chrg%site_flux(chrg%Nx,chrg%Ny,chrg%Nz))    
    allocate(chrg%R(maxN),chrg%iE(maxN),chrg%exists(maxN))
    allocate(chrg%E0(0:chrg%NE))
    chrg%exists = .false.
    chrg%iE     = 0
    chrg%E0     = Ediff(0:chrg%NE)
    chrg%N      = 0
    endsubroutine Eindex_init
!--------------------------------------------------------------------------------------------!
! check_charge_matrix (remampped to the "%check_mat" procedure for charge data types):       !
!     chacks that the number of charge particles is consistent across all arrays             !
! usage, for charge type varbiable a                                                         !
!     call a%check_mat                                                                       !
!--------------------------------------------------------------------------------------------!
    subroutine check_charge_matrix(chrg,iteration,me_ip,text,i)
    implicit none
    class(charge)       :: chrg
    integer, intent(in) :: iteration, me_ip, i
    character(4)        :: text
    ! Local variables:
    integer :: ix, iy, iz, check
     check = 0
     do iz=1, chrg%Nz
      do iy=1, chrg%Ny
       do ix=1, chrg%Nx
        if(chrg%distri(ix,iy,iz)) check = check + 1
       enddo
      enddo
     enddo
     if(check .ne. chrg%N) then
        write(*,'(2x,a)') '*****************************************************'
        write(*,'(2x,a,2x,a)') 'WARNING: check_matrix.f90!', text
        write(*,'(2x,a,2x,i8)') 'Number of chrg from matrix: ', check
        write(*,'(2x,a,2x,i8)') 'Number of chrg from inputs: ', chrg%N
        write(*,'(2x,a,2x,i8)') 'At iteration              : ', iteration
        write(*,'(2x,a,2x,i8)') 'IP                        : ', me_ip
        write(*,'(2x,a,2x,i8)') 'i                         : ', i
        do ix=1,chrg%N; write(*,*)chrg%R(ix)%x, chrg%R(ix)%y, chrg%R(ix)%z; enddo
        write(*,'(2x,a)') '*****************************************************'
        flush(6)
     end if
    end subroutine check_charge_matrix
!--------------------------------------------------------------------------------------------!
! update_all (remampped to the "%update" procedure for charge data types):                   !
!     updates a bunch of stuff in the object ot account for a hopping event                  !
! usage, for charge type varbiable a                                                         !
!     call a%update                                                                          !
!--------------------------------------------------------------------------------------------!
    subroutine update_all(chrg,me,dir)
    implicit none
    class(charge)                      :: chrg
    integer,                intent(in) :: me
    double precision,       intent(in) :: dir
    chrg%flux = chrg%flux + dir 
    if(chrg%distri(chrg%d%x,chrg%d%y,chrg%d%z)) chrg%relax_count = chrg%relax_count + 1
    chrg%R(me) = chrg%d
    call chrg%update_flux(chrg%o)
    chrg%distri(chrg%o%x,chrg%o%y,chrg%o%z) = .false.
    chrg%distri(chrg%d%x,chrg%d%y,chrg%d%z) = .true.
    end subroutine update_all

    subroutine site_flux_update(chrg,R)
    implicit none
    class(charge) :: chrg
    type(indexes) :: R
    chrg%site_flux(R%x,R%y,R%z) = chrg%site_flux(R%x,R%y,R%z) + 1
    end subroutine site_flux_update
!--------------------------------------------------------------------------------------------!
! indexes_find (remampped to the ".find." operator for indexes data types):                  !
!     finds the charge particle with the same index                                          !
! usage, for indexes type varbiables a, a charge type variable b and a integer variable c    !
!     c = a.find.b                                                                           !
!--------------------------------------------------------------------------------------------!
    function charge_indexes_find(chrg,a) result(i)
    implicit none
    class(charge), intent(in) :: chrg
    type(indexes), intent(in) :: a
    integer                   :: i
    integer                   :: j
    i = 0
    do j=1,chrg%N
      if(a%x.eq.chrg%R(j)%x .and. a%y.eq.chrg%R(j)%y .and. a%z.eq.chrg%R(j)%z) i=j
    enddo
    endfunction charge_indexes_find
!--------------------------------------------------------------------------------------------!
! hops_indexdirection (remampped to the "%hopid" procedure for charge data types):           !
!     determines the direction of the hopping event for use in current and mobility analysis !
! usage, for charge type varbiable a and indexes data type j                                 !
!   j = a%hopid()                                                                            !
!--------------------------------------------------------------------------------------------!   
    function hops_indexdirection(chrg) result(j) 
    implicit none
    class(charge), intent(in) :: chrg
    type(indexes)             :: j
    j = chrg%d - chrg%o
    if(j%x.gt.1) then; j%x = -1; elseif(j%x.lt.-1) then; j%x = 1; endif
    if(j%y.gt.1) then; j%y = -1; elseif(j%y.lt.-1) then; j%y = 1; endif
    if(j%z.gt.1) then; j%z = -1; elseif(j%z.lt.-1) then; j%z = 1; endif
    endfunction hops_indexdirection
!--------------------------------------------------------------------------------------------!
! charge_index_check (remapped to the ".chk." operator for indexes and charge data types):   !
!     gives the distribution value for a charge data type at an indexes data type grid point !
! usage, for indexes type varbiable a, charge type variable b, and logical variable c        !
!     c = a.chk.b                                                                            !
!--------------------------------------------------------------------------------------------!
    function charge_index_check(chrg,i) result(l)
    implicit none
    class(charge), intent(in) :: chrg
    type(indexes), intent(in) :: i
    logical                   :: l
    l = chrg%distri(i%x,i%y,i%z)
    endfunction charge_index_check
!--------------------------------------------------------------------------------------------!
endmodule type_charge

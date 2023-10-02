!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% module coulomb.f90: ...description...
!%
!% Copyright (C) 2020 Yvelin Giret & Thomas Pope
!% \author Yvelin Giret & Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module coulomb
    use kinds
    use constants
    use pbc
    use param
    use type_indexes
    use my_mpi
    implicit none
    integer                                 :: Ncx,Ncy,Ncz !indexes
    double precision, allocatable, dimension(:,:,:) :: U           !coulomb matrix
    contains

    subroutine init_coloumb_matrix()
    implicit none
    integer :: i,j,k,di,dj,dk

    if(me_ip.eq.0) then
     write(screen,fmdlm1) blank1
     write(screen,fmchr1) "Initializing the Coulomb Matrix:"
     write(screen,fmdlm2) delim2, blank1
     flush(screen)
    endif
    
    Ncx=int(Nx/2); if(Ncx*2.ne.Nx) Ncx=Ncx+1
    Ncy=int(Ny/2); if(Ncy*2.ne.Ny) Ncy=Ncy+1
    Ncz=int(Nz/2); if(Ncz*2.ne.Nz) Ncz=Ncz+1
    
    allocate(U(Nx,Ny,Nz)); U=0.0d0 !for reasons, this array has to be a bit bigger than needed. Don't "fix" it!
    
    do i=1,Nz;     di = i - 1
      do j=1,Ny;   dj = j - 1
        do k=1,Nx; dk = k - 1
          U(i,j,k) = coulomb_factor / sqrt(dble(di**2 + dj**2 + dk**2))
          ! * include_neighbours(Nx,Ny,Nz,di,dj,dk) 
        enddo
      enddo
    enddo
    U(1,1,1) = coulomb_factor
    if(i_am_root) then
     open(100,file="OUTPUT_coulomb.dat")
     do i=1,Ncx; do j=1,Ncy; do k=1,Ncz
       write(100,*) i-1,j-1,k-1,U(i,j,k)/coulomb_factor
     enddo; enddo; enddo
     close(100)
     write(screen,10) Ncx, Ncy, Ncz
     write(screen,fmdlm2) blank1, delim1
     flush(screen)
    endif    
10  format(t2,"|",t5,2(i0,"X"),i0," grid point distance matrix has been constructed",t100,"|")
    endsubroutine init_coloumb_matrix
    
    function include_neighbours(Nx,Ny,Nz,ix,iy,iz) result(d)
    integer  :: Nx,Ny,Nz,ix,iy,iz
    integer  :: ux,uy,uz
    double precision :: d
    d = 0.0d0
    do ux=-1,1; do uy=-1,1; do uz=-1,1
      if(ux.eq.0.and.uy.eq.0.and.uz.eq.0.and.ix.eq.0.and.iy.eq.0.and.iz.eq.0) then
        d = d + 1
      else
        d = d + 1 / sqrt(dble((ux*Nx+ix)**2 + (uy*Ny+iy)**2 + (uz*Nz+iz)**2))
      endif
    enddo; enddo; enddo
    endfunction include_neighbours

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%% routine for calculating differnce in coulomb energy before and after a hopping %%%%%%%
!%%%%%% event in a single charge type system                                           %%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function coulomb_difference_single_sparse(N,o,d,R,me) result(Ec_diff)
       implicit none
       ! Input Variables:
       integer,       intent(in)               :: N        ! number of charge paticles 
       integer,       intent(in)               :: me       ! index of the active charge particle
       type(indexes), intent(in)               :: o,d      ! origin/destination point
       type(indexes), intent(in), dimension(N) :: R        ! location of all charges
       ! Output Variables:
       double precision                        :: Ec_diff  ! difference in coulomb energy
       ! Local variables:
       type(indexes)                            :: v        ! direction vector
       double precision                        :: oEc, dEc ! coloumb energy at origin and destination
       integer                                 :: i        ! generic index
       oEc = 0.0d0
       dEc = 0.0d0
       
       do i=1,N
         if(i.ne.me) then
           v   = (R(i).vec.o)
           oEc = oEc + U(v%x,v%y,v%z)           ! add coulomb energy for origin
           v   = (R(i).vec.d)
           dEc = dEc + U(v%x,v%y,v%z)            ! add coulomb energy for destination
!           oEc = oEc + ((R(i).vec.o).of.U)            ! add coulomb energy for origin
!           dEc = dEc + ((R(i).vec.d).of.U)            ! add coulomb energy for destination
         endif
       enddo
       Ec_diff = dEc - oEc
    endfunction coulomb_difference_single_sparse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%% routine for calculating differnce in coulomb energy before and after a hopping %%%%%%%
!%%%%%% event in a two charge type system                                              %%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function coulomb_difference_double_sparse(N1,N2,o,d,R1,R2,me) result(Ec_diff)
       implicit none
       ! Input Variables:
       integer,       intent(in)                :: N1, N2   ! number of charge paticles 
       integer,       intent(in)                :: me       ! index of the active charge particle
       type(indexes), intent(in)                :: o,d      ! origin/destination point
       type(indexes), intent(in), dimension(N1) :: R1       ! location of alike charges
       type(indexes), intent(in), dimension(N2) :: R2       ! location of other charges
       ! Output Variables:
       double precision                         :: Ec_diff  ! difference in coulomb energy
       ! Local variables:
       type(indexes)                            :: v        ! direction vector
       double precision                         :: oEc, dEc ! coloumb energy at origin and destination
       integer                                  :: i        ! generic index
       oEc = 0.0d0
       dEc = 0.0d0
       
       do i=1,N1                                       ! Loop over alike particles
         if(i.ne.me) then
           v   = (R1(i).vec.o)
           oEc = oEc + U(v%x,v%y,v%z)           ! add coulomb energy for origin
           v   = (R1(i).vec.d)
           dEc = dEc + U(v%x,v%y,v%z)            ! add coulomb energy for destination
!           oEc = oEc - ((R1(i).vec.o).of.U)              ! add coulomb energy for origin
!           dEc = dEc - ((R1(i).vec.d).of.U)              ! add coulomb energy for destination
         endif
       enddo
       do i=1,N2                                       ! Loop over other particles
         v   = (R2(i).vec.o)
         oEc = oEc - U(v%x,v%y,v%z)           ! add coulomb energy for origin
         v   = (R2(i).vec.d)
         dEc = dEc - U(v%x,v%y,v%z)            ! add coulomb energy for destination
!         oEc = oEc - ((R2(i).vec.o).of.U)              ! add coulomb energy for origin
!         dEc = dEc - ((R2(i).vec.d).of.U)              ! add coulomb energy for destination
       enddo
       Ec_diff = dEc - oEc
    endfunction coulomb_difference_double_sparse

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%% routine for calculating difference in coulomb energy in a one charge type system %%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function coulomb_energy_single(N,x,R,me) result(Ec)
       implicit none
       ! Input Variables:
       integer,       intent(in)               :: N  ! number of charge paticles 
       integer,       intent(in)               :: me ! index of the active charge particle
       type(indexes), intent(in)               :: x  ! charge location
       type(indexes), intent(in), dimension(N) :: R  ! location of alike charges
       ! Output Variables:
       double precision                        :: Ec ! coloumb energy at charge location
       ! Local variables:
       type(indexes)                           :: v  ! direction vector
       integer                                 :: i  ! generic index
       Ec = 0.0d0       
       do i=1,N                                       ! Loop over alike particles
        if(i.ne.me) then
         v  = (R(i).vec.x)
         Ec = Ec + U(v%x,v%y,v%z)  ! add coulomb energy for origin
        endif
       enddo
    endfunction coulomb_energy_single

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%% routine for calculating difference in coulomb energy in a two charge type system %%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function coulomb_energy_double(N1,N2,x,R1,R2,me) result(Ec)
       implicit none
       ! Input Variables:
       integer,       intent(in)                :: N1, N2   ! number of charge paticles 
       integer,       intent(in)                :: me       ! index of the active charge particle
       type(indexes), intent(in)                :: x        ! charge location
       type(indexes), intent(in), dimension(N1) :: R1       ! location of alike charges
       type(indexes), intent(in), dimension(N2) :: R2       ! location of other charges
       ! Output Variables:
       double precision                         :: Ec       ! coloumb energy at charge location
       ! Local variables:
       type(indexes)                            :: v        ! direction vector
       integer                                  :: i        ! generic index
       Ec = 0.0d0       
       do i=1,N1                                       ! Loop over alike particles
        if(i.ne.me) then
         v  = (R1(i).vec.x)
         Ec = Ec + U(v%x,v%y,v%z)  ! add coulomb energy for origin
        endif
       enddo
       do i=1,N2                   ! Loop over other particles
        v  = (R2(i).vec.x)
        Ec = Ec - U(v%x,v%y,v%z) ! add coulomb energy for origin
       enddo
    endfunction coulomb_energy_double


end module coulomb




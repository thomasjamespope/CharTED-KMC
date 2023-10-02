!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% module random.f90: ...description...
!%
!% Copyright (C) 2020 Yvelin Giret & Thomas Pope
!% \author Yvelin Giret & Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module random
    use kinds
    use constants
    use printings
    implicit none
    contains
!#####################################################################!
!#### Check the distributions and print the results ##################!
!#####################################################################!
    subroutine check_distribution(Nx,Ny,Nz,E)
    use my_mpi
    implicit none
    !input
    integer,                               intent(in) :: Nx, Ny, Nz
    double precision, dimension(Nx,Ny,Nz), intent(in) :: E
    !routine parameters
    double precision                   :: sigma, ave_E
    integer                            :: ix, iy, iz, i, ip
    integer(long_integer)              :: Ntot, long_Nx, long_Ny, long_Nz, long_Np
    integer,      dimension(-15:15)    :: n, np
    double precision, dimension(-3:3)  :: dist
    double precision, dimension(0:15)  :: box
    double precision, dimension(we_np) :: sigma_all, ave_E_all
    double precision                   :: ave_sigma, ave_ave_E
    double precision                   :: sig_sigma, sig_ave_E
    logical                            :: looking
    character(1), dimension(-15:15,9)  :: graph
    double precision, dimension(9)     :: ybox
    character(1), parameter            :: grph_top="x"
    character(1), parameter            :: grph_col="|"
    logical                            :: grph_go
    graph = " "
    n     = 0
    ave_E = sum(E) / (Nx*Ny*Nz)
    sigma = sqrt(sum((E-ave_E)**2) / (Nx*Ny*Nz))
    forall(i=0:15) box(i) = 0.2 * sigma * i
    do ix=1,Nz
      do iy=1,Ny
        do iz=1,Nx
          looking = .true.; i=1
          if(E(ix,iy,iz) - ave_E.lt.box(0)) then
            do while (i.le.15.and.looking)
              if(E(ix,iy,iz) - ave_E.ge.-box(i)) then
                looking = .false.
                n(-i) = n(-i) + 1
              endif
              i = i + 1
            enddo
          else
            do while (i.le.15.and.looking)
              if(E(ix,iy,iz) - ave_E.le.box(i)) then
                looking = .false.
                n(i) = n(i) + 1
              endif
              i = i + 1
            enddo
          endif
          if(looking) n(0) = n(0) + 1
        enddo
      enddo
    enddo
    if(i_am_root) then
     do ip=1,we_np-1
      call mpi_recv(ave_E_all(ip),1,we_real,ip,0,we_comm,istat,ierr)
      call mpi_recv(sigma_all(ip),1,we_real,ip,0,we_comm,istat,ierr)
      call mpi_recv(np(-15:15),31,mpi_integer,ip,0,we_comm,istat,ierr)
      n = n + np
     enddo
     ave_E_all(we_np) = ave_E
     ave_ave_E = sum(ave_E_all) / we_np
     sig_ave_E = sqrt(sum((ave_E_all - ave_ave_E)**2)/(we_np-1))
     sigma_all(we_np) = sigma
     ave_sigma = sum(sigma_all) / we_np
     sig_sigma = sqrt(sum((sigma_all - ave_sigma)**2)/(we_np-1))
     long_Nx   = Nx
     long_Ny   = Ny
     long_Nz   = Nz
     long_Np   = we_np
     Ntot = long_Nx * long_Ny * long_Nz * long_Np
     dist(-3) = 100.0d0 * dble(sum(n(-15:-11))) / dble(Ntot)
     dist(-2) = 100.0d0 * dble(sum(n(-10:-6)))  / dble(Ntot)
     dist(-1) = 100.0d0 * dble(sum(n(-5:-1)))   / dble(Ntot)
     dist(0)  = 100.0d0 * dble(n(0))            / dble(Ntot)
     dist(1)  = 100.0d0 * dble(sum(n(1:5)))     / dble(Ntot)
     dist(2)  = 100.0d0 * dble(sum(n(6:10)))    / dble(Ntot)
     dist(3)  = 100.0d0 * dble(sum(n(11:15)))   / dble(Ntot)
     forall (i=1:9) ybox(i) = i * maxval(n) / 9
     do ix=-15,15
       grph_go = .true.
       iy      = 1
       do while(grph_go.and.iy.le.9)
         if(n(ix).le.ybox(iy)) then
           graph(ix,iy)     = grph_top
           graph(ix,1:iy-1) = grph_col
           grph_go          = .false.
         endif
         iy = iy + 1
       enddo 
     enddo
     write(screen,17) 'Average Energy     [eV]       :', ave_ave_E * Eh_to_eV, sig_ave_E * Eh_to_eV
     write(screen,17) 'Width of the DOS   [eV]       :', ave_sigma * Eh_to_eV, sig_sigma * Eh_to_eV
     write(screen,16) 'Total number of energies      :', Ntot
     write(screen,15) graph(-15:-1,9), graph(1:15,9)
     write(screen,14) "Distribution:", graph(-15:-1,8), graph(1:15,8)
     write(screen,10) "Energy Range",      "Ideal", "Actual", graph(-15:-1,7), graph(1:15,7)
     write(screen,11) "-3*sigma to -2*sigma", 2.14, dist(-3), graph(-15:-1,6), graph(1:15,6)
     write(screen,11) "-2*sigma to -sigma",  13.59, dist(-2), graph(-15:-1,5), graph(1:15,5)
     write(screen,11) "-sigma to 0",         34.13, dist(-1), graph(-15:-1,4), graph(1:15,4)
     write(screen,11) "0 to sigma",          34.13, dist(1),  graph(-15:-1,3), graph(1:15,3)
     write(screen,11) "sigma to 2*sigma",    13.59, dist(2),  graph(-15:-1,2), graph(1:15,2)
     write(screen,11) "2*sigma to 3*sigma",   2.14, dist(3),  graph(-15:-1,1), graph(1:15,1)
     write(screen,12) "Out of range",         0.28, dist(0),  "_________________|________________"
     write(screen,13) "-3*s -2*s  -s    0    s   2*s  3*s"
     write(screen,fmdlm1) blank1       
     write(screen,fmdlm1) blank1           
    else
     call mpi_ssend(ave_E,1,we_real,0,0,we_comm,ierr)
     call mpi_ssend(sigma,1,we_real,0,0,we_comm,ierr)
     call mpi_ssend(n(-15:15),31,mpi_integer,0,0,we_comm,ierr)
    endif
10  format(t2,"|",t5,a,t28,a,t43,a,t59,15(a1),"|",15(a1),t100,"|")
11  format(t2,"|",t5,a,t28,f5.2,"%",t43,f5.2,"%",t59,15(a1),"|",15(a1),t100,"|")    
12  format(t2,"|",t5,a,t28,f5.2,"%",t43,f5.2,"%",t57,a,t100,"|")
13  format(t2,"|",t57,a,t100,"|")
14  format(t2,"|",t5,a,t59,15(a1),"|",15(a1),t100,"|")
15  format(t2,"|",t59,15(a1),"|",15(a1),t100,"|")
16  format(t2,"|",t5,a,2x,i0,t74,"|",t100,"|")
17  format(t2,"|",t5,a,2x,f9.4," (",f5.4,")",t100,"|")
    endsubroutine check_distribution
!#####################################################################!
!#### Generate a random list of sites ################################!
!#####################################################################!
    function random_pairs(Nx,Ny,Nz,rho) result(L)
    use general_functions
    integer,          intent(in)       :: Nx,Ny,Nz
    double precision, intent(in)       :: rho
    logical, dimension(Nx,Ny,Nz,2)     :: L
    logical, dimension(Nx,Ny,Nz)       :: Lp
    integer                            :: ix, iy, iz, k, Np
    integer                            :: ixp, ixm, iyp, iym, izp, izm
    logical                            :: redo
    L(:,:,:,1)  = random_sites(Nx,Ny,Nz,rho); L(:,:,:,2)=.false.; Lp = L(:,:,:,1)
    Np = nint(Nx * Ny * Nz * rho) * 3
    do iz=1,Nz; do iy=1,Ny; do ix=1,Nx
      if(L(ix,iy,iz,1)) then
        redo   = .true.
        ixp = ix + 1; if(ixp.gt.Nx) ixp = 1
        ixm = ix - 1; if(ixm.lt.1)  ixm = Nx
        iyp = iy + 1; if(iyp.gt.Ny) iyp = 1
        iym = iy - 1; if(iym.lt.1)  iym = Ny
        izp = iz + 1; if(izp.gt.Nz) izp = 1
        izm = iz - 1; if(izm.lt.1)  izm = Nz   
        if(Lp(ixp,iy,iz).and.Lp(ix,iyp,iz).and.Lp(ix,iy,izp).and. &
        &  Lp(ixm,iy,iz).and.Lp(ix,iym,iz).and.Lp(ix,iy,izm)) then
         L(ix,iy,iz,1)  = .false.
         Lp(ix,iy,iz) = .false.
        else
          do while(redo)
            k   = floor(myrandom() * 6) + 1
            if(    k.eq.1.and..not.Lp(ixp,iy,iz)) then
              Lp(ixp,iy,iz) = .true.
              redo          = .false.
            elseif(k.eq.2.and..not.Lp(ixm,iy,iz)) then
              Lp(ixm,iy,iz) = .true.
              redo          = .false.
            elseif(k.eq.3.and..not.Lp(ix,iyp,iz)) then
              Lp(ix,iyp,iz) = .true.
              redo         = .false.
            elseif(k.eq.4.and..not.Lp(ix,iym,iz)) then
              Lp(ix,iym,iz) = .true.
              redo          = .false.
            elseif(k.eq.5.and..not.Lp(ix,iy,izp)) then
              Lp(ix,iy,izp) = .true.
              redo          = .false.
            elseif(k.eq.6.and..not.Lp(ix,iy,izm)) then
              Lp(ix,iy,izm) = .true.
              redo          = .false.
            endif
          enddo
        endif
      endif
    enddo; enddo; enddo
    do iz=1,Nz; do iy=1,Ny; do ix=1,Nx
     if(Lp(ix,iy,iz).and..not.L(ix,iy,iz,1)) L(ix,iy,iz,2) = Lp(ix,iy,iz)
    enddo; enddo; enddo
    endfunction random_pairs 
!#####################################################################!
!#### Uncorrelated Disorder (Ziggurat algorithm) #####################!
!#####################################################################!
    function uncorrelated(Nx,Ny,Nz,sigma) result(E)
    use general_functions
    implicit none
    integer,          intent(in)          :: Nx, Ny, Nz
    double precision, intent(in)          :: sigma
    double precision, dimension(Nx,Ny,Nz) :: E
    ! Local variables:
    integer,          parameter           :: N      = 256
    double precision, parameter           :: x_max  = 3.6541528853610088d0
    double precision, parameter           :: area   = 0.00492867323399d0
    double precision, dimension(0:256)    :: x, f, k
    integer                               :: ix, iy, iz, ib
    logical                               :: accepted 
    double precision                      :: r0, rx, ry, rs 
    E    = 0.0d0
    x(0) = 0.0d0
    f(0) = 1.0d0
    x(N) = x_max
    f(N) = ziggurat_function(x(N))
    do ix=N-1, 1, -1
      x(ix) = ziggurat_inverse_function(f(ix+1) + area / x(ix+1))
      f(ix) = ziggurat_function(x(ix))
    enddo
    k(0) = x(N) * f(N) / area
    forall(ix=1:N) k(ix) = x(ix) / x(ix-1)
    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          accepted = .false.
          do while(.not. accepted)
            ib = 1 + floor(myrandom() * dble(N)) ! chose a random box between 0 and N_red
            r0 = 2.0d0 * myrandom() - 1.0d0      ! generate random number between -1 and 1
            rs = r0 * x(ib)                      ! pick a random point in the box
            if(abs(r0).lt.k(ib)) then            ! 
              accepted = .true.                  ! leave the loop
            elseif(myrandom() * (f(ib-1) - f(ib)) .lt. (ziggurat_function(rs) - f(ib))) then
              accepted = .true.                  ! leave the loop              
            elseif(ib.eq.0) then                 ! 
              ry = 0.0d0; rx = 1.0d0             ! pick a random point from the tail using
              do while(ry.le.rx**2)              ! the method detailed in:
                rx = -log(myrandom()) / x(N)     ! 1964-Marsaglia-Technometrics-6-101
                ry = -2.0d0 * log(myrandom())    ! 
              enddo                              ! 
              rs = x(N) + rx                     ! 
              accepted = .true.                  ! ... and leave the loop
            endif
          enddo
          E(ix,iy,iz) = sigma * rs
        enddo
      enddo
    enddo
    call check_distribution(Nx,Ny,Nz,E)
    endfunction uncorrelated

!#####################################################################!
!#### Spactially Correlated Disorder #################################!
!#####################################################################!
    function correlated(NN,Nx,Ny,Nz,dipole,lattice,epsilon_r,order_parameter) result(E)
    use type_coor
    use pbc
    use my_mpi, only: i_am_root
    use general_functions
    implicit none
    integer,                       intent(in)   :: NN, Nx, Ny, Nz
    double precision,              intent(in)   :: dipole, lattice, epsilon_r, order_parameter
    double precision, dimension(Nx,Ny,Nz)       :: E
    ! Local Variables
    type(coor), dimension(Nx,Ny,Nz)             :: u
    type(coor), dimension(-NN:NN,-NN:NN,-NN:NN) :: Rij
    type(coor)                                  :: p1, p2, p3, p
    double precision                            :: Rmag
    double precision                            :: crnd, phi, costht
    double precision                            :: factor
    integer                                     :: ix, iy, iz
    integer                                     :: jx, jy, jz
    double precision, dimension(-NN:NN)         :: latt
    integer, dimension(Nx,-NN:NN)               :: kx
    integer, dimension(Ny,-NN:NN)               :: ky
    integer, dimension(Nz,-NN:NN)               :: kz
    !initialize the variables
    E       = 0.0d0
    factor  = - dipole / epsilon_r
    !initialize index arrays
    do ix=1,Nx; do jx=-NN,NN; kx(ix,jx) = indx_periodic_1D(ix+jx,Nx); enddo; enddo
    do iy=1,Ny; do jy=-NN,NN; ky(iy,jy) = indx_periodic_1D(iy+jy,Ny); enddo; enddo
    do iz=1,Nz; do jz=-NN,NN; kz(iz,jz) = indx_periodic_1D(iz+jz,Nz); enddo; enddo
    do ix = -NN, NN; latt(ix) = lattice * dble(ix); enddo
    !set up the antiferromagnetic reference frame
    p1%x = (myrandom() - 0.5d0) * 2.0d0
    p1%y = (myrandom() - 0.5d0) * 2.0d0
    p1%z = (myrandom() - 0.5d0) * 2.0d0
    p1   = p1%norm()
    p2   = p1%perp()    
    p3   = p1.cross.p2  
    !generate the dipole vectors
    if(order_parameter.lt.1.0d0) then
      crnd = 2.0d0 * pi * (1.0d0 - order_parameter)
      do iz = 1, Nz; do iy = 1, Ny; do ix = 1, Nx
       phi           = crnd * myrandom()
       costht        = (myrandom() - 0.5d0) * 2.0d0
       p%x           = sin(phi) * costht
       p%y           = sin(phi) * sqrt(1-costht**2)
       p%z           = cos(phi)
       u(ix,iy,iz)%x = p%x * p1%x + p%y * p2%x + p%z * p3%x
       u(ix,iy,iz)%y = p%x * p1%y + p%y * p2%y + p%z * p3%y
       u(ix,iy,iz)%z = p%x * p1%z + p%y * p2%z + p%z * p3%z
      enddo; enddo; enddo 
    else
       u = p1
    endif
    !apply AF alternation
    do iz = 1, Nz; do iy = 1, Ny; do ix = 1, Nx
      if(mod(ix + iy + iz,2) == 0) then
        u(ix,iy,iz) = u(ix,iy,iz)%neg()
      else
        u(ix,iy,iz) = u(ix,iy,iz)
      end if
    enddo; enddo; enddo
    !set up the grid vectors
    do ix = -NN, NN; do iy = -NN, NN; do iz = -NN, NN
      Rij(ix,iy,iz)%x = latt(ix); Rij(ix,iy,iz)%y = latt(iy); Rij(ix,iy,iz)%z = latt(iz)
      Rmag            = 1.0d0 / Rij(ix,iy,iz)%mag()**3
      Rij(ix,iy,iz)%x = Rij(ix,iy,iz)%x * Rmag
      Rij(ix,iy,iz)%y = Rij(ix,iy,iz)%y * Rmag
      Rij(ix,iy,iz)%z = Rij(ix,iy,iz)%z * Rmag
    enddo; enddo; enddo
    !do the dot products
    do jz= -NN, NN; do jy= -NN, NN; do jx= -NN, NN;
     if(jx.ne.0.or.jy.ne.0.or.jz.ne.0) then
      do iz = 1, Nz; do iy = 1, Ny; do ix = 1, Nx
       E(ix,iy,iz) = E(ix,iy,iz) + (u(kx(ix,jx),ky(iy,jy),kz(iz,jz)) * Rij(jx,jy,jz))
      enddo; enddo; enddo
     endif
    enddo; enddo; enddo
    E     = factor * E
    if(i_am_root) then
      if(order_parameter.eq.1.0) then
        write(screen,fmchr1) "Dipoles uniformly distributed antiferromagnetically"
      elseif(order_parameter.eq.0.0) then
        write(screen,fmchr1) "Dipoles randomly distributed"
      else
        write(screen,fmchfl) "Dipoles have an order paramter of", order_parameter
      endif
      write(screen,11) "Dipole Magnetude [D]          :", dipole * ea0_to_D
      call print_3D_matrix(Nx,Ny,Nz,E,'CHECK_E_dipole2.dat')
    endif
    call check_distribution(Nx,Ny,Nz,E)
11  format(t2,"|",t5,a,2x,f9.4,t100,"|")
    endfunction correlated

!#####################################################################!
!#### Spactially Correlated Disorder with Random Dimers ##############!
!#####################################################################!
    function correlated_dimers(NN,Nx,Ny,Nz,dipole,lattice,epsilon_r,dimer_dipole,dimer_rho) result(E)
    use type_coor
    use pbc
    use my_mpi
    use general_functions
    implicit none
    integer,                       intent(in)   :: NN, Nx, Ny, Nz
    double precision,              intent(in)   :: dipole, lattice, epsilon_r, dimer_dipole, dimer_rho
    double precision, dimension(Nx,Ny,Nz)       :: E, dE
    ! Local Variables
    type(coor), dimension(Nx,Ny,Nz)             :: u
    type(coor), dimension(-NN:NN,-NN:NN,-NN:NN) :: Rij
    type(coor)                                  :: p1, p2, p3, p
    logical,    dimension(Nx,Ny,Nz,2)           :: dimer_list
    double precision                            :: dimer_swap
    double precision                            :: Rmag
    double precision                            :: crnd, phi, costht
    double precision                            :: factor
    integer                                     :: ix, iy, iz, ip
    integer                                     :: jx, jy, jz
    double precision, dimension(-NN:NN)         :: latt
    integer, dimension(Nx,-NN:NN)               :: kx
    integer, dimension(Ny,-NN:NN)               :: ky
    integer, dimension(Nz,-NN:NN)               :: kz
    double precision                            :: new_rho, all_rho(0:we_np-1), std_rho
    !initialize the variables
    E          = 0.0d0
    factor     = - 1.d0 / epsilon_r
    new_rho    = dimer_rho / 2.d0
    dimer_swap = dimer_dipole / dipole
    !initialize index arrays
    do ix=1,Nx; do jx=-NN,NN; kx(ix,jx) = indx_periodic_1D(ix+jx,Nx); enddo; enddo
    do iy=1,Ny; do jy=-NN,NN; ky(iy,jy) = indx_periodic_1D(iy+jy,Ny); enddo; enddo
    do iz=1,Nz; do jz=-NN,NN; kz(iz,jz) = indx_periodic_1D(iz+jz,Nz); enddo; enddo
    do ix = -NN, NN; latt(ix) = lattice * dble(ix); enddo
    !set up the antiferromagnetic reference frame
    p1%x = (myrandom() - 0.5d0) * 2.0d0
    p1%y = (myrandom() - 0.5d0) * 2.0d0
    p1%z = (myrandom() - 0.5d0) * 2.0d0
    p1   = p1%norm()
    p2   = p1%perp()    
    p3   = p1.cross.p2  
    !generate the dipole vectors
    crnd = 2.0d0 * pi
    do iz = 1, Nz; do iy = 1, Ny; do ix = 1, Nx
     phi           = crnd * myrandom()
     costht        = (myrandom() - 0.5d0) * 2.0d0
     p%x           = sin(phi) * costht
     p%y           = sin(phi) * sqrt(1-costht**2)
     p%z           = cos(phi)
     u(ix,iy,iz)%x = (p%x * p1%x + p%y * p2%x + p%z * p3%x) * dipole
     u(ix,iy,iz)%y = (p%x * p1%y + p%y * p2%y + p%z * p3%y) * dipole
     u(ix,iy,iz)%z = (p%x * p1%z + p%y * p2%z + p%z * p3%z) * dipole
    enddo; enddo; enddo 
    !apply AF alternation
    do iz = 1, Nz; do iy = 1, Ny; do ix = 1, Nx
      if(mod(ix + iy + iz,2) == 0) then
        u(ix,iy,iz) = u(ix,iy,iz)%neg()
      else
        u(ix,iy,iz) = u(ix,iy,iz)
      end if
    enddo; enddo; enddo
    !generate the dimers
    dimer_list = random_pairs(Nx,Ny,Nz,new_rho)
    new_rho    = dble(count(dimer_list)) / dble(Nx*Ny*Nz)
    do iz = 1, Nz; do iy = 1, Ny; do ix = 1, Nx
      if(dimer_list(ix,iy,iz,1)) then
        u(ix,iy,iz)%x = u(ix,iy,iz)%x * dimer_swap
        u(ix,iy,iz)%y = u(ix,iy,iz)%y * dimer_swap 
        u(ix,iy,iz)%z = u(ix,iy,iz)%z * dimer_swap
      elseif(dimer_list(ix,iy,iz,2)) then
        u(ix,iy,iz)%x = 0.d0
        u(ix,iy,iz)%y = 0.d0
        u(ix,iy,iz)%z = 0.d0
      endif
    enddo; enddo; enddo 
    !set up the grid vectors
    do ix = -NN, NN; do iy = -NN, NN; do iz = -NN, NN
      Rij(ix,iy,iz)%x = latt(ix); Rij(ix,iy,iz)%y = latt(iy); Rij(ix,iy,iz)%z = latt(iz)
      Rmag            = 1.0d0 / Rij(ix,iy,iz)%mag()**3
      Rij(ix,iy,iz)%x = Rij(ix,iy,iz)%x * Rmag
      Rij(ix,iy,iz)%y = Rij(ix,iy,iz)%y * Rmag
      Rij(ix,iy,iz)%z = Rij(ix,iy,iz)%z * Rmag
    enddo; enddo; enddo
    !do the dot products !! This is the expensive bit, but I don't think it can go any faster :(
    do jz= -NN, NN; do jy= -NN, NN; do jx= -NN, NN;
     if(jx.ne.0.or.jy.ne.0.or.jz.ne.0) then
      do iz = 1, Nz; do iy = 1, Ny; do ix = 1, Nx
       E(ix,iy,iz) = E(ix,iy,iz) + (u(kx(ix,jx),ky(iy,jy),kz(iz,jz)) * Rij(jx,jy,jz))
      enddo; enddo; enddo
     endif
    enddo; enddo; enddo
    E     = factor * E
    if(i_am_root) then
      write(screen,fmchr1) "Dimers randomly distributed"
      write(screen,11) "Monomer Dipole Magnetude [D]  :", dipole * ea0_to_D
      write(screen,11) "Dimer Dipole Magnetude [D]    :", dimer_dipole * ea0_to_D
      all_rho(0) = new_rho
      do ip=1,we_np-1
       call mpi_recv(all_rho(ip),1,we_real,ip,0,we_comm,istat,ierr)
      enddo
      new_rho = sum(all_rho) / we_np   
      std_rho = sqrt(sum((all_rho - new_rho)**2)/(we_np-1))
      write(screen,10) "System Dimer Density [/nm]    :", new_rho * lattice * a0_to_nm, std_rho * lattice * a0_to_nm
      call print_3D_matrix(Nx,Ny,Nz,E,'CHECK_E_dipole2.dat')
    else
      call mpi_ssend(new_rho,1,we_real,0,0,we_comm,ierr)
    endif
    call check_distribution(Nx,Ny,Nz,E)
    !temporary site analysis
    if(i_am_root) write(screen,fmchr1) "Energy Difference Analysis:"
    dE = 0.0d0
    do jz= -1, 1; do jy= -1, 1; do jx= -1, 1;
     if(jx.ne.0.or.jy.ne.0.or.jz.ne.0) then
      do iz = 1, Nz; do iy = 1, Ny; do ix = 1, Nx
       dE(ix,iy,iz) = dE(ix,iy,iz) + (E(ix,iy,iz) - E(kx(ix,jx),ky(iy,jy),kz(iz,jz)))
      enddo; enddo; enddo
     endif
    enddo; enddo; enddo
    dE = dE / 26.d0
    call check_distribution(Nx,Ny,Nz,dE)    
10  format(t2,"|",t5,a,2x,f9.4," (",f5.4,")",t100,"|")
11  format(t2,"|",t5,a,2x,f9.4,t100,"|")
    endfunction correlated_dimers
!#####################################################################################################

end module random

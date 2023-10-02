module general_functions
contains
!---------------------------------------------------------------------!
! findeta:                                                            !
!      search for the first number above a threshold in a cummulative !
!      sum array using a binary search                                !
!---------------------------------------------------------------------!
    function findeta(s,eta) result(i)
    implicit none
    double precision, dimension(:), intent(in) :: s
    double precision,               intent(in) :: eta
    integer                                    :: i
    integer                                    :: u,l,m       
    if(eta.lt.s(2)) then
     i = 1
    elseif(eta.gt.s(size(s)-1)) then
     i = size(s)
    else
     l=0; u=size(s)
     do while(u-l.gt.1) 
      m = (u + l) / 2
      if(s(m).lt.eta) then
       l = m
      else
       u = m
      endif
     enddo
     if(s(l).gt.eta) then
      i = l
     else
      i = u
     endif
    endif
    endfunction findeta
!---------------------------------------------------------------------!
! myrandom:                                                           !
!      intrinsic random number routine repackaged the into a function !
!---------------------------------------------------------------------!
    function myrandom() result(r)
    double precision :: r
    call random_number(r)
    endfunction myrandom
!---------------------------------------------------------------------!
! ziggurat_function:                                                  !
!      return the result of the ziggurat function                     !
!---------------------------------------------------------------------!
    function ziggurat_function(x) result(f)
    implicit none
    double precision, intent(in) :: x
    double precision             :: f
    f = exp(-0.5d0 * x**2)
    endfunction ziggurat_function
!---------------------------------------------------------------------!
! ziggurat_inverse_function:                                          !
!      return the result of the inverse ziggurat function             !
!---------------------------------------------------------------------!
    function ziggurat_inverse_function(f) result(x)
    implicit none
    double precision, intent(in) :: f
    double precision             :: x
    x = sqrt(-2.0d0 * log(f))
    endfunction ziggurat_inverse_function
!---------------------------------------------------------------------!
! cumulative_sum:                                                     !
!      create a cumulative sum array from the input array             !
!---------------------------------------------------------------------!    
    function cumulative_sum(p) result(s)
    implicit none
    double precision, dimension(:), intent(in) :: p
    double precision, dimension(size(p))       :: s
    integer                                    :: i
    s(1) = p(1)
    do i=2, size(p)
     s(i) = s(i-1) + p(i)
    enddo
    s = s / s(size(p))
    endfunction cumulative_sum
!---------------------------------------------------------------------!
! random_sites:                                                       !
!      create a random list of sites                                  !
!---------------------------------------------------------------------!
    function random_sites(Nx,Ny,Nz,rho) result(L)
    integer,          intent(in)       :: Nx,Ny,Nz
    double precision, intent(in)       :: rho
    logical, dimension(Nx,Ny,Nz)       :: L
    integer                            :: N, Np, i, j
    logical, allocatable, dimension(:) :: L1D
    logical                            :: redo
    N  = Nx * Ny * Nz
    Np = nint(N * rho)
    allocate(L1D(N)); L1D=.false.
    do i=1,Np
     redo = .true.
     do while(redo)
       j = floor(myrandom() * N) + 1
       if(.not.L1D(j)) then
        L1D(j) = .true.
        redo   = .false.
       endif
     enddo
    enddo
    L = reshape(L1D,(/Nx,Ny,Nz/))   
    endfunction random_sites 
!---------------------------------------------------------------------!
! random_noise:                                                       !
!      create a 1D array of random gaussian noise of specified width  !
!---------------------------------------------------------------------!
    function random_noise(NA,sigma) result(A)
    implicit none
    integer,          intent(in)     :: NA
    double precision, intent(in)     :: sigma
    double precision, dimension(NA)  :: A
    ! Local variables:
    integer,          parameter      :: N      = 256
    double precision, parameter      :: x_max  = 3.6541528853610088d0
    double precision, parameter      :: area   = 0.00492867323399d0
    double precision, dimension(0:N) :: x, f, k
    integer                          :: ia, ib
    logical                          :: accepted 
    double precision                 :: r0, rx, ry, rs 
    A    = 0.0d0
    x(0) = 0.0d0
    f(0) = 1.0d0
    x(N) = x_max
    f(N) = ziggurat_function(x(N))
    do ia=N-1, 1, -1
      x(ia) = ziggurat_inverse_function(f(ia+1) + area / x(ia+1))
      f(ia) = ziggurat_function(x(ia))
    enddo
    k(0) = x(N) * f(N) / area
    forall(ia=1:N) k(ia) = x(ia) / x(ia-1)
    do ia = 1, NA
     accepted = .false.
     do while(.not. accepted)
      ib = 1 + floor(myrandom() * dble(N)) ! chose a random box between 0 and N_red
      r0 = 2.0d0 * myrandom() - 1.0d0      ! generate random number between -1 and 1
      rs = r0 * x(ib)                      ! pick a random point in the box
      if(abs(r0).lt.k(ib)) then            ! 
       accepted = .true.                   ! leave the loop
      elseif(myrandom() * (f(ib-1) - f(ib)) .lt. (ziggurat_function(rs) - f(ib))) then
       accepted = .true.                   ! leave the loop              
      elseif(ib.eq.0) then                 ! 
       ry = 0.0d0; rx = 1.0d0              ! pick a random point from the tail using
       do while(ry.le.rx**2)               ! the method detailed in:
        rx = -log(myrandom()) / x(N)       ! 1964-Marsaglia-Technometrics-6-101
        ry = -2.0d0 * log(myrandom())      ! 
       enddo                               ! 
       rs = x(N) + rx                      ! 
       accepted = .true.                   ! ... and leave the loop
      endif
     enddo
     A(ia) = sigma * rs
    enddo
    endfunction random_noise
!---------------------------------------------------------------------!
endmodule general_functions


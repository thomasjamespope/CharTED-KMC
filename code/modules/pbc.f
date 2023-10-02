!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% module pbc.f90: ...description...
!%
!% Copyright (C) 2020 Yvelin Giret & Thomas Pope
!% \author Yvelin Giret & Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    module pbc
    implicit none
    contains
!---------------------------------------------------------------------!
! indx_periodic_1D:                                                   !
!      Adjust an index, i, so that it falls between 1:N assuming a    !
!      periodic structure                                             !
!---------------------------------------------------------------------!
    function indx_periodic_1D(i,N) result(j)
    implicit none
    integer, intent(in) :: i, N
    integer             :: j
    j = i
    if(i.lt.1) then
      if(i.gt.-N) then
        j = i + N 
      elseif(i.gt.-2*N) then
        j = i + 2 * N
      elseif(i.gt.-3*N) then
        j = i + 3 * N
      endif
    elseif(i.gt.N) then
      if(i.lt.2*N) then
        j = i - N 
      elseif(i.lt.3*N) then
        j = i - 2 * N 
      elseif(i.lt.4*N) then
        j = i - 3 * N 
      endif
    endif
    endfunction indx_periodic_1D
!---------------------------------------------------------------------!
! indx_periodic_1D:                                                   !
!      Adjust an index data type, i, so that it falls on a grid point !
!      in a box defined by the dimnesions Nx,Ny,Nz, assuming the box  !
!      is periodic                                                    !
!---------------------------------------------------------------------!
    function indx_periodic(i) result(j)
    use type_indexes
    use param, only: Nx, Ny, Nz
    implicit none
    type(indexes), intent(in)  :: i
    type(indexes)              :: j
    j%x = indx_periodic_1D(i%x,Nx)
    j%y = indx_periodic_1D(i%y,Ny)
    j%z = indx_periodic_1D(i%z,Nz)
    endfunction indx_periodic

    endmodule pbc

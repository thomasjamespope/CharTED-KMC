!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% type_indexes.f
!%
!% module type_indexes: 
!%  Definition of the indexes derived data type
!%  and its associated procedures
!%
!% Copyright (C) 2020 Yvelin Giret & Thomas Pope
!% \author Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module type_indexes
    use kinds
!********************************************************************************************!
! indexes data type:                                                                         !
! contains three integers (x,y,z), which reference a point on a 3-dimensional grid           !
!********************************************************************************************!
    type indexes
      integer   :: x,y,z
    contains
      procedure :: mag => indexes_mag
      procedure :: indexes_sum
      procedure :: indexes_min
      procedure :: indexes_vec
      procedure :: index_matrix
      procedure :: indexes_dot
      procedure :: indexes_nearto
      procedure :: indexes_scalar_projection
      generic   :: operator(+)        => indexes_sum
      generic   :: operator(-)        => indexes_min
      generic   :: operator(.vec.)    => indexes_vec
      generic   :: operator(.of.)     => index_matrix
      generic   :: operator(.dot.)    => indexes_dot
      generic   :: operator(.nearto.) => indexes_nearto
      generic   :: operator(.sp.)     => indexes_scalar_projection
    endtype indexes
    contains
!--------------------------------------------------------------------------------------------!
! indexes_mag (remampped to the "%mag" operator for indexes data types):                     !
!     gives the magnetude of a vector defines by the x, y and z components of an index       !
! usage, for indexes tybe variable a                                                         !
!     a%mag                                                                                  !
!--------------------------------------------------------------------------------------------!
    function indexes_mag(a) result(r)
    implicit none
    class(indexes), intent(in)  :: a
    double precision                    :: r
    r = sqrt(dble(a%x**2 + a%y**2 + a%z**2))
    endfunction indexes_mag
!--------------------------------------------------------------------------------------------!
! indexes_sum (remampped to the "+" operator for indexes data types):                        !
!     sums the x,y and z components of two indexes type variables                            !
! usage, for indexes type varbiables a,b and c                                               !
!     c = a + b                                                                              !
!--------------------------------------------------------------------------------------------!
    function indexes_sum(a,b) result(c)
    implicit none
    class(indexes), intent(in)  :: a, b
    type(indexes)               :: c
    c%x = a%x + b%x
    c%y = a%y + b%y
    c%z = a%z + b%z
    endfunction indexes_sum
!--------------------------------------------------------------------------------------------!
! indexes_min (remampped to the "-" operator for indexes data types):                        !
!     takes the difference of the x,y and z components of two indexes type variables         !
! usage, for indexes type varbiables a,b and c                                               !
!     c = a - b                                                                              !
!--------------------------------------------------------------------------------------------!
    function indexes_min(a,b) result(c)
    implicit none
    class(indexes), intent(in)  :: a, b
    type(indexes)               :: c
    c%x = a%x - b%x
    c%y = a%y - b%y
    c%z = a%z - b%z
    endfunction indexes_min
!--------------------------------------------------------------------------------------------!
! indexes_vec (remampped to the ".vec." operator for indexes data types):                    !
!     takes the difference of the x,y and z components of two indexes type variables         !
! usage, for indexes type varbiables a,b and c                                               !
!     c = a.vec.b                                                                            !
!--------------------------------------------------------------------------------------------!
    function indexes_vec(a,b) result(c)
    use param, only: Nx, Ny, Nz
    implicit none
    class(indexes), intent(in)  :: a, b
    type(indexes)               :: c
    c%x = a%x - b%x; c%x = min(abs(c%x), abs(c%x - Nx), abs(c%x + Nx)) + 1 
    c%y = a%y - b%y; c%y = min(abs(c%y), abs(c%y - Ny), abs(c%y + Ny)) + 1
    c%z = a%z - b%z; c%z = min(abs(c%z), abs(c%z - Nz), abs(c%z + Nz)) + 1
    endfunction indexes_vec
!--------------------------------------------------------------------------------------------!
! index_matrix (remampped to the ".of." operator for indexes and charge data types):         !
!     gives the energy value for a charge data type at an indexes data type grid point       !
! usage, for indexes type varbiable a, real 3d array variable b, and real variable c         !
!     c = a.of.b                                                                             !
!--------------------------------------------------------------------------------------------!
    function index_matrix(i,M) result(r)
    use param, only: Nx, Ny, Nz
    implicit none
    class(indexes),                        intent(in) :: i
    double precision, dimension(Nx,Ny,Nz), intent(in) :: M
    double precision                                  :: r
    r = M(i%x,i%y,i%z)
    endfunction index_matrix
!--------------------------------------------------------------------------------------------!
! indexes_dot (remampped to the ".dot." operator for indexes data types):                    !
!     finds the dot product of two indexes                                                   !
! usage, for indexes type varbiables a and b, and a real variable c                          !
!     c = a.dot.b                                                                            !
!--------------------------------------------------------------------------------------------!
    function indexes_dot(a,b) result(c)
    implicit none
    class(indexes), intent(in) :: a, b
    double precision           :: c
    c = a%x * b%x + a%y * b%y + a%z * b%z
    endfunction indexes_dot
!--------------------------------------------------------------------------------------------!
! indexes_distance (remampped to the ".nearto." operator for indexes data types):            !
!     finds the distance between two indexes                                                 !
! usage, for indexes type varbiables a and b, and a real variable c                          !
!     c = a.nearto.b                                                                         !
!--------------------------------------------------------------------------------------------!
    function indexes_nearto(a,b) result(near)
    use param, only: coul_cutoff, Nx, Ny, Nz
    implicit none
    class(indexes), intent(in) :: a, b
    type(indexes)              :: c
    logical                    :: near
    c%x = a%x - b%x; c%x = min(abs(c%x), abs(c%x - Nx), abs(c%x + Nx))
    if(c%x.lt.coul_cutoff) then
     c%y = a%y - b%y; c%y = min(abs(c%y), abs(c%y - Ny), abs(c%y + Ny))
     if(c%y.lt.coul_cutoff) then
      c%z = a%z - b%z; c%z = min(abs(c%z), abs(c%z - Nz), abs(c%z + Nz))
      if(c%z.lt.coul_cutoff) then
       if(c%mag().lt.coul_cutoff) then
        near=.true.
       else
        near=.false.
       endif
      else
       near=.false.
      endif
     else
      near=.false.
     endif
    else
     near=.false.
    endif
    endfunction indexes_nearto
!--------------------------------------------------------------------------------------------!
! indexes_scalar_projection (remampped to the ".sp." operator for indexes data types):       !
!     finds the scalar projection of one index onto another                                  !
! usage, for indexes type varbiables a and b, and a real variable c                          !
!     c = a.sp.b                                                                             !
!--------------------------------------------------------------------------------------------!
    function indexes_scalar_projection(a,b) result(c)
    implicit none
    class(indexes), intent(in) :: a, b
    double precision                   :: c
    c = (a.dot.b) / b%mag()
    endfunction indexes_scalar_projection
endmodule type_indexes

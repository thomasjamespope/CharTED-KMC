!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% type_coor.f
!%
!% module type_coor: 
!%  Definition of the coor derived data type
!%  and its associated procedures
!%
!% Copyright (C) 2020 Yvelin Giret & Thomas Pope
!% \author Thomas Pope
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module type_coor
    use kinds
!********************************************************************************************!
! coor data type:                                                                            !
! contains three reals (x,y,z), which reference a vector in 3-dimensional space              !
!********************************************************************************************!
    type coor
      double precision  :: x,y,z
    contains
      procedure :: mag  => coor_mag
      procedure :: norm => coor_normalize
      procedure :: neg  => coor_negative
      procedure :: rot  => coor_rotate
      procedure :: perp => coor_perpendicular
      procedure :: coor_dot
      procedure :: coor_scalar_projection
      procedure :: coor_sum
      procedure :: coor_min
      procedure :: coor_cross
      generic   :: operator(+)       => coor_sum
      generic   :: operator(-)       => coor_min
      generic   :: operator(*)       => coor_dot
      generic   :: operator(.cross.) => coor_cross
      generic   :: operator(.sp.)    => coor_scalar_projection
    endtype coor
    contains
!--------------------------------------------------------------------------------------------!
! coor_mag (remampped to the "%mag" operator for coor data types):                           !
!     gives the magnetude of a vector defines by the x, y and z components of an coor        !
! usage, for coor tybe variable a                                                            !
!     a%mag                                                                                  !
!--------------------------------------------------------------------------------------------!
    function coor_mag(a) result(r)
    implicit none
    class(coor), intent(in)  :: a
    double precision                 :: r
    r = sqrt(a%x**2 + a%y**2 + a%z**2)
    endfunction coor_mag
!--------------------------------------------------------------------------------------------!
! coor_normalize (remampped to the "%norm" operator for coor data types):                    !
!     normalizes the vector defined by the x, y and z components of an coor type             !
! usage, for coor type variables a and b                                                     !
!     b=a%norm()                                                                             !
!--------------------------------------------------------------------------------------------!
    function coor_normalize(a) result(b)
    implicit none
    class(coor), intent(in)  :: a
    type(coor)               :: b
    double precision         :: nrm
    nrm = a%mag()
    b%x = a%x / nrm
    b%y = a%y / nrm
    b%z = a%z / nrm
    endfunction coor_normalize
!--------------------------------------------------------------------------------------------!
! coor_sum (remampped to the "+" operator for coor data types):                              !
!     sums the x,y and z components of two coor type variables                               !
! usage, for coor type varbiables a,b and c                                                  !
!     c = a + b                                                                              !
!--------------------------------------------------------------------------------------------!
    function coor_sum(a,b) result(c)
    implicit none
    class(coor), intent(in)  :: a, b
    type(coor)               :: c
    c%x = a%x + b%x
    c%y = a%y + b%y
    c%z = a%z + b%z
    end function coor_sum
!--------------------------------------------------------------------------------------------!
! coor_min (remampped to the "-" operator for coor data types):                              !
!     takes the difference of the x,y and z components of two coor type variables            !
! usage, for coor type varbiables a,b and c                                                  !
!     c = a - b                                                                              !
!--------------------------------------------------------------------------------------------!
    function coor_min(a,b) result(c)
    implicit none
    class(coor), intent(in)  :: a, b
    type(coor)               :: c
    c%x = a%x - b%x
    c%y = a%y - b%y
    c%z = a%z - b%z
    endfunction coor_min
!--------------------------------------------------------------------------------------------!
! coor_dot (remampped to the "*" operator for coor data types):                              !
!     finds the dot product of two coor                                                      !
! usage, for coor type varbiables a and b, and a real variable c                             !
!     c = a * b                                                                              !
!--------------------------------------------------------------------------------------------!
    function coor_dot(a,b) result(c)
    implicit none
    class(coor), intent(in) :: a, b
    double precision                   :: c
    c = a%x * b%x + a%y * b%y + a%z * b%z
    endfunction coor_dot
!--------------------------------------------------------------------------------------------!
! coor_cross (remampped to the ".cross." operator for coor data types):                      !
!     finds the cross product of two coor                                                    !
! usage, for coor type varbiables a, b, and c                                                !
!     c = a.cross.b                                                                          !
!--------------------------------------------------------------------------------------------!
    function coor_cross(a,b) result(c)
    implicit none
    class(coor), intent(in) :: a, b
    type(coor)              :: c
    c%x = a%y * b%z - a%z * b%y
    c%y = a%z * b%x - a%x * b%z
    c%z = a%x * b%y - a%y * b%x
    c   = c%norm()
    endfunction coor_cross
!--------------------------------------------------------------------------------------------!
! coor_perpendicular (remampped to the "%perp" operator for coor data types):                !
!     finds an arbitrary vector perpendicular to an input vector                             !
! usage, for coor type varbiables a, and b                                                   !
!     c = a%perp()                                                                           !
!--------------------------------------------------------------------------------------------!
    function coor_perpendicular(a) result(b)
    implicit none
    class(coor), intent(in) :: a
    type(coor)              :: b,c
    if(a%x.eq.1.0d0.and.a%y.eq.0.0d0.and.a%z.eq.0.0d0) then
      b%x = 0.0d0; b%y = 1.0d0; b%z = 0.0d0
    else
      c%x = 1.0d0; c%y = 0.0d0; c%z = 0.0d0
      b = c.cross.a
    endif
    endfunction coor_perpendicular
!--------------------------------------------------------------------------------------------!
! coor_rotate (remampped to the "%rot" operator for coor data types):                        !
!    returns a vector rotated through an angle theta around the axis vector ax               !
! usage, for coor type varbiables a, b and c and real variable theta                         !
!     c = a%rot(b,theta)                                                                       !
!--------------------------------------------------------------------------------------------!
    function coor_rotate(v0,ax,theta) result(v)
    implicit none
    class(coor), intent(in)  :: v0
    type(coor),  intent(in)  :: ax
    double precision,    intent(in)  :: theta 
    double precision                 :: cth, sth
    double precision, dimension(3,3) :: R
    type(coor)               :: v
    cth    = cos(theta)
    sth    = sin(theta)
    R(1,1) = cth + ax%x**2 * (1 - cth)
    R(1,2) =   ax%x * ax%y * (1 - cth) - ax%z * sth
    R(1,3) =   ax%x * ax%z * (1 - cth) + ax%y * sth
    R(2,1) =   ax%x * ax%y * (1 - cth) + ax%z * sth
    R(2,2) = cth + ax%y**2 * (1 - cth)
    R(2,3) =   ax%y * ax%z * (1 - cth) - ax%x * sth
    R(3,1) =   ax%x * ax%z * (1 - cth) - ax%y * sth
    R(3,2) =   ax%y * ax%z * (1 - cth) + ax%x * sth
    R(3,3) = cth + ax%z**2 * (1 - cth)
    v      = coormul(R,v0)
    endfunction coor_rotate
!--------------------------------------------------------------------------------------------!
! coor_negative (remampped to the "%neg" operator for coor data types):                      !
!    returns the negative vector of a coor type variable                                     !
! usage, for coor type varbiables a and b                                                    !
!     b = a%neg()                                                                            !
!--------------------------------------------------------------------------------------------!
    function coor_negative(a) result(b)
    implicit none
    class(coor), intent(in)  :: a
    type(coor)               :: b
    b%x = -a%x
    b%y = -a%y
    b%z = -a%z
    endfunction coor_negative
!--------------------------------------------------------------------------------------------!
! coor_scalar_projection (remampped to the ".sp." operator for coor data types):             !
!     finds the scalar projection of one index onto another                                  !
! usage, for coor type varbiables a and b, and a real variable c                             !
!     c = a.sp.b                                                                             !
!--------------------------------------------------------------------------------------------!
    function coor_scalar_projection(a,b) result(c)
    implicit none
    class(coor), intent(in) :: a, b
    double precision                   :: c
    c = (a * b) / b%mag()
    endfunction coor_scalar_projection
!--------------------------------------------------------------------------------------------!
! coor_R_perturbation (remampped to the ".sR." operator for coor data types):                !
!     finds the perturbation parameter for a shifts index point                              !
! usage, for coor type varbiables a, an indexes type b, and a real variable c                !
!     c = a.sp.b                                                                             !
!--------------------------------------------------------------------------------------------!
    function coor_R_perturbation(d,I) result(p)
    use param, only: lattice, local
    use type_indexes
    implicit none
    class(coor), intent(in)   :: d
    type(indexes), intent(in) :: I
    type(coor)                :: R,Rp
    double precision                  :: p
    R%x = I%x * lattice; R%y = I%y * lattice; R%z = I%z * lattice
    Rp  = R - d
    p   = exp( - local * Rp%mag() )
    endfunction coor_R_perturbation
!--------------------------------------------------------------------------------------------!
! coormul                                                                                    !
!     multiplies a vector in the coor format by a 3x3 matrix                                 !
! usage, for coor type varbiables a and b, and a 3x3 array M                                 !
!     b = coormul(M,a)                                                                       !
!--------------------------------------------------------------------------------------------!
    function coormul(M,vi) result(vo)
    implicit none
    double precision, dimension(3,3), intent(in) :: M
    type(coor),               intent(in) :: vi
    type(coor)                           :: vo
    vo%x = M(1,1) * vi%x + M(1,2) * vi%y + M(1,3) * vi%z
    vo%y = M(2,1) * vi%x + M(2,2) * vi%y + M(2,3) * vi%z
    vo%z = M(3,1) * vi%x + M(3,2) * vi%y + M(3,3) * vi%z
    endfunction coormul
endmodule type_coor

!######################################################################################################
! Point2DModule
! type Point2D
! membres :
!       real(kind=8) :: x
!       real(kind=8) :: y
!
! subroutine :
!       copy_Point2D : permet de définir l'opérateur d'affection pour le type
!######################################################################################################


module Point2DModule

    implicit none
    type, public :: Point2D
        real(kind=8) :: x
        real(kind=8) :: y
    contains
    end type Point2D

    interface assignment(=)
        module procedure copy_Point2D
    end interface
contains
  subroutine copy_Point2D(point_left, point_right)
        type(Point2D), intent(out) :: point_left
        type(Point2D), intent(in) :: point_right

        point_left%x = point_right%x
        point_left%y = point_right%y
    end subroutine copy_Point2D

end module Point2DModule

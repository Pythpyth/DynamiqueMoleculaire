!######################################################################################################
! ForceCalculatorModule
! type ForceCalculator
! membres :
!       real(kind=8) :: box_size        : la taille de la boite
!       real(kind=8) :: half_box_size   : la moitié de la taille de la boite
!       real(kind=8) :: r_truncated     : le rayon de coupure
!
! subroutine :
!       initialize : initialise les membres
!
! fonctions :
!
! compute_sum_forces_on_reference_particle : calcule la somme des forces sur la particule d'indice
!                                            position_reference_index dans le vecteur des particules
!
! compute_force_on_reference_particle : calcule la force exercée par la particule en position_2
!                                       sur la particule en position_reference
!
!######################################################################################################


module ForceCalculatorModule
    use Point2DModule
    use DistanceCalculatorModule
    implicit none

    type, public :: ForceCalculator
        private
        real(kind=8) :: box_size
        real(kind=8) :: half_box_size
        real(kind=8) :: r_truncated
    contains
        procedure, public :: initialize, compute_sum_forces_on_reference_particle
    end type ForceCalculator

contains

    subroutine initialize(this, box_size, r_truncated)
        class(ForceCalculator), intent(inout) :: this
        real(kind=8), intent(in):: box_size, r_truncated

        this%box_size = box_size
        this%half_box_size = box_size / 2.0d0
        this%r_truncated = r_truncated
    end subroutine


    function compute_sum_forces_on_reference_particle(this, position_reference, position_reference_index, &
                                                      positions, nb_particle)

        class(ForceCalculator), intent(in) :: this
        type(Point2D),  intent(in):: position_reference
        integer(kind=4), intent(in) :: nb_particle, position_reference_index
        type(Point2D), dimension(nb_particle), intent(in):: positions
        type(Point2D) :: compute_sum_forces_on_reference_particle
        type(Point2D) :: force_i, sum_forces
        integer(kind=4) :: i

        sum_forces%x = 0.0d0
        sum_forces%y = 0.0d0

        do i = 1, nb_particle

            if ( i/= position_reference_index ) then
                force_i = compute_force_on_reference_particle(this, position_reference, positions(i))

                sum_forces%x = sum_forces%x + force_i%x
                sum_forces%y = sum_forces%y + force_i%y
            end if

        end do

        compute_sum_forces_on_reference_particle = sum_forces

    end function


    function compute_force_on_reference_particle(this, position_reference, position_2)
        class(ForceCalculator), intent(in) :: this
        type(Point2D),  intent(in):: position_reference, position_2
        type(Point2D) :: compute_force_on_reference_particle
        type(Point2D) :: closest_image_position_2
        real(kind=8) :: distance, common_component, x_force, y_force


        closest_image_position_2 = compute_closest_image_position_2(&
        position_reference, position_2, this%half_box_size, this%box_size)

        distance = compute_distance(position_reference, closest_image_position_2)

        if (distance.gt.this%r_truncated) then
            compute_force_on_reference_particle = Point2D(0.0d0, 0.0d0)
        else
            common_component = 48.0d0 * (distance ** (-14)) - 24.0d0 * (distance ** (-8))

            x_force = common_component * (position_reference%x - closest_image_position_2%x)
            y_force = common_component * (position_reference%y - closest_image_position_2%y)

            compute_force_on_reference_particle = Point2D(x_force, y_force)

        end if

    end function

end module ForceCalculatorModule




!######################################################################################################
! VerletPropagatorModule
! type VerletPropagator
!        real(kind=8) :: box_size       : taille de la boite
!        real(kind=8) :: half_box_size  : la moitié de la taille de la boite
!        real(kind=8) :: delta_t        : le pas de temps
!
! subroutine :
!       initialize : initialise les membres
!
! fonction :
!       apply_position_boundary_conditions : on considère une position après une évolution sur un pas de temps
!                                            dans le cas ou la position est sortie de la boite,
!                                            on retourne la position retourne la position image dans la boite
!                                            sinon on retourne la position
!
!       evolve_position_1_dimension : évolution 1 dimension de la position avec l'algorithme Velocity-Verlet
!       evolve_position : évolution de la position avec l'algorithme Velocity-Verlet
!
!       evolve_velocity_1_dimension : évolution de la vitesse avec l'algorithme Velocity-Verlet
!       evolve_velocity : évolution de la vitesse avec l'algorithme Velocity-Verlet
!######################################################################################################


module VelocityVerletPropagatorModule
    use Point2DModule
    implicit none

    type, public :: VelocityVerletPropagator
        private
        real(kind=8) :: box_size
        real(kind=8) :: half_box_size
        real(kind=8) :: delta_t
    contains
        procedure, public :: initialize, evolve_position, evolve_velocity
    end type VelocityVerletPropagator

contains
    subroutine initialize(this, box_size, delta_t)
        class(VelocityVerletPropagator), intent(inout) :: this
        real(kind=8), intent(in):: box_size, delta_t

        this%box_size = box_size
        this%half_box_size = box_size / 2.0d0
        this%delta_t = delta_t
    end subroutine

    function apply_position_boundary_conditions(this, position_i)
        class(VelocityVerletPropagator), intent(in) :: this
        real(kind=8), intent(in):: position_i
        real(kind=8) :: apply_position_boundary_conditions

        if(position_i.gt.(this%half_box_size)) then
            apply_position_boundary_conditions = position_i - this%box_size
        else if (position_i.lt.(-this%half_box_size)) then
            apply_position_boundary_conditions = position_i + this%box_size
        else
            apply_position_boundary_conditions = position_i
        end if
    end function

    function evolve_position_1_dimension(this, position_i, velocity_i, force_i)
        class(VelocityVerletPropagator), intent(in) :: this
        real(kind=8), intent(in):: position_i, velocity_i, force_i
        real(kind=8) :: evolve_position_1_dimension

        evolve_position_1_dimension = position_i + velocity_i * this%delta_t + 0.5d0 * force_i * (this%delta_t **2)
        evolve_position_1_dimension = apply_position_boundary_conditions(this, evolve_position_1_dimension)
    end function

    function evolve_position(this, position_i, velocity_i, force_i)
        class(VelocityVerletPropagator), intent(in) :: this
        type(Point2D), intent(in):: position_i, velocity_i, force_i
        real(kind=8) ::  evolve_x, evolve_y
        type(Point2D) :: evolve_position

        evolve_x = evolve_position_1_dimension(this, position_i%x, velocity_i%x, force_i%x)
        evolve_y = evolve_position_1_dimension(this, position_i%y, velocity_i%y, force_i%y)
        evolve_position = Point2D(evolve_x, evolve_y)
    end function

    function evolve_velocity_1_dimension(this, velocity_i, force_i, force_i_plus)
        class(VelocityVerletPropagator), intent(in) :: this
        real(kind=8), intent(in):: velocity_i, force_i, force_i_plus
        real(kind=8) :: evolve_velocity_1_dimension

        evolve_velocity_1_dimension = velocity_i + 0.5d0 * (force_i + force_i_plus ) * this%delta_t
    end function

    function evolve_velocity(this, velocity_i, force_i, force_i_plus_1)
        class(VelocityVerletPropagator), intent(in) :: this
        type(Point2D), intent(in):: velocity_i, force_i, force_i_plus_1
        real(kind=8) ::  evolve_x, evolve_y
        type(Point2D) :: evolve_velocity

        evolve_x = evolve_velocity_1_dimension(this, velocity_i%x, force_i%x, force_i_plus_1%x)
        evolve_y = evolve_velocity_1_dimension(this, velocity_i%y, force_i%y, force_i_plus_1%y)
        evolve_velocity = Point2D(evolve_x, evolve_y)
    end function

end module VelocityVerletPropagatorModule




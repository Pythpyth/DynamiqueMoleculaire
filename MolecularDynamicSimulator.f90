!######################################################################################################
! MolecularDynamicSimulatorModule
! type MolecularDynamicSimulator
! membres :
!       real(kind=8) :: box_size        : la taille de la boite
!       real(kind=8) :: r_truncated     : le rayon de coupure
!
! subroutine :
!       initialize : initialiser les membres
!
! fonctions :
!
!   run : effectue la simulation en bouclant sur chaque pas de temps
!         retourne pour le moment la liste des énergies totales à chaque temps, à regarder pour
!         voir ce dont on aura besoin
!
!
!######################################################################################################


module MolecularDynamicSimulatorModule
    use Point2DModule
    Use EnergyCalculatorModule
    use ForceCalculatorModule, only : ForceCalculator
    use VelocityVerletPropagatorModule, only : VelocityVerletPropagator
    implicit none

    type, public :: MolecularDynamicSimulator
        private
        real(kind=8) :: box_size
        real(kind=8) :: r_truncated
    contains
        procedure, public :: initialize, run
    end type MolecularDynamicSimulator

contains

    subroutine initialize(this, box_size, r_truncated)
        class(MolecularDynamicSimulator), intent(inout) :: this
        real(kind=8), intent(in) :: box_size, r_truncated
        this%box_size = box_size
        this%r_truncated = r_truncated
    end subroutine


    function run(this, initial_positions, initial_velocities, nb_particle, delta_t, nb_time_step)
        class(MolecularDynamicSimulator), intent(in) :: this
        integer(kind=4), intent(in) :: nb_particle,  nb_time_step
        type(Point2D), dimension(nb_particle), intent(in):: initial_positions, initial_velocities
        real(kind=8), intent(in) :: delta_t
        real(kind=8), dimension(nb_time_step + 1) :: run
        real(kind=8), dimension(nb_time_step + 1) :: kinetic_energy, potential_energy, total_energy
        integer(kind=4) i, j

        type(Point2D), dimension(nb_particle):: positions_i, velocities_i
        type(Point2D), dimension(nb_particle) :: positions_i_plus_1, velocities_i_plus_1
        type(Point2D), dimension(nb_particle) :: forces_i
        type(Point2D) :: forces_i_plus_1
        type(VelocityVerletPropagator) :: propagator
        type(ForceCalculator) :: force_calculator
        real(kind=8) :: half_box_size, potential_shift

        positions_i = initial_positions
        velocities_i = initial_velocities
        half_box_size = this%box_size / 2.0d0

        !calcul du shift pour assurer la continuité du potentiel en r_truncated
        potential_shift = 4.0d0 * ( (this%r_truncated ** (-12)) - (this%r_truncated ** (-6)) )

        call propagator%initialize(this%box_size, delta_t)
        call force_calculator%initialize(this%box_size, this%r_truncated)

        kinetic_energy(1) = compute_sum_kinetic_energy(velocities_i, nb_particle)
        potential_energy(1) = compute_sum_potential_energy(positions_i, nb_particle, this%box_size, half_box_size, &
                                                           this%r_truncated, potential_shift)
        total_energy(1) = kinetic_energy(1) + potential_energy(1)

        do i = 1, nb_time_step

            do j = 1, nb_particle
                forces_i(j) = force_calculator%compute_sum_forces_on_reference_particle(positions_i(j), j, positions_i , &
                                                                                        nb_particle)

                positions_i_plus_1(j) = propagator%evolve_position(positions_i(j), velocities_i(j), forces_i(j))

            end do

             do j = 1, nb_particle
                forces_i_plus_1 = force_calculator%compute_sum_forces_on_reference_particle(positions_i_plus_1(j), &
                                                                            j, positions_i_plus_1 , nb_particle)

                velocities_i_plus_1(j) = propagator%evolve_velocity(velocities_i(j), forces_i(j), forces_i_plus_1)

            end do

            kinetic_energy(i+1) = compute_sum_kinetic_energy(velocities_i_plus_1, nb_particle)

            potential_energy(i+1) = compute_sum_potential_energy(positions_i_plus_1, nb_particle, this%box_size, &
                                                                 half_box_size, this%r_truncated, potential_shift )

            total_energy(i+1) = kinetic_energy(i+1) + potential_energy(i+1)

            positions_i =  positions_i_plus_1
            velocities_i = velocities_i_plus_1

        end do

        print*, total_energy(nb_time_step+1) - total_energy(1)
        print*, kinetic_energy(nb_time_step+1) - kinetic_energy(1)

        run = total_energy

    end function

end module MolecularDynamicSimulatorModule




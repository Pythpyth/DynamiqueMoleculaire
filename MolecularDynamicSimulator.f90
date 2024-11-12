!######################################################################################################
! MolecularDynamicSimulatorModule
! type MolecularDynamicSimulator
! membres :
!       real(kind=8) :: box_size        : la taille de la boite
!       real(kind=8) :: r_truncated     : le rayon de coupure
!       real(kind=8) :: potential_shift : le shift du potentiel
!
! subroutine :
!       initialize : initialiser les membres
!
! fonctions :
!
!   run : effectue la simulation en bouclant sur chaque pas de temps
!         retourne un type(simulationresult) qui contient les observales de notre simulation
!
!
!######################################################################################################


module MolecularDynamicSimulatorModule
    use Point2DModule
    Use EnergyCalculatorModule
    use ForceCalculatorModule, only : ForceCalculator
    use VelocityVerletPropagatorModule, only : VelocityVerletPropagator
    use SimulationResultModule
    use GridPrinterModule
    implicit none

    type, public :: MolecularDynamicSimulator
        private
        real(kind=8) :: box_size
        real(kind=8) :: r_truncated
        real(kind=8) :: potential_shift
    contains
        procedure, public :: initialize, run
    end type MolecularDynamicSimulator

contains

    subroutine initialize(this, box_size, r_truncated, potential_shift)
        class(MolecularDynamicSimulator), intent(inout) :: this
        real(kind=8), intent(in) :: box_size, r_truncated, potential_shift
        this%box_size = box_size
        this%r_truncated = r_truncated
        this%potential_shift = potential_shift
    end subroutine

    function run(this, initial_positions, initial_velocities, nb_particle, delta_t, nb_time_step)
        class(MolecularDynamicSimulator), intent(in) :: this
        integer(kind=4), intent(in) :: nb_particle,  nb_time_step
        type(Point2D), dimension(nb_particle), intent(in):: initial_positions, initial_velocities
        real(kind=8), intent(in) :: delta_t
        type(SimulationResult) :: run
        real(kind=8), dimension(nb_time_step + 1) :: kinetic_energy, potential_energy,&
                                                     total_energy, temperature, time, pressure
        integer(kind=4) i, j, r

        type(Point2D), dimension(nb_particle):: positions_i, velocities_i
        type(Point2D), dimension(nb_particle) :: positions_i_plus_1, velocities_i_plus_1
        type(Point2D), dimension(nb_particle) :: forces_i
        type(Point2D) :: forces_i_plus_1
        type(VelocityVerletPropagator) :: propagator
        type(ForceCalculator) :: force_calculator
        real(kind=8) :: half_box_size
        character(len=100), allocatable :: file_name_pos
        logical :: is_all_particle_in_box_boolean
        integer, parameter :: distribution_size = 10
        real(kind=8), dimension(distribution_size) :: radial_distribution

        positions_i = initial_positions
        velocities_i = initial_velocities
        half_box_size = this%box_size / 2.0d0

        call propagator%initialize(this%box_size, delta_t)
        call force_calculator%initialize(this%box_size, this%r_truncated)


        time(1) = 0.0d0
        kinetic_energy(1) = compute_sum_kinetic_energy(velocities_i, nb_particle)
        potential_energy(1) = compute_sum_potential_energy(positions_i, nb_particle, this%box_size, half_box_size, &
                                                            this%r_truncated, this%potential_shift)

        total_energy(1) = kinetic_energy(1) + potential_energy(1)
        temperature(1) = compute_temperature(kinetic_energy(1), nb_particle)
        pressure(1) = compute_pressure(positions_i, nb_particle, this%box_size, temperature(1), force_calculator)

        file_name_pos  = ''

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

            time(i+1) = delta_t * real(i,8)
            kinetic_energy(i+1) = compute_sum_kinetic_energy(velocities_i_plus_1, nb_particle)

            potential_energy(i+1) = compute_sum_potential_energy(positions_i_plus_1, nb_particle, this%box_size, half_box_size, &
                                                                   this%r_truncated, this%potential_shift)

            total_energy(i+1) = kinetic_energy(i+1) + potential_energy(i+1)
            temperature(i+1) = compute_temperature(kinetic_energy(i+1), nb_particle)
            pressure(i+1) = compute_pressure(positions_i_plus_1, nb_particle, this%box_size, temperature(i+1), &
                                             force_calculator)

            positions_i =  positions_i_plus_1
            velocities_i = velocities_i_plus_1

            if (mod(i,1000) ==1) then
                write(file_name_pos, '(A,I7.7,A)') '..\position\\Position_dt', i,'.dat'
                open(i*113, file = trim(file_name_pos), status = 'replace')
                do r=1, size(positions_i)
                    write(i*113,*) positions_i(r)
                end do

                close(i*113)
                print*, i

                is_all_particle_in_box_boolean = is_all_particle_in_box(positions_i_plus_1,nb_particle, half_box_size)

                if (.not. is_all_particle_in_box_boolean) then
                    print* , 'a particule is out of the box'
                end if

                radial_distribution = &
                compute_radial_distribution(positions_i, nb_particle, half_box_size, this%box_size, distribution_size)

                write(file_name_pos, '(A,I7.7,A)') '..\Distribution\\distribution_dt', i,'.dat'
                open(i*112, file = trim(file_name_pos), status = 'replace')

                do r=1, distribution_size
                    write(i*112,*) radial_distribution(r)
                end do

                close(i*112)

            end if


        end do

        run%nb_time = nb_time_step + 1
        run%time = time
        run%potential_energy = potential_energy
        run%kinetic_energy = kinetic_energy
        run%total_energy = total_energy
        run%temperature = temperature
        run%pressure = pressure
    end function

end module MolecularDynamicSimulatorModule




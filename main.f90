!######################################################################################################
! programe MolecularDynamicSimulation
!
! On considère 121 particules dans une boite 11 x 11 :
! on a 11 particules par dimension. La distance dans chaque dimension entre chaque particule vaut 1
! Les particules les plus proches du bord de la boite sont à une distance de 0.5 du bord de la boite
! dans chaque dimension. On a ainsi une distance de 1 dans chaque dimension avec l'image dans l'image
! périodique d'une boite en contact.
! On initialise les vitesses selon une loi normale centrée et d'écart-type 0.01 à regarder
!
!######################################################################################################

program MolecularDynamicSimulation
    use MolecularDynamicSimulatorModule
    use Point2DModule
    use GridPrinterModule
    use UniformRandomGeneratorModule
    use UnitConverterModule
    implicit none

    type(MolecularDynamicSimulator) molecular_dynamic_simulator
    ! t_star : temps utilisé pour tracer la ln|E(t*)-Et(0)| = f(dt)
    ! time_end : temps de la fin pour la simulation complète
    real(kind=8) :: t_star, time_end, delta_t
    real(kind=8), parameter :: r_truncated = 2.5d0
    real(kind=8), parameter :: box_size = 11.0d0
    integer(kind=4), parameter :: nb_particle_per_dim = 11
    integer(kind=4), parameter  :: nb_particle = 121
    real(kind=8), parameter:: pi = 3.1415926535897932d0

    integer(kind=4), parameter :: nb_time_step_simulation = 100000
    integer(kind=4), parameter :: nb_time_step_for_validation = 5
    integer(kind=4), dimension(nb_time_step_for_validation) ::  &
    nb_time_steps = (/5000, 10000, 20000, 40000, 80000/)

    real(kind=8), dimension(nb_time_step_for_validation) :: ln_delta_t, ln_energy_diff
    type(Point2D), dimension(nb_particle) :: initial_positions, initial_velocities
    character(len=:), allocatable :: result_file
    integer(kind=4) i, j
    real(kind=8) :: u1, u2

    !température du système fixée à t=0
    real(kind=8), parameter :: initial_temperature_K = 298.15
    real(kind=8), parameter :: particle_mass_atomic_unit = 21874.6618353931 !masse du carbone en unité atomique
    !paramètres du potentiel de Lennard-Jones en unité atomique (sigma_LJ =3.5 Å et epsilon_LJ : 0.00286 eV)
    real(kind=8) :: sigma_LJ_atomic_unit = 6.61404144067724, epsilon_LJ_atomic_unit = 0.0001051030614224d0
    !Boltzmann constant atomic unit by K
    real(kind=8), parameter :: kb_atomic_unit_by_K = 3.1668115634564d-6



    ! paramètres du générateur implémenté en TP avec la plus grande période
    integer(kind=8), parameter :: k = 7 **5, l = 0, m= 2**31 -1, x_0= 12
    type(UniformRandomGenerator) :: uniform_random_generator

    real(kind=8) :: kinetic_energy_random, temperature_random, &
                    temperature_random_K, adjustment_velocity

    type(SimulationResult) :: simulation_result
    real(kind=8), dimension(:), allocatable :: kinetic_energy, potential_energy, temperature, &
                                                total_energy, time_list


    ! initialisation du générateur aléatoire uniforme
    call uniform_random_generator%initialize(k, l, m, x_0)


    kinetic_energy_random = 0
    !initialisation des conditions initiales
    do i = 1, nb_particle_per_dim
        do j = 1, nb_particle_per_dim
            initial_positions(j + (i-1) * nb_particle_per_dim)%x = -6.0d0 + real(i,8) * 1.0d0
            initial_positions(j + (i-1) * nb_particle_per_dim)%y = -6.0d0 + real(j,8) * 1.0d0

            u1 =  uniform_random_generator%generate()
            u2 =  uniform_random_generator%generate()

            !A partir de deux nombres aléatoires u1 et u2 distribués de manière uniforme sur l’intervalle [0, 1[,
            !l'algorithme de Box-Müller fournit deux nombres g1, g2 distribué selon une loi normale centrée réduite

            initial_velocities(j + (i-1) *nb_particle_per_dim)%x = sqrt(- 2.0d0 * log(u1)) * cos(2.0d0 * pi * u2)
            initial_velocities(j + (i-1) *nb_particle_per_dim)%y = sqrt(- 2.0d0 * log(u1)) * sin(2.0d0 * pi * u2)
        end do
    end do

    !On calcule l'énergie cinétique te la température en unité réduite
    kinetic_energy_random = compute_sum_kinetic_energy(initial_velocities, nb_particle)
    temperature_random = compute_temperature(kinetic_energy_random, nb_particle)

    !on convertit en K
    temperature_random_K = convert_reduced_unit_to_kelvin(temperature_random, epsilon_LJ_atomic_unit, kb_atomic_unit_by_K)

    !on ajuste les vitesses pour obtenir la température initiale
    adjustment_velocity = sqrt(initial_temperature_K /temperature_random_K)

    do i = 1, nb_particle_per_dim
        do j = 1, nb_particle_per_dim

            initial_velocities(j + (i-1) *nb_particle_per_dim)%x = &
            initial_velocities(j + (i-1) *nb_particle_per_dim)%x * adjustment_velocity
            initial_velocities(j + (i-1) *nb_particle_per_dim)%y = &
            initial_velocities(j + (i-1) *nb_particle_per_dim)%y  * adjustment_velocity

        end do
    end do



    call molecular_dynamic_simulator%initialize(box_size, r_truncated)

    !Partie pour valider la relation ln|E(t*)-E(t0)| = a * ln(delta_t) + b
    !Avec Velocity-Verlet, on devrait obtenir a = 2
    t_star = 0.02d0

    !Boucle pour différents pas de temps
    do i = 1, nb_time_step_for_validation

        delta_t = t_star / real(nb_time_steps(i), 8)

        simulation_result = molecular_dynamic_simulator%run(initial_positions, initial_velocities, nb_particle, &
                                                            delta_t, nb_time_steps(i))


        total_energy = convert_reduced_unit_to_energy_unit(simulation_result%total_energy, nb_time_steps(i)+1, &
                                                            epsilon_LJ_atomic_unit)

        ln_delta_t(i) = log(delta_t)
        ln_energy_diff(i) = log(abs(simulation_result%total_energy(nb_time_steps(i)+ 1) - simulation_result%total_energy(1)))

        deallocate(total_energy)
        call simulation_result%terminate()
    end do

    result_file = "..\\Results\ln_energy_diff_f_ln_delta_t.dat"
    call print_energy(result_file, ln_delta_t, ln_energy_diff, nb_time_step_for_validation)


    time_end = 20.0d0
    delta_t = time_end / real(nb_time_step_simulation, 8)

    simulation_result = molecular_dynamic_simulator%run(initial_positions, initial_velocities, nb_particle, &
                                                        delta_t, nb_time_step_simulation)


    kinetic_energy = convert_reduced_unit_to_energy_unit(simulation_result%kinetic_energy, nb_time_step_simulation +1, &
                                                             epsilon_LJ_atomic_unit)
    potential_energy = convert_reduced_unit_to_energy_unit(simulation_result%potential_energy, nb_time_step_simulation+1, &
                                                             epsilon_LJ_atomic_unit)

    total_energy = convert_reduced_unit_to_energy_unit(simulation_result%total_energy, nb_time_step_simulation+1, &
                                                            epsilon_LJ_atomic_unit)

    temperature = convert_list_reduced_unit_to_kelvin(simulation_result%temperature, nb_time_step_simulation+1, &
                                                          epsilon_LJ_atomic_unit, kb_atomic_unit_by_K)

    time_list = convert_list_reduced_unit_to_second(simulation_result%time, nb_time_step_simulation+1, &
                                                    particle_mass_atomic_unit, &
                                                    epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)

    print* ,temperature(1)
    print* ,temperature(nb_time_step_simulation + 1)

    result_file = "..\\Results\simulation_result.dat"
    call print_simulation_result(result_file, time_list, kinetic_energy, &
                                 potential_energy, total_energy, temperature, nb_time_step_simulation)


    deallocate(kinetic_energy)
    deallocate(potential_energy)
    deallocate(total_energy)
    deallocate(temperature)
    deallocate(time_list)

end program


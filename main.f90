!######################################################################################################
! programe MolecularDynamicSimulation
!
! On considère 121 particules dans une boite 11 x 11 :
! on a 11 particules par dimension. La distance dans chaque dimension entre chaque particule vaut 1
! Les particules les plus proches du bord de la boite sont à une distance de 0.5 du bord de la boite
! dans chaque dimension. On a ainsi une distance de 1 dans chaque dimension avec l'image dans l'image
! périodique d'une boite en contact.
! On initialise les vitesses selon une loi normale centrée réduite puis ajustée pour obtenir la température T0
!
! Il faut créer un dossier '..\\Results', '..\Position', '..\Distribution'
!######################################################################################################

program MolecularDynamicSimulation
    use MolecularDynamicSimulatorModule
    use EnergyCalculatorModule
    use Point2DModule
    use GridPrinterModule
    use UniformRandomGeneratorModule
    use UnitConverterModule
    Use ObservableCalculatorModule
    implicit none

    type(MolecularDynamicSimulator) molecular_dynamic_simulator
    ! t_star : temps utilisé pour tracer la ln|E(t*)-Et(0)| = f(dt)
    ! time_end : temps de la fin pour la simulation complète
    ! t_ref : temps à partir duquel on calcul la moyenne temporelle
    real(kind=8) :: t_star, t_ref, time_end, delta_t
    real(kind=8), parameter :: r_truncated = 2.5d0
    real(kind=8), parameter :: box_size = 10.0d0
    integer(kind=4), parameter :: nb_particle_per_dim = 10
    integer(kind=4), parameter  :: nb_particle = 100

    integer(kind=4), parameter :: nb_time_step_simulation = 20000
    integer(kind=4), parameter :: nb_time_step_for_validation = 5
    integer(kind=4), dimension(nb_time_step_for_validation) ::  &
    nb_time_steps = (/50, 100, 200, 400, 800/)

    real(kind=8), dimension(nb_time_step_for_validation) :: ln_delta_t, ln_energy_diff
    type(Point2D), dimension(nb_particle) :: initial_positions, initial_velocities
    character(len=:), allocatable :: result_file
    integer(kind=4) i, j
    real(kind=8) :: u1, u2

    !température du système fixée à t=0
    real(kind=8), parameter :: initial_temperature_K = 298.15
    real(kind=8), parameter :: particle_mass_atomic_unit = 21874.6618353931 !masse du carbone en unité atomique
    !paramètres du potentiel de Lennard-Jones en unité atomique (sigma_LJ =3.5 Å et epsilon_LJ : 0.00286 eV)
    real(kind=8) :: sigma_LJ_atomic_unit = 6.61404144067724, epsilon_LJ_atomic_unit = 0.0001051030614224d0, &
                    epsilon_LJ_eV = 0.00286d0 !pour le calcul des énergies
    !Boltzmann constant atomic unit by K
    real(kind=8), parameter :: kb_atomic_unit_by_K = 3.1668115634564d-6
    real(kind=8) :: potential_shift, potential_energy_tail_correction, potential_energy_shift_correction, &
                    pressure_tail_correction

    real(kind=8) :: mean_kinetic_energy, mean_potential_energy, mean_total_energy, mean_temperature, mean_pressure
    real(kind=8) :: stdev_kinetic_energy, stdev_potential_energy, stdev_total_energy, stdev_temperature, stdev_pressure

    ! paramètres du générateur implémenté en TP avec la plus grande période
    integer, parameter :: int64 = selected_int_kind(18)
    integer(kind=8), parameter :: k = 7 ** 5, l = 0, m = int(2,int64)**31 - 1, x_0 = 12

    type(UniformRandomGenerator) :: uniform_random_generator

    real(kind=8) :: kinetic_energy_random, temperature_random, &
                    temperature_random_K, adjustment_velocity

    type(SimulationResult) :: simulation_result
    real(kind=8), dimension(:), allocatable :: kinetic_energy, potential_energy, temperature, &
                                                total_energy, time_list, pressure
    integer :: number_observable_from_t_ref
    logical :: is_validation_step

    ! initialisation du générateur aléatoire uniforme
    call uniform_random_generator%initialize(k, l, m, x_0)

    !kinetic_energy_random = 0.0d0
    !initialisation des conditions initiales
    do i = 1, nb_particle_per_dim
        do j = 1, nb_particle_per_dim
            initial_positions(j + (i-1) * nb_particle_per_dim)%x = -5.5d0 + real(i,8) * 1.0d0
            initial_positions(j + (i-1) * nb_particle_per_dim)%y = -5.5d0 + real(j,8) * 1.0d0

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

    !calcul du shift pour assurer la continuité du potentiel en r_truncated
    potential_shift = 4.0d0 * ( (r_truncated ** (-12)) - (r_truncated ** (-6)) )
    call molecular_dynamic_simulator%initialize(box_size, r_truncated, potential_shift)


    is_validation_step = .true.

    if (is_validation_step) then
        !Partie pour valider la relation ln|E(t*)-E(t0)| = a * ln(delta_t) + b
        !Avec Velocity-Verlet, on devrait obtenir a = 2
        !premier max de l'énergie totale à 70fs
        t_star = 70.0d0 / convert_reduced_unit_to_fs(1.0d0, particle_mass_atomic_unit, epsilon_LJ_atomic_unit, &
                                                 sigma_LJ_atomic_unit)

        print* , 't* : ', t_star
        !Boucle pour différents pas de temps
        do i = 1, nb_time_step_for_validation

            delta_t = t_star / real(nb_time_steps(i), 8)

            simulation_result = molecular_dynamic_simulator%run(initial_positions, initial_velocities, nb_particle, &
                                                            delta_t, nb_time_steps(i))

            ln_delta_t(i) = log(delta_t)
            ln_energy_diff(i) = log(abs(simulation_result%total_energy(nb_time_steps(i)+ 1) - simulation_result%total_energy(1)))


            call simulation_result%terminate()
        end do

        result_file = "..\\Results\ln_energy_diff_f_ln_delta_t.dat"
        call print_energy(result_file, ln_delta_t, ln_energy_diff, nb_time_step_for_validation)

    end if


    time_end = 4.0d0
    delta_t = time_end / real(nb_time_step_simulation, 8)


    simulation_result = molecular_dynamic_simulator%run(initial_positions, initial_velocities, nb_particle, &
                                                        delta_t, nb_time_step_simulation)


    kinetic_energy = convert_reduced_unit_to_energy_unit(simulation_result%kinetic_energy, nb_time_step_simulation +1, &
                                                             epsilon_LJ_eV)


    potential_energy = convert_reduced_unit_to_energy_unit(simulation_result%potential_energy, nb_time_step_simulation+1, &
                                                             epsilon_LJ_eV)

    total_energy = convert_reduced_unit_to_energy_unit(simulation_result%total_energy, nb_time_step_simulation+1, &
                                                            epsilon_LJ_eV)

    temperature = convert_list_reduced_unit_to_kelvin(simulation_result%temperature, nb_time_step_simulation+1, &
                                                          epsilon_LJ_atomic_unit, kb_atomic_unit_by_K)

    time_list = convert_list_reduced_unit_to_fs(simulation_result%time, nb_time_step_simulation+1, &
                                                    particle_mass_atomic_unit, &
                                                    epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)

    pressure = convert_list_reduced_unit_to_newton_by_meter(simulation_result%pressure, nb_time_step_simulation+1 , &
                                                            epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)


    !correction de l'énergie potentielle de queue
    potential_energy_tail_correction = compute_potential_energy_tail_correction(nb_particle, box_size, &
                                                                                r_truncated, epsilon_LJ_eV)

    !correction liée au shift de l'énergie potentielle
    potential_energy_shift_correction = compute_potential_energy_shift_correction(nb_particle, box_size, r_truncated, &
                                                       potential_shift, epsilon_LJ_eV)


    !correction liée à la queue de pression
    pressure_tail_correction = compute_pressure_tail_correction(nb_particle, box_size, r_truncated, &
                                                                epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)

    print*,'potential energy tail correction ' , potential_energy_tail_correction, ' ev'
    print*, 'potential energy shift correction ' , potential_energy_shift_correction, ' ev'
    print*, 'pressure tail correction', pressure_tail_correction, ' N/m'

    do i = 1, nb_time_step_simulation + 1
        potential_energy(i) = potential_energy(i) + potential_energy_tail_correction + potential_energy_shift_correction
        total_energy(i) = total_energy(i) + potential_energy_tail_correction + potential_energy_shift_correction
        pressure(i) = pressure(i) + pressure_tail_correction
    end do


    t_ref = 500.0d0 !en fs

    mean_kinetic_energy = compute_mean_observable_from_t_ref(t_ref, time_list, kinetic_energy, nb_time_step_simulation + 1)
    mean_potential_energy = compute_mean_observable_from_t_ref(t_ref, time_list, potential_energy, nb_time_step_simulation + 1)
    mean_total_energy = compute_mean_observable_from_t_ref(t_ref, time_list, total_energy, nb_time_step_simulation + 1)



    mean_temperature = compute_mean_observable_from_t_ref(t_ref, time_list, temperature, nb_time_step_simulation + 1)
    mean_pressure = compute_mean_observable_from_t_ref(t_ref, time_list, pressure, nb_time_step_simulation + 1)


    stdev_kinetic_energy = compute_stdev_observable_from_t_ref(t_ref, time_list, kinetic_energy, nb_time_step_simulation + 1)
    stdev_potential_energy = compute_stdev_observable_from_t_ref(t_ref, time_list, potential_energy, nb_time_step_simulation + 1)
    stdev_total_energy = compute_stdev_observable_from_t_ref(t_ref, time_list, total_energy, nb_time_step_simulation + 1)



    stdev_temperature = compute_stdev_observable_from_t_ref(t_ref, time_list, temperature, nb_time_step_simulation + 1)
    stdev_pressure = compute_stdev_observable_from_t_ref(t_ref, time_list, pressure, nb_time_step_simulation + 1)



    print*, 'mean kinetic energy : ' , mean_kinetic_energy, ' eV'
    print*, 'mean potential energy : ' , mean_potential_energy, ' eV'
    print*, 'mean total energy : ' , mean_total_energy, ' eV'
    print*, 'mean temperature : ' , mean_temperature, ' K'
    print*, 'mean pressure : ' , mean_pressure, ' N/m'

    print*, 'stdev kinetic energy : ' , stdev_kinetic_energy
    print*, 'stdev potential energy : ' , stdev_potential_energy
    print*, 'stdev total energy : ' , stdev_total_energy
    print*, 'stdev pressure : ' , stdev_pressure
    print*, 'stdev temperature : ' , stdev_temperature


    number_observable_from_t_ref = compute_number_observable_from_t_ref(t_ref, time_list, &
                                                                        nb_time_step_simulation + 1)

    print*, 'number_observable_from_t_ref  : ' , number_observable_from_t_ref

    !les écarts sont calculés sans prendre en compte l'autocorrélation
    !3 écart-types => intervalle de confiance à 99.7%
    print*, 'kinetic energy confidence interval  ' , mean_kinetic_energy - &
                                                3.0d0 * stdev_kinetic_energy / sqrt(real(number_observable_from_t_ref-1,8)) &
                                             , mean_kinetic_energy &
                                             + 3.0d0 * stdev_kinetic_energy / sqrt(real(number_observable_from_t_ref-1,8))

    print*, 'potential energy confidence interval  ' , mean_potential_energy - &
                                                3.0d0 * stdev_potential_energy / sqrt(real(number_observable_from_t_ref-1,8)) &
                                             , mean_potential_energy &
                                             + 3.0d0 * stdev_potential_energy / sqrt(real(number_observable_from_t_ref-1,8))


    print*, 'total energy confidence interval  ' , mean_total_energy - &
                                                3.0d0 * stdev_total_energy / sqrt(real(number_observable_from_t_ref-1,8)) &
                                             , mean_total_energy &
                                             + 3.0d0 * stdev_total_energy / sqrt(real(number_observable_from_t_ref-1,8))


    print*, 'pressure confidence interval  ' , mean_pressure - &
                                                3.0d0 * stdev_pressure / sqrt(real(number_observable_from_t_ref-1,8)) &
                                             , mean_pressure &
                                             + 3.0d0 * stdev_pressure / sqrt(real(number_observable_from_t_ref-1,8))


    print*, 'temperature confidence interval  ' , mean_temperature - &
                                                3.0d0 * stdev_temperature / sqrt(real(number_observable_from_t_ref-1,8)) &
                                             , mean_temperature &
                                             + 3.0d0 * stdev_temperature / sqrt(real(number_observable_from_t_ref-1,8))

    print* ,'temperature t0 ', temperature(1)
    print* ,'temperature tend ', temperature(nb_time_step_simulation + 1)
    print* ,'pressure t0 ', pressure(1)
    print* ,'pressure tend ',pressure(nb_time_step_simulation + 1)

    result_file = "..\\Results\simulation_result.dat"
    call print_simulation_result(result_file, time_list, kinetic_energy, &
                                 potential_energy, total_energy, temperature, pressure, nb_time_step_simulation)


    deallocate(kinetic_energy)
    deallocate(potential_energy)
    deallocate(total_energy)
    deallocate(temperature)
    deallocate(time_list)
    deallocate(pressure)

end program


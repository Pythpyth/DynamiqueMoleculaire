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
    implicit none

    type(MolecularDynamicSimulator) molecular_dynamic_simulator
    real(kind=8) :: t_star, delta_t
    real(kind=8), dimension(:), allocatable :: energy_list
    real(kind=8), parameter :: r_truncated = 2.5d0
    real(kind=8), parameter :: box_size = 11.0d0
    integer(kind=4), parameter :: nb_particle_per_dim = 11
    integer(kind=4), parameter  :: nb_particle = 121
    real(kind=8), parameter:: pi = 3.1415926535897932d0

    integer(kind=4), parameter :: nb_time_step_for_validation = 5
    integer(kind=4), dimension(nb_time_step_for_validation) ::  &
    nb_time_steps = (/5000, 10000, 20000, 40000, 80000/)


    real(kind=8), dimension(nb_time_step_for_validation) :: ln_delta_t, ln_energy_diff
    type(Point2D), dimension(nb_particle) :: initial_positions, initial_velocities
    character(len=:), allocatable :: energy_file
    integer(kind=4) i, j
    real(kind=8) :: u1, u2

    ! paramètres du générateur implémenté en TP avec la plus grande période
    integer(kind=8), parameter :: k = 7 **5, l = 0, m= 2**31 -1, x_0= 12
    type(UniformRandomGenerator) :: uniform_random_generator

    ! initialisation du générateur aléatoire uniforme
    call uniform_random_generator%initialize(k, l, m, x_0)

    !initialisation des conditions initiales
    do i = 1, nb_particle_per_dim
        do j = 1, nb_particle_per_dim
            initial_positions(j + (i-1) * nb_particle_per_dim)%x = -6.0d0 + real(i,8) * 1.0d0
            initial_positions(j + (i-1) * nb_particle_per_dim)%y = -6.0d0 + real(j,8) * 1.0d0

            u1 =  uniform_random_generator%generate()
            u2 =  uniform_random_generator%generate()

            !A partir de deux nombres aléatoires u1 et u2 distribués de manière uniforme sur l’intervalle [0, 1[,
            !l'algorithme de Box-Müller fournit deux nombres g1, g2 distribué selon une loi normale centrée réduite
            !On obtient une loi normale centrée d'écart-type sigma par g1' = sigma * g1 et g'2 =sigma * g2

            initial_velocities(j + (i-1) *nb_particle_per_dim)%x = 0.01d0 * sqrt(- 2.0d0 * log(u1)) * cos(2.0d0 * pi * u2)
            initial_velocities(j + (i-1) *nb_particle_per_dim)%y = 0.01d0 * sqrt(- 2.0d0 * log(u1)) * sin(2.0d0 * pi * u2)

        end do
    end do

    call molecular_dynamic_simulator%initialize(box_size, r_truncated)

    !Partie pour valider la relation ln|E(t*)-E(t0)| = a * ln(delta_t) + b
    !Avec Velocity-Verlet, on devrait obtenir a = 2
    t_star = 0.5d0

    !Boucle pour différents pas de temps
    do i = 1, nb_time_step_for_validation

        delta_t = t_star / real(nb_time_steps(i), 8)
        allocate(energy_list(nb_time_steps(i)+ 1))

        energy_list = molecular_dynamic_simulator%run(initial_positions, initial_velocities, nb_particle, &
                                                      delta_t, nb_time_steps(i))


        ln_delta_t(i) = log(delta_t)
        ln_energy_diff(i) = log(abs(energy_list(nb_time_steps(i)+ 1) - energy_list(1)))

        deallocate(energy_list)

    end do

    energy_file = "..\\Results\ln_energy_diff_f_ln_delta_t.dat"
    call print_energy(energy_file, ln_delta_t, ln_energy_diff, nb_time_step_for_validation)

end program


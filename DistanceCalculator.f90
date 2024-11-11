!######################################################################################################
! Module DistanceCalculatorModule
! Fonctions :
!   compute_distance : calcule la distance entre 2 particules
!
!   compute_closest_image_position_2 : on considère une position de référence et une position 2.
!                                      la fonction retourne l’image périodique de la position 2 la plus
!                                      proche de la position référence
!
!    is_all_particle_in_box : vérifie que toutes les particules sont dans la boite
!
!    compute_radial_distribution : calcul de la distribution radiale
!
!######################################################################################################


module DistanceCalculatorModule
    use Point2DModule
    implicit none

contains
    function compute_distance(position_1, position_2)
        type(Point2D), intent(in):: position_1, position_2
        real(kind=8) :: compute_distance
        compute_distance = sqrt( (position_1%x - position_2%x)**2 + (position_1%y - position_2%y)**2 )
    end function

    function compute_closest_image_position_2(position_reference, position_2, half_box_size, box_size)
        type(Point2D),  intent(in):: position_reference, position_2
        real(kind=8) ::  half_box_size, box_size
        type(Point2D) :: compute_closest_image_position_2
        real(kind=8) :: distance_x, distance_y

        distance_x = position_reference%x - position_2%x
        distance_y = position_reference%y - position_2%y

        compute_closest_image_position_2 = position_2

        if (distance_x.gt.(half_box_size)) then
            compute_closest_image_position_2%x = compute_closest_image_position_2%x + box_size
        else if (distance_x.lt.(-half_box_size)) then
             compute_closest_image_position_2%x = compute_closest_image_position_2%x - box_size
        end if

        if (distance_y.gt.(half_box_size)) then
            compute_closest_image_position_2%y = compute_closest_image_position_2%y + box_size
        elseif (distance_y.lt.(-half_box_size)) then
            compute_closest_image_position_2%y = compute_closest_image_position_2%y - box_size
        end if

    end function


    function is_all_particle_in_box(positions, nb_particle, half_box_size)
        integer(kind=4), intent(in) :: nb_particle
        type(Point2D), dimension(nb_particle) ::  positions
        real(kind=8), intent(in):: half_box_size
        logical :: is_all_particle_in_box
        integer :: i
        is_all_particle_in_box = .true.

        do i =1, nb_particle

            if( abs(positions(i)%x) > half_box_size ) then
                is_all_particle_in_box = .false.
                exit
            end if

            if( abs(positions(i)%y) > half_box_size ) then
                is_all_particle_in_box = .false.
                exit
            end if

        end do

    end function


    function compute_radial_distribution(positions, nb_particle, half_box_size, box_size, distribution_size)

        integer(kind=4), intent(in) :: nb_particle
        real(kind=8), intent(in):: box_size, half_box_size
        integer, intent(in) :: distribution_size
        type(Point2D), dimension(nb_particle), intent(in) :: positions
        real(kind=8) :: compute_sum_potential_energy
        real(kind=8), dimension(nb_particle) :: compute_radial_distribution
        real(kind=8) :: delta_r, distance
        integer(kind=4) :: i,j, domaine_index
        type(Point2D) :: closest_image_position_j

        compute_sum_potential_energy = 0.0d0
        delta_r = half_box_size / real(distribution_size, 8)

        do i = 1, distribution_size
            compute_radial_distribution(i) = 0.0d0
        end do

        do i = 1, nb_particle

             do j = i + 1, nb_particle ! Compter la paire (i, j) pas besoin de compter (j,i) car identique

                closest_image_position_j = compute_closest_image_position_2(positions(i), positions(j), &
                                                                             half_box_size, box_size)

                distance = compute_distance(positions(i), closest_image_position_j)

                if (distance < half_box_size) then
                    ! Calculer l'indice du domaine correspondant
                    domaine_index = int(distance / delta_r) + 1
                    compute_radial_distribution(domaine_index) = &
                    compute_radial_distribution(domaine_index) + 1.0d0
                end if

             end do

        end do

    end function

end module DistanceCalculatorModule


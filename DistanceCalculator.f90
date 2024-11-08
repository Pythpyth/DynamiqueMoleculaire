!######################################################################################################
! Module DistanceCalculatorModule
! Fonctions :
!   compute_distance : calcule la distance entre 2 particules
!   compute_closest_image_position_2 : on considère une position de référence et une position 2.
!                                      la fonction retourne l’image périodique de la position 2 la plus
!                                      proche de la position référence
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

end module DistanceCalculatorModule


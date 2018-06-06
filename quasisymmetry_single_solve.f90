subroutine quasisymmetry_single_solve

  use quasisymmetry_variables

  implicit none

  integer :: iteration, new_N_phi

  iota_tolerance_achieved = .false.
  elongation_tolerance_achieved = .false.
  iteration = 0
  N_phi = N_phi_original
  do 
     iteration = iteration + 1
     print "(a)"," -------------------------------------------------------"
     print "(a,i4)"," Solving system using N_phi=",N_phi
     
     call quasisymmetry_init_phi()

     call quasisymmetry_init_axis()

     call quasisymmetry_solve()

     if (trim(resolution_option) == resolution_option_fixed) exit

     if (iteration > 1) then
        print "(a,es10.3)"," abs(iota - last_iota)                     =",abs(iota - last_iota)
        print "(a,es10.3)"," abs(max_elongation - last_max_elongation) =",abs(max_elongation - last_max_elongation)
        if (abs(iota - last_iota) <= iota_tolerance) then
           print *,"iota_tolerance achieved."
           iota_tolerance_achieved = .true.
        end if
        if (abs(max_elongation - last_max_elongation) <= elongation_tolerance) then
           print *,"elongation_tolerance achieved."
           elongation_tolerance_achieved = .true.
        end if
        if (iota_tolerance_achieved .and. elongation_tolerance_achieved) exit
     end if

     last_iota = iota
     last_max_elongation = max_elongation

     new_N_phi = N_phi * 2 + 1
     if (new_N_phi > max_N_phi) then
        print *,"Stopping N_phi refinement since max_N_phi exceeded."
        exit
     end if
     N_phi = new_N_phi

  end do

end subroutine quasisymmetry_single_solve

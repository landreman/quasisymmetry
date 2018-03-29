subroutine quasisymmetry_single_solve

  use quasisymmetry_variables

  implicit none

  integer :: j_N_phi

  do j_N_phi = 1:N_N_phi

     N_phi = N_phis(j_N_phi)
     print *,"Handling N_phi=",N_phi
     
     call quasisymmetry_init_phi()

     call quasisymmetry_init_axis()

     call quasisymmetry_solve()

  end do

end subroutine quasisymmetry_single_solve

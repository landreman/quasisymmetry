subroutine quasisymmetry_Jacobian

  use quasisymmetry_variables

  implicit none

  integer :: j

  Jacobian = 0

  ! d (Riccati equation) / d iota:
  Jacobian(1:N_phi, N_phi+1) = B1Squared_over_curvatureSquared * B1Squared_over_curvatureSquared + 1 + sigma * sigma

  ! d (Riccati equation) / d sigma:
  Jacobian(1:N_phi, 1:N_phi) = d_d_zeta
  do j = 1,N_phi
     Jacobian(j,j) = Jacobian(j,j) + iota * 2 * sigma(j)
  end do

  ! d (Z(0)) / d sigma:
  Jacobian(N_phi+1, 1) = &
       (-RZ_to_XY_b(1) * (0 + sign_G * curvature(1) * (0 + B1c_over_B0 * 1) / (B1c_over_B0*B1c_over_B0 + B1c_over_B0*B1s_over_B0) * binormal_cylindrical(1,1)) &
       + RZ_to_XY_a(1) * (0 + sign_G * curvature(1) * (0 + B1c_over_B0 * 1) / (B1c_over_B0*B1c_over_B0 + B1c_over_B0*B1s_over_B0) * binormal_cylindrical(1,3)))

end subroutine quasisymmetry_Jacobian

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

  ! At zeta=0, we want theta=0 to correspond to the direction e_R.
  ! Therefore we impose Z1c(1) = 0.
  ! We use
  ! Z1c = ( binormal_cylindrical(:,1) .* X1c - normal_cylindrical(:,1) .* Y1c) .* d_l_d_phi ./ R0
  ! and
  ! X1c = B1c_over_B0 ./ curvature
  ! and 
  ! Y1c = sign_G * curvature .* (-B1s_over_B0 + B1c_over_B0 * sigma) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0);
  Jacobian(matrix_size, 1) = (- normal_cylindrical(1,1) * sign_G * curvature(1) * B1c_over_B0) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0) * d_l_d_phi(1) / R0(1)

end subroutine quasisymmetry_Jacobian

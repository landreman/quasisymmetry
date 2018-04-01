subroutine quasisymmetry_residual

  use quasisymmetry_variables

  implicit none

  integer :: j

  residual(1:N_phi) = matmul(d_d_zeta, sigma) ! Could use BLAS for speed if needed

  residual(1:N_phi) = residual(1:N_phi) &
       + iota * ( B1Squared_over_curvatureSquared * B1Squared_over_curvatureSquared + 1 + sigma * sigma) &
       - 2 * B1Squared_over_curvatureSquared * (torsion + I2_over_B0) / B0_over_abs_G0

  ! At zeta=0, we want theta=0 to correspond to the direction e_R.
  ! Therefore we impose Z1c(1) = 0.
  ! We use
  ! Z1c = ( binormal_cylindrical(:,1) .* X1c - normal_cylindrical(:,1) .* Y1c) .* d_l_d_phi ./ R0
  ! and
  ! X1c = B1c_over_B0 ./ curvature
  ! and 
  ! Y1c = sign_G * curvature .* (-B1s_over_B0 + B1c_over_B0 * sigma) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0);
  residual(matrix_size) = ( binormal_cylindrical(1,1) * B1c_over_B0 / curvature(1) &
       - normal_cylindrical(1,1) * sign_G * curvature(1) * (-B1s_over_B0 + B1c_over_B0 * sigma(1)) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0)) * d_l_d_phi(1) / R0(1)

end subroutine quasisymmetry_residual

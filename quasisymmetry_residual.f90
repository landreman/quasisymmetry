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
  ! Z1c = (-c .* (X1c .* normal_cylindrical(:,1) + Y1c .* binormal_cylindrical(:,1)) ...
  !    + a .* (X1c .* normal_cylindrical(:,3) + Y1c .* binormal_cylindrical(:,3))) ...
  !    ./ (a .* d - c .* b);
  ! and
  ! Y1c = sign_G * curvature .* (-B1s_over_B0 + B1c_over_B0 * W) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0);
  ! so
  ! Z1c(1) \propto (-c(1) * (X1c(1) * normal_cylindrical(1,1) + sign_G * curvature(1) * (-B1s_over_B0 + B1c_over_B0 * W(1)) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0) * binormal_cylindrical(1,1)) ...
  !                + a(1) * (X1c(1) * normal_cylindrical(1,3) + sign_G * curvature(1) * (-B1s_over_B0 + B1c_over_B0 * W(1)) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0) * binormal_cylindrical(1,3)))
  residual(matrix_size) = &
       (-RZ_to_XY_b(1) * (X1c(1) * normal_cylindrical(1,1) + sign_G * curvature(1) * (-B1s_over_B0 + B1c_over_B0 * sigma(1)) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0) * binormal_cylindrical(1,1)) &
       + RZ_to_XY_a(1) * (X1c(1) * normal_cylindrical(1,3) + sign_G * curvature(1) * (-B1s_over_B0 + B1c_over_B0 * sigma(1)) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0) * binormal_cylindrical(1,3)))

end subroutine quasisymmetry_residual

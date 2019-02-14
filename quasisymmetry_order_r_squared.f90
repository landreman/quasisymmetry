subroutine quasisymmetry_order_r_squared

  use quasisymmetry_variables

  implicit none

  real(dp), dimension(:), allocatable :: V1, V2, V3, qs, qc, rs, rc
  real(dp), dimension(:), allocatable :: Y2s_from_X20, Y2s_inhomogeneous, Y2c_from_X20, Y2c_inhomogeneous
  real(dp), dimension(:), allocatable :: fX0_from_X20, fX0_from_Y20, fX0_inhomogeneous
  real(dp), dimension(:), allocatable :: fXs_from_X20, fXs_from_Y20, fXs_inhomogeneous
  real(dp), dimension(:), allocatable :: fXc_from_X20, fXc_from_Y20, fXc_inhomogeneous
  real(dp), dimension(:), allocatable :: fY0_from_X20, fY0_from_Y20, fY0_inhomogeneous
  real(dp), dimension(:), allocatable :: fYs_from_X20, fYs_from_Y20, fYs_inhomogeneous
  real(dp), dimension(:), allocatable :: fYc_from_X20, fYc_from_Y20, fYc_inhomogeneous
  real(dp) :: factor, iota_N, abs_G0_over_B0, beta_1s, normalizer
  real(dp), dimension(:,:), allocatable :: matrix
  real(dp), dimension(:), allocatable :: right_hand_side
  integer :: j, iunit=20
  ! Variables needed by LAPACK:
  integer :: INFO
  integer, dimension(:), allocatable :: IPIV

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  iota_N = iota + axis_helicity*nfp
  abs_G0_over_B0 = 1 / B0_over_abs_G0

  if ((verbose) .and. abs(iota_N) < 1e-8) print "(a,es20.14)","Warning: |iota_N| is very small so O(r^2) solve will be poorly conditioned. iota_N=",iota_N

  if (allocated(X20)) deallocate(X20)
  if (allocated(Y20)) deallocate(Y20)
  if (allocated(Z20)) deallocate(Z20)
  if (allocated(X2s)) deallocate(X2s)
  if (allocated(Y2s)) deallocate(Y2s)
  if (allocated(Z2s)) deallocate(Z2s)
  if (allocated(X2c)) deallocate(X2c)
  if (allocated(Y2c)) deallocate(Y2c)
  if (allocated(Z2c)) deallocate(Z2c)

  if (allocated(X20_untwisted)) deallocate(X20_untwisted)
  if (allocated(Y20_untwisted)) deallocate(Y20_untwisted)
  if (allocated(Z20_untwisted)) deallocate(Z20_untwisted)
  if (allocated(X2s_untwisted)) deallocate(X2s_untwisted)
  if (allocated(Y2s_untwisted)) deallocate(Y2s_untwisted)
  if (allocated(Z2s_untwisted)) deallocate(Z2s_untwisted)
  if (allocated(X2c_untwisted)) deallocate(X2c_untwisted)
  if (allocated(Y2c_untwisted)) deallocate(Y2c_untwisted)
  if (allocated(Z2c_untwisted)) deallocate(Z2c_untwisted)

  allocate(X20(N_phi))
  allocate(X2s(N_phi))
  allocate(X2c(N_phi))
  allocate(Y20(N_phi))
  allocate(Y2s(N_phi))
  allocate(Y2c(N_phi))
  allocate(Z20(N_phi))
  allocate(Z2s(N_phi))
  allocate(Z2c(N_phi))

  allocate(X20_untwisted(N_phi))
  allocate(X2s_untwisted(N_phi))
  allocate(X2c_untwisted(N_phi))
  allocate(Y20_untwisted(N_phi))
  allocate(Y2s_untwisted(N_phi))
  allocate(Y2c_untwisted(N_phi))
  allocate(Z20_untwisted(N_phi))
  allocate(Z2s_untwisted(N_phi))
  allocate(Z2c_untwisted(N_phi))

  allocate(V1(N_phi))
  allocate(V2(N_phi))
  allocate(V3(N_phi))
  allocate(qs(N_phi))
  allocate(qc(N_phi))
  allocate(rs(N_phi))
  allocate(rc(N_phi))

  V1 = X1c * X1c + Y1c * Y1c + Y1s * Y1s
  V2 = 2 * Y1s * Y1c
  V3 = X1c * X1c + Y1c * Y1c - Y1s * Y1s

  ! The "matmul"s that follow could be sped up using BLAS2
  factor = - B0_over_abs_G0 / 8;
  Z20 = factor*matmul(d_d_zeta,V1)
  Z2s = factor*(matmul(d_d_zeta,V2) - 2 * iota_N * V3)
  Z2c = factor*(matmul(d_d_zeta,V3) + 2 * iota_N * V2)

  qs = -iota_N * X1c - Y1s * torsion * abs_G0_over_B0
  qc = matmul(d_d_zeta,X1c) - Y1c * torsion * abs_G0_over_B0
  rs = matmul(d_d_zeta,Y1s) - iota_N * Y1c
  rc = matmul(d_d_zeta,Y1c) + iota_N * Y1s + X1c * torsion * abs_G0_over_B0

  X2s = B0_over_abs_G0 * (matmul(d_d_zeta,Z2s) - 2*iota_N*Z2c + B0_over_abs_G0 * ( abs_G0_over_B0*abs_G0_over_B0*B2s/B0 + (qc * qs + rc * rs)/2)) / curvature

  X2c = B0_over_abs_G0 * (matmul(d_d_zeta,Z2c) + 2*iota_N*Z2s - B0_over_abs_G0 * (-abs_G0_over_B0*abs_G0_over_B0*B2c/B0 &
       + abs_G0_over_B0*abs_G0_over_B0*eta_bar*eta_bar/2 - (qc * qc - qs * qs + rc * rc - rs * rs)/4)) / curvature

  beta_1s = -4 * sign_psi * sign_G * mu0 * p2 * eta_bar * abs_G0_over_B0 / (iota_N * B0 * B0)

  allocate(Y2s_from_X20(N_phi))
  allocate(Y2s_inhomogeneous(N_phi))
  allocate(Y2c_from_X20(N_phi))
  allocate(Y2c_inhomogeneous(N_phi))

  allocate(fX0_from_X20(N_phi))
  allocate(fX0_from_Y20(N_phi))
  allocate(fX0_inhomogeneous(N_phi))

  allocate(fXs_from_X20(N_phi))
  allocate(fXs_from_Y20(N_phi))
  allocate(fXs_inhomogeneous(N_phi))

  allocate(fXc_from_X20(N_phi))
  allocate(fXc_from_Y20(N_phi))
  allocate(fXc_inhomogeneous(N_phi))

  allocate(fY0_from_X20(N_phi))
  allocate(fY0_from_Y20(N_phi))
  allocate(fY0_inhomogeneous(N_phi))

  allocate(fYs_from_X20(N_phi))
  allocate(fYs_from_Y20(N_phi))
  allocate(fYs_inhomogeneous(N_phi))

  allocate(fYc_from_X20(N_phi))
  allocate(fYc_from_Y20(N_phi))
  allocate(fYc_inhomogeneous(N_phi))


  Y2s_from_X20 = -sign_G * sign_psi * curvature * curvature / (eta_bar * eta_bar)
  Y2s_inhomogeneous = sign_G * sign_psi * (-curvature/2 + curvature*curvature/(eta_bar*eta_bar)*(-X2c + X2s * sigma))

  Y2c_from_X20 = -sign_G * sign_psi * curvature * curvature * sigma / (eta_bar * eta_bar)
  Y2c_inhomogeneous = sign_G * sign_psi * curvature * curvature / (eta_bar * eta_bar) * (X2s + X2c * sigma)

  ! Note: in the fX* and fY* quantities below, I've omitted the contributions from X20 and Y20 to the d/dzeta terms. These contributions are
  ! handled later when we assemble the large matrix.

  fX0_from_X20 = -4 * sign_G * sign_psi * abs_G0_over_B0 * (Y2c_from_X20 * Z2s - Y2s_from_X20 * Z2c)
  fX0_from_Y20 = -torsion * abs_G0_over_B0 - 4 * sign_G * sign_psi * abs_G0_over_B0 * (Z2s) &
       - sign_psi * I2_over_B0 * (-2) * abs_G0_over_B0
  fX0_inhomogeneous = curvature * abs_G0_over_B0 * Z20 - 4 * sign_G * sign_psi * abs_G0_over_B0 * (Y2c_inhomogeneous * Z2s - Y2s_inhomogeneous * Z2c) &
       - sign_psi * I2_over_B0 * (0.5d+0 * curvature * sign_G * sign_psi) * abs_G0_over_B0 + beta_1s * abs_G0_over_B0 / 2 * Y1c

  fXs_from_X20 = -torsion * abs_G0_over_B0 * Y2s_from_X20 - 4 * sign_psi * sign_G * abs_G0_over_B0 * (Y2c_from_X20 * Z20) &
       - sign_psi * I2_over_B0 * (- 2 * Y2s_from_X20) * abs_G0_over_B0
  fXs_from_Y20 = - 4 * sign_psi * sign_G * abs_G0_over_B0 * (-Z2c + Z20)
  fXs_inhomogeneous = matmul(d_d_zeta,X2s) - 2 * iota_N * X2c - torsion * abs_G0_over_B0 * Y2s_inhomogeneous + curvature * abs_G0_over_B0 * Z2s &
       - 4 * sign_psi * sign_G * abs_G0_over_B0 * (Y2c_inhomogeneous * Z20) &
       - sign_psi * I2_over_B0 * (0.5d+0 * curvature * sign_psi * sign_G - 2 * Y2s_inhomogeneous) * abs_G0_over_B0 &
       - (0.5d+0) * abs_G0_over_B0 * beta_1s * Y1s

  fXc_from_X20 = - torsion * abs_G0_over_B0 * Y2c_from_X20 - 4 * sign_psi * sign_G * abs_G0_over_B0 * (-Y2s_from_X20 * Z20) &
       - sign_psi * I2_over_B0 * (- 2 * Y2c_from_X20) * abs_G0_over_B0
  fXc_from_Y20 = - torsion * abs_G0_over_B0 - 4 * sign_psi * sign_G * abs_G0_over_B0 * (Z2s) &
       - sign_psi * I2_over_B0 * (-2) * abs_G0_over_B0
  fXc_inhomogeneous = matmul(d_d_zeta,X2c) + 2 * iota_N * X2s - torsion * abs_G0_over_B0 * Y2c_inhomogeneous + curvature * abs_G0_over_B0 * Z2c &
       - 4 * sign_psi * sign_G * abs_G0_over_B0 * (-Y2s_inhomogeneous * Z20) &
       - sign_psi * I2_over_B0 * (0.5d+0 * curvature * sign_G * sign_psi - 2 * Y2c_inhomogeneous) * abs_G0_over_B0 &
       - (0.5d+0) * abs_G0_over_B0 * beta_1s * Y1c

  fY0_from_X20 = torsion * abs_G0_over_B0 - sign_psi * I2_over_B0 * (2) * abs_G0_over_B0
  fY0_from_Y20 = 0
  fY0_inhomogeneous = -4 * sign_psi * sign_G * abs_G0_over_B0 * (X2s * Z2c - X2c * Z2s) &
       - sign_psi * I2_over_B0 * (-0.5d+0 * curvature * X1c * X1c) * abs_G0_over_B0 - (0.5d+0) * abs_G0_over_B0 * beta_1s * X1c

  fYs_from_X20 = -2 * iota_N * Y2c_from_X20 - 4 * sign_psi * sign_G * abs_G0_over_B0 * (Z2c)
  fYs_from_Y20 = -2 * iota_N
  fYs_inhomogeneous = matmul(d_d_zeta,Y2s_inhomogeneous) - 2 * iota_N * Y2c_inhomogeneous + torsion * abs_G0_over_B0 * X2s &
       - 4 * sign_psi * sign_G * abs_G0_over_B0 * (-X2c * Z20) - 2 * sign_psi * I2_over_B0 * X2s * abs_G0_over_B0

  fYc_from_X20 = 2 * iota_N * Y2s_from_X20 - 4 * sign_psi * sign_G * abs_G0_over_B0 * (-Z2s)
  fYc_from_Y20 = 0
  fYc_inhomogeneous = matmul(d_d_zeta,Y2c_inhomogeneous) + 2 * iota_N * Y2s_inhomogeneous + torsion * abs_G0_over_B0 * X2c &
       - 4 * sign_psi * sign_G * abs_G0_over_B0 * (X2s * Z20) &
       - sign_psi * I2_over_B0 * (-0.5d+0 * curvature * X1c * X1c + 2 * X2c) * abs_G0_over_B0 + 0.5d+0 * abs_G0_over_B0 * beta_1s * X1c

  allocate(matrix(2*N_phi,2*N_phi))
  allocate(right_hand_side(2*N_phi))

  matrix = 0
  do j = 1, N_phi
     ! Handle the terms involving d X_0 / d zeta and d Y_0 / d zeta:
     ! ----------------------------------------------------------------

     ! Equation 1, terms involving X0:
     ! Contributions arise from Y1c * fYs - Y1s * fYc.
     matrix(j,1:N_phi) = Y1c(j) * d_d_zeta(j,:) * Y2s_from_X20 - Y1s(j) * d_d_zeta(j,:) * Y2c_from_X20

     ! Equation 1, terms involving Y0:
     ! Contributions arise from -Y1s * fY0 - Y1s * fYc, and they happen to be equal.
     matrix(j,(N_phi+1):(2*N_phi)) = -2 * Y1s(j) * d_d_zeta(j,:)

     ! Equation 2, terms involving X0:
     ! Contributions arise from -X1c * fX0 + Y1s * fYs + Y1c * fYc
     matrix(j+N_phi,1:N_phi) = -X1c(j) * d_d_zeta(j,:) + Y1s(j) * d_d_zeta(j,:) * Y2s_from_X20 + Y1c(j) * d_d_zeta(j,:) * Y2c_from_X20

     ! Equation 2, terms involving Y0:
     ! Contributions arise from -Y1c * fY0 + Y1c * fYc, but they happen to cancel.

     ! Now handle the terms involving X_0 and Y_0 without d/dzeta derivatives:
     ! ----------------------------------------------------------------

     matrix(j,j      ) = matrix(j,j      ) + X1c(j) * fXs_from_X20(j) - Y1s(j) * fY0_from_X20(j) + Y1c(j) * fYs_from_X20(j) - Y1s(j) * fYc_from_X20(j)
     matrix(j,j+N_phi) = matrix(j,j+N_phi) + X1c(j) * fXs_from_Y20(j) - Y1s(j) * fY0_from_Y20(j) + Y1c(j) * fYs_from_Y20(j) - Y1s(j) * fYc_from_Y20(j)

     matrix(j+N_phi,j      ) = matrix(j+N_phi,j      ) - X1c(j) * fX0_from_X20(j) + X1c(j) * fXc_from_X20(j) - Y1c(j) * fY0_from_X20(j) + Y1s(j) * fYs_from_X20(j) + Y1c(j) * fYc_from_X20(j)
     matrix(j+N_phi,j+N_phi) = matrix(j+N_phi,j+N_phi) - X1c(j) * fX0_from_Y20(j) + X1c(j) * fXc_from_Y20(j) - Y1c(j) * fY0_from_Y20(j) + Y1s(j) * fYs_from_Y20(j) + Y1c(j) * fYc_from_Y20(j)
  end do

  right_hand_side(1:N_phi) = -(X1c * fXs_inhomogeneous - Y1s * fY0_inhomogeneous + Y1c * fYs_inhomogeneous - Y1s * fYc_inhomogeneous)
  right_hand_side((N_phi+1):(2*N_phi)) = -(- X1c * fX0_inhomogeneous + X1c * fXc_inhomogeneous - Y1c * fY0_inhomogeneous + Y1s * fYs_inhomogeneous + Y1c * fYc_inhomogeneous)

  print *,"Here comes abs_G0_over_B0:",abs_G0_over_B0
  print *,"Here comes X1c:"
  print *,X1c
  print *,"Here comes X1s:"
  print *,X1s
  print *,"Here comes Y1c:"
  print *,Y1c
  print *,"Here comes Y1s:"
  print *,Y1s
  print *,"Here comes torsion:"
  print *,torsion

  print *,"Here comes Z20:"
  print *,Z20
  print *,"Here comes Z2s:"
  print *,Z2s
  print *,"Here comes Z2c:"
  print *,Z2c

  print *,"Here comes qs:"
  print *,qs
  print *,"Here comes qc:"
  print *,qc
  print *,"Here comes rs:"
  print *,rs
  print *,"Here comes rc:"
  print *,rc

  print *,"Here comes X2c:"
  print *,X2c
  print *,"Here comes X2s:"
  print *,X2s
  print *,"Here comes Y2c_inhomogeneous:"
  print *,Y2c_inhomogeneous
  print *,"Here comes Y2s_inhomogeneous:"
  print *,Y2s_inhomogeneous
  print *," "
  print *,"Here comes fX0_inhomogeneous:"
  print *,fX0_inhomogeneous
  print *,"Here comes fXs_inhomogeneous:"
  print *,fXs_inhomogeneous
  print *,"Here comes fXc_inhomogeneous:"
  print *,fXc_inhomogeneous
  print *,"Here comes fY0_inhomogeneous:"
  print *,fY0_inhomogeneous
  print *,"Here comes fYs_inhomogeneous:"
  print *,fYs_inhomogeneous
  print *,"Here comes fYc_inhomogeneous:"
  print *,fYc_inhomogeneous

  print *,"Here comes right_hand_side:"
  print *,right_hand_side

  open(unit=iunit,file="matrix.dat")
  do j = 1, N_phi*2
     write(iunit, "(*(es24.15))") matrix(j,:)
  end do
  close(iunit)

  ! We will use the LAPACK subroutine DGESV to solve a general (asymmetric) linear system
  ! solution = matrix \ right_hand_side
  ! Note that LAPACK will over-write "right_hand_side" with the solution, and over-write "matrix" with the LU factorization.
  allocate(IPIV(2*N_phi))
  call DGESV(2*N_phi, 1, matrix, 2*N_phi, IPIV, right_hand_side, 2*N_phi, INFO)
  deallocate(IPIV)
  if (INFO /= 0) then
     print *, "Error in LAPACK call DGESV: info = ", INFO
     stop
  end if

  X20 = right_hand_side(1:N_phi)
  Y20 = right_hand_side((N_phi+1):(2*N_phi))

  B20 = B0 * (curvature * X20 - B0_over_abs_G0 * matmul(d_d_zeta,Z20) + (0.5d+0) * eta_bar * eta_bar - mu0 * p2 / (B0 * B0) &
       - (0.25d+0) * B0_over_abs_G0 * B0_over_abs_G0 * (qc * qc + qs * qs + rc * rc + rs * rs))

  print *,"Here comes X20:"
  print *,X20
  print *,"Here comes Y20:"
  print *,Y20
  print *,"Here comes B20:"
  print *,B20

  normalizer = 1 / sum(d_l_d_phi)
  B20_mean = sum(B20 * d_l_d_phi) * normalizer
  B20_residual = sqrt(sum((B20 - B20_mean) * (B20 - B20_mean) * d_l_d_phi) * normalizer)

  deallocate(V1, V2, V3, qs, qc, rs, rc)
  deallocate(matrix, right_hand_side)
  deallocate(Y2s_from_X20, Y2s_inhomogeneous, Y2c_from_X20, Y2c_inhomogeneous)
  deallocate(fX0_from_X20, fX0_from_Y20, fX0_inhomogeneous)
  deallocate(fXs_from_X20, fXs_from_Y20, fXs_inhomogeneous)
  deallocate(fXc_from_X20, fXc_from_Y20, fXc_inhomogeneous)
  deallocate(fY0_from_X20, fY0_from_Y20, fY0_inhomogeneous)
  deallocate(fYs_from_X20, fYs_from_Y20, fYs_inhomogeneous)
  deallocate(fYc_from_X20, fYc_from_Y20, fYc_inhomogeneous)

end subroutine quasisymmetry_order_r_squared
! Equations for this subroutine are derived in
! 20190727-01 Error field for O(r1) GarrenBoozer.docx
! and
! 20200120-01 Computing B2 from r1 quasisymmetry construction.nb
subroutine quasisymmetry_compute_B2_for_r1

  use quasisymmetry_variables

  implicit none

  real(dp), dimension(:), allocatable :: V1, V2, V3, qs, qc, rs, rc
  real(dp), dimension(:), allocatable :: fX0_from_X20, fX0_from_Y20, fX0_inhomogeneous
  real(dp), dimension(:), allocatable :: fXs_from_X20, fXs_from_Y20, fXs_inhomogeneous
  real(dp), dimension(:), allocatable :: fXc_from_X20, fXc_from_Y20, fXc_inhomogeneous
  real(dp), dimension(:), allocatable :: fY0_from_X20, fY0_from_Y20, fY0_inhomogeneous
  real(dp), dimension(:), allocatable :: fYs_from_X20, fYs_from_Y20, fYs_inhomogeneous
  real(dp), dimension(:), allocatable :: fYc_from_X20, fYc_from_Y20, fYc_inhomogeneous
  real(dp), allocatable, dimension(:) :: fX0, fXs, fXc, fY0, fYs, fYc, eq1residual, eq2residual
  real(dp) :: factor, iota_N, beta_1s, max_eq1residual, max_eq2residual, temp
  integer :: j, k, row
  real(dp), allocatable :: factors(:), residuals(:), matrix(:,:), right_hand_side(:)
  real(dp), allocatable :: matrix2(:,:), right_hand_side2(:)
  integer :: INFO
  integer, dimension(:), allocatable :: IPIV
  real(dp), allocatable, dimension(:) :: X2s_from_X20, X2s_from_Y20
  real(dp), allocatable, dimension(:) :: X2c_from_X20, X2c_inhomogeneous
  real(dp), allocatable, dimension(:) :: Y2s_from_X20, Y2s_from_Y20, Y2s_inhomogeneous
  real(dp), allocatable, dimension(:) :: Y2c_from_X20, Y2c_from_Y20, Y2c_inhomogeneous
  real(dp), allocatable, dimension(:) :: X2s_test, X2c_test
  real(dp), allocatable, dimension(:) :: Q, A, q_tilde
  real(dp) :: I2, G0


  if (verbose) print "(a)", " Hello from quasisymmetry_compute_B2_for_r1"

  if (allocated(X20)) deallocate(X20)
  if (allocated(X2s)) deallocate(X2s)
  if (allocated(X2c)) deallocate(X2c)
  if (allocated(Y20)) deallocate(Y20)
  if (allocated(Y2s)) deallocate(Y2s)
  if (allocated(Y2c)) deallocate(Y2c)
  if (allocated(Z20)) deallocate(Z20)
  if (allocated(Z2s)) deallocate(Z2s)
  if (allocated(Z2c)) deallocate(Z2c)

  allocate(X20(N_phi))
  allocate(X2s(N_phi))
  allocate(X2c(N_phi))
  allocate(Y20(N_phi))
  allocate(Y2s(N_phi))
  allocate(Y2c(N_phi))
  allocate(Z20(N_phi))
  allocate(Z2s(N_phi))
  allocate(Z2c(N_phi))

  allocate(V1(N_phi))
  allocate(V2(N_phi))
  allocate(V3(N_phi))
  allocate(qs(N_phi))
  allocate(qc(N_phi))
  allocate(rs(N_phi))
  allocate(rc(N_phi))

  if (allocated(B02)) deallocate(B02)
  allocate(B02(N_phi))

  allocate(factors(N_phi))

  iota_N = iota + axis_helicity*nfp
  abs_G0_over_B0 = 1 / B0_over_abs_G0

  V1 = X1c * X1c + Y1c * Y1c + Y1s * Y1s
  V2 = 2 * Y1s * Y1c
  V3 = X1c * X1c + Y1c * Y1c - Y1s * Y1s

  ! The "matmul"s that follow could be sped up using BLAS2
  factor = - B0_over_abs_G0 / 8;
  Z20 = factor*matmul(d_d_zeta,V1)
  Z2s = factor*(matmul(d_d_zeta,V2) - 2 * iota_N * V3)
  Z2c = factor*(matmul(d_d_zeta,V3) + 2 * iota_N * V2)

  allocate(d_Z20_d_zeta(N_phi))
  allocate(d_Z2s_d_zeta(N_phi))
  allocate(d_Z2c_d_zeta(N_phi))
  d_Z20_d_zeta = matmul(d_d_zeta, Z20)
  d_Z2s_d_zeta = matmul(d_d_zeta, Z2s)
  d_Z2c_d_zeta = matmul(d_d_zeta, Z2c)

  qs = -iota_N * X1c - Y1s * torsion * abs_G0_over_B0
  qc = matmul(d_d_zeta,X1c) - Y1c * torsion * abs_G0_over_B0
  rs = matmul(d_d_zeta,Y1s) - iota_N * Y1c
  rc = matmul(d_d_zeta,Y1c) + iota_N * Y1s + X1c * torsion * abs_G0_over_B0

  beta_1s = -4 * sign_psi * sign_G * mu0 * p2 * eta_bar * abs_G0_over_B0 / (iota_N * B0 * B0)

!  print *,"iota_N:",iota_N,", abs_G0_over_B0:",abs_G0_over_B0,", beta_1s:",beta_1s,", I2_over_B0:", I2_over_B0, ", B0:", B0,", sign_G:", sign_G,", sign_psi:",sign_psi

  ! If true, use the direct method involving a 6N * 6N linear system.
  ! If false, use the more complicated method involving only a 2N * 2N linear system.
  if (.false.) then
     ! Method with a 6N * 6N linear system.

     allocate(matrix(N_phi*6, N_phi*6))
     allocate(right_hand_side(N_phi*6))
     matrix = 0
     right_hand_side = 0

     ! Order of columns in the matrix, and of rows in the vector of unknowns: X20, X2s, X2c, Y20, Y2s, Y2c

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! First new equation, specific to the tilde equations when the construction is truncated at O(r^1)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     factors = 1 / (Y1s * Y1s + Y1c * Y1c)
     do j = 1, N_phi
        row = j
        ! X2c tilde term
        matrix(row, j + 2*N_phi) = -1 / X1c(j)
        ! Y2s tilde term
        matrix(row, j + 4*N_phi) = factors(j) * Y1s(j)
        ! Y2c tilde term
        matrix(row, j + 5*N_phi) = factors(j) * Y1c(j)
     end do
     ! There are no inhomogeneous terms for this equation.

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Second new equation, specific to the tilde equations when the construction is truncated at O(r^1)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 1, N_phi
        row = j + N_phi
        ! X2s tilde term
        matrix(row, j + 1*N_phi) = -1 / X1c(j)
        ! Y2s tilde term
        matrix(row, j + 4*N_phi) = factors(j) * Y1c(j)
        ! Y2c tilde term
        matrix(row, j + 5*N_phi) =-factors(j) * Y1s(j)
     end do
     ! There are no inhomogeneous terms for this equation.

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Eq (A32) of Landreman & Sengupta JPP (2019)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! -X1c * Y2c + X1c * Y20 + X2s * Y1s + X2c * Y1c - X20 * Y1c = 0
     do j = 1, N_phi
        row = j + 2 * N_phi
        ! X20 term
        matrix(row, j + 0*N_phi) = -Y1c(j)
        ! X2s term
        matrix(row, j + 1*N_phi) =  Y1s(j)
        ! X2c term
        matrix(row, j + 2*N_phi) =  Y1c(j)
        ! Y20 term
        matrix(row, j + 3*N_phi) =  X1c(j)
        ! Y2s term
        !matrix(row, j + 4*N_phi) = 0
        ! Y2c term
        matrix(row, j + 5*N_phi) = -X1c(j)
     end do
     ! There are no inhomogeneous terms for this equation.

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Eq (A33) of Landreman & Sengupta JPP (2019)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! X1c * Y2s + X2c * Y1s - X2s * Y1c + X20 * Y1s + sign_G * sign_psi * X1c * curvature / 2 = 0
     do j = 1, N_phi
        row = j + 3 * N_phi
        ! X20 term
        matrix(row, j + 0*N_phi) =  Y1s(j)
        ! X2s term
        matrix(row, j + 1*N_phi) = -Y1c(j)
        ! X2c term
        matrix(row, j + 2*N_phi) =  Y1s(j)
        ! Y20 term
        !matrix(row, j + 3*N_phi) = 0
        ! Y2s term
        matrix(row, j + 4*N_phi) =  X1c(j)
        ! Y2c term
        !matrix(row, j + 5*N_phi) = 0
     end do
     right_hand_side(3*N_phi+1 : 4*N_phi) = -sign_G * sign_psi * X1c * curvature * 0.5d+0
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Eq (A33) of Landreman & Sengupta JPP (2019)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     ! We will use the LAPACK subroutine DGESV to solve a general (asymmetric) linear system
     ! solution = matrix \ right_hand_side
     ! Note that LAPACK will over-write "right_hand_side" with the solution, and over-write "matrix" with the LU factorization.
     allocate(IPIV(6*N_phi))
     call DGESV(6*N_phi, 1, matrix, 6*N_phi, IPIV, right_hand_side, 6*N_phi, INFO)
     deallocate(IPIV)
     if (INFO /= 0) then
        print *, "Error in LAPACK call DGESV: info = ", INFO
        stop
     end if

     ! Unpack the solution vector
     X20 = right_hand_side(0*N_phi+1 : 1*N_phi)
     X2s = right_hand_side(1*N_phi+1 : 2*N_phi)
     X2c = right_hand_side(2*N_phi+1 : 3*N_phi)
     Y20 = right_hand_side(3*N_phi+1 : 4*N_phi)
     Y2s = right_hand_side(4*N_phi+1 : 5*N_phi)
     Y2c = right_hand_side(5*N_phi+1 : 6*N_phi)

  else
     ! Method with a 2N * 2N linear system.

     ! Order of columns in the matrix, and of rows in the vector of unknowns: [X20; Y20]
     allocate(matrix(N_phi*2, N_phi*2))
     allocate(right_hand_side(N_phi*2))
     matrix = 0
     right_hand_side = 0

     allocate(X2s_from_X20(N_phi))
     allocate(X2s_from_Y20(N_phi))

     allocate(X2c_from_X20(N_phi))
     allocate(X2c_inhomogeneous(N_phi))

     allocate(Y2s_from_X20(N_phi))
     allocate(Y2s_from_Y20(N_phi))
     allocate(Y2s_inhomogeneous(N_phi))

     allocate(Y2c_from_X20(N_phi))
     allocate(Y2c_from_Y20(N_phi))
     allocate(Y2c_inhomogeneous(N_phi))
     
     X2s_from_X20 =  Y1c / (2 * Y1s)
     X2s_from_Y20 = -X1c / (2 * Y1s)

     X2c_from_X20 = -0.5d+0
     X2c_inhomogeneous = -sign_G * sign_psi * X1c * curvature / (4 * Y1s)

     Y2s_from_X20 = (Y1c * Y1c - Y1s * Y1s) / (2 * X1c * Y1s)
     Y2s_from_Y20 = -Y1c / (2 * Y1s)
     Y2s_inhomogeneous = -sign_G * sign_psi * curvature * 0.25d+0

     Y2c_from_X20 = -Y1c / X1c
     Y2c_from_Y20 = 0.5d+0
     Y2c_inhomogeneous = -sign_G * sign_psi * Y1c * curvature / (4 * Y1s)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Eq (A41) of Landreman & Sengupta JPP (2019)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Part of this equation that depends on X20
     do j = 1, N_phi
        do k = 1, N_phi
           ! d/dphi terms
           temp = (X1c(k) * X2s_from_X20(j) + Y1c(k) * Y2s_from_X20(j) - Y1s(k) * Y2c_from_X20(j)) * d_d_zeta(k,j)
           ! Terms that do not involve d/dphi
           if (j==k) then
              temp = temp &
                   + 0 & ! X1c * fXs
                   + X1c(k) * (-2 * iota_N * X2c_from_X20(k) - torsion(k) * abs_G0_over_B0 * Y2s_from_X20(k) -4 * sign_G * sign_psi * abs_G0_over_B0 * (Y2c_from_X20(k) * Z20(k))&
                   - I2_over_B0 * sign_psi * (-2 * Y2s_from_X20(k)) * abs_G0_over_B0) &
                   + 0 & ! - Y1s * fY0
                   - Y1s(k) * (torsion(k) * abs_G0_over_B0 - 4 * sign_G * sign_psi * abs_G0_over_B0 * (X2s_from_X20(k) * Z2c(k) - X2c_from_X20(k) * Z2s(k)) &
                   - I2_over_B0 * sign_psi * 2 * abs_G0_over_B0) &
                   + 0 & ! + Y1c * fYs
                   + Y1c(k) * (-2 * iota_N * Y2c_from_X20(k) + torsion(k) * abs_G0_over_B0 * X2s_from_X20(k) - 4 * sign_G * sign_psi * abs_G0_over_B0 * (Z2c(k) - X2c_from_X20(k) * Z20(k)) &
                   - I2_over_B0 * sign_psi * 2 * X2s_from_X20(k) * abs_G0_over_B0) &
                   + 0 & ! - Y1s * fYc
                   - Y1s(k) * (2 * iota_N * Y2s_from_X20(k) + torsion(k) * abs_G0_over_B0 * X2c_from_X20(k) - 4 * sign_G * sign_psi * abs_G0_over_B0 * (X2s_from_X20(k) * Z20(k) - Z2s(k)) &
                   - I2_over_B0 * sign_psi * 2 * X2c_from_X20(k) * abs_G0_over_B0)
           end if
           matrix(k,j) = temp
        end do
     end do

     ! Part of this equation that depends on Y20. Note X2c does not depend on Y20.
     do j = 1, N_phi
        do k = 1, N_phi
           ! d/dphi terms
           temp = (X1c(k) * X2s_from_Y20(j) - Y1s(k) + Y1c(k) * Y2s_from_Y20(j) - Y1s(k) * Y2c_from_Y20(j)) * d_d_zeta(k,j)
           ! Terms that do not involve d/dphi
           if (j==k) then
              temp = temp &
                   + 0 & ! X1c * fXs
                   + X1c(k) * (- torsion(k) * abs_G0_over_B0 * Y2s_from_Y20(k) - 4 * sign_G * sign_psi * abs_G0_over_B0 * (-Z2c(k) + Y2c_from_Y20(k) * Z20(k)) &
                   - I2_over_B0 * sign_psi * (-2 * Y2s_from_Y20(k)) * abs_G0_over_B0) &
                   + 0 & ! - Y1s * fY0
                   - Y1s(k) * ( - 4 * sign_G * sign_psi * abs_G0_over_B0 * (X2s_from_Y20(k) * Z2c(k) )) &
                   + 0 & ! + Y1c * fYs
                   + Y1c(k) * (-2 * iota_N * Y2c_from_Y20(k) + torsion(k) * abs_G0_over_B0 * X2s_from_Y20(k) &
                   - I2_over_B0 * sign_psi * 2 * X2s_from_Y20(k) * abs_G0_over_B0) &
                   + 0 & ! - Y1s * fYc
                   - Y1s(k) * (2 * iota_N * Y2s_from_Y20(k) - 4 * sign_G * sign_psi * abs_G0_over_B0 * (X2s_from_Y20(k) * Z20(k)))
           end if
           matrix(k,j + N_phi) = temp
        end do
     end do

     ! "Inhomogeneous" part of this term, i.e. the part that is independent of X20 and Y20
     ! Note X2s has no inhomogeneous term. The matmuls below could perhaps be sped up.
     right_hand_side(1:N_phi) = -(&
          0 & !  X1c * fXs
          + X1c * (-2 * iota_N * X2c_inhomogeneous - torsion * abs_G0_over_B0 * Y2s_inhomogeneous + curvature * abs_G0_over_B0 * Z2s &
          - 4 * sign_G * sign_psi * abs_G0_over_B0 * (Y2c_inhomogeneous * Z20) - I2_over_B0 * sign_psi * (curvature * 0.5d+0 * X1c * Y1s - 2 * Y2s_inhomogeneous) * abs_G0_over_B0 &
           - 0.5d+0 * abs_G0_over_B0 * beta_1s * Y1s) &
          + 0 & ! -Y1s * fY0
          - Y1s * (- 4 * sign_G * sign_psi * abs_G0_over_B0 * (-X2c_inhomogeneous * Z2s) - I2_over_B0 * sign_psi * (-curvature * 0.5d+0 * X1c * X1c) * abs_G0_over_B0 &
          - 0.5d+0 * abs_G0_over_B0 * beta_1s * X1c) &
          + 0 & ! +Y1c * fYs
          + Y1c * (matmul(d_d_zeta,Y2s_inhomogeneous) - 2 * iota_N * Y2c_inhomogeneous &
          - 4 * sign_G * sign_psi * abs_G0_over_B0 * (-X2c_inhomogeneous * Z20)) &
          + 0 & ! -Y1s * fYc
          - Y1s * (matmul(d_d_zeta,Y2c_inhomogeneous) + 2 * iota_N * Y2s_inhomogeneous + torsion * abs_G0_over_B0 * X2c_inhomogeneous &
          - I2_over_B0 * sign_psi * (-curvature * 0.5d+0 * X1c * X1c + 2 * X2c_inhomogeneous) * abs_G0_over_B0 + 0.5d+0 * abs_G0_over_B0 * beta_1s * X1c) &
          )

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Eq (A42) of Landreman & Sengupta JPP (2019)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Part of this equation that depends on X20
     do j = 1, N_phi
        do k = 1, N_phi
           ! d/dphi terms
           temp = (-X1c(k) + X1c(k) * X2c_from_X20(j) + Y1s(k) * Y2s_from_X20(j) + Y1c(k) * Y2c_from_X20(j)) * d_d_zeta(k,j)
           ! Terms that do not involve d/dphi
           if (j==k) then
              temp = temp &
                   + 0 & ! -X1c * fX0
                   - X1c(k) * (- 4 * sign_G * sign_psi * abs_G0_over_B0 * (Y2c_from_X20(k) * Z2s(k) - Y2s_from_X20(k) * Z2c(k))) &
                   + 0 & ! +X1c * fXc
                   + X1c(k) * (2 * iota_N * X2s_from_X20(k) - torsion(k) * abs_G0_over_B0 * Y2c_from_X20(k) - 4 * sign_G * sign_psi * abs_G0_over_B0 * (-Y2s_from_X20(k) * Z20(k)) &
                   - I2_over_B0 * sign_psi * (-2 * Y2c_from_X20(k)) * abs_G0_over_B0) &
                   + 0 & ! -Y1c * fY0
                   - Y1c(k) * (torsion(k) * abs_G0_over_B0 - 4 * sign_G * sign_psi * abs_G0_over_B0 * (X2s_from_X20(k) * Z2c(k) - X2c_from_X20(k) * Z2s(k)) &
                   - I2_over_B0 * sign_psi * 2 * abs_G0_over_B0) &
                   + 0 & ! +Y1s * fYs
                   + Y1s(k) * (-2 * iota_N * Y2c_from_X20(k) + torsion(k) * abs_G0_over_B0 * X2s_from_X20(k) - 4 * sign_G * sign_psi * abs_G0_over_B0 * (Z2c(k) - X2c_from_X20(k) * Z20(k)) &
                   - I2_over_B0 * sign_psi * 2 * X2s_from_X20(k) * abs_G0_over_B0) &
                   + 0 & ! +Y1c * fYc
                   + Y1c(k) * (2 * iota_N * Y2s_from_X20(k) + torsion(k) * abs_G0_over_B0 * X2c_from_X20(k) - 4 * sign_G * sign_psi * abs_G0_over_B0 * (X2s_from_X20(k) * Z20(k) - Z2s(k)) &
                   - I2_over_B0 * sign_psi * 2 * X2c_from_X20(k) * abs_G0_over_B0)
           end if
           matrix(k + N_phi, j) = temp
        end do
     end do

     ! Part of this equation that depends on Y20. Note X2c does not depend on Y20.
     do j = 1, N_phi
        do k = 1, N_phi
           ! d/dphi terms
           temp = ( - Y1c(k) + Y1s(k) * Y2s_from_Y20(j) + Y1c(k) * Y2c_from_Y20(j)) * d_d_zeta(k,j)
           ! Terms that do not involve d/dphi
           if (j==k) then
              temp = temp &
                   + 0 & ! -X1c * fX0
                   - X1c(k) * (- torsion(k) * abs_G0_over_B0 - 4 * sign_G * sign_psi * abs_G0_over_B0 * (Y2c_from_Y20(k) * Z2s(k) - Y2s_from_Y20(k) * Z2c(k)) &
                   - I2_over_B0 * sign_psi * (-2) * abs_G0_over_B0) &
                   + 0 & ! +X1c * fXc
                   + X1c(k) * (2 * iota_N * X2s_from_Y20(k) - torsion(k) * abs_G0_over_B0 * Y2c_from_Y20(k) - 4 * sign_G * sign_psi * abs_G0_over_B0 * (Z2s(k) - Y2s_from_Y20(k) * Z20(k)) &
                   - I2_over_B0 * sign_psi * (-2 * Y2c_from_Y20(k)) * abs_G0_over_B0) &
                   + 0 & ! -Y1c * fY0
                   - Y1c(k) * ( - 4 * sign_G * sign_psi * abs_G0_over_B0 * (X2s_from_Y20(k) * Z2c(k) )) &
                   + 0 & ! +Y1s * fYs
                   + Y1s(k) * (-2 * iota_N * Y2c_from_Y20(k) + torsion(k) * abs_G0_over_B0 * X2s_from_Y20(k) &
                   - I2_over_B0 * sign_psi * 2 * X2s_from_Y20(k) * abs_G0_over_B0) &
                   + 0 & ! +Y1c * fYc
                   + Y1c(k) * (2 * iota_N * Y2s_from_Y20(k) - 4 * sign_G * sign_psi * abs_G0_over_B0 * (X2s_from_Y20(k) * Z20(k)))
           end if
           matrix(k + N_phi, j + N_phi) = temp
        end do
     end do

     ! "Inhomogeneous" part of this term, i.e. the part that is independent of X20 and Y20
     ! Note X2s has no inhomogeneous term. The matmuls below could perhaps be sped up.
     right_hand_side(N_phi + 1 : 2*N_phi) = -(&
          0 & ! -X1c * fX0
          - X1c * (curvature * abs_G0_over_B0 * Z20 - 4 * sign_G * sign_psi * abs_G0_over_B0 * (Y2c_inhomogeneous * Z2s - Y2s_inhomogeneous * Z2c) &
          - I2_over_B0 * sign_psi * (curvature * 0.5d+0 * X1c * Y1c) * abs_G0_over_B0 - 0.5d+0 * abs_G0_over_B0 * (-beta_1s * Y1c)) &
          + 0 & ! +X1c * fXc
          + X1c * (matmul(d_d_zeta, X2c_inhomogeneous) - torsion * abs_G0_over_B0 * Y2c_inhomogeneous + curvature * abs_G0_over_B0 * Z2c &
          - 4 * sign_G * sign_psi * abs_G0_over_B0 * (-Y2s_inhomogeneous * Z20) &
          - I2_over_B0 * sign_psi * (curvature * 0.5d+0 * X1c * Y1c - 2 * Y2c_inhomogeneous) * abs_G0_over_B0 - 0.5d+0 * abs_G0_over_B0 * beta_1s * Y1c) &
          + 0 & ! -Y1c * fY0
          - Y1c * (- 4 * sign_G * sign_psi * abs_G0_over_B0 * (-X2c_inhomogeneous * Z2s) - I2_over_B0 * sign_psi * (-curvature * 0.5d+0 * X1c * X1c) * abs_G0_over_B0 &
          - 0.5d+0 * abs_G0_over_B0 * beta_1s * X1c) &
          + 0 & ! +Y1s * fYs
          + Y1s * (matmul(d_d_zeta,Y2s_inhomogeneous) - 2 * iota_N * Y2c_inhomogeneous &
          - 4 * sign_G * sign_psi * abs_G0_over_B0 * (-X2c_inhomogeneous * Z20)) &
          + 0 & ! +Y1c * fYc
          + Y1c * (matmul(d_d_zeta,Y2c_inhomogeneous) + 2 * iota_N * Y2s_inhomogeneous + torsion * abs_G0_over_B0 * X2c_inhomogeneous &
          - I2_over_B0 * sign_psi * (-curvature * 0.5d+0 * X1c * X1c + 2 * X2c_inhomogeneous) * abs_G0_over_B0 + 0.5d+0 * abs_G0_over_B0 * beta_1s * X1c) &
          )

     ! Done assembling the matrix and right-hand side.

     ! 2nd pass typing up the equations
     if (debug) then
        allocate(right_hand_side2(2*N_phi))
        right_hand_side2 = 0
        ! Note: X2s has no inhomogeneous term
        right_hand_side2(1:N_phi) = -( &
             0   & !  X1c * fXs
             + X1c * (-2 * iota_N * X2c_inhomogeneous - torsion * abs_G0_over_B0 * Y2s_inhomogeneous + curvature * abs_G0_over_B0 * Z2s &
             - 4 * sign_G * sign_psi * abs_G0_over_B0 * (Y2c_inhomogeneous * Z20) - I2_over_B0 * sign_psi * (curvature * 0.5d+0 * X1c * Y1s - 2 * Y2s_inhomogeneous) * abs_G0_over_B0 &
             - 0.5d+0 * abs_G0_over_B0 * beta_1s * Y1s) &
             + 0 & ! -Y1s * fY0
             - Y1s * (- 4 *sign_G * sign_psi * abs_G0_over_B0 * (-X2c_inhomogeneous * Z2s) - I2_over_B0 * sign_psi * (-curvature * 0.5d+0 * X1c * X1c) * abs_G0_over_B0 &
             - 0.5d+0 * abs_G0_over_B0 * beta_1s * X1c) &
             + 0 & ! +Y1c * fYs
             + Y1c * (matmul(d_d_zeta, Y2s_inhomogeneous) - 2 * iota_N * Y2c_inhomogeneous - 4 *sign_G * sign_psi * abs_G0_over_B0 * (-X2c_inhomogeneous * Z20)) &
             + 0 & ! -Y1s * fYc
             - Y1s * (matmul(d_d_zeta, Y2c_inhomogeneous) + 2 * iota_N * Y2s_inhomogeneous + torsion * abs_G0_over_B0 * X2c_inhomogeneous &
             - I2_over_B0 * sign_psi * (curvature * 0.5d+0 * (-X1c * X1c) + 2 * X2c_inhomogeneous) * abs_G0_over_B0 + 0.5d+0 * abs_G0_over_B0 * beta_1s * X1c) &
             )
        
        right_hand_side2((N_phi+1):(2*N_phi)) = -( &
             0   & ! -X1c * fX0
             - X1c * (curvature * abs_G0_over_B0 * Z20 - 4 * sign_G * sign_psi * abs_G0_over_B0 * (Y2c_inhomogeneous * Z2s - Y2s_inhomogeneous * Z2c) &
             - I2_over_B0 * sign_psi * (curvature * 0.5d+0 * X1c * Y1c) * abs_G0_over_B0 - 0.5d+0 * abs_G0_over_B0 * (-beta_1s * Y1c)) &
             + 0 & ! +X1c * fXc
             + X1c * (matmul(d_d_zeta, X2c_inhomogeneous) - torsion * abs_G0_over_B0 * Y2c_inhomogeneous + curvature * abs_G0_over_B0 * Z2c &
             - 4 * sign_G * sign_psi * abs_G0_over_B0 * (-Y2s_inhomogeneous * Z20) - I2_over_B0 * sign_psi * (curvature * 0.5d+0 * X1c * Y1c - 2 * Y2c_inhomogeneous) * abs_G0_over_B0 &
             - 0.5d+0 * abs_G0_over_B0 * beta_1s * Y1c) &
             + 0 & ! -Y1c * fY0
             - Y1c * (- 4 *sign_G * sign_psi * abs_G0_over_B0 * (-X2c_inhomogeneous * Z2s) - I2_over_B0 * sign_psi * (-curvature * 0.5d+0 * X1c * X1c) * abs_G0_over_B0 &
             - 0.5d+0 * abs_G0_over_B0 * beta_1s * X1c) &
             + 0 & ! +Y1s * fYs
             + Y1s * (matmul(d_d_zeta, Y2s_inhomogeneous) - 2 * iota_N * Y2c_inhomogeneous - 4 *sign_G * sign_psi * abs_G0_over_B0 * (-X2c_inhomogeneous * Z20)) &
             + 0 & ! +Y1c * fYc
             + Y1c * (matmul(d_d_zeta, Y2c_inhomogeneous) + 2 * iota_N * Y2s_inhomogeneous + torsion * abs_G0_over_B0 * X2c_inhomogeneous &
             - I2_over_B0 * sign_psi * (curvature * 0.5d+0 * (-X1c * X1c) + 2 * X2c_inhomogeneous) * abs_G0_over_B0 + 0.5d+0 * abs_G0_over_B0 * beta_1s * X1c) &
             )
             
        ! print *,"Difference in RHS:",maxval(abs(right_hand_side(1:N_phi) - right_hand_side2(1:N_phi)))
        if (verbose) print *,"Difference in RHS:",maxval(abs(right_hand_side - right_hand_side2))
        deallocate(right_hand_side2)
     end if

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

     ! Unpack the solution vector
     X20 = right_hand_side(0*N_phi+1 : 1*N_phi)
     Y20 = right_hand_side(1*N_phi+1 : 2*N_phi)

     X2s = X2s_from_X20 * X20 + X2s_from_Y20 * Y20
     X2c = X2c_from_X20 * X20 + X2c_inhomogeneous
     Y2s = Y2s_from_X20 * X20 + Y2s_from_Y20 * Y20 + Y2s_inhomogeneous
     Y2c = Y2c_from_X20 * X20 + Y2c_from_Y20 * Y20 + Y2c_inhomogeneous

     deallocate(X2s_from_X20, X2s_from_Y20)
     deallocate(X2c_from_X20, X2c_inhomogeneous)
     deallocate(Y2s_from_X20, Y2s_from_Y20, Y2s_inhomogeneous)
     deallocate(Y2c_from_X20, Y2c_from_Y20, Y2c_inhomogeneous)

  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! At this point, we have completed solving for the tilde versions of
  ! X20, X2s, X2c, Y20, Y2s, Y2c.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (debug) then
     ! Verify that residuals in each equation are zero.
     allocate(eq1residual(N_phi))
     allocate(eq2residual(N_phi))

     ! Pair of new equations from the tilde analysis
     eq1residual = (X2s * X1s + X2c * X1c) / (X1s * X1s + X1c * X1c) - (Y2s * Y1s + Y2c * Y1c) / (Y1s * Y1s + Y1c * Y1c)
     eq2residual = (-X2c * X1s + X2s * X1c) / (X1s * X1s + X1c * X1c) - (-Y2c * Y1s + Y2s * Y1c) / (Y1s * Y1s + Y1c * Y1c)
     max_eq1residual = maxval(abs(eq1residual))
     max_eq2residual = maxval(abs(eq2residual))
     if (verbose) print *,"max(abs(residual)) for 1st equation from tilde analysis:",max_eq1residual
     if (verbose) print *,"max(abs(residual)) for 2nd equation from tilde analysis:",max_eq2residual
     if (max_eq1residual > 1e-12) stop "Large residual in 1st equation."
     if (max_eq2residual > 1e-12) stop "Large residual in 2nd equation."

     allocate(fX0(N_phi))
     allocate(fXs(N_phi))
     allocate(fXc(N_phi))
     allocate(fY0(N_phi))
     allocate(fYs(N_phi))
     allocate(fYc(N_phi))

     fX0 = matmul(d_d_zeta,X20) - torsion * abs_G0_over_B0 * Y20 + curvature * abs_G0_over_B0 * Z20 &
          -4*sign_G*sign_psi*abs_G0_over_B0*(Y2c * Z2s - Y2s * Z2c) &
          - sign_psi * I2_over_B0 * (curvature/2 * X1c * Y1c - 2 * Y20) * abs_G0_over_B0 + abs_G0_over_B0 * beta_1s * Y1c / 2

     fXs = matmul(d_d_zeta,X2s) - 2 * iota_N * X2c - torsion * abs_G0_over_B0 * Y2s + curvature * abs_G0_over_B0 * Z2s &
          -4*sign_G*sign_psi*abs_G0_over_B0*(-Y20 * Z2c + Y2c * Z20) &
          -sign_psi * I2_over_B0 * (curvature/2 * X1c * Y1s - 2 * Y2s) * abs_G0_over_B0 - abs_G0_over_B0 * beta_1s * Y1s / 2

     fXc = matmul(d_d_zeta,X2c) + 2 * iota_N * X2s - torsion * abs_G0_over_B0 * Y2c + curvature * abs_G0_over_B0 * Z2c &
          -4*sign_G*sign_psi*abs_G0_over_B0*(Y20 * Z2s - Y2s * Z20) &
          -sign_psi * I2_over_B0 * (curvature/2 * X1c * Y1c - 2 * Y2c) * abs_G0_over_B0 - abs_G0_over_B0 * beta_1s * Y1c / 2

     fY0 = matmul(d_d_zeta,Y20) + torsion * abs_G0_over_B0 * X20 - 4*sign_G*sign_psi*abs_G0_over_B0*(X2s * Z2c - X2c * Z2s) &
          -sign_psi * I2_over_B0 * (-curvature/2*X1c*X1c + 2*X20) * abs_G0_over_B0 - abs_G0_over_B0 * beta_1s * X1c / 2

     fYs = matmul(d_d_zeta,Y2s) - 2 * iota_N * Y2c + torsion * abs_G0_over_B0 * X2s &
          -4*sign_G*sign_psi*abs_G0_over_B0*(X20 * Z2c - X2c * Z20) - 2*sign_psi* I2_over_B0 * X2s * abs_G0_over_B0

     fYc = matmul(d_d_zeta,Y2c) + 2 * iota_N * Y2s + torsion * abs_G0_over_B0 * X2c &
          -4*sign_G*sign_psi*abs_G0_over_B0*(X2s * Z20 - X20 * Z2s) &
          -sign_psi * I2_over_B0 * (-curvature/2 * X1c * X1c + 2 * X2c) * abs_G0_over_B0 + abs_G0_over_B0 * beta_1s * X1c / 2


     eq1residual = X1c * fXs - Y1s * fY0 + Y1c * fYs - Y1s * fYc

     eq2residual = -X1c * fX0 + X1c * fXc - Y1c * fY0 + Y1s * fYs + Y1c * fYc

     max_eq1residual = maxval(abs(eq1residual))
     max_eq2residual = maxval(abs(eq2residual))
     if (verbose) print *,"max(abs(eq1residual)):",max_eq1residual
     if (verbose) print *,"max(abs(eq2residual)):",max_eq2residual

     if (max_eq1residual > 1e-6) stop "Equation 1 residual is large !!!"
     if (max_eq2residual > 1e-6) stop "Equation 2 residual is large !!!"

     ! Now check the two equations that were used to determine Y2s and Y2c:

     eq1residual = -X1c * Y2c + X1c * Y20 + X2s * Y1s + X2c * Y1c - X20 * Y1c

     eq2residual = X1c * Y2s + X2c * Y1s - X2s * Y1c + X20 * Y1s + sign_G * sign_psi * X1c * curvature / 2

     max_eq1residual = maxval(abs(eq1residual))
     max_eq2residual = maxval(abs(eq2residual))
     if (verbose) print *,"max(abs(Y2c eq residual)):",max_eq1residual
     if (verbose) print *,"max(abs(Y2s eq residual)):",max_eq2residual

     if (max_eq1residual > 1e-6) stop "Y2c equation residual is large !!!"
     if (max_eq2residual > 1e-6) stop "Y2s equation residual is large !!!"

     deallocate(fX0, fXs, fXc, fY0, fYs, fYc, eq1residual, eq2residual)
  end if

  ! Now that we know the tilde versions of X20, X2s, and X2c, we can compute the tilde versions of B20, B2s, and B2c

  B20       = B0 * (curvature * X20 - B0_over_abs_G0 * d_Z20_d_zeta + (0.5d+0) * eta_bar * eta_bar - mu0 * p2 / (B0 * B0) &
       - (0.25d+0) * B0_over_abs_G0 * B0_over_abs_G0 * (qc * qc + qs * qs + rc * rc + rs * rs))

  B2s_array = B0 * (curvature * X2s - B0_over_abs_G0 * (d_Z2s_d_zeta - 2 * iota_N * Z2c) &
       - (0.25d+0) * B0_over_abs_G0 * B0_over_abs_G0 * 2 * (qc * qs + rc * rs))

  B2c_array = B0 * (curvature * X2c - B0_over_abs_G0 * (d_Z2c_d_zeta + 2 * iota_N * Z2s) + (0.5d+0) * eta_bar * eta_bar &
       - (0.25d+0) * B0_over_abs_G0 * B0_over_abs_G0 * (qc * qc - qs * qs + rc * rc - rs * rs))


  if (debug) then
     ! Verify the B2s and B2c equations
     allocate(X2s_test(N_phi))
     allocate(X2c_test(N_phi))
     X2s_test = B0_over_abs_G0 * (d_Z2s_d_zeta - 2*iota_N*Z2c + B0_over_abs_G0 * ( abs_G0_over_B0*abs_G0_over_B0*B2s_array/B0 + (qc * qs + rc * rs)/2)) / curvature

     X2c_test = B0_over_abs_G0 * (d_Z2c_d_zeta + 2*iota_N*Z2s - B0_over_abs_G0 * (-abs_G0_over_B0*abs_G0_over_B0*B2c_array/B0 &
          + abs_G0_over_B0*abs_G0_over_B0*eta_bar*eta_bar/2 - (qc * qc - qs * qs + rc * rc - rs * rs)/4)) / curvature

     max_eq1residual = maxval(abs(X2s - X2s_test))
     max_eq2residual = maxval(abs(X2c - X2c_test))
     if (verbose) print *,"max(abs(residual)) in X2s equation:",max_eq1residual
     if (verbose) print *,"max(abs(residual)) in X2c equation:",max_eq2residual

     deallocate(X2s_test, X2c_test)
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now compute \tilde{B}_0^{(2)}
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(Q(N_phi))
  allocate(A(N_phi))
  allocate(q_tilde(N_phi))

  I2 = I2_over_B0 * B0
  G0 = sign_G * abs_G0_over_B0 * B0

  ! The expression below is derived in the O(r^2) paper, and in "20190318-01 Wrick's streamlined Garren-Boozer method, MHD.nb" in the section "Not assuming quasisymmetry".
  ! Note Q = (1/2) * (XYEquation0 without X3 and Y3 terms) where XYEquation0 is the quantity in the above notebook.
  Q = -sign_psi * B0 * abs_G0_over_B0 / (2*G0*G0) * (iota_N * I2 + mu0 * p2 * G0 / (B0 * B0)) + 2 * (X2c * Y2s - X2s * Y2c) &
       + sign_psi * B0 / (2*G0) * (abs_G0_over_B0 * X20 * curvature - d_Z20_d_zeta) &
       + I2 / (4 * G0) * (-abs_G0_over_B0 * torsion * (X1c*X1c + Y1s*Y1s + Y1c*Y1c) + Y1c * matmul(d_d_zeta,X1c) - X1c * matmul(d_d_zeta,Y1c))

  ! See "20200318-01 Checking equations in Error field in O(r1) GarrenBoozer.docx"
  ! A = torsion * (-Z2c * 2 * Y1s * Y1c - Z2s * (-X1c*X1c + Y1c*Y1c - Y1s*Y1s)) &
  A = torsion * (-Z2c * 2 * Y1s * Y1c - Z2s * (-X1c*X1c - Y1c*Y1c + Y1s*Y1s)) &
       -(Z2s / abs_G0_over_B0) * (-X1c * matmul(d_d_zeta,Y1c) + Y1c * matmul(d_d_zeta,X1c)) - (Z2c / abs_G0_over_B0) * (X1c * matmul(d_d_zeta,Y1s) - Y1s * matmul(d_d_zeta,X1c))

  q_tilde = Q + A/2 + 2 * (X2s*X2s + X2c*X2c) * (-Y1s)/X1c

  B02 = -sign_G * sign_psi * B0 * q_tilde

  max_B2tilde = max( &
       maxval(abs(B2s_array)), &
       maxval(abs(B2c_array)), &
       (maxval(B20) - minval(B20))/2, &
       (maxval(B02) - minval(B02))/2)

  if (verbose) print *,"max_B2tilde:",max_B2tilde

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute difference between the tilde and non-tilde poloidal angles:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(t1s(N_phi))
  allocate(t1c(N_phi))

  t1c =  2 * X2s / X1c
  t1s = -2 * X2c / X1c

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Clean up
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  deallocate(Q,A,q_tilde)

  deallocate(matrix, right_hand_side, factors)
  deallocate(V1, V2, V3, qs, qc, rs, rc)
  deallocate(d_Z20_d_zeta, d_Z2s_d_zeta, d_Z2c_d_zeta)

end subroutine quasisymmetry_compute_B2_for_r1

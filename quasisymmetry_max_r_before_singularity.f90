subroutine quasisymmetry_max_r_before_singularity(d_Z20_d_zeta, d_Z2s_d_zeta, d_Z2c_d_zeta)

  use quasisymmetry_variables, only: dp, N_phi, abs_G0_over_B0, X1c, Y1s, Y1c, X20, X2s, X2c, Y20, Y2s, Y2c, Z20, Z2s, Z2c, &
       curvature, torsion, r_singularity, r_singularity_vs_zeta, d_X1c_d_zeta, d_Y1s_d_zeta, d_Y1c_d_zeta, verbose, &
       general_option, general_option_single, &
       eta_bar, B2c, R0c, Z0s, axis_nmax

  implicit none

  real(dp), intent(in) :: d_Z20_d_zeta(N_phi), d_Z2s_d_zeta(N_phi), d_Z2c_d_zeta(N_phi)
  integer :: j, jr, varsigma, sign_quadratic
  real(dp) :: g0, g1c, g20, g2s, g2c, sigma_denominator, abs_costheta, denominator, r, residual, rc
  real(dp) :: K0, K2s, K2c, K4s, K4c, sintheta, costheta, sin2theta, cos2theta, lp
  real(dp) :: coefficients(5), real_parts(4), imag_parts(4)
  real(dp) :: quadratic_A, quadratic_B, quadratic_C, radical, sigma
  real(dp) :: start_time, end_time
  integer :: max_j_phi
  logical :: local_verbose = .false.

  call cpu_time(start_time)

  lp = abs_G0_over_B0 ! Shorthand

  max_j_phi = N_phi
  !if (local_verbose) max_j_phi = 1
  do j = 1, max_j_phi

     ! Write sqrt(g) = r * [g0 + r*g1c*cos(theta) + (r^2)*(g20 + g2s*sin(2*theta) + g2c*cos(2*theta) + ...]
     ! The coefficients are evaluated in "20200322-02 Max r for Garren Boozer.nb", in the section "Order r^2 construction, expanding"

     g0 = lp * X1c(j) * Y1s(j)

     !g1s = -2*X20(j)*Y1c(j) + 2*X2c(j)*Y1c(j) + 2*X2s(j)*Y1s(j) + 2*X1c(j)*Y20(j) - 2*X1c(j)*Y2c(j)

     g1c = -2*X2s(j)*Y1c(j) + 2*X20(j)*Y1s(j) + 2*X2c(j)*Y1s(j) + 2*X1c(j)*Y2s(j) - X1c(j)**2*Y1s(j)*curvature(j)

     g20 = -4*lp*X2s(j)*Y2c(j) + 4*lp*X2c(j)*Y2s(j) + lp*X1c(j)*X2s(j)*Y1c(j)*curvature(j) - &
          2*lp*X1c(j)*X20(j)*Y1s(j)*curvature(j) - lp*X1c(j)*X2c(j)*Y1s(j)*curvature(j) - &
          lp*X1c(j)**2*Y2s(j)*curvature(j) + 2*lp*Y1c(j)*Y1s(j)*Z2c(j)*torsion(j) - &
          lp*X1c(j)**2*Z2s(j)*torsion(j) - lp*Y1c(j)**2*Z2s(j)*torsion(j) + lp*Y1s(j)**2*Z2s(j)*torsion(j) - &
          Y1s(j)*Z20(j)*d_X1c_d_zeta(j) - Y1s(j)*Z2c(j)*d_X1c_d_zeta(j) + &
          Y1c(j)*Z2s(j)*d_X1c_d_zeta(j) - X1c(j)*Z2s(j)*d_Y1c_d_zeta(j) - &
          X1c(j)*Z20(j)*d_Y1s_d_zeta(j) + X1c(j)*Z2c(j)*d_Y1s_d_zeta(j) + &
          X1c(j)*Y1s(j)*d_Z20_d_zeta(j)

     g2c = -4*lp*X2s(j)*Y20(j) + 4*lp*X20(j)*Y2s(j) + &
          lp*X1c(j)*X2s(j)*Y1c(j)*curvature(j) - lp*X1c(j)*X20(j)*Y1s(j)*curvature(j) - &
          2*lp*X1c(j)*X2c(j)*Y1s(j)*curvature(j) - lp*X1c(j)**2*Y2s(j)*curvature(j) + &
          2*lp*Y1c(j)*Y1s(j)*Z20(j)*torsion(j) - lp*X1c(j)**2*Z2s(j)*torsion(j) - &
          lp*Y1c(j)**2*Z2s(j)*torsion(j) - lp*Y1s(j)**2*Z2s(j)*torsion(j) - &
          Y1s(j)*Z20(j)*d_X1c_d_zeta(j) - Y1s(j)*Z2c(j)*d_X1c_d_zeta(j) + &
          Y1c(j)*Z2s(j)*d_X1c_d_zeta(j) - X1c(j)*Z2s(j)*d_Y1c_d_zeta(j) + &
          X1c(j)*Z20(j)*d_Y1s_d_zeta(j) - X1c(j)*Z2c(j)*d_Y1s_d_zeta(j) + &
          X1c(j)*Y1s(j)*d_Z2c_d_zeta(j)

     g2s = 4*lp*X2c(j)*Y20(j) - 4*lp*X20(j)*Y2c(j) + &
          lp*X1c(j)*X20(j)*Y1c(j)*curvature(j) - lp*X1c(j)*X2c(j)*Y1c(j)*curvature(j) - &
          2*lp*X1c(j)*X2s(j)*Y1s(j)*curvature(j) - lp*X1c(j)**2*Y20(j)*curvature(j) + &
          lp*X1c(j)**2*Y2c(j)*curvature(j) - lp*X1c(j)**2*Z20(j)*torsion(j) - &
          lp*Y1c(j)**2*Z20(j)*torsion(j) + lp*Y1s(j)**2*Z20(j)*torsion(j) + &
          lp*X1c(j)**2*Z2c(j)*torsion(j) + lp*Y1c(j)**2*Z2c(j)*torsion(j) + &
          lp*Y1s(j)**2*Z2c(j)*torsion(j) + Y1c(j)*Z20(j)*d_X1c_d_zeta(j) - &
          Y1c(j)*Z2c(j)*d_X1c_d_zeta(j) - Y1s(j)*Z2s(j)*d_X1c_d_zeta(j) - &
          X1c(j)*Z20(j)*d_Y1c_d_zeta(j) + X1c(j)*Z2c(j)*d_Y1c_d_zeta(j) - &
          X1c(j)*Z2s(j)*d_Y1s_d_zeta(j) + X1c(j)*Y1s(j)*d_Z2s_d_zeta(j)

     ! We consider the system sqrt(g) = 0 and
     ! d (sqrtg) / d theta = 0.
     ! We algebraically eliminate r in "20200322-02 Max r for Garren Boozer.nb", in the section
     ! "Keeping first 3 orders in the Jacobian".
     ! We end up with the form in "20200322-01 Max r for GarrenBoozer.docx":
     ! K0 + K2s*sin(2*theta) + K2c*cos(2*theta) + K4s*sin(4*theta) + K4c*cos(4*theta) = 0.

     K0 = 2*g1c**2*g20 - 3*g1c**2*g2c + 8*g0*g2c**2 + 8*g0*g2s**2

     K2s = 2*g1c**2*g2s

     K2c = -2*g1c**2*g20 + 2*g1c**2*g2c

     K4s = g1c**2*g2s - 16*g0*g2c*g2s

     K4c = g1c**2*g2c - 8*g0*g2c**2 + 8*g0*g2s**2

     coefficients(1) = 4*(K4c**2 + K4s**2)

     coefficients(2) = 4*(K4s*K2c - K2s*K4c)

     coefficients(3) = K2s**2 + K2c**2 - 4*K0*K4c - 4*K4c**2 - 4*K4s**2

     coefficients(4) = 2*K0*K2s + 2*K4c*K2s - 4*K4s*K2c

     coefficients(5) = (K0 + K4c)**2 - K2c**2

     call quasisymmetry_quartic_roots(coefficients, real_parts, imag_parts)

     if (local_verbose) then
        print *,"g0:",g0,"  g1c:",g1c
        print *,"g20:",g20,"  g2s:",g2s,"  g2c:",g2c
        print *,"K0:",K0,"  K2s:",K2s,"  K2c:",K2c
        print *,"K4s:",K4s,"  K4c:",K4c
        print *,"coefficients:",coefficients
        print *,"real_parts:",real_parts
        print *,"imag_parts:",imag_parts
     end if

     rc = 1d+100 ! This huge number indicates a true solution has not yet been found.

     do jr = 1, 4 ! Loop over the roots of the equation for w.

        ! If root is not purely real, skip it.
        if (abs(imag_parts(jr)) > 1e-7) then
           if (local_verbose) print *,"Skipping root with jr=",jr," since imag part is",imag_parts(jr)
           cycle
        end if

        sin2theta = real_parts(jr)

        ! Discard any roots that have magnitude larger than 1. (I'm not sure this ever happens, but check to be sure.)
        if (abs(sin2theta) > 1) then
           if (local_verbose) print *,"Skipping root with jr=",jr," since sin2theta=",sin2theta
           cycle
        end if

        sigma_denominator = ((K4s*2*sin2theta + K2c) * sqrt(1 - sin2theta*sin2theta))
        if (abs(sigma_denominator) < 1e-8) print *,"WARNING!!! sigma_denominator=",sigma_denominator
        sigma = -(K0 + K2s * sin2theta + K4c*(1 - 2*sin2theta*sin2theta)) / sigma_denominator
        if (abs(sigma*sigma-1) > 1e-3) then
           print *,"WARNING!!! abs(sigma*sigma-1) =",abs(sigma*sigma-1)
        end if
        if (abs(sigma*sigma-1) > 1e-1) then
           print "(3(a,es24.15),2(a,3(es24.15)))","WARNING!!! abs(sigma*sigma-1) =",abs(sigma*sigma-1)," etabar=",eta_bar," B2c=",B2c," R0c=",R0c(2:4)," Z0s=",Z0s(2:4)
        end if
        sigma = nint(sigma) ! Ensure sigma is exactly either +1 or -1.
        cos2theta = sigma * sqrt(1 - sin2theta*sin2theta)
        abs_costheta = sqrt(0.5*(1 + cos2theta))
        if (local_verbose) print *,"  jr=",jr,"  sin2theta=",sin2theta,"  cos2theta=",cos2theta
        do varsigma = -1, 1, 2 ! so varsigma will be either -1 or +1.
           costheta = varsigma * abs_costheta
           sintheta = sin2theta / (2 * costheta)
           if (local_verbose) print *,"    varsigma=",varsigma,"  costheta=",costheta,"  sintheta=",sintheta

           ! Sanity test
           if (abs(costheta*costheta + sintheta*sintheta - 1) > 1e-4) then
              print *,"Error! sintheta=",sintheta,"  costheta=",costheta
              print *,"j=",j,"  jr=",jr,"  sin2theta=",sin2theta,"  cos2theta=",cos2theta,"  sigma_denominator=",sigma_denominator
              print *,"abs(costheta*costheta + sintheta*sintheta - 1):",abs(costheta*costheta + sintheta*sintheta - 1)
              if (trim(general_option)==general_option_single) stop
           end if

           ! Try to get r using the simpler method, the equation that is linear in r.
           denominator = 2*(g2s*cos2theta - g2c*sin2theta)
           if (abs(denominator) > 1e-8) then ! This method cannot be used if we would need to divide by 0
              r = g1c*sintheta / denominator
              residual = g0 + r*g1c*costheta + r*r*(g20 + g2s*sin2theta + g2c*cos2theta) ! Residual in the equation sqrt(g)=0.
              if (local_verbose) print *,"    Linear method: r=",r,"  residual=",residual
              if ((r>0) .and. (abs(residual) < 1e-5)) then
                 rc = min(rc, r)
                 if (local_verbose) print *,"      > rc =",rc
              end if
           else
              ! Use the more complicated method to determine r by solving a quadratic equation.
              quadratic_A = g20 + g2s * sin2theta + g2c * cos2theta
              quadratic_B = costheta * g1c
              quadratic_C = g0
              radical = sqrt(quadratic_B * quadratic_B - 4 * quadratic_A * quadratic_C)
              do sign_quadratic = -1, 1, 2 ! So sign_quadratic = +1 or -1
                 r = (-quadratic_B + sign_quadratic * radical) / (2 * quadratic_A) ! This is the quadratic formula.
                 residual = -g1c*sintheta + 2*r*(g2s*cos2theta - g2c*sin2theta) ! Residual in the equation d sqrt(g) / d theta = 0.
                 if (local_verbose) print *,"    Quadratic method: r=",r,"  residual=",residual
                 if ((r>0) .and. (abs(residual) < 1e-5)) then
                    rc = min(rc, r)
                    if (local_verbose) print *,"      > rc =",rc
                 end if
              end do
           end if
        end do
     end do
     r_singularity_vs_zeta(j) = rc

  end do
  r_singularity = minval(r_singularity_vs_zeta)

  call cpu_time(end_time)
  if (verbose) print "(a,es11.4,a,es10.3,a)"," r_singularity:",r_singularity,"  Time to compute:",end_time - start_time," sec."

end subroutine quasisymmetry_max_r_before_singularity

! ---------------------------------------------------

!> Compute the roots of a quartic polynomial by computing eigenvalues of a companion matrix
!>
!> This subroutine follows the same algorithm as Matlab's 'roots' routine.
!> @params coefficients Ordered the same way as in matlab, from the coefficient of x^4 to the coefficient of x^0.
!> @params real_parts Real parts of the roots
!> @params imag_parts Imaginary parts of the roots
subroutine quasisymmetry_quartic_roots(coefficients, real_parts, imag_parts)

  use quasisymmetry_variables, only: dp

  implicit none
  real(dp), intent(in) :: coefficients(5)
  real(dp), intent(out) :: real_parts(4), imag_parts(4)

  real(dp) :: matrix(4,4)
  integer, parameter :: LWORK = 100 ! Work array for LAPACK
  real(dp) :: WORK(LWORK), VL(1), VR(1)
  integer :: INFO

  matrix(1,1) = -coefficients(2) / coefficients(1)
  matrix(1,2) = -coefficients(3) / coefficients(1)
  matrix(1,3) = -coefficients(4) / coefficients(1)
  matrix(1,4) = -coefficients(5) / coefficients(1)

  matrix(2:4,:) = 0.0_dp

  matrix(2,1) = 1.0_dp
  matrix(3,2) = 1.0_dp
  matrix(4,3) = 1.0_dp

  call dgeev('N', 'N', 4, matrix, 4, real_parts, imag_parts, VL, 1, VR, 1, WORK, LWORK, INFO)
  if (INFO .ne. 0) then
     print *,"Error in DGEEV: info=",INFO
     stop
  end if

end subroutine quasisymmetry_quartic_roots

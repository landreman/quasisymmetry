subroutine quasisymmetry_validate_input

  use quasisymmetry_variables

  implicit none

  integer :: j

  if (nfp<1) stop "Error! nfp must be positive."

  if (sign_G .ne. 1 .and. sign_G .ne. -1) stop "sign_G must be +1 or -1."

  if (N_phi < 0) stop "Error! N_phi must be positive"
  
  ! Ensure N_phi is always odd.
  if (mod(N_phi,2) ==0) N_phi = N_phi + 1

  N_phi_original = N_phi

  select case (trim(resolution_option))
  case (resolution_option_fixed)
  case (resolution_option_adaptive)
  case default
     print *,"Error! Invalid resolution_option:",resolution_option
     stop
  end select

  select case (trim(general_option))
  case (general_option_single)
  case (general_option_scan)
  case (general_option_random)
  case default
     print *,"Error! Invalid general_option:",general_option
     stop
  end select

  select case (trim(verbose_option))
  case (verbose_option_all)
  case (verbose_option_proc0)
  case (verbose_option_summary)
  case default
     print *,"Error! Invalid verbose_option:",verbose_option
     stop
  end select

  select case (trim(eta_bar_scan_option))
  case (eta_bar_scan_option_linear)
  case (eta_bar_scan_option_log)
  case (eta_bar_scan_option_2_sided_log)
  case default
     print *,"Error! Invalid eta_bar_scan_option:",eta_bar_scan_option
     stop
  end select

  select case (trim(sigma_initial_scan_option))
  case (sigma_initial_scan_option_linear)
  case (sigma_initial_scan_option_log)
  case (sigma_initial_scan_option_2_sided_log)
  case default
     print *,"Error! Invalid sigma_initial_scan_option:",sigma_initial_scan_option
     stop
  end select

  select case (trim(B2s_scan_option))
  case (B2s_scan_option_linear)
  case (B2s_scan_option_log)
  case (B2s_scan_option_2_sided_log)
  case default
     print *,"Error! Invalid B2s_scan_option:",B2s_scan_option
     stop
  end select

  select case (trim(B2c_scan_option))
  case (B2c_scan_option_linear)
  case (B2c_scan_option_log)
  case (B2c_scan_option_2_sided_log)
  case default
     print *,"Error! Invalid B2c_scan_option:",B2c_scan_option
     stop
  end select

  select case (trim(Fourier_scan_option))
  case (Fourier_scan_option_linear)
  case (Fourier_scan_option_2_sided_log)
  case (Fourier_scan_option_2_sided_log_except_Z0s1)
  case default
     print *,"Error! Invalid Fourier_scan_option:",Fourier_scan_option
     stop
  end select

  select case (trim(finite_r_option))
  case (finite_r_option_linear)
  case (finite_r_option_nonlinear)
  case default
     print *,"Error! Invalid finite_r_option:",finite_r_option
     stop
  end select

  select case (trim(order_r_option))
  case (order_r_option_r1)
  case (order_r_option_r1_compute_B2)
  case (order_r_option_r2)
  case (order_r_option_r3_simplified)
  case (order_r_option_r3_simplified_with_Z3)
  case (order_r_option_r3_flux_constraint)
  case (order_r_option_r3_flux_constraint_const_B20)
  case (order_r_option_r3_B3)
  case (order_r_option_r3_X3s3_X3c3)
  !case (order_r_option_r3_X3s3_Y3s3)
  !case (order_r_option_r3_X3c3_Y3c3)
  !case (order_r_option_r3_Y3s3_Y3c3)
  case default
     print *,"Error! Invalid order_r_option:",order_r_option
     stop
  end select

!!$  if (order_r_squared .and. trim(finite_r_option)==finite_r_option_linear) then
!!$     print "(a)"," NOTE: Since order_r_squared==.true., finite_r_option is being set to 'nonlinear'."
!!$     finite_r_option = finite_r_option_nonlinear
!!$  end if

  ! Ensuring that "max" values are not less than "min" values simplifies the logic in quasisymmetry_random()
  if (eta_bar_max < eta_bar_min) stop "eta_bar_max < eta_bar_min"
  if (sigma_initial_max < sigma_initial_min) stop "sigma_initial_max < sigma_initial_min"
  if (B2s_max < B2s_min) stop "B2s_max < B2s_min"
  if (B2c_max < B2c_min) stop "B2c_max < B2c_min"
  do j = 1, axis_nmax+1
     if (R0s_max(j) < R0s_min(j)) stop "An entry of R0s_max is < the corresponding entry in R0s_min"
     if (R0c_max(j) < R0c_min(j)) stop "An entry of R0c_max is < the corresponding entry in R0c_min"
     if (Z0s_max(j) < Z0s_min(j)) stop "An entry of Z0s_max is < the corresponding entry in Z0s_min"
     if (Z0c_max(j) < Z0c_min(j)) stop "An entry of Z0c_max is < the corresponding entry in Z0c_min"
  end do

  if (min_R0_to_keep <= 0) stop "min_R0_to_keep should be positive."

end subroutine quasisymmetry_validate_input

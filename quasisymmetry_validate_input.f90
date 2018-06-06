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
  end select

  select case (trim(general_option))
  case (general_option_single)
  case (general_option_scan)
  case default
     print *,"Error! Invalid general_option:",general_option
  end select

  select case (trim(constraint_option))
  case (constraint_option_no_Z_component)
  case (constraint_option_sigma_initial)
  case default
     print *,"Error! Invalid constraint_option:",constraint_option
  end select

end subroutine quasisymmetry_validate_input

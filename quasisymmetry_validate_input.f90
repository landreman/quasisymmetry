subroutine quasisymmetry_validate_input

  use quasisymmetry_variables

  implicit none

  integer :: j

  if (nfp<1) stop "Error! nfp must be positive."

  if (sign_G .ne. 1 .and. sign_G .ne. -1) stop "sign_G must be +1 or -1."

  ! Ensure N_phi is always odd.
  do j = 1,N_N_phi
     if (N_phis(j) < 0) stop "Error! N_phi must be positive"
     if (mod(N_phis(j),2) ==0) N_phis(j) = N_phis(j) + 1
  end do

end subroutine quasisymmetry_validate_input

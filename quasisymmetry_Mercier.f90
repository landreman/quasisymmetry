subroutine quasisymmetry_Mercier

  use quasisymmetry_variables

  implicit none

  real(dp), allocatable :: integrand(:)
  real(dp) :: integral, iota_N0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  iota_N0 = iota + axis_helicity*nfp

  allocate(integrand(N_phi))

  ! See Overleaf note "Mercier criterion near the magnetic axis- detailed notes".
  ! See also "20200604-02 Checking sign in Mercier DGeod near axis.docx"

  !integrand = d_l_d_phi * (Y1c * Y1c + X1c * (X1c + Y1s)) / (Y1c * Y1c + (X1c + Y1s) * (X1c + Y1s))
  integrand = d_l_d_phi * (eta_bar*eta_bar*eta_bar*eta_bar + curvature*curvature*curvature*curvature*sigma*sigma + eta_bar*eta_bar*curvature*curvature) &
       / (eta_bar*eta_bar*eta_bar*eta_bar + curvature*curvature*curvature*curvature*(1+sigma*sigma) + 2*eta_bar*eta_bar*curvature*curvature)

  integral = sum(integrand) * d_phi * nfp * 2 * pi / axis_length

  !DGeod_times_r2 = -(2 * sign_G * sign_psi * mu0 * mu0 * p2 * p2 * G0 * G0 * G0 * G0 * eta_bar * eta_bar &
  DGeod_times_r2 = -(2 * mu0 * mu0 * p2 * p2 * G0 * G0 * G0 * G0 * eta_bar * eta_bar &
       / (pi * pi * pi * B0 * B0 * B0 * B0 * B0 * B0 * B0 * B0 * B0 * B0 * iota_N0 * iota_N0)) &
       * integral

  if (trim(order_r_option) .ne. order_r_option_r1 .and. trim(order_r_option) .ne. order_r_option_r1_compute_B2) then

     d2_volume_d_psi2 = 4*pi*pi*abs(G0)/(B0*B0*B0)*(3*eta_bar*eta_bar - 4*B20_mean/B0 + 2*(G2+iota*I2)/G0)

     DWell_times_r2 = (mu0 * p2 * abs(G0) / (8 * pi * pi * pi * pi * B0 * B0 * B0)) * (d2_volume_d_psi2 - 8 * pi * pi * mu0 * p2 * abs(G0) / (B0 * B0 * B0 * B0 * B0))

     DMerc_times_r2 = DWell_times_r2 + DGeod_times_r2

  end if

end subroutine quasisymmetry_Mercier

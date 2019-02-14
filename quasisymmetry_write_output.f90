subroutine quasisymmetry_write_output

  use quasisymmetry_variables
  use ezcdf
  use vmec_input, only: lasym, ntor, RBC, RBS, ZBC, ZBS

  implicit none

  real :: time1, end_time
  integer :: ierr, ncid

  ! Same convention as in VMEC:
  ! Prefix vn_ indicates the variable name used in the .nc file.

  ! Scalars:
  character(len=*), parameter :: &
       vn_general_option = "general_option", &
       vn_nfp = "nfp", &
       vn_sign_G = "sign_G", &
       vn_sign_psi = "sign_psi", &
       vn_resolution_option = "resolution_option", &
       vn_max_precise_elongation = "max_precise_elongation", &
       vn_sigma_initial_min = "sigma_initial_min", &
       vn_sigma_initial_max = "sigma_initial_max", &
       vn_sigma_initial_N_scan = "sigma_initial_N_scan", &
       vn_sigma_initial_scan_option = "sigma_initial_scan_option", &
       vn_eta_bar_min = "eta_bar_min", &
       vn_eta_bar_max = "eta_bar_max", &
       vn_eta_bar_N_scan = "eta_bar_N_scan", &
       vn_eta_bar_scan_option = "eta_bar_scan_option", &
       vn_Fourier_scan_option = "Fourier_scan_option", &
       vn_max_elongation_to_keep = "max_elongation_to_keep", &
       vn_max_max_curvature_to_keep = "max_max_curvature_to_keep", &
       vn_min_iota_to_keep = "min_iota_to_keep", &
       vn_eta_bar = "eta_bar", &
       vn_sigma_initial = "sigma_initial", &
       vn_I2_over_B0 = "I2_over_B0", &
       vn_iota = "iota", &
       vn_max_elongation = "max_elongation", &
       vn_rms_curvature = "rms_curvature", &
       vn_max_curvature = "max_curvature", &
       vn_axis_length = "axis_length", &
       vn_standard_deviation_of_R = "standard_deviation_of_R", &
       vn_standard_deviation_of_Z = "standard_deviation_of_Z", &
       vn_axis_helicity = "axis_helicity", &
       vn_B_helicity = "B_helicity", &
       vn_effective_nfp = "effective_nfp", &
       vn_Newton_tolerance_achieved = "Newton_tolerance_achieved", &
       vn_iota_tolerance_achieved = "iota_tolerance_achieved", &
       vn_elongation_tolerance_achieved = "elongation_tolerance_achieved", &
       vn_finite_r_option = "finite_r_option", &
       vn_N_phi = "N_phi", &
       vn_r = "r", &
       vn_mpol = "mpol", &
       vn_ntor = "ntor", &
       vn_lasym = "lasym"

  ! Arrays with dimension 1
  character(len=*), parameter :: &
       vn_iotas = "iotas", &
       vn_max_elongations = "max_elongations", &
       vn_rms_curvatures = "rms_curvatures", &
       vn_max_curvatures = "max_curvatures", &
       vn_axis_lengths = "axis_lengths", &
       vn_standard_deviations_of_R = "standard_deviations_of_R", &
       vn_standard_deviations_of_Z = "standard_deviations_of_Z", &
       vn_axis_helicities = "axis_helicities", &
       vn_B_helicities = "B_helicities", &
       vn_effective_nfps = "effective_nfps", &
       vn_Newton_tolerance_achieveds = "Newton_tolerance_achieveds", &
       vn_iota_tolerance_achieveds = "iota_tolerance_achieveds", &
       vn_elongation_tolerance_achieveds = "elongation_tolerance_achieveds", &
       vn_R0s_min = "R0s_min", &
       vn_R0s_max = "R0s_max", &
       vn_R0s_N_scan = "R0s_N_scan", &
       vn_R0c_min = "R0c_min", &
       vn_R0c_max = "R0c_max", &
       vn_R0c_N_scan = "R0c_N_scan", &
       vn_Z0s_min = "Z0s_min", &
       vn_Z0s_max = "Z0s_max", &
       vn_Z0s_N_scan = "Z0s_N_scan", &
       vn_Z0c_min = "Z0c_min", &
       vn_Z0c_max = "Z0c_max", &
       vn_Z0c_N_scan = "Z0c_N_scan", &
       vn_scan_sigma_initial = "scan_sigma_initial", &
       vn_sigma_initial_values = "sigma_initial_values", &
       vn_scan_eta_bar = "scan_eta_bar", &
       vn_eta_bar_values = "eta_bar_values", &
       vn_R0c = "R0c", &
       vn_R0s = "R0s", &
       vn_Z0c = "Z0c", &
       vn_Z0s = "Z0s", &
       vn_phi = "phi", &
       vn_Boozer_toroidal_angle = "Boozer_toriodal_angle", &
       vn_R0 = "R0", &
       vn_z0 = "z0", &
       vn_curvature = "curvature", &
       vn_torsion = "torsion", &
       vn_sigma = "sigma", &
       vn_X1c = "X1c", &
       vn_Y1c = "Y1c", &
       vn_Y1s = "Y1s", &
       vn_X1c_untwisted = "X1c_untwisted", &
       vn_X1s_untwisted = "X1s_untwisted", &
       vn_Y1c_untwisted = "Y1c_untwisted", &
       vn_Y1s_untwisted = "Y1s_untwisted", &
       vn_elongation = "elongation"

  ! Arrays with dimension 2
  character(len=*), parameter :: &
       vn_N_scan_array  = "N_scan_array", &
       vn_scan_R0c  = "scan_R0c", &
       vn_scan_R0s  = "scan_R0s", &
       vn_scan_Z0c  = "scan_Z0c", &
       vn_scan_Z0s  = "scan_Z0s", &
       vn_RBC = "RBC", &
       vn_RBS = "RBS", &
       vn_ZBC = "ZBC", &
       vn_ZBS = "ZBS"

!!$  ! Arrays with dimension 3
!!$  character(len=*), parameter :: &
!!$       vn_r_plasma  = "r_plasma"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now create variables that name the dimensions.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Arrays with dimension 1:
  character(len=*), parameter, dimension(1) :: &
       N_scan_dim = (/'N_scan'/), &
       axis_nmax_plus_1_dim = (/'axis_nmax_plus_1'/), &
       sigma_initial_N_scan_dim = (/'sigma_initial_N_scan'/), &
       eta_bar_N_scan_dim = (/'eta_bar_N_scan'/), &
       N_phi_dim = (/'N_phi'/)

  ! Arrays with dimension 2:
  ! The form of the array declarations here is inspired by
  ! http://stackoverflow.com/questions/21552430/gfortran-does-not-allow-character-arrays-with-varying-component-lengths
  character(len=*), parameter, dimension(2) :: &
       axis_nmax_plus_1_4_dim = (/ character(len=50) :: 'axis_nmax_plus_1','4'/), &
       N_scan_axis_nmax_plus_1_dim = (/ character(len=50) :: 'N_scan','axis_nmax_plus_1'/), &
       ntor_mpol_dim = (/ character(len=50) :: 'ntor_times_2_plus_1','mpol_plus_1' /)

!!$  ! Arrays with dimension 3:
!!$  character(len=*), parameter, dimension(3) :: &
!!$       xyz_ntheta_nzetal_plasma_dim = (/ character(len=50) :: 'xyz','ntheta_plasma','nzetal_plasma'/), &

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Only proc 0 writes.
  if (.not. proc0) return

  call cpu_time(time1)

  call cdf_open(ncid,output_filename,'w',ierr)
  IF (ierr .ne. 0) then
     print *,"Error opening output file ",output_filename
     stop
  end IF

  ! Scalars

  call cdf_define(ncid, vn_general_option, general_option)

  call cdf_define(ncid, vn_nfp, nfp)
  call cdf_setatt(ncid, vn_nfp, 'Number of field periods, i.e. the number of identical toroidal segments, 5 for W7-X, 4 for HSX, etc. ' // &
       'Equivalent to the VMEC variable of the same name.')

  call cdf_define(ncid, vn_sign_G, sign_G)
  call cdf_define(ncid, vn_sign_psi, sign_psi)

  call cdf_define(ncid, vn_resolution_option, resolution_option)
  !call cdf_setatt(ncid, vn_resolution_option, 'Method used to define the geometry of the plasma surface.' // input_parameter_text)

  call cdf_define(ncid, vn_max_precise_elongation, max_precise_elongation)

  select case (trim(general_option))
  case (general_option_single)
     call cdf_define(ncid, vn_eta_bar, eta_bar)
     call cdf_define(ncid, vn_sigma_initial, sigma_initial)
     call cdf_define(ncid, vn_I2_over_B0, I2_over_B0)
     call cdf_define(ncid, vn_iota, iota)
     call cdf_define(ncid, vn_max_elongation, max_elongation)
     call cdf_define(ncid, vn_rms_curvature, rms_curvature)
     call cdf_define(ncid, vn_max_curvature, max_curvature)
     call cdf_define(ncid, vn_axis_length, axis_length)
     call cdf_define(ncid, vn_standard_deviation_of_R, standard_deviation_of_R)
     call cdf_define(ncid, vn_standard_deviation_of_Z, standard_deviation_of_Z)
     call cdf_define(ncid, vn_axis_helicity, axis_helicity)
     call cdf_define(ncid, vn_B_helicity, B_helicity)
     call cdf_define(ncid, vn_effective_nfp, effective_nfp)
     call cdf_define(ncid, vn_Newton_tolerance_achieved, Newton_tolerance_achieved)
     call cdf_define(ncid, vn_iota_tolerance_achieved, iota_tolerance_achieved)
     call cdf_define(ncid, vn_elongation_tolerance_achieved, elongation_tolerance_achieved)
     call cdf_define(ncid, vn_finite_r_option, finite_r_option)
     call cdf_define(ncid, vn_N_phi, N_phi)
     call cdf_define(ncid, vn_r, r)
     call cdf_define(ncid, vn_mpol, mpol_nonzero)
     call cdf_define(ncid, vn_ntor, ntor)
     call cdf_define(ncid, vn_lasym, lasym)
  case (general_option_scan)
     call cdf_define(ncid, vn_sigma_initial_min, sigma_initial_min)
     call cdf_define(ncid, vn_sigma_initial_max, sigma_initial_max)
     call cdf_define(ncid, vn_sigma_initial_N_scan, sigma_initial_N_scan)
     call cdf_define(ncid, vn_sigma_initial_scan_option, sigma_initial_scan_option)
     call cdf_define(ncid, vn_eta_bar_min, eta_bar_min)
     call cdf_define(ncid, vn_eta_bar_max, eta_bar_max)
     call cdf_define(ncid, vn_eta_bar_N_scan, eta_bar_N_scan)
     call cdf_define(ncid, vn_eta_bar_scan_option, eta_bar_scan_option)
     call cdf_define(ncid, vn_Fourier_scan_option, Fourier_scan_option)
     call cdf_define(ncid, vn_max_elongation_to_keep, max_elongation_to_keep)
     call cdf_define(ncid, vn_max_max_curvature_to_keep, max_max_curvature_to_keep)
     call cdf_define(ncid, vn_min_iota_to_keep, min_iota_to_keep)
  end select

  ! Arrays with dimension 1

  select case (trim(general_option))
  case (general_option_single)
     call cdf_define(ncid, vn_R0c, R0c(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_R0s, R0s(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0c, Z0c(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0s, Z0s(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_phi, phi, dimname=N_phi_dim)
     call cdf_define(ncid, vn_Boozer_toroidal_angle, Boozer_toroidal_angle, dimname=N_phi_dim)
     call cdf_define(ncid, vn_R0, R0, dimname=N_phi_dim)
     call cdf_define(ncid, vn_z0, z0, dimname=N_phi_dim)
     call cdf_define(ncid, vn_curvature, curvature, dimname=N_phi_dim)
     call cdf_define(ncid, vn_torsion, torsion, dimname=N_phi_dim)
     call cdf_define(ncid, vn_sigma, sigma, dimname=N_phi_dim)
     call cdf_define(ncid, vn_X1c, X1c, dimname=N_phi_dim)
     call cdf_define(ncid, vn_Y1c, Y1c, dimname=N_phi_dim)
     call cdf_define(ncid, vn_Y1s, Y1s, dimname=N_phi_dim)
     call cdf_define(ncid, vn_X1c_untwisted, X1c_untwisted, dimname=N_phi_dim)
     call cdf_define(ncid, vn_X1s_untwisted, X1s_untwisted, dimname=N_phi_dim)
     call cdf_define(ncid, vn_Y1c_untwisted, Y1c_untwisted, dimname=N_phi_dim)
     call cdf_define(ncid, vn_Y1s_untwisted, Y1s_untwisted, dimname=N_phi_dim)
     call cdf_define(ncid, vn_elongation, elongation, dimname=N_phi_dim)
  case (general_option_scan)
     call cdf_define(ncid, vn_iotas, iotas, dimname=N_scan_dim)
     call cdf_define(ncid, vn_max_elongations, max_elongations, dimname=N_scan_dim)
     call cdf_define(ncid, vn_rms_curvatures, rms_curvatures, dimname=N_scan_dim)
     call cdf_define(ncid, vn_max_curvatures, max_curvatures, dimname=N_scan_dim)
     call cdf_define(ncid, vn_axis_lengths, axis_lengths, dimname=N_scan_dim)
     call cdf_define(ncid, vn_standard_deviations_of_R, standard_deviations_of_R, dimname=N_scan_dim)
     call cdf_define(ncid, vn_standard_deviations_of_Z, standard_deviations_of_Z, dimname=N_scan_dim)
     call cdf_define(ncid, vn_axis_helicities, axis_helicities, dimname=N_scan_dim)
     call cdf_define(ncid, vn_B_helicities, B_helicities, dimname=N_scan_dim)
     call cdf_define(ncid, vn_effective_nfps, effective_nfps, dimname=N_scan_dim)
     call cdf_define(ncid, vn_Newton_tolerance_achieveds, Newton_tolerance_achieveds, dimname=N_scan_dim)
     call cdf_define(ncid, vn_iota_tolerance_achieveds, iota_tolerance_achieveds, dimname=N_scan_dim)
     call cdf_define(ncid, vn_elongation_tolerance_achieveds, elongation_tolerance_achieveds, dimname=N_scan_dim)
     call cdf_define(ncid, vn_R0s_min, R0s_min(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_R0s_max, R0s_max(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_R0s_N_scan, R0s_N_scan(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_R0c_min, R0c_min(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_R0c_max, R0c_max(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_R0c_N_scan, R0c_N_scan(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0s_min, Z0s_min(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0s_max, Z0s_max(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0s_N_scan, Z0s_N_scan(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0c_min, Z0c_min(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0c_max, Z0c_max(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_Z0c_N_scan, Z0c_N_scan(1:axis_nmax+1), dimname=axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_scan_sigma_initial, scan_sigma_initial, dimname=N_scan_dim)
     call cdf_define(ncid, vn_sigma_initial_values, sigma_initial_values, dimname=sigma_initial_N_scan_dim)
     call cdf_define(ncid, vn_scan_eta_bar, scan_eta_bar, dimname=N_scan_dim)
     call cdf_define(ncid, vn_eta_bar_values, eta_bar_values, dimname=eta_bar_N_scan_dim)
  end select

  ! Arrays with dimension 2

  select case (trim(general_option))
  case (general_option_single)
     call cdf_define(ncid, vn_RBC, RBC(-ntor:ntor, 0:mpol_nonzero), dimname=ntor_mpol_dim)
     call cdf_define(ncid, vn_RBS, RBS(-ntor:ntor, 0:mpol_nonzero), dimname=ntor_mpol_dim)
     call cdf_define(ncid, vn_ZBC, ZBC(-ntor:ntor, 0:mpol_nonzero), dimname=ntor_mpol_dim)
     call cdf_define(ncid, vn_ZBS, ZBS(-ntor:ntor, 0:mpol_nonzero), dimname=ntor_mpol_dim)
  case (general_option_scan)
     call cdf_define(ncid, vn_N_scan_array,  N_scan_array(1:axis_nmax+1,:), dimname=axis_nmax_plus_1_4_dim)
     call cdf_define(ncid, vn_scan_R0c,  scan_R0c, dimname=N_scan_axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_scan_R0s,  scan_R0s, dimname=N_scan_axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_scan_Z0c,  scan_Z0c, dimname=N_scan_axis_nmax_plus_1_dim)
     call cdf_define(ncid, vn_scan_Z0s,  scan_Z0s, dimname=N_scan_axis_nmax_plus_1_dim)
  end select

  ! Arrays with dimension 3

  !call cdf_define(ncid, vn_r_plasma,  r_plasma,  dimname=xyz_ntheta_nzetal_plasma_dim)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! Done with cdf_define calls. Now write the data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  ! Scalars

  call cdf_write(ncid, vn_general_option, general_option)
  call cdf_write(ncid, vn_nfp, nfp)
  call cdf_write(ncid, vn_sign_G, sign_G)
  call cdf_write(ncid, vn_sign_psi, sign_psi)
  call cdf_write(ncid, vn_resolution_option, resolution_option)
  call cdf_write(ncid, vn_max_precise_elongation, max_precise_elongation)

  select case (trim(general_option))
  case (general_option_single)
     call cdf_write(ncid, vn_eta_bar, eta_bar)
     call cdf_write(ncid, vn_sigma_initial, sigma_initial)
     call cdf_write(ncid, vn_I2_over_B0, I2_over_B0)
     call cdf_write(ncid, vn_iota, iota)
     call cdf_write(ncid, vn_max_elongation, max_elongation)
     call cdf_write(ncid, vn_rms_curvature, rms_curvature)
     call cdf_write(ncid, vn_max_curvature, max_curvature)
     call cdf_write(ncid, vn_axis_length, axis_length)
     call cdf_write(ncid, vn_standard_deviation_of_R, standard_deviation_of_R)
     call cdf_write(ncid, vn_standard_deviation_of_Z, standard_deviation_of_Z)
     call cdf_write(ncid, vn_axis_helicity, axis_helicity)
     call cdf_write(ncid, vn_B_helicity, B_helicity)
     call cdf_write(ncid, vn_effective_nfp, effective_nfp)
     call cdf_write(ncid, vn_Newton_tolerance_achieved, Newton_tolerance_achieved)
     call cdf_write(ncid, vn_iota_tolerance_achieved, iota_tolerance_achieved)
     call cdf_write(ncid, vn_elongation_tolerance_achieved, elongation_tolerance_achieved)
     call cdf_write(ncid, vn_finite_r_option, finite_r_option)
     call cdf_write(ncid, vn_N_phi, N_phi)
     call cdf_write(ncid, vn_r, r)
     call cdf_write(ncid, vn_mpol, mpol_nonzero)
     call cdf_write(ncid, vn_ntor, ntor)
     call cdf_write(ncid, vn_lasym, lasym)
  case (general_option_scan)
     call cdf_write(ncid, vn_sigma_initial_min, sigma_initial_min)
     call cdf_write(ncid, vn_sigma_initial_max, sigma_initial_max)
     call cdf_write(ncid, vn_sigma_initial_N_scan, sigma_initial_N_scan)
     call cdf_write(ncid, vn_eta_bar_min, eta_bar_min)
     call cdf_write(ncid, vn_eta_bar_max, eta_bar_max)
     call cdf_write(ncid, vn_eta_bar_N_scan, eta_bar_N_scan)
     call cdf_write(ncid, vn_eta_bar_scan_option, eta_bar_scan_option)
     call cdf_write(ncid, vn_Fourier_scan_option, Fourier_scan_option)
     call cdf_write(ncid, vn_sigma_initial_scan_option, sigma_initial_scan_option)
     call cdf_write(ncid, vn_max_elongation_to_keep, max_elongation_to_keep)
     call cdf_write(ncid, vn_max_max_curvature_to_keep, max_max_curvature_to_keep)
     call cdf_write(ncid, vn_min_iota_to_keep, min_iota_to_keep)
  end select

  ! Arrays with dimension 1

  select case (trim(general_option))
  case (general_option_single)
     call cdf_write(ncid, vn_R0c, R0c(1:axis_nmax+1))
     call cdf_write(ncid, vn_R0s, R0s(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0c, Z0c(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0s, Z0s(1:axis_nmax+1))
     call cdf_write(ncid, vn_phi, phi)
     call cdf_write(ncid, vn_Boozer_toroidal_angle, Boozer_toroidal_angle)
     call cdf_write(ncid, vn_R0, R0)
     call cdf_write(ncid, vn_z0, z0)
     call cdf_write(ncid, vn_curvature, curvature)
     call cdf_write(ncid, vn_torsion, torsion)
     call cdf_write(ncid, vn_sigma, sigma)
     call cdf_write(ncid, vn_X1c, X1c)
     call cdf_write(ncid, vn_Y1c, Y1c)
     call cdf_write(ncid, vn_Y1s, Y1s)
     call cdf_write(ncid, vn_X1c_untwisted, X1c_untwisted)
     call cdf_write(ncid, vn_X1s_untwisted, X1s_untwisted)
     call cdf_write(ncid, vn_Y1c_untwisted, Y1c_untwisted)
     call cdf_write(ncid, vn_Y1s_untwisted, Y1s_untwisted)
     call cdf_write(ncid, vn_elongation, elongation)
  case (general_option_scan)
     call cdf_write(ncid, vn_iotas, iotas)
     call cdf_write(ncid, vn_max_elongations, max_elongations)
     call cdf_write(ncid, vn_rms_curvatures, rms_curvatures)
     call cdf_write(ncid, vn_max_curvatures, max_curvatures)
     call cdf_write(ncid, vn_axis_lengths, axis_lengths)
     call cdf_write(ncid, vn_standard_deviations_of_R, standard_deviations_of_R)
     call cdf_write(ncid, vn_standard_deviations_of_Z, standard_deviations_of_Z)
     call cdf_write(ncid, vn_axis_helicities, axis_helicities)
     call cdf_write(ncid, vn_B_helicities, B_helicities)
     call cdf_write(ncid, vn_effective_nfps, effective_nfps)
     call cdf_write(ncid, vn_Newton_tolerance_achieveds, Newton_tolerance_achieveds)
     call cdf_write(ncid, vn_iota_tolerance_achieveds, iota_tolerance_achieveds)
     call cdf_write(ncid, vn_elongation_tolerance_achieveds, elongation_tolerance_achieveds)
     call cdf_write(ncid, vn_R0s_min, R0s_min(1:axis_nmax+1))
     call cdf_write(ncid, vn_R0s_max, R0s_max(1:axis_nmax+1))
     call cdf_write(ncid, vn_R0s_N_scan, R0s_N_scan(1:axis_nmax+1))
     call cdf_write(ncid, vn_R0c_min, R0c_min(1:axis_nmax+1))
     call cdf_write(ncid, vn_R0c_max, R0c_max(1:axis_nmax+1))
     call cdf_write(ncid, vn_R0c_N_scan, R0c_N_scan(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0s_min, Z0s_min(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0s_max, Z0s_max(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0s_N_scan, Z0s_N_scan(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0c_min, Z0c_min(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0c_max, Z0c_max(1:axis_nmax+1))
     call cdf_write(ncid, vn_Z0c_N_scan, Z0c_N_scan(1:axis_nmax+1))
     call cdf_write(ncid, vn_scan_sigma_initial, scan_sigma_initial)
     call cdf_write(ncid, vn_sigma_initial_values, sigma_initial_values)
     call cdf_write(ncid, vn_scan_eta_bar, scan_eta_bar)
     call cdf_write(ncid, vn_eta_bar_values, eta_bar_values)
  end select

  ! Arrays with dimension 2

  select case (trim(general_option))
  case (general_option_single)
     call cdf_write(ncid, vn_RBC, RBC(-ntor:ntor, 0:mpol_nonzero))
     call cdf_write(ncid, vn_RBS, RBS(-ntor:ntor, 0:mpol_nonzero))
     call cdf_write(ncid, vn_ZBC, ZBC(-ntor:ntor, 0:mpol_nonzero))
     call cdf_write(ncid, vn_ZBS, ZBS(-ntor:ntor, 0:mpol_nonzero))
  case (general_option_scan)
     call cdf_write(ncid, vn_N_scan_array,  N_scan_array(1:axis_nmax+1,:))
     call cdf_write(ncid, vn_scan_R0c,  scan_R0c)
     call cdf_write(ncid, vn_scan_R0s,  scan_R0s)
     call cdf_write(ncid, vn_scan_Z0c,  scan_Z0c)
     call cdf_write(ncid, vn_scan_Z0s,  scan_Z0s)
  end select

  ! Arrays with dimension 3

  !call cdf_write(ncid, vn_r_plasma, r_plasma)

  ! Finish up:
  call cdf_close(ncid)

  call cpu_time(end_time)
  print "(a,es10.3,a)"," Time to write output file:",end_time - time1," sec"

end subroutine quasisymmetry_write_output

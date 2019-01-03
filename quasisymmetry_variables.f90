module quasisymmetry_variables

  use stel_kinds

  implicit none

  real(dp), parameter :: pi = 3.14159265358979d+0

  character(len=*), parameter :: &
       resolution_option_fixed = "fixed", &
       resolution_option_adaptive = "adaptive"
  character(len=50) :: resolution_option = resolution_option_fixed
  ! "fixed"    = Run using the specified N_phi.
  ! "adaptive" = Keep doubling N_phi (approximately, so N_phi remains odd) until iota_tolerance is achieved, or N_phi > max_N_phi.

  character(len=*), parameter :: &
       general_option_single = "single", &
       general_option_scan = "scan"
  character(len=50) :: general_option = general_option_single

  character(len=*), parameter :: &
       verbose_option_all = "all", &
       verbose_option_proc0 = "proc0", &
       verbose_option_summary = "summary"
  character(len=50) :: verbose_option = verbose_option_all

  character(len=*), parameter :: &
       eta_bar_scan_option_linear = "linear", &
       eta_bar_scan_option_log = "log"
  character(len=50) :: eta_bar_scan_option = eta_bar_scan_option_linear

  real(dp) :: sigma_initial = 0

  integer :: nfp = 3

  integer :: sign_G = 1
  integer :: sign_psi = 1
  real(dp) :: I2_over_B0 = 0

  integer :: N_iterations = 20
  integer :: N_line_search = 10
  real(dp) :: Newton_tolerance = 1.0d-12
  real(dp) :: iota_tolerance = 1.0d-6
  real(dp) :: elongation_tolerance = 1.0d-2

  integer :: N_phi = 15
  integer :: N_phi_original
  integer :: max_N_phi = 100

  integer, parameter :: max_axis_nmax = 5
  integer :: axis_nmax = 1
  real(dp), dimension(max_axis_nmax + 1) :: R0s, R0c, Z0s, Z0c ! Fourier coefficients for the magnetic axis
  real(dp) :: eta_bar

  integer :: matrix_size, helicity
  real(dp) :: last_iota, last_max_elongation, d_phi
  real(dp), dimension(:,:), allocatable :: d_d_phi, d_d_zeta
  real(dp), dimension(:), allocatable :: phi_extended, R0_extended, Z0_extended
  real(dp), dimension(:), allocatable :: phi, R0, Z0, R0p, Z0p, R0pp, Z0pp, R0ppp, Z0ppp
  real(dp), dimension(:), allocatable :: d_l_d_phi, curvature, torsion, B1Squared_over_curvatureSquared
  real(dp), dimension(:,:), allocatable :: tangent_cylindrical, normal_cylindrical, binormal_cylindrical
  real(dp), dimension(:,:), allocatable :: tangent_Cartesian, normal_Cartesian, binormal_Cartesian
  real(dp), dimension(:), allocatable :: sigma, X1s, X1c, Y1s, Y1c, R1s, R1c, Z1s, Z1c, elongation
  real(dp) :: B0_over_abs_G0, iota, max_elongation, rms_curvature, max_curvature, axis_length
  real(dp) :: max_precise_elongation = 10 ! Above this value, we won't do a precise solve, just take maxval over the phi grid.
  real(dp) :: max_elongation_to_keep = 10 ! Discard solutions with max(elongation) higher than this value. Set to e.g. 1.0e200 to keep all solutions.
  real(dp), dimension(:,:), allocatable :: Jacobian
  real(dp), dimension(:), allocatable :: residual, step_direction

  integer :: dimension_Fourier = 0
  real(dp), dimension(:,:), allocatable :: sin_n_phi, cos_n_phi

  character(len=200) :: output_filename
  character(len=200) :: vmec_template_filename = ''
  character(len=200) :: new_vmec_filename

  real(dp) :: total_time

  real(dp), dimension(max_axis_nmax+1) :: R0s_min, R0s_max, R0c_min, R0c_max, Z0s_min, Z0s_max, Z0c_min, Z0c_max
  real(dp) :: eta_bar_min = 1, eta_bar_max = 1, sigma_initial_min = 0, sigma_initial_max = 0
  integer, dimension(max_axis_nmax+1) :: R0s_N_scan=0, R0c_N_scan=0, Z0s_N_scan=0, Z0c_N_scan=0
  integer :: eta_bar_N_scan=0, sigma_initial_N_scan=0
  integer :: N_scan
  real(dp), dimension(:), allocatable :: iotas, max_elongations, rms_curvatures, max_curvatures, axis_lengths, eta_bar_values, sigma_initial_values
  integer, dimension(:), allocatable :: helicities
  logical, dimension(:), allocatable :: iota_tolerance_achieveds, elongation_tolerance_achieveds, Newton_tolerance_achieveds
  logical :: iota_tolerance_achieved, elongation_tolerance_achieved, Newton_tolerance_achieved
  integer, dimension(max_axis_nmax+1, 4) :: N_scan_array

  real(dp), dimension(:), allocatable :: scan_eta_bar, scan_sigma_initial
  real(dp), dimension(:,:), allocatable :: scan_R0c, scan_R0s, scan_Z0c, scan_Z0s
  real(dp) :: r = 0.1d+0

  integer :: N_procs, mpi_rank
  logical :: proc0, verbose = .true.

  namelist / quasisymmetry / resolution_option, general_option, verbose_option, nfp, sign_G, sign_psi, I2_over_B0, vmec_template_filename, r, &
       N_iterations, N_line_search, Newton_tolerance, iota_tolerance, elongation_tolerance, max_precise_elongation, max_elongation_to_keep, N_phi, max_N_phi, &
       R0s, R0c, Z0s, Z0c, eta_bar, sigma_initial, eta_bar_scan_option, &
       R0s_min, R0s_max, R0s_N_scan, R0c_min, R0c_max, R0c_N_scan, Z0s_min, Z0s_max, Z0s_N_scan, Z0c_min, Z0c_max, Z0c_N_scan, &
       eta_bar_min, eta_bar_max, eta_bar_N_scan, sigma_initial_min, sigma_initial_max, sigma_initial_N_scan

end module quasisymmetry_variables


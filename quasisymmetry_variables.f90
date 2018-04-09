module quasisymmetry_variables

  use stel_kinds

  implicit none

  real(dp), parameter :: pi = 3.14159265358979d+0

  integer :: resolution_option=1
  ! 1 = Run using the specified N_phi.
  ! 2 = Keep doubling N_phi (approximately, so N_phi remains odd) until iota_tolerance is achieved, or N_phi > max_N_phi.

  integer :: general_option=1

  integer :: nfp = 3

  integer :: sign_G = 1
  real(dp) :: I2_over_B0 = 0

  integer :: N_iterations = 20
  integer :: N_line_search = 10
  real(dp) :: Newton_tolerance = 1.0d-12
  real(dp) :: iota_tolerance = 1.0d-6
  real(dp) :: elongation_tolerance = 1.0d-2

  integer :: N_phi = 15
  integer :: N_phi_original
  integer :: max_N_phi = 100

  integer, parameter :: max_axis_nmax = 2
  integer :: axis_nmax = 1
  real(dp), dimension(max_axis_nmax + 1) :: R0s, R0c, Z0s, Z0c ! Fourier coefficients for the magnetic axis
  real(dp) :: B1s_over_B0, B1c_over_B0

  integer :: matrix_size, helicity
  real(dp) :: last_iota, last_max_elongation, d_phi
  real(dp), dimension(:,:), allocatable :: d_d_phi, d_d_zeta
  real(dp), dimension(:), allocatable :: phi_extended, R0_extended, Z0_extended
  real(dp), dimension(:), allocatable :: phi, R0, Z0, R0p, Z0p, R0pp, Z0pp, R0ppp, Z0ppp
  real(dp), dimension(:), allocatable :: d_l_d_phi, curvature, torsion, B1Squared_over_curvatureSquared
  real(dp), dimension(:,:), allocatable :: tangent_cylindrical, normal_cylindrical, binormal_cylindrical
  real(dp), dimension(:,:), allocatable :: tangent_Cartesian, normal_Cartesian, binormal_Cartesian
  real(dp), dimension(:), allocatable :: sigma, X1s, X1c, Y1s, Y1c, R1s, R1c, Z1s, Z1c, elongation
  real(dp) :: B0_over_abs_G0, iota, max_elongation
  real(dp), dimension(:,:), allocatable :: Jacobian
  real(dp), dimension(:), allocatable :: residual, step_direction

  integer :: dimension_Fourier = 0
  real(dp), dimension(:,:), allocatable :: sin_n_phi, cos_n_phi

  character(len=200) :: output_filename

  real(dp) :: total_time

  real(dp), dimension(max_axis_nmax+1) :: R0s_min, R0s_max, R0c_min, R0c_max, Z0s_min, Z0s_max, Z0c_min, Z0c_max
  real(dp) :: B1s_min, B1s_max, B1c_min, B1c_max
  integer, dimension(max_axis_nmax+1) :: R0s_N_scan=0, R0c_N_scan=0, Z0s_N_scan=0, Z0c_N_scan=0
  integer :: B1s_N_scan=0, B1c_N_scan=0
  integer :: N_scan
  real(dp), dimension(:), allocatable :: iotas, max_elongations
  integer, dimension(:), allocatable :: helicities
  logical, dimension(:), allocatable :: iota_tolerance_achieveds, elongation_tolerance_achieveds, Newton_tolerance_achieveds
  logical :: iota_tolerance_achieved, elongation_tolerance_achieved, Newton_tolerance_achieved
  integer, dimension(max_axis_nmax+1, 4) :: N_scan_array

  real(dp), dimension(:), allocatable :: scan_B1c, scan_B1s
  real(dp), dimension(:,:), allocatable :: scan_R0c, scan_R0s, scan_Z0c, scan_Z0s

end module quasisymmetry_variables


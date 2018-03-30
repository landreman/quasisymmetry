module quasisymmetry_variables

  use stel_kinds

  implicit none

  real(dp), parameter :: pi = 3.14159265358979d+0

  integer :: general_option=1

  integer :: nfp = 3

  integer :: sign_G = 1
  real(dp) :: I2_over_B0 = 0

  integer :: N_iterations = 20
  integer :: N_line_search = 10
  real(dp) :: Newton_tolerance = 1.0d-12

  integer :: N_N_phis = 2
  integer, parameter :: max_N_N_phis = 10
  integer, dimension(max_N_N_phis) :: N_phis = (/ 15, 31, 0, 0, 0, 0, 0, 0, 0, 0 /)
  integer :: N_phi = 0

  integer, parameter :: max_axis_nmax = 1
  integer :: axis_nmax = 1
  real(dp), dimension(max_axis_nmax + 1) :: R0s, R0c, Z0s, Z0c ! Fourier coefficients for the magnetic axis
  real(dp) :: B1s_over_B0, B1c_over_B0

  integer :: matrix_size
  real(dp), dimension(max_N_N_phis) :: iotas, max_elongations
  real(dp), dimension(:,:), allocatable :: d_d_phi, d_d_zeta
  real(dp), dimension(:), allocatable :: phi_extended, R0_extended, Z0_extended
  real(dp), dimension(:), allocatable :: phi, R0, Z0, R0p, Z0p, R0pp, Z0pp, R0ppp, Z0ppp
  real(dp), dimension(:), allocatable :: d_l_d_phi, curvature, torsion, B1Squared_over_curvatureSquared
  real(dp), dimension(:), allocatable :: RZ_to_XY_a, RZ_to_XY_b, RZ_to_XY_d
  real(dp), dimension(:,:), allocatable :: tangent_cylindrical, normal_cylindrical, binormal_cylindrical
  real(dp), dimension(:,:), allocatable :: tangent_Cartesian, normal_Cartesian, binormal_Cartesian
  real(dp), dimension(:), allocatable :: sigma, X1s, X1c, Y1s, Y1c, R1s, R1c, Z1s, Z1c, elongation
  real(dp) :: B0_over_abs_G0, iota, max_elongation
  real(dp), dimension(:,:), allocatable :: Jacobian
  real(dp), dimension(:), allocatable :: residual, step_direction

  character(len=200) :: output_filename

  real(dp) :: total_time

end module quasisymmetry_variables


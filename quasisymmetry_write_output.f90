subroutine quasisymmetry_write_output

  use quasisymmetry_variables
  use ezcdf

  implicit none

  integer :: ierr, ncid

  ! Same convention as in VMEC:
  ! Prefix vn_ indicates the variable name used in the .nc file.

  ! Scalars:
  character(len=*), parameter :: &
       vn_nfp = "nfp", &
       vn_resolution_option = "resolution_option", &
       vn_B1s_min = "B1s_min", &
       vn_B1s_max = "B1s_max", &
       vn_B1s_N_scan = "B1s_N_scan", &
       vn_B1c_min = "B1c_min", &
       vn_B1c_max = "B1c_max", &
       vn_B1c_N_scan = "B1c_N_scan"

  ! Arrays with dimension 1
  character(len=*), parameter :: &
       vn_iotas = "iotas", &
       vn_max_elongations = "max_elongations", &
       vn_helicities = "helicities", &
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
       vn_Z0c_N_scan = "Z0c_N_scan"

  ! Arrays with dimension 2
  character(len=*), parameter :: &
       vn_N_scan_array  = "N_scan_array"

!!$  ! Arrays with dimension 3
!!$  character(len=*), parameter :: &
!!$       vn_r_plasma  = "r_plasma"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now create variables that name the dimensions.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Arrays with dimension 1:
  character(len=*), parameter, dimension(1) :: &
       N_scan_dim = (/'N_scan'/), &
       max_axis_nmax_plus_1_dim = (/'max_axis_nmax_plus_1'/)

  ! Arrays with dimension 2:
  ! The form of the array declarations here is inspired by
  ! http://stackoverflow.com/questions/21552430/gfortran-does-not-allow-character-arrays-with-varying-component-lengths
  character(len=*), parameter, dimension(2) :: &
       max_axis_nmax_plus_1_4_dim = (/ character(len=50) :: 'max_axis_nmax_plus_1','4'/)

!!$  ! Arrays with dimension 3:
!!$  character(len=*), parameter, dimension(3) :: &
!!$       xyz_ntheta_nzetal_plasma_dim = (/ character(len=50) :: 'xyz','ntheta_plasma','nzetal_plasma'/), &

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call cdf_open(ncid,output_filename,'w',ierr)
  IF (ierr .ne. 0) then
     print *,"Error opening output file ",output_filename
     stop
  end IF

  ! Scalars

  call cdf_define(ncid, vn_nfp, nfp)
  call cdf_setatt(ncid, vn_nfp, 'Number of field periods, i.e. the number of identical toroidal segments, 5 for W7-X, 4 for HSX, etc. ' // &
       'Equivalent to the VMEC variable of the same name.')

  call cdf_define(ncid, vn_resolution_option, resolution_option)
  !call cdf_setatt(ncid, vn_resolution_option, 'Method used to define the geometry of the plasma surface.' // input_parameter_text)

  call cdf_define(ncid, vn_B1s_min, B1s_min)
  call cdf_define(ncid, vn_B1s_max, B1s_max)
  call cdf_define(ncid, vn_B1s_N_scan, B1s_N_scan)
  call cdf_define(ncid, vn_B1c_min, B1c_min)
  call cdf_define(ncid, vn_B1c_max, B1c_max)
  call cdf_define(ncid, vn_B1c_N_scan, B1c_N_scan)

  ! Arrays with dimension 1

  call cdf_define(ncid, vn_iotas, iotas, dimname=N_scan_dim)
  call cdf_define(ncid, vn_max_elongations, max_elongations, dimname=N_scan_dim)
  call cdf_define(ncid, vn_helicities, helicities, dimname=N_scan_dim)
  call cdf_define(ncid, vn_Newton_tolerance_achieveds, Newton_tolerance_achieveds, dimname=N_scan_dim)
  call cdf_define(ncid, vn_iota_tolerance_achieveds, iota_tolerance_achieveds, dimname=N_scan_dim)
  call cdf_define(ncid, vn_elongation_tolerance_achieveds, elongation_tolerance_achieveds, dimname=N_scan_dim)
  call cdf_define(ncid, vn_R0s_min, R0s_min, dimname=max_axis_nmax_plus_1_dim)
  call cdf_define(ncid, vn_R0s_max, R0s_max, dimname=max_axis_nmax_plus_1_dim)
  call cdf_define(ncid, vn_R0s_N_scan, R0s_N_scan, dimname=max_axis_nmax_plus_1_dim)
  call cdf_define(ncid, vn_R0c_min, R0c_min, dimname=max_axis_nmax_plus_1_dim)
  call cdf_define(ncid, vn_R0c_max, R0c_max, dimname=max_axis_nmax_plus_1_dim)
  call cdf_define(ncid, vn_R0c_N_scan, R0c_N_scan, dimname=max_axis_nmax_plus_1_dim)
  call cdf_define(ncid, vn_Z0s_min, Z0s_min, dimname=max_axis_nmax_plus_1_dim)
  call cdf_define(ncid, vn_Z0s_max, Z0s_max, dimname=max_axis_nmax_plus_1_dim)
  call cdf_define(ncid, vn_Z0s_N_scan, Z0s_N_scan, dimname=max_axis_nmax_plus_1_dim)
  call cdf_define(ncid, vn_Z0c_min, Z0c_min, dimname=max_axis_nmax_plus_1_dim)
  call cdf_define(ncid, vn_Z0c_max, Z0c_max, dimname=max_axis_nmax_plus_1_dim)
  call cdf_define(ncid, vn_Z0c_N_scan, Z0c_N_scan, dimname=max_axis_nmax_plus_1_dim)


  ! Arrays with dimension 2

  call cdf_define(ncid, vn_N_scan_array,  N_scan_array, dimname=max_axis_nmax_plus_1_4_dim)

  ! Arrays with dimension 3

  !call cdf_define(ncid, vn_r_plasma,  r_plasma,  dimname=xyz_ntheta_nzetal_plasma_dim)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ! Done with cdf_define calls. Now write the data.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  ! Scalars

  call cdf_write(ncid, vn_nfp, nfp)
  call cdf_write(ncid, vn_resolution_option, resolution_option)
  call cdf_write(ncid, vn_B1s_min, B1s_min)
  call cdf_write(ncid, vn_B1s_max, B1s_max)
  call cdf_write(ncid, vn_B1s_N_scan, B1s_N_scan)
  call cdf_write(ncid, vn_B1c_min, B1c_min)
  call cdf_write(ncid, vn_B1c_max, B1c_max)
  call cdf_write(ncid, vn_B1c_N_scan, B1c_N_scan)

  ! Arrays with dimension 1

  call cdf_write(ncid, vn_iotas, iotas)
  call cdf_write(ncid, vn_max_elongations, max_elongations)
  call cdf_write(ncid, vn_helicities, helicities)
  call cdf_write(ncid, vn_Newton_tolerance_achieveds, Newton_tolerance_achieveds)
  call cdf_write(ncid, vn_iota_tolerance_achieveds, iota_tolerance_achieveds)
  call cdf_write(ncid, vn_elongation_tolerance_achieveds, elongation_tolerance_achieveds)
  call cdf_write(ncid, vn_R0s_min, R0s_min)
  call cdf_write(ncid, vn_R0s_max, R0s_max)
  call cdf_write(ncid, vn_R0s_N_scan, R0s_N_scan)
  call cdf_write(ncid, vn_R0c_min, R0c_min)
  call cdf_write(ncid, vn_R0c_max, R0c_max)
  call cdf_write(ncid, vn_R0c_N_scan, R0c_N_scan)
  call cdf_write(ncid, vn_Z0s_min, Z0s_min)
  call cdf_write(ncid, vn_Z0s_max, Z0s_max)
  call cdf_write(ncid, vn_Z0s_N_scan, Z0s_N_scan)
  call cdf_write(ncid, vn_Z0c_min, Z0c_min)
  call cdf_write(ncid, vn_Z0c_max, Z0c_max)
  call cdf_write(ncid, vn_Z0c_N_scan, Z0c_N_scan)

  ! Arrays with dimension 2

  call cdf_write(ncid, vn_N_scan_array,  N_scan_array)

  ! Arrays with dimension 3

  !call cdf_write(ncid, vn_r_plasma, r_plasma)

  ! Finish up:
  call cdf_close(ncid)

end subroutine quasisymmetry_write_output

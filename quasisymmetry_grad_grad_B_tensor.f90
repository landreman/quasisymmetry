subroutine quasisymmetry_grad_grad_B_tensor

  use quasisymmetry_variables

  real(dp), allocatable :: grad_grad_B(:,:,:,:)
  real(dp), allocatable, dimension(:) :: d_curvature_d_zeta, d_torsion_d_zeta, d_X20_d_zeta, d_X2s_d_zeta, d_X2c_d_zeta, d_Y20_d_zeta, d_Y2s_d_zeta, d_Y2c_d_zeta
  real(dp), allocatable, dimension(:) :: d_Z20_d_zeta, d_Z2s_d_zeta, d_Z2c_d_zeta, d2_X1c_d_zeta2, d2_Y1c_d_zeta2, d2_Y1s_d_zeta2, difference
  real(dp) :: G0, iota_N0, lp, G2, I2, norm_difference
  integer :: i,j,k

  ! Only proceed if we have O(r^2) information
  if ((trim(order_r_option)==order_r_option_r1) .or. (trim(order_r_option)==order_r_option_r1_compute_B2)) return

  allocate(grad_grad_B(N_phi,3,3,3))
  allocate(difference(N_phi))
  allocate(d_curvature_d_zeta(N_phi))
  allocate(d_torsion_d_zeta(N_phi))
  allocate(d_X20_d_zeta(N_phi))
  allocate(d_X2s_d_zeta(N_phi))
  allocate(d_X2c_d_zeta(N_phi))
  allocate(d_Y20_d_zeta(N_phi))
  allocate(d_Y2s_d_zeta(N_phi))
  allocate(d_Y2c_d_zeta(N_phi))
  allocate(d_Z20_d_zeta(N_phi))
  allocate(d_Z2s_d_zeta(N_phi))
  allocate(d_Z2c_d_zeta(N_phi))
  allocate(d2_X1c_d_zeta2(N_phi))
  allocate(d2_Y1c_d_zeta2(N_phi))
  allocate(d2_Y1s_d_zeta2(N_phi))

  iota_N0 = iota + axis_helicity*nfp
  lp = abs_G0_over_B0
  I2 = I2_over_B0 * B0
  G0 = sign_G * abs_G0_over_B0 * B0
  G2 = -mu0 * p2 * G0 / (B0 * B0) - iota * I2

  d_curvature_d_zeta = matmul(d_d_zeta, curvature)
  d_torsion_d_zeta = matmul(d_d_zeta, torsion)
  d_X20_d_zeta = matmul(d_d_zeta, X20)
  d_X2s_d_zeta = matmul(d_d_zeta, X2s)
  d_X2c_d_zeta = matmul(d_d_zeta, X2c)
  d_Y20_d_zeta = matmul(d_d_zeta, Y20)
  d_Y2s_d_zeta = matmul(d_d_zeta, Y2s)
  d_Y2c_d_zeta = matmul(d_d_zeta, Y2c)
  d_Z20_d_zeta = matmul(d_d_zeta, Z20)
  d_Z2s_d_zeta = matmul(d_d_zeta, Z2s)
  d_Z2c_d_zeta = matmul(d_d_zeta, Z2c)
  d2_X1c_d_zeta2 = matmul(d_d_zeta, d_X1c_d_zeta)
  d2_Y1c_d_zeta2 = matmul(d_d_zeta, d_Y1c_d_zeta)
  d2_Y1s_d_zeta2 = matmul(d_d_zeta, d_Y1s_d_zeta)

  ! The elements that follow are computed in the Mathematica notebook "20200407-01 Grad grad B tensor near axis"
  ! and then formatted for fortran by the python script process_grad_grad_B_tensor_code

  ! Element 111
  grad_grad_B(:,1,1,1) =(B0*B0*B0*B0*lp*lp*(8*iota_N0*X2c*Y1c*&
       Y1s + 4*iota_N0*X2s*&
       (-Y1c*Y1c + Y1s*Y1s) + &
       2*iota_N0*X1c*Y1s*Y20 + &
       2*iota_N0*X1c*Y1s*Y2c - &
       2*iota_N0*X1c*Y1c*Y2s + &
       5*iota_N0*X1c*X1c*Y1c*Y1s*&
       curvature - &
       2*Y1c*Y20*d_X1c_d_zeta + &
       2*Y1c*Y2c*d_X1c_d_zeta + &
       2*Y1s*Y2s*d_X1c_d_zeta + &
       5*X1c*Y1s*Y1s*curvature*&
       d_X1c_d_zeta + &
       2*Y1c*Y1c*d_X20_d_zeta + &
       2*Y1s*Y1s*d_X20_d_zeta - &
       2*Y1c*Y1c*d_X2c_d_zeta + &
       2*Y1s*Y1s*d_X2c_d_zeta - &
       4*Y1c*Y1s*d_X2s_d_zeta))/&
       (G0*G0*G0)
  
  ! Element 112
  grad_grad_B(:,1,1,2) =(B0*B0*B0*B0*lp*lp*(Y1c*Y1c*&
       (-6*iota_N0*Y2s + &
       5*iota_N0*X1c*Y1s*&
       curvature + &
       2*(lp*X20*torsion - &
       lp*X2c*torsion + &
       d_Y20_d_zeta - &
       d_Y2c_d_zeta)) + &
       Y1s*(5*iota_N0*X1c*Y1s*Y1s*&
       curvature + &
       2*(lp*X1c*Y2s*torsion + &
       Y2s*d_Y1c_d_zeta - &
       (Y20 + Y2c)*&
       d_Y1s_d_zeta) + &
       Y1s*(6*iota_N0*Y2s + &
       2*lp*X20*torsion + &
       2*lp*X2c*torsion + &
       5*lp*X1c*X1c*curvature*&
       torsion + &
       5*X1c*curvature*&
       d_Y1c_d_zeta + &
       2*d_Y20_d_zeta + &
       2*d_Y2c_d_zeta)) + &
       Y1c*(2*(lp*X1c*&
       (-Y20 + Y2c)*torsion - &
       Y20*d_Y1c_d_zeta + &
       Y2c*d_Y1c_d_zeta + &
       Y2s*d_Y1s_d_zeta) + &
       Y1s*(12*iota_N0*Y2c - &
       4*lp*X2s*torsion - &
       5*X1c*curvature*&
       d_Y1s_d_zeta - &
       4*d_Y2s_d_zeta))))/(G0*G0*G0)
  
  ! Element 113
  grad_grad_B(:,1,1,3) =-((B0*B0*B0*lp*lp*(2*Y1c*Y1c*&
       (2*B2c*G0*lp + B0*G2*lp + B0*I2*lp*iota - &
       2*G0*lp*B20 + 2*B0*G0*iota_N0*Z2s + &
       B0*G0*lp*X20*curvature - &
       B0*G0*lp*X2c*curvature - &
       B0*G0*d_Z20_d_zeta + &
       B0*G0*d_Z2c_d_zeta) + &
       Y1s*(-2*B0*G0*lp*X1c*Y2s*&
       curvature + &
       Y1s*(-4*B2c*G0*lp + 2*B0*G2*lp + &
       2*B0*I2*lp*iota - 4*G0*lp*B20 - &
       4*B0*G0*iota_N0*Z2s + &
       2*B0*G0*lp*X20*curvature + &
       2*B0*G0*lp*X2c*curvature + &
       B0*G0*lp*X1c*X1c*curvature*curvature - &
       2*B0*G0*d_Z20_d_zeta - &
       2*B0*G0*d_Z2c_d_zeta)) + &
       2*G0*Y1c*(B0*lp*X1c*&
       (Y20 - Y2c)*curvature + &
       2*Y1s*(2*B2s*lp - 2*B0*iota_N0*Z2c - &
       B0*lp*X2s*curvature + &
       B0*d_Z2s_d_zeta))))/(G0*G0*G0*G0))
  
  ! Element 121
  grad_grad_B(:,1,2,1) =-((B0*B0*B0*B0*lp*lp*(3*iota_N0*X1c*X1c*X1c*Y1s*&
       curvature + &
       3*lp*X1c*X1c*Y1s*Y1s*curvature*&
       torsion + &
       2*(X2s*Y1s*&
       (-2*lp*Y1c*torsion + &
       d_X1c_d_zeta) + &
       X20*(lp*Y1c*Y1c*torsion + &
       lp*Y1s*Y1s*torsion - &
       Y1c*d_X1c_d_zeta) + &
       X2c*(-(lp*Y1c*Y1c*&
       torsion) + &
       lp*Y1s*Y1s*torsion + &
       Y1c*d_X1c_d_zeta)) - &
       2*X1c*(3*iota_N0*X2s*Y1c - &
       iota_N0*X20*Y1s - &
       3*iota_N0*X2c*Y1s + &
       lp*Y1c*Y20*torsion - &
       lp*Y1c*Y2c*torsion - &
       lp*Y1s*Y2s*torsion - &
       Y1c*d_X20_d_zeta + &
       Y1c*d_X2c_d_zeta + &
       Y1s*d_X2s_d_zeta)))/&
       (G0*G0*G0))
  
  ! Element 122
  grad_grad_B(:,1,2,2) =(B0*B0*B0*B0*lp*lp*(-4*iota_N0*X1c*Y1s*&
       Y2c + 4*iota_N0*X1c*Y1c*&
       Y2s - 3*iota_N0*X1c*X1c*Y1c*&
       Y1s*curvature + &
       2*X20*Y1c*d_Y1c_d_zeta + &
       2*X20*Y1s*d_Y1s_d_zeta + &
       3*X1c*X1c*Y1s*curvature*&
       d_Y1s_d_zeta + &
       2*X2s*(iota_N0*Y1c*Y1c - &
       Y1s*(iota_N0*Y1s + &
       d_Y1c_d_zeta) - &
       Y1c*d_Y1s_d_zeta) - &
       2*X2c*(Y1c*&
       (2*iota_N0*Y1s + d_Y1c_d_zeta) & ! Removed \
       - Y1s*d_Y1s_d_zeta) - &
       2*X1c*Y1c*d_Y20_d_zeta + &
       2*X1c*Y1c*d_Y2c_d_zeta + &
       2*X1c*Y1s*d_Y2s_d_zeta))/&
       (G0*G0*G0)
!       (2*iota_N0*Y1s + d_Y1c_d_zeta) \&
  
  ! Element 123
  grad_grad_B(:,1,2,3) =(2*B0*B0*B0*lp*lp*X1c*&
       (Y1c*(2*B2c*G0*lp + B0*G2*lp + B0*I2*lp*iota - &
       2*G0*lp*B20 + 2*B0*G0*iota_N0*Z2s + &
       2*B0*G0*lp*X20*curvature - &
       2*B0*G0*lp*X2c*curvature - &
       B0*G0*d_Z20_d_zeta + &
       B0*G0*d_Z2c_d_zeta) + &
       G0*Y1s*(2*B2s*lp - 2*B0*iota_N0*Z2c - &
       2*B0*lp*X2s*curvature + &
       B0*d_Z2s_d_zeta)))/(G0*G0*G0*G0)
  
  ! Element 131
  grad_grad_B(:,1,3,1) =(B0*B0*B0*B0*lp*(-4*lp*lp*X2s*Y1c*Y1s*&
       curvature + &
       2*lp*lp*X2c*(-Y1c*Y1c + Y1s*Y1s)*&
       curvature + &
       2*lp*lp*X20*(Y1c*Y1c + Y1s*Y1s)*&
       curvature - &
       2*lp*lp*X1c*Y1c*Y20*&
       curvature + &
       2*lp*lp*X1c*Y1c*Y2c*&
       curvature + &
       2*lp*lp*X1c*Y1s*Y2s*&
       curvature + &
       3*lp*lp*X1c*X1c*Y1s*Y1s*&
       curvature*curvature + &
       lp*iota_N0*X1c*X1c*X1c*Y1s*&
       torsion - lp*iota_N0*X1c*&
       Y1c*Y1c*Y1s*torsion - &
       lp*iota_N0*X1c*Y1s*Y1s*Y1s*&
       torsion - Y1s*Y1s*&
       d_X1c_d_zeta*d_X1c_d_zeta + &
       iota_N0*X1c*X1c*Y1s*&
       d_Y1c_d_zeta - &
       lp*X1c*Y1s*Y1s*torsion*&
       d_Y1c_d_zeta - &
       iota_N0*X1c*X1c*Y1c*&
       d_Y1s_d_zeta + &
       lp*X1c*Y1c*Y1s*&
       torsion*d_Y1s_d_zeta + &
       X1c*Y1s*Y1s*d2_X1c_d_zeta2))/&
       (G0*G0*G0)
  
  ! Element 132
  grad_grad_B(:,1,3,2) =(B0*B0*B0*B0*lp*(-(Y1s*d_X1c_d_zeta*&
       (iota_N0*Y1c*Y1c + &
       Y1s*(iota_N0*Y1s + &
       d_Y1c_d_zeta) - &
       Y1c*d_Y1s_d_zeta)) + &
       lp*X1c*X1c*Y1s*&
       (2*iota_N0*Y1c*torsion - &
       torsion*d_Y1s_d_zeta + &
       Y1s*d_torsion_d_zeta) + &
       X1c*(Y1c*d_Y1s_d_zeta*&
       (-(iota_N0*Y1c) + d_Y1s_d_zeta) & ! Removed \
       + Y1s*Y1s*(lp*torsion*&
       d_X1c_d_zeta + &
       iota_N0*d_Y1s_d_zeta + &
       d2_Y1c_d_zeta2) - &
       Y1s*(d_Y1c_d_zeta*&
       d_Y1s_d_zeta + &
       Y1c*(-2*iota_N0*d_Y1c_d_zeta + &
       d2_Y1s_d_zeta2)))))/(G0*G0*G0)
!       (-(iota_N0*Y1c) + d_Y1s_d_zeta) \&
  
  ! Element 133
  grad_grad_B(:,1,3,3) =(B0*B0*B0*B0*lp*lp*X1c*Y1s*&
       (-(Y1s*curvature*&
       d_X1c_d_zeta) + &
       X1c*(-(iota_N0*Y1c*&
       curvature) + &
       Y1s*d_curvature_d_zeta)))/&
       (G0*G0*G0)
  
  ! Element 211
  grad_grad_B(:,2,1,1) =(-2*B0*B0*B0*B0*lp*lp*X1c*&
       (-2*iota_N0*X2s*Y1c + &
       2*iota_N0*X2c*Y1s - &
       iota_N0*X1c*Y2s + &
       iota_N0*X1c*X1c*Y1s*curvature + &
       lp*X1c*Y1s*Y1s*curvature*&
       torsion - Y20*&
       d_X1c_d_zeta + &
       Y2c*d_X1c_d_zeta + &
       Y1c*d_X20_d_zeta - &
       Y1c*d_X2c_d_zeta - &
       Y1s*d_X2s_d_zeta))/(G0*G0*G0)
  
  ! Element 212
  grad_grad_B(:,2,1,2) =(2*B0*B0*B0*B0*lp*lp*X1c*&
       (lp*X1c*Y20*torsion - &
       lp*X1c*Y2c*torsion + &
       Y20*d_Y1c_d_zeta - &
       Y2c*d_Y1c_d_zeta - &
       Y2s*d_Y1s_d_zeta + &
       Y1c*(3*iota_N0*Y2s - &
       lp*X20*torsion + &
       lp*X2c*torsion - &
       d_Y20_d_zeta + d_Y2c_d_zeta) & ! Removed \ at end
       + Y1s*(iota_N0*Y20 - &
       3*iota_N0*Y2c - &
       iota_N0*X1c*Y1c*curvature + &
       lp*X2s*torsion + &
       X1c*curvature*&
       d_Y1s_d_zeta + d_Y2s_d_zeta))& ! Removed \ at end
       )/(G0*G0*G0)
!       d_Y20_d_zeta + d_Y2c_d_zeta) \&
!       d_Y1s_d_zeta + d_Y2s_d_zeta))\&
  
  ! Element 213
  grad_grad_B(:,2,1,3) =(2*B0*B0*B0*lp*lp*X1c*&
       (Y1c*(2*B2c*G0*lp + B0*G2*lp + B0*I2*lp*iota - &
       2*G0*lp*B20 + 2*B0*G0*iota_N0*Z2s + &
       B0*G0*lp*X20*curvature - &
       B0*G0*lp*X2c*curvature - &
       B0*G0*d_Z20_d_zeta + &
       B0*G0*d_Z2c_d_zeta) + &
       G0*(B0*lp*X1c*(Y20 - Y2c)*&
       curvature + &
       Y1s*(2*B2s*lp - 2*B0*iota_N0*Z2c - &
       B0*lp*X2s*curvature + &
       B0*d_Z2s_d_zeta))))/(G0*G0*G0*G0)
  
  ! Element 221
  grad_grad_B(:,2,2,1) =(-2*B0*B0*B0*B0*lp*lp*X1c*&
       (lp*X2c*Y1c*torsion + &
       lp*X2s*Y1s*torsion - &
       X2c*d_X1c_d_zeta + &
       X20*(-(lp*Y1c*torsion) + &
       d_X1c_d_zeta) + &
       X1c*(3*iota_N0*X2s + &
       lp*Y20*torsion - &
       lp*Y2c*torsion - &
       d_X20_d_zeta + d_X2c_d_zeta)))/&
       (G0*G0*G0)
  
  ! Element 222
  grad_grad_B(:,2,2,2) =(-2*B0*B0*B0*B0*lp*lp*X1c*&
       (-(iota_N0*X2c*Y1s) + &
       2*iota_N0*X1c*Y2s - &
       X2c*d_Y1c_d_zeta + &
       X20*(iota_N0*Y1s + &
       d_Y1c_d_zeta) + &
       X2s*(iota_N0*Y1c - &
       d_Y1s_d_zeta) - &
       X1c*d_Y20_d_zeta + &
       X1c*d_Y2c_d_zeta))/(G0*G0*G0)
  
  ! Element 223
  grad_grad_B(:,2,2,3) =(-2*B0*B0*B0*lp*lp*X1c*X1c*&
       (2*B2c*G0*lp + B0*G2*lp + B0*I2*lp*iota - 2*G0*lp*B20 + &
       2*B0*G0*iota_N0*Z2s + &
       2*B0*G0*lp*X20*curvature - &
       2*B0*G0*lp*X2c*curvature - &
       B0*G0*d_Z20_d_zeta + &
       B0*G0*d_Z2c_d_zeta))/(G0*G0*G0*G0)
  
  ! Element 231
  grad_grad_B(:,2,3,1) =(B0*B0*B0*B0*lp*X1c*(-2*lp*lp*X20*Y1c*&
       curvature + &
       2*lp*lp*X2c*Y1c*curvature + &
       2*lp*lp*X2s*Y1s*curvature + &
       2*lp*lp*X1c*Y20*curvature - &
       2*lp*lp*X1c*Y2c*curvature + &
       2*lp*iota_N0*X1c*Y1c*Y1s*&
       torsion - iota_N0*X1c*Y1s*&
       d_X1c_d_zeta + &
       lp*Y1s*Y1s*torsion*&
       d_X1c_d_zeta + &
       iota_N0*X1c*X1c*d_Y1s_d_zeta - &
       lp*X1c*Y1s*torsion*&
       d_Y1s_d_zeta - &
       lp*X1c*Y1s*Y1s*&
       d_torsion_d_zeta))/(G0*G0*G0)
  
  ! Element 232
  grad_grad_B(:,2,3,2) =(B0*B0*B0*B0*lp*X1c*(-(lp*iota_N0*X1c*X1c*&
       Y1s*torsion) + &
       lp*Y1s*torsion*&
       (iota_N0*Y1c*Y1c + &
       Y1s*(iota_N0*Y1s + &
       d_Y1c_d_zeta) - &
       Y1c*d_Y1s_d_zeta) + &
       X1c*((iota_N0*Y1c - &
       d_Y1s_d_zeta)*d_Y1s_d_zeta & ! Removed \ at end
       + Y1s*(-(iota_N0*d_Y1c_d_zeta) + &
       d2_Y1s_d_zeta2))))/(G0*G0*G0)
!       d_Y1s_d_zeta)*d_Y1s_d_zeta \&
  
  ! Element 233
  grad_grad_B(:,2,3,3) =(B0*B0*B0*B0*lp*lp*X1c*X1c*Y1s*curvature*&
       (iota_N0*X1c + 2*lp*Y1s*torsion))/&
       (G0*G0*G0)
  
  ! Element 311
  grad_grad_B(:,3,1,1) =(B0*B0*B0*B0*lp*X1c*Y1s*&
       (lp*iota_N0*X1c*X1c*torsion - &
       lp*iota_N0*Y1c*Y1c*torsion - &
       lp*iota_N0*Y1s*Y1s*torsion - &
       lp*Y1s*torsion*&
       d_Y1c_d_zeta + &
       X1c*(2*lp*lp*Y1s*curvature*curvature + &
       iota_N0*d_Y1c_d_zeta) + &
       d_X1c_d_zeta*d_Y1s_d_zeta + &
       Y1c*(iota_N0*d_X1c_d_zeta + &
       lp*torsion*d_Y1s_d_zeta) + &
       Y1s*d2_X1c_d_zeta2))/(G0*G0*G0)
  
  ! Element 312
  grad_grad_B(:,3,1,2) =(B0*B0*B0*B0*lp*X1c*Y1s*&
       (lp*X1c*(2*iota_N0*Y1c*&
       torsion + &
       Y1s*d_torsion_d_zeta) + &
       Y1s*(2*lp*torsion*&
       d_X1c_d_zeta + &
       2*iota_N0*d_Y1s_d_zeta + &
       d2_Y1c_d_zeta2) + &
       Y1c*(2*iota_N0*d_Y1c_d_zeta - &
       d2_Y1s_d_zeta2)))/(G0*G0*G0)
  
  ! Element 313
  grad_grad_B(:,3,1,3) =(B0*B0*B0*B0*lp*lp*X1c*X1c*Y1s*&
       (-(iota_N0*Y1c*curvature) + &
       curvature*d_Y1s_d_zeta + &
       Y1s*d_curvature_d_zeta))/&
       (G0*G0*G0)
  
  ! Element 321
  grad_grad_B(:,3,2,1) =-((B0*B0*B0*B0*lp*X1c*X1c*Y1s*&
       (-2*lp*iota_N0*Y1c*torsion + &
       2*iota_N0*d_X1c_d_zeta + &
       2*lp*torsion*d_Y1s_d_zeta + &
       lp*Y1s*d_torsion_d_zeta))/&
       (G0*G0*G0))
  
  ! Element 322
  grad_grad_B(:,3,2,2) =-((B0*B0*B0*B0*lp*X1c*Y1s*&
       (lp*iota_N0*X1c*X1c*torsion - &
       lp*iota_N0*Y1c*Y1c*torsion - &
       lp*iota_N0*Y1s*Y1s*torsion - &
       lp*Y1s*torsion*&
       d_Y1c_d_zeta - &
       d_X1c_d_zeta*d_Y1s_d_zeta + &
       Y1c*(iota_N0*d_X1c_d_zeta + &
       lp*torsion*d_Y1s_d_zeta) + &
       X1c*(iota_N0*d_Y1c_d_zeta - &
       d2_Y1s_d_zeta2)))/(G0*G0*G0))
  
  ! Element 323
  grad_grad_B(:,3,2,3) =(B0*B0*B0*B0*lp*lp*X1c*X1c*Y1s*curvature*&
       (iota_N0*X1c + 2*lp*Y1s*torsion))/&
       (G0*G0*G0)
  
  ! Element 331
  grad_grad_B(:,3,3,1) =(B0*B0*B0*B0*lp*lp*X1c*X1c*Y1s*&
       (-(iota_N0*Y1c*curvature) + &
       curvature*d_Y1s_d_zeta + &
       Y1s*d_curvature_d_zeta))/&
       (G0*G0*G0)
  
  ! Element 332
  grad_grad_B(:,3,3,2) =-((B0*B0*B0*B0*lp*lp*X1c*Y1s*curvature*&
       (iota_N0*Y1c*Y1c + &
       Y1s*(iota_N0*Y1s + &
       d_Y1c_d_zeta) - &
       Y1c*d_Y1s_d_zeta))/(G0*G0*G0))
  
  ! Element 333
  grad_grad_B(:,3,3,3) =(-2*B0*B0*B0*B0*lp*lp*lp*X1c*X1c*Y1s*Y1s*&
       curvature*curvature)/(G0*G0*G0)
  
  if (debug) then
     ! Verify symmetry in the first 2 indices
     do i = 1, 3
        do j = 1, i-1
           do k = 1, 3
              difference = grad_grad_B(:,i,j,k) - grad_grad_B(:,j,i,k)
              norm_difference = maxval(abs(difference))
              if (verbose) print "(a,6(i1,a),es24.14)", " grad_grad_B(",i,",",j,",",k,") - grad_grad_B(",j,",",i,",",k,") = ",norm_difference
              if (norm_difference > 1e-8) then
                 print "(a,6(i1,a),es24.14)", "ERROR!!! grad_grad_B(",i,",",j,",",k,") - grad_grad_B(",j,",",i,",",k,") = ",norm_difference
                 stop
              end if
           end do
        end do
     end do
     
     ! If curl-free, verify symmetry in last 2 indices
     if (abs(I2) < 1e-13) then
        do i = 1, 3
           do j = 1, 3
              do k = 1, j-1
                 difference = grad_grad_B(:,i,j,k) - grad_grad_B(:,i,k,j)
                 norm_difference = maxval(abs(difference))
                 if (verbose) print "(a,6(i1,a),es24.14)", " grad_grad_B(",i,",",j,",",k,") - grad_grad_B(",i,",",k,",",j,") = ",norm_difference
                 if (norm_difference > 1e-8) then
                    print "(a,6(i1,a),es24.14)", "ERROR!! grad_grad_B(",i,",",j,",",k,") - grad_grad_B(",i,",",k,",",j,") = ",norm_difference
                    stop
                 end if
              end do
           end do
        end do
     end if
  end if

  deallocate(grad_grad_B)
  deallocate(d_curvature_d_zeta, d_torsion_d_zeta, d_X20_d_zeta, d_X2s_d_zeta, d_X2c_d_zeta, d_Y20_d_zeta, d_Y2s_d_zeta, d_Y2c_d_zeta)
  deallocate(d_Z20_d_zeta, d_Z2s_d_zeta, d_Z2c_d_zeta, d2_X1c_d_zeta2, d2_Y1c_d_zeta2, d2_Y1s_d_zeta2, difference)
  
end subroutine quasisymmetry_grad_grad_B_tensor

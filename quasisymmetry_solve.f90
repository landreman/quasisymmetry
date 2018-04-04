subroutine quasisymmetry_solve

  use quasisymmetry_variables

  implicit none

  integer :: iteration, j_line_search
  real(dp) :: residual_norm, last_residual_norm, step_scale
  real(dp), dimension(:), allocatable :: state, state0

  ! Variables needed by LAPACK:                                                                                            
  integer :: INFO
  integer, dimension(:), allocatable :: IPIV

  allocate(state(matrix_size))
  allocate(state0(matrix_size))
  allocate(IPIV(matrix_size))

  ! Initialize state
  sigma = 1
  iota = 0

  state(1:N_phi) = sigma
  state(matrix_size) = iota

  call quasisymmetry_residual()
  residual_norm = sum(residual * residual)
  print "(a,es10.3)","                 Initial residual norm:",residual_norm

  ! Here is the main Newton iteration:
  Newton: do iteration = 1, N_iterations
     last_residual_norm = residual_norm
     if (residual_norm < Newton_tolerance) exit Newton
     
     call quasisymmetry_Jacobian()

     state0 = state
     print "(a,i3)","  Newton iteration ",iteration
     ! We will use the LAPACK subroutine DGESV to solve a general (asymmetric) linear system
     ! step_direction = - matrix \ residual
     step_direction = -residual ! Note that LAPACK will over-write step_direction and with the solution, and over-write Jacobian with the LU factorization.
     call DGESV(matrix_size, 1, Jacobian, matrix_size, IPIV, step_direction, matrix_size, INFO)
     if (INFO /= 0) then
        print *, "Error in LAPACK call DGESV: info = ", INFO
        stop
     end if

     step_scale = 1
     line_search: do j_line_search = 1, N_line_search
        state = state0 + step_scale * step_direction

        sigma = state(1:N_phi)
        iota = state(matrix_size)
        call quasisymmetry_residual()
        residual_norm = sum(residual * residual)
        print "(a,i3,a,es10.3,a,es23.15)","    Line search step",j_line_search,"  Residual norm:",residual_norm,"  iota:",iota
        if (residual_norm < last_residual_norm) exit line_search

        step_scale = step_scale / 2
     end do line_search
  end do Newton
  ! End of Newton solve.
  ! Now compute quantities that are derived from the solution:

  Y1s = sign_G * curvature * ( B1c_over_B0 + B1s_over_B0 * sigma) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0);
  Y1c = sign_G * curvature * (-B1s_over_B0 + B1c_over_B0 * sigma) / (B1c_over_B0*B1c_over_B0 + B1s_over_B0*B1s_over_B0);
    
!!$  R1c = (RZ_to_XY_d * (X1c * normal_cylindrical(:,1) + Y1c * binormal_cylindrical(:,1)) &
!!$       - RZ_to_XY_b * (X1c * normal_cylindrical(:,3) + Y1c * binormal_cylindrical(:,3))) &
!!$       / (RZ_to_XY_a * RZ_to_XY_d - RZ_to_XY_b * RZ_to_XY_b)
!!$    
!!$  Z1c = (-RZ_to_XY_b * (X1c * normal_cylindrical(:,1) + Y1c * binormal_cylindrical(:,1)) &
!!$       + RZ_to_XY_a * (X1c * normal_cylindrical(:,3) + Y1c * binormal_cylindrical(:,3))) &
!!$       / (RZ_to_XY_a * RZ_to_XY_d - RZ_to_XY_b * RZ_to_XY_b)
!!$    
!!$  R1s = (RZ_to_XY_d * (X1s * normal_cylindrical(:,1) + Y1s * binormal_cylindrical(:,1)) &
!!$       - RZ_to_XY_b * (X1s * normal_cylindrical(:,3) + Y1s * binormal_cylindrical(:,3))) &
!!$       / (RZ_to_XY_a * RZ_to_XY_d - RZ_to_XY_b * RZ_to_XY_b)
!!$
!!$  Z1s = (-RZ_to_XY_b * (X1s * normal_cylindrical(:,1) + Y1s * binormal_cylindrical(:,1)) &
!!$       + RZ_to_XY_a * (X1s * normal_cylindrical(:,3) + Y1s * binormal_cylindrical(:,3))) &
!!$       / (RZ_to_XY_a * RZ_to_XY_d - RZ_to_XY_b * RZ_to_XY_b)

  R1c = (-binormal_cylindrical(:,3) * X1c + normal_cylindrical(:,3) * Y1c) * d_l_d_phi / R0
  R1s = (-binormal_cylindrical(:,3) * X1s + normal_cylindrical(:,3) * Y1s) * d_l_d_phi / R0
  Z1c = ( binormal_cylindrical(:,1) * X1c - normal_cylindrical(:,1) * Y1c) * d_l_d_phi / R0
  Z1s = ( binormal_cylindrical(:,1) * X1s - normal_cylindrical(:,1) * Y1s) * d_l_d_phi / R0

  call quasisymmetry_elongation()

  call quasisymmetry_determine_helicity()

  print *,"elongation:"
  print *,elongation
  print "(a,es23.15,a,es23.15)", " Final iota:",iota,"  max elongation:",max_elongation

  deallocate(state, state0, IPIV)

  ! Print R1c in MATLAB format, for comparing with the matlab version:
!!$  print *,"R1c_fortran=[",R1c(1),";"
!!$  do iteration = 2,N_phi-1
!!$     print *,R1c(iteration),";"
!!$  end do
!!$  print *,R1c(N_phi),"];"

end subroutine quasisymmetry_solve

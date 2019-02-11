module quasisymmetry_Frenet_to_cylindrical_mod

  ! The MATLAB version of this algorithm can be found in
  ! m20180501_01_quasisymmetricEquilibria_Frenet_altTransformation.m

  implicit none

  private

  public :: quasisymmetry_Frenet_to_cylindrical

contains

  subroutine quasisymmetry_Frenet_to_cylindrical

    use quasisymmetry_variables
    use vmec_input, only: vmec_nfp => nfp, lasym, ntor, RBC, RBS, ZBC, ZBS
    use vparams, only: ntord
    use quasisymmetry_splines

    implicit none

    integer :: N_theta, j_theta, N_phi_conversion, j_phi, i
    real(dp) :: costheta, sintheta, final_R, final_z
    real(dp), dimension(:,:), allocatable :: R_2D, z_2D, phi0_2D
    real(dp), dimension(:), allocatable :: theta, phi_conversion
    real(dp), dimension(:), allocatable :: X1, Y1
    type(periodic_spline) :: normal_R_spline, normal_phi_spline, normal_z_spline, binormal_R_spline, binormal_phi_spline, binormal_z_spline
    type(periodic_spline) :: X1_spline, Y1_spline, R0_spline, z0_spline
    real(dp) :: rootSolve_abserr, rootSolve_relerr, phi0_rootSolve_min, phi0_rootSolve_max
    real(dp) :: phi0_solution, phi_target
    integer :: fzeroFlag

    !----------------------------------------------

    N_theta = 28
    N_phi_conversion = N_phi
    allocate(theta(N_theta))
    allocate(phi_conversion(N_phi_conversion))
    allocate(phi0_2D(N_theta,N_phi_conversion))
    allocate(R_2D(N_theta,N_phi_conversion))
    allocate(z_2D(N_theta,N_phi_conversion))
    theta = [( 2*pi*i/N_theta, i=0,N_theta-1 )]
    phi_conversion = [( 2*pi*i/(nfp*N_phi_conversion), i=0,N_phi_conversion-1 )]

    call new_periodic_spline(N_phi, phi, R0, 2*pi/nfp, R0_spline)
    call new_periodic_spline(N_phi, phi, z0, 2*pi/nfp, z0_spline)
    call new_periodic_spline(N_phi, phi, normal_cylindrical(:,1), 2*pi/nfp, normal_R_spline)
    call new_periodic_spline(N_phi, phi, normal_cylindrical(:,2), 2*pi/nfp, normal_phi_spline)
    call new_periodic_spline(N_phi, phi, normal_cylindrical(:,3), 2*pi/nfp, normal_z_spline)
    call new_periodic_spline(N_phi, phi, binormal_cylindrical(:,1), 2*pi/nfp, binormal_R_spline)
    call new_periodic_spline(N_phi, phi, binormal_cylindrical(:,2), 2*pi/nfp, binormal_phi_spline)
    call new_periodic_spline(N_phi, phi, binormal_cylindrical(:,3), 2*pi/nfp, binormal_z_spline)

    rootSolve_abserr = 1.0e-10_dp
    rootSolve_relerr = 1.0e-10_dp

    do j_theta = 1, N_theta
       costheta = cos(theta(j_theta))
       sintheta = sin(theta(j_theta))
       X1 = X1c_untwisted * costheta + X1s_untwisted * sintheta
       Y1 = Y1c_untwisted * costheta + Y1s_untwisted * sintheta
       call new_periodic_spline(N_phi, phi, X1, 2*pi/nfp, X1_spline)
       call new_periodic_spline(N_phi, phi, Y1, 2*pi/nfp, Y1_spline)
       do j_phi = 1, N_phi_conversion
          ! Solve for the phi0 such that r0 + X1 n + Y1 b has the desired phi

          phi_target = phi_conversion(j_phi)
          phi0_rootSolve_min = phi_target - 0.3
          phi0_rootSolve_max = phi_target + 0.3
          call quasisymmetry_fzero(Frenet_to_cylindrical_residual, phi0_rootSolve_min, phi0_rootSolve_max, phi_target, &
               rootSolve_relerr, rootSolve_abserr, fzeroFlag)
          ! Note: fzero returns its answer in phi0_rootSolve_min
          phi0_solution = phi0_rootSolve_min
          if (fzeroFlag == 4) then
             print *,"j_theta=",j_theta," j_phi=",j_phi
             stop "ERROR: fzero returned error 4: no sign change in residual"
          else if (fzeroFlag > 2) then
             print *,"WARNING in cosm: fzero returned an error code:",fzeroFlag
          end if

          call Frenet_to_cylindrical_1_point(phi0_solution, final_R, final_z)
          R_2D(j_theta,j_phi) = final_R
          z_2D(j_theta,j_phi) = final_z
          phi0_2D(j_theta,j_phi) = phi0_solution
       end do
       call delete_periodic_spline(X1_spline)
       call delete_periodic_spline(Y1_spline)
    end do
    
!!$    print *,"N_theta:"
!!$    print *,N_theta
!!$    print *,"N_zeta:"
!!$    print *,N_phi_conversion
!!$    print *,"phi0_2D:"
!!$    do j_theta = 1,N_theta
!!$       !print "(*(f7.3))",phi0_2D(j_theta,:)
!!$       print "(*(f20.15))",phi0_2D(j_theta,:)
!!$    end do
!!$    print *,"R_2D:"
!!$    do j_theta = 1,N_theta
!!$       !print "(*(f7.3))",R_2D(j_theta,:)
!!$       print "(*(f20.15))",R_2D(j_theta,:)
!!$    end do
!!$    print *,"z_2D:"
!!$    do j_theta = 1,N_theta
!!$       !print "(*(f7.3))",Z_2D(j_theta,:)
!!$       print "(*(f20.15))",Z_2D(j_theta,:)
!!$    end do

    deallocate(theta,phi_conversion,R_2D,z_2D)
    call delete_periodic_spline(R0_spline)
    call delete_periodic_spline(z0_spline)
    call delete_periodic_spline(normal_R_spline)
    call delete_periodic_spline(normal_phi_spline)
    call delete_periodic_spline(normal_z_spline)
    call delete_periodic_spline(binormal_R_spline)
    call delete_periodic_spline(binormal_phi_spline)
    call delete_periodic_spline(binormal_z_spline)

  contains

    function Frenet_to_cylindrical_residual(phi0)
      ! Given a point on the axis with toroidal angle phi0, compute phi for the associated point at r>0,
      ! and find the difference between this phi and the target.

      implicit none

      real(dp) :: Frenet_to_cylindrical_residual, phi0
      real(dp) :: total_x, total_y, R0_at_phi0, X1_at_phi0, Y1_at_phi0
      real(dp) :: cosphi0, sinphi0
      real(dp) :: normal_R, normal_phi, normal_x, normal_y
      real(dp) :: binormal_R, binormal_phi, binormal_x, binormal_y
      
      sinphi0 = sin(phi0)
      cosphi0 = cos(phi0)
      R0_at_phi0   = periodic_splint(phi0,R0_spline)
      X1_at_phi0   = periodic_splint(phi0,X1_spline)
      Y1_at_phi0   = periodic_splint(phi0,Y1_spline)
      normal_R     = periodic_splint(phi0,normal_R_spline)
      normal_phi   = periodic_splint(phi0,normal_phi_spline)
      binormal_R   = periodic_splint(phi0,binormal_R_spline)
      binormal_phi = periodic_splint(phi0,binormal_phi_spline)

      normal_x   =   normal_R * cosphi0 -   normal_phi * sinphi0
      normal_y   =   normal_R * sinphi0 +   normal_phi * cosphi0
      binormal_x = binormal_R * cosphi0 - binormal_phi * sinphi0
      binormal_y = binormal_R * sinphi0 + binormal_phi * cosphi0

      total_x = R0_at_phi0 * cosphi0 + r * (X1_at_phi0 * normal_x + Y1_at_phi0 * binormal_x)
      total_y = R0_at_phi0 * sinphi0 + r * (X1_at_phi0 * normal_y + Y1_at_phi0 * binormal_y)
      Frenet_to_cylindrical_residual = atan2(total_y, total_x) - phi_target

    end function Frenet_to_cylindrical_residual

    ! -----------------------------------------------

    subroutine Frenet_to_cylindrical_1_point(phi0, total_R, total_z)
      ! Given a point on the axis with toroidal angle phi0, compute R and z for the associated point at r>0.

      implicit none

      real(dp), intent(in) :: phi0
      real(dp), intent(out) :: total_R, total_z
      real(dp) :: total_x, total_y, R0_at_phi0, z0_at_phi0, X1_at_phi0, Y1_at_phi0
      real(dp) :: cosphi0, sinphi0
      real(dp) :: normal_R, normal_phi, normal_x, normal_y, normal_z
      real(dp) :: binormal_R, binormal_phi, binormal_x, binormal_y, binormal_z
      
      sinphi0 = sin(phi0)
      cosphi0 = cos(phi0)
      R0_at_phi0   = periodic_splint(phi0,R0_spline)
      z0_at_phi0   = periodic_splint(phi0,z0_spline)
      X1_at_phi0   = periodic_splint(phi0,X1_spline)
      Y1_at_phi0   = periodic_splint(phi0,Y1_spline)
      normal_R     = periodic_splint(phi0,normal_R_spline)
      normal_phi   = periodic_splint(phi0,normal_phi_spline)
      normal_z     = periodic_splint(phi0,normal_z_spline)
      binormal_R   = periodic_splint(phi0,binormal_R_spline)
      binormal_phi = periodic_splint(phi0,binormal_phi_spline)
      binormal_z   = periodic_splint(phi0,binormal_z_spline)

      normal_x   =   normal_R * cosphi0 -   normal_phi * sinphi0
      normal_y   =   normal_R * sinphi0 +   normal_phi * cosphi0
      binormal_x = binormal_R * cosphi0 - binormal_phi * sinphi0
      binormal_y = binormal_R * sinphi0 + binormal_phi * cosphi0

      total_x = R0_at_phi0 * cosphi0 + r * (X1_at_phi0 * normal_x + Y1_at_phi0 * binormal_x)
      total_y = R0_at_phi0 * sinphi0 + r * (X1_at_phi0 * normal_y + Y1_at_phi0 * binormal_y)
      total_z = z0_at_phi0           + r * (X1_at_phi0 * normal_z + Y1_at_phi0 * binormal_z)
      total_R = sqrt(total_x * total_x + total_y * total_y)

    end subroutine Frenet_to_cylindrical_1_point

  end subroutine quasisymmetry_Frenet_to_cylindrical

end module quasisymmetry_Frenet_to_cylindrical_mod

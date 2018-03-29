subroutine quasisymmetry_init_phi

  use quasisymmetry_variables

  implicit none

  real(prec), dimension(:,:), allocatable :: temp_matrix
  real(prec), dimension(:,:), allocatable :: temp_vector
  integer :: option, quadrature_option

  if (allocated(phi)) deallocate(phi)
  if (allocated(phi_extended)) deallocate(phi_extended)
  if (allocated(d_d_phi)) deallocate(d_d_phi)

  if (allocated(X1s)) deallocate(X1s)
  if (allocated(X1c)) deallocate(X1c)
  if (allocated(Y1s)) deallocate(Y1s)
  if (allocated(Y1c)) deallocate(Y1c)
  if (allocated(sigma)) deallocate(sigma)
  if (allocated(elongation)) deallocate(elongation)
  
  allocate(phi(N_phi))
  allocate(d_d_phi(N_phi,N_phi))
  allocate(temp_matrix(N_phi,N_phi))
  allocate(temp_vector(N_phi))

  allocate(X1s(N_phi))
  allocate(X1c(N_phi))
  allocate(Y1s(N_phi))
  allocate(Y1c(N_phi))
  allocate(sigma(N_phi))
  allocate(elongation(N_phi))

  option = 20
  quadrature_option = 0
  call uniformDiffMatrices(N_phi,0_prec, 2*pi/nfp, option, quadrature_option, phi, temp_vector, d_d_phi, temp_matrix)

  phi_extended = [( 2*pi*i/(N_phi*nfp), i=0,N_phi*nfp-1 )]

  deallocate(temp_matrix, temp_vector)

end subroutine quasisymmetry_init_phi

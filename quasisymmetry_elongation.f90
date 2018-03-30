subroutine quasisymmetry_elongation

  use quasisymmetry_variables

  implicit none

  real(dp), dimension(:), allocatable :: p, q

  allocate(q(N_phi))
  allocate(p(N_phi))

  ! See my note 20180329-03 for derivation of the formula below for elongation:
  p = R1s*R1s + R1c*R1c + Z1s*Z1s + Z1c*Z1c
  q = R1s*Z1c - R1c*Z1s
  elongation = 2*abs(q) / (p - sqrt(p*p-4*q*q))

  ! Search for maximum using Fourier interpolation...

  deallocate(p,q)

end subroutine quasisymmetry_elongation

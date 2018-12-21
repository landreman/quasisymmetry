subroutine quasisymmetry_write_vmec_input

  use quasisymmetry_variables
  use safe_open_mod
  use vmec_input, only: vmec_nfp => nfp, lasym, ntor, raxis_cc, raxis_cs, zaxis_cc, zaxis_cs, &
       read_indata_namelist, write_indata_namelist, lfreeb, RBC, RBS, ZBC, ZBS
  use vparams, only: ntord

  implicit none

  integer :: iunit = 40, istat
  integer :: max_n, n
  real(dp) :: half_sum, half_difference

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read in the VMEC template namelist:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  istat = 0
  call safe_open(iunit,istat,trim(vmec_template_filename),'old','formatted')
  if (istat /= 0) then
     print "(a,a)"," Error opening vmec_template_filename: ",trim(vmec_template_filename)
     stop
  end if

  call read_indata_namelist(iunit, istat)
  if (istat /= 0) then
     print "(a,a)"," Error reading &indata namelist from vmec_template_filename: ",trim(vmec_template_filename)
     stop
  end if

  close(iunit)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Over-write a few vmec parameters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  vmec_nfp = nfp
  lfreeb = .false.

  ! The output is not stellarator-symmetric if (1) R0s is nonzero, (2) Z0c is nonzero, or (3) sigma_initial is nonzero:
  lasym = (maxval(abs(R0s))>0 .or. maxval(abs(Z0c)) > 0 .or. abs(sigma_initial) > 0)

  ! We should be able to resolve (N_phi-1)/2 modes (note integer division!), but in case N_phi is very large, don't attempt more than the vmec arrays can handle.
  ntor = min((N_phi - 1) / 2, ntord)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Over-write vmec axis shape
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  raxis_cc = 0
  raxis_cs = 0
  zaxis_cc = 0
  zaxis_cs = 0

  max_n = min(ntord,max_axis_nmax)

  ! To convert sin(...) modes to vmec, we introduce a minus sign. This is because in vmec,
  ! R and Z ~ sin(m theta - n phi), which for m=0 is sin(-n phi) = -sin(n phi).

  raxis_cc(0:max_n) = R0c(1:max_n+1)
  raxis_cs(0:max_n) = -R0s(1:max_n+1)
  zaxis_cc(0:max_n) = Z0c(1:max_n+1)
  zaxis_cs(0:max_n) = -Z0s(1:max_n+1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute new RBC, RBS, ZBC, ZBS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RBC = 0
  RBS = 0
  ZBC = 0
  ZBS = 0

  ! Handle the m=0 modes of the boundary, which are the same as the axis shape:
  RBC(0:max_n,0) = raxis_cc(0:max_n)
  RBS(0:max_n,0) = raxis_cs(0:max_n)
  ZBC(0:max_n,0) = zaxis_cc(0:max_n)
  ZBS(0:max_n,0) = zaxis_cs(0:max_n)

  ! Handle the n=0 m=1 modes:
  RBC(0,1) = r * sum(R1c) / N_phi
  RBS(0,1) = r * sum(R1s) / N_phi
  ZBC(0,1) = r * sum(Z1c) / N_phi
  ZBS(0,1) = r * sum(Z1s) / N_phi

  ! Handle the m=1 modes that have nonzero n:
  do n = 1, ntor
     ! RBC:
     half_sum        =  r * sum(R1c * cos_n_phi(:,n+1)) / N_phi
     half_difference =  r * sum(R1s * sin_n_phi(:,n+1)) / N_phi
     RBC( n,1) = half_sum + half_difference
     RBC(-n,1) = half_sum - half_difference

     ! ZBC:
     half_sum        =  r * sum(Z1c * cos_n_phi(:,n+1)) / N_phi
     half_difference =  r * sum(Z1s * sin_n_phi(:,n+1)) / N_phi
     ZBC( n,1) = half_sum + half_difference
     ZBC(-n,1) = half_sum - half_difference

     ! RBS:
     half_sum        =  r * sum(R1s * cos_n_phi(:,n+1)) / N_phi
     half_difference = -r * sum(R1c * sin_n_phi(:,n+1)) / N_phi
     RBS( n,1) = half_sum + half_difference
     RBS(-n,1) = half_sum - half_difference

     ! ZBS:
     half_sum        =  r * sum(Z1s * cos_n_phi(:,n+1)) / N_phi
     half_difference = -r * sum(Z1c * sin_n_phi(:,n+1)) / N_phi
     ZBS( n,1) = half_sum + half_difference
     ZBS(-n,1) = half_sum - half_difference
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  istat = 0
  call safe_open(iunit,istat,trim(new_vmec_filename),'unknown','formatted')
  if (istat /= 0) then
     print "(a,a)"," Error opening new_vmec_filename: ",trim(new_vmec_filename)
     stop
  end if

  write(iunit,'(a)') "! This &INDATA namelist was generated by quasisymmetry.f90"
  write(iunit,'(a,a)') "! Based on template file ",trim(vmec_template_filename)
  write(iunit,'(a,es22.12)') "! r =",r

  call write_indata_namelist(iunit, istat)
  if (istat /= 0) then
     print "(a,a)"," Error writing &indata namelist to new_vmec_filename: ",trim(new_vmec_filename)
     stop
  end if

  close(iunit)

end subroutine quasisymmetry_write_vmec_input

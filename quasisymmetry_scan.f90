subroutine quasisymmetry_scan

  use quasisymmetry_variables

  implicit none

  integer :: j_scan, j_Fourier_scan, j, k, N_Fourier_scan, j_B1s, j_B1c
  integer, dimension(max_axis_nmax+1, 4) :: scan_state, N_scan_array

  ! Clean up scan arrays
  do j = 1,max_axis_nmax+1
     if (R0s_N_scan(j)<1) R0s_N_scan(j) = 1
     if (R0c_N_scan(j)<1) R0c_N_scan(j) = 1
     if (Z0s_N_scan(j)<1) Z0s_N_scan(j) = 1
     if (Z0c_N_scan(j)<1) Z0c_N_scan(j) = 1

     N_scan_array(j,1) = max(R0s_N_scan(j),1)
     N_scan_array(j,2) = max(R0c_N_scan(j),1)
     N_scan_array(j,3) = max(Z0s_N_scan(j),1)
     N_scan_array(j,4) = max(Z0c_N_scan(j),1)
  end do
  if (B1s_N_scan<1) B1s_N_scan = 1
  if (B1c_N_scan<1) B1c_N_scan = 1

  !N_scan = product(R0s_N_scan)*product(R0c_N_scan)*product(Z0s_N_scan)*product(Z0c_N_scan)*product(B1s_N_scan)*product(B1c_N_scan)
  N_Fourier_scan = product(N_scan_array)
  N_scan = N_Fourier_scan * B1s_N_scan * B1c_N_scan
  print *,"N_Fourier_scan:",N_Fourier_scan,"N_scan:",N_scan
  allocate(iotas(N_scan))
  allocate(max_elongations(N_scan))

  scan_state = 1
  j_scan = 0
  do j_Fourier_scan = 1, N_Fourier_scan
!!$     print *,"scan_state:"
!!$     do k = 1,4
!!$        print *,scan_state(:,k)
!!$     end do

     ! Set Fourier amplitudes for axis:
     do j = 1, max_axis_nmax+1
        R0s(j) = R0s_min(j) + (R0s_max(j) - R0s_min(j)) * (scan_state(j,1) - 1) / max(N_scan_array(j,1) - 1, 1)
        R0c(j) = R0c_min(j) + (R0c_max(j) - R0c_min(j)) * (scan_state(j,2) - 1) / max(N_scan_array(j,2) - 1, 1)
        Z0s(j) = Z0s_min(j) + (Z0s_max(j) - Z0s_min(j)) * (scan_state(j,3) - 1) / max(N_scan_array(j,3) - 1, 1)
        Z0c(j) = Z0c_min(j) + (Z0c_max(j) - Z0c_min(j)) * (scan_state(j,4) - 1) / max(N_scan_array(j,4) - 1, 1)
     end do

     do j_B1s = 1, B1s_N_scan
        B1s_over_B0 = B1s_min + (B1s_max - B1s_min) * (j_B1s - 1) / max(B1s_N_scan - 1, 1)
        do j_B1c = 1, B1c_N_scan
           B1c_over_B0 = B1c_min + (B1c_max - B1c_min) * (j_B1c - 1) / max(B1c_N_scan - 1, 1)
           j_scan = j_scan + 1
           print *,"###################################################################################"
           print "(a,i5,a,i5)"," Scan case",j_scan," of",N_scan
           print *,"B1s =", B1s_over_B0, "   B1c =",B1c_over_B0
           print *,"R0s:",R0s
           print *,"R0c:",R0c
           print *,"Z0s:",Z0s
           print *,"Z0c:",Z0c
        end do
     end do

     ! Update scan state for the next solve
     k = 1
     !do k = 1,4
     do while (k <= 4)
        j = 1
        !do j = 1,max_axis_nmax+1
        do while (j <= max_axis_nmax+1)
           if (scan_state(j,k) < N_scan_array(j,k)) then
              scan_state(j,k) = scan_state(j,k) + 1
              ! exit both j and k loops
              k = 4
              j = max_axis_nmax+1
           else
              scan_state(j,k) = 1
           end if
           j = j + 1
        end do
        k = k + 1
     end do
  end do

end subroutine quasisymmetry_scan

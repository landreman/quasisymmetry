subroutine quasisymmetry_random

  ! N_random = maximum # of successful configurations to keep on each MPI processor.
  ! random_time = maximum number of seconds to run.
  ! This subroutine will exit when either of the above limits is reached.

  use quasisymmetry_variables

  implicit none

  include 'mpif.h'

  integer :: j, k
  integer*8 :: j_scan, index
  integer :: ierr, tag
  integer*8 :: dummy_long(1), print_summary_stride
  integer :: mpi_status(MPI_STATUS_SIZE)
  real(dp), dimension(:), allocatable :: iotas_local, max_elongations_local, mean_elongations_local, rms_curvatures_local, max_curvatures_local, axis_lengths_local
  real(dp), dimension(:), allocatable :: standard_deviations_of_R_local, standard_deviations_of_Z_local, max_modBinv_sqrt_half_grad_B_colon_grad_Bs_local
  integer, dimension(:), allocatable :: axis_helicities_local
  logical, dimension(:), allocatable :: iota_tolerance_achieveds_local, elongation_tolerance_achieveds_local, Newton_tolerance_achieveds_local
  integer*8, dimension(:), allocatable :: N_solves_kept
  real(dp), dimension(:), allocatable :: scan_eta_bar_local, scan_sigma_initial_local, scan_B2s_local, scan_B2c_local
  real(dp), dimension(:,:), allocatable :: scan_R0c_local, scan_R0s_local, scan_Z0c_local, scan_Z0s_local
  real(dp), dimension(:), allocatable :: max_B2tildes_local, r_singularities_local, d2_volume_d_psi2s_local, B20_variations_local
  logical :: keep_going
  real(dp) :: thresh, time1, time2
  real(dp) :: rand
  real(dp) :: log_eta_bar_min, log_eta_bar_max, log_sigma_initial_min, log_sigma_initial_max
  real(dp), dimension(max_axis_nmax+1) :: log_R0s_min, log_R0s_max, log_R0c_min, log_R0c_max, log_Z0s_min, log_Z0s_max, log_Z0c_min, log_Z0c_max
  real(dp) :: log_B2s_min, log_B2s_max, log_B2c_min, log_B2c_max

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! So initial lines printed by proc0 are sure to come first.

  ! Allocate arrays to store the scan results:
  allocate(iotas_local(N_random))
  allocate(max_elongations_local(N_random))
  allocate(mean_elongations_local(N_random))
  allocate(rms_curvatures_local(N_random))
  allocate(max_curvatures_local(N_random))
  allocate(max_modBinv_sqrt_half_grad_B_colon_grad_Bs_local(N_random))
  allocate(axis_lengths_local(N_random))
  allocate(standard_deviations_of_R_local(N_random))
  allocate(standard_deviations_of_Z_local(N_random))
  allocate(axis_helicities_local(N_random))
  allocate(iota_tolerance_achieveds_local(N_random))
  allocate(elongation_tolerance_achieveds_local(N_random))
  allocate(Newton_tolerance_achieveds_local(N_random))
  if (trim(order_r_option) == order_r_option_r1_compute_B2) then
     allocate(max_B2tildes_local(N_random))
  end if
  if (trim(order_r_option) == order_r_option_r2) then
     allocate(r_singularities_local(N_random))
     allocate(d2_volume_d_psi2s_local(N_random))
     allocate(B20_variations_local(N_random))
  end if

  ! Allocate arrays to store the input parameters that correspond to saved results:
  allocate(scan_eta_bar_local(N_random))
  allocate(scan_sigma_initial_local(N_random))
  allocate(scan_R0c_local(N_random,axis_nmax+1))
  allocate(scan_R0s_local(N_random,axis_nmax+1))
  allocate(scan_Z0c_local(N_random,axis_nmax+1))
  allocate(scan_Z0s_local(N_random,axis_nmax+1))
  if (trim(order_r_option) == order_r_option_r2) then
     allocate(scan_B2s_local(N_random))
     allocate(scan_B2c_local(N_random))
  end if

  log_eta_bar_min = log(abs(eta_bar_min))
  log_eta_bar_max = log(abs(eta_bar_max))
  log_sigma_initial_min = log(abs(sigma_initial_min))
  log_sigma_initial_max = log(abs(sigma_initial_max))
  log_R0s_min = log(abs(R0s_min))
  log_R0s_max = log(abs(R0s_max))
  log_R0c_min = log(abs(R0c_min))
  log_R0c_max = log(abs(R0c_max))
  log_Z0s_min = log(abs(Z0s_min))
  log_Z0s_max = log(abs(Z0s_max))
  log_Z0c_min = log(abs(Z0c_min))
  log_Z0c_max = log(abs(Z0c_max))
  log_B2s_min = log(abs(B2s_min))
  log_B2s_max = log(abs(B2s_max))
  log_B2c_min = log(abs(B2c_min))
  log_B2c_max = log(abs(B2c_max))

  j_scan = 0
  keep_going = .true.
  do while (keep_going)

     ! If needed, pick a random value for eta_bar.
     if (eta_bar_max > eta_bar_min) then
        call random_number(rand)
        select case (trim(eta_bar_scan_option))
        case (eta_bar_scan_option_linear)
           eta_bar = eta_bar_min + rand * (eta_bar_max - eta_bar_min)
        case (eta_bar_scan_option_log)
           eta_bar = exp(log_eta_bar_min + rand * (log_eta_bar_max - log_eta_bar_min))
        case (eta_bar_scan_option_2_sided_log)
           eta_bar = exp(log_eta_bar_min + rand * (log_eta_bar_max - log_eta_bar_min))
           call random_number(rand)
           if (rand > 0.5) eta_bar = -eta_bar ! Flip sign with 50% probability.
        case default
           stop "Error! Unrecognized eta_bar_scan_option"
        end select
     else
        eta_bar = eta_bar_max
     end if

     ! If needed, pick a random value for sigma_initial.
     if (sigma_initial_max > sigma_initial_min) then
        call random_number(rand)
        select case (trim(sigma_initial_scan_option))
        case (sigma_initial_scan_option_linear)
           sigma_initial = sigma_initial_min + rand * (sigma_initial_max - sigma_initial_min)
        case (sigma_initial_scan_option_log)
           sigma_initial = exp(log_sigma_initial_min + rand * (log_sigma_initial_max - log_sigma_initial_min))
        case (sigma_initial_scan_option_2_sided_log)
           sigma_initial = exp(log_sigma_initial_min + rand * (log_sigma_initial_max - log_sigma_initial_min))
           call random_number(rand)
           if (rand > 0.5) sigma_initial = -sigma_initial ! Flip sign with 50% probability.
        case default
           stop "Error! Unrecognized sigma_initial_scan_option"
        end select
     else
        sigma_initial = sigma_initial_max
     end if

     if (trim(order_r_option) == order_r_option_r2) then
        ! If needed, pick a random value for B2s.
        if (B2s_max > B2s_min) then
           call random_number(rand)
           select case (trim(B2s_scan_option))
           case (B2s_scan_option_linear)
              B2s = B2s_min + rand * (B2s_max - B2s_min)
           case (B2s_scan_option_log)
              B2s = exp(log_B2s_min + rand * (log_B2s_max - log_B2s_min))
           case (B2s_scan_option_2_sided_log)
              B2s = exp(log_B2s_min + rand * (log_B2s_max - log_B2s_min))
              call random_number(rand)
              if (rand > 0.5) B2s = -B2s ! Flip sign with 50% probability.
           case default
              stop "Error! Unrecognized B2s_scan_option"
           end select
        else
           B2s = B2s_max
        end if
        
        ! If needed, pick a random value for B2c.
        if (B2c_max > B2c_min) then
           call random_number(rand)
           select case (trim(B2c_scan_option))
           case (B2c_scan_option_linear)
              B2c = B2c_min + rand * (B2c_max - B2c_min)
           case (B2c_scan_option_log)
              B2c = exp(log_B2c_min + rand * (log_B2c_max - log_B2c_min))
           case (B2c_scan_option_2_sided_log)
              B2c = exp(log_B2c_min + rand * (log_B2c_max - log_B2c_min))
              call random_number(rand)
              if (rand > 0.5) B2c = -B2c ! Flip sign with 50% probability.
           case default
              stop "Error! Unrecognized B2c_scan_option"
           end select
        else
           B2c = B2c_max
        end if
     end if

     ! Set Fourier amplitudes for axis:
     select case (trim(Fourier_scan_option))
     case (Fourier_scan_option_linear)
        do j = 1, axis_nmax+1
           if (R0s_max(j) > R0s_min(j)) then
              call random_number(rand)
              R0s(j) = R0s_min(j) + (R0s_max(j) - R0s_min(j)) * rand
           else
              R0s(j) = R0s_max(j)
           end if
           if (R0c_max(j) > R0c_min(j)) then
              call random_number(rand)
              R0c(j) = R0c_min(j) + (R0c_max(j) - R0c_min(j)) * rand
           else
              R0c(j) = R0c_max(j)
           end if
           if (Z0s_max(j) > Z0s_min(j)) then
              call random_number(rand)
              Z0s(j) = Z0s_min(j) + (Z0s_max(j) - Z0s_min(j)) * rand
           else
              Z0s(j) = Z0s_max(j)
           end if
           if (Z0c_max(j) > Z0c_min(j)) then
              call random_number(rand)
              Z0c(j) = Z0c_min(j) + (Z0c_max(j) - Z0c_min(j)) * rand
           else
              Z0c(j) = Z0c_max(j)
           end if
        end do

     case (Fourier_scan_option_2_sided_log, Fourier_scan_option_2_sided_log_except_Z0s1)
        do j = 1, axis_nmax+1
           if (R0s_max(j) > R0s_min(j)) then
              call random_number(rand)
              R0s(j) = exp(log_R0s_min(j) + (log_R0s_max(j) - log_R0s_min(j)) * rand)
              call random_number(rand)
              if (rand > 0.5) R0s(j) = -R0s(j) ! Flip sign with 50% probability.
           else
              R0s(j) = R0s_max(j)
           end if
           if (R0c_max(j) > R0c_min(j)) then
              call random_number(rand)
              R0c(j) = exp(log_R0c_min(j) + (log_R0c_max(j) - log_R0c_min(j)) * rand)
              call random_number(rand)
              if (rand > 0.5) R0c(j) = -R0c(j) ! Flip sign with 50% probability.
           else
              R0c(j) = R0c_max(j)
           end if
           if (Z0s_max(j) > Z0s_min(j)) then
              call random_number(rand)
              Z0s(j) = exp(log_Z0s_min(j) + (log_Z0s_max(j) - log_Z0s_min(j)) * rand)
              if (j>2 .or. trim(Fourier_scan_option)==Fourier_scan_option_2_sided_log) then
                 call random_number(rand)
                 if (rand > 0.5) Z0s(j) = -Z0s(j) ! Flip sign with 50% probability.
              end if
           else
              Z0s(j) = Z0s_max(j)
           end if
           if (Z0c_max(j) > Z0c_min(j)) then
              call random_number(rand)
              Z0c(j) = exp(log_Z0c_min(j) + (log_Z0c_max(j) - log_Z0c_min(j)) * rand)
              call random_number(rand)
              if (rand > 0.5) Z0c(j) = -Z0c(j) ! Flip sign with 50% probability.
           else
              Z0c(j) = Z0c_max(j)
           end if
        end do
     case default
        stop "Error! Unrecognized Fourier_scan_option"
     end select

     if (verbose) then
        print "(a)"," ###################################################################################"
        print "(a)"," Random scan: trying the following parameters"
        print *,"eta_bar =",eta_bar
        print *,"sigma_initial = ",sigma_initial
        print *,"R0s:",R0s(1:axis_nmax+1)
        print *,"R0c:",R0c(1:axis_nmax+1)
        print *,"Z0s:",Z0s(1:axis_nmax+1)
        print *,"Z0c:",Z0c(1:axis_nmax+1)
        if (trim(order_r_option) == order_r_option_r2) then
           print *,"B2s =",B2s,", B2c = ",B2c
        end if
     end if
           
     call quasisymmetry_single_solve()

     ! Check whether one of the stopping criteria has been reached.
     ! We do this before all the 'cycles' in case most of the parameters are rejected
     call cpu_time(time1)
     if (time1 - start_time >= random_time) keep_going = .false.

     ! If we can reject this solution, go back and pick new input parameters.
     if (skipped_solve) cycle ! In case R0 <= 0 or some other reason caused quasisymmetry_single_solve to exit prematurely.
     if (max_elongation > max_elongation_to_keep) cycle
     if (abs(iota) < min_iota_to_keep) cycle
     if (max_modBinv_sqrt_half_grad_B_colon_grad_B > max_max_modBinv_sqrt_half_grad_B_colon_grad_B_to_keep) cycle
     if (trim(order_r_option) == order_r_option_r1_compute_B2) then
        if (max_B2tilde > max_B2tilde_to_keep) cycle
     end if
     if (trim(order_r_option) == order_r_option_r2) then
        if (r_singularity < min_r_singularity_to_keep) cycle
        if (B20_variation > max_B20_variation_to_keep) cycle
     end if
           
     ! If we made it this far, then record the results
     j_scan = j_scan + 1
     
     scan_eta_bar_local(j_scan) = eta_bar
     scan_sigma_initial_local(j_scan) = sigma_initial
     scan_R0c_local(j_scan,:) = R0c(1:axis_nmax+1)
     scan_R0s_local(j_scan,:) = R0s(1:axis_nmax+1)
     scan_Z0c_local(j_scan,:) = Z0c(1:axis_nmax+1)
     scan_Z0s_local(j_scan,:) = Z0s(1:axis_nmax+1)
     if (trim(order_r_option) == order_r_option_r2) then
        scan_B2s_local(j_scan) = B2s
        scan_B2c_local(j_scan) = B2c
     end if

     iotas_local(j_scan) = iota
     max_elongations_local(j_scan) = max_elongation
     mean_elongations_local(j_scan) = mean_elongation
     rms_curvatures_local(j_scan) = rms_curvature
     max_curvatures_local(j_scan) = max_curvature
     max_modBinv_sqrt_half_grad_B_colon_grad_Bs_local(j_scan) = max_modBinv_sqrt_half_grad_B_colon_grad_B
     axis_lengths_local(j_scan) = axis_length
     standard_deviations_of_R_local(j_scan) = standard_deviation_of_R
     standard_deviations_of_Z_local(j_scan) = standard_deviation_of_Z
     axis_helicities_local(j_scan) = axis_helicity
     iota_tolerance_achieveds_local(j_scan) = iota_tolerance_achieved
     elongation_tolerance_achieveds_local(j_scan) = elongation_tolerance_achieved
     Newton_tolerance_achieveds_local(j_scan) = Newton_tolerance_achieved
     if (trim(order_r_option) == order_r_option_r1_compute_B2) then
        max_B2tildes_local(j_scan) = max_B2tilde
     end if
     if (trim(order_r_option) == order_r_option_r2) then
        r_singularities_local(j_scan) = r_singularity
        d2_volume_d_psi2s_local(j_scan) = d2_volume_d_psi2
        B20_variations_local(j_scan) = B20_variation
     end if

     ! Check whether one of the stopping criteria has been reached
     if (j_scan >= N_random) keep_going = .false.
  end do

  call cpu_time(time1)
  print "(a,i5,a,es10.3,a)"," Proc",mpi_rank," finished after",time1 - start_time," sec."

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! So proc0 does not start printing results until all procs have finished.

  call cpu_time(time1)

  ! Finally, send all results to proc 0.
  if (proc0) then
     if (N_procs > 1) then
        print "(a)"," ###################################################################################"
        print *,"Transferring results to proc 0"
     end if

     allocate(N_solves_kept(N_procs))
     N_solves_kept(1) = j_scan
     do j = 1, N_procs-1
        ! Use mpi_rank as the tag
        call mpi_recv(dummy_long,1,MPI_INTEGER8,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        N_solves_kept(j+1) = dummy_long(1)
     end do
     print *,"# solves kept on each proc:",N_solves_kept
     N_scan = sum(N_solves_kept)
     print *,"Total # of solves kept:",N_scan

     ! Now that we know the total number of runs that were kept, we can allocate the arrays for the final results:
     allocate(iotas(N_scan))
     allocate(max_elongations(N_scan))
     allocate(mean_elongations(N_scan))
     allocate(rms_curvatures(N_scan))
     allocate(max_curvatures(N_scan))
     allocate(max_modBinv_sqrt_half_grad_B_colon_grad_Bs(N_scan))
     allocate(axis_lengths(N_scan))
     allocate(standard_deviations_of_R(N_scan))
     allocate(standard_deviations_of_Z(N_scan))
     allocate(axis_helicities(N_scan))
     allocate(iota_tolerance_achieveds(N_scan))
     allocate(elongation_tolerance_achieveds(N_scan))
     allocate(Newton_tolerance_achieveds(N_scan))
     if (trim(order_r_option) == order_r_option_r1_compute_B2) then
        allocate(max_B2tildes(N_scan))
     end if
     if (trim(order_r_option) == order_r_option_r2) then
        allocate(r_singularities(N_scan))
        allocate(d2_volume_d_psi2s(N_scan))
        allocate(B20_variations(N_scan))
     end if

     allocate(scan_eta_bar(N_scan))
     allocate(scan_sigma_initial(N_scan))
     allocate(scan_R0c(N_scan,axis_nmax+1))
     allocate(scan_R0s(N_scan,axis_nmax+1))
     allocate(scan_Z0c(N_scan,axis_nmax+1))
     allocate(scan_Z0s(N_scan,axis_nmax+1))
     if (trim(order_r_option) == order_r_option_r2) then
        allocate(scan_B2s(N_scan))
        allocate(scan_B2c(N_scan))
     end if

     ! Store results from proc0 in the final arrays:
     iotas(1:N_solves_kept(1)) = iotas_local(1:N_solves_kept(1))
     max_elongations(1:N_solves_kept(1)) = max_elongations_local(1:N_solves_kept(1))
     mean_elongations(1:N_solves_kept(1)) = mean_elongations_local(1:N_solves_kept(1))
     rms_curvatures(1:N_solves_kept(1)) = rms_curvatures_local(1:N_solves_kept(1))
     max_curvatures(1:N_solves_kept(1)) = max_curvatures_local(1:N_solves_kept(1))
     max_modBinv_sqrt_half_grad_B_colon_grad_Bs(1:N_solves_kept(1)) = max_modBinv_sqrt_half_grad_B_colon_grad_Bs_local(1:N_solves_kept(1))
     axis_lengths(1:N_solves_kept(1)) = axis_lengths_local(1:N_solves_kept(1))
     standard_deviations_of_R(1:N_solves_kept(1)) = standard_deviations_of_R_local(1:N_solves_kept(1))
     standard_deviations_of_Z(1:N_solves_kept(1)) = standard_deviations_of_Z_local(1:N_solves_kept(1))
     axis_helicities(1:N_solves_kept(1)) = axis_helicities_local(1:N_solves_kept(1))
     iota_tolerance_achieveds(1:N_solves_kept(1)) = iota_tolerance_achieveds_local(1:N_solves_kept(1))
     elongation_tolerance_achieveds(1:N_solves_kept(1)) = elongation_tolerance_achieveds_local(1:N_solves_kept(1))
     Newton_tolerance_achieveds(1:N_solves_kept(1)) = Newton_tolerance_achieveds_local(1:N_solves_kept(1))
     if (trim(order_r_option) == order_r_option_r1_compute_B2) then
        max_B2tildes(1:N_solves_kept(1)) = max_B2tildes_local(1:N_solves_kept(1))
     end if
     if (trim(order_r_option) == order_r_option_r2) then
        r_singularities(1:N_solves_kept(1)) = r_singularities_local(1:N_solves_kept(1))
        d2_volume_d_psi2s(1:N_solves_kept(1)) = d2_volume_d_psi2s_local(1:N_solves_kept(1))
        B20_variations(1:N_solves_kept(1)) = B20_variations_local(1:N_solves_kept(1))
     end if

     scan_eta_bar(1:N_solves_kept(1)) = scan_eta_bar_local(1:N_solves_kept(1))
     scan_sigma_initial(1:N_solves_kept(1)) = scan_sigma_initial_local(1:N_solves_kept(1))
     scan_R0c(1:N_solves_kept(1),:) = scan_R0c_local(1:N_solves_kept(1),:)
     scan_R0s(1:N_solves_kept(1),:) = scan_R0s_local(1:N_solves_kept(1),:)
     scan_Z0c(1:N_solves_kept(1),:) = scan_Z0c_local(1:N_solves_kept(1),:)
     scan_Z0s(1:N_solves_kept(1),:) = scan_Z0s_local(1:N_solves_kept(1),:)
     if (trim(order_r_option) == order_r_option_r2) then
        scan_B2s(1:N_solves_kept(1)) = scan_B2s_local(1:N_solves_kept(1))
        scan_B2c(1:N_solves_kept(1)) = scan_B2c_local(1:N_solves_kept(1))
     end if

     index = N_solves_kept(1) + 1
     do j = 1, N_procs-1
        print "(a,i20,a,i4)"," Proc 0 is receiving results from ",N_solves_kept(j+1)," solves on proc",j
        call mpi_recv(iotas(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(max_elongations(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(mean_elongations(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(rms_curvatures(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(max_curvatures(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(max_modBinv_sqrt_half_grad_B_colon_grad_Bs(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(axis_lengths(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(standard_deviations_of_R(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(standard_deviations_of_Z(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(axis_helicities(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(Newton_tolerance_achieveds(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(iota_tolerance_achieveds(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(elongation_tolerance_achieveds(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_INT,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        if (trim(order_r_option) == order_r_option_r1_compute_B2) then
           call mpi_recv(max_B2tildes(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        end if
        if (trim(order_r_option) == order_r_option_r2) then
           call mpi_recv(r_singularities(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
           call mpi_recv(d2_volume_d_psi2s(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
           call mpi_recv(B20_variations(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        end if

        call mpi_recv(scan_eta_bar(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_sigma_initial(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_R0c(index:index+N_solves_kept(j+1)-1,:),N_solves_kept(j+1)*(axis_nmax+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_R0s(index:index+N_solves_kept(j+1)-1,:),N_solves_kept(j+1)*(axis_nmax+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_Z0c(index:index+N_solves_kept(j+1)-1,:),N_solves_kept(j+1)*(axis_nmax+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        call mpi_recv(scan_Z0s(index:index+N_solves_kept(j+1)-1,:),N_solves_kept(j+1)*(axis_nmax+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        if (trim(order_r_option) == order_r_option_r2) then
           call mpi_recv(scan_B2s(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
           call mpi_recv(scan_B2c(index:index+N_solves_kept(j+1)-1),N_solves_kept(j+1),MPI_DOUBLE,j,j,MPI_COMM_WORLD,mpi_status,ierr)
        end if

        index = index + N_solves_kept(j+1)
     end do


  else
     ! Send the number of runs this proc kept:
     dummy_long(1) = j_scan
     ! Use mpi_rank as the tag
     call mpi_send(dummy_long,1,MPI_INTEGER8,0,mpi_rank,MPI_COMM_WORLD,ierr)

     ! Send the other results:
     call mpi_send(iotas_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(max_elongations_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(mean_elongations_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(rms_curvatures_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(max_curvatures_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(max_modBinv_sqrt_half_grad_B_colon_grad_Bs_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(axis_lengths_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(standard_deviations_of_R_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(standard_deviations_of_Z_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(axis_helicities_local(1:j_scan),j_scan,MPI_INT,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(Newton_tolerance_achieveds_local(1:j_scan),j_scan,MPI_LOGICAL,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(iota_tolerance_achieveds_local(1:j_scan),j_scan,MPI_LOGICAL,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(elongation_tolerance_achieveds_local(1:j_scan),j_scan,MPI_LOGICAL,0,mpi_rank,MPI_COMM_WORLD,ierr)
     if (trim(order_r_option) == order_r_option_r1_compute_B2) then
        call mpi_send(max_B2tildes_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     end if
     if (trim(order_r_option) == order_r_option_r2) then
        call mpi_send(r_singularities_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
        call mpi_send(d2_volume_d_psi2s_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
        call mpi_send(B20_variations_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     end if

     call mpi_send(scan_eta_bar_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_sigma_initial_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_R0c_local(1:j_scan,:),j_scan*(axis_nmax+1),MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_R0s_local(1:j_scan,:),j_scan*(axis_nmax+1),MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_Z0c_local(1:j_scan,:),j_scan*(axis_nmax+1),MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     call mpi_send(scan_Z0s_local(1:j_scan,:),j_scan*(axis_nmax+1),MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     if (trim(order_r_option) == order_r_option_r2) then
        call mpi_send(scan_B2s_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
        call mpi_send(scan_B2c_local(1:j_scan),j_scan,MPI_DOUBLE,0,mpi_rank,MPI_COMM_WORLD,ierr)
     end if
  end if

  if (proc0) then
     call cpu_time(time2)
     print "(a,es10.3,a)"," Time for communication:",time2 - time1," sec."
     print "(a)"," ###################################################################################"
     print *,"Scan complete."
     if (N_scan < 5000) then
        print "(a,99999(f8.2))"," iotas:",iotas
        print *," "
        print "(a,99999(f8.2))"," eta_bar:",scan_eta_bar
        print *," "
        print "(a,99999(f8.1))"," elongations:",max_elongations
        print *," "
        print "(a,99999(f8.2))"," rms_curvatures:",rms_curvatures
        print *," "
        print "(a,99999(f8.2))"," max_curvatures:",max_curvatures
        print *," "
        print "(a,99999(f8.2))"," max_modBinv_sqrt_half_grad_B_colon_grad_Bs:",max_modBinv_sqrt_half_grad_B_colon_grad_Bs
        print *," "
        if (trim(order_r_option) == order_r_option_r1_compute_B2) then
           print "(a,99999(f8.2))"," max_B2tildes:",max_B2tildes
           print *," "
        end if
        if (trim(order_r_option) == order_r_option_r2) then
           print "(a,99999(f8.2))"," r_singularities:",r_singularities
           print *," "
           print "(a,99999(f8.2))"," B20_variations:",B20_variations
           print *," "
        end if
        print "(a,99999(f8.2))"," axis_lengths:",axis_lengths
        print *," "
        print "(a,99999(i2))","                axis_helicities:",axis_helicities
        print *,"    Newton_tolerance_achieveds:",Newton_tolerance_achieveds
        print *,"      iota_tolerance_achieveds:",iota_tolerance_achieveds
        print *,"elongation_tolerance_achieveds:",elongation_tolerance_achieveds
     end if
  end if

end subroutine quasisymmetry_random


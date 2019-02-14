! Main program

program quasisymmetry

  use quasisymmetry_variables, only: start_time, total_time, general_option, general_option_single, general_option_scan, &
       N_procs, mpi_rank, proc0

  implicit none

  include 'mpif.h'

  integer :: tic, toc, countrate, ierr
  real :: end_time


  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
  proc0 = (mpi_rank==0)

  if (proc0) then
     print "(a)"," -------------------------------------------------------------"
     print *,"Quasisymmetry solver"
  end if
  !call system_clock(tic,countrate)
  call cpu_time(start_time)

  call quasisymmetry_read_input()
  call quasisymmetry_validate_input()

  select case (trim(general_option))
  case (general_option_single)
     call quasisymmetry_single_solve()
     call quasisymmetry_write_vmec_input()
  case (general_option_scan)
     call quasisymmetry_scan()
  case default
     print *,"Invalid general_option:",general_option
     stop
  end select

  !call system_clock(toc)
  !total_time = real(toc-tic)/countrate
  call cpu_time(end_time)
  total_time = end_time - start_time

  call quasisymmetry_write_output()

  if (proc0) then
     print "(a)"," -------------------------------------------------------------"
     print "(a,es10.3,a)","Quasisymmetry solver is complete. Total time=",total_time," sec."
  end if

  call mpi_finalize(ierr)

end program quasisymmetry

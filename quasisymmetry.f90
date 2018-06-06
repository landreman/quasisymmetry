! Main program

program quasisymmetry

  use quasisymmetry_variables, only: total_time, general_option, general_option_single, general_option_scan

  implicit none

  integer :: tic, toc, countrate
  real :: start_time, end_time

  print "(a)"," -------------------------------------------------------------"
  print *,"Quasisymmetry solver"
  !call system_clock(tic,countrate)
  call cpu_time(start_time)

  call quasisymmetry_read_input()
  call quasisymmetry_validate_input()

  select case (trim(general_option))
  case (general_option_single)
     call quasisymmetry_single_solve()
  case (general_option_scan)
     call quasisymmetry_scan()
     call quasisymmetry_write_output()
  case default
     print *,"Invalid general_option:",general_option
     stop
  end select

  !call system_clock(toc)
  !total_time = real(toc-tic)/countrate
  call cpu_time(end_time)
  total_time = end_time - start_time

  !call write_output()

  print "(a)"," -------------------------------------------------------------"
  print "(a,es10.3,a)","Quasisymmetry solver is complete. Total time=",total_time," sec."

end program quasisymmetry

! Main program

program quasisymmetry

  use quasisymmetry_variables, only: total_time, general_option

  implicit none

  integer :: tic, toc, countrate

  print *,"Quasisymmetry solver"
  call system_clock(tic,countrate)

  call quasisymmetry_read_input()
  call validate_input()

  select case (general_option)
  case (1)
     call quasisymmetry_single_solve()
  !case (2)
  !   call quasisymmetry_scan()
  case default
     print *,"Invalid general_option:",general_option
     stop
  end select

  call system_clock(toc)
  total_time = real(toc-tic)/countrate

  call write_output()
 
  print *,"quasisymmetry solver is complete. Total time=",total_time,"sec."

end program quasisymmetry

subroutine quasisymmetry_read_input

  use quasisymmetry_variables

  implicit none

  integer :: numargs
  character(len=200) :: input_filename
  integer :: fileUnit, didFileAccessWork, i
  integer, parameter :: uninitialized = -9999
  real(dp) :: threshold

  R0s = 0
  R0c = 0
  Z0s = 0
  Z0c = 0

  ! getcarg is in LIBSTELL
  call getcarg(1, input_filename, numargs)

  if (numargs<1) then
     stop "One argument is required: the input namelist file, which must be named quasisymmetry_in.XXXXX"
  end if
  if (numargs>1) then
     print *,"WARNING: Arguments after the first will be ignored."
  end if
  if (input_filename(1:17) .ne. "quasisymmetry_in.") then
     stop "Input file must be named quasisymmetry_in.XXX for some extension XXX"
  end if

  output_filename = "quasisymmetry_out" // trim(input_filename(17:)) // ".nc"
  new_vmec_filename = "input" // trim(input_filename(17:))

  fileUnit=11
  open(unit=fileUnit, file=input_filename, action="read", status="old", iostat=didFileAccessWork)
  if (didFileAccessWork /= 0) then
     print *,"Error opening input file ", trim(input_filename)
     stop
  else
     read(fileUnit, nml=quasisymmetry, iostat=didFileAccessWork)
     if (didFileAccessWork /= 0) then
        print *,"Error!  I was able to open the file ", trim(input_filename), &
               " but not read data from the quasisymmetry namelist in it."
        if (didFileAccessWork==-1) then
           print *,"Make sure there is a carriage return after the / at the end of the namelist!"
        end if
        stop
     end if
     if (proc0) print *,"Successfully read parameters from quasisymmetry namelist in ", trim(input_filename), "."
  end if
  close(unit = fileUnit)

  verbose = (trim(verbose_option)==verbose_option_all .or. (proc0 .and. trim(verbose_option)==verbose_option_proc0))

  if (proc0) print *,"Number of field periods (nfp):",nfp

!!$  ! Count the number of nonzero entries at the beginning of N_phis:
!!$  N_N_phis = 0
!!$  do i=1,max_N_N_phis
!!$     if (N_phis(i) == 0) then
!!$        N_N_phis = i-1
!!$        exit
!!$     end if
!!$  end do

  ! Find how many Fourier modes we are keeping in the axis shape:
  threshold = 1.0d-14
  if (trim(general_option)==general_option_single) then
     ! Set axis_nmax using the single-run options.
     do i = max_axis_nmax+1,1,-1
        !print *,'********* i=',i
        !if (R0c(i) .ne. 0 .or. R0s(i) .ne. 0 .or. Z0s(i) .ne. 0 .or. Z0c(i) .ne. 0) then
        if (abs(R0c(i)) > threshold .or. abs(R0s(i)) > threshold .or. abs(Z0s(i)) > threshold .or. abs(Z0c(i)) > threshold) then
           !print *,"R0c(i):",R0c(i)
           !print *,"R0s(i):",R0s(i)
           !print *,"Z0c(i):",Z0c(i)
           !print *,"Z0s(i):",Z0s(i)
           axis_nmax = i-1
           exit
        end if
     end do

     if (proc0) then
        print *,"R0c:", R0c(1:axis_nmax+1)
        print *,"R0s:", R0s(1:axis_nmax+1)
        print *,"Z0c:", Z0c(1:axis_nmax+1)
        print *,"Z0s:", Z0s(1:axis_nmax+1)
     end if
  else
     ! Set axis_nmax using the scan options
     do i = max_axis_nmax+1,1,-1
        if (abs(R0c_min(i)) > threshold .or. abs(R0s_min(i)) > threshold .or. abs(Z0s_min(i)) > threshold .or. abs(Z0c_min(i)) > threshold &
             .or. abs(R0c_max(i)) > threshold .or. abs(R0s_max(i)) > threshold .or. abs(Z0s_max(i)) > threshold .or. abs(Z0c_max(i)) > threshold) then
           axis_nmax = i-1
           exit
        end if
     end do
  end if

  if (proc0) print *,"axis_nmax:",axis_nmax

end subroutine quasisymmetry_read_input

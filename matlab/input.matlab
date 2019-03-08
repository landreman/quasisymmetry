! This &INDATA namelist was generated by quasisymmetry.f90
! Based on template file /Users/mattland/quasisymmetry/input.li383_vacuum
! r =    1.050000000000E-01
!----- Runtime Parameters -----
&INDATA
  DELT =    9.000000000000E-001
  NITER = 2500
  NSTEP = 200
  TCON0 =    2.000000000000E+000
  NS_ARRAY =                16            49
  FTOL_ARRAY =    1.000000E-06  1.000000E-11
  PRECON_TYPE = 'NONE'
  PREC2D_THRESHOLD =   1.000000E-30
!----- Grid Parameters -----
  LASYM = F
  NFP = 0003
  MPOL = 0009
  NTOR = 0015
  PHIEDGE =    3.463605900583E-002
!----- Free Boundary Parameters -----
  LFREEB = F
!----- Pressure Parameters -----
  GAMMA =    0.000000000000E+000
  BLOAT =    1.000000000000E+000
  SPRES_PED =    1.000000000000E+000
  PRES_SCALE =    1.000000000000E+000
  PMASS_TYPE = 'power_series'
  AM =  -0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00
!----- ANI/FLOW Parameters -----
  BCRIT =    1.000000000000E+000
  PT_TYPE = 'power_series'
  AT =   1.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00
  PH_TYPE = 'power_series'
  AH =   0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00
!----- Current/Iota Parameters -----
  CURTOR =    0.000000000000E+000
  NCURR = 1
  PIOTA_TYPE = 'power_series'
  AI =   0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00
  PCURR_TYPE = 'power_series'
  AC =    1.000000000000E+000   0.000000000000E+000   0.000000000000E+000   0.000000000000E+000
   0.000000000000E+000   0.000000000000E+000   0.000000000000E+000   0.000000000000E+000
   0.000000000000E+000   0.000000000000E+000   0.000000000000E+000   0.000000000000E+000
   0.000000000000E+000   0.000000000000E+000   0.000000000000E+000   0.000000000000E+000
   0.000000000000E+000   0.000000000000E+000   0.000000000000E+000   0.000000000000E+000
   0.000000000000E+000
!----- Axis Parameters ----- 
  RAXIS_CC =   1.00000000000000E+00  2.84313600000000E-01  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
  ZAXIS_CS =   -0.000000000000E+000  -3.000000000000E-001  -0.000000000000E+000  -0.000000000000E+000
  -0.000000000000E+000  -0.000000000000E+000  -0.000000000000E+000  -0.000000000000E+000
  -0.000000000000E+000  -0.000000000000E+000  -0.000000000000E+000  -0.000000000000E+000
  -0.000000000000E+000  -0.000000000000E+000  -0.000000000000E+000  -0.000000000000E+000
!----- Boundary Parameters -----
  RBC( 000,000) =    1.000000000000E+000    ZBS( 000,000) =   -0.000000000000E+000
  RBC( 001,000) =    2.843136000000E-001    ZBS( 001,000) =   -3.000000000000E-001
  RBC(-015,001) =    1.217957572736E-006    ZBS(-015,001) =    1.678150521751E-006
  RBC(-014,001) =   -3.651839132292E-006    ZBS(-014,001) =   -3.187091955802E-006
  RBC(-013,001) =    7.499262593206E-006    ZBS(-013,001) =    7.836192048548E-006
  RBC(-012,001) =   -1.615530031466E-005    ZBS(-012,001) =   -1.800626213596E-005
  RBC(-011,001) =    3.570787426753E-005    ZBS(-011,001) =    4.035484388849E-005
  RBC(-010,001) =   -7.946629261034E-005    ZBS(-010,001) =   -9.001718640851E-005
  RBC(-009,001) =    1.772313954843E-004    ZBS(-009,001) =    2.007710251596E-004
  RBC(-008,001) =   -3.955936473114E-004    ZBS(-008,001) =   -4.479786396947E-004
  RBC(-007,001) =    8.835472528931E-004    ZBS(-007,001) =    1.000256036749E-003
  RBC(-006,001) =   -1.983445482204E-003    ZBS(-006,001) =   -2.244016364867E-003
  RBC(-005,001) =    4.538274846928E-003    ZBS(-005,001) =    5.118866983197E-003
  RBC(-004,001) =   -1.065434021071E-002    ZBS(-004,001) =   -1.192383028403E-002
  RBC(-003,001) =    2.392519270621E-002    ZBS(-003,001) =    2.655464027549E-002
  RBC(-002,001) =   -4.071322629120E-002    ZBS(-002,001) =   -4.634415721375E-002
  RBC(-001,001) =   -4.803895688381E-002    ZBS(-001,001) =   -2.968344296476E-002
  RBC( 000,001) =    8.462362077270E-002    ZBS( 000,001) =    2.306728177445E-002
  RBC( 001,001) =   -1.556149847767E-001    ZBS( 001,001) =    1.306358425660E-001
  RBC( 002,001) =   -9.072914557461E-003    ZBS( 002,001) =    1.483303462462E-002
  RBC( 003,001) =    2.023128581140E-002    ZBS( 003,001) =   -2.264276060112E-002
  RBC( 004,001) =   -1.076881496867E-002    ZBS( 004,001) =    1.198068101934E-002
  RBC( 005,001) =    4.640822392636E-003    ZBS( 005,001) =   -5.214690676986E-003
  RBC( 006,001) =   -1.994874566023E-003    ZBS( 006,001) =    2.255469635889E-003
  RBC( 007,001) =    8.823057912935E-004    ZBS( 007,001) =   -9.991504346405E-004
  RBC( 008,001) =   -3.950203536812E-004    ZBS( 008,001) =    4.474218199416E-004
  RBC( 009,001) =    1.771584042845E-004    ZBS( 009,001) =   -2.006964697253E-004
  RBC( 010,001) =   -7.946749155228E-005    ZBS( 010,001) =    9.001757551736E-005
  RBC( 011,001) =    3.570972845132E-005    ZBS( 011,001) =   -4.035658309450E-005
  RBC( 012,001) =   -1.615554723883E-005    ZBS( 012,001) =    1.800650740644E-005
  RBC( 013,001) =    7.499252693321E-006    ZBS( 013,001) =   -7.836184318366E-006
  RBC( 014,001) =   -3.651830816698E-006    ZBS( 014,001) =    3.187083928395E-006
  RBC( 015,001) =    1.217956370439E-006    ZBS( 015,001) =   -1.678149327980E-006
/

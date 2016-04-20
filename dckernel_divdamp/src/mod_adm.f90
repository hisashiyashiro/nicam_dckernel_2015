module mod_adm
  implicit none

  public :: ADM_MPItime
  !
  !------ Character length of system control
  integer, public, parameter :: ADM_NSYS = 32
  !
  !------ Maximum length of file name
  integer, public, parameter :: ADM_MAXFNAME = 128
  !
  !====== Basic definition & information ======
  !
  !------ Log file ID & Control file ID
  integer, public, save      :: ADM_LOG_FID = 6 ! default is STDOUT
  integer, public, parameter :: ADM_CTL_FID = 35
  !
  !------ Identifier for single computation or parallel computation
  integer, public, parameter :: ADM_SINGLE_PRC = 0
  integer, public, parameter :: ADM_MULTI_PRC  = 1
  !
  !------ Identifiers of directions of region edges
  integer, public, parameter :: ADM_SW = 1
  integer, public, parameter :: ADM_NW = 2
  integer, public, parameter :: ADM_NE = 3
  integer, public, parameter :: ADM_SE = 4
  !
  !------ Identifiers of directions of region vertices
  integer, public, parameter :: ADM_W = 1
  integer, public, parameter :: ADM_N = 2
  integer, public, parameter :: ADM_E = 3
  integer, public, parameter :: ADM_S = 4
  !
  !--- Identifier of triangle element (i-axis-side or j-axis side)
  integer, public, parameter :: ADM_TI = 1
  integer, public, parameter :: ADM_TJ = 2
  !
  !--- Identifier of line element (i-axis-side, ij-axis side, or j-axis side)
  integer, public, parameter :: ADM_AI  = 1
  integer, public, parameter :: ADM_AIJ = 2
  integer, public, parameter :: ADM_AJ  = 3
  !
  !------ Identifier of 1 variable
  integer, public, parameter :: ADM_KNONE = 1
  integer, public, parameter :: ADM_VNONE = 1
  !
  !------ Identifier of poles (north pole or south pole)
  integer, public, parameter :: ADM_NPL = 1
  integer, public, parameter :: ADM_SPL = 2
  !
  !------ Fist colomn on the table for region and direction
  integer, public, parameter :: ADM_RID = 1
  integer, public, parameter :: ADM_DIR = 2
  !
  real(8), public, parameter :: ADM_VMISS = 1.D0

  !
  !====== Information for processes ======
  !
  !------ Communication world for NICAM
  integer, public, save      :: ADM_COMM_WORLD
  logical, public, save      :: ADM_MPI_alive = .false.
  !
  !------ Master process
  integer, public, parameter :: ADM_prc_run_master = 1
  !
  !------ Total number of process
  integer, public, save      :: ADM_prc_all
  !
  !------ My process ID
  integer, public, save      :: ADM_prc_me
  !
  !------ Process ID which manages the pole regions.
  integer, public, save      :: ADM_prc_pl
  !
  !------ Process ID which have the pole regions.
  integer, public, save      :: ADM_prc_npl
  integer, public, save      :: ADM_prc_spl
  integer, public, save      :: ADM_prc_nspl(ADM_NPL:ADM_SPL)
  logical, public, save      :: ADM_have_pl

  !
  !====== Information for processes-region relationship ======
  !

  !
  !====== Information for regions ======
  !
  !------ Region division level
  integer, public, save      :: ADM_rlevel
  !
  !------ Total number of regular regions managed by all process
  integer, public, save      :: ADM_rgn_nmax
  !
  !------ Maximum number of pole regions
  integer, public, parameter :: ADM_rgn_nmax_pl = 2
  !
  !------ Local region number
  integer, public, save      :: ADM_lall
  !
  !------ Local region number for poles
  integer, public, save      :: ADM_lall_pl = ADM_rgn_nmax_pl
  !
  !------ Present Local region number ! 2010.4.26 M.Satoh
  integer, public, save      :: ADM_l_me

  logical, public, allocatable, save :: ADM_have_sgp(:) ! region have singlar point?

  !
  !====== Grid resolution informations  ======
  !
  !------ Grid division level
  integer, public, save      :: ADM_glevel
  !
  !------ Horizontal grid numbers
  integer, public, save      :: ADM_gmin
  integer, public, save      :: ADM_gmax
  integer, public, save      :: ADM_gall_1d
  integer, public, save      :: ADM_gall
  !
  !----- grid number of inner region in the diamond
  integer, public, save      :: ADM_gall_in
  !
  !------ Identifiers of grid points around poles.
  integer, public, parameter :: ADM_gslf_pl = 1
  integer, public, parameter :: ADM_gmin_pl = 2
  integer, public, save      :: ADM_gmax_pl     ! [mod] S.Iga 100607
  integer, public, save      :: ADM_gall_pl     ! [mod] S.Iga 100607
  !
  !------ Vertica grid numbers
  integer, public, save      :: ADM_vlayer
  integer, public, save      :: ADM_kmin
  integer, public, save      :: ADM_kmax
  integer, public, save      :: ADM_kall

contains

  function ADM_MPItime() result(time)
    implicit none
    real(8) :: time
       !cx  time = real(MPI_WTIME(), kind=8)
       call cpu_time(time)
  end function ADM_MPItime

end module mod_adm

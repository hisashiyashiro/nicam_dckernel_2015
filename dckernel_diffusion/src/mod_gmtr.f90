module mod_gmtr
  implicit none
  public :: GMTR_setup

  integer, public, parameter :: GMTR_T_nmax_var = 7

  integer, public, parameter :: GMTR_T_AREA  = 1
  integer, public, parameter :: GMTR_T_RAREA = 2
  integer, public, parameter :: GMTR_T_W1    = 3
  integer, public, parameter :: GMTR_T_W2    = 4
  integer, public, parameter :: GMTR_T_W3    = 5
  integer, public, parameter :: GMTR_T_LAT   = 6
  integer, public, parameter :: GMTR_T_LON   = 7

  integer, public, parameter :: GMTR_P_nmax_var = 10

  integer, public, parameter :: GMTR_P_AREA  = 1
  integer, public, parameter :: GMTR_P_RAREA = 2
  integer, public, parameter :: GMTR_P_IX    = 3
  integer, public, parameter :: GMTR_P_IY    = 4
  integer, public, parameter :: GMTR_P_IZ    = 5
  integer, public, parameter :: GMTR_P_JX    = 6
  integer, public, parameter :: GMTR_P_JY    = 7
  integer, public, parameter :: GMTR_P_JZ    = 8
  integer, public, parameter :: GMTR_P_LAT   = 9
  integer, public, parameter :: GMTR_P_LON   = 10

  integer, public, parameter :: GMTR_A_nmax_var    = 12
  integer, public, parameter :: GMTR_A_nmax_var_pl = 18

  integer, public, parameter :: GMTR_A_HNX  = 1
  integer, public, parameter :: GMTR_A_HNY  = 2
  integer, public, parameter :: GMTR_A_HNZ  = 3
  integer, public, parameter :: GMTR_A_HTX  = 4
  integer, public, parameter :: GMTR_A_HTY  = 5
  integer, public, parameter :: GMTR_A_HTZ  = 6
  integer, public, parameter :: GMTR_A_TNX  = 7
  integer, public, parameter :: GMTR_A_TNY  = 8
  integer, public, parameter :: GMTR_A_TNZ  = 9
  integer, public, parameter :: GMTR_A_TTX  = 10
  integer, public, parameter :: GMTR_A_TTY  = 11
  integer, public, parameter :: GMTR_A_TTZ  = 12

  integer, public, parameter :: GMTR_A_TN2X = 13
  integer, public, parameter :: GMTR_A_TN2Y = 14
  integer, public, parameter :: GMTR_A_TN2Z = 15
  integer, public, parameter :: GMTR_A_TT2X = 16
  integer, public, parameter :: GMTR_A_TT2Y = 17
  integer, public, parameter :: GMTR_A_TT2Z = 18

  !cx   real(8), public, allocatable, save :: GMTR_P_var   (:,:,:,:)
  !cx   real(8), public, allocatable, save :: GMTR_T_var   (:,:,:,:,:)
  !cx   real(8), public, allocatable, save :: GMTR_A_var   (:,:,:,:,:)

  real(8), public, allocatable, save :: GMTR_P_var_pl(:,:,:,:)
  real(8), public, allocatable, save :: GMTR_T_var_pl(:,:,:,:)
  real(8), public, allocatable, save :: GMTR_A_var_pl(:,:,:,:)
contains

  subroutine GMTR_setup
    use mod_adm, only: &
       ADM_TI,           &
       ADM_TJ,           &
       ADM_AI,           &
       ADM_AJ,           &
       ADM_gmin,         &
       ADM_gmax,         &
       ADM_gall,         &
       ADM_gall_1d,      &
       ADM_gall_pl,      &
       ADM_lall,         &
       ADM_lall_pl,      &
       K0 => ADM_KNONE
    implicit none

    write(0,'(a,10i8)') "** To Do: <GMTR_setup> should first touch the arrays."

    ! --- setup triangle data
    allocate( GMTR_T_var_pl(ADM_gall_pl,K0,ADM_lall_pl,GMTR_T_nmax_var) )
    GMTR_T_var_pl(:,:,:,:)   = 0.D0
    !--- setup point data
    allocate( GMTR_P_var_pl(ADM_gall_pl,K0,ADM_lall_pl,GMTR_P_nmax_var) )
    GMTR_P_var_pl(:,:,:,:) = 0.D0
    !--- setup arc data
    allocate( GMTR_A_var_pl(ADM_gall_pl,K0,ADM_lall_pl,GMTR_A_nmax_var_pl) )
    GMTR_A_var_pl(:,:,:,:)   = 0.D0

    !cx allocate(GMTR_T_var(ADM_gall,K0,ADM_lall,ADM_TI:ADM_TJ,GMTR_T_nmax_var))
    !cx allocate(GMTR_P_var(ADM_gall,K0,ADM_lall,GMTR_P_nmax_var))
    !cx allocate(GMTR_A_var(ADM_gall,K0,ADM_lall,ADM_AI:ADM_AJ,GMTR_A_nmax_var))
    !cx GMTR_T_var   (:,:,:,:,:) = 0.D0
    !cx GMTR_P_var   (:,:,:,:) = 0.D0
    !cx GMTR_A_var   (:,:,:,:,:) = 0.D0

    return
  end subroutine GMTR_setup

end module mod_gmtr


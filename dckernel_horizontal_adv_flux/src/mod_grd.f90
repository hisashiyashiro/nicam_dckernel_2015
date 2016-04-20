module mod_grd
  use mod_precision
  implicit none
  public :: GRD_setup

  !------ Indentifiers for the directions in the Cartesian coordinate.
  integer, public, parameter :: GRD_XDIR=1
  integer, public, parameter :: GRD_YDIR=2
  integer, public, parameter :: GRD_ZDIR=3

  !------ Grid points ( CELL ARC )
  real(RP), public, allocatable, save :: GRD_xr   (:,:,:,:,:)
  real(RP), public, allocatable, save :: GRD_xr_pl(:,:,:,:)

contains

  subroutine GRD_setup
    use mod_adm, only: &
        ADM_gall,    &
        ADM_gall_pl, &
        ADM_lall,    &
        ADM_lall_pl, &
        ADM_kall
    use mod_adm, only: &
        AI  => ADM_AI,  &
        AJ  => ADM_AJ,  &
        K0  => ADM_KNONE
    implicit none

    allocate(GRD_xr(ADM_gall ,K0, ADM_lall ,AI:AJ, GRD_XDIR:GRD_ZDIR))
    allocate(GRD_xr_pl(ADM_gall_pl, K0, ADM_lall_pl, GRD_XDIR:GRD_ZDIR))
    return
  end subroutine GRD_setup

end module mod_grd

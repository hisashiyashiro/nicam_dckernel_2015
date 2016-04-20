module mod_grd
  use mod_precision
  implicit none
  public :: GRD_setup
  real(RP), public, allocatable, save :: GRD_rdgzh(:)
  real(RP), public, allocatable, save :: GRD_afac(:)
  real(RP), public, allocatable, save :: GRD_bfac(:)

contains

  subroutine GRD_setup
    use mod_adm, only: &
        ADM_kall
    implicit none
        allocate( GRD_rdgzh(ADM_kall) )
        allocate( GRD_afac(ADM_kall) )
        allocate( GRD_bfac(ADM_kall) )
    return
  end subroutine GRD_setup

end module mod_grd

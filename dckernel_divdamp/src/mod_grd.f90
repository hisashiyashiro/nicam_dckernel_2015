module mod_grd
  use mod_precision
  public :: GRD_setup

  real(RP), public, allocatable, save ::  GRD_rdgz (:)
contains
  subroutine GRD_setup
    use mod_adm, only :  ADM_kall
    implicit none
    allocate( GRD_rdgz (ADM_kall) )
    return
  end subroutine GRD_setup
end module mod_grd


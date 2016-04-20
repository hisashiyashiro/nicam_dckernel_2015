module mod_oprt
   use mod_precision

  implicit none

  !++ Public parameters & variables
  integer, public, save :: OPRT_nstart
  integer, public, save :: OPRT_nend

  ! < for diffusion operator >
  real(RP), public, allocatable, save :: cinterp_TN (:,:,:,:)
  real(RP), public, allocatable, save :: cinterp_HN (:,:,:,:)
  real(RP), public, allocatable, save :: cinterp_TRA(:,:,:)
  real(RP), public, allocatable, save :: cinterp_PRA(:,:)

end module mod_oprt
!-------------------------------------------------------------------------------


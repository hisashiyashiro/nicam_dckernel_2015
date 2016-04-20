module mod_vmtr
  use mod_precision
  implicit none
  private
  !++ Public procedure
  public :: VMTR_setup

  !++ Public parameters & variables
  ! index for VMTR_C2Wfact, W2Cfact
  integer, public, parameter :: I_a = 1
  integer, public, parameter :: I_b = 2

  integer, public, parameter :: I_c = 1
  integer, public, parameter :: I_d = 2

  integer, public, parameter :: I_a_GZXH = 1
  integer, public, parameter :: I_b_GZXH = 2
  integer, public, parameter :: I_a_GZYH = 3
  integer, public, parameter :: I_b_GZYH = 4
  integer, public, parameter :: I_a_GZZH = 5
  integer, public, parameter :: I_b_GZZH = 6

  !--- Gamma^2 at the full level
  real(RP), public, allocatable :: VMTR_GAM2   (:,:,:)
  real(RP), public, allocatable :: VMTR_GAM2_pl(:,:,:)

  !--- Gamma^2 at the half level
  real(RP), public, allocatable :: VMTR_GAM2H   (:,:,:)
  real(RP), public, allocatable :: VMTR_GAM2H_pl(:,:,:)

  !--- G^1/2 X Gamma^2 at the full level
  real(RP), public, allocatable :: VMTR_GSGAM2   (:,:,:)
  real(RP), public, allocatable :: VMTR_GSGAM2_pl(:,:,:)

  !--- G^1/2 X Gamma^2 at the half level
  real(RP), public, allocatable :: VMTR_GSGAM2H   (:,:,:)
  real(RP), public, allocatable :: VMTR_GSGAM2H_pl(:,:,:)

  !--- 1 / G^1/2 at the half level
  real(RP), public, allocatable :: VMTR_RGSQRTH   (:,:,:)
  real(RP), public, allocatable :: VMTR_RGSQRTH_pl(:,:,:)

  !--- 1 / Gamma at the integer level
  real(RP), public, allocatable :: VMTR_RGAM   (:,:,:)
  real(RP), public, allocatable :: VMTR_RGAM_pl(:,:,:)

  !--- 1 / Gamma at the half level
  real(RP), public, allocatable :: VMTR_RGAMH   (:,:,:)
  real(RP), public, allocatable :: VMTR_RGAMH_pl(:,:,:)

  !--- 1 / (G^1/2 X Gamma^2) at the full level
  real(RP), public, allocatable :: VMTR_RGSGAM2   (:,:,:)
  real(RP), public, allocatable :: VMTR_RGSGAM2_pl(:,:,:)

  !--- 1 / (G^1/2 X Gamma^2) at the half level
  real(RP), public, allocatable :: VMTR_RGSGAM2H   (:,:,:)
  real(RP), public, allocatable :: VMTR_RGSGAM2H_pl(:,:,:)

  !--- volume at the full level
  real(RP), public, allocatable :: VMTR_VOLUME   (:,:,:)
  real(RP), public, allocatable :: VMTR_VOLUME_pl(:,:,:)

  !--- geopotential at the full level
  real(RP), public, allocatable :: VMTR_PHI   (:,:,:)
  real(RP), public, allocatable :: VMTR_PHI_pl(:,:,:)

  !--- factor for half to full level
  real(RP), public, allocatable :: VMTR_W2Cfact   (:,:,:,:)
  real(RP), public, allocatable :: VMTR_W2Cfact_pl(:,:,:,:)

  !--- factor for full to half level
  real(RP), public, allocatable :: VMTR_C2Wfact   (:,:,:,:)
  real(RP), public, allocatable :: VMTR_C2Wfact_pl(:,:,:,:)

  !--- factor for full to half level with Gz
  real(RP), public, allocatable :: VMTR_C2WfactGz   (:,:,:,:)
  real(RP), public, allocatable :: VMTR_C2WfactGz_pl(:,:,:,:)


contains
  !-----------------------------------------------------------------------------
  !>
  !> Setup the vertical metrics
  !>
  subroutine VMTR_setup
    use mod_adm, only: &
       ADM_gall,      &
       ADM_gall_pl,   &
       ADM_lall,      &
       ADM_lall_pl,   &
       ADM_kall,      &
       ADM_kmin,      &
       ADM_kmax,      &
       ADM_KNONE,     &
       ADM_gmin,      &
       ADM_gall_1d

    implicit none
    integer, parameter :: var_max = 6

    real(RP) :: var      (ADM_gall,   ADM_kall,ADM_lall,   var_max)
    real(RP) :: var_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,var_max)

    !--- G^1/2
    real(RP) :: GSQRT    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GSQRT_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: GSQRTH   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GSQRTH_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !--- Gamma factor
    real(RP) :: GAM      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GAM_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: GAMH     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GAMH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !--- vector G^z at the full level
    real(RP) :: GZX      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GZX_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: GZY      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GZY_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: GZZ      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GZZ_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !--- vector G^z at the half level
    real(RP) :: GZXH     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GZXH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: GZYH     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GZYH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: GZZH     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: GZZH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)


    !--- initialization
    allocate( VMTR_GAM2        (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GAM2_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_GAM2H       (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GAM2H_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_GSGAM2      (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GSGAM2_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_GSGAM2H     (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_GSGAM2H_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_RGSQRTH     (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSQRTH_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGAM        (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGAM_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGAMH       (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGAMH_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGSGAM2     (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSGAM2_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGSGAM2H    (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSGAM2H_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_VOLUME      (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_VOLUME_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_PHI         (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_PHI_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_W2Cfact     (2,ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_W2Cfact_pl  (2,ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_C2Wfact     (2,ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_C2Wfact_pl  (2,ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_C2WfactGz   (6,ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( VMTR_C2WfactGz_pl(6,ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    return
  end subroutine VMTR_setup

end module mod_vmtr
!-------------------------------------------------------------------------------

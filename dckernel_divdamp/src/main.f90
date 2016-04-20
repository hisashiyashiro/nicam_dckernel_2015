program main
!    use mpi
    use mod_debug, only: &
       DEBUG_rapreport
    use mod_snapshot_for_kernel, only : &
        snapshot_seq_bin_open,    &
        snapshot_seq_bin_close,    &
        ADM_snap_read,    &
        ADM_ascii_write
    use mod_oprt3d, only : &
        OPRT3D_snap_write, &
        OPRT3D_snap_read, &
        OPRT3D_divdamp
    use mod_gmtr, only : GMTR_setup
    use mod_vmtr, only : VMTR_setup
    use mod_grd, only : GRD_setup

    integer :: my_id=1
    logical :: test_with_snapshot = .false.
!    integer :: ierr
    !---------------------------------------------------------------------------

!    call mpi_init(ierr)

    if (test_with_snapshot) then
    ! for testing with realistic variables from snapshot.OPRT3D_divdamp.* file
        call snapshot_seq_bin_open ("OPRT3D_divdamp", my_id)
    else
    ! for testing without snapshot file
        call snapshot_seq_bin_open ("", my_id)
    endif

    call ADM_snap_read
    call ADM_ascii_write
    call GRD_setup
    call GMTR_setup
    call VMTR_setup

    call dynamics_step

    call DEBUG_rapreport

!    call mpi_finalize(ierr)

    stop
end program main


  subroutine dynamics_step
    use mod_precision
    use mod_oprt3d, only : &
        OPRT3D_snap_write, &
        OPRT3D_snap_read, &
        OPRT3D_divdamp
    use mod_adm, only: &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    implicit none
    integer :: i

    real(RP) :: ddivdx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!cx Better To Do: first touch above arrays

    call OPRT3D_snap_read( rhogvx, rhogvy, rhogvz, rhogw )

  !cx call fapp_start( 'OPRT3D_divdamp', 1, 1 )
    do i=1,14
    call OPRT3D_divdamp( ddivdx, ddivdy, ddivdz, &
                            rhogvx, rhogvy, rhogvz, rhogw  )

    end do
    !cx call fapp_stop( 'OPRT3D_divdamp', 1, 1 )

    return
  end subroutine dynamics_step


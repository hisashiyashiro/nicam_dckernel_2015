program main
!    use mpi
    use mod_debug, only: &
       DEBUG_rapreport
    use mod_snapshot_for_kernel, only : &
        snapshot_seq_bin_open,    &
        snapshot_seq_bin_close,    &
        ADM_snap_read,    &
        ADM_ascii_write
    use mod_src_tracer, only : &
        V_Adv_limiter_snap_write, &
        V_Adv_limiter_snap_read, &
        vertical_limiter_thuburn
#ifdef DO_WE_NEED_EXTRA_MODULES
    use mod_gmtr, only : GMTR_setup
    use mod_vmtr, only : VMTR_setup
    use mod_grd, only : GRD_setup
#endif

    integer :: my_id=1
    logical :: test_with_snapshot = .true.
!    integer :: ierr
    !---------------------------------------------------------------------------

!    call mpi_init(ierr)

    if (test_with_snapshot) then
        call snapshot_seq_bin_open ("Vertical_Adv_limiter", my_id)
    else ! for testing without snapshot file
        call snapshot_seq_bin_open ("", my_id)
    endif

    call ADM_snap_read
    call ADM_ascii_write
#ifdef DO_WE_NEED_EXTRA_MODULES
    call GMTR_setup
    call VMTR_setup
    call GRD_setup
#endif

    call dynamics_step

    call DEBUG_rapreport

!    call mpi_finalize(ierr)

    stop
end program main


  subroutine dynamics_step
    use mod_precision
    use mod_src_tracer, only : &
        V_Adv_limiter_snap_write, &
        V_Adv_limiter_snap_read, &
        vertical_limiter_thuburn
    use mod_adm, only: &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    implicit none
    integer :: i

    real(RP) :: q_h   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: q_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: q     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: q_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: d     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: d_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ck    (ADM_gall,   ADM_kall,ADM_lall   ,2)
    real(RP) :: ck_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)

    call V_Adv_limiter_snap_read (q_h, q_h_pl, q, q_pl, d, d_pl, ck, ck_pl)
    !cx call fapp_start( '_Vertical_Adv_limiter', 1, 1 )
    do i=1,60
    call vertical_limiter_thuburn (q_h, q_h_pl, q, q_pl, d, d_pl, ck, ck_pl)
    end do
    !cx call fapp_stop( '_Vertical_Adv_limiter', 1, 1 )

    return
  end subroutine dynamics_step


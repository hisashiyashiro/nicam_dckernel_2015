program main
!    use mpi
    use mod_debug, only: &
       DEBUG_rapreport
    use mod_snapshot_for_kernel, only : &
        snapshot_seq_bin_open,    &
        snapshot_seq_bin_close,    &
        ADM_snap_read,    &
        ADM_ascii_write

    integer :: my_id=1
    logical :: test_with_snapshot = .true.
!    integer :: ierr
    !---------------------------------------------------------------------------

!    call mpi_init(ierr)

    if (test_with_snapshot) then
        call snapshot_seq_bin_open ("Horizontal_Adv_limiter", my_id)
    else ! for testing without snapshot file
        call snapshot_seq_bin_open ("", my_id)
    endif
    call ADM_snap_read

    call dynamics_step

    call DEBUG_rapreport

!    call mpi_finalize(ierr)

    stop
end program main


  subroutine dynamics_step
    use mod_precision
    use mod_src_tracer, only : &
        H_Adv_limiter_snap_write, &
        H_Adv_limiter_snap_read, &
        horizontal_limiter_thuburn
    use mod_adm, only: &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    implicit none
    integer :: i

    real(RP) q_a     (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) q_a_pl  (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) q       (  ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) q_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) d       (  ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) d_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) ch      (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) ch_pl   (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
!fj    real(RP) cmask   (6,ADM_gall   ,ADM_kall,ADM_lall   )
!fj>
    real(RP) cmask   (3,ADM_gall   ,ADM_kall,ADM_lall   )
!fj<
    real(RP) cmask_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl)

    call H_Adv_limiter_snap_read &
                (q_a, q_a_pl, q, q_pl, d, d_pl, ch, ch_pl, cmask, cmask_pl)


    !cx call fapp_start( '____Horizontal_Adv_limiter', 1, 1 )
    do i=1,30
    call horizontal_limiter_thuburn &
                (q_a, q_a_pl, q, q_pl, d, d_pl, ch, ch_pl, cmask, cmask_pl)
    end do
    !cx call fapp_stop( '____Horizontal_Adv_limiter', 1, 1 )

    return
  end subroutine dynamics_step


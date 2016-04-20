program main
!    use mpi
    use mod_debug, only: &
       DEBUG_rapreport
    use mod_snapshot_for_kernel, only : &
        snapshot_seq_bin_open,    &
        snapshot_seq_bin_close,    &
        ADM_snap_read,    &
        ADM_ascii_write
    use mod_gmtr, only : GMTR_setup
    use mod_cnst, only : CNST_setup
    use mod_grd, only : GRD_setup
    use mod_vmtr, only : VMTR_setup
    use mod_vi, only : vi_setup

    integer :: my_id=1
    logical :: test_with_snapshot = .false.
!    integer :: ierr
    !---------------------------------------------------------------------------

!    call mpi_init(ierr)

    if (test_with_snapshot) then
        call snapshot_seq_bin_open ("vi_rhow_solver", my_id)
    else ! for testing without snapshot file
        call snapshot_seq_bin_open ("", my_id)
    endif

    call ADM_snap_read
    call GMTR_setup
    call CNST_setup
    call GRD_setup
    call VMTR_setup
    call vi_setup

    call dynamics_step

    call DEBUG_rapreport

!    call mpi_finalize(ierr)

    stop
end program main


  subroutine dynamics_step
    use mod_precision
    use mod_vi, only : &
        vi_rhow_snap_read, &
        vi_rhow_solver
    use mod_adm, only: &
        ADM_gall,    &
        ADM_gall_pl, &
        ADM_lall,    &
        ADM_lall_pl, &
        ADM_kall

    implicit none
    integer :: i

    real(RP) rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) rhogw0   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) rhogw0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) preg0    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) preg0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) rhog0    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) rhog0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) Sr       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) Sr_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) Sw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) Sw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) Sp       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) Sp_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) dt

    call vi_rhow_snap_read ( &
       rhogw,  rhogw_pl,  &
       rhogw0, rhogw0_pl, &
       preg0,  preg0_pl,  &
       rhog0,  rhog0_pl,  &
       Sr,     Sr_pl,     &
       Sw,     Sw_pl,     &
       Sp,     Sp_pl,     &
       dt                 )

    !cx call fapp_start( '____vi_rhow_solver', 1, 1 )
    do i=1,11
    call vi_rhow_solver( &
       rhogw,  rhogw_pl,  &
       rhogw0, rhogw0_pl, &
       preg0,  preg0_pl,  &
       rhog0,  rhog0_pl,  &
       Sr,     Sr_pl,     &
       Sw,     Sw_pl,     &
       Sp,     Sp_pl,     &
       dt                 )
    end do
 !cx call fapp_stop( '____vi_rhow_solver', 1, 1 )

    return
  end subroutine dynamics_step


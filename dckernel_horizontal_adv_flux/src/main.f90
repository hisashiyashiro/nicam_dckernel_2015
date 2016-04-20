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
    use mod_grd, only : GRD_setup

    integer :: my_id=1
    logical :: test_with_snapshot = .false.
!    integer :: ierr
    !---------------------------------------------------------------------------

!    call mpi_init(ierr)

    if (test_with_snapshot) then
        call snapshot_seq_bin_open ("Horizontal_Adv_flux", my_id)
    else ! for testing without snapshot file
        call snapshot_seq_bin_open ("", my_id)
    endif

    call ADM_snap_read
    call ADM_ascii_write
    call GRD_setup
    call GMTR_setup

    call dynamics_step

    call DEBUG_rapreport

!    call mpi_finalize(ierr)

    stop
end program main


  subroutine dynamics_step
    use mod_precision
    use mod_src_tracer, only : &
        H_Adv_flux_snap_read, &
        horizontal_flux
    use mod_adm, only: &
        ADM_gall,    &
        ADM_gall_pl, &
        ADM_lall,    &
        ADM_lall_pl, &
        ADM_kall
    use mod_adm, only: AI=>ADM_AI, AJ=>ADM_AJ, K0=>ADM_KNONE
    use mod_grd, only: XDIR=>GRD_XDIR, YDIR=>GRD_YDIR, ZDIR=>GRD_ZDIR

    implicit none
    integer :: i

    real(RP) flx_h    (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) flx_h_pl (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) GRD_xc   (ADM_gall   ,ADM_kall,ADM_lall   ,AI:AJ,XDIR:ZDIR)
    real(RP) GRD_xc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,      XDIR:ZDIR)
    real(RP) rho      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) rhovx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) rhovx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) rhovy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) rhovy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) rhovz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) rhovz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) dt

    call H_Adv_flux_snap_read ( &
       rho,    rho_pl,    &
       rhovx,  rhovx_pl,  &
       rhovy,  rhovy_pl,  &
       rhovz,  rhovz_pl,  &
       dt                 )

    !cx call fapp_start( '____Horizontal_Adv_flux', 1, 1 )
    do i=1,1
    call horizontal_flux( &
       flx_h,  flx_h_pl,  &
       GRD_xc, GRD_xc_pl, &
       rho,    rho_pl,    &
       rhovx,  rhovx_pl,  &
       rhovy,  rhovy_pl,  &
       rhovz,  rhovz_pl,  &
       dt                 )
    end do
    !cx call fapp_stop( '____Horizontal_Adv_flux', 1, 1 )

    return
  end subroutine dynamics_step





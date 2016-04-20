program main
!    use mpi
    use mod_precision
    use mod_debug, only: &
       DEBUG_rapreport
    use mod_snapshot_for_kernel, only : &
       snapshot_seq_bin_open, &
       ADM_snap_read
    use mod_adm, only: &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    use mod_gmtr, only: &
       GMTR_setup
    use mod_oprt, only : &
      OPRT_snap_read,   &
      OPRT_diffusion
    implicit none

    real(RP), allocatable :: dscl   (:,:,:)
    real(RP), allocatable :: dscl_pl(:,:,:)
    real(RP), allocatable :: scl    (:,:,:)
    real(RP), allocatable :: scl_pl (:,:,:)
    real(RP), allocatable :: kh     (:,:,:)
    real(RP), allocatable :: kh_pl  (:,:,:)
    real(RP) :: mfact

    integer :: my_id=1
    logical :: test_with_snapshot = .false.
!    integer :: ierr

    integer :: i
    !---------------------------------------------------------------------------

!    call mpi_init(ierr)

    if (test_with_snapshot) then
    ! for testing with realistic variables from snapshot.OPRT_diffusion.* file
        call snapshot_seq_bin_open ("OPRT_diffusion", my_id)
    else
    ! for testing without snapshot file
        call snapshot_seq_bin_open ("", my_id)
    endif

    write(0,'(a)') "Note: the arrays should be first touched in"
    write(0,'(a)') " <GMTR_setup> and <OPRT_snap_read>"

    call ADM_snap_read
    call GMTR_setup

    allocate( dscl   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( scl    (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( kh     (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    call OPRT_snap_read( scl, scl_pl, kh, kh_pl, mfact )

    do i = 1, 6
       call OPRT_diffusion( dscl, dscl_pl, scl, scl_pl, kh, kh_pl, mfact=1.0_RP)
    enddo

    call DEBUG_rapreport

!    call mpi_finalize(ierr)

    stop
end program main

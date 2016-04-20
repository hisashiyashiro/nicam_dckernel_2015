!-------------------------------------------------------------------------------
!>
!! Tracer advection module
!!
!! @par Description
!!         This module contains subroutines for tracer advection
!!
!! @li  chosen for performance evaluation targetting post-K
!!
!<
!-------------------------------------------------------------------------------

module mod_src_tracer
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID,    &
     TI  => ADM_TI,  &
     TJ  => ADM_TJ,  &
     AI  => ADM_AI,  &
     AIJ => ADM_AIJ, &
     AJ  => ADM_AJ,  &
     K0  => ADM_KNONE

    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_kall,     &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl,  &
       ADM_kmin,     &
       ADM_kmax
    use mod_snapshot_for_kernel


    use mod_adm, only: &
       ADM_have_pl,    &
       ADM_have_sgp,   &
       ADM_lall,       &
       ADM_lall_pl,    &
       ADM_gall,       &
       ADM_gall_pl,    &
       ADM_kall,       &
       ADM_gall_1d,    &
       ADM_gmin,       &
       ADM_gmax,       &
       ADM_gslf_pl,    &
       ADM_gmin_pl,    &
       ADM_gmax_pl
    use mod_cnst, only: &
       CNST_MAX_REAL, &
       CNST_EPS_ZERO
    use mod_snapshot_for_kernel
    use mod_adm, only: ADM_prc_me
  implicit none
  private

  public :: horizontal_limiter_thuburn
  public :: H_Adv_limiter_snap_write
  public :: H_Adv_limiter_snap_read


contains


  !-----------------------------------------------------------------------------
  !> Miura(2004)'s scheme with Thuburn(1996) limiter
  subroutine horizontal_limiter_thuburn( &
       q_a,    q_a_pl,  &
!fj       q,      q_pl,    &
       q_tmp,      q_pl,    &
       d,      d_pl,    &
       ch,     ch_pl,   &
       cmask,  cmask_pl )

!cx    use mod_comm, only: &
!cx       COMM_data_transfer

    implicit none

    real(RP), intent(inout) :: q_a     (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: q_a_pl  (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
!fj    real(RP), intent(in)    :: q       (  ADM_gall   ,ADM_kall,ADM_lall   )
!fj>
    real(RP), intent(in)    :: q_tmp   (  ADM_gall   ,ADM_kall,ADM_lall   )
!fj<
    real(RP), intent(in)    :: q_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: d       (  ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: d_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: ch      (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: ch_pl   (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
!fj    real(RP), intent(in)    :: cmask   (6,ADM_gall   ,ADM_kall,ADM_lall   )
!fj>
    real(RP), intent(in)    :: cmask   (3,ADM_gall   ,ADM_kall,ADM_lall   )
!fj<
    real(RP), intent(in)    :: cmask_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: q_min_AI, q_min_AIJ, q_min_AJ, q_min_pl
    real(RP) :: q_max_AI, q_max_AIJ, q_max_AJ, q_max_pl

!fj    real(8) :: qnext_min   (ADM_gall), qnext_min_pl
!fj    real(8) :: qnext_max   (ADM_gall), qnext_max_pl
!fj    real(8) :: Cin_sum     (ADM_gall), Cin_sum_pl
!fj    real(8) :: Cout_sum    (ADM_gall), Cout_sum_pl
!fj    real(8) :: CQin_max_sum(ADM_gall), CQin_max_sum_pl
!fj    real(8) :: CQin_min_sum(ADM_gall), CQin_min_sum_pl
    !fj>
    real(RP) :: qnext_min   , qnext_min_pl
    real(RP) :: qnext_max   , qnext_max_pl
    real(RP) :: Cin_sum     , Cin_sum_pl
    real(RP) :: Cout_sum    , Cout_sum_pl
    real(RP) :: CQin_max_sum, CQin_max_sum_pl
    real(RP) :: CQin_min_sum, CQin_min_sum_pl
    real(RP) :: q       (  -ADM_gall_1d+1:ADM_gall   ,ADM_kall,ADM_lall   )
    !fj<

    integer, parameter :: I_min = 1
    integer, parameter :: I_max = 2
!fj    real(8) :: Qin    (6,ADM_gall   ,ADM_kall,ADM_lall   ,2)
!fj>
    real(RP) :: Qin1    (ADM_gall      ,2)
    real(RP) :: Qin2    (ADM_gall      ,2)
    real(RP) :: Qin3    (ADM_gall      ,2)
    real(RP) :: Qin4    (ADM_gall      ,2)
    real(RP) :: Qin5    (ADM_gall      ,2)
    real(RP) :: Qin6    (ADM_gall      ,2)
!fj<
    real(RP) :: Qin_pl (2,ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    real(RP) :: Qout   (  ADM_gall   ,ADM_kall,ADM_lall   ,2)
    real(RP) :: Qout_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl,2)

    real(RP) :: zerosw
    !fj>
    real(RP) :: tmp,cmaskch1,cmaskch2,cmaskch3,cmaskch4,cmaskch5,cmaskch6
    !fj<

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1, ip2jp1
    integer :: im1j, ijm1

    integer :: nstart, nend
    integer :: n, k, l, v

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____Horizontal_Adv_limiter')

!cx    if (mark_for_snapshots(ip_Horizontal_Adv_limiter).ne.0) then
!cx    call snapshot_seq_bin_open ("Horizontal_Adv_limiter", ADM_prc_me)
!cx    call H_Adv_limiter_snap_write &
!cx                (q_a, q_a_pl, q, q_pl, d, d_pl, ch, ch_pl, cmask, cmask_pl)
!cx    call snapshot_seq_bin_close
!cx    mark_for_snapshots(ip_Horizontal_Adv_limiter) = 0
!cx    endif


    !---< (i) define inflow bounds, eq.(32)&(33) >---
    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = -ADM_gall_1d+1, 0
    q(n,k,l) = 0.0_RP
    enddo
    do n = 1, ADM_gall
    q(n,k,l) = q_tmp(n,k,l)
    enddo
    enddo
    enddo

    do l = 1, ADM_lall
    do k = 1, ADM_kall
!fj       nstart = suf(ADM_gmin,ADM_gmin)
       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax,ADM_gmax)
       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d

          q_min_AI  = min( q(ij,k,l), q(ijm1,k,l), q(ip1j,k,l), q(ip1jp1,k,l) )
          q_max_AI  = max( q(ij,k,l), q(ijm1,k,l), q(ip1j,k,l), q(ip1jp1,k,l) )
          q_min_AIJ = min( q(ij,k,l), q(ip1j,k,l), q(ip1jp1,k,l), q(ijp1,k,l) )
          q_max_AIJ = max( q(ij,k,l), q(ip1j,k,l), q(ip1jp1,k,l), q(ijp1,k,l) )
          q_min_AJ  = min( q(ij,k,l), q(ip1jp1,k,l), q(ijp1,k,l), q(im1j,k,l) )
          q_max_AJ  = max( q(ij,k,l), q(ip1jp1,k,l), q(ijp1,k,l), q(im1j,k,l) )

!fj          Qin(1,ij,    k,l,I_min) = (      cmask(1,n,k,l) ) * q_min_AI         &
!fj                                  + ( 1.D0-cmask(1,n,k,l) ) * CNST_MAX_REAL
!fj          Qin(4,ip1j,  k,l,I_min) = (      cmask(1,n,k,l) ) * CNST_MAX_REAL    &
!fj                                  + ( 1.D0-cmask(1,n,k,l) ) * q_min_AI
!fj          Qin(1,ij,    k,l,I_max) = (      cmask(1,n,k,l) ) * q_max_AI         &
!fj                                  + ( 1.0_RP-cmask(1,n,k,l) ) * (-CNST_MAX_REAL)
!fj          Qin(4,ip1j,  k,l,I_max) = (      cmask(1,n,k,l) ) * (-CNST_MAX_REAL) &
!fj                                  + ( 1.0_RP-cmask(1,n,k,l) ) * q_max_AI
!fj
!fj          Qin(2,ij,    k,l,I_min) = (      cmask(2,n,k,l) ) * q_min_AIJ        &
!fj                                  + ( 1.0_RP-cmask(2,n,k,l) ) * CNST_MAX_REAL
!fj          Qin(5,ip1jp1,k,l,I_min) = (      cmask(2,n,k,l) ) * CNST_MAX_REAL    &
!fj                                  + ( 1.0_RP-cmask(2,n,k,l) ) * q_min_AIJ
!fj          Qin(2,ij,    k,l,I_max) = (      cmask(2,n,k,l) ) * q_max_AIJ        &
!fj                                  + ( 1.0_RP-cmask(2,n,k,l) ) * (-CNST_MAX_REAL)
!fj          Qin(5,ip1jp1,k,l,I_max) = (      cmask(2,n,k,l) ) * (-CNST_MAX_REAL) &
!fj                                  + ( 1.0_RP-cmask(2,n,k,l) ) * q_max_AIJ
!fj
!fj          Qin(3,ij,    k,l,I_min) = (      cmask(3,n,k,l) ) * q_min_AJ         &
!fj                                  + ( 1.0_RP-cmask(3,n,k,l) ) * CNST_MAX_REAL
!fj          Qin(6,ijp1,  k,l,I_min) = (      cmask(3,n,k,l) ) * CNST_MAX_REAL    &
!fj                                  + ( 1.0_RP-cmask(3,n,k,l) ) * q_min_AJ
!fj          Qin(3,ij,    k,l,I_max) = (      cmask(3,n,k,l) ) * q_max_AJ         &
!fj                                  + ( 1.0_RP-cmask(3,n,k,l) ) * (-CNST_MAX_REAL)
!fj          Qin(6,ijp1,  k,l,I_max) = (      cmask(3,n,k,l) ) * (-CNST_MAX_REAL) &
!fj                                  + ( 1.0_RP-cmask(3,n,k,l) ) * q_max_AJ
!fj>
          Qin1(ij,    I_min) = (      cmask(1,n,k,l) ) * q_min_AI         &
                             + ( 1.0_RP-cmask(1,n,k,l) ) * CNST_MAX_REAL
          Qin4(ip1j,    I_min) = (      cmask(1,n,k,l) ) * CNST_MAX_REAL    &
                             + ( 1.0_RP-cmask(1,n,k,l) ) * q_min_AI
          Qin1(ij,    I_max) = (      cmask(1,n,k,l) ) * q_max_AI         &
                             + ( 1.0_RP-cmask(1,n,k,l) ) * (-CNST_MAX_REAL)
          Qin4(ip1j,    I_max) = (      cmask(1,n,k,l) ) * (-CNST_MAX_REAL) &
                             + ( 1.0_RP-cmask(1,n,k,l) ) * q_max_AI

          Qin2(ij,    I_min) = (      cmask(2,n,k,l) ) * q_min_AIJ        &
                             + ( 1.0_RP-cmask(2,n,k,l) ) * CNST_MAX_REAL
          Qin5(ip1jp1,    I_min) = (      cmask(2,n,k,l) ) * CNST_MAX_REAL    &
                             + ( 1.0_RP-cmask(2,n,k,l) ) * q_min_AIJ
          Qin2(ij,    I_max) = (      cmask(2,n,k,l) ) * q_max_AIJ        &
                             + ( 1.0_RP-cmask(2,n,k,l) ) * (-CNST_MAX_REAL)
          Qin5(ip1jp1,    I_max) = (      cmask(2,n,k,l) ) * (-CNST_MAX_REAL) &
                             + ( 1.0_RP-cmask(2,n,k,l) ) * q_max_AIJ

          Qin3(ij,    I_min) = (      cmask(3,n,k,l) ) * q_min_AJ         &
                             + ( 1.0_RP-cmask(3,n,k,l) ) * CNST_MAX_REAL
          Qin6(ijp1,    I_min) = (      cmask(3,n,k,l) ) * CNST_MAX_REAL    &
                             + ( 1.0_RP-cmask(3,n,k,l) ) * q_min_AJ
          Qin3(ij,    I_max) = (      cmask(3,n,k,l) ) * q_max_AJ         &
                             + ( 1.0_RP-cmask(3,n,k,l) ) * (-CNST_MAX_REAL)
          Qin6(ijp1,    I_max) = (      cmask(3,n,k,l) ) * (-CNST_MAX_REAL) &
                             + ( 1.0_RP-cmask(3,n,k,l) ) * q_max_AJ
!fj<
       enddo

       ! peeling
!fj       nstart = suf(ADM_gmin-1,ADM_gmin  )
!fj       nend   = suf(ADM_gmin-1,ADM_gmin  )
!fj
!fj       do n = nstart, nend
!fj          ij     = n
!fj          ip1j   = n + 1
!fj          ip1jp1 = n + 1 + ADM_gall_1d
!fj          ijm1   = n     - ADM_gall_1d
!fj
!fj          q_min_AI  = min( q(ij,k,l), q(ip1j,k,l), q(ip1jp1,k,l), q(ijm1,k,l) )
!fj          q_max_AI  = max( q(ij,k,l), q(ip1j,k,l), q(ip1jp1,k,l), q(ijm1,k,l) )

!fj          Qin(1,ij,    k,l,I_min) = (      cmask(1,n,k,l) ) * q_min_AI         &
!fj                                  + ( 1.0_RP-cmask(1,n,k,l) ) * CNST_MAX_REAL
!fj          Qin(4,ip1j,  k,l,I_min) = (      cmask(1,n,k,l) ) * CNST_MAX_REAL    &
!fj                                  + ( 1.0_RP-cmask(1,n,k,l) ) * q_min_AI
!fj          Qin(1,ij,    k,l,I_max) = (      cmask(1,n,k,l) ) * q_max_AI         &
!fj                                  + ( 1.0_RP-cmask(1,n,k,l) ) * (-CNST_MAX_REAL)
!fj          Qin(4,ip1j,  k,l,I_max) = (      cmask(1,n,k,l) ) * (-CNST_MAX_REAL) &
!fj                                  + ( 1.0_RP-cmask(1,n,k,l) ) * q_max_AI
!fj>
!fj          Qin1(ij,    I_min) = (      cmask(1,n,k,l) ) * q_min_AI         &
!fj                             + ( 1.0_RP-cmask(1,n,k,l) ) * CNST_MAX_REAL
!fj          Qin4(ip1j,    I_min) = (      cmask(1,n,k,l) ) * CNST_MAX_REAL    &
!fj                             + ( 1.0_RP-cmask(1,n,k,l) ) * q_min_AI
!fj          Qin1(ij,    I_max) = (      cmask(1,n,k,l) ) * q_max_AI         &
!fj                             + ( 1.0_RP-cmask(1,n,k,l) ) * (-CNST_MAX_REAL)
!fj          Qin4(ip1j,    I_max) = (      cmask(1,n,k,l) ) * (-CNST_MAX_REAL) &
!fj                             + ( 1.0_RP-cmask(1,n,k,l) ) * q_max_AI
!fj<
!fj       enddo
!fj
!fj       ! peeling
!fj       nstart = suf(ADM_gmin-1,ADM_gmin-1)
!fj       nend   = suf(ADM_gmin-1,ADM_gmin  )
!fj
!fj       do n = nstart, nend
!fj          ij     = n
!fj          ip1j   = n + 1
!fj          ijp1   = n     + ADM_gall_1d
!fj          ip1jp1 = n + 1 + ADM_gall_1d
!fj
!fj          q_min_AIJ = min( q(ij,k,l), q(ip1jp1,k,l), q(ip1j,k,l), q(ijp1,k,l) )
!fj          q_max_AIJ = max( q(ij,k,l), q(ip1jp1,k,l), q(ip1j,k,l), q(ijp1,k,l) )
!fj          Qin(2,ij,    k,l,I_min) = (      cmask(2,n,k,l) ) * q_min_AIJ        &
!fj                                  + ( 1.0_RP-cmask(2,n,k,l) ) * CNST_MAX_REAL
!fj          Qin(5,ip1jp1,k,l,I_min) = (      cmask(2,n,k,l) ) * CNST_MAX_REAL    &
!fj                                  + ( 1.0_RP-cmask(2,n,k,l) ) * q_min_AIJ
!fj          Qin(2,ij,    k,l,I_max) = (      cmask(2,n,k,l) ) * q_max_AIJ        &
!fj                                  + ( 1.0_RP-cmask(2,n,k,l) ) * (-CNST_MAX_REAL)
!fj          Qin(5,ip1jp1,k,l,I_max) = (      cmask(2,n,k,l) ) * (-CNST_MAX_REAL) &
!fj                                  + ( 1.0_RP-cmask(2,n,k,l) ) * q_max_AIJ
!fj          Qin2(ij,    I_min) = (      cmask(2,n,k,l) ) * q_min_AIJ        &
!fj                             + ( 1.0_RP-cmask(2,n,k,l) ) * CNST_MAX_REAL
!fj          Qin5(ip1jp1,    I_min) = (      cmask(2,n,k,l) ) * CNST_MAX_REAL    &
!fj                             + ( 1.0_RP-cmask(2,n,k,l) ) * q_min_AIJ
!fj          Qin2(ij,    I_max) = (      cmask(2,n,k,l) ) * q_max_AIJ        &
!fj                             + ( 1.0_RP-cmask(2,n,k,l) ) * (-CNST_MAX_REAL)
!fj          Qin5(ip1jp1,    I_max) = (      cmask(2,n,k,l) ) * (-CNST_MAX_REAL) &
!fj                                  + ( 1.0_RP-cmask(2,n,k,l) ) * q_max_AIJ
!fj       enddo

       ! peeling
!fj       nstart = suf(ADM_gmin,  ADM_gmin-1)
!fj       nend   = suf(ADM_gmin-1,ADM_gmin  )
!fj
!fj       do n = nstart, nend
!fj          ij     = n
!fj          ijp1   = n     + ADM_gall_1d
!fj          ip1jp1 = n + 1 + ADM_gall_1d
!fj          im1j   = n - 1
!fj
!fj          q_min_AJ  = min( q(ij,k,l), q(ijp1,k,l), q(ip1jp1,k,l), q(im1j,k,l) )
!fj          q_max_AJ  = max( q(ij,k,l), q(ijp1,k,l), q(ip1jp1,k,l), q(im1j,k,l) )

!fj          Qin(3,ij,    k,l,I_min) = (      cmask(3,n,k,l) ) * q_min_AJ         &
!fj                                  + ( 1.0_RP-cmask(3,n,k,l) ) * CNST_MAX_REAL
!fj          Qin(6,ijp1,  k,l,I_min) = (      cmask(3,n,k,l) ) * CNST_MAX_REAL    &
!fj                                  + ( 1.0_RP-cmask(3,n,k,l) ) * q_min_AJ
!fj          Qin(3,ij,    k,l,I_max) = (      cmask(3,n,k,l) ) * q_max_AJ         &
!fj                                  + ( 1.0_RP-cmask(3,n,k,l) ) * (-CNST_MAX_REAL)
!fj          Qin(6,ijp1,  k,l,I_max) = (      cmask(3,n,k,l) ) * (-CNST_MAX_REAL) &
!fj                                  + ( 1.0_RP-cmask(3,n,k,l) ) * q_max_AJ
!fj>
!fj          Qin3(ij,    I_min) = (      cmask(3,n,k,l) ) * q_min_AJ         &
!fj                             + ( 1.D0-cmask(3,n,k,l) ) * CNST_MAX_REAL
!fj          Qin6(ijp1,    I_min) = (      cmask(3,n,k,l) ) * CNST_MAX_REAL    &
!fj                             + ( 1.D0-cmask(3,n,k,l) ) * q_min_AJ
!fj          Qin3(ij,    I_max) = (      cmask(3,n,k,l) ) * q_max_AJ         &
!fj                             + ( 1.D0-cmask(3,n,k,l) ) * (-CNST_MAX_REAL)
!fj          Qin6(ijp1,    I_max) = (      cmask(3,n,k,l) ) * (-CNST_MAX_REAL) &
!fj                                  + ( 1.D0-cmask(3,n,k,l) ) * q_max_AJ
!fj       enddo

       if ( ADM_have_sgp(l) ) then
          n = suf(ADM_gmin-1,ADM_gmin-1)
          ij     = n
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          ip2jp1 = n + 2 + ADM_gall_1d

          q_min_AIJ = min( q(ij,k,l), q(ip1jp1,k,l), q(ip2jp1,k,l), q(ijp1,k,l) )
          q_max_AIJ = max( q(ij,k,l), q(ip1jp1,k,l), q(ip2jp1,k,l), q(ijp1,k,l) )

!fj          Qin(2,ij,    k,l,I_min) = (      cmask(2,n,k,l) ) * q_min_AIJ        &
!fj                                  + ( 1.D0-cmask(2,n,k,l) ) * CNST_MAX_REAL
!fj          Qin(5,ip1jp1,k,l,I_min) = (      cmask(2,n,k,l) ) * CNST_MAX_REAL    &
!fj                                  + ( 1.D0-cmask(2,n,k,l) ) * q_min_AIJ
!fj          Qin(2,ij,    k,l,I_max) = (      cmask(2,n,k,l) ) * q_max_AIJ        &
!fj                                  + ( 1.D0-cmask(2,n,k,l) ) * (-CNST_MAX_REAL)
!fj          Qin(5,ip1jp1,k,l,I_max) = (      cmask(2,n,k,l) ) * (-CNST_MAX_REAL) &
!fj                                  + ( 1.D0-cmask(2,n,k,l) ) * q_max_AIJ
!fj>
          Qin2(ij,    I_min) = (      cmask(2,n,k,l) ) * q_min_AIJ        &
                             + ( 1.0_RP-cmask(2,n,k,l) ) * CNST_MAX_REAL
          Qin5(ip1jp1,    I_min) = (      cmask(2,n,k,l) ) * CNST_MAX_REAL    &
                             + ( 1.0_RP-cmask(2,n,k,l) ) * q_min_AIJ
          Qin2(ij,    I_max) = (      cmask(2,n,k,l) ) * q_max_AIJ        &
                             + ( 1.0_RP-cmask(2,n,k,l) ) * (-CNST_MAX_REAL)
          Qin5(ip1jp1,    I_max) = (      cmask(2,n,k,l) ) * (-CNST_MAX_REAL) &
                             + ( 1.0_RP-cmask(2,n,k,l) ) * q_max_AIJ
!fj<
       endif

!fj    enddo
!fj    enddo


    !---< (iii) define allowable range of q at next step, eq.(42)&(43) >---
    nstart = suf(ADM_gmin, ADM_gmin )
    nend   = suf(ADM_gmax ,ADM_gmax )

!fj    do l = 1, ADM_lall
!fj    do k = 1, ADM_kall


!ocl simd
       do n = nstart, nend
       !fj>
       cmaskch1=min(0.0_RP,ch(1,n,k,l))
       cmaskch2=min(0.0_RP,ch(2,n,k,l))
       cmaskch3=min(0.0_RP,ch(3,n,k,l))
       cmaskch4=min(0.0_RP,ch(4,n,k,l))
       cmaskch5=min(0.0_RP,ch(5,n,k,l))
       cmaskch6=min(0.0_RP,ch(6,n,k,l))
       !fj<
!fj          Cin_sum(n)  = (      cmask(1,n,k,l) ) * ch(1,n,k,l) &
!fj                      + (      cmask(2,n,k,l) ) * ch(2,n,k,l) &
!fj                      + (      cmask(3,n,k,l) ) * ch(3,n,k,l) &
!fj                      + (      cmask(4,n,k,l) ) * ch(4,n,k,l) &
!fj                      + (      cmask(5,n,k,l) ) * ch(5,n,k,l) &
!fj                      + (      cmask(6,n,k,l) ) * ch(6,n,k,l)
!fj>
          Cin_sum     =       cmaskch1 &
                      +       cmaskch2 &
                      +       cmaskch3 &
                      +       cmaskch4 &
                      +       cmaskch5 &
                      +       cmaskch6
!fj<

!fj          Cout_sum(n) = ( 1.D0-cmask(1,n,k,l) ) * ch(1,n,k,l) &
!fj                      + ( 1.D0-cmask(2,n,k,l) ) * ch(2,n,k,l) &
!fj                      + ( 1.D0-cmask(3,n,k,l) ) * ch(3,n,k,l) &
!fj                      + ( 1.D0-cmask(4,n,k,l) ) * ch(4,n,k,l) &
!fj                      + ( 1.D0-cmask(5,n,k,l) ) * ch(5,n,k,l) &
!fj                      + ( 1.D0-cmask(6,n,k,l) ) * ch(6,n,k,l)
          Cout_sum    =  ch(1,n,k,l)-cmaskch1 &
                      +  ch(2,n,k,l)-cmaskch2 &
                      +  ch(3,n,k,l)-cmaskch3 &
                      +  ch(4,n,k,l)-cmaskch4 &
                      +  ch(5,n,k,l)-cmaskch5 &
                      +  ch(6,n,k,l)-cmaskch6
!fj       enddo

!fj       do n = nstart, nend
!fj          CQin_min_sum(n) = cmask(1,n,k,l) * ch(1,n,k,l) * Qin(1,n,k,l,I_min) &
!fj                          + cmask(2,n,k,l) * ch(2,n,k,l) * Qin(2,n,k,l,I_min) &
!fj                          + cmask(3,n,k,l) * ch(3,n,k,l) * Qin(3,n,k,l,I_min) &
!fj                          + cmask(4,n,k,l) * ch(4,n,k,l) * Qin(4,n,k,l,I_min) &
!fj                          + cmask(5,n,k,l) * ch(5,n,k,l) * Qin(5,n,k,l,I_min) &
!fj                          + cmask(6,n,k,l) * ch(6,n,k,l) * Qin(6,n,k,l,I_min)
!fj
!fj          CQin_max_sum(n) = cmask(1,n,k,l) * ch(1,n,k,l) * Qin(1,n,k,l,I_max) &
!fj                          + cmask(2,n,k,l) * ch(2,n,k,l) * Qin(2,n,k,l,I_max) &
!fj                          + cmask(3,n,k,l) * ch(3,n,k,l) * Qin(3,n,k,l,I_max) &
!fj                          + cmask(4,n,k,l) * ch(4,n,k,l) * Qin(4,n,k,l,I_max) &
!fj                          + cmask(5,n,k,l) * ch(5,n,k,l) * Qin(5,n,k,l,I_max) &
!fj                          + cmask(6,n,k,l) * ch(6,n,k,l) * Qin(6,n,k,l,I_max)
!fj>
          CQin_min_sum    = cmaskch1 * Qin1(n,I_min) &
                          + cmaskch2 * Qin2(n,I_min) &
                          + cmaskch3 * Qin3(n,I_min) &
                          + cmaskch4 * Qin4(n,I_min) &
                          + cmaskch5 * Qin5(n,I_min) &
                          + cmaskch6 * Qin6(n,I_min)

          CQin_max_sum    = cmaskch1 * Qin1(n,I_max) &
                          + cmaskch2 * Qin2(n,I_max) &
                          + cmaskch3 * Qin3(n,I_max) &
                          + cmaskch4 * Qin4(n,I_max) &
                          + cmaskch5 * Qin5(n,I_max) &
                          + cmaskch6 * Qin6(n,I_max)
          qnext_min    = min(Qin1(n,I_min),Qin2(n,I_min),Qin3(n,I_min),Qin4(n,I_min),&
                             Qin5(n,I_min),Qin6(n,I_min) )
          qnext_max    = max(Qin1(n,I_max),Qin2(n,I_max),Qin3(n,I_max),Qin4(n,I_max),&
                             Qin5(n,I_max),Qin6(n,I_max) )
!fj       enddo
          tmp=qnext_min
          qnext_min=q(n,k,l)
          if( tmp .ne.  CNST_MAX_REAL ) qnext_min = tmp
          tmp=qnext_max
          qnext_max=q(n,k,l)
          if( tmp .ne. -CNST_MAX_REAL ) qnext_max = tmp
!fj       do n = nstart, nend
!fj          zerosw = 0.5D0 - sign(0.5D0,abs(Cout_sum(n))-CNST_EPS_ZERO) ! if Cout_sum = 0, sw = 1
          zerosw = 0.5_RP - sign(0.5_RP,abs(Cout_sum)-CNST_EPS_ZERO) ! if Cout_sum = 0, sw = 1

!fj          Qout(n,k,l,I_min) = ( q(n,k,l) - CQin_max_sum(n) - qnext_max(n)*(1.0_RP-Cin_sum(n)-Cout_sum(n)+d(n,k,l)) ) &
!fj                            / ( Cout_sum(n) + zerosw ) * ( 1.0_RP - zerosw )                                         &
!fj                            + q(n,k,l) * zerosw
!fj          Qout(n,k,l,I_max) = ( q(n,k,l) - CQin_min_sum(n) - qnext_min(n)*(1.0_RP-Cin_sum(n)-Cout_sum(n)+d(n,k,l)) ) &
!fj                            / ( Cout_sum(n) + zerosw ) * ( 1.0_RP - zerosw )                                         &
!fj                            + q(n,k,l) * zerosw
          Qout(n,k,l,I_min) = ( q(n,k,l) - CQin_max_sum - qnext_max*(1.0_RP-Cin_sum-Cout_sum+d(n,k,l)) ) &
                        / ( Cout_sum + zerosw ) * ( 1.0_RP - zerosw )                                         &
                        + q(n,k,l) * zerosw
          Qout(n,k,l,I_max) = ( q(n,k,l) - CQin_min_sum - qnext_min*(1.0_RP-Cin_sum-Cout_sum+d(n,k,l)) ) &
                        / ( Cout_sum + zerosw ) * ( 1.0_RP - zerosw )                                         &
                        + q(n,k,l) * zerosw
       enddo ! n loop
!fj    enddo ! k loop
!fj    enddo ! l loop

!fj    Qout(     1:nstart-1,:,:,I_min) = q(     1:nstart-1,:,:)
!fj    Qout(nend+1:ADM_gall,:,:,I_min) = q(nend+1:ADM_gall,:,:)
!fj    Qout(     1:nstart-1,:,:,I_max) = q(     1:nstart-1,:,:)
!fj    Qout(nend+1:ADM_gall,:,:,I_max) = q(nend+1:ADM_gall,:,:)

!cx    call COMM_data_transfer( Qout(:,:,:,:), Qout_pl(:,:,:,:) )

    !---- apply inflow/outflow limiter
    nstart = suf(ADM_gmin-1,ADM_gmin-1)
    nend   = suf(ADM_gmax  ,ADM_gmax  )

!fj    do l = 1, ADM_lall
!fj    do k = 1, ADM_kall
       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d

!fj          q_a(1,n,k,l) = (      cmask(1,n,k,l) ) * min(max(q_a(1,n,k,l), Qin (1,ij    ,k,l,I_min)), Qin (1,ij    ,k,l,I_max)) &
!fj                       + ( 1.0_RP-cmask(1,n,k,l) ) * min(max(q_a(1,n,k,l), Qin (4,ip1j  ,k,l,I_min)), Qin (4,ip1j  ,k,l,I_max))
!fj          q_a(1,n,k,l) = (      cmask(1,n,k,l) ) * max(min(q_a(1,n,k,l), Qout(  ip1j  ,k,l,I_max)), Qout(  ip1j  ,k,l,I_min)) &
!fj                       + ( 1.0_RP-cmask(1,n,k,l) ) * max(min(q_a(1,n,k,l), Qout(  ij    ,k,l,I_max)), Qout(  ij    ,k,l,I_min))
!fj          q_a(4,ip1j,k,l) = q_a(1,n,k,l)
          q_a(1,n,k,l) = (      cmask(1,n,k,l) ) * min(max(q_a(1,n,k,l), Qin1(ij   ,I_min)), Qin1(ij    ,I_max)) &
                       + ( 1.0_RP-cmask(1,n,k,l) ) * min(max(q_a(1,n,k,l), Qin4(ip1j   ,I_min)), Qin4(ip1j    ,I_max))

          q_a(2,n,k,l) = (      cmask(2,n,k,l) ) * min(max(q_a(2,n,k,l), Qin2(ij,I_min)), Qin2(ij,I_max)) &
                       + ( 1.0_RP-cmask(2,n,k,l) ) * min(max(q_a(2,n,k,l), Qin5(ip1jp1,I_min)), Qin5(ip1jp1,I_max))
!fj          q_a(5,ip1jp1,k,l) = q_a(2,n,k,l)

!fj          q_a(3,n,k,l) = (      cmask(3,n,k,l) ) * min(max(q_a(3,n,k,l), Qin (3,ij    ,k,l,I_min)), Qin (3,ij    ,k,l,I_max)) &
!fj                       + ( 1.0_RP-cmask(3,n,k,l) ) * min(max(q_a(3,n,k,l), Qin (6,ijp1  ,k,l,I_min)), Qin (6,ijp1  ,k,l,I_max))
!fj          q_a(3,n,k,l) = (      cmask(3,n,k,l) ) * max(min(q_a(3,n,k,l), Qout(  ijp1  ,k,l,I_max)), Qout(  ijp1  ,k,l,I_min)) &
!fj                       + ( 1.0_RP-cmask(3,n,k,l) ) * max(min(q_a(3,n,k,l), Qout(  ij    ,k,l,I_max)), Qout(  ij    ,k,l,I_min))
          q_a(3,n,k,l) = (      cmask(3,n,k,l) ) * min(max(q_a(3,n,k,l), Qin3(ij    ,I_min)), Qin3(ij    ,I_max)) &
                       + ( 1.0_RP-cmask(3,n,k,l) ) * min(max(q_a(3,n,k,l), Qin6(ijp1    ,I_min)), Qin6(ijp1    ,I_max))
!fj          q_a(6,ijp1,k,l) = q_a(3,n,k,l)
       enddo
    enddo
    enddo
!fj>
    do l = 1, ADM_lall
    do k = 1, ADM_kall
       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d

          q_a(1,n,k,l) = (      cmask(1,n,k,l) ) * max(min(q_a(1,n,k,l), Qout(  ip1j  ,k,l,I_max)), Qout(  ip1j  ,k,l,I_min)) &
                       + ( 1.0_RP-cmask(1,n,k,l) ) * max(min(q_a(1,n,k,l), Qout(  ij    ,k,l,I_max)), Qout(  ij    ,k,l,I_min))
          q_a(4,ip1j,k,l) = q_a(1,n,k,l)

          q_a(2,n,k,l) = (      cmask(2,n,k,l) ) * max(min(q_a(2,n,k,l), Qout(  ip1jp1,k,l,I_max)), Qout(  ip1jp1,k,l,I_min)) &
                       + ( 1.0_RP-cmask(2,n,k,l) ) * max(min(q_a(2,n,k,l), Qout(  ij    ,k,l,I_max)), Qout(  ij    ,k,l,I_min))
          q_a(5,ip1jp1,k,l) = q_a(2,n,k,l)

          q_a(3,n,k,l) = (      cmask(3,n,k,l) ) * max(min(q_a(3,n,k,l), Qout(  ijp1  ,k,l,I_max)), Qout(  ijp1  ,k,l,I_min)) &
                       + ( 1.0_RP-cmask(3,n,k,l) ) * max(min(q_a(3,n,k,l), Qout(  ij    ,k,l,I_max)), Qout(  ij    ,k,l,I_min))
          q_a(6,ijp1,k,l) = q_a(3,n,k,l)
       enddo
    enddo
    enddo
 !fj<

    call DEBUG_rapend  ('____Horizontal_Adv_limiter')

    return
  end subroutine horizontal_limiter_thuburn


  subroutine H_Adv_limiter_snap_read &
            (q_a, q_a_pl, q, q_pl, d, d_pl, ch, ch_pl, cmask, cmask_pl)

    implicit none

    real(RP), intent(inout) :: q_a     (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: q_a_pl  (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: q       (  ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: q_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: d       (  ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: d_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: ch      (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: ch_pl   (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
!fj    real(8), intent(inout) :: cmask   (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: cmask_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !fj>
    real(RP), intent(inout) :: cmask   (3,ADM_gall   ,ADM_kall,ADM_lall   )
    real(8):: tmp_q_a     (6,ADM_gall   ,42,ADM_lall   )
    real(8):: tmp_q_a_pl  (  ADM_gall_pl,42,ADM_lall_pl)
    real(8):: tmp_q       (  ADM_gall   ,42,ADM_lall   )
    real(8):: tmp_q_pl    (  ADM_gall_pl,42,ADM_lall_pl)
    real(8):: tmp_d       (  ADM_gall   ,42,ADM_lall   )
    real(8):: tmp_d_pl    (  ADM_gall_pl,42,ADM_lall_pl)
    real(8):: tmp_ch      (6,ADM_gall   ,42,ADM_lall   )
    real(8):: tmp_ch_pl   (  ADM_gall_pl,42,ADM_lall_pl)
    real(8):: tmp_cmask   (6,ADM_gall   ,42,ADM_lall   )
    real(8):: tmp_cmask_pl(  ADM_gall_pl,42,ADM_lall_pl)
    !fj<
    integer i,l,k,n

    if (my_snapshot.eq.0) then
!fj        q_a     (6,1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
!fj        q_a     (1:6,1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
       q_a_pl  (  1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
!fj        q       (  1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        q_pl    (  1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
!fj        d       (  1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        d_pl    (  1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
!fj        ch      (6,1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
!fj        ch      (1:6,1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        ch_pl   (  1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
!fj        cmask   (6,1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
!fj        cmask   (1:6,1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        cmask_pl(  1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        call random_seed
        call random_number(q_a)
        call random_number(q)
        call random_number(d)
        call random_number(ch)
        call random_number(cmask)
    do l = 1, ADM_lall
    do k = 1, ADM_kall
       do n = 1,ADM_gall
       do i= 1,3
       if(cmask(i,n,k,l)>0.5_RP) then
          cmask(i,n,k,l)=1.0_RP
       else
          cmask(i,n,k,l)=0.0_RP
       endif
       enddo
       enddo
    enddo
    enddo

        return
    endif

    rewind(my_snapshot)
    do i=1,16
    read (my_snapshot) ! skip ADM records
    end do

!fj    read (my_snapshot) q_a
!fj    read (my_snapshot) q_a_pl
!fj    read (my_snapshot) q
!fj    read (my_snapshot) q_pl
!fj    read (my_snapshot) d
!fj    read (my_snapshot) d_pl
!fj    read (my_snapshot) ch
!fj    read (my_snapshot) ch_pl
!fj    read (my_snapshot) cmask
!fj    read (my_snapshot) cmask_pl
    read (my_snapshot) tmp_q_a
    read (my_snapshot) tmp_q_a_pl
    read (my_snapshot) tmp_q
    read (my_snapshot) tmp_q_pl
    read (my_snapshot) tmp_d
    read (my_snapshot) tmp_d_pl
    read (my_snapshot) tmp_ch
    read (my_snapshot) tmp_ch_pl
    read (my_snapshot) tmp_cmask
    read (my_snapshot) tmp_cmask_pl
    q_a(1:6,1:ADM_gall,1:42,1)= tmp_q_a(1:6,1:ADM_gall,1:42,1)
    q_a_pl  (  1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
    q_pl    (  1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
    d_pl    (  1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
    ch_pl   (  1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
    cmask_pl(  1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
!    q(1:ADM_gall,1:42,1)= tmp_q(1:ADM_gall,1:42,1)
!    d(1:ADM_gall,1:42,1)= tmp_d(1:ADM_gall,1:42,1)
!    ch(1:6,1:ADM_gall,1:42,1)= tmp_ch(1:6,1:ADM_gall,1:42,1)
!    cmask(1:6,1:ADM_gall,1:42,1)= tmp_cmask(1:6,1:ADM_gall,1:42,1)
!    q_a(1:6,1:ADM_gall,43:84,1)= tmp_q_a(1:6,1:ADM_gall,1:42,1)
!    q(1:ADM_gall,43:84,1)= tmp_q(1:ADM_gall,1:42,1)
!    d(1:ADM_gall,43:84,1)= tmp_d(1:ADM_gall,1:42,1)
!    ch(1:6,1:ADM_gall,43:84,1)= tmp_ch(1:6,1:ADM_gall,1:42,1)
!    cmask(1:6,1:ADM_gall,43:84,1)= tmp_cmask(1:6,1:ADM_gall,1:42,1)
!    q_a(1:6,1:ADM_gall,43:84,1)= tmp_q_a(1:6,1:ADM_gall,1:42,1)
!    q(1:ADM_gall,85:96,1)= tmp_q(1:ADM_gall,1:12,1)
!    d(1:ADM_gall,85:96,1)= tmp_d(1:ADM_gall,1:12,1)
!    ch(1:6,1:ADM_gall,85:96,1)= tmp_ch(1:6,1:ADM_gall,1:12,1)
!    cmask(1:6,1:ADM_gall,85:96,1)= tmp_cmask(1:6,1:ADM_gall,1:12,1)
    do l = 1, ADM_lall
    do k = 1, ADM_kall
       do n = 1,ADM_gall
       do i= 1,3
       if(cmask(i,n,k,l)>0.5_RP) then
          cmask(i,n,k,l)=1.0_RP
       else
          cmask(i,n,k,l)=0.0_RP
       endif
       enddo
       enddo
    enddo
    enddo
    do k= 1,ADM_kall
    if ( k <= 42) then
    l=k
    elseif( k <= 84) then
    l=k-42
    else
    l=k-84
    endif
    q_a(:,:,k,1)= tmp_q_a(:,:,l,1)
    q(:,k,1)= tmp_q(:,l,1)
    d(:,k,1)= tmp_d(:,l,1)
    ch(:,:,k,1)= tmp_ch(:,:,l,1)
    cmask(:,:,k,1)= tmp_cmask(1:3,:,l,1)
    enddo



    write(0,'(a,10i8)') "<H_Adv_limiter_snap_read>"
    return
  end subroutine H_Adv_limiter_snap_read


  subroutine H_Adv_limiter_snap_write &
            (q_a, q_a_pl, q, q_pl, d, d_pl, ch, ch_pl, cmask, cmask_pl)

    implicit none

    real(RP), intent(in) :: q_a     (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in) :: q_a_pl  (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in) :: q       (  ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in) :: q_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in) :: d       (  ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in) :: d_pl    (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in) :: ch      (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in) :: ch_pl   (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in) :: cmask   (6,ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in) :: cmask_pl(  ADM_gall_pl,ADM_kall,ADM_lall_pl)

    write(0,'(a,10i8)') "<H_Adv_limiter_snap_write>"

    write(my_snapshot) TI, TJ, AI, AIJ, AJ, K0
    write(my_snapshot) ADM_have_pl
    write(my_snapshot) ADM_lall
    write(my_snapshot) ADM_lall_pl
    write(my_snapshot) ADM_gall
    write(my_snapshot) ADM_gall_pl
    write(my_snapshot) ADM_kall
    write(my_snapshot) ADM_gall_1d
    write(my_snapshot) ADM_gmin
    write(my_snapshot) ADM_gmax
    write(my_snapshot) ADM_gslf_pl
    write(my_snapshot) ADM_gmin_pl
    write(my_snapshot) ADM_gmax_pl
    write(my_snapshot) ADM_kmin
    write(my_snapshot) ADM_kmax
    write(my_snapshot) ADM_have_sgp(1:ADM_lall)

    write(my_snapshot) q_a
    write(my_snapshot) q_a_pl
    write(my_snapshot) q
    write(my_snapshot) q_pl
    write(my_snapshot) d
    write(my_snapshot) d_pl
    write(my_snapshot) ch
    write(my_snapshot) ch_pl
    write(my_snapshot) cmask
    write(my_snapshot) cmask_pl

    return
  end subroutine H_Adv_limiter_snap_write

end module mod_src_tracer

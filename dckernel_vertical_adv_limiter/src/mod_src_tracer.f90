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
#ifdef DO_WE_NEED_EXTRA_MODULES
  use mod_grd, only: &
     XDIR => GRD_XDIR, &
     YDIR => GRD_YDIR, &
     ZDIR => GRD_ZDIR
  use mod_gmtr, only: &
     P_RAREA => GMTR_P_RAREA, &
     T_RAREA => GMTR_T_RAREA, &
     W1      => GMTR_T_W1,    &
     W2      => GMTR_T_W2,    &
     W3      => GMTR_T_W3,    &
     HNX     => GMTR_A_HNX,   &
     HNY     => GMTR_A_HNY,   &
     HNZ     => GMTR_A_HNZ,   &
     HTX     => GMTR_A_HTX,   &
     HTY     => GMTR_A_HTY,   &
     HTZ     => GMTR_A_HTZ,   &
     TNX     => GMTR_A_TNX,   &
     TNY     => GMTR_A_TNY,   &
     TNZ     => GMTR_A_TNZ,   &
     TN2X    => GMTR_A_TN2X,  &
     TN2Y    => GMTR_A_TN2Y,  &
     TN2Z    => GMTR_A_TN2Z
#endif

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

  implicit none
  private


  public :: vertical_limiter_thuburn
  public :: V_Adv_limiter_snap_write
  public :: V_Adv_limiter_snap_read


contains


  subroutine vertical_limiter_thuburn( &
       q_h, q_h_pl, &
       q,   q_pl,   &
       d,   d_pl,   &
       ck,  ck_pl   )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CNST_MAX_REAL, &
       CNST_EPS_ZERO
    use mod_snapshot_for_kernel
    use mod_adm, only: ADM_prc_me

    implicit none

    real(RP), intent(inout) :: q_h   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: q_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: q     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: q_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: d     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: d_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: ck    (ADM_gall,   ADM_kall,ADM_lall   ,2)
    real(RP), intent(in)    :: ck_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)

!fj    real(RP) :: Qin_min    (ADM_gall,   ADM_kall,2)
    real(RP) :: Qin_min_pl (ADM_gall_pl,ADM_kall,2)
!fj    real(RP) :: Qin_max    (ADM_gall,   ADM_kall,2)
    real(RP) :: Qin_max_pl (ADM_gall_pl,ADM_kall,2)

    real(RP) :: Qout_min   (ADM_gall,   ADM_kall)
    real(RP) :: Qout_min_pl(ADM_gall_pl,ADM_kall)
    real(RP) :: Qout_max   (ADM_gall,   ADM_kall)
    real(RP) :: Qout_max_pl(ADM_gall_pl,ADM_kall)

    real(RP) :: qnext_min, qnext_min_pl
    real(RP) :: qnext_max, qnext_max_pl
    real(RP) :: Cin,       Cin_pl
    real(RP) :: Cout,      Cout_pl
    real(RP) :: CQin_min,  CQin_min_pl
    real(RP) :: CQin_max,  CQin_max_pl

    real(RP) :: mask, mask1, mask2
    !fj>
    real(RP) :: mask3,qin_min1,qin_min2,qin_max1,qin_max2
    !fj<
    real(RP) :: zerosw

    integer :: n, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____Vertical_Adv_limiter')
!cx    if (mark_for_snapshots(ip_Vertical_Adv_limiter).ne.0) then
!cx    call snapshot_seq_bin_open ("Vertical_Adv_limiter", ADM_prc_me)
!cx    call V_Adv_limiter_snap_write (q_h, q_h_pl, q, q_pl, d, d_pl, ck, ck_pl)
!cx    call snapshot_seq_bin_close
!cx    mark_for_snapshots(ip_Vertical_Adv_limiter) = 0
!cx    endif


    do l = 1, ADM_lall
!fj       Qin_min(:,:,:) =  CNST_MAX_REAL
!fj       Qin_max(:,:,:) = -CNST_MAX_REAL
!fj
!fj       do k = ADM_kmin, ADM_kmax+1
!fj       do n = 1,ADM_gall
!fj          if ( ck(n,k,l,1) > 0.D0 ) then ! incoming
!fj             Qin_min(n,k-1,2) = min( q(n,k,l), q(n,k-1,l) )
!fj             Qin_max(n,k-1,2) = max( q(n,k,l), q(n,k-1,l) )
!fj          else ! outgoing
!fj             Qin_min(n,k  ,1) = min( q(n,k,l), q(n,k-1,l) )
!fj             Qin_max(n,k  ,1) = max( q(n,k,l), q(n,k-1,l) )
!fj          endif
!fj       enddo
!fj       enddo
!fj
!fj       do k = 1, ADM_kall
!fj>
       do k = 2, ADM_kall-1
!fj<
       do n = 1, ADM_gall
!fj>
          mask1 = 0.5_RP - sign(0.5_RP,ck(n,k,l,1)-CNST_EPS_ZERO)
          mask3 = 0.5_RP - sign(0.5_RP,ck(n,k+1,l,1)-CNST_EPS_ZERO)
          qin_min1 = (1.0_RP-mask1)*CNST_MAX_REAL +min(q(n,k,l),q(n,k-1,l))
          qin_min2 = mask3*CNST_MAX_REAL        +min(q(n,k,l),q(n,k+1,l))
          qin_max1 = -(1.0_RP-mask1)*CNST_MAX_REAL+max(q(n,k,l),q(n,k-1,l))
          qin_max2 = -mask3*CNST_MAX_REAL       +max(q(n,k,l),q(n,k+1,l))
!fj<

!fj          qnext_min = min( Qin_min(n,k,1), Qin_min(n,k,2) )
!fj>
          qnext_min = min(qin_min1,qin_min2,q(n,k,l))
          qnext_max = max(qin_max1,qin_max2,q(n,k,l))
!fj<
!fj          if( qnext_min ==  CNST_MAX_REAL ) qnext_min = q(n,k,l)

!fj          qnext_max = max( Qin_max(n,k,1), Qin_max(n,k,2) )
!fj          if( qnext_max == -CNST_MAX_REAL ) qnext_max = q(n,k,l)

!fj          mask1 = 0.5_RP - sign(0.5_RP,ck(n,k,l,1)-CNST_EPS_ZERO)
          mask2 = 0.5_RP - sign(0.5_RP,ck(n,k,l,2)-CNST_EPS_ZERO)

          Cin  = (      mask1 ) * ck(n,k,l,1) &
               + (      mask2 ) * ck(n,k,l,2)
          Cout = ( 1.0_RP-mask1 ) * ck(n,k,l,1) &
               + ( 1.0_RP-mask2 ) * ck(n,k,l,2)

!fj          CQin_max = mask1 * ( ck(n,k,l,1) * Qin_max(n,k,1) ) &
!fj                   + mask2 * ( ck(n,k,l,2) * Qin_max(n,k,2) )
!fj          CQin_min = mask1 * ( ck(n,k,l,1) * Qin_min(n,k,1) ) &
!fj                   + mask2 * ( ck(n,k,l,2) * Qin_min(n,k,2) )
!fj>
          CQin_max = mask1 * ( ck(n,k,l,1) * qin_max1 ) &
                   + mask2 * ( ck(n,k,l,2) * qin_max2 )
          CQin_min = mask1 * ( ck(n,k,l,1) * qin_min1 ) &
                   + mask2 * ( ck(n,k,l,2) * qin_min2 )
!fj<

          zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-CNST_EPS_ZERO) ! if Cout = 0, sw = 1

          Qout_min(n,k) = ( q(n,k,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d(n,k,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q(n,k,l) * zerosw
          Qout_max(n,k) = ( q(n,k,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d(n,k,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q(n,k,l) * zerosw

       enddo
       enddo
!fj>
       k=1
       do n = 1, ADM_gall
          mask1 = 0.5_RP - sign(0.5_RP,ck(n,k,l,1)-CNST_EPS_ZERO)
          mask3 = 0.5_RP - sign(0.5_RP,ck(n,k+1,l,1)-CNST_EPS_ZERO)
          qin_min2 = mask3*CNST_MAX_REAL        +min(q(n,k,l),q(n,k+1,l))
          qin_max2 = -mask3*CNST_MAX_REAL       +max(q(n,k,l),q(n,k+1,l))
          qnext_min = min(qin_min2,q(n,k,l))
          qnext_max = max(qin_max2,q(n,k,l))
          mask2 = 0.5_RP - sign(0.5_RP,ck(n,k,l,2)-CNST_EPS_ZERO)
          Cin  = (      mask1 ) * ck(n,k,l,1) &
               + (      mask2 ) * ck(n,k,l,2)
          Cout = ( 1.0_RP-mask1 ) * ck(n,k,l,1) &
               + ( 1.0_RP-mask2 ) * ck(n,k,l,2)
          CQin_max = mask1 * ( ck(n,k,l,1) * qin_max1 ) &
                   + mask2 * ( ck(n,k,l,2) * qin_max2 )
          CQin_min = mask1 * ( ck(n,k,l,1) * qin_min1 ) &
                   + mask2 * ( ck(n,k,l,2) * qin_min2 )
          zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-CNST_EPS_ZERO) ! if Cout = 0, sw = 1
          Qout_min(n,k) = ( q(n,k,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d(n,k,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q(n,k,l) * zerosw
          Qout_max(n,k) = ( q(n,k,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d(n,k,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q(n,k,l) * zerosw
       enddo
       k=ADM_kmax+1
       do n = 1, ADM_gall
          mask1 = 0.5_RP - sign(0.5_RP,ck(n,k,l,1)-CNST_EPS_ZERO)
          qin_min1 = (1.0_RP-mask1)*CNST_MAX_REAL +min(q(n,k,l),q(n,k-1,l))
          qin_max1 = -(1.0_RP-mask1)*CNST_MAX_REAL+max(q(n,k,l),q(n,k-1,l))
          qnext_min = min(qin_min1,q(n,k,l))
          qnext_max = max(qin_max1,q(n,k,l))
          mask2 = 0.5_RP - sign(0.5_RP,ck(n,k,l,2)-CNST_EPS_ZERO)
          Cin  = (      mask1 ) * ck(n,k,l,1) &
               + (      mask2 ) * ck(n,k,l,2)
          Cout = ( 1.0_RP-mask1 ) * ck(n,k,l,1) &
               + ( 1.0_RP-mask2 ) * ck(n,k,l,2)
          CQin_max = mask1 * ( ck(n,k,l,1) * qin_max1 ) &
                   + mask2 * ( ck(n,k,l,2) * qin_max2 )
          CQin_min = mask1 * ( ck(n,k,l,1) * qin_min1 ) &
                   + mask2 * ( ck(n,k,l,2) * qin_min2 )
          zerosw = 0.5_RP - sign(0.5_RP,abs(Cout)-CNST_EPS_ZERO) ! if Cout = 0, sw = 1
          Qout_min(n,k) = ( q(n,k,l) - CQin_max - qnext_max*(1.0_RP-Cin-Cout+d(n,k,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q(n,k,l) * zerosw
          Qout_max(n,k) = ( q(n,k,l) - CQin_min - qnext_min*(1.0_RP-Cin-Cout+d(n,k,l)) ) &
                        / ( Cout + zerosw ) * ( 1.0_RP - zerosw )                        &
                        + q(n,k,l) * zerosw
       enddo
!fj<

       do k = ADM_kmin, ADM_kmax+1
       do n = 1, ADM_gall
          mask = 0.5_RP - sign(0.5_RP,ck(n,k,l,1)-CNST_EPS_ZERO)

!fj          q_h(n,k,l) = (      mask ) * min( max( q_h(n,k,l), Qin_min(n,k  ,1) ), Qin_max(n,k  ,1) ) &
!fj                     + ( 1._RP-mask ) * min( max( q_h(n,k,l), Qin_min(n,k-1,2) ), Qin_max(n,k-1,2) )
          q_h(n,k,l) = min( max( q_h(n,k,l), min(q(n,k,l),q(n,k-1,l)) ), max(q(n,k,l),q(n,k-1,l)) )

          q_h(n,k,l) = (      mask ) * max( min( q_h(n,k,l), Qout_max(n,k-1) ), Qout_min(n,k-1) ) &
                     + ( 1.0_RP-mask ) * max( min( q_h(n,k,l), Qout_max(n,k  ) ), Qout_min(n,k  ) )
       enddo
       enddo

    enddo

    call DEBUG_rapend  ('____Vertical_Adv_limiter')

    return
  end subroutine vertical_limiter_thuburn


  subroutine V_Adv_limiter_snap_read (q_h, q_h_pl, q, q_pl, d, d_pl, ck, ck_pl)
    implicit none

    real(RP), intent(inout) :: q_h   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: q_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: q     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: q_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: d     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: d_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: ck    (ADM_gall,   ADM_kall,ADM_lall   ,2)
    real(RP), intent(inout) :: ck_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)
    !fj>
    real(8) :: tmp_q_h   (ADM_gall,   42,ADM_lall   )
    real(8) :: tmp_q_h_pl(ADM_gall_pl,42,ADM_lall_pl)
    real(8) :: tmp_q     (ADM_gall,   42,ADM_lall   )
    real(8) :: tmp_q_pl  (ADM_gall_pl,42,ADM_lall_pl)
    real(8) :: tmp_d     (ADM_gall,   42,ADM_lall   )
    real(8) :: tmp_d_pl  (ADM_gall_pl,42,ADM_lall_pl)
    real(8) :: tmp_ck    (ADM_gall,   42,ADM_lall   ,2)
    real(8) :: tmp_ck_pl (ADM_gall_pl,42,ADM_lall_pl,2)
    !fj<
    integer i,j,k,n,l

    if (my_snapshot.eq.0) then
!fj       do l = 1, ADM_lall
!fj       do k = 1, ADM_kall
!fj       do n = 1, ADM_gall
!fj        q_h   (1:ADM_gall,   1:ADM_kall,1:ADM_lall   ) = drand(0)-0.5d0
!fj        q     (1:ADM_gall,   1:ADM_kall,1:ADM_lall   ) = drand(0)-0.5d0
!fj        d     (1:ADM_gall,   1:ADM_kall,1:ADM_lall   ) = drand(0)-0.5d0
!fj        ck    (1:ADM_gall,   1:ADM_kall,1:ADM_lall   ,1) = drand(0)-0.5d0
!fj        ck    (1:ADM_gall,   1:ADM_kall,1:ADM_lall   ,2) = drand(0)-0.5d0
!fj        enddo
!fj        enddo
!fj        enddo
        call random_seed
        call random_number(q_h)
        call random_number(q)
        call random_number(d)
        call random_number(ck)
        return
    endif

    rewind(my_snapshot)
    do i=1,16
    read (my_snapshot) ! skip ADM records
    end do

!fj    read (my_snapshot) q_h
!fj    read (my_snapshot) q_h_pl
!fj    read (my_snapshot) q
!fj    read (my_snapshot) q_pl
!fj    read (my_snapshot) d
!fj    read (my_snapshot) d_pl
!fj    read (my_snapshot) ck
!fj    read (my_snapshot) ck_pl
    read (my_snapshot) tmp_q_h
    read (my_snapshot) tmp_q_h_pl
    read (my_snapshot) tmp_q
    read (my_snapshot) tmp_q_pl
    read (my_snapshot) tmp_d
    read (my_snapshot) tmp_d_pl
    read (my_snapshot) tmp_ck
    read (my_snapshot) tmp_ck_pl
    do k= 1,ADM_kall
    if ( k <= 42) then
    l=k
    elseif( k <= 84) then
    l=k-42
    else
    l=k-84
    endif
    q_h(:,k,1)= tmp_q_h(:,l,1)
    q(:,k,1)= tmp_q(:,l,1)
    d(:,k,1)= tmp_d(:,l,1)
    ck(:,k,1,1)= tmp_ck(:,l,1,1)
    ck(:,k,1,2)= tmp_ck(:,l,1,2)
    enddo

    write (0,'(a,10i8)') "<V_Adv_limiter_snap_read>"
    return
  end subroutine V_Adv_limiter_snap_read


  subroutine V_Adv_limiter_snap_write (q_h, q_h_pl, q, q_pl, d, d_pl, ck, ck_pl)

    implicit none
    real(RP), intent(inout) :: q_h   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: q_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: q     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: q_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: d     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: d_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: ck    (ADM_gall,   ADM_kall,ADM_lall   ,2)
    real(RP), intent(in)    :: ck_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,2)

    write(0,'(a,10i8)') "<V_Adv_limiter_snap_write>"

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

    write(my_snapshot) q_h
    write(my_snapshot) q_h_pl
    write(my_snapshot) q
    write(my_snapshot) q_pl
    write(my_snapshot) d
    write(my_snapshot) d_pl
    write(my_snapshot) ck
    write(my_snapshot) ck_pl

    write (0,'(a,10i8)') "<V_Adv_limiter_snap_write>"
    return
  end subroutine V_Adv_limiter_snap_write

end module mod_src_tracer

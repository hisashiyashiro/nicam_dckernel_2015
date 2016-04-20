!  This module is for the vertical implicit scheme of non-hydorostatic model.
module mod_vi

    use mod_precision
    use mod_adm, only: &
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
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_pl
    use mod_grd, only: &
       GRD_rdgzh, &
       GRD_afac,  &
       GRD_bfac
    use mod_vmtr, only: &
       VMTR_RGSGAM2,     &
       VMTR_RGSGAM2_pl,  &
       VMTR_RGSGAM2H,    &
       VMTR_RGSGAM2H_pl, &
       VMTR_RGAMH,       &
       VMTR_RGAMH_pl,    &
       VMTR_RGAM,        &
       VMTR_RGAM_pl,     &
       VMTR_GSGAM2H,     &
       VMTR_GSGAM2H_pl
    use mod_cnst, only: &
       GRAV  => CNST_EGRAV, &
       Rdry  => CNST_RAIR,  &
       CVdry => CNST_CV
    use mod_snapshot_for_kernel, only : my_snapshot
    use mod_debug

  implicit none
  private

  !++ Public procedure
  public :: vi_setup
  public :: vi_rhow_solver
  public :: vi_rhow_snap_write
  public :: vi_rhow_snap_read

  !++ Private parameters & variables
  real(RP), private, allocatable, save :: Mc   (:,:,:)
  real(RP), private, allocatable, save :: Mc_pl(:,:,:)
  real(RP), private, allocatable, save :: Ml   (:,:,:)
  real(RP), private, allocatable, save :: Ml_pl(:,:,:)
  real(RP), private, allocatable, save :: Mu   (:,:,:)
  real(RP), private, allocatable, save :: Mu_pl(:,:,:)

contains


  subroutine vi_setup
    use mod_adm, only: &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    implicit none
    allocate( Mc   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Mc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Mu   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Mu_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Ml   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Ml_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    return
  end subroutine vi_setup


  subroutine vi_rhow_solver( &
       rhogw,  rhogw_pl,  & !--- [INOUT]
       rhogw0, rhogw0_pl, & !--- [IN]
       preg0,  preg0_pl,  & !--- [IN]
       rhog0,  rhog0_pl,  & !--- [IN]
       Sr,     Sr_pl,     & !--- [IN]
       Sw,     Sw_pl,     & !--- [IN]
       Sp,     Sp_pl,     & !--- [IN]
       dt                 )

    implicit none
    real(RP), intent(inout) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 ), n+1
    real(RP), intent(inout) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)    :: rhogw0   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhogw0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: preg0    (ADM_gall   ,ADM_kall,ADM_lall   ) ! pressure prime ( G^1/2 x gam2 )
    real(RP), intent(in)    :: preg0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhog0    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho            ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhog0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sr       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for rho  at the full level
    real(RP), intent(in)    :: Sr_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sw       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for rhow at the half level
    real(RP), intent(in)    :: Sw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sp       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for pres at the full level
    real(RP), intent(in)    :: Sp_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: dt

    real(RP) :: Sall    (ADM_gall,   ADM_kall)
    real(RP) :: Sall_pl (ADM_gall_pl,ADM_kall)
    real(RP) :: beta    (ADM_gall   )
    real(RP) :: beta_pl (ADM_gall_pl)
    real(RP) :: gamma   (ADM_gall,   ADM_kall)
    real(RP) :: gamma_pl(ADM_gall_pl,ADM_kall)

    real(RP) :: alfa
    real(RP) :: CVovRt2   ! Cv / R / dt**2

    integer :: g, k, l
    !---------------------------------------------------------------------------
!cx    if (mark_for_snapshots(ip_vi_rhow_solver).ne.0) then
!cx    call snapshot_seq_bin_open ("vi_rhow_solver", ADM_prc_me)
!cx    call vi_rhow_snap_write ( &
!cx       rhogw,  rhogw_pl,  & !--- [INOUT]
!cx       rhogw0, rhogw0_pl, & !--- [IN]
!cx       preg0,  preg0_pl,  & !--- [IN]
!cx       rhog0,  rhog0_pl,  & !--- [IN]
!cx       Sr,     Sr_pl,     & !--- [IN]
!cx       Sw,     Sw_pl,     & !--- [IN]
!cx       Sp,     Sp_pl,     & !--- [IN]
!cx       dt                 )
!cx    call snapshot_seq_bin_close
!cx    mark_for_snapshots(ip_vi_rhow_solver) = 0
!cx    endif

    call DEBUG_rapstart('____vi_rhow_solver')

    alfa = 1.0_RP ! if hydrostatic, alfa=0
    CVovRt2 = CVdry / Rdry / (dt*dt)


    do l = 1, ADM_lall
       !--- < calc Sall > ---
       do k  = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          Sall(g,k) = (   ( rhogw0(g,k,  l)*alfa + dt * Sw(g,k,  l) ) * VMTR_RGAMH  (g,k,  l)**2            &
                      - ( ( preg0 (g,k,  l)      + dt * Sp(g,k,  l) ) * VMTR_RGSGAM2(g,k,  l)               &
                        - ( preg0 (g,k-1,l)      + dt * Sp(g,k-1,l) ) * VMTR_RGSGAM2(g,k-1,l)               &
                        ) * dt * GRD_rdgzh(k)                                                               &
                      - ( ( rhog0 (g,k,  l)      + dt * Sr(g,k,  l) ) * VMTR_RGAM(g,k,  l)**2 * GRD_afac(k) &
                        + ( rhog0 (g,k-1,l)      + dt * Sr(g,k-1,l) ) * VMTR_RGAM(g,k-1,l)**2 * GRD_bfac(k) &
                        ) * dt * 0.5_RP * GRAV                                                               &
                      ) * CVovRt2
       enddo
       enddo

       !--- boundary conditions
       do g = 1, ADM_gall
          rhogw(g,ADM_kmin,  l) = rhogw(g,ADM_kmin,  l) * VMTR_RGSGAM2H(g,ADM_kmin,  l)
          rhogw(g,ADM_kmax+1,l) = rhogw(g,ADM_kmax+1,l) * VMTR_RGSGAM2H(g,ADM_kmax+1,l)
          Sall (g,ADM_kmin+1)   = Sall (g,ADM_kmin+1) - Ml(g,ADM_kmin+1,l) * rhogw(g,ADM_kmin,  l)
          Sall (g,ADM_kmax  )   = Sall (g,ADM_kmax  ) - Mu(g,ADM_kmax,  l) * rhogw(g,ADM_kmax+1,l)
       enddo

       !--- < solve tri-daigonal matrix > ---

       ! condition at ADM_kmin+1
       k = ADM_kmin+1
       do g = 1, ADM_gall
          beta (g)     = Mc(g,k,l)
          rhogw(g,k,l) = Sall(g,k) / beta(g)
       enddo

       !--- forward
       do k = ADM_kmin+2, ADM_kmax
       do g = 1, ADM_gall
          gamma(g,k)   = Mu(g,k-1,l) / beta(g)
          beta (g)     = Mc(g,k,l) - Ml(g,k,l) * gamma(g,k) ! update beta
          rhogw(g,k,l) = ( Sall(g,k) - Ml(g,k,l) * rhogw(g,k-1,l) ) / beta(g)
       enddo
       enddo

       !--- backward
       do k = ADM_kmax-1, ADM_kmin+1, -1
       do g = 1, ADM_gall
          rhogw(g,k,l) = rhogw(g,k,l) - gamma(g,k+1) * rhogw(g,k+1,l)
          rhogw(g,k+1,l) = rhogw(g,k+1,l) * VMTR_GSGAM2H(g,k+1,l)
       enddo
       enddo

       !--- return value ( G^1/2 x gam2 )
!fj       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall
!fj          rhogw(g,k,l) = rhogw(g,k,l) * VMTR_GSGAM2H(g,k,l)
!fj>
          rhogw(g,adm_kmax+1,l) = rhogw(g,adm_kmax+1,l) * VMTR_GSGAM2H(g,adm_kmax+1,l)
          rhogw(g,adm_kmin,l) = rhogw(g,adm_kmin,l) * VMTR_GSGAM2H(g,adm_kmin,l)
          rhogw(g,adm_kmin+1,l) = rhogw(g,adm_kmin+1,l) * VMTR_GSGAM2H(g,adm_kmin+1,l)
!fj       enddo
!fj<
       enddo
    enddo

    call DEBUG_rapend('____vi_rhow_solver')

    return
  end subroutine vi_rhow_solver


  subroutine vi_rhow_snap_read ( &
       rhogw,  rhogw_pl,  & !--- [INOUT]
       rhogw0, rhogw0_pl, & !--- [IN]
       preg0,  preg0_pl,  & !--- [IN]
       rhog0,  rhog0_pl,  & !--- [IN]
       Sr,     Sr_pl,     & !--- [IN]
       Sw,     Sw_pl,     & !--- [IN]
       Sp,     Sp_pl,     & !--- [IN]
       dt                 )

    implicit none
    real(RP), intent(inout)    :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)    :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)    :: rhogw0   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)    :: rhogw0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)    :: preg0    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)    :: preg0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)    :: rhog0    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)    :: rhog0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)    :: Sr       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)    :: Sr_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)    :: Sw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)    :: Sw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)    :: Sp       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)    :: Sp_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)    :: dt
    real(RP) :: alfa
    integer i

    if (my_snapshot.eq.0) then
        rhogw    (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        rhogw_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        rhogw0   (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        rhogw0_pl(1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        preg0    (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        preg0_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        rhog0    (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        rhog0_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        Sr       (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        Sr_pl    (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        Sw       (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        Sw_pl    (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        Sp       (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        Sp_pl    (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
!fj>
       VMTR_RGSGAM2(:,:,:)=1.0
       VMTR_RGSGAM2H(:,:,:)=1.0
       VMTR_RGAMH(:,:,:)=1.0
       VMTR_RGAM(:,:,:)=1.0
       VMTR_GSGAM2H(:,:,:)=1.0
       alfa=1.0
       Mu(:,:,:)=1.0
       Mc(:,:,:)=1.0
       ML(:,:,:)=1.0
!fj<
        dt =1.0
        return
    endif

    rewind(my_snapshot)
    do i=1,16
    read (my_snapshot) ! skip ADM records
    end do

    read (my_snapshot) rhogw
    read (my_snapshot) rhogw_pl
    read (my_snapshot) rhogw0
    read (my_snapshot) rhogw0_pl
    read (my_snapshot) preg0
    read (my_snapshot) preg0_pl
    read (my_snapshot) rhog0
    read (my_snapshot) rhog0_pl
    read (my_snapshot) Sr
    read (my_snapshot) Sr_pl
    read (my_snapshot) Sw
    read (my_snapshot) Sw_pl
    read (my_snapshot) Sp
    read (my_snapshot) Sp_pl
    read (my_snapshot) dt

    read (my_snapshot) Mc
    read (my_snapshot) Mc_pl
    read (my_snapshot) Mu
    read (my_snapshot) Mu_pl
    read (my_snapshot) Ml
    read (my_snapshot) Ml_pl

    read (my_snapshot) GRD_rdgzh
    read (my_snapshot) GRD_afac
    read (my_snapshot) GRD_bfac

    read (my_snapshot) VMTR_RGSGAM2
    read (my_snapshot) VMTR_RGSGAM2_pl
    read (my_snapshot) VMTR_RGSGAM2H
    read (my_snapshot) VMTR_RGSGAM2H_pl
    read (my_snapshot) VMTR_RGAMH
    read (my_snapshot) VMTR_RGAMH_pl
    read (my_snapshot) VMTR_RGAM
    read (my_snapshot) VMTR_RGAM_pl
    read (my_snapshot) VMTR_GSGAM2H
    read (my_snapshot) VMTR_GSGAM2H_pl
    read (my_snapshot) alfa

    write(0,'(a,10i8)') "<vi_rhow_snap_read>"
        write(0,*) "shape(rhogw)=", shape(rhogw)
        write(0,*) "shape(Mc)=", shape(Mc)
        write(0,*) "shape(GRD_rdgzh)=", shape(GRD_rdgzh)
        write(0,*) "shape(VMTR_RGSGAM2)=", shape(VMTR_RGSGAM2)
    return
  end subroutine vi_rhow_snap_read

  subroutine vi_rhow_snap_write ( &
       rhogw,  rhogw_pl,  & !--- [INOUT]
       rhogw0, rhogw0_pl, & !--- [IN]
       preg0,  preg0_pl,  & !--- [IN]
       rhog0,  rhog0_pl,  & !--- [IN]
       Sr,     Sr_pl,     & !--- [IN]
       Sw,     Sw_pl,     & !--- [IN]
       Sp,     Sp_pl,     & !--- [IN]
       dt                 )

    implicit none
    real(RP), intent(in) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhogw0   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: rhogw0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: preg0    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: preg0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhog0    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: rhog0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sr       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: Sr_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: Sw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sp       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: Sp_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: dt
    real(RP) :: alfa


    write(0,'(a,10i8)') "<vi_rhow_snap_write>"

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

    write(my_snapshot) rhogw
    write(my_snapshot) rhogw_pl
    write(my_snapshot) rhogw0
    write(my_snapshot) rhogw0_pl
    write(my_snapshot) preg0
    write(my_snapshot) preg0_pl
    write(my_snapshot) rhog0
    write(my_snapshot) rhog0_pl
    write(my_snapshot) Sr
    write(my_snapshot) Sr_pl
    write(my_snapshot) Sw
    write(my_snapshot) Sw_pl
    write(my_snapshot) Sp
    write(my_snapshot) Sp_pl
    write(my_snapshot) dt

    !cx allocate( Mc   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    !cx allocate( Mc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    !cx allocate( Mu   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    !cx allocate( Mu_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    !cx allocate( Ml   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    !cx allocate( Ml_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )

    write(my_snapshot) Mc
    write(my_snapshot) Mc_pl
    write(my_snapshot) Mu
    write(my_snapshot) Mu_pl
    write(my_snapshot) Ml
    write(my_snapshot) Ml_pl

    write(my_snapshot) GRD_rdgzh
    write(my_snapshot) GRD_afac
    write(my_snapshot) GRD_bfac

    write(my_snapshot) VMTR_RGSGAM2
    write(my_snapshot) VMTR_RGSGAM2_pl
    write(my_snapshot) VMTR_RGSGAM2H
    write(my_snapshot) VMTR_RGSGAM2H_pl
    write(my_snapshot) VMTR_RGAMH
    write(my_snapshot) VMTR_RGAMH_pl
    write(my_snapshot) VMTR_RGAM
    write(my_snapshot) VMTR_RGAM_pl
    write(my_snapshot) VMTR_GSGAM2H
    write(my_snapshot) VMTR_GSGAM2H_pl
    write(my_snapshot) alfa

    return
  end subroutine vi_rhow_snap_write

end module mod_vi

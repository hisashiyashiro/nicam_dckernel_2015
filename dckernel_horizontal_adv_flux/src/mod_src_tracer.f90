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
module mod_src_tracer
  !-----------------------------------------------------------------------------

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
       ADM_gmax_pl,  &
       ADM_kmin,     &
       ADM_kmax
    use mod_grd, only: &
       GRD_xr,   &
       GRD_xr_pl
    use mod_gmtr, only: &
       GMTR_P_var_pl, &
       GMTR_T_var,    &
       GMTR_T_var_pl, &
       GMTR_A_var_pl
    use mod_oprt, only: &
       cinterp_HN,  &
       cinterp_PRA
    use mod_cnst, only: &
       CNST_EPS_ZERO
    use mod_snapshot_for_kernel
    use mod_adm, only: ADM_prc_me

  implicit none
  private

  public :: horizontal_flux
  public :: H_Adv_flux_snap_write
  public :: H_Adv_flux_snap_read


contains


  !-----------------------------------------------------------------------------
  !> prepare horizontal advection trem: mass flux, GRD_xc
  subroutine horizontal_flux( &
       flx_h,  flx_h_pl,  &
       GRD_xc, GRD_xc_pl, &
       rho,    rho_pl,    &
       rhovx,  rhovx_pl,  &
       rhovy,  rhovy_pl,  &
       rhovz,  rhovz_pl,  &
       dt                 )

    implicit none

    real(RP), intent(out) :: flx_h    (6,ADM_gall   ,ADM_kall,ADM_lall   )               ! horizontal mass flux
    real(RP), intent(out) :: flx_h_pl (  ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: GRD_xc   (ADM_gall   ,ADM_kall,ADM_lall   ,AI:AJ,XDIR:ZDIR) ! mass centroid position
    real(RP), intent(out) :: GRD_xc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,      XDIR:ZDIR)
    real(RP), intent(in)  :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )                 ! rho at cell center
    real(RP), intent(in)  :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: dt

!fj    real(8) :: rhot     (ADM_gall   ,TI:TJ) ! rho at cell vertex
    real(RP) :: rhot1     (ADM_gall ) ! rho at cell vertex
    real(RP) :: rhot2     (ADM_gall ) ! rho at cell vertex
    real(RP) :: rhot_pl  (ADM_gall_pl)
!fj    real(8) :: rhovxt   (ADM_gall   ,TI:TJ)
    real(RP) :: rhovxt1   (ADM_gall)
    real(RP) :: rhovxt2   (ADM_gall)
    real(RP) :: rhovxt_pl(ADM_gall_pl)
!fj    real(8) :: rhovyt   (ADM_gall   ,TI:TJ)
    real(RP) :: rhovyt1   (ADM_gall)
    real(RP) :: rhovyt2   (ADM_gall)
    real(RP) :: rhovyt_pl(ADM_gall_pl)
!fj    real(8) :: rhovzt   (ADM_gall   ,TI:TJ)
    real(RP) :: rhovzt1   (ADM_gall)
    real(RP) :: rhovzt2   (ADM_gall)
    real(RP) :: rhovzt_pl(ADM_gall_pl)
    !fj>
    real(RP) :: tmp_rhovxt
    real(RP) :: tmp_rhovyt
    real(RP) :: tmp_rhovzt
    !fj<

    real(RP) :: flux
    real(RP) :: rrhoa2

    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1

    integer :: nstart,nend
    integer :: n, k, l, v

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

!cx    if (mark_for_snapshots(ip_Horizontal_Adv_flux).ne.0) then
!cx    call snapshot_seq_bin_open ("Horizontal_Adv_flux", ADM_prc_me)
!cx    call H_Adv_flux_snap_write ( &
!cx       rho,    rho_pl,    &
!cx       rhovx,  rhovx_pl,  &
!cx       rhovy,  rhovy_pl,  &
!cx       rhovz,  rhovz_pl,  &
!cx       dt                 )
!cx    call snapshot_seq_bin_close
!cx    mark_for_snapshots(ip_Horizontal_Adv_flux) = 0
!cx    endif

    call DEBUG_rapstart('____Horizontal_Adv_flux')

    do l = 1, ADM_lall
    do k = 1, ADM_kall

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
!fj          ip1jp1 = n + 1 + ADM_gall_1d

!fj          rhot  (n,TI) = rho  (ij    ,k,l) * GMTR_T_var(n,K0,l,TI,W1) &
!fj                       + rho  (ip1j  ,k,l) * GMTR_T_var(n,K0,l,TI,W2) &
          rhot1  (n)   = rho  (ij    ,k,l) * GMTR_T_var(n,K0,l,TI,W1) &
                       + rho  (ip1j  ,k,l) * GMTR_T_var(n,K0,l,TI,W2) 
!fj                       + rho  (ip1jp1,k,l) * GMTR_T_var(n,K0,l,TI,W3)
!fj          rhovxt(n,TI) = rhovx(ij    ,k,l) * GMTR_T_var(n,K0,l,TI,W1) &
!fj                       + rhovx(ip1j  ,k,l) * GMTR_T_var(n,K0,l,TI,W2) &
          rhovxt1(n)   = rhovx(ij    ,k,l) * GMTR_T_var(n,K0,l,TI,W1) &
                       + rhovx(ip1j  ,k,l) * GMTR_T_var(n,K0,l,TI,W2) 
!fj                       + rhovx(ip1jp1,k,l) * GMTR_T_var(n,K0,l,TI,W3)
!fj          rhovyt(n,TI) = rhovy(ij    ,k,l) * GMTR_T_var(n,K0,l,TI,W1) &
!fj                       + rhovy(ip1j  ,k,l) * GMTR_T_var(n,K0,l,TI,W2) &
          rhovyt1(n)   = rhovy(ij    ,k,l) * GMTR_T_var(n,K0,l,TI,W1) &
                       + rhovy(ip1j  ,k,l) * GMTR_T_var(n,K0,l,TI,W2) 
!fj                       + rhovy(ip1jp1,k,l) * GMTR_T_var(n,K0,l,TI,W3)
!fj          rhovzt(n,TI) = rhovz(ij    ,k,l) * GMTR_T_var(n,K0,l,TI,W1) &
!fj                       + rhovz(ip1j  ,k,l) * GMTR_T_var(n,K0,l,TI,W2) &
          rhovzt1(n)   = rhovz(ij    ,k,l) * GMTR_T_var(n,K0,l,TI,W1) &
                       + rhovz(ip1j  ,k,l) * GMTR_T_var(n,K0,l,TI,W2) 
!fj                       + rhovz(ip1jp1,k,l) * GMTR_T_var(n,K0,l,TI,W3)
!fj>
          rhot2  (n) = rho  (ij    ,k,l) * GMTR_T_var(n,K0,l,TJ,W1) 
          rhovxt2(n) = rhovx(ij    ,k,l) * GMTR_T_var(n,K0,l,TJ,W1) 
          rhovyt2(n) = rhovy(ij    ,k,l) * GMTR_T_var(n,K0,l,TJ,W1) 
          rhovzt2(n) = rhovz(ij    ,k,l) * GMTR_T_var(n,K0,l,TJ,W1) 
!fj<
       enddo

       do n = nstart, nend
!fj          ij     = n
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d

!fj          rhot  (n,TJ) = rho  (ij    ,k,l) * GMTR_T_var(n,K0,l,TJ,W1) &
!fj                       + rho  (ip1jp1,k,l) * GMTR_T_var(n,K0,l,TJ,W2) &
!fj                       + rho  (ijp1  ,k,l) * GMTR_T_var(n,K0,l,TJ,W3)
!fj          rhovxt(n,TJ) = rhovx(ij    ,k,l) * GMTR_T_var(n,K0,l,TJ,W1) &
!fj                       + rhovx(ip1jp1,k,l) * GMTR_T_var(n,K0,l,TJ,W2) &
!fj                       + rhovx(ijp1  ,k,l) * GMTR_T_var(n,K0,l,TJ,W3)
!fj          rhovyt(n,TJ) = rhovy(ij    ,k,l) * GMTR_T_var(n,K0,l,TJ,W1) &
!fj                       + rhovy(ip1jp1,k,l) * GMTR_T_var(n,K0,l,TJ,W2) &
!fj                       + rhovy(ijp1  ,k,l) * GMTR_T_var(n,K0,l,TJ,W3)
!fj          rhovzt(n,TJ) = rhovz(ij    ,k,l) * GMTR_T_var(n,K0,l,TJ,W1) &
!fj                       + rhovz(ip1jp1,k,l) * GMTR_T_var(n,K0,l,TJ,W2) &
!fj                       + rhovz(ijp1  ,k,l) * GMTR_T_var(n,K0,l,TJ,W3)
!fj>
        rhot1(n)     = rhot1  (n) + rho(ip1jp1,k,l)   * GMTR_T_var(n,K0,l,TI,W3)
        rhovxt1(n)   = rhovxt1(n) + rhovx(ip1jp1,k,l) * GMTR_T_var(n,K0,l,TI,W3)
        rhovyt1(n)   = rhovyt1(n) + rhovy(ip1jp1,k,l) * GMTR_T_var(n,K0,l,TI,W3)
        rhovzt1(n)   = rhovzt1(n) + rhovz(ip1jp1,k,l) * GMTR_T_var(n,K0,l,TI,W3)
        rhot2 (n )   = rhot2 (n)  + rho  (ip1jp1,k,l) * GMTR_T_var(n,K0,l,TJ,W2) & 
                                  + rho  (ijp1,k,l)   * GMTR_T_var(n,K0,l,TJ,W3) 
        rhovxt2(n)   = rhovxt2(n) + rhovx(ip1jp1,k,l) * GMTR_T_var(n,K0,l,TJ,W2) & 
                                  + rhovx(ijp1,k,l)   * GMTR_T_var(n,K0,l,TJ,W3) 
        rhovyt2(n)   = rhovyt2(n) + rhovy(ip1jp1,k,l) * GMTR_T_var(n,K0,l,TJ,W2) & 
                                  + rhovy(ijp1,k,l)   * GMTR_T_var(n,K0,l,TJ,W3) 
        rhovzt2(n)   = rhovzt2(n) + rhovz(ip1jp1,k,l) * GMTR_T_var(n,K0,l,TJ,W2) &
                                  + rhovz(ijp1,k,l)   * GMTR_T_var(n,K0,l,TJ,W3) 
!fj<
       enddo

       if ( ADM_have_sgp(l) ) then
!fj          rhot  (suf(ADM_gmin-1,ADM_gmin-1),TI) = rhot  (suf(ADM_gmin,ADM_gmin-1),TJ)
!fj          rhovxt(suf(ADM_gmin-1,ADM_gmin-1),TI) = rhovxt(suf(ADM_gmin,ADM_gmin-1),TJ)
!fj          rhovyt(suf(ADM_gmin-1,ADM_gmin-1),TI) = rhovyt(suf(ADM_gmin,ADM_gmin-1),TJ)
!fj          rhovzt(suf(ADM_gmin-1,ADM_gmin-1),TI) = rhovzt(suf(ADM_gmin,ADM_gmin-1),TJ)
!fj>
          rhot1  (suf(ADM_gmin-1,ADM_gmin-1)) = rhot2  (suf(ADM_gmin,ADM_gmin-1))
          rhovxt1(suf(ADM_gmin-1,ADM_gmin-1)) = rhovxt2(suf(ADM_gmin,ADM_gmin-1))
          rhovyt1(suf(ADM_gmin-1,ADM_gmin-1)) = rhovyt2(suf(ADM_gmin,ADM_gmin-1))
          rhovzt1(suf(ADM_gmin-1,ADM_gmin-1)) = rhovzt2(suf(ADM_gmin,ADM_gmin-1))
!fj<
       endif

       !--- calculate flux and mass centroid position

       nstart = suf(ADM_gmin-1,ADM_gmin  )
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ijm1   = n     - ADM_gall_1d
          ip1j   = n + 1

!fj>
          tmp_rhovxt=rhovxt2(ijm1)+rhovxt1(ij)
          tmp_rhovyt=rhovyt2(ijm1)+rhovyt1(ij)
          tmp_rhovzt=rhovzt2(ijm1)+rhovzt1(ij)
!fj<
!fj          flux = 0.5D0 * ( (rhovxt(ij,TJ)+rhovxt(ij,TI)) * cinterp_HN(AI,1,ij,l) &
!fj                         + (rhovyt(ij,TJ)+rhovyt(ij,TI)) * cinterp_HN(AI,2,ij,l) &
!fj                         + (rhovzt(ij,TJ)+rhovzt(ij,TI)) * cinterp_HN(AI,3,ij,l) )
          flux = 0.5_RP * ( tmp_rhovxt * cinterp_HN(ij,l,AI ,1) &
                         + tmp_rhovyt * cinterp_HN(ij,l,AI ,2) &
                         + tmp_rhovzt * cinterp_HN(ij,l,AI ,3) )

          flx_h(1,ij  ,k,l) =  flux * cinterp_PRA(ij  ,l) * dt
          flx_h(4,ip1j,k,l) = -flux * cinterp_PRA(ip1j,l) * dt
!fj       enddo

!fj       do n = nstart, nend
!fj          ij     = n
!fj          ijm1   = n     - ADM_gall_1d

!fj          rrhoa2 = 1.D0 / max( rhot(ijm1,TJ) + rhot(ij,TI), CNST_EPS_ZERO ) ! doubled
          rrhoa2 = 1.0_RP / max( rhot2(ijm1) + rhot1(ij), CNST_EPS_ZERO ) ! doubled

!fj          GRD_xc(n,k,l,AI,XDIR) = GRD_xr(n,K0,l,AI,XDIR) - (rhovxt(ijm1,TJ)+rhovxt(ij,TI)) * rrhoa2 * dt * 0.5D0
!fj          GRD_xc(n,k,l,AI,YDIR) = GRD_xr(n,K0,l,AI,YDIR) - (rhovyt(ijm1,TJ)+rhovyt(ij,TI)) * rrhoa2 * dt * 0.5D0
!fj          GRD_xc(n,k,l,AI,ZDIR) = GRD_xr(n,K0,l,AI,ZDIR) - (rhovzt(ijm1,TJ)+rhovzt(ij,TI)) * rrhoa2 * dt * 0.5D0
          GRD_xc(n,k,l,AI,XDIR) = GRD_xr(n,K0,l,AI,XDIR) - tmp_rhovxt * rrhoa2 * dt * 0.5_RP
          GRD_xc(n,k,l,AI,YDIR) = GRD_xr(n,K0,l,AI,YDIR) - tmp_rhovyt * rrhoa2 * dt * 0.5_RP
          GRD_xc(n,k,l,AI,ZDIR) = GRD_xr(n,K0,l,AI,ZDIR) - tmp_rhovzt * rrhoa2 * dt * 0.5_RP
       enddo

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          ip1jp1 = n + 1 + ADM_gall_1d
          !fj>
          tmp_rhovxt=rhovxt1(ij)+rhovxt2(ij)
          tmp_rhovyt=rhovyt1(ij)+rhovyt2(ij)
          tmp_rhovzt=rhovzt1(ij)+rhovzt2(ij)
          !fj<

!fj          flux = 0.5D0 * ( (rhovxt(ij,TI)+rhovxt(ij,TJ)) * cinterp_HN(AIJ,1,ij,l) &
!fj                         + (rhovyt(ij,TI)+rhovyt(ij,TJ)) * cinterp_HN(AIJ,2,ij,l) &
!fj                         + (rhovzt(ij,TI)+rhovzt(ij,TJ)) * cinterp_HN(AIJ,3,ij,l) )
!fj>
          flux = 0.5_RP * ( tmp_rhovxt * cinterp_HN(ij,l,AIJ,1) &
                         + tmp_rhovyt * cinterp_HN(ij,l,AIJ,2) &
                         + tmp_rhovzt * cinterp_HN(ij,l,AIJ,3) )
!fj<

          flx_h(2,ij    ,k,l) =  flux * cinterp_PRA(ij    ,l) * dt
          flx_h(5,ip1jp1,k,l) = -flux * cinterp_PRA(ip1jp1,l) * dt
!fj       enddo

!fj       do n = nstart, nend
!fj          ij     = n

!fj          rrhoa2 = 1.D0 / max( rhot(ij,TI) + rhot(ij,TJ), CNST_EPS_ZERO ) ! doubled
          rrhoa2 = 1.0_RP / max( rhot1(ij) + rhot2(ij), CNST_EPS_ZERO ) ! doubled

!fj          GRD_xc(n,k,l,AIJ,XDIR) = GRD_xr(n,K0,l,AIJ,XDIR) - (rhovxt(ij,TI)+rhovxt(ij,TJ)) * rrhoa2 * dt * 0.5D0
!fj          GRD_xc(n,k,l,AIJ,YDIR) = GRD_xr(n,K0,l,AIJ,YDIR) - (rhovyt(ij,TI)+rhovyt(ij,TJ)) * rrhoa2 * dt * 0.5D0
!fj          GRD_xc(n,k,l,AIJ,ZDIR) = GRD_xr(n,K0,l,AIJ,ZDIR) - (rhovzt(ij,TI)+rhovzt(ij,TJ)) * rrhoa2 * dt * 0.5D0
          GRD_xc(n,k,l,AIJ,XDIR) = GRD_xr(n,K0,l,AIJ,XDIR) - tmp_rhovxt * rrhoa2 * dt * 0.5_RP
          GRD_xc(n,k,l,AIJ,YDIR) = GRD_xr(n,K0,l,AIJ,YDIR) - tmp_rhovyt * rrhoa2 * dt * 0.5_RP
          GRD_xc(n,k,l,AIJ,ZDIR) = GRD_xr(n,K0,l,AIJ,ZDIR) - tmp_rhovzt * rrhoa2 * dt * 0.5_RP
       enddo

       nstart = suf(ADM_gmin  ,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do n = nstart, nend
          ij     = n
          im1j   = n - 1
          ijp1   = n     + ADM_gall_1d
          !fj>
          tmp_rhovxt=rhovxt2(ij)+rhovxt1(im1j)
          tmp_rhovyt=rhovyt2(ij)+rhovyt1(im1j)
          tmp_rhovzt=rhovzt2(ij)+rhovzt1(im1j)
          !fj<

!fj          flux = 0.5D0 * ( (rhovxt(ij,TJ)+rhovxt(im1j,TI)) * cinterp_HN(AJ ,1,ij,l) &
!fj                         + (rhovyt(ij,TJ)+rhovyt(im1j,TI)) * cinterp_HN(AJ ,2,ij,l) &
!fj                         + (rhovzt(ij,TJ)+rhovzt(im1j,TI)) * cinterp_HN(AJ ,3,ij,l) )
          flux = 0.5_RP * ( tmp_rhovxt * cinterp_HN(ij,l,AJ ,1) &
                         + tmp_rhovyt * cinterp_HN(ij,l,AJ ,2) &
                         + tmp_rhovzt * cinterp_HN(ij,l,AJ ,3) )

          flx_h(3,ij  ,k,l) =  flux * cinterp_PRA(ij  ,l) * dt
          flx_h(6,ijp1,k,l) = -flux * cinterp_PRA(ijp1,l) * dt
!fj       enddo

!fj       do n = nstart, nend
!fj          ij     = n
!fj          im1j   = n - 1

!fj          rrhoa2 = 1.0_RP / max( rhot(ij,TJ) + rhot(im1j,TI), CNST_EPS_ZERO ) ! doubled
          rrhoa2 = 1.0_RP / max( rhot2(ij) + rhot1(im1j), CNST_EPS_ZERO ) ! doubled

!fj          GRD_xc(n,k,l,AJ,XDIR) = GRD_xr(n,K0,l,AJ,XDIR) - (rhovxt(ij,TJ)+rhovxt(im1j,TI)) * rrhoa2 * dt * 0.5_RP
!fj          GRD_xc(n,k,l,AJ,YDIR) = GRD_xr(n,K0,l,AJ,YDIR) - (rhovyt(ij,TJ)+rhovyt(im1j,TI)) * rrhoa2 * dt * 0.5_RP
!fj          GRD_xc(n,k,l,AJ,ZDIR) = GRD_xr(n,K0,l,AJ,ZDIR) - (rhovzt(ij,TJ)+rhovzt(im1j,TI)) * rrhoa2 * dt * 0.5_RP
          GRD_xc(n,k,l,AJ,XDIR) = GRD_xr(n,K0,l,AJ,XDIR) - tmp_rhovxt * rrhoa2 * dt * 0.5_RP
          GRD_xc(n,k,l,AJ,YDIR) = GRD_xr(n,K0,l,AJ,YDIR) - tmp_rhovyt * rrhoa2 * dt * 0.5_RP
          GRD_xc(n,k,l,AJ,ZDIR) = GRD_xr(n,K0,l,AJ,ZDIR) - tmp_rhovzt * rrhoa2 * dt * 0.5_RP
       enddo

       if ( ADM_have_sgp(l) ) then
          flx_h(6,suf(ADM_gmin,ADM_gmin),k,l) = 0.0_RP
       endif

    enddo
    enddo


    call DEBUG_rapend  ('____Horizontal_Adv_flux')

    return
  end subroutine horizontal_flux


  subroutine H_Adv_flux_snap_read ( &
       rho,    rho_pl,    &
       rhovx,  rhovx_pl,  &
       rhovy,  rhovy_pl,  &
       rhovz,  rhovz_pl,  &
       dt                 )

    implicit none
    real(RP), intent(inout)  :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)  :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)  :: rhovx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)  :: rhovx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)  :: rhovy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)  :: rhovy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)  :: rhovz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout)  :: rhovz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout)  :: dt
    integer i

!cx Better To Do: first touch the following arrays
    !cx allocate( cinterp_TN(AI:AJ,1:3,ADM_gall,ADM_lall) )
!fj    allocate( cinterp_HN(AI:AJ,1:3,ADM_gall,ADM_lall) )
    allocate( cinterp_HN(ADM_gall,ADM_lall,AI:AJ,1:3) )
    !cx allocate( cinterp_TRA(TI:TJ,ADM_gall,ADM_lall) )
    allocate( cinterp_PRA(      ADM_gall,ADM_lall) )


    if (my_snapshot.eq.0) then
        rho      (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        rho_pl   (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        rhovx    (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        rhovx_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        rhovy    (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        rhovy_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        rhovz    (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        rhovz_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        dt = 1.0
        !fj> add
        GRD_xr(:,:,:,:,:)=1.0
        cinterp_HN(:,:,:,:)=1.0
        cinterp_PRA(:,:)=1.0
        !fj< add
        return
    endif

    rewind(my_snapshot)
    do i=1,16
    read (my_snapshot) ! skip ADM records
    end do

    read (my_snapshot) rho
    read (my_snapshot) rho_pl
    read (my_snapshot) rhovx
    read (my_snapshot) rhovx_pl
    read (my_snapshot) rhovy
    read (my_snapshot) rhovy_pl
    read (my_snapshot) rhovz
    read (my_snapshot) rhovz_pl
    read (my_snapshot) dt

    read (my_snapshot) cinterp_HN
    read (my_snapshot) cinterp_PRA
    read (my_snapshot) GMTR_T_var
    read (my_snapshot) GMTR_A_var_pl
    read (my_snapshot) GMTR_P_var_pl
    read (my_snapshot) GMTR_T_var_pl
    read (my_snapshot) GRD_xr
    read (my_snapshot) GRD_xr_pl

    write(0,'(a,10i8)') "<H_Adv_flux_snap_read>"
    return
  end subroutine H_Adv_flux_snap_read


  subroutine H_Adv_flux_snap_write ( &
       rho,    rho_pl,    &
       rhovx,  rhovx_pl,  &
       rhovy,  rhovy_pl,  &
       rhovz,  rhovz_pl,  &
       dt                 )

    implicit none
    real(RP), intent(in)  :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhovz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhovz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: dt

    write(0,'(a,10i8)') "<H_Adv_flux_snap_write>"

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

    write(my_snapshot) rho
    write(my_snapshot) rho_pl
    write(my_snapshot) rhovx
    write(my_snapshot) rhovx_pl
    write(my_snapshot) rhovy
    write(my_snapshot) rhovy_pl
    write(my_snapshot) rhovz
    write(my_snapshot) rhovz_pl
    write(my_snapshot) dt

    write(my_snapshot) cinterp_HN
    write(my_snapshot) cinterp_PRA
    write(my_snapshot) GMTR_T_var
    write(my_snapshot) GMTR_A_var_pl
    write(my_snapshot) GMTR_P_var_pl
    write(my_snapshot) GMTR_T_var_pl
    write(my_snapshot) GRD_xr
    write(my_snapshot) GRD_xr_pl

    return
  end subroutine H_Adv_flux_snap_write

end module mod_src_tracer

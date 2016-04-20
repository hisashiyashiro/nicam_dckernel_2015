!-------------------------------------------------------------------------------
!>
!! 3D Operator module
!!
!! @par Description
!!         This module contains the subroutines for differential oeprators using vertical metrics.
!!
!! @li  chosen for performance evaluation targetting post-K
!!
!<
module mod_oprt3d
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
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
  use  mod_snapshot_for_kernel, only : &
    my_snapshot

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
    use mod_oprt, only: &
       OPRT_nstart, &
       OPRT_nend,   &
       cinterp_TN,  &
       cinterp_HN,  &
       cinterp_TRA, &
       cinterp_PRA
    use mod_grd, only: &
       GRD_rdgz
    use mod_vmtr, only: &
       VMTR_RGAM,        &
       VMTR_RGAMH,       &
       VMTR_RGSQRTH,     &
       VMTR_C2WfactGz

  implicit none
  private

  !++ Public procedure
  public :: OPRT3D_snap_read
  public :: OPRT3D_snap_write
  public :: OPRT3D_divdamp

  !++ Public parameters & variables
  !++ Private procedures
  !++ Private parameters & variables

contains

  subroutine OPRT3D_divdamp( ddivdx, ddivdy, ddivdz, &
                            rhogvx, rhogvy, rhogvz, rhogw  )

!cx Please remark that the number of arguments are reduced from the
!cx original NICAM-DC source which has the following API.
!cx    subroutine OPRT3D_divdamp( &
!cx       ddivdx, ddivdx_pl, &
!cx       ddivdy, ddivdy_pl, &
!cx       ddivdz, ddivdz_pl, &
!cx       rhogvx, rhogvx_pl, &
!cx       rhogvy, rhogvy_pl, &
!cx       rhogvz, rhogvz_pl, &
!cx       rhogw,  rhogw_pl   )

    use mod_snapshot_for_kernel
    use mod_adm, only: ADM_prc_me

    implicit none
    real(RP), intent(out) :: ddivdx   (ADM_gall   ,ADM_kall,ADM_lall   ) ! tendency
    real(RP), intent(out) :: ddivdy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*vx { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*vy { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*vz { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w  { gam2 x G^1/2 }

!fj    real(RP) :: sclt         (ADM_gall   ,ADM_kall,TI:TJ) ! scalar on the hexagon vertex
    real(RP) :: sclt         (ADM_gall   ,TI:TJ) ! scalar on the hexagon vertex
    real(RP) :: sclt_rhogw

!fj    real(8) :: rhogvx_vm   (ADM_gall   ,ADM_kall) ! rho*vx / vertical metrics
!fj    real(8) :: rhogvy_vm   (ADM_gall   ,ADM_kall) ! rho*vy / vertical metrics
!fj    real(8) :: rhogvz_vm   (ADM_gall   ,ADM_kall) ! rho*vz / vertical metrics
    real(RP) :: rhogvx_vm   (ADM_gall) ! rho*vx / vertical metrics
    real(RP) :: rhogvy_vm   (ADM_gall) ! rho*vy / vertical metrics
    real(RP) :: rhogvz_vm   (ADM_gall) ! rho*vz / vertical metrics
    real(RP) :: rhogw_vm    (ADM_gall,   ADM_kall) ! rho*w  / vertical metrics

    integer :: nstart, nend
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: g, k, l, v, n

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------
!cx    if (mark_for_snapshots(ip_OPRT3D_divdamp).ne.0) then
!cx    call snapshot_seq_bin_open ("OPRT3D_divdamp", ADM_prc_me)
!cx    call OPRT3D_snap_write( rhogvx, rhogvy, rhogvz, rhogw )
!cx    mark_for_snapshots(ip_OPRT3D_divdamp) = 0
!cx    endif
!fj>add
    OPRT_nstart = suf(ADM_gmin,ADM_gmin)
    OPRT_nend   = suf(ADM_gmax,ADM_gmax)
!fj<add


    call DEBUG_rapstart('OPRT3D_divdamp')

    do l = 1, ADM_lall
       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          rhogw_vm(g,k) = ( VMTR_C2WfactGz(1,g,k,l) * rhogvx(g,k  ,l) &
                          + VMTR_C2WfactGz(2,g,k,l) * rhogvx(g,k-1,l) &
                          + VMTR_C2WfactGz(3,g,k,l) * rhogvy(g,k  ,l) &
                          + VMTR_C2WfactGz(4,g,k,l) * rhogvy(g,k-1,l) &
                          + VMTR_C2WfactGz(5,g,k,l) * rhogvz(g,k  ,l) &
                          + VMTR_C2WfactGz(6,g,k,l) * rhogvz(g,k-1,l) &
                          ) * VMTR_RGAMH(g,k,l)                       & ! horizontal contribution
                        + rhogw(g,k,l) * VMTR_RGSQRTH(g,k,l)            ! vertical   contribution
       enddo
       enddo
       do g = 1, ADM_gall
          rhogw_vm(g,ADM_kmin  ) = 0.0_RP
          rhogw_vm(g,ADM_kmax+1) = 0.0_RP
       enddo
!fj>
       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax,  ADM_gmax  )
!fj<
!$omp parallel do private(rhogvx_vm,rhogvy_vm,rhogvz_vm,ij,ip1j,ip1jp1,sclt_rhogw,ijp1,im1j,sclt,ijm1,im1jm1)
       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall
!fj          rhogvx_vm(g,k) = rhogvx(g,k,l) * VMTR_RGAM(g,k,l)
!fj          rhogvy_vm(g,k) = rhogvy(g,k,l) * VMTR_RGAM(g,k,l)
!fj          rhogvz_vm(g,k) = rhogvz(g,k,l) * VMTR_RGAM(g,k,l)
          rhogvx_vm(g) = rhogvx(g,k,l) * VMTR_RGAM(g,k,l)
          rhogvy_vm(g) = rhogvy(g,k,l) * VMTR_RGAM(g,k,l)
          rhogvz_vm(g) = rhogvz(g,k,l) * VMTR_RGAM(g,k,l)
       enddo
!fj       enddo

!fj       nstart = suf(ADM_gmin-1,ADM_gmin-1)
!fj       nend   = suf(ADM_gmax,  ADM_gmax  )

!fj       do k = ADM_kmin, ADM_kmax
       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ip1jp1 = n + 1 + ADM_gall_1d

          sclt_rhogw = ( ( rhogw_vm(ij,k+1) + rhogw_vm(ip1j,k+1) + rhogw_vm(ip1jp1,k+1) ) &
                       - ( rhogw_vm(ij,k  ) + rhogw_vm(ip1j,k  ) + rhogw_vm(ip1jp1,k  ) ) &
                       ) / 3.0_RP * GRD_rdgz(k)

          sclt(n,TI) =   ( - (rhogvx_vm(ij    )+rhogvx_vm(ip1j  )) * cinterp_TN(AI ,1,ij  ,l) &
                           - (rhogvx_vm(ip1j  )+rhogvx_vm(ip1jp1)) * cinterp_TN(AJ ,1,ip1j,l) &
                           + (rhogvx_vm(ip1jp1)+rhogvx_vm(ij    )) * cinterp_TN(AIJ,1,ij  ,l) &
                           - (rhogvy_vm(ij    )+rhogvy_vm(ip1j  )) * cinterp_TN(AI ,2,ij  ,l) &
                           - (rhogvy_vm(ip1j  )+rhogvy_vm(ip1jp1)) * cinterp_TN(AJ ,2,ip1j,l) &
                           + (rhogvy_vm(ip1jp1)+rhogvy_vm(ij    )) * cinterp_TN(AIJ,2,ij  ,l) &
                           - (rhogvz_vm(ij    )+rhogvz_vm(ip1j  )) * cinterp_TN(AI ,3,ij  ,l) &
                           - (rhogvz_vm(ip1j  )+rhogvz_vm(ip1jp1)) * cinterp_TN(AJ ,3,ip1j,l) &
                           + (rhogvz_vm(ip1jp1)+rhogvz_vm(ij    )) * cinterp_TN(AIJ,3,ij  ,l) &
                         ) * 0.5_RP * cinterp_TRA(TI,ij,l) &
                       + sclt_rhogw
       enddo
!fj       enddo

!fj       do k = ADM_kmin, ADM_kmax
       do n = nstart, nend
          ij     = n
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d

          sclt_rhogw = ( ( rhogw_vm(ij,k+1) + rhogw_vm(ijp1,k+1) + rhogw_vm(ip1jp1,k+1) ) &
                       - ( rhogw_vm(ij,k  ) + rhogw_vm(ijp1,k  ) + rhogw_vm(ip1jp1,k  ) ) &
                       ) / 3.0_RP * GRD_rdgz(k)

          sclt(n,TJ) =   ( - (rhogvx_vm(ij    )+rhogvx_vm(ip1jp1)) * cinterp_TN(AIJ,1,ij  ,l) &
                           + (rhogvx_vm(ip1jp1)+rhogvx_vm(ijp1  )) * cinterp_TN(AI ,1,ijp1,l) &
                           + (rhogvx_vm(ijp1  )+rhogvx_vm(ij    )) * cinterp_TN(AJ ,1,ij  ,l) &
                           - (rhogvy_vm(ij    )+rhogvy_vm(ip1jp1)) * cinterp_TN(AIJ,2,ij  ,l) &
                           + (rhogvy_vm(ip1jp1)+rhogvy_vm(ijp1  )) * cinterp_TN(AI ,2,ijp1,l) &
                           + (rhogvy_vm(ijp1  )+rhogvy_vm(ij    )) * cinterp_TN(AJ ,2,ij  ,l) &
                           - (rhogvz_vm(ij    )+rhogvz_vm(ip1jp1)) * cinterp_TN(AIJ,3,ij  ,l) &
                           + (rhogvz_vm(ip1jp1)+rhogvz_vm(ijp1  )) * cinterp_TN(AI ,3,ijp1,l) &
                           + (rhogvz_vm(ijp1  )+rhogvz_vm(ij    )) * cinterp_TN(AJ ,3,ij  ,l) &
                         ) * 0.5_RP * cinterp_TRA(TJ,ij,l) &
                       + sclt_rhogw
       enddo
!fj       enddo

!fj       do k = ADM_kmin, ADM_kmax
       do n = OPRT_nstart, OPRT_nend
          ij     = n
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d
          im1jm1 = n - 1 - ADM_gall_1d

          ddivdx(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(AI ,1,ij,    l) &
                            + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(AIJ,1,ij,    l) &
                            + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(AJ ,1,ij,    l) &
                            - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(AI ,1,im1j,  l) &
                            - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(AIJ,1,im1jm1,l) &
                            - ( sclt(ijm1  ,TJ) + sclt(im1jm1,TI) ) * cinterp_HN(AJ ,1,ijm1,  l) &
                          ) * 0.5_RP * cinterp_PRA(ij,l)

          ddivdy(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(AI ,2,ij,    l) &
                            + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(AIJ,2,ij,    l) &
                            + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(AJ ,2,ij,    l) &
                            - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(AI ,2,im1j,  l) &
                            - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(AIJ,2,im1jm1,l) &
                            - ( sclt(ijm1  ,TJ) + sclt(im1jm1,TI) ) * cinterp_HN(AJ ,2,ijm1,  l) &
                          ) * 0.5_RP * cinterp_PRA(ij,l)

          ddivdz(n,k,l) = ( + ( sclt(ijm1,  TJ) + sclt(ij,    TI) ) * cinterp_HN(AI ,3,ij,    l) &
                            + ( sclt(ij,    TI) + sclt(ij,    TJ) ) * cinterp_HN(AIJ,3,ij,    l) &
                            + ( sclt(ij,    TJ) + sclt(im1j,  TI) ) * cinterp_HN(AJ ,3,ij,    l) &
                            - ( sclt(im1jm1,TJ) + sclt(im1j,  TI) ) * cinterp_HN(AI ,3,im1j,  l) &
                            - ( sclt(im1jm1,TI) + sclt(im1jm1,TJ) ) * cinterp_HN(AIJ,3,im1jm1,l) &
                            - ( sclt(ijm1  ,TJ) + sclt(im1jm1,TI) ) * cinterp_HN(AJ ,3,ijm1,  l) &
                          ) * 0.5_RP * cinterp_PRA(ij,l)
       enddo
       enddo

       do g = 1, ADM_gall
          ddivdx(g,ADM_kmin-1,l) = 0.0_RP
          ddivdx(g,ADM_kmax+1,l) = 0.0_RP
          ddivdy(g,ADM_kmin-1,l) = 0.0_RP
          ddivdy(g,ADM_kmax+1,l) = 0.0_RP
          ddivdz(g,ADM_kmin-1,l) = 0.0_RP
          ddivdz(g,ADM_kmax+1,l) = 0.0_RP
       enddo
    enddo

!cx cut off other loops

    call DEBUG_rapend('OPRT3D_divdamp')

    return
  end subroutine OPRT3D_divdamp


  subroutine OPRT3D_snap_read( rhogvx, rhogvy, rhogvz, rhogw )

    implicit none
    real(RP), intent(out)  :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    integer :: i

!cx Better To Do: first touch the following arrays
    allocate( cinterp_TN(AI:AJ,1:3,ADM_gall,ADM_lall) )
    allocate( cinterp_HN(AI:AJ,1:3,ADM_gall,ADM_lall) )
    allocate( cinterp_TRA(TI:TJ,ADM_gall,ADM_lall) )
    allocate( cinterp_PRA(      ADM_gall,ADM_lall) )


    if (my_snapshot.eq.0) then
        cinterp_TN(AI:AJ,1:3,1:ADM_gall,1:ADM_lall) = 1.0_RP
        cinterp_HN(AI:AJ,1:3,1:ADM_gall,1:ADM_lall) = 1.0_RP
        cinterp_TRA(TI:TJ,1:ADM_gall,1:ADM_lall) = 1.0_RP
        cinterp_PRA(      1:ADM_gall,1:ADM_lall) = 1.0_RP
        rhogvx (1:ADM_gall ,1:ADM_kall, 1:ADM_lall) = 1.0_RP
        rhogvy (1:ADM_gall ,1:ADM_kall, 1:ADM_lall) = 1.0_RP
        rhogvz (1:ADM_gall ,1:ADM_kall, 1:ADM_lall) = 1.0_RP
        rhogw  (1:ADM_gall ,1:ADM_kall, 1:ADM_lall) = 1.0_RP
        !fj> add
        VMTR_C2WfactGz(:,:,:,:)=1.0_RP
        VMTR_RGAM(:,:,:)=1.0_RP
        VMTR_RGAMH(:,:,:)=1.0_RP
        VMTR_RGSQRTH(:,:,:)=1.0_RP
        !fj< add
        return
    endif

    rewind(my_snapshot)
    do i=1,16
    read (my_snapshot) ! skip ADM records
    end do

    read (my_snapshot) cinterp_TN
    read (my_snapshot) cinterp_HN
    read (my_snapshot) cinterp_TRA
    read (my_snapshot) cinterp_PRA
    read (my_snapshot) rhogvx
    read (my_snapshot)
    read (my_snapshot) rhogvy
    read (my_snapshot)
    read (my_snapshot) rhogvz
    read (my_snapshot)
    read (my_snapshot) rhogw
    read (my_snapshot)

    read (my_snapshot)
    read (my_snapshot)
    read (my_snapshot)

    read (my_snapshot) GRD_rdgz
    read (my_snapshot) OPRT_nstart
    read (my_snapshot) OPRT_nend
    read (my_snapshot) VMTR_C2WfactGz
    read (my_snapshot)
    read (my_snapshot) VMTR_RGAM
    read (my_snapshot)
    read (my_snapshot) VMTR_RGAMH
    read (my_snapshot)
    read (my_snapshot) VMTR_RGSQRTH
    read (my_snapshot)

    write(0,'(a,10i8)') "<OPRT3D_snap_read>"
    return
  end subroutine OPRT3D_snap_read


  subroutine OPRT3D_snap_write( rhogvx, rhogvy, rhogvz, rhogw )

    use mod_grd, only: &
       GRD_rdgz
    use mod_oprt, only: &
       OPRT_nstart, &
       OPRT_nend,   &
       cinterp_TN,  &
       cinterp_HN,  &
       cinterp_TRA, &
       cinterp_PRA
    use mod_vmtr, only: &
       VMTR_RGAM,        &
       VMTR_RGAMH,       &
       VMTR_RGSQRTH,     &
       VMTR_C2WfactGz
    use mod_adm, only: ADM_prc_me
    use mod_snapshot_for_kernel, only : my_snapshot

    implicit none
    real(RP), intent(in)  :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )

    write(0,'(a,10i8)') "<OPRT3D_snap_write>"

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
    write(my_snapshot) cinterp_TN
    write(my_snapshot) cinterp_HN
    write(my_snapshot) cinterp_TRA
    write(my_snapshot) cinterp_PRA
    write(my_snapshot) rhogvx
    write(my_snapshot) rhogvy
    write(my_snapshot) rhogvz
    write(my_snapshot) rhogw

    write(my_snapshot) GRD_rdgz
    write(my_snapshot) OPRT_nstart
    write(my_snapshot) OPRT_nend
    write(my_snapshot) VMTR_C2WfactGz
    write(my_snapshot) VMTR_RGAM
    write(my_snapshot) VMTR_RGAMH
    write(my_snapshot) VMTR_RGSQRTH

    return
  end subroutine OPRT3D_snap_write

end module mod_oprt3d
!-------------------------------------------------------------------------------

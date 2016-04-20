!-------------------------------------------------------------------------------
!>
!! Operator module
!!
!! @par Description
!!         This module contains the subroutines for differential oeprators.
!!
!! @li  chosen for performance evaluation targetting post-K
!!
!<
!-------------------------------------------------------------------------------

module mod_oprt
  use mod_precision
  use mod_adm, only: &
    ADM_have_pl,    &   !cx logical
    ADM_have_sgp,   &   !cx logical array
    ADM_lall,       &
    ADM_lall_pl,    &
    ADM_gall,       &
    ADM_gall_pl,    &
    ADM_kall,       &
    ADM_gall_1d,    &
    ADM_gmin,       &
    ADM_gmax,       &
    ADM_gslf_pl,    &   !cx parameter
    ADM_gmin_pl,    &   !cx parameter
    ADM_gmax_pl,    &
    TI  => ADM_TI,  &
    TJ  => ADM_TJ,  &
    AI  => ADM_AI,  &
    AIJ => ADM_AIJ, &
    AJ  => ADM_AJ,  &
    K0  => ADM_KNONE

  use mod_adm, only: my_id => ADM_prc_me

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

  implicit none

  !++ Public procedure
  public :: OPRT_snap_read
  public :: OPRT_snap_write
  public :: OPRT_diffusion

  !++ Public parameters & variables
  integer, public, save :: OPRT_nstart
  integer, public, save :: OPRT_nend

  ! < for diffusion operator >
  real(RP), public, allocatable, save :: cinterp_TN (:,:,:,:)
  real(RP), public, allocatable, save :: cinterp_HN (:,:,:,:)
  real(RP), public, allocatable, save :: cinterp_TRA(:,:,:)
  real(RP), public, allocatable, save :: cinterp_PRA(:,:)

contains


  subroutine OPRT_diffusion( &
       dscl, dscl_pl, &
       scl,  scl_pl,  &
       kh,   kh_pl,   &
       mfact          )
    use mod_debug, only: DEBUG_rapstart, DEBUG_rapend
    implicit none

    real(RP), intent(out) :: dscl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: kh     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: mfact

    real(RP)  :: vxt    (ADM_gall   ,ADM_kall,TI:TJ)
    real(RP)  :: vyt    (ADM_gall   ,ADM_kall,TI:TJ)
    real(RP)  :: vzt    (ADM_gall   ,ADM_kall,TI:TJ)
    real(RP)  :: flux   (ADM_gall   ,ADM_kall,AI:AJ)
    real(RP)  :: vxt_pl (ADM_gall_pl,ADM_kall)
    real(RP)  :: vyt_pl (ADM_gall_pl,ADM_kall)
    real(RP)  :: vzt_pl (ADM_gall_pl,ADM_kall)
    real(RP)  :: flux_pl(ADM_gall_pl,ADM_kall)

    real(RP) :: u1, u2, u3, smean

    integer :: nstart, nend
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, k, l, v

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------
    !cx  write(0,'(a,10i8)') "<OPRT_diffusion>"

    call DEBUG_rapstart('OPRT_diffusion')

    do l = 1, ADM_lall
       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do k = 1, ADM_kall
       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ip1jp1 = n + 1 + ADM_gall_1d

          smean = ( scl(ij,k,l) + scl(ip1j,k,l) + scl(ip1jp1,k,l) ) / 3.D0

          u1 = 0.5D0 * (scl(ij    ,k,l)+scl(ip1j  ,k,l)) - smean
          u2 = 0.5D0 * (scl(ip1j  ,k,l)+scl(ip1jp1,k,l)) - smean
          u3 = 0.5D0 * (scl(ip1jp1,k,l)+scl(ij    ,k,l)) - smean

          vxt(n,k,TI) = ( - u1 * cinterp_TN(AI ,1,ij  ,l) &
                          - u2 * cinterp_TN(AJ ,1,ip1j,l) &
                          + u3 * cinterp_TN(AIJ,1,ij  ,l) ) * cinterp_TRA(TI,ij,l)
          vyt(n,k,TI) = ( - u1 * cinterp_TN(AI ,2,ij  ,l) &
                          - u2 * cinterp_TN(AJ ,2,ip1j,l) &
                          + u3 * cinterp_TN(AIJ,2,ij  ,l) ) * cinterp_TRA(TI,ij,l)
          vzt(n,k,TI) = ( - u1 * cinterp_TN(AI ,3,ij  ,l) &
                          - u2 * cinterp_TN(AJ ,3,ip1j,l) &
                          + u3 * cinterp_TN(AIJ,3,ij  ,l) ) * cinterp_TRA(TI,ij,l)
       enddo
       enddo

       do k = 1, ADM_kall
       do n = nstart, nend
          ij     = n
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d

          smean = ( scl(ij,k,l) + scl(ip1jp1,k,l) + scl(ijp1,k,l) ) / 3.D0

          u1 = 0.5D0 * (scl(ij    ,k,l)+scl(ip1jp1,k,l)) - smean
          u2 = 0.5D0 * (scl(ip1jp1,k,l)+scl(ijp1  ,k,l)) - smean
          u3 = 0.5D0 * (scl(ijp1  ,k,l)+scl(ij    ,k,l)) - smean

          vxt(n,k,TJ) = ( - u1 * cinterp_TN(AIJ,1,ij  ,l) &
                          + u2 * cinterp_TN(AI ,1,ijp1,l) &
                          + u3 * cinterp_TN(AJ ,1,ij  ,l) ) * cinterp_TRA(TJ,ij,l)
          vyt(n,k,TJ) = ( - u1 * cinterp_TN(AIJ,2,ij  ,l) &
                          + u2 * cinterp_TN(AI ,2,ijp1,l) &
                          + u3 * cinterp_TN(AJ ,2,ij  ,l) ) * cinterp_TRA(TJ,ij,l)
          vzt(n,k,TJ) = ( - u1 * cinterp_TN(AIJ,3,ij  ,l) &
                          + u2 * cinterp_TN(AI ,3,ijp1,l) &
                          + u3 * cinterp_TN(AJ ,3,ij  ,l) ) * cinterp_TRA(TJ,ij,l)
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          do k = 1, ADM_kall
             vxt(suf(ADM_gmin-1,ADM_gmin-1),k,TI) = vxt(suf(ADM_gmin,ADM_gmin-1),k,TJ)
             vyt(suf(ADM_gmin-1,ADM_gmin-1),k,TI) = vyt(suf(ADM_gmin,ADM_gmin-1),k,TJ)
             vzt(suf(ADM_gmin-1,ADM_gmin-1),k,TI) = vzt(suf(ADM_gmin,ADM_gmin-1),k,TJ)
          enddo
       endif

       nstart = suf(ADM_gmin-1,ADM_gmin  )
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do k = 1, ADM_kall
       do n = nstart, nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d

          flux(n,k,AI ) = 0.25D0 * ( (vxt(ijm1,k,TJ)+vxt(ij  ,k,TI)) * cinterp_HN(AI ,1,ij,l) &
                                   + (vyt(ijm1,k,TJ)+vyt(ij  ,k,TI)) * cinterp_HN(AI ,2,ij,l) &
                                   + (vzt(ijm1,k,TJ)+vzt(ij  ,k,TI)) * cinterp_HN(AI ,3,ij,l) &
                                   ) * (kh(ij,k,l)+kh(ip1j  ,k,l))
       enddo
       enddo

       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do k = 1, ADM_kall
       do n = nstart, nend
          ij     = n
          ip1jp1 = n + 1 + ADM_gall_1d

          flux(n,k,AIJ) = 0.25D0 * ( (vxt(ij  ,k,TI)+vxt(ij  ,k,TJ)) * cinterp_HN(AIJ,1,ij,l) &
                                   + (vyt(ij  ,k,TI)+vyt(ij  ,k,TJ)) * cinterp_HN(AIJ,2,ij,l) &
                                   + (vzt(ij  ,k,TI)+vzt(ij  ,k,TJ)) * cinterp_HN(AIJ,3,ij,l) &
                                   ) * (kh(ij,k,l)+kh(ip1jp1,k,l))
       enddo
       enddo

       nstart = suf(ADM_gmin  ,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

       do k = 1, ADM_kall
       do n = nstart, nend
          ij     = n
          ijp1   = n     + ADM_gall_1d
          im1j   = n - 1

          flux(n,k,AJ ) = 0.25D0 * ( (vxt(ij  ,k,TJ)+vxt(im1j,k,TI)) * cinterp_HN(AJ ,1,ij,l) &
                                   + (vyt(ij  ,k,TJ)+vyt(im1j,k,TI)) * cinterp_HN(AJ ,2,ij,l) &
                                   + (vzt(ij  ,k,TJ)+vzt(im1j,k,TI)) * cinterp_HN(AJ ,3,ij,l) &
                                   ) * (kh(ij,k,l)+kh(ijp1  ,k,l))
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          do k = 1, ADM_kall
             flux(suf(ADM_gmin,ADM_gmin-1),k,AJ) = 0.D0
          enddo
       endif

       do k = 1, ADM_kall
       do n = OPRT_nstart, OPRT_nend
          ij     = n
          im1j   = n - 1
          im1jm1 = n - 1 - ADM_gall_1d
          ijm1   = n     - ADM_gall_1d

          dscl(n,k,l) = ( flux(ij,k,AI ) - flux(im1j  ,k,AI ) &
                        + flux(ij,k,AIJ) - flux(im1jm1,k,AIJ) &
                        + flux(ij,k,AJ ) - flux(ijm1  ,k,AJ ) ) * cinterp_PRA(ij,l) * mfact
       enddo
       enddo

    enddo

!cx cut off other loops

    call DEBUG_rapend('OPRT_diffusion')

    return
  end subroutine OPRT_diffusion


  subroutine OPRT_snap_read ( scl, scl_pl, kh, kh_pl, mfact )
    implicit none
    integer :: i

    real(RP), intent(out)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out)  :: kh     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out)  :: mfact

!cx Better To Do: first touch the following arrays
    allocate( cinterp_TN(AI:AJ,1:3,ADM_gall,ADM_lall) )
    allocate( cinterp_HN(AI:AJ,1:3,ADM_gall,ADM_lall) )
    allocate( cinterp_TRA(TI:TJ,ADM_gall,ADM_lall) )
    allocate( cinterp_PRA(      ADM_gall,ADM_lall) )

    if (my_snapshot.eq.0) then
        cinterp_TN(AI:AJ,1:3,1:ADM_gall,1:ADM_lall) = 1.0
        cinterp_HN(AI:AJ,1:3,1:ADM_gall,1:ADM_lall) = 1.0
        cinterp_TRA(TI:TJ,1:ADM_gall,1:ADM_lall) = 1.0
        cinterp_PRA(      1:ADM_gall,1:ADM_lall) = 1.0
        scl    (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        scl_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        kh     (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0
        kh_pl  (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0
        mfact = 1.0
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
    read (my_snapshot) scl
    read (my_snapshot) scl_pl
    read (my_snapshot) kh
    read (my_snapshot) kh_pl
    read (my_snapshot) mfact

    write(0,'(a,10i8)') "<OPRT_snap_read> snapshot file was read"
    return
  end subroutine OPRT_snap_read


  subroutine OPRT_snap_write ( scl, scl_pl, kh, kh_pl, mfact )
    implicit none
    real(RP), intent(in)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: kh     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: mfact

    if (my_snapshot.eq.0) then
        return
    endif
    write(0,'(a,10i8)') "<OPRT_snap_write>"

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

    write(my_snapshot) ADM_have_sgp
    write(my_snapshot) cinterp_TN
    write(my_snapshot) cinterp_HN
    write(my_snapshot) cinterp_TRA
    write(my_snapshot) cinterp_PRA
    write(my_snapshot) scl
    write(my_snapshot) scl_pl
    write(my_snapshot) kh
    write(my_snapshot) kh_pl
    write(my_snapshot) mfact
    return
  end subroutine OPRT_snap_write

end module mod_oprt


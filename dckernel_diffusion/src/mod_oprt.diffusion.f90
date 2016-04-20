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


  subroutine OPRT_snap_read ( scl, scl_pl, kh, kh_pl, mfact )
    implicit none
    integer :: i

    real(RP), intent(out)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out)  :: kh     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out)  :: mfact

!cx Better To Do: first touch the following arrays
!fj    allocate( cinterp_TN(AI:AJ,1:3,ADM_gall,ADM_lall) )
!fj    allocate( cinterp_HN(AI:AJ,1:3,ADM_gall,ADM_lall) )
!fj    allocate( cinterp_TRA(TI:TJ,ADM_gall,ADM_lall) )
!fj    allocate( cinterp_PRA(      ADM_gall,ADM_lall) )
!fj>
    allocate( cinterp_TN(ADM_gall,ADM_lall,AI:AJ,1:3) )
    allocate( cinterp_HN(ADM_gall,ADM_lall,AI:AJ,1:3) )
    allocate( cinterp_TRA(ADM_gall,ADM_lall,TI:TJ) )
    allocate( cinterp_PRA(      ADM_gall,ADM_lall) )
!fj<


    if (my_snapshot.eq.0) then
!fj        cinterp_TN(AI:AJ,1:3,1:ADM_gall,1:ADM_lall) = 1.0
!fj        cinterp_HN(AI:AJ,1:3,1:ADM_gall,1:ADM_lall) = 1.0
!fj        cinterp_TRA(TI:TJ,1:ADM_gall,1:ADM_lall) = 1.0
!fj        cinterp_PRA(      1:ADM_gall,1:ADM_lall) = 1.0
!fj>
        cinterp_TN(1:ADM_gall,1:ADM_lall,AI:AJ,1:3) = 1.0_RP
        cinterp_HN(1:ADM_gall,1:ADM_lall,AI:AJ,1:3) = 1.0_RP
        cinterp_TRA(1:ADM_gall,1:ADM_lall,TI:TJ) = 1.0_RP
        cinterp_PRA(      1:ADM_gall,1:ADM_lall) = 1.0_RP
!fj<
        scl    (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0_RP
        scl_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0_RP
        kh     (1:ADM_gall   ,1:ADM_kall,1:ADM_lall   ) = 1.0_RP
        kh_pl  (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0_RP
        mfact = 1.0_RP
        return
    endif

    rewind(my_snapshot)
    do i=1,14
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


  subroutine OPRT_diffusion( &
       dscl, dscl_pl, &
       scl,  scl_pl,  &
       kh,   kh_pl,   &
       mfact          )
    use mod_gmtr, only: &
       GMTR_P_var_pl, &
       GMTR_T_var_pl, &
       GMTR_A_var_pl
    use mod_debug, only: DEBUG_rapstart, DEBUG_rapend
    implicit none

    real(RP), intent(out) :: dscl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: kh     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: mfact

!fj    real(RP)  :: vxt    (ADM_gall   ,ADM_kall,TI:TJ)
!fj    real(RP)  :: vyt    (ADM_gall   ,ADM_kall,TI:TJ)
!fj    real(RP)  :: vzt    (ADM_gall   ,ADM_kall,TI:TJ)
!fj    real(RP)  :: flux   (ADM_gall   ,ADM_kall,AI:AJ)
    real(RP)  :: vxt    (ADM_gall   ,TI:TJ)
    real(RP)  :: vyt    (ADM_gall   ,TI:TJ)
    real(RP)  :: vzt    (ADM_gall   ,TI:TJ)
    real(RP)  :: flux   (ADM_gall   ,AI:AJ)
    real(RP)  :: vxt_pl (ADM_gall_pl,ADM_kall)
    real(RP)  :: vyt_pl (ADM_gall_pl,ADM_kall)
    real(RP)  :: vzt_pl (ADM_gall_pl,ADM_kall)
    real(RP)  :: flux_pl(ADM_gall_pl,ADM_kall)

    real(RP) :: u1, u2, u3, smean

    integer :: nstart1, nend
    integer :: nstart2, nstart3
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, k, l, v

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !fj add >
    OPRT_nstart = suf(ADM_gmin,ADM_gmin)
    OPRT_nend   = suf(ADM_gmax,ADM_gmax)
    !fj add<
    !---------------------------------------------------------------------------
    !cx  write(0,'(a,10i8)') "<OPRT_diffusion>"

    call DEBUG_rapstart('OPRT_diffusion')

       do l = 1, ADM_lall
       nstart1 = suf(ADM_gmin-1,ADM_gmin-1)
       nstart2 = suf(ADM_gmin-1,ADM_gmin)
       nstart3 = suf(ADM_gmin,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )

!$omp  parallel default(none),private(k,ij,ip1j,ip1jp1,smean,u1,u2,&
!$omp  u3,ijp1,im1j,ijm1,im1jm1,n),&
!$omp  shared(nstart1,nstart2,nstart3,nend,scl,ADM_gall_1d,vxt,&
!$omp  vyt,vzt,cinterp_TN,cinterp_TRA,kh,l,cinterp_HN,ADM_gmin, &
!$omp  ADM_have_sgp,flux,cinterp_PRA,mfact,OPRT_nstart,dscl,OPRT_nend,ADM_kall)
       do k = 1, ADM_kall
!$omp do
       do n = nstart1, nend
          ij     = n
          ip1j   = n + 1
          ip1jp1 = n + 1 + ADM_gall_1d

          smean = ( scl(ij,k,l) + scl(ip1j,k,l) + scl(ip1jp1,k,l) ) / 3._RP

          u1 = 0.5_RP * (scl(ij    ,k,l)+scl(ip1j  ,k,l)) - smean
          u2 = 0.5_RP * (scl(ip1j  ,k,l)+scl(ip1jp1,k,l)) - smean
          u3 = 0.5_RP * (scl(ip1jp1,k,l)+scl(ij    ,k,l)) - smean

!fj          vxt(n,k,TI) = ( - u1 * cinterp_TN(AI ,1,ij  ,l) &
!fj                          - u2 * cinterp_TN(AJ ,1,ip1j,l) &
!fj                          + u3 * cinterp_TN(AIJ,1,ij  ,l) ) * cinterp_TRA(TI,ij,l)
!fj          vyt(n,k,TI) = ( - u1 * cinterp_TN(AI ,2,ij  ,l) &
!fj                          - u2 * cinterp_TN(AJ ,2,ip1j,l) &
!fj                          + u3 * cinterp_TN(AIJ,2,ij  ,l) ) * cinterp_TRA(TI,ij,l)
!fj          vzt(n,k,TI) = ( - u1 * cinterp_TN(AI ,3,ij  ,l) &
!fj                          - u2 * cinterp_TN(AJ ,3,ip1j,l) &
!fj                          + u3 * cinterp_TN(AIJ,3,ij  ,l) ) * cinterp_TRA(TI,ij,l)
!fj>
          vxt(n,TI) = ( - u1 * cinterp_TN(ij  ,l,AI ,1) &
                        - u2 * cinterp_TN(ip1j,l,AJ ,1) &
                        + u3 * cinterp_TN(ij  ,l,AIJ,1) ) * cinterp_TRA(ij,l,TI)
          vyt(n,TI) = ( - u1 * cinterp_TN(ij  ,l,AI ,2) &                      
                        - u2 * cinterp_TN(ip1j,l,AJ ,2) &                      
                        + u3 * cinterp_TN(ij  ,l,AIJ,2) ) * cinterp_TRA(ij,l,TI)
          vzt(n,TI) = ( - u1 * cinterp_TN(ij  ,l,AI ,3) &                      
                          - u2 * cinterp_TN(ip1j,l,AJ ,3) &                      
                          + u3 * cinterp_TN(ij  ,l,AIJ,3) ) * cinterp_TRA(ij,l,TI)
!fj<
       enddo
!$omp enddo nowait
!fj       enddo

!fj       do k = 1, ADM_kall
!$omp do
       do n = nstart1, nend
          ij     = n
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d

          smean = ( scl(ij,k,l) + scl(ip1jp1,k,l) + scl(ijp1,k,l) ) / 3.0_RP

          u1 = 0.5_RP * (scl(ij    ,k,l)+scl(ip1jp1,k,l)) - smean
          u2 = 0.5_RP * (scl(ip1jp1,k,l)+scl(ijp1  ,k,l)) - smean
          u3 = 0.5_RP * (scl(ijp1  ,k,l)+scl(ij    ,k,l)) - smean

          vxt(n,TJ) = ( - u1 * cinterp_TN(ij  ,l,AIJ,1) &
                        + u2 * cinterp_TN(ijp1,l,AI ,1) &
                        + u3 * cinterp_TN(ij  ,l,AJ ,1) ) * cinterp_TRA(ij,l,TJ)
          vyt(n,TJ) = ( - u1 * cinterp_TN(ij  ,l,AIJ,2) &                      
                        + u2 * cinterp_TN(ijp1,l,AI ,2) &                      
                        + u3 * cinterp_TN(ij  ,l,AJ ,2) ) * cinterp_TRA(ij,l,TJ)
          vzt(n,TJ) = ( - u1 * cinterp_TN(ij  ,l,AIJ,3) &                      
                          + u2 * cinterp_TN(ijp1,l,AI ,3) &                      
                          + u3 * cinterp_TN(ij  ,l,AJ ,3) ) * cinterp_TRA(ij,l,TJ)
       enddo
!fj       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
!fj          do k = 1, ADM_kall
!$omp master
             vxt(suf(ADM_gmin-1,ADM_gmin-1),TI) = vxt(suf(ADM_gmin,ADM_gmin-1),TJ)
             vyt(suf(ADM_gmin-1,ADM_gmin-1),TI) = vyt(suf(ADM_gmin,ADM_gmin-1),TJ)
             vzt(suf(ADM_gmin-1,ADM_gmin-1),TI) = vzt(suf(ADM_gmin,ADM_gmin-1),TJ)
!fj          enddo
!$omp end master
!$omp barrier
       endif

!fj       nstart = suf(ADM_gmin-1,ADM_gmin  )

!fj       do k = 1, ADM_kall
!$omp do
       do n = nstart2, nend
          ij     = n
          ip1j   = n + 1
          ijp1   = n     + ADM_gall_1d
          ip1jp1 = n + 1 + ADM_gall_1d
          im1j   = n - 1
          ijm1   = n     - ADM_gall_1d

          flux(n,AI ) = 0.25_RP *   ( (vxt(ijm1,TJ)+vxt(ij  ,TI)) * cinterp_HN(ij,l,AI,1) &
                                   + (vyt(ijm1,TJ)+vyt(ij  ,TI)) * cinterp_HN(ij,l,AI,2) &
                                   + (vzt(ijm1,TJ)+vzt(ij  ,TI)) * cinterp_HN(ij,l,AI,3) &
                                   ) * (kh(ij,k,l)+kh(ip1j  ,k,l))
       enddo
!$omp enddo nowait
!fj       enddo

!fj       nstart3 = suf(ADM_gmin-1,ADM_gmin-1)
!fj       nend   = suf(ADM_gmax  ,ADM_gmax  )

 !fj      do k = 1, ADM_kall
!$omp do
       do n = nstart1, nend
          ij     = n
          ip1jp1 = n + 1 + ADM_gall_1d

          flux(n,AIJ) =   0.25_RP * ( (vxt(ij  ,TI)+vxt(ij  ,TJ)) * cinterp_HN(ij,l,AIJ,1) &
                                   + (vyt(ij  ,TI)+vyt(ij  ,TJ)) * cinterp_HN(ij,l,AIJ,2) &
                                   + (vzt(ij  ,TI)+vzt(ij  ,TJ)) * cinterp_HN(ij,l,AIJ,3) &
                                   ) * (kh(ij,k,l)+kh(ip1jp1,k,l))
       enddo
!$omp enddo nowait
!fj       enddo

!fj       nstart = suf(ADM_gmin  ,ADM_gmin-1)
!fj       nend   = suf(ADM_gmax  ,ADM_gmax  )

!fj       do k = 1, ADM_kall
!$omp do
       do n = nstart3, nend
          ij     = n
          ijp1   = n     + ADM_gall_1d
          im1j   = n - 1

          flux(n,AJ ) =   0.25_RP * ( (vxt(ij  ,TJ)+vxt(im1j,TI)) * cinterp_HN(ij,l,AJ ,1) &
                                   + (vyt(ij  ,TJ)+vyt(im1j,TI)) * cinterp_HN(ij,l,AJ ,2) &
                                   + (vzt(ij  ,TJ)+vzt(im1j,TI)) * cinterp_HN(ij,l,AJ ,3) &
                                   ) * (kh(ij,k,l)+kh(ijp1  ,k,l))
       enddo
!fj       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
!fj          do k = 1, ADM_kall
!$omp master
             flux(suf(ADM_gmin,ADM_gmin-1),AJ) = 0._RP
!$omp end master 
!$omp barrier
!fj          enddo
       endif

!fj       do k = 1, ADM_kall
!$omp do 
       do n = OPRT_nstart, OPRT_nend
          ij     = n
          im1j   = n - 1
          im1jm1 = n - 1 - ADM_gall_1d
          ijm1   = n     - ADM_gall_1d

          dscl(n,k,l) = ( flux(ij,AI ) - flux(im1j  ,AI ) &
                        + flux(ij,AIJ) - flux(im1jm1,AIJ) &
                        + flux(ij,AJ ) - flux(ijm1  ,AJ ) ) * cinterp_PRA(ij,l) * mfact
       enddo
!$omp enddo nowait
       enddo
!$omp end parallel

    enddo
!    call stop_colleciton('OPRT_diffusion_loop1')

!cx cut off other loops

    call DEBUG_rapend('OPRT_diffusion')

    return
  end subroutine OPRT_diffusion
end module mod_oprt
!-------------------------------------------------------------------------------


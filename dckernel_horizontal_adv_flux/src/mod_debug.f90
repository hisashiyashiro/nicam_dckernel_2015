module mod_debug
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_adm, only: &
     ADM_LOG_FID, &
     ADM_NSYS,    &
     ADM_MAXFNAME
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: DEBUG_rapstart
  public :: DEBUG_rapend
  public :: DEBUG_rapreport
  !
  !++ Private procedure
  !
  private :: DEBUG_rapid

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                 private, parameter :: DEBUG_rapnlimit = 100
  integer,                 private,      save :: DEBUG_rapnmax   = 0
  character(len=ADM_NSYS), private,      save :: DEBUG_rapname(DEBUG_rapnlimit)
  real(8),                 private,      save :: DEBUG_raptstr(DEBUG_rapnlimit)
  real(8),                 private,      save :: DEBUG_rapttot(DEBUG_rapnlimit)
  integer,                 private,      save :: DEBUG_rapnstr(DEBUG_rapnlimit)
  integer,                 private,      save :: DEBUG_rapnend(DEBUG_rapnlimit)

contains


  !-----------------------------------------------------------------------------
  function DEBUG_rapid( rapname ) result(id)
    implicit none

    character(len=*), intent(in) :: rapname

    integer :: id
    !---------------------------------------------------------------------------

    if ( DEBUG_rapnmax >= 1 ) then
       do id = 1, DEBUG_rapnmax
          if( trim(rapname) == trim(DEBUG_rapname(id)) ) return
       enddo
    endif

    DEBUG_rapnmax     = DEBUG_rapnmax + 1
    id                = DEBUG_rapnmax
    DEBUG_rapname(id) = trim(rapname)
    DEBUG_raptstr(id) = 0.D0
    DEBUG_rapttot(id) = 0.D0
    DEBUG_rapnstr(id) = 0
    DEBUG_rapnend(id) = 0

  end function DEBUG_rapid

  !-----------------------------------------------------------------------------
  subroutine DEBUG_rapstart( rapname )
    use mod_adm, only: &
       ADM_MPItime
    implicit none

    character(len=*), intent(in) :: rapname

    integer :: id
    !---------------------------------------------------------------------------

    id = DEBUG_rapid( rapname )

    DEBUG_raptstr(id) = ADM_MPItime()
    DEBUG_rapnstr(id) = DEBUG_rapnstr(id) + 1

    !write(ADM_LOG_FID,*) rapname, DEBUG_rapnstr(id)

#ifdef _FAPP_
    call fapp_start( rapname, id, 1 )
    !call START_COLLECTION( rapname )
#endif

    return
  end subroutine DEBUG_rapstart

  !-----------------------------------------------------------------------------
  subroutine DEBUG_rapend( rapname )
    use mod_adm, only: &
       ADM_MPItime
    implicit none

    character(len=*), intent(in) :: rapname

    integer :: id
    !---------------------------------------------------------------------------

    id = DEBUG_rapid( rapname )

    DEBUG_rapttot(id) = DEBUG_rapttot(id) + ( ADM_MPItime()-DEBUG_raptstr(id) )
    DEBUG_rapnend(id) = DEBUG_rapnend(id) + 1

#ifdef _FAPP_
    call fapp_stop( rapname, id, 1 )
    !call STOP_COLLECTION( rapname )
#endif

    return
  end subroutine DEBUG_rapend

  !-----------------------------------------------------------------------------
  subroutine DEBUG_rapreport
    implicit none

    integer :: id
    !---------------------------------------------------------------------------

    if ( DEBUG_rapnmax >= 1 ) then

       do id = 1, DEBUG_rapnmax
          if ( DEBUG_rapnstr(id) /= DEBUG_rapnend(id) ) then
              write(*,*) '*** Mismatch Report',id,DEBUG_rapname(id),DEBUG_rapnstr(id),DEBUG_rapnend(id)
          endif
       enddo

       write(*,*)
       write(*,*) '*** Computational Time Report'

       do id = 1, DEBUG_rapnmax
          write(*,'(1x,A,I3.3,A,A,A,F10.3,A,I7)') &
          '*** ID=',id,' : ',DEBUG_rapname(id),' T=',DEBUG_rapttot(id),' N=',DEBUG_rapnstr(id)
       enddo

    else
       write(*,*)
       write(*,*) '*** Computational Time Report: NO item.'
    endif

    return
  end subroutine DEBUG_rapreport

end module mod_debug


module mod_snapshot_for_kernel

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
    ADM_kmin,    &
    ADM_kmax,    &
    TI  => ADM_TI,  &
    TJ  => ADM_TJ,  &
    AI  => ADM_AI,  &
    AIJ => ADM_AIJ, &
    AJ  => ADM_AJ,  &
    K0  => ADM_KNONE

  implicit none
  public :: snapshot_seq_bin_open
  public :: snapshot_seq_bin_close
  public :: ADM_snap_read
  public :: ADM_ascii_write
  integer, public, save :: my_snapshot
  integer, private, save :: my_id_saved
contains


! Open the snapshot file.
! The first argument is used to identify the snapshot file name
! If the first argument is NULL "", then snapshot file will not be read
! The second argument should correspond to ADM_prc_me, i.e. 1, 2, ..
  subroutine snapshot_seq_bin_open (prefix, my_id)
    implicit none
    !cx use mod_adm, only: my_id => ADM_prc_me
    character*(*), intent(in) :: prefix
    integer, intent(in) :: my_id
    integer :: i
    character*(1) :: string_null=""
    character(LEN=80) :: fname
    logical :: i_opened

    if(trim(prefix).eq.string_null) then
        my_snapshot = 0
        my_id_saved = my_id
        return
    endif

    write(0,'(a,10i8)') "<snapshot_seq_bin_open> my_id=", my_id
    do i=7,99
       inquire (i, opened=i_opened)
       if(.not.i_opened) exit
    enddo
    my_snapshot = i ! if i==99, the inquire logic is somewhat questionable
    my_id_saved = my_id
    write(fname,'(a,i3.3)') "snapshot."//trim(prefix)//".", my_id
    open( unit=my_snapshot, file=trim(fname), form='unformatted')

    return
  end subroutine snapshot_seq_bin_open


  subroutine snapshot_seq_bin_close
    implicit none
    write(0,'(a,10i8)') "<snapshot_seq_bin_close> my_id_saved=", my_id_saved
    if (my_snapshot.ne.0) then
    close( unit=my_snapshot)
    endif
    return
  end subroutine snapshot_seq_bin_close

  subroutine ADM_snap_read
    implicit none

    if (my_snapshot.eq.0) then
       write(0,'(a,10i8)') "<ADM_snap_read> cold start"
       ! parameters :: TI=1, TJ=2, AI=1, AIJ=2, AJ=3, K0=1
       ADM_have_pl=.true.
       ADM_lall=1
       ADM_lall_pl=2
       ADM_gall=16900
       ADM_gall_pl=6
       ADM_kall= 96
       ADM_gall_1d=130
       ADM_gmin=2
       ADM_gmax=129
       ADM_gmax_pl=6
       ADM_kmin=2
       ADM_kmax=95
       allocate( ADM_have_sgp(ADM_lall) )
       ADM_have_sgp(1:ADM_lall)=.true.
    else

        write(0,'(a,10i8)') "<ADM_snap_read> restart from snapshot file"
        rewind(my_snapshot)
        read (my_snapshot) ! parameters :: TI=1, TJ=2, AI=1, AIJ=2, AJ=3, K0=1
        read (my_snapshot) ADM_have_pl
        read (my_snapshot) ADM_lall
        read (my_snapshot) ADM_lall_pl
        read (my_snapshot) ADM_gall
        read (my_snapshot) ADM_gall_pl
        read (my_snapshot) ADM_kall
        read (my_snapshot) ADM_gall_1d
        read (my_snapshot) ADM_gmin
        read (my_snapshot) ADM_gmax
        read (my_snapshot) ! ADM_gslf_pl
        read (my_snapshot) ! ADM_gmin_pl
        read (my_snapshot) ADM_gmax_pl
        read (my_snapshot) ADM_kmin
        read (my_snapshot) ADM_kmax
        allocate( ADM_have_sgp(ADM_lall) )
        read (my_snapshot) ADM_have_sgp

    end if

    return
  end subroutine ADM_snap_read


  subroutine ADM_ascii_write
    implicit none
    integer :: suf,i,j, l, nstart, nend
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i) ! statement function !
    character(LEN=80) :: fname
    integer :: my_ascii=77

    write(fname,'(a,i3.3)') "snapshot.ADM_ascii.", my_id_saved
    open( unit=my_ascii, file=trim(fname), form='formatted')

    write(my_ascii,'(a,10i8)') "<ADM_ascii_write> for my_id:", my_id_saved
    write(my_ascii,'(a,10i8)') "TI,TJ,AI,AIJ,AJ,K0:", TI,TJ,AI,AIJ,AJ,K0
    write(my_ascii,'(a,L2)')   "ADM_have_pl:", ADM_have_pl
    write(my_ascii,'(a,10i8)') "ADM_lall   :", ADM_lall
    write(my_ascii,'(a,10i8)') "ADM_lall_pl:", ADM_lall_pl
    write(my_ascii,'(a,10i8)') "ADM_gall   :", ADM_gall
    write(my_ascii,'(a,10i8)') "ADM_gall_pl:", ADM_gall_pl
    write(my_ascii,'(a,10i8)') "ADM_kall:", ADM_kall
    write(my_ascii,'(a,10i8)') "ADM_gall_1d:", ADM_gall_1d
    write(my_ascii,'(a,10i8)') "ADM_gslf_pl:", ADM_gslf_pl
    write(my_ascii,'(a,10i8)') "ADM_gmin,    ADM_gmax:   ", ADM_gmin, ADM_gmax
    write(my_ascii,'(a,10i8)') "ADM_gmin_pl, ADM_gmax_pl:", ADM_gmin_pl, ADM_gmax_pl
    write(my_ascii,'(a,10i8)') "ADM_kmin,    ADM_kmax:   ", ADM_kmin, ADM_kmax
    write(my_ascii,'(a,(2x,10L2/))') "ADM_have_sgp(1:*):", (ADM_have_sgp(l),l=1, ADM_lall)
    do l = 1, ADM_lall
       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax  ,ADM_gmax  )
       write(my_ascii,'(a,10i8)') "(1:ADM_lall), nstart, nend=", nstart, nend
    enddo
    close(my_ascii)
    return
  end subroutine ADM_ascii_write

end module mod_snapshot_for_kernel


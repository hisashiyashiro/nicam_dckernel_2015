module mod_cnst
  use mod_precision
  implicit none
  private
  !++ Public procedure
  public :: CNST_setup

  !++ Public parameters & variables
  real(RP), public, save :: CNST_ERADIUS = 6.37122E6_RP ! Radius of the Earth [m]
  real(RP), public, save :: CNST_EOHM    = 7.292E-5_RP   ! Angular velocity of the Earth [/s]
  real(RP), public, save :: CNST_EGRAV   = 9.80616_RP  ! Gravitational accerlaration of the Earth [m/s2]

  real(RP), public, save :: CNST_RAIR    =  287.0_RP   ! Gas constant of air
  real(RP), public, save :: CNST_RVAP    =  461.5_RP   ! Gas constant of vapor

  real(RP), public, save :: CNST_CP      = 1004.5_RP   ! Specific heat of air (consant pressure)
  real(RP), public, save :: CNST_CV                   ! Specific heat of air (consant volume)

  real(RP), public, save :: CNST_CPV     = 1846.0_RP   ! Specific heat of vapor (consant pressure)
  real(RP), public, save :: CNST_CVV                  ! Specific heat of vapor (consant volume)
  real(RP), public, save :: CNST_CL      = 4218.0_RP   ! Specific heat of water
  real(RP), public, save :: CNST_CI      = 2006.0_RP   ! Specific heat of ice

  !------ cp/cv
  real(RP), public, save :: CNST_GAMMA
  !<----- calculated in sub[CNST_setup].
  !
  !------ R/cp
  real(RP), public, save :: CNST_KAPPA
  !<----- calculated in sub[CNST_setup].
  !
  !------ dry lapse rate [K/m]
  real(RP), public, save :: CNST_LAPS
  !<----- calculated in sub[CNST_setup].
  !
  !------ molecular weight ( water/air )
  real(RP), public, save :: CNST_EPSV
  !<----- calculated in sub[CNST_setup].
  !
  !------ 1/epsv-1
  real(RP), public, save :: CNST_EPSVT
  !<----- calculated in sub[CNST_setup].
  !
  !------ Density of water
  real(RP), public, save :: CNST_DWATR = 1000.0_RP
  !
  !------ Saturate pressure of water vapor at 0C
  real(RP), public, save :: CNST_PSAT0 = 610.7_RP
  !<----- unit : [Pa]
  !
  !------ Latent heat of vaporizaion at 0C
!  real(RP), public, save :: CNST_LH0   = 2.5008D+6 [mod] 20120704 H.Yashiro
  real(RP), public, save :: CNST_LH0   = 2.501E+6_RP
  !
  !------ Latent heat of vaporizaion at 0K
  real(RP), public, save :: CNST_LH00
  !<----- calculated in sub[CNST_setup].
  !
  !------ Latent heat of sublimation at 0C
!  real(RP), public, save :: CNST_LHS0  = 2.8342E+6 [mod] 20120704 H.Yashiro
  real(RP), public, save :: CNST_LHS0  = 2.834E+6_RP
  !
  !------ Latent heat of sublimation at 0K
  real(RP), public, save :: CNST_LHS00
  !<----- calculated in sub[CNST_setup].
  !
  !------ Latent heat of fusion at 0C
  real(RP), public, save :: CNST_LHF0
  !<----- calculated in sub[CNST_setup].
  !
  !------ latent heat of fusion at 0K
  real(RP), public, save :: CNST_LHF00
  !<----- calculated in sub[CNST_setup].
  !
  !------ Latent heat of melting
  real(RP), public, save :: CNST_EMELT = 3.40E+5_RP
  !
  !------ Melting temperature of water
  real(RP), public, save :: CNST_TMELT = 273.15_RP
  !
  !------ Freeze point of sea
  real(RP), public, save :: CNST_TFRZS  = 271.35_RP
  !
  !------ Wet-bulb temp. rain/snow
  real(RP), public, save :: CNST_TQICE = 273.15_RP
  !
  !------ Stefan-Boltzman constant
  real(RP), public, save :: CNST_STB   = 5.67E-8_RP
  !
  !------ Karman constant
  real(RP), public, save :: CNST_KARMAN = 0.4_RP
  !
  !------ Surface pressure
  real(RP), public, save :: CNST_PRES0    = 101325.0_RP
  !
  !------ Surface temperature
  real(RP), public, save :: CNST_TEMS0    = 300.0_RP
  !
  !------ Standard pressure
  real(RP), public, save :: CNST_PRE00    = 1.0E+5_RP
  !
  !------ Standard temperature
  real(RP), public, save :: CNST_TEM00    = 273.15_RP
  !
  !------ Standard density
  real(RP), public, save :: CNST_RHO00
  !<----- calculated in sub[CNST_setup].
  !
  !====== Misc. constants ======
  !
  !------ Definition of PI
  real(RP), public, save :: CNST_PI = 3.14159265358979323846_RP

  real(RP), public, save :: CNST_D2R

  !------ Allowable minimum value
  real(RP), public :: CNST_EPS_ZERO = 1.E-16_RP
  !
  !------ Allowable maximum value
  real(RP), public, parameter :: CNST_MAX_REAL = 1.E+30_RP
  !
  !------ Missing value
  real(RP), public, parameter :: CNST_VMISS    = 0.0_RP
  !
  !------ Undefined value
  real(DP), public, parameter :: CNST_UNDEF    = -99.9D+33
  !
  !------ Undefined value
  real(4), public, parameter :: CNST_UNDEF4   = -99.9E+33
  !
  !------ Undefined value
  integer(4), public, parameter :: CNST_UNDEF2   = -32768

  !-----------------------------------------------------------------------------

contains

  subroutine CNST_setup
    implicit none
    real(RP) :: earth_radius               ! Earth radius
    real(RP) :: earth_angvel               ! Anguler velocity of the earth
    real(RP) :: small_planet_factor = 1.0_RP ! small planet factor
    real(RP) :: earth_gravity              ! Gravitational accelaration
    real(RP) :: gas_cnst                   ! Gas constant of dry air
    real(RP) :: gas_cnst_vap               ! Gas constant of water vapour
    real(RP) :: specific_heat_pre          ! Specific heat of air( const pre )
    real(RP) :: specific_heat_pre_vap      ! Specific heat of water vapour ( const pre )
    real(RP) :: latent_heat_vap            ! latent heat of vaporization LH0 ( 0 deg )
    real(RP) :: latent_heat_sub            ! latent heat of sublimation LHS0 ( 0 deg )

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- initialization of controled parameters
    earth_radius = CNST_ERADIUS
    earth_angvel = CNST_EOHM
    earth_gravity   = CNST_EGRAV
    gas_cnst    = CNST_RAIR
    gas_cnst_vap = CNST_RVAP
    specific_heat_pre      = CNST_CP
    specific_heat_pre_vap  = CNST_CPV
    latent_heat_vap = CNST_LH0
    latent_heat_sub = CNST_LHS0


    CNST_ERADIUS = earth_radius / small_planet_factor
    CNST_EOHM    = earth_angvel * small_planet_factor
    CNST_EGRAV   = earth_gravity
    CNST_RAIR    = gas_cnst
    CNST_RVAP    = gas_cnst_vap
    CNST_CP      = specific_heat_pre
    CNST_CPV     = specific_heat_pre_vap
    CNST_LH0     = latent_heat_vap
    CNST_LHS0    = latent_heat_sub

    !--- calculate other parameters
    CNST_PI    = 4.0_RP * atan( 1.0_RP )
    CNST_D2R   = CNST_PI / 180.0_RP

    CNST_CV    = CNST_CP - CNST_RAIR
    CNST_GAMMA = CNST_CP / CNST_CV
    CNST_KAPPA = CNST_RAIR / CNST_CP
    CNST_RHO00 = CNST_PRE00 / CNST_RAIR / CNST_TEM00
    CNST_LAPS  = CNST_EGRAV / CNST_CP

    CNST_CVV   = CNST_CPV - CNST_RVAP
    CNST_EPSV  = CNST_RAIR / CNST_RVAP
    CNST_EPSVT = 1.0_RP/CNST_EPSV - 1.0_RP

    CNST_LH00  = CNST_LH0  - ( CNST_CPV - CNST_CL ) * CNST_TEM00
    CNST_LHS00 = CNST_LHS0 - ( CNST_CPV - CNST_CI ) * CNST_TEM00
    CNST_LHF0  = CNST_LHS0 - CNST_LH0
    CNST_LHF00 = CNST_LHF0 - ( CNST_CL - CNST_CI ) * CNST_TEM00 ! bugfix: CNST_LHS0 -> CNST_LHF0

    return
  end subroutine CNST_setup

end module mod_cnst
!-------------------------------------------------------------------------------

#include "rundeck_opts.h"
module CHEM_DRV
  !================================================================================================
  ! Module CHEM_DRV is a module that enables ModelE to drive the GEOS-Chem
  ! chemistry-transport model. Initial version Jul 12, 2020.
  ! 
  ! Author: Lee T. Murray (lee.murray@rochester.edu)
  !===============================================================================================
  USE QUSDEF,      ONLY : nmom
  USE RESOLUTION,  ONLY : im, jm, lm
  USE ERRCODE_MOD, ONLY : GC_SUCCESS 
  USE CHEM_COM,    ONLY : IsAdvected, NTM, TrM, TrMom, TrName

  USE Input_Opt_Mod,           ONLY: OptInput
  USE State_Chm_Mod,           ONLY: ChmState
  USE State_Grid_Mod,          ONLY: GrdState
  USE State_Met_Mod,           ONLY: MetState
  USE State_Diag_Mod,          ONLY: DgnState
  USE DiagList_Mod,            ONLY: DgnList, Init_DiagList, Print_DiagList
  USE TaggedDiagList_Mod,      ONLY: TaggedDgnList, Init_TaggedDiagList, Print_TaggedDiagList
  USE HCO_Types_Mod,           ONLY: ConfigObj
  USE Precision_Mod,           ONLY: f8, fp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: INIT_CHEM
  PUBLIC :: DO_CHEM
  PUBLIC :: IO_CHEM
  PUBLIC :: TrDYNAM
  PUBLIC :: accumGCsubdd
  PUBLIC :: NYMDb
  PUBLIC :: NHMSb
  PUBLIC :: NYMDe
  PUBLIC :: NHMSe

  SAVE

  TYPE(OptInput)             :: Input_Opt  ! Input Options (same for all domains)
  TYPE(MetState)             :: State_Met  ! Meteorology state
  TYPE(ChmState)             :: State_Chm  ! Chemistry state
  TYPE(DgnState)             :: State_Diag ! Diagnostics state
  TYPE(DgnList)              :: Diag_List  ! Diagnostics state
  TYPE(TaggedDgnList)        :: TaggedDiag_List ! Diagnostics state
  TYPE(GrdState)             :: State_Grid ! Grid state


  ! Start, stop and size of main grid
  INTEGER                                     :: J_1, J_0, I_1, I_0, J_0H, J_1H, NI, NJ    
  ! Start and end date and time
  INTEGER   :: NYMDb,   NHMSb,    NYMDe,   NHMSe

  ! Default flags for GEOS-Chem operators to use; may be overwritten by rundeck
  LOGICAL                               :: DoGCConv     = .true.
  LOGICAL                               :: DoGCEmis     = .true.
  LOGICAL                               :: DoGCTend     = .true.
  LOGICAL                               :: DoGCTurb     = .true.
  LOGICAL                               :: DoGCChem     = .true.
  LOGICAL                               :: DoGCDryDep   = .true.
  LOGICAL                               :: DoGCWetDep   = .true.

  LOGICAL                               :: first_chem = .true.

  !-----------------------------------------------------------------
  ! 40-level GISS grid
  !-----------------------------------------------------------------

  ! Ap [hPa] for 40 levels (41 edges)
  REAL(fp), PARAMETER :: AP(41) = (/                  &
        0.000000,   3.597122,   7.553957,  12.050360, &
       16.906475,  22.302158,  28.597122,  35.791367, &
       43.884892,  52.517986,  61.510791,  70.683453, &
       80.035971,  89.028777,  97.661871, 105.755396, &
      113.309353, 120.143885, 126.258993, 131.834532, &
      136.870504, 141.546763, 145.863309, 150.000000, &
      128.000000, 108.000000,  90.000000,  73.000000, &
       57.000000,  43.000000,  31.000000,  20.000000, &
       10.000000,   5.620000,   3.160000,   1.780000, &
        1.000000,   0.562000,   0.316000,   0.178000, &
        0.100000                                       /)

  ! Bp [unitless] for 40 levels (41 edges)
  REAL(fp), PARAMETER :: BP(41) = (/                   &
       1.00000000, 0.97601918, 0.94964029, 0.91966427, &
       0.88729017, 0.85131894, 0.80935252, 0.76139089, &
       0.70743405, 0.64988010, 0.58992806, 0.52877698, &
       0.46642686, 0.40647482, 0.34892086, 0.29496403, &
       0.24460432, 0.19904077, 0.15827338, 0.12110312, &
       0.08752998, 0.05635492, 0.02757794, 0.00000000, &
       0.00000000, 0.00000000, 0.00000000, 0.00000000, &
       0.00000000, 0.00000000, 0.00000000, 0.00000000, &
       0.00000000, 0.00000000, 0.00000000, 0.00000000, &
       0.00000000, 0.00000000, 0.00000000, 0.00000000, &
       0.00000000                                       /)

CONTAINS

  !==========================================================================================================

  SUBROUTINE DO_CHEM

    ! ModelE modules
    USE DOMAIN_DECOMP_ATM, ONLY : AM_I_ROOT, GRID, getDomainBounds, hasnorthpole, hassouthpole
    USE DOMAIN_DECOMP_1D,  ONLY : HALO_UPDATE, SOUTH, NORTH
    USE DYNAMICS,          ONLY : SIGE
    USE MODEL_COM,         ONLY : modelEclock, itime, ItimeI, DTsrc
    USE ATM_COM,           ONLY : gz, mma, mws, pedn, pk, pmid, ptropo, q, qci, qcl, t, ualij, &
                                  valij, zatmo
#ifdef CALC_MERRA2_LIKE_DIAGS
    USE CLOUDS_COM,        ONLY : tauss, taumc, cldmc, cldss, cldss3d, pficu, pflcu, pfilsan, pfllsan
    USE CLOUDS_COM,        ONLY : dtrain, dqrcu, dqrlsan, reevapcn, reevapls, cmfmc
#endif
    USE FLUXES,            ONLY : atmsrf, atmlnd, prec, precss, focean, fland, flice
    USE GEOM,              ONLY : axyp , byaxyp
    USE GHY_COM,           ONLY : fearth, wearth, aiearth, wfcs
#ifdef CALC_MERRA2_LIKE_DIAGS
    USE GHY_COM,           ONLY : lai_save, z0m_save
#endif
    USE LAKES_COM,         ONLY : flake
#ifdef CALC_MERRA2_LIKE_DIAGS
    USE O3mod,             ONLY : save_to3
    USE RAD_COM,           ONLY : save_alb, taui3d, tauw3d
#endif
    USE RAD_COM,           ONLY : cfrac, srdn, fsrdir, srvissurf, cosz1, save_cosz2
    USE SEAICE_COM,        ONLY : si_atm, si_ocn
    USE CONSTANT,          ONLY : avog, bygasc, bygrav, lhe, tf, teeny
    USE LIGHTNING,         ONLY : cth_save, flash_dens

    ! GEOS-Chem modules
    USE HCO_Interface_Common, ONLY : SetHcoTime
    USE Time_Mod,          ONLY : Accept_External_Date_Time
    USE Emissions_Mod,     ONLY : Emissions_Run
    USE State_Chm_Mod,     ONLY : Ind_
    USE Calc_Met_Mod,      ONLY : AirQnt
    USE Pressure_Mod,      ONLY : Set_Floating_Pressures
    USE Pressure_Mod,      ONLY : Accept_External_Pedge
    USE Calc_Met_Mod,      ONLY : Set_Dry_Surface_Pressure
    USE Calc_Met_Mod,      ONLY : GCHP_Cap_Tropopause_Prs
    USE PBL_Mix_Mod,       ONLY : Compute_PBL_Height
    USE VDIFF_Mod,         ONLY : Max_PblHt_For_Vdiff
    USE ERROR_MOD,         ONLY : Safe_Div, IT_IS_NAN, ERROR_STOP
    USE UnitConv_Mod,      ONLY : Convert_Spc_Units, KG_SPECIES, &
                                  MOLES_SPECIES_PER_MOLES_DRY_AIR

    USE Photolysis_Mod,  ONLY : Init_Photolysis

    IMPLICIT NONE

    INTEGER   :: NYMD,    NHMS,     YEAR,    MONTH,    DAY
    INTEGER   :: DOY,     HOUR,     MINUTE,  SECOND
    INTEGER   :: I,       J,        L,       K,        N,       M
    INTEGER   :: II,      JJ,       III,     JJJ,      RC
    REAL*4    :: MINUTES, hElapsed, UTC
    REAL*8    :: sElapsed
    LOGICAL   :: IsChemTime, IsRadTime

    CHARACTER(LEN=256)       :: ThisLoc
    CHARACTER(LEN=512)       :: ErrMsg

    ! External functions (from shared/Utilities.F90)
    REAL*8 SLP
    REAL*8 QSAT

    LOGICAL, SAVE :: FIRST_CHEM = .true.

    ! Assume initial success
    RC = GC_SUCCESS

    !====================================
    ! Get time information
    !====================================
    YEAR    = modelEclock%getYear()
    MONTH   = modelEclock%getMonth()
    DAY     = modelEclock%getDate()
    HOUR    = modelEclock%getHour()
    MINUTE  = DTsrc*ITIME ! This works because model must start at top of hour
    MINUTES = modulo( MINUTE, 3600 ) / 60d0
    MINUTE  = floor( MINUTES )
    SECOND  = 0

    NYMD    = year*10000 + month*100 + day
    NHMS    = hour*10000 + minute*100 + second

    UTC     = HOUR + MINUTE / 60.0

    hElapsed = (DTsrc*(ITIME-ItimeI)) / 3600d0 ! Hours elapsed
    sElapsed = (DTsrc*(ITIME-ItimeI))          ! Seconds elapsed

    ! Is it time for chemistry?
    IF ( ( sElapsed / DTsrc ) == FLOOR( sElapsed / DTsrc ) ) THEN
       IsChemTime = .TRUE.
    ELSE
       IsChemTime = .FALSE.
    ENDIF
    
    IsRadTime = .FALSE. ! Hardwire RRTMG off

    DO JJJ = J_0, J_1
       DO III = I_0, I_1

          ! GEOS-Chem local index
          II = III - I_0 + 1
          JJ = JJJ - J_0 + 1

          ! GISS meteorology index (GISS only has one polar box)
          I = III
          J = JJJ
          if(hassouthpole(grid) .and. JJJ .eq. J_0 ) I = 1
          if(hasnorthpole(grid) .and. JJJ .eq. J_1 ) I = 1

          !----------------------------------------------------------------------
          ! Surface fields
          !----------------------------------------------------------------------

#ifdef CALC_MERRA2_LIKE_DIAGS
          ! Visible surface albedo [1]
          State_Met%ALBD        (II,JJ) = save_alb(i,j)
#endif

          ! Grid box surface area [cm2]
          State_Met%AREA_M2     (II,JJ) = axyp(i,j)

          ! Chemistry grid level [1]
          State_Met%ChemGridLev (II,JJ) = LM

          ! Column cloud fraction [1]
          State_Met%CLDFRC      (II,JJ) = cfrac(i,j)

#ifdef CALC_MERRA2_LIKE_DIAGS
          ! Max cloud top height [levels]
          State_Met%CLDTOPS(II,JJ) = 1
          DO K = LM, 1, -1
             IF ( State_Met%CMFMC(II,JJ,K) > 0d0 ) THEN
                State_Met%CLDTOPS(II,JJ) = K + 1
                EXIT
             ENDIF
          ENDDO
#endif

#ifdef MODEL_GEOS
          ! Convective fraction [1] (only used by GEOS)
          State_Met%CNV_FRC     (II,JJ) = 0.0
#endif

          ! Convective cloud depth [m]
          State_Met%CONV_DEPTH  (II,JJ) = cth_save(i,j)

          ! Latent heat flux [W/m2]
          State_Met%EFLUX       (II,JJ) = -atmsrf%latht(i,j)/dtsrc

          ! Lightning flash density [km-2 s-1]
          State_Met%FLASH_DENS  (II,JJ) = flash_dens(i,j)

          ! Olson land fraction [1]
          State_Met%FRCLND      (II,JJ) = fland(i,j)

          ! Fraction of lake [1]
          State_Met%FRLAKE      (II,JJ) = flake(i,j)

          ! Fraction of land [1]
          State_Met%FRLAND      (II,JJ) = fland(i,j)

          ! Fraction of land ice [1]
          State_Met%FRLANDICE   (II,JJ) = flice(i,j)

          ! Fraction of ocean [1]
          State_Met%FROCEAN     (II,JJ) = focean(i,j)

          ! Sfc sea ice fraction [1]
          State_Met%FRSEAICE    (II,JJ) = si_atm%RSI(i,j)*focean(i,j)

          ! Surface snow fraction [1]
          State_Met%FRSNOW       (II,JJ) = 0.0
          if ( si_ocn%snowi(i,j) > 0. ) &
               State_Met%FRSNOW(II,JJ) = si_atm%rsi(i,j)*flake(i,j)
          if ( atmlnd%SNOWE(i,j) > 0. ) &
               State_Met%FRSNOW(II,JJ) = State_Met%FRSNOW(II,JJ) + atmlnd%snowfr(i,j)*fearth(i,j)
          State_Met%FRSNOW(II,JJ) = min( 1.0, State_Met%FRSNOW(II,JJ) )

          ! Root soil wetness [1]
          State_Met%GWETROOT    (II,JJ) = 0.0
          if ( fearth(i,j) .gt. 0 ) then
             State_Met%GWETROOT  (II,JJ) = (wearth(i,j)+aiearth(i,j))/(wfcs(i,j)+1e-20)
          else
             State_Met%GWETROOT  (II,JJ) = 1 ! Set to 1 over oceans to match MERRA-2
          end if

          ! Top soil moisture [1] (assume same as GWETROOT for now)
          State_Met%GWETTOP     (II,JJ) = 0.0
          if ( fearth(i,j) .gt. 0 ) then
             State_Met%GWETTOP  (II,JJ) = (wearth(i,j)+aiearth(i,j))/(wfcs(i,j)+1e-20)
          else
             State_Met%GWETTOP  (II,JJ) = 1 ! Set to 1 over oceans to match MERRA-2
          end if

          ! Sensible heat flux [W/m2]
          State_Met%HFLUX       (II,JJ) = -atmsrf%sensht(i,j)/dtsrc

#ifdef CALC_MERRA2_LIKE_DIAGS
          ! Leaf area index [m2/m2] (online)
          State_Met%LAI         (II,JJ) = lai_save(i,j)
#endif

          ! Direct photsynthetically active radiation [W/m2]
          State_Met%PARDR       (II,JJ) = 0.82*srvissurf(i,j)*(fsrdir(i,j))*cosz1(i,j)

          ! Diffuse photsynthetically active radiation [W/m2]
          State_Met%PARDF       (II,JJ) = 0.82*srvissurf(i,j)*(1d0-fsrdir(i,j))*cosz1(i,j)

          ! PBL height [m] PBL top layer [1]
          State_Met%PBLH        (II,JJ) = atmsrf%dblavg(i,j)

          ! Surface geopotential height [m]
          State_Met%PHIS        (II,JJ) = zatmo(i,j)

          ! Anvil previp @ ground [kg/m2/s] -> mm/d
          State_Met%PRECANV     (II,JJ) = 0.0

          ! Conv  precip @ ground [kg/m2/s] -> mm/d
          State_Met%PRECCON     (II,JJ) = 86400d0*max(0.,prec(i,j)-precss(i,j))/dtsrc

          ! LS precip @ ground [kg/m2/s] -> mm/d
          State_Met%PRECLSC     (II,JJ) = 86400d0*precss(i,j)/dtsrc

          ! Total precip @ ground [kg/m2/s] -> mm/d
          State_Met%PRECTOT     (II,JJ) = 86400d0*prec(i,j)/dtsrc

          ! Wet surface pressure at start of timestep [hPa]
          State_Met%PS1_WET     (II,JJ) = pedn(1,i,j)

          ! Wet surface pressure at end of timestep [hPa]
          State_Met%PS2_WET     (II,JJ) = pedn(1,i,j)

          ! Wet interpolated surface pressure [hPa]
          State_Met%PSC2_WET    (II,JJ) = pedn(1,i,j)

          ! Dry surface pressure at start of timestep [hPa]
          State_Met%PS1_DRY     (II,JJ) = pedn(1,i,j)

          ! Dry surface pressure at end of timestep [hPa]
          State_Met%PS2_DRY     (II,JJ) = pedn(1,i,j)

          ! Dry interpolated surface pressure [hPa]
          State_Met%PSC2_DRY    (II,JJ) = pedn(1,i,j)

          ! Sea ice coverage 00-10% to 90-100% (only used by Hg)
          State_Met%SEAICE00    (II,JJ) = 0.0
          State_Met%SEAICE10    (II,JJ) = 0.0
          State_Met%SEAICE20    (II,JJ) = 0.0
          State_Met%SEAICE30    (II,JJ) = 0.0
          State_Met%SEAICE40    (II,JJ) = 0.0
          State_Met%SEAICE50    (II,JJ) = 0.0
          State_Met%SEAICE60    (II,JJ) = 0.0
          State_Met%SEAICE70    (II,JJ) = 0.0
          State_Met%SEAICE80    (II,JJ) = 0.0
          State_Met%SEAICE90    (II,JJ) = 0.0

          ! Sea level pressure [hPa]
          State_Met%SLP         (II,JJ) = slp(pedn(1,i,j),atmsrf%tsavg(i,j),bygrav*zatmo(i,j))*100.

          ! Snow depth [m]
          State_Met%SNODP       (II,JJ) = atmsrf%SNOWDP(i,j) * ( 1d0 - flice(i,j) )

          ! Snow mass [kg/m2]
          State_Met%SNOMAS      (II,JJ) = atmsrf%SNOW(i,j)

          ! COS(solar zenith angle) at current time
          State_Met%SUNCOS      (II,JJ) = cosz1(i,j)

          ! COS(solar zenith angle) at midpoint of chem timestep
          State_Met%SUNCOSmid   (II,JJ) = save_cosz2(i,j)

          ! Incident radiation @ ground [W/m2]
          State_Met%SWGDN       (II,JJ) = srdn(i,j)*save_cosz2(i,j)

#ifdef CALC_MERRA2_LIKE_DIAGS
          ! Total overhead O3 column [DU]
          State_Met%TO3         (II,JJ) = save_to3(i,j)
#endif

          ! Tropopause level [1]
          State_Met%TropLev(II,JJ) = 1
          DO K = LM, 1, -1
             IF ( pedn(k,i,j) >= ptropo(i,j) ) THEN
                State_Met%TropLev(II,JJ) = K
                EXIT
             ENDIF
          ENDDO

          ! Tropopause height [km]
          State_Met%TropHt      (II,JJ) = 0
          DO K = 1, State_Met%TropLev(II,JJ)-1
             State_Met%TropHt(II,JJ) = State_Met%TropHt(II,JJ) + State_Met%BXHEIGHT(II,JJ,K) * 1d-3
          ENDDO
          State_Met%TropHt(II,JJ) = State_Met%TropHt(II,JJ) + &
               State_Met%BXHEIGHT(II,JJ,State_Met%TropLev(II,JJ)) * 0.5d-3

          ! Tropopause pressure [hPa]
          State_Met%TROPP       (II,JJ) = ptropo(i,j)

          ! Surface temperature [K]
          State_Met%TS          (II,JJ) = atmsrf%tsavg(i,j) - tf + 273.15

          ! Surface skin temperature [K]
          State_Met%TSKIN       (II,JJ) = atmsrf%gtempr(i,j)

          ! E/W wind speed @ 10m ht [m/s]
          State_Met%U10M        (II,JJ) = atmsrf%usavg(i,j)

          ! Friction velocity [m/s]
          State_Met%USTAR       (II,JJ) = atmsrf%ustar_pbl(i,j)

#ifdef CALC_MERRA2_LIKE_DIAGS
          ! UV surface albedo [1]
          State_Met%UVALBEDO    (II,JJ) = save_alb(i,j)
#endif

          ! N/S wind speed @ 10m ht [m/s]
          State_Met%V10M        (II,JJ) = atmsrf%vsavg(i,j)

#ifdef CALC_MERRA2_LIKE_DIAGS
          ! Surface roughness height [m]
          State_Met%Z0          (II,JJ) = z0m_save(i,j)
#endif

       ENDDO
    ENDDO

    DO K=1,LM
       DO JJJ=J_0,J_1
          DO III=I_0,I_1

             ! GEOS-Chem local index
             II = III - I_0 + 1
             JJ = JJJ - J_0 + 1

             ! GISS meteorology index (GISS only has one polar box)
             I = III
             J = JJJ
             if(hassouthpole(grid) .and. JJJ .eq. J_0 ) I = 1
             if(hasnorthpole(grid) .and. JJJ .eq. J_1 ) I = 1

             ! Dry air mass [kg]
             State_Met%AD          (II,JJ,K) = mma(i,j,k)

             ! Dry air density [kg m-3]
             ! NOTE: Formula from DIFFG function definition in model/DRYDEF.f with pressure PMID
             State_Met%AIRDEN      (II,JJ,K) = PMID(II,JJ,K)*avog*bygasc/t(i,j,k)*pk(k,i,j)

             ! Volume of grid box [m3]
             ! TODO: State_Met%AIRVOL      (II,JJ,K) = ???
             ! NOTE: We could calculate this using State_Grid%XEdge and State_Grid%YEdge but what
             !       about the edge cases?

             ! Water vapor mixing ratio (w/r/t dry air)
             ! NOTE: QSAT is a function defined in model/shared/Utilities.F90
             State_Met%AVGW        (II,JJ,K) = QSAT(t(i,j,k)*pk(k,i,j),LHE,pmid(k,i,j))

             ! Grid box height [m]
             ! NOTE: This is the geopotential height from Model E. Do we need to deduce the heights
             !       from the GEOS-Chem State_Grid derived type?
             State_Met%BXHEIGHT    (II,JJ,K) = gz(i,j,k)

#ifdef CALC_MERRA2_LIKE_DIAGS
             ! 3-D cloud fraction [1]
             State_Met%CLDF        (II,JJ,K) = min(1.0,CLDSS3D(k,i,j) + CLDMC(k,i,j))

             ! Cloud mass flux [kg/m2/s]
             State_Met%CMFMC       (II,JJ,K) = cmfmc(i,j,k)

#endif

             ! Delta-pressure across grid box(wet air) [hPa]
             ! NOTE: Based on EdgePressure_GISS and DeltPressure_GISS from model/FV_UTILS.f
             ! NOTE: Some of the conditional logic related to the RESOLUTION module was dropped
             !       here because we aren't using GISS grids.
             ! NOTE: Surely the calculation is flawed at the top level?
             ! NOTE: Uses PMID for pressure, which probably isn't right.
             State_Met%DELP        (II,JJ,K) = (SIGE(K)*PMID(II,JJ,K) - SIGE(K+1)*PMID(II,JJ,K))

             ! Delta-pressure across grid box (dry air) [hPa]
             ! TODO: State_Met%DELPDRY     (II,JJ,K) = ???

#ifdef CALC_MERRA2_LIKE_DIAGS

             ! Conv precip production rate [kg/kg/s] (assume per dry air)
             State_Met%DQRCU       (II,JJ,K) = dqrcu(i,j,k)

             ! LS precip prod rate [kg/kg/s] (assume per dry air)
             State_Met%DQRLSAN     (II,JJ,K) = dqrlsan(i,j,k)

             ! Detrainment flux [kg/m2/s]
             State_Met%DTRAIN      (II,JJ,K) = dtrain(i,j,k)

#endif

             ! Vertical pressure velocity [Pa/s]
             State_Met%OMEGA       (II,JJ,K) = MWs(i,j,k)*byaxyp(i,j)*100.0/dtsrc

#ifdef CALC_MERRA2_LIKE_DIAGS
             ! Visible optical depth [1]
             State_Met%OPTD        (II,JJ,K) = (cldss(k,i,j)*TAUSS(k,i,j) + &
                  cldmc(k,i,j)*TAUMC(k,i,j)  ) / ( cldss(k,i,j) + cldmc(k,i,j) + teeny )
#endif

             ! Pressure (w/r/t moist air) at level edges [hPa]
             State_Met%PEDGE       (II,JJ,K) = pedn(k,i,j)

             ! Pressure (w/r/t dry air) at level edges [hPa]
             ! TODO: State_Met%PEDGE_DRY   (II,JJ,K) = ???

#ifdef CALC_MERRA2_LIKE_DIAGS
             ! Dwn flux ice prec:conv [kg/m2/s]
             State_Met%PFICU       (II,JJ,K) = pficu(i,j,k)

             ! Dwn flux ice prec:LS+anv [kg/m2/s]
             State_Met%PFILSAN     (II,JJ,K) = pfilsan(i,j,k)

             ! Dwn flux liq prec:conv [kg/m2/s]
             State_Met%PFLCU       (II,JJ,K) = pflcu(i,j,k)

             ! Dwn flux ice prec:LS+anv [kg/m2/s]
             State_Met%PFLLSAN     (II,JJ,K) = pfllsan(i,j,k)
#endif

             ! Ice mixing ratio [kg/kg dry air]
             State_Met%QI          (II,JJ,K) = qci(i,j,k)

             ! Water mixing ratio [kg/kg dry air]
             State_Met%QL          (II,JJ,K) = qcl(i,j,k)

#ifdef CALC_MERRA2_LIKE_DIAGS
             ! Evap of precip conv [kg/kg/s] (assume per dry air)
             State_Met%REEVAPCN    (II,JJ,K) = reevapcn(i,j,k)

             ! Evap of precip LS+anvil [kg/kg/s] (assume per dry air)
             State_Met%REEVAPLS    (II,JJ,K) = reevapls(i,j,k)
#endif

             ! Relative humidity [%]
             State_Met%RH          (II,JJ,K) = 100.*q(i,j,k)/State_Met%AVGW(II,JJ,K)
             IF ( IT_IS_NAN( State_Met%RH(II,JJ,K) ) ) THEN
                WRITE(6,*) II,JJ,K, q(i,j,k), State_Met%AVGW(II,JJ,K), &
                     t(i,j,k), pk(k,i,j), LHE, pmid(k,i,j)
                CALL STOP_MODEL("Bad RH",255)
             ENDIF

             ! Specific humidity [g H2O/kg tot air]
             State_Met%SPHU        (II,JJ,K) = q(i,j,k)

             ! Specific humidity at start of timestep [g/kg]
             State_Met%SPHU1       (II,JJ,K) = q(i,j,k)

             ! Specific humidity at end of timestep [g/kg]  
             State_Met%SPHU2       (II,JJ,K) = q(i,j,k)

             ! Temperature [K]
             State_Met%T           (II,JJ,K) = t(i,j,k)*pk(k,i,j)

             ! Potential temperature [K]
             State_Met%THETA       (II,JJ,K) = t(i,j,k)

#ifdef CALC_MERRA2_LIKE_DIAGS
             ! Optical depth of ice clouds [1]
             State_Met%TAUCLI      (II,JJ,K) = taui3d(i,j,k)

             ! Optical depth of H2O clouds [1]
             State_Met%TAUCLW      (II,JJ,K) = tauw3d(i,j,k)
#endif

             ! Temperature at start of timestep [K]
             State_Met%TMPU1       (II,JJ,K) = t(i,j,k)*pk(k,i,j)

             ! Temperature at end of timestep [K]
             State_Met%TMPU2       (II,JJ,K) = t(i,j,k)*pk(k,i,j)

             ! E/W component of wind [m s-1]
             State_Met%U           (II,JJ,K) = ualij(k,i,j)

#ifdef MODEL_GEOS
             ! Updraft vertical velocity [hPa/s] (only used by GEOS)
             State_Met%UPDVVEL     (II,JJ,K) = 0d0
#endif

             ! N/S component of wind [m s-1]
             State_Met%V           (II,JJ,K) = valij(k,i,j)

          ENDDO
       ENDDO
    ENDDO

    ! Model top
    DO J=J_0,J_1
       DO I=I_0,I_1

          II = I - I_0 + 1
          JJ = J - J_0 + 1

          State_Met%PEDGE       (II,JJ,LM+1) =  pedn(LM+1,i,j)
          ! TODO: State_Met%PEDGE_DRY   (II,JJ,LM+1) =  ???
#ifdef CALC_MERRA2_LIKE_DIAGS
          State_Met%CMFMC       (II,JJ,LM+1) =  cmfmc(i,j,LM+1)
          State_Met%PFICU       (II,JJ,LM+1) =  pficu(i,j,LM+1)
          State_Met%PFILSAN     (II,JJ,LM+1) =  pfilsan(i,j,LM+1)
          State_Met%PFLCU       (II,JJ,LM+1) =  pflcu(i,j,LM+1)
          State_Met%PFLLSAN     (II,JJ,LM+1) =  pfllsan(i,j,LM+1)
#endif

       ENDDO
    ENDDO

    DO K=1,LM
       DO JJJ=J_0,J_1
          DO III=I_0,I_1

             ! GEOS-Chem local index
             II = III - I_0 + 1
             JJ = JJJ - J_0 + 1

             ! GISS meteorology index (GISS only has one polar box)
             I = III
             J = JJJ
             if(hassouthpole(grid) .and. JJJ .eq. J_0 ) I = 1
             if(hasnorthpole(grid) .and. JJJ .eq. J_1 ) I = 1

             ! TODO: Prod_?PRD? - NOTE: Loop case
             
             ! Production of hydrophilic black carbon from hydrophobic black carbon [kg]
             ProdLoss%ProdBCPIfromBCPO             (II,JJ,K) = prodbcpifrombcpo(i,j,k)
        
             ! Production of hydrophilic organic carbon from hydrophobic organic carbon [kg]
             ProdLoss%ProdOCPIfromOCPO             (II,JJ,K) = prodocpifromocpo(i,j,k)

             ! Production of HMS from aqueous reaction of SO2 and HCHO in clouds [kg S s-1]
             ProdLoss%ProdHMSfromSO2andHCHOinCloud             (II,JJ,K) = prodhmsfromso2andhchoincloud(i,j,k)

             ! Production of SO2 and HCHO from aqueous reaction of HS and OH- in clouds [kg S s-1]
             ProdLoss%ProdSO2andHCHOfromHMSinCloud             (II,JJ,K) = prodso2andhchofromhmsincloud(i,j,k)

             ! Production of SO4 from aqueous oxidation of O3 in clouds [kg S  s-1]
             ProdLoss%ProdSO4fromHMSinCloud             (II,JJ,K) = prodso4fromhmsincloud(i,j,k)

             ! Production of SO4 from aqueous oxidation of H2O2 in clouds [kg S s-1]
             ProdLoss%ProdSO4fromH2O2inCloud             (II,JJ,K) = prodso4fromh2o2incloud(i,j,k)

             ! Production of SO4 from aqueous oxidation of O2 from metals in cloud [kg S]
             ProdLoss%ProdSO4fromO2inCloudMetal             (II,JJ,K) = prodso4fromo2incloudmetal(i,j,k)

             ! Production of SO4 from aqueus oxidation of O3 in cloud [kg S s-1]
             ProdLoss%ProdSO4fromO3inCloud             (II,JJ,K) = prodso4fromo3incloud(i,j,k)

             ! Production of SO4 from O3 in sea salt [kg S s-1]
             ProdLoss%ProdSO4fromO3inSeaSalt             (II,JJ,K) = prodso4fromo3inseasalt(i,j,k)

             ! Production of SO4 from aqueus oxidation of HOBr in clouds [kg S s-1]
             ProdLoss%ProdSO4fromHOBrInCloud             (II,JJ,K) = prodso4fromhobrincloud(i,j,k)

             ! Production of SO4 from sulphur production rate of O3 [kg S s-1]
             ProdLoss%ProdSO4fromSRO3             (II,JJ,K) = prodso4fromsro3(i,j,k)

             ! Production of SO4 from sulphur production rate of HOBr + O3 [kg S s-1]
             ProdLoss%ProdSO4fromSRHObr             (II,JJ,K) = prodso4fromsrhobr(i,j,k)

             ! Production of SO4 from aqueous phase SO3 loss by OH [kg S s-1]
             ProdLoss%ProdSO4fromO3s             (II,JJ,K) = prodso4fromo3s(i,j,k)

             ! TODO: Loss_?LOS? - NOTE: Loop case

             ! Loss of HNO3 on sea salt aerosols [kg s-1]
             ProdLoss%LossHNO3onSeaSalt             (II,JJ,K) = losshno3onseasalt(i,j,k)

             ! Production of CO from CH4 [kg s-1]
             ProdLoss%ProdCOfromCH4             (II,JJ,K) = prodcofromch4(i,j,k)

             ! Production of CO from NMVOCs SO3 - loss by OH [kg s-1]
             ProdLoss%ProdCOfromNMVOC             (II,JJ,K) = prodcofromnmvoc(i,j,k)
          ENDDO
       ENDDO
    ENDDO

    ! Set the pressure at level edges [hPa] from the GCM
    CALL Accept_External_Pedge( State_Met  = State_Met,   &
         State_Grid = State_Grid,  &
         RC         = RC          )

    ! Set dry surface pressure (PS1_DRY) from State_Met%PS1_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 1 )

    ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 2 )

    ! Initialize surface pressures to match the post-advection pressures
    State_Met%PSC2_WET = State_Met%PS1_WET
    State_Met%PSC2_DRY = State_Met%PS1_DRY
    CALL SET_FLOATING_PRESSURES( State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Define airmass and related quantities
    CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, RC, .FALSE. )    
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "AirQnt", 255 )

    ! Cap the polar tropopause pressures at 200 hPa, in order to avoid
    ! tropospheric chemistry from happening too high up (cf. J. Logan)
    CALL GCHP_Cap_Tropopause_Prs( Input_Opt      = Input_Opt,  &
         State_Grid     = State_Grid, &
         State_Met      = State_Met,  &
         RC             = RC         )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "GCHP_Cap_Tropopause_Prs", 255 )

    ! Call PBL quantities. Those are always needed
    CALL Compute_Pbl_Height( Input_Opt, State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "COMPUTE_PBL_HEIGHT", 255 )

    IF ( am_I_Root() ) THEN
       WRITE(6,"(I4.4,A,I2.2,A,I2.2,X,I2.2,A,I2.2,A,I2.2)") &
            YEAR, '-', MONTH, '-', DAY, HOUR, ':', MINUTE, ':', SECOND
    ENDIF

    IF ( FIRST_CHEM ) THEN

       ! Set species units to kg to be put into TrM
       DO N=1, State_Chm%nSpecies
         State_Chm%Species(N)%Units = KG_SPECIES
       ENDDO

       ! Put State_Chm back in TrM
       DO N=1,NTM
          DO L=1,LM
             DO J=J_0,J_1
                DO I=I_0,I_1
                   II = I - I_0 + 1
                   JJ = J - J_0 + 1
                   TrM( I, J, L, N ) = State_Chm%Species(N)%Conc(II,JJ,L)
                ENDDO
             ENDDO
          ENDDO
          CALL HALO_UPDATE( GRID, TrM(:,:,:,N) )
          DO M=1,NMOM
             CALL HALO_UPDATE( GRID, TrMom(M,:,:,:,N) )
          ENDDO
       ENDDO

       ! Initialize PBL quantities from the initial met fields
       CALL Compute_Pbl_Height( Input_Opt, State_Grid, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "COMPUTE_PBL_HEIGHT" at initialization!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       ! Once the initial met fields have been read in, we need to find
       ! the maximum PBL level for the non-local mixing algorithm.
       CALL Max_PblHt_For_Vdiff( Input_Opt, State_Grid, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Max_PblHt_for_Vdiff"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       ! Initialize photolysis, including reading files for optical properties
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. &
            Input_Opt%ITS_AN_AEROSOL_SIM .or. &
            Input_Opt%ITS_A_MERCURY_SIM  ) THEN
          CALL Init_Photolysis( Input_Opt, State_Grid, State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Init_Photolysis"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF
       ENDIF

       FIRST_CHEM = .FALSE.
    ELSE
      ! Convert to kg
      CALL Convert_Spc_Units(                                                  &
          Input_Opt  = Input_Opt,                                              &
          State_Chm  = State_Chm,                                              &
          State_Grid = State_Grid,                                             &
          State_Met  = State_Met,                                              &
          new_units  = KG_SPECIES,                                             &
          RC         = RC                                                    )
      IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Convert_Spc_Units", 255 )
    ENDIF

    ! Copy TrM into State_Chm
    DO N=1,NTM
       DO L=1,LM
          DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                State_Chm%Species(N)%Conc(II,JJ,L) = TrM( I, J, L, N )
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! Convert to v/v dry
    CALL Convert_Spc_Units(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         new_units  = MOLES_SPECIES_PER_MOLES_DRY_AIR,                       &
         RC         = RC                                                    )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Convert_Spc_Units", 255 )

    !=====================
    ! Call GEOS-Chem
    !=====================
    CALL CHEM_CHUNK_RUN(                                 &
         nymd,       nhms,       year,       month,      &
         day,        doy,        hour,       minute,     &
         second,     utc,        hElapsed,   Input_Opt,  &
         State_Chm,  State_Diag, State_Grid, State_Met,  &
         -1, IsChemTime, IsRadTime, RC                         )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Chem_Chunk_Run", 255 )

    ! Convert back to kg for GCM advection
    CALL Convert_Spc_Units(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         new_units  = KG_SPECIES,                                            &
         RC         = RC                                                    )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Convert_Spc_Units", 255 )

    ! Copy State_Chm back in TrM
    DO N=1,NTM
       DO L=1,LM
          DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                TrM( I, J, L, N ) = State_Chm%Species(N)%Conc(II,JJ,L)
             ENDDO
          ENDDO
       ENDDO
       CALL HALO_UPDATE( GRID, TrM(:,:,:,N) )
       DO M=1,NMOM
          CALL HALO_UPDATE( GRID, TrMom(M,:,:,:,N) )
       ENDDO
    ENDDO

    ! IF ( AM_I_ROOT() ) THEN
    !    WRITE(6,*) State_Chm%Species(182)%Conc(1,:,1)
    ! ENDIF
    
    !IF ( AM_I_ROOT() ) THEN
    !   WRITE(6,*) ""
    !   WRITE(6,*) 'O3:', TrM(1,:,1,182)
    !   WRITE(6,*) ""
    !   WRITE(6,*) "OH:", TrM(1,:,1,306)
    !   WRITE(6,*) ""
    !ENDIF

    RETURN

  END SUBROUTINE DO_CHEM

  !==========================================================================================================

  SUBROUTINE TrDYNAM

    USE DOMAIN_DECOMP_ATM, ONLY : AM_I_ROOT, GRID, HALO_UPDATE
    USE TRACER_ADV,        ONLY : AADVQ, sfbm, sfcm

    IMPLICIT NONE

    INTEGER N, M

    IF ( .not. ALLOCATED( sfbm ) ) THEN 
       WRITE(6,*) 'Not allocated yet'
       CALL FLUSH(6)
       STOP
    ENDIF

    ! Uses the fluxes MUs,MVs,MWs from DYNAM and QDYNAM
    DO N=1,NTM

       IF ( IsAdvected(N) ) THEN
          
          CALL HALO_UPDATE( GRID, TrM(:,:,:,n) )
          DO M=1,NMOM
             CALL HALO_UPDATE( GRID, TrMom(M,:,:,:,N) )
          ENDDO
          
          CALL AADVQ( TrM(:,:,:,n), TrMom(:,:,:,:,n), .true., TrName(n) )
       
          CALL HALO_UPDATE( GRID, TrM(:,:,:,n) )
          DO M=1,NMOM
             CALL HALO_UPDATE( GRID, TrMom(M,:,:,:,N) )
          ENDDO

       ENDIF
          
    ENDDO

    RETURN

  END SUBROUTINE TrDYNAM

  !==========================================================================================================

  SUBROUTINE CHEM_Chunk_Run( nymd,       nhms,       year,       month,      &
                             day,        dayOfYr,    hour,       minute,     &
                             second,     utc,        hElapsed,   Input_Opt,  &
                             State_Chm,  State_Diag, State_Grid, State_Met,  &
                             Phase,      IsChemTime, IsRadTime,              &
                             RC )
  

    ! Based on GIGC_Chunk_Run

    ! GEOS-Chem state objects
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState

    ! GEOS-Chem components
    USE Chemistry_Mod,      ONLY : Do_Chemistry, Recompute_OD
    USE Convection_Mod,     ONLY : Do_Convection
    USE DryDep_Mod,         ONLY : Do_DryDep
    USE Emissions_Mod,      ONLY : Emissions_Run
    USE Mixing_Mod,         ONLY : Do_Tend, Do_Mixing
    USE WetScav_Mod,        ONLY : Setup_WetScav, Do_WetDep

    ! HEMCO components (eventually moved to a separate GridComp?)
    USE HCO_State_GC_Mod,   ONLY : HcoState, ExtState
    USE HCO_Interface_Common, ONLY : SetHcoTime
    USE HCO_Interface_GC_Mod, ONLY : Compute_Sflx_For_Vdiff

    ! Specialized subroutines
    USE Calc_Met_Mod,       ONLY : AirQnt
    USE Calc_Met_Mod,       ONLY : Set_Dry_Surface_Pressure
    USE Calc_Met_Mod,       ONLY : Set_Clock_Tracer
    USE Calc_Met_Mod,       ONLY : GCHP_Cap_Tropopause_Prs
    USE Set_Global_CH4_Mod, ONLY : Set_CH4
    USE MODIS_LAI_Mod,      ONLY : Compute_XLAI
    USE PBL_Mix_Mod,        ONLY : Compute_PBL_Height
    USE Pressure_Mod,       ONLY : Set_Floating_Pressures
    USE TOMS_Mod,           ONLY : Compute_Overhead_O3
    USE UCX_Mod,            ONLY : Set_H2O_Trac
    USE Vdiff_Mod,          ONLY : Max_PblHt_for_Vdiff

    ! Utilities
    USE Pressure_Mod,       ONLY : Accept_External_Pedge
    USE State_Chm_Mod,      ONLY : IND_
    USE Time_Mod,           ONLY : Accept_External_Date_Time
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units, KG_SPECIES_PER_KG_DRY_AIR

    ! Diagnostics
    USE Diagnostics_Mod,    ONLY : Zero_Diagnostics_StartofTimestep
    USE Diagnostics_Mod,    ONLY : Set_Diagnostics_EndofTimestep
    USE Diagnostics_Mod,    ONLY : Set_AerMass_Diagnostic

    USE Calc_Met_Mod,       ONLY : GET_COSINE_SZA
    USE Species_Mod,        ONLY : Species

!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: nymd        ! YYYY/MM/DD @ current time
    INTEGER,        INTENT(IN)    :: nhms        ! hh:mm:ss   @ current time
    INTEGER,        INTENT(IN)    :: year        ! UTC year
    INTEGER,        INTENT(IN)    :: month       ! UTC month
    INTEGER,        INTENT(IN)    :: day         ! UTC day
    INTEGER,        INTENT(IN)    :: dayOfYr     ! UTC day of year
    INTEGER,        INTENT(IN)    :: hour        ! UTC hour
    INTEGER,        INTENT(IN)    :: minute      ! UTC minute
    INTEGER,        INTENT(IN)    :: second      ! UTC second
    REAL*4,         INTENT(IN)    :: utc         ! UTC time [hrs]
    REAL*4,         INTENT(IN)    :: hElapsed    ! Elapsed hours
    INTEGER,        INTENT(IN)    :: Phase       ! Run phase (-1, 1 or 2)
    LOGICAL,        INTENT(IN)    :: IsChemTime  ! Time for chemistry?
    LOGICAL,        INTENT(IN)    :: IsRadTime   ! Time for RRTMG?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt   ! Input Options obj
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm   ! Chemistry State obj
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag  ! Diagnostics State obj
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid  ! Grid State obj
    TYPE(MetState),      INTENT(INOUT) :: State_Met   ! Meteorology State obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Jul 2011 - M. Long     - Initial Version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!    TYPE(ESMF_STATE)               :: INTSTATE
!    TYPE(MAPL_MetaComp), POINTER   :: STATE
!    TYPE(ESMF_VM)                  :: VM            ! ESMF VM object
!    TYPE(ESMF_Field)               :: IntField
    REAL*8                         :: DT
    CHARACTER(LEN=512)             :: Iam
    INTEGER                        :: HCO_PHASE, previous_units

    ! Local logicals to turn on/off individual components
    ! The parts to be executed are based on the input options,
    ! the time step and the phase.
    LOGICAL                        :: DoConv
    LOGICAL                        :: DoDryDep
    LOGICAL                        :: DoEmis
    LOGICAL                        :: DoTend
    LOGICAL                        :: DoTurb
    LOGICAL                        :: DoChem
    LOGICAL                        :: DoWetDep
    LOGICAL                        :: DoRad

    ! First call?
    LOGICAL, SAVE                  :: FIRST    = .TRUE.

    ! # of times this routine has been called. Only temporary for printing
    ! processes on the first 10 calls.
    INTEGER, SAVE                  :: NCALLS = 0

    ! Strat. H2O settings
    LOGICAL                        :: SetStratH2O

    ! Whether to scale mixing ratio with meteorology update in AirQnt
    LOGICAL, SAVE                  :: scaleMR = .FALSE.

    ! Debug variables
    INTEGER, parameter             :: I_DBG = 6, J_DBG = 5, L_DBG=1

    !=======================================================================
    ! CHEM_CHUNK_RUN begins here
    !=======================================================================        

    ! Error trap
    Iam = 'GCHP_CHUNK_RUN (gchp_chunk_mod.F90)'
    
    ! Assume success
    RC = GC_SUCCESS

    !=======================================================================
    ! Define processes to be covered in this phase
    !
    ! In the standard GEOS-Chem, the following operator sequence is used:
    ! 1. DryDep (kg)
    ! 2. Emissions (kg)
    ! 3. Turbulence (v/v)
    ! 4. Convection (v/v)
    ! 5. Chemistry (kg)
    ! 6. Wetdep (kg)
    !
    ! The GEOS-5 operator sequence is:
    ! 1. Gravity wave drag
    ! 2. Moist (convection)
    ! 3. Chemistry 1 (drydep and emissions)
    ! 4. Surface 1
    ! 5. Turbulence 1
    ! 6. Surface 2
    ! 7. Turbulence 2
    ! 8. Chemistry 2 (chemistry and wet deposition)
    ! 9. Radiation
    !
    ! Here, we use the following operator sequence:
    !
    ! 1.  Convection (v/v) --> Phase 1
    ! 2.  DryDep (kg)      --> Phase 1
    ! 3.  Emissions (kg)   --> Phase 1
    ! 4a. Tendencies (v/v) --> Phase 1
    ! -------------------------------
    ! 4b. Turbulence (v/v) --> Phase 2
    ! 5.  Chemistry (kg)   --> Phase 2
    ! 6.  WetDep (kg)      --> Phase 2
    !
    ! Any of the listed processes is only executed if the corresponding switch
    ! in the geoschem_config.yml file is enabled. If the physics component
    ! already covers convection or turbulence, they should not be applied here!
    ! The tendencies are only applied if turbulence is not done within
    ! GEOS-Chem (ckeller, 10/14/14).
    !
    ! The standard number of phases in GCHP is 1, set in GCHP.rc, which
    ! results in Phase -1 in gchp_chunk_run. This results in executing
    ! all GEOS-Chem components in a single run rather than splitting up
    ! across two runs as is done in GEOS-5. (ewl, 10/26/18)
    !=======================================================================

    ! By default, do processes as defined in geoschem_config.yml. DoTend
    ! defined below.
    !DoConv   = Input_Opt%LCONV                    ! dynamic time step
    !DoDryDep = Input_Opt%LDRYD .AND. IsChemTime   ! chemistry time step
    !DoEmis   = IsChemTime                         ! chemistry time step
    !DoTurb   = Input_Opt%LTURB                    ! dynamic time step
    !DoChem   = Input_Opt%LCHEM .AND. IsChemTime   ! chemistry time step
    !DoWetDep = Input_Opt%LWETD                    ! dynamic time step
    !DoRad    = Input_Opt%LRAD  .AND. IsRadTime    ! radiation time step

    ! By default, do processes as defined in rundeck
    DoConv   = DoGCConv                          ! dynamic time step
    DoDryDep = DOGCDryDep .AND. IsChemTime       ! chemistry time step
    DoEmis   = DoGCEmis                          ! chemistry time step
    DoTurb   = DoGCTurb                          ! dynamic time step
    DoChem   = DoGCChem .AND. IsChemTime         ! chemistry time step
    DoWetDep = DoGCWetDep                        ! dynamic time step
    DoRad    = .false.

    IF ( Input_Opt%ITS_A_CARBON_SIM ) THEN
       DoDryDep = .false.
       DoWetDep = .false.
    ENDIF
    
    IF ( Input_Opt%AmIRoot .and. NCALLS < 10 ) THEN
       write(6,*) 'DoConv   : ', DoConv
       write(6,*) 'DoDryDep : ', DoDryDep
       write(6,*) 'DoEmis   : ', DoEmis
       write(6,*) 'DoTurb   : ', DoTurb
       write(6,*) 'DoChem   : ', DoChem
       write(6,*) 'DoWetDep : ', DoWetDep
       write(6,*) ' '
    ENDIF

    ! If Phase is not -1, only do selected processes for given phases:
    ! Phase 1: disable turbulence, chemistry and wet deposition.
    IF ( Phase == 1 ) THEN
       DoTurb   = .FALSE.
       DoChem   = .FALSE.
       DoWetDep = .FALSE.

    ! Phase 2: disable convection, drydep and emissions.
    ELSEIF ( Phase == 2 ) THEN
       DoConv   = .FALSE.
       DoDryDep = .FALSE.
       DoEmis   = .FALSE.
    ENDIF

    ! Check if tendencies need be applied. The drydep and emission calls
    ! only calculates the emission / drydep rates, but do not apply the
    ! tendencies to the tracer array yet. If turbulence is done as part of
    ! GEOS-5, we need to make sure that these tendencies are applied to the
    ! tracer array. If turbulence is explicitly covered by GEOS-Chem,
    ! however, the tendencies become automatically applied within the PBL
    ! mixing routines (DO_MIXING), so we should never apply the tendencies
    ! in this case.
    DoTend = ( DoEmis .OR. DoDryDep ) .AND. .NOT. Input_Opt%LTURB

    !-------------------------------------------------------------------------
    ! Pre-Run assignments
    !-------------------------------------------------------------------------

    ! Zero out certain State_Diag arrays. This should not be done in a phase 2
    ! call since this can erase diagnostics filled during phase 1 (e.g., drydep)
    ! (ckeller, 1/21/2022).
    IF ( Phase /= 2 ) THEN
       CALL Zero_Diagnostics_StartOfTimestep( Input_Opt, State_Diag, RC )
    ENDIF
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Zero_Diagnostics_StartOfTimestep", 255 )

    ! Pass time values obtained from the ESMF environment to GEOS-Chem
    CALL Accept_External_Date_Time( value_NYMD     = nymd,       &
                                    value_NHMS     = nhms,       &
                                    value_YEAR     = year,       &
                                    value_MONTH    = month,      &
                                    value_DAY      = day,        &
                                    value_DAYOFYR  = dayOfYr,    &
                                    value_HOUR     = hour,       &
                                    value_MINUTE   = minute,     &
                                    value_HELAPSED = hElapsed,   &
                                    value_UTC      = utc,        &
                                    RC             = RC         )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Accept_External_Date_Time", 255 )
    
    ! Pass time values obtained from the ESMF environment to HEMCO
    CALL SetHcoTime ( HcoState,   ExtState,   year,    month,   day,   &
                      dayOfYr,    hour,       minute,  second,  DoEmis,  RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "SetHcoTime", 255 )
    
    ! Calculate MODIS leaf area indexes needed for dry deposition
    CALL Compute_XLAI( Input_Opt, State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Compute_XLAI", 255 )
    
    ! Set the pressure at level edges [hPa] from the ESMF environment
    CALL Accept_External_Pedge( State_Met  = State_Met,   &
                                State_Grid = State_Grid,  &
                                RC         = RC          )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Accept_External_Pedge", 255 )
    
    ! Set dry surface pressure (PS1_DRY) from State_Met%PS1_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 1 )
    
    ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 2 )

    ! Initialize surface pressures to match the post-advection pressures
    State_Met%PSC2_WET = State_Met%PS1_WET
    State_Met%PSC2_DRY = State_Met%PS1_DRY
    CALL SET_FLOATING_PRESSURES( State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "SET_FLOATING_PRESSURES", 255 )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Define airmass and related quantities
    ! Scale mixing ratio with changing met only if FV advection is off.
    ! Only do this the first timestep if DELP_DRY found in restart.
    IF ( FIRST .and. .not. Input_Opt%LTRAN ) THEN       
       CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, RC, scaleMR )
       scaleMR = .TRUE.
    ELSE
       CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, RC, scaleMR )
    ENDIF
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "AirQnt", 255 )

    ! Initialize/reset wetdep after air quantities computed
    IF ( DoConv .OR. DoChem .OR. DoWetDep ) THEN
       CALL SETUP_WETSCAV( Input_Opt, State_Chm, State_Grid, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "SETUP_WETSCAV", 255 )
    ENDIF

    ! Cap the polar tropopause pressures at 200 hPa, in order to avoid
    ! tropospheric chemistry from happening too high up (cf. J. Logan)
    CALL GCHP_Cap_Tropopause_Prs( Input_Opt      = Input_Opt,  &
                                  State_Grid     = State_Grid, &
                                  State_Met      = State_Met,  &
                                  RC             = RC         )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "GCHP_Cap_Tropopause_Prs", 255 )
    
    ! Update clock tracer if relevant
    IF (  IND_('CLOCK','A') > 0 ) THEN
       CALL Set_Clock_Tracer( State_Chm, State_Grid )
    ENDIF

    ! Call PBL quantities. Those are always needed
    CALL Compute_Pbl_Height( Input_Opt, State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "COMPUTE_PBL_HEIGHT", 255 )

    ! Convert to dry mixing ratio
    CALL Convert_Spc_Units(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         new_units  = KG_SPECIES_PER_KG_DRY_AIR,                             &
         previous_units   = previous_units,                                  &
         RC         = RC                                                    )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "CONVERT_SPC_UNITS", 255 )

    !=======================================================================
    ! Always prescribe H2O in both the stratosphere and troposhere in GEOS.
    ! This is now done right after passing the species from the internal
    ! state to State_Chm (in Chem_GridCompMod.F90). It is important to do it
    ! there to make sure that any H2O tendencies are properly calculated
    ! cakelle2, 2023/10/14
    !=======================================================================
#if !defined( MODEL_GEOS )
    ! SDE 05/28/13: Set H2O to STT if relevant
    IF ( IND_('H2O','A') > 0 ) THEN
       SetStratH2O = .FALSE.
       IF ( Input_Opt%LSETH2O ) THEN
          SetStratH2O = .TRUE.
       ENDIF
       CALL SET_H2O_TRAC( SetStratH2O, Input_Opt, State_Chm, &
                          State_Grid,  State_Met, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "SET_H2O_TRAC", 255 )

      ! Only force strat once
       IF ( Input_Opt%LSETH2O ) Input_Opt%LSETH2O = .FALSE.
    ENDIF
#endif

    !=======================================================================
    ! EMISSIONS. Pass HEMCO Phase 1 which only updates the HEMCO clock
    ! and the HEMCO data list. Should be called every time to make sure
    ! that the HEMCO clock and the HEMCO data list are up to date.
    !=======================================================================
    HCO_PHASE = 1
    CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, DoEmis, HCO_PHASE, RC  )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "EMISSIONS_RUN", 255 )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                                PHASE 1 or -1                           !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 1. Convection
    !
    ! Call GEOS-Chem internal convection routines if convection is enabled
    ! in geoschem_config.yml. This should only be done if convection is not
    ! covered by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================
    IF ( DoConv ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do convection now'

       CALL DO_CONVECTION ( Input_Opt, State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "DO_CONVECTION", 255 )

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Convection done!'
    ENDIF

    !=======================================================================
    ! 2. Dry deposition
    !
    ! Calculates the deposition rates in [s-1].
    !=======================================================================
    IF ( DoDryDep ) THEN

       if(Input_Opt%AmIRoot.and.NCALLS<10) THEN
          write(*,*) ' --- Do drydep now'
          write(*,*) '     Use FULL PBL: ', Input_Opt%PBL_DRYDEP
       endif

       ! Do dry deposition
       CALL Do_DryDep ( Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Do_DryDep", 255 )

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Drydep done!'
    ENDIF

    !=======================================================================
    ! 3. Emissions (HEMCO)
    !
    ! HEMCO must be called on first time step to make sure that the HEMCO
    ! data lists are all properly set up.
    !=======================================================================
    IF ( DoEmis ) THEN

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do emissions now'

       ! Do emissions. Pass HEMCO Phase 2 which performs the emissions
       ! calculations.
       HCO_PHASE = 2
       CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
                           State_Grid, State_Met, DoEmis, HCO_PHASE, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "EMISSIONS_RUN - 2", 255 )

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Emissions done!'

    ENDIF

    !=======================================================================
    ! If physics covers turbulence, simply add the emission and dry
    ! deposition fluxes calculated above to the tracer array, without caring
    ! about the vertical distribution. The tracer tendencies are only added
    ! to the tracers array after emissions, drydep. So we need to use the
    ! emissions time step here.
    !=======================================================================
    IF ( DoTend ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*)   &
                           ' --- Add emissions and drydep to tracers'

       ! Get emission time step [s].
       !_ASSERT(ASSOCIATED(HcoState), 'Error: HcoState not associated')
       DT = HcoState%TS_EMIS

       ! Apply tendencies over entire PBL. Use emission time step.
       CALL DO_TEND( Input_Opt, State_Chm, State_Diag, &
                     State_Grid, State_Met, .FALSE., RC, DT=DT )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "DO_TEND", 255 )

       ! testing only
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*)   &
                                 '     Tendency time step [s]: ', DT

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*)   &
                                 ' --- Fluxes applied to tracers!'
    ENDIF ! Tendencies

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                              PHASE 2 or -1                             !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 4. Turbulence
    !
    ! Call GEOS-Chem internal turbulence routines if turbulence is enabled
    ! in geoschem_config.yml. This should only be done if turbulence is not
    ! covered by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================
    IF ( DoTurb ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do turbulence now'

       ! Only do the following for the non-local PBL mixing
       IF ( Input_Opt%LNLPBL ) THEN

          ! Once the initial met fields have been read in, we need to find
          ! the maximum PBL level for the non-local mixing algorithm.
          ! This only has to be done once. (bmy, 5/28/20)
          IF ( FIRST ) THEN
             CALL Max_PblHt_For_Vdiff( Input_Opt, State_Grid, State_Met, RC )
             IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "MAX_PBLHT_FOR_VDIFF", 255 )
          ENDIF

          ! Compute the surface flux for the non-local mixing,
          ! (which means getting emissions & drydep from HEMCO)
          ! and store it in State_Chm%Surface_Flux
          CALL Compute_Sflx_For_Vdiff( Input_Opt,  State_Chm, State_Diag,    &
                                       State_Grid, State_Met, RC            )
          IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "COMPUTE_SFLX_FOR_VDIFF", 255 )
       ENDIF

       ! Do mixing and apply tendencies. This will use the dynamic time step,
       ! which is fine since this call will be executed on every time step.
       CALL DO_MIXING ( Input_Opt, State_Chm, State_Diag,                    &
                        State_Grid, State_Met, RC                           )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "DO_MIXING", 255 )

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Turbulence done!'
    ENDIF

    ! Set tropospheric CH4 concentrations and fill species array with
    ! current values.
    IF ( Phase /= 2 .AND. Input_Opt%ITS_A_FULLCHEM_SIM  &
         .AND. IND_('CH4','A') > 0 ) THEN

       CALL SET_CH4 ( Input_Opt, State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "SET_CH4", 255 )
    ENDIF

    !=======================================================================
    ! 5. Chemistry
    !=======================================================================
    IF ( DoChem ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do chemistry now'

       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
          ! Calculate TOMS O3 overhead. For now, always use it from the
          ! Met field. State_Met%TO3 is imported from PCHEM (ckeller, 10/21/2014).
          CALL COMPUTE_OVERHEAD_O3( Input_Opt, State_Grid, State_Chm, DAY, &
                                    .TRUE., State_Met%TO3, RC )
       ENDIF

#if !defined( MODEL_GEOS )
       ! Set H2O to species value if H2O is advected
       IF ( IND_('H2O','A') > 0 ) THEN
          CALL SET_H2O_TRAC( .FALSE., Input_Opt, &
                             State_Chm, State_Grid, State_Met, RC )
       ENDIF
#endif

       ! Do chemistry
       CALL Do_Chemistry( Input_Opt, State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Do_Chemistry", 255 )

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Chemistry done!'

    ENDIF

    !=======================================================================
    ! 6. Wet deposition
    !=======================================================================
    IF ( DoWetDep ) THEN

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do wetdep now'

       ! Do wet deposition
       CALL DO_WETDEP( Input_Opt, State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "DO_WETDEP", 255 )

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Wetdep done!'
    ENDIF

    !=======================================================================
    ! Diagnostics
    !=======================================================================

    !==============================================================
    !      ***** U P D A T E  O P T I C A L  D E P T H *****
    !==============================================================
    ! Recalculate the optical depth at the wavelength(s) specified
    ! in the Radiation Menu. This must be done before the call to any
    ! diagnostic and only on a chemistry timestep.
    ! (skim, 02/05/11)
    IF ( DoChem ) THEN
       CALL RECOMPUTE_OD ( Input_Opt, State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "RECOMPUTE_OD", 255 )
    ENDIF

    if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do diagnostics now'

    ! Set certain diagnostics dependent on state at end of step. This
    ! includes species concentration and dry deposition flux.
    ! For GEOS, this is now done in Chem_GridCompMod.F90. This makes sure
    ! that the diagnostics include any post-run updates (e.g., if assimilation
    ! increments are being applied (ckeller, 2/7/22).
#if !defined( MODEL_GEOS )
    CALL Set_Diagnostics_EndofTimestep( Input_Opt,  State_Chm, State_Diag, &
                                        State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Set_Diagnostics_EndofTimestep", 255 )
#endif

    ! Archive aerosol mass and PM2.5 diagnostics
    IF ( State_Diag%Archive_AerMass ) THEN
       CALL Set_AerMass_Diagnostic( Input_Opt,  State_Chm, State_Diag, &
                                    State_Grid, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Set_AerMass_Diagnostic", 255 )
    ENDIF

    if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Diagnostics done!'

    !=======================================================================
    ! Convert State_Chm%Species units
    !=======================================================================
    CALL Convert_Spc_Units(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         new_units    = previous_units,                                      &
         RC         = RC                                                    )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "CONVERT_SPC_UNITS", 255 )

    !=======================================================================
    ! Clean up
    !=======================================================================

    ! testing only
    IF ( HCO_PHASE /= 1 .AND. NCALLS < 10 ) NCALLS = NCALLS + 1

    ! First call is done
    FIRST = .FALSE.

    ! Return success
    RC = GC_SUCCESS

  END SUBROUTINE CHEM_CHUNK_RUN

  !==========================================================================================================

  SUBROUTINE INIT_CHEM( grid, is_coldstart )

    USE DOMAIN_DECOMP_1D,        ONLY : getMpiCommunicator 
    USE DOMAIN_DECOMP_ATM,       ONLY : DIST_GRID, Am_I_Root, getDomainBounds
    USE GEOM,                    ONLY : axyp, lat2d_dg, lon2d_dg
    USE CONSTANT,                ONLY : Pi
    USE MODEL_COM,               ONLY : DTsrc, rsf_file_name
    USE Dictionary_mod,          ONLY : sync_param
    USE CHEM_COM,                ONLY : SpcChmID_to_TrID, t_qlimit, TrFullName, TrID_to_SpcChmID, NSP, SpName
    USE ERROR_MOD,               ONLY : Debug_Msg, Error_Stop, Init_Error

    USE GC_Environment_Mod,      ONLY : GC_Allocate_All
    USE State_Grid_Mod,          ONLY : Init_State_Grid
    USE Input_Opt_Mod,           ONLY : Set_Input_Opt
    USE Input_Mod,               ONLY : Read_Input_File
    USE Time_Mod,                ONLY : GET_NHMS, GET_NHMSb, GET_NYMD, GET_NYMDb, GET_TAU, &
                                        GET_TAUb, SET_TIMESTEPS
    USE TIMERS_MOD,              ONLY : Timer_End, Timer_Start
    USE grid_registry_mod,       ONLY : Init_Grid_Registry
    USE LINOZ_MOD,               ONLY : Linoz_Read
    USE HISTORY_MOD,             ONLY : History_Init
    USE OLSON_LANDMAP_MOD,       ONLY : Compute_Olson_Landmap, Init_LandTypeFrac

    USE Emissions_Mod,           ONLY : Emissions_Init, Emissions_Run
    USE GC_Environment_Mod,      ONLY : GC_Init_StateObj, GC_Init_Extra, GC_Init_Grid
    USE GC_Grid_Mod,             ONLY : SetGridFromCtr
    USE Pressure_Mod,            ONLY : Init_Pressure, Accept_External_ApBp
    USE UCX_MOD,                 ONLY : Init_UCX
    USE PhysConstants,           ONLY : PI_180
    USE State_Chm_Mod,           ONLY : Ind_

    USE LINEAR_CHEM_MOD,         ONLY : Init_Linear_Chem
    USE Photolysis_Mod,          ONLY : Init_Photolysis
    USE Vdiff_Mod,               ONLY : Max_PblHt_for_Vdiff

    USE pario,                   ONLY : par_open, par_close
    USE UnitConv_Mod,            ONLY : Convert_Spc_Units, KG_SPECIES, &
                                        MOLES_SPECIES_PER_MOLES_DRY_AIR

    IMPLICIT NONE

    TYPE (DIST_GRID), INTENT(IN) :: grid
    LOGICAL, INTENT(IN)          :: is_coldstart

    LOGICAL   :: isRoot, prtDebug, TimeForEmis
    INTEGER   :: RC, previous_units

    INTEGER   :: NYMD, NHMS
    REAL*8    :: DT

    INTEGER   :: I, J, L, N, NN, II, JJ, I_0H, I_1H
    INTEGER   :: NYMDb, NHMSb, NHMSe
    INTEGER   :: id_H2O, id_CH4, id_CLOCK

    INTEGER   :: TAU, TAUb

    INTEGER   :: fid
    INTEGER   :: KDISK

    CHARACTER(LEN=255)       :: ThisLoc, historyConfigFile
    CHARACTER(LEN=512)       :: ErrMsg, Instr

    !--------------------------------------------------------------------------
    ! Read the user-defined options for the simulation, etc.
    !--------------------------------------------------------------------------

    isRoot = am_I_root()

    ! Initialize fields of the Input Options object (including amIRoot)
    CALL Set_Input_Opt( isRoot, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Set_Input_Opt"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Initialize fields of the Grid State object
    CALL Init_State_Grid( Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Set_Grid_State"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Read GEOS-Chem input file at very beginning of simulation
    CALL Read_Input_File( Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Read_Input_File"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    Input_Opt%numCPUs = grid%npes_world          ! Number of MPI procs
    Input_Opt%thisCPU = grid%rank                ! Local MPI process handle
    Input_Opt%MPIComm = getMpiCommunicator(grid) ! MPI Communicator Handle
    Input_Opt%isMPI   = .true.                   ! Is this an MPI sim?
    Input_Opt%amIRoot = am_I_root()              ! Is this the root cpu?

    Input_Opt%LTRAN   = .false.                  ! Do not use GEOS-Chem for transport

    CALL sync_param( "DoGCConv",   DoGCConv   )
    CALL sync_param( "DoGCEmis",   DoGCEmis   )
    CALL sync_param( "DoGCTend",   DoGCTend   )
    CALL sync_param( "DoGCTurb",   DoGCTurb   )
    CALL sync_param( "DoGCChem",   DoGCChem   )
    CALL sync_param( "DoGCDryDep", DoGCDryDep )
    CALL sync_param( "DoGCWetDep", DoGCWetDep )

    !================================================================
    ! Specify local domain
    !================================================================

    call getDomainBounds( grid, I_STRT = I_0, I_STOP = I_1, &
         J_STRT = J_0, J_STOP = J_1, &
         I_STRT_HALO = I_0H, I_STOP_HALO = I_1H, &
         J_STRT_HALO = J_0H, J_STOP_HALO = J_1H )

    NI = I_1 - I_0 + 1
    NJ = J_1 - J_0 + 1

    State_Grid%DX           = 2.5e+0_fp
    State_Grid%DY           = 2.0e+0_fp
    State_Grid%XMin         = lon2d_dg(i_0,lbound(lon2d_dg,dim=2))
    State_Grid%XMax         = lon2d_dg(i_1,lbound(lon2d_dg,dim=2))
    State_Grid%YMin         = max( lat2d_dg(1,j_0), -89.0_fp )
    State_Grid%YMax         = min( lat2d_dg(1,j_1),  89.0_fp )

    State_Grid%NX           = NI
    State_Grid%NY           = NJ
    State_Grid%NZ           = LM
    State_Grid%HalfPolar    = .FALSE.
    State_Grid%NestedGrid   = .FALSE.
    State_Grid%NorthBuffer  = 0
    State_Grid%SouthBuffer  = 0
    State_Grid%EastBuffer   = 0
    State_Grid%WestBuffer   = 0

    State_Grid%GlobalNX     = IM
    State_Grid%GlobalNY     = JM
    State_Grid%NativeNZ     = LM
    State_Grid%MaxChemLev   = LM
    State_Grid%MaxStratLev  = LM
    State_Grid%MaxTropLev   = LM
    State_Grid%XMinOffset   = 0
    State_Grid%XMaxOffset   = 0
    State_Grid%YMinOffset   = 0
    State_Grid%YMaxOffset   = 0

    State_Grid%GlobalXMid   => NULL()
    State_Grid%GlobalYMid   => NULL()
    State_Grid%XMid         => NULL()
    State_Grid%XEdge        => NULL()
    State_Grid%YMid         => NULL()
    State_Grid%YEdge        => NULL()
    State_Grid%YMid_R       => NULL()
    State_Grid%YEdge_R      => NULL()
    State_Grid%YSIN         => NULL()
    State_Grid%Area_M2      => NULL()

    ! Initialize GEOS-Chem horizontal grid structure
    CALL GC_Init_Grid( Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error in "GC_Init_Grid"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Call the routine GC_Allocate_All (located in module file
    ! GeosCore/gc_environment_mod.F90) to allocate all lat/lon
    ! allocatable arrays used by GEOS-Chem.
    CALL GC_Allocate_All( Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "GC_Allocate_All"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Store shadow copies of am_I_Root, Input_Opt in error_mod.F
    CALL Init_Error(Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_Error"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Set grid based on passed mid-points
    WRITE(6,*) "Manually Set Grid"

    ! Compute number of grid boxes on global grid
    State_Grid%GlobalNX =   360.0_fp / State_Grid%DX
    if ( State_Grid%HalfPolar ) then
       State_Grid%GlobalNY = ( 180.0_fp / State_Grid%DY ) + 1
    else
       State_Grid%GlobalNY = ( 180.0_fp / State_Grid%DY )
    endif

    !----------------------------------------------------------------------
    ! Calculate grid box centers on global grid
    !----------------------------------------------------------------------

    ! Allocate arrays
    ALLOCATE( State_Grid%GlobalXMid(State_Grid%GlobalNX,State_Grid%GlobalNY), STAT=RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Allocate GlobalXMid", 255 )
    State_Grid%GlobalXMid = 0e+0_fp

    ALLOCATE( State_Grid%GlobalYMid(State_Grid%GlobalNX,State_Grid%GlobalNY), STAT=RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Allocate GlobalYMid", 255 )
    State_Grid%GlobalYMid = 0e+0_fp

    ! Loop over horizontal grid
    DO J = 1, State_Grid%GlobalNY
       DO I = 1, State_Grid%GlobalNX

          !--------------------------------
          ! Longitude centers [degrees]
          !--------------------------------
          State_Grid%GlobalXMid(I,J) = ( State_Grid%DX * (I-1) ) - &
               180e+0_fp + State_Grid%DX / 2d0

          !--------------------------------
          ! Latitude centers [degrees]
          !--------------------------------
          IF ( State_Grid%HalfPolar ) THEN
             IF ( J == 1)  THEN
                ! South Pole
                State_Grid%GlobalYMid(I,J) = -90e+0_fp + (0.25e+0_fp * State_Grid%DY)
             ELSEIF ( J == State_Grid%GlobalNY ) THEN
                ! North Pole
                State_Grid%GlobalYMid(I,J) = +90e+0_fp - (0.25e+0_fp * State_Grid%DY)
             ELSE
                State_Grid%GlobalYMid(I,J) = ( State_Grid%DY * (J-1) ) - 90e+0_fp
             ENDIF
          ELSE                  
             State_Grid%GlobalYMid(I,J) = ( State_Grid%DY * (J-1) ) - 89e+0_fp
          ENDIF

       ENDDO
    ENDDO

    !======================================================================
    ! User-defined Horizontal Grid
    !======================================================================

    ! Determine X offsets based on global grid
    DO I = 1, State_Grid%GlobalNX
       IF ( State_Grid%GlobalXMid(I,1) >= State_Grid%XMin ) THEN
          State_Grid%XMinOffset = I-1
          EXIT
       ENDIF
    ENDDO
    DO I = 1, State_Grid%GlobalNX
       IF ( State_Grid%GlobalXMid(I,1)+State_Grid%DX >= State_Grid%XMax ) THEN
          State_Grid%XMaxOffset = I
          EXIT
       ENDIF
    ENDDO

    ! Determine Y offsets based on global grid
    DO J = 1, State_Grid%GlobalNY
       IF ( State_Grid%GlobalYMid(1,J) >= State_Grid%YMin ) THEN
          State_Grid%YMinOffset = J-1
          EXIT
       ENDIF
    ENDDO
    DO J = 1, State_Grid%GlobalNY
       IF ( State_Grid%GlobalYMid(1,J)+State_Grid%DY >= State_Grid%YMax ) THEN
          State_Grid%YMaxOffset = J
          EXIT
       ENDIF
    ENDDO

    !----------------------------------------------------------------------
    ! Calculate grid box centers and edges on local grid
    !----------------------------------------------------------------------

    DO J = J_0, J_1
       DO I = I_0, I_1

          JJ = J - J_0 + 1
          II = I - I_0 + 1

          State_Grid%XMid(II,JJ)     = lon2d_dg(i,j)                      ! Longitude at center [degE]
          State_Grid%YMid(II,JJ)     = lat2d_dg(i,j)                      ! Latitude at center [degN]

          ! Fix polar boxes, which GISS gives wrong lat value
          IF ( State_Grid%YMid(II,JJ) >  89.0 ) State_Grid%YMid   (II,  JJ) =  89.0
          IF ( State_Grid%YMid(II,JJ) < -89.0 ) State_Grid%YMid   (II,  JJ) = -89.0

          State_Grid%XEdge(II,JJ) = &
               State_Grid%XMid(II,JJ) - State_Grid%DX/2d0  ! Longitude at edge [degE]

          State_Grid%YEdge(II,JJ) = &
               State_Grid%YMid(II,JJ) - State_Grid%DY/2d0  ! Latitude at edge [degN]

          State_Grid%YMid_R(II,JJ)   = State_Grid%YMid(II,JJ)  * PI_180   ! Latitude at center [rad]
          State_Grid%YEdge_R(II,JJ)  = State_Grid%YEdge(II,JJ) * PI_180   ! Latitude at edge [rad]
          State_Grid%YSIN(II,JJ)     = SIN( State_Grid%YEdge_R(II,JJ) )   ! sin(lat) at edge

          State_Grid%Area_M2(II,JJ)  = axyp(i,j)                          ! Grid box area [m2]

          IF ( J .eq. J_1 ) THEN
             State_Grid%YEdge(II,JJ+1) = &
                  State_Grid%YMid(II,JJ) + State_Grid%DY/2d0 ! Latitude at edge [degN]
             State_Grid%YEdge_R(II,JJ+1)  = State_Grid%YEdge(II,JJ+1) * PI_180
             State_Grid%YSIN(II,JJ+1)     = SIN( State_Grid%YEdge_R(II,JJ+1) )
          ENDIF

          IF ( I .eq. IM ) THEN
             State_Grid%XEdge(II+1,JJ) = lon2d_dg(i,j) + State_Grid%DX/2d0 ! Longitude at edge [degE]
          ENDIF

          ! Keep latitudes to -90 to 90 range
          IF ( State_Grid%YEdge(II,JJ+1)   >  90.0     ) State_Grid%YEdge  (II,JJ+1) =  90.0
          IF ( State_Grid%YEdge(II,JJ)     < -90.0     ) State_Grid%YEdge  (II,  JJ) = -90.0
          IF ( State_Grid%YEdge_R(II,JJ+1) > ( Pi/2d0) ) State_Grid%YEdge_R(II,JJ+1) =  Pi / 2d0
          IF ( State_Grid%YEdge_R(II,JJ)   < (-Pi/2d0) ) State_Grid%YEdge_R(II,  JJ) = -Pi / 2d0
          IF ( State_Grid%YSIN(II,JJ+1)    >  1d0      ) State_Grid%YSIN   (II,JJ+1) =  1d0
          IF ( State_Grid%YSIN(II,JJ)      < -1d0      ) State_Grid%YSIN   (II,  JJ) = -1d0 

       ENDDO
    ENDDO

    !======================================================================
    ! Echo info to stdout
    !======================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a)' )
       WRITE( 6, '(''%%%%%%%%%%%%%%% GLOBAL GRID %%%%%%%%%%%%%%%'')' )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%GlobalXMid(I,1), &
            I=1,State_Grid%GlobalNX )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%GlobalYMid(1,J), &
            J=1,State_Grid%GlobalNY )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''%%%%%%%%%%%% USER-DEFINED GRID %%%%%%%%%%%%'')' )
       WRITE( 6, '(a)' )
       WRITE( 6, * ) ' XMinOffset : ', State_Grid%XMinOffset
       WRITE( 6, * ) ' XMaxOffset : ', State_Grid%XMaxOffset
       WRITE( 6, * ) ' YMinOffset : ', State_Grid%YMinOffset
       WRITE( 6, * ) ' YMaxOffset : ', State_Grid%YMaxOffset
       WRITE( 6, '(a)' )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%XMid(I,1), I=1,State_Grid%NX )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box longitude edges [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%XEdge(I,1), I=1,State_Grid%NX+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%YMid(1,J), J=1,State_Grid%NY )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''Grid box latitude edges [degrees]: '')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%YEdge(1,J), J=1,State_Grid%NY+1 )
       WRITE( 6, '(a)' )
       WRITE( 6, '(''SIN( grid box latitude edges )'')' )
       WRITE( 6, '(8(f8.3,1x))' ) ( State_Grid%YSIN(1,J), J=1,State_Grid%NY+1 )
    ENDIF

    ! Set a flag to denote if we should print debug output
    prtDebug = am_I_Root()

    ! Debug output
    IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a READ_INPUT_FILE' )

    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Start( "All diagnostics",           RC )
       CALL Timer_Start( "=> History (netCDF diags)", RC )
    ENDIF

    ! Initialize the Diag_List (list of all diagnostics)
    historyConfigFile = 'HISTORY.rc'
    CALL Init_DiagList( Input_Opt%amIroot, historyConfigFile, Diag_List, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_DiagList"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Initialize the TaggedDiag_List (list of wildcards/tags per diagnostic)
    CALL Init_TaggedDiagList( Input_Opt%amIroot, Diag_List,  &
         TaggedDiag_List,   RC         )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_TaggedDiagList"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    IF ( prtDebug ) THEN
       CALL Print_DiagList( Input_Opt%amIRoot, Diag_List, RC )
       CALL Print_TaggedDiagList( Input_Opt%amIRoot, TaggedDiag_List, RC )
    ENDIF

    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "All diagnostics",           RC )
       CALL Timer_End( "=> History (netCDF diags)", RC )
    ENDIF

    !--------------------------------------------------------------------------
    ! %%%% REPLICATING GCHP FUNCTIONALITY IN EXISTING GEOS-CHEM %%%%
    !
    ! To replicate the functionality of the ESMF interface, we must
    ! initialize the Meteorology State (i.e. State_Met) and the
    ! Chemistry State (i.e. State_Chm) objects.  These objects hold
    ! several individual data fields that need to be passed as
    ! inputs to the chemistry routines.
    !
    ! The Meteorology State has replaced all of the individual
    ! met field arrays contained in module dao_mod.F. Likewise,
    ! the Chemistry State has replaced the STT tracer array
    ! and CSPEC chemical species array.
    !
    ! The Chemistry and Meteorology State objects facilitate using
    ! GEOS-Chem directly from the ESMF interface.  This is the main
    ! reason we are migrating towards used of these objects instead
    ! of the existing ALLOCATABLE module arrays. (bmy, 10/25/12)
    !--------------------------------------------------------------------------

    ! Initialize State_Met, State_Chm, and State_Diag objects
    CALL GC_Init_StateObj( Diag_List       = Diag_List,                       &
         TaggedDiag_List = TaggedDiag_List,                 &
         Input_Opt       = Input_Opt,                       &
         State_Chm       = State_Chm,                       &
         State_Diag      = State_Diag,                      &
         State_Grid      = State_Grid,                      &
         State_Met       = State_Met,                       &
         RC              = RC                              )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "GC_Init_StateObj!"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Copy to State_Met%AREA_M2 to avoid breaking GCHP benchmarks,
    ! which require the AREA_M2 field saved out to the StateMet
    ! diagnostic collection for computing emission totals.
    State_Met%Area_M2 = State_Grid%Area_M2

    CALL sync_param( "DTsrc", DTsrc )   ! GISS chemistry timestep [sec]
    CALL sync_param( "DT",    DT    )   ! GISS dynamic timestep [sec]
    Input_Opt%TS_CHEM = INT( DTsrc  )   
    Input_Opt%TS_EMIS = INT( DTsrc  )   
    Input_Opt%TS_DYN  = INT( DTsrc  )   
    Input_Opt%TS_CONV = INT( DTsrc  )   
    Input_Opt%TS_RAD  = INT( DTsrc  )

    ! Set start and finish time from rundeck
    Input_Opt%NYMDb   = 20141201 ! nymdB
    Input_Opt%NHMSb   = 000000   ! nhmsB
    Input_Opt%NYMDe   = 20141202 ! nymdE
    Input_Opt%NHMSe   = 000000   ! nhmsE

    ! Set GEOS-Chem timesteps on all CPUs
    WRITE(6,*) "Calling SET_TIMESTEPS"
    CALL SET_TIMESTEPS( Input_Opt,                                       &
         Chemistry  = Input_Opt%TS_CHEM,                  &
         Convection = Input_Opt%TS_CONV,                  &
         Dynamics   = Input_Opt%TS_DYN,                   &
         Emission   = Input_Opt%TS_EMIS,                  &
         Radiation  = Input_Opt%TS_RAD,                   &
         Unit_Conv  = MAX( Input_Opt%TS_DYN,              &
         Input_Opt%TS_CONV ),           &
         Diagnos    = Input_Opt%TS_CHEM         )

    !--------------------------------------------------------------------------
    ! For regular simulations, initialize various module arrays etc.
    ! This removes the init calls from the run-stage, which cannot
    ! happen when connecting GEOS-Chem to external ESMs.
    !--------------------------------------------------------------------------
    CALL GC_Init_Extra( Diag_List,  Input_Opt,  State_Chm, &
         State_Diag, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "GC_Init_Extra"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Define advected species ID flags for use below
    id_H2O   = Ind_('H2O',   'A')
    id_CH4   = Ind_('CH4',   'A')
    id_CLOCK = Ind_('CLOCK', 'A')

    !-----------------------------------------------------------------------
    ! OBSPACK Diagnostics: Get information from the species
    ! database for all requested ObsPack output species
    !-----------------------------------------------------------------------
    !    IF ( Input_Opt%Do_ObsPack ) THEN
    !       CALL ObsPack_SpeciesMap_Init( Input_Opt, State_Chm, State_Diag, RC )
    !       IF ( RC /= GC_SUCCESS ) THEN
    !          ErrMsg = 'Error encountered in "ObsPack_SpeciesMap_Init"!'
    !          CALL Error_Stop( ErrMsg, ThisLoc )
    !       ENDIF
    !    ENDIF

    ! LTM: Skipping RRTMG initialization code
    ! LTM: Skipping APM initialization code
    ! LTM: Skipping BPCH_DIAG initialization code

    ! Initialize the GEOS-Chem pressure module (set Ap & Bp)
    CALL Init_Pressure( Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Init_Pressure", 255 )

    ! Set Ap and Bp
    CALL Accept_External_ApBp( State_Grid, Ap, Bp, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Accept_External_ApBp", 255 )

    !--------------------------------------------------------------------------
    ! Register the horizontal and vertical grid information so that
    ! the History component can use it for netCDF metadata
    !--------------------------------------------------------------------------
    CALL Init_Grid_Registry( Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_Grid_Registry"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    !--------------------------------------------------------------------------
    ! Added to read input file for Linoz O3
    !--------------------------------------------------------------------------
    IF ( Input_Opt%LLINOZ ) THEN
       CALL Linoz_Read( Input_Opt, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Linoz_Read"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF
    ENDIF

    ! Define time variables for use below
    NHMS  = GET_NHMS()
    NHMSb = GET_NHMSb()
    NYMD  = GET_NYMD()
    NYMDb = GET_NYMDb()
    TAU   = GET_TAU()
    TAUb  = GET_TAUb()

    !--------------------------------------------------------------------------
    !         ***** H I S T O R Y   I N I T I A L I Z A T I O N *****
    !--------------------------------------------------------------------------
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Start( "All diagnostics",           RC )
       CALL Timer_Start( "=> History (netCDF diags)", RC )
    ENDIF

    ! For now, just hardwire the input file for the History component
    Input_Opt%HistoryInputFile = './HISTORY.rc'

    ! LTM: This is still broken. May need to be called through to set up State_Diag
    ! Initialize the GEOS-Chem history component
    !CALL History_Init( Input_Opt,  State_Met,  State_Chm,                      &
    !     State_Diag, State_Grid, RC                             )

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "History_Init"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "All diagnostics",           RC )
       CALL Timer_End( "=> History (netCDF diags)", RC )
    ENDIF

    !--------------------------------------------------------------------------
    !            ***** I N I T I A L I Z A T I O N  continued *****
    !--------------------------------------------------------------------------

    ! To enable FlexGrid, need to initialize HEMCO and run phase 1
    ! before reading initial metfields.
    ! (Jiawei Zhuang 2017/6)

    ! Initialize HEMCO. This reads the HEMCO configuration file
    ! and creates entries for all data files needed for emission
    ! calculation.
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Start( "HEMCO", RC )
    ENDIF

    CALL Emissions_Init( Input_Opt, State_Chm, State_Grid, State_Met, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Emissions_Init"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Run HEMCO phase 0 as simplfied phase 1 to get initial met fields
    ! and restart file fields
    TimeForEmis = .FALSE.
    CALL Emissions_Run( Input_Opt, State_Chm,   State_Diag, State_Grid,  &
         State_Met, TimeForEmis, 0,          RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Emissions_Run", Phase 0'
       Instr  = 'This error can indicate a missing file. Please check '// &
            'the HEMCO log file for additional error messages! '
       CALL Error_Stop( ErrMsg, ThisLoc, Instr )
    ENDIF

    ! In the case of a cold restart, initialise GEOS-Chem from its restart file
    IF (is_coldstart) THEN
      CALL Get_GC_Restart( Input_Opt, State_Chm, State_Grid, State_Met, RC )

      ! IF ( AM_I_ROOT() ) THEN
      !   WRITE(6,*) State_Chm%Species(182)%Conc(1,:,1)
      ! ENDIF
  
      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Get_GC_Restart"'
        Instr  = ''
        CALL Error_Stop( ErrMsg, ThisLoc, Instr )
      ENDIF
    ENDIF

    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End ( "HEMCO", RC )
    ENDIF

    
    ! Populate the State_Met%LandTypeFrac field with data from HEMCO
    CALL Init_LandTypeFrac( Input_Opt, State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_LandTypeFrac"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF
    
    ! Compute the Olson landmap fields of State_Met
    ! (e.g. State_Met%IREG, State_Met%ILAND, etc.)
    CALL Compute_Olson_Landmap( Input_Opt, State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Compute_Olson_Landmap"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF
  
    !==========================================================================
    !            *****  I N I T I A L   C O N D I T I O N S *****
    !==========================================================================
    
    ! Initialize the UCX routines
    CALL INIT_UCX( Input_Opt, State_Chm, State_Diag, State_Grid )
    IF ( isRoot ) CALL DEBUG_MSG( '### MAIN: a INIT_UCX' )

    ! Capture initial state of atmosphere for STE flux calc (ltm, 06/10/12)
    IF ( Input_Opt%LINEAR_CHEM ) THEN
       CALL Init_Linear_Chem( Input_Opt,  State_Chm, State_Met, State_Grid, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Linear_Chem"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF
    ENDIF

    !-----------------------------------------------------------------------------
    
    NSP = State_Chm%nSpecies
    NTM = State_Chm%nAdvect
    
    ALLOCATE( TrID_to_SpcChmID(NTM) )
    ALLOCATE( SpcChmID_to_TrID(NSP) )
    TrID_to_SpcChmID = 0
    SpcChmID_to_TrID = 0

    ALLOCATE( SpName(NSP) )
    ALLOCATE( TrName(NTM) )
    ALLOCATE( TrFullName(NTM) )
    ALLOCATE( IsAdvected(NTM) )
    ALLOCATE( t_qlimit(NTM) )
    ALLOCATE(     TrM(       I_0H:I_1H, J_0H:J_1H, LM, NTM  ) )
    ALLOCATE(   TrMom( NMOM, I_0H:I_1H, J_0H:J_1H, LM, NTM ) )

    NN=1
    DO N = 1, NSP
       SpName(N) = TRIM( State_Chm%SpcData(N)%Info%Name )
       IF ( State_Chm%SpcData(N)%Info%Is_Advected ) THEN
 
          TrName(NN) =     TRIM( State_Chm%SpcData(N)%Info%Name )
          TrFullName(NN) = TRIM( State_Chm%SpcData(N)%Info%FullName ) // " (" // &
               TRIM( State_Chm%SpcData(N)%Info%Formula  ) // ")"
          IsAdvected(NN) = State_Chm%SpcData(N)%Info%Is_Advected
          IF ( Am_I_Root() ) WRITE(6,*) NN, N, TrName(NN)

          TrID_to_SpcChmID(NN) = N
          SpcChmID_to_TrID(N) = NN

          NN = NN + 1
       ENDIF
    ENDDO
    t_qlimit(:) = .true.

    !-----------------------------------------------------------------------------

    TrM    = 0d0
    TrMom  = 0d0
    IF (is_coldstart) THEN
      !------------------------------------------------------------------------
      ! In the case of a cold restart, copy State_Chm into TrM with the
      ! appropriate units
      !------------------------------------------------------------------------
      DO N=1,NTM
         DO L=1,LM
            DO J=J_0,J_1
               DO I=I_0,I_1
                  II = I - I_0 + 1
                  JJ = J - J_0 + 1
                  TrM( I, J, L, N ) = State_Chm%Species(N)%Conc(II,JJ,L)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
    ELSE
      !------------------------------------------------------------------------
      ! In the case of a non-cold-restart, initialise GEOS-Chem from the Model
      ! E restart file
      !------------------------------------------------------------------------

      ! Determine which was the latest restart file to be written to
      call find_later_rsf(KDISK)

      ! NOTE: Tried reading with io_rsf rather than the manual code below but it gave an MPI abort
      ! USE MODEL_COM, only : ioread, Itime
      ! INTEGER :: ioerr
      ! call io_rsf(rsf_file_name(KDISK),Itime,ioread,ioerr)

      ! Read the TrM and TrMom values from the restart file in parallel
      fid = par_open( grid, trim(rsf_file_name(KDISK))//'.nc', 'read' )
      CALL IO_CHEM( fid, 'read_dist' )
      call par_close( grid, fid )

      ! Species are read from restart file in units of kg
      DO N=1, State_Chm%nSpecies
        State_Chm%Species(N)%Units = KG_SPECIES
      ENDDO
      DO N=1,NTM
         DO L=1,LM
            DO J=J_0,J_1
               DO I=I_0,I_1
                  II = I - I_0 + 1
                  JJ = J - J_0 + 1
                  State_Chm%Species(N)%Conc(II,JJ,L) = TrM( I, J, L, N )
               ENDDO
            ENDDO
         ENDDO
      ENDDO
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    RETURN
  END SUBROUTINE INIT_CHEM

  !==========================================================================================================

#ifdef CACHED_SUBDD
  SUBROUTINE accumGCsubdd

    use domain_decomp_atm, only : grid, am_i_root
    use subdd_mod, only : subdd_groups,subdd_type,subdd_ngroups, &
         inc_subdd,find_groups, LmaxSUBDD
    ! use geom, only : byaxyp
    ! use atm_com, only : byma
    USE UnitConv_Mod

    implicit none

    integer :: igrp,ngroups,grpids(subdd_ngroups)
    type(subdd_type), pointer :: subdd
    integer :: L, n, k, RC, i, j, ii, jj
    real*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
         grid%j_strt_halo:grid%j_stop_halo) :: sddarr2d
    real*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
         grid%j_strt_halo:grid%j_stop_halo, &
         LM                               ) :: sddarr3d
    ! real*8 :: convert
    integer :: previous_units
    
!    ! 3-D diagnostics of advected tracers on model levels
!    call find_groups('taijlh',grpids,ngroups)
!    do igrp=1,ngroups
!       subdd => subdd_groups(grpids(igrp))
!       do k=1,subdd%ndiags
!          ntm_loop: do n=1,ntm
!             ! tracer 3D mixing ratios (SUBDD names are just tracer name):
!             if( trim(trname(n)) .eq. trim(subdd%name(k)) ) then
!                convert = 1d9 * 28.97d0 / State_Chm%SpcData( TrID_to_SpcChmID(N) )%Info%MW_g ! kg/kg -> ppbv
!                do L=1,LmaxSUBDD
!                   sddarr3d(:,:,L) = &
!                        trm(:,:,L,n)*convert*byaxyp(:,:)*byma(L,:,:)
!                end do
!                call inc_subdd(subdd,k,sddarr3d)
!                exit ntm_loop
!             end if
!          end do ntm_loop
!       enddo ! k
!    enddo ! igroup
    
    ! Convert from kg to dry mixing ratio
    CALL Convert_Spc_Units(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         new_units  = MOLES_SPECIES_PER_MOLES_DRY_AIR,                       &
         previous_units   = previous_units,                                  &
         RC         = RC                                                    )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "CONVERT_SPC_UNITS", 255 )
    
    ! 3-D diagnostics of all tracers (advected and non-advected) on model levels
    call find_groups('taijlh',grpids,ngroups)
    do igrp=1,ngroups
       subdd => subdd_groups(grpids(igrp))
       do k=1,subdd%ndiags
          ntm_loop: do n=1,State_Chm%nSpecies
             ! tracer 3D mixing ratios (SUBDD names are just tracer name):
             if( trim( State_Chm%SpcData(N)%Info%Name ) .eq. trim(subdd%name(k)) ) then
                DO L=1,LmaxSUBDD
                DO J=J_0,J_1
                DO I=I_0,I_1
                   II = I - I_0 + 1
                   JJ = J - J_0 + 1
                   sddarr3d(I,J,L) = State_Chm%Species(N)%Conc(II,JJ,L)
                ENDDO
                ENDDO
                ENDDO
                call inc_subdd(subdd,k,sddarr3d)
                exit ntm_loop
             end if
          end do ntm_loop
       enddo ! k
    enddo ! igroup

    do igrp=1,ngroups
       subdd => subdd_groups(grpids(igrp))
       do k=1,subdd%ndiags
          if ( trim(subdd%name(k)) == "StateMet_AD" ) then
            DO L=1,LmaxSUBDD
            DO J=J_0,J_1
            DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%AD(II,JJ,L)
            ENDDO
            ENDDO
            ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_AIRDEN" ) then
            DO L=1,LmaxSUBDD
            DO J=J_0,J_1
            DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%AIRDEN(II,JJ,L)
            ENDDO
            ENDDO
            ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_AIRVOL" ) then
            DO L=1,LmaxSUBDD
            DO J=J_0,J_1
            DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%AIRVOL(II,JJ,L)
            ENDDO
            ENDDO
            ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_AVGW" ) then
            DO L=1,LmaxSUBDD
            DO J=J_0,J_1
            DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%AVGW(II,JJ,L)
            ENDDO
            ENDDO
            ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_BXHEIGHT" ) then
            DO L=1,LmaxSUBDD
            DO J=J_0,J_1
            DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%BXHEIGHT(II,JJ,L)
            ENDDO
            ENDDO
            ENDDO
          endif

#ifdef CALC_MERRA2_LIKE_DIAGS

          if ( trim(subdd%name(k)) == "StateMet_CLDF" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%CLDF(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_CMFMC" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%CMFMC(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

#endif

          if ( trim(subdd%name(k)) == "StateMet_DELP" ) then
            DO L=1,LmaxSUBDD
            DO J=J_0,J_1
            DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%DELP(II,JJ,L)
            ENDDO
            ENDDO
            ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_DELPDRY" ) then
            DO L=1,LmaxSUBDD
            DO J=J_0,J_1
            DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%DELP_DRY(II,JJ,L)
            ENDDO
            ENDDO
            ENDDO
          endif

#ifdef CALC_MERRA2_LIKE_DIAGS

          if ( trim(subdd%name(k)) == "StateMet_DQRCU" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%DQRCU(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_DQRLSAN" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%DQRLSAN(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_DTRAIN" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%DTRAIN(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

#endif

          if ( trim(subdd%name(k)) == "StateMet_OMEGA" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%OMEGA(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

#ifdef CALC_MERRA2_LIKE_DIAGS

          if ( trim(subdd%name(k)) == "StateMet_OPTD" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%OPTD(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

#endif

          if ( trim(subdd%name(k)) == "StateMet_PEDGE" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%PEDGE(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PEDGEDRY" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%PEDGE_DRY(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

#ifdef CALC_MERRA2_LIKE_DIAGS

          if ( trim(subdd%name(k)) == "StateMet_PFICU" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%PFICU(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PFILSAN" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%PFILSAN(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PFLCU" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%PFLCU(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PFLLSAN" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%PFLLSAN(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

#endif

          if ( trim(subdd%name(k)) == "StateMet_QI" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%QI(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_QL" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%QL(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

#ifdef CALC_MERRA2_LIKE_DIAGS

          if ( trim(subdd%name(k)) == "StateMet_REEVAPCN" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%REEVAPCN(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_REEVAPLS" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%REEVAPLS(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

#endif

          if ( trim(subdd%name(k)) == "StateMet_RH" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%RH(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SPHU" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%SPHU(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SPHU1" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%SPHU1(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SPHU2" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%SPHU2(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_T" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%T(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_THETA" ) then
            DO L=1,LmaxSUBDD
            DO J=J_0,J_1
            DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%THETA(II,JJ,L)
            ENDDO
            ENDDO
            ENDDO
          endif

#ifdef CALC_MERRA2_LIKE_DIAGS

          if ( trim(subdd%name(k)) == "StateMet_TAUCLI" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%TAUCLI(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_TAUCLW" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%TAUCLW(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

#endif

          if ( trim(subdd%name(k)) == "StateMet_TMPU1" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%TMPU1(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_TMPU2" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%TMPU2(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_U" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%U(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

#ifdef MODEL_GEOS

          if ( trim(subdd%name(k)) == "StateMet_UPDVVEL" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%UPDVVEL(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

#endif

          if ( trim(subdd%name(k)) == "StateMet_V" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = State_Met%V(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          ! TODO: Prod_?PRD? - NOTE: Loop case
          
         
          if ( trim(subdd%name(k)) == "ProdLoss_ProdBCPIfromBCPO" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdBCPIfromBCPO(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdOCPIfromOCPO" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdOCPIfromOCPO(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdHMSfromSO2andHCHOinCloud" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdHMSfromSO2andHCHOinCloud(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdSO2andHCHOfromHMSinCloud" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdSO2andHCHOfromHMSinCloud(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdSO4fromHMSinCloud" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdSO4fromHMSinCloud(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdSO4fromH2O2inCloud" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdSO4fromH2O2inCloud(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdSO4fromO2inCloudMetal" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdSO4fromO2inCloudMetal(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdSO4fromO3inCloud" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdSO4fromO3inCloud(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdSO4fromO3inSeaSalt" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdSO4fromO3inSeaSalt(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdSO4fromHOBrInCloud" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdSO4fromHOBrInCloud(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdSO4fromSRO3" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdSO4fromSRO3(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdSO4fromSRHObr" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdSO4fromSRHObr(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdSO4fromO3s" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdSO4fromO3s(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          ! TODO: Loss_?LOS? - NOTE: Loop case
          
          
          if ( trim(subdd%name(k)) == "ProdLoss_LossHNO3onSeaSalt" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%LossHNO3onSeaSalt(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdCOfromCH4" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdCOfromCH4(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif
          
          if ( trim(subdd%name(k)) == "ProdLoss_ProdCOfromNMVOC" ) then
             DO L=1,LmaxSUBDD
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr3d(I,J,L) = ProdLoss%ProdCOfromNMVOC(II,JJ,L)
             ENDDO
             ENDDO
             ENDDO
          endif

          call inc_subdd(subdd,k,sddarr3d)
       enddo ! k
    enddo ! igroup

    ! Convert back to kg species
    CALL Convert_Spc_Units(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         new_units    = previous_units,                                      &
         RC         = RC                                                    )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "CONVERT_SPC_UNITS", 255 )
        
    ! 2-D diagnostics
    call find_groups('taijh',grpids,ngroups)
    do igrp=1,ngroups
       subdd => subdd_groups(grpids(igrp))
       do k=1,subdd%ndiags

#ifdef CALC_MERRA2_LIKE_DIAGS
          if ( trim(subdd%name(k)) == "StateMet_ALBD" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%ALBD(II,JJ)
             ENDDO
             ENDDO
          endif
#endif

          ! NOTE: Not in HISTORY file
          if ( trim(subdd%name(k)) == "StateMet_AREAM2" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%AREA_M2(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_ChemGridLev" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%ChemGridLev(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_CLDFRC" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%CLDFRC(II,JJ)
             ENDDO
             ENDDO
          endif

#ifdef CALC_MERRA2_LIKE_DIAGS
          if ( trim(subdd%name(k)) == "StateMet_CLDTOPS" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%CLDTOPS(II,JJ)
             ENDDO
             ENDDO
          endif
#endif

#ifdef MODEL_GEOS
          if ( trim(subdd%name(k)) == "StateMet_CNVFRC" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%CNV_FRC(II,JJ)
             ENDDO
             ENDDO
          endif
#endif

          if ( trim(subdd%name(k)) == "StateMet_CONVDEPTH" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%CONV_DEPTH(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_EFLUX" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%EFLUX(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_FLASHDENS" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%FLASH_DENS(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_FRCLND" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%FRCLND(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_FRLAKE" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%FRLAKE(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_FRLAND" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%FRLAND(II,JJ)
             ENDDO
             ENDDO
          endif

          ! NOTE: FRLANDIC in HISTORY file
          if ( trim(subdd%name(k)) == "StateMet_FRLANDICE" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%FRLANDICE(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_FROCEAN" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%FROCEAN(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_FRSEAICE" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%FRSEAICE(II,JJ)
             ENDDO
             ENDDO
          endif

          ! NOTE: FRSNO in HISTORY file
          if ( trim(subdd%name(k)) == "StateMet_FRSNOW" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%FRSNOW(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_GWETROOT" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%GWETROOT(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_GWETTOP" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%GWETTOP(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_HFLUX" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%HFLUX(II,JJ)
             ENDDO
             ENDDO
          endif

#ifdef CALC_MERRA2_LIKE_DIAGS
          if ( trim(subdd%name(k)) == "StateMet_LAI" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%LAI(II,JJ)
             ENDDO
             ENDDO
          endif
#endif

          if ( trim(subdd%name(k)) == "StateMet_PARDR" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PARDR(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PARDF" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PARDF(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PBLH" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PBLH(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PHIS" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PHIS(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PRECANV" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PRECANV(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PRECCON" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PRECCON(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PRECLSC" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PRECLSC(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PRECTOT" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PRECTOT(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PS1WET" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PS1_WET(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PS2WET" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PS2_WET(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PSC2WET" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PSC2_WET(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PS1DRY" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PS1_DRY(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_PS2DRY" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%PS2_DRY(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SEAICE00" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SEAICE00(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SEAICE10" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SEAICE10(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SEAICE20" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SEAICE20(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SEAICE30" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SEAICE30(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SEAICE40" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SEAICE40(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SEAICE50" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SEAICE50(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SEAICE60" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SEAICE60(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SEAICE70" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SEAICE70(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SEAICE80" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SEAICE80(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SEAICE90" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SEAICE90(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SLP" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SLP(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SNODP" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SNODP(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SNOMAS" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SNOMAS(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SUNCOS" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SUNCOS(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SUNCOSmid" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SUNCOSmid(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_SWGDN" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%SWGDN(II,JJ)
             ENDDO
             ENDDO
          endif

#ifdef CALC_MERRA2_LIKE_DIAGS

          if ( trim(subdd%name(k)) == "StateMet_TO3" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%TO3(II,JJ)
             ENDDO
             ENDDO
          endif
#endif

          if ( trim(subdd%name(k)) == "StateMet_TropLev" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%TropLev(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_TropHt" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%TropHt(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_TROPP" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%TROPP(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_TS" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%TS(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_TSKIN" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%TSKIN(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_U10M" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%U10M(II,JJ)
             ENDDO
             ENDDO
          endif

          if ( trim(subdd%name(k)) == "StateMet_USTAR" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%USTAR(II,JJ)
             ENDDO
             ENDDO
          endif

#ifdef CALC_MERRA2_LIKE_DIAGS

          if ( trim(subdd%name(k)) == "StateMet_UVALBEDO" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%UVALBEDO(II,JJ)
             ENDDO
             ENDDO
          endif
#endif

          if ( trim(subdd%name(k)) == "StateMet_V10M" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%V10M(II,JJ)
             ENDDO
             ENDDO
          endif

#ifdef CALC_MERRA2_LIKE_DIAGS
          if ( trim(subdd%name(k)) == "StateMet_Z0" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Met%Z0(II,JJ)
             ENDDO
             ENDDO
          endif
#endif

          if ( trim(subdd%name(k)) == "lat2d" ) then
             DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                sddarr2d(I,J) = State_Grid%YMID(II,JJ)
             ENDDO
             ENDDO
          endif

          call inc_subdd(subdd,k,sddarr2d)
       enddo ! k
    enddo ! igroup

  END SUBROUTINE accumGCsubdd

#endif

  !==========================================================================================================     

  SUBROUTINE IO_CHEM( fid, action ) 

    USE ParallelIo_mod,    ONLY : doVar, ParallelIo
    USE domain_decomp_atm, ONLY : grid

    implicit none

    integer, intent(in) :: fid
    character(len=*), intent(in) :: action
    type (ParallelIo) :: handle
    integer :: n

    handle = ParallelIo( grid, fid )

    do n=1,NTM
       call doVar( handle, action,     TrM(:,:,:,n),   'trm_' // trim(TrName(n)) //      '(dist_im,dist_jm,lm)'         )
       call doVar( handle, action, TrMom(:,:,:,:,n), 'trmom_' // trim(TrName(n)) // '(nmom,dist_im,dist_jm,lm)', jdim=3 )
    enddo

    RETURN
  END SUBROUTINE IO_CHEM

  ! Read restart file and put into State_Chm
  SUBROUTINE Get_GC_Restart( Input_Opt, State_Chm, State_Grid, State_Met, RC )
    !
    ! !USES:
    !
    USE CMN_SIZE_Mod,      ONLY : NDUST
    USE Error_Mod,         ONLY : Debug_Msg
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_GetPtr
    USE PhysConstants,     ONLY : AIRMW
    USE Input_Opt_Mod,     ONLY : OptInput
    USE Species_Mod,       ONLY : Species, SpcConc
    USE State_Chm_Mod,     ONLY : ChmState
    USE State_Grid_Mod,    ONLY : GrdState
    USE State_Met_Mod,     ONLY : MetState
    USE Time_Mod,          ONLY : Expand_Date
    USE UnitConv_Mod,      ONLY : KG_SPECIES_PER_KG_DRY_AIR, MOLECULES_SPECIES_PER_CM3, &
                                  Convert_Spc_Units
#ifdef APM
    USE APM_Init_Mod,      ONLY : APMIDS
#endif
    !
    ! !INPUT PARAMETERS:
    !
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid ! Grid State object
    !
    ! !INPUT/OUTPUT PARAMETERS:
    !
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
    !
    ! !OUTPUT PARAMETERS:
    !
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?
    !
    ! !REVISION HISTORY:
    !
    !  09 Feb 2016 - E. Lundgren - Initial version
    !  See https://github.com/geoschem/geos-chem for complete history
    !EOP
    !------------------------------------------------------------------------------
    !BOC
    !
    ! !LOCAL VARIABLES:
    !
    INTEGER                   :: I, J, L, M, N      ! lon, lat, lev, indexes
    INTEGER                   :: previous_units
    LOGICAL                   :: FOUND              ! Found in restart file?
    CHARACTER(LEN=60)         :: Prefix             ! utility string
    CHARACTER(LEN=255)        :: LOC                ! routine location
    CHARACTER(LEN=255)        :: MSG                ! message
    CHARACTER(LEN=255)        :: v_name             ! variable name
    REAL(fp)                  :: MW_g               ! species molecular weight
    REAL(fp)                  :: SMALL_NUM          ! small number threshold

    ! Temporary arrays and pointers
    REAL*4,  TARGET           :: Temp3D(State_Grid%NX,State_Grid%NY, &
         State_Grid%NZ)
    REAL*4,  POINTER          :: Ptr2D(:,:  )
    REAL*4,  POINTER          :: Ptr3D(:,:,:)

    ! For Hg simulation
    CHARACTER(LEN=60)         :: HgSpc

    ! Objects
    TYPE(SpcConc),    POINTER :: Spc(:)
    TYPE(Species),    POINTER :: SpcInfo

    !=================================================================
    ! READ_GC_RESTART begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Initialize pointers
    Ptr2D       => NULL()
    Ptr3D       => NULL()
    SpcInfo     => NULL()

    ! Name of this routine
    LOC = ' -> at Get_GC_Restart (in model/CHEM_DRV.F90)'

    ! Set minimum value threshold for [mol/mol]
    SMALL_NUM = 1.0e-30_fp

    ! Set pointer to species concentrations
    Spc => State_Chm%Species

    !=================================================================
    ! Open GEOS-Chem restart file
    !=================================================================

    ! Write read message to log
    WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
    WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'

    !=================================================================
    ! Read species concentrations from NetCDF or use default
    ! background [mol/mol]; store in State_Chm%Species%Conc in [kg/kg dry]
    !=================================================================

    ! Print header for min/max concentration to log
    WRITE( 6, 110 )
110 FORMAT( 'Min and Max of each species in restart file [mol/mol]:' )

    ! Loop over species
    DO N = 1, State_Chm%nSpecies

       ! Initialize species concentration to all zeroes
       Spc(N)%Conc = 0.e+0_fp

       ! Get info about this species from the species database
       SpcInfo => State_Chm%SpcData(N)%Info
       MW_g    =  SpcInfo%MW_g

       ! Define variable name
       v_name = 'SPC_' // TRIM( SpcInfo%Name )

       ! Initialize temporary array for this species and point to it
       Temp3D = 0.0_fp
       Ptr3D => Temp3D

       ! Get variable from HEMCO and store in local array
       CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM(v_name), &
            Ptr3D,     RC,         FOUND=FOUND )

       ! Check if species data is in file
       IF ( FOUND ) THEN
          SpcInfo%Is_InRestart = .TRUE.
       ELSE
          SpcInfo%Is_InRestart = .FALSE.
       ENDIF

       ! If data is in file, read in as [mol/mol] and convert to
       ! [kg/kg dry]. Otherwise, set to background value [mol/mol]
       ! either stored in species database (advected species all levels and
       ! non-advected species levels in the chemistry grid) or a small number
       ! (non-advected species levels above the chemistry grid) converted to
       ! [kg/kg dry]
       IF ( SpcInfo%Is_InRestart ) THEN

          ! Print the min & max of each species as it is read from
          ! the restart file in mol/mol
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 120 ) N, TRIM( SpcInfo%Name ), &
                  MINVAL( Ptr3D ), MAXVAL( Ptr3D ), SUM ( Ptr3D(:,:,1:State_Grid%NZ) )
120          FORMAT( 'Species ', i3, ', ', a8, ': Min = ', es15.9, &
                  '  Max = ',es15.9, '  Sum = ',es15.9)
          ENDIF

          ! Convert file value [mol/mol] to [kg/kg dry] for storage
          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( I, J, L )
          DO L = 1, State_Grid%NZ
             DO J = 1, State_Grid%NY
                DO I = 1, State_Grid%NX
                   Spc(N)%Conc(I,J,L) = Ptr3D(I,J,L) * MW_g / AIRMW
                ENDDO
             ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ELSE

          ! Set species to the background value converted to [kg/kg dry]
          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( I, J, L )
          ! Loop over all grid boxes
          DO L = 1, State_Grid%NZ
             DO J = 1, State_Grid%NY
                DO I = 1, State_Grid%NX

                   ! For non-advected species at levels above chemistry grid,
                   ! use a small number for background
                   IF ( L > State_Grid%MaxChemLev .and. &
                        .NOT. SpcInfo%Is_Advected ) THEN

                      Spc(N)%Conc(I,J,L) = SMALL_NUM * MW_g / AIRMW

                      ! For all other cases, use the background value
                      ! stored in the species database
                   ELSE

                      Spc(N)%Conc(I,J,L) = SpcInfo%BackgroundVV &
                           * MW_g / AIRMW

                      ! Print to log if debugging is on
                      IF ( Input_Opt%amIRoot .AND. &
                           I == 1 .AND. J == 1 .AND. L == 1 ) THEN
                         WRITE( 6, 140 ) N, TRIM( SpcInfo%Name ), SpcInfo%BackgroundVV
140                      FORMAT('Species ', i3, ', ', a9, &
                              ': Use background = ', es15.9)
                      ENDIF


                   ENDIF

                ENDDO
             ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ENDIF

       ! Free pointer
       SpcInfo => NULL()

    ENDDO

    ! Set species units
    DO N=1, State_Chm%nSpecies
       State_Chm%Species(N)%Units = KG_SPECIES_PER_KG_DRY_AIR
    ENDDO
    
    ! If in debug mode, print out species min and max in [molec/cm3]
    IF ( .false. ) THEN ! Input_Opt%Verbose ) THEN

       ! Convert units
       PRINT *, " "
       PRINT *, "Species min and max in molec/cm3"

       CALL Convert_Spc_Units(                                                &
            Input_Opt  = Input_Opt,                                           &
            State_Chm  = State_Chm,                                           &
            State_Grid = State_Grid,                                          &
            State_Met  = State_Met,                                           &
            new_units  = MOLECULES_SPECIES_PER_CM3,                           &
            previous_units = previous_units,                                  &
            RC         = RC                                                  )

       ! Trap error
       IF ( RC /= GC_SUCCESS ) THEN
          Msg = 'Error returned from Convert_Spc_Units, call #1!'
          CALL GC_Error( Msg, RC, Loc )
          RETURN
       ENDIF

       ! Print values
       DO N = 1, State_Chm%nSpecies
          SpcInfo => State_Chm%SpcData(N)%Info
          WRITE(6,150) N, TRIM( SpcInfo%Name ),         &
               MINVAL( Spc(N)%Conc(:,:,:) ), &
               MAXVAL( Spc(N)%Conc(:,:,:) )
150       FORMAT( 'Species ', i3, ', ', a9,             &
               ': Min = ', es15.9, ', Max = ', es15.9 )
          SpcInfo => NULL()
       ENDDO

       !      ! Convert units back
       !      CALL Convert_Spc_Units(                                                &
       !           Input_Opt  = Input_Opt,                                           &
       !           State_Chm  = State_Chm,                                           &
       !           State_Grid = State_Grid,                                          &
       !           State_Met  = State_Met,                                           &
       !           new_units    = previous_units,                                            &
       !           RC         = RC                                                  )
       !
       !      ! Trap error
       !      IF ( RC /= GC_SUCCESS ) THEN
       !         Msg = 'Error returned from Convert_Spc_Units, call #2!'
       !         CALL GC_Error( Msg, RC, Loc )
       !         RETURN
       !      ENDIF

    ENDIF

    !=========================================================================
    ! Get variables for KPP mechanisms (right now just fullchem and Hg)
    !=========================================================================
    IF ( ( Input_Opt%ITS_A_FULLCHEM_SIM .or.                                  &
         Input_Opt%ITS_A_MERCURY_SIM        ) .and. Input_Opt%LCHEM ) THEN

       !----------------------------------------------------------------------
       ! KPP_HVALUE (count of internal timesteps at each grid box)
       !----------------------------------------------------------------------
       v_name = 'KPP_HVALUE'

       ! Get variable from HEMCO and store in local array
       CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
            Ptr3D,     RC,         FOUND=FOUND )

       ! Check if variable is in file
       IF ( FOUND ) THEN
          State_Chm%KPPHvalue = Ptr3D
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 510 ) ADJUSTL( v_name              ),                  &
                  MINVAL(  State_Chm%KPPHvalue ),                  &
                  MAXVAL(  State_Chm%KPPHvalue ),                  &
                  SUM(     State_Chm%KPPHvalue )
          ENDIF
       ELSE
          State_Chm%KPPHvalue = 0.0_fp
          IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
       ENDIF

       ! Nullify pointer
       Ptr3D => NULL()

       ! FORMAT strings
500    FORMAT( a                                                              )
510    FORMAT( a21, ': Min = ', es15.9, '  Max = ', es15.9, '  Sum = ',es15.9 )
520    FORMAT( a21, ': not found in restart, set to zero'                      )

    ENDIF

    !=========================================================================
    ! Get variables for Soil NOx emissions
    !=========================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       !----------------------------------------------------------------------
       ! WETDEP_N (wet-deposited nitrogen)
       !----------------------------------------------------------------------
       v_name = 'WETDEP_N'

       ! Get variable from HEMCO and store in local array
       CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
            Ptr2D,     RC,         FOUND=FOUND )

       ! Check if variable is in file
       IF ( FOUND ) THEN
          State_Chm%WetDepNitrogen = Ptr2D
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 510 ) ADJUSTL( v_name                   ),             &
                  MINVAL(  State_Chm%WetDepNitrogen ),             &
                  MAXVAL(  State_Chm%WetDepNitrogen ),             &
                  SUM(     State_Chm%WetDepNitrogen )
          ENDIF
       ELSE
          State_Chm%WetDepNitrogen = 0.0_fp
          IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) TRIM( v_name )
       ENDIF

       ! Nullify pointer
       Ptr2D => NULL()

       !----------------------------------------------------------------------
       ! DRYDEP_N (dry-deposited nitrogen)
       !----------------------------------------------------------------------
       v_name = 'DRYDEP_N'

       ! Get variable from HEMCO and store in local array
       CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
            Ptr2D,     RC,         FOUND=FOUND )

       ! Check if variable is in file
       IF ( FOUND ) THEN
          State_Chm%DryDepNitrogen = Ptr2D
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 510 ) ADJUSTL( v_name                   ),             &
                  MINVAL(  State_Chm%DryDepNitrogen ),             &
                  MAXVAL(  State_Chm%DryDepNitrogen ),             &
                  SUM(     State_Chm%DryDepNitrogen )
          ENDIF
       ELSE
          State_Chm%DryDepNitrogen = 0.0_fp
          IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
       ENDIF

       ! Nullify pointer
       Ptr2D => NULL()

    ENDIF

    !=========================================================================
    ! Read variables for sulfate chemistry and aerosols
    !=========================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. &
         Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       !----------------------------------------------------------------------
       ! H2O2_AFTERCHEM
       !----------------------------------------------------------------------
       v_name = 'H2O2_AFTERCHEM'

       ! Get variable from HEMCO and store in local array
       CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
            Ptr3D,     RC,         FOUND=FOUND                )

       ! Check if variable is in file
       IF ( FOUND ) THEN
          State_Chm%H2O2AfterChem = Ptr3D
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 510 ) ADJUSTL( v_name                  ),              &
                  MINVAL(  State_Chm%H2O2AfterChem ),              &
                  MAXVAL(  State_Chm%H2O2AfterChem ),              &
                  SUM(     State_Chm%H2O2AfterChem )
          ENDIF
       ELSE
          State_Chm%H2O2AfterChem = 0.0_fp
          IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
       ENDIF

       ! Nullify pointer
       Ptr3D => NULL()

       !----------------------------------------------------------------------
       ! SO2_AFTERCHEM
       !----------------------------------------------------------------------
       v_name = 'SO2_AFTERCHEM'

       ! Get variable from HEMCO and store in local array
       CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
            Ptr3D,     RC,         FOUND=FOUND                )

       ! Check if variable is in file
       IF ( FOUND ) THEN
          State_Chm%SO2AfterChem = Ptr3D
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 510 ) ADJUSTL( v_name                 ),               &
                  MINVAL(  State_Chm%SO2AfterChem ),               &
                  MAXVAL(  State_Chm%SO2AfterChem ),               &
                  SUM(     State_Chm%SO2AfterChem )
          ENDIF
       ELSE
          State_Chm%SO2AfterChem = 0.0_fp
          IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
       ENDIF

       ! Nullify pointer
       Ptr3D => NULL()

       !----------------------------------------------------------------------
       ! AeroH2O_SNA
       !----------------------------------------------------------------------
       v_name = 'AEROH2O_SNA'

       ! Get variable from HEMCO and store in local array
       CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
            Ptr3D,     RC,         FOUND=FOUND                )

       ! Check if variable is in file
       IF ( FOUND ) THEN
          State_Chm%AeroH2O(:,:,:,NDUST+1) = Ptr3D
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 510 ) ADJUSTL( v_name                           ),     &
                  MINVAL(  State_Chm%AeroH2O(:,:,:,NDUST+1) ),     &
                  MAXVAL(  State_Chm%AeroH2O(:,:,:,NDUST+1) ),     &
                  SUM(     State_Chm%AeroH2O(:,:,:,NDUST+1) )
          ENDIF
       ELSE
          State_Chm%AeroH2O(:,:,:,NDUST+1) = 0.0_fp
          IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
       ENDIF

       ! Nullify pointer
       Ptr3D => NULL()

       !----------------------------------------------------------------------
       ! ORVCsesq
       !----------------------------------------------------------------------
       IF ( Input_Opt%LCARB .AND. Input_Opt%LSOA ) THEN

          v_name = 'ORVCSESQ'

          ! Get variable from HEMCO and store in local array
          CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
               Ptr3D,     RC,         FOUND=FOUND                )

          ! Check if variable is in file
          IF ( FOUND ) THEN
             State_Chm%ORVCsesq(:,:,:) = Ptr3D
             IF ( Input_Opt%amIRoot ) THEN
                WRITE( 6, 510 ) ADJUSTL( v_name                           ),     &
                     MINVAL(  State_Chm%ORVCsesq(:,:,:) ),     &
                     MAXVAL(  State_Chm%ORVCsesq(:,:,:) ),     &
                     SUM(     State_Chm%ORVCsesq(:,:,:) )
             ENDIF
          ELSE
             State_Chm%ORVCsesq(:,:,:) = 0.0_fp
             IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
          ENDIF

          ! Nullify pointer
          Ptr3D => NULL()

       ENDIF

    ENDIF

    !=========================================================================
    ! Read variables for UCX and the HEMCO PARANOx extension
    !=========================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       !----------------------------------------------------------------------
       ! STATE_PSC (needed to initialize UCX)
       !----------------------------------------------------------------------
       v_name = 'STATE_PSC'

       ! Get variable from HEMCO and store in local array
       CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
            Ptr3D,     RC,         FOUND=FOUND                )

       ! Check if variable is in file
       IF ( FOUND ) THEN
          State_Chm%STATE_PSC = Ptr3D
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 510 ) ADJUSTL( v_name              ),                  &
                  MINVAL(  State_Chm%STATE_PSC ),                  &
                  MAXVAL(  State_Chm%STATE_PSC ),                  &
                  SUM(     State_Chm%STATE_PSC )

          ENDIF
       ELSE
          IF ( Input_Opt%amIRoot ) THEN
#ifdef ESMF_
             ! ExtData and HEMCO behave ambiguously - if the file was found
             ! but was full of zeros throughout the domain of interest, it
             ! will result in the same output from ExtData as if the field
             ! was missing from the file. As such, HEMCO cannot distinguish
             ! between a missing file and a field of zeros
             WRITE(6,*) 'PSC restart either all zeros in the '
             WRITE(6,*) 'root domain, or the restart file did '
             WRITE(6,*) 'not contain STATE_PSC. Root domain '
             WRITE(6,*) 'will be initialized PSC-free'
          ENDIF
#else
          WRITE( 6, 500 ) &
               'STATE_PSC not found in restart, initialize PSC-free'
       ENDIF
#endif
    ENDIF

    ! Nullify pointer
    Ptr3D => NULL()

    !----------------------------------------------------------------------
    ! JOH (needed to initialize PARANOx)
    !----------------------------------------------------------------------
    v_name = 'JOH'

    ! Get variable from HEMCO and store in local array
    CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
         Ptr2D,     RC,         FOUND=FOUND                )

    ! Check if variable is in file
    IF ( FOUND ) THEN
       State_Chm%JOH = Ptr2D
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 510 ) ADJUSTL( v_name       ),                         &
               MINVAL(  State_Chm%JOH ),                        &
               MAXVAL(  State_Chm%JOH ),                        &
               SUM(     State_Chm%JOH )
       ENDIF
    ELSE
       State_Chm%JOH = 0.0_fp
       IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
    ENDIF

    ! Nullify pointer
    Ptr2D => NULL()

    !----------------------------------------------------------------------
    ! JNO2 (needed to initialize PARANOx)
    !----------------------------------------------------------------------
    v_name = 'JNO2'

    ! Get variable from HEMCO and store in local array
    CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
         Ptr2D,     RC,         FOUND=FOUND                )

    ! Check if variable is in file
    IF ( FOUND ) THEN
       State_Chm%JNO2 = Ptr2D
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 510 ) ADJUSTL( v_name         ),                       &
               MINVAL(  State_Chm%JNO2 ),                       &
               MAXVAL(  State_Chm%JNO2 ),                       &
               SUM(     State_Chm%JNO2 )
       ENDIF
    ELSE
       State_Chm%JNO2 = 0.0_fp
       IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
    ENDIF
    ! Nullify pointer
    Ptr2D => NULL()

 ENDIF

 !=========================================================================
 ! Read ocean mercury variables
 !=========================================================================
 IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

    ! Print total mass to log
    WRITE( 6, 220 )
220 FORMAT(/, 'Total mass of each ocean and snow Hg species:')

    !----------------------------------------------------------------------
    ! Total Hg in ocean
    !----------------------------------------------------------------------
    DO M = 1, 3

       ! Define variable name
       SELECT CASE( M )
       CASE ( 1 )
          HgSpc    = 'Hg0'
       CASE ( 2 )
          HgSpc    = 'Hg2'
       CASE ( 3 )
          HgSpc    = 'HgP'
       END SELECT
       v_name = 'OCEAN_' // TRIM( HgSpc )

       ! Get variable from HEMCO and store in local array
       CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),          &
            Ptr2D,     RC,         FOUND=FOUND             )

       ! Check if variable is in file
       IF ( FOUND ) THEN

          ! Check for negative concentrations (jaf, 7/6/11)
          DO I = 1, State_Grid%NX
             DO J = 1, State_Grid%NY
                IF ( Ptr2D(I,J) < 0.0d4 ) THEN
                   Ptr2D(I,J) = 0.0d4
                ENDIF
             ENDDO
          ENDDO

          ! Assign ocean mercury data and write total mass to log file
          SELECT CASE( M )
          CASE ( 1 )
             State_Chm%OceanHg0 = Ptr2D
             WRITE( 6, 240 ) TRIM( v_name            ),                 &
                  SUM( State_Chm%OceanHg0 ), 'kg'
          CASE ( 2 )
             State_Chm%OceanHg2 = Ptr2D
             WRITE( 6, 240 ) TRIM( v_name            ),                 &
                  SUM( State_Chm%OceanHg2 ), 'kg'
          CASE ( 3 )
             State_Chm%OceanHgP = Ptr2D
             WRITE( 6, 240 ) TRIM( v_name            ),                 &
                  SUM( State_Chm%OceanHgP ), 'kg'
          END SELECT

       ELSE
          WRITE( 6, 230 ) TRIM( v_name )
       ENDIF

       ! Nullify pointer
       Ptr2D => NULL()

    ENDDO

    !--------------------------------------------------------------
    ! Hg snowpack on land and ocean
    !--------------------------------------------------------------
    DO M = 1, 4

       ! Define variable name prefix
       SELECT CASE( M )
       CASE ( 1 )
          Prefix = 'SNOW_HG_OCEAN'        ! Reducible on ocean
       CASE ( 2 )
          Prefix = 'SNOW_HG_OCEAN_STORED' ! Non-reducible on ocean
       CASE ( 3 )
          Prefix = 'SNOW_HG_LAND'         ! Reducible on land
       CASE ( 4 )
          Prefix = 'SNOW_HG_LAND_STORED'  ! Non-reducible on land
       END SELECT

       v_name = TRIM( Prefix )

       ! Get variable from HEMCO and store in local array
       CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),          &
            Ptr2D,     RC,         FOUND=FOUND             )

       ! Check if variable is in file
       IF ( FOUND ) THEN

          ! Assign ocean mercury data and write total mass to file
          SELECT CASE( M )
          CASE ( 1 )
             State_Chm%SnowHgOcean = Ptr2D
             WRITE( 6, 240 ) TRIM( v_name                     ),        &
                  SUM( State_Chm%SnowHgOcean       ), 'kg'
          CASE ( 2 )
             State_Chm%SnowHgOceanStored = Ptr2D
             WRITE( 6, 240 ) TRIM( v_name                     ),        &
                  SUM( State_Chm%SnowHgOceanStored ),'kg'
          CASE ( 3 )
             State_Chm%SnowHgLand = Ptr2D
             WRITE( 6, 240 ) TRIM( v_name                     ),        &
                  SUM( State_Chm%SnowHgLand        ), 'kg'
          CASE ( 4 )
             State_Chm%SnowHgLandStored = Ptr2D
             WRITE( 6, 240 ) TRIM( v_name                     ),        &
                  SUM( State_Chm%SnowHgLandStored  ), 'kg'
          END SELECT

       ELSE
          WRITE( 6, 230 ) TRIM( v_name )
       ENDIF

    ENDDO

    ! Format strings
230 FORMAT( a24, ' not found in restart file, set to zero')
240 FORMAT( a24, ':   ', es15.9, 1x, a4)

 ENDIF

 !=================================================================
 ! Clean up
 !=================================================================

 ! Free pointer
 Spc => NULL()

 ! Mark end of section in log
 IF ( Input_Opt%Verbose .AND. Input_Opt%amIRoot ) THEN
    CALL DEBUG_MSG('### DONE GET_GC_RESTART')
 ENDIF
 WRITE( 6, '(a)' ) REPEAT( '=', 79 )

END SUBROUTINE Get_GC_Restart


END MODULE CHEM_DRV
!==========================================================================================================

#ifdef CACHED_SUBDD

SUBROUTINE tijlh_defs(arr,nmax,decl_count)
  ! Needs to be outside the module to prevent circular dependencies with SUBDD
  ! 3D tracer outputs (model horizontal grid and layers).
use subdd_mod, only : info_type
! info_type_ is a homemade structure constructor for older compilers
use subdd_mod, only : info_type_
use chem_com, only : ntm, trname, nsp, spname
implicit none
integer :: nmax,decl_count
integer :: n
type(info_type) :: arr(nmax)

decl_count = 0

! First, diagnostics available for all tracers:
do n=1,nsp
   ! 3D mixing ratios (SUBDD string is just tracer name):
   arr(next()) = info_type_(                     &
       sname = trim(spname(n)),                  &
       lname = trim(spname(n))//' mixing ratio', &
       units = 'mol mol-1'                       &
       )
end do ! tracers loop

arr(next()) = info_type_(                        &
    sname = 'StateMet_AD',                       &
    lname = 'StateMet_AD',                       &
    units = 'kg'                                 &
    )

arr(next()) = info_type_(                        &
    sname = 'StateMet_AIRDEN',                   &
    lname = 'StateMet_AIRDEN',                   &
    units = 'kg m-3'                             &
    )

arr(next()) = info_type_(                        &
    sname = 'StateMet_AIRVOL',                   &
    lname = 'StateMet_AIRVOL',                   &
    units = 'm3'                                 &
    )

arr(next()) = info_type_(                        &
    sname = 'StateMet_AVGW',                     &
    lname = 'StateMet_AVGW',                     &
    units = 'vol vol-1'                          &
    )

arr(next()) = info_type_(                        &
    sname = 'StateMet_BXHEIGHT',                 &
    lname = 'StateMet_BXHEIGHT',                 &
    units = 'm'                                  &
    )

#ifdef CALC_MERRA2_LIKE_DIAGS

arr(next()) = info_type_(                        &
     sname = 'StateMet_CLDF',                    &
     lname = 'StateMet_CLDF',                    &
     units = '1'                                 &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_CMFMC',                   &
     lname = 'StateMet_CMFMC',                   &
     units = 'kg m-2 s-1'                        &
     )

#endif

arr(next()) = info_type_(                        &
    sname = 'StateMet_DELP',                     &
    lname = 'StateMet_DELP',                     &
    units = 'hPa'                                &
    )

arr(next()) = info_type_(                        &
    sname = 'StateMet_DELPDRY',                  &
    lname = 'StateMet_DELPDRY',                  &
    units = 'hPa'                                &
    )

#ifdef CALC_MERRA2_LIKE_DIAGS

arr(next()) = info_type_(                        &
     sname = 'StateMet_DQRCU',                   &
     lname = 'StateMet_DQRCU',                   &
     units = 'kg kg-1 s-1'                       &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_DQRLSAN',                 &
     lname = 'StateMet_DQRLSAN',                 &
     units = 'kg kg-1 s-1'                       &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_DTRAIN',                  &
     lname = 'StateMet_DTRAIN',                  &
     units = 'kg / m2 / s'                       &
     )

#endif

arr(next()) = info_type_(                        &
     sname = 'StateMet_OMEGA',                   &
     lname = 'StateMet_OMEGA',                   &
     units = 'Pa s-1'                            &
     )

#ifdef CALC_MERRA2_LIKE_DIAGS

arr(next()) = info_type_(                        &
     sname = 'StateMet_OPTD',                    &
     lname = 'StateMet_OPTD',                    &
     units = '1'                                 &
     )

#endif

arr(next()) = info_type_(                        &
     sname = 'StateMet_PEDGE',                   &
     lname = 'StateMet_PEDGE',                   &
     units = 'hPa'                               &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_PEDGEDRY',                &
     lname = 'StateMet_PEDGEDRY',                &
     units = 'hPa'                               &
     )

#ifdef CALC_MERRA2_LIKE_DIAGS

arr(next()) = info_type_(                        &
     sname = 'StateMet_PFICU',                   &
     lname = 'StateMet_PFICU',                   &
     units = 'kg m-2 s-1'                        &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_PFILSAN',                 &
     lname = 'StateMet_PFILSAN',                 &
     units = 'kg m-2 s-1'                        &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_PFLCU',                   &
     lname = 'StateMet_PFLCU',                   &
     units = 'kg m-2 s-1'                        &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_PFLLSAN',                 &
     lname = 'StateMet_PFLLSAN',                 &
     units = 'kg m-2 s-1'                        &
     )

#endif

arr(next()) = info_type_(                        &
     sname = 'StateMet_QI',                      &
     lname = 'StateMet_QI',                      &
     units = 'kg kg-1'                           &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_QL',                      &
     lname = 'StateMet_QL',                      &
     units = 'kg kg-1'                           &
     )

#ifdef CALC_MERRA2_LIKE_DIAGS

arr(next()) = info_type_(                        &
     sname = 'StateMet_REEVAPCN',                &
     lname = 'StateMet_REEVAPCN',                &
     units = 'kg kg-1 s-1'                       &
     )

! NOTE: In GEOS-Chem/Header/state_met_mod.F90 the units are 'kg '. Typo?
arr(next()) = info_type_(                        &
     sname = 'StateMet_REEVAPLS',                &
     lname = 'StateMet_REEVAPLS',                &
     units = 'kg kg-1 s-1'                       &
     )

#endif

arr(next()) = info_type_(                        &
     sname = 'StateMet_RH',                      &
     lname = 'StateMet_RH',                      &
     units = '%'                                 &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_SPHU',                    &
     lname = 'StateMet_SPHU',                    &
     units = 'g kg-1'                            &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_SPHU1',                   &
     lname = 'StateMet_SPHU1',                   &
     units = 'g kg-1'                            &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_SPHU2',                   &
     lname = 'StateMet_SPHU2',                   &
     units = 'g kg-1'                            &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_T',                       &
     lname = 'StateMet_T',                       &
     units = 'K'                                 &
     )

arr(next()) = info_type_(                        &
    sname = 'StateMet_THETA',                    &
    lname = 'StateMet_THETA',                    &
    units = 'K'                                  &
    )

#ifdef CALC_MERRA2_LIKE_DIAGS

arr(next()) = info_type_(                        &
     sname = 'StateMet_TAUCLI',                  &
     lname = 'StateMet_TAUCLI',                  &
     units = '1'                                 &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_TAUCLW',                  &
     lname = 'StateMet_TAUCLW',                  &
     units = '1'                                 &
     )

#endif

arr(next()) = info_type_(                        &
     sname = 'StateMet_TMPU1',                   &
     lname = 'StateMet_TMPU1',                   &
     units = 'K'                                 &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_TMPU2',                   &
     lname = 'StateMet_TMPU2',                   &
     units = 'K'                                 &
     )

arr(next()) = info_type_(                        &
     sname = 'StateMet_U',                       &
     lname = 'StateMet_U',                       &
     units = 'm s-1'                             &
     )

#ifdef MODEL_GEOS

arr(next()) = info_type_(                        &
     sname = 'StateMet_UPDVVEL',                 &
     lname = 'StateMet_UPDVVEL',                 &
     units = 'hPa s-1'                           &
     )

#endif

arr(next()) = info_type_(                        &
     sname = 'StateMet_V',                       &
     lname = 'StateMet_V',                       &
     units = 'm s-1'                             &
     )

! ProdLoss Collection
! TODO: Prod_?PRD? - NOTE: Loop case

arr(next()) = info_type_(                        &
     sname = 'ProdBCPIfromBCPO',                 &
     lname = 'ProdBCPIfromBCPO',                 &
     units = 'kg'                                &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdBCPIfromOCPO',                 &
     lname = 'ProdBCPIfromOCPO',                 &
     units = 'kg'                                &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdHMSfromSO2andHCHOinCloud',     &
     lname = 'ProdHMSfromSO2andHCHOinCloud',     &
     units = 'kg S s-1'                          &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdSO2andHCHOfromHMSinCloud',     &
     lname = 'ProdSO2andHCHOfromHMSinCloud',     &
     units = 'kg S s-1'                          &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdSO4fromHMSinCloud',            &
     lname = 'ProdSO4fromHMSinCloud',            &
     units = 'kg S s-1'                          &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdSO4fromH2O2inCloud',           &
     lname = 'ProdSO4fromH2O2inCloud',           &
     units = 'kg S s-1'                          &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdSO4fromO2inCloudMetal',        &
     lname = 'ProdSO4fromO2inCloudMetal',        &
     units = 'kg S'                              &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdSO4fromO3inCloud',             &
     lname = 'ProdSO4fromO3inCloud',             &
     units = 'kg S s-1'                          &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdSO4fromO3inSeaSalt',           &
     lname = 'ProdSO4fromO3inSeaSalt',           &
     units = 'kg S s-1'                          &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdSO4fromHOBrInCloud',           &
     lname = 'ProdSO4fromHOBrInCloud',           &
     units = 'kg S s-1'                          &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdSO4fromSRO3',                  &
     lname = 'ProdSO4fromSRO3',                  &
     units = 'kg S s-1'                          &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdSO4fromSRHObr',                 &
     lname = 'ProdSO4fromSRHObr',                &
     units = 'kg S s-1'                          &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdSO4fromO3s',                   &
     lname = 'ProdSO4fromO3s',                   &
     units = 'kg S s-1'                          &
     )
! TODO: Loss_?LOS? - NOTE: Loop case

arr(next()) = info_type_(                        &
     sname = 'LossHNO3onSeaSalt',                &
     lname = 'LossHNO3onSeaSalt',                &
     units = 'kg s-1'                            &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdCOfromCH4',                    &
     lname = 'ProdCOfromCH4',                    &
     units = 'kg s-1'                            &
     )

arr(next()) = info_type_(                        &
     sname = 'ProdCOfromNMVOC',                  &
     lname = 'ProdCOfromNMVOC',                  &
     units = 'kg s-1'                            &
     )
return
contains
integer function next()
 decl_count = decl_count + 1
 next = decl_count
end function next
END SUBROUTINE tijlh_defs

SUBROUTINE tijh_defs(arr,nmax,decl_count)
  ! Needs to be outside the module to prevent circular dependencies with SUBDD
  ! 3D tracer outputs (model horizontal grid and layers).
use subdd_mod, only : info_type
! info_type_ is a homemade structure constructor for older compilers
use subdd_mod, only : info_type_
implicit none
integer :: nmax,decl_count
integer :: n
character*80 :: unitString
type(info_type) :: arr(nmax)

decl_count = 0

!do n=1,1

#ifdef CALC_MERRA2_LIKE_DIAGS
  arr(next()) = info_type_(                      &
       sname = 'StateMet_ALBD',                  &
       lname = 'StateMet_ALBD',                  &
       units = '1'                               & 
       )
#endif

  arr(next()) = info_type_(                      &
       sname = 'StateMet_AREAM2',                &
       lname = 'StateMet_AREAM2',                &
       units = 'm2'                              &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_ChemGridLev',           &
       lname = 'StateMet_ChemGridLev',           &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_CLDFRC',                &
       lname = 'StateMet_CLDFRC',                &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_CLDFRC',                &
       lname = 'StateMet_CLDFRC',                &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_CLDTOPS',               &
       lname = 'StateMet_CLDTOPS',               &
       units = 'level'                           &
       )

#ifdef MODEL_GEOS
  arr(next()) = info_type_(                      &
       sname = 'StateMet_CNVFRC',                &
       lname = 'StateMet_CNVFRC',                &
       units = '1'                               &
       )
#endif

  arr(next()) = info_type_(                      &
       sname = 'StateMet_CONVDEPTH',             &
       lname = 'StateMet_CONVDEPTH',             &
       units = 'm'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_EFLUX',                 &
       lname = 'StateMet_EFLUX',                 &
       units = 'W m-2'                           &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_FLASHDENS',             &
       lname = 'StateMet_FLASHDENS',             &
       units = 'km-2 s-1'                        &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_FRCLND',                &
       lname = 'StateMet_FRCLND',                &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_FRLAKE',                &
       lname = 'StateMet_FRLAKE',                &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_FRLAND',                &
       lname = 'StateMet_FRLAND',                &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_FRLANDICE',             &
       lname = 'StateMet_FRLANDICE',             &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_FROCEAN',               &
       lname = 'StateMet_FROCEAN',               &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_FRSEAICE',              &
       lname = 'StateMet_FRSEAICE',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_FRSNOW',                &
       lname = 'StateMet_FRSNOW',                &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_GWETROOT',              &
       lname = 'StateMet_GWETROOT',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_GWETTOP',               &
       lname = 'StateMet_GWETTOP',               &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_HFLUX',                 &
       lname = 'StateMet_HFLUX',                 &
       units = 'W m-2'                           &
       )

#ifdef CALC_MERRA2_LIKE_DIAGS
  arr(next()) = info_type_(                      &
       sname = 'StateMet_LAI',                   &
       lname = 'StateMet_LAI',                   &
       units = 'm2 m-2'                          &
       )
#endif

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PARDR',                 &
       lname = 'StateMet_PARDR',                 &
       units = 'W m-2'                           &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PARDF',                 &
       lname = 'StateMet_PARDF',                 &
       units = 'W m-2'                           &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PBLH',                  &
       lname = 'StateMet_PBLH',                  &
       units = 'm'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PHIS',                  &
       lname = 'StateMet_PHIS',                  &
       units = 'm'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PRECANV',               &
       lname = 'StateMet_PRECANV',               &
       units = 'kg m-2 s-1'                      &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PRECCON',               &
       lname = 'StateMet_PRECCON',               &
       units = 'mm day-1'                        &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PRECLSC',               &
       lname = 'StateMet_PRECLSC',               &
       units = 'kg m-2 s-1'                      &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PRECTOT',               &
       lname = 'StateMet_PRECTOT',               &
       units = 'mm day-1'                        &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PS1WET',                &
       lname = 'StateMet_PS1WET',                &
       units = 'hPa'                             &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PS2WET',                &
       lname = 'StateMet_PS2WET',                &
       units = 'hPa'                             &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PSC2WET',               &
       lname = 'StateMet_PSC2WET',               &
       units = 'hPa'                             &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PS1DRY',                &
       lname = 'StateMet_PS1DRY',                &
       units = 'hPa'                             &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_PS2DRY',                &
       lname = 'StateMet_PS2DRY',                &
       units = 'hPa'                             &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SEAICE00',              &
       lname = 'StateMet_SEAICE00',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SEAICE10',              &
       lname = 'StateMet_SEAICE10',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SEAICE20',              &
       lname = 'StateMet_SEAICE20',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SEAICE30',              &
       lname = 'StateMet_SEAICE30',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SEAICE40',              &
       lname = 'StateMet_SEAICE40',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SEAICE50',              &
       lname = 'StateMet_SEAICE50',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SEAICE60',              &
       lname = 'StateMet_SEAICE60',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SEAICE70',              &
       lname = 'StateMet_SEAICE70',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SEAICE80',              &
       lname = 'StateMet_SEAICE80',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SEAICE90',              &
       lname = 'StateMet_SEAICE90',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SLP',                   &
       lname = 'StateMet_SLP',                   &
       units = 'hPa'                             &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SNODP',                 &
       lname = 'StateMet_SNODP',                 &
       units = 'm'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SNOMAS',                &
       lname = 'StateMet_SNOMAS',                &
       units = 'kg m-2'                          &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SUNCOS',                &
       lname = 'StateMet_SUNCOS',                &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SUNCOSmid',             &
       lname = 'StateMet_SUNCOSmid',             &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_SWGDN',                 &
       lname = 'StateMet_SWGDN',                 &
       units = 'W m-2'                           &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_TO3',                   &
       lname = 'StateMet_TO3',                   &
       units = 'dobsons'                         &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_TROPLEV',               &
       lname = 'StateMet_TROPLEV',               &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_TROPHT',                &
       lname = 'StateMet_TROPHT',                &
       units = 'km'                              &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_TROPP',                 &
       lname = 'StateMet_TROPP',                 &
       units = 'hPa'                             &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_TS',                    &
       lname = 'StateMet_TS',                    &
       units = 'K'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_TSKIN',                 &
       lname = 'StateMet_TSKIN',                 &
       units = 'K'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_U10M',                  &
       lname = 'StateMet_U10M',                  &
       units = 'm s-1'                           &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_USTAR',                 &
       lname = 'StateMet_USTAR',                 &
       units = 'm s-1'                           &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_UVALBEDO',              &
       lname = 'StateMet_UVALBEDO',              &
       units = '1'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_V10M',                  &
       lname = 'StateMet_V10M',                  &
       units = 'm s-1'                           &
       )

  arr(next()) = info_type_(                      &
       sname = 'StateMet_Z0',                    &
       lname = 'StateMet_Z0',                    &
       units = 'm'                               &
       )

  arr(next()) = info_type_(                      &
       sname = 'lat2d',                          &
       lname = 'lat2d',                          &
       units = 'degrees_north'                   &
       )

!end do ! tracers loop

return
contains
integer function next()
 decl_count = decl_count + 1
 next = decl_count
end function next
END SUBROUTINE tijh_defs

#endif

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
  USE DIAG_COM
  USE CHEM_COM

  USE Input_Opt_Mod,     ONLY : OptInput
  USE State_Chm_Mod,     ONLY : ChmState
  USE State_Grid_Mod,    ONLY : GrdState 
  USE State_Met_Mod,     ONLY : MetState
  USE State_Diag_Mod,    ONLY : DgnState
  USE DiagList_Mod,      ONLY : DgnList
  USE HCO_Types_Mod,     ONLY : ConfigObj
  USE Precision_Mod

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
  TYPE(GrdState)             :: State_Grid ! Grid state
  TYPE(ConfigObj), POINTER   :: HcoConfig
  
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
  LOGICAL                               :: DoGCDiagn    = .false.
  LOGICAL                               :: coupled_chem = .false.

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
    USE MODEL_COM,         ONLY : modelEclock, itime, ItimeI, DTsrc
    USE ATM_COM,           ONLY : pedn, pmid, pk, ptropo, zatmo, mws, t, q, ualij, valij, qci, qcl
    USE CLOUDS_COM,        ONLY : tauss, taumc, cldmc, cldss, cldss3d, pficu, pflcu, pfilsan, pfllsan
    USE CLOUDS_COM,        ONLY : dtrain, dqrcu, dqrlsan, reevapcn, reevapls, cmfmc
    USE FLUXES,            ONLY : atmsrf, atmlnd, prec, precss, focean, fland, flice
    USE GEOM,              ONLY : axyp, byaxyp
    USE GHY_COM,           ONLY : fearth, wearth, aiearth, wfcs, lai_save, z0m_save
    USE LAKES_COM,         ONLY : flake
    USE O3mod,             ONLY : save_to3
    USE RAD_COM,           ONLY : save_alb, cfrac, srdn, fsrdir, srvissurf, cosz1, save_cosz2
    USE RAD_COM,           ONLY : taui3d, tauw3d
    USE SEAICE_COM,        ONLY : si_atm, si_ocn
    USE CONSTANT,          ONLY : bygrav, lhe, tf, teeny

    ! GEOS-Chem modules
    USE HCO_Interface_Mod, ONLY : SetHcoTime
    USE Time_Mod,          ONLY : Accept_External_Date_Time
    USE Emissions_Mod,     ONLY : Emissions_Run         
    USE State_Chm_Mod,     ONLY : Ind_
    USE Calc_Met_Mod,      ONLY : AirQnt
    USE UnitConv_Mod,      ONLY : Convert_Spc_Units
    USE Pressure_Mod,      ONLY : Set_Floating_Pressures
    USE Pressure_Mod,      ONLY : Accept_External_Pedge
    USE Calc_Met_Mod,      ONLY : Set_Dry_Surface_Pressure
    USE Calc_Met_Mod,      ONLY : GIGC_Cap_Tropopause_Prs
    USE PBL_Mix_Mod,       ONLY : Compute_PBL_Height
    USE ERROR_MOD,         ONLY : Safe_Div, IT_IS_NAN

    IMPLICIT NONE

    INTEGER   :: NYMD,    NHMS,     YEAR,    MONTH,    DAY
    INTEGER   :: DOY,     HOUR,     MINUTE,  SECOND
    INTEGER   :: I,       J,        L,       K,        N
    INTEGER   :: II,      JJ,       III,     JJJ,      RC
    INTEGER   :: STATUS
    REAL*4    :: MINUTES, hElapsed, UTC
    REAL*8    :: sElapsed
    LOGICAL   :: IsChemTime
    CHARACTER(LEN=63)  :: OrigUnit

    ! External functions (from shared/Utilities.F90)
    REAL*8 SLP
    REAL*8 QSAT

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

          ! Visible surface albedo [1]
          State_Met%ALBD        (II,JJ) = save_alb(i,j)  

          ! Grid box surface area [cm2]
          State_Met%AREA_M2     (II,JJ) = axyp(i,j)      

          ! Chemistry grid level [1]
          State_Met%ChemGridLev (II,JJ) = LM             

          ! Column cloud fraction [1]
          State_Met%CLDFRC      (II,JJ) = cfrac(i,j)     

          ! Max cloud top height [levels]
          State_Met%CLDTOPS(II,JJ) = 1                   
          DO K = LM, 1, -1
             IF ( State_Met%CMFMC(II,JJ,K) > 0d0 ) THEN
                State_Met%CLDTOPS(II,JJ) = K + 1
                EXIT
             ENDIF
          ENDDO

          ! Latent heat flux [W/m2]
          State_Met%EFLUX       (II,JJ) = -atmsrf%latht(i,j)/dtsrc    

          ! Olson land fraction [1]
          State_Met%FRCLND      (II,JJ) = fland(i,j)                  

          ! Fraction of lake [1]
          State_Met%FRLAKE      (II,JJ) = flake(i,j)                  

          ! Fraction of land [1]
          State_Met%FRLAND      (II,JJ) = fland(i,j)                  

          ! Fraction of land ice [1]
          State_Met%FRLANDIC    (II,JJ) = flice(i,j)                  

          ! Fraction of ocean [1]
          State_Met%FROCEAN     (II,JJ) = focean(i,j)                 

          ! Sfc sea ice fraction [1]
          State_Met%FRSEAICE    (II,JJ) = si_atm%RSI(i,j)*focean(i,j) 

          ! Surface snow fraction [1]
          State_Met%FRSNO       (II,JJ) = 0.0                         
          if ( si_ocn%snowi(i,j) > 0. ) &
               State_Met%FRSNO(II,JJ) = si_atm%rsi(i,j)*flake(i,j)
          if ( atmlnd%SNOWE(i,j) > 0. ) &
               State_Met%FRSNO(II,JJ) = State_Met%FRSNO(II,JJ) + atmlnd%snowfr(i,j)*fearth(i,j)
          State_Met%FRSNO(II,JJ) = min( 1.0, State_Met%FRSNO(II,JJ) )

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
             State_Met%GWETROOT  (II,JJ) = (wearth(i,j)+aiearth(i,j))/(wfcs(i,j)+1e-20)
          else
             State_Met%GWETROOT  (II,JJ) = 1 ! Set to 1 over oceans to match MERRA-2
          end if

          ! Sensible heat flux [W/m2]
          State_Met%HFLUX       (II,JJ) = -atmsrf%sensht(i,j)/dtsrc   

          ! Leaf area index [m2/m2] (online)
          State_Met%LAI         (II,JJ) = lai_save(i,j)               

          ! Land/water/ice indices [1]
          State_Met%LWI         (II,JJ) = 1                        
          if ( focean(i,j) > fearth(i,j) ) State_Met%LWI(II,JJ) = 0
          if ( si_atm%rsi(i,j)*focean(i,j) > 0.5 ) State_Met%LWI(II,JJ) = 2

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

          ! Total precip @ ground [kg/m2/s] -> mm/d
          State_Met%PRECTOT     (II,JJ) = 86400d0*prec(i,j)/dtsrc                                   

          ! LS precip @ ground [kg/m2/s] -> mm/d
          State_Met%PRECLSC     (II,JJ) = 86400d0*precss(i,j)/dtsrc                                 

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

          ! Total overhead O3 column [DU]
          State_Met%TO3         (II,JJ) = save_to3(i,j)                                             

          ! Tropopause pressure [hPa]     
          State_Met%TROPP       (II,JJ) = ptropo(i,j)                                               

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
          
          ! Surface temperature [K]
          State_Met%TS          (II,JJ) = atmsrf%tsavg(i,j) - tf + 273.15                               

          ! Surface skin temperature [K]
          State_Met%TSKIN       (II,JJ) = atmsrf%gtempr(i,j)                                        
          
          ! E/W wind speed @ 10m ht [m/s]
          State_Met%U10M        (II,JJ) = atmsrf%usavg(i,j)                                         

          ! Friction velocity [m/s]
          State_Met%USTAR       (II,JJ) = atmsrf%ustar_pbl(i,j)                                     

          ! UV surface albedo [1]
          State_Met%UVALBEDO    (II,JJ) = save_alb(i,j)                                             

          ! N/S wind speed @ 10m ht [m/s]
          State_Met%V10M        (II,JJ) = atmsrf%vsavg(i,j)                                         

          ! Surface roughness height [m]
          State_Met%Z0          (II,JJ) = z0m_save(i,j)                                             

          ! Convective fraction [1] (only used by GEOS)
          State_Met%CNV_FRC     (II,JJ) = 0.0                                                       

       ENDDO
    ENDDO

    !IF ( am_I_Root() ) WRITE(6,*) State_Met%SUNCOSmid(:,40)

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

             ! 3-D cloud fraction [1]
             State_Met%CLDF        (II,JJ,K) = min(1.0,CLDSS3D(k,i,j) + CLDMC(k,i,j))                   

             ! Cloud mass flux [kg/m2/s]
             State_Met%CMFMC       (II,JJ,K) = cmfmc(i,j,k)                                             

             ! Conv precip production rate [kg/kg/s] (assume per dry air)
             State_Met%DQRCU       (II,JJ,K) = dqrcu(i,j,k)                                             

             ! LS precip prod rate [kg/kg/s] (assume per dry air)
             State_Met%DQRLSAN     (II,JJ,K) = dqrlsan(i,j,k)                                           

             ! Detrainment flux [kg/m2/s]
             State_Met%DTRAIN      (II,JJ,K) = dtrain(i,j,k)                                            

             ! Vertical pressure velocity [Pa/s]
             State_Met%OMEGA       (II,JJ,K) = MWs(i,j,k)*byaxyp(i,j)*100.0/dtsrc                       

             ! Visible optical depth [1]
             State_Met%OPTD        (II,JJ,K) = (cldss(k,i,j)*TAUSS(k,i,j) + &                           
                  cldmc(k,i,j)*TAUMC(k,i,j)  ) / ( cldss(k,i,j) + cldmc(k,i,j) + teeny )

             ! Wet air press @ level edges [hPa]
             State_Met%PEDGE       (II,JJ,K) = pedn(k,i,j)                                              

             ! Dwn flux ice prec:conv [kg/m2/s]
             State_Met%PFICU       (II,JJ,K) = pficu(i,j,k)                                             

             ! Dwn flux ice prec:LS+anv [kg/m2/s]
             State_Met%PFILSAN     (II,JJ,K) = pfilsan(i,j,k)                                           

             ! Dwn flux liq prec:conv [kg/m2/s]
             State_Met%PFLCU       (II,JJ,K) = pflcu(i,j,k)                                             

             ! Dwn flux ice prec:LS+anv [kg/m2/s]
             State_Met%PFLLSAN     (II,JJ,K) = pfllsan(i,j,k)                                           

             ! Ice mixing ratio [kg/kg dry air]
             State_Met%QI          (II,JJ,K) = qci(i,j,k)                                               

             ! Water mixing ratio [kg/kg dry air]
             State_Met%QL          (II,JJ,K) = qcl(i,j,k)                                               

             ! Evap of precip conv [kg/kg/s] (assume per dry air)
             State_Met%REEVAPCN    (II,JJ,K) = reevapcn(i,j,k)                                          

             ! Evap of precip LS+anvil [kg/kg/s] (assume per dry air)
             State_Met%REEVAPLS    (II,JJ,K) = reevapls(i,j,k)                                          

             ! Relative humidity [%]
             State_Met%RH          (II,JJ,K) = 100.*q(i,j,k)/QSAT(t(i,j,k)*pk(k,i,j),LHE,pmid(k,i,j))   
             IF ( IT_IS_NAN( State_Met%RH(II,JJ,K) ) ) THEN
                WRITE(6,*) II,JJ,K, q(i,j,k), QSAT(t(i,j,k)*pk(k,i,j),LHE,pmid(k,i,j)), &
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

             ! Optical depth of ice clouds [1]
             State_Met%TAUCLI      (II,JJ,K) = taui3d(i,j,k)                                            

             ! Optical depth of H2O clouds [1]
             State_Met%TAUCLW      (II,JJ,K) = tauw3d(i,j,k)

             ! Temperature at start of timestep [K]
             State_Met%TMPU1       (II,JJ,K) = t(i,j,k)*pk(k,i,j)                                       

             ! Temperature at end of timestep [K]
             State_Met%TMPU2       (II,JJ,K) = t(i,j,k)*pk(k,i,j)                                       

             ! E/W component of wind [m s-1]
             State_Met%U           (II,JJ,K) = ualij(k,i,j)                                             

             ! Updraft vertical velocity [hPa/s] (only used by GEOS)
             State_Met%UPDVVEL     (II,JJ,K) = 0d0                                                      

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
          State_Met%CMFMC       (II,JJ,LM+1) =  cmfmc(i,j,LM+1)
          State_Met%PFICU       (II,JJ,LM+1) =  pficu(i,j,LM+1)
          State_Met%PFILSAN     (II,JJ,LM+1) =  pfilsan(i,j,LM+1)
          State_Met%PFLCU       (II,JJ,LM+1) =  pflcu(i,j,LM+1)
          State_Met%PFLLSAN     (II,JJ,LM+1) =  pfllsan(i,j,LM+1)

       ENDDO
    ENDDO

    ! Set the pressure at level edges [hPa] from the GCM
    CALL Accept_External_Pedge( State_Met = State_Met,  &
         RC        = RC         )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Accept_External_Pedge", 255 )
    
    ! Set dry surface pressure (PS1_DRY) from State_Met%PS1_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 1 )
    
    ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 2 )
    
    ! Initialize surface pressures to match the post-advection pressures
    CALL SET_FLOATING_PRESSURES( State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "SET_FLOATING_PRESSURES", 255 )
    
    ! Define airmass and related quantities
    CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, RC, .FALSE. )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "AirQnt", 255 )
    
    ! Cap the polar tropopause pressures at 200 hPa, in order to avoid
    ! tropospheric chemistry from happening too high up (cf. J. Logan)
    CALL GIGC_Cap_Tropopause_Prs( Input_Opt      = Input_Opt,  &
                                  State_Grid     = State_Grid, &
                                  State_Met      = State_Met,  &
                                  RC             = RC         )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "GIGC_Cap_Tropopause_Prs", 255 )

    ! Call PBL quantities. Those are always needed
    CALL COMPUTE_PBL_HEIGHT( State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "COMPUTE_PBL_HEIGHT", 255 )

    IF ( FIRST_CHEM ) THEN

       IF ( am_I_Root() ) THEN
       WRITE(6,*) "Chemistry initialization ..."
       ENDIF

       ! Pass time values to GEOS-Chem
       CALL Accept_External_Date_Time( value_NYMD     = nymd,       &
            value_NHMS     = nhms,       &
            value_YEAR     = year,       &
            value_MONTH    = month,      &
            value_DAY      = day,        &
            value_DAYOFYR  = DOY,        &
            value_HOUR     = hour,       &
            value_MINUTE   = minute,     &
            value_HELAPSED = hElapsed,   &
            value_UTC      = utc,        &
            RC             = RC            )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Accept_External_Date_Time", 255 )

       ! Set initial HEMCO time
       CALL SetHcoTime ( .true., RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "SetHcoTime", 255 )

       !=======================================================================
       ! Read Initial Conditions
       !=======================================================================

       CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
            State_Grid, State_Met, .True.,  0,  RC  )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "EMISSIONS_RUN - Phase 0", 255 )

       CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, 'kg', RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Convert_Spc_Units", 255 )

       ! Copy initial conditions into TrM if absolute beginning
       IF ( hElapsed .eq. 0 ) THEN
         DO N=1,NTM
         DO L=1,LM
         DO J=J_0,J_1
         DO I=I_0,I_1
           II = I - I_0 + 1
           JJ = J - J_0 + 1
           TrM( I, J, L, N )      = State_Chm%Species(II,JJ,L,N)
           !TrMom( I, J, L, :, N ) = 0d0 ! Assuming initially uniformly distributed 
           ! Hack until HALO_UPDATE works for initial conditions
           if ( J_0 .ne. J_0H .and. J .eq. J_0 ) TrM( I, J_0H, L, N ) = TrM(I,J,L,N)
           if ( J_1 .ne. J_1H .and. J .eq. J_1 ) TrM( I, J_1H, L, N ) = TrM(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
       ENDIF

       ! Update halos
       !DO N=1,NTM
       !  if ( .not. hasnorthpole(grid) ) CALL HALO_UPDATE( grid, TrM(:,:,:,N),
       !  FROM=NORTH )
       !  if ( .not. hassouthpole(grid) ) CALL HALO_UPDATE( grid, TrM(:,:,:,N),
       !  FROM=SOUTH )
       !ENDDO

       ! Should this be _0H and _1H?
       ALLOCATE ( GCAIJL_out( I_0:I_1,J_0:J_1,LM,NTM),stat=status)

       IF ( am_I_Root() ) THEN
       WRITE(6,*) "                      ... complete!"
       ENDIF

    ENDIF ! FIRST_CHEM

    IF ( am_I_Root() ) THEN
       WRITE(6,"(I4.4,A,I2.2,A,I2.2,X,I2.2,A,I2.2,A,I2.2)") &
            YEAR, '-', MONTH, '-', DAY, HOUR, ':', MINUTE, ':', SECOND
    ENDIF

    ! Copy TrM into State_Chm
    DO N=1,NTM
    DO L=1,LM
    DO J=J_0,J_1
    DO I=I_0,I_1
       II = I - I_0 + 1
       JJ = J - J_0 + 1
       State_Chm%Species(II,JJ,L,N) = TrM( I, J, L, N )
    ENDDO
    ENDDO
    ENDDO
    ENDDO   
    State_Chm%Spc_Units = 'kg'
       
    ! Convert to v/v dry
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
         'v/v dry', RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Convert_Spc_Units", 255 )

    IF ( Am_I_Root() .and. .false. ) THEN
    print*,NTM     ! 278 in 12.9.0 (199 advected)
    print*,i_H2O   ! 61
    print*,i_CO2   ! 200
    print*,i_O3    ! 164
    print*,i_O2    ! 277
    print*,i_NO2   ! 161
    print*,i_N2O   ! 154
    print*,i_CH4   ! 154
    print*,i_CFC11 ! 36
    print*,i_CFC12 ! 21
    print*,i_N2    ! 276
    print*,i_CFCY  ! 0
    print*,i_CFCZ  ! 0
    print*,i_SO2   ! 191
    CALL FLUSH(6)
    ENDIF
    
    !=====================
    ! Call GEOS-Chem
    !=====================
    CALL CHEM_CHUNK_RUN(                                 &
         nymd,       nhms,       year,       month,      &
         day,        doy,        hour,       minute,     &
         second,     utc,        hElapsed,   Input_Opt,  &
         State_Chm,  State_Diag, State_Grid, State_Met,  &
         IsChemTime, RC                                  )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Chem_Chunk_Run", 255 )

    ! Convert back to kg for GCM advection
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
         'kg', RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Convert_Spc_Units", 255 )

    ! Copy State_Chm back in TrM
    DO N=1,NTM
       DO L=1,LM
          DO J=J_0,J_1
             DO I=I_0,I_1
                II = I - I_0 + 1
                JJ = J - J_0 + 1
                TrM( I, J, L, N ) = State_Chm%Species(II,JJ,L,N)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    FIRST_CHEM = .false.

    RETURN

  END SUBROUTINE DO_CHEM
  
  !==========================================================================================================

  SUBROUTINE TrDYNAM

    USE DOMAIN_DECOMP_ATM, ONLY : AM_I_ROOT
    USE TRACER_ADV,     only : AADVQ, sfbm, sbm, sbf, sfcm, scm, scf, safv, sbfv

    IMPLICIT NONE

    INTEGER N

    IF ( .not. ALLOCATED( sfbm ) ) THEN 
       WRITE(6,*) 'Not allocated yet'
       CALL FLUSH(6)
       STOP
    ENDIF

    ! Uses the fluxes MUs,MVs,MWs from DYNAM and QDYNAM
    DO N=1,NTM
       !IF ( Am_I_Root() ) WRITE(6,*) TrName(N), IsAdvected(N)
       IF ( IsAdvected(N) ) THEN
         CALL AADVQ( TrM(:,:,:,n), TrMom(:,:,:,:,n), .true., TrName(n) )
       ENDIF
    ENDDO

    RETURN

  END SUBROUTINE TrDYNAM

  !==========================================================================================================

  SUBROUTINE CHEM_CHUNK_RUN( nymd,       nhms,       year,       month,      &
                             day,        dayOfYr,    hour,       minute,     &
                             second,     utc,        hElapsed,   Input_Opt,  &
                             State_Chm,  State_Diag, State_Grid, State_Met,  &
                             IsChemTime, RC                                  )

    ! Based on GIGC_Chunk_Run

    ! GEOS-Chem state objects
    USE HCO_Interface_Mod,  ONLY : HcoState
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState

    ! GEOS-Chem components
    USE Chemistry_Mod,      ONLY : Do_Chemistry, Recompute_OD
    USE Convection_Mod,     ONLY : Do_Convection
    USE DryDep_Mod,         ONLY : Do_DryDep
    USE Emissions_Mod,      ONLY : Emissions_Run
    USE Mixing_Mod,         ONLY : Do_Tend, Do_Mixing
    USE WetScav_Mod,        ONLY : Setup_WetScav, Do_WetDep

    ! Specialized subroutines
    USE Calc_Met_Mod,       ONLY : AirQnt, Set_Dry_Surface_Pressure
    USE Calc_Met_Mod,       ONLY : GIGC_Cap_Tropopause_Prs
    USE Set_Global_CH4_Mod, ONLY : Set_CH4
    USE MODIS_LAI_Mod,      ONLY : Compute_XLAI
    USE PBL_Mix_Mod,        ONLY : Compute_PBL_Height
    USE Pressure_Mod,       ONLY : Set_Floating_Pressures
    USE TOMS_Mod,           ONLY : Compute_Overhead_O3
    USE UCX_Mod,            ONLY : Set_H2O_Trac

    ! Utilities
    USE ErrCode_Mod
    USE HCO_Error_Mod
    USE Pressure_Mod,       ONLY : Accept_External_Pedge
    USE State_Chm_Mod,      ONLY : IND_
    USE Time_Mod,           ONLY : Accept_External_Date_Time
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
    USE HCO_Interface_Mod,  ONLY : SetHcoTime

    ! Diagnostics
    USE Diagnostics_Mod,    ONLY : Set_Diagnostics_EndofTimestep
    USE Aerosol_Mod,        ONLY : Set_AerMass_Diagnostic
    USE CHEM_COM,           ONLY : gcaijl_out !, ijlt_vmr !, gcaijs=>gcaijs_loc
    USE DOMAIN_DECOMP_ATM,  ONLY : GRID, getDomainBounds, am_I_Root
  
    IMPLICIT NONE
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
    LOGICAL,        INTENT(IN)    :: IsChemTime  ! Time for chemistry?
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

    ! First call?
    LOGICAL, SAVE                      :: FIRST = .TRUE.

    CHARACTER(LEN=255)                 :: Iam, OrigUnit
    INTEGER                            :: STATUS, HCO_PHASE, RST

    ! Local logicals to turn on/off individual components
    ! The parts to be executed are based on the input options,
    ! the time step and the phase.
    LOGICAL                            :: DoConv
    LOGICAL                            :: DoDryDep
    LOGICAL                            :: DoEmis
    LOGICAL                            :: DoTurb
    LOGICAL                            :: DoChem
    LOGICAL                            :: DoWetDep

    ! # of times this routine has been called. Only temporary for printing
    ! processes on the first 10 calls.
    INTEGER, SAVE                      :: NCALLS = 0

    ! Strat. H2O settings
    LOGICAL                            :: SetStratH2O

    ! Whether to scale mixing ratio with meteorology update in AirQnt
    LOGICAL, SAVE                      :: scaleMR = .FALSE.

    INTEGER :: N, I, J, K, L, II, JJ

    !=======================================================================
    ! GIGC_CHUNK_RUN begins here
    !=======================================================================        

    ! By default, do processes as defined in rundeck
    DoConv   = DoGCConv                          ! dynamic time step
    DoDryDep = DOGCDryDep .AND. IsChemTime       ! chemistry time step
    DoEmis   = DoGCEmis                          ! chemistry time step
    DoTurb   = DoGCTurb                          ! dynamic time step
    DoChem   = DoGCChem .AND. IsChemTime         ! chemistry time step
    DoWetDep = DoGCWetDep                        ! dynamic time step

    ! Assume success
    RC = GC_SUCCESS

    IF ( Input_Opt%AmIRoot .and. NCALLS < 10 ) THEN
       write(6,*) 'DoConv   : ', DoConv
       write(6,*) 'DoDryDep : ', DoDryDep
       write(6,*) 'DoEmis   : ', DoEmis
       write(6,*) 'DoTurb   : ', DoTurb
       write(6,*) 'DoChem   : ', DoChem
       write(6,*) 'DoWetDep : ', DoWetDep
       write(6,*) ' '
    ENDIF

    !-------------------------------------------------------------------------
    ! Pre-Run assignments
    !-------------------------------------------------------------------------

    ! Eventually initialize/reset wetdep
    IF ( DoConv .OR. DoChem .OR. DoWetDep ) THEN
       CALL SETUP_WETSCAV( Input_Opt, State_Chm, State_Grid, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "SETUP_WETSCAV", 255 )
    ENDIF

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

    ! Set HEMCO time
    CALL SetHcoTime ( DoEmis, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "SetHcoTime", 255 )

    ! Calculate MODIS leaf area indexes needed for dry deposition
    CALL Compute_XLAI( Input_Opt, State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Compute_XLAI", 255 )

    ! Convert to dry mixing ratio
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
         'kg/kg dry', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Convert_Spc_Units", 255 )     

    ! SDE 05/28/13: Set H2O to STT if relevant
    IF ( IND_('H2O','A') > 0 ) THEN
       
       SetStratH2O = .FALSE.
       IF ( Input_Opt%LSETH2O .OR. .NOT. Input_Opt%LUCX ) THEN
          SetStratH2O = .TRUE.
       ENDIF
       CALL SET_H2O_TRAC( SetStratH2O, Input_Opt, State_Chm, &
            State_Grid,  State_Met, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "SET_H2O_TRAC", 255 )   

       ! Only force strat once if using UCX
       IF (Input_Opt%LSETH2O) Input_Opt%LSETH2O = .FALSE.

    ENDIF

    ! LTM 
    !IF ( State_Grid%NY .eq. 2 ) THEN
    !WRITE(6,'(A6,2(F8.2),8X,3X,A8,2(F8.2),8X,3X,A8,3(F8.2),8X)') "YMID: ", State_Grid%YMID(1,:), "YMID_R: ", State_Grid%YMID_R(1,:), "YEDGE", State_Grid%YEDGE(1,:)
    !ENDIF

    !IF ( State_Grid%NY .eq. 3 ) THEN
    !WRITE(6,'(A6,3(F8.2),3X,A8,3(F8.2),3X,A8,4(F8.2))') "YMID: ", State_Grid%YMID(1,:), "YMID_R: ", State_Grid%YMID_R(1,:), "YEDGE:", State_Grid%YEDGE(1,:)
    !ENDIF

    !CALL FLUSH(6)

    !=======================================================================
    ! EMISSIONS. Pass HEMCO Phase 1 which only updates the HEMCO clock
    ! and the HEMCO data list. Should be called every time to make sure
    ! that the HEMCO clock and the HEMCO data list are up to date.
    !=======================================================================
    HCO_PHASE = 1
    CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
         State_Grid, State_Met, DoEmis, HCO_PHASE, RC  )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "EMISSIONS_RUN - Phase 1", 255 )

    !=======================================================================
    ! 1. Convection
    !
    ! Call GEOS-Chem internal convection routines if convection is enabled
    ! in input.geos. This should only be done if convection is not covered
    ! by another gridded component and/or the GC species are not made
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
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "DO_DryDep", 255 )

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
       ! calculations. Note that this does not apply the emissions.
       HCO_PHASE = 2
       CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
            State_Grid, State_Met, DoEmis, HCO_PHASE, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "EMISSIONS_RUN - Phase 2", 255 )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Emissions done!'

    ENDIF

    !=======================================================================
    ! If physics covers turbulence, simply add the emission and dry
    ! deposition fluxes calculated above to the tracer array, without caring
    ! about the vertical distribution. The tracer tendencies are only added
    ! to the tracers array after emissions, drydep. So we need to use the
    ! emissions time step here.
    !=======================================================================
    ! Not implemented   

    !=======================================================================
    ! 4. Turbulence
    !
    ! Call GEOS-Chem internal turbulence routines if turbulence is enabled
    ! in input.geos. This should only be done if turbulence is not covered
    ! by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================
    IF ( DoTurb ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do turbulence now'

       ! Do mixing and apply tendencies. This will use the dynamic time step,
       ! which is fine since this call will be executed on every time step.
       !CALL DO_MIXING ( Input_Opt, State_Chm, State_Diag, &
       !     State_Grid, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "DO_MIXING", 255 )

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Turbulence done!'
    ENDIF

    ! Set tropospheric CH4 concentrations and fill species array with
    ! current values.
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM  &
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

       ! Calculate TOMS O3 overhead. For now, always use it from the
       ! Met field. State_Met%TO3 is imported from PCHEM (ckeller, 10/21/2014).
       CALL COMPUTE_OVERHEAD_O3( Input_Opt, State_Grid, State_Chm, DAY, .TRUE., &
            State_Met%TO3, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "COMPUTE_OVERHEAD_O3", 255 )

       ! Set H2O to species value if H2O is advected
       IF ( IND_('H2O','A') > 0 ) THEN
          CALL SET_H2O_TRAC( (.not. Input_Opt%LUCX), Input_Opt, &
               State_Chm, State_Grid, State_Met, RC )
          IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "SET_H2O_TRAC", 255 )
       ENDIF
       
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
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Do_WetDep", 255 )

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
    CALL Set_Diagnostics_EndofTimestep( Input_Opt,  State_Chm, State_Diag, &
         State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Set_Diagnostics_EndofTimestep", 255 )

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
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
         OrigUnit, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Convert_Spc_Units @ end of CHEM_CHUNK_RUN", 255 )    

    if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Unit conversion done!'

    !=======================================================================
    ! Accumulate diagnostics
    !=======================================================================

    !IF ( .not. FIRST ) THEN
    call getDomainBounds( grid, I_STRT = I_0, I_STOP = I_1, J_STRT = J_0, J_STOP = J_1 )
    DO N=1,NTM
    !IF ( ijlt_vmr(N) == 0 ) CALL STOP_MODEL( "Error with ijlt_vmr(N)", 255 )
    DO L=1,LM
    DO J=J_0,J_1
    DO I=I_0,I_1
       II = I - I_0 + 1
       JJ = J - J_0 + 1
       ! Archive volumetric mixing ratio in ppbv
       gcaijl_out( I, J, L, N ) = gcaijl_out( I, J, L, N ) + &
                                  State_Chm%Species(II,JJ,L,N) * 1d9 
       !IF ( N .eq. 1 .and. I .eq. 1 .and. L .eq. 1 ) WRITE(6,*) 'LTM:', n, i, j, l, gcaijl_out(i,j,l,n)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !ENDIF

    if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Diagnostic accumulation done!'

    !=======================================================================
    ! Clean up
    !=======================================================================

    ! testing only
    IF ( NCALLS < 10 ) NCALLS = NCALLS + 1

    ! First call is done
    FIRST = .FALSE.

    ! Return success
    RC = GC_SUCCESS

  END SUBROUTINE CHEM_CHUNK_RUN

  !==========================================================================================================

  SUBROUTINE INIT_CHEM( grid )

    USE DOMAIN_DECOMP_1D,        ONLY : getMpiCommunicator 
    USE DOMAIN_DECOMP_ATM,       ONLY : DIST_GRID, Am_I_Root, getDomainBounds
    USE GEOM,                    ONLY : axyp, lat2d_dg, lon2d_dg
    USE CONSTANT,                ONLY : Pi
    USE MODEL_COM,               ONLY : modelEclock, itime, ItimeI, DTsrc
    USE Dictionary_mod,          ONLY : sync_param
    USE CHEM_COM

    USE GC_Environment_Mod,      ONLY : GC_Allocate_All
    USE State_Grid_Mod,          ONLY : Init_State_Grid
    USE Input_Opt_Mod,           ONLY : Set_Input_Opt
    USE Input_Mod,               ONLY : Read_Input_File
    USE Time_Mod,                ONLY : Set_Timesteps

    USE Chemistry_Mod,           ONLY : Init_Chemistry
    USE Emissions_Mod,           ONLY : Emissions_Init
    USE GC_Environment_Mod
    USE GC_Grid_Mod,             ONLY : SetGridFromCtr
    USE PBL_Mix_Mod,             ONLY : Init_PBL_Mix
    USE Pressure_Mod,            ONLY : Init_Pressure, Accept_External_ApBp
    USE UCX_MOD,                 ONLY : Init_UCX
    USE UnitConv_Mod,            ONLY : Convert_Spc_Units
    USE PhysConstants,           ONLY : PI_180
    USE State_Chm_Mod,           ONLY : Ind_

    IMPLICIT NONE

    TYPE (DIST_GRID), INTENT(IN) :: grid

    LOGICAL   :: isRoot
    INTEGER   :: myPET, NPES, RC

    INTEGER   :: NYMD, NHMS, YEAR, MONTH, DAY, DOY, HOUR, MINUTE, SECOND
    REAL*4    :: MINUTES, hElapsed, UTC
    REAL*8    :: DT

    INTEGER   :: I, J, L, N, NN, II, JJ, I_0H, I_1H
    INTEGER   :: NYMDb, NHMSb, NYMDe, NHMSe, NSP, PLC

    ! Fill Input_Opt with input.geos
    WRITE(6,*) "Calling Set_Input_Opt"
    CALL Set_Input_Opt( am_I_Root(), Input_Opt, RC )

    WRITE(6,*) "Calling Init_State_Grid"
    CALL Init_State_Grid( Input_Opt, State_Grid, RC )

    WRITE(6,*) "Calling Read_Input_File"

    CALL Read_Input_File( Input_Opt, State_Grid, RC )

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
    
    ! LTM: Debug Overrides
    !DoGCConv   = .true.  ! Works
    !DoGCEmis   = .true.  ! Works (make sure 3-D emissions on correct vert grid)
    !DoGCTend   = .false. 
    !DoGCTurb   = .true.  ! Works
    !DoGCChem   = .true.  ! Works (make sure to turn off linear strat)
    !DoGCDryDep = .true.  ! Works
    !DoGCWetDep = .true.  ! Works
    
    !================================================================
    ! Specify local domain
    !================================================================

    call getDomainBounds( grid, I_STRT = I_0, I_STOP = I_1, &
         J_STRT = J_0, J_STOP = J_1, &
         I_STRT_HALO = I_0H, I_STOP_HALO = I_1H, &
         J_STRT_HALO = J_0H, J_STOP_HALO = J_1H )

    NI = I_1 - I_0 + 1
    NJ = J_1 - J_0 + 1

    !WRITE(6,*) i_0, i_1, NI
    !WRITE(6,*) j_0, j_1, NJ
    !WRITE(6,*) SHAPE(lon2d_dg)
    !CALL FLUSH(6)
    
    State_Grid%DX           = 2.5e+0_fp
    State_Grid%DY           = 2.0e+0_fp
    State_Grid%XMin         = lon2d_dg(i_0,j_0)
    State_Grid%XMax         = lon2d_dg(i_1,j_0)
    State_Grid%YMin         = max( lat2d_dg(i_0,j_0), -89.0_fp )
    State_Grid%YMax         = min( lat2d_dg(i_0,j_1),  89.0_fp )
    
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

    WRITE(6,*) "Calling GC_Init_Grid"
    CALL GC_Init_Grid( Input_Opt, State_Grid, RC )

    WRITE(6,*) "Calling GC_Allocate_All"
    CALL GC_Allocate_All( Input_Opt, State_Grid, RC )

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

    ! Initialize State_Met, State_Chm, and State_Diag objects
    WRITE(6,*) "Calling GC_Init_StateObj"
    CALL GC_Init_StateObj( Diag_List,  Input_Opt,  State_Chm, &
         State_Diag, State_Grid, State_Met, RC )
    
    CALL sync_param( "DTsrc", DTsrc )
    CALL sync_param( "DT",    DT    )         
    Input_Opt%TS_CHEM = INT( DTsrc  )   ! Chemistry timestep [sec]
    Input_Opt%TS_EMIS = INT( DTsrc  )   ! Chemistry timestep [sec]
    Input_Opt%TS_DYN  = INT( DT     )   ! Dynamic   timestep [sec]
    Input_Opt%TS_CONV = INT( DT     )   ! Dynamic   timestep [sec]
    Input_Opt%TS_RAD  = INT( DT     )
    Input_Opt%TS_DIAG = INT( DTsrc  )

    ! Set start and finish time from rundeck
    Input_Opt%NYMDb   = 20160701 ! nymdB
    Input_Opt%NHMSb   = 000000   ! nhmsB
    Input_Opt%NYMDe   = 20160801 ! nymdE
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
                        Diagnos    = Input_Opt%TS_DIAG         )

    ! Initialize other GEOS-Chem modules
    WRITE(6,*) "Calling GC_Init_Extra"
    CALL GC_Init_Extra( Diag_List, Input_Opt,    &
         State_Chm, State_Diag, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "GC_Init_Extra", 255 )

    ! Initialize the GEOS-Chem pressure module (set Ap & Bp)
    CALL Init_Pressure( Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Init_Pressure", 255 )

    ! Set Ap and Bp
    CALL Accept_External_ApBp( State_Grid, Ap, Bp, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Accept_External_ApBp", 255 )

    ! Initialize the PBL mixing module
    CALL INIT_PBL_MIX( Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Init_PBL_Mix", 255 )

    ! Initialize chemistry mechanism
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .OR. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL INIT_CHEMISTRY ( Input_Opt,  State_Chm, State_Diag, &
            State_Grid, RC )
       IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Init_Chemistry", 255 )
    ENDIF

    ! Initialize HEMCO
    CALL EMISSIONS_INIT( Input_Opt, State_Chm, State_Grid, State_Met, RC, &
         HcoConfig=HcoConfig )
    IF ( RC /= GC_SUCCESS ) CALL STOP_MODEL( "Emissions_Init", 255 )

    ! Stratosphere - can't be initialized without HEMCO because of STATE_PSC
    IF ( Input_Opt%LUCX ) THEN
       
       ! Initialize stratospheric routines
       CALL INIT_UCX( Input_Opt, State_Chm, State_Diag, State_Grid )
       
    ENDIF

    NSP = State_Chm%nSpecies
    NTM = State_Chm%nAdvect + 1

    ALLOCATE( TrName(NTM) )
    ALLOCATE( TrFullName(NTM) )
    ALLOCATE( IsAdvected(NTM) )    
    ALLOCATE( t_qlimit(NTM) )
    ALLOCATE(     TrM(       I_0H:I_1H, J_0H:J_1H, LM, NTM  ) )
    ALLOCATE(   TrMom( NMOM, I_0H:I_1H, J_0H:J_1H, LM, NTM ) )
    
    NN=1
    DO N = 1, NSP
       IF ( State_Chm%SpcData(N)%Info%Is_Advected .or. &
            TRIM( State_Chm%SpcData(N)%Info%Name ) .eq. "OH" ) THEN
 
       TrName(NN) = TRIM( State_Chm%SpcData(N)%Info%Name ) 
       TrFullName(NN) = TRIM( State_Chm%SpcData(N)%Info%FullName ) // " (" // &
                        TRIM( State_Chm%SpcData(N)%Info%Formula  ) // ")"
       IsAdvected(NN) = State_Chm%SpcData(N)%Info%Is_Advected
       SELECT CASE ( TrName(NN) )
         CASE('H2O')
           i_H2O   = N
         CASE('CO2')
           i_CO2   = N
         CASE('O3')
           i_O3    = N
         CASE('O2')
           i_O2    = N
         CASE('NO2')
           i_NO2   = N
         CASE('N2O')
           i_N2O   = N
         CASE('CH4')
           i_CH4   = N
         CASE('CFC11')
           i_CFC11 = N
         CASE('CFC12')
           i_CFC12 = N
         CASE('N2')
           i_N2    = N
         CASE('CFCY')
           i_CFCY  = N
         CASE('CFCZ')
           i_CFCZ  = N
         CASE('SO2')
           i_SO2   = N
         CASE DEFAULT
           plc = 0
       END SELECT

       IF ( Am_I_Root() ) WRITE(6,*) NN, N, TrName(NN)

       NN = NN + 1

       ENDIF    
    ENDDO
    t_qlimit(:) = .true.
    TrM    = 0d0
    TrMom  = 0d0

    ! Return success
    RC = GC_SUCCESS

    RETURN
  END SUBROUTINE INIT_CHEM

  !==========================================================================================================

#ifdef CACHED_SUBDD
  SUBROUTINE accumGCsubdd

    use domain_decomp_atm, only : grid
    use subdd_mod, only : subdd_groups,subdd_type,subdd_ngroups, &
         inc_subdd,find_groups, LmaxSUBDD
    use geom, only : byaxyp
    use atm_com, only : byma

    implicit none

    integer :: igrp,ngroups,grpids(subdd_ngroups)
    type(subdd_type), pointer :: subdd
    integer :: L, n, k
    real*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
         grid%j_strt_halo:grid%j_stop_halo) :: sddarr2d
    real*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
         grid%j_strt_halo:grid%j_stop_halo, &
         LM                               ) :: sddarr3d
    real*8 :: convert

    ! 3-D diagnostics on model levels
    call find_groups('taijlh',grpids,ngroups)
    do igrp=1,ngroups
       subdd => subdd_groups(grpids(igrp))
       do k=1,subdd%ndiags
          ntm_loop: do n=1,ntm
             ! tracer 3D mixing ratios (SUBDD names are just tracer name):
             if(trim(trname(n)).eq.trim(subdd%name(k))) then
                convert = 1d9 * 28.97d0 / State_Chm%SpcData(N)%Info%MW_g ! kg/kg -> ppbv
                do L=1,LmaxSUBDD
                   sddarr3d(:,:,L) = &
                        trm(:,:,L,n)*convert*byaxyp(:,:)*byma(L,:,:)
                end do
                call inc_subdd(subdd,k,sddarr3d)
                exit ntm_loop
             end if
          end do ntm_loop
       enddo ! k
    enddo ! igroup

!!$       ! 2-D I-J diags
!!$       call find_groups('taijh',grpids,ngroups)
!!$       do igrp=1,ngroups
!!$          subdd => subdd_groups(grpids(igrp))
!!$          diag_loop: do k=1,subdd%ndiags
!!$             ntm_loop3: do n=1,ntm
!!$                
!!$                ! tracer surface mixing ratios:
!!$                if(trim(trname(n))//'sm'.eq.trim(subdd%name(k))) then
!!$                   if (to_volume_MixRat(n) == 1) then
!!$                      sddarr2d(:,:)=trcsurf(:,:,n)*mass2vol(n)
!!$                   else
!!$                      sddarr2d(:,:)=trcsurf(:,:,n)
!!$                   endif
!!$                   call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop
!!$                end if
!!$                
!!$                ! tracer surface concentrations:
!!$                if(trim(trname(n))//'sc'.eq.trim(subdd%name(k))) then
!!$                   sddarr2d(:,:)=trcSurfByVol(:,:,n)
!!$                   call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop
!!$                end if
!!$                
!!$                ! tracer L=1 mixing ratios:
!!$                if(trim(trname(n))//'l1m'.eq.trim(subdd%name(k))) then
!!$                   if (to_volume_MixRat(n) == 1) then
!!$                      sddarr2d(:,:)= &
!!$                           trm(:,:,1,n)*mass2vol(n)*byaxyp(:,:)*byma(1,:,:)
!!$                   else
!!$                      sddarr2d(:,:)=trm(:,:,1,n)*byaxyp(:,:)*byma(1,:,:)
!!$                   endif
!!$                   call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop 
!!$                end if
!!$                
!!$             enddo ntm_loop3   
!!$          enddo diag_loop
!!$       enddo ! igroup

  END SUBROUTINE accumGCsubdd

#endif

  !==========================================================================================================     

  SUBROUTINE IO_CHEM( fid, action ) 

    use ParallelIo_mod
    use domain_decomp_atm, only : grid

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

END MODULE CHEM_DRV
!==========================================================================================================

#ifdef CACHED_SUBDD

SUBROUTINE tijlh_defs(arr,nmax,decl_count)
  ! Needs to be outside the module to prevent circular dependencies with SUBDD
  ! 3D tracer outputs (model horizontal grid and layers).
  use subdd_mod, only : info_type
  ! info_type_ is a homemade structure constructor for older compilers
  use subdd_mod, only : info_type_
  use chem_com, only : ntm, trname
  implicit none
  integer :: nmax,decl_count
  integer :: n
  character*80 :: unitString
  type(info_type) :: arr(nmax)

  decl_count = 0

  ! First, diagnostics available for all tracers:
  do n=1,ntm

     ! 3D mixing ratios (SUBDD string is just tracer name):
     unitString='ppbv'
     arr(next()) = info_type_(                      &
          sname = trim(trname(n)),                  &
          lname = trim(trname(n))//' mixing ratio', &
          units = trim(unitString)                  &
          )

  end do ! tracers loop

  ! Special diagnostics

  !      arr(next()) = info_type_(
  !     &  sname = 'OH_conc', ! because not a tracer
  !     &  lname = 'OH concentration',
  !     &  units = 'molecules cm-3'
  !     &  )

  return
contains
  integer function next()
    decl_count = decl_count + 1
    next = decl_count
  end function next
END SUBROUTINE tijlh_defs

#endif

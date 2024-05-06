#include "rundeck_opts.h"

#ifdef SKIP_TRACERS_RAD
#undef TRACERS_ON
#endif

!@sum RAD_DRV contains drivers for the radiation related routines
!@ver  2009/05/11
!@cont init_RAD, RADIA
!**** semi-random cloud overlap (computed opt.d+diagn)
!**** to be used with R99E or later radiation  routines.  carbon/2
!****

SUBROUTINE CALC_ZENITH_ANGLE
  !@sum calculate zenith angle for current time step
  !@auth Gavin Schmidt (from RADIA)
  USE CONSTANT,          ONLY : twopi
  USE MODEL_COM,         ONLY : itime, nday, dtsrc, calendar
  USE TIMECONSTANTS_MOD, ONLY : SECONDS_PER_DAY
  USE RAD_COM,           ONLY : cosz1
  USE RAD_COSZ0,         ONLY : COSZT
  USE TIMEINTERVAL_MOD

  IMPLICIT NONE

  INTEGER JTIME
  REAL*8 ROT1, ROT2
  TYPE (TIMEINTERVAL) :: sPerDay

  JTIME   = MOD(ITIME,NDAY)
  ROT1    = (TWOPI*JTIME)/NDAY
  sPerDay = calendar%GETSECONDSPERDAY()
  ROT2    = ROT1 + TWOPI*DTsrc/REAL(sPerDay)
  
  CALL COSZT(ROT1,ROT2,COSZ1)

END SUBROUTINE CALC_ZENITH_ANGLE

SUBROUTINE INIT_RAD( istart )
  !@sum  init_RAD initialises radiation code
  !@auth Original Development Team
  !@calls RADPAR : RCOMP1, ORBPAR
  USE FILEMANAGER
  USE RUNTIMECONTROLS_MOD, ONLY : tracers_minerals
  USE DICTIONARY_MOD
  USE CONSTANT,            ONLY : GRAV, BYSHA, TWOPI, planet_name
  USE RESOLUTION,          ONLY : jm, lm, psf
  USE ATM_COM,             ONLY : t, pk, kradia, lm_req
  USE MODEL_COM,           ONLY : DTSRC, IYEAR1, MODELECLOCK, master_yr
  USE MODEL_COM,           ONLY : orbit
  USE ATM_COM,             ONLY : pednl00
  USE DOMAIN_DECOMP_ATM,   ONLY : grid, WRITE_PARALLEL, AM_I_ROOT,        &
                                  READT_PARALLEL, GETDOMAINBOUNDS
#ifndef CUBED_SPHERE
  USE GEOM,                ONLY : lat_dg
#endif
  
  USE RADPAR,              ONLY : PTLISO, KTREND, LMR => NL, PLB,         &
                                  LS1_loc, planck_tmin, planck_tmax,      &
                                  transmission_corrections, KCLDEM,       &
                                  KSIALB, KSOLAR, SHL, snoage_fac_max,    &
                                  KZSNOW, KYEARS, KJDAYS, MADLUV,         &
                                  KYEARG, KJDAYG, MADGHG, KYEARO,         &
                                  KJDAYO, MADO3M, KYEARA, KJDAYA,         &
                                  MADAER, KYEARD, KJDAYD, MADDST,         &
                                  KYEARV, KJDAYV, MADVOL, KYEARE,         &
                                  KJDAYE, MADEPS, KYEARR, KJDAYR,         &
                                  ITR, nraero_aod => NTRACE, FS8OPX,      &
                                  FT8OPX, TRRDRY, KRHTRA, TRADEN,         &
                                  REFDRY, RCOMP1, WRITER, WRITET,         &
                                  FSTASC, FTTASC

  ! turning on options for extra aerosols
#ifdef ALTER_RADF_BY_LAT
  USE RADPAR,              ONLY : FS8OPX_orig, FT8OPX_orig
#endif

#ifdef TRACERS_SPECIAL_Shindell
  USE PHOTOLYSIS,          ONLY : aer2, miedx2, nbfastj
#endif

#ifdef TRACERS_ON
  USE RAD_COM,             ONLY : nraero_rf, nraero_seasalt, nraero_koch, &
                                  nraero_nitrate, nraero_dust,            &
																																		nraero_OMA, nraero_AMP, nraero_TOMAS
#endif
  
  USE RAD_COM,             ONLY : rqt, s0x, co2x, n2ox, ch4x, cfc11x,     &
                                  cfc12x, xGHGx, o2x, no2x, n2cx, yGHGx,  &
                                  so2x, CH4X_RADoverCHEM, snoage_def,     &
                                  s0_yr, s0_day, ghg_yr, ghg_day, volc_yr,&
                                  volc_day, aero_yr, dust_yr, O3_yr,      &
                                  H2ObyCH4, dH2O, h2ostratx, O3x, RHfix,  &
                                  CLDx, ref_mult, COSZ1, OBLIQ, ECCN,     &
                                  OMEGT, OBLIQ_DEF, ECCN_DEF, OMEGT_DEF,  &
                                  CC_cdncx, OD_cdncx, cdncl, pcdnc, vcdnc,&
                                  cloud_rad_forc, TAero_aod_diag,         &
                                  aer_rad_forc, PLB0, SHL0, albsn_yr,     &
                                  dALBsnX, nradfrc, rad_interact_aer,     &
                                  clim_interact_chem, rad_forc_lev,       &
                                  ntrix_aod, ntrix_rf, wttr,              &
                                  variable_orb_par, orb_par_year_bp,      &
                                  orb_par, nrad, RADIATIONSETORBIT,       &
                                  chl_from_obio, chl_from_seawifs

#ifdef TRACERS_ON
  USE RAD_COM,             ONLY : njaero, nraero_aod_rsf, nraero_rf_rsf,  &
		                                tau_as, tau_cs, tau_dry
#ifdef CACHED_SUBDD
  USE RAD_COM,             ONLY : abstau_as, abstau_cs, abstau_dry,       &
                                  swfrc, lwfrc
#endif
#endif

  USE RAD_COSZ0,           ONLY : COSZ_INIT
  USE CLOUDS_COM,          ONLY : llow
  USE DIAG_COM,            ONLY : IWRITE, JWRITE, ITWRITE																																		

#ifdef ALTER_RADF_BY_LAT
  USE RAD_COM,             ONLY : FULGAS_lat, FS8OPX_lat, FT8OPX_lat
#endif

#ifdef TRACERS_ON
  USE DIAG_COM, ONLY : save3dAOD
  USE TRACER_COM,          ONLY : NTM
  USE TRACER_COM,          ONLY : n_BCIA, n_BCB, n_NO3p
  USE TRACER_COM,          ONLY : n_Clay, n_Silt1, n_Silt2, n_Silt3,      &
                                  n_Silt4, n_Silt5
  USE TRACER_COM,          ONLY : n_SO4, n_Seasalt1, n_Seasalt2
  USE TRACER_COM,          ONLY : n_OCB, n_OCIA, n_Isopp1a, n_SO4
  USE TRACER_COM,          ONLY : n_vbsAm2
  USE RAD_COM,             ONLY : diag_fc
  USE TRDIAG_COM,          ONLY : save_dry_aod
#ifdef TRACERS_TOMAS       
  USE TRACER_COM,          ONLY : N_ASO4,  N_ANACL, N_AECOB, N_AECIL,     &
		                                N_AOCOB, N_AOCIL, N_ADUST
#endif
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
  USE OLDTRACER_MOD,      ONLY : TRPDENS
  USE TRDUST_MOD,         ONLY : imDust, nSubClays, dryEffRadMinerals,    &
                                 SUBCLAYWEIGHTS
  USE TRDUST_DRV,         ONLY : CALCSUBCLAYWEIGHTS
  USE TRACER_COM,         ONLY : ntm_clay, ntm_sil1, ntm_sil2, ntm_sil3,  &
                                 ntm_sil4, ntm_sil5, N_SOILDUST
  USE RAD_COM,            ONLY : nr_soildust
#endif

#ifdef TRACERS_MINERALS
  USE TRACER_COM,         ONLY : n_clayilli, n_claykaol, n_claysmec,      &
                                 n_claycalc, n_clayquar, n_clayfeld,      &
                                 n_clayhema, n_claygyps, n_clayilhe,      &
                                 n_claykahe, n_claysmhe, n_claycahe,      &
                                 n_clayquhe, n_clayfehe, n_claygyhe,      &
                                 n_sil1quar, n_sil1feld, n_sil1calc,      &
                                 n_sil1illi, n_sil1kaol, n_sil1smec,      &
                                 n_sil1hema, n_sil1gyps, n_sil1quhe,      &
                                 n_sil1fehe, n_sil1cahe, n_sil1gyhe,      &
                                 n_sil1ilhe, n_sil1kahe, n_sil1smhe,      &
                                 n_sil2quar, n_sil2feld, n_sil2calc,      &
                                 n_sil2hema, n_sil2gyps, n_sil2illi,      &
                                 n_sil2kaol, n_sil2smec, n_sil2quhe,      &
                                 n_sil2fehe, n_sil2cahe, n_sil2gyhe,      &
                                 n_sil2ilhe, n_sil2kahe, n_sil2smhe,      &
                                 n_sil3quar, n_sil3feld, n_sil3calc,      &
                                 n_sil3hema, n_sil3gyps, n_sil3illi,      &
                                 n_sil3kaol, n_sil3smec, n_sil3quhe,      &
                                 n_sil3fehe, n_sil3cahe, n_sil3gyhe,      &
                                 n_sil3ilhe, n_sil3kahe, n_sil3smhe,      &
                                 n_sil4quar, n_sil4feld, n_sil4calc,      &
                                 n_sil4hema, n_sil4gyps, n_sil4illi,      &
                                 n_sil4kaol, n_sil4smec, n_sil4quhe,      &
                                 n_sil4fehe, n_sil4cahe, n_sil4gyhe,      &
                                 n_sil4ilhe, n_sil4kahe, n_sil4smhe,      &
                                 n_sil5quar, n_sil5feld, n_sil5calc,      &
                                 n_sil5hema, n_sil5gyps, n_sil5illi,      &
                                 n_sil5kaol, n_sil5smec, n_sil5quhe,      &
                                 n_sil5fehe, n_sil5cahe, n_sil5gyhe,      &
                                 n_sil5ilhe, n_sil5kahe, n_sil5smhe
#endif

#ifdef TRACERS_AMP
  USE AERO_CONFIG,        ONLY : nmodes
  USE TRACER_COM,         ONLY : n_N_AKK_1, n_N_ACC_1, n_N_DD1_1,         &
                                 n_N_DS1_1, n_N_DD2_1, n_N_DS2_1,         &
                                 n_N_SSA_1, n_N_SSC_1, n_N_OCC_1,         &
                                 n_N_BC1_1, n_N_BC2_1, n_N_BC3_1,         &
                                 n_N_DBC_1, n_N_BOC_1, n_N_BCS_1,         &
                                 n_N_MXX_1
#endif

#ifdef TRACERS_TOMAS
  USE TOMAS_AEROSOL,      ONLY : icomp
#endif

  USE AERPARAM_MOD,       ONLY : aermix

#ifdef OLD_BCdalbsn
  USE AERPARAM_MOD,       ONLY : depoBC, depoBC_1990
#endif

  USE ABSTRACTORBIT_MOD,  ONLY : ABSTRACTORBIT

  ! begin section for radiation-only SCM
  USE CONSTANT,           ONLY : gasc, tf, mair, mwat, pi, lhe, lhs,      &
                                 mb2kg, kg2mb, kapa
  USE ATM_COM,            ONLY : q, p, PMID, pedn, PDSIG, pek, MA, BYMA,  &
                                 ltropo
  USE ATM_COM,            ONLY : AML00, BYAML00, req_fac, kradia, lm_req
  USE RESOLUTION,         ONLY : im, plbot, ls1 => LS1_NOMINAL
  USE RESOLUTION,         ONLY : MFIX, MFRAC

#ifndef STDHYB
  USE RESOLUTION,         ONLY : mfixs, mtop
#endif

  USE RAD_COM,            ONLY : modrd
  USE RADPAR,             ONLY : u0gas, ulgas, set_gases_internally
  USE RADPAR,             ONLY : set_aerosols_internally, sraext, srasct, &
                                 sragcb, srdext, srdsct, srdgcb, srvext,  &
                                 srvsct, srvgcb, srbext, srbsct, srbgcb,  &
                                 traalk, trdalk, trvalk, trbalk
  USE RADPAR,             ONLY : keepal, srbalb, srxalb, FSTOPX, FTTOPX
  USE RADPAR,             ONLY : skip_AOD_in_rad
  USE PARIO,              ONLY : PAR_OPEN, PAR_CLOSE, READ_DATA, READ_DIST_DATA
  USE FLUXES,             ONLY : atmsrf, ASFLX4, focean, fland, flice
  USE FLUXES,             ONLY : atmocn, atmice, atmgla, atmlnd
  USE GHY_COM,            ONLY : fearth
  USE LAKES_COM,          ONLY : flake
  USE SEAICE_COM,         ONLY : si_atm
  USE CLOUDS_COM,         ONLY : SVLHX, SVLAT, RHSAV
  ! end section for radiation-only SCM
		
#ifdef GCAP
  USE RAD_COM,            ONLY : SAVE_COSZ2
#endif

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: istart

  INTEGER L, LR, n1, n, nn, iu2
  REAL*8 PLBx(LM+1), pyear

  !@var NRFUN indices of unit numbers for radiation routines
  INTEGER NRFUN(14), IU, DONOTREAD

  !@var RUNSTR names of files for radiation routines
  CHARACTER*5  ::  RUNSTR(14) = &
                        (/ "RADN1", "RADN2", "RADN3", "RADN4", "RADN5",  &
                           "RADN6", "RADN7", "RADN8", "RADN9", "RADNA",  &
                        			"RADNB", "RADNC", "RADND", "RADNE"           /)

  !@var QBIN true if files for radiation input files are binary
  LOGICAL      ::  QBIN(14) = &
                        (/  .TRUE.,  .TRUE., .FALSE.,  .TRUE.,  .TRUE.,  &
                            .TRUE.,  .TRUE.,  .TRUE., .FALSE.,  .TRUE.,  &
                        				.TRUE.,  .TRUE.,  .TRUE.,  .TRUE.           /)

#ifdef TRACERS_MINERALS
  REAL(KIND=8)  ::  densclay(4*ntm_clay), denssil1(ntm_sil1),            &
                      denssil2(ntm_sil2), denssil3(ntm_sil3),            &
                      denssil4(ntm_sil4), denssil5(ntm_sil5)
#endif

  CHARACTER(LEN=300)  ::  out_line
  CHARACTER*6         ::  skip

  ! begin section for radiation-only SCM
  REAL*8            ::  cosz_const, mvar
  CHARACTER(LEN=6)  ::  gasnames(13)
  INTEGER           ::  fid, igas
  REAL*8            ::  szadeg, s0cosz, s0_tmp, cosz_tmp, tloc
  INTEGER           ::  rad_scm_int
  LOGICAL           ::  rad_scm = .FALSE.

  REAL*8 QSAT ! external function
  ! end section for radiation-only SCM

  INTEGER           ::  I, J
  INTEGER           ::  I_0, I_1, J_0, J_1
  INTEGER           ::  I_0H, I_1H
  INTEGER           ::  J_0H, J_1H

  !**** sync radiation parameters from input
  CALL SYNC_PARAM("NRAD",NRAD)

  IF ( IS_SET_PARAM("variable_orb_par") ) THEN
     CALL GET_PARAM("variable_orb_par",variable_orb_par)
  ELSEIF ( master_yr==0 ) THEN
     variable_orb_par = 1
  ELSE
     variable_orb_par = 0
  ENDIF

  IF ( IS_SET_PARAM("orb_par_year_bp") ) THEN
     CALL GET_PARAM("orb_par_year_bp",orb_par_year_bp)
  ELSEIF ( master_yr==0 ) THEN
     orb_par_year_bp = 0
  ELSE
     orb_par_year_bp = 1950 - master_yr
  ENDIF

  CALL SYNC_PARAM("orb_par",orb_par,3)
  CALL SYNC_PARAM("S0X",S0X)
  CALL SYNC_PARAM("CO2X",CO2X)              ! fulgas(2)
  CALL SYNC_PARAM("O2X",O2X)                ! fulgas(4)
  CALL SYNC_PARAM("NO2X",NO2X)              ! fulgas(5)
  CALL SYNC_PARAM("N2OX",N2OX)              ! fulgas(6)
  CALL SYNC_PARAM("CH4X",CH4X)              ! fulgas(7)
  CALL SYNC_PARAM("CH4X_RADoverCHEM",CH4X_RADoverCHEM)
  CALL SYNC_PARAM("CFC11X",CFC11X)          ! fulgas(8)
  CALL SYNC_PARAM("CFC12X",CFC12X)          ! fulgas(9)
  CALL SYNC_PARAM("N2CX",N2CX)              ! fulgas(10)
  CALL SYNC_PARAM("XGHGX",XGHGX)            ! fulgas(11)
  CALL SYNC_PARAM("YGHGX",YGHGX)            ! fulgas(12)
  CALL SYNC_PARAM("SO2X",SO2X)              ! fulgas(13)
  CALL SYNC_PARAM("H2OstratX",H2OstratX)    ! fulgas(1)
  CALL SYNC_PARAM("O3X",O3X)                ! fulgas(3)
  CALL SYNC_PARAM("CLDX",CLDX)
  CALL SYNC_PARAM("H2ObyCH4",H2ObyCH4)
  CALL GET_PARAM("S0_yr",S0_yr,DEFAULT=master_yr)

  IF ( IS_SET_PARAM("S0_day") ) THEN
     CALL GET_PARAM("S0_day",S0_day)
  ELSE
     IF ( s0_yr==0 ) s0_day = 0
     ! else use default value
  ENDIF

  CALL GET_PARAM("ghg_yr",ghg_yr,DEFAULT=master_yr)
  IF ( IS_SET_PARAM("ghg_day") ) THEN
     CALL GET_PARAM("ghg_day",ghg_day)
  ELSE
     IF ( ghg_yr==0 ) ghg_day = 0
     ! else use default value
  ENDIF

  CALL GET_PARAM("volc_yr",volc_yr,DEFAULT=master_yr)
  IF ( IS_SET_PARAM("volc_day") ) THEN
     CALL GET_PARAM("volc_day",volc_day)
  ELSE
     IF ( volc_yr==0 ) volc_day = 0
     ! else use default value
  ENDIF

  CALL  GET_PARAM("aero_yr",aero_yr,DEFAULT=master_yr)
  CALL  GET_PARAM("dust_yr",dust_yr,DEFAULT=master_yr)
  CALL SYNC_PARAM("dALBsnX",dALBsnX)
  CALL  GET_PARAM("albsn_yr",albsn_yr,DEFAULT=master_yr)
  CALL SYNC_PARAM("aermix",aermix,13)
  CALL SYNC_PARAM("REFdry",REFdry,8)
  CALL SYNC_PARAM("FS8OPX",FS8OPX,8)
  CALL SYNC_PARAM("FT8OPX",FT8OPX,8)
  CALL SYNC_PARAM("RHfix",RHfix)
  CALL SYNC_PARAM("CC_cdncx",CC_cdncx)
  CALL SYNC_PARAM("OD_cdncx",OD_cdncx)
  CALL  GET_PARAM("O3_yr",O3_yr,DEFAULT=master_yr)

  IF ( planet_name/='Earth' ) PTLISO = .015D0*psf
  ! reasonable default

  CALL SYNC_PARAM("PTLISO",PTLISO)
  CALL SYNC_PARAM("KSOLAR",KSOLAR)
  CALL SYNC_PARAM("KSIALB",KSIALB)
  CALL SYNC_PARAM("KZSNOW",KZSNOW)
  CALL SYNC_PARAM("snoage_def",snoage_def)
  CALL SYNC_PARAM("snoage_fac_max",snoage_fac_max)
  CALL SYNC_PARAM("nradfrc",nradfrc)
  IF ( snoage_fac_max<0. .OR. snoage_fac_max>1. ) THEN
     WRITE (out_line,*) 'set 0<snoage_fac_max<1, not',              &
          snoage_fac_max
     CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
     CALL STOP_MODEL('init_RAD :  snoage_fac_max out of range',255)
  ENDIF

  CALL SYNC_PARAM("rad_interact_aer",rad_interact_aer)
  CALL SYNC_PARAM("clim_interact_chem",clim_interact_chem)
  CALL SYNC_PARAM("rad_forc_lev",rad_forc_lev)
  CALL SYNC_PARAM("cloud_rad_forc",cloud_rad_forc)
  CALL SYNC_PARAM("TAero_aod_diag",TAero_aod_diag)
  CALL SYNC_PARAM("aer_rad_forc",aer_rad_forc)
  CALL SYNC_PARAM("ref_mult",ref_mult)

#ifdef TRACERS_ON
  CALL SYNC_PARAM("save3dAOD",save3dAOD)
  CALL SYNC_PARAM("save_dry_aod",save_dry_aod)
#endif

  REFdry = REFdry*ref_mult

  IF ( IS_SET_PARAM('planck_tmin') )                                &
       CALL GET_PARAM('planck_tmin',planck_tmin)
  IF ( IS_SET_PARAM('planck_tmax') )                                &
       CALL GET_PARAM('planck_tmax',planck_tmax)

  CALL GETDOMAINBOUNDS(grid,I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,       &
       J_STOP=J_1)
  CALL GETDOMAINBOUNDS(grid,J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)

  I_0H = grid%I_STRT_HALO
  I_1H = grid%I_STOP_HALO

  CALL RADIATIONSETORBIT(orbit)

  IF ( IS_SET_PARAM('rad_scm') ) THEN
     CALL GET_PARAM('rad_scm',rad_scm_int)
     rad_scm = im==1 .AND. jm==1 .AND. rad_scm_int==1
  ENDIF

  IF ( rad_scm ) THEN
     ! Some of the initializations in this block will be
     ! moved elsewhere once refactorings in the rest of
     ! the GCM are complete.

     modrd = 0
     ! just in case
     i = 1
     j = 1

     IF ( FILE_EXISTS('TOPO') ) THEN
        fid = PAR_OPEN(grid,'TOPO','read')
        CALL READ_DIST_DATA(grid,fid,'focean',focean)
        CALL READ_DIST_DATA(grid,fid,'fgice',flice)
        CALL PAR_CLOSE(grid,fid)
     ELSE
        CALL GET_PARAM('focean',focean(1,1))
        CALL GET_PARAM('flice',flice(1,1))
        CALL GET_PARAM('rsi',si_atm%RSI(1,1))
     ENDIF

     fland  = 1D0 - focean
     fearth = 1D0 - focean - flice
     flake  = 0.

     pednl00(1:LM+1) = plbot
     ! default
     pednl00(lm+2 : lm+lm_req) = req_fac(1:LM_req-1)*pednl00(lm+1)
     pednl00(lm+lm_req+1) = 0.

     IF ( FILE_EXISTS('TEMP1D') ) THEN
        fid = PAR_OPEN(grid,'TEMP1D','read')
        CALL READ_DATA(grid,fid,'t',t(i,j, : ))
        CALL PAR_CLOSE(grid,fid)
     ELSEIF ( FILE_EXISTS('AIC') ) THEN
        fid = PAR_OPEN(grid,'AIC','read')
        CALL READ_DIST_DATA(grid,fid,'t',t)
        CALL READ_DIST_DATA(grid,fid,'q',q)
        CALL READ_DIST_DATA(grid,fid,'p',p)
        !  surface pressure (mb)
        CALL READ_DIST_DATA(grid,fid,'tsurf',atmsrf%GTEMPR)
        !****     Rescaling :  pednl00(1) = ps (= p)
#ifdef STDHYB
        mvar = p(1,1)*mb2kg
#else
        mvar = p(1,1)*mb2kg - mfixs - mtop
#endif
        pednl00(lm+1) = plbot(lm+1)
        !mtop*kg2mb
        DO l = lm, 1, -1
           pednl00(l) = pednl00(l+1) + (MFIX(l)+mvar*MFRAC(l))*kg2mb
        ENDDO
        CALL PAR_CLOSE(grid,fid)
     ELSEIF ( IS_SET_PARAM('temp1d') ) THEN
        CALL GET_PARAM('temp1d',t(i,j, : ),lm)
     ENDIF

     DO l = 1, lm + lm_req
        AML00(l) = mb2kg*(pednl00(l)-pednl00(l+1))
        BYAML00(l) = 1D0/AML00(l)
     ENDDO

     pedn(1:LM+1,i,j) = pednl00(1:LM+1)

     DO l = 1, lm
        PMID(l,i,j) = .5D0*(pedn(l,i,j)+pedn(l+1,i,j))
        PDSIG(l,i,j) = (pedn(l,i,j)-pedn(l+1,i,j))
        MA(l,i,j) = PDSIG(l,i,j)*mb2kg
        BYMA(l,i,j) = 1D0/MA(l,i,j)
     ENDDO

     pk = PMID**kapa
     pek = pedn**kapa
     p(i,j) = pedn(1,i,j) - pedn(ls1,i,j)

     DO l = 1, lm
        t(i,j,l) = t(i,j,l)/pk(l,i,j)
     ENDDO

     IF ( FILE_EXISTS('GTEMPR') ) THEN
        fid = PAR_OPEN(grid,'GTEMPR','read')
        CALL READ_DATA(grid,fid,'gtempr',atmsrf%GTEMPR)
        CALL PAR_CLOSE(grid,fid)
     ELSEIF ( IS_SET_PARAM('gtempr') ) THEN
        CALL GET_PARAM('gtempr',atmsrf%GTEMPR(1,1))
     ENDIF

     CALL GET_PARAM('wsavg',atmsrf%WSAVG(1,1),DEFAULT=7D0)
     CALL GET_PARAM('bare_soil_wetness',                            &
          atmlnd%BARE_SOIL_WETNESS(1,1),DEFAULT=1D0)
     CALL GET_PARAM('snow',atmsrf%SNOW(1,1),DEFAULT=0D0)
#ifdef RESET_SURFACE_FRACTIONS_ON_SOUTH_POLE
     CALL GET_PARAM('vdata',vdata(1,1, : ),SIZE(vdata,3),             &
          DEFAULT=(/0D0,0D0,1D0,0D0,0D0,0D0,0D0,0D0,0D0,  &
          0D0,0D0,0D0/))
#endif

     IF ( FILE_EXISTS('SUN') ) THEN
        fid = PAR_OPEN(grid,'SUN','read')
        CALL READ_DATA(grid,fid,'szadeg',szadeg)
        CALL READ_DATA(grid,fid,'s0cosz',s0cosz)
        CALL PAR_CLOSE(grid,fid)
        cosz_tmp = COS(szadeg*pi/180D0)
        s0_tmp = s0cosz/cosz_tmp
        CALL SYNC_PARAM('cosz',cosz_tmp)
        CALL SYNC_PARAM('s0',s0_tmp)
     ENDIF

     LTROPO = 1

     atmsrf%TSAVG = atmsrf%GTEMPR
     atmocn%GTEMPR = atmsrf%GTEMPR
     atmice%GTEMPR = atmsrf%GTEMPR
     atmgla%GTEMPR = atmsrf%GTEMPR
     atmlnd%GTEMPR = atmsrf%GTEMPR
     DO n = 1, 4
        ASFLX4(n)%GTEMPR = atmsrf%GTEMPR
     ENDDO
     atmgla%SNOW = atmsrf%SNOW

  ENDIF

  IF ( IS_SET_PARAM('srxalb') ) THEN
     keepal = 1
     CALL GET_PARAM('srxalb',srxalb,SIZE(srxalb))
     CALL GET_PARAM('srbalb',srbalb,SIZE(srbalb))
  ENDIF

  IF ( IS_SET_PARAM('cosz') ) THEN
     CALL GET_PARAM('cosz',cosz_const)
     CALL COSZ_INIT(cosz_const=cosz_const)
  ELSE
     CALL COSZ_INIT
  ENDIF

  CALL SYNC_PARAM("chl_from_obio",chl_from_obio)
  CALL SYNC_PARAM("chl_from_seawifs",chl_from_seawifs)
  IF ( chl_from_obio>0 .AND. chl_from_seawifs>0 )                   &
       CALL STOP_MODEL("Make your mind which chl to use",255)

  IF ( istart == 2 ) THEN
     ! replace with cold vs warm start logic
     !**** SET RADIATION EQUILIBRIUM TEMPERATURES FROM LAYER LM TEMPERATURE
     DO J = J_0, J_1
        DO I = I_0, I_1
           RQT( : ,I,J) = T(I,J,LM)*PK(LM,I,J)
        ENDDO
     ENDDO
  ENDIF


  !****
  !**** SET THE CONTROL PARAMETERS FOR THE RADIATION (need mean pressures)
  !****
  LMR = LM + LM_REQ
  PLB(1:LMR+1) = PEDNL00(1:LMR+1)
  DO L = 1, LM
     PLBx(L) = PLB(L)        ! needed for CH4 prod. H2O
  ENDDO
  PLBx(LM+1) = 0.
  DO LR = LM + 1, LMR
     PLB0(LR-LM) = PLB(LR+1)
  ENDDO
  cdncl = 0
  CALL RETERP(vcdnc,pcdnc,7,cdncl,plb,llow+2)

  KTREND = 1 !  GHgas trends are determined by input file
  !note KTREND=0 is a possible but virtually obsolete option
  !****
  !             Model Add-on Data of Extended Climatology Enable Parameter
  !     MADO3M  = -1   Reads                      Ozone data the GCM way
  !     MADAER  =  1   Reads   Tropospheric Aerosol climatology 1850-2050
  !     MADAER  =  3   uses Koch,Bauer 2008 aerosol climatology 1890-2000
  !     MADDST  =  1   Reads   Dust-windblown mineral climatology   RFILE6
  !     MADVOL  =  1   Reads   Volcanic 1950-00 aerosol climatology RFILE7
  !     MADEPS  =  1   Reads   Epsilon cloud heterogeniety data     RFILE8
  !     MADLUV  =  1   Reads   Lean''s SolarUV 1882-1998 variability RFILE9
  !**** Radiative forcings are either constant = obs.value at given yr/day
  !****    or time dependent (year=0); if day=0 an annual cycle is used
  !****                                         even if the year is fixed
  KYEARS = s0_yr
  KJDAYS = s0_day
  MADLUV = 1                                   ! solar 'constant'
  KYEARG = ghg_yr
  KJDAYG = ghg_day                             ! well-mixed GHGases

#ifndef ALTER_RADF_BY_LAT
  IF ( ghg_yr>0 ) MADGHG = 0                   ! skip GHG-updating
#endif

  KYEARO = O3_yr
  KJDAYO = 0
  MADO3M = -1                                  ! ozone (ann.cycle)
  IF ( KYEARO>0 ) KYEARO = -KYEARO            ! use ONLY KYEARO-data
  KYEARA = Aero_yr
  KJDAYA = 0                ! MADAER=1 or 3, trop.aeros (ann.cycle)
  IF ( KYEARA>0 ) KYEARA = -KYEARA            ! use ONLY KYEARA-data
  IF ( FILE_EXISTS('TAero_SSA') ) MADAER = 3
  ! one of the TAero_XXX set
  KYEARD = Dust_yr
  IF ( KYEARD>0 ) KYEARD = -KYEARD            ! use ONLY KYEARD-data
  KYEARV = Volc_yr
  KJDAYV = Volc_day
  IF ( FILE_EXISTS('RADN7') ) MADVOL = 1
  ! Volc. Aerosols
  CALL SYNC_PARAM("MADVOL",MADVOL)
  !***  KYEARV=0  :  use current year
  !***  KYEARV<0  :  use long term mean stratospheric aerosols (use -1)
  !     Hack :  KYEARV= -2000 and -2010 were used for 2 specific runs that
  !           ended in 2100 and repeated some 20th century volcanos
  !***  KYEARV=-2000 :  use volcanos from 100 yrs ago after 2000
  !***  KYEARV=-2010 :  repeat 2nd half, then first half of 20th century
  IF ( KYEARV<=-2000 ) KYEARV = 0
  ! use current year (before 2000)
  !**** NO time history (yet), except for ann.cycle, for forcings below;
  !****  if KJDAY?=day0 (1->365), data from that day are used all year
  KYEARE = 0
  KJDAYE = 0
  KYEARR = 0
  KJDAYR = 0                          ! surf.reflectance (ann.cycle)
  KCLDEM = 1                ! 0 : old 1 : new LW cloud scattering scheme

  IF ( FILE_EXISTS('DUSTaer') ) MADDST = 1
  ! Desert dust
  IF ( FILE_EXISTS('RADN8') ) MADEPS = 1  ! cloud Epsln - KCLDEP
  transmission_corrections = FILE_EXISTS('RADN4')

  !**** Aerosols : 
  !**** Currently there are five different default aerosol controls
  !****   1 : total 2 : background+tracer 3 : Climatology 4 : dust 5 : volcanic
  !**** By adjusting FSXAER,FTXAER you can remove the default
  !**** aerosols and replace them with your version if required
  !**** (through TRACER in RADIA).
  !**** FSXAER is for the shortwave,    FTXAER is for the longwave effects
  !aer  FSXAER = (/ 1.,1.,1.,1.,1. /) ; FTXAER = (/ 1.,1.,1.,1.,1. /)

  !**** climatology aerosols are grouped into 6 types from 13 sources : 
  !****  Pre-Industrial+Natural 1850 Level  Industrial Process  BioMBurn
  !****  ---------------------------------  ------------------  --------
  !****   1    2    3    4    5    6    7    8    9   10   11   12   13
  !****  SNP  SBP  SSP  ANP  ONP  OBP  BBP  SUI  ANI  OCI  BCI  OCB  BCB
  !**** using the following default scaling/tuning factors  AERMIX(1-13)
  !****  1.0, 1.0, .26, 1.0, 2.5, 2.5, 1.9, 1.0, 1.0, 2.5, 1.9, 2.5, 1.9
  !**** The 8 groups are (adding dust and volcanic aerosols as 7. and 8.)
  !**** 1. Sulfates (industr and natural), 2. Sea Salt, 3. Nitrates
  !**** 4. Organic Carbons, 5. industr Black Carbons(BC), 6. Biomass BC
  !**** 7. Dust aerosols, 8. Volcanic aerosols
  !**** use FS8OPX and FT8OPX to enhance the optical effect; defaults : 
  !aer  FS8OPX = (/1., 1., 1., 1., 2., 2.,    1.   ,   1./)     solar
  !aer  FT8OPX = (/1., 1., 1., 1., 1., 1.,    1.3d0,   1./)     thermal
!!!!! Note :  FS|T8OPX(7-8) makes FS|TXAER(4-5) redundant.
  !**** Particle sizes of the first 4 groups have RelHum dependence

  !**** To add up to 8 further aerosols : 
  !****  1) set nraero_aod to the number of extra aerosol fields
  !****  2) ITR defines which set of Mie parameters get used, choose
  !****     from the following : 
  !****     1 SO4,  2 seasalt, 3 nitrate, 4 OCX organic carbons
  !****     5 BCI,  6 BCB,     7 dust,    8 H2SO4 volc
  !****  2b) set up the indexing array ntrix_aod to map the RADIATION tracers
  !****      to the main model tracers
  !****  2c) set up the weighting array WTTR to weight main model tracers,
  !****      if needed (default value is 1).
  !****
  !****  3) Use FSTOPX/FTTOPX(1 : nraero_aod) to scale them in RADIA
  !****  4) Set TRRDRY to dry radius
  !****  5) Set KRHTRA=1 if aerosol has RH dependence, 0 if not
  !**** Note :  whereas FSXAER/FTXAER are global (shared), FSTOPX/FTTOPX
  !****       have to be reset for each grid box to allow for the way it
  !****       is used in RADIA (TRACERS_AEROSOLS_Koch)
  !aer   nraero_aod = 0
  !aer   ITR = (/ 0,0,0,0, 0,0,0,0 /)
  !aer   TRRDRY=(/ .1d0, .1d0, .1d0, .1d0, .1d0, .1d0, .1d0, .1d0/)
  !aer   KRHTRA=(/1,1,1,1,1,1,1,1/)

#if defined( TRACERS_ON )

#if defined( TRACERS_AMP )
  nraero_AMP = nmodes
  IF ( diag_fc==2 ) THEN
     nraero_rf = nraero_rf + nraero_AMP
  ELSEIF ( diag_fc==1 ) THEN
     IF ( nraero_AMP>0 ) nraero_rf = nraero_rf + 1
  ENDIF
#elif defined( TRACERS_TOMAS )
  !TOMAS does not include NO3 AND VOL, which use its default radiation.
#ifndef TRACERS_NITRATE
  nraero_TOMAS = icomp - 2
#else
  nraero_TOMAS = icomp - 1
#endif
  IF ( diag_fc==2 ) THEN
     nraero_rf = nraero_rf + nraero_TOMAS
  ELSEIF ( diag_fc==1 ) THEN
     IF ( nraero_TOMAS>0 ) nraero_rf = nraero_rf + 1
  ENDIF
#else
  nraero_OMA = nraero_seasalt + nraero_koch + nraero_nitrate +      &
       nraero_dust
  IF ( diag_fc==2 ) THEN
     nraero_rf = nraero_rf + nraero_OMA
  ELSEIF ( diag_fc==1 ) THEN
     IF ( nraero_OMA>0 ) nraero_rf = nraero_rf + 1
  ENDIF
#endif

  nraero_aod = nraero_OMA + nraero_AMP + nraero_TOMAS

  IF ( nraero_aod_rsf>0 ) THEN
     IF ( nraero_aod_rsf/=nraero_aod )                              &
          CALL STOP_MODEL('nraero_aod_rsf /= nraero_aod',255)
  ENDIF

  IF ( nraero_rf_rsf>0 ) THEN
     IF ( nraero_rf_rsf/=nraero_rf )                                &
          CALL STOP_MODEL('nraero_rf_rsf /= nraero_rf',255)
  ENDIF

  IF ( nraero_aod>0 ) THEN
     ALLOCATE (ntrix_aod(nraero_aod))
     ntrix_aod = 0
     IF ( nraero_rf>0 ) ALLOCATE (ntrix_rf(nraero_rf))
     ntrix_rf = 0
     ALLOCATE (wttr(nraero_aod))
     wttr = 1.

     IF ( .NOT.ALLOCATED(tau_as) ) THEN
        ALLOCATE (tau_as(I_0H : I_1H,J_0H : J_1H,lm,nraero_aod))
        ALLOCATE (tau_cs(I_0H : I_1H,J_0H : J_1H,lm,nraero_aod))
        tau_as = 0.D0
        tau_cs = 0.D0
        IF ( save_dry_aod>0 ) THEN
           ALLOCATE (tau_dry(I_0H : I_1H,J_0H : J_1H,lm,nraero_aod))
           tau_dry = 0.D0
        ENDIF
#ifdef CACHED_SUBDD
        ALLOCATE (abstau_as(I_0H : I_1H,J_0H : J_1H,lm,nraero_aod))
        ALLOCATE (abstau_cs(I_0H : I_1H,J_0H : J_1H,lm,nraero_aod))
        abstau_as = 0.D0
        abstau_cs = 0.D0
        IF ( save_dry_aod>0 ) THEN
           ALLOCATE (abstau_dry(I_0H : I_1H,J_0H : J_1H,lm,nraero_aod))
           abstau_dry = 0.D0
        ENDIF
        IF ( nraero_rf>0 ) THEN
           ALLOCATE (swfrc(I_0H : I_1H,J_0H : J_1H,nraero_rf))
           ALLOCATE (lwfrc(I_0H : I_1H,J_0H : J_1H,nraero_rf))
           swfrc = 0.D0
           lwfrc = 0.D0
        ENDIF
#endif  /* CACHED_SUBDD */
     ENDIF
  ENDIF
#ifdef TRACERS_SPECIAL_Shindell
#if (! defined(TRACERS_AMP)) &  (! defined(TRACERS_TOMAS))
  njaero = nraero_aod + 2
#else
  njaero = 2
#endif
  ALLOCATE (miedx2(nbfastj,njaero))
  ALLOCATE (aer2(nbfastj,njaero))
#endif  /* TRACERS_SPECIAL_Shindell */

  !=======================================================================
  ! Define indices to map model aerosol tracer arrays to radiation arrays
  ! and other radiation-related aerosol properties
  !=======================================================================
  n = 0
  !-----------------------------------------------------------------------
#ifdef TRACERS_AEROSOLS_SEASALT
  IF ( nraero_seasalt>0 ) THEN
     IF ( rad_interact_aer>0 ) THEN
        FS8OPX(2) = 0.D0
        FT8OPX(2) = 0.D0
     ENDIF
     ntrix_aod(n+1 : n+nraero_seasalt) = (/n_seasalt1,n_seasalt2/)
     trrdry(n+1 : n+nraero_seasalt) = (/0.44D0,1.7D0/)
     itr(n+1 : n+nraero_seasalt) = (/2,2/)
  ENDIF
  n = n + nraero_seasalt
#endif  /* TRACERS_AEROSOLS_SEASALT */
  !-----------------------------------------------------------------------
#ifdef TRACERS_AEROSOLS_Koch
  IF ( nraero_koch>0 ) THEN
     IF ( rad_interact_aer>0 ) THEN ! if BC''s sol.effect are doubled : 
        FS8OPX(1) = 0.D0
        FT8OPX(1) = 0.D0
#ifndef SULF_ONLY_AEROSOLS
        FS8OPX(4 : 6) = 0.D0
        FT8OPX(4 : 6) = 0.D0
#endif
     ENDIF
     ntrix_aod(n+1) = n_SO4
     trrdry(n+1) = 0.15D0
     itr(n+1) = 1

#ifndef SULF_ONLY_AEROSOLS

#if defined( TRACERS_AEROSOLS_VBS ) &  defined( TRACERS_AEROSOLS_SOA )

     ntrix_aod(n+2 : n+nraero_koch) = (/n_vbsAm2,n_isopp1a,n_BCIA,    &
          n_BCB/)
     trrdry(n+2 : n+nraero_koch) = (/0.2D0,0.2D0,0.08D0,0.08D0/)
     itr(n+2 : n+nraero_koch) = (/4,5,6/)
     krhtra(n+2 : n+nraero_koch) = (/1,0,0/)
     !       Augment BC by 50 %
     fstasc(n+2 : n+nraero_koch) = (/1.D0,1.5D0,1.5D0/)

#elif defined( TRACERS_AEROSOLS_VBS )

     ntrix_aod(n+2 : n+nraero_koch) = (/n_vbsAm2,n_BCIA,n_BCB/)
     trrdry(n+2 : n+nraero_koch) = (/0.2D0,0.08D0,0.08D0/)
     itr(n+2 : n+nraero_koch) = (/4,5,6/)
     krhtra(n+2 : n+nraero_koch) = (/1,0,0/)
     !       Augment BC by 50 %
     fstasc(n+2 : n+nraero_koch) = (/1.D0,1.5D0,1.5D0/)

#elif defined( TRACERS_AEROSOLS_SOA )

     ntrix_aod(n+2 : n+nraero_koch) = (/n_OCIA,n_OCB,n_isopp1a,n_BCIA,&
          n_BCB/)
     trrdry(n+2 : n+nraero_koch) = (/0.2D0,0.2D0,0.2D0,0.08D0,0.08D0/)
     itr(n+2 : n+nraero_koch) = (/4,4,5,6/)
     krhtra(n+2 : n+nraero_koch) = (/1,1,0,0/)
     !       Augment BC by 50 %
     fstasc(n+2 : n+nraero_koch) = (/1.D0,1.D0,1.5D0,1.5D0/)

#else

     ntrix_aod(n+2 : n+nraero_koch) = (/n_OCIA,n_OCB,n_BCIA,n_BCB/)
     trrdry(n+2 : n+nraero_koch) = (/0.2D0,0.2D0,0.08D0,0.08D0/)
     itr(n+2 : n+nraero_koch) = (/4,4,5,6/)
     krhtra(n+2 : n+nraero_koch) = (/1,1,0,0/)
     !       Augment BC by 50 %
     fstasc(n+2 : n+nraero_koch) = (/1.D0,1.D0,1.5D0,1.5D0/)

#endif

#endif  /* SULF_ONLY_AEROSOLS */
  ENDIF
  n = n + nraero_koch
#endif  /* TRACERS_AEROSOLS_Koch */
  !-----------------------------------------------------------------------
#ifdef TRACERS_NITRATE
  IF ( nraero_nitrate>0 ) THEN
#ifdef SULF_ONLY_AEROSOLS
     CALL STOP_MODEL('SULF_ONLY_AEROSOLS and TRACERS_NITRATE on',   &
          255)
#endif /* OFF :  SULF_ONLY_AEROSOLS */
     IF ( rad_interact_aer>0 ) THEN
        ! turn off default nitrate
        FS8OPX(3) = 0.D0
        FT8OPX(3) = 0.D0
     ENDIF
     ntrix_aod(n+1 : n+nraero_nitrate) = (/n_NO3p/)
     trrdry(n+1 : n+nraero_nitrate) = (/0.15D0/)
     itr(n+1 : n+nraero_nitrate) = (/3/)
  ENDIF
  n = n + nraero_nitrate
#endif  /* TRACERS_NITRATE */
  !-----------------------------------------------------------------------
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
  IF ( nraero_dust>0 ) THEN
     IF ( rad_interact_aer>0 ) THEN
        ! turn off default dust
        FS8OPX(7) = 0.D0
        FT8OPX(7) = 0.D0
     ENDIF
     nr_soildust = n + 1

#ifdef TRACERS_MINERALS

#if defined( TRACERS_DUST_Silt4 ) &  defined( TRACERS_DUST_Silt5 )

     ! Adjust if number of dust tracers changes.
     ntrix_aod(n+1 : n+nraero_dust) = (/(n_clayilli,i=1,nSubClays),(  &
          n_claykaol,i=1,nSubClays),      &
          (n_claysmec,i=1,nSubClays),     &
          (n_claycalc,i=1,nSubClays),     &
          (n_clayquar,i=1,nSubClays),     &
          (n_clayfeld,i=1,nSubClays),     &
          (n_clayhema,i=1,nSubClays),     &
          (n_claygyps,i=1,nSubClays),     &
          (n_clayilhe,i=1,nSubClays),     &
          (n_claykahe,i=1,nSubClays),     &
          (n_claysmhe,i=1,nSubClays),     &
          (n_claycahe,i=1,nSubClays),     &
          (n_clayquhe,i=1,nSubClays),     &
          (n_clayfehe,i=1,nSubClays),     &
          (n_claygyhe,i=1,nSubClays),     &
          n_sil1illi,n_sil1kaol,          &
          n_sil1smec,n_sil1calc,          &
          n_sil1quar,n_sil1feld,          &
          n_sil1hema,n_sil1gyps,          &
          n_sil1ilhe,n_sil1kahe,          &
          n_sil1smhe,n_sil1cahe,          &
          n_sil1quhe,n_sil1fehe,          &
          n_sil1gyhe,n_sil2illi,          &
          n_sil2kaol,n_sil2smec,          &
          n_sil2calc,n_sil2quar,          &
          n_sil2feld,n_sil2hema,          &
          n_sil2gyps,n_sil2ilhe,          &
          n_sil2kahe,n_sil2smhe,          &
          n_sil2cahe,n_sil2quhe,          &
          n_sil2fehe,n_sil2gyhe,          &
          n_sil3illi,n_sil3kaol,          &
          n_sil3smec,n_sil3calc,          &
          n_sil3quar,n_sil3feld,          &
          n_sil3hema,n_sil3gyps,          &
          n_sil3ilhe,n_sil3kahe,          &
          n_sil3smhe,n_sil3cahe,          &
          n_sil3quhe,n_sil3fehe,          &
          n_sil3gyhe,n_sil4illi,          &
          n_sil4kaol,n_sil4smec,          &
          n_sil4calc,n_sil4quar,          &
          n_sil4feld,n_sil4hema,          &
          n_sil4gyps,n_sil4ilhe,          &
          n_sil4kahe,n_sil4smhe,          &
          n_sil4cahe,n_sil4quhe,          &
          n_sil4fehe,n_sil4gyhe,          &
          n_sil5illi,n_sil5kaol,          &
          n_sil5smec,n_sil5calc,          &
          n_sil5quar,n_sil5feld,          &
          n_sil5hema,n_sil5gyps,          &
          n_sil5ilhe,n_sil5kahe,          &
          n_sil5smhe,n_sil5cahe,          &
          n_sil5quhe,n_sil5fehe,          &
          n_sil5gyhe/)

     trrdry(n+1 : n+nraero_dust) = (/(dryEffRadMinerals(1 : nSubClays),i&
          =1,ntm_clay),                      &
          (dryEffRadMinerals(5),i=1,ntm_sil1)&
          ,                                  &
          (dryEffRadMinerals(6),i=1,ntm_sil2)&
          ,                                  &
          (dryEffRadMinerals(7),i=1,ntm_sil3)&
          ,                                  &
          (dryEffRadMinerals(8),i=1,ntm_sil4)&
          ,                                  &
          (dryEffRadMinerals(9),i=1,ntm_sil5)&
          /)


#elif defined( TRACERS_DUST_Silt5 )

     ! Adjust if number of dust tracers changes.
     ntrix_aod(n+1 : n+nraero_dust) = (/(n_clayilli,i=1,nSubClays),(  &
          n_claykaol,i=1,nSubClays),      &
          (n_claysmec,i=1,nSubClays),     &
          (n_claycalc,i=1,nSubClays),     &
          (n_clayquar,i=1,nSubClays),     &
          (n_clayfeld,i=1,nSubClays),     &
          (n_clayhema,i=1,nSubClays),     &
          (n_claygyps,i=1,nSubClays),     &
          (n_clayilhe,i=1,nSubClays),     &
          (n_claykahe,i=1,nSubClays),     &
          (n_claysmhe,i=1,nSubClays),     &
          (n_claycahe,i=1,nSubClays),     &
          (n_clayquhe,i=1,nSubClays),     &
          (n_clayfehe,i=1,nSubClays),     &
          (n_claygyhe,i=1,nSubClays),     &
          n_sil1illi,n_sil1kaol,          &
          n_sil1smec,n_sil1calc,          &
          n_sil1quar,n_sil1feld,          &
          n_sil1hema,n_sil1gyps,          &
          n_sil1ilhe,n_sil1kahe,          &
          n_sil1smhe,n_sil1cahe,          &
          n_sil1quhe,n_sil1fehe,          &
          n_sil1gyhe,n_sil2illi,          &
          n_sil2kaol,n_sil2smec,          &
          n_sil2calc,n_sil2quar,          &
          n_sil2feld,n_sil2hema,          &
          n_sil2gyps,n_sil2ilhe,          &
          n_sil2kahe,n_sil2smhe,          &
          n_sil2cahe,n_sil2quhe,          &
          n_sil2fehe,n_sil2gyhe,          &
          n_sil3illi,n_sil3kaol,          &
          n_sil3smec,n_sil3calc,          &
          n_sil3quar,n_sil3feld,          &
          n_sil3hema,n_sil3gyps,          &
          n_sil3ilhe,n_sil3kahe,          &
          n_sil3smhe,n_sil3cahe,          &
          n_sil3quhe,n_sil3fehe,          &
          n_sil3gyhe,n_sil5illi,          &
          n_sil5kaol,n_sil5smec,          &
          n_sil5calc,n_sil5quar,          &
          n_sil5feld,n_sil5hema,          &
          n_sil5gyps,n_sil5ilhe,          &
          n_sil5kahe,n_sil5smhe,          &
          n_sil5cahe,n_sil5quhe,          &
          n_sil5fehe,n_sil5gyhe/)

     trrdry(n+1 : n+nraero_dust) = (/(dryEffRadMinerals(1 : nSubClays),i&
          =1,ntm_clay),                      &
          (dryEffRadMinerals(5),i=1,ntm_sil1)&
          ,                                  &
          (dryEffRadMinerals(6),i=1,ntm_sil2)&
          ,                                  &
          (dryEffRadMinerals(7),i=1,ntm_sil3)&
          ,                                  &
          (dryEffRadMinerals(9),i=1,ntm_sil5)&
          /)

     !
#elif defined( TRACERS_DUST_Silt4 )

     ! Adjust if number of dust tracers changes.
     ntrix_aod(n+1 : n+nraero_dust) = (/(n_clayilli,i=1,nSubClays),(  &
          n_claykaol,i=1,nSubClays),      &
          (n_claysmec,i=1,nSubClays),     &
          (n_claycalc,i=1,nSubClays),     &
          (n_clayquar,i=1,nSubClays),     &
          (n_clayfeld,i=1,nSubClays),     &
          (n_clayhema,i=1,nSubClays),     &
          (n_claygyps,i=1,nSubClays),     &
          (n_clayilhe,i=1,nSubClays),     &
          (n_claykahe,i=1,nSubClays),     &
          (n_claysmhe,i=1,nSubClays),     &
          (n_claycahe,i=1,nSubClays),     &
          (n_clayquhe,i=1,nSubClays),     &
          (n_clayfehe,i=1,nSubClays),     &
          (n_claygyhe,i=1,nSubClays),     &
          n_sil1illi,n_sil1kaol,          &
          n_sil1smec,n_sil1calc,          &
          n_sil1quar,n_sil1feld,          &
          n_sil1hema,n_sil1gyps,          &
          n_sil1ilhe,n_sil1kahe,          &
          n_sil1smhe,n_sil1cahe,          &
          n_sil1quhe,n_sil1fehe,          &
          n_sil1gyhe,n_sil2illi,          &
          n_sil2kaol,n_sil2smec,          &
          n_sil2calc,n_sil2quar,          &
          n_sil2feld,n_sil2hema,          &
          n_sil2gyps,n_sil2ilhe,          &
          n_sil2kahe,n_sil2smhe,          &
          n_sil2cahe,n_sil2quhe,          &
          n_sil2fehe,n_sil2gyhe,          &
          n_sil3illi,n_sil3kaol,          &
          n_sil3smec,n_sil3calc,          &
          n_sil3quar,n_sil3feld,          &
          n_sil3hema,n_sil3gyps,          &
          n_sil3ilhe,n_sil3kahe,          &
          n_sil3smhe,n_sil3cahe,          &
          n_sil3quhe,n_sil3fehe,          &
          n_sil3gyhe,n_sil4illi,          &
          n_sil4kaol,n_sil4smec,          &
          n_sil4calc,n_sil4quar,          &
          n_sil4feld,n_sil4hema,          &
          n_sil4gyps,n_sil4ilhe,          &
          n_sil4kahe,n_sil4smhe,          &
          n_sil4cahe,n_sil4quhe,          &
          n_sil4fehe,n_sil4gyhe/)

     trrdry(n+1 : n+nraero_dust) = (/(dryEffRadMinerals(1 : nSubClays),i&
          =1,ntm_clay),                      &
          (dryEffRadMinerals(5),i=1,ntm_sil1)&
          ,                                  &
          (dryEffRadMinerals(6),i=1,ntm_sil2)&
          ,                                  &
          (dryEffRadMinerals(7),i=1,ntm_sil3)&
          ,                                  &
          (dryEffRadMinerals(8),i=1,ntm_sil4)&
          /)
     !
#else

     ! Adjust if number of dust tracers changes.
     ntrix_aod(n+1 : n+nraero_dust) = (/(n_clayilli,i=1,nSubClays),(  &
          n_claykaol,i=1,nSubClays),      &
          (n_claysmec,i=1,nSubClays),     &
          (n_claycalc,i=1,nSubClays),     &
          (n_clayquar,i=1,nSubClays),     &
          (n_clayfeld,i=1,nSubClays),     &
          (n_clayhema,i=1,nSubClays),     &
          (n_claygyps,i=1,nSubClays),     &
          (n_clayilhe,i=1,nSubClays),     &
          (n_claykahe,i=1,nSubClays),     &
          (n_claysmhe,i=1,nSubClays),     &
          (n_claycahe,i=1,nSubClays),     &
          (n_clayquhe,i=1,nSubClays),     &
          (n_clayfehe,i=1,nSubClays),     &
          (n_claygyhe,i=1,nSubClays),     &
          n_sil1illi,n_sil1kaol,          &
          n_sil1smec,n_sil1calc,          &
          n_sil1quar,n_sil1feld,          &
          n_sil1hema,n_sil1gyps,          &
          n_sil1ilhe,n_sil1kahe,          &
          n_sil1smhe,n_sil1cahe,          &
          n_sil1quhe,n_sil1fehe,          &
          n_sil1gyhe,n_sil2illi,          &
          n_sil2kaol,n_sil2smec,          &
          n_sil2calc,n_sil2quar,          &
          n_sil2feld,n_sil2hema,          &
          n_sil2gyps,n_sil2ilhe,          &
          n_sil2kahe,n_sil2smhe,          &
          n_sil2cahe,n_sil2quhe,          &
          n_sil2fehe,n_sil2gyhe,          &
          n_sil3illi,n_sil3kaol,          &
          n_sil3smec,n_sil3calc,          &
          n_sil3quar,n_sil3feld,          &
          n_sil3hema,n_sil3gyps,          &
          n_sil3ilhe,n_sil3kahe,          &
          n_sil3smhe,n_sil3cahe,          &
          n_sil3quhe,n_sil3fehe,          &
          n_sil3gyhe/)

     trrdry(n+1 : n+nraero_dust) = (/(dryEffRadMinerals(1 : nSubClays),i&
          =1,ntm_clay),                      &
          (dryEffRadMinerals(5),i=1,ntm_sil1)&
          ,                                  &
          (dryEffRadMinerals(6),i=1,ntm_sil2)&
          ,                                  &
          (dryEffRadMinerals(7),i=1,ntm_sil3)&
          /)

#endif

     IF ( tracers_minerals ) CALL CALCSUBCLAYWEIGHTS

     wttr(n+1 : n+nraero_dust) = (/((SUBCLAYWEIGHTS(i,j),i=1,nSubClays&
          ),j=1,ntm_clay),                     &
          (1.D0,i=1,ntm_sil1+ntm_sil2+ntm_sil3+&
          ntm_sil4+ntm_sil5)/)

     densclay = (/(TRPDENS(n_clayilli),i=1,nSubClays),              &
          (TRPDENS(n_claykaol),i=1,nSubClays),                &
          (TRPDENS(n_claysmec),i=1,nSubClays),                &
          (TRPDENS(n_claycalc),i=1,nSubClays),                &
          (TRPDENS(n_clayquar),i=1,nSubClays),                &
          (TRPDENS(n_clayfeld),i=1,nSubClays),                &
          (TRPDENS(n_clayhema),i=1,nSubClays),                &
          (TRPDENS(n_claygyps),i=1,nSubClays),                &
          (TRPDENS(n_clayilhe),i=1,nSubClays),                &
          (TRPDENS(n_claykahe),i=1,nSubClays),                &
          (TRPDENS(n_claysmhe),i=1,nSubClays),                &
          (TRPDENS(n_claycahe),i=1,nSubClays),                &
          (TRPDENS(n_clayquhe),i=1,nSubClays),                &
          (TRPDENS(n_clayfehe),i=1,nSubClays),                &
          (TRPDENS(n_claygyhe),i=1,nSubClays)/)
     denssil1 = (/TRPDENS(n_sil1illi),TRPDENS(n_sil1kaol),          &
          TRPDENS(n_sil1smec),TRPDENS(n_sil1calc),            &
          TRPDENS(n_sil1quar),TRPDENS(n_sil1feld),            &
          TRPDENS(n_sil1hema),TRPDENS(n_sil1gyps),            &
          TRPDENS(n_sil1ilhe),TRPDENS(n_sil1kahe),            &
          TRPDENS(n_sil1smhe),TRPDENS(n_sil1cahe),            &
          TRPDENS(n_sil1quhe),TRPDENS(n_sil1fehe),            &
          TRPDENS(n_sil1gyhe)/)
     denssil2 = denssil1
     denssil3 = denssil1
#ifdef TRACERS_DUST_Silt4
     denssil4 = denssil1
#endif  /* TRACERS_DUST_Silt4 */
#ifdef TRACERS_DUST_Silt5
     denssil5 = denssil1
#endif  /* TRACERS_DUST_Silt5 */


#if defined( TRACERS_DUST_Silt4 ) &  defined( TRACERS_DUST_Silt5 )

     traden(n+1 : n+nraero_dust) = (/densclay( : ),denssil1( : ),denssil2(&
          : ),denssil3( : ),denssil4( : ),        &
          denssil5( : )/)*1D-3
     ! Convert from kg/m^3 to g/cm^3

#elif defined( TRACERS_DUST_Silt5 )

     traden(n+1 : n+nraero_dust) = (/densclay( : ),denssil1( : ),denssil2(&
          : ),denssil3( : ),denssil5( : )/)*1D-3
     ! Convert from kg/m^3 to g/cm^3

#elif defined( TRACERS_DUST_Silt4 )

     traden(n+1 : n+nraero_dust) = (/densclay( : ),denssil1( : ),denssil2(&
          : ),denssil3( : ),denssil4( : )/)*1D-3
     ! Convert from kg/m^3 to g/cm^3

#else
     !
     traden(n+1 : n+nraero_dust) = (/densclay( : ),denssil1( : ),denssil2(&
          : ),denssil3( : )/)*1D-3
     ! Convert from kg/m^3 to g/cm^3

#endif

#else  /* not TRACERS_MINERALS */


#if defined( TRACERS_DUST_Silt4 ) &  defined( TRACERS_DUST_Silt5 )

     ntrix_aod(n+1 : n+nraero_dust) = (/(n_clay,i=1,nSubClays),n_silt1&
          ,n_silt2,n_silt3,n_silt4,       &
          n_silt5/)

     trrdry(n+1 : n+nraero_dust) = (/(dryEffRadMinerals(1 : nSubClays),i&
          =1,ntm_clay),                      &
          (dryEffRadMinerals(5),i=1,ntm_sil1)&
          ,                                  &
          (dryEffRadMinerals(6),i=1,ntm_sil2)&
          ,                                  &
          (dryEffRadMinerals(7),i=1,ntm_sil3)&
          ,                                  &
          (dryEffRadMinerals(8),i=1,ntm_sil4)&
          ,                                  &
          (dryEffRadMinerals(9),i=1,ntm_sil5)&
          /)

#elif defined( TRACERS_DUST_Silt5 )

     ntrix_aod(n+1 : n+nraero_dust) = (/(n_clay,i=1,nSubClays),n_silt1&
          ,n_silt2,n_silt3,n_silt5/)

     trrdry(n+1 : n+nraero_dust) = (/(dryEffRadMinerals(1 : nSubClays),i&
          =1,ntm_clay),                      &
          (dryEffRadMinerals(5),i=1,ntm_sil1)&
          ,                                  &
          (dryEffRadMinerals(6),i=1,ntm_sil2)&
          ,                                  &
          (dryEffRadMinerals(7),i=1,ntm_sil3)&
          ,                                  &
          (dryEffRadMinerals(9),i=1,ntm_sil5)&
          /)
     !
#elif defined( TRACERS_DUST_Silt4 )

     ntrix_aod(n+1 : n+nraero_dust) = (/(n_clay,i=1,nSubClays),n_silt1&
          ,n_silt2,n_silt3,n_silt4/)

     trrdry(n+1 : n+nraero_dust) = (/(dryEffRadMinerals(1 : nSubClays),i&
          =1,ntm_clay),                      &
          (dryEffRadMinerals(5),i=1,ntm_sil1)&
          ,                                  &
          (dryEffRadMinerals(6),i=1,ntm_sil2)&
          ,                                  &
          (dryEffRadMinerals(7),i=1,ntm_sil3)&
          ,                                  &
          (dryEffRadMinerals(8),i=1,ntm_sil4)&
          /)
     !
#else

     ntrix_aod(n+1 : n+nraero_dust) = (/(n_clay,i=1,nSubClays),n_silt1&
          ,n_silt2,n_silt3/)

     trrdry(n+1 : n+nraero_dust) = (/(dryEffRadMinerals(1 : nSubClays),i&
          =1,ntm_clay),                      &
          (dryEffRadMinerals(5),i=1,ntm_sil1)&
          ,                                  &
          (dryEffRadMinerals(6),i=1,ntm_sil2)&
          ,                                  &
          (dryEffRadMinerals(7),i=1,ntm_sil3)&
          /)
     !
#endif

     IF ( imDust>=4 ) CALL CALCSUBCLAYWEIGHTS

     wttr(n+1 : n+nraero_dust) = (/((SUBCLAYWEIGHTS(i,j),i=1,nSubClays&
          ),j=1,ntm_clay),                     &
          (1.D0,i=1,ntm_sil1+ntm_sil2+ntm_sil3+&
          ntm_sil4+ntm_sil5)/)

     ! Particle density of dust
#if defined( TRACERS_DUST_Silt4 ) &  defined( TRACERS_DUST_Silt5 )

     traden(n+1 : n+nraero_dust) = (/(TRPDENS(n_clay),i=1,nSubClays), &
          TRPDENS(n_silt1),TRPDENS(n_silt2), &
          TRPDENS(n_silt3),TRPDENS(n_silt4), &
          TRPDENS(n_silt5)/)*1D-3
     ! Convert from kg/m^3 to g/cm^3

#elif defined( TRACERS_DUST_Silt5 )

     traden(n+1 : n+nraero_dust) = (/(TRPDENS(n_clay),i=1,nSubClays), &
          TRPDENS(n_silt1),TRPDENS(n_silt2), &
          TRPDENS(n_silt3),TRPDENS(n_silt5)/)&
          *1D-3   ! Convert from kg/m^3 to g/cm^3
     !
#elif defined( TRACERS_DUST_Silt4 )

     traden(n+1 : n+nraero_dust) = (/(TRPDENS(n_clay),i=1,nSubClays), &
          TRPDENS(n_silt1),TRPDENS(n_silt2), &
          TRPDENS(n_silt3),TRPDENS(n_silt4)/)&
          *1D-3   ! Convert from kg/m^3 to g/cm^3
     !
#else

     traden(n+1 : n+nraero_dust) = (/(TRPDENS(n_clay),i=1,nSubClays), &
          TRPDENS(n_silt1),TRPDENS(n_silt2), &
          TRPDENS(n_silt3)/)*1D-3
     ! Convert from kg/m^3 to g/cm^3

#endif


#endif  /* TRACERS_MINERALS */

     itr(n+1 : n+nraero_dust) = 7
     ! all dust cases, outside ifdefs
     krhtra(n+1 : n+nraero_dust) = 0
     ! no deliq for dust or minerals
     fttasc(n+1 : n+nraero_dust) = 1.3D0
     ! increase dust AOD by 1.3 in LW
  ENDIF
  n = n + nraero_dust
#endif  /* (defined TRACERS_DUST) || (defined TRACERS_MINERALS) */
  !-----------------------------------------------------------------------
  !define ntrix_rf, based on the OMA tracers above
  IF ( n>0 ) THEN
     IF ( diag_fc==2 ) THEN
        ntrix_rf(1 : nraero_OMA) = ntrix_aod(1 : nraero_OMA)
     ELSEIF ( diag_fc==1 ) THEN
        ntrix_rf(1) = ntrix_aod(1)
     ENDIF
  ENDIF
  !-----------------------------------------------------------------------
#if (defined TRACERS_AMP) || (defined TRACERS_AMP_M1)
  IF ( nraero_AMP>0 ) THEN
     IF ( rad_interact_aer>0 ) THEN
        FS8OPX(1 : 7) = 0.D0
        FT8OPX(1 : 7) = 0.D0
     ENDIF
     ntrix_aod(n+1 : n+nraero_AMP) = (/n_N_AKK_1,n_N_ACC_1,n_N_DD1_1, &
          n_N_DS1_1,n_N_DD2_1,n_N_DS2_1,   &
          n_N_SSA_1,n_N_SSC_1,n_N_OCC_1,   &
          n_N_BC1_1,n_N_BC2_1,n_N_BC3_1,   &
          n_N_DBC_1,n_N_BOC_1,n_N_BCS_1,   &
          n_N_MXX_1/)
     IF ( diag_fc==2 ) THEN
        ntrix_rf(n+1 : n+nraero_AMP) = ntrix_aod(n+1 : n+nraero_AMP)
     ELSEIF ( diag_fc==1 ) THEN
        ntrix_rf(n+1) = ntrix_aod(n+1)
     ENDIF
  ENDIF
  n = n + nraero_AMP
#endif  /* (defined TRACERS_AMP) || (defined TRACERS_AMP_M1) */
  !-----------------------------------------------------------------------
#ifdef TRACERS_TOMAS
  IF ( nraero_TOMAS>0 ) THEN
     IF ( rad_interact_aer>0 ) THEN
        FS8OPX(1 : 2) = 0.D0
        FS8OPX(4 : 7) = 0.D0
        FT8OPX(1 : 2) = 0.D0
        FT8OPX(4 : 7) = 0.D0
#ifdef TRACERS_NITRATE
        FS8OPX(3) = 0.D0
        FT8OPX(3) = 0.D0
#endif  /* TRACERS_NITRATE */
     ENDIF

     ntrix_aod(n+1 : n+nraero_TOMAS)                                  &
          = (/N_ASO4(1),N_ANACL(1),N_AECOB(1),N_AECIL(1),N_AOCOB(1),  &
          N_AOCIL(1),N_ADUST(1)/)
     itr(n+1 : n+nraero_TOMAS) = (/1,2,6,5,4,4,7/)
     krhtra(n+1 : n+nraero_TOMAS) = 0
     ! ANUM(1) for internal-mixing case. Others(ncomp-1) for external-mixing case.
     IF ( diag_fc==2 ) THEN
        ntrix_rf(n+1 : n+nraero_TOMAS) = ntrix_aod(n+1 : n+nraero_TOMAS)
     ELSEIF ( diag_fc==1 ) THEN
        ntrix_rf(n+1) = ntrix_aod(n+1)
     ENDIF
  ENDIF
  n = n + nraero_TOMAS
#endif
  !=======================================================================
  !=======================================================================
#endif  /* TRACERS_ON */

  ! set default FSTOPX and FTTOPX values
  IF ( rad_interact_aer>0 ) THEN
     FSTOPX( : ) = 1.D0
     FTTOPX( : ) = 1.D0
  ELSE
     FSTOPX( : ) = 0.D0
     FTTOPX( : ) = 0.D0
  ENDIF
  skip_AOD_in_rad = rad_interact_aer>0

  IF ( ktrend/=0 ) THEN
     !****   Read in time history of well-mixed greenhouse gases
     CALL OPENUNIT('GHG',iu,.FALSE.,.TRUE.)
     CALL GHGHST(iu)
     CALL CLOSEUNIT(iu)
     IF ( FILE_EXISTS('dH2O') .AND. H2ObyCH4/=0. .AND. Kradia<=0 )  &
          THEN
        !****     Read in dH2O :  H2O prod.rate in kg/m^2 per day and ppm_CH4
        CALL OPENUNIT('dH2O',iu,.FALSE.,.TRUE.)
#if defined(CUBED_SPHERE)
        CALL READ_QMA(iu,plbx)
#else
        CALL GETQMA(iu,lat_dg,plbx,dh2o,lm,jm)
#endif
        CALL CLOSEUNIT(iu)
     ELSE
        H2ObyCH4 = 0.
     ENDIF
  ENDIF
#ifdef OLD_BCdalbsn
  IF ( dalbsnX/=0. ) THEN
     CALL UPDBCD(1990)
     depoBC_1990 = depoBC
  ENDIF
#endif
  !**** set up unit numbers for 14 more radiation input files
  donotread = -9999
  nrfun( : ) = 0 ! green light
  nrfun(12 : 13) = donotread ! not used in GCM
  nrfun(10 : 11) = donotread ! obsolete O3 data
  nrfun(6) = donotread     ! dust read externally now
  IF ( .NOT.transmission_corrections ) nrfun(4) = donotread
  IF ( madvol==0 ) nrfun(7) = donotread
  IF ( madeps==0 ) nrfun(8) = donotread
  !      if(ksolar < 0)  nrfun(9) = donotread
  nrfun(9) = donotread     ! open/read RADN9 inside RCOMP1
  DO IU = 1, 14
     IF ( nrfun(iu)==donotread ) CYCLE
     CALL OPENUNIT(RUNSTR(IU),NRFUN(IU),QBIN(IU),.TRUE.)
  ENDDO

  LS1_loc = 1
  ! default
  !***********************************************************************
  !     Main Radiative Initializations
  !     ------------------------------------------------------------------
  CALL RCOMP1(NRFUN)
  IF ( AM_I_ROOT() ) CALL WRITER(6,0)
  ! print rad. control parameters
  !***********************************************************************
  DO IU = 1, 14
     IF ( nrfun(iu)==donotread ) CYCLE
     CALL CLOSEUNIT(NRFUN(IU))
  ENDDO
  !**** Save initial (currently permanent and global) Q in rad.layers
  DO LR = 1, LM_REQ
     SHL0(LR) = SHL(LM+LR)
  ENDDO
  WRITE (out_line,*) 'spec.hum in rad.equ.layers : ', SHL0
  CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)

#ifdef ALTER_RADF_BY_LAT
  !**** Save initial rad forcing alterations : 
  FS8OPX_orig( : ) = FS8OPX( : )
  FT8OPX_orig( : ) = FT8OPX( : )                         ! aerosols

  !**** Read in the factors used for alterations : 
  CALL OPENUNIT('ALT_GHG_LAT',iu2,.FALSE.,.TRUE.)
  READ (iu2,*)
  ! skip first line
  DO n = 1, 46
     READ (iu2,'(a6,13D8.3)') skip, (FULGAS_lat(nn,n),nn=1,13)
  ENDDO
  CALL CLOSEUNIT(iu2)
  CALL OPENUNIT('ALT_AER_LAT',iu2,.FALSE.,.TRUE.)
  READ (iu2,*)
  ! skip first line
  DO n = 1, 46
     READ (iu2,'(a6,8D8.3)') skip, (FS8OPX_lat(nn,n),nn=1,8)
  ENDDO
  READ (iu2,*)
  ! skip first line
  DO n = 1, 46
     READ (iu2,'(a6,8D8.3)') skip, (FT8OPX_lat(nn,n),nn=1,8)
  ENDDO
  CALL CLOSEUNIT(iu2)
#endif

  ! transplanted from main().  needs reviving
  !      USE RAD_COM, only  :  dimrad_sv
  !      CHARACTER aDATE*14
  !      if (Kradia.ne.0 .and. Kradia<10) then
  !        write(aDATE(1 : 7),'(a3,I4.4)') aMON(1 : 3),Jyear
  !        if (Kradia.gt.0) aDATE(4 : 7)='    '
  !        call openunit(trim('RAD'//aDATE(1 : 7)),iu_RAD,.true.,.false.)
  !        if (Kradia.lt.0) call io_POS(iu_RAD,Itime-1,2*dimrad_sv,Nrad)
  !      end if

  IF ( rad_scm ) THEN
     IF ( FILE_EXISTS('GASES') ) THEN
        !     GAS NUMBER    1         2    3      4    5         6           7
        !                 H2O       CO2   O3     O2  NO2       N2O         CH4
        !     GAS NUMBER    8         9   10        11          12          13
        !              CCL3F1    CCL2F2   N2     CFC-Y       CFC-Z         SO2

        gasnames = (/ 'h2o   ','co2   ','o3    ','o2    ','no2   ', &
             'n2o   ','ch4   ','cfc11 ','cfc12 ','n2    ', &
             'cfc-y ','cfc-z ','so2   ' /)
        set_gases_internally = .FALSE.
        u0gas = 0.
        fid = PAR_OPEN(grid,'GASES','read')
        DO igas = 1, SIZE(gasnames)
           CALL READ_DATA(grid,fid,TRIM(gasnames(igas)),            &
                u0gas(1:LM, igas))
           u0gas(lm+1 : ,igas) = u0gas(lm,igas)
           ! fill lm+1:LM+lm_req
           IF ( TRIM(gasnames(igas))=='h2o' ) q(1,1, : )              &
                = u0gas(1:LM, igas)*(mwat/mair)    ! vol. ratio -> sp. hum.
           u0gas( : ,igas) = AML00*u0gas( : ,igas)                      &
                *((1D5/mair)*(gasc*tf/101325D0))
           ! vol. ratio -> cm-atm
        ENDDO
        CALL PAR_CLOSE(grid,fid)

        !fulgas = 1. ! needed?

        ulgas = u0gas

        ! Multiply gas amounts by rundeck scaling factors.
        ! Looping not an option since fulgas array does not yet
        ! contain the factors.

        !ulgas( : , 1) = ulgas( : , 1)*H2OstratX
        ulgas( : ,2) = ulgas( : ,2)*CO2X
        ulgas( : ,3) = ulgas( : ,3)*O3X
        ulgas( : ,4) = ulgas( : ,4)*O2X
        ulgas( : ,5) = ulgas( : ,5)*NO2X
        ulgas( : ,6) = ulgas( : ,6)*N2OX
        ulgas( : ,7) = ulgas( : ,7)*CH4X
        ulgas( : ,8) = ulgas( : ,8)*CFC11X
        ulgas( : ,9) = ulgas( : ,9)*CFC12X
        ulgas( : ,10) = ulgas( : ,10)*N2CX
        ulgas( : ,11) = ulgas( : ,11)*XGHGX
        ulgas( : ,12) = ulgas( : ,12)*YGHGX
        ulgas( : ,13) = ulgas( : ,13)*SO2X

     ENDIF
     IF ( FILE_EXISTS('VISAODangstr') ) THEN
        set_aerosols_internally = .FALSE.
        !          fid = par_open(grid,'VISAODangstr','read')
        !          not needed for initial CIRC cases which have zero aerosol
        !          todo :   read optical depths and scale with Angstrom exponent
        !          weighted by solar flux
        !          ....
        !          call par_close(grid,fid)
        sraext = 0.
        srasct = 0.
        sragcb = 0.
        srdext = 0.
        srdsct = 0.
        srdgcb = 0.
        srvext = 0.
        srvsct = 0.
        srvgcb = 0.
        srbext = 0.
        srbsct = 0.
        srbgcb = 0.
        traalk = 0.
        trdalk = 0.
        trvalk = 0.
        trbalk = 0.
     ENDIF
     i = 1
     j = 1
     DO l = 1, lm
        tloc = t(i,j,l)*pk(l,i,j)
        IF ( tloc>=tf ) THEN
           SVLHX(l,i,j) = lhe
        ELSE
           SVLHX(l,i,j) = lhs
        ENDIF
        SVLAT(l,i,j) = SVLHX(l,i,j)
        RHSAV(l,i,j) = q(i,j,l)/QSAT(tloc,SVLHX(l,i,j),PMID(l,i,j))
     ENDDO
     !llow=1; lmid=2; lhi=3
  ENDIF

END SUBROUTINE INIT_RAD

SUBROUTINE SETATM             ! dummy routine in gcm
END SUBROUTINE SETATM

SUBROUTINE GETVEG(LONR,LATR)  ! dummy routine in gcm
  INTEGER LONR, LATR
END SUBROUTINE GETVEG

SUBROUTINE DAILY_RAD(end_of_day)
  !@sum  daily_RAD sets radiation parameters that change every day
  !@auth G. Schmidt
  !@calls RADPAR : RCOMPT
  USE DOMAIN_DECOMP_ATM, ONLY : AM_I_ROOT
  USE DOMAIN_DECOMP_ATM, ONLY : GRID, GETDOMAINBOUNDS
  USE MODEL_COM, ONLY : MODELECLOCK
  USE RADPAR, ONLY : FULGAS, JYEARR => JYEAR, JDAYR => JDAY, XREF,    &
       KYEARV
#ifdef ALTER_RADF_BY_LAT
  USE RADPAR, ONLY : FULGAS_orig
#endif
  USE RADPAR, ONLY : RCOMPT, WRITET
  USE RAD_COM, ONLY : co2x, n2ox, ch4x, cfc11x, cfc12x, xGHGx,        &
       h2ostratx, o2x, no2x, n2cx, yghgx, so2x, o3x, o3_yr, ghg_yr,  &
       co2ppm, Volc_yr, albsn_yr, dalbsnX, SNOAGE, snoage_def,       &
       chl_from_seawifs
  USE DIAG_COM, ONLY : iwrite, jwrite, itwrite, TDIURN
  USE GEOM, ONLY : IMAXJ
  IMPLICIT NONE
  LOGICAL, INTENT(IN)  ::  end_of_day
  INTEGER  ::  year, dayOfYear
  INTEGER  ::  i, j, i_0, i_1, j_0, j_1, itype

  CALL MODELECLOCK%GET(year=year,dayOfYear=dayOfYear)

  !**** Update time dependent radiative parameters each day
  !     Get black carbon deposition data for the appropriate year
  !     (does nothing except at a restart or the beginning of a new year)
  IF ( dalbsnX/=0. ) THEN
     IF ( albsn_yr==0 ) THEN
#ifdef OLD_BCdalbsn
        CALL UPDBCD(year)
#else
        CALL UPDBCDALBSN(year,dayofyear)
#endif
     ELSE
#ifdef OLD_BCdalbsn
        CALL UPDBCD(albsn_yr)
#else
        ! as per radiation-code convention, pass -albsn_yr to indicate
        ! perpetual-year mode
        CALL UPDBCDALBSN(-albsn_yr,dayofyear)
#endif
     ENDIF
  ENDIF
  !     Hack :  2 specific volc. eruption scenarios for 2000-2100 period
  IF ( volc_yr==-2010 ) THEN             ! repeat some old volcanos
     KYEARV = YEAR
     IF ( YEAR>2010 ) KYEARV = YEAR - 100
     ! go back 100 years
  ENDIF
  IF ( volc_yr==-2000 ) THEN
     KYEARV = YEAR
     IF ( YEAR>2000 ) KYEARV = YEAR - 50
     ! go back 50 years til 2050
     IF ( YEAR>2050 ) KYEARV = YEAR - 150
     ! then go back 150 years
  ENDIF

  JDAYR = dayOfYear
  JYEARR = YEAR
  CALL RCOMPT
  !     FULGAS(2 : ) is set only in the first call to RCOMPT unless ghg_yr=0
  !     Optional scaling of the observed value only in case it was (re)set
  IF ( .NOT.end_of_day .AND. H2OstratX>=0. ) FULGAS(1) = FULGAS(1)  &
       *H2OstratX
  IF ( .NOT.end_of_day .OR. O3_yr==0. ) FULGAS(3) = FULGAS(3)*O3X
  IF ( ghg_yr==0 .OR. .NOT.end_of_day ) THEN
     FULGAS(2) = FULGAS(2)*CO2X
     FULGAS(6) = FULGAS(6)*N2OX
     FULGAS(7) = FULGAS(7)*CH4X
     FULGAS(8) = FULGAS(8)*CFC11X
     FULGAS(9) = FULGAS(9)*CFC12X
     FULGAS(11) = FULGAS(11)*XGHGX
     FULGAS(12) = FULGAS(12)*YGHGX
  ENDIF
  IF ( .NOT.end_of_day ) THEN
     FULGAS(4) = FULGAS(4)*O2X
     FULGAS(5) = FULGAS(5)*NO2X
     FULGAS(10) = FULGAS(10)*N2CX
     FULGAS(13) = FULGAS(13)*SO2X
     ! no effect since FULGAS(13)=0.
  ENDIF

  !**** write trend table for forcing 'itwrite' for years iwrite->jwrite
  !**** itwrite :  1-2=GHG 3=So 4-5=O3 6-9=aerosols :  Trop,DesDust,Volc,Total
  IF ( AM_I_ROOT() .AND. jwrite>1500 )                              &
       CALL WRITET(6,itwrite,iwrite,jwrite,1,0)

#ifdef ALTER_RADF_BY_LAT
  !**** Save initial rad forcing alterations : 
  FULGAS_orig( : ) = FULGAS( : )
  ! GHGs
#endif

  !**** Define CO2 (ppm) for rest of model
  co2ppm = FULGAS(2)*XREF(1)

  IF ( chl_from_seawifs>0 ) CALL GET_CHL_FROM_SEAWIFS

  IF ( end_of_day ) THEN

     CALL GETDOMAINBOUNDS(grid,J_STRT=J_0,J_STOP=J_1)
     CALL GETDOMAINBOUNDS(grid,I_STRT=I_0,I_STOP=I_1)

     DO j = J_0, J_1
        DO i = I_0, IMAXJ(j)
           !****
           !**** increase snow age depending on snoage_def
           !****
           IF ( snoage_def==0 ) THEN
              ! update indep. of ts
              DO itype = 1, 3
                 SNOAGE(itype,i,j) = 1. + .98D0*SNOAGE(itype,i,j)
              ENDDO
           ELSEIF ( snoage_def==1 ) THEN
              ! update if max T>0
              IF ( TDIURN(i,j,7)>0 ) SNOAGE(1,i,j)                  &
                   = 1. + .98D0*SNOAGE(1,i,j)
              ! ocean ice (not currently used)
              IF ( TDIURN(i,j,8)>0 ) SNOAGE(2,i,j)                  &
                   = 1. + .98D0*SNOAGE(2,i,j)
              ! land ice
              IF ( TDIURN(i,j,2)>0 ) SNOAGE(3,i,j)                  &
                   = 1. + .98D0*SNOAGE(3,i,j)
              ! land
           ELSE
              WRITE (6,*) "This snoage_def is not defined :  ",       &
                   snoage_def
              WRITE (6,*) "Please use :  0 (update indep of T)"
              WRITE (6,*) "            1 (update if T>0)"
              CALL STOP_MODEL('stopped in RAD_DRV.f',255)
           ENDIF
        ENDDO
     ENDDO
  ENDIF

END SUBROUTINE DAILY_RAD

SUBROUTINE GET_CHL_FROM_SEAWIFS

  USE DOMAIN_DECOMP_ATM, ONLY : GRID, REWIND_PARALLEL, READT_PARALLEL,&
       GETDOMAINBOUNDS
  USE FLUXES, ONLY : FOCEAN, atmocn
  USE CONSTANT, ONLY : by12
  USE MODEL_COM, ONLY : MODELECLOCK, calendar
  USE RESOLUTION, ONLY : im, jm
  USE FILEMANAGER, ONLY : NAMEUNIT
  USE GEOM, ONLY : IMAXJ
  USE FILEMANAGER, ONLY : OPENUNIT
  USE CALENDARMONTH_MOD
  IMPLICIT NONE

  REAL*8  ::  TEMP_LOCAL(grid%I_STRT_HALO:grid%I_STOP_HALO,         &
       grid%J_STRT_HALO:grid%J_STOP_HALO,2)
  INTEGER  ::  month, date, year
  LOGICAL  ::  HAVE_NORTH_POLE, HAVE_SOUTH_POLE
  INTEGER  ::  LSTMON, I, J, J_0, J_1, I_0, I_1
  INTEGER, SAVE  ::  IMON0 = 0
  INTEGER, SAVE  ::  iu_chl = -1
  !@var ACHL,ECHL1,ECHL0,BCHL,CCHL arrays for the reading in chlorophyll
  REAL*8, ALLOCATABLE, DIMENSION( : , : ), SAVE  ::  ACHL, ECHL1, ECHL0,&
       BCHL, CCHL
  REAL*8  ::  TIME
  INTEGER  ::  I_0H, I_1H, J_0H, J_1H
  TYPE (CALENDARMONTH) :: cMonth

  I_0H = grid%I_STRT_HALO
  I_1H = grid%I_STOP_HALO
  J_0H = grid%J_STRT_HALO
  J_1H = grid%J_STOP_HALO
  IF ( iu_chl<0 ) THEN
     CALL OPENUNIT("CHL_DATA",iu_CHL,.TRUE.,.TRUE.)
     ALLOCATE (ACHL(I_0H : I_1H,J_0H : J_1H),ECHL1(I_0H : I_1H,J_0H : J_1H),&
          ECHL0(I_0H : I_1H,J_0H : J_1H),BCHL(I_0H : I_1H,J_0H : J_1H),&
          CCHL(I_0H : I_1H,J_0H : J_1H))
  ENDIF
  CALL MODELECLOCK%GET(month=month,date=date)
  CALL GETDOMAINBOUNDS(GRID,J_STRT=J_0,J_STOP=J_1,                  &
       HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,             &
       HAVE_NORTH_POLE=HAVE_NORTH_POLE)
  I_0 = grid%I_STRT
  I_1 = grid%I_STOP

  !**** Read in Seawifs files here
  IF ( month/=IMON0 ) THEN
     IF ( IMON0==0 ) THEN
        !**** READ IN LAST MONTH''S END-OF-MONTH DATA
        LSTMON = month - 1
        IF ( lstmon==0 ) lstmon = 12
        CALL READT_PARALLEL(grid,iu_CHL,NAMEUNIT(iu_CHL),TEMP_LOCAL,&
             LSTMON)
        ECHL0 = TEMP_LOCAL( : , : ,2)
     ELSE
        !**** COPY END-OF-OLD-MONTH DATA TO START-OF-NEW-MONTH DATA
        ECHL0 = ECHL1
     ENDIF
     !**** READ IN CURRENT MONTHS DATA :  MEAN AND END-OF-MONTH
     IMON0 = month
     IF ( month==1 ) CALL REWIND_PARALLEL(iu_CHL)
     CALL READT_PARALLEL(grid,iu_CHL,NAMEUNIT(iu_CHL),TEMP_LOCAL,1)
     ACHL = TEMP_LOCAL( : , : ,1)
     ECHL1 = TEMP_LOCAL( : , : ,2)

     !**** FIND INTERPOLATION COEFFICIENTS (LINEAR/QUADRATIC FIT)
     DO J = J_0, J_1
        DO I = I_0, IMAXJ(J)
           BCHL(I,J) = ECHL1(I,J) - ECHL0(I,J)
           CCHL(I,J) = 3.*(ECHL1(I,J)+ECHL0(I,J)) - 6.*ACHL(I,J)
        ENDDO
     ENDDO
  ENDIF
  !**** Calculate CHL for current day
  cMonth = calendar%GETCALENDARMONTH(month,year)
  TIME = (DATE-.5)/cMonth%DAYSINMONTH - .5
  ! -.5<TIME<.5
  DO J = J_0, J_1
     DO I = I_0, IMAXJ(J)
        IF ( FOCEAN(I,J)>0 ) THEN
           !**** CHL always uses quadratic fit
           atmocn%CHL(I,J) = ACHL(I,J) + BCHL(I,J)*TIME + CCHL(I,J) &
                *(TIME**2-BY12)
           IF ( atmocn%CHL(I,J)<0 ) atmocn%CHL(I,J) = 0.
           ! just in case
        ENDIF
     ENDDO
  ENDDO
  !**** REPLICATE VALUES AT POLE
  IF ( HAVE_NORTH_POLE ) THEN
     IF ( FOCEAN(1,JM)>0 ) atmocn%CHL(2 : IM,JM) = atmocn%CHL(1,JM)
  ENDIF
  IF ( HAVE_SOUTH_POLE ) THEN
     IF ( FOCEAN(1,1)>0 ) atmocn%CHL(2 : IM,1) = atmocn%CHL(1,1)
  ENDIF
  atmocn%CHL_DEFINED = .TRUE.

END SUBROUTINE GET_CHL_FROM_SEAWIFS

SUBROUTINE DAILY_ORBIT(end_of_day)
  !@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
  !@auth Original Development Team
  !@calls constant : orbit
  USE MODEL_COM, ONLY : MODELECLOCK
  USE RAD_COM, ONLY : RSDIST, COSD, SIND, COSZ_day, SUNSET,           &
       VARIABLE_ORB_PAR, ORB_PAR_YEAR_BP, USEORBIT => ORBIT
  USE DOMAIN_DECOMP_ATM, ONLY : AM_I_ROOT
  USE RAD_COSZ0, ONLY : DAILY_COSZ
  USE BASETIME_MOD
  USE TIMEINTERVAL_MOD
  USE RATIONAL_MOD
  IMPLICIT NONE
  REAL*8  ::  SUNLON, SUNLAT, LAM, EDPY, VEDAY, PYEAR
  LOGICAL, INTENT(IN)  ::  end_of_day
  INTEGER  ::  year, dayOfYear
  TYPE (BASETIME) :: t
  REAL*8  ::  declinationAngle
  TYPE (TIMEINTERVAL) :: halfDay

  CALL MODELECLOCK%GET(year=year,dayOfYear=dayOfYear)

  !**** CALCULATE SOLAR ANGLES AND ORBIT POSITION
  !**** This is for noon (GMT) for new day.

  !**** The orbital calculation will need to vary depending on the kind
  !**** of calendar adopted (i.e. a generic 365 day year, or a transient
  !**** calendar including leap years etc.).  For transient calendars the
  !**** dayOfYear passed to orbit needs to be adjusted to represent the number
  !**** of days from Jan 1 2000AD.
  !      EDPY=365.2425d0, VEDAY=79.3125d0  ! YR 2000AD
  !      dayOfYear => dayOfYear + 365 * (YEAR-2000) + appropriate number of leaps
  !**** Default calculation (no leap, VE=Mar 21 hr 0)
  !      EDPY=365d0 ; VEDAY=79d0           ! Generic year
  !**** PMIP calculation (no leap, VE=Mar 21 hr 12)
  EDPY = 365D0
  VEDAY = 79.5D0                      ! Generic year
  !**** Update orbital parameters at start of year
  IF ( dayOfYear==1 ) CALL USEORBIT%SETYEAR(REAL(year,KIND=8))

  ! Use time for the _middle_ of the day to compute
  ! zenith angle : 

  halfDay = TIMEINTERVAL(USEORBIT%GETMEANDAY()/2)
  t = NEWBASETIME(MODELECLOCK%GETTIMEATBEGINNINGOFCURRENTDAY()      &
       +halfDay)

  sinD = USEORBIT%GETSINDECLINATIONANGLE(t)
  cosD = SQRT(1-sinD**2)
  rsdist = USEORBIT%GETDISTANCE(t)**2

  CALL DAILY_COSZ(sind,cosd,cosz_day,sunset)

END SUBROUTINE DAILY_ORBIT

SUBROUTINE DAILY_CH4OX(end_of_day)
  !@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
  !@vers 2013/03/27
  !@auth Original Development Team
  !@calls constant : orbit
  USE RESOLUTION, ONLY : im, jm, lm
  USE ATM_COM, ONLY : Q
  USE MODEL_COM, ONLY : MODELECLOCK
  USE MODEL_COM, ONLY : itime
  USE GEOM, ONLY : AXYP, IMAXJ, LAT2D
  USE ATM_COM, ONLY : BYMA
  USE RADPAR, ONLY : GHGAM, ghgyr2, ghgyr1
  USE RAD_COM, ONLY : DH2O, H2ObyCH4, ghg_yr
#ifdef TRACERS_WATER
  USE OLDTRACER_MOD, ONLY : TR_WD_TYPE, NWATER, TR_H2OBYCH4, ITIME_TR0
  USE TRACER_COM, ONLY : TRM, NTM
#endif
  USE DIAG_COM, ONLY : FTYPE, ntype, AIJ => AIJ_LOC
  USE DIAG_COM_RAD, ONLY : j_h2och4, ij_h2och4
  USE DOMAIN_DECOMP_ATM, ONLY : grid, GETDOMAINBOUNDS, AM_I_ROOT
  IMPLICIT NONE
  REAL*8  ::  xCH4, xdH2O
  INTEGER i, j, l, iy, it
  LOGICAL, INTENT(IN)  ::  end_of_day
#ifdef TRACERS_WATER
  INTEGER n
#endif
  !**** Extract domain decomposition info
  INTEGER  ::  J_0, J_1, I_0, I_1
  LOGICAL  ::  HAVE_SOUTH_POLE, HAVE_NORTH_POLE
  INTEGER  ::  year, month

  CALL MODELECLOCK%GET(year=year,month=month)

  CALL GETDOMAINBOUNDS(grid,J_STRT=J_0,J_STOP=J_1,                  &
       HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,             &
       HAVE_NORTH_POLE=HAVE_NORTH_POLE)
  I_0 = grid%I_STRT
  I_1 = grid%I_STOP

  IF ( .NOT.end_of_day ) RETURN

  !**** Tasks to be done at end of day only
  IF ( H2ObyCH4>0 ) THEN
     !****   Add obs. H2O generated by CH4(*H2ObyCH4) using a 2 year lag
     iy = year - 2 - ghgyr1 + 1
     IF ( ghg_yr>0 ) iy = ghg_yr - 2 - ghgyr1 + 1
     IF ( iy<1 ) iy = 1
     IF ( iy>ghgyr2-ghgyr1+1 ) iy = ghgyr2 - ghgyr1 + 1
     xCH4 = GHGAM(3,iy)*H2ObyCH4
     !        If (AM_I_ROOT())
     !          write(6,*) 'add in stratosphere :  H2O gen. by CH4(ppm)=',xCH4

     DO l = 1, lm
        DO j = J_0, J_1
           DO i = I_0, IMAXJ(j)
#ifdef CUBED_SPHERE
              CALL LAT_INTERP_QMA(LAT2D(i,j),l,month,xdH2O)
#else
              xdH2O = DH2O(j,l,month)
#endif
              Q(i,j,l) = Q(i,j,l) + xCH4*xdH2O*BYMA(l,i,j)
#ifdef TRACERS_WATER
              !**** Add water to relevant tracers as well
              DO n = 1, ntm
                 IF ( ITIME_TR0(n)<=itime ) THEN
                    SELECT CASE (TR_WD_TYPE(n))
                    CASE (NWATER)
                       ! water :  add CH4-sourced water to tracers
                       TRM(i,j,l,n) = TRM(i,j,l,n) + TR_H2OBYCH4(n) &
                            *xCH4*xdH2O*AXYP(i,j)
                    ENDSELECT
                 ENDIF
              ENDDO
#endif
              DO it = 1, ntype
                 CALL INC_AJ(i,j,it,j_h2och4,                       &
                      xCH4*xdH2O*FTYPE(it,i,j))
              ENDDO
              AIJ(i,j,ij_h2och4) = AIJ(i,j,ij_h2och4) + xCH4*xdH2O
           ENDDO
        ENDDO
        IF ( HAVE_NORTH_POLE ) Q(2 : im,jm,l) = Q(1,jm,l)
        IF ( HAVE_SOUTH_POLE ) Q(2 : im,1,l) = Q(1,1,l)
#ifdef TRACERS_WATER
        DO n = 1, ntm
           IF ( HAVE_SOUTH_POLE ) TRM(2 : im,1,l,n) = TRM(1,1,l,n)
           IF ( HAVE_NORTH_POLE ) TRM(2 : im,jm,l,n) = TRM(1,jm,l,n)
        ENDDO
#endif
     ENDDO
  ENDIF

END SUBROUTINE DAILY_CH4OX

SUBROUTINE RADIA
  !@sum  RADIA adds the radiation heating to the temperatures
  !@vers 2013/03/27
  !@auth Original Development Team
  !@calls tropwmo,coszs,coszt, RADPAR : rcompx ! writer,writet
  USE CONSTANT, ONLY : lhe, lhs, twopi, tf, stbo, rhow, mair, grav,   &
       bysha, pi, radian, areag
  USE RESOLUTION, ONLY : pmtop
  USE RESOLUTION, ONLY : im, jm, lm
#ifdef TRACERS_SPECIAL_Shindell
  USE RESOLUTION, ONLY : LS1 => LS1_NOMINAL
#endif
  USE ATM_COM, ONLY : kradia, lm_req, p, t, Q, iu_rad, req_fac_d
  USE MODEL_COM
  USE TIMECONSTANTS_MOD, ONLY : SECONDS_PER_DAY, INT_DAYS_PER_YEAR
  USE ATM_COM, ONLY : BYAML00
  USE GEOM, ONLY : IMAXJ, AXYP, BYAXYP, LAT2D, LON2D
  ! for threadprivate copyin common block
  !     INPUT DATA         ! not (i,j) dependent
  USE RADPAR, ONLY : LX, tauwc0, tauic0, WRITER, RCOMPX, UPDGHG,      &
       S00WM2, RATLS0, S0, JYEARR => JYEAR, JDAYR => JDAY, FULGAS,   &
       use_tracer_chem, FS8OPX, FT8OPX, use_o3_ref, KYEARG, KJDAYG,  &
       planck_tmin, planck_tmax
  ! set in radpar block data
#ifdef ALTER_RADF_BY_LAT
  USE RADPAR, ONLY : FS8OPX_orig, FT8OPX_orig, FULGAS_orig
#endif
  !     INPUT DATA  (i,j) dependent
  USE RADPAR, ONLY : JLAT46 => JLAT, ILON72 => ILON, JGCM, IGCM, L1,  &
       LMR => NL, PLB, TLB, TLM, SHL, RHL, ltopcl, TAUWC, TAUIC,     &
       SIZEWC, SIZEIC, kdeliq, POCEAN, PEARTH, POICE, PLICE, PLAKE,  &
       COSZ, PVT, TGO, TGE, TGOI, TGLI, TSL, WMAG, WEARTH, AGESN,    &
       SNOWD, SNOWOI, SNOWLI, dALBsn, ZSNWOI, ZOICE, zmp, fmp, flags,&
       LS1_loc, snow_frac, zlake, TRACER, FSTOPX, FTTOPX, chem_IN,   &
       nraero_aod => NTRACE, FTAUC, LOC_CHL, FSTASC, FTTASC
#ifdef HEALY_LM_DIAGS
  USE RADPAR, ONLY : VTAULAT
#endif
#ifdef GCC_COUPLE_RAD
  USE RADPAR, ONLY : GCCco2_IN, use_tracer_GCCco2, GCCCO2_OUT
#endif

  !     OUTPUT DATA
  USE RADPAR, ONLY : TRDFLB, TRNFLB, TRUFLB, TRFCRL, chem_out, SRDFLB,&
       SRNFLB, SRUFLB, SRFHRL, PLAVIS, PLANIR, ALBVIS, ALBNIR,       &
       FSRNFG, SRRVIS, SRRNIR, SRAVIS, SRANIR, SRXVIS, SRDVIS,       &
       BTEMPW, SRAEXT, SRASCT, SRAGCB, SRDEXT, SRDSCT, SRDGCB,       &
       SRVEXT, SRVSCT, SRVGCB, aesqex, aesqsc, aesqcb, CO2outCol,    &
       aesqex_dry, aesqsc_dry, aesqcb_dry, SRXNIR, SRDNIR
  USE RAD_COM, ONLY : modrd, nrad
  USE RAD_COM, ONLY : rqt, SRHR, TRHR, FSF, COSZ1, s0x, rsdist,       &
       nradfrc, CH4X_RADoverCHEM, snoage, PLB0, SHL0, TCHG, ALB,     &
       FSRDIR, SRVISSURF, SRDN, cfrac, RCLD, chem_tracer_save,       &
       rad_interact_aer, kliq, RHfix, CLDx, GHG_YR, CO2X, N2OX, CH4X,&
       CFC11X, CFC12X, XGHGX, rad_forc_lev, NTRIX_AOD, NTRIX_RF,     &
       WTTR, cloud_rad_forc, CC_cdncx, OD_cdncx, cdncl, dALBsnX,     &
       rad_to_chem, TRSURF, DIRVIS, FSRDIF, DIRNIR, DIFNIR,          &
       aer_rad_forc, clim_interact_chem, TAUSUMW, TAUSUMI,           &
       TAero_aod_diag, chl_from_obio, chl_from_seawifs
#ifdef GCC_COUPLE_RAD
  USE RAD_COM, ONLY : GCCCO2_TRACER_SAVE, GCCCO2RAD_TO_CHEM
#endif
#ifdef GCAP
  USE RAD_COM, ONLY : TAUW3D, TAUI3D
#endif
#ifdef mjo_subdd
  USE RAD_COM, ONLY : SWHR, LWHR, SWHR_cnt, LWHR_cnt, OLR_ACC,        &
       OLR_cnt, SWU_AVG, swu_cnt
#endif
#ifdef ALTER_RADF_BY_LAT
  USE RAD_COM, ONLY : FULGAS_lat, FS8OPX_lat, FT8OPX_lat
#endif
#ifdef TRACERS_DUST
  USE RAD_COM, ONLY : srnflb_save, trnflb_save
#endif
#if (defined SHINDELL_STRAT_EXTRA) &(defined ACCMIP_LIKE_DIAGS)
  USE RAD_COM, ONLY : STRATO3_TRACER_SAVE
#endif
#ifdef TRACERS_ON
  USE RAD_COM, ONLY : tau_as, tau_cs, tau_dry, nraero_rf
#ifdef CACHED_SUBDD
  USE CONSTANT, ONLY : grav, Rgas
  USE RAD_COM, ONLY : abstau_as, abstau_cs, abstau_dry, swfrc, lwfrc
  USE RUNTIMECONTROLS_MOD, ONLY : tracers_amp, tracers_tomas
#endif  /* CACHED_SUBDD */
#endif
  USE RANDOM
  USE CLOUDS_COM, ONLY : TAUSS, TAUMC, SVLHX, RHSAV, SVLAT, CLDSAV,   &
       CLDMC, CLDSS, CSIZMC, CSIZSS, llow, lmid, lhi, FSS, TAUSSIP,  &
       CSIZSSIP, QLSS, QISS, QLMC, QIMC, GET_CLD_OVERLAP
  !  subroutine
#ifdef GCAP
  USE CLOUDS_COM, ONLY : CLDSS3D
  USE CONSTANT, ONLY : teeny
#endif
  USE DIAG_COM, ONLY : ia_rad, JREG, AIJ => AIJ_LOC, AIJL => AIJL_LOC,&
       ntype, FTYPE, itocean, itlake, itearth, itlandi, itoice,      &
       itlkice, ADIURN => ADIURN_LOC, ndiuvar, ia_rad_frc
#ifdef USE_HDIURN
  USE DIAG_COM, ONLY : HDIURN => HDIURN_LOC
#endif
  USE DIAG_COM, ONLY : iwrite, jwrite, itwrite, ndiupt, IJDD, AFLX_ST,&
       hr_in_day, hr_in_month
  USE DIAG_COM_RAD
#ifdef TRACERS_ON
  USE DIAG_COM, ONLY : adiurn_dust, SAVE3DAOD
  USE RAD_COM, ONLY : diag_fc
#endif
  USE ATM_COM, ONLY : PK, PEDN, PMID, PDSIG, ltropo, MA, BYMA
  USE SEAICE_COM, ONLY : si_atm
  USE GHY_COM, ONLY : FEARTH, snowd_ij => snowd
  USE ENT_COM, ONLY : ENTCELLS
  USE ENT_MOD, ONLY : ENT_GET_EXPORTS, N_COVERTYPES
  !YKIM-temp hack
  USE ENT_DRV, ONLY : MAP_ENT2GISS    !YKIM-temp hack
  USE LAKES_COM, ONLY : flake, dlake !,mwl
  USE FLUXES, ONLY : ASFLX4, atmocn, atmice, atmgla, atmlnd, atmsrf,  &
       FLICE, FLAND, FOCEAN
  USE DOMAIN_DECOMP_ATM, ONLY : grid, WRITE_PARALLEL
  USE DOMAIN_DECOMP_ATM, ONLY : GLOBALSUM, GETDOMAINBOUNDS
  USE RAD_COSZ0, ONLY : COSZT, COSZS

#ifdef TRACERS_ON
  USE OLDTRACER_MOD, ONLY : TRNAME, TRPDENS
  USE TRACER_COM, ONLY : NTM, n_Ox, TRM, n_OCB, n_BCII, n_BCIA,       &
       n_OCIA, N_OCII, N_SO4_D2, N_SO4_D3, N_SO4, n_stratOx,         &
       N_N_AKK_1
#ifdef TRACERS_NITRATE
  USE OLDTRACER_MOD, ONLY : TR_MM
  USE TRACER_COM, ONLY : n_NH4, n_NO3p
#endif
#ifdef TRACERS_AEROSOLS_SOA
  USE TRACER_COM, ONLY : n_isopp1a, n_isopp2a
#ifdef TRACERS_TERP
  USE TRACER_COM, ONLY : n_apinp1a, n_apinp2a
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_OCEAN
  USE TRACER_COM, ONLY : n_ococean
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifdef GCC_COUPLE_RAD
  USE TRACER_COM, ONLY : n_CO2n
  USE CONSTANT, ONLY : avog
  USE OLDTRACER_MOD, ONLY : TR_MM
#endif
#ifdef TRACERS_AEROSOLS_VBS
  USE TRACERS_VBS, ONLY : vbs_tr
#endif
  USE TRDIAG_COM, ONLY : TAIJS => TAIJS_LOC, taijls => TAIJLS_LOC,    &
       IJTS_FC, IJTS_TAU, IJTS_TAUSUB, IJTS_FCSUB, IJLT_3DTAU,       &
       IJLT_3DAAOD, IJLT_3DTAUCS, IJLT_3DAAODCS, IJLT_3DTAUDRY,      &
       IJLT_3DAAODDRY, IJTS_SQEX, IJTS_SQEXSUB, IJTS_SQSC,           &
       IJTS_SQSCSUB, IJTS_SQCB, IJTS_SQCBSUB, diag_rad, diag_aod_3d, &
       save_dry_aod
#ifdef AUXILIARY_OX_RADF
  USE TRDIAG_COM, ONLY : IJTS_AUXFC
#endif /* AUXILIARY_OX_RADF */
#ifdef BC_ALB
  USE TRDIAG_COM, ONLY : IJTS_ALB, ijts_sunlit_snow
#endif  /* BC_ALB */
#ifdef TRACERS_SPECIAL_Shindell
  USE TRCHEM_SHINDELL_COM, ONLY : Lmax_rad_O3, Lmax_rad_CH4
#endif /* TRACERS_SPECIAL_Shindell */
#ifdef TRACERS_TOMAS
  USE TOMAS_AEROSOL, ONLY : icomp
#endif
#endif /* TRACERS_ON */
  USE AERPARAM_MOD, ONLY : DCDNC_EST
#ifdef OLD_BCdalbsn
  USE AERPARAM_MOD, ONLY : DEPOBC, DEPOBC_1990
#else
  USE AERPARAM_MOD, ONLY : BCDALBSN
#endif
  USE TIMERPACKAGE_MOD, ONLY : STARTTIMER => START, STOPTIMER => STOP
  USE DICTIONARY_MOD, ONLY : GET_PARAM, IS_SET_PARAM
#ifdef CACHED_SUBDD
  USE SUBDD_MOD, ONLY : sched_rad, SUBDD_GROUPS, SUBDD_TYPE,          &
       subdd_ngroups, INC_SUBDD, FIND_GROUPS, lmaxsubdd
#endif
#ifdef SCM
  USE SCM_COM, ONLY : SCMopt, SCMin
  USE CONSTANT, ONLY : SHA
  USE ATM_COM, ONLY : QCL
#endif
  USE DIAG_COM, ONLY : IJ_NINTAEREXT, IJ_NINTAERSCA, IJ_NINTAERASY
  USE RADPAR, ONLY : nintaerext, nintaersca, nintaerasy

#ifdef GCAP
  USE RAD_COM,  ONLY : SAVE_ALB, save_cosz2
  USE O3MOD,    ONLY : SAVE_TO3
#endif
#ifdef TRACERS_GC
  USE CHEM_COM,          ONLY : TrM, i_O3, i_CH4
  USE DICTIONARY_MOD,    ONLY : SYNC_PARAM
  USE DOMAIN_DECOMP_ATM, ONLY : AM_I_ROOT
  USE RAD_COM,           ONLY : SAVE_RF, SAVE_RF_TP, SAVE_RF_3D
#endif

  IMPLICIT NONE
  REAL*8 dz, rho
  !
  !@var wtrtau,icetau per-layer opacity for cloud water,ice
  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
       grid%J_STRT_HALO:grid%J_STOP_HALO,lm)           &
       ::  wtrtau, icetau
#ifdef SCM
  REAL*8 q_above(LM+1), q_below(LM+1), Frad(LM+1)
#endif
  !     INPUT DATA   partly (i,j) dependent, partly global
  REAL*8 U0GAS, taulim
#ifdef OLD_BCdalbsn
  REAL*8 xdalbs, sumda, tauda, fsnow
  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
       grid%J_STRT_HALO:grid%J_STOP_HALO)              &
       ::  sumda_psum, tauda_psum
#endif
  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
       grid%J_STRT_HALO:grid%J_STOP_HALO)  ::  COSZ2,  &
       COSZA, TRINCG, BTMPW, WSOIL, fmp_com
  REAL*8, DIMENSION(4,grid%I_STRT_HALO:grid%I_STOP_HALO,            &
                      grid%J_STRT_HALO:grid%J_STOP_HALO)  ::  SNFS, TNFS
  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
       grid%J_STRT_HALO:grid%J_STOP_HALO)  ::  SNFSCRF,&
       TNFSCRF, SNFSCRF2, TNFSCRF2, LWDNCS,            &
       SNFS_AS_noA, TNFS_AS_noA, SNFS_CS_noA,          &
       TNFS_CS_noA, SWUS, CTT, CTP, WTRCLD, ICECLD
  REAL*8, DIMENSION(18,grid%I_STRT_HALO:grid%I_STOP_HALO,           &
       grid%J_STRT_HALO:grid%J_STOP_HALO)              &
       ::  SNFSAERRF, TNFSAERRF
#ifdef CFMIP3_SUBDD
  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
       grid%J_STRT_HALO:grid%J_STOP_HALO)  ::  swut,   &
       swutcs, cfmip_twp, swdcls, swucls, swdt
  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
       grid%J_STRT_HALO:grid%J_STOP_HALO,lm)           &
       ::  cfmip_cf, cfmip_qci, cfmip_qcl
#endif
#ifdef CACHED_SUBDD
  INTEGER  ::  igrp, ngroups, grpids(subdd_ngroups)
  TYPE (SUBDD_TYPE), POINTER :: subdd
  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
       grid%J_STRT_HALO:grid%J_STOP_HALO)  ::  SDDARR
  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
       grid%J_STRT_HALO:grid%J_STOP_HALO,lm)           &
       ::  SDDARR3D
#ifdef TRACERS_ON
  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
       grid%J_STRT_HALO:grid%J_STOP_HALO,nraero_rf)    &
       ::  sddarr3drf
  INTEGER  ::  f
#endif  /* TRACERS_ON */
#ifdef SCM
  !     radiative flux profiles for sub-daily output, generalized
  !     for GCM grid but currently limited to SCM use
  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
       grid%J_STRT_HALO:grid%J_STOP_HALO,lm)           &
       ::  TRDFLB_prof, TRUFLB_prof, SRDFLB_prof,     &
       SRUFLB_prof
#endif
#ifdef TRACERS_ON
  ! types of aods to be saved
  ! The name will be any combination of {,TRNAME}{as,cs}{,a}aod
  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
       grid%J_STRT_HALO:grid%J_STOP_HALO,lm,nraero_aod)&
       ::  sddarr4d
  CHARACTER(LEN=10), DIMENSION(2)                                   &
       ::  sgroups = (/'taijh ','taijlh'/)
  CHARACTER(LEN=10), DIMENSION(3)  ::  ssky = (/'as ','cs ','dry'/)
  CHARACTER(LEN=10), DIMENSION(2)  ::  sabs = (/' ','a'/)
  CHARACTER(LEN=10), DIMENSION(2)  ::  sfrc = (/'swf','lwf'/)
  CHARACTER(LEN=10)  ::  spcname
  CHARACTER(LEN=50)  ::  sname
  INTEGER  ::  g, s, a
#endif  /* TRACERS_ON */
  !@var CO2out for holding 3D CO2 from rad code for SUBDD
  REAL*8, DIMENSION(LM,grid%I_STRT_HALO:grid%I_STOP_HALO,           &
       grid%J_STRT_HALO:grid%J_STOP_HALO)  ::  CO2out
#endif /* CACHED_SUBDD */
#if (defined ACCMIP_LIKE_DIAGS)
#ifndef SKIP_ACCMIP_GHG_RADF_DIAGS
  !@var snfs_ghg,tnfs_ghg like SNFS/TNFS but with reference GHG for
  !@+   radiative forcing calculations. TOA only.
  !@+   index 1=CH4, 2=N2O, 3=CFC11, 4=CFC12
  REAL*8, DIMENSION(4,grid%I_STRT_HALO:grid%I_STOP_HALO,            &
       grid%J_STRT_HALO:grid%J_STOP_HALO)              &
       ::  snfs_ghg, tnfs_ghg
  REAL*8, DIMENSION(4)  ::  sv_fulgas_ref, sv_fulgas_now
  INTEGER  ::  nf, GFrefY, GFrefD, GFnowY, GFnowD
  !@var nfghg fulgas( ) index of radf diag ghgs : 
  INTEGER, DIMENSION(4)  ::  nfghg = (/7,6,8,9/)
#endif
#endif

#if defined ( TRACERS_GC )
  !@var snfs_ghg,tnfs_ghg like SNFS/TNFS but with reference GHG for
  !@+   radiative forcing calculations. TOA only.
  !@+   index 1=CH4, 2=N2O, 3=CFC11, 4=CFC12
  REAL*8, DIMENSION(4,grid%I_STRT_HALO:grid%I_STOP_HALO,            &
                      grid%J_STRT_HALO:grid%J_STOP_HALO)            &
                        ::  snfs_ghg,    tnfs_ghg,                  &
                            snfs_ghg_tp, tnfs_ghg_tp                     
																																																																							
  REAL*8, DIMENSION(4)  ::  sv_fulgas_ref, sv_fulgas_now
  INTEGER               ::  nf, GFrefY, GFrefD, GFnowY, GFnowD
  !@var nfghg fulgas( ) index of radf diag ghgs : 
  INTEGER, DIMENSION(4) ::  nfghg = (/7,6,8,9/)
  ! For ozone
  REAL*8, DIMENSION(5,grid%I_STRT_HALO:grid%I_STOP_HALO,            &
                      grid%J_STRT_HALO:grid%J_STOP_HALO)            &
                                   :: SNFST_o3ref, TNFST_o3ref
  ! For 3D fluxes
		! 20 radiatively-active species (1:13 gases + 12:20 aerosol particle types + 21:21 clouds)
  REAL*8, DIMENSION( grid%I_STRT_HALO:grid%I_STOP_HALO,             &
                     grid%J_STRT_HALO:grid%J_STOP_HALO,             &
																					LM+LM_REQ+1, 21 ) :: SNFS_3D_pert, TNFS_3D_pert
  ! Baseline 3-D radiation fluxes 
  REAL*8, DIMENSION( grid%I_STRT_HALO:grid%I_STOP_HALO,             &
                     grid%J_STRT_HALO:grid%J_STOP_HALO,             &
  																			LM+LM_REQ+1 )     :: SNFS_3D,      TNFS_3D																																						

  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
                    grid%J_STRT_HALO:grid%J_STOP_HALO,LM+LM_REQ+1)  &
                                       ::  SDDARRFLX
#endif

  ! variables for running uncoupled concentration-driven GCC
#ifdef GCC_UNCOUPLE_RAD_CONCEN
  REAL*8  ::  GCCco2_fulgas_ref, GCCco2_fulgas_now
  INTEGER  ::  GCCco2nowY, GCCco2nowD
#endif

#ifdef HEALY_LM_DIAGS
  !  GHG Effective forcing relative to 1850
  REAL*8  ::  ghg_totforc, CO2I = 285.2, N2OI = .2754, CH4I = .791
  ! 1850  GHG's
  REAL*8  ::  CO2R = 337.9, N2OR = .3012, CH4R = 1.547     ! RAD's 1979 Reference values
  REAL*8  ::  FCO2, FN2O, FCH4                             ! Current Model GHG
  REAL*8  ::  FE
  !! Function
#endif
#ifdef TRACERS_ON
  !@var SNFST,TNFST like SNFS/TNFS but with/without specific tracers for
  !@+   radiative forcing calculations
  REAL*8, DIMENSION(2,nraero_rf,grid%I_STRT_HALO:grid%I_STOP_HALO,  &
       grid%J_STRT_HALO:grid%J_STOP_HALO)  ::  SNFST,  &
       TNFST
  !@var SNFST_o3ref,TNFST_o3ref like snfst,tnfst for special case ozone for
  !@+   which nraero_rf fields are not defined. Indicies are  : 
  !@+   1=LTROPO,reference, 2=TOA,reference; not saving surface forcing.
  !@+   3=LTROPO or LS1-1,auxiliary, 4=TOA,auxiliary; 5=LS1-1,reference
#if (defined SHINDELL_STRAT_EXTRA) &(defined ACCMIP_LIKE_DIAGS)
  REAL*8, DIMENSION(5,grid%I_STRT_HALO:grid%I_STOP_HALO,            &
       grid%J_STRT_HALO:grid%J_STOP_HALO)              &
       ::  SNFST_o3ref, TNFST_o3ref, snfst_stratOx,   &
       tnfst_stratOx
#endif /* SHINDELL_STRAT_EXTRA &ACCMIP_LIKE_DIAGS */
#ifdef BC_ALB
  REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,              &
       grid%J_STRT_HALO:grid%J_STOP_HALO)  ::  ALBNBC, &
       NFSNBC, dALBsnBC
  ! not to be confused with BCdalbsn from an input file
  LOGICAL, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,             &
       grid%J_STRT_HALO:grid%J_STOP_HALO)             &
       ::  bc_snow_present
#endif /* BC_ALB */
#endif /* TRACERS_ON */
  REAL*8, DIMENSION(LM_REQ,grid%I_STRT_HALO:grid%I_STOP_HALO,       &
       grid%J_STRT_HALO:grid%J_STOP_HALO)  ::  TRHRS,  &
       SRHRS
  REAL*8, DIMENSION(0 : LM+LM_REQ,grid%I_STRT_HALO:grid%I_STOP_HALO,  &
       grid%J_STRT_HALO:grid%J_STOP_HALO)  ::  TRHRA,  &
       SRHRA
  ! for adj.frc
  REAL*8, DIMENSION(LM)  ::  TOTCLD, SS_CLD, dcc_cdncl, dod_cdncl
  INTEGER I, J, L, K, KR, LR, JR, IH, IHM, INCH, JK, IT, iy, iend,  &
       N, onoff_aer, onoff_chem, LFRC, JTIME, n1, moddrf
  REAL*8 ROT1, ROT2, PLAND, CSS, CMC, DEPTH, QSS, TAUSSL, TAUSSLIP, &
       TAUMCL, ELHX, CLDCV, X, OPNSKY, CSZ2, tauup, taudn,        &
       ptype4(4), taucl, wtlin, MSTRAT, STRATQ, STRJ, MSTJ, optdw,&
       optdi, rsign_aer, rsign_chem, tauex5, tauex6, tausct,      &
       taugcb, dcdnc,                                             &
       QR(LM,grid%I_STRT_HALO:grid%I_STOP_HALO,grid%J_STRT_HALO :   &
       grid%J_STOP_HALO),                                         &
       CLDinfo(LM,3,grid%I_STRT_HALO:grid%I_STOP_HALO,            &
       grid%J_STRT_HALO:grid%J_STOP_HALO)
  REAL*8 tmpS(8), tmpT(8)
  REAL*8 QSAT
#ifdef BC_ALB
  REAL*8 dALBsn1
#endif
  LOGICAL set_clayilli, set_claykaol, set_claysmec, set_claycalc,   &
       set_clayquar
  !
  REAL*8 RDSS(LM,grid%I_STRT_HALO:grid%I_STOP_HALO,                 &
       grid%J_STRT_HALO:grid%J_STOP_HALO),                        &
       RDMC(grid%I_STRT_HALO:grid%I_STOP_HALO,                    &
       grid%J_STRT_HALO:grid%J_STOP_HALO)

  REAL*8  ::  TMP(NDIUVAR)
  INTEGER, PARAMETER  ::  NLOC_DIU_VAR = 8
  INTEGER  ::  idx(NLOC_DIU_VAR)
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
  INTEGER, PARAMETER  ::  NLOC_DIU_VARB = 5
#else
  INTEGER, PARAMETER  ::  NLOC_DIU_VARB = 3
#endif
  INTEGER  ::  idxb(NLOC_DIU_VARB)

  INTEGER  ::  aj_alb_inds(8)
  REAL*8, DIMENSION(lm_req)  ::  bydpreq

  !     INTEGER ICKERR,JCKERR,KCKERR
  INTEGER  ::  J_0, J_1, I_0, I_1
  INTEGER  ::  J_0S, J_1S
  LOGICAL  ::  HAVE_SOUTH_POLE, HAVE_NORTH_POLE
  CHARACTER(LEN=300)  ::  out_line

  INTEGER  ::  NIJ_BEFORE_J0, NIJ_AFTER_J1, NIJ_AFTER_I1
  INTEGER  ::  initial_GHG_setup

  REAL*8  ::  PVT0(N_COVERTYPES), HVT0(N_COVERTYPES)
#ifdef TRACERS_NITRATE
  REAL*8  ::  nh4_on_no3
#endif
#ifdef TRACERS_TOMAS
  REAL*8  ::  qcb_col(6,ICOMP-2), qcb_col_dry(6,ICOMP-2)
#endif

  REAL*8, DIMENSION( : , : ), POINTER  ::  RSI, ZSI, SNOWI, POND_MELT
  LOGICAL, DIMENSION( : , : ), POINTER  ::  FLAG_DSWS
  REAL*8  ::  rhodz
  ! air density times layer thickness (kg/m2
  INTEGER  ::  year, dayOfYear, hour, date

#ifdef TRACERS_ON
  !@var nsub_ntrix  array of index counters for sub classes of tracers
  INTEGER, DIMENSION(ntm)  ::  nsub_ntrix
#endif

#ifdef GCC_COUPLE_RAD
  INTEGER  ::  Lmax_rad_CO2 = LM
#endif

  CALL MODELECLOCK%GET(year=year,dayOfYear=dayOfYear,hour=hour,     &
       date=date)

  RSI => SI_ATM%RSI
  ZSI => SI_ATM%ZSI
  SNOWI => SI_ATM%SNOWI
  POND_MELT => SI_ATM%POND_MELT
  FLAG_DSWS => SI_ATM%FLAG_DSWS

  !
  !****
  CALL STARTTIMER('RADIA()')

  idx = (/(IDD_CL7+i-1,i=1,7),IDD_CCV/)
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
  idxb = (/IDD_PALB,IDD_GALB,IDD_ABSA,idd_aot,idd_aot2/)
#else
  idxb = (/IDD_PALB,IDD_GALB,IDD_ABSA/)
#endif
  CALL GETDOMAINBOUNDS(grid,HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,        &
       HAVE_NORTH_POLE=HAVE_NORTH_POLE)
  I_0 = grid%I_STRT
  I_1 = grid%I_STOP
  J_0 = grid%J_STRT
  J_1 = grid%J_STOP
  J_0S = grid%J_STRT_SKP
  J_1S = grid%J_STOP_SKP


  !****
  !**** FLAND     LAND COVERAGE (1)
  !**** FLICE     LAND ICE COVERAGE (1)
  !****
  !**** GTEMPR RADIATIVE TEMPERATURE ARRAY OVER ALL SURFACE TYPES (K)
  !****   RSI  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
  !****
  !**** VDATA  1-11 RATIOS FOR THE 11 VEGETATION TYPES (1)
  !****

  !**** limit optical cloud depth from below :  taulim
  taulim = MIN(tauwc0,tauic0)
  ! currently both .001
  tauwc0 = taulim
  tauic0 = taulim
  !**** Calculate mean cosine of zenith angle for the current physics step
  JTIME = MOD(ITIME,NDAY)
  ROT1 = (TWOPI*JTIME)/NDAY
  !      ROT2=ROT1+TWOPI*DTsrc/SECONDS_PER_DAY
  !      CALL COSZT (ROT1,ROT2,COSZ1)
  CALL CALC_ZENITH_ANGLE   ! moved to main loop

  IF ( kradia>0 ) THEN     ! read in all rad. input data (frc.runs)
     iend = 1
     it = itime - 1        ! make sure, at least 1 record is read
     DO WHILE ( MOD(itime-it,NDAY*INT_DAYS_PER_YEAR)/=0 )
        !****   input data :           WARNINGS
        !****        1 - any changes here also go in later (look for 'iu_rad')
        !****        2 - keep "dimrad_sv" up-to-date :          dimrad_sv=IM*JM*{
        ! LM+LM_REQ+1+
        !     *     ,(((GTEMPR(k,i,j),k=1,4),i=1,im),j=1,jm)    ! (4+)
        ! LM+1+3*LM+1+1+
        ! 1+1+1+1+1+
        ! 3+1+.5+.5+
        !****   output data :  really needed only if kradia=2
        ! 2+1+1
        !****   total :  dimrad_sv= IM*JM*(7*LM + 3*LM_REQ + 24 (+4)) => RAD_COM.f
        READ (iu_rad,END=10,ERR=10) it, T, RQT, atmsrf%TSAVG, QR, P,&
             CLDinfo, rsi, zsi, wsoil,       &
             atmsrf%WSAVG, snowi,            &
             atmgla%SNOW, atmlnd%SNOWE,      &
             snoage, fmp_com, flag_dsws,     &
             ltropo, atmlnd%FR_SNOW_RAD,     &
             dlake, flake, srhra, trhra, iy
        ! 2(LM+LM_REQ+1)}
        IF ( qcheck ) THEN
           WRITE (out_line,*) 'reading RADfile at Itime', Itime, it,&
                iy
           CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
        ENDIF
     ENDDO
     iend = 0
10   IF ( it/=iy .OR. iend==1 ) THEN
        WRITE (out_line,*) 'RAD input file bad or too short : ',      &
             itime, it, iy, iend
        CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
        CALL STOP_MODEL('RADIA :  input file bad or too short',255)
     ENDIF
  ENDIF

  IF ( MODRD==0 ) THEN
     IDACC(ia_rad) = IDACC(ia_rad) + 1
     moddrf = 1
     ! skip rad.forcing diags if nradfrc.le.0
     IF ( nradfrc>0 ) moddrf = MOD(itime-itimei,nrad*nradfrc)
     !****
     IF ( moddrf==0 ) IDACC(ia_rad_frc) = IDACC(ia_rad_frc) + 1
     !**** Interface with radiation routines, done only every NRAD time steps
     !****
     !**** Calculate mean cosine of zenith angle for the full radiation step
     ROT2 = ROT1 + TWOPI*NRAD*DTsrc/SECONDS_PER_DAY
     CALL COSZS(ROT1,ROT2,COSZ2,COSZA)
#ifdef GCAP
     save_COSZ2 = COSZ2
#endif
     JDAYR = dayOfYear
     JYEARR = YEAR

     IF ( IS_SET_PARAM('s0') ) THEN
        ! typically only used for SCM
        CALL GET_PARAM('s0',s0)
        s00wm2 = s0
        ! just in case
     ELSE
        S0 = S0X*S00WM2*RATLS0/RSDIST
     ENDIF

#ifdef OLD_BCdalbsn
     !**** find scaling factors for surface albedo reduction
     ! LTM :  Fix, testing equality of reals is not reliable
     IF ( dalbsnX/=0 ) THEN
        IF ( HAVE_SOUTH_POLE ) THEN
           sumda_psum( : ,1) = AXYP(1,1)
           tauda_psum( : ,1) = AXYP(1,1)*DEPOBC_1990(1,1)
        ENDIF
        DO j = J_0S, J_1S
           DO i = I_0, I_1
              ! ilon72, jlat46 are indices w.r.t 72x46 grid
              !      JLAT46=INT(1.+(J-1.)*0.25*DLAT_DG+.5)   ! slightly more general
              !      ILON72=INT(.5+(I-.5)*72./IM+.5)
              ilon72 = 1 + INT(72D0*LON2D(i,j)/twopi)
              jlat46 = 1 + INT(45D0*(LAT2D(i,j)+92D0*radian)/pi)
              fsnow = FLICE(i,j) + rsi(i,j)*(1-FLAND(i,j))
              IF ( atmlnd%SNOWE(I,J)>0. ) fsnow = fsnow +           &
                   FEARTH(i,j)
              sumda_psum(i,j) = AXYP(i,j)*fsnow
              tauda_psum(i,j) = AXYP(i,j)*fsnow*DEPOBC_1990(i,j)
           ENDDO
        ENDDO
        IF ( HAVE_NORTH_POLE ) THEN
           sumda_psum( : ,JM) = AXYP(1,jm)*rsi(1,jm)
           tauda_psum( : ,JM) = AXYP(1,jm)*rsi(1,jm)*DEPOBC_1990(1,jm)
        ENDIF
        CALL GLOBALSUM(grid,sumda_psum,sumda,ALL=.TRUE.)
        CALL GLOBALSUM(grid,tauda_psum,tauda,ALL=.TRUE.)

        xdalbs = -dalbsnX*sumda/tauda
        IF ( QCHECK ) WRITE (6,*) 'coeff. for snow alb reduction',  &
             xdalbs
     ENDIF
     ! dalbsnX not zero
#endif

     IF ( kradia<=0 ) THEN
        IF ( QCHECK ) THEN
           !****   Calculate mean strat water conc
           STRATQ = 0.
           MSTRAT = 0.
           DO J = J_0, J_1
              STRJ = 0.
              MSTJ = 0.
              DO I = I_0, IMAXJ(J)
                 DO L = LTROPO(I,J) + 1, LM
                    STRJ = STRJ + Q(I,J,L)*MA(L,I,J)*AXYP(I,J)
                    MSTJ = MSTJ + MA(L,I,J)*AXYP(I,J)
                 ENDDO
              ENDDO
              IF ( J==1 .OR. J==JM ) THEN
                 STRJ = STRJ*IM
                 MSTJ = MSTJ*IM
              ENDIF
              STRATQ = STRATQ + STRJ
              MSTRAT = MSTRAT + MSTJ
           ENDDO
           PRINT *, "Strat water vapour (ppmv), mass (mb)",         &
                1D6*STRATQ*mair/(18.*MSTRAT),                      &
                PMTOP + 1D-2*GRAV*MSTRAT/AREAG
        ENDIF

        !**** Get the random numbers outside openMP parallel regions
        !**** but keep MC calculation separate from SS clouds
        !**** To get parallel consistency also with mpi, force each process
        !**** to generate random numbers for all latitudes (using BURN_RANDOM)

        !**** MC clouds are considered as a block for each I,J grid point

        CALL BURN_RANDOM(NIJ_BEFORE_J0(J_0))

        DO J = J_0, J_1           ! complete overlap
           CALL BURN_RANDOM((I_0-1))
           DO I = I_0, IMAXJ(J)
              RDMC(I,J) = RANDU(X)
              ! 1 random number per column
           ENDDO
           CALL BURN_RANDOM(NIJ_AFTER_I1(I_1))
        ENDDO

        CALL BURN_RANDOM((NIJ_AFTER_J1(J_1)))

        !**** SS clouds are considered as a block for each continuous cloud
        CALL BURN_RANDOM(NIJ_BEFORE_J0(j_0)*LM)

        DO J = J_0, J_1           ! semi-random overlap
           CALL BURN_RANDOM((I_0-1)*LM)
           DO I = I_0, IMAXJ(J)
              ! reverse loop kept only for consistency with previous version
              DO L = LM, 1, -1
                 !   better :   1,LM
                 IF ( TAUSS(L,I,J)<=taulim ) CLDSS(L,I,J) = 0.
                 IF ( TAUMC(L,I,J)<=taulim ) CLDMC(L,I,J) = 0.
                 RDSS(L,I,J) = RANDU(X)
              ENDDO
           ENDDO
           CALL BURN_RANDOM(NIJ_AFTER_I1(I_1)*LM)
        ENDDO

        CALL BURN_RANDOM(NIJ_AFTER_J1(j_1)*LM)

     ENDIF                  ! kradia le 0

#if (defined ACCMIP_LIKE_DIAGS)
#ifndef SKIP_ACCMIP_GHG_RADF_DIAGS
     ! because of additional updghg calls, these factors will not apply : 
     ! LTM :  This needs to be fixed, testing equality of reals is not reliable
     IF ( CO2X/=1. ) CALL STOP_MODEL('CO2x.ne.1 accmip diags',255)
     IF ( N2OX/=1. ) CALL STOP_MODEL('N2Ox.ne.1 accmip diags',255)
     IF ( CH4X/=1. ) CALL STOP_MODEL('CH4x.ne.1 accmip diags',255)
     IF ( CFC11X/=1. ) CALL STOP_MODEL('CFC11x.ne.1 accmip diags',  &
          255)
     IF ( CFC12X/=1. ) CALL STOP_MODEL('CFC12x.ne.1 accmip diags',  &
          255)
     IF ( XGHGX/=1. ) CALL STOP_MODEL('XGHGx.ne.1 accmip diags',255)
     GFrefY = 1850
     GFrefD = 182             ! ghg forcing refrnce year, day
     GFnowY = JyearR
     GFnowD = JdayR           ! ghg current desired year, day
     IF ( KJDAYG>0 ) GFnowD = KJDAYG
     ! unless presribed in deck
     IF ( KYEARG>0 ) GFnowY = KYEARG
     !
     CALL UPDGHG(GFrefY,GFrefD)
     sv_fulgas_ref(1:4) = FULGAS(nfghg(1:4))
     CALL UPDGHG(GFnowY,GFnowD)
     sv_fulgas_now(1:4) = FULGAS(nfghg(1:4))
#endif
#endif

#if defined ( TRACERS_GC )

     ! because of additional updghg calls, these factors will not apply : 
     ! LTM :  This needs to be fixed, testing equality of reals is not reliable
     IF ( CO2X/=1. ) CALL STOP_MODEL('CO2x.ne.1 accmip diags',255)
     IF ( N2OX/=1. ) CALL STOP_MODEL('N2Ox.ne.1 accmip diags',255)
     IF ( CH4X/=1. ) CALL STOP_MODEL('CH4x.ne.1 accmip diags',255)
     IF ( CFC11X/=1. ) CALL STOP_MODEL('CFC11x.ne.1 accmip diags',  &
          255)
     IF ( CFC12X/=1. ) CALL STOP_MODEL('CFC12x.ne.1 accmip diags',  &
          255)
     IF ( XGHGX/=1. ) CALL STOP_MODEL('XGHGx.ne.1 accmip diags',255)
     GFrefY = 1850
     GFrefD = 182             ! ghg forcing refrnce year, day
     GFnowY = JyearR
     GFnowD = JdayR           ! ghg current desired year, day
     IF ( KJDAYG>0 ) GFnowD = KJDAYG
     ! unless presribed in deck
     IF ( KYEARG>0 ) GFnowY = KYEARG
     !

     CALL UPDGHG(GFrefY,GFrefD)
     sv_fulgas_ref(1:4) = FULGAS(nfghg(1:4))

     CALL UPDGHG(GFnowY,GFnowD)
     sv_fulgas_now(1:4) = FULGAS(nfghg(1:4))

#endif



     ! Set variables used in storing reference CO2 for uncoupled runs
#ifdef GCC_UNCOUPLE_RAD_CONCEN
     IF ( KJDAYG>0 ) GCCco2nowD = KJDAYG
     ! unless presribed in deck
     IF ( KYEARG>0 ) GCCco2nowY = KYEARG
     !
     CALL UPDGHG(1850,182)
     GCCco2_fulgas_ref = FULGAS(2)
     CALL UPDGHG(GCCco2nowD,GCCco2nowD)
     GCCco2_fulgas_now = FULGAS(2)
#endif

#ifdef HEALY_LM_DIAGS
     FCO2 = FULGAS(2)*CO2R
     FN2O = FULGAS(6)*N2OR
     FCH4 = FULGAS(7)*CH4R
     !
     !      write(6,*) 'RJH :  GHG :  CONC=',
     !     * FCO2,FN2O,FCH4
     ghg_totforc = 5.35D0*LOG(FCO2/CO2I)                            &
          + 0.036D0*(SQRT(FCH4)-SQRT(CH4I))                &
          - (FE(FCH4,N2OI)-FE(CH4I,N2OI))                  &
          + 0.12D0*(SQRT(FN2O)-SQRT(N2OI))                 &
          - (FE(CH4I,FN2O)-FE(CH4I,N2OI))
     !      write(6,*) 'RJH :  GHG :  FORC=',ghg_totforc
#endif

     aj_alb_inds = (/J_PLAVIS,J_PLANIR,J_ALBVIS,J_ALBNIR,J_SRRVIS,  &
          J_SRRNIR,J_SRAVIS,J_SRANIR/)

     cfrac = 0.
     wtrcld = 0.
     icecld = 0.
     tausumw = 0.
     tausumi = 0.
     ctp = 0.
     ctt = 0.
     swus = 0.
     wtrtau = 0.
     icetau = 0.
#ifdef CFMIP3_SUBDD
     swut = 0.
     swutcs = 0.
     cfmip_twp = 0.
     swdcls = 0.
     swucls = 0.
     swdt = 0.
     cfmip_cf = 0.
     cfmip_qci = 0.
     cfmip_qcl = 0.
#endif
#ifdef GCAP
     tauw3d = 0.
     taui3d = 0.
#endif
#ifdef TRACERS_GC
     save_rf = 0.
					save_rf_tp = 0.
					save_rf_3D = 0.
					SNFS_3D_pert = 0.
					TNFS_3D_pert = 0.
     SNFS_3D = 0.
					TNFS_3D = 0.
#endif
					
     !****
     !**** MAIN J LOOP
     !****
     DO J = J_0, J_1

        !     ICKERR=0
        !     JCKERR=0
        !     KCKERR=0

        !****
        !**** MAIN I LOOP
        !****
        DO I = I_0, IMAXJ(J)
           !**** Radiation input files use a 72x46 grid independent of IM and JM
           !**** (ilon72,jlat46) is the 4x5 box containing the center of box (i,j)
           !      JLAT46=INT(1.+(J-1.)*45./(JM-1.)+.5)  !  lat_index w.r.to 72x46 grid
           !      JLAT46=INT(1.+(J-1.)*0.25*DLAT_DG+.5) ! slightly more general
           !      ILON72=INT(.5+(I-.5)*72./IM+.5)  ! lon_index w.r.to 72x46 grid
           igcm = i
           jgcm = j
           ilon72 = 1 + INT(72D0*LON2D(i,j)/twopi)
           jlat46 = 1 + INT(45D0*(LAT2D(i,j)+92D0*radian)/pi)
#ifdef ALTER_RADF_BY_LAT
           FULGAS( : ) = FULGAS_orig( : )*FULGAS_lat( : ,JLAT46)
           FS8OPX( : ) = FS8OPX_orig( : )*FS8OPX_lat( : ,JLAT46)
           FT8OPX( : ) = FT8OPX_orig( : )*FT8OPX_lat( : ,JLAT46)
#endif
           L1 = 1                ! lowest layer above ground
           LMR = LM + LM_REQ     ! radiation allows var. # of layers
           JR = JREG(I,J)
           !**** DETERMINE FRACTIONS FOR SURFACE TYPES AND COLUMN PRESSURE
           PLAND = FLAND(I,J)
           POICE = RSI(I,J)*(1.-PLAND)
           POCEAN = (1.-PLAND) - POICE
           PLAKE = FLAKE(I,J)
           PLICE = FLICE(I,J)
           PEARTH = FEARTH(I,J)
           ptype4(1) = pocean
           ! open ocean and open lake
           ptype4(2) = poice
           ! ocean/lake ice
           ptype4(3) = plice
           ! glacial ice
           ptype4(4) = pearth
           ! non glacial ice covered soil

           !**** CHECK SURFACE TEMPERATURES
           DO IT = 1, 4
              IF ( ptype4(IT)>0. ) THEN
                 !CC         STOP 'In Radia :  Grnd Temp out of range'
                 !           ICKERR=ICKERR+1
                 IF ( INT(ASFLX4(it)%GTEMPR(I,J))<planck_tmin .OR.  &
                      INT(ASFLX4(it)%GTEMPR(I,J))>=planck_tmax )    &
                      WRITE (6,*) 'In Radia :  Time,I,J,IT,TG1',      &
                      ITime, I, J, IT, ASFLX4(it)       &
                      %GTEMPR(I,J)
              ENDIF
           ENDDO

           !**** Set Chlorophyll concentration
           IF ( POCEAN>0 ) THEN
              IF ( (chl_from_seawifs>0 .OR. chl_from_obio>0) .AND.  &
                   atmocn%CHL_DEFINED ) THEN
                 LOC_CHL = atmocn%CHL(I,J)
                 IF ( ij_chl>0 ) AIJ(I,J,IJ_CHL) = AIJ(I,J,IJ_CHL)  &
                      + atmocn%CHL(I,J)*FOCEAN(I,J)
                 !         write(*,'(a,3i5,e12.4)')'RAD_DRV : ',
                 !    .    itime,i,j,chl(i,j)
              ELSE
                 LOC_CHL = -1.D30
              ENDIF
           ENDIF

           LS1_loc = LTROPO(I,J) + 1
           ! define stratosphere for radiation
           !**** kradia>1 :  adjusted forcing, i.e. T adjusts in L=LS1_loc->LM+3
           IF ( kradia>1 ) LS1_loc = LS1_loc + 2 - kradia
           ! favorite : kradia=3
           IF ( kradia>3 ) LS1_loc = 1      ! favorite : kradia=3
           kdeliq = 0
           ! initialize mainly for L>LM
           IF ( kradia>0 ) THEN
              ! rad forcing model
              DO l = 1, lm
                 TLM(l) = T(i,j,l)*PK(l,i,j)
                 SHL(l) = QR(l,i,j)
                 IF ( SHL(l)<0 ) SHL(l) = 0
                 TAUWC(l) = cldx*CLDinfo(l,1,i,j)
                 TAUIC(l) = cldx*CLDinfo(l,2,i,j)
                 SIZEWC(L) = CLDinfo(l,3,i,j)
                 SIZEIC(L) = SIZEWC(L)
              ENDDO
           ELSE             ! full model
              !****
              !**** DETERMINE CLOUDS (AND THEIR OPTICAL DEPTHS) SEEN BY RADIATION
              !****
              CSS = 0.
              CMC = 0.
              CLDCV = 0.
              DEPTH = 0.
              OPTDW = 0.
              OPTDI = 0.
              ! LTM :  Fix, testing equality of reals is not reliable.
              IF ( cc_cdncx/=0. .OR. od_cdncx/=0. ) THEN
                 CALL DCDNC_EST(i,j,pland,dCDNC)
              ELSE
                 dCDNC = 0.
              ENDIF
              dCC_CDNCL = CC_cdncx*dCDNC*CDNCL
              dOD_CDNCL = OD_cdncx*dCDNC*CDNCL

              !**** Adjust RDSS for semi-random overlap
              CALL GET_CLD_OVERLAP(lm,CLDSS( : ,i,j),                 &
                   RANDSS=rdss( : ,i,j))

              DO L = 1, LM
                 IF ( Q(i,j,l)<0 ) THEN
                    WRITE (6,*) 'In Radia :  Time,I,J,L,Q<0', ITime,  &
                         I, J, L, Q, '->0'
                    Q(I,J,L) = 0.
                 ENDIF
                 QSS = Q(I,J,L)/(RHSAV(L,I,J)+1.D-20)
                 SHL(L) = QSS
                 IF ( FSS(L,I,J)*CLDSAV(L,I,J)<1. ) SHL(L)          &
                      = (Q(I,J,L)-QSS*FSS(L,I,J)*CLDSAV(L,I,J))     &
                      /(1.-FSS(L,I,J)*CLDSAV(L,I,J))
                 TLM(L) = T(I,J,L)*PK(L,I,J)
                 rhodz = PDSIG(l,i,j)*100/grav
                 TAUSSL = 0.
                 TAUSSLIP = 0.
                 TAUMCL = 0.
                 TAUWC(L) = 0.
                 TAUIC(L) = 0.
                 SIZEWC(L) = 0.
                 SIZEIC(L) = 0.
                 TOTCLD(L) = 0.
                 SS_CLD(L) = 0.
                 !**** Determine large scale and moist convective cloud cover for radia
                 IF ( CLDSS(L,I,J)*(1.+dcc_cdncl(l))>RDSS(L,I,J) )  &
                      THEN
                    TAUSSL = TAUSS(L,I,J)*(1.+dod_cdncl(l))
                    ! tausslip is tau of ice precip in a supercooled water cloud
                    TAUSSLIP = TAUSSIP(L,I,J)*(1.+dod_cdncl(l))
                    SHL(L) = QSS
                    CSS = 1.
                    CALL INC_AJL(i,j,l,jl_sscld,css)
#ifdef CFMIP3_SUBDD
                    ! LS Cloud
                    cfmip_cf(i,j,l) = cfmip_cf(i,j,l) + 1.
#endif
                 ENDIF
                 IF ( CLDMC(L,I,J)>RDMC(I,J) ) THEN
                    CMC = 1.
                    CALL INC_AJL(i,j,l,jl_mccld,cmc)
#ifdef CFMIP3_SUBDD
                    ! MC Cloud
                    cfmip_cf(i,j,l) = MIN(cfmip_cf(i,j,l)+1.,1.)
#endif
                    DEPTH = DEPTH + PDSIG(L,I,J)
                    IF ( TAUMC(L,I,J)>TAUSSL+TAUSSLIP ) THEN
                       TAUMCL = TAUMC(L,I,J)
                       ELHX = LHE
                       IF ( TLM(L)<=TF ) ELHX = LHS
                       SHL(L) = QSAT(TLM(L),ELHX,PMID(L,I,J))
                    ENDIF
                 ENDIF
                 IF ( TAUSSL+TAUSSLIP+TAUMCL>0. ) THEN
                    CLDCV = 1.
                    TOTCLD(L) = 1.
                    CALL INC_AJL(i,j,l,jl_totcld,1D0)
                    !**** save 3D cloud fraction as seen by radiation
                    IF ( cldx>0 ) AIJL(I,J,L,IJL_CF)                &
                         = AIJL(I,J,L,IJL_CF) + 1.
                    IF ( TAUMCL>TAUSSL+TAUSSLIP ) THEN
                       SIZEWC(L) = CSIZMC(L,I,J)
                       SIZEIC(L) = CSIZMC(L,I,J)
                       IF ( SVLAT(L,I,J)==LHE ) THEN
                          TAUWC(L) = cldx*TAUMCL
                          OPTDW = OPTDW + TAUWC(L)
#ifdef GCAP
                          TAUW3D(I,J,L) = TAUW3D(I,J,L) + TAUMCL
                          ! in-cloud vs. in-cell TAUWC(L)
#endif
                          CALL INC_AJL(i,j,l,jl_wcld,1D0)
                          CALL INC_AJL(i,j,l,jl_wcldwt,PDSIG(l,i,j))
                          AIJ(i,j,ij_lwprad) = AIJ(i,j,ij_lwprad)   &
                               + QLMC(l,i,j)*rhodz/CLDMC(l,i,j)
                          AIJL(i,j,l,ijl_QLrad)                     &
                               = AIJL(i,j,l,ijl_QLrad) + QLMC(l,i,j)  &
                               *PDSIG(l,i,j)/CLDMC(l,i,j)
#ifdef CFMIP3_SUBDD
                          ! MC Cloud Liquid
                          cfmip_twp(i,j) = cfmip_twp(i,j)           &
                               + QLMC(l,i,j)*rhodz/CLDMC(l,i,j)
                          cfmip_qcl(i,j,l) = QLMC(l,i,j)            &
                               /CLDMC(l,i,j)
#endif
                       ELSE
                          TAUIC(L) = cldx*TAUMCL
                          OPTDI = OPTDI + TAUIC(L)
#ifdef GCAP
                          TAUI3D(I,J,L) = TAUI3D(I,J,L) + TAUMCL
                          ! in-cloud vs. in-cell TAUIC(L)
#endif
                          CALL INC_AJL(i,j,l,jl_icld,1D0)
                          CALL INC_AJL(i,j,l,jl_icldwt,PDSIG(l,i,j))
                          AIJ(i,j,ij_iwprad) = AIJ(i,j,ij_iwprad)   &
                               + QIMC(l,i,j)*rhodz/CLDMC(l,i,j)
                          AIJL(i,j,l,ijl_QIrad)                     &
                               = AIJL(i,j,l,ijl_QIrad) + QIMC(l,i,j)  &
                               *PDSIG(l,i,j)/CLDMC(l,i,j)
#ifdef CFMIP3_SUBDD
                          ! MC Cloud Ice
                          cfmip_twp(i,j) = cfmip_twp(i,j)           &
                               + QIMC(l,i,j)*rhodz/CLDMC(l,i,j)
                          cfmip_qci(i,j,l) = QIMC(l,i,j)            &
                               /CLDMC(l,i,j)
#endif
                       ENDIF
                    ELSE
                       SS_CLD(L) = 1.
                       SIZEWC(L) = CSIZSS(L,I,J)
                       SIZEIC(L) = CSIZSS(L,I,J)
                       IF ( SVLHX(L,I,J)==LHE ) THEN
                          TAUWC(L) = cldx*TAUSSL
                          OPTDW = OPTDW + TAUWC(L)
#ifdef GCAP
                          TAUW3D(I,J,L) = TAUW3D(I,J,L) + TAUSSL
                          ! in-cloud vs. in-cell TAUWC(L)
#endif
                          CALL INC_AJL(i,j,l,jl_wcld,1D0)
                          CALL INC_AJL(i,j,l,jl_wcldwt,PDSIG(l,i,j))
                          AIJ(i,j,ij_lwprad) = AIJ(i,j,ij_lwprad)   &
                               + QLSS(l,i,j)*rhodz/CLDSS(l,i,j)
                          AIJL(i,j,l,ijl_QLrad)                     &
                               = AIJL(i,j,l,ijl_QLrad) + QLSS(l,i,j)  &
                               *PDSIG(l,i,j)/CLDSS(l,i,j)
#ifdef CFMIP3_SUBDD
                          ! LS Cloud Liquid
                          cfmip_twp(i,j) = cfmip_twp(i,j)           &
                               + QLSS(l,i,j)*rhodz/CLDSS(l,i,j)
                          cfmip_qcl(i,j,l) = QLSS(l,i,j)            &
                               /CLDSS(l,i,j)
#endif
                          IF ( tausslip>0. ) THEN
                             SIZEIC(L) = CSIZSSIP(L,I,J)
                             TAUIC(L) = cldx*TAUSSLIP
                             OPTDI = OPTDI + TAUIC(L)
#ifdef GCAP
                             TAUI3D(I,J,L) = TAUI3D(I,J,L)          &
                                  + TAUSSLIP         ! in-cloud vs. in-cell TAUIC(L)
#endif
                             CALL INC_AJL(i,j,l,jl_icld,1D0)
                             CALL INC_AJL(i,j,l,jl_icldwt,          &
                                  PDSIG(l,i,j))
                             AIJ(i,j,ij_iwprad) = AIJ(i,j,ij_iwprad)&
                                  + QISS(l,i,j)*rhodz/CLDSS(l,i,j)
                             AIJL(i,j,l,ijl_QIrad)                  &
                                  = AIJL(i,j,l,ijl_QIrad)             &
                                  + QISS(l,i,j)*PDSIG(l,i,j)          &
                                  /CLDSS(l,i,j)
#ifdef CFMIP3_SUBDD
                             ! LS Snow in supercooled liquid
                             cfmip_twp(i,j) = cfmip_twp(i,j)        &
                                  + QISS(l,i,j)*rhodz/CLDSS(l,i,j)
#endif
                          ENDIF
                       ELSE
                          TAUIC(L) = cldx*TAUSSL
                          OPTDI = OPTDI + TAUIC(L)
#ifdef GCAP
                          TAUI3D(I,J,L) = TAUI3D(I,J,L) + TAUSSL
                          ! in-cloud vs. in-cell TAUIC(L)
#endif
                          CALL INC_AJL(i,j,l,jl_icld,1D0)
                          CALL INC_AJL(i,j,l,jl_icldwt,PDSIG(l,i,j))
                          AIJ(i,j,ij_iwprad) = AIJ(i,j,ij_iwprad)   &
                               + QISS(l,i,j)*rhodz/CLDSS(l,i,j)
                          AIJL(i,j,l,ijl_QIrad)                     &
                               = AIJL(i,j,l,ijl_QIrad) + QISS(l,i,j)  &
                               *PDSIG(l,i,j)/CLDSS(l,i,j)
#ifdef CFMIP3_SUBDD
                          ! LS Cloud Ice
                          cfmip_twp(i,j) = cfmip_twp(i,j)           &
                               + QISS(l,i,j)*rhodz/CLDSS(l,i,j)
                          cfmip_qci(i,j,l) = QISS(l,i,j)            &
                               /CLDSS(l,i,j)
#endif
                       ENDIF
                    ENDIF
                    CALL INC_AJL(i,j,l,jl_wcod,TAUWC(l))
                    CALL INC_AJL(i,j,l,jl_icod,TAUIC(l))
                    CALL INC_AJL(i,j,l,jl_wcsiz,SIZEWC(l)*TAUWC(l))
                    CALL INC_AJL(i,j,l,jl_icsiz,SIZEIC(l)*TAUIC(l))
                    AIJL(i,j,l,ijl_wtrtau) = AIJL(i,j,l,ijl_wtrtau) &
                         + TAUWC(l)
                    AIJL(i,j,l,ijl_icetau) = AIJL(i,j,l,ijl_icetau) &
                         + TAUIC(l)
                    wtrtau(i,j,l) = TAUWC(l)
                    icetau(i,j,l) = TAUIC(l)
                 ENDIF
                 !**** save some radiation/cloud fields for wider use
                 RCLD(L,I,J) = TAUWC(L) + TAUIC(L)
              ENDDO
              CFRAC(I,J) = CLDCV
              ! cloud fraction consistent with radiation
              !**** effective cloud cover diagnostics
              OPNSKY = 1. - CLDCV
              DO IT = 1, NTYPE
                 CALL INC_AJ(i,j,it,J_PCLDSS,CSS*FTYPE(IT,I,J))
                 CALL INC_AJ(i,j,it,J_PCLDMC,CMC*FTYPE(IT,I,J))
                 CALL INC_AJ(i,j,it,J_CLDDEP,DEPTH*FTYPE(IT,I,J))
                 CALL INC_AJ(i,j,it,J_PCLD,CLDCV*FTYPE(IT,I,J))
              ENDDO
              CALL INC_AREG(i,j,jr,J_PCLDSS,CSS)
              CALL INC_AREG(i,j,jr,J_PCLDMC,CMC)
              CALL INC_AREG(i,j,jr,J_CLDDEP,DEPTH)
              CALL INC_AREG(i,j,jr,J_PCLD,CLDCV)
              AIJ(I,J,IJ_PMCCLD) = AIJ(I,J,IJ_PMCCLD) + CMC
              AIJ(I,J,IJ_CLDCV) = AIJ(I,J,IJ_CLDCV) + CLDCV
              DO L = 1, LLOW
                 ! LTM :  Fix, testing equality of reals is not reliable
                 IF ( TOTCLD(L)/=1. ) CYCLE
                 AIJ(I,J,IJ_PCLDL) = AIJ(I,J,IJ_PCLDL) + 1.
                 EXIT
              ENDDO
              DO L = LLOW + 1, LMID
                 ! LTM :  Fix, testing equality of reals is not reliable
                 IF ( TOTCLD(L)/=1. ) CYCLE
                 AIJ(I,J,IJ_PCLDM) = AIJ(I,J,IJ_PCLDM) + 1.
                 EXIT
              ENDDO
              DO L = LMID + 1, LHI
                 ! LTM :  Fix, testing equality of reals is not reliable
                 IF ( TOTCLD(L)/=1. ) CYCLE
                 AIJ(I,J,IJ_PCLDH) = AIJ(I,J,IJ_PCLDH) + 1.
                 EXIT
              ENDDO
              DO L = 1, LLOW
                 ! LTM :  Fix, testing equality of reals is not reliable
                 IF ( SS_CLD(L)/=1. ) CYCLE
                 AIJ(I,J,IJ_PCLDL_SS) = AIJ(I,J,IJ_PCLDL_SS) + 1.
                 EXIT
              ENDDO

              TAUSUMW(I,J) = OPTDW
              TAUSUMI(I,J) = OPTDI
              IF ( optdw>0. ) THEN
                 AIJ(I,J,IJ_optdw) = AIJ(I,J,IJ_optdw) + optdw
                 AIJ(I,J,IJ_wtrcld) = AIJ(I,J,IJ_wtrcld) + 1.
                 WTRCLD(I,J) = 1.
              ENDIF
              IF ( optdi>0. ) THEN
                 AIJ(I,J,IJ_optdi) = AIJ(I,J,IJ_optdi) + optdi
                 AIJ(I,J,IJ_icecld) = AIJ(I,J,IJ_icecld) + 1.
                 ICECLD(I,J) = 1.
              ENDIF

              DO KR = 1, NDIUPT
                 IF ( I==IJDD(1,KR) .AND. J==IJDD(2,KR) ) THEN
                    !**** Warning :  this replication may give inaccurate results for hours
                    !****          1->(NRAD-1)*DTsrc (ADIURN) or skip them (HDIURN)
                    TMP(IDD_CL7 : IDD_CL7+6) = TOTCLD(1 : 7)
                    TMP(IDD_CCV) = CLDCV
                    DO INCH = 1, NRAD
                       IHM = 1 + (JTIME+INCH-1)*HR_IN_DAY/NDAY
                       IH = IHM
                       IF ( IH>HR_IN_DAY ) IH = IH - HR_IN_DAY
                       ADIURN(IDX( : ),KR,IH) = ADIURN(IDX( : ),KR,IH)  &
                            + TMP(IDX( : ))
#ifdef USE_HDIURN
                       IHM = IHM + (DATE-1)*HR_IN_DAY
                       IF ( IHM<=HR_IN_MONTH ) HDIURN(IDX( : ),KR,IHM)&
                            = HDIURN(IDX( : ),KR,IHM) + TMP(IDX( : ))
#endif
                    ENDDO
                 ENDIF
              ENDDO
           ENDIF
           ! kradia le 0 (full model)
           !****
           !**** SET UP VERTICAL ARRAYS OMITTING THE I AND J INDICES
           !****
           !**** EVEN PRESSURES
#ifdef TRACERS_TOMAS
           aesqex( : , : , : ) = 0.0
           aesqsc( : , : , : ) = 0.0
           aesqcb( : , : , : ) = 0.0
           aesqex_dry( : , : , : ) = 0.0
           aesqsc_dry( : , : , : ) = 0.0
           aesqcb_dry( : , : , : ) = 0.0
#endif
           PLB(LM+1) = PEDN(LM+1,I,J)
           DO L = 1, LM
              PLB(L) = PEDN(L,I,J)
              !**** TEMPERATURES
              !---- TLm(L)=T(I,J,L)*PK(L,I,J)     ! already defined
              IF ( INT(TLM(L))<planck_tmin .OR. INT(TLM(L))         &
                   >=planck_tmax ) THEN
                 WRITE (6,*) 'In Radia :  Time,I,J,L,TL', ITime, I, J,&
                      L, TLM(L)
                 WRITE (6,*) 'GTEMPR : ', ASFLX4(1)%GTEMPR(I,J),      &
                      ASFLX4(2)%GTEMPR(I,J), ASFLX4(3)       &
                      %GTEMPR(I,J), ASFLX4(4)%GTEMPR(I,J)
                 !CC       STOP 'In Radia :  Temperature out of range'
                 !         ICKERR=ICKERR+1
              ENDIF
              !**** MOISTURE VARIABLES
              !---- shl(L)=Q(I,J,L)        ! already defined and reset to 0 if <0
              !       if(shl(l).lt.0.) then
              !         WRITE(0,*)'In Radia :  Time,I,J,L,QL<0',ITime,I,J,L,shl(L),'->0'
              !         KCKERR=KCKERR+1
              !         shl(l)=0.
              !       end if
              RHL(L) = SHL(L)/QSAT(TLM(L),LHE,PMID(L,I,J))
              IF ( RHfix>=0. ) RHL(L) = RHfix
              !**** Extra aerosol data
              !**** For up to nraero_aod aerosols, define the aerosol amount to
              !**** be used (kg/m^2)
              !**** Only define TRACER if individual tracer is actually defined.
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
              (defined TRACERS_MINERALS) || (defined TRACERS_AEROSOLS_SEASALT)
              !**** loop over tracers that are passed to radiation.
              !**** Some special cases for black carbon, organic carbon, SOAs where
              !**** more than one tracer is lumped together for radiation purposes
              DO n = 1, nraero_aod
                 SELECT CASE (TRNAME(NTRIX_AOD(n)))
                 CASE ("OCIA","vbsAm2")
#ifdef TRACERS_AEROSOLS_VBS
                    TRACER(L,n) = SUM(TRM(i,j,l,vbs_tr%IAER))
#else
                    TRACER(L,n) = TRM(i,j,l,n_OCII)                 &
                         + TRM(i,j,l,n_OCIA)
#endif  /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_OCEAN
                    TRACER(L,n) = TRACER(L,n) + TRM(i,j,l,n_ococean)
#endif  /* TRACERS_AEROSOLS_OCEAN */
                    TRACER(L,n) = TRACER(L,n)*BYAXYP(I,J)
                 CASE ("OCB")
#ifdef TRACERS_AEROSOLS_VBS
                    TRACER(L,n) = 0.D0
#else
                    TRACER(L,n) = TRM(i,j,l,n_OCB)*BYAXYP(I,J)
#endif  /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_SOA
                 CASE ("isopp1a")
                    TRACER(L,n) = TRM(i,j,l,n_isopp1a)              &
                         + TRM(i,j,l,n_isopp2a)
#ifdef TRACERS_TERP
                    TRACER(L,n) = TRACER(L,n) + TRM(i,j,l,n_apinp1a)&
                         + TRM(i,j,l,n_apinp2a)
#endif /* TRACERS_TERP */
                    TRACER(L,n) = TRACER(L,n)*BYAXYP(I,J)
#endif /* TRACERS_AEROSOLS_SOA */
                 CASE ("BCIA")
                    TRACER(L,n) = (TRM(i,j,l,n_BCII)+TRM(i,j,l,     &
                         n_BCIA))*BYAXYP(I,J)
                 CASE DEFAULT
#ifdef TRACERS_NITRATE
                    ! assume full neutralization of NO3p, if NH4 suffice
                    SELECT CASE (TRNAME(NTRIX_AOD(n)))
                    CASE ("NO3p")
                       IF ( TRM(i,j,l,NTRIX_AOD(n))>0.D0 ) THEN
                          nh4_on_no3 = MIN(TRM(i,j,l,n_NO3p)        &
                               *(TR_MM(n_NO3p)+TR_MM(n_NH4))          &
                               /TR_MM(n_NO3p)-TRM(i,j,l,n_NO3p),      &
                               TRM(i,j,l,n_NH4))
                          WTTR(n) = (nh4_on_no3+TRM(i,j,l,NTRIX_AOD(&
                               n)))/TRM(i,j,l,NTRIX_AOD(n))
                       ENDIF
                    CASE ("SO4")
                       IF ( TRM(i,j,l,NTRIX_AOD(n))>0.D0 ) THEN
                          nh4_on_no3 = MIN(TRM(i,j,l,n_NO3p)        &
                               *(TR_MM(n_NO3p)+TR_MM(n_NH4))          &
                               /TR_MM(n_NO3p)-TRM(i,j,l,n_NO3p),      &
                               TRM(i,j,l,n_NH4))
                          WTTR(n) = (TRM(i,j,l,n_NH4)-nh4_on_no3+TRM&
                               (i,j,l,NTRIX_AOD(n)))           &
                               /TRM(i,j,l,NTRIX_AOD(n))
                       ENDIF
                    ENDSELECT
#endif
                    TRACER(L,n) = WTTR(n)*TRM(i,j,l,NTRIX_AOD(n))   &
                         *BYAXYP(I,J)
                 ENDSELECT
              ENDDO
#endif /* TRACERS_AEROSOLS_Koch/DUST/MINERALS/SEASALT */

#ifdef TRACERS_AMP
              CALL SETAMP_LEV(i,j,l)
#endif
#ifdef TRACERS_TOMAS
              CALL SETTOMAS_LEV(i,j,l)
#endif
           ENDDO
           !**** Radiative Equilibrium Layer data
           DO K = 1, LM_REQ
              !CC     STOP 'In Radia :  RQT out of range'
              !       JCKERR=JCKERR+1
              IF ( INT(RQT(K,I,J))<planck_tmin .OR. INT(RQT(K,I,J)) &
                   >=planck_tmax ) WRITE (6,*)                      &
                   'In RADIA :  Time,I,J,L,TL', ITime, I, J, LM + K, &
                   RQT(K,I,J)
              TLM(LM+K) = RQT(K,I,J)
              PLB(LM+k+1) = PLB0(k)
              SHL(LM+k) = SHL0(k)
              RHL(LM+k) = SHL(LM+k)                                 &
                   /QSAT(TLM(LM+k),LHE,.5D0*(PLB(LM+k)       &
                   +PLB(LM+k+1)))
              TAUWC(LM+k) = 0.
              TAUIC(LM+k) = 0.
              SIZEWC(LM+k) = 0.
              SIZEIC(LM+k) = 0.
#ifdef TRACERS_ON
              !**** set radiative equilibrium extra tracer amount to zero
              IF ( nraero_aod>0 ) TRACER(LM+k,1 : nraero_aod) = 0.
#endif
           ENDDO
           IF ( kradia>1 ) THEN
              DO l = 1, lm + lm_req
                 TLM(l) = TLM(l) + TCHG(l,i,j)
                 AFLX_ST(L,I,J,5) = AFLX_ST(L,I,J,5) + TCHG(L,I,J)
              ENDDO
           ENDIF
           !**** Zenith angle and GROUND/SURFACE parameters
           COSZ = COSZA(I,J)
           TGO = atmocn%GTEMPR(I,J)
           TGOI = atmice%GTEMPR(I,J)
           TGLI = atmgla%GTEMPR(I,J)
           TGE = atmlnd%GTEMPR(I,J)
           TSL = atmsrf%TSAVG(I,J)
           SNOWOI = SNOWI(I,J)
           SNOWLI = atmgla%SNOW(I,J)
           !SNOWE=atmlnd%SNOWE(I,J)                    ! snow depth (kg/m**2)
           SNOWD( : ) = snowd_ij( : ,I,J)
           snow_frac( : ) = atmlnd%FR_SNOW_RAD( : ,i,j)
           ! snow cover (1)
           AGESN(1) = SNOAGE(3,I,J)
           ! land         ! ? why are these numbers
           AGESN(2) = SNOAGE(1,I,J)
           ! ocean ice        so confusing ?
           AGESN(3) = SNOAGE(2,I,J)
           ! land ice
           !      print*,"snowage",i,j,SNOAGE(1,I,J)
           !**** set up parameters for new sea ice and snow albedo
           zsnwoi = atmice%ZSNOWI(I,J)
           ! LTM :  Fix, testing equality of reals is not reliable
           IF ( dalbsnX/=0. ) THEN
#ifdef OLD_BCdalbsn
              dALBsn = xdalbs*DEPOBC(i,j)
#else
              dALBsn = dalbsnX*BCDALBSN(i,j)
#endif
           ELSE
              dALBsn = 0.
           ENDIF

           ! to use on-line tracer albedo impact, set dALBsnX=0. in rundeck
#ifdef BC_ALB
           CALL GET_BC_DALBEDO(i,j,dALBsn1,bc_snow_present(i,j))
           IF ( rad_interact_aer>0 ) dALBsn = dALBsn1
           dALBsnBC(I,J) = dALBsn1
#endif  /* BC_ALB */
           IF ( poice>0. ) THEN
              zoice = ZSI(i,j)
              flags = flag_dsws(i,j)
              IF ( kradia<=0 ) THEN
                 fmp = MIN(1.6D0*SQRT(pond_melt(i,j)/rhow),1D0)
                 AIJ(I,J,IJ_FRMP) = AIJ(I,J,IJ_FRMP) + fmp*POICE
              ELSE
                 fmp = fmp_com(i,j)
              ENDIF
              zmp = MIN(0.8D0*fmp,0.9D0*zoice)
           ELSE
              zoice = 0.
              flags = .FALSE.
              fmp = 0.
              zmp = 0.
           ENDIF
           !**** set up new lake depth parameter to incr. albedo for shallow lakes
           !      zlake=0.
           !      if (plake.gt.0) then
           !        zlake = MWL(I,J)/(RHOW*PLAKE*AXYP(I,J))
           !      end if
           zlake = dlake(i,j)
           !****
           IF ( kradia<=0 ) THEN
              !WEARTH=(WEARTH_COM(I,J)+AIEARTH(I,J))/(WFCS(I,J)+1.D-20)
              WEARTH = atmlnd%BARE_SOIL_WETNESS(i,j)
              IF ( wearth>1. ) wearth = 1.
           ELSE                   ! rad.frc. model
              wearth = wsoil(i,j)
           ENDIF
           IF ( FEARTH(i,j)>0.D0 ) THEN
              CALL ENT_GET_EXPORTS(ENTCELLS(i,j),                   &
                   VEGETATION_FRACTIONS=PVT0,       &
                   VEGETATION_HEIGHTS=HVT0)
              CALL MAP_ENT2GISS(PVT0,HVT0,PVT)
              !temp hack :  ent pfts->giss veg
           ELSE
              PVT( : ) = 0.D0
              ! actually PVT is not supposed to be used in this case
           ENDIF
           WMAG = atmsrf%WSAVG(I,J)
           !****
           !**** Radiative interaction and forcing diagnostics : 
           !**** If no radiatively active tracers are defined, nothing changes.
           !**** Currently this works for aerosols and ozone but should be extended
           !**** to cope with all trace gases.
           !****
           FTAUC = 1.
           ! deflt (clouds on)
           use_tracer_chem(:) = 0
           ! by default use climatological ozone/ch4/co2
           !**** Set level for inst. rad. forc. calcs for aerosols/trace gases
           !**** This is set from the rundeck.
           LFRC = LM + LM_REQ + 1
           ! TOA
           IF ( rad_forc_lev>0 ) LFRC = LTROPO(I,J)
           ! TROPOPAUSE
#ifdef ACCMIP_LIKE_DIAGS
           IF ( rad_forc_lev>0 ) CALL STOP_MODEL(                   &
                &'ACCMIP_LIKE_DIAGS desires TOA RADF diags',255)
#endif
           !**** The calculation of the forcing is slightly different.
           !**** depending on whether full radiative interaction is turned on
           !**** or not.
           onoff_aer = 0
           onoff_chem = 0
           IF ( rad_interact_aer   > 0 ) onoff_aer  = 1
           IF ( clim_interact_chem > 0 ) onoff_chem = 1
           use_o3_ref = 0

#ifdef TRACERS_GC

           IF ( SUM( TrM(I,J,1:LM,i_CH4) ) .lt. 1d-20 ) THEN
             ! If methane is not initialized yet assume 1.5 ppmv everywhere
             CHEM_IN(2,1:LM) = MA(1:LM,I,J) * 1.5e-6 * 16.04 / 28.97
											ELSE
												 !IF ( i_CH4 > 0 )
													IF ( Am_I_Root() .and. I .eq. 1 .and. J .eq. 1 ) &
													     WRITE(6,*) 'Updating methane in radiation code...'  
													CHEM_IN(2,1:LM) = TrM(I,J,1:LM,i_CH4) * BYAXYP(I,J)
											ENDIF
												
											IF ( SUM( TrM(I,J,1:LM,i_O3) ) .lt. 1d-20 ) THEN	
             ! If ozone is not initialized yet assume 5 ppmv everywhere
             CHEM_IN(1,1:LM) = MA(1:LM,I,J) * 5.0e-6 * 47.997 / 28.97
											ELSE
             !IF (  i_O3 > 0 ) 
													IF ( Am_I_Root() .and. I .eq. 1 .and. J .eq. 1 ) &
                  WRITE(6,*) 'Updating ozone in radiation code...'
													CHEM_IN(1,1:LM) = TrM(I,J,1:LM, i_O3) * BYAXYP(I,J)	  						
											ENDIF
											
           IF ( clim_interact_chem > 0 ) THEN
              use_tracer_chem(1) = LM ! Lmax_rad_O3  ! O3
              use_tracer_chem(2) = LM ! Lmax_rad_CH4 ! CH4
           ENDIF
#endif

#ifdef TRACERS_SPECIAL_Shindell
           !**** Ozone and Methane : 
           CHEM_IN(1,1:LM) = chem_tracer_save(1,1:LM, I,J)
           CHEM_IN(2,1:LM) = chem_tracer_save(2,1:LM, I,J)           &
                *CH4X_RADoverCHEM
           IF ( clim_interact_chem>0 ) THEN
              use_tracer_chem(1) = Lmax_rad_O3
              ! O3
              use_tracer_chem(2) = Lmax_rad_CH4
              ! CH4
           ENDIF
#if (defined SHINDELL_STRAT_EXTRA) &(defined ACCMIP_LIKE_DIAGS)
           IF ( clim_interact_chem<=0 ) CALL STOP_MODEL(            &
                &"stratOx RADF on, clim_interact_chem<=0",255)
#endif /* SHINDELL_STRAT_EXTRA &ACCMIP_LIKE_DIAGS */
#endif /* TRACERS_SPECIAL_Shindell */

           !**** CO2
           ! Update GCCco2_IN and GCCco2_tracer_save with trm to pass on CO2
           ! information to RADIATION
#ifdef GCC_COUPLE_RAD
           DO L = 1, LM
              GCCCO2_TRACER_SAVE(L,i,j) = (TRM(i,j,L,n_CO2n))       &
                   *BYAXYP(i,j)*avog/(TR_MM(n_CO2n)*2.69E20)
           ENDDO
           GCCco2_IN(1:LM) = GCCCO2_TRACER_SAVE(1:LM, I,J)*CO2X
           use_tracer_GCCco2 = Lmax_rad_CO2
           ! CO2
#endif /* GCC_COUPLE_RAD */

           IF ( moddrf==0 ) THEN
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
              (defined TRACERS_MINERALS) || (defined TRACERS_AEROSOLS_SEASALT) ||\
              (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
              !**** Aerosols (OMA, MATRIX, TOMAS) : 
              DO n = 1, nraero_rf
                 IF ( TRNAME(NTRIX_RF(n))=="seasalt2" ) CYCLE
                 ! not for seasalt2
                 IF ( diag_fc==2 ) THEN
                    FSTOPX(n) = 1 - onoff_aer
                    !turns off online tracer
                    FTTOPX(n) = 1 - onoff_aer
                    !
                    !**** Warning :  small bit of hardcoding assumes that seasalt2 immediately
                    !****          succeeds seasalt1 in nraero_rf array
                    IF ( TRNAME(NTRIX_RF(n))=="seasalt1" ) THEN
                       !add seasalt2
                       FSTOPX(n+1) = 1 - onoff_aer
                       FTTOPX(n+1) = 1 - onoff_aer        !to seasalt1
                    ENDIF
                 ELSEIF ( diag_fc==1 ) THEN
                    FSTOPX(1 : nraero_aod) = 1 - onoff_aer
                    !turns off online tracer
                    FTTOPX(1 : nraero_aod) = 1 - onoff_aer
                    !
                 ENDIF
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
                 (defined TRACERS_MINERALS) || (defined TRACERS_AEROSOLS_SEASALT)
                 kdeliq(1:LM, 1:4) = kliq(1:LM, 1:4,i,j)
#endif
                 CALL RCOMPX
                 ! tr.aero.Koch/dust/miner./seasalt
                 SNFST(1,n,I,J) = SRNFLB(1)
                 ! surface forcing
                 TNFST(1,n,I,J) = TRNFLB(1)
                 SNFST(2,n,I,J) = SRNFLB(LFRC)
                 ! Tropopause forcing
                 TNFST(2,n,I,J) = TRNFLB(LFRC)
                 IF ( diag_fc==2 ) THEN
                    FSTOPX(n) = onoff_aer
                    !turns on online tracer
                    FTTOPX(n) = onoff_aer
                    !
                    IF ( TRNAME(NTRIX_RF(n))=="seasalt1" ) THEN
                       ! also for seasalt2
                       FSTOPX(n+1) = onoff_aer
                       FTTOPX(n+1) = onoff_aer
                    ENDIF
                 ELSEIF ( diag_fc==1 ) THEN
                    FSTOPX(1 : nraero_aod) = onoff_aer
                    !turns on online tracer
                    FTTOPX(1 : nraero_aod) = onoff_aer
                    !
                 ENDIF
              ENDDO
#endif

#ifdef TRACERS_GC
              ! Get flux values minus ozone at various heights
														! Use constant reference year for first call as with Shindell tracers
              !use_o3_ref = 1
														! Do not use constant reference year
              use_o3_ref = 0
              use_tracer_chem(1) = 0
              kdeliq( 1:LM, 1:4 ) = kliq( 1:LM, 1:4, i, j ) 
              CALL RCOMPX
              ! Meteorological tropopause
              SNFST_o3ref(1,I,J)    = SRNFLB(LTROPO(I,J))
              TNFST_o3ref(1,I,J)    = TRNFLB(LTROPO(I,J))
              ! Top of the atmosphere
              SNFST_o3ref(2,I,J)    = SRNFLB(LM+LM_REQ+1)
              TNFST_o3ref(2,I,J)    = TRNFLB(LM+LM_REQ+1)
              ! Whole atmosphere
														SNFS_3D_pert(I,J,:,5) = SRNFLB
														TNFS_3D_pert(I,J,:,5) = TRNFLB

              use_o3_ref = 0
              use_tracer_chem(1) = onoff_chem * LM !Lmax_rad_O3

              IF ( SUM( TrM(I,J,1:LM,i_CH4) ) .lt. 1d-20 ) THEN
                ! If methane is not initialized yet assume 1.5 ppmv everywhere
                CHEM_IN(2,1:LM) = MA(1:LM,I,J) * 1.5e-6 * 16.04 / 28.97
              ELSE
              	 !IF ( i_CH4 > 0 )
              		CHEM_IN(2,1:LM) = TrM(I,J,1:LM,i_CH4) * BYAXYP(I,J)
              ENDIF
              	
              IF ( SUM( TrM(I,J,1:LM,i_O3) ) .lt. 1d-20 ) THEN	
                ! If ozone is not initialized yet assume 5 ppmv everywhere
                CHEM_IN(1,1:LM) = MA(1:LM,I,J) * 5.0e-6 * 47.997 / 28.97
              ELSE
                !IF (  i_O3 > 0 ) 
              		CHEM_IN(1,1:LM) = TrM(I,J,1:LM, i_O3) * BYAXYP(I,J)	  						
              ENDIF
																			
#endif

#ifdef TRACERS_SPECIAL_Shindell
              !**** Ozone : 
              ! ozone rad forcing diags now use a constant reference year
              ! for this first call. And no tracer values...
              use_o3_ref = 1
              use_tracer_chem(1) = 0
              kdeliq(1:LM, 1:4) = kliq(1:LM, 1:4,i,j)
              CALL RCOMPX
              ! tr_Shindell Ox tracer
              SNFST_o3ref(1,I,J) = SRNFLB(LTROPO(I,J))
              ! meteorological tropopause
              TNFST_o3ref(1,I,J) = TRNFLB(LTROPO(I,J))
              SNFST_o3ref(2,I,J) = SRNFLB(LM+LM_REQ+1)
              ! T.O.A.
              TNFST_o3ref(2,I,J) = TRNFLB(LM+LM_REQ+1)
              SNFST_o3ref(5,I,J) = SRNFLB(LS1-1)
              ! fixed tropopause
              TNFST_o3ref(5,I,J) = TRNFLB(LS1-1)

#ifdef AUXILIARY_OX_RADF
              ! if needed, also save the auxiliary ozone field (i.e. climatology
              ! if tracer is used in final call, tracers if climatology is used.)
#ifdef AUX_OX_RADF_TROP
              ! forces use of tracer from L=1,LS1-1 and reference above that :
              use_o3_ref = 1
              use_tracer_chem(1) = LS1 - 1
#else
              ! use tracer or climatology, whichever won''t be used in final call :
              use_o3_ref = 0
              use_tracer_chem(1) = (1-onoff_chem)*Lmax_rad_O3
#endif
              kdeliq(1:LM, 1:4) = kliq(1:LM, 1:4,i,j)
              CALL RCOMPX
              ! tr_Shindell Ox tracer
#ifdef AUX_OX_RADF_TROP
              SNFST_o3ref(3,I,J) = SRNFLB(LS1-1)
              ! fixed tropopause
              TNFST_o3ref(3,I,J) = TRNFLB(LS1-1)
#else
              SNFST_o3ref(3,I,J) = SRNFLB(LTROPO(I,J))
              ! meteorological tropopause
              TNFST_o3ref(3,I,J) = TRNFLB(LTROPO(I,J))
#endif
              SNFST_o3ref(4,I,J) = SRNFLB(LM+LM_REQ+1)
              ! T.O.A.
              TNFST_o3ref(4,I,J) = TRNFLB(LM+LM_REQ+1)
#endif /* AUXILIARY_OX_RADF */
              ! After AUX call, use either climatological or tracer O3 : 
              use_o3_ref = 0
              use_tracer_chem(1) = onoff_chem*Lmax_rad_O3
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
              ! Optional intermediate call with stratOx tracer : 
              !NEED CHEM_IN(1,1:LM)=stratO3_tracer_save(1:LM, I,J)
              !NEED kdeliq(1:LM, 1:4)=kliq(1:LM, 1:4,i,j)
              !NEED CALL RCOMPX        ! stratOx diag tracer
              ! Tropopause
              SNFST_stratOx(1,I,J) = SRNFLB(LTROPO(I,J))
              TNFST_stratOx(1,I,J) = TRNFLB(LTROPO(I,J))
              ! T.O.A.
              SNFST_stratOx(2,I,J) = SRNFLB(LM+LM_REQ+1)
              TNFST_stratOx(2,I,J) = TRNFLB(LM+LM_REQ+1)
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS */
              CHEM_IN(1,1:LM) = chem_tracer_save(1,1:LM, I,J)                    ! Ozone
              CHEM_IN(2,1:LM) = chem_tracer_save(2,1:LM, I,J) * CH4X_RADoverCHEM ! Methane
#if (defined ACCMIP_LIKE_DIAGS)
#ifndef SKIP_ACCMIP_GHG_RADF_DIAGS
              ! TOA GHG rad forcing :  nf=1,4 are CH4, N2O, CFC11, CFC12 : 
              ! Initial calls are reference year/day : 
              DO nf = 1, 4
                 IF ( nf==1 ) THEN
                    ! CH4 reference call must not use tracer
                    use_tracer_chem(2) = 0
                 ELSE
                    ! N2O and CFC call's CH4 should match final call
                    use_tracer_chem(2) = onoff_chem*Lmax_rad_CH4
                 ENDIF
                 FULGAS(nfghg(nf)) = sv_fulgas_ref(nf)
                 kdeliq(1:LM, 1:4) = kliq(1:LM, 1:4,i,j)
                 CALL RCOMPX
                 SNFS_ghg(nf,I,J) = SRNFLB(LM+LM_REQ+1)
                 TNFS_ghg(nf,I,J) = TRNFLB(LM+LM_REQ+1)
                 FULGAS(nfghg(nf)) = sv_fulgas_now(nf)
              ENDDO
#endif /* NOT DEFINED SKIP_ACCMIP_GHG_RADF_DIAGS */
#endif /* ACCMIP_LIKE_DIAGS */
#endif /* TRACERS_SPECIAL_Shindell */

#if defined ( TRACERS_GC )

              IF ( SUM( TrM(I,J,1:LM,i_CH4) ) .lt. 1d-20 ) THEN
                ! If methane is not initialized yet assume 1.5 ppmv everywhere
                CHEM_IN(2,1:LM) = MA(1:LM,I,J) * 1.5e-6 * 16.04 / 28.97
              ELSE
              	 !IF ( i_CH4 > 0 )
              		CHEM_IN(2,1:LM) = TrM(I,J,1:LM,i_CH4) * BYAXYP(I,J)
              ENDIF
              	
              IF ( SUM( TrM(I,J,1:LM,i_O3) ) .lt. 1d-20 ) THEN	
                ! If ozone is not initialized yet assume 5 ppmv everywhere
                CHEM_IN(1,1:LM) = MA(1:LM,I,J) * 5.0e-6 * 47.997 / 28.97
              ELSE
                !IF (  i_O3 > 0 ) 
              		CHEM_IN(1,1:LM) = TrM(I,J,1:LM, i_O3) * BYAXYP(I,J)	  						
              ENDIF

              ! TOA GHG rad forcing :  nf=1,4 are CH4, N2O, CFC11, CFC12 : 
              ! Initial calls are reference year/day : 
              DO nf = 1, 4
                 IF ( nf==1 ) THEN
                    ! CH4 reference call must not use tracer
                    use_tracer_chem(2) = 0
                 ELSE
                    ! N2O and CFC call's CH4 should match final call
                    use_tracer_chem(2) = onoff_chem*LM ! Lmax_rad_CH4
                 ENDIF
                 FULGAS(nfghg(nf)) = sv_fulgas_ref(nf)
                 kdeliq(1:LM, 1:4) = kliq(1:LM, 1:4,i,j)
                 CALL RCOMPX
																	! TOA
                 SNFS_ghg(nf,I,J)       = SRNFLB(LM+LM_REQ+1)
                 TNFS_ghg(nf,I,J)       = TRNFLB(LM+LM_REQ+1)
                 ! Tropopause
                 SNFS_ghg_tp(nf,I,J)    = SRNFLB(LTROPO(I,J))
                 TNFS_ghg_tp(nf,I,J)    = TRNFLB(LTROPO(I,J))
                 ! Whole atmosphere
																	SNFS_3D_pert(I,J,:,nf) = SRNFLB
																	TNFS_3D_pert(I,J,:,nf) = TRNFLB
                 FULGAS(nfghg(nf))      = sv_fulgas_now(nf)
              ENDDO

#endif

           ENDIF
           ! moddrf=0
#if (defined GCC_COUPLE_RAD)
           ! final (main) RCOMPX call can use tracer co2 (or not) : 
           use_tracer_GCCco2 = Lmax_rad_CO2
           IF ( IS_SET_PARAM('initial_GHG_setup') ) THEN
              CALL GET_PARAM('initial_GHG_setup',initial_GHG_setup)
              IF ( initial_GHG_setup==1 .AND. itime==itimeI )       &
                   use_tracer_GCCco2 = 0
              ! special case; model outputs climatology
           ENDIF
#endif /* GCC_COUPLE_RAD */
#if (defined TRACERS_SPECIAL_Shindell)
           ! final (main) RCOMPX call can use tracer methane (or not) : 
           use_tracer_chem(2) = onoff_chem*Lmax_rad_CH4
           IF ( IS_SET_PARAM('initial_GHG_setup') ) THEN
              CALL GET_PARAM('initial_GHG_setup',initial_GHG_setup)
              IF ( initial_GHG_setup==1 .AND. itime==itimeI )       &
                   use_tracer_chem(2) = 0
              ! special case; model outputs climatology
           ENDIF
#endif /* TRACERS_SPECIAL_Shindell */

#ifdef GCC_UNCOUPLE_RAD_CONCEN
           ! Use reference year CO2 for uncoupling radiation
           FULGAS(2) = GCCco2_fulgas_ref
#endif

           IF ( moddrf==0 ) THEN
#ifdef BC_ALB
              IF ( rad_interact_aer>0 ) dalbsn = 0.D0
              CALL RCOMPX
              NFSNBC(I,J) = SRNFLB(LM+LM_REQ+1)
              !       NFSNBC(I,J)=SRNFLB(LFRC)
              ALBNBC(I,J) = SRNFLB(1)/(SRDFLB(1)+1.D-20)
              ! set for BC-albedo effect
              IF ( rad_interact_aer>0 ) dALBsn = dALBsn1
#endif
              !**** Optional calculation of CRF using a clear sky calc.
              IF ( cloud_rad_forc>0 ) THEN
                 FTAUC = 0.
                 ! turn off cloud tau (tauic +tauwc)
                 kdeliq(1:LM, 1:4) = kliq(1:LM, 1:4,i,j)
                 CALL RCOMPX
                 ! cloud_rad_forc>0  :  clr sky
                 SNFSCRF(I,J) = SRNFLB(LM+LM_REQ+1)
                 ! always TOA
                 TNFSCRF(I,J) = TRNFLB(LM+LM_REQ+1)
                 ! always TOA
                 LWDNCS(I,J) = STBO*(POCEAN*atmocn%GTEMPR(I,J)**4+  &
                      POICE*atmice%GTEMPR(I,J)             &
                      **4+PLICE*atmgla%GTEMPR(I,J)         &
                      **4+PEARTH*atmlnd%GTEMPR(I,J)**4)    &
                      - TRNFLB(1)
                 ! clr sky trhr(0)
                 !         BEGIN AMIP
                 AIJ(I,J,IJ_SWDCLS) = AIJ(I,J,IJ_SWDCLS) + SRDFLB(1)&
                      *COSZ2(I,J)
                 AIJ(I,J,IJ_SWNCLS) = AIJ(I,J,IJ_SWNCLS) + SRNFLB(1)&
                      *COSZ2(I,J)
                 AIJ(I,J,IJ_LWDCLS) = AIJ(I,J,IJ_LWDCLS) + TRDFLB(1)
                 AIJ(I,J,IJ_SWNCLT) = AIJ(I,J,IJ_SWNCLT)            &
                      + SRNFLB(LM+LM_REQ+1)*COSZ2(I,J)
                 AIJ(I,J,IJ_LWNCLT) = AIJ(I,J,IJ_LWNCLT)            &
                      + TRNFLB(LM+LM_REQ+1)
                 !       END AMIP
#ifdef CFMIP3_SUBDD
                 ! SW upward flux at TOA, Csky
                 !swutcs(i,j)=sruflb(lm)*csz2
                 swutcs(i,j) = SRUFLB(lm)*cosz2(i,j)
                 ! SW downward flux at SFC, Csky
                 swdcls(i,j) = SRDFLB(1)*cosz2(i,j)
                 ! SW upward flux at SFC, Csky
                 swucls(i,j) = SRUFLB(1)*cosz2(i,j)
#endif
              ENDIF
              FTAUC = 1.
              ! default :  turn on cloud tau


              !**** 2nd Optional calculation of CRF using a clear sky calc. without aerosols and Ox
              IF ( cloud_rad_forc==2 ) THEN
                 FTAUC = 0.
                 ! turn off cloud tau (tauic +tauwc)
                 kdeliq(1:LM, 1:4) = kliq(1:LM, 1:4,i,j)
                 ! Including turn off of aerosols and Ox during crf calc.+++++++++++++++++++
#ifdef TRACERS_SPECIAL_Shindell
                 use_o3_ref = 1
                 use_tracer_chem(1) = 0 !turns off ozone
#endif
                 FSTOPX( : ) = 0
                 !turns off aerosol tracers
                 FTTOPX( : ) = 0
                 CALL RCOMPX
                 ! cloud_rad_forc=2  :  clr sky
                 FSTOPX( : ) = onoff_aer
                 !turns on aerosol tracers, if requested
                 FTTOPX( : ) = onoff_aer
                 !
#ifdef TRACERS_SPECIAL_Shindell
                 use_o3_ref = 0
                 use_tracer_chem(1) = onoff_chem*Lmax_rad_O3
                 ! turns on ozone tracers
#endif
                 SNFSCRF2(I,J) = SRNFLB(LM+LM_REQ+1)
                 ! always TOA
                 TNFSCRF2(I,J) = TRNFLB(LM+LM_REQ+1)
                 ! always TOA
              ENDIF
              FTAUC = 1.
              ! default :  turn on cloud tau

              IF ( cloud_rad_forc>0 ) THEN
                 !**** all sky calc. without aerosol
                 kdeliq(1:LM, 1:4) = kliq(1:LM, 1:4,i,j)
                 FSTOPX( : ) = 0
                 !turns off aerosol tracers
                 FTTOPX( : ) = 0
                 CALL RCOMPX
                 !  all sky
                 FSTOPX( : ) = onoff_aer
                 !turns on aerosol tracers, if requested
                 FTTOPX( : ) = onoff_aer
                 !

                 SNFS_AS_noA(I,J) = SRNFLB(LM+LM_REQ+1)
                 ! always TOA
                 TNFS_AS_noA(I,J) = TRNFLB(LM+LM_REQ+1)
                 ! always TOA

                 !**** clear sky calc. without aerosol
                 FTAUC = 0.
                 ! turn off cloud tau (tauic +tauwc)
                 kdeliq(1:LM, 1:4) = kliq(1:LM, 1:4,i,j)
                 FSTOPX( : ) = 0
                 !turns off aerosol tracers
                 FTTOPX( : ) = 0
                 CALL RCOMPX
                 !  clr sky
                 FSTOPX( : ) = onoff_aer
                 !turns on aerosol tracers, if requested
                 FTTOPX( : ) = onoff_aer
                 !

                 SNFS_CS_noA(I,J) = SRNFLB(LM+LM_REQ+1)
                 ! always TOA
                 TNFS_CS_noA(I,J) = TRNFLB(LM+LM_REQ+1)
                 ! always TOA
                 FTAUC = 1.
                 ! default :  turn on cloud tau
              ENDIF

              !**** Optional calculation of the impact of NINT aerosols
              IF ( aer_rad_forc>0 ) THEN
                 !**** first, separate aerosols
                 DO N = 1, 8
                    tmpS(N) = FS8OPX(N)
                    tmpT(N) = FT8OPX(N)
                    FS8OPX(N) = 0.
                    FT8OPX(N) = 0.
                    kdeliq(1:LM, 1:4) = kliq(1:LM, 1:4,i,j)
                    CALL RCOMPX
                    ! aer_rad_forc>0  :  no aerosol N
                    SNFSAERRF(N,I,J) = SRNFLB(LM+LM_REQ+1)
                    ! TOA
                    TNFSAERRF(N,I,J) = TRNFLB(LM+LM_REQ+1)
                    ! TOA
                    SNFSAERRF(N+8,I,J) = SRNFLB(1)
                    ! SURF
                    TNFSAERRF(N+8,I,J) = TRNFLB(1)
                    ! SURF
                    FS8OPX(N) = tmpS(N)
                    FT8OPX(N) = tmpT(N)
                 ENDDO
                 !**** second, net aerosols
                 tmpS( : ) = FS8OPX( : )
                 tmpT( : ) = FT8OPX( : )
                 FS8OPX( : ) = 0.
                 FT8OPX( : ) = 0.
                 kdeliq(1:LM, 1:4) = kliq(1:LM, 1:4,i,j)
                 CALL RCOMPX  ! aer_rad_forc>0  :  no aerosols
                 SNFSAERRF(17,I,J) = SRNFLB(LM+LM_REQ+1)
                 ! TOA
                 TNFSAERRF(17,I,J) = TRNFLB(LM+LM_REQ+1)
                 ! TOA
                 SNFSAERRF(18,I,J) = SRNFLB(1)
                 ! SURF
                 TNFSAERRF(18,I,J) = TRNFLB(1)
                 ! SURF
                 FS8OPX( : ) = tmpS( : )
                 FT8OPX( : ) = tmpT( : )
              ENDIF
           ENDIF
           ! moddrf=0

           !**** End of initial computations for optional forcing diagnostics

           !**** Localize fields that are modified by RCOMPX
           kdeliq(1:LM, 1:4) = kliq(1:LM, 1:4,i,j)

           !*****************************************************
           !     Main RADIATIVE computations, SOLAR and THERM(A)L
           CALL RCOMPX
           !*****************************************************

#ifdef CACHED_SUBDD
           CO2out(1:LM, i,j) = CO2outCol(1:LM)
#endif
#ifdef GCC_UNCOUPLE_RAD_CONCEN
           ! Put back the actual amount of CO2 to fulgas
           FULGAS(2) = GCCco2_fulgas_now
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
           (defined TRACERS_MINERALS) || (defined TRACERS_AMP) ||\
           (defined TRACERS_TOMAS) || (defined TRACERS_AEROSOLS_SEASALT)

           !**** Save optical depth diags
           nsub_ntrix = 0
           DO n = 1, nraero_aod
              SELECT CASE (TRNAME(NTRIX_AOD(n)))
              CASE ('Clay','ClayIlli','ClayKaol','ClaySmec',        &
                   &'ClayCalc','ClayQuar','ClayFeld','ClayHema',    &
                   &'ClayGyps','ClayIlHe','ClayKaHe','ClaySmHe',    &
                   &'ClayCaHe','ClayQuHe','ClayFeHe','ClayGyHe')
                 nsub_ntrix(NTRIX_AOD(n)) = nsub_ntrix(NTRIX_AOD(n))&
                      + 1

                 ! 3d aod
                 IF ( diag_aod_3d>0 .AND. diag_aod_3d<5 ) THEN
                    ! valid values are 1-4
                    IF ( IJLT_3DAAOD(n)>0 )                         &
                         taijls(i,j,1:LM, IJLT_3DAAOD(n))            &
                         = taijls(i,j,1:LM, IJLT_3DAAOD(n))          &
                         + (aesqex(1:LM, 6,n)-aesqsc(1:LM, 6,n))
                    IF ( IJLT_3DAAODCS(n)>0 )                       &
                         taijls(i,j,1:LM, IJLT_3DAAODCS(n))          &
                         = taijls(i,j,1:LM, IJLT_3DAAODCS(n))        &
                         + (aesqex(1:LM, 6,n)-aesqsc(1:LM, 6,n))      &
                         *OPNSKY
                    IF ( IJLT_3DAAODDRY(n)>0 )                      &
                         taijls(i,j,1:LM, IJLT_3DAAODDRY(n))         &
                         = taijls(i,j,1:LM, IJLT_3DAAODDRY(n))       &
                         + (aesqex_dry(1:LM, 6,n)                    &
                         -aesqsc_dry(1:LM, 6,n))
                    IF ( IJLT_3DTAU(n)>0 )                          &
                         taijls(i,j,1:LM, IJLT_3DTAU(n))             &
                         = taijls(i,j,1:LM, IJLT_3DTAU(n))           &
                         + aesqex(1:LM, 6,n)
                    IF ( IJLT_3DTAUCS(n)>0 )                        &
                         taijls(i,j,1:LM, IJLT_3DTAUCS(n))           &
                         = taijls(i,j,1:LM, IJLT_3DTAUCS(n))         &
                         + aesqex(1:LM, 6,n)*OPNSKY
                    IF ( IJLT_3DTAUDRY(n)>0 )                       &
                         taijls(i,j,1:LM, IJLT_3DTAUDRY(n))          &
                         = taijls(i,j,1:LM, IJLT_3DTAUDRY(n))        &
                         + aesqex_dry(1:LM, 6,n)
                 ELSEIF ( diag_aod_3d<0 .AND. diag_aod_3d>-5 ) THEN
                    ! if negative, save total
                    IF ( IJLT_3DAAOD(1)>0 )                         &
                         taijls(i,j,1:LM, IJLT_3DAAOD(1))            &
                         = taijls(i,j,1:LM, IJLT_3DAAOD(1))          &
                         + (aesqex(1:LM, 6,n)-aesqsc(1:LM, 6,n))
                    IF ( IJLT_3DAAODCS(1)>0 )                       &
                         taijls(i,j,1:LM, IJLT_3DAAODCS(1))          &
                         = taijls(i,j,1:LM, IJLT_3DAAODCS(1))        &
                         + (aesqex(1:LM, 6,n)-aesqsc(1:LM, 6,n))      &
                         *OPNSKY
                    IF ( IJLT_3DAAODDRY(1)>0 )                      &
                         taijls(i,j,1:LM, IJLT_3DAAODDRY(1))         &
                         = taijls(i,j,1:LM, IJLT_3DAAODDRY(1))       &
                         + (aesqex_dry(1:LM, 6,n)                    &
                         -aesqsc_dry(1:LM, 6,n))
                    IF ( IJLT_3DTAU(1)>0 )                          &
                         taijls(i,j,1:LM, IJLT_3DTAU(1))             &
                         = taijls(i,j,1:LM, IJLT_3DTAU(1))           &
                         + aesqex(1:LM, 6,n)
                    IF ( IJLT_3DTAUCS(1)>0 )                        &
                         taijls(i,j,1:LM, IJLT_3DTAUCS(1))           &
                         = taijls(i,j,1:LM, IJLT_3DTAUCS(1))         &
                         + aesqex(1:LM, 6,n)*OPNSKY
                    IF ( IJLT_3DTAUDRY(1)>0 )                       &
                         taijls(i,j,1:LM, IJLT_3DTAUDRY(1))          &
                         = taijls(i,j,1:LM, IJLT_3DTAUDRY(1))        &
                         + aesqex_dry(1:LM, 6,n)
                 ENDIF
                 ! 0<diag_aod_3d<5 or 0>diag_aod_3d>-5

                 ! 2d aod, per band or just band6, depending on diag_rad
                 IF ( diag_rad/=1 ) THEN
                    IF ( IJTS_TAUSUB(1,NTRIX_AOD(n),nsub_ntrix(     &
                         NTRIX_AOD(n)))>0 )                         &
                         TAIJS(i,j,IJTS_TAUSUB(1,NTRIX_AOD(n),      &
                         nsub_ntrix(NTRIX_AOD(n))))                 &
                         = TAIJS(i,j,IJTS_TAUSUB(1,NTRIX_AOD(n),    &
                         nsub_ntrix(NTRIX_AOD(n))))                 &
                         + SUM(aesqex(1:LM, 6,n))
                    IF ( IJTS_TAUSUB(2,NTRIX_AOD(n),nsub_ntrix(     &
                         NTRIX_AOD(n)))>0 )                         &
                         TAIJS(i,j,IJTS_TAUSUB(2,NTRIX_AOD(n),      &
                         nsub_ntrix(NTRIX_AOD(n))))                 &
                         = TAIJS(i,j,IJTS_TAUSUB(2,NTRIX_AOD(n),    &
                         nsub_ntrix(NTRIX_AOD(n))))                 &
                         + SUM(aesqex(1:LM, 6,n))*OPNSKY
                    IF ( IJTS_TAUSUB(3,NTRIX_AOD(n),nsub_ntrix(     &
                         NTRIX_AOD(n)))>0 )                         &
                         TAIJS(i,j,IJTS_TAUSUB(3,NTRIX_AOD(n),      &
                         nsub_ntrix(NTRIX_AOD(n))))                 &
                         = TAIJS(i,j,IJTS_TAUSUB(3,NTRIX_AOD(n),    &
                         nsub_ntrix(NTRIX_AOD(n))))                 &
                         + SUM(aesqex_dry(1:LM, 6,n))
                 ELSE
                    DO kr = 1, 6
                       IF ( IJTS_SQEXSUB(1,kr,NTRIX_AOD(n),         &
                            nsub_ntrix(NTRIX_AOD(n)))>0 )           &
                            TAIJS(i,j,IJTS_SQEXSUB(1,kr,NTRIX_AOD(n)&
                            ,nsub_ntrix(NTRIX_AOD(n))))             &
                            = TAIJS(i,j,IJTS_SQEXSUB(1,kr,          &
                            NTRIX_AOD(n),nsub_ntrix(NTRIX_AOD(n)))) &
                            + SUM(aesqex(1:LM, kr,n))
                       IF ( IJTS_SQEXSUB(2,kr,NTRIX_AOD(n),         &
                            nsub_ntrix(NTRIX_AOD(n)))>0 )           &
                            TAIJS(i,j,IJTS_SQEXSUB(2,kr,NTRIX_AOD(n)&
                            ,nsub_ntrix(NTRIX_AOD(n))))             &
                            = TAIJS(i,j,IJTS_SQEXSUB(2,kr,          &
                            NTRIX_AOD(n),nsub_ntrix(NTRIX_AOD(n)))) &
                            + SUM(aesqex(1:LM, kr,n))*OPNSKY
                       IF ( IJTS_SQEXSUB(3,kr,NTRIX_AOD(n),         &
                            nsub_ntrix(NTRIX_AOD(n)))>0 )           &
                            TAIJS(i,j,IJTS_SQEXSUB(3,kr,NTRIX_AOD(n)&
                            ,nsub_ntrix(NTRIX_AOD(n))))             &
                            = TAIJS(i,j,IJTS_SQEXSUB(3,kr,          &
                            NTRIX_AOD(n),nsub_ntrix(NTRIX_AOD(n)))) &
                            + SUM(aesqex_dry(1:LM, kr,n))
                       IF ( IJTS_SQSCSUB(1,kr,NTRIX_AOD(n),         &
                            nsub_ntrix(NTRIX_AOD(n)))>0 )           &
                            TAIJS(i,j,IJTS_SQSCSUB(1,kr,NTRIX_AOD(n)&
                            ,nsub_ntrix(NTRIX_AOD(n))))             &
                            = TAIJS(i,j,IJTS_SQSCSUB(1,kr,          &
                            NTRIX_AOD(n),nsub_ntrix(NTRIX_AOD(n)))) &
                            + SUM(aesqsc(1:LM, kr,n))
                       IF ( IJTS_SQSCSUB(2,kr,NTRIX_AOD(n),         &
                            nsub_ntrix(NTRIX_AOD(n)))>0 )           &
                            TAIJS(i,j,IJTS_SQSCSUB(2,kr,NTRIX_AOD(n)&
                            ,nsub_ntrix(NTRIX_AOD(n))))             &
                            = TAIJS(i,j,IJTS_SQSCSUB(2,kr,          &
                            NTRIX_AOD(n),nsub_ntrix(NTRIX_AOD(n)))) &
                            + SUM(aesqsc(1:LM, kr,n))*OPNSKY
                       IF ( IJTS_SQSCSUB(3,kr,NTRIX_AOD(n),         &
                            nsub_ntrix(NTRIX_AOD(n)))>0 )           &
                            TAIJS(i,j,IJTS_SQSCSUB(3,kr,NTRIX_AOD(n)&
                            ,nsub_ntrix(NTRIX_AOD(n))))             &
                            = TAIJS(i,j,IJTS_SQSCSUB(3,kr,          &
                            NTRIX_AOD(n),nsub_ntrix(NTRIX_AOD(n)))) &
                            + SUM(aesqsc_dry(1:LM, kr,n))
                       IF ( IJTS_SQCBSUB(1,kr,NTRIX_AOD(n),         &
                            nsub_ntrix(NTRIX_AOD(n)))>0 )           &
                            TAIJS(i,j,IJTS_SQCBSUB(1,kr,NTRIX_AOD(n)&
                            ,nsub_ntrix(NTRIX_AOD(n))))             &
                            = TAIJS(i,j,IJTS_SQCBSUB(1,kr,          &
                            NTRIX_AOD(n),nsub_ntrix(NTRIX_AOD(n)))) &
                            + SUM(aesqcb(1:LM, kr,n))                &
                            /(SUM(aesqsc(1:LM, kr,n))+1.D-10)
                       IF ( IJTS_SQCBSUB(2,kr,NTRIX_AOD(n),         &
                            nsub_ntrix(NTRIX_AOD(n)))>0 )           &
                            TAIJS(i,j,IJTS_SQCBSUB(2,kr,NTRIX_AOD(n)&
                            ,nsub_ntrix(NTRIX_AOD(n))))             &
                            = TAIJS(i,j,IJTS_SQCBSUB(2,kr,          &
                            NTRIX_AOD(n),nsub_ntrix(NTRIX_AOD(n)))) &
                            + SUM(aesqcb(1:LM, kr,n))                &
                            /(SUM(aesqsc(1:LM, kr,n))+1.D-10)*OPNSKY
                       IF ( IJTS_SQCBSUB(3,kr,NTRIX_AOD(n),         &
                            nsub_ntrix(NTRIX_AOD(n)))>0 )           &
                            TAIJS(i,j,IJTS_SQCBSUB(3,kr,NTRIX_AOD(n)&
                            ,nsub_ntrix(NTRIX_AOD(n))))             &
                            = TAIJS(i,j,IJTS_SQCBSUB(3,kr,          &
                            NTRIX_AOD(n),nsub_ntrix(NTRIX_AOD(n)))) &
                            + SUM(aesqcb_dry(1:LM, kr,n))            &
                            /(SUM(aesqsc_dry(1:LM, kr,n))+1.D-10)
                    ENDDO
                 ENDIF
              CASE DEFAULT

                 ! 3d aod
                 IF ( diag_aod_3d>0 .AND. diag_aod_3d<5 ) THEN
                    ! valid values are 1-4
                    IF ( IJLT_3DAAOD(n)>0 )                         &
                         taijls(i,j,1:LM, IJLT_3DAAOD(n))            &
                         = taijls(i,j,1:LM, IJLT_3DAAOD(n))          &
                         + (aesqex(1:LM, 6,n)-aesqsc(1:LM, 6,n))
                    IF ( IJLT_3DAAODCS(n)>0 )                       &
                         taijls(i,j,1:LM, IJLT_3DAAODCS(n))          &
                         = taijls(i,j,1:LM, IJLT_3DAAODCS(n))        &
                         + (aesqex(1:LM, 6,n)-aesqsc(1:LM, 6,n))      &
                         *OPNSKY
                    IF ( IJLT_3DAAODDRY(n)>0 )                      &
                         taijls(i,j,1:LM, IJLT_3DAAODDRY(n))         &
                         = taijls(i,j,1:LM, IJLT_3DAAODDRY(n))       &
                         + (aesqex_dry(1:LM, 6,n)-aesqsc(1:LM, 6,n))
                    IF ( IJLT_3DTAU(n)>0 )                          &
                         taijls(i,j,1:LM, IJLT_3DTAU(n))             &
                         = taijls(i,j,1:LM, IJLT_3DTAU(n))           &
                         + aesqex(1:LM, 6,n)
                    IF ( IJLT_3DTAUCS(n)>0 )                        &
                         taijls(i,j,1:LM, IJLT_3DTAUCS(n))           &
                         = taijls(i,j,1:LM, IJLT_3DTAUCS(n))         &
                         + aesqex(1:LM, 6,n)*OPNSKY
                    IF ( IJLT_3DTAUDRY(n)>0 )                       &
                         taijls(i,j,1:LM, IJLT_3DTAUDRY(n))          &
                         = taijls(i,j,1:LM, IJLT_3DTAUDRY(n))        &
                         + aesqex_dry(1:LM, 6,n)
                 ELSEIF ( diag_aod_3d<0 .AND. diag_aod_3d>-5 ) THEN
                    ! if negative, save total
                    IF ( IJLT_3DAAOD(1)>0 )                         &
                         taijls(i,j,1:LM, IJLT_3DAAOD(1))            &
                         = taijls(i,j,1:LM, IJLT_3DAAOD(1))          &
                         + (aesqex(1:LM, 6,n)-aesqsc(1:LM, 6,n))
                    IF ( IJLT_3DAAODCS(1)>0 )                       &
                         taijls(i,j,1:LM, IJLT_3DAAODCS(1))          &
                         = taijls(i,j,1:LM, IJLT_3DAAODCS(1))        &
                         + (aesqex(1:LM, 6,n)-aesqsc(1:LM, 6,n))      &
                         *OPNSKY
                    IF ( IJLT_3DAAODDRY(1)>0 )                      &
                         taijls(i,j,1:LM, IJLT_3DAAODDRY(1))         &
                         = taijls(i,j,1:LM, IJLT_3DAAODDRY(1))       &
                         + (aesqex_dry(1:LM, 6,n)                    &
                         -aesqsc_dry(1:LM, 6,n))
                    IF ( IJLT_3DTAU(1)>0 )                          &
                         taijls(i,j,1:LM, IJLT_3DTAU(1))             &
                         = taijls(i,j,1:LM, IJLT_3DTAU(1))           &
                         + aesqex(1:LM, 6,n)
                    IF ( IJLT_3DTAUCS(1)>0 )                        &
                         taijls(i,j,1:LM, IJLT_3DTAUCS(1))           &
                         = taijls(i,j,1:LM, IJLT_3DTAUCS(1))         &
                         + aesqex(1:LM, 6,n)*OPNSKY
                    IF ( IJLT_3DTAUDRY(1)>0 )                       &
                         taijls(i,j,1:LM, IJLT_3DTAUDRY(1))          &
                         = taijls(i,j,1:LM, IJLT_3DTAUDRY(1))        &
                         + aesqex_dry(1:LM, 6,n)
                 ENDIF
                 ! 0<diag_aod_3d<5 or 0>diag_aod_3d>-5

                 ! 2d aod, per band or just band6, depending on diag_rad
                 IF ( diag_rad/=1 ) THEN
                    IF ( IJTS_TAU(1,NTRIX_AOD(n))>0 )               &
                         TAIJS(i,j,IJTS_TAU(1,NTRIX_AOD(n)))        &
                         = TAIJS(i,j,IJTS_TAU(1,NTRIX_AOD(n)))      &
                         + SUM(aesqex(1:LM, 6,n))
                    IF ( IJTS_TAU(2,NTRIX_AOD(n))>0 )               &
                         TAIJS(i,j,IJTS_TAU(2,NTRIX_AOD(n)))        &
                         = TAIJS(i,j,IJTS_TAU(2,NTRIX_AOD(n)))      &
                         + SUM(aesqex(1:LM, 6,n))*OPNSKY
                    IF ( IJTS_TAU(3,NTRIX_AOD(n))>0 )               &
                         TAIJS(i,j,IJTS_TAU(3,NTRIX_AOD(n)))        &
                         = TAIJS(i,j,IJTS_TAU(3,NTRIX_AOD(n)))      &
                         + SUM(aesqex_dry(1:LM, 6,n))
                 ELSE
                    DO kr = 1, 6
                       !               print*,'SUSA  diag',SUM(aesqex(1:LM, kr,n))
                       IF ( IJTS_SQEX(1,kr,NTRIX_AOD(n))>0 )        &
                            TAIJS(i,j,IJTS_SQEX(1,kr,NTRIX_AOD(n))) &
                            = TAIJS(i,j,IJTS_SQEX(1,kr,NTRIX_AOD(n))&
                            ) + SUM(aesqex(1:LM, kr,n))
                       IF ( IJTS_SQEX(2,kr,NTRIX_AOD(n))>0 )        &
                            TAIJS(i,j,IJTS_SQEX(2,kr,NTRIX_AOD(n))) &
                            = TAIJS(i,j,IJTS_SQEX(2,kr,NTRIX_AOD(n))&
                            ) + SUM(aesqex(1:LM, kr,n))*OPNSKY
                       IF ( IJTS_SQEX(3,kr,NTRIX_AOD(n))>0 )        &
                            TAIJS(i,j,IJTS_SQEX(3,kr,NTRIX_AOD(n))) &
                            = TAIJS(i,j,IJTS_SQEX(3,kr,NTRIX_AOD(n))&
                            ) + SUM(aesqex_dry(1:LM, kr,n))
                       IF ( IJTS_SQSC(1,kr,NTRIX_AOD(n))>0 )        &
                            TAIJS(i,j,IJTS_SQSC(1,kr,NTRIX_AOD(n))) &
                            = TAIJS(i,j,IJTS_SQSC(1,kr,NTRIX_AOD(n))&
                            ) + SUM(aesqsc(1:LM, kr,n))
                       IF ( IJTS_SQSC(2,kr,NTRIX_AOD(n))>0 )        &
                            TAIJS(i,j,IJTS_SQSC(2,kr,NTRIX_AOD(n))) &
                            = TAIJS(i,j,IJTS_SQSC(2,kr,NTRIX_AOD(n))&
                            ) + SUM(aesqsc(1:LM, kr,n))*OPNSKY
                       IF ( IJTS_SQSC(3,kr,NTRIX_AOD(n))>0 )        &
                            TAIJS(i,j,IJTS_SQSC(3,kr,NTRIX_AOD(n))) &
                            = TAIJS(i,j,IJTS_SQSC(3,kr,NTRIX_AOD(n))&
                            ) + SUM(aesqsc_dry(1:LM, kr,n))
#ifndef TRACERS_TOMAS
                       IF ( IJTS_SQCB(1,kr,NTRIX_AOD(n))>0 )        &
                            TAIJS(i,j,IJTS_SQCB(1,kr,NTRIX_AOD(n))) &
                            = TAIJS(i,j,IJTS_SQCB(1,kr,NTRIX_AOD(n))&
                            ) + SUM(aesqcb(1:LM, kr,n))            &
                            /(SUM(aesqsc(1:LM, kr,n))+1.D-10)
                       IF ( IJTS_SQCB(2,kr,NTRIX_AOD(n))>0 )        &
                            TAIJS(i,j,IJTS_SQCB(2,kr,NTRIX_AOD(n))) &
                            = TAIJS(i,j,IJTS_SQCB(2,kr,NTRIX_AOD(n))&
                            ) + SUM(aesqcb(1:LM, kr,n))            &
                            /(SUM(aesqsc(1:LM, kr,n))+1.D-10)      &
                            *OPNSKY
                       IF ( IJTS_SQCB(3,kr,NTRIX_AOD(n))>0 )        &
                            TAIJS(i,j,IJTS_SQCB(3,kr,NTRIX_AOD(n))) &
                            = TAIJS(i,j,IJTS_SQCB(3,kr,NTRIX_AOD(n))&
                            ) + SUM(aesqcb_dry(1:LM, kr,n))        &
                            /(SUM(aesqsc_dry(1:LM, kr,n))+1.D-10)
#else
                       qcb_col(kr,n) = 0.D0
                       qcb_col_dry(kr,n) = 0.D0
                       DO l = 1, lm
                          qcb_col(kr,n) = qcb_col(kr,n)             &
                               + aesqcb(l,kr,n)*aesqsc(l,kr,n)
                          qcb_col_dry(kr,n) = qcb_col(kr,n)         &
                               + aesqcb_dry(l,kr,n)*aesqsc_dry(l,kr,n)
                       ENDDO

                       IF ( IJTS_SQCB(1,kr,NTRIX_AOD(n))>0 )        &
                            TAIJS(i,j,IJTS_SQCB(1,kr,NTRIX_AOD(n))) &
                            = TAIJS(i,j,IJTS_SQCB(1,kr,NTRIX_AOD(n))&
                            ) + qcb_col(kr,n)                     &
                            /(SUM(aesqsc(1:LM, kr,n))+1.D-10)
                       IF ( IJTS_SQCB(2,kr,NTRIX_AOD(n))>0 )        &
                            TAIJS(i,j,IJTS_SQCB(2,kr,NTRIX_AOD(n))) &
                            = TAIJS(i,j,IJTS_SQCB(2,kr,NTRIX_AOD(n))&
                            ) + qcb_col(kr,n)                     &
                            /(SUM(aesqsc(1:LM, kr,n))+1.D-10)      &
                            *OPNSKY
                       IF ( IJTS_SQCB(3,kr,NTRIX_AOD(n))>0 )        &
                            TAIJS(i,j,IJTS_SQCB(3,kr,NTRIX_AOD(n))) &
                            = TAIJS(i,j,IJTS_SQCB(3,kr,NTRIX_AOD(n))&
                            ) + qcb_col_dry(kr,n)                 &
                            /(SUM(aesqsc_dry(1:LM, kr,n))+1.D-10)
#endif
                    ENDDO
                    ! kr
                 ENDIF
                 ! diag_rad
              ENDSELECT
              ! clay or not
           ENDDO
           ! nraero_aod

#endif  /* Koch||DUST||MINERALS||AMP||TOMAS||SEASALT */

           IF ( TAero_aod_diag>0 ) THEN
              DO n = 1, 8
                 ! 8 radiatively active aerosol tracers
                 DO kr = 1, 6
                    ! 6 bands in the shortwave
                    IF ( TAero_aod_diag==2 .AND. kr/=6 ) CYCLE
                    ! only save band6
                    AIJ(i,j,IJ_NINTAEREXT(kr,n))                    &
                         = AIJ(i,j,IJ_NINTAEREXT(kr,n))               &
                         + SUM(nintaerext(1:LM, kr,n))
                    AIJ(i,j,IJ_NINTAERSCA(kr,n))                    &
                         = AIJ(i,j,IJ_NINTAERSCA(kr,n))               &
                         + SUM(nintaersca(1:LM, kr,n))
                    AIJ(i,j,IJ_NINTAERASY(kr,n))                    &
                         = AIJ(i,j,IJ_NINTAERASY(kr,n))               &
                         + SUM(nintaerasy(1:LM, kr,n)                  &
                         *nintaersca(1:LM, kr,n))                      &
                         /(SUM(nintaersca(1:LM, kr,n))+1.D-10)
                 ENDDO
                 ! kr
              ENDDO
              ! n
           ENDIF


#ifdef TRACERS_ON
           IF ( nraero_aod>0 ) THEN
              tau_as(i,j,1:LM, 1 : nraero_aod)                         &
                   = aesqex(1:LM, 6,1 : nraero_aod)
              tau_cs(i,j,1:LM, 1 : nraero_aod)                         &
                   = aesqex(1:LM, 6,1 : nraero_aod)*OPNSKY
              IF ( save_dry_aod>0 ) tau_dry(i,j,1:LM, 1 : nraero_aod)  &
                   = aesqex_dry(1:LM, 6,1 : nraero_aod)
#ifdef CACHED_SUBDD
              abstau_as(i,j,1:LM, 1 : nraero_aod)                      &
                   = (aesqex(1:LM, 6,1 : nraero_aod)                     &
                   -aesqsc(1:LM, 6,1 : nraero_aod))
              abstau_cs(i,j,1:LM, 1 : nraero_aod)                      &
                   = (aesqex(1:LM, 6,1 : nraero_aod)                     &
                   -aesqsc(1:LM, 6,1 : nraero_aod))*OPNSKY
              IF ( save_dry_aod>0 )                                 &
                   abstau_dry(i,j,1:LM, 1 : nraero_aod)                &
                   = (aesqex_dry(1:LM, 6,1 : nraero_aod)               &
                   -aesqsc_dry(1:LM, 6,1 : nraero_aod))
#endif  /* CACHED_SUBDD */
           ENDIF
#endif /* TRACERS_ON */

           IF ( I==IWRITE .AND. J==JWRITE ) CALL WRITER(6,ITWRITE)
           CSZ2 = COSZ2(I,J)
           DO L = 1, LM
#ifdef GCC_COUPLE_RAD
              GCCCO2RAD_TO_CHEM(L,i,j) = GCCCO2_OUT(L)
#endif
              rad_to_chem( : ,L,i,j) = chem_out(L, : )
              rad_to_chem(4,L,i,j) = chem_out(L,4)/CH4X_RADoverCHEM
              DO k = 1, 4
                 kliq(L,k,i,j) = kdeliq(L,k)
                 ! save updated flags
              ENDDO
           ENDDO
           IF ( kradia>0 ) THEN
              ! rad. forc. model; acc diagn
              DO L = 1, LM + LM_REQ + 1
                 AFLX_ST(L,I,J,1) = AFLX_ST(L,I,J,1) + SRUFLB(L)    &
                      *CSZ2
                 AFLX_ST(L,I,J,2) = AFLX_ST(L,I,J,2) + SRDFLB(L)    &
                      *CSZ2
                 AFLX_ST(L,I,J,3) = AFLX_ST(L,I,J,3) + TRUFLB(L)
                 AFLX_ST(L,I,J,4) = AFLX_ST(L,I,J,4) + TRDFLB(L)
              ENDDO
              IF ( kradia==1 ) THEN
                 tauex6 = 0.
                 tauex5 = 0.
                 tausct = 0.
                 taugcb = 0.
                 DO L = 1, LM
                    AFLX_ST(L,I,J,5) = AFLX_ST(L,I,J,5)             &
                         + 1.D2*RHL(L)
                    tauex6 = tauex6 + SRAEXT(L,6) + SRDEXT(L,6)     &
                         + SRVEXT(L,6)
                    tauex5 = tauex5 + SRAEXT(L,5) + SRDEXT(L,5)     &
                         + SRVEXT(L,5)
                    tausct = tausct + SRASCT(L,6) + SRDSCT(L,6)     &
                         + SRVSCT(L,6)
                    taugcb = taugcb + SRASCT(L,6)*SRAGCB(L,6)       &
                         + SRDSCT(L,6)*SRDGCB(L,6) + SRVSCT(L,6)&
                         *SRVGCB(L,6)
                 ENDDO
                 AFLX_ST(LM+1,I,J,5) = AFLX_ST(LM+1,I,J,5) + tauex5
                 AFLX_ST(LM+2,I,J,5) = AFLX_ST(LM+2,I,J,5) + tauex6
                 AFLX_ST(LM+3,I,J,5) = AFLX_ST(LM+3,I,J,5) + tausct
                 AFLX_ST(LM+4,I,J,5) = AFLX_ST(LM+4,I,J,5) + taugcb
                 CYCLE
              ENDIF
              DO l = LS1_loc, lm
                 TCHG(l,i,j) = TCHG(l,i,j)                          &
                      + (SRFHRL(l)*csz2-srhra(l,i,j)       &
                      +(-TRFCRL(l)-trhra(l,i,j)))          &
                      *nrad*DTsrc*bysha*BYMA(l,i,j)
              ENDDO
              DO l = lm + 1, lm + lm_req
                 TCHG(l,i,j) = TCHG(l,i,j)                          &
                      + (SRFHRL(l)*csz2-srhra(l,i,j)       &
                      +(-TRFCRL(l)-trhra(l,i,j)))          &
                      *nrad*DTsrc*bysha*BYAML00(l)
              ENDDO
              CYCLE
           ELSEIF ( kradia<0 ) THEN
              ! save i/o data for frc.runs
              fmp_com(i,j) = fmp        ! input data
              wsoil(i,j) = wearth
              DO L = 1, LM
                 QR(L,I,J) = SHL(L)
                 CLDinfo(L,1,I,J) = TAUWC(L)
                 CLDinfo(L,2,I,J) = TAUIC(L)
                 CLDinfo(L,3,I,J) = SIZEIC(L)
                 ! sizeic=sizewc currently
              ENDDO
              SRHRA(0,I,J) = SRNFLB(1)*CSZ2
              ! output data (for adj frc)
              TRHRA(0,I,J) = -TRNFLB(1)
              DO L = 1, LM + LM_REQ
                 SRHRA(L,I,J) = SRFHRL(L)*CSZ2
                 TRHRA(L,I,J) = -TRFCRL(L)
              ENDDO
           ENDIF
           !****
           !**** Save relevant output in model arrays
           !****
           !**** (some generalisation and coherence needed in the rad surf type calc)
           FSF(1,I,J) = FSRNFG(1)
           !  ocean
           FSF(2,I,J) = FSRNFG(3)
           !  ocean ice
           FSF(3,I,J) = FSRNFG(4)
           !  land ice
           FSF(4,I,J) = FSRNFG(2)
           !  soil
           SRHR(0,I,J) = SRNFLB(1)
           TRHR(0,I,J) = STBO*(POCEAN*atmocn%GTEMPR(I,J)**4+POICE*  &
                atmice%GTEMPR(I,J)                         &
                **4+PLICE*atmgla%GTEMPR(I,J)               &
                **4+PEARTH*atmlnd%GTEMPR(I,J)**4)          &
                - TRNFLB(1)
           TRSURF(1,I,J) = STBO*atmocn%GTEMPR(I,J)**4
           !  ocean
           TRSURF(2,I,J) = STBO*atmice%GTEMPR(I,J)**4
           !  ocean ice
           TRSURF(3,I,J) = STBO*atmgla%GTEMPR(I,J)**4
           !  land ice
           TRSURF(4,I,J) = STBO*atmlnd%GTEMPR(I,J)**4
           !  soil
           DO L = 1, LM
              SRHR(L,I,J) = SRFHRL(L)
              TRHR(L,I,J) = -TRFCRL(L)
           ENDDO
           DO LR = 1, LM_REQ
              SRHRS(LR,I,J) = SRFHRL(LM+LR)
              TRHRS(LR,I,J) = -TRFCRL(LM+LR)
           ENDDO
#ifdef SCM
           !**** possibly turn off radiative heating in atmosphere
           !**** and use specified profile for thermal heating rate
           !**** converting units from K/s to W/m2
           IF ( SCMopt%QRAD ) THEN
              SRHR(1:LM, I,J) = 0.
              TRHR(1:LM, I,J) = 0.
              SRHRS(1:LM_REQ,I,J) = 0.
              TRHRS(1:LM_REQ,I,J) = 0.
              TRHR(1:LM, I,J) = SCMin%QRAD(1:LM)*SHA*MA(1:LM, I,J)
           ENDIF
           !**** possibly turn off radiative heating in atmosphere
           !**** and use Beers Law for thermal heating rate as
           !**** difference of net flux over layer
           IF ( SCMopt%BEERSLAW ) THEN
              SRHR(1:LM, I,J) = 0.
              TRHR(1:LM, I,J) = 0.
              SRHRS(1:LM_REQ,I,J) = 0.
              TRHRS(1:LM_REQ,I,J) = 0.
              ! cumulative cloud water paths * extinction coefficient
              q_above(LM+1) = 0.
              DO L = LM, 1, -1
                 q_above(L) = q_above(L+1)                          &
                      + SCMin%BEERSLAW_KAPPA*MA(L,i,j)      &
                      *QCL(i,j,L)
              ENDDO
              q_below(1) = 0.
              DO L = 1, LM
                 q_below(L+1) = q_below(L)                          &
                      + SCMin%BEERSLAW_KAPPA*MA(L,i,j)    &
                      *QCL(i,j,L)
              ENDDO
              ! net upward radiative flux at layer edges
              Frad( : ) = SCMin%BEERSLAW_F0*EXP(-q_above( : ))          &
                   + SCMin%BEERSLAW_F1*EXP(-q_below( : ))
              ! radiative flux difference over each layer
              TRHR(1:LM, I,J) = Frad(1:LM) - Frad(2 : LM+1)
           ENDIF
           !**** save radiative flux profiles for sub-daily output
#ifdef CACHED_SUBDD
           TRDFLB_prof(I,J,1:LM) = TRDFLB(1:LM)
           TRUFLB_prof(I,J,1:LM) = TRUFLB(1:LM)
           SRDFLB_prof(I,J,1:LM) = SRDFLB(1:LM)
           SRUFLB_prof(I,J,1:LM) = SRUFLB(1:LM)
#endif
#endif
           !**** Save fluxes at four levels surface, P0, P1, LTROPO
           ! Surface
           SNFS(1,I,J) = SRNFLB(1)
           TNFS(1,I,J) = TRNFLB(1)
           ! P1
           SNFS(2,I,J) = SRNFLB(LM+1)
           TNFS(2,I,J) = TRNFLB(LM+1)
           ! P0 = TOA
           SNFS(3,I,J) = SRNFLB(LM+LM_REQ+1)
           TNFS(3,I,J) = TRNFLB(LM+LM_REQ+1)
           ! LTROPO
           SNFS(4,I,J) = SRNFLB(LTROPO(I,J))
           TNFS(4,I,J) = TRNFLB(LTROPO(I,J))

#ifdef TRACERS_GC
           SNFS_3D(I,J,:) = SRNFLB(:)
           TNFS_3D(I,J,:) = TRNFLB(:)
											! Archive total flux for diagnostis
           SAVE_RF_3D(I,J,:,0,1) = SRNFLB(:) 
           SAVE_RF_3D(I,J,:,0,2) = TRNFLB(:)
#endif

           !****
           TRINCG(I,J) = TRDFLB(1)
           BTMPW(I,J) = BTEMPW - TF
           ALB(I,J,1) = SRNFLB(1)/(SRDFLB(1)+1.D-20)
           ALB(I,J,2) = PLAVIS
           ALB(I,J,3) = PLANIR
           ALB(I,J,4) = ALBVIS
           ALB(I,J,5) = ALBNIR
           ALB(I,J,6) = SRRVIS
           ALB(I,J,7) = SRRNIR
           ALB(I,J,8) = SRAVIS
           ALB(I,J,9) = SRANIR

#ifdef TRACERS_DUST
           IF ( adiurn_dust==1 ) THEN
              srnflb_save(i,j,1:LM) = SRNFLB(1:LM)
              trnflb_save(i,j,1:LM) = TRNFLB(1:LM)
           ENDIF
#endif
#ifdef mjo_subdd
           SWU_AVG(I,J) = SWU_AVG(I,J) + SRUFLB(1)*CSZ2
#endif

           SWUS(I,J) = SRUFLB(1)*CSZ2
#ifdef CFMIP3_SUBDD
           ! SW upward flux at TOA
           swut(i,j) = SRUFLB(lm)*csz2
           ! SW downward flux at TOA
           swdt(i,j) = SRDFLB(lm)*csz2
#endif
           SRDN(I,J) = SRDFLB(1)
           ! save total solar flux at surface
           !**** SALB(I,J)=ALB(I,J,1)      ! save surface albedo (pointer)
           FSRDIR(I,J) = SRXVIS
           ! direct visible solar at surface **coefficient
           SRVISSURF(I,J) = SRDVIS
           ! total visible solar at surface
           DIRVIS(I,J) = SRXVIS*SRDVIS
           ! direct visible solar at surface
           FSRDIF(I,J) = SRDVIS*(1-SRXVIS)
           ! diffuse visible solar at surface

           DIRNIR(I,J) = SRXNIR*SRDNIR
           ! direct beam nir solar at surface
           DIFNIR(I,J) = SRDNIR*(1-SRXNIR)
           ! diffuse     nir solar at surface

           !diag write(*,'(a,2i5,6e12.4)')'RAD_DRV :  ',
           !diag.    I,J,FSRDIR(I,J),SRVISSURF(I,J),FSRDIF(I,J),
           !diag.        DIRNIR(I,J),SRDNIR,DIFNIR(I,J)
           !**** Save clear sky/tropopause diagnostics here
           AIJ(I,J,IJ_CLR_SRINCG) = AIJ(I,J,IJ_CLR_SRINCG)          &
                + OPNSKY*SRDFLB(1)*CSZ2
           AIJ(I,J,IJ_CLR_SRNFG) = AIJ(I,J,IJ_CLR_SRNFG)            &
                + OPNSKY*SRNFLB(1)*CSZ2
           AIJ(I,J,IJ_CLR_TRDNG) = AIJ(I,J,IJ_CLR_TRDNG)            &
                + OPNSKY*TRHR(0,I,J)
           AIJ(I,J,IJ_CLR_SRUPTOA) = AIJ(I,J,IJ_CLR_SRUPTOA)        &
                + OPNSKY*SRUFLB(LM+LM_REQ+1)*CSZ2
           AIJ(I,J,IJ_CLR_TRUPTOA) = AIJ(I,J,IJ_CLR_TRUPTOA)        &
                + OPNSKY*TRUFLB(LM+LM_REQ+1)
           AIJ(I,J,IJ_CLR_SRNTP) = AIJ(I,J,IJ_CLR_SRNTP)            &
                + OPNSKY*SRNFLB(LTROPO(I,J))*CSZ2
           AIJ(I,J,IJ_CLR_TRNTP) = AIJ(I,J,IJ_CLR_TRNTP)            &
                + OPNSKY*TRNFLB(LTROPO(I,J))
           AIJ(I,J,IJ_SRNTP) = AIJ(I,J,IJ_SRNTP)                    &
                + SRNFLB(LTROPO(I,J))*CSZ2
           AIJ(I,J,IJ_TRNTP) = AIJ(I,J,IJ_TRNTP)                    &
                + TRNFLB(LTROPO(I,J))
           AIJ(I,J,IJ_SISWD) = AIJ(I,J,IJ_SISWD) + POICE*SRDFLB(1)  &
                *CSZ2
           AIJ(I,J,IJ_SISWU) = AIJ(I,J,IJ_SISWU)                    &
                + POICE*(SRDFLB(1)-FSRNFG(3))*CSZ2

           DO IT = 1, NTYPE
              CALL INC_AJ(i,j,it,J_CLRTOA,                          &
                   OPNSKY*(SRNFLB(LM+LM_REQ+1)               &
                   *CSZ2-TRNFLB(LM+LM_REQ+1))*FTYPE(IT,I,J))
              CALL INC_AJ(i,j,it,J_CLRTRP,                          &
                   OPNSKY*(SRNFLB(LTROPO(I,J))               &
                   *CSZ2-TRNFLB(LTROPO(I,J)))*FTYPE(IT,I,J))
              CALL INC_AJ(i,j,it,J_TOTTRP,                          &
                   (SRNFLB(LTROPO(I,J))*CSZ2-TRNFLB          &
                   (LTROPO(I,J)))*FTYPE(IT,I,J))
           ENDDO
           CALL INC_AREG(i,j,jr,J_CLRTOA,                           &
                OPNSKY*(SRNFLB(LM+LM_REQ+1)*CSZ2-          &
                TRNFLB(LM+LM_REQ+1)))
           CALL INC_AREG(i,j,jr,J_CLRTRP,                           &
                OPNSKY*(SRNFLB(LTROPO(I,J))*CSZ2-          &
                TRNFLB(LTROPO(I,J))))
           CALL INC_AREG(i,j,jr,J_TOTTRP,                           &
                (SRNFLB(LTROPO(I,J))*CSZ2-TRNFLB           &
                (LTROPO(I,J))))
           !**** Save cloud top diagnostics here
           IF ( CLDCV>0. ) THEN
              AIJ(I,J,IJ_CLDTPPR) = AIJ(I,J,IJ_CLDTPPR)             &
                   + PLB(ltopcl+1)
              AIJ(I,J,IJ_CLDTPT) = AIJ(I,J,IJ_CLDTPT)               &
                   + (TLB(ltopcl+1)-tf)
              CTT(i,j) = (TLB(ltopcl+1)-tf)
              CTP(i,j) = PLB(ltopcl+1)
              !**** Save cloud tau=1 related diagnostics here (opt.depth=1 level)
              tauup = 0.
              DO L = LM, 1, -1
                 taucl = TAUWC(l) + TAUIC(l)
                 taudn = tauup + taucl
                 IF ( taudn>1. ) THEN
                    AIJ(i,j,ij_cldcv1) = AIJ(i,j,ij_cldcv1) + 1.
                    wtlin = (1.-tauup)/taucl
                    AIJ(i,j,ij_cldt1t) = AIJ(i,j,ij_cldt1t)         &
                         + (TLB(l+1)-tf+(TLB(l)-TLB(l+1))*wtlin)
                    AIJ(i,j,ij_cldt1p) = AIJ(i,j,ij_cldt1p)         &
                         + (PLB(l+1)+(PLB(l)-PLB(l+1))*wtlin)
                    EXIT
                 ENDIF
                 tauup = taudn
              ENDDO
           ENDIF

        ENDDO
        !****
        !**** END OF MAIN LOOP FOR I INDEX
        !****

     ENDDO
     !****
     !**** END OF MAIN LOOP FOR J INDEX
     !****

#ifdef mjo_subdd
     swu_cnt = swu_cnt + 1.
#endif

     IF ( kradia>0 ) THEN
        CALL STOPTIMER('RADIA()')
        RETURN
     ENDIF
     !**** Stop if temperatures were out of range
     !**** Now only warning messages are printed for T,Q errors
     !     IF(ICKERR.GT.0)
     !         call stop_model('In Radia :  Temperature out of range',11)
     !     IF(JCKERR.GT.0)  call stop_model('In Radia :  RQT out of range',11)
     !     IF(KCKERR.GT.0)  call stop_model('In Radia :  Q<0',255)
     !**** save all input data to disk if kradia<0
     ! LM+LM_REQ+1+
     !         ,(((GTEMPR(k,i,j),k=1,4),i=1,im),j=1,jm) ! (4+)
     ! LM+1+3*LM+1+1+
     ! 1+1+1+1+1+
     ! 3+1+.5+.5+
     !**** output data :  really needed only if kradia=2
     ! 2+1+1
     IF ( kradia<0 ) WRITE (iu_rad) itime, T, RQT, atmsrf%TSAVG, QR,&
          P, CLDinfo, rsi, zsi, wsoil,    &
          atmsrf%WSAVG, snowi,            &
          atmgla%SNOW, atmlnd%SNOWE,      &
          snoage, fmp_com, flag_dsws,     &
          ltropo, atmlnd%FR_SNOW_RAD,     &
          dlake, flake, srhra, trhra,     &
          itime
     ! 2(LM+LM_REQ+1)
     !****
     !**** ACCUMULATE THE RADIATION DIAGNOSTICS
     !****
     bydpreq( : ) = 1D0/(req_fac_d( : )*pmtop)
     DO J = J_0, J_1
        DO I = I_0, IMAXJ(J)
           DO l = 1, lm
              CALL INC_AJL(i,j,l,jl_srhr,SRHR(L,I,J)*COSZ2(I,J))
              CALL INC_AJL(i,j,l,jl_trcr,TRHR(L,I,J))
           ENDDO
           CSZ2 = COSZ2(I,J)
           JR = JREG(I,J)
           DO LR = 1, LM_REQ
              CALL INC_ASJL(i,j,lr,3,bydpreq(lr)*SRHRS(LR,I,J)*CSZ2)
              CALL INC_ASJL(i,j,lr,4,bydpreq(lr)*TRHRS(LR,I,J))
           ENDDO
           DO KR = 1, NDIUPT
              IF ( I==IJDD(1,KR) .AND. J==IJDD(2,KR) ) THEN
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
                 TMP(idd_aot) = SUM(aesqex(1:LM, 6,1 : nraero_aod))
                 !*OPNSKY
                 TMP(idd_aot2) = SUM(aesqsc(1:LM, 6,1 : nraero_aod))
                 !*OPNSKY
#endif
                 TMP(IDD_PALB) = (1.-SNFS(3,I,J)/S0)
                 TMP(IDD_GALB) = (1.-ALB(I,J,1))
                 TMP(IDD_ABSA) = (SNFS(3,I,J)-SRHR(0,I,J))*CSZ2
                 DO INCH = 1, NRAD
                    IHM = 1 + (JTIME+INCH-1)*HR_IN_DAY/NDAY
                    IH = IHM
                    IF ( IH>HR_IN_DAY ) IH = IH - HR_IN_DAY
                    ADIURN(IDXB( : ),KR,IH) = ADIURN(IDXB( : ),KR,IH)   &
                         + TMP(IDXB( : ))
#ifdef USE_HDIURN
                    IHM = IHM + (DATE-1)*HR_IN_DAY
                    IF ( IHM<=HR_IN_MONTH ) HDIURN(IDXB( : ),KR,IHM)  &
                         = HDIURN(IDXB( : ),KR,IHM) + TMP(IDXB( : ))
#endif
                 ENDDO
              ENDIF
           ENDDO

           DO IT = 1, NTYPE
              CALL INC_AJ(I,J,IT,J_SRINCP0,(S0*CSZ2)*FTYPE(IT,I,J))
              CALL INC_AJ(I,J,IT,J_SRNFP0,(SNFS(3,I,J)*CSZ2)        &
                   *FTYPE(IT,I,J))
              CALL INC_AJ(I,J,IT,J_SRINCG,                          &
                   (SRHR(0,I,J)*CSZ2/(ALB(I,J,1)+1.D-20))    &
                   *FTYPE(IT,I,J))
              CALL INC_AJ(I,J,IT,J_BRTEMP,BTMPW(I,J)*FTYPE(IT,I,J))
              CALL INC_AJ(I,J,IT,J_TRINCG,TRINCG(I,J)*FTYPE(IT,I,J))
              CALL INC_AJ(I,J,IT,J_HSURF,-(TNFS(3,I,J)-TNFS(1,I,J)) &
                   *FTYPE(IT,I,J))
              CALL INC_AJ(I,J,IT,J_TRNFP0,-TNFS(3,I,J)*FTYPE(IT,I,J)&
                   )
              CALL INC_AJ(I,J,IT,J_TRNFP1,-TNFS(2,I,J)*FTYPE(IT,I,J)&
                   )
              CALL INC_AJ(I,J,IT,J_SRNFP1,SNFS(2,I,J)               &
                   *CSZ2*FTYPE(IT,I,J))
              CALL INC_AJ(I,J,IT,J_HATM,-(TNFS(2,I,J)-TNFS(1,I,J))  &
                   *FTYPE(IT,I,J))
#ifdef HEALY_LM_DIAGS
              CALL INC_AJ(I,J,IT,j_vtau,10.*-20*VTAULAT(J)          &
                   *FTYPE(IT,I,J))
              CALL INC_AJ(I,J,IT,j_ghg,10.*ghg_totforc*FTYPE(IT,I,J)&
                   )
#endif

           ENDDO
           !**** Note :  confusing because the types for radiation are a subset
           CALL INC_AJ(I,J,ITOCEAN,J_SRNFG,(FSF(1,I,J)*CSZ2)        &
                *FOCEAN(I,J)*(1.-RSI(I,J)))
           CALL INC_AJ(I,J,ITLAKE,J_SRNFG,(FSF(1,I,J)*CSZ2)         &
                *FLAKE(I,J)*(1.-RSI(I,J)))
           CALL INC_AJ(I,J,ITEARTH,J_SRNFG,(FSF(4,I,J)*CSZ2)        &
                *FEARTH(I,J))
           CALL INC_AJ(I,J,ITLANDI,J_SRNFG,(FSF(3,I,J)*CSZ2)        &
                *FLICE(I,J))
           CALL INC_AJ(I,J,ITOICE,J_SRNFG,(FSF(2,I,J)*CSZ2)         &
                *FOCEAN(I,J)*RSI(I,J))
           CALL INC_AJ(I,J,ITLKICE,J_SRNFG,(FSF(2,I,J)*CSZ2)        &
                *FLAKE(I,J)*RSI(I,J))
           !****
           CALL INC_AREG(I,J,JR,J_SRINCP0,(S0*CSZ2))
           CALL INC_AREG(I,J,JR,J_SRNFP0,(SNFS(3,I,J)*CSZ2))
           CALL INC_AREG(I,J,JR,J_SRNFP1,(SNFS(2,I,J)*CSZ2))
           CALL INC_AREG(I,J,JR,J_SRINCG,                           &
                (SRHR(0,I,J)*CSZ2/(ALB(I,J,1)+1.D-20)))
           CALL INC_AREG(I,J,JR,J_HATM,-(TNFS(2,I,J)-TNFS(1,I,J)))
           CALL INC_AREG(I,J,JR,J_SRNFG,(SRHR(0,I,J)*CSZ2))
           CALL INC_AREG(I,J,JR,J_HSURF,-(TNFS(3,I,J)-TNFS(1,I,J)))
           CALL INC_AREG(I,J,JR,J_BRTEMP,BTMPW(I,J))
           CALL INC_AREG(I,J,JR,J_TRINCG,TRINCG(I,J))
           CALL INC_AREG(I,J,JR,J_TRNFP0,-TNFS(3,I,J))
           CALL INC_AREG(I,J,JR,J_TRNFP1,-TNFS(2,I,J))
           DO K = 2, 9
              JK = AJ_ALB_INDS(K-1)
              ! accumulate 8 radiation diags.
              DO IT = 1, NTYPE
                 CALL INC_AJ(I,J,IT,JK,(S0*CSZ2)*ALB(I,J,K)         &
                      *FTYPE(IT,I,J))
              ENDDO
              CALL INC_AREG(I,J,JR,JK,(S0*CSZ2)*ALB(I,J,K))
           ENDDO
           AIJ(I,J,IJ_SRINCG) = AIJ(I,J,IJ_SRINCG)                  &
                + (SRHR(0,I,J)*CSZ2/(ALB(I,J,1)     &
                +1.D-20))
           AIJ(I,J,IJ_SRNFG) = AIJ(I,J,IJ_SRNFG)                    &
                + (SRHR(0,I,J)*CSZ2)
           AIJ(I,J,IJ_BTMPW) = AIJ(I,J,IJ_BTMPW) + BTMPW(I,J)
           AIJ(I,J,IJ_SRREF) = AIJ(I,J,IJ_SRREF)                    &
                + S0*CSZ2*ALB(I,J,2)
           AIJ(I,J,IJ_SRVIS) = AIJ(I,J,IJ_SRVIS)                    &
                + S0*CSZ2*ALB(I,J,4)
           AIJ(I,J,IJ_TRNFP0) = AIJ(I,J,IJ_TRNFP0) - TNFS(3,I,J)
           AIJ(I,J,IJ_SRNFP0) = AIJ(I,J,IJ_SRNFP0)                  &
                + (SNFS(3,I,J)*CSZ2)
           AIJ(I,J,IJ_RNFP1) = AIJ(I,J,IJ_RNFP1)                    &
                + (SNFS(2,I,J)*CSZ2-TNFS(2,I,J))
           AIJ(I,J,ij_srvdir) = AIJ(I,J,ij_srvdir) + FSRDIR(I,J)    &
                *SRVISSURF(I,J)
           AIJ(I,J,IJ_SRVISSURF) = AIJ(I,J,IJ_SRVISSURF)            &
                + SRVISSURF(I,J)
#ifdef mjo_subdd
           OLR_ACC(I,J) = OLR_ACC(I,J) - TNFS(3,I,J)
#endif
           !**** CRF diags if required
           IF ( moddrf==0 ) THEN
              IF ( cloud_rad_forc>0 ) THEN
                 !    CRF diagnostics
                 AIJ(I,J,IJ_SWCRF) = AIJ(I,J,IJ_SWCRF)              &
                      + (SNFS(3,I,J)-SNFSCRF(I,J))*CSZ2
                 AIJ(I,J,IJ_LWCRF) = AIJ(I,J,IJ_LWCRF)              &
                      - (TNFS(3,I,J)-TNFSCRF(I,J))
              ENDIF
              IF ( cloud_rad_forc==2 ) THEN
                 !    CRF diagnostics without aerosols and Ox
                 AIJ(I,J,IJ_SWCRF2) = AIJ(I,J,IJ_SWCRF2)            &
                      + (SNFS(3,I,J)-SNFSCRF2(I,J))*CSZ2
                 AIJ(I,J,IJ_LWCRF2) = AIJ(I,J,IJ_LWCRF2)            &
                      - (TNFS(3,I,J)-TNFSCRF2(I,J))
              ENDIF

              !**** AERRF diags if required
              IF ( aer_rad_forc>0 ) THEN
                 DO N = 1, 8
                    AIJ(I,J,IJ_SWAERRF+N-1)                         &
                         = AIJ(I,J,IJ_SWAERRF+N-1)                    &
                         + (SNFS(3,I,J)-SNFSAERRF(N,I,J))*CSZ2
                    AIJ(I,J,IJ_LWAERRF+N-1)                         &
                         = AIJ(I,J,IJ_LWAERRF+N-1)                    &
                         - (TNFS(3,I,J)-TNFSAERRF(N,I,J))
                    AIJ(I,J,IJ_SWAERSRF+N-1)                        &
                         = AIJ(I,J,IJ_SWAERSRF+N-1)                   &
                         + (SNFS(1,I,J)-SNFSAERRF(N+8,I,J))*CSZ2
                    AIJ(I,J,IJ_LWAERSRF+N-1)                        &
                         = AIJ(I,J,IJ_LWAERSRF+N-1)                   &
                         - (TNFS(1,I,J)-TNFSAERRF(N+8,I,J))
                 ENDDO
                 AIJ(I,J,IJ_SWAERRFNT) = AIJ(I,J,IJ_SWAERRFNT)      &
                      + (SNFS(3,I,J)-SNFSAERRF(17,I,J))*CSZ2
                 AIJ(I,J,IJ_LWAERRFNT) = AIJ(I,J,IJ_LWAERRFNT)      &
                      - (TNFS(3,I,J)-TNFSAERRF(17,I,J))
                 AIJ(I,J,IJ_SWAERSRFNT) = AIJ(I,J,IJ_SWAERSRFNT)    &
                      + (SNFS(1,I,J)-SNFSAERRF(18,I,J))*CSZ2
                 AIJ(I,J,IJ_LWAERSRFNT) = AIJ(I,J,IJ_LWAERSRFNT)    &
                      - (TNFS(1,I,J)-TNFSAERRF(18,I,J))
              ENDIF

              !***** Clear Sky and All Sky TOA Forcing without aerosol
              AIJ(I,J,IJ_SW_AS_noA) = AIJ(I,J,IJ_SW_AS_noA)         &
                   + (SNFS(3,I,J)-SNFS_AS_noA(I,J))*CSZ2
              AIJ(I,J,IJ_LW_AS_noA) = AIJ(I,J,IJ_LW_AS_noA)         &
                   - (TNFS(3,I,J)-TNFS_AS_noA(I,J))
              AIJ(I,J,IJ_SW_CS_noA) = AIJ(I,J,IJ_SW_CS_noA)         &
                   + (SNFS(3,I,J)-SNFS_CS_noA(I,J))*CSZ2
              AIJ(I,J,IJ_LW_CS_noA) = AIJ(I,J,IJ_LW_CS_noA)         &
                   - (TNFS(3,I,J)-TNFS_CS_noA(I,J))


#ifdef TRACERS_GC
              !**** Generic diagnostics for radiative forcing calculations
              !**** Depending on whether tracers radiative interaction is turned on,
              !**** diagnostic sign changes (for aerosols)
!              rsign_aer = 1.
              rsign_chem = -1.
!               IF ( rad_interact_aer>0 ) rsign_aer = -1.
!               !**** define SNFS/TNFS level (TOA/TROPO) for calculating forcing
!               LFRC = 3      ! TOA
!               IF ( rad_forc_lev>0 ) LFRC = 4
																			
              ! TOA 
		            SAVE_RF(I,J,5,1)    =  rsign_chem*(SNFST_o3ref(2,I,J)-SNFS(3,I,J))*CSZ2 ! SW
		            SAVE_RF(I,J,5,2)    = -rsign_chem*(TNFST_o3ref(2,I,J)-TNFS(3,I,J))      ! LW
														! Topopause
		            SAVE_RF_TP(I,J,5,1) =  rsign_chem*(SNFST_o3ref(1,I,J)-SNFS(4,I,J))*CSZ2 ! SW
		            SAVE_RF_TP(I,J,5,2) = -rsign_chem*(TNFST_o3ref(1,I,J)-TNFS(4,I,J))      ! LW
														! Whole atmosphere
		            SAVE_RF_3D(I,J,:,5,1) =  rsign_chem*(SNFS_3D_pert(I,J,:,5)-SNFS_3D(I,J,:))*CSZ2 ! SW
		            SAVE_RF_3D(I,J,:,5,2) = -rsign_chem*(TNFS_3D_pert(I,J,:,5)-TNFS_3D(I,J,:))      ! LW

#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
              (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_MINERALS) ||\
              (defined TRACERS_AMP) || (defined TRACERS_TOMAS) ||\
              (defined TRACERS_AEROSOLS_SEASALT)
              !**** Generic diagnostics for radiative forcing calculations
              !**** Depending on whether tracers radiative interaction is turned on,
              !**** diagnostic sign changes (for aerosols)
              rsign_aer = 1.
              rsign_chem = -1.
              IF ( rad_interact_aer>0 ) rsign_aer = -1.
              !**** define SNFS/TNFS level (TOA/TROPO) for calculating forcing
              LFRC = 3      ! TOA
              IF ( rad_forc_lev>0 ) LFRC = 4
              ! TROPOPAUSE
#ifdef BC_ALB
              IF ( IJTS_ALB(1)>0 .AND. bc_snow_present(i,j) .AND.   &
                   csz2>0. ) THEN
                 TAIJS(I,J,ijts_sunlit_snow)                        &
                      = TAIJS(I,J,ijts_sunlit_snow) + 1.
                 TAIJS(I,J,IJTS_ALB(1)) = TAIJS(I,J,IJTS_ALB(1))    &
                      + dALBsnBC(I,J)
                 !           + 100.d0*(ALBNBC(I,J)-ALB(I,J,1))
              ENDIF
              IF ( IJTS_ALB(2)>0 ) TAIJS(i,j,IJTS_ALB(2))           &
                   = TAIJS(i,j,IJTS_ALB(2))                         &
                   + (SNFS(3,I,J)-NFSNBC(I,J))*CSZ2
#endif /* BC_ALB */
              !     ..........
              !     accumulation of forcings for tracers for which nraero_rf fields are
              !     defined
              !     ..........
              nsub_ntrix = 0
              DO n = 1, nraero_rf
                 SELECT CASE (TRNAME(NTRIX_RF(n)))
                 CASE ('Clay','ClayIlli','ClayKaol','ClaySmec',     &
                      &'ClayCalc','ClayQuar','ClayFeld','ClayHema', &
                      &'ClayGyps','ClayIlHe','ClayKaHe','ClaySmHe', &
                      &'ClayCaHe','ClayQuHe','ClayFeHe','ClayGyHe')
                    nsub_ntrix(NTRIX_RF(n))                         &
                         = nsub_ntrix(NTRIX_RF(n)) + 1
                    ! shortwave forcing (TOA or TROPO) of Clay sub size classes
                    IF ( IJTS_FCSUB(1,NTRIX_RF(n),nsub_ntrix(       &
                         NTRIX_RF(n)))>0 )                          &
                         TAIJS(i,j,IJTS_FCSUB(1,NTRIX_RF(n),        &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         = TAIJS(i,j,IJTS_FCSUB(1,NTRIX_RF(n),      &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         + rsign_aer*(snfst(2,n,i,j)-snfs(lfrc,i,j))&
                         *csz2
                    ! longwave forcing  (TOA or TROPO) of Clay size sub classes
                    IF ( IJTS_FCSUB(2,NTRIX_RF(n),nsub_ntrix(       &
                         NTRIX_RF(n)))>0 )                          &
                         TAIJS(i,j,IJTS_FCSUB(2,NTRIX_RF(n),        &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         = TAIJS(i,j,IJTS_FCSUB(2,NTRIX_RF(n),      &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         - rsign_aer*(tnfst(2,n,i,j)-tnfs(lfrc,i,j))
                    ! shortwave forcing (TOA or TROPO) clear sky of Clay sub size classes
                    IF ( IJTS_FCSUB(5,NTRIX_RF(n),nsub_ntrix(       &
                         NTRIX_RF(n)))>0 )                          &
                         TAIJS(i,j,IJTS_FCSUB(5,NTRIX_RF(n),        &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         = TAIJS(i,j,IJTS_FCSUB(5,NTRIX_RF(n),      &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         + rsign_aer*(snfst(2,n,i,j)-snfs(lfrc,i,j))&
                         *csz2*(1.D0-cfrac(i,j))
                    ! longwave forcing  (TOA or TROPO) clear sky of Clay sub size classes
                    IF ( IJTS_FCSUB(6,NTRIX_RF(n),nsub_ntrix(       &
                         NTRIX_RF(n)))>0 )                          &
                         TAIJS(i,j,IJTS_FCSUB(6,NTRIX_RF(n),        &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         = TAIJS(i,j,IJTS_FCSUB(6,NTRIX_RF(n),      &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         - rsign_aer*(tnfst(2,n,i,j)-tnfs(lfrc,i,j))&
                         *(1.D0-cfrac(i,j))
                    ! shortwave forcing at surface (if required) of Clay sub size classes
                    IF ( IJTS_FCSUB(3,NTRIX_RF(n),nsub_ntrix(       &
                         NTRIX_RF(n)))>0 )                          &
                         TAIJS(i,j,IJTS_FCSUB(3,NTRIX_RF(n),        &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         = TAIJS(i,j,IJTS_FCSUB(3,NTRIX_RF(n),      &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         + rsign_aer*(snfst(1,n,i,j)-snfs(1,i,j))   &
                         *csz2
                    ! longwave forcing at surface (if required) of Clay sub size classes
                    IF ( IJTS_FCSUB(4,NTRIX_RF(n),nsub_ntrix(       &
                         NTRIX_RF(n)))>0 )                          &
                         TAIJS(i,j,IJTS_FCSUB(4,NTRIX_RF(n),        &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         = TAIJS(i,j,IJTS_FCSUB(4,NTRIX_RF(n),      &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         - rsign_aer*(tnfst(1,n,i,j)-tnfs(1,i,j))
                    ! shortwave forcing at surface clear sky (if required) of Clay sub size classes
                    IF ( IJTS_FCSUB(7,NTRIX_RF(n),nsub_ntrix(       &
                         NTRIX_RF(n)))>0 )                          &
                         TAIJS(i,j,IJTS_FCSUB(7,NTRIX_RF(n),        &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         = TAIJS(i,j,IJTS_FCSUB(7,NTRIX_RF(n),      &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         + rsign_aer*(snfst(1,n,i,j)-snfs(1,i,j))   &
                         *csz2*(1.D0-cfrac(i,j))
                    ! longwave forcing at surface clear sky (if required) of Clay sub size classes
                    IF ( IJTS_FCSUB(8,NTRIX_RF(n),nsub_ntrix(       &
                         NTRIX_RF(n)))>0 )                          &
                         TAIJS(i,j,IJTS_FCSUB(8,NTRIX_RF(n),        &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         = TAIJS(i,j,IJTS_FCSUB(8,NTRIX_RF(n),      &
                         nsub_ntrix(NTRIX_RF(n))))                  &
                         - rsign_aer*(tnfst(1,n,i,j)-tnfs(1,i,j))   &
                         *(1.D0-cfrac(i,j))
                 CASE DEFAULT
                    SELECT CASE (TRNAME(NTRIX_RF(n)))
                    CASE ('seasalt2')
                       CYCLE
                    ENDSELECT
                    ! shortwave forcing (TOA or TROPO)
                    IF ( IJTS_FC(1,NTRIX_RF(n))>0 )                 &
                         TAIJS(i,j,IJTS_FC(1,NTRIX_RF(n)))          &
                         = TAIJS(i,j,IJTS_FC(1,NTRIX_RF(n)))        &
                         + rsign_aer*(SNFST(2,N,I,J)-SNFS(LFRC,I,J))&
                         *CSZ2
                    ! longwave forcing  (TOA or TROPO)
                    IF ( IJTS_FC(2,NTRIX_RF(n))>0 )                 &
                         TAIJS(i,j,IJTS_FC(2,NTRIX_RF(n)))          &
                         = TAIJS(i,j,IJTS_FC(2,NTRIX_RF(n)))        &
                         - rsign_aer*(TNFST(2,N,I,J)-TNFS(LFRC,I,J))
                    ! shortwave forcing (TOA or TROPO) clear sky
                    IF ( IJTS_FC(5,NTRIX_RF(n))>0 )                 &
                         TAIJS(i,j,IJTS_FC(5,NTRIX_RF(n)))          &
                         = TAIJS(i,j,IJTS_FC(5,NTRIX_RF(n)))        &
                         + rsign_aer*(SNFST(2,N,I,J)-SNFS(LFRC,I,J))&
                         *CSZ2*(1.D0-CFRAC(I,J))
                    ! longwave forcing  (TOA or TROPO) clear sky
                    IF ( IJTS_FC(6,NTRIX_RF(n))>0 )                 &
                         TAIJS(i,j,IJTS_FC(6,NTRIX_RF(n)))          &
                         = TAIJS(i,j,IJTS_FC(6,NTRIX_RF(n)))        &
                         - rsign_aer*(TNFST(2,N,I,J)-TNFS(LFRC,I,J))&
                         *(1.D0-CFRAC(I,J))
                    ! shortwave forcing at surface (if required)
                    IF ( IJTS_FC(3,NTRIX_RF(n))>0 )                 &
                         TAIJS(i,j,IJTS_FC(3,NTRIX_RF(n)))          &
                         = TAIJS(i,j,IJTS_FC(3,NTRIX_RF(n)))        &
                         + rsign_aer*(SNFST(1,N,I,J)-SNFS(1,I,J))   &
                         *CSZ2
                    ! longwave forcing at surface (if required)
                    IF ( IJTS_FC(4,NTRIX_RF(n))>0 )                 &
                         TAIJS(i,j,IJTS_FC(4,NTRIX_RF(n)))          &
                         = TAIJS(i,j,IJTS_FC(4,NTRIX_RF(n)))        &
                         - rsign_aer*(TNFST(1,N,I,J)-TNFS(1,I,J))
                    ! shortwave forcing at surface clear sky (if required)
                    IF ( IJTS_FC(7,NTRIX_RF(n))>0 )                 &
                         TAIJS(i,j,IJTS_FC(7,NTRIX_RF(n)))          &
                         = TAIJS(i,j,IJTS_FC(7,NTRIX_RF(n)))        &
                         + rsign_aer*(SNFST(1,N,I,J)-SNFS(1,I,J))   &
                         *CSZ2*(1.D0-CFRAC(I,J))
                    ! longwave forcing at surface clear sky (if required)
                    IF ( IJTS_FC(8,NTRIX_RF(n))>0 )                 &
                         TAIJS(i,j,IJTS_FC(8,NTRIX_RF(n)))          &
                         = TAIJS(i,j,IJTS_FC(8,NTRIX_RF(n)))        &
                         - rsign_aer*(TNFST(1,N,I,J)-TNFS(1,I,J))   &
                         *(1.D0-CFRAC(I,J))
                 ENDSELECT
              ENDDO
              ! n=1,nraero_rf

              ! ..........
              ! accumulation of forcings for special case ozone (nraero_rf fields
              ! not defined) Warning :  indicies used differently, since we don't
              ! need CS or Surface, but are doing both TOA and Ltropo : 
              ! ..........
              IF ( n_Ox>0 ) THEN
                 ! ------ main Ox tracer -------
                 ! shortwave forcing at tropopause
                 IF ( IJTS_FC(1,n_Ox)>0 ) TAIJS(i,j,IJTS_FC(1,n_Ox))&
                      = TAIJS(i,j,IJTS_FC(1,n_Ox))                  &
                      + rsign_chem*(SNFST_o3ref(1,I,J)-SNFS(4,I,J)&
                      )*CSZ2
                 ! longwave forcing at tropopause
                 IF ( IJTS_FC(2,n_Ox)>0 ) TAIJS(i,j,IJTS_FC(2,n_Ox))&
                      = TAIJS(i,j,IJTS_FC(2,n_Ox))                  &
                      - rsign_chem*(TNFST_o3ref(1,I,J)-TNFS(4,I,J)&
                      )
                 ! shortwave forcing at TOA
                 IF ( IJTS_FC(3,n_Ox)>0 ) TAIJS(i,j,IJTS_FC(3,n_Ox))&
                      = TAIJS(i,j,IJTS_FC(3,n_Ox))                  &
                      + rsign_chem*(SNFST_o3ref(2,I,J)-SNFS(3,I,J)&
                      )*CSZ2
                 ! longwave forcing at TOA
                 IF ( IJTS_FC(4,n_Ox)>0 ) TAIJS(i,j,IJTS_FC(4,n_Ox))&
                      = TAIJS(i,j,IJTS_FC(4,n_Ox))                  &
                      - rsign_chem*(TNFST_o3ref(2,I,J)-TNFS(3,I,J)&
                      )
              ENDIF
#ifdef AUXILIARY_OX_RADF
              ! shortwave forcing at tropopause

#ifdef AUX_OX_RADF_TROP
              IF ( IJTS_AUXFC(1)>0 ) TAIJS(i,j,IJTS_AUXFC(1))       &
                   = TAIJS(i,j,IJTS_AUXFC(1))                       &
                   + rsign_chem*(SNFST_o3ref(5,I,J)                 &
                   -SNFST_o3ref(3,I,J))*CSZ2
#else
              IF ( IJTS_AUXFC(1)>0 ) TAIJS(i,j,IJTS_AUXFC(1))       &
                   = TAIJS(i,j,IJTS_AUXFC(1))                       &
                   + rsign_chem*(SNFST_o3ref(1,I,J)                 &
                   -SNFST_o3ref(3,I,J))*CSZ2
#endif
              ! longwave forcing at tropopause
#ifdef AUX_OX_RADF_TROP
              IF ( IJTS_AUXFC(2)>0 ) TAIJS(i,j,IJTS_AUXFC(2))       &
                   = TAIJS(i,j,IJTS_AUXFC(2))                       &
                   - rsign_chem*(TNFST_o3ref(5,I,J)                 &
                   -TNFST_o3ref(3,I,J))
#else
              IF ( IJTS_AUXFC(2)>0 ) TAIJS(i,j,IJTS_AUXFC(2))       &
                   = TAIJS(i,j,IJTS_AUXFC(2))                       &
                   - rsign_chem*(TNFST_o3ref(1,I,J)                 &
                   -TNFST_o3ref(3,I,J))
#endif
              ! shortwave forcing at TOA
              IF ( IJTS_AUXFC(3)>0 ) TAIJS(i,j,IJTS_AUXFC(3))       &
                   = TAIJS(i,j,IJTS_AUXFC(3))                       &
                   + rsign_chem*(SNFST_o3ref(2,I,J)                 &
                   -SNFST_o3ref(4,I,J))*CSZ2
              ! longwave forcing at TOA
              IF ( IJTS_AUXFC(4)>0 ) TAIJS(i,j,IJTS_AUXFC(4))       &
                   = TAIJS(i,j,IJTS_AUXFC(4))                       &
                   - rsign_chem*(TNFST_o3ref(2,I,J)                 &
                   -TNFST_o3ref(4,I,J))
#endif /* AUXILIARY_OX_RADF */
#if (defined SHINDELL_STRAT_EXTRA) &(defined ACCMIP_LIKE_DIAGS)
              ! ------ diag stratOx tracer -------
              ! note for now for this diag, there is a failsafe that stops model
              ! if clim_interact_chem .le. 0 when the below would be wrong : 
              ! shortwave forcing at tropopause
              IF ( IJTS_FC(1,n_stratOx)>0 )                         &
                   TAIJS(i,j,IJTS_FC(1,n_stratOx))                  &
                   = TAIJS(i,j,IJTS_FC(1,n_stratOx))                &
                   + rsign_chem*(SNFST_o3ref(1,I,J)                 &
                   -SNFST_stratOx(1,I,J))*CSZ2
              ! longwave forcing at tropopause
              IF ( IJTS_FC(2,n_stratOx)>0 )                         &
                   TAIJS(i,j,IJTS_FC(2,n_stratOx))                  &
                   = TAIJS(i,j,IJTS_FC(2,n_stratOx))                &
                   - rsign_chem*(TNFST_o3ref(1,I,J)                 &
                   -TNFST_stratOx(1,I,J))
              ! shortwave forcing at TOA
              IF ( IJTS_FC(3,n_stratOx)>0 )                         &
                   TAIJS(i,j,IJTS_FC(3,n_stratOx))                  &
                   = TAIJS(i,j,IJTS_FC(3,n_stratOx))                &
                   + rsign_chem*(SNFST_o3ref(2,I,J)                 &
                   -SNFST_stratOx(2,I,J))*CSZ2
              ! longwave forcing at TOA
              IF ( IJTS_FC(4,n_stratOx)>0 )                         &
                   TAIJS(i,j,IJTS_FC(4,n_stratOx))                  &
                   = TAIJS(i,j,IJTS_FC(4,n_stratOx))                &
                   - rsign_chem*(TNFST_o3ref(2,I,J)                 &
                   -TNFST_stratOx(2,I,J))
#endif /* SHINDELL_STRAT_EXTRA &ACCMIP_LIKE_DIAGS*/
#endif /* any of various tracer groups defined */

#ifdef TRACERS_GC
              !============================================
              ! Methane, N2O, CFC11 and CFC12  
              !============================================
														DO nf = 1, 4
															
                 ! TOA
                 SAVE_RF(I,J,nf,1) = (SNFS(3,I,J)-SNFS_ghg(nf,I,J))*CSZ2 ! SW
                 SAVE_RF(I,J,nf,2) = (TNFS_ghg(nf,I,J)-TNFS(3,I,J))      ! LW

                 IF ( IJ_FCGHG(1,nf)>0 )                                  &
																	       AIJ(i,j,IJ_FCGHG(1,nf)) = AIJ(i,j,IJ_FCGHG(1,nf)) &
                      + (SNFS(3,I,J)-SNFS_ghg(nf,I,J))*CSZ2
																	IF ( IJ_FCGHG(2,nf)>0 )                                  &
																	       AIJ(i,j,IJ_FCGHG(2,nf)) = AIJ(i,j,IJ_FCGHG(2,nf)) &
                      + (TNFS_ghg(nf,I,J)-TNFS(3,I,J))

                 ! Tropopause
                 SAVE_RF_TP(I,J,nf,1) = (SNFS(4,I,J)-SNFS_ghg_tp(nf,I,J))*CSZ2         ! SW
                 SAVE_RF_TP(I,J,nf,2) = (TNFS_ghg_tp(nf,I,J)-TNFS(4,I,J))              ! LW

                 ! Whole Atmosphere
                 SAVE_RF_3D(I,J,:,nf,1) = (SNFS_3D(I,J,:)-SNFS_3D_pert(I,J,:,nf))*CSZ2 ! SW
                 SAVE_RF_3D(I,J,:,nf,2) = (TNFS_3D_pert(I,J,:,nf)-TNFS_3D(I,J,:))      ! LW
																	
              ENDDO

#endif

#ifdef ACCMIP_LIKE_DIAGS
#ifndef SKIP_ACCMIP_GHG_RADF_DIAGS
              DO nf = 1, 4
                 ! CH4, N2O, CFC11, and CFC12 : 
                 ! shortwave GHG forcing at TOA
                 IF ( IJ_FCGHG(1,nf)>0 ) AIJ(i,j,IJ_FCGHG(1,nf))    &
                      = AIJ(i,j,IJ_FCGHG(1,nf))                     &
                      + (SNFS(3,I,J)-SNFS_ghg(nf,I,J))*CSZ2
                 ! longwave GHG forcing at TOA
                 IF ( IJ_FCGHG(2,nf)>0 ) AIJ(i,j,IJ_FCGHG(2,nf))    &
                      = AIJ(i,j,IJ_FCGHG(2,nf))                     &
                      + (TNFS_ghg(nf,I,J)-TNFS(3,I,J))
              ENDDO
#endif /* NOT DEFINED SKIP_ACCMIP_GHG_RADF_DIAGS */
#endif /* ACCMIP_LIKE_DIAGS */

#ifdef CACHED_SUBDD
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
              (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_MINERALS) ||\
              (defined TRACERS_AMP) || (defined TRACERS_TOMAS) ||\
              (defined TRACERS_AEROSOLS_SEASALT)
              IF ( nraero_rf>0 ) THEN
                 swfrc(i,j,1 : nraero_rf)                             &
                      = rsign_aer*(SNFST(2,1 : nraero_rf,I,J)           &
                      -SNFS(LFRC,I,J))*CSZ2
                 lwfrc(i,j,1 : nraero_rf)                             &
                      = -rsign_aer*(TNFST(2,1 : nraero_rf,I,J)          &
                      -TNFS(LFRC,I,J))
              ENDIF
#endif /* any of various tracer groups defined */
#endif  /* CACHED_SUBDD */

           ENDIF
        ENDDO
     ENDDO

#ifdef mjo_subdd
     OLR_cnt = OLR_cnt + 1.
#endif

     DO J = J_0, J_1
        DO I = I_0, I_1
           DO L = 1, LM
              AIJL(i,j,l,IJL_RC) = AIJL(i,j,l,IJL_RC)               &
                   + (SRHR(L,I,J)*COSZ2(I,J)        &
                   +TRHR(L,I,J))
#ifdef mjo_subdd
              SWHR(I,J,L) = SWHR(I,J,L) + SRHR(L,I,J)*COSZ2(I,J)    &
                   *bysha*BYMA(L,I,J)
              LWHR(I,J,L) = LWHR(I,J,L) + TRHR(L,I,J)               &
                   *bysha*BYMA(L,I,J)
#endif
           ENDDO
        ENDDO
     ENDDO
#ifdef mjo_subdd
     SWHR_cnt = SWHR_cnt + 1
     LWHR_cnt = LWHR_cnt + 1
#endif

#ifdef CACHED_SUBDD
     DO k = 1, subdd_ngroups
        subdd => SUBDD_GROUPS(k)
        subdd%NACC(subdd%SUBDD_PERIOD,sched_rad)                    &
             = subdd%NACC(subdd%SUBDD_PERIOD,sched_rad) + 1
     ENDDO
     !****
     !**** Collect some high-frequency outputs
     !****
     CALL FIND_GROUPS('rijh',grpids,ngroups)
     DO igrp = 1, ngroups
        subdd => SUBDD_GROUPS(grpids(igrp))
        DO k = 1, subdd%NDIAGS
           SELECT CASE (subdd%NAME(k))
           CASE ('olrrad')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = tnfs(3,i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('olrcs')
              IF ( cloud_rad_forc<=0. ) CALL STOP_MODEL(            &
                   &'diagnostic olrcs needs cloud_rad_forc>0',255)
              CALL INC_SUBDD(subdd,k,TNFSCRF)
           CASE ('lwds')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = TRHR(0,i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('lwdscs')
              IF ( cloud_rad_forc<=0. ) CALL STOP_MODEL(            &
                   &'diagnostic lwdscs needs cloud_rad_forc>0',255)
              CALL INC_SUBDD(subdd,k,lwdncs)
           CASE ('lwus')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = TRHR(0,i,j) + tnfs(1,i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('swus')
              CALL INC_SUBDD(subdd,k,SWUS)
           CASE ('swds')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SRDN(i,j)*cosz2(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('swdf')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = FSRDIF(i,j) + DIFNIR(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('swtoa')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = snfs(3,i,j)*cosz2(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
              !Net solar flux at surface : 
           CASE ('swns')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = snfs(1,i,j)*cosz2(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
              !Net Longwave flux at surface : 
           CASE ('lwns')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = tnfs(1,i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('totcld')
              CALL INC_SUBDD(subdd,k,cfrac)
           CASE ('totcld_diag')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    CALL GET_CLD_OVERLAP(lm,CLDSS( : ,i,j),           &
                         CLDMCL=CLDMC( : ,i,j),CLDTOT=sddarr(i,j))
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('cldss_2d')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    CALL GET_CLD_OVERLAP(lm,CLDSS( : ,i,j),           &
                         CLDSS=sddarr(i,j))
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('cldmc_2d')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    CALL GET_CLD_OVERLAP(lm,CLDSS( : ,i,j),           &
                         CLDMCL=CLDMC( : ,i,j),CLDMC=sddarr(i,j))
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('wtrcld')
              CALL INC_SUBDD(subdd,k,WTRCLD)
           CASE ('icecld')
              CALL INC_SUBDD(subdd,k,ICECLD)
           CASE ('cod')
              CALL INC_SUBDD(subdd,k,TAUSUMW)
           CASE ('cid')
              CALL INC_SUBDD(subdd,k,TAUSUMI)
           CASE ('ctp')
              CALL INC_SUBDD(subdd,k,CTP)
           CASE ('ctt')
              CALL INC_SUBDD(subdd,k,CTT)
#ifdef CFMIP3_SUBDD
           CASE ('rtmt')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = (snfs(3,i,j)*cosz2(i,j))          &
                         - tnfs(2,i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('swut')
              CALL INC_SUBDD(subdd,k,swut)
           CASE ('swutcs')
              CALL INC_SUBDD(subdd,k,swutcs)
           CASE ('clwvi')
              CALL INC_SUBDD(subdd,k,cfmip_twp)
           CASE ('swdcls')
              CALL INC_SUBDD(subdd,k,swdcls)
           CASE ('swucls')
              CALL INC_SUBDD(subdd,k,swucls)
           CASE ('swdt')
              CALL INC_SUBDD(subdd,k,swdt)
#endif
           ENDSELECT

        ENDDO
     ENDDO

#ifdef GCAP
     CALL FIND_GROUPS('aijlh',grpids,ngroups)
     DO igrp = 1, ngroups
        subdd => SUBDD_GROUPS(grpids(igrp))
        DO k = 1, subdd%NDIAGS
           SELECT CASE (subdd%NAME(k))
           CASE ('OPTDEPTH')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    DO l = 1, lmaxsubdd
                       ! Weight mean cloud optical thickness by relative 2-D area fractions
                       sddarr3d(i,j,l)                              &
                            = (CLDSS(l,i,j)*TAUSS(l,i,j)+CLDMC(l,i,j) &
                            *TAUMC(l,i,j))                            &
                            /(CLDSS(l,i,j)+CLDMC(l,i,j)+teeny)
                    ENDDO
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr3d)
           CASE ('CLOUD')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    DO l = 1, lmaxsubdd
                       sddarr3d(i,j,l)                              &
                            = MIN(1.0,CLDSS3D(l,i,j)+CLDMC(l,i,j))
                    ENDDO
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr3d)
           CASE ('TAUCLI')
              CALL INC_SUBDD(subdd,k,taui3d)
           CASE ('TAUCLW')
              CALL INC_SUBDD(subdd,k,tauw3d)
           ENDSELECT
        ENDDO
     ENDDO

     CALL FIND_GROUPS('aijh',grpids,ngroups)
     DO igrp = 1, ngroups
        subdd => SUBDD_GROUPS(grpids(igrp))
        DO k = 1, subdd%NDIAGS
           SELECT CASE (subdd%NAME(k))
           CASE ('PARDF')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = 0.82*SRVISSURF(i,j)               &
                         *(1D0-FSRDIR(i,j))*COSZ1(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('PARDR')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = 0.82*SRVISSURF(i,j)*(FSRDIR(i,j)) &
                         *COSZ1(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('ALBEDO')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    ! Saves non-zero albedos for nighttime
                    IF ( SRDN(i,j)>0 ) SAVE_ALB(i,j)                &
                         = 1D0 - ALB(i,j,1)
                    sddarr(i,j) = SAVE_ALB(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('CLDTOT')
              ! totcld in standard model
              CALL INC_SUBDD(subdd,k,cfrac)
           CASE ('SWGDN')
              ! swds in standard model
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SRDN(i,j)*cosz2(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('TO3')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_TO3(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)										
           ENDSELECT
        ENDDO
     ENDDO
#endif

#ifdef TRACERS_GC
     CALL FIND_GROUPS('aijh',grpids,ngroups)
     DO igrp = 1, ngroups
        subdd => SUBDD_GROUPS(grpids(igrp))
        DO k = 1, subdd%NDIAGS
           SELECT CASE (subdd%NAME(k))
           
											!==================================================
           ! Top of the Atmosphere
           !==================================================
											CASE ('SW_CH4')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF(i,j,1,1)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('LW_CH4')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF(i,j,1,2)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('SW_N2O')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF(i,j,2,1)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('LW_N2O')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF(i,j,2,2)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('SW_CFC11')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF(i,j,3,1)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('LW_CFC11')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF(i,j,3,2)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('SW_CFC12')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF(i,j,4,1)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('LW_CFC12')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF(i,j,4,2)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)													
           CASE ('SW_O3')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF(i,j,5,1)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('LW_O3')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF(i,j,5,2)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)							

           !==================================================
           ! Radiative Forcing @ Tropopause
           !==================================================
           CASE ('SW_CH4_TP')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF_TP(i,j,1,1)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('LW_CH4_TP')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF_TP(i,j,1,2)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('SW_N2O_TP')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF_TP(i,j,2,1)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('LW_N2O_TP')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF_TP(i,j,2,2)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('SW_CFC11_TP')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF_TP(i,j,3,1)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('LW_CFC11_TP')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF_TP(i,j,3,2)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('SW_CFC12_TP')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF_TP(i,j,4,1)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('LW_CFC12_TP')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF_TP(i,j,4,2)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)													
           CASE ('SW_O3_TP')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF_TP(i,j,5,1)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)														
           CASE ('LW_O3_TP')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_RF_TP(i,j,5,2)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)							

           ENDSELECT
        ENDDO
     ENDDO

     CALL FIND_GROUPS('rijleh',grpids,ngroups)
     DO igrp = 1, ngroups
        subdd => SUBDD_GROUPS(grpids(igrp))
        DO k = 1, subdd%NDIAGS
           SELECT CASE (subdd%NAME(k))
           
											!==================================================
           ! Whole Atmosphere Radiative Forcing
           !==================================================
											CASE ('SW_FLUX')
											   SDDARRFLX = 0.
											   DO l =   1, LM+LM_REQ+1
              DO j = j_0, j_1
              DO i = i_0, IMAXJ(j)
                 SDDARRFLX(i,j,l) = SAVE_RF_3D(i,j,l,0,1)
														ENDDO
              ENDDO
														ENDDO
              CALL INC_SUBDD(subdd,k,SDDARRFLX)
           CASE ('LW_FLUX')
								      SDDARRFLX = 0.
								      DO l =   1, LM+LM_REQ+1
              DO j = j_0, j_1
              DO i = i_0, IMAXJ(j)
                 SDDARRFLX(i,j,l) = SAVE_RF_3D(i,j,l,0,2)
											   ENDDO
              ENDDO
											   ENDDO
              CALL INC_SUBDD(subdd,k,SDDARRFLX)
											CASE ('SW_CH4_3D')
											   SDDARRFLX = 0.
											   DO l =   1, LM+LM_REQ+1
			           DO j = j_0, j_1
			           DO i = i_0, IMAXJ(j)
			              SDDARRFLX(i,j,l) = SAVE_RF_3D(i,j,l,1,1)
														ENDDO
			           ENDDO
														ENDDO
			           CALL INC_SUBDD(subdd,k,SDDARRFLX)
			        CASE ('LW_CH4_3D')
											   SDDARRFLX = 0.
											   DO l =   1, LM+LM_REQ+1
			           DO j = j_0, j_1
			           DO i = i_0, IMAXJ(j)
			              SDDARRFLX(i,j,l) = SAVE_RF_3D(i,j,l,1,2)
											   ENDDO
			           ENDDO
											   ENDDO
			           CALL INC_SUBDD(subdd,k,SDDARRFLX)
           CASE ('SW_N2O_3D')
              SDDARRFLX = 0.
              DO l =   1, LM+LM_REQ+1
              DO j = j_0, j_1
              DO i = i_0, IMAXJ(j)
                 SDDARRFLX(i,j,l) = SAVE_RF_3D(i,j,l,2,1)
           			ENDDO
              ENDDO
           			ENDDO
              CALL INC_SUBDD(subdd,k,SDDARRFLX)
           CASE ('LW_N2O_3D')
              SDDARRFLX = 0.
              DO l =   1, LM+LM_REQ+1
              DO j = j_0, j_1
              DO i = i_0, IMAXJ(j)
                 SDDARRFLX(i,j,l) = SAVE_RF_3D(i,j,l,2,2)
              ENDDO
              ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,SDDARRFLX)														
			        CASE ('SW_CFC11_3D')
			           SDDARRFLX = 0.
			           DO l =   1, LM+LM_REQ+1
			           DO j = j_0, j_1
			           DO i = i_0, IMAXJ(j)
			              SDDARRFLX(i,j,l) = SAVE_RF_3D(i,j,l,3,1)
			        			ENDDO
			           ENDDO
			        			ENDDO
			           CALL INC_SUBDD(subdd,k,SDDARRFLX)
			        CASE ('LW_CFC11_3D')
			           SDDARRFLX = 0.
			           DO l =   1, LM+LM_REQ+1
			           DO j = j_0, j_1
			           DO i = i_0, IMAXJ(j)
			              SDDARRFLX(i,j,l) = SAVE_RF_3D(i,j,l,3,2)
			           ENDDO
			           ENDDO
			           ENDDO
			           CALL INC_SUBDD(subdd,k,SDDARRFLX)							
						     CASE ('SW_CFC12_3D')
						        SDDARRFLX = 0.
						        DO l =   1, LM+LM_REQ+1
						        DO j = j_0, j_1
						        DO i = i_0, IMAXJ(j)
						           SDDARRFLX(i,j,l) = SAVE_RF_3D(i,j,l,4,1)
						     			ENDDO
						        ENDDO
						     			ENDDO
						        CALL INC_SUBDD(subdd,k,SDDARRFLX)
						     CASE ('LW_CFC12_3D')
						        SDDARRFLX = 0.
						        DO l =   1, LM+LM_REQ+1
						        DO j = j_0, j_1
						        DO i = i_0, IMAXJ(j)
						           SDDARRFLX(i,j,l) = SAVE_RF_3D(i,j,l,4,2)
						        ENDDO
						        ENDDO
						        ENDDO
						        CALL INC_SUBDD(subdd,k,SDDARRFLX)							
           CASE ('SW_O3_3D')
              SDDARRFLX = 0.
              DO l =   1, LM+LM_REQ+1
              DO j = j_0, j_1
              DO i = i_0, IMAXJ(j)
                 SDDARRFLX(i,j,l) = SAVE_RF_3D(i,j,l,5,1)
           			ENDDO
              ENDDO
           			ENDDO
              CALL INC_SUBDD(subdd,k,SDDARRFLX)
           CASE ('LW_O3_3D')
              SDDARRFLX = 0.
              DO l =   1, LM+LM_REQ+1
              DO j = j_0, j_1
              DO i = i_0, IMAXJ(j)
                 SDDARRFLX(i,j,l) = SAVE_RF_3D(i,j,l,5,2)
              ENDDO
              ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,SDDARRFLX)
           ENDSELECT
        ENDDO
     ENDDO
					
#endif
					
#ifdef GCAP
					
     CALL FIND_GROUPS('aijh',grpids,ngroups)
     DO igrp = 1, ngroups
        subdd => SUBDD_GROUPS(grpids(igrp))
        DO k = 1, subdd%NDIAGS
           SELECT CASE (subdd%NAME(k))
           CASE ('PARDF')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = 0.82*SRVISSURF(i,j)               &
                         *(1D0-FSRDIR(i,j))*COSZ1(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('PARDR')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = 0.82*SRVISSURF(i,j)*(FSRDIR(i,j)) &
                         *COSZ1(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('ALBEDO')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    ! Saves non-zero albedos for nighttime
                    IF ( SRDN(i,j)>0 ) SAVE_ALB(i,j)                &
                         = 1D0 - ALB(i,j,1)
                    sddarr(i,j) = SAVE_ALB(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('CLDTOT')
              ! totcld in standard model
              CALL INC_SUBDD(subdd,k,cfrac)
           CASE ('SWGDN')
              ! swds in standard model
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SRDN(i,j)*cosz2(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)
           CASE ('TO3')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    sddarr(i,j) = SAVE_TO3(i,j)
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr)										
           ENDSELECT
        ENDDO
     ENDDO

#endif

     CALL FIND_GROUPS('rijlh',grpids,ngroups)
     DO igrp = 1, ngroups
        subdd => SUBDD_GROUPS(grpids(igrp))
        DO k = 1, subdd%NDIAGS
           SELECT CASE (subdd%NAME(k))
           CASE ('MRCO2rad')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    DO l = 1, lmaxsubdd
                       sddarr3d(i,j,l) = CO2out(l,i,j)
                    ENDDO
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr3d)
           CASE ('wtrtau')
              CALL INC_SUBDD(subdd,k,wtrtau)
           CASE ('icetau')
              CALL INC_SUBDD(subdd,k,icetau)
           ENDSELECT
        ENDDO
     ENDDO

#ifdef SCM
     CALL FIND_GROUPS('rijlh',grpids,ngroups)
     DO igrp = 1, ngroups
        subdd => SUBDD_GROUPS(grpids(igrp))
        DO k = 1, subdd%NDIAGS
           SELECT CASE (subdd%NAME(k))
           CASE ('dth_sw')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    DO l = 1, lmaxsubdd
                       sddarr3d(i,j,l) = SRHR(L,I,J)                &
                            *bysha*BYMA(L,I,J)*COSZ2(I,J)/PK(L,I,J)
                    ENDDO
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr3d)
           CASE ('dth_lw')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    DO l = 1, lmaxsubdd
                       sddarr3d(i,j,l) = TRHR(L,I,J)                &
                            *bysha*BYMA(L,I,J)/PK(L,I,J)
                    ENDDO
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr3d)
           CASE ('dth_rad')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    DO l = 1, lmaxsubdd
                       sddarr3d(i,j,l)                              &
                            = (SRHR(L,I,J)*COSZ2(I,J)+TRHR(L,I,J))    &
                            *bysha*BYMA(L,I,J)/PK(L,I,J)
                    ENDDO
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr3d)
           CASE ('lwdp')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    DO l = 1, lmaxsubdd
                       sddarr3d(i,j,l) = TRDFLB_prof(i,j,l)
                    ENDDO
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr3d)
           CASE ('lwup')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    DO l = 1, lmaxsubdd
                       sddarr3d(i,j,l) = TRUFLB_prof(i,j,l)
                    ENDDO
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr3d)
           CASE ('swdp')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    DO l = 1, lmaxsubdd
                       sddarr3d(i,j,l) = SRDFLB_prof(i,j,l)         &
                            *COSZ2(I,J)
                    ENDDO
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr3d)
           CASE ('swup')
              DO j = j_0, j_1
                 DO i = i_0, IMAXJ(j)
                    DO l = 1, lmaxsubdd
                       sddarr3d(i,j,l) = SRUFLB_prof(i,j,l)         &
                            *COSZ2(I,J)
                    ENDDO
                 ENDDO
              ENDDO
              CALL INC_SUBDD(subdd,k,sddarr3d)
           ENDSELECT
        ENDDO
     ENDDO
#endif
#ifdef CFMIP3_SUBDD
     CALL FIND_GROUPS('rijlh',grpids,ngroups)
     DO igrp = 1, ngroups
        subdd => SUBDD_GROUPS(grpids(igrp))
        DO k = 1, subdd%NDIAGS
           SELECT CASE (subdd%NAME(k))
           CASE ('cf')
              CALL INC_SUBDD(subdd,k,cfmip_cf)
           CASE ('qcirad')
              CALL INC_SUBDD(subdd,k,cfmip_qci)
           CASE ('qclrad')
              CALL INC_SUBDD(subdd,k,cfmip_qcl)
           ENDSELECT
        ENDDO
     ENDDO
#endif
#ifdef TRACERS_ON

     ! aod
     DO g = 1, SIZE(sgroups)
        CALL FIND_GROUPS(sgroups(g),grpids,ngroups)
        DO igrp = 1, ngroups
           subdd => SUBDD_GROUPS(grpids(igrp))
           DO k = 1, subdd%NDIAGS
              DO s = 1, SIZE(ssky)
                 IF ( TRIM(ssky(s))=='dry' .AND. save_dry_aod==0 )  &
                      CYCLE
                 DO a = 1, SIZE(sabs)
                    SELECT CASE (TRIM(ssky(s))//TRIM(sabs(a)))
                    CASE ('as')
                       sddarr4d = tau_as
                    CASE ('cs')
                       sddarr4d = tau_cs
                    CASE ('dry')
                       sddarr4d = tau_dry
                    CASE ('asa')
                       sddarr4d = abstau_as
                    CASE ('csa')
                       sddarr4d = abstau_cs
                    CASE ('drya')
                       sddarr4d = abstau_dry
                    CASE DEFAULT
                       CYCLE
                       ! not implemented, silently ignore
                    ENDSELECT
                    DO n = 1, nraero_aod + 1
                       ! +1 for total
                       IF ( n<=nraero_aod ) THEN
                          spcname = TRIM(TRNAME(NTRIX_AOD(n)))
                       ELSE
                          spcname = ''
                       ENDIF
                       !aod
                       sname = TRIM(spcname)//TRIM(ssky(s))         &
                            //TRIM(sabs(a))//'aod'
                       IF ( TRIM(sgroups(g))=='taijlh' )            &
                            sname = TRIM(sname)//'3d'
                       IF ( TRIM(sname)==TRIM(subdd%NAME(k)) ) THEN
                          ! not select case here
                          IF ( n<=nraero_aod ) THEN
                             sddarr3d = sddarr4d( : , : , : ,n)
                          ELSE
                             sddarr3d = SUM(sddarr4d,DIM=4)
                          ENDIF
                          SELECT CASE (TRIM(sgroups(g)))
                          CASE ('taijh')
                             sddarr = SUM(sddarr3d,DIM=3)
                             CALL INC_SUBDD(subdd,k,sddarr)
                          CASE ('taijlh')
                             CALL INC_SUBDD(subdd,k,sddarr3d)
                          ENDSELECT
                       ENDIF
                       !bext (bcoef) or babs (abcoef)
                       IF ( TRIM(sgroups(g))=='taijlh' ) THEN
                          sname = TRIM(spcname)//TRIM(ssky(s))      &
                               //TRIM(sabs(a))//'bcoef3d'
                          IF ( TRIM(sname)==TRIM(subdd%NAME(k)) )   &
                               THEN                ! not select case here
                             IF ( n<=nraero_aod ) THEN
                                sddarr3d = sddarr4d( : , : , : ,n)
                             ELSE
                                sddarr3d = SUM(sddarr4d,DIM=4)
                             ENDIF
                             DO j = j_0, j_1
                                DO i = i_0, IMAXJ(j)
                                   DO l = 1, lm
                                      TLM(l) = T(i,j,l)*PK(l,i,j)
                                      rho = PMID(l,i,j)              &
                                           *100./(Rgas*TLM(l))
                                      dz = MA(l,i,j)/rho
                                      sddarr3d(i,j,l)                &
                                           = sddarr3d(i,j,l)/dz
                                   ENDDO
                                ENDDO
                             ENDDO
                             CALL INC_SUBDD(subdd,k,sddarr3d)
                          ENDIF
                       ENDIF
                    ENDDO
                    ! n
                 ENDDO
                 ! a
              ENDDO
              ! s
           ENDDO
           ! k
        ENDDO
        ! igrp
     ENDDO
     ! g

     ! rf
     CALL FIND_GROUPS('taijh',grpids,ngroups)
     DO igrp = 1, ngroups
        subdd => SUBDD_GROUPS(grpids(igrp))
        DO k = 1, subdd%NDIAGS
           DO f = 1, SIZE(sfrc)
              SELECT CASE (TRIM(sfrc(f)))
              CASE ('swf')
                 sddarr3drf = swfrc
              CASE ('lwf')
                 sddarr3drf = lwfrc
              CASE DEFAULT
                 CYCLE
                 ! not implemented, silently ignore
              ENDSELECT
              DO n = 1, nraero_rf
                 IF ( diag_fc==2 ) THEN
                    spcname = TRIM(TRNAME(NTRIX_RF(n)))
                 ELSEIF ( diag_fc==1 ) THEN
                    IF ( tracers_amp ) THEN
                       spcname = 'AMP'
                    ELSEIF ( tracers_tomas ) THEN
                       spcname = 'TOMAS'
                    ELSE
                       spcname = 'OMA'
                    ENDIF
                 ENDIF
                 sname = TRIM(sfrc(f))//'_'//TRIM(spcname)
                 IF ( TRIM(sname)==TRIM(subdd%NAME(k)) )            &
                      CALL INC_SUBDD(subdd,k,sddarr3drf( : , : ,n))
                 ! not select case here
              ENDDO
              ! n
           ENDDO
           ! f
        ENDDO
        ! k
     ENDDO
     ! igrp

#endif  /* TRACERS_ON */

#endif  /* CACHED_SUBDD */

     !****
     !**** Update radiative equilibrium temperatures
     !****
     DO J = J_0, J_1
        DO I = I_0, IMAXJ(J)
           DO LR = 1, LM_REQ
              RQT(LR,I,J) = RQT(LR,I,J)                             &
                   + (SRHRS(LR,I,J)*COSZ2(I,J)+TRHRS(LR,I, &
                   J))*NRAD*DTsrc*bysha*BYAML00(lr+lm)
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  !****
  !**** Update other temperatures every physics time step
  !****
  DO J = J_0, J_1
     DO I = I_0, IMAXJ(J)
        DO L = 1, LM
           T(I,J,L) = T(I,J,L)                                      &
                + (SRHR(L,I,J)*COSZ1(I,J)+TRHR(L,I,J))        &
                *DTsrc*bysha*BYMA(l,i,j)/PK(L,I,J)
        ENDDO
        AIJ(I,J,IJ_SRINCP0) = AIJ(I,J,IJ_SRINCP0) + (S0*COSZ1(I,J))
     ENDDO
  ENDDO

  !**** daily diagnostics
  IH = 1 + MODELECLOCK%GETHOUR()
  IHM = IH + (MODELECLOCK%GETDATE()-1)*24
  DO KR = 1, NDIUPT
     I = IJDD(1,KR)
     J = IJDD(2,KR)
     IF ( (J>=J_0) .AND. (J<=J_1) .AND. (I>=I_0) .AND. (I<=I_1) )   &
          THEN
        ADIURN(IDD_ISW,KR,IH) = ADIURN(IDD_ISW,KR,IH)               &
             + S0*COSZ1(I,J)
#ifdef USE_HDIURN
        HDIURN(IDD_ISW,KR,IHM) = HDIURN(IDD_ISW,KR,IHM)             &
             + S0*COSZ1(I,J)
#endif
     ENDIF
  ENDDO

  CALL STOPTIMER('RADIA()')
END SUBROUTINE RADIA

SUBROUTINE RESET_SURF_FLUXES(I,J,ITYPE_OLD,ITYPE_NEW,FTYPE_ORIG,  &
     FTYPE_NOW)
  !@sum set incident solar and upward thermal fluxes appropriately
  !@+   as fractions change to conserve energy, prevent restart problems
  !@auth Gavin Schmidt
  USE RAD_COM, ONLY : FSF, TRSURF
  IMPLICIT NONE
  !@var itype_old, itype_new indices for the old type turning to new type
  INTEGER, INTENT(IN)  ::  i, j, itype_old, itype_new
  !@var ftype_orig, ftype_now original and current fracs of the 'new' type
  REAL*8, INTENT(IN)  ::  ftype_orig, ftype_now
  REAL*8  ::  delf
  ! change in fraction from old to new

  IF ( ((ITYPE_OLD==1 .AND. ITYPE_NEW==2) .OR.                      &
       (ITYPE_OLD==2 .AND. ITYPE_NEW==1)) .AND.                     &
       (FTYPE_NOW<=0. .OR. FTYPE_NOW>1.) ) THEN
     WRITE (6,*) &
          'RESET_SURF_FLUXES : I, J, ITYPE_OLD, ITYPE_NEW, FTYPE_ORIG, FTYPE_NOW = ', &
          I, J, ITYPE_OLD, ITYPE_NEW, FTYPE_ORIG, FTYPE_NOW
     CALL STOP_MODEL('RESET_SURF_FLUXES :  INCORRECT RESET',255)
  ENDIF

  delf = FTYPE_NOW - FTYPE_ORIG
  !**** Constrain fsf_1*ftype_1+fsf_2*ftype_2 to be constant
  FSF(ITYPE_NEW,I,J) = (FSF(ITYPE_NEW,I,J)*FTYPE_ORIG+FSF(ITYPE_OLD,&
       I,J)*DELF)/FTYPE_NOW

  !**** Same for upward thermal
  TRSURF(ITYPE_NEW,I,J) = (TRSURF(ITYPE_NEW,I,J)*FTYPE_ORIG+TRSURF( &
       ITYPE_OLD,I,J)*DELF)/FTYPE_NOW

END SUBROUTINE RESET_SURF_FLUXES

SUBROUTINE GHGHST(iu)
  !@sum  reads history for nghg well-mixed greenhouse gases
  !@auth R. Ruedy

  USE DOMAIN_DECOMP_ATM, ONLY : WRITE_PARALLEL
  USE RADPAR, ONLY : nghg, ghgyr1, ghgyr2, ghgam
  USE RAD_COM, ONLY : ghg_yr
  IMPLICIT NONE
  INTEGER  ::  iu, n, k, nhead = 4, iyr
  CHARACTER*80 title
  CHARACTER(LEN=300)  ::  out_line

  WRITE (out_line,*) ! print header lines and first data line
  CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
  DO n = 1, nhead + 1
     READ (iu,'(a)') title
     WRITE (out_line,'(1x,a80)') title
     CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
  ENDDO
  IF ( title(1 : 2)=='--' ) THEN                ! older format
     READ (iu,'(a)') title
     WRITE (out_line,'(1x,a80)') title
     CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
     nhead = 5
  ENDIF

  !**** find range of table :  ghgyr1 - ghgyr2
  READ (title,*) ghgyr1
  DO
     READ (iu,'(a)',END=20) title
  ENDDO
20 READ (title,*) ghgyr2
  REWIND iu  !   position to data lines
  DO n = 1, nhead
     READ (iu,'(a)')
  ENDDO

  ALLOCATE (ghgam(nghg,ghgyr2-ghgyr1+1))
  DO n = 1, ghgyr2 - ghgyr1 + 1
     READ (iu,*) iyr, (ghgam(k,n),k=1,nghg)
     DO k = 1, nghg
        ! replace -999. by reasonable numbers
        IF ( ghgam(k,n)<0. ) ghgam(k,n) = ghgam(k,n-1)
     ENDDO
     IF ( ghg_yr>0 .AND. ABS(ghg_yr-iyr)<=1 ) THEN
        WRITE (out_line,'(i5,6f10.4)') iyr, (ghgam(k,n),k=1,nghg)
        CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
     ENDIF
  ENDDO
  WRITE (out_line,*) 'read GHG table for years', ghgyr1, ' - ',     &
       ghgyr2
  CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
END SUBROUTINE GHGHST

#if defined(CUBED_SPHERE)
SUBROUTINE READ_QMA(iu,plb)
  !@sum  reads H2O production rates induced by CH4 (Tim Hall)
  !@auth R. Ruedy
  USE DOMAIN_DECOMP_ATM, ONLY : WRITE_PARALLEL
  USE RAD_COM, ONLY : DH2O, jma => JM_DH2O, lat_dh2o
  USE RESOLUTION, ONLY : lm
  USE CONSTANT, ONLY : radian
  USE TIMECONSTANTS_MOD, ONLY : DAYS_PER_YEAR
  IMPLICIT NONE
  INTEGER, PARAMETER  ::  LMA = 24
  INTEGER m, iu, j, l, ll, ldn(lm), lup(lm)
  REAL*8  ::  plb(lm+1)
  REAL*4 pb(0 : LMA+1), h2o(jma,0 : LMA), z(LMA), dz(0 : LMA)
  CHARACTER*100 title
  REAL*4 pdn, pup, dh, fracl
  CHARACTER(LEN=300)  ::  out_line

  !**** read headers/latitudes
  READ (iu,'(a)') title
  WRITE (out_line,'(''0'',a100)') title
  CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
  READ (iu,'(a)') title
  WRITE (out_line,'(1x,a100)') title
  CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
  READ (iu,'(a)') title
  !      write(6,'(1x,a100)') title
  READ (title(10 : 100),*) (lat_dh2o(j),j=1,jma)
  lat_dh2o( : ) = lat_dh2o( : )*radian

  !**** read heights z(km) and data (kg/km^3/year)
  DO m = 1, 12
     READ (iu,'(a)') title
     WRITE (out_line,'(1x,a100)') title
     CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
     DO l = LMA, 1, -1
        READ (iu,'(a)') title
        !          write(6,'(1x,a100)') title
        READ (title,*) z(l), (H2O(j,l),j=1,jma)
     ENDDO
     DO j = 1, jma
        h2o(j,0) = 0.
     ENDDO

     !**** Find edge heights and pressures
     dz(0) = 0.
     dz(1) = z(2) - z(1)
     DO l = 2, LMA - 1
        dz(l) = .5*(z(l+1)-z(l-1))
     ENDDO
     dz(LMA) = z(LMA) - z(LMA-1)

     pb(0) = plb(1)
     DO l = 1, LMA
        Pb(l) = 1000.*10.**(-(z(l)-.5*dz(l))/16.)
     ENDDO
     !**** extend both systems vertically to p=0
     pb(LMA+1) = 0.
     plb(lm+1) = 0.

     !**** Interpolate vertical resolution to model layers
     ldn( : ) = 0
     DO l = 1, lm
        DO WHILE ( pb(ldn(l)+1)>=plb(l) .AND. ldn(l)<LMA )
           ldn(l) = ldn(l) + 1
        ENDDO
        lup(l) = ldn(l)
        DO WHILE ( pb(lup(l)+1)>plb(l+1) .AND. lup(l)<LMA )
           lup(l) = lup(l) + 1
        ENDDO
     ENDDO

     !**** Interpolate (extrapolate) vertically
     DO j = 1, jma
        DO l = 1, lm
           dh = 0.
           pdn = plb(l)
           IF ( lup(l)>0 ) THEN
              DO ll = ldn(l), lup(l)
                 pup = MAX(REAL(pb(ll+1),KIND=8),plb(l+1))
                 fracl = (pdn-pup)/(pb(ll)-pb(ll+1))
                 dh = dh + h2o(j,ll)*fracl*dz(ll)
                 pdn = pup
              ENDDO
           ENDIF
           DH2O(j,l,m) = 1.D-6*dh/1.74D0/DAYS_PER_YEAR
           !->(kg/m^2/ppm_CH4/day)
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE READ_QMA

SUBROUTINE LAT_INTERP_QMA(rlat,lev,mon,dh2o_interp)
  !@sum  interpolate CH4->H2O production rates in latitude
  !@auth R. Ruedy
  USE RAD_COM, ONLY : jma => JM_DH2O, XLAT => LAT_DH2O, DH2O
  IMPLICIT NONE
  REAL*8  ::  rlat
  ! input latitude (radians)
  INTEGER  ::  lev, mon
  ! input level, month
  REAL*8  ::  dh2o_interp
  ! output
  REAL*8 w1, w2
  INTEGER  ::  j1, j2

  !**** Interpolate (extrapolate) horizontally
  j2 = 2 + (jma-1)*(rlat-XLAT(1))/(XLAT(jma)-XLAT(1))
  ! first guess
  j2 = MIN(MAX(2,j2),jma)
  j1 = j2 - 1
  IF ( rlat>XLAT(j2) ) THEN ! j guess was too low
     DO WHILE ( j2<jma .AND. rlat>XLAT(j2) )
        j2 = j2 + 1
     ENDDO
     j1 = j2 - 1
  ELSEIF ( rlat<XLAT(j1) ) THEN ! j guess was too high
     DO WHILE ( j1>1 .AND. rlat<XLAT(j1) )
        j1 = j1 - 1
     ENDDO
     j2 = j1 + 1
  ENDIF
  !**** coeff. for latitudinal linear inter/extrapolation
  w1 = (XLAT(j2)-rlat)/(XLAT(j2)-XLAT(j1))
  !**** for extrapolations, only use half the slope
  IF ( w1>1. ) w1 = .5 + .5*w1
  IF ( w1<0. ) w1 = .5*w1
  w2 = 1. - w1
  dh2o_interp = w1*DH2O(j1,lev,mon) + w2*DH2O(j2,lev,mon)
END SUBROUTINE LAT_INTERP_QMA

#endif /* CUBED_SPHERE */

SUBROUTINE GETQMA(iu,dglat,plb,dh2o,lm,jm)
  !@sum  reads H2O production rates induced by CH4 (Tim Hall)
  !@auth R. Ruedy
  USE DOMAIN_DECOMP_ATM, ONLY : grid, GETDOMAINBOUNDS, WRITE_PARALLEL
  USE TIMECONSTANTS_MOD, ONLY : DAYS_PER_YEAR
  IMPLICIT NONE
  INTEGER, PARAMETER  ::  JMA = 18, LMA = 24
  INTEGER m, iu, jm, lm, j, j1, j2, l, ll, ldn(lm), lup(lm)
  REAL*8 PLB(lm+1), dH2O(grid%J_STRT_HALO:grid%J_STOP_HALO,lm,12),  &
       dglat(jm)
  REAL*4 pb(0 : LMA+1), h2o(JMA,0 : LMA), xlat(JMA), z(LMA), dz(0 : LMA)
  CHARACTER*100 title
  REAL*4 pdn, pup, w1, w2, dh, fracl
  INTEGER  ::  j_0, j_1
  CHARACTER(LEN=300)  ::  out_line
  CALL GETDOMAINBOUNDS(grid,J_STRT=J_0,J_STOP=J_1)

  !**** read headers/latitudes
  READ (iu,'(a)') title
  WRITE (out_line,'(''0'',a100)') title
  CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
  READ (iu,'(a)') title
  WRITE (out_line,'(1x,a100)') title
  CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
  READ (iu,'(a)') title
  !      write(6,'(1x,a100)') title
  READ (title(10 : 100),*) (xlat(j),j=1,JMA)

  !**** read heights z(km) and data (kg/km^3/year)
  DO m = 1, 12
     READ (iu,'(a)') title
     WRITE (out_line,'(1x,a100)') title
     CALL WRITE_PARALLEL(TRIM(out_line),UNIT=6)
     DO l = LMA, 1, -1
        READ (iu,'(a)') title
        !          write(6,'(1x,a100)') title
        READ (title,*) z(l), (H2O(j,l),j=1,JMA)
     ENDDO
     DO j = 1, JMA
        h2o(j,0) = 0.
     ENDDO

     !**** Find edge heights and pressures
     dz(0) = 0.
     dz(1) = z(2) - z(1)
     DO l = 2, LMA - 1
        dz(l) = .5*(z(l+1)-z(l-1))
     ENDDO
     dz(LMA) = z(LMA) - z(LMA-1)

     pb(0) = plb(1)
     DO l = 1, LMA
        Pb(l) = 1000.*10.**(-(z(l)-.5*dz(l))/16.)
     ENDDO
     !**** extend both systems vertically to p=0
     pb(LMA+1) = 0.
     plb(lm+1) = 0.

     !**** Interpolate vertical resolution to model layers
     ldn( : ) = 0
     DO l = 1, lm
        DO WHILE ( pb(ldn(l)+1)>=plb(l) .AND. ldn(l)<LMA )
           ldn(l) = ldn(l) + 1
        ENDDO
        lup(l) = ldn(l)
        DO WHILE ( pb(lup(l)+1)>plb(l+1) .AND. lup(l)<LMA )
           lup(l) = lup(l) + 1
        ENDDO
     ENDDO

     !**** Interpolate (extrapolate) horizontally and vertically
     j2 = 2
     DO j = j_0, j_1
        !**** coeff. for latitudinal linear inter/extrapolation
        DO WHILE ( j2<JMA .AND. dglat(j)>xlat(j2) )
           j2 = j2 + 1
        ENDDO
        j1 = j2 - 1
        w1 = (xlat(j2)-dglat(j))/(xlat(j2)-xlat(j1))
        !**** for extrapolations, only use half the slope
        IF ( w1>1. ) w1 = .5 + .5*w1
        IF ( w1<0. ) w1 = .5*w1
        w2 = 1. - w1
        DO l = 1, lm
           dh = 0.
           pdn = plb(l)
           IF ( lup(l)>0 ) THEN
              DO ll = ldn(l), lup(l)
                 pup = MAX(REAL(pb(ll+1),KIND=8),plb(l+1))
                 fracl = (pdn-pup)/(pb(ll)-pb(ll+1))
                 dh = dh + (w1*h2o(j1,ll)+w2*h2o(j2,ll))            &
                      *fracl*dz(ll)
                 pdn = pup
              ENDDO
           ENDIF
           dh2o(j,l,m) = 1.D-6*dh/1.74D0/DAYS_PER_YEAR
           !->(kg/m^2/ppm_CH4/day)
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE GETQMA

SUBROUTINE ORBIT(DOBLIQ,ECCEN,DOMEGVP,VEDAY,EDPY,DAY,SDIST,SIND,  &
     COSD,SUNLON,SUNLAT,EQTIME)
  !****
  !**** ORBIT receives orbital parameters and time of year, and returns
  !**** distance from Sun, declination angle, and Sun's overhead position.
  !**** Reference for following caculations is :   V.M.Blanco and
  !**** S.W.McCuskey, 1961, "Basic Physics of the Solar System", pages
  !**** 135 - 151.  Existence of Moon and heavenly bodies other than
  !**** Earth and Sun are ignored.  Earth is assumed to be spherical.
  !****
  !**** Program author :  Gary L. Russell 2004/11/16
  !**** Angles, longitude and latitude are measured in radians.
  !****
  !**** Input :  ECCEN  = eccentricity of the orbital ellipse
  !****        OBLIQ  = latitude of Tropic of Cancer
  !****        OMEGVP = longitude of perihelion (sometimes Pi is added) =
  !****               = spatial angle from vernal equinox to perihelion
  !****                 with Sun as angle vertex
  !****        DAY    = days measured since 2000 January 1, hour 0
  !****
  !****        EDPY  = Earth days per year
  !****                tropical year = 365.2425 (Gregorgian Calendar)
  !****                tropical year = 365      (Generic Year)
  !****        VEDAY = Vernal equinox
  !****                79.0 (Generic year Mar 21 hour 0)
  !****                79.5 (Generic year Mar 21 hour 12 - PMIP standard)
  !****                79.3125d0 for days from 2000 January 1, hour 0 till vernal
  !****                     equinox of year 2000 = 31 + 29 + 19 + 7.5/24
  !****
  !**** Intermediate quantities : 
  !****    BSEMI = semi minor axis in units of semi major axis
  !****   PERIHE = perihelion in days since 2000 January 1, hour 0
  !****            in its annual revolution about Sun
  !****       TA = true anomaly = spatial angle from perihelion to
  !****            current location with Sun as angle vertex
  !****       EA = eccentric anomaly = spatial angle measured along
  !****            eccentric circle (that circumscribes Earth's orbit)
  !****            from perihelion to point above (or below) Earth's
  !****            absisca (where absisca is directed from center of
  !****            eccentric circle to perihelion)
  !****       MA = mean anomaly = temporal angle from perihelion to
  !****            current time in units of 2*Pi per tropical year
  !****   TAofVE = TA(VE) = true anomaly of vernal equinox = - OMEGVP
  !****   EAofVE = EA(VE) = eccentric anomaly of vernal equinox
  !****   MAofVE = MA(VE) = mean anomaly of vernal equinox
  !****   SLNORO = longitude of Sun in Earth's nonrotating reference frame
  !****   VEQLON = longitude of Greenwich Meridion in Earth's nonrotating
  !****            reference frame at vernal equinox
  !****   ROTATE = change in longitude in Earth's nonrotating reference
  !****            frame from point's location on vernal equinox to its
  !****            current location where point is fixed on rotating Earth
  !****   SLMEAN = longitude of fictitious mean Sun in Earth's rotating
  !****            reference frame (normal longitude and latitude)
  !****
  !**** Output :  SIND = sine of declination angle = sin(SUNLAT)
  !****         COSD = cosine of the declination angle = cos(SUNLAT)
  !****       SUNDIS = distance to Sun in units of semi major axis
  !****       SUNLON = longitude of point on Earth directly beneath Sun
  !****       SUNLAT = latitude of point on Earth directly beneath Sun
  !****       EQTIME = Equation of Time =
  !****              = longitude of fictitious mean Sun minus SUNLON
  !****
  !**** From the above reference : 
  !**** (4-54) :  [1 - ECCEN*cos(EA)]*[1 + ECCEN*cos(TA)] = (1 - ECCEN^2)
  !**** (4-55) :  tan(TA/2) = sqrt[(1+ECCEN)/(1-ECCEN)]*tan(EA/2)
  !**** Yield :   tan(EA) = sin(TA)*sqrt(1-ECCEN^2) / [cos(TA) + ECCEN]
  !****    or :   tan(TA) = sin(EA)*sqrt(1-ECCEN^2) / [cos(EA) - ECCEN]
  !****
  USE CONSTANT, ONLY : twopi, pi, radian
  IMPLICIT NONE
  REAL*8, INTENT(IN)  ::  DOBLIQ, ECCEN, DOMEGVP, DAY, VEDAY, EDPY
  REAL*8, INTENT(OUT)  ::  SIND, COSD, SDIST, SUNLON, SUNLAT, EQTIME

  REAL*8 MA, OMEGVP, OBLIQ, EA, DEA, BSEMI, TAofVE, EAofVE, MAofVE, &
       SUNDIS, TA, SUNX, SUNY, SLNORO, VEQLON, ROTATE, SLMEAN
  !      REAL*8, PARAMETER :: EDAYzY=365.2425d0, VE2000=79.3125d0
  !      REAL*8, PARAMETER :: EDAYzY=365d0, VE2000=79d0  ! original parameters
  REAL*8 EDAYzY, VE2000
  !****
  VE2000 = VEDAY
  EDAYzY = EDPY
  OMEGVP = DOMEGVP*radian
  OBLIQ = DOBLIQ*radian
  !**** Determine EAofVE from geometry :  tan(EA) = b*sin(TA) / [e+cos(TA)]
  !**** Determine MAofVE from Kepler's equation :  MA = EA - e*sin(EA)
  !**** Determine MA knowing time from vernal equinox to current day
  !****
  BSEMI = SQRT(1-ECCEN*ECCEN)
  TAofVE = -OMEGVP
  EAofVE = ATAN2(BSEMI*SIN(TAofVE),ECCEN+COS(TAofVE))
  MAofVE = EAofVE - ECCEN*SIN(EAofVE)
  !     PERIHE = VE2000 - MAofVE*EDAYzY/TWOPI
  MA = MODULO(TWOPI*(DAY-VE2000)/EDAYzY+MAofVE,TWOPI)
  !****
  !**** Numerically invert Kepler's equation :  MA = EA - e*sin(EA)
  !****
  EA = MA + ECCEN*(SIN(MA)+ECCEN*SIN(2*MA)/2)
  DO
     dEA = (MA-EA+ECCEN*SIN(EA))/(1-ECCEN*COS(EA))
     EA = EA + dEA
     IF ( ABS(dEA)<=1D-10 ) THEN
        !****
        !**** Calculate distance to Sun and true anomaly
        !****
        SUNDIS = 1 - ECCEN*COS(EA)
        TA = ATAN2(BSEMI*SIN(EA),COS(EA)-ECCEN)
        SDIST = SUNDIS*SUNDIS
        ! added for compatiblity
        !****
        !**** Change reference frame to be nonrotating reference frame, angles
        !**** fixed according to stars, with Earth at center and positive x
        !**** axis be ray from Earth to Sun were Earth at vernal equinox, and
        !**** x-y plane be Earth's equatorial plane.  Distance from current Sun
        !**** to this x axis is SUNDIS sin(TA-TAofVE).  At vernal equinox, Sun
        !**** is located at (SUNDIS,0,0).  At other times, Sun is located at : 
        !****
        !**** SUN = (SUNDIS cos(TA-TAofVE),
        !****        SUNDIS sin(TA-TAofVE) cos(OBLIQ),
        !****        SUNDIS sin(TA-TAofVE) sin(OBLIQ))
        !****
        SIND = SIN(TA-TAofVE)*SIN(OBLIQ)
        COSD = SQRT(1-SIND*SIND)
        SUNX = COS(TA-TAofVE)
        SUNY = SIN(TA-TAofVE)*COS(OBLIQ)
        SLNORO = ATAN2(SUNY,SUNX)
        !****
        !**** Determine Sun location in Earth's rotating reference frame
        !**** (normal longitude and latitude)
        !****
        VEQLON = TWOPI*VE2000 - PI + MAofVE - TAofVE
        !  modulo 2*Pi
        ROTATE = TWOPI*(DAY-VE2000)*(EDAYzY+1)/EDAYzY
        SUNLON = MODULO(SLNORO-ROTATE-VEQLON,TWOPI)
        IF ( SUNLON>PI ) SUNLON = SUNLON - TWOPI
        SUNLAT = ASIN(SIN(TA-TAofVE)*SIN(OBLIQ))
        !****
        !**** Determine longitude of fictitious mean Sun
        !**** Calculate Equation of Time
        !****
        SLMEAN = PI - TWOPI*(DAY-FLOOR(DAY))
        EQTIME = MODULO(SLMEAN-SUNLON,TWOPI)
        IF ( EQTIME>PI ) EQTIME = EQTIME - TWOPI
        EXIT
     ENDIF
  ENDDO
  !****
END SUBROUTINE ORBIT

#ifdef HEALY_LM_DIAGS
REAL*8 FUNCTION FE(M,N)
  REAL*8 M, N

  FE = 0.47D0*LOG(1.+2.01D-5*(M*N)**(0.75)+5.31D-15*M*(M*N)**(1.52))
END FUNCTION FE
#endif

#ifdef CACHED_SUBDD
SUBROUTINE RIJH_DEFS(arr,nmax,decl_count)
  !
  ! 2D outputs
  !
  USE SUBDD_MOD, ONLY : INFO_TYPE, sched_rad
  ! info_type_ is a homemade structure constructor for older compilers
  USE SUBDD_MOD, ONLY : INFO_TYPE_
  IMPLICIT NONE
  INTEGER  ::  nmax, decl_count
  TYPE (INFO_TYPE) :: arr(nmax)
  !
  ! note :  next() is a locally declared function to increment decl_count
  !

  decl_count = 0

  !
  arr(NEXT()) = INFO_TYPE_(SNAME='olrrad',LNAME=                    &
       &'OUTGOING LW RADIATION at TOA (in RADIA)',          &
       &UNITS='W/m^2',SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='olrcs',LNAME=                     &
       &'OUTGOING LW RADIATION at TOA, CLEAR-SKY',          &
       &UNITS='W/m^2',SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='lwds',LNAME=                      &
       &'LONGWAVE DOWNWARD FLUX at SURFACE',UNITS='W/m^2',  &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='lwdscs',LNAME=                    &
       &'LONGWAVE DOWNWARD FLUX at SURFACE, CLEAR-SKY',     &
       UNITS='W/m^2',SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='lwus',LNAME=                      &
       &'LONGWAVE UPWARD FLUX at SURFACE',UNITS='W/m^2',    &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='totcld',LNAME=                    &
       &'Total Cloud Cover (as seen by rad)',UNITS='%',     &
       SCALE=1D2,SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='totcld_diag',LNAME=               &
       &'Total Cloud Cover (continuous, not seen by rad)',  &
       UNITS='%',SCALE=1D2,SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='cldss_2d',LNAME=                  &
       &'Stratiform Cloud Cover',UNITS='%',SCALE=1D2,       &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='cldmc_2d',LNAME=                  &
       &'Convective Cloud Cover',UNITS='%',SCALE=1D2,       &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='cod',LNAME=                       &
       &'Cloud optical depth warm clouds',UNITS='-',        &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='cid',LNAME=                       &
       &'Cloud optical depth ice clouds',UNITS='-',         &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='wtrcld',LNAME=                    &
       &'Water cloud frequency',UNITS='-',SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='icecld',LNAME=                    &
       &'Ice cloud frequency',UNITS='-',SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='ctt',LNAME='Cloud top temperature'&
       ,UNITS='C',SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='ctp',LNAME='Cloud top pressure',  &
       UNITS='hPa',SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swds',LNAME=                      &
       &'SOLAR DOWNWARD FLUX at SURFACE',UNITS='W/m^2',     &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swus',LNAME=                      &
       &'SOLAR UPWARD FLUX at SURFACE',UNITS='W/m^2',       &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swdf',LNAME=                      &
       &'SOLAR DOWNWARD DIFFUSE FLUX at SURFACE',           &
       &UNITS='W/m^2',SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swtoa',LNAME='SOLAR NET FLUX, TOA'&
       ,UNITS='W/m^2',SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swns',LNAME=                      &
       &'Solar net flux at surface',UNITS='W/m^2',          &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='lwns',LNAME=                      &
       &'Longwave net flux at surface',UNITS='W/m^2',       &
       SCHED=sched_rad)
  !
#ifdef CFMIP3_SUBDD  /* CFMIP3_SUBDD */
  arr(NEXT()) = INFO_TYPE_(SNAME='rtmt',LNAME=                      &
       &'Net downward radiative flux, TOA',UNITS='W/m^2',   &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swut',LNAME='TOA outgoing SW',    &
       UNITS='W/m^2',SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swutcs',LNAME=                    &
       &'TOA outgoing SW, CSKY',UNITS='W/m^2',              &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='clwvi',LNAME='Total water path',  &
       UNITS='kg/m^2',SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swdcls',LNAME=                    &
       &'SFC downward radiative flux, CSKY',UNITS='kg/m^2', &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swucls',LNAME=                    &
       &'SFC upward radiative flux, CSKY',UNITS='kg/m^2',   &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swus',LNAME=                      &
       &'SFC upward radiative flux',UNITS='kg/m^2',         &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swdt',LNAME='TOA incoming SW',    &
       UNITS='W/m^2',SCHED=sched_rad)
#endif  /* CFMIP3_SUBDD */
  RETURN
CONTAINS
  INTEGER FUNCTION NEXT()
    decl_count = decl_count + 1
    NEXT = decl_count
  END FUNCTION NEXT
END SUBROUTINE RIJH_DEFS

SUBROUTINE RIJLH_DEFS(arr,nmax,decl_count)
  !
  ! 3D outputs
  !
  USE SUBDD_MOD, ONLY : INFO_TYPE, sched_rad
  ! info_type_ is a homemade structure constructor for older compilers
  USE SUBDD_MOD, ONLY : INFO_TYPE_
  USE CONSTANT, ONLY : kapa
  USE TIMECONSTANTS_MOD, ONLY : SECONDS_PER_DAY
  IMPLICIT NONE
  INTEGER  ::  nmax, decl_count
  TYPE (INFO_TYPE) :: arr(nmax)
  !
  ! note :  next() is a locally declared function to increment decl_count
  !
  decl_count = 0
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='dth_sw',LNAME=                    &
       &'theta tendency from shortwave radiative heating',  &
       UNITS='K/day',SCALE=1000.**kapa*SECONDS_PER_DAY,    &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='dth_lw',LNAME=                    &
       &'theta tendency from longwave radiative heating',   &
       UNITS='K/day',SCALE=1000.**kapa*SECONDS_PER_DAY,    &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='dth_rad',LNAME=                   &
       &'theta tendency from radiative heating',            &
       &UNITS='K/day',SCALE=1000.**kapa*SECONDS_PER_DAY,    &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='lwdp',LNAME=                      &
       &'LONGWAVE DOWNWARD FLUX profile',UNITS='W/m^2',     &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='lwup',LNAME=                      &
       &'LONGWAVE UPWARD FLUX profile',UNITS='W/m^2',       &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swdp',LNAME=                      &
       &'SHORTWAVE DOWNWARD FLUX profile',UNITS='W/m^2',    &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='swup',LNAME=                      &
       &'SHORTWAVE UPWARD FLUX profile',UNITS='W/m^2',      &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='MRCO2rad',LNAME=                  &
       &'radiation code CO2 volume mixing ratio',           &
       &UNITS='mole species / mole air',SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='wtrtau',LNAME=                    &
       &'Cloud Water Opacity Seen by Radiation',UNITS='-',  &
       SCHED=sched_rad)
  !
  arr(NEXT()) = INFO_TYPE_(SNAME='icetau',LNAME=                    &
       &'Cloud Ice Opacity Seen by Radiation',UNITS='-',    &
       SCHED=sched_rad)
  !
#ifdef CFMIP3_SUBDD
  arr(NEXT()) = INFO_TYPE_(SNAME='cf',LNAME='Cloud Fraction',       &
       &UNITS='%',SCALE=1D2,SCHED=sched_rad)
  arr(NEXT()) = INFO_TYPE_(SNAME='qcirad',LNAME=                    &
       &'Ice Water Mass Mixing Ratio Seen by Radiation',    &
       UNITS='kg/kg',SCHED=sched_rad)
  arr(NEXT()) = INFO_TYPE_(SNAME='qclrad',LNAME=                    &
       &'Liquid Water Mass Mixing Ratio Seen by Radiation', &
       UNITS='kg/kg',SCHED=sched_rad)
#endif
  !
  RETURN
CONTAINS
  INTEGER FUNCTION NEXT()
    decl_count = decl_count + 1
    NEXT = decl_count
  END FUNCTION NEXT
END SUBROUTINE RIJLH_DEFS

#endif

SUBROUTINE READIFILE(IFile)
  ! Consolidated duplicate of MODELE.f code snippets that read the
  ! I-file containing the parameter database and INPUTZ namelist.
  ! Currently used by radiation-only configuration; to be moved
  ! to MODELE.f and used by all configurations once full testing
  ! is completed.
  ! Note that INPUTZ contains fewer variables in this version.
  USE FILEMANAGER, ONLY : OPENUNIT, CLOSEUNIT
  USE PARSER_MOD
  USE MODEL_COM, ONLY : xlabel, lrunid
  USE MODEL_COM, ONLY : HOURI, DATEI, MONTHI, YEARI, IYEAR1
  USE DIAG_COM, ONLY : itwrite
  IMPLICIT NONE
  !**** Command line options
  CHARACTER(LEN=*), INTENT(IN)  ::  IFile

  INTEGER  ::  iu_IFILE
  !@nlparam IHRI,TIMEE,IHOURE   end of model run
  !@var  IHRI,IHOURE start and end of run in hours (from 1/1/IYEAR1 hr 0)
  CHARACTER NLREC*80, RLABEL*132
  !****    List of parameters that are disregarded at restarts
  NAMELIST /INPUTZ/ ITWRITE, HOURI, DATEI, MONTHI, YEARI
  !****    List of parameters that are disregarded at restarts
  NAMELIST /INPUTZ_COLD/ ITWRITE, HOURI, DATEI, MONTHI, YEARI
  CHARACTER*132  ::  bufs
  INTEGER, PARAMETER  ::  MAXLEN_RUNID = 32
  INTEGER  ::  lid1, lid2, fid, noff


  !****
  !**** Reading rundeck (I-file) options
  !****
  CALL OPENUNIT(TRIM(ifile),iu_IFILE,.FALSE.,.TRUE.)
  CALL PARSE_PARAMS(iu_IFILE)
  CALL CLOSEUNIT(iu_IFILE)

  !****
  !**** Print Header and Label (2 lines) from rundeck
  !****
  CALL OPENUNIT(TRIM(ifile),iu_IFILE,.FALSE.,.TRUE.)
  !if (AM_I_ROOT())
  WRITE (6,'(A,40X,A/)') '0', 'GISS CLIMATE MODEL'
  READ (iu_IFILE,'(A80)') XLABEL(1 : 80), NLREC
  NOFF = 0
  IF ( XLABEL(73 : 80)=='        ' ) NOFF = 8 ! for 72-column rundecks
  XLABEL(81-NOFF : 132) = NLREC(1 : 52+NOFF)
  !if (AM_I_ROOT())
  WRITE (6,'(A,A/)') '0', XLABEL
  RLABEL = XLABEL !@var RLABEL rundeck-label

  lid1 = INDEX(XLABEL,'(') - 1
  IF ( lid1<1 ) lid1 = MAXLEN_RUNID + 1
  lid2 = INDEX(XLABEL,' ') - 1
  IF ( lid2<1 ) lid2 = MAXLEN_RUNID + 1
  LRUNID = MIN(lid1,lid2)
  IF ( LRUNID>MAXLEN_RUNID ) CALL STOP_MODEL(                       &
       &'INPUT :  Rundeck name too long. Shorten to 32 char or less',  &
       255)

  !****
  !**** Read parameters from the rundeck to database and namelist
  !****
  DO
     READ (iu_IFILE,*,ERR=910,END=910) bufs
     ! achar(38) is an ampersand
     IF ( bufs == achar(38)//achar(38)//'END_PARAMETERS' ) EXIT
  ENDDO

  READ (iu_IFILE,NML=INPUTZ,ERR=900)

  CALL CLOSEUNIT(iu_IFILE)

  IF ( YearI<0 ) THEN
     WRITE (6,*) 'Please choose a proper start year yearI, not',    &
          yearI
     CALL STOP_MODEL('INPUT :  yearI not provided',255)
  ENDIF

  RETURN
  !****
  !**** TERMINATE BECAUSE OF IMPROPER PICK-UP
  !****
900 WRITE (6,*) 'Error in NAMELIST parameters'
  CALL STOP_MODEL('Error in NAMELIST parameters',255)
910 WRITE (6,*) 'Error readin I-file'
  CALL STOP_MODEL('Error reading I-file',255)

END SUBROUTINE READIFILE

SUBROUTINE RUN_RADONLY(IFile)
  !@sum Call single-column radiation-only model once
  USE DICTIONARY_MOD
  USE DOMAIN_DECOMP_1D, ONLY : INIT_APP
  USE MODEL_COM, ONLY : itime, itimeE, master_yr, xlabel, lrunid
  USE MODEL_COM, ONLY : YEARI, IYEAR1
  USE DOMAIN_DECOMP_ATM, ONLY : grid, INIT_GRID
#ifdef CACHED_SUBDD
  USE DIAG_COM
  USE GEOM, ONLY : lon_dg, lat_dg
  USE CDL_MOD
#endif
  USE TIMERPACKAGE_MOD,                                             &
       ONLY : INITIALIZETIMERPACKAGE_MOD => INITIALIZE
  IMPLICIT NONE
  !**** Command line options
  CHARACTER(LEN=*), INTENT(IN)  ::  IFile
  !
  INTEGER  ::  i, j, l, n
  CHARACTER(LEN=80)  ::  filenm

  CALL INIT_APP()

  CALL INITIALIZETIMERPACKAGE_MOD() ! avoid probs when RADIA calls timers

  CALL READIFILE(IFile)

  IF ( IS_SET_PARAM("master_yr") ) THEN
     CALL GET_PARAM("master_yr",master_yr)
  ELSE
     CALL STOP_MODEL('Please define master_yr in the rundeck.',255)
  ENDIF

  Iyear1 = yearI

  CALL SUNDIAL

  itimeE = itime + 1
  ! for length-1 nominal time axis in diags

  CALL INIT_GRID(grid,1,1,1,WIDTH=0)

  !call alloc_clouds_com(grid)
  CALL ALLOC_RAD_COM(grid)
  !call alloc_veg_com(grid)

  CALL GEOM_1PT

  CALL INIT_RAD(2) ! istart=2
  CALL DAILY_ORBIT(.FALSE.)             ! not end_of_day
  CALL DAILY_RAD(.FALSE.)

  CALL PRINT_PARAM(6)

#ifdef CACHED_SUBDD
  ! Initialize diagnostics framework
  CALL INIT_CDL_TYPE('cdl_aij',cdl_ij_template)
  CALL ADD_COORD(cdl_ij_template,'lon',1,UNITS='degrees_east',      &
       COORDVALUES=lon_dg( : ,1))
  CALL ADD_COORD(cdl_ij_template,'lat',1,UNITS='degrees_north',     &
       COORDVALUES=lat_dg( : ,1))
  CALL PARSE_SUBDD
  CALL RESET_CACHED_SUBDD
  CALL GET_SUBDD_VINTERP_COEFFS
  CALL SET_SUBDD_PERIOD()
#endif

  CALL CALC_ZENITH_ANGLE
  CALL RADIA

#ifdef CACHED_SUBDD
  filenm = 'allsteps.subdd'//XLABEL(1 : LRUNID)
  CALL WRITE_SUBDD_ACCFILE(filenm)
#endif

  CALL STOP_MODEL('Radiation calculations completed.',13)

CONTAINS

  SUBROUTINE GEOM_1PT
    USE GEOM
    USE CONSTANT, ONLY : pi, twopi, radian
    USE DICTIONARY_MOD, ONLY : GET_PARAM, SYNC_PARAM
    IMPLICIT NONE
    REAL*8  ::  lon_targ, lat_targ

    ! mandatory rundeck parameters :  lon and lat of target point
    CALL GET_PARAM('lon_targ',lon_targ)
    CALL GET_PARAM('lat_targ',lat_targ)

    IF ( ABS(lon_targ)>180D0 .OR. ABS(lat_targ)>90D0 )                &
         CALL STOP_MODEL(                                             &
         &'geom_atm :  invalid lon_targ,lat_targ in rundeck',255)

    LON2D_DG(1,1) = lon_targ
    LAT2D_DG(1,1) = lat_targ

    AXYP(1,1) = 1.

    BYAXYP(1,1) = 1D0/AXYP(1,1)

    LON2D(1,1) = LON2D_DG(1,1)*radian
    LAT2D(1,1) = LAT2D_DG(1,1)*radian

    SINLAT2D(1,1) = SIN(LAT2D(1,1))
    COSLAT2D(1,1) = COS(LAT2D(1,1))
    LON2D(1,1) = LON2D(1,1) + pi ! IDL has a value of zero
    IF ( LON2D(1,1)<0. ) LON2D(1,1) = LON2D(1,1) + twopi

    imaxj = 1

    lon_dg = LON2D_DG
    lat_dg = LAT2D_DG

  END SUBROUTINE GEOM_1PT

  SUBROUTINE SUNDIAL
    ! Duplicate of relevant snippets of clock initialization in MODELE.f.
    ! Currently used by radiation-only configuration; will disappear
    ! once the clock initialization in MODELE.f has been cleanly isolated
    ! from other intialization activities.
    USE DICTIONARY_MOD
    USE MODEL_COM, ONLY : nday, dtsrc, itime, itimei, HOURI, DATEI,     &
         MONTHI, YEARI
    USE MODEL_COM, ONLY : modelEclock, calendar
    USE MODELCLOCK_MOD, ONLY : MODELCLOCK
#ifdef TRACERS_GC
    USE TEMPUS_MOD
#else
    USE TIME_MOD
#endif
    USE BASETIME_MOD
    USE RATIONAL_MOD
    USE TIMEINTERVAL_MOD
    IMPLICIT NONE

    TYPE (TIME) :: MODELETIME0
    TYPE (TIME) :: MODELETIME
    TYPE (TIMEINTERVAL) :: dtSrcUsed
    TYPE (TIMEINTERVAL) :: secsPerDay

    !**** Get those parameters which are needed in this subroutine
    CALL GET_PARAM("DTsrc",DTsrc)

    !@var NDAY=(1 day)/DTsrc  :  even integer; adjust DTsrc to be commensurate
    secsPerDay = calendar%GETSECONDSPERDAY()
    NDAY = 2*NINT((secsPerDay/(DTsrc*2)))
    dtSrcUsed = TIMEINTERVAL(secsPerDay/NDAY)
    DTsrc = REAL(dtSrcUsed)

    MODELETIME0 = NEWTIME(calendar)
    MODELETIME = NEWTIME(calendar)

    CALL MODELETIME%SETBYDATE(yearI,monthI,dateI,hourI)
    CALL MODELETIME0%SETBYDATE(yearI,MONTH=1,DATE=1,HOUR=0)

    ITimeI = NINT((MODELETIME-MODELETIME0)/dtSrcUsed)
    Itime = ItimeI

    modelEclock = MODELCLOCK(MODELETIME,dtSrcUsed,itime)

    CALL DAILY_CAL(.FALSE.)   ! not end_of_day

  END SUBROUTINE SUNDIAL

END SUBROUTINE RUN_RADONLY

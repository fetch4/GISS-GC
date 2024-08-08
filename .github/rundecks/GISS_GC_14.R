GISS_GC_14.R GISS ModelE Lat-Lon Atmosphere Model, 1850 atm./ocean

! GISS_GC_14 is based on E6F40 with updated aerosol/ozone input files for CMIP6
! simulations
!
! It uses GEOS-Chem at version 14.3.1
!
! Lat-lon: 2x2.5 degree horizontal resolution
! F40: 40 vertical layers with standard hybrid coordinate, top at .1 mb
! Atmospheric composition for year 1850
! Ocean climatology prescribed from years 1876-1885, CMIP6
! Uses turbulence scheme (no dry conv), grav.wave drag
! timesteps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
! Filters: U,V in E-W and N-S direction (after every physics timestep)
!          U,V in E-W direction near poles (after every dynamics timestep)
!          sea level pressure (after every physics timestep)

Preprocessor Options
#define STDHYB                   ! standard hybrid vertical coordinate
#define ATM_LAYERING L40         ! 40 layers, top at .1 mb
#define NEW_IO                   ! new I/O (netcdf) on
#define IRRIGATION_ON
#define MODIS_LAI
#define NEW_BCdalbsn
#define NEW_IO_SUBDD
#define CACHED_SUBDD
#define GCAP
#define CALCULATE_LIGHTNING
#define TRACERS_GC               ! tracers using GISS-GC coupling
#define MERRA_NUDGING
#define NUDGE_ON
End Preprocessor Options

Object modules:

! resolution-specific source codes
Atm144x90                           ! horizontal resolution 144x90 -> 2x2.5deg
AtmLayering                         ! vertical resolution
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform

IO_DRV                              ! new i/o

! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)
NUDGE

! latitude-longitude grid specific source codes
AtmRes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT                       ! diagn/post-processing output
MODEL_COM                           ! calendar, timing variables
MODELE_DRV                          ! ModelE cap
MODELE                              ! initialization and main loop
ATM_COM                             ! main atmospheric variables
ATM_DRV                             ! driver for atmosphere-grid components
ATMDYN_COM                          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
CHEM_DRV CHEM_COM                   ! GEOS-Chem
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES      ! surface calculation and fluxes
GHY_COM GHY_DRV                     ! land surface and soils + snow model
VEG_DRV                             ! vegetation
ENT_DRV ENT_COM                     ! new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
IRRIGMOD                            ! irrigation module
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV     ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO ocalbedo ! radiation and albedo
DIAG_COM DIAG DEFACC                ! diagnostics
OCN_DRV                             ! driver for ocean-grid components
OCEAN OCNML                         ! ocean modules

lightning
SUBDD

Components:
shared MPI_Support solvers giss_LSM
dd2d
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT

Data input files:
! Atmospheric initial conditions for automatic relayering to model vertical grid
AIC=NCARIC.144x90.D7712010_ext.nc
! Ground initial condition
GIC=GIC.144X90.DEC01.1.ext_1.nc

OSST=OST_144x90.1876-1885avg.CMIP6.nc       ! climatological ocean temperature
SICE=SICE_144x90.1876-1885avg.CMIP6.nc      ! climatological sea ice cover
ZSIFAC=ZSIfac_144x90.1876-1885avg.CMIP6.nc  ! climatological sea ice thickness
TOPO=Z2HX2fromZ1QX1N.BS1.nc                 ! ocean frac. and surface topography
RVR=RD_Fd.nc                                ! river direction file
NAMERVR=RD_Fd.names.txt                     ! named river outlets

CDN=CD144X90.ext.nc
VEG=V144x90_EntMM16_lc_max_trimmed_scaled_nocrops.ext.nc
LAIMAX=V144x90_EntMM16_lai_max_trimmed_scaled_ext.nc
HITEent=V144x90_EntMM16_height_trimmed_scaled_ext.nc
LAI=V144x90_EntMM16_lai_trimmed_scaled_ext.nc
CROPS=CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp.nc
IRRIG=Irrig144x90_1848to2100_FixedFuture_v3.nc
SOIL=S144X900098M.ext.nc
TOP_INDEX=top_index_144x90_a.ij.ext.nc
ZVAR=ZVAR2X25A.nc ! topographic variation for gravity wave drag

! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_2x2.5.nc
GLMELT=GLMELT_144X90_gas.OCN.nc
RADN1=sgpgxg.table8                           ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800  ! rad.tables and history files
RADN4=LWCorrTables33k                         ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
RADN3=miescatpar.abcdv2

RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-2014_CMIP6_hdr                  ! needs MADVOL=2
RADN8=cloud.epsilon4.72x46
RADN9=solar.CMIP6official.ann1850-2299_with_E3_fastJ.nc ! needs KSOLAR=2
RADNE=topcld.trscat8

ISCCP=ISCCP.tautables
GHG=GHG.CMIP6.1-2014.txt  !  GreenHouse Gases for CMIP6 runs up to 2014
CO2profile=CO2profile.Jul16-2017.txt ! scaling of CO2 in stratosphere
dH2O=dH2O_by_CH4_monthly

! NINT E2.1 input files
BCdalbsn=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/BCdalbsn
DUSTaer=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/DUST
TAero_SUL=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/SUL
TAero_SSA=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/SSA
TAero_NIT=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/NIT
TAero_OCA=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/OCA
TAero_BCA=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/BCA
TAero_BCB=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/BCB
u2014.nc=nudging/merra2/uwnd.2014.MERRA2onGISSE2.nc4
v2014.nc=nudging/merra2/vwnd.2014.MERRA2onGISSE2.nc4
u2015.nc=nudging/merra2/uwnd.2015.MERRA2onGISSE2.nc4
v2015.nc=nudging/merra2/vwnd.2015.MERRA2onGISSE2.nc4
u2016.nc=nudging/merra2/uwnd.2016.MERRA2onGISSE2.nc4
v2016.nc=nudging/merra2/vwnd.2016.MERRA2onGISSE2.nc4
O3file=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/O3
Ox_ref=o3_2010_shindell_144x90x49_April1850.nc

MSU_wts=MSU_SSU_RSS_weights.txt   ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
GISS_GC_14 (LLF40 + updated aerosol/ozone input files for CMIP6 simulations, 1850 atm/ocean)

&&PARAMETERS
! parameters set for choice of ocean model:
KOCEAN=0        ! ocean is prescribed
Kvflxo=0        ! usually set to 1 only during a prescr.ocn run by editing "I"
variable_lk=1   ! variable lakes

! drag params if gravity wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33   34   35   36  37  38   39 40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,0.00,0.00,0.00,0.3,0.6,0.83,1.

! Gravity wave parameters
PBREAK = 200.       ! the level for GW breaking above.
DEFTHRESH=0.000055  ! threshold (1/s) for triggering deformation waves
PCONPEN=400.        ! penetrating convection defn for GWDRAG
CMC = 0.0000002     ! parameter for GW Moist Convective drag
CSHEAR=10.          ! shear drag coefficient
CMTN=0.1            ! default is 0.5
CDEF=1.6            ! tuning factor for deformation -> momentum flux
XCDNST=400.,10000.  ! strat. gw drag parameters
QGWMTN=1            ! mountain waves ON
QGWDEF=1            ! deformation waves ON
QGWSHR=0            ! shear drag OFF
QGWCNV=0            ! convective drag OFF

! following two lines are only used when aerosol/radiation interactions are off
FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,1.
FT8OPX=1.,1.,1.,1.,1.,1.,1.3,1.

! increasing U00a decreases the high cloud cover (tune first)
U00a=0.655 ! above 850mb w/o MC region; tune to get 30-35% high clouds
! increasing U00b decreases net rad at TOA (tune last)
U00b=1.00  ! below 850mb and MC regions; tune this to get radiative balance
WMUI_multiplier=2.
use_vmp=1
radius_multiplier=1.1

PTLISO=0.   ! pressure(mb) above which radiation assumes isothermal layers
H2ObyCH4=1. ! activates stratospheric H2O generated by CH4 without interactive chemistry
KSOLAR=2    ! use long annual mean file

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
master_yr=0        ! transient run
volc_yr=-1
od_cdncx=0.        ! do not include 1st indirect effect
cc_cdncx=0.        ! do not include 2nd indirect effect (used 0.0036)
dalbsnX=1.

MADVOL=2

! GEOS-Chem Operators (1=true, 0=false) ! TODO: Turn these on
! DoGCConv=1
! DoGCEmis=1
! DoGCTend=0
! DoGCTurb=1
! DoGCChem=1
! DoGCDryDep=1
! DoGCWetDep=1

DTsrc=1800.      ! physics timestep (cannot be changed after a run starts)
DT=225.          ! advection timestep

! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! surface interaction computed NIsurf times per source timestep
NRAD=1           ! radiation computed NRAD times per source timestep
! parameters that affect at most diagn. output:  standard if DTsrc=1800. (sec)
TAero_aod_diag=2 ! save band6 only
aer_rad_forc=0   ! turn off aerosol radiative forcing diagnostics
cloud_rad_forc=1 ! turn on cloud radiative forcing diagnostics

! diagnostics
! SUBDD='OH:4 NO:4 O3:4 NO2:4 CO:4 CH4:4 PS:4' ! TODO: Turn these on
SUBDD='SAT:6i'
NSUBDD=1         ! saving sub-daily diags every NSUBDD-th physics timestep
DAYS_PER_FILE=1
KCOPY=1          ! save accumulated diagnostics files
KRSF=12          ! save restart file at the beginning of every 12 months
isccp_diags=1    ! include all key diagnostics
nda5d=13
nda5s=13
ndaa=13
nda5k=13
nda4=48
Nssw=2
Ndisk=960        ! write fort.1.nc or fort.2.nc every NDISK source timestep
&&END_PARAMETERS

&INPUTZ
 YEARI=2016,MONTHI=7,DATEI=1,HOURI=0,
 YEARE=2016,MONTHE=8,DATEE=1,HOURE=0,     KDIAG=13*0,
 ISTART=2,IRANDI=0, YEARE=2016,MONTHE=7,DATEE=2,HOURE=0,
/

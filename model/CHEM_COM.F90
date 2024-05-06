module CHEM_COM

  !================================================================================================
  ! Module CHEM_COM is a module that holds diagnostic arrays for
  ! GEOS-Chem
  ! 
  ! Author: Lee T. Murray (lee.murray@rochester.edu)
  !===============================================================================================
  USE MDIAG_COM,   ONLY : sname_strlen,lname_strlen,units_strlen
  USE CDL_MOD

  IMPLICIT NONE
  SAVE

  !PUBLIC :: Init_Chem_Diagnostics

  INTEGER, PUBLIC                                       :: NTM        ! Number of tracers to advect
  REAL*8,  PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:)      :: TrM        ! Tracer array (kg)
  REAL*8,  PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:,:)    :: TrMom      ! Second order moments for tracers (kg)

  CHARACTER(LEN=8), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: TrName     ! Species name
  CHARACTER(LEN=163), PUBLIC, ALLOCATABLE, DIMENSION(:) :: TrFullName ! Full name
  LOGICAL, PUBLIC, ALLOCATABLE, DIMENSION(:)            :: IsAdvected ! Advect this tracer?
  LOGICAL, PUBLIC, ALLOCATABLE, DIMENSION(:)            :: t_qlimit   ! Limit fluxes in QUS?    

  !====================
  ! Diagnostics
  !====================

      INTEGER, PARAMETER :: kgcaijl=500

!@var SNAME_IJLT: Names of 3D tracer IJL diagnostics
      character(len=sname_strlen), dimension(kgcaijl) :: sname_ijlt
!@var DNAME_IJLT, DENOM_IJLT: Short names, indices of gcaijls denominators.
!@+   Currently, dname is specified along with the standard metadata and
!@+   the denom indices are looked up afterward.
      character(len=sname_strlen), dimension(kgcaijl) :: dname_ijlt=''
      integer, dimension(kgcaijl) :: denom_ijlt=0
!@var LNAME_IJLT,UNITS_IJLT: descriptions/units of 3D tracer diagnostics
      character(len=lname_strlen), dimension(kgcaijl) :: lname_ijlt = 'unused'
      character(len=units_strlen), dimension(kgcaijl) :: units_ijlt
!@var SCALE_IJLT: printout scaling factor for 3D tracer diagnostics
      REAL*8, dimension(kgcaijl) :: scale_ijlt
!@var IR_IJLT: range index of IJL diagnostics
      integer, dimension(kgcaijl) :: ir_ijlt
!@var IA_IJLT: accumulation index for IJL diagnostics
      integer, dimension(kgcaijl) :: ia_ijlt
!@var ijlt_power: power of 10 used for tracer IJL 3D diags
      INTEGER, DIMENSION(kgcaijl) :: ijlt_power

      ! Indices to radiatively active species
      INTEGER :: i_H2O   = 0
      INTEGER :: i_CO2   = 0
      INTEGER :: i_O3    = 0
      INTEGER :: i_O2    = 0
      INTEGER :: i_NO2   = 0
      INTEGER :: i_N2O   = 0
      INTEGER :: i_CH4   = 0
      INTEGER :: i_CFC11 = 0
      INTEGER :: i_CFC12 = 0
      INTEGER :: i_N2    = 0
      INTEGER :: i_CFCY  = 0
      INTEGER :: i_CFCZ  = 0
      INTEGER :: i_SO2   = 0
      INTEGER :: i_SO4   = 0
      INTEGER :: i_NIT   = 0
      INTEGER :: i_OCPI  = 0
      INTEGER :: i_OCPO  = 0
      INTEGER :: i_BCPI  = 0
      INTEGER :: i_BCPO  = 0
      INTEGER :: i_SOAS  = 0
      
! This section declares arrays used to write tracer
! diagnostics accumulations in a format suitable for offline
! postprocessing.  The size of these arrays cannot be known
! precisely a priori due to the large number of outputs not
! declared in advance.  Thus, metadata arrays are declared with
! a larger-than-needed size, and acc arrays are re-allocated
! to the correct size after the size is determined.
! The "n" and "s" instances of diagnostics classes are
! merged as follows:
!@var gcaijl_out combines gcaijln and gcaijls.  Distributed.
!@var gcaij_out  combines gcaijn  and gcaijs.   Distributed.
!@var tajl_out  combines tajln  and tajls.
!@var tconsrv_out combines the ktcon/ntmxcon dims of tconsrv
!@+               and omits its unused elements.
! All contents of txxx_out arrays are PER UNIT AREA.
! Denominator information is stored in denom_xxx.
! Metadata for scaled outputs is consolidated in cdl_xxx.
! 
      ! This was NTM, but NTM is now dynamic. Introducing MAXNTM instead.
      !integer, parameter :: MAXNTM = 1000 ! rather than making all kgcaij_ arrays allocatable - TLC
      !integer, parameter :: kgcaij_ = (kgcaij*MAXNTM+kgcaijs)*3/2  ! make 50% larger for denoms and extra specials
      !integer :: kgcaij_out ! actual number of qtys in gcaij_out
      !real*8, dimension(:,:,:), allocatable :: gcaij_out
      !integer, dimension(kgcaij_) :: ir_gcaij,ia_gcaij,denom_gcaij
      !character(len=lname_strlen), dimension(kgcaij_) :: lname_gcaij
      !character(len=sname_strlen), dimension(kgcaij_) :: sname_gcaij
      !character(len=units_strlen), dimension(kgcaij_) :: units_gcaij
      !real*8, dimension(kgcaij_) :: scale_gcaij
      !type(cdl_type) :: cdl_gcaij,cdl_gcaij_latlon
      !real*8, dimension(:,:,:), allocatable :: hemis_gcaij

      integer :: kgcaijl_
      integer :: kgcaijl_out ! actual number of qtys in gcaijl_out
      real*8, dimension(:,:,:,:), allocatable :: gcaijl_out
      integer, allocatable, dimension(:) ::ir_gcaijl,ia_gcaijl,denom_gcaijl,ijlt_vmr
      character(len=lname_strlen),allocatable,dimension(:) ::lname_gcaijl
      character(len=sname_strlen),allocatable,dimension(:) ::sname_gcaijl
      character(len=units_strlen),allocatable,dimension(:) ::units_gcaijl
      real*8, allocatable, dimension(:) :: scale_gcaijl
      type(cdl_type) :: cdl_gcaijl,cdl_gcaijl_latlon

      !integer :: ktajl_
      !integer :: ktajl_out ! actual number of qtys in tajl_out
      !real*8, dimension(:,:,:), allocatable :: tajl_out
      !integer, allocatable, dimension(:) :: pow_tajl,ia_tajl,denom_tajl,&
      !   jgrid_tajl,lgrid_tajl,ltop_tajl
      !character(len=lname_strlen),allocatable,dimension(:) :: lname_tajl
      !character(len=sname_strlen),allocatable,dimension(:) :: sname_tajl
      !character(len=units_strlen),allocatable,dimension(:) :: units_tajl
      !real*8, allocatable, dimension(:) :: scale_tajl
      !type(cdl_type) :: cdl_tajl
      !real*8, dimension(:,:,:), allocatable :: hemis_tajl,vmean_tajl

  CONTAINS

end module CHEM_COM

      subroutine def_rsf_trdiag(fid,r4_on_disk)
!@sum  def_rsf_trdiag defines tracer diag array structure in restart+acc files
!@auth M. Kelley
!@ver  beta
      use CHEM_COM, only : GCAIJL=>GCAIJL_out ! , &
!                            GCAIJLN=>GCAIJLN_loc, &
!                            GCAIJLS=>GCAIJLS_loc, &
!                            GCAIJN=>GCAIJN_loc, &
!                            GCAIJS=>GCAIJS_loc, &
!                            TAJLN, &
!                            TAJLS, &
!                            GCAIJ=>GCAIJ_out, &
!                            TAJL=>TAJL_out
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid           !@var fid file id
      logical :: r4_on_disk !@var r4_on_disk if true, real*8 stored as real*4

      if(r4_on_disk) then ! acc file
        call defvar(grid,fid,gcaijl,'gcaijl(dist_im,dist_jm,lm,kgcaijl)',r4_on_disk=.true.)
!         call defvar(grid,fid,gcaij,
!         'gcaij(dist_im,dist_jm,kgcaij)',r4_on_disk=.true.)
!         call defvar(grid,fid,tajl, 'tajl(jm_budg,lm,ktajl)',r4_on_disk=.true.)
!        call write_src_dist_data(fid, .true.)
!       else
!         call defvar(grid,fid,gcaijln,'gcaijln(dist_im,dist_jm,lm,ntm)')
!         call defvar(grid,fid,gcaijls,'gcaijls(dist_im,dist_jm,lm,kgcaijl)')
!         call defvar(grid,fid,gcaijs,'gcaijs(dist_im,dist_jm,kgcaijs)')
!         call defvar(grid,fid,gcaijn,'gcaijn(dist_im,dist_jm,kgcaij,ntm)')
!         call defvar(grid,fid,tajln,'tajln(jm_budg,lm,ktajlx,ntm)')
!         call defvar(grid,fid,tajls,'tajls(jm_budg,lm,ktajls)')
      endif

!      call def_rsf_tcons(fid,r4_on_disk)

      return
      end subroutine def_rsf_trdiag

      subroutine new_io_trdiag(fid,iaction)
!@sum  new_io_trdiag read/write tracer acc arrays from/to restart+acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite,iowrite_single
      use CHEM_COM, only : GCAIJL=>GCAIJL_out ! , &
!                            GCAIJLN=>GCAIJLN_loc, &
!                            GCAIJLS=>GCAIJLS_loc, &
!                            GCAIJN=>GCAIJN_loc, &
!                            GCAIJS=>GCAIJS_loc, &
!                            TAJLN, &
!                            TAJLS, &
!                            GCAIJ=>GCAIJ_out, &
!                            TAJL=>TAJL_out
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data,write_data,read_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite_single)     ! output to acc file
        call write_dist_data(grid,fid,'gcaijl',gcaijl)
        !call write_dist_data(grid,fid,'gcaij',gcaij)
        !call write_data(grid,fid,'tajl',tajl)
        !call write_src_dist_data(fid, .false.)
      case (iowrite)            ! output to restart file
        !call gather_zonal_trdiag
        !call write_dist_data(grid,fid,'gcaijln',gcaijln)
        !call write_dist_data(grid,fid,'gcaijls',gcaijls)
        !call write_dist_data(grid,fid,'gcaijs',gcaijs)
        !call write_dist_data(grid,fid,'gcaijn',gcaijn)
        !call write_data(grid,fid,'tajln',tajln)
        !call write_data(grid,fid,'tajls',tajls)
      case (ioread)            ! input from restart file
        !call read_dist_data(grid,fid,'gcaijln',gcaijln)
        !call read_dist_data(grid,fid,'gcaijls',gcaijls)
        !call read_dist_data(grid,fid,'gcaijs',gcaijs)
        !call read_dist_data(grid,fid,'gcaijn',gcaijn)
        !call read_data(grid,fid,'tajln',tajln)
        !call read_data(grid,fid,'tajls',tajls)
        !call scatter_zonal_trdiag
      end select

      !call new_io_tcons(fid,iaction)

      return
      end subroutine new_io_trdiag

      subroutine def_meta_trdiag(fid)
!@sum  def_meta_trdiag defines tracer metadata in acc files
!@auth M. Kelley
!@ver  beta
      use CHEM_COM
      use pario, only : defvar,write_attr
      use domain_decomp_atm, only : grid
      use cdl_mod, only : defvar_cdl
      implicit none
      integer :: fid         !@var fid file id

!       call write_attr(grid,fid,'gcaij','reduction','sum')
!       call write_attr(grid,fid,'gcaij','split_dim',3)
!       call defvar(grid,fid,ia_gcaij(1:kgcaij_out), 'ia_gcaij(kgcaij)')
!       call defvar(grid,fid,denom_gcaij(1:kgcaij_out), 'denom_gcaij(kgcaij)')
!       call defvar(grid,fid,scale_gcaij(1:kgcaij_out), 'scale_gcaij(kgcaij)')
!       call defvar(grid,fid,sname_gcaij(1:kgcaij_out),
!       'sname_gcaij(sname_strlen,kgcaij)')
!       call defvar(grid,fid,hemis_gcaij,'hemis_gcaij(one,shnhgm,kgcaij)',
!       r4_on_disk=.true.)
!       call write_attr(grid,fid,'hemis_gcaij','reduction','sum')
!       call defvar_cdl(grid,fid,cdl_gcaij, 'cdl_gcaij(cdl_strlen,kcdl_gcaij)')
! #ifdef CUBED_SPHERE
!       call defvar_cdl(grid,fid,cdl_gcaij_latlon,
!       'cdl_gcaij_latlon(cdl_strlen,kcdl_gcaij_latlon)')
! #endif

      call write_attr(grid,fid,'gcaijl','reduction','sum')
      call write_attr(grid,fid,'gcaijl','split_dim',4)
      call defvar(grid,fid,ia_gcaijl(1:kgcaijl_out), 'ia_gcaijl(kgcaijl)')
      call defvar(grid,fid,denom_gcaijl(1:kgcaijl_out), 'denom_gcaijl(kgcaijl)')
      call defvar(grid,fid,scale_gcaijl(1:kgcaijl_out), 'scale_gcaijl(kgcaijl)')
      call defvar(grid,fid,sname_gcaijl(1:kgcaijl_out), 'sname_gcaijl(sname_strlen,kgcaijl)')
      call defvar_cdl(grid,fid,cdl_gcaijl, 'cdl_gcaijl(cdl_strlen,kcdl_gcaijl)')
#ifdef CUBED_SPHERE
      call defvar_cdl(grid,fid,cdl_gcaijl_latlon, 'cdl_gcaijl_latlon(cdl_strlen,kcdl_gcaijl_latlon)')
#endif

!       call write_attr(grid,fid,'tajl','reduction','sum')
!       call write_attr(grid,fid,'tajl','split_dim',3)
!       call defvar(grid,fid,ia_tajl(1:ktajl_out), 'ia_tajl(ktajl)')
!       call defvar(grid,fid,denom_tajl(1:ktajl_out), 'denom_tajl(ktajl)')
!       call defvar(grid,fid,scale_tajl(1:ktajl_out), 'scale_tajl(ktajl)')
!       call defvar(grid,fid,sname_tajl(1:ktajl_out),
!       'sname_tajl(sname_strlen,ktajl)')
!       call defvar_cdl(grid,fid,cdl_tajl, 'cdl_tajl(cdl_strlen,kcdl_tajl)')
!       call defvar(grid,fid,hemis_tajl,'hemis_tajl(shnhgm,lm,ktajl)',
!       r4_on_disk=.true.)
!       call write_attr(grid,fid,'hemis_tajl','reduction','sum')
!       call defvar(grid,fid,vmean_tajl,
!       'vmean_tajl(jm_budg_plus3,one,ktajl)',r4_on_disk=.true.)
!       call write_attr(grid,fid,'vmean_tajl','reduction','sum')

      return
      end subroutine def_meta_trdiag

      subroutine write_meta_trdiag(fid)
!@sum  write_meta_trdiag write tracer accumulation metadata to file
!@auth M. Kelley
      use CHEM_COM
      use pario, only : write_dist_data,write_data
      use domain_decomp_atm, only : grid
      use cdl_mod, only : write_cdl
      implicit none
      integer :: fid         !@var fid file id

!       call write_data(grid,fid,'hemis_gcaij',hemis_gcaij)
!       call write_data(grid,fid,'ia_gcaij',ia_gcaij(1:kgcaij_out))
!       call write_data(grid,fid,'denom_gcaij',denom_gcaij(1:kgcaij_out))
!       call write_data(grid,fid,'scale_gcaij',scale_gcaij(1:kgcaij_out))
!       call write_data(grid,fid,'sname_gcaij',sname_gcaij(1:kgcaij_out))
!       call write_cdl(grid,fid,'cdl_gcaij',cdl_gcaij)
! #ifdef CUBED_SPHERE
!       call write_cdl(grid,fid,'cdl_gcaij_latlon',cdl_gcaij_latlon)
! #endif

      call write_data(grid,fid,'ia_gcaijl',ia_gcaijl(1:kgcaijl_out))
      call write_data(grid,fid,'denom_gcaijl',denom_gcaijl(1:kgcaijl_out))
      call write_data(grid,fid,'scale_gcaijl',scale_gcaijl(1:kgcaijl_out))
      call write_data(grid,fid,'sname_gcaijl',sname_gcaijl(1:kgcaijl_out))
      call write_cdl(grid,fid,'cdl_gcaijl',cdl_gcaijl)
#ifdef CUBED_SPHERE
      call write_cdl(grid,fid,'cdl_gcaijl_latlon',cdl_gcaijl_latlon)
#endif

!       call write_data(grid,fid,'hemis_tajl',hemis_tajl)
!       call write_data(grid,fid,'vmean_tajl',vmean_tajl)
!       call write_data(grid,fid,'ia_tajl',ia_tajl(1:ktajl_out))
!       call write_data(grid,fid,'denom_tajl',denom_tajl(1:ktajl_out))
!       call write_data(grid,fid,'scale_tajl',scale_tajl(1:ktajl_out))
!       call write_data(grid,fid,'sname_tajl',sname_tajl(1:ktajl_out))
!       call write_cdl(grid,fid,'cdl_tajl',cdl_tajl)

      return
      end subroutine write_meta_trdiag

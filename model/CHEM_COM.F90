module CHEM_COM

  !================================================================================================
  ! Module CHEM_COM is a module that holds diagnostic arrays for
  ! GEOS-Chem
  ! 
  ! Author: Lee T. Murray (lee.murray@rochester.edu)
  !===============================================================================================
  USE RESOLUTION,  ONLY : im, jm, lm
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

  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:)            :: TrID_to_SpcChmID
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:)            :: SpcChmID_to_TrID

  INTEGER :: ijlt_OHvmr
  
  !====================
  ! Diagnostics
  !====================

  INTEGER, PARAMETER :: kgcijl=500

  !@var SNAME_IJLT: Names of 3D tracer IJL diagnostics
  character(len=sname_strlen), dimension(kgcijl) :: sname_ijlt
  !@var DNAME_IJLT, DENOM_IJLT: Short names, indices of gcaijls denominators.
  !@+   Currently, dname is specified along with the standard metadata and
  !@+   the denom indices are looked up afterward.
  character(len=sname_strlen), dimension(kgcijl) :: dname_ijlt=''
  integer, dimension(kgcijl) :: denom_ijlt=0
  !@var LNAME_IJLT,UNITS_IJLT: descriptions/units of 3D tracer diagnostics
  character(len=lname_strlen), dimension(kgcijl) :: lname_ijlt = 'unused'
  character(len=units_strlen), dimension(kgcijl) :: units_ijlt
  !@var SCALE_IJLT: printout scaling factor for 3D tracer diagnostics
  REAL*8, dimension(kgcijl) :: scale_ijlt
  !@var IR_IJLT: range index of IJL diagnostics
  integer, dimension(kgcijl) :: ir_ijlt
  !@var IA_IJLT: accumulation index for IJL diagnostics
  integer, dimension(kgcijl) :: ia_ijlt
  !@var ijlt_power: power of 10 used for tracer IJL 3D diags
  INTEGER, DIMENSION(kgcijl) :: ijlt_power

  integer :: kgcijl_
  integer :: kgcijl_out ! actual number of qtys in gcijl_out
  real*8, dimension(:,:,:,:), allocatable :: gcijl_out
  integer, allocatable, dimension(:) ::ir_gcijl,ia_gcijl,denom_gcijl,ijlt_vmr
  character(len=lname_strlen),allocatable,dimension(:) ::lname_gcijl
  character(len=sname_strlen),allocatable,dimension(:) ::sname_gcijl
  character(len=units_strlen),allocatable,dimension(:) ::units_gcijl
  real*8, allocatable, dimension(:) :: scale_gcijl
  type(cdl_type) :: cdl_gcijl,cdl_gcijl_latlon

  !@var GCIJLN 3D tracer diagnostics (all tracers)
  real*8, allocatable, dimension(:,:,:,:) :: gcijln
  real*8, allocatable, dimension(:,:,:,:) :: gcijln_loc

  !@var GCIJLS 3D tracer diagnostics
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: GCIJLS
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: GCIJLS_loc

CONTAINS

end module CHEM_COM

!      subroutine write_src_dist_data(fid, def)
!      use pario, only: defvar, write_dist_data, write_data
!      use domain_decomp_atm, only : grid
!      use chem_com, only: taijln=>taijln_loc
!      !use chem_com, only: taijn=>taijn_loc
!      !use chem_com, only: tij_prec
!      use oldtracer_mod, only: src_dist_base, src_dist_index
!      use tracer_com, only: ntm, xyztr, ntm_sph, ntm_reg
!      implicit none
!      integer, intent(in) :: fid
!      logical, intent(in) :: def
!      real*8, dimension(:, :, :, :), allocatable, save :: tr2
!      real*8, dimension(:, :, :, :, :), allocatable, save :: tr3
!      integer, dimension(ntm) :: lst
!      integer :: i, n, ndist, nindex
!
!      if (def) then
!        ndist=0
!        nindex=0
!        do n=1, ntm
!          if (src_dist_index(n)==1) then
!            ndist=ndist+1
!            lst(n)=ndist
!          end if
!          nindex=max(nindex, src_dist_index(n))
!        end do
!        if (ndist==0) return
!        allocate(tr2(size(taijn, 1), size(taijn, 2), ndist, nindex))
!        allocate(tr3(size(taijln, 1), size(taijln, 2), size(taijln, 3),ndist, nindex))
!        do n=1, ntm
!          if (src_dist_index(n)==1) then
!            do i=n, ntm
!              if (src_dist_base(i)==src_dist_base(n)) then
!                tr2(:, :, lst(n), src_dist_index(i))=taijn(:, :, tij_prec, i)
!                tr3(:, :, :, lst(n), src_dist_index(i))=taijln(:, :, :, i)
!              end if
!            end do
!          end if
!        end do
!        call defvar(grid,fid,tr2,'src_dist2(dist_im,dist_jm,ndist,nbasis)')
!        call defvar(grid,fid,tr3,'src_dist3(dist_im,dist_jm,lm,ndist,nbasis)')
!        call defvar(grid,fid,xyztr,'src_dist_basis(nbasis,dist_im,dist_jm)')
!        call defvar(grid,fid,ntm_sph,'ntm_sph')
!        call defvar(grid,fid,ntm_reg,'ntm_reg')
!      else
!        if (allocated(tr2).and.allocated(tr3)) then
!          call write_dist_data(grid,fid,'src_dist2',tr2)
!          call write_dist_data(grid,fid,'src_dist3',tr3)
!          call write_data(grid,fid,'ntm_sph',ntm_sph)
!          call write_data(grid,fid,'ntm_reg',ntm_reg)
!          deallocate(tr2)
!          deallocate(tr3)
!        end if
!        if (allocated(xyztr)) &
!         call write_dist_data(grid,fid,'src_dist_basis',xyztr,jdim=3)
!      endif
!      end subroutine write_src_dist_data

subroutine def_rsf_gcdiag(fid,r4_on_disk)
  !@sum  def_rsf_gcdiag defines tracer diag array structure in restart+acc files
  !@auth M. Kelley
  !@ver  beta
  use chem_com, only : GCIJLN=>GCIJLN_loc
  use chem_com, only : GCIJLS=>GCIJLS_loc
  use chem_com, only : GCIJL=>GCIJL_out
  !     &     GCIJN=>GCIJN_loc,
  !     &     GCIJS=>GCIJS_loc,
  !     &     GCJLN,
  !     &     GCJLS,
  !     &     GCIJ=>GCIJ_out,
  !     &     GCJL=>GCJL_out
  use domain_decomp_atm, only : grid
  use pario, only : defvar
  implicit none  
  integer fid           !@var fid file id
  logical :: r4_on_disk !@var r4_on_disk if true, real*8 stored as real*4
  if(r4_on_disk) then ! acc file
     call defvar(grid,fid,gcijl,'gcijl(dist_im,dist_jm,lm,kgcijl)',r4_on_disk=.true.)
     !        call defvar(grid,fid,gcij,
     !     &       'gcij(dist_im,dist_jm,kgcij)',r4_on_disk=.true.)
     !        call defvar(grid,fid,gcjl,
     !     &       'gcjl(jm_budg,lm,ktajl)',r4_on_disk=.true.)
     !        call write_src_dist_data(fid, .true.)
  else
     call defvar(grid,fid,gcijln,'gcijln(dist_im,dist_jm,lm,ntm)')
     call defvar(grid,fid,gcijls,'gcijls(dist_im,dist_jm,lm,kgcijl)')
     !        call defvar(grid,fid,gcijs,'gcijs(dist_im,dist_jm,kgcijs)')
     !        call defvar(grid,fid,gcijn,'gcijn(dist_im,dist_jm,kgcij,ntm)')
     !        call defvar(grid,fid,gcjln,'gcjln(jm_budg,lm,ktajlx,ntm)')
     !        call defvar(grid,fid,gcjls,'gcjls(jm_budg,lm,ktajls)')
  endif

  ! call def_rsf_tcons(fid,r4_on_disk)

  return
end subroutine def_rsf_gcdiag

subroutine new_io_gcdiag(fid,iaction)
  !@sum  new_io_gcdiag read/write tracer acc arrays from/to restart+acc files
  !@auth M. Kelley
  !@ver  beta new_ prefix avoids name clash with the default version
  use model_com, only : ioread,iowrite,iowrite_single
  use chem_com, only : GCIJLN=>GCIJLN_loc
  use chem_com, only : GCIJLS=>GCIJLS_loc
  use chem_com, only : GCIJL=>GCIJL_out
  !     &     GCIJN=>GCIJN_loc,
  !     &     GCIJS=>GCIJS_loc,
  !     &     GCJLN,
  !     &     GCJLS,
  !     &     GCIJ=>GCIJ_out,
  !     &     GCJL=>GCJL_out
  use domain_decomp_atm, only : grid
  use pario, only : write_dist_data,read_dist_data,write_data,read_data
  implicit none
  integer fid   !@var fid unit number of read/write
  integer iaction !@var iaction flag for reading or writing to file
  select case (iaction)
  case (iowrite_single)     ! output to acc file
     call write_dist_data(grid,fid,'gcijl',gcijl)
     !call write_dist_data(grid,fid,'gcij',gcij)
     !call write_data(grid,fid,'gcjl',gcjl)
     !call write_src_dist_data(fid, .false.)
  case (iowrite)            ! output to restart file
     !call gather_zonal_gcdiag
     call write_dist_data(grid,fid,'gcijln',gcijln)
     call write_dist_data(grid,fid,'gcijls',gcijls)
     !call write_dist_data(grid,fid,'gcijs',gcijs)
     !call write_dist_data(grid,fid,'gcijn',gcijn)
     !call write_data(grid,fid,'gcjln',gcjln)
     !call write_data(grid,fid,'gcjls',gcjls)
  case (ioread)            ! input from restart file
     call read_dist_data(grid,fid,'gcijln',gcijln)
     call read_dist_data(grid,fid,'gcijls',gcijls)
     !call read_dist_data(grid,fid,'gcijs',gcijs)
     !call read_dist_data(grid,fid,'gcijn',gcijn)
     !call read_data(grid,fid,'gcjln',gcjln)
     !call read_data(grid,fid,'gcjls',gcjls)
     !call scatter_zonal_gcdiag
  end select

  !call new_io_tcons(fid,iaction)

  return
end subroutine new_io_gcdiag

subroutine def_meta_gcdiag(fid)
  !@sum  def_meta_gcdiag defines tracer metadata in acc files
  !@auth M. Kelley
  !@ver  beta
  use chem_com
  use pario, only : defvar,write_attr
  use domain_decomp_atm, only : grid
  use cdl_mod, only : defvar_cdl
  implicit none
  integer :: fid         !@var fid file id
  
  !      call write_attr(grid,fid,'gcij','reduction','sum')
  !      call write_attr(grid,fid,'gcij','split_dim',3)
  !      call defvar(grid,fid,ia_gcij(1:kgcij_out),
  !     &     'ia_gcij(kgcij)')
  !      call defvar(grid,fid,denom_gcij(1:kgcij_out),
  !     &     'denom_gcij(kgcij)')
  !      call defvar(grid,fid,scale_gcij(1:kgcij_out),
  !     &     'scale_gcij(kgcij)')
  !      call defvar(grid,fid,sname_gcij(1:kgcij_out),
  !     &     'sname_gcij(sname_strlen,kgcij)')
  !      call defvar(grid,fid,hemis_gcij,'hemis_gcij(one,shnhgm,kgcij)',
  !     &     r4_on_disk=.true.)
  !      call write_attr(grid,fid,'hemis_gcij','reduction','sum')
  !      call defvar_cdl(grid,fid,cdl_gcij,
  !     &     'cdl_gcij(cdl_strlen,kcdl_gcij)')
  !#ifdef CUBED_SPHERE
  !      call defvar_cdl(grid,fid,cdl_gcij_latlon,
  !     &     'cdl_gcij_latlon(cdl_strlen,kcdl_gcij_latlon)')
  !#endif

  call write_attr(grid,fid,'gcijl','reduction','sum')
  call write_attr(grid,fid,'gcijl','split_dim',4)
  call defvar(grid,fid,ia_gcijl(1:kgcijl_out),'ia_gcijl(kgcijl)')
  call defvar(grid,fid,denom_gcijl(1:kgcijl_out),'denom_gcijl(kgcijl)')
  call defvar(grid,fid,scale_gcijl(1:kgcijl_out),'scale_gcijl(kgcijl)')
  call defvar(grid,fid,sname_gcijl(1:kgcijl_out),'sname_gcijl(sname_strlen,kgcijl)')
  call defvar_cdl(grid,fid,cdl_gcijl,'cdl_gcijl(cdl_strlen,kcdl_gcijl)')
#ifdef CUBED_SPHERE
  call defvar_cdl(grid,fid,cdl_gcijl_latlon,'cdl_gcijl_latlon(cdl_strlen,kcdl_gcijl_latlon)')
#endif

  !     call write_attr(grid,fid,'gcjl','reduction','sum')
  !     call write_attr(grid,fid,'gcjl','split_dim',3)
  !     call defvar(grid,fid,ia_gcjl(1:ktajl_out),
  !    &     'ia_gcjl(ktajl)')
  !     call defvar(grid,fid,denom_gcjl(1:ktajl_out),
  !    &     'denom_gcjl(ktajl)')
  !     call defvar(grid,fid,scale_gcjl(1:ktajl_out),
  !    &     'scale_gcjl(ktajl)')
  !     call defvar(grid,fid,sname_gcjl(1:ktajl_out),
  !    &     'sname_gcjl(sname_strlen,ktajl)')
  !     call defvar_cdl(grid,fid,cdl_gcjl,
  !    &     'cdl_gcjl(cdl_strlen,kcdl_gcjl)')
  !     call defvar(grid,fid,hemis_gcjl,'hemis_gcjl(shnhgm,lm,ktajl)',
  !    &     r4_on_disk=.true.)
  !     call write_attr(grid,fid,'hemis_gcjl','reduction','sum')
  !     call defvar(grid,fid,vmean_gcjl,
  !    &     'vmean_gcjl(jm_budg_plus3,one,ktajl)',r4_on_disk=.true.)
  !     call write_attr(grid,fid,'vmean_gcjl','reduction','sum')

  !      call def_meta_tcons(fid)

  return
end subroutine def_meta_gcdiag

subroutine write_meta_gcdiag(fid)
  !@sum  write_meta_gcdiag write tracer accumulation metadata to file
  !@auth M. Kelley
  use chem_com
  use pario, only : write_dist_data,write_data
  use domain_decomp_atm, only : grid
  use cdl_mod, only : write_cdl
  implicit none
  integer :: fid         !@var fid file id

  !      call write_data(grid,fid,'hemis_gcij',hemis_gcij)
  !      call write_data(grid,fid,'ia_gcij',ia_gcij(1:kgcij_out))
  !      call write_data(grid,fid,'denom_gcij',denom_gcij(1:kgcij_out))
  !      call write_data(grid,fid,'scale_gcij',scale_gcij(1:kgcij_out))
  !      call write_data(grid,fid,'sname_gcij',sname_gcij(1:kgcij_out))
  !      call write_cdl(grid,fid,'cdl_gcij',cdl_gcij)
  !#ifdef CUBED_SPHERE
  !      call write_cdl(grid,fid,'cdl_gcij_latlon',cdl_gcij_latlon)
  !#endif

  call write_data(grid,fid,'ia_gcijl',ia_gcijl(1:kgcijl_out))
  call write_data(grid,fid,'denom_gcijl',denom_gcijl(1:kgcijl_out))
  call write_data(grid,fid,'scale_gcijl',scale_gcijl(1:kgcijl_out))
  call write_data(grid,fid,'sname_gcijl',sname_gcijl(1:kgcijl_out))
  call write_cdl(grid,fid,'cdl_gcijl',cdl_gcijl)
#ifdef CUBED_SPHERE
  call write_cdl(grid,fid,'cdl_gcijl_latlon',cdl_gcijl_latlon)
#endif

  !      call write_data(grid,fid,'hemis_gcjl',hemis_gcjl)
  !      call write_data(grid,fid,'vmean_gcjl',vmean_gcjl)
  !      call write_data(grid,fid,'ia_gcjl',ia_gcjl(1:ktajl_out))
  !      call write_data(grid,fid,'denom_gcjl',denom_gcjl(1:ktajl_out))
  !      call write_data(grid,fid,'scale_gcjl',scale_gcjl(1:ktajl_out))
  !      call write_data(grid,fid,'sname_gcjl',sname_gcjl(1:ktajl_out))
  !      call write_cdl(grid,fid,'cdl_gcjl',cdl_gcjl)

  !      call write_meta_tcons(fid)

  return
end subroutine write_meta_gcdiag

!      subroutine def_rsf_tcons(fid,r4_on_disk)
!!@sum  def_rsf_tcons defines tracer diag array structure in restart+acc files
!!@auth M. Kelley
!!@ver  beta
!      use chem_com, only :
!     &     TCONSRV,
!     &     TCONSRV_out
!      use domain_decomp_atm, only : grid
!      use pario, only : defvar
!      implicit none
!      integer fid           !@var fid file id
!      logical :: r4_on_disk !@var r4_on_disk if true, real*8 stored as real*4
!      if(r4_on_disk) then ! acc file
!        call defvar(grid,fid,tconsrv_out,
!     &       'tconsrv(jm_budg,ktcon)',r4_on_disk=.true.)
!      else
!        call defvar(grid,fid,tconsrv,
!     &       'tconsrv(jm_budg,ktcon,ntmxcon)',r4_on_disk=r4_on_disk)
!      endif
!      return
!      end subroutine def_rsf_tcons
!
!      subroutine new_io_tcons(fid,iaction)
!!@sum  new_io_tcons read/write tconsrv arrays from/to restart+acc files
!!@auth M. Kelley
!!@ver  beta new_ prefix avoids name clash with the default version
!      use model_com, only : ioread,iowrite,iowrite_single
!      use chem_com, only :
!     &     TCONSRV,
!     &     TCONSRV_out
!      use domain_decomp_atm, only : grid
!      use pario, only : write_dist_data,read_dist_data,
!     &     write_data,read_data
!      implicit none
!      integer fid   !@var fid unit number of read/write
!      integer iaction !@var iaction flag for reading or writing to file
!      select case (iaction)
!      case (iowrite_single)     ! output to acc file
!        call write_data(grid,fid,'tconsrv',tconsrv_out)
!      case (iowrite)            ! output to restart file
!        call gather_zonal_tcons
!        call write_data(grid,fid,'tconsrv',tconsrv)
!      case (ioread)            ! input from restart file
!        call read_data(grid,fid,'tconsrv',tconsrv)
!        call scatter_zonal_tcons
!      end select
!      return
!      end subroutine new_io_tcons
!
!      subroutine def_meta_tcons(fid)
!!@sum  def_meta_gcdiag defines tconsrv metadata in acc files
!!@auth M. Kelley
!!@ver  beta
!      use chem_com
!      use pario, only : defvar,write_attr
!      use domain_decomp_atm, only : grid
!      use cdl_mod, only : defvar_cdl
!      implicit none
!      integer :: fid         !@var fid file id
!
!      call write_attr(grid,fid,'tconsrv','reduction','sum')
!      call write_attr(grid,fid,'tconsrv','split_dim',2)
!      call defvar(grid,fid,hemis_tconsrv,'hemis_tconsrv(shnhgm,ktcon)',
!     &     r4_on_disk=.true.)
!      call write_attr(grid,fid,'hemis_tconsrv','reduction','sum')
!      call defvar(grid,fid,ia_tcon_out,'ia_tconsrv(ktcon)')
!      call defvar(grid,fid,scale_tcon_out,'scale_tconsrv(ktcon)')
!      call defvar(grid,fid,sname_tconsrv_out,
!     &     'sname_tconsrv(sname_strlen,ktcon)')
!      call defvar_cdl(grid,fid,cdl_tconsrv,
!     &     'cdl_tconsrv(cdl_strlen,kcdl_tconsrv)')
!
!      return
!      end subroutine def_meta_tcons
!
!      subroutine write_meta_tcons(fid)
!!@sum  write_meta_tcons write tconsrv accumulation metadata to file
!!@auth M. Kelley
!      use chem_com
!      use pario, only : write_dist_data,write_data
!      use domain_decomp_atm, only : grid
!      use cdl_mod, only : write_cdl
!      implicit none
!      integer :: fid         !@var fid file id
!
!      call write_data(grid,fid,'hemis_tconsrv',hemis_tconsrv)
!      call write_data(grid,fid,'ia_tconsrv',ia_tcon_out)
!      call write_data(grid,fid,'scale_tconsrv',scale_tcon_out)
!      call write_data(grid,fid,'sname_tconsrv',sname_tconsrv_out)
!      call write_cdl(grid,fid,'cdl_tconsrv',cdl_tconsrv)
!
!      return
!      end subroutine write_meta_tcons

SUBROUTINE ALLOC_CHEM_COM
  USE DIAG_COM, only : jm_budg
  USE CHEM_COM
  USE DOMAIN_DECOMP_ATM, only : getDomainBounds, AM_I_ROOT, GRID
  use diag_zonal, only : get_alloc_bounds
  use fluxes, only : atmocn
  implicit none
  INTEGER :: J_0H,J_1H, I_0H,I_1H
  INTEGER :: status
  integer :: j_0budg,j_1budg
  integer :: img, jmg

  call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
  I_0H=GRID%I_STRT_HALO
  I_1H=GRID%I_STOP_HALO
  
  call get_alloc_bounds(grid,j_strt_budg=j_0budg,j_stop_budg=j_1budg)

  ALLOCATE ( GCIJLN_loc(I_0H:I_1H,J_0H:J_1H,LM,ntm), stat=status )
  ALLOCATE ( GCIJLS_loc(I_0H:I_1H,J_0H:J_1H,LM,kgcijl), stat=status)
  if(am_i_root()) then
     img = IM
     jmg = JM
  else
     img = 1
     jmg = 1
  end if
  ALLOCATE ( GCIJLN(img,jmg,LM,ntm), stat=status )
  ALLOCATE ( GCIJLS(img,jmg,LM,kgcijl), stat=status )

  kgcijl_ = (ntm + kgcijl) !*3/2  ! make 50% larger for denoms and extra specials
  kgcijl_out = kgcijl_
  
  allocate(ir_gcijl(kgcijl_))
  ir_gcijl = 0
  allocate(ia_gcijl(kgcijl_))
  ia_gcijl = 0
  allocate(denom_gcijl(kgcijl_))
  denom_gcijl = 0
  allocate(lname_gcijl(kgcijl_))
  allocate(sname_gcijl(kgcijl_))
  allocate(units_gcijl(kgcijl_))
  allocate(scale_gcijl(kgcijl_))

  ALLOCATE ( GCIJL_out( I_0H:I_1H,J_0H:J_1H,LM,kgcijl_),stat=status)

  RETURN
END SUBROUTINE ALLOC_CHEM_COM

!      subroutine reset_tcons
!        USE CHEM_COM, only: TCONSRV_loc, TCONSRV
!        USE DOMAIN_DECOMP_ATM, only : am_i_root
!        implicit none
!        TCONSRV_loc=0.
!        if(am_i_root()) TCONSRV=0.
!        return
!      end subroutine reset_tcons
!
!      subroutine gather_zonal_tcons
!        USE DIAG_COM, only : ia_inst,jm_budg
!        USE CHEM_COM, only: TCONSRV_loc, TCONSRV, ia_tcon,ntmxcon,ktcon
!        USE DOMAIN_DECOMP_ATM, only : GRID,am_i_root,sumxpe
!        implicit none
!        real*8, dimension(:,:,:), allocatable :: tconsrv_sv
!        integer :: n,k
!
!        ! TCONSRV is a mixture of accumulations and instantaneous values, hence the
!        ! complicated logic
!        ! set inst qtys to zero in the "global" array before summing over PEs
!        if(am_i_root()) then
!           allocate(tconsrv_sv(jm_budg,ktcon,ntmxcon)); tconsrv_sv=tconsrv
!           do n=1,ntmxcon
!              do k=1,ktcon
!                 if(ia_tcon(k,n).eq.ia_inst) tconsrv(:,k,n)=0.
!              enddo
!           enddo
!        endif
!
!        call sumxpe (TCONSRV_loc,TCONSRV, increment=.true. )
!
!        do n=1,ntmxcon ! keep instantaneous values
!           do k=1,ktcon
!              if(ia_tcon(k,n).ne.ia_inst) tconsrv_loc(:,k,n)=0. 
!           enddo
!        enddo
!
!        if(am_i_root()) then
!           do n=1,ntmxcon
!              do k=1,ktcon
!                 if(ia_tcon(k,n).eq.ia_inst .and. all(tconsrv(:,k,n)==0.)) &
!                          tconsrv(:,k,n) = tconsrv_sv(:,k,n)
!              enddo
!           enddo
!           deallocate(tconsrv_sv)
!        endif
!
!        return
!      end subroutine gather_zonal_tcons
!
!      subroutine scatter_zonal_tcons
!        USE CHEM_COM, only: TCONSRV_loc, TCONSRV
!        USE DOMAIN_DECOMP_ATM, only : GRID
!        USE DIAG_ZONAL, only : unpack_lc
!        implicit none
!        !call unpack_lc   (grid, TCONSRV, TCONSRV_loc)
!        tconsrv_loc = 0
!        return
!      end subroutine scatter_zonal_tcons

subroutine reset_gcdiag
  USE CHEM_COM, only: GCIJLN_loc, GCIJLS_loc ! GCIJN_loc, GCIJS_loc, GCJLN_loc, GCJLS_loc, GCJLN, GCJLS
  USE DOMAIN_DECOMP_ATM, only : am_i_root
  implicit none

  !TAJLN_loc=0.
  !TAJLS_loc=0. 
  GCIJLN_loc=0.
  GCIJLS_loc=0.
  !TAIJN_loc=0.
  !TAIJS_loc=0.
  !if(am_i_root()) then
  !   GCJLN=0. ; GCJLS=0. 
  !endif
  !call reset_tcons

  return
end subroutine reset_gcdiag

subroutine gather_gcdiag
  USE CHEM_COM, only : GCIJLN, GCIJLN_loc, GCIJLS, GCIJLS_loc !, GCIJN, GCIJN_loc,GCIJS, GCIJS_loc
  USE DOMAIN_DECOMP_ATM, only : grid
  USE DOMAIN_DECOMP_1D, ONLY : PACK_DATA
  implicit none

  CALL PACK_DATA (GRID, GCIJLN_loc, GCIJLN)
  CALL PACK_DATA (GRID, GCIJLS_loc, GCIJLS)
  !CALL PACK_DATA (GRID, GCIJN_loc , GCIJN)
  !CALL PACK_DATA (GRID, GCIJS_loc , GCIJS)
  !call gather_zonal_gcdiag
  !call gather_zonal_tcons

  return
end subroutine gather_gcdiag

!subroutine gather_zonal_gcdiag
!  USE CHEM_COM, only : GCJLN , GCJLN_loc, GCJLS, GCJLS_loc
!  USE DOMAIN_DECOMP_ATM, ONLY : GRID
!  USE DIAG_ZONAL, only : pack_lc
!  implicit none
!  call pack_lc   (grid, GCJLN_loc,  GCJLN )
!  call pack_lc   (grid, GCJLS_loc,  GCJLS )
!  return
!end subroutine gather_zonal_gcdiag

subroutine scatter_gcdiag
  USE CHEM_COM, only : GCIJLN, GCIJLN_loc, GCIJLS, GCIJLS_loc !, GCIJN, GCIJN_loc, GCIJS, GCIJS_loc
  USE DOMAIN_DECOMP_ATM, only : grid
  USE DOMAIN_DECOMP_1D, ONLY : UNPACK_DATA
  implicit none
  CALL UNPACK_DATA (GRID, GCIJLN, GCIJLN_loc)
  CALL UNPACK_DATA (GRID, GCIJLS, GCIJLS_loc)
  !CALL UNPACK_DATA (GRID, GCIJN , GCIJN_loc)
  !CALL UNPACK_DATA (GRID, GCIJS , GCIJS_loc)
  !call scatter_zonal_gcdiag
  !call scatter_zonal_tcons
  return
end subroutine scatter_gcdiag

!subroutine scatter_zonal_gcdiag
!  USE CHEM_COM, only : GCJLN , GCJLN_loc, GCJLS, GCJLS_loc
!  USE DOMAIN_DECOMP_ATM, ONLY : GRID
!  USE DIAG_ZONAL, only : unpack_lc
!  implicit none
!  call unpack_lc (grid, GCJLN,  GCJLN_loc )
!  call unpack_lc (grid, GCJLS,  GCJLS_loc )
!  return
!end subroutine scatter_zonal_gcdiag

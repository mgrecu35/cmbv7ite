      !This file contains various functions to read and interpolate parameters from look-up tables.
      !These are called internally and should not have to be modified unless the LUT dimensions
      !or file locations are modified. 
      
      !Last modified: Joe Munchak 7/9/2014
!Define LUT type

module LUT_def

  type LUTType
    integer :: nch
    
!     real :: tervels(60)
! 
    integer :: gas_nP, gas_nT, gas_nV
    real, dimension(:), allocatable :: gas_press, gas_temp, gas_wv
    real, dimension(:,:,:,:), allocatable :: gas_kext

    
    integer :: water_emis_nT, water_emis_nW, water_emis_nA
    real, dimension(:), allocatable :: water_emis_temp, water_emis_wind
    real, dimension(:,:), allocatable :: water_emis_angle
    real, dimension(:,:,:,:,:), allocatable :: water_emis
    real, dimension(:,:,:,:,:,:), allocatable :: water_eharm
    real, dimension(:,:,:), allocatable :: water_ebar
    
	integer :: clw_nT
    real, dimension(:), allocatable :: clw_temp
    real, dimension(:,:), allocatable :: clw_kext
!     
!     integer :: clwPR_nT
!     real, dimension(:), allocatable :: clwPR_temp
!     real, dimension(:), allocatable :: clwPR_kext
    
    integer :: water_sigma0_nA, water_sigma0_nW
    real, dimension(:), allocatable :: water_sigma0_wind
    real, dimension(:,:), allocatable :: water_sigma0_Ku, water_sigma0_Ku_std
    real, dimension(:,:), allocatable :: water_sigma0_Ka, water_sigma0_Ka_std
    real, dimension(:,:), allocatable :: water_sigma0_corr, water_sigma0_diff, water_sigma0_dstd
    real, dimension(:,:,:), allocatable :: water_sigma0_Ku_Harm, water_sigma0_Ka_harm
    
    integer :: land_nclass, land_class_nx, land_class_ny, land_nch, land_nrayKu, land_nrayKa
    integer(1), dimension(:,:), allocatable :: land_class_map, land_class_map_bare, land_class_map_snow
    real, dimension(:,:), allocatable :: land_class_emis_mean
    real, dimension(:,:), allocatable :: land_class_sigma0Ku_mean, land_class_sigma0Ku_std
    real, dimension(:,:), allocatable :: land_class_sigma0Ka_mean, land_class_sigma0Ka_std
    real, dimension(:,:,:), allocatable :: land_class_emis_eofs_LF !low frequency only EOFs - to be used when TPW > 10 mm
    real, dimension(:,:,:), allocatable :: land_class_emis_eofs_AF !all frequency EOFS - to be used when TPW <= 10 mm
    
    
    integer :: sfc_nclass, sfc_class_nx, sfc_class_ny, sfc_nch
    integer, dimension(:,:), allocatable :: bare_class_map, snow_class_map, seaice_class_map

    integer :: emis_map_nx, emis_map_ny, emis_map_nch
    real, dimension(:,:,:), allocatable :: emis_map_bare, emis_map_snow, emis_map_seaice, emis_map_all
    real, dimension(:,:,:), allocatable :: estd_map_bare, estd_map_snow, estd_map_seaice, estd_map_all
    
    !class emissivity correlation matrices
    real, dimension(:,:,:), allocatable :: emis_class_corr

    
    !integer :: s0Ku_map_nx, s0Ku_map_ny, s0Ku_map_nray
    !real, dimension (:,:,:), allocatable :: s0Ku_map_bare, s0Ku_map_snow, s0Ku_map_all
    !real, dimension (:,:,:), allocatable :: s0Ku_std_bare, s0Ku_std_snow, s0Ku_std_all
    
    !integer :: s0Ka_map_nx, s0Ka_map_ny, s0Ka_map_nray
    !real, dimension (:,:,:), allocatable :: s0Ka_map_bare, s0Ka_map_snow, s0Ka_map_all
    !real, dimension (:,:,:), allocatable :: s0Ka_std_bare, s0Ka_std_snow, s0Ka_std_all
    
    !Atmosphere profile EOFs
	integer :: nprof, prof_nts, prof_nst, prof_nlev
    integer :: prof_mints, prof_minst, prof_minlev
    integer, dimension(:,:,:), allocatable :: prof_index
    real, dimension(:,:), allocatable :: prof_W, prof_avg, prof_std, pc_std
    real, dimension(:,:,:), allocatable :: prof_eofs
    real, dimension(:), allocatable :: prof_tpw_mean, prof_tpw_rms, prof_tskin_mean, prof_tskin_rms
    
    !error covariance matrix from EOF representation of temp/wv profiles
    real, dimension(:,:,:), allocatable :: eof_covmat_land
    real, dimension(:,:,:), allocatable :: eof_covmat_water
    
    !RSS LUTs
    real :: acoef(5,2,6), bcoef(5,2,4,6), rcoef(3,0:2,6), scatterm(91,50,26,13,2)

    
  end type LUTType
  
  Type(LUTType) :: LUT
end module LUT_def

subroutine read_err_eofs(dirname)

  use LUT_def
  use profile_def
  implicit none
  include 'parametersMERRA.inc'
  
  character*90 :: dirname, fname
  integer :: i, j
  
  !there are only 2 model error profile eofs, one for land and another for ocean
  LUT%prof_nts = 1
  LUT%prof_nst = 2
  LUT%prof_nlev= 1
  LUT%nprof = 2
  allocate(LUT%prof_index(LUT%prof_nts,LUT%prof_nst,LUT%prof_nlev))
  allocate(LUT%prof_W(LUT%nprof,1+MERRA_NLEV*2),LUT%pc_std(LUT%nprof,1+MERRA_NLEV*2), LUT%prof_std(LUT%nprof,1+MERRA_NLEV*2))
  allocate(LUT%prof_eofs(LUT%nprof,1+MERRA_NLEV*2,1+MERRA_NLEV*2))
  
  !read in error EOFS
  i=1 !ocean
  write(fname, '(A,A41)') trim(dirname), 'raob_err_eofs.ocean.wv_wgt2.reordered.bin'
  open(unit=1,file=fname,form='unformatted',status='old',recl=((MERRA_NLEV*2+1)),access='direct')
  read(1,rec=1) LUT%prof_W(i,:)
  read(1,rec=2) LUT%prof_std(i,:)
  read(1,rec=3) LUT%pc_std(i,:)
!   print '(10F6.3)', LUT%prof_W(i,:)
!   print '(10F6.3)', LUT%prof_std(i,:)
!   print '(10F6.3)', LUT%pc_std(i,:)
!   stop
  do j=1,MERRA_NLEV*2+1
    read(1,rec=3+j) LUT%prof_eofs(i,j,:)
    !print*, j, LUT%prof_eofs(i,j,:)
  end do
  close(1)
  i=2 !land
  write(fname, '(A,A40)') trim(dirname), 'raob_err_eofs.land.wv_wgt2.reordered.bin'
  open(unit=1,file=fname,form='unformatted',status='old',recl=((MERRA_NLEV*2+1)),access='direct')
  read(1,rec=1) LUT%prof_W(i,:)
  read(1,rec=2) LUT%prof_std(i,:)
  read(1,rec=3) LUT%pc_std(i,:)
  !print '(10F6.3)', LUT%prof_W(i,:)
  !print '(10F6.3)', LUT%prof_std(i,:)
  !print '(10F6.3)', LUT%pc_std(i,:)
  do j=1,MERRA_NLEV*2+1
    read(1,rec=3+j) LUT%prof_eofs(i,j,:)
    !print*, j, LUT%prof_eofs(i,j,:)
  end do
  close(1)
  
  allocate(LUT%eof_covmat_land(MERRA_NLEV*2+1,13,13))
  allocate(LUT%eof_covmat_water(MERRA_NLEV*2+1,13,13))
  !read in associated GMI error covariances
  i=1 !ocean
  write(fname, '(A,A43)') trim(dirname), 'err_covmat_raob.ocean.wv_wgt2.reordered.dat'
  open(unit=1,file=fname,form='unformatted',status='old',recl=(((MERRA_NLEV*2+1)*13*13)),access='direct')
  read(1,rec=1) LUT%eof_covmat_water
  close(1)
  i=2 !land
  write(fname, '(A,A42)') trim(dirname), 'err_covmat_raob.land.wv_wgt2.reordered.dat'
  open(unit=1,file=fname,form='unformatted',status='old',recl=(((MERRA_NLEV*2+1)*13*13)),access='direct')
  read(1,rec=1) LUT%eof_covmat_land
  close(1)
  
  !print '(13F5.2)', LUT%eof_covmat_land(5,:,:)
  !stop
  
!   do i=1,MERRA_NLEV*2+1
!     print '(i5,13F8.5)', i, sqrt(LUT%eof_covmat_land(i,1,1)), &
!                          sqrt(LUT%eof_covmat_land(i,2,2)), &
!                          sqrt(LUT%eof_covmat_land(i,3,3)), &
!                          sqrt(LUT%eof_covmat_land(i,4,4)), &
!                          sqrt(LUT%eof_covmat_land(i,5,5)), &
!                          sqrt(LUT%eof_covmat_land(i,6,6)), &
!                          sqrt(LUT%eof_covmat_land(i,7,7)), &
!                          sqrt(LUT%eof_covmat_land(i,8,8)), &
!                          sqrt(LUT%eof_covmat_land(i,9,9)), &
!                          sqrt(LUT%eof_covmat_land(i,10,10)), &
!                          sqrt(LUT%eof_covmat_land(i,11,11)), &
!                          sqrt(LUT%eof_covmat_land(i,12,12)), &
!                          sqrt(LUT%eof_covmat_land(i,13,13))
! 
!   end do
!   stop
  !copy land to ocean for now
!   LUT%prof_W(1,:) = LUT%prof_W(2,:)
!   LUT%prof_std(1,:) = LUT%prof_std(2,:)
!   LUT%pc_std(1,:) = LUT%pc_std(2,:)
!   LUT%prof_eofs(1,:,:) = LUT%prof_eofs(2,:,:)
!   LUT%eof_covmat_water = LUT%eof_covmat_land
  
end subroutine read_err_eofs

subroutine read_LUTclw(filename)

use LUT_def
implicit none
      
integer :: i,j,k,l,ch
character*90 :: filename

!open LUT file
open(unit = 1, file = filename,status = 'old',form='unformatted')
     
!Check that LUT dimensions are what we think they are
read(1) LUT%nch,LUT%clw_nT
! print*, LUT%nch,LUT%clw_nT
allocate(LUT%clw_temp(LUT%clw_nT))
allocate(LUT%clw_kext(LUT%nch,LUT%clw_nT))

!Read LUT
do i = 1, LUT%nch
  do k = 1, LUT%clw_nT
    read(1) ch,LUT%clw_temp(k),LUT%clw_kext(i,k)
  end do
end do
close(1)

! print '(15F8.2)', LUT%clw_temp
! print '(10F8.2)', LUT%clw_kext
! stop
      
end      

subroutine intplte_gas(nf,Tavg,press,wv,kext)
!nf- freq number in LUT
!Tavg - temp in K
!Press - pressure in mb
!wv - vapor pressure in mb
!kext - extinction coef in km^-1
use LUT_def
implicit   none

real       Tavg, T, press,P, wv,V, dT, dP, dV
integer    nf, iT, iP, iV
real	   weights(2,2,2)
real	   kext

      
! print*, nf, Tavg, press, wv
P = press
if ( P .lt. LUT%gas_press(1) ) P = LUT%gas_press(1)
if ( P .gt. LUT%gas_press(LUT%gas_nP) ) P = LUT%gas_press(LUT%gas_nP)
T = Tavg
! print*, T
if ( T .lt. LUT%gas_temp(1) ) then
  !print*, Tavg
  T = LUT%gas_temp(1)
  !print*, T
endif
if ( T .gt. LUT%gas_temp(LUT%gas_nT) ) then
  !print*, Tavg
  T = LUT%gas_temp(LUT%gas_nT)
  !print*, T
endif
! print*, T
V = wv
if ( V .lt. LUT%gas_wv(1) ) V = LUT%gas_wv(1)
if ( V .gt. LUT%gas_wv(LUT%gas_nV) ) V = LUT%gas_wv(LUT%gas_nV)

!print*, Tavg,T

!print*, wv, V
! stop
!--Simple bilinear interpolation

call get_indx(P, LUT%gas_nP, LUT%gas_press, 117, iP)
call get_indx(T, LUT%gas_nT, LUT%gas_temp, 117, iT)
call get_indx(V, LUT%gas_nV, LUT%gas_wv, 117, iV)
! print*, iP, LUT%gas_press(iP)
! print*, iT, LUT%gas_temp(iT)
! print*, iV, LUT%gas_wv(iV)
dP = (P - LUT%gas_press(iP-1))/(LUT%gas_press(iP) - LUT%gas_press(iP-1))
dT = (T - LUT%gas_temp(iT-1))/(LUT%gas_temp(iT) - LUT%gas_temp(iT-1))
dV = (V - LUT%gas_wv(iV-1))/(LUT%gas_wv(iV) - LUT%gas_wv(iV-1))

! print*, dT, dP, dV

weights(1,1,1) = (1.-dP)*(1.-dT)*(1.-dV)
weights(1,1,2) = (1.-dP)*(1.-dT)*(   dV)
weights(1,2,1) = (1.-dP)*(   dT)*(1.-dV)
weights(1,2,2) = (1.-dP)*(   dT)*(   dV)
weights(2,1,1) = (   dP)*(1.-dT)*(1.-dV)
weights(2,1,2) = (   dP)*(1.-dT)*(   dV)
weights(2,2,1) = (   dP)*(   dT)*(1.-dV)
weights(2,2,2) = (   dP)*(   dT)*(   dV)

kext=sum(weights*LUT%gas_kext(nf,iP-1:iP,iT-1:iT,iV-1:iV))
      
return 
end

subroutine intplte_clw(nf,temp,kext)

use LUT_def
implicit   none

real       temp, T, dT
integer    nf, iT
real	   weights(2)
real       kext

      

T = temp
if ( T .lt. LUT%clw_temp(1) ) T = LUT%clw_temp(1)
if ( T .gt. LUT%clw_temp(LUT%clw_nT) ) T = LUT%clw_temp(LUT%clw_nT)

! print*, temp,T
! stop
!--Simple bilinear interpolation

call get_indx(T, LUT%clw_nT, LUT%clw_temp, 117, iT)

dT = (T - LUT%clw_temp(iT-1))/(LUT%clw_temp(iT) - LUT%clw_temp(iT-1))
! print*, dT, dW
! stop
weights(1) = (1.-dT)
weights(2) = (   dT)

kext=sum(weights*LUT%clw_kext(nf,iT-1:iT))

      
return 
end

subroutine read_RSS_tables(dirname)
  
  use LUT_def
  
  implicit none
  
  character*90, intent(in) :: dirname
  character(len=200) :: file_coeffs_wind_isotropic
  ! [MW 2012] Tables 3 + 4
  character(len=200) :: file_coeffs_wind_direction 
  ! [MW 2012] Tables 8, 9, 10, 11, 12
  character(len=200) :: file_coeffs_sctterm 
  ! [MW 2012] eq. (15) for fast computation
  character(len=200) :: file_em0_ref_freq_sst 
  
  write(file_coeffs_wind_isotropic, '(A,A)') trim(dirname), 'finetune_emiss_wind.dat'
  open(unit=3,file=file_coeffs_wind_isotropic,form='binary',status='old',action='read')
  read(3) LUT%acoef
  close(3)
  
  write(file_coeffs_wind_direction, '(A,A)') trim(dirname), 'fit_emiss_phir_wind.dat'
  open(unit=3,file=file_coeffs_wind_direction,form='binary',status='old',action='read')
  read(3) LUT%bcoef
  close(3) 
  
  write(file_em0_ref_freq_sst, '(A,A)') trim(dirname), 'em0_ref_freq_sst.dat'
  open(unit=3,file=file_em0_ref_freq_sst,status='old',form='binary',action='read')
  read(3) LUT%rcoef
  close(3)
  
  write(file_coeffs_sctterm, '(A,A)') trim(dirname), 'mk_scatterm_table_all.dat'
  open(unit=3,file=file_coeffs_sctterm,status='old',form='binary',action='read')
  read(3) LUT%scatterm
  close(3)
  !stop
  
end subroutine read_RSS_tables

subroutine read_surfmap(lutfile,day)
  
  use LUT_def
  implicit none
  
  character*90 :: lutfile
  
  integer :: d, day
  integer :: header(5)
  integer*1 :: cdata(5760,2880)
  
  print*, lutfile, day
  
  open(unit=1,file=lutfile,status='old',form='unformatted',access='stream')
  do d=1,day
    read(1) header, cdata
    !print*, header
  end do
  close(1)
  LUT%land_class_nx = header(1)
  LUT%land_class_ny = header(2)
  allocate(LUT%land_class_map(LUT%land_class_nx,LUT%land_class_ny))
  LUT%land_class_map = cdata
   

end subroutine read_surfmap

subroutine read_surfmap_month(lutfile)
  
  use LUT_def
  implicit none
  
  character*90 :: lutfile
  
  integer :: header(5)
  integer*1 :: cdata(5760,2880)
  
  print*, lutfile
  
  open(unit=1,file=lutfile,status='old',form='unformatted',access='stream')
  read(1) header, cdata

  LUT%land_class_nx = header(1)
  LUT%land_class_ny = header(2)
  allocate(LUT%land_class_map_bare(LUT%land_class_nx,LUT%land_class_ny))
  allocate(LUT%land_class_map_snow(LUT%land_class_nx,LUT%land_class_ny))
  LUT%land_class_map_bare = cdata
  read(1) header, cdata
  LUT%land_class_map_snow = cdata

end subroutine read_surfmap_month

subroutine read_emis_s0_map(month)

use LUT_def
implicit none

integer :: month, i
character*90 :: lutfile
real :: fdata(1440,560)


LUT%emis_map_nx = 1440
LUT%emis_map_ny = 560
LUT%emis_map_nch = 10
allocate(LUT%emis_map_all(LUT%emis_map_nx,LUT%emis_map_ny,LUT%emis_map_nch))
!allocate(LUT%emis_map_bare(LUT%emis_map_nx,LUT%emis_map_ny,LUT%emis_map_nch))
!allocate(LUT%emis_map_snow(LUT%emis_map_nx,LUT%emis_map_ny,LUT%emis_map_nch))
!allocate(LUT%emis_map_bare(LUT%emis_map_nch,LUT%emis_map_ny,LUT%emis_map_nx))
!allocate(LUT%emis_map_snow(LUT%emis_map_nch,LUT%emis_map_ny,LUT%emis_map_nx))

write(lutfile, '(A16,I2.2,A11)') 'Emiss/emis_mean.',month,'.ITE101.dat'
open(unit = 1, file = lutfile,status = 'old',form='unformatted',access='stream')
do i=1,LUT%emis_map_nch
  read(1) LUT%emis_map_all(:,:,i)
  read(1) fdata!LUT%emis_map_snow(:,:,i)
  read(1) fdata!LUT%emis_map_bare(:,:,i)
end do

!begin  SJM 2/14/17 initiation of new emissivity std deviation variables
allocate(LUT%estd_map_all(LUT%emis_map_nx,LUT%emis_map_ny,LUT%emis_map_nch))
!allocate(LUT%estd_map_bare(LUT%emis_map_nx,LUT%emis_map_ny,LUT%emis_map_nch))
!allocate(LUT%estd_map_snow(LUT%emis_map_nx,LUT%emis_map_ny,LUT%emis_map_nch))
write(lutfile, '(A16,I2.2,A11)') 'Emiss/emis_std.',month,'.ITE101.dat'
open(unit = 1, file = lutfile,status = 'old',form='unformatted',access='stream')
do i=1,LUT%emis_map_nch
  read(1) LUT%estd_map_all(:,:,i)
  read(1) fdata!LUT%estd_map_snow(:,:,i)
  read(1) fdata!LUT%estd_map_bare(:,:,i)
end do
!end    SJM 2/14/17
return
!do i=1,560
!  print*, i/4.-70., maxval(LUT%emis_map_bare(:,i,1)), maxval(LUT%emis_map_bare(:,i,2))
!end do

!print '(10F8.3)', LUT%emis_map_bare(4*(180+77),4*(70+20),:)
!print '(10F8.3)', LUT%emis_map_snow(4*(180+77),4*(70+20),:)
!print '(10F8.3)', LUT%emis_map_bare(:,293,734)
!print '(10F8.3)', LUT%emis_map_snow(:,293,734)
!print*, maxval(LUT%emis_map_bare), maxloc(LUT%emis_map_bare,1)

! LUT%s0Ku_map_nx = 720
! LUT%s0Ku_map_ny = 280
! LUT%s0Ku_map_nray = 49
! allocate(LUT%s0Ku_map_bare(LUT%s0Ku_map_nx,LUT%s0Ku_map_ny,LUT%s0Ku_map_nray))
! allocate(LUT%s0Ku_map_snow(LUT%s0Ku_map_nx,LUT%s0Ku_map_ny,LUT%s0Ku_map_nray))
! allocate(LUT%s0Ku_map_all(LUT%s0Ku_map_nx,LUT%s0Ku_map_ny,LUT%s0Ku_map_nray))
! 
! write(lutfile, '(A20,I2.2,A4)') 'Emiss/s0Ku_mean.0.5.',month,'.dat'
! open(unit = 1, file = lutfile,status = 'old',form='unformatted',access='stream')
! do i=1,LUT%s0Ku_map_nray
!   read(1) LUT%s0Ku_map_bare(:,:,i)
!   read(1) LUT%s0Ku_map_snow(:,:,i)
!   read(1) LUT%s0Ku_map_all(:,:,i)
! end do
! close(1)
! 
! allocate(LUT%s0Ku_std_bare(LUT%s0Ku_map_nx,LUT%s0Ku_map_ny,LUT%s0Ku_map_nray))
! allocate(LUT%s0Ku_std_snow(LUT%s0Ku_map_nx,LUT%s0Ku_map_ny,LUT%s0Ku_map_nray))
! allocate(LUT%s0Ku_std_all(LUT%s0Ku_map_nx,LUT%s0Ku_map_ny,LUT%s0Ku_map_nray))
! 
! write(lutfile, '(A20,I2.2,A4)') 'Emiss/s0Ku_std.0.5.',month,'.dat'
! open(unit = 1, file = lutfile,status = 'old',form='unformatted',access='stream')
! do i=1,LUT%s0Ku_map_nray
!   read(1) LUT%s0Ku_std_bare(:,:,i)
!   read(1) LUT%s0Ku_std_snow(:,:,i)
!   read(1) LUT%s0Ku_std_all(:,:,i)
! end do
! close(1)

!print '(49F5.1)', LUT%s0Ku_map_bare(1*(180+77),1*(70+20),:)
!print '(49F5.1)', LUT%s0Ku_map_snow(1*(180+77),1*(70+20),:)

! LUT%s0Ka_map_nx = 360*2
! LUT%s0Ka_map_ny = 140*2
! LUT%s0Ka_map_nray = 25
! allocate(LUT%s0Ka_map_bare(LUT%s0Ka_map_nx,LUT%s0Ka_map_ny,LUT%s0Ka_map_nray))
! allocate(LUT%s0Ka_map_snow(LUT%s0Ka_map_nx,LUT%s0Ka_map_ny,LUT%s0Ka_map_nray))
! allocate(LUT%s0Ka_map_all(LUT%s0Ka_map_nx,LUT%s0Ka_map_ny,LUT%s0Ka_map_nray))
! write(lutfile, '(A20,I2.2,A4)') 'Emiss/s0Ka_mean.0.5.',month,'.dat'
! open(unit = 1, file = lutfile,status = 'old',form='unformatted',access='stream')
! do i=1,LUT%s0Ka_map_nray
!   read(1) LUT%s0Ka_map_bare(:,:,i)
!   read(1) LUT%s0Ka_map_snow(:,:,i)
!   read(1) LUT%s0Ka_map_all(:,:,i)
! end do
! close(1)
! allocate(LUT%s0Ka_std_bare(LUT%s0Ka_map_nx,LUT%s0Ka_map_ny,LUT%s0Ka_map_nray))
! allocate(LUT%s0Ka_std_snow(LUT%s0Ka_map_nx,LUT%s0Ka_map_ny,LUT%s0Ka_map_nray))
! allocate(LUT%s0Ka_std_all(LUT%s0Ka_map_nx,LUT%s0Ka_map_ny,LUT%s0Ka_map_nray))
! write(lutfile, '(A20,I2.2,A4)') 'Emiss/s0Ka_std.0.5.',month,'.dat'
! open(unit = 1, file = lutfile,status = 'old',form='unformatted',access='stream')
! do i=1,LUT%s0Ka_map_nray
!   read(1) LUT%s0Ka_std_bare(:,:,i)
!   read(1) LUT%s0Ka_std_snow(:,:,i)
!   read(1) LUT%s0Ka_std_all(:,:,i)
! end do
! close(1)

end subroutine read_emis_s0_map

subroutine read_LUTlandclass(filename)

use LUT_def
implicit none
      
integer :: i,j,k,l,ch,recl,n
character*90 :: filename, lutfile
real*8 :: ddata(10)
real*4 :: fdata(49*2)

!hard code these for now
LUT%land_nclass = 13
LUT%land_nch = 10
LUT%land_nrayKu = 49
LUT%land_nrayKa = 49
allocate(LUT%land_class_emis_mean(LUT%land_nclass,LUT%land_nch))
allocate(LUT%land_class_sigma0Ku_mean(LUT%land_nclass,LUT%land_nrayKu))
allocate(LUT%land_class_sigma0Ka_mean(LUT%land_nclass,LUT%land_nrayKa))
allocate(LUT%land_class_sigma0Ku_std(LUT%land_nclass,LUT%land_nrayKu))
allocate(LUT%land_class_sigma0Ka_std(LUT%land_nclass,LUT%land_nrayKa))
allocate(LUT%land_class_emis_eofs_LF(LUT%land_nclass,8+2,10+49+49))
LUT%land_class_emis_eofs_LF = 0.
allocate(LUT%land_class_emis_eofs_AF(LUT%land_nclass,10+2,10+49+49))
LUT%land_class_emis_eofs_AF = 0.
!allocate(LUT%land_class_sigma0Ka_mean(LUT%land_nclass,LUT%land_nrayKa))
!allocate(LUT%land_class_emis_eofs_MS(LUT%land_nclass,LUT%land_nch,LUT%land_nch+50))
!print*, trim(filename)
do i=2,14
  write(lutfile, '(A9,I2.2,A22,A,A4)') 'Emiss/sfc',i,'_emis_s0_eofs.NS.08ch.',trim(filename),'.dat'
  open(unit = 1, file = lutfile,status = 'old',form='unformatted',access='stream')
  write(lutfile, '(A9,I2.2,A22,A,A4)') 'Emiss/sfc',i,'_emis_s0_eofs.NS.10ch.',trim(filename),'.dat'
  open(unit = 2, file = lutfile,status = 'old',form='unformatted',access='stream')
  !print*, lutfile
  !read mean emissivity and sigma_zero for each class
  read(1) ddata(1:8), fdata
  !print '(10F10.3)', ddata
  !print '(25F5.1)', fdata(1:25)
  !print '(25F5.1)', fdata(25:49)
  !print '(25F5.1)', fdata(62:86)
  
  LUT%land_class_emis_mean(i-1,1:8) = ddata(1:8)
  LUT%land_class_sigma0Ku_mean(i-1,:) = fdata(1:49)
  LUT%land_class_sigma0Ka_mean(i-1,:) = fdata(50:98)
  read(2) ddata(1:10), fdata
  !print '(10F10.3)', ddata
  LUT%land_class_emis_mean(i-1,9:10) = ddata(9:10)
  !read standard deviation of emissivity and sigma_zero for each class. Currently not doing anything with this information.
  read(1) ddata(1:8), fdata
  !LUT%land_class_emis_mean(i-1,1:8) = ddata(1:8)
  LUT%land_class_sigma0Ku_std(i-1,:) = fdata(1:49)
  LUT%land_class_sigma0Ka_std(i-1,:) = fdata(50:98)
  !print '(10F10.3)', ddata
  !print '(25F5.1)', fdata(1:25)
  !print '(25F5.1)', fdata(25:49)
  !print '(25F5.1)', fdata(62:86)
  read(2) ddata(1:10), fdata
  !print '(10F10.3)', ddata
  !print '(25F5.1)', fdata(1:25)
  !print '(25F5.1)', fdata(25:49)
  !print '(25F5.1)', fdata(62:86)
  !In 8-channel EOFs, set 166 Ghz EOF8 to standard deviation (shouldn't matter, but at least this provides some surface variability)
  LUT%land_class_emis_eofs_LF(i-1,8,9:10) = ddata(9:10)
  !read in EOFs and regression to sigma_zero
  do j=1,10
    if(j .le. 8) then
      read(1) ddata(1:8), fdata
      !print '(8F7.3)', ddata(1:8)
      !print '(25F5.1)', fdata(1:25)
      !print '(25F5.1)', fdata(25:49)
      !print '(25F5.1)', fdata(62:86)
      LUT%land_class_emis_eofs_LF(i-1,j,1:8) = ddata(1:8)
      LUT%land_class_emis_eofs_LF(i-1,j,11:59) = fdata(1:49)
      LUT%land_class_emis_eofs_LF(i-1,j,60:108) = fdata(50:98)
    endif
    read(2) ddata(1:10), fdata
    !print '(10F7.3)', ddata(1:10)
    !print '(25F5.1)', fdata(1:25)
    !print '(25F5.1)', fdata(25:49)
    !print '(25F5.1)', fdata(62:86)
    LUT%land_class_emis_eofs_AF(i-1,j,1:10) = ddata(1:10)
    LUT%land_class_emis_eofs_AF(i-1,j,11:59) = fdata(1:49)
    LUT%land_class_emis_eofs_AF(i-1,j,60:108) = fdata(50:98)
    
  end do
  
  read(1) fdata
  !print '(25F5.1)', fdata(1:25)
  !print '(25F5.1)', fdata(25:49)
  !print '(25F5.1)', fdata(62:86)
  LUT%land_class_emis_eofs_LF(i-1,9,11:59) = fdata(1:49)
  LUT%land_class_emis_eofs_LF(i-1,9,60:108) = fdata(50:98)
  read(1) fdata
  !print '(25F5.1)', fdata(1:25)
  !print '(25F5.1)', fdata(25:49)
  !print '(25F5.1)', fdata(62:86)
  LUT%land_class_emis_eofs_LF(i-1,10,11:59) = fdata(1:49)
  LUT%land_class_emis_eofs_LF(i-1,10,60:108) = fdata(50:98)
  close(1)
  
  read(2) fdata
  !print '(25F5.1)', fdata(1:25)
  !print '(25F5.1)', fdata(25:49)
  !print '(25F5.1)', fdata(62:86)
  LUT%land_class_emis_eofs_AF(i-1,11,11:59) = fdata(1:49)
  LUT%land_class_emis_eofs_AF(i-1,11,60:108) = fdata(50:98)
  read(2) fdata
  !print '(25F5.1)', fdata(1:25)
  !print '(25F5.1)', fdata(25:49)
  !print '(25F5.1)', fdata(62:86)
  LUT%land_class_emis_eofs_AF(i-1,12,11:59) = fdata(1:49)
  LUT%land_class_emis_eofs_AF(i-1,12,60:108) = fdata(50:98)
  close(2)
  
end do
!stop
end    

!read routine for the new classes for all-sky retrieval
!read in emissivity covariance matrices, emissivity mean and std deviation maps for this month, and class maps for this month
subroutine read_LUT_sfc_class(imonth)

use LUT_def
implicit none
      
integer, intent(in) :: imonth

integer :: i,j,k,l,ch,recl,n
character*90 :: filename
real*8 :: ddata(10)
real*4 :: fdata(11,11)

!emissivity covariance matrices
!hard code these for now
LUT%sfc_nclass = 50 !1-20: bare land, 21-40: snow, 41-50: sea ice
LUT%sfc_nch = 11
allocate(LUT%emis_class_corr(LUT%sfc_nclass,LUT%sfc_nch,LUT%sfc_nch))

!print*, trim(filename)
do i=1,LUT%sfc_nclass
  write(filename, '(A16,I2.2,A14)') 'Emiss/sfc_class.',i,'.emis_corr.bin'
  open(unit = 1, file = filename,status = 'old',form='unformatted',access='stream')
  !print*, lutfile
  !read mean emissivity and sigma_zero for each class
  read(1) fdata
  close(1)
  LUT%emis_class_corr(i,:,:) = fdata

  
end do

!emissivity mean and std maps
print*,imonth
LUT%emis_map_nx = 1440
LUT%emis_map_ny = 560
LUT%emis_map_nch = 11
allocate(LUT%emis_map_bare(LUT%emis_map_nch,LUT%emis_map_nx,LUT%emis_map_ny))
allocate(LUT%emis_map_snow(LUT%emis_map_nch,LUT%emis_map_nx,LUT%emis_map_ny))
allocate(LUT%emis_map_seaice(LUT%emis_map_nch,LUT%emis_map_nx,LUT%emis_map_ny))
allocate(LUT%estd_map_bare(LUT%emis_map_nch,LUT%emis_map_nx,LUT%emis_map_ny))
allocate(LUT%estd_map_snow(LUT%emis_map_nch,LUT%emis_map_nx,LUT%emis_map_ny))
allocate(LUT%estd_map_seaice(LUT%emis_map_nch,LUT%emis_map_nx,LUT%emis_map_ny))

write(filename, '(A25,I2.2,A4)') 'Emiss/gmi_emis_mean.bare.',imonth,'.bin'
open(unit=1,file = filename, status='old', form='unformatted', access='stream')
read(1) LUT%emis_map_bare
close(1)
!print*, LUT%emis_map_bare(:,329,389)
!stop
write(filename, '(A25,I2.2,A4)') 'Emiss/gmi_emis_mean.snow.',imonth,'.bin'
open(unit=1,file = filename, status='old', form='unformatted', access='stream')
read(1) LUT%emis_map_snow
close(1)

write(filename, '(A27,I2.2,A4)') 'Emiss/gmi_emis_mean.seaice.',imonth,'.bin'
open(unit=1,file = filename, status='old', form='unformatted', access='stream')
read(1) LUT%emis_map_seaice
close(1)

write(filename, '(A25,I2.2,A4)') 'Emiss/gmi_emis_std.bare.',imonth,'.bin'
open(unit=1,file = filename, status='old', form='unformatted', access='stream')
read(1) LUT%estd_map_bare
close(1)

write(filename, '(A25,I2.2,A4)') 'Emiss/gmi_emis_std.snow.',imonth,'.bin'
open(unit=1,file = filename, status='old', form='unformatted', access='stream')
read(1) LUT%estd_map_snow
close(1)

write(filename, '(A27,I2.2,A4)') 'Emiss/gmi_emis_std.seaice.',imonth,'.bin'
open(unit=1,file = filename, status='old', form='unformatted', access='stream')
read(1) LUT%estd_map_seaice
close(1)

!class maps
LUT%sfc_class_nx = 1440
LUT%sfc_class_ny = 560
allocate(LUT%bare_class_map(LUT%sfc_class_nx,LUT%sfc_class_ny))
allocate(LUT%snow_class_map(LUT%sfc_class_nx,LUT%sfc_class_ny))
allocate(LUT%seaice_class_map(LUT%sfc_class_nx,LUT%sfc_class_ny))

write(filename, '(A17,I2.2,A4)') 'Emiss/land_class.',imonth,'.bin'
open(unit=1,file = filename, status='old', form='unformatted', access='stream')
read(1) LUT%bare_class_map
close(1)

write(filename, '(A17,I2.2,A4)') 'Emiss/snow_class.',imonth,'.bin'
open(unit=1,file = filename, status='old', form='unformatted', access='stream')
read(1) LUT%snow_class_map
close(1)

write(filename, '(A19,I2.2,A4)') 'Emiss/seaice_class.',imonth,'.bin'
open(unit=1,file = filename, status='old', form='unformatted', access='stream')
read(1) LUT%seaice_class_map
close(1)

!stop
end    

subroutine get_gmi_sfc_class(lat, lon, sfc_type, ice_type)

use LUT_def

implicit none

real, intent(in) :: lat, lon

integer, intent(out) :: sfc_type, ice_type
integer :: snow_type, seaice_type
integer :: x, y

x = mod(floor((lon+180.)*LUT%sfc_class_nx/360),LUT%sfc_class_nx)+1
y = floor((lat+70.)*LUT%sfc_class_ny)/140+1

sfc_type = LUT%bare_class_map(x,y)
snow_type = LUT%snow_class_map(x,y)
seaice_type = LUT%seaice_class_map(x,y)
if(snow_type .gt. 0) then
	ice_type = snow_type+20
else if(seaice_type .gt. 0) then !sea ice supersedes land ice
	ice_type = seaice_type+30
else
	ice_type = 0
endif


end subroutine get_gmi_sfc_class

subroutine get_gmi_emis_std(lat, lon, stype, emis, estd)

use LUT_def

implicit none

real, intent(in) :: lat, lon
integer, intent(in) :: stype
real, intent(out) :: emis(LUT%sfc_nch), estd(LUT%sfc_nch)
integer :: x, y

x = mod(floor((lon+180.)*LUT%sfc_class_nx/360),LUT%sfc_class_nx)+1
y = floor((lat+70.)*LUT%sfc_class_ny)/140+1
emis(:) = 0
estd(:) = 0
!print '(2F8.2,3I5,11F6.2)', lat, lon, stype, x, y, LUT%emis_map_bare(:,x,y)
if((stype .ge. 1) .and. (stype .le. 20)) then !bare land
  !print*, 'bare'
  emis(:) = LUT%emis_map_bare(:,x,y)
  estd(:) = LUT%estd_map_bare(:,x,y)
else if ((stype .ge. 21) .and. (stype .le. 40)) then !snow
  !print*, 'snow'
  emis(:) = LUT%emis_map_snow(:,x,y)
  estd(:) = LUT%estd_map_snow(:,x,y)
else if((stype .ge. 41) .and. (stype .le. 50)) then !sea ice
  !print*, 'sea ice'
  emis(:) = LUT%emis_map_seaice(:,x,y)
  estd(:) = LUT%estd_map_seaice(:,x,y)
endif
!print '(2F8.2,3I5,11F6.2)', lat, lon, stype, x, y, estd
!extrapolate 89 to 166 if 166 is not in database
if((emis(10) .le. 0) .and. (emis(8) .gt. 0)) then
  emis(10) = emis(8)
  estd(10) = estd(8)
endif
if((emis(11) .le. 0) .and. (emis(9) .gt. 0)) then
  emis(11) = emis(9)
  estd(11) = estd(9)
endif

end subroutine get_gmi_emis_std

subroutine read_LUTgas_monortm(filename)

use LUT_def
implicit none
      
integer :: i,j,k,l,ch,nchan
character*90 :: filename
real, allocatable :: freq(:), kext(:,:,:,:)
integer, allocatable :: ifreq(:), ipol(:)
!print*, filename

open(unit=10, file=trim(filename), access='stream')

read(10) LUT%nch
! print*, LUT%nch
allocate(freq(LUT%nch))
read(10) freq
! print*, freq
read(10) nchan
! print*, nchan
allocate(ifreq(nchan))
allocate(ipol(nchan))
read(10) ifreq
! print*, ifreq
read(10) ipol
! print*, ipol

read(10) LUT%gas_nP
! print*, LUT%gas_nP
allocate(LUT%gas_press(LUT%gas_nP))
read(10) LUT%gas_press
! print*,  LUT%gas_press
read(10) LUT%gas_nT
! print*, LUT%gas_nT
allocate(LUT%gas_temp(LUT%gas_nT))
read(10) LUT%gas_temp
! print*, LUT%gas_temp
read(10) LUT%gas_nV
! print*, LUT%gas_nV
allocate(LUT%gas_wv(LUT%gas_nV))
read(10) LUT%gas_wv
! print*, LUT%gas_wv
LUT%gas_wv = 1e-3*LUT%gas_wv
allocate(kext(LUT%gas_nP,LUT%gas_nT,LUT%gas_nV,LUT%nch))
read(10) kext
close(10)
! stop
allocate(LUT%gas_kext(LUT%nch+2,LUT%gas_nP,LUT%gas_nT,LUT%gas_nV))
do i=1,LUT%nch
  LUT%gas_kext(i,:,:,:) = kext(:,:,:,i)
end do



!read in DPR LUT
deallocate(freq)
deallocate(ifreq)
deallocate(ipol)
deallocate(kext)
open(unit=10, file='AncData/berg/MonoRTM_v5.3-DPR.tbl', access='stream')

read(10) LUT%nch
! print*, LUT%nch
allocate(freq(LUT%nch))
read(10) freq
! print*, freq
read(10) nchan
! print*, nchan
allocate(ifreq(nchan))
allocate(ipol(nchan))
read(10) ifreq
! print*, ifreq
read(10) ipol
! print*, ipol

read(10) LUT%gas_nP
! print*, LUT%gas_nP
! allocate(LUT%gas_press(LUT%gas_nP))
read(10) LUT%gas_press
! print*,  LUT%gas_press
read(10) LUT%gas_nT
! print*, LUT%gas_nT
! allocate(LUT%gas_temp(LUT%gas_nT))
read(10) LUT%gas_temp
! print*, LUT%gas_temp
read(10) LUT%gas_nV
! print*, LUT%gas_nV
! allocate(LUT%gas_wv(LUT%gas_nV))
read(10) LUT%gas_wv
!print*, LUT%gas_wv
LUT%gas_wv = 1e-3*LUT%gas_wv
allocate(kext(LUT%gas_nP,LUT%gas_nT,LUT%gas_nV,LUT%nch))
read(10) kext
close(10)

do i=14,15
  LUT%gas_kext(i,:,:,:) = kext(:,:,:,i-13)
end do

! print '(15F8.2)', LUT%gas_press
! print '(15F8.2)', LUT%gas_temp
! print '(15F8.4)', LUT%gas_wv
! stop

end      

! subroutine get_s0(lat, lon, ray, stype, s0Ku, s0Ka, s0Kustd, s0Kastd)
! 
! use LUT_def
! 
! implicit none
! 
! real :: lon, lat, s0Ku, s0Ka, s0Kustd, s0Kastd
! 
! integer :: x, y, ray, stype
! 
! s0Ku = -99.
! s0Kustd = -99.
! s0Ka = -99.
! s0Kastd = -99.
! 
! if(stype .eq. 1) return
! 
! x = mod(floor(lon+180.)*LUT%s0Ku_map_nx/360,LUT%s0Ku_map_nx)+1
! y = floor(lat+70.)*LUT%s0Ku_map_ny/140.+1
! 
! if(stype .eq. 2 .or. (stype .ge. 8 .and. stype .le. 11) .or. stype .eq. 14) then
!   s0Ku = LUT%s0Ku_map_snow(x,y,ray)
!   s0Kustd = LUT%s0Ku_std_snow(x,y,ray)
!   if(ray .gt. 12 .and. ray .lt. 38) then
!     s0Ka = LUT%s0Ka_map_snow(x,y,ray-12)
!     s0Kastd = LUT%s0Ka_map_snow(x,y,ray-12)
!   endif
! else
!   s0Ku = LUT%s0Ku_map_bare(x,y,ray)
!   s0Kustd = LUT%s0Ku_std_bare(x,y,ray)
!   if(ray .gt. 12 .and. ray .lt. 38) then
!     s0Ka = LUT%s0Ka_map_bare(x,y,ray-12)
!     s0Kastd = LUT%s0Ka_std_bare(x,y,ray-12)
!   endif
! endif
! 
! if(s0Ku .eq. -99.) then
!   s0Ku = LUT%s0Ku_map_all(x,y,ray)
!   s0Kustd = LUT%s0Ku_std_all(x,y,ray)
! endif
! if(s0Ku .eq. -99.) then
!   s0Ku = LUT%land_class_sigma0Ku_mean(stype-1,ray)
!   s0Kustd = LUT%land_class_sigma0Ku_std(stype-1,ray)
! endif
! if(s0Kustd .lt. 0.) s0Kustd = LUT%land_class_sigma0Ku_std(stype-1,ray)
! 
! if(s0Ka .eq. -99. .and. ray .gt. 12 .and. ray .lt. 38) then
!   s0Ka = LUT%s0Ka_map_all(x,y,ray-12)
!   s0Kastd = LUT%s0Ka_std_all(x,y,ray-12)
! endif
! if(s0Ka .eq. -99. .and. ray .gt. 12 .and. ray .lt. 38) then
!   s0Ka = LUT%land_class_sigma0Ka_mean(stype-1,ray)
!   s0Kastd = LUT%land_class_sigma0Ka_std(stype-1,ray)
! endif
! if(s0Kastd .lt. 0.) s0Kastd = LUT%land_class_sigma0Ka_std(stype-1,ray)
! 
! end subroutine get_s0

subroutine read_LUTwatemis(filename)

use LUT_def

implicit none
      
integer :: i,j,k,l,ch,junk
character*90 :: filename

!open LUT file
open(unit = 1, file = filename,status = 'old',form='unformatted')
     
!Check that LUT dimensions are what we think they are
read(1) LUT%nch,LUT%water_emis_nT,LUT%water_emis_nW, LUT%water_emis_nA
! print*, LUT%nch,LUT%water_emis_nT,LUT%water_emis_nW, LUT%water_emis_nA

allocate(LUT%water_emis_temp(LUT%water_emis_nT), LUT%water_emis_wind(LUT%water_emis_nW))
allocate(LUT%water_emis_angle(LUT%nch,LUT%water_emis_nA))
allocate(LUT%water_emis(LUT%nch,LUT%water_emis_nT,LUT%water_emis_nW,LUT%water_emis_nA,2))
allocate(LUT%water_eharm(LUT%nch,LUT%water_emis_nT,LUT%water_emis_nW,LUT%water_emis_nA,2,2))
allocate(LUT%water_ebar(LUT%nch,LUT%water_emis_nT,LUT%water_emis_nW))

!Read LUT
do ch = 1, LUT%nch
  do k = 1, LUT%water_emis_nT
    do l = 1, LUT%water_emis_nW
      do i=1, LUT%water_emis_nA
        read(1) junk,LUT%water_emis_temp(k),LUT%water_emis_wind(l),LUT%water_emis_angle(ch,i), LUT%water_emis(ch,k,l,i,:),LUT%water_ebar(ch,k,l),LUT%water_eharm(ch,k,l,i,1,1),LUT%water_eharm(ch,k,l,i,1,2),LUT%water_eharm(ch,k,l,i,2,1),LUT%water_eharm(ch,k,l,i,2,2)
      end do
    end do
  end do
end do
close(1)
! 
! print '(15F8.2)', LUT%water_emis_temp
! print '(15F8.2)', LUT%water_emis_wind
! print '(10F8.2)', LUT%water_emis_angle
! print '(10F8.4)', LUT%water_emis(:,5,1,:,2)
! print '(10F8.4)', LUT%water_eharm(1,5,:,1,2,2)
! stop
end      

subroutine read_lutwatsigma0(filename)

  use LUT_def
  
  implicit none
  
  character*90 :: filename
  integer :: irec, nr, nw, junk
  real :: fdata(12)

  open(unit = 1, file = filename,status = 'old',form='unformatted', recl=12,access='direct')
  read(1,rec=1) LUT%water_sigma0_nA, LUT%water_sigma0_nW
  !print*, LUT%water_sigma0_nA, LUT%water_sigma0_nW
  
  allocate(LUT%water_sigma0_wind(LUT%water_sigma0_nW))
  allocate(LUT%water_sigma0_Ku(LUT%water_sigma0_nA, LUT%water_sigma0_nW))
  allocate(LUT%water_sigma0_Ku_std(LUT%water_sigma0_nA, LUT%water_sigma0_nW))
  allocate(LUT%water_sigma0_Ku_harm(LUT%water_sigma0_nA, LUT%water_sigma0_nW,2))
  allocate(LUT%water_sigma0_Ka(LUT%water_sigma0_nA, LUT%water_sigma0_nW))
  allocate(LUT%water_sigma0_Ka_std(LUT%water_sigma0_nA, LUT%water_sigma0_nW))
  allocate(LUT%water_sigma0_Ka_harm(LUT%water_sigma0_nA, LUT%water_sigma0_nW,2))
  allocate(LUT%water_sigma0_corr(LUT%water_sigma0_nA, LUT%water_sigma0_nW))
  allocate(LUT%water_sigma0_diff(LUT%water_sigma0_nA, LUT%water_sigma0_nW))
  allocate(LUT%water_sigma0_dstd(LUT%water_sigma0_nA, LUT%water_sigma0_nW))
  
  irec=1
  do nr = 1, LUT%water_sigma0_nA
    do nW = 1, LUT%water_sigma0_nW
      irec=irec+1
      read(1,rec=irec) fdata
      !print '(2I5,12F10.3)', nr, nw, fdata
      LUT%water_sigma0_wind(nw) = fdata(1)
      LUT%water_sigma0_Ku(nr,nw) = fdata(2)
      LUT%water_sigma0_Ku_std(nr,nw) = fdata(3)
      LUT%water_sigma0_Ku_harm(nr,nw,:) = fdata(4:5)
      LUT%water_sigma0_Ka(nr,nw) = fdata(6)
      LUT%water_sigma0_Ka_std(nr,nw) = fdata(7)
      LUT%water_sigma0_Ka_harm(nr,nw,:) = fdata(8:9)
      LUT%water_sigma0_corr(nr,nw) = fdata(10)
      LUT%water_sigma0_diff(nr,nw) = fdata(11)
      LUT%water_sigma0_dstd(nr,nw) = fdata(12)
      !if(nw .eq. 7) print*,nr, LUT%water_sigma0_wind(nw), LUT%water_sigma0_Ku(nr,nw), LUT%water_sigma0_Ka(nr,nw)
    end do
    
  end do
  close(1)
  !stop
!   
!   print*, LUT%water_sigma0_wind
!   stop
end subroutine read_lutwatsigma0


subroutine get_indx(Var, LUT_nVar, LUT_Var, nerror, idx_Var)
       
!--Code determines the index (entry just above the variable of
!--interest in the LUT).  LUT(index-1) < Variable < LUT(index)

implicit   none
        
real       LUT_Var(LUT_nvar)
real 	   Var
integer    i, nerror, LUT_nVar, idx_Var
logical    check_bdry                  
      
 check_bdry = .true.
 i = 1
 do while (check_bdry)
   i = i + 1
   if ( i .gt. LUT_nVar) then
     print*, Var
     call err(nerror)
   endif
   if (Var .le. LUT_Var(i)) check_bdry = .false.
 end do
 !if ( i .gt. LUT_nVar) then
 !  call err(nerror)
 !endif
 idx_Var =  i

 return
end

subroutine intplte_emis(nf,np,temp,wind,relAz,eia,emis,ebar)

use LUT_def
implicit   none

real       temp, T, wind,W, dT, dW, relAz, eia, A, dA
integer    nf, np, iT, iW, iA
real	   weights(2,2,2)
real       emis,refl,ebar, eharm(2)
real	   e0(4), ewind(2), eharmMW(2,4)
      

T = temp
if ( T .lt. LUT%water_emis_temp(1) ) T = LUT%water_emis_temp(1)
if ( T .gt. LUT%water_emis_temp(LUT%water_emis_nT) ) T = LUT%water_emis_temp(LUT%water_emis_nT)
W = wind
if ( W .lt. LUT%water_emis_wind(1) ) W = LUT%water_emis_wind(1)
if ( W .gt. LUT%water_emis_wind(LUT%water_emis_nW) ) W = LUT%water_emis_wind(LUT%water_emis_nW)
A = eia
if ( A .lt. LUT%water_emis_angle(nf,1) ) A = LUT%water_emis_angle(nf,1)
if ( A .gt. LUT%water_emis_angle(nf,LUT%water_emis_nA) ) A = LUT%water_emis_angle(nf,LUT%water_emis_nA)

! print*, temp,T
! print*, wind, W
! stop
!--Simple bilinear interpolation

call get_indx(T, LUT%water_emis_nT, LUT%water_emis_temp, 115, iT)
call get_indx(W, LUT%water_emis_nW, LUT%water_emis_wind, 116, iW)
call get_indx(A, LUT%water_emis_nA, LUT%water_emis_angle(nf,:), 117, iA)

dT = (T - LUT%water_emis_temp(iT-1))/(LUT%water_emis_temp(iT) - LUT%water_emis_temp(iT-1))
dW = (W - LUT%water_emis_wind(iW-1))/(LUT%water_emis_wind(iW) - LUT%water_emis_wind(iW-1))
dA = (A - LUT%water_emis_angle(nf,iA-1))/(LUT%water_emis_angle(nf,iA) - LUT%water_emis_angle(nf,iA-1))
! print*, dT, dW
! stop
weights(1,1,1) = (1.-dT)*(1.-dW)*(1.-dA)
weights(1,2,1) = (1.-dT)*(   dW)*(1.-dA)
weights(2,1,1) = (   dT)*(1.-dW)*(1.-dA)
weights(2,2,1) = (   dT)*(   dW)*(1.-dA)
weights(1,1,2) = (1.-dT)*(1.-dW)*(   dA)
weights(1,2,2) = (1.-dT)*(   dW)*(   dA)
weights(2,1,2) = (   dT)*(1.-dW)*(   dA)
weights(2,2,2) = (   dT)*(   dW)*(   dA)

emis=sum(weights*LUT%water_emis(nf,iT-1:iT,iW-1:iW,iA-1:iA,np+1))
ebar=sum(weights(:,:,1)*LUT%water_ebar(nf,iT-1:iT,iW-1:iW))+sum(weights(:,:,2)*LUT%water_ebar(nf,iT-1:iT,iW-1:iW))
eharm(1)=sum(weights*LUT%water_eharm(nf,iT-1:iT,iW-1:iW,iA-1:iA,np+1,1))
eharm(2)=sum(weights*LUT%water_eharm(nf,iT-1:iT,iW-1:iW,iA-1:iA,np+1,2))

emis = emis+eharm(1)*cos(relAz)+eharm(2)*cos(2.*relAz)
      
return 
end


subroutine intplte_water_sigma0(ray,wind,relAz,sigma0_Ku,sigma0_Ka,std_Ku,std_Ka,corr)
  use LUT_def
  implicit   none

  real :: wind, w, dW
  integer :: ray,r, iW
  real :: weights(2), a_Ku(2), a_Ka(2), relAz
  real :: sigma0_ku,sigma0_Ka,std_Ku,std_Ka,corr
      
  r = ray!abs(25-ray)+1
  w = wind
  if (w .lt. LUT%water_sigma0_wind(1) ) w = LUT%water_sigma0_wind(1)
  if (w .gt. LUT%water_sigma0_wind(LUT%water_sigma0_nW) ) w = LUT%water_sigma0_wind(LUT%water_sigma0_nW)

  call get_indx(w, LUT%water_sigma0_nW, LUT%water_sigma0_wind, 118, iW)
   
  dW = (w - LUT%water_sigma0_wind(iW-1))/(LUT%water_sigma0_wind(iW) - LUT%water_sigma0_wind(iW-1))

  weights(1) = (1.-dW)
  weights(2) = (   dW)
      
  sigma0_Ku=sum(weights*LUT%water_sigma0_Ku(r,iW-1:iW))
  sigma0_Ka=sum(weights*LUT%water_sigma0_Ka(r,iW-1:iW))
  a_Ku(1) = sum(weights*LUT%water_sigma0_Ku_harm(r,iW-1:iW,1))
  a_Ku(2) = sum(weights*LUT%water_sigma0_Ku_harm(r,iW-1:iW,2))
  !print*, ray, wind, relAz, a_Ku
  !sigma0_Ku = sigma0_Ku + a_Ku(1)*cos(relAz) + a_Ku(2)*cos(2*relAz)
  if(ray .ge. 13 .and. ray .le. 37) then
    a_Ka(1) = sum(weights*LUT%water_sigma0_Ka_harm(r,iW-1:iW,1))
    a_Ka(2) = sum(weights*LUT%water_sigma0_Ka_harm(r,iW-1:iW,2))
    sigma0_Ka = sigma0_Ka + a_Ka(1)*cos(relAz) + a_Ka(2)*cos(2*relAz)
  endif
  !print*, ray, wind, relAz, a_Ku
  std_Ku=sum(weights*LUT%water_sigma0_Ku_std(r,iW-1:iW))
  std_Ka=sum(weights*LUT%water_sigma0_Ka_std(r,iW-1:iW))
  corr=sum(weights*LUT%water_sigma0_corr(r,iW-1:iW))
  return 
end subroutine intplte_water_sigma0

subroutine err(n)
integer   n      
if (n .eq. 101) then
  write(*,*) ' Error in rain LUT frequencies '
  stop
elseif (n .eq. 102) then 
  write(*,*) ' Error in snow LUT frequencies '
  stop
elseif (n .eq. 103) then         
  write(*,*) ' Error in graupel LUT frequencies '
  stop
elseif (n .eq. 104) then         
  write(*,*) ' Error in emissivity LUT frequencies '
  stop
elseif (n .eq. 105) then         
  write(*,*) ' Error in drizzle LUT frequencies '
  stop
elseif (n .eq. 115) then         
  write(*,*) ' Temperature is out of bounds in water emissivity LUT '
  stop
elseif (n .eq. 116) then         
  write(*,*) ' Wind is out of bounds in water emissivity LUT '
  stop
elseif (n .eq. 117) then         
  write(*,*) ' incidence Angle is out of bounds in water emissivity LUT '
  stop
elseif (n .eq. 118) then         
  write(*,*) ' Wind is out of bounds in sigma_zero LUT '
  stop
elseif (n .eq. 141) then         
  write(*,*) ' Tskin is outside of Ts range in emisivity LUT'
  stop
elseif (n .eq. 142) then         
  write(*,*) ' sfc wind is outside wind range in emissvity LUT'
  stop
elseif (n .eq. 151) then         
  write(*,*) ' Tb is outside physical range'
  stop
else
  write(*,*) ' Unspecified error'
  stop
endif
end

subroutine calc_relAz(lon0, lat0, lon1, lat1, u, v, relAz)

implicit none

real, intent(in) :: lon0, lat0, lon1, lat1, u, v
real, intent(out) :: relAz
real :: x0,y0,x1,y1, dotprod, norm

relAz = 0.
y0=lat0
y1=lat1
x0=lon0
x1=lon1
if(abs(lon0-lon1) .gt. 90.) then
  if(lon0 .lt. 0. .and. lon1 .gt. 0.) x0 = lon0+360.
  if(lon1 .lt. 0. .and. lon0 .gt. 0.) x1 = lon1+360. 
endif
  !unfold longitude by adding 360 to negative

dotprod = (x0-x1)*u+(y0-y1)*v
norm = sqrt((x0-x1)**2+(y0-y1)**2)*sqrt(u**2+v**2)
    !print*, dotprod, norm

if(abs(dotprod/norm) .lt. 1.) relAz = acos(dotprod/norm)

end

subroutine gmiretsub(icL,ichunk,orbNumb,ialg,idir)
  use globalData
  use f90DataTypes
  use f90Types
  use profile_def
  use LUT_def
  use ocean_ret_raob
  use missingMod
  
  implicit none
  
  include 'parametersMERRA.inc'
  include 'constants_sjm.inc'
  
  integer*4 :: ichunk
  integer :: ialg, icL, orbNumb, idir
  
  integer :: i,j,kp,kg,i2,j2
  integer :: minGMIscan(2), maxGMIscan(2),minGMIpix(2),maxGMIpix(2) !S1 and S2
  
  integer nearest_gmiS1(3,2,2), nearest_gmiS2(3,2,2) !3 rays, 2 scan ends, and 2 indices
  
  integer :: dpri(3), dprj(2), S2i(2), idpr,jdpr, gloc(2)
  
  real :: gmi_dist(221,500)
  real :: dp, wf, rho, maxdist, env_wgt, clwp
  type(profile) :: gmiprof(221), gmidata_temp
  type(ret_data) :: gmi_ret_ocean, gmi_ret_land, gmi_ret_ice
  type(ret_data), dimension(:,:), allocatable :: gmi_ret
  real, dimension(:,:), allocatable :: gmi_dist2, gmi_wgt, gmi_mask
  real :: gmi_temp_prof(MERRA_NLEV), gmi_qv_prof(MERRA_NLEV), gmi_cloud_prof(MERRA_NLEV)
  
  !!temporary storage for variables interpolated to DPR footprint
  !real, dimension(:,:), allocatable :: dprSknTemp, dprSfcTemp, dprSfcWind
  
  
  print*, orbNumb,ichunk,icL
  print*, gmidata%n1b11
  
  print*, size(gmidata%S1lon,1), size(gmidata%S1lon,2)

  !Find the range of gmi scans and pixels to run retrieval for this chunk
  !set min and max scan based on chunk number
 

  print*, ndpr, dprdata%n1c21
  minGMIscan(1) = icL/float(ndpr)*gmidata%n1b11-200
  if(minGMIscan(1) < 1) minGMIscan(1) = 1
  maxGMIscan(1) = minGMIscan(1)+500
  if(maxGMIscan(1) > gmidata%n1b11) then
  	maxGMIscan(1) = gmidata%n1b11
  	minGMIscan(1) = maxGMIscan(1)-500
  endif
  print*,minGMIscan(1),maxGMIscan(1)
  !find matches to scan dpr scan center and ends at first and last scan in this chunk
  print '(4F8.2)', dprdata%xlon(1,1), dprdata%xlat(1,1), dprdata%xlon(1,dprdata%n1c21), dprdata%xlat(1,dprdata%n1c21)
  print '(4F8.2)', dprdata%xlon(25,1), dprdata%xlat(25,1), dprdata%xlon(25,dprdata%n1c21), dprdata%xlat(25,dprdata%n1c21)
  print '(4F8.2)', dprdata%xlon(49,1), dprdata%xlat(49,1), dprdata%xlon(49,dprdata%n1c21), dprdata%xlat(49,dprdata%n1c21)
  dpri = (/1,25,49/)
  dprj(1) = 1
  dprj(2) = dprdata%n1c21
  do i=1,3
    do j=1,2
      print '(2I5,2F8.2)', dpri(i), dprj(j), dprdata%xlon(dpri(i),dprj(j)), dprdata%xlat(dpri(i),dprj(j))
      gmi_dist = (dprdata%xlon(dpri(i),dprj(j))-gmidata%S1lon(:,minGMIscan(1):maxGMIscan(1)))**2+(dprdata%xlat(dpri(i),dprj(j))-gmidata%S1lat(:,minGMIscan(1):maxGMIscan(1)))**2
      nearest_gmiS1(i,j,:) = minloc(gmi_dist)
      nearest_gmiS1(i,j,2) = nearest_gmiS1(i,j,2)+minGMIscan(1)-1
      gmi_dist = (dprdata%xlon(dpri(i),dprj(j))-gmidata%S2lon(:,minGMIscan(1):maxGMIscan(1)))**2+(dprdata%xlat(dpri(i),dprj(j))-gmidata%S2lat(:,minGMIscan(1):maxGMIscan(1)))**2
      nearest_gmiS2(i,j,:) = minloc(gmi_dist)
      nearest_gmiS2(i,j,2) = nearest_gmiS2(i,j,2)+minGMIscan(1)-1
	  print '(2I5,2F8.2,2I5,2F8.2)', nearest_gmiS1(i,j,:), gmidata%S1lon(nearest_gmiS1(i,j,1),nearest_gmiS1(i,j,2)), &
  	                                           gmidata%S1lat(nearest_gmiS1(i,j,1),nearest_gmiS1(i,j,2)), &
  	                                           nearest_gmiS2(i,j,:), gmidata%S2lon(nearest_gmiS2(i,j,1),nearest_gmiS2(i,j,2)), &
  	                                           gmidata%S2lat(nearest_gmiS2(i,j,1),nearest_gmiS2(i,j,2)) 
    end do                                    
  end do
  
  !Run GMI retrieval over identified scan range (add 1 pixel buffer for interpolation to DPR)
  minGMIscan(1) = minval(nearest_GMIS1(:,:,2))-1
  maxGMIscan(1) = maxval(nearest_GMIS1(:,:,2))+1
  minGMIpix(1) = minval(nearest_GMIS1(:,:,1))-1
  maxGMIpix(1) = maxval(nearest_GMIS1(:,:,1))+1
  minGMIscan(2) = minval(nearest_GMIS2(:,:,2))-1
  maxGMIscan(2) = maxval(nearest_GMIS2(:,:,2))+1
  minGMIpix(2) = minval(nearest_GMIS2(:,:,1))-1
  maxGMIpix(2) = maxval(nearest_GMIS2(:,:,1))+1
  
  !Bounds check the min/max GMI pix - this is only needed when yaw flips mess up the scan geometry
  if(minGMIpix(1) .lt. 50) minGMIpix(1) = 50
  if(maxGMIpix(1) .gt. 170) maxGMIpix(1) = 170
  if(minGMIpix(2) .lt. 50) minGMIpix(2) = 50
  if(maxGMIpix(2) .gt. 170) maxGMIpix(2) = 170
  print*, minGMIscan, maxGMIscan, minGMIpix, maxGMIpix
  !allocate GMI ret array
  allocate(gmi_ret(maxGMIpix(1)-minGMIpix(1)+1,maxGMIscan(1)-minGMIscan(1)+1))
  allocate(gmi_dist2(maxGMIpix(1)-minGMIpix(1)+1,maxGMIscan(1)-minGMIscan(1)+1))
  allocate(gmi_wgt(maxGMIpix(1)-minGMIpix(1)+1,maxGMIscan(1)-minGMIscan(1)+1))
  allocate(gmi_mask(maxGMIpix(1)-minGMIpix(1)+1,maxGMIscan(1)-minGMIscan(1)+1))
  
  !open(unit=12,file='temp/input.bin',form='unformatted',status='unknown')
  print*, 'Running GMI 1DVAR'
  do i = minGMIscan(1),maxGMIscan(1)
    !set up GMI data structures
    !print*, i
    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP DO PRIVATE (j,gmi_dist,gmiprof,S2i,idpr,jdpr,kg,kp,dp,rho,wf,gmi_ret_ocean,gmi_ret_land,gmi_ret_ice,gmidata_temp,i2,j2)
    do j=minGMIpix(1),maxGMIpix(1)
      !print*,i,j
      gmiprof(j)%sclat = gmidata%sclat(i)
      gmiprof(j)%sclon = gmidata%sclon(i)
      gmiprof(j)%scalt = gmidata%scalt(i)
      gmiprof(j)%scyaw = gmidata%scyaw(i)
      gmiprof(j)%lat = gmidata%S1lat(j,i)
      gmiprof(j)%lon = gmidata%S1lon(j,i)
      gmiprof(j)%scan = i
      gmiprof(j)%pos = j
      gmiprof(j)%gmi_tb(1:9) = gmidata%gmiS1(:,j,i)
      
      !get nearest neighbor Tbs from S2
      gmi_dist(minGMIpix(2):maxGMIpix(2),1:maxGMIscan(2)-minGMIscan(2)+1) = &
      (gmiprof(j)%lon-gmidata%S2lon(minGMIpix(2):maxGMIpix(2),minGMIscan(2):maxGMIscan(2)))**2 + &
      (gmiprof(j)%lat-gmidata%S2lat(minGMIpix(2):maxGMIpix(2),minGMIscan(2):maxGMIscan(2)))**2
      S2i = minloc(gmi_dist(minGMIpix(2):maxGMIpix(2),1:maxGMIscan(2)-minGMIscan(2)+1))
      S2i(1) = S2i(1)+ minGMIpix(2)-1
      S2i(2) = S2i(2)+ minGMIscan(2)-1
      !print '(2I8, 2F8.2)', S2i, gmidata%S2lat(S2i(1),S2i(2)),gmidata%S2lon(S2i(1),S2i(2))
      !write(12) gmiprof(j)%lat, gmiprof(j)%lon, gmidata%S2lat(S2i(1),S2i(2)),gmidata%S2lon(S2i(1),S2i(2))
      gmiprof(j)%gmi_tb(10:13) = gmidata%gmiS2(:,S2i(1),S2i(2))
      gmiprof(j)%eia(1) = gmidata%S1eia(j,i)
      gmiprof(j)%eia(2) = gmidata%S2eia(S2i(1),S2i(2))
      gmiprof(j)%sga = gmidata%S1sga(j,i)
      !print '(2F8.2, F7.2, I3, 2F8.2, 2I5, 13F6.1, 3F7.2)', gmiprof(j)%sclat, gmiprof(j)%sclon, gmiprof(j)%scalt, gmiprof(j)%scyaw, gmiprof(j)%lat, gmiprof(j)%lon, &
      !                                     gmiprof(j)%scan, gmiprof(j)%pos, gmiprof(j)%gmi_tb, gmiprof(j)%eia, gmiprof(j)%sga
      
      !get ancillary data (a priori profiles, tskin) form DPR Env data & interpolate to MERRA pressure levels
      !Find nearest DPR index - doesn't have to be exact match, since ancillary data are coarse anyway
      idpr = int((i-nearest_GMIS1(1,1,2))/float(nearest_GMIS1(1,2,2)-nearest_GMIS1(1,1,2))*dprdata%n1c21)
      jdpr = int((j-nearest_GMIS1(1,1,1))/float(nearest_GMIS1(3,1,1)-nearest_GMIS1(1,1,1))*49)
      if(idpr < 1) idpr = 1
      if(idpr > dprdata%n1c21) idpr = dprdata%n1c21
      if(jdpr < 1) jdpr = 1
      if(jdpr > 49) jdpr = 49
      !print*, idpr,jdpr
      
      gmiprof(j)%anc_tskin = dprdata%envSknTemp(jdpr,idpr)
      gmiprof(j)%anc_t2m = dprdata%envSfcTemp(jdpr,idpr)
      gmiprof(j)%anc_psfc = dprdata%envSfcPress(jdpr,idpr)
      gmiprof(j)%anc_u10m = dprdata%envSfcWindU(jdpr,idpr)
      gmiprof(j)%anc_v10m = dprdata%envSfcWindV(jdpr,idpr)
      gmiprof(j)%sfc_elev = 1e-3*dprdata%elevation(jdpr,idpr)
      gmiprof(j)%anc_clwp = 0.25*cos(dtor*dprdata%localzenithangle(jdpr,idpr))*sum(dprdata%envCloud(:,jdpr,idpr))
      gmiprof(j)%tskin_index = nint(gmiprof(j)%anc_tskin)
      
      gmiprof(j)%tskin = -99.
      gmiprof(j)%t2m = -99.
      gmiprof(j)%psfc = -99.
      gmiprof(j)%w10m = -99.

      
      !interpolate DPR ENV profiles to MERRA2 pressure levels
      gmiprof(j)%press_lev = MERRA_PRESSLEV
      gmiprof(j)%anc_temp_lev(:) = -99
      gmiprof(j)%anc_clw_lev(:) = 0.
      !print*, dprdata%ngates
      !print*, gmiprof(j)%anc_psfc, gmiprof(j)%anc_sfc_elev, dprdata%binRealSurface(jdpr,idpr)
      do kg=1,MERRA_NLEV
      	do kp=1,dprdata%ngates-1
      	  if((dprdata%envPress(kp,jdpr,idpr) .le. gmiprof(j)%press_lev(kg)) .and. (dprdata%envPress(kp+1,jdpr,idpr) .gt. gmiprof(j)%press_lev(kg))) then
      	  	dp = (gmiprof(j)%press_lev(kg)-dprdata%envPress(kp,jdpr,idpr))/(dprdata%envPress(kp+1,jdpr,idpr)-dprdata%envPress(kp,jdpr,idpr))
      	  	
      	  	gmiprof(j)%anc_temp_lev(kg) = (1-dp)*dprdata%envTemp(kp,jdpr,idpr)+dp*dprdata%envTemp(kp+1,jdpr,idpr)
      	  	!convert kg/m^3 to kg/kg mixing ratio
      	  	rho = 100*gmiprof(j)%press_lev(kg)/(Rd*gmiprof(j)%anc_temp_lev(kg))
      	  	gmiprof(j)%anc_qv_lev(kg) = 1e-3*((1-dp)*dprdata%envQv(kp,jdpr,idpr)+dp*dprdata%envQv(kp+1,jdpr,idpr))/rho
      	  	gmiprof(j)%anc_clw_lev(kg) = 1e-3*((1-dp)*dprdata%envCloud(kp,jdpr,idpr)+dp*dprdata%envCloud(kp+1,jdpr,idpr))/rho
      	  	!gmiprof(j)%anc_hgt_lev(kg) = 250*cos(DTOR*dprdata%localZenithAngle(jdpr,idpr))*(dprdata%binRealSurface(jdpr,idpr)-kp-dp+1)+dprdata%elevation(jdpr,idpr)
      	  	!print '(2I5,2F8.2,2F8.5,F10.1)', kg, kp,gmiprof(j)%press_lev(kg),gmiprof(j)%anc_temp_lev(kg), gmiprof(j)%anc_qv_lev(kg), gmiprof(j)%anc_clw_lev(kg)!, gmiprof(j)%anc_hgt_lev(kg)
      	  endif
      	end do
      	!if level has not been interpolated because it is above the top of the DPR range, use nearest temp and QV and extrapolate height hydrostatically
      	if((gmiprof(j)%anc_temp_lev(kg) .eq. -99) .and. (gmiprof(j)%press_lev(kg) .lt. dprdata%envPress(1,jdpr,idpr))) then
      	  gmiprof(j)%anc_temp_lev(kg) = dprdata%envTemp(1,jdpr,idpr)
      	  rho = 100*gmiprof(j)%press_lev(kg)/(Rd*gmiprof(j)%anc_temp_lev(kg))
      	  gmiprof(j)%anc_qv_lev(kg) = 1e-3*dprdata%envQv(1,jdpr,idpr)/rho
      	  gmiprof(j)%anc_clw_lev(kg) = 1e-3*dprdata%envCloud(1,jdpr,idpr)/rho
      	  !gmiprof(j)%anc_hgt_lev(kg) = gmiprof(j)%anc_hgt_lev(kg-1) + (gmiprof(j)%press_lev(kg-1)-gmiprof(j)%press_lev(kg))*287*gmiprof(j)%anc_temp_lev(kg)/(gmiprof(j)%press_lev(kg)*9.81)
      	  !print '(2I5,2F8.2,2F8.5,F10.1)', kg, kp,gmiprof(j)%press_lev(kg),gmiprof(j)%anc_temp_lev(kg), gmiprof(j)%anc_qv_lev(kg), gmiprof(j)%anc_clw_lev(kg)!, gmiprof(j)%anc_hgt_lev(kg)
      	 
      	endif
      end do
      
      !print*, MERRA_PRESSLEV
      !print*, dprdata%envPress(:,jdpr,idpr)
      call getwfraction(gmiprof(j)%lat, gmiprof(j)%lon, wf)
      !get sfc type from my new atlas
      !print*, LUT%emis_map_bare(:,329,389)
  
      call get_gmi_sfc_class(gmiprof(j)%lat, gmiprof(j)%lon, gmiprof(j)%sfc_type, gmiprof(j)%ice_type)
      
      !Consider up to 3 sfc types based on wfract, temp, and sfc type
      !initilize gmi_ret chisq since these will be used to judge best retrieval
      gmi_ret_ocean%chi_squared = 1e7
      gmi_ret_land%chi_squared = 1e7
      gmi_ret_ice%chi_squared = 1e7
      !Conditions for ocean: wfract = 100 or sfc_type = 0
      if(((wf .gt. 99.99) .or. (gmiprof(j)%sfc_type .eq. 0))) then
        gmidata_temp = gmiprof(j)
        gmidata_temp%sfc_type=0 !set flag for RTM
        call gmi_ocean_ret_noclw(gmidata_temp, gmi_ret_ocean)
        if(gmi_ret_ocean%chi_squared .gt. 1.) then
		  gmidata_temp = gmiprof(j)
		  gmidata_temp%sfc_type=0 !set flag for RTM
          call gmi_ocean_ret_clw1d(gmidata_temp, gmi_ret_ocean)
		  !gmi_ret(i) = gmi_ret2(i)
        endif
      endif
      !conditions for bare ground: wfract < 100 or sfc type != 0 (use class 1 - coast - if sfc type=0)
      if((wf .lt. 99.99) .or. (gmiprof(j)%sfc_type .ne. 0)) then
      	gmidata_temp = gmiprof(j)
      	if(gmidata_temp%sfc_type .eq. 0) gmidata_temp%sfc_type=1
      	!print*, 'land', gmidata_temp%sfc_type, gmidata_temp%ice_type
      	call gmi_land_ret(gmidata_temp, gmi_ret_land)
      	!print*, LUT%emis_map_bare(:,329,389)
      	!call get_gmi_emis_std(gmidata_temp%lat, gmidata_temp%lon, gmidata_temp%sfc_type, gmi_ret_land%emis(1:11), gmi_ret_land%semis(1:11))
      	!print*, gmidata_temp%lat, gmidata_temp%lon, gmidata_temp%sfc_type, gmi_ret_land%emis(1:11), gmi_ret_land%semis(1:11)
      	!gmi_ret_land%chi_squared = 1
      endif
      !conditions for snow/sea ice: ice type !=0 and tskin < 275
      if((gmiprof(j)%ice_type .ne. 0) .and. (gmiprof(j)%tskin .lt. 275)) then
        gmidata_temp = gmiprof(j)
        !print*, 'ice', gmidata_temp%sfc_type, gmidata_temp%ice_type
        gmidata_temp%sfc_type = gmidata_temp%ice_type
      	call gmi_land_ret(gmidata_temp, gmi_ret_ice)
      endif
      i2 = i-minGMIscan(1)+1
      j2 = j-minGMIpix(1)+1
      !find best solution and save to gmi_ret_array
      !print*, i2,j2, gmi_ret_ocean%chi_squared, gmi_ret_land%chi_squared, gmi_ret_ice%chi_squared
	  if((gmi_ret_ocean%chi_squared .le. gmi_ret_land%chi_squared) .and. (gmi_ret_ocean%chi_squared .le. gmi_ret_ice%chi_squared)) then
	    gmi_ret(j2,i2) = gmi_ret_ocean
	    gmiprof(j)%sfc_type = 0
	    gmi_ret(j2,i2)%sfc_type=0
	  else if((gmi_ret_land%chi_squared .lt. gmi_ret_ocean%chi_squared) .and. (gmi_ret_land%chi_squared .le. gmi_ret_ice%chi_squared)) then
	  	gmi_ret(j2,i2) = gmi_ret_land
	  	gmi_ret(j2,i2)%sfc_type = gmiprof(j)%sfc_type
	  else
	    gmi_ret(j2,i2) = gmi_ret_ice
	    gmi_ret(j2,i2)%sfc_type = gmiprof(j)%ice_type
	  endif
      gmi_ret(j2,i2)%chisq_ocean = gmi_ret_ocean%chi_squared
      gmi_ret(j2,i2)%chisq_land = gmi_ret_land%chi_squared
      gmi_ret(j2,i2)%chisq_ice = gmi_ret_ice%chi_squared
    end do
    !$OMP END DO
    !$OMP END PARALLEL  
   
    do j=minGMIpix(1),maxGMIpix(1)
      i2 = i-minGMIscan(1)+1
      j2 = j-minGMIpix(1)+1
      gmi_wgt(j2,i2) = 1-alog10(gmi_ret(j2,i2)%chi_squared)
      if(gmi_wgt(j2,i2) .gt. 1-alog10(0.5)) gmi_wgt(j2,i2) = 1-alog10(0.5)
      if(gmi_wgt(j2,i2) .lt. 0) gmi_wgt(j2,i2) = 0.

      !write(12), i,j,gmiprof(j)%lat, gmiprof(j)%lon, gmiprof(j)%anc_tskin, gmiprof(j)%anc_t2m, &
      !           gmiprof(j)%anc_psfc, gmiprof(j)%anc_u10m, gmiprof(j)%anc_v10m, gmiprof(j)%sfc_elev, &
      !           gmiprof(j)%anc_clwp, wf, gmiprof(j)%sfc_type, gmiprof(j)%ice_type, gmi_ret(j2,i2)%chisq_ocean, gmi_ret(j2,i2)%chisq_land, gmi_ret(j2,i2)%chisq_ice, &
      !           gmiprof(j)%gmi_tb, gmi_ret(j2,i2)%sim_tb0, gmi_ret(j2,i2)%sim_tb, &
      !           gmi_ret(j2,i2)%w10m_init, gmi_ret(j2,i2)%relAz_init, gmi_ret(j2,i2)%tpw_init, gmi_ret(j2,i2)%tskin_init, gmi_ret(j2,i2)%clwp_init, &
      !           gmi_ret(j2,i2)%w10m_ret, gmi_ret(j2,i2)%relAz_ret, gmi_ret(j2,i2)%tpw_ret, gmi_ret(j2,i2)%tskin_ret, gmi_ret(j2,i2)%clwp_ret, &
      !           gmi_ret(j2,i2)%emis, gmi_ret(j2,i2)%semis, gmi_ret(j2,i2)%aemis, gmi_wgt(j2,i2), gmi_ret(j2,i2)%stpw_ret
            
    end do

  end do
   
  

  !close(12)
  !stop
  !interplate gmi retrievals to DPR and over-write ENV profiles
  
  !parameters to replace (create alternate versions of these):
  !dprdata%envSknTemp(jdpr,idpr) - done
  !dprdata%envSfcTemp(jdpr,idpr) - done
  !dprdata%envSfcPress(jdpr,idpr) - don't change
  !dprdata%envSfcWindU(jdpr,idpr) - just scale (done)
  !dprdata%envSfcWindV(jdpr,idpr) - scale (done)
  !dprdata%envTemp(kp,jdpr,idpr) - interpolate vertically (done)
  !dprdata%envQv(kp,jdpr,idpr) - interpolate and convert units (done)
  !dprdata%envCloud(kp,jdpr,idpr) - interpolate (done)
  !additional output - copy to appropriate fields and ensure output is written in DPR routines
  !NS/MS/FS airTemperature (done)
  !NS/MS/FS cloudLiqWaterCont (done)
  !NS/MS/FS cloudLiqSigma (done)
  !NS/MS/FS columnVaporSigma (done)
  !NS/MS/FS errorOfDataFit (chisq) - this will be a strict nearest neighbor value so Sarah can use to identify light precip. (done)
  !NS/MS/FS pia (done)
  !NS/MS/FS simulatedBrightTemp (done, but needs separate structure to avoid conflict from radar retrieval)
  !NS/MS/FS skinTempSigma
  !NS/MS/FS skinTemperature (done)
  !NS/MS/FS surfEmissSigma
  !NS/MS/FS surfEmissivity
  !NS/MS/FS surfaceAirTemperature
  !NS/MS/FS surfaceVaporDensity
  !NS/MS/FS tenMeterWindSigma
  !NS/MS/FS tenMeterWindSpeed
  !NS/MS/FS vaporDensity
  
  !parameters for Bill to add
  !NS/MS/FS columnWaterVapor
  !emissivityAvgKernel
  !total liqWaterPath
  !total IceWaterPath
  
  !retrieved emissivity, observed and simulated Tb
  !Use 10 km radius, weight by 1-log10(chisq), cap at 1-log10(0.5)
  !if no good obs found, double search radius but weight env data linearly, with weight decreasing from 1@10km to 0@40km
  print*, 'Interpolating GMI 1DVAR'
  !open(unit=12,file='temp/dprout.bin',form='unformatted',status='unknown')
  do idpr=1,dprdata%n1c21
    !print*, idpr
    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP DO PRIVATE (jdpr,gmi_dist2, gloc, maxdist, gmi_mask, env_wgt, kg, gmi_temp_prof, gmi_qv_prof, gmi_cloud_prof, kp, dp, rho, clwp)
  	do jdpr=1,49
  	  !calc dist to all gmi pixels in this chunk
  	  gmi_dist2 = 6371.*ACOS(SIN(dprdata%xlat(jdpr,idpr)*DTOR)*SIN(gmidata%S1lat(minGMIpix(1):maxGMIpix(1),minGMIscan(1):maxGMIscan(1))*DTOR) + &
  	  COS(dprdata%xlat(jdpr,idpr)*DTOR)*COS(gmidata%S1lat(minGMIpix(1):maxGMIpix(1),minGMIscan(1):maxGMIscan(1))*DTOR)*COS(dprdata%xlon(jdpr,idpr)*DTOR-gmidata%S1lon(minGMIpix(1):maxGMIpix(1),minGMIscan(1):maxGMIscan(1))*DTOR))
  	  gloc = minloc(gmi_dist2)
  	  !print*, idpr, jdpr, gloc
  	  !Copy nearest GMI to optEst variables w/o filtering or interpolation
  	  !Three dimensional fields
  	  do kp=1,dprdata%ngates
  	    do kg=2,MERRA_NLEV
      	  if((MERRA_PRESSLEV(kg) .le. dprdata%envPress(kp,jdpr,idpr)) .and. (MERRA_PRESSLEV(kg-1) .gt. dprdata%envPress(kp,jdpr,idpr))) then
      	  	dp = (dprdata%envPress(kp,jdpr,idpr)-MERRA_PRESSLEV(kg-1))/(MERRA_PRESSLEV(kg)-MERRA_PRESSLEV(kg-1))
      	  	!print '(4F8.2)', dprdata%envPress(kp,jdpr,idpr), MERRA_PRESSLEV(kg-1:kg), dp
      	  	dprdata%OETemp(kp,jdpr,idpr) = (1-dp)*gmi_ret(gloc(1),gloc(2))%temp_lev(kg-1)+dp*gmi_ret(gloc(1),gloc(2))%temp_lev(kg)
      	  	!convert kg/m^3 to kg/kg mixing ratio
      	  	rho = 100*dprdata%envPress(kp,jdpr,idpr)/(Rd*dprdata%OETemp(kp,jdpr,idpr))
      	  	dprdata%OEQv(kp,jdpr,idpr) = 1e3*rho*((1-dp)*gmi_ret(gloc(1),gloc(2))%qv_lev(kg-1)+dp*gmi_ret(gloc(1),gloc(2))%qv_lev(kg))
      	  	dprdata%OECloud(kp,jdpr,idpr) = ((1-dp)*gmi_ret(gloc(1),gloc(2))%cloud_water(kg-1)+dp*gmi_ret(gloc(1),gloc(2))%cloud_water(kg))
      	    !print '(2I5,3F8.2,4F12.5)', kg, kp,dprdata%envPress(kp,jdpr,idpr), dprdata%modTemp(kp,jdpr,idpr), dprdata%envTemp(kp,jdpr,idpr), &
      	    !  dprdata%modQv(kp,jdpr,idpr), dprdata%envQv(kp,jdpr,idpr), dprdata%modCloud(kp,jdpr,idpr), dprdata%envCloud(kp,jdpr,idpr)
      	  endif
      	end do
  	  end do
  	  !calc sfc qv from t2m and lowest sfc Q
  	  rho = 100*dprdata%envSfcPress(jdpr,idpr)/(Rd*gmi_ret(gloc(1), gloc(2))%t2m_ret)
  	  !print*, rho, 1e3*rho*gmi_ret(gloc(1),gloc(2))%qv_lev(1)
  	  
  	  !Two dimensional fields
  	  dprdata%OEcloudLiqPath(jdpr,idpr) = gmi_ret(gloc(1),gloc(2))%clwp_ret
  	  if(clwp .gt. 0) then
  	  	dprdata%OEcloudLiqSigma(jdpr,idpr) = exp(gmi_ret(gloc(1),gloc(2))%sclwp_ret+alog(dprdata%OEcloudLiqPath(jdpr,idpr)))-dprdata%OEcloudLiqPath(jdpr,idpr)!/clwp-1
  	  else
  	  	dprdata%OEcloudLiqSigma(jdpr,idpr) = missing_r4
  	  endif
  	  dprdata%OEcloudLiqSigma(jdpr,idpr) = gmi_ret(gloc(1), gloc(2))%sclwp_ret
      dprdata%OESfcWind(jdpr,idpr) = gmi_ret(gloc(1), gloc(2))%w10m_ret
      dprdata%OESfcWindSigma(jdpr,idpr) = gmi_ret(gloc(1), gloc(2))%sw10m_ret
      dprdata%OESknTemp(jdpr,idpr) = gmi_ret(gloc(1), gloc(2))%tskin_ret
      dprdata%OEskinTempSigma(jdpr,idpr) = gmi_ret(gloc(1), gloc(2))%stskin_ret
      dprdata%OESfcTemp(jdpr,idpr) = gmi_ret(gloc(1), gloc(2))%t2m_ret
      dprdata%OESfcQv(jdpr,idpr) = 1e3*rho*gmi_ret(gloc(1),gloc(2))%qv_lev(1)
	  dprdata%OEtpw(jdpr,idpr) = gmi_ret(gloc(1), gloc(2))%tpw_ret
	  dprdata%OEtpwSigma(jdpr,idpr) = gmi_ret(gloc(1), gloc(2))%stpw_ret
  	  dprdata%OEchiSq(jdpr,idpr) = gmi_ret(gloc(1), gloc(2))%chi_squared
  	  dprdata%OEstype(jdpr,idpr) = gmi_ret(gloc(1), gloc(2))%sfc_type
  	  
  	  !if(dprdata%skintempSigma(jdpr,idpr) < 0) dprdata%skintempSigma(jdpr,idpr) = missing_r4
  	  !Simulated Observations
  	  do kg=1,13
  	  	dprdata%OEsimTbNonRain(kg,jdpr,idpr) = gmi_ret(gloc(1),gloc(2))%sim_tb(kg)
  	  	dprdata%OEemis(kg,jdpr,idpr) = gmi_ret(gloc(1),gloc(2))%emis(kg)
  	  	dprdata%OEemisSigma(kg,jdpr,idpr) = gmi_ret(gloc(1),gloc(2))%semis(kg)
  	  	dprdata%OEemisA(kg,jdpr,idpr) = gmi_ret(gloc(1),gloc(2))%aemis(kg)
  	  end do
  	  
  	  dprdata%OEpiaNonRain(1,jdpr,idpr) = gmi_ret(gloc(1),gloc(2))%PIAKu
  	  dprdata%OEpiaNonRain(2,jdpr,idpr) = gmi_ret(gloc(1),gloc(2))%PIAKa

  	  !set up gmi mask for interpolation
  	  !print*, gloc,dprdata%chiSq(jdpr,idpr) 
  	  maxdist = 0
  	  gmi_mask(:,:) = 0.
  	  do while((maxval(gmi_mask*gmi_wgt) .eq. 0) .and. (maxdist .lt. 50))
  	    maxdist=maxdist+7.5
  	    where(gmi_dist2 .lt. maxdist)
  	      gmi_mask = 1
  	    elsewhere
  	      gmi_mask = 0
  	    end where
  	    
  	  end do
  	  
  	  !Do filtering/interpolation for modified environmental variables
  	  
  	  !print '(2I5,2F8.2)', idpr, jdpr, maxdist, maxval(gmi_wgt)
  	  env_wgt = (maxdist-10)/40
  	  if(env_wgt .lt. 0) env_wgt = 0
  	  if(env_wgt .ge. 1) then
  	    env_wgt = 1
  	    dprdata%modSknTemp(jdpr,idpr) = dprdata%envSknTemp(jdpr,idpr)
  	  	dprdata%modSfcTemp(jdpr,idpr) = dprdata%envSfcTemp(jdpr,idpr)
  	  	dprdata%modSfcWind(jdpr,idpr) = dprdata%envSfcWind(jdpr,idpr)
  	  else
  	  	dprdata%modSknTemp(jdpr,idpr) = (1-env_wgt)*sum(gmi_ret%tskin_ret*gmi_mask*gmi_wgt)/sum(gmi_mask*gmi_wgt)+env_wgt*dprdata%envSknTemp(jdpr,idpr)
  	  	dprdata%modSfcTemp(jdpr,idpr) = (1-env_wgt)*sum(gmi_ret%t2m_ret*gmi_mask*gmi_wgt)/sum(gmi_mask*gmi_wgt)+env_wgt*dprdata%envSfcTemp(jdpr,idpr)
  	  	dprdata%modSfcWind(jdpr,idpr) = (1-env_wgt)*sum(gmi_ret%w10m_ret*gmi_mask*gmi_wgt)/sum(gmi_mask*gmi_wgt)+env_wgt*dprdata%envSfcWind(jdpr,idpr)
  	  endif
  	  
  	  	
  	  !Just scale 10m wind for now. azimuth correction or U/V conversion can be done when combined GMI/DPR retrieval implemented.
  	  if(dprdata%envSfcWind(jdpr,idpr) .gt. 0) then
  	    dprdata%modSfcWindU(jdpr,idpr) = dprData%envSfcWindU(jdpr,idpr)*dprdata%modSfcWind(jdpr,idpr)/dprdata%envSfcWind(jdpr,idpr)
  	    dprdata%modSfcWindV(jdpr,idpr) = dprData%envSfcWindV(jdpr,idpr)*dprdata%modSfcWind(jdpr,idpr)/dprdata%envSfcWind(jdpr,idpr)
  	  else
  	    dprdata%modSfcWindU(idpr,jdpr) = dprdata%modSfcWind(jdpr,idpr)
  	  endif
  	  !interpolate temp, cloud and water vapor profiles to pressure levels of original dpr data
  	  dprdata%modQv(:,jdpr,idpr) = dprdata%envQv(:,jdpr,idpr)
      dprdata%modTemp(:,jdpr,idpr) = dprdata%envTemp(:,jdpr,idpr)
      dprdata%modCloud(:,jdpr,idpr) = 0
      !calculate weighted average profiles from retrievals
      do kg=1,MERRA_NLEV
        gmi_temp_prof(kg) = sum(gmi_ret%temp_lev(kg)*gmi_mask*gmi_wgt)/sum(gmi_mask*gmi_wgt)
        gmi_qv_prof(kg) = sum(gmi_ret%qv_lev(kg)*gmi_mask*gmi_wgt)/sum(gmi_mask*gmi_wgt)
        gmi_cloud_prof(kg) = sum(gmi_ret%cloud_water(kg)*gmi_mask*gmi_wgt)/sum(gmi_mask*gmi_wgt)
        !print '(I5,2F8.2,2F12.5)', kg, MERRA_PRESSLEV(kg), gmi_temp_prof(kg), gmi_qv_prof(kg), gmi_cloud_prof(kg)
      end do

	  do kp=1,dprdata%ngates
  	    do kg=2,MERRA_NLEV
      	  if((MERRA_PRESSLEV(kg) .le. dprdata%envPress(kp,jdpr,idpr)) .and. (MERRA_PRESSLEV(kg-1) .gt. dprdata%envPress(kp,jdpr,idpr))) then
      	  	dp = (dprdata%envPress(kp,jdpr,idpr)-MERRA_PRESSLEV(kg-1))/(MERRA_PRESSLEV(kg)-MERRA_PRESSLEV(kg-1))
      	  	!print '(4F8.2)', dprdata%envPress(kp,jdpr,idpr), MERRA_PRESSLEV(kg-1:kg), dp
      	  	dprdata%modTemp(kp,jdpr,idpr) = (1-dp)*gmi_temp_prof(kg-1)+dp*gmi_temp_prof(kg)
      	  	!convert kg/m^3 to kg/kg mixing ratio
      	  	rho = 100*dprdata%envPress(kp,jdpr,idpr)/(Rd*dprdata%modTemp(kp,jdpr,idpr))
      	  	dprdata%modQv(kp,jdpr,idpr) = 1e3*rho*((1-dp)*gmi_qv_prof(kg-1)+dp*gmi_qv_prof(kg))
      	  	dprdata%modCloud(kp,jdpr,idpr) = ((1-dp)*gmi_cloud_prof(kg-1)+dp*gmi_cloud_prof(kg))
      	    !print '(2I5,3F8.2,4F12.5)', kg, kp,dprdata%envPress(kp,jdpr,idpr), dprdata%modTemp(kp,jdpr,idpr), dprdata%envTemp(kp,jdpr,idpr), &
      	    !  dprdata%modQv(kp,jdpr,idpr), dprdata%envQv(kp,jdpr,idpr), dprdata%modCloud(kp,jdpr,idpr), dprdata%envCloud(kp,jdpr,idpr)
      	  endif
      	end do
  	  end do
  	  !now weight by distance to good retrieval
  	  
  	  !parameters in retdata structure
  	  if(env_wgt .lt. 1) then
  	    clwp = sum(gmi_ret%clwp_ret*gmi_mask*gmi_wgt)/sum(gmi_mask*gmi_wgt)
  	    dprdata%modTemp(:,jdpr,idpr) = (1-env_wgt)*dprdata%modTemp(:,jdpr,idpr) + env_wgt*dprdata%envTemp(:,jdpr,idpr)
  	  	dprdata%modQv(:,jdpr,idpr) = (1-env_wgt)*dprdata%modQv(:,jdpr,idpr) + env_wgt*dprdata%envQv(:,jdpr,idpr)
  	  	dprdata%modCloud(:,jdpr,idpr) = (1-env_wgt)*dprdata%modCloud(:,jdpr,idpr) + env_wgt*dprdata%envCloud(:,jdpr,idpr)
  	  else
  	    dprdata%modTemp(:,jdpr,idpr) = dprdata%envTemp(:,jdpr,idpr)
  	  	dprdata%modQv(:,jdpr,idpr) = dprdata%envQv(:,jdpr,idpr)
  	  	dprdata%modCloud(:,jdpr,idpr) = dprdata%envCloud(:,jdpr,idpr)
  	  endif
  	 
  	  
  	end do
  	!$OMP END DO
    !$OMP END PARALLEL  
    
  	!do jdpr=1,49
  	!  write(12), idpr, jdpr, dprdata%xlon(jdpr,idpr), dprdata%xlat(jdpr,idpr), maxdist, sum(gmi_mask*gmi_wgt), dprdata%modSknTemp(jdpr,idpr), &
  	!                        dprdata%modSfcTemp(jdpr,idpr), dprdata%modSfcWind(jdpr,idpr), 0.25*cos(dtor*dprdata%localzenithangle(jdpr,idpr))*sum(dprdata%envCloud(:,jdpr,idpr)), &
  	!                        dprdata%cloudLiqSigma(jdpr,idpr), dprdata%tpwSigma(jdpr,idpr)
  	!end do
  end do
  !close(12)
  deallocate(gmi_ret)
  deallocate(gmi_dist2)
  deallocate(gmi_wgt)
  
  !deallocate(dprSknTemp,dprSfcTemp,dprSfcWind)
  
end subroutine gmiretsub

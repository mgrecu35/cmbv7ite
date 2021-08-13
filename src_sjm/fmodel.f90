

subroutine get_state_land_rh(nvar,neof,nemis,emis_ch,nlev,x,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata)

use profile_def

implicit none

include 'parametersMERRA.inc'
include 'constants_sjm.inc'

integer, intent(in) :: nvar,neof,nemis,nlev,emis_ch(11)
real, intent(in) :: x(nvar), xbar(nvar), xstd(nvar)
real, intent(in) :: prof_avg(1+2*MERRA_NLEV), prof_std(1+2*MERRA_NLEV), prof_W(1+2*MERRA_NLEV)
real, intent(in) :: prof_eofs(1+2*MERRA_NLEV,1+2*MERRA_NLEV)
type(profile), intent(inout) :: gdata

integer :: i, j, i0
real :: Tavg, Qavg, dz, R, clwp, tpw, wgt, rhov, rhod, es, pv!, sclT, sclQ
real :: rh_lev(MERRA_NLEV)

!print*, nlev
!initialize with average profile
!print*, prof_std



!if sfc pressure is < 500mb, just use partial profiles
i0=0
! print*, gdata%psfc
do i=1,MERRA_NLEV
  if(gdata%psfc .ge. MERRA_PRESSLEV(i)) exit
  i0=i0+1
end do

if(gdata%anc_tskin .gt. 0.) then
  gdata%tskin = gdata%anc_tskin
else
  print*, 'ancillary data missing for tskin'
  stop
  !gdata%tskin = prof_avg(1)
endif

if(gdata%anc_t2m .gt. 0.) then
  gdata%t2m = gdata%anc_t2m
else
  print*, 'ancillary data missing for t2m'
  stop
endif

do i=1,nlev
  if(gdata%anc_temp_lev(i+MERRA_NLEV-nlev) .ge. 0.) then
    gdata%temp_lev(i) = gdata%anc_temp_lev(i+MERRA_NLEV-nlev)
    !gdata%temp_lev(i) = prof_avg(i+1+i0)
  else if(i .eq. 1) then
    gdata%temp_lev(i) = gdata%t2m
  else if(i .gt. 1) then
    gdata%temp_lev(i) = gdata%temp_lev(i-1)
  else
    print*, 'ancillary data missing for temp lev',i
    stop
  endif
  gdata%press_lev(i) = MERRA_presslev(i+MERRA_NLEV-nlev)
  if(gdata%anc_qv_lev(i+MERRA_NLEV-nlev) .ge. 0.) then
    gdata%qv_lev(i) = gdata%anc_qv_lev(i+MERRA_NLEV-nlev)
    !calculate MERRA2 RH
  else if(gdata%anc_qv_lev(i+MERRA_NLEV-nlev+1) .ge. 0) then
    gdata%anc_qv_lev(i+MERRA_NLEV-nlev) = gdata%anc_qv_lev(i+MERRA_NLEV-nlev+1)
  else 
     print*, 'ancillary data missing for Qv lev',i, gdata%anc_qv_lev(:), gdata%anc_temp_lev(:)
    stop
  endif
  call get_es(gdata%temp_lev(i),es)
  pv = (gdata%qv_lev(i)*100*gdata%press_lev(i)*Rv)/(Rd*(1+(gdata%qv_lev(i)*Rv)/Rd))
  !rhov = (100*es)/(Rv*gdata%temp_lev(i))
  !rhod = (100*MERRA_PRESSLEV(i)-100*es)/(Rd*gdata%temp_lev(i))
  rh_lev(i) = pv/(100*es)!gdata%qv_lev(i)/(rhov/rhod)
end do
!print*, nlev, i0
!print*, gdata%temp_lev(1:nlev)
!print '(10F8.4)', prof_std(:)
!print '(10F8.4)', prof_std(2+i0:1+i0+nlev)
!print '(10F8.4)', prof_std(2+i0+MERRA_NLEV:1+i0+MERRA_NLEV+nlev)
!stop
!if(i0+1+MERRA_NLEV+nlev .gt. 55) print*, i0, nlev, gdata%psfc
!add eofs
do i=1,neof
  !print*, i, x(i), xstd(i), xbar(i)
  gdata%t2m = gdata%t2m+prof_eofs(i,1)*(x(i)*xstd(i)+xbar(i))*prof_std(1)
  gdata%temp_lev(1:nlev) = gdata%temp_lev(1:nlev)+prof_eofs(i,2+i0:2+i0+nlev)*(x(i)*xstd(i)+xbar(i))*prof_std(2+i0:2+i0+nlev)
  rh_lev(1:nlev) = rh_lev(1:nlev)+prof_eofs(i,2+i0+MERRA_NLEV:1+i0+MERRA_NLEV+nlev)*(x(i)*xstd(i)+xbar(i))*prof_std(2+i0+MERRA_NLEV:1+i0+MERRA_NLEV+nlev)
end do

!to do: check saturated/negative RH and convert back to Qv
do i=1,nlev
  if(rh_lev(i) .lt. 0.) rh_lev(i) = 0.
  if(rh_lev(i) .gt. 1.01) rh_lev(i) = 1.
  call get_es(gdata%temp_lev(i),es)
  rhov = (rh_lev(i)*es)/(Rv*gdata%temp_lev(i))
  rhod = (gdata%press_lev(i)-rh_lev(i)*es)/(Rd*gdata%temp_lev(i))
  gdata%qv_lev(i) = rhov/rhod
  !print*, gdata%temp_lev(i), rh_lev(i), gdata%qv_lev(i), rhov/rhod
end do
!stop

!print*, gdata%qv_lev
!calculate height levels
Tavg = 0.5*(gdata%t2m+gdata%temp_lev(1))
Qavg = gdata%qv_lev(1)
R = (1.-Qavg)*Rd+Qavg*Rv
dz = R*Tavg/9.81*alog(gdata%psfc/gdata%press_lev(1))
gdata%hgt_lev(1) = gdata%sfc_elev+0.001*dz
!print*, 1, Tavg, Qavg, gdata%press_lev(1), gdata%hgt_lev(1)
do i=2,nlev
  Tavg = 0.5*(gdata%temp_lev(i-1)+gdata%temp_lev(i))
  Qavg = 0.5*(gdata%qv_lev(i-1)+gdata%qv_lev(i-1))
  R = (1.-Qavg)*Rd+Qavg*Rv
  dz = R*Tavg/9.81*alog(gdata%press_lev(i-1)/gdata%press_lev(i))
  gdata%hgt_lev(i) = gdata%hgt_lev(i-1)+0.001*dz
  !print*, i, Tavg, Qavg, gdata%press_lev(i), gdata%hgt_lev(i), gdata%anc_hgt_lev(i+MERRA_NLEV-nlev)
end do

!emissivity
gdata%emis = -99.

do i=1,nemis
  !need to know which channel corresponds to which variable
  j = nvar-(nemis+1)+i
  if(emis_ch(i) .eq. 1) gdata%emis(1) = x(j)*xstd(j)+xbar(j) !10V
  if(emis_ch(i) .eq. 2) gdata%emis(2) = x(j)*xstd(j)+xbar(j) !10H
  if(emis_ch(i) .eq. 3) gdata%emis(3) = x(j)*xstd(j)+xbar(j) !18V
  if(emis_ch(i) .eq. 4) gdata%emis(4) = x(j)*xstd(j)+xbar(j) !18H
  if(emis_ch(i) .eq. 5) gdata%emis(5) = x(j)*xstd(j)+xbar(j) !23V
  if(emis_ch(i) .eq. 6) gdata%emis(6) = x(j)*xstd(j)+xbar(j) !36V
  if(emis_ch(i) .eq. 7) gdata%emis(7) = x(j)*xstd(j)+xbar(j) !36H
  if(emis_ch(i) .eq. 8) gdata%emis(8) = x(j)*xstd(j)+xbar(j) !89V
  if(emis_ch(i) .eq. 9) gdata%emis(9) = x(j)*xstd(j)+xbar(j) !89H
  if(emis_ch(i) .eq. 10) gdata%emis(10) = x(j)*xstd(j)+xbar(j) !166V
  if(emis_ch(i) .eq. 11) gdata%emis(11) = x(j)*xstd(j)+xbar(j) !166H
end do
!if(gdata%emis(6) .gt. 0. .and. gdata%emis(3) .gt. 0. .and. gdata%emis(5) .eq. -99.) gdata%emis(5) = 0.28*gdata%emis(6)+(1.-0.28)*gdata%emis(3) !interpolate 23V
wgt = x(nvar)*xstd(nvar)+xbar(nvar)
!print*, wgt
gdata%emis(5) = wgt*gdata%emis(6)+(1.-wgt)*gdata%emis(3) !interpolate 23V
if(gdata%emis(10) .gt. 0.) gdata%emis(12:13) = gdata%emis(10) !183 +/-3 and 183 +/-7 same as 166V

gdata%cloud_water = 0.
gdata%rain_water = 0.
gdata%cloud_ice = 0.
!clwp = exp(x(nvar)*xstd(nvar)+xbar(nvar))
!print*, x(nvar), xstd(nvar), xbar(nvar), clwp
!call get_clw_profile(gdata,nlev,clwp)

end subroutine get_state_land_rh


subroutine get_state_noclw_rh(nvar,neof,nlev,x,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata)

  use profile_def

  implicit none

  include 'parametersMERRA.inc'
  include 'constants_sjm.inc'

  integer, intent(in) :: nvar,nlev,neof
  real, intent(in) :: x(nvar), xbar(nvar), xstd(nvar)
  real, intent(in) :: prof_avg(1+2*MERRA_NLEV), prof_W(1+2*MERRA_NLEV), prof_std(1+2*MERRA_NLEV)
  real, intent(in) :: prof_eofs(1+2*MERRA_NLEV,1+2*MERRA_NLEV)
  type(profile), intent(inout) :: gdata

  integer :: i, j, i0
  real :: Tavg, Qavg, dz, R, clwp, tpw, es, pv, rhov, rhod, rh_lev(MERRA_NLEV)

  i0=0

  do i=1,MERRA_NLEV
    if(gdata%psfc .ge. MERRA_PRESSLEV(i)) exit
    i0=i0+1
  end do

  gdata%tskin = x(nvar)*xstd(nvar)+xbar(nvar)
  !print*, x(nvar), xstd(nvar), xbar(nvar), gdata%tskin
  
  if(gdata%anc_t2m .gt. 0.) then
    gdata%t2m = gdata%anc_t2m
    !gdata%t2m = prof_avg(1)
  else! if(prof_avg(1) .gt. 0.) then
    !gdata%t2m = prof_avg(1)
    print*, 'ancillary data missing for t2m'
    stop
  endif

  
  do i=1,nlev
    if(gdata%anc_temp_lev(i+MERRA_NLEV-nlev) .ge. 0.) then
      gdata%temp_lev(i) = gdata%anc_temp_lev(i+MERRA_NLEV-nlev)
      !gdata%temp_lev(i) = prof_avg(i+1+i0)
    else if(i .eq. 1) then
      gdata%temp_lev(i) = gdata%t2m
    else if(i .gt. 1) then
      gdata%temp_lev(i) = gdata%temp_lev(i-1)
    else
      print*, 'ancillary data missing for temp lev',i
      stop
    endif
    gdata%press_lev(i) = MERRA_presslev(i+MERRA_NLEV-nlev)
    if(gdata%anc_qv_lev(i+MERRA_NLEV-nlev) .ge. 0.) then
      gdata%qv_lev(i) = gdata%anc_qv_lev(i+MERRA_NLEV-nlev)
      !calculate MERRA2 RH
    else if(gdata%anc_qv_lev(i+MERRA_NLEV-nlev+1) .ge. 0) then
      gdata%anc_qv_lev(i+MERRA_NLEV-nlev) = gdata%anc_qv_lev(i+MERRA_NLEV-nlev+1)
    else 
      print*, 'ancillary data missing for Qv lev',i, gdata%anc_qv_lev(:), gdata%anc_temp_lev(:)
      stop
    endif
    call get_es(gdata%temp_lev(i),es)
    pv = (gdata%qv_lev(i)*100*gdata%press_lev(i)*Rv)/(Rd*(1+(gdata%qv_lev(i)*Rv)/Rd))
    !rhov = (100*es)/(Rv*gdata%temp_lev(i))
    !rhod = (100*MERRA_PRESSLEV(i)-100*es)/(Rd*gdata%temp_lev(i))
    rh_lev(i) = pv/(100*es)!gdata%qv_lev(i)/(rhov/rhod)
    !print*, i, gdata%temp_lev(i), gdata%qv_lev(i), pv, es, rh_lev(i)
  end do
  !print*, nlev, i0
  !print*, gdata%temp_lev(1:nlev)
  !print '(10F8.4)', prof_std(:)
  !print '(10F8.4)', prof_std(2+i0:1+i0+nlev)
  !print '(10F8.4)', prof_std(2+i0+MERRA_NLEV:1+i0+MERRA_NLEV+nlev)
  !stop
  !if(i0+1+MERRA_NLEV+nlev .gt. 55) print*, i0, nlev, gdata%psfc
  !add eofs
  do i=1,neof
    !print*, i, x(i), xstd(i), xbar(i)
    gdata%t2m = gdata%t2m+prof_eofs(i,1)*(x(i)*xstd(i)+xbar(i))*prof_std(1)
    gdata%temp_lev(1:nlev) = gdata%temp_lev(1:nlev)+prof_eofs(i,2+i0:2+i0+nlev)*(x(i)*xstd(i)+xbar(i))*prof_std(2+i0:2+i0+nlev)
    rh_lev(1:nlev) = rh_lev(1:nlev)+prof_eofs(i,2+i0+MERRA_NLEV:1+i0+MERRA_NLEV+nlev)*(x(i)*xstd(i)+xbar(i))*prof_std(2+i0+MERRA_NLEV:1+i0+MERRA_NLEV+nlev)
  end do

  !to do: check saturated/negative RH and convert back to Qv
  do i=1,nlev
    if(rh_lev(i) .lt. 0.) rh_lev(i) = 0.
    if(rh_lev(i) .gt. 1.01) rh_lev(i) = 1.
    call get_es(gdata%temp_lev(i),es)
    rhov = (rh_lev(i)*es)/(Rv*gdata%temp_lev(i))
    rhod = (gdata%press_lev(i)-rh_lev(i)*es)/(Rd*gdata%temp_lev(i))
    gdata%qv_lev(i) = rhov/rhod
    !print*, gdata%temp_lev(i), rh_lev(i), gdata%qv_lev(i), rhov/rhod
  end do
  !do i=1,nlev
  !  print*, i, gdata%press_lev(i), gdata%temp_lev(i), prof_avg(i+2+i0), gdata%qv_lev(i), prof_avg(2+nlev+i+i0)
  !end do
  !stop
  !print '(40F6.3)', gdata%qv_lev(1:nlev)
  !print*, gdata%qv_lev
  !calculate height levels
  Tavg = 0.5*(gdata%t2m+gdata%temp_lev(1))
  Qavg = gdata%qv_lev(1)
  R = (1.-Qavg)*Rd+Qavg*Rv
  dz = R*Tavg/9.81*alog(gdata%psfc/gdata%press_lev(1))
  gdata%hgt_lev(1) = gdata%sfc_elev+0.001*dz
  !print*, 1, Tavg, Qavg, gdata%press_lev(1), gdata%hgt_lev(1)
  do i=2,nlev
    Tavg = 0.5*(gdata%temp_lev(i-1)+gdata%temp_lev(i))
    Qavg = 0.5*(gdata%qv_lev(i-1)+gdata%qv_lev(i-1))
    R = (1.-Qavg)*Rd+Qavg*Rv
    dz = R*Tavg/9.81*alog(gdata%press_lev(i-1)/gdata%press_lev(i))
    gdata%hgt_lev(i) = gdata%hgt_lev(i-1)+0.001*dz
    !print*, i, Tavg, Qavg, gdata%press_lev(i), gdata%hgt_lev(i)
  end do

  !add wind and cloud water
  gdata%w10m = x(nvar-2)*xstd(nvar-2)+xbar(nvar-2)!gdata%iws
  if(gdata%w10m .lt. 0.) gdata%w10m = 0.

  !print*, gdata%w10m
  !relative azimuth
  gdata%relAz = x(nvar-1)*xstd(nvar-1)+xbar(nvar-1)
  
  gdata%cloud_water = 0.
  !print '(4F8.2)', gdata%psfc, gdata%tskin, gdata%t2m, gdata%w10m
  !print '(31F7.2)', gdata%temp_lev
  !print '(31F7.4)', gdata%qv_lev
  !print '(31F7.4)', gdata%cloud_water
  !print '(31F7.2)', gdata%hgt_lev
end subroutine get_state_noclw_rh

subroutine get_state_clw1d_rh(nvar,neof,nlev,x,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata)

  use profile_def

  implicit none

  include 'parametersMERRA.inc'
  include 'constants_sjm.inc'

  integer, intent(in) :: nvar,nlev,neof
  real, intent(in) :: x(nvar), xbar(nvar), xstd(nvar)
  real, intent(in) :: prof_avg(1+2*MERRA_NLEV), prof_W(1+2*MERRA_NLEV), prof_std(1+2*MERRA_NLEV)
  real, intent(in) :: prof_eofs(1+2*MERRA_NLEV,1+2*MERRA_NLEV)
  type(profile) :: gdata

  integer :: i, j, i0, istatus
  real :: Tavg, Qavg, dz, R, clwp, tpw, es, pv, rhov, rhod, rh_lev(MERRA_NLEV)

  i0=0

  do i=1,MERRA_NLEV
    if(gdata%psfc .ge. MERRA_PRESSLEV(i)) exit
    i0=i0+1
  end do

  gdata%tskin = x(neof+3)*xstd(neof+3)+xbar(neof+3)
  !print*, x(nvar), xstd(nvar), xbar(nvar), gdata%tskin
  
  if(gdata%anc_t2m .gt. 0.) then
    gdata%t2m = gdata%anc_t2m
    !gdata%t2m = prof_avg(1)
  else! if(prof_avg(1) .gt. 0.) then
    !gdata%t2m = prof_avg(1)
    print*, 'ancillary data missing for t2m'
    stop
  endif

  
  do i=1,nlev
    if(gdata%anc_temp_lev(i+MERRA_NLEV-nlev) .ge. 0.) then
      gdata%temp_lev(i) = gdata%anc_temp_lev(i+MERRA_NLEV-nlev)
      !gdata%temp_lev(i) = prof_avg(i+1+i0)
    else if(i .eq. 1) then
      gdata%temp_lev(i) = gdata%t2m
    else if(i .gt. 1) then
      gdata%temp_lev(i) = gdata%temp_lev(i-1)
    else
      print*, 'ancillary data missing for temp lev',i
      stop
    endif
    gdata%press_lev(i) = MERRA_presslev(i+MERRA_NLEV-nlev)
    if(gdata%anc_qv_lev(i+MERRA_NLEV-nlev) .ge. 0.) then
      gdata%qv_lev(i) = gdata%anc_qv_lev(i+MERRA_NLEV-nlev)
      !calculate MERRA2 RH
    else if(gdata%anc_qv_lev(i+MERRA_NLEV-nlev+1) .ge. 0) then
      gdata%anc_qv_lev(i+MERRA_NLEV-nlev) = gdata%anc_qv_lev(i+MERRA_NLEV-nlev+1)
    else 
      print*, 'ancillary data missing for Qv lev',i, gdata%anc_qv_lev(:), gdata%anc_temp_lev(:)
      stop
    endif
    call get_es(gdata%temp_lev(i),es)
    pv = (gdata%qv_lev(i)*100*gdata%press_lev(i)*Rv)/(Rd*(1+(gdata%qv_lev(i)*Rv)/Rd))
    !rhov = (100*es)/(Rv*gdata%temp_lev(i))
    !rhod = (100*MERRA_PRESSLEV(i)-100*es)/(Rd*gdata%temp_lev(i))
    rh_lev(i) = pv/(100*es)!gdata%qv_lev(i)/(rhov/rhod)
  end do
  !print*, nlev, i0
  !print*, gdata%temp_lev(1:nlev)
  !print '(10F8.4)', prof_std(:)
  !print '(10F8.4)', prof_std(2+i0:1+i0+nlev)
  !print '(10F8.4)', prof_std(2+i0+MERRA_NLEV:1+i0+MERRA_NLEV+nlev)
  !stop
  !if(i0+1+MERRA_NLEV+nlev .gt. 55) print*, i0, nlev, gdata%psfc
  !add eofs
  do i=1,neof
    !print*, i, x(i), xstd(i), xbar(i)
    gdata%t2m = gdata%t2m+prof_eofs(i,1)*(x(i)*xstd(i)+xbar(i))*prof_std(1)
    gdata%temp_lev(1:nlev) = gdata%temp_lev(1:nlev)+prof_eofs(i,2+i0:2+i0+nlev)*(x(i)*xstd(i)+xbar(i))*prof_std(2+i0:2+i0+nlev)
    rh_lev(1:nlev) = rh_lev(1:nlev)+prof_eofs(i,2+i0+MERRA_NLEV:1+i0+MERRA_NLEV+nlev)*(x(i)*xstd(i)+xbar(i))*prof_std(2+i0+MERRA_NLEV:1+i0+MERRA_NLEV+nlev)
  end do

  !to do: check saturated/negative RH and convert back to Qv
  do i=1,nlev
    if(rh_lev(i) .lt. 0.) rh_lev(i) = 0.
    if(rh_lev(i) .gt. 1.01) rh_lev(i) = 1.
    call get_es(gdata%temp_lev(i),es)
    rhov = (rh_lev(i)*es)/(Rv*gdata%temp_lev(i))
    rhod = (gdata%press_lev(i)-rh_lev(i)*es)/(Rd*gdata%temp_lev(i))
    gdata%qv_lev(i) = rhov/rhod
    !print*, gdata%temp_lev(i), rh_lev(i), gdata%qv_lev(i), rhov/rhod
  end do
  !do i=1,nlev
  !  print*, i, gdata%press_lev(i), gdata%temp_lev(i), prof_avg(i+2+i0), gdata%qv_lev(i), prof_avg(2+nlev+i+i0)
  !end do
  !stop
  !print '(40F6.3)', gdata%qv_lev(1:nlev)
  !print*, gdata%qv_lev
  !calculate height levels
  Tavg = 0.5*(gdata%t2m+gdata%temp_lev(1))
  Qavg = gdata%qv_lev(1)
  R = (1.-Qavg)*Rd+Qavg*Rv
  dz = R*Tavg/9.81*alog(gdata%psfc/gdata%press_lev(1))
  gdata%hgt_lev(1) = gdata%sfc_elev+0.001*dz
  !print*, 1, Tavg, Qavg, gdata%press_lev(1), gdata%hgt_lev(1)
  do i=2,nlev
    Tavg = 0.5*(gdata%temp_lev(i-1)+gdata%temp_lev(i))
    Qavg = 0.5*(gdata%qv_lev(i-1)+gdata%qv_lev(i-1))
    R = (1.-Qavg)*Rd+Qavg*Rv
    dz = R*Tavg/9.81*alog(gdata%press_lev(i-1)/gdata%press_lev(i))
    gdata%hgt_lev(i) = gdata%hgt_lev(i-1)+0.001*dz
    !print*, i, Tavg, Qavg, gdata%press_lev(i), gdata%hgt_lev(i)
  end do

  !add wind and cloud water
  gdata%w10m = x(neof+1)*xstd(neof+1)+xbar(neof+1)!gdata%iws
  if(gdata%w10m .lt. 0.) gdata%w10m = 0.

  !print*, gdata%w10m
  !relative azimuth
  gdata%relAz = x(neof+2)*xstd(neof+2)+xbar(neof+2)
  
  clwp = exp(x(neof+4)*xstd(neof+4)+xbar(neof+4))
  !print*, clwp
  gdata%cloud_water = 0.
  if(gdata%anc_clwp .gt. 0.) then
    do i=1,nlev
      if(not(isnan(gdata%anc_clw_lev(i+MERRA_NLEV-nlev)))) gdata%cloud_water(i) = 1e3*gdata%anc_clw_lev(i+MERRA_NLEV-nlev)/gdata%anc_clwp*clwp
      !print*, i, gdata%temp_lev(i), gdata%hgt_lev(i), gdata%press_lev(i), gdata%anc_clw_lev(i), gdata%cloud_water(i)
    end do
  else
     call check_nan(gdata%hgt_lev,gdata%temp_lev,istatus)
     if(istatus.eq.1) then
        print*, xbar
        print*, xstd
        print*, x
     endif
     call get_clw_profile(gdata,nlev,clwp) 
  endif
  
  !stop
  !gdata%cloud_water = 0.
  !print '(4F8.2)', gdata%psfc, gdata%tskin, gdata%t2m, gdata%w10m
  !print '(31F7.2)', gdata%temp_lev
  !print '(31F7.4)', gdata%qv_lev
  !print '(31F7.4)', gdata%cloud_water
  !print '(31F7.2)', gdata%hgt_lev
end subroutine get_state_clw1d_rh

subroutine get_state_clw2d_rh(nvar,neof,nclw,nlev,x,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,clwprofs,gdata)

  use profile_def

  implicit none

  include 'parametersMERRA.inc'
  include 'constants_sjm.inc'

  integer, intent(in) :: nvar,nlev,neof,nclw
  real, intent(in) :: x(nvar), xbar(nvar), xstd(nvar)
  real, intent(in) :: prof_avg(1+2*MERRA_NLEV), prof_W(1+2*MERRA_NLEV), prof_std(1+2*MERRA_NLEV)
  real, intent(in) :: prof_eofs(1+2*MERRA_NLEV,1+2*MERRA_NLEV)
  type(profile), intent(inout) :: gdata
  real, intent(out) :: clwprofs(nclw, MERRA_NLEV)

  integer :: i, j, i0
  real :: Tavg, Qavg, dz, R, clwp, tpw, es, pv, rhov, rhod, rh_lev(MERRA_NLEV)

  i0=0

  do i=1,MERRA_NLEV
    if(gdata%psfc .ge. MERRA_PRESSLEV(i)) exit
    i0=i0+1
  end do

  gdata%tskin = x(neof+3)*xstd(neof+3)+xbar(neof+3)
  !print*, x(nvar), xstd(nvar), xbar(nvar), gdata%tskin
  
  if(gdata%anc_t2m .gt. 0.) then
    gdata%t2m = gdata%anc_t2m
    !gdata%t2m = prof_avg(1)
  else! if(prof_avg(1) .gt. 0.) then
    !gdata%t2m = prof_avg(1)
    print*, 'ancillary data missing for t2m'
    stop
  endif

  
  do i=1,nlev
    if(gdata%anc_temp_lev(i+MERRA_NLEV-nlev) .ge. 0.) then
      gdata%temp_lev(i) = gdata%anc_temp_lev(i+MERRA_NLEV-nlev)
      !gdata%temp_lev(i) = prof_avg(i+1+i0)
    else if(i .eq. 1) then
      gdata%temp_lev(i) = gdata%t2m
    else if(i .gt. 1) then
      gdata%temp_lev(i) = gdata%temp_lev(i-1)
    else
      print*, 'ancillary data missing for temp lev',i
      stop
    endif
    gdata%press_lev(i) = MERRA_presslev(i+MERRA_NLEV-nlev)
    if(gdata%anc_qv_lev(i+MERRA_NLEV-nlev) .ge. 0.) then
      gdata%qv_lev(i) = gdata%anc_qv_lev(i+MERRA_NLEV-nlev)
      !calculate MERRA2 RH
    else if(gdata%anc_qv_lev(i+MERRA_NLEV-nlev+1) .ge. 0) then
      gdata%anc_qv_lev(i+MERRA_NLEV-nlev) = gdata%anc_qv_lev(i+MERRA_NLEV-nlev+1)
    else 
      print*, 'ancillary data missing for Qv lev',i, gdata%anc_qv_lev(:), gdata%anc_temp_lev(:)
      stop
    endif
    call get_es(gdata%temp_lev(i),es)
    pv = (gdata%qv_lev(i)*100*gdata%press_lev(i)*Rv)/(Rd*(1+(gdata%qv_lev(i)*Rv)/Rd))
    !rhov = (100*es)/(Rv*gdata%temp_lev(i))
    !rhod = (100*MERRA_PRESSLEV(i)-100*es)/(Rd*gdata%temp_lev(i))
    rh_lev(i) = pv/(100*es)!gdata%qv_lev(i)/(rhov/rhod)
  end do
  !print*, nlev, i0
  !print*, gdata%temp_lev(1:nlev)
  !print '(10F8.4)', prof_std(:)
  !print '(10F8.4)', prof_std(2+i0:1+i0+nlev)
  !print '(10F8.4)', prof_std(2+i0+MERRA_NLEV:1+i0+MERRA_NLEV+nlev)
  !stop
  !if(i0+1+MERRA_NLEV+nlev .gt. 55) print*, i0, nlev, gdata%psfc
  !add eofs
  do i=1,neof
    !print*, i, x(i), xstd(i), xbar(i)
    gdata%t2m = gdata%t2m+prof_eofs(i,1)*(x(i)*xstd(i)+xbar(i))*prof_std(1)
    gdata%temp_lev(1:nlev) = gdata%temp_lev(1:nlev)+prof_eofs(i,2+i0:2+i0+nlev)*(x(i)*xstd(i)+xbar(i))*prof_std(2+i0:2+i0+nlev)
    rh_lev(1:nlev) = rh_lev(1:nlev)+prof_eofs(i,2+i0+MERRA_NLEV:1+i0+MERRA_NLEV+nlev)*(x(i)*xstd(i)+xbar(i))*prof_std(2+i0+MERRA_NLEV:1+i0+MERRA_NLEV+nlev)
  end do

  !to do: check saturated/negative RH and convert back to Qv
  do i=1,nlev
    if(rh_lev(i) .lt. 0.) rh_lev(i) = 0.
    if(rh_lev(i) .gt. 1.01) rh_lev(i) = 1.
    call get_es(gdata%temp_lev(i),es)
    rhov = (rh_lev(i)*es)/(Rv*gdata%temp_lev(i))
    rhod = (gdata%press_lev(i)-rh_lev(i)*es)/(Rd*gdata%temp_lev(i))
    gdata%qv_lev(i) = rhov/rhod
    !print*, gdata%temp_lev(i), rh_lev(i), gdata%qv_lev(i), rhov/rhod
  end do
  !do i=1,nlev
  !  print*, i, gdata%press_lev(i), gdata%temp_lev(i), prof_avg(i+2+i0), gdata%qv_lev(i), prof_avg(2+nlev+i+i0)
  !end do
  !stop
  !print '(40F6.3)', gdata%qv_lev(1:nlev)
  !print*, gdata%qv_lev
  !calculate height levels
  Tavg = 0.5*(gdata%t2m+gdata%temp_lev(1))
  Qavg = gdata%qv_lev(1)
  R = (1.-Qavg)*Rd+Qavg*Rv
  dz = R*Tavg/9.81*alog(gdata%psfc/gdata%press_lev(1))
  gdata%hgt_lev(1) = gdata%sfc_elev+0.001*dz
  !print*, 1, Tavg, Qavg, gdata%press_lev(1), gdata%hgt_lev(1)
  do i=2,nlev
    Tavg = 0.5*(gdata%temp_lev(i-1)+gdata%temp_lev(i))
    Qavg = 0.5*(gdata%qv_lev(i-1)+gdata%qv_lev(i-1))
    R = (1.-Qavg)*Rd+Qavg*Rv
    dz = R*Tavg/9.81*alog(gdata%press_lev(i-1)/gdata%press_lev(i))
    gdata%hgt_lev(i) = gdata%hgt_lev(i-1)+0.001*dz
    !print*, i, Tavg, Qavg, gdata%press_lev(i), gdata%hgt_lev(i)
  end do

  !add wind and cloud water
  gdata%w10m = x(neof+1)*xstd(neof+1)+xbar(neof+1)!gdata%iws
  if(gdata%w10m .lt. 0.) gdata%w10m = 0.

  !print*, gdata%w10m
  !relative azimuth
  gdata%relAz = x(neof+2)*xstd(neof+2)+xbar(neof+2)
  
  do i=1,nclw
    clwp = exp(x(neof+3+i)*xstd(neof+3+1)+xbar(neof+3+i))
    !print '(I5,F8.3)', i, clwp
    if(gdata%anc_clwp .gt. 0.) then
      do j=1,nlev
        if(not(isnan(gdata%anc_clw_lev(j+MERRA_NLEV-nlev)))) clwprofs(i,j) = 1e3*gdata%anc_clw_lev(j+MERRA_NLEV-nlev)/gdata%anc_clwp*clwp
        !print*, i, gdata%temp_lev(i), gdata%anc_clw_lev(i), gdata%anc_clw_lev(i)/gdata%anc_clwp*clwp
      end do
    else
      call get_clw_profile(gdata,nlev,clwp)
      clwprofs(i,:) = gdata%cloud_water
    endif
  end do
  
  !gdata%cloud_water = 0.
  !print '(4F8.2)', gdata%psfc, gdata%tskin, gdata%t2m, gdata%w10m
  !print '(31F7.2)', gdata%temp_lev
  !print '(31F7.4)', gdata%qv_lev
  !print '(31F7.4)', gdata%cloud_water
  !print '(31F7.2)', gdata%hgt_lev
end subroutine get_state_clw2d_rh

subroutine get_clw_profile(atm,nlev,clwp) 
  !This routine creates a cloud by lifting a surface parcel and condensing until the desired cloud LWP is reached.
  
  use profile_def
  implicit none
  include 'constants_sjm.inc'
  
  type(profile) :: atm
  real :: clwp, satQv, es, clwp_sum, clwp_add, cloudQv1,cloudQv2,Tavg, Pavg
  integer :: nlev,i, i0, LCL
  
  atm%cloud_water = 0.
  if(clwp .eq. 0.) return
  !Find LCL
  do i=1,nlev
    if(atm%hgt_lev(i) .ge. atm%sfc_elev .and. atm%temp_lev(i) .gt. 0.) then
      i0=i
      exit
    endif
  end do
  LCL=i0

!diagnostic
!  if(LCL .le. 0) then
!    do i=1,nlev
!      write(6, '("i: ", i10, "  sfc_lev: ", f12.4, "  hgt_lev: ", f12.4, "  temp_lev: ", f12.4)') &
!       i, atm%sfc_elev, atm%hgt_lev(i), atm%temp_lev(i)
!    end do
!    write(*, '("LCL: ", i10, "  atmheight: ", f12.4, "  atmtemp: ", f12.4)') LCL, atm%hgt_lev(LCL), atm%temp_lev(LCL)
!  endif
!end diagnostic

  do i=i0,nlev
    call get_es(atm%temp_lev(i), es)
    satQv = es/(es+(atm%press_lev(i)-es)*Rv/Rd)
    !print*, atm%Qv_lev(i0), satQv
    if(atm%Qv_lev(i0) .gt. satQv) then
      LCL = i
      exit
    endif
  end do

!diagnostic
!  if(LCL .le. 0) then
!    write(*, '("LCL: ", i10)') LCL 
!  endif
!end diagnostic
  
  !print*, LCL, atm%hgt_lev(LCL), atm%temp_lev(LCL)
  !print*, clwp
  !stop
  !Now create the cloud. The cloud water content at each level should be the excess of the surface mixing ratio minus the saturation value.
  !Convert mixing ratio to mass, and add to cloud water path. If CLWP exceeds prescribed value, then reduce top level so that they are equal, then exit
  cloudQv1 = 0.
  clwp_sum = 0.
  do i=LCL+1, nlev
    if(atm%temp_lev(i) .lt. 240.) exit !do not allow extreme supercooling or clouds in the stratosphere - 
    call get_es(atm%temp_lev(i), es)
    satQv = es/(es+(atm%press_lev(i)-es)*Rv/Rd)
    cloudQv2 = atm%Qv_lev(i0)-satQv
    if(cloudQv2 .lt. 0.) cloudQv2 = 0.
    !print*, cloudQv2
    
    Tavg = 0.5*(atm%temp_lev(i-1)+atm%temp_lev(i))
    Pavg = (atm%press_lev(i)-atm%press_lev(i-1)) / log(dble(atm%press_lev(i))/dble(atm%press_lev(i-1)))
    atm%cloud_water(i) = 1000.*(0.5*(cloudQv1+cloudQv2))*(100.*es/(Rv*Tavg)+100.*(Pavg-es)/(Rd*Tavg)) !g/m^3

    clwp_add = (atm%hgt_lev(i)-atm%hgt_lev(i-1))*atm%cloud_water(i)
    if(clwp_add+clwp_sum .le. clwp) then
      clwp_sum = clwp_sum+clwp_add
      !print*, i, atm%hgt_lev(i), atm%cloud_water(i), clwp_sum
    else
      atm%cloud_water(i)=atm%cloud_water(i)*(clwp-clwp_sum)/clwp_add
      !print*, i, atm%hgt_lev(i), atm%cloud_water(i), clwp_sum+(atm%hgt_lev(i)-atm%hgt_lev(i-1))*atm%cloud_water(i)
      exit
    endif
    cloudQv1=cloudQv2
  end do
  
  !convert excess cloud water to rain water
!   if(clwp .gt. 0.2) then
!     atm%cloud_water = atm%cloud_water*0.2/clwp
!     atm%rain_water(1:nlev-1) = (clwp-0.2)/(atm%hgt_lev(nlev-1)-atm%sfc_elev)
!   endif
  
  
end subroutine get_clw_profile

subroutine check_saturation(atm,nlev)
  !This routine  checks if the water vapor mixing ratio exceeds saturation at any level.
  !If so, reduce to saturation value
  
  use profile_def
  
  implicit none
  
  include 'constants_sjm.inc'
  
  type(profile) :: atm
  integer :: i, nlev
  real :: es, satQv
  
  do i=1,nlev
    if(atm%temp_lev(i) .gt. 0.) then
      call get_es(atm%temp_lev(i), es)
      satQv = es/(es+(atm%press_lev(i)-es)*Rv/Rd)
      !print*, es, satQv, atm%Qv_lev(i)

      if(atm%Qv_lev(i) .gt. satQv .and. es .gt. 0.) atm%Qv_lev(i) = satQv
      if(atm%Qv_lev(i) .lt. 0.) atm%Qv_lev(i) = 0.
    endif
  end do

  
end subroutine check_saturation


subroutine get_es(temp, es)
  !Just calculates the saturation vapor pressure for a given temperature (in K)
  !Based on formula in Pruppacher and Klett 1995
  real :: a0, a1, a2, a3, a4, a5, a6, tc 
  real,intent(in) ::   temp
  real,intent(out) ::  es

  data a0 / 6.11176750 /
  data a1 / 4.43986062e-1 /
  data a2 / 1.43053301e-2 /
  data a3 / 2.65027242e-4 /
  data a4 / 3.02246994e-6 /
  data a5 / 2.03886313e-8 /
  data a6 / 6.38780966e-11 /
      
  tc = temp-273.16
  !es = a0+tc*(a1+tc*(a2+tc*(a3+tc*(a4+tc*(a5+a6*tc)))))
  !if(es .lt. 0.) es = 0.
  !if(tc .lt. -50.) es = 0.
  es = 6.1070*exp(17.15*(tc)/(tc+273.16-38.25))
  
end subroutine get_es

subroutine calc_tpw(atm,nlev,tpw)
  !Integrate water vapor to total precipitable water (mm)
  use profile_def
  
  implicit none
  
  include 'constants_sjm.inc'
  
  integer :: i0,i
  integer, intent(in) :: nlev
  type(profile), intent(in) :: atm
  real, intent(out) :: tpw
  
  real :: rhov_bot, rhov_top, Tv, Qv_sfc, Tavg, Pavg
  
  real :: delp
  
  !calculate rhov in sfc to first above-surface layer
  Tavg = 0.5*(atm%t2m+atm%temp_lev(1))
!begin WSO 8/27/21; prevent psfc <= press_lev(1)
  delp = 0.1
  if(atm%psfc-atm%press_lev(1) > delp) then
     Pavg = (atm%psfc-atm%press_lev(1)) / log(dble(atm%psfc)/dble(atm%press_lev(1)))
  else
     Pavg = delp / log(dble(atm%press_lev(1) + delp)/dble(atm%press_lev(1)))
!     write(*, '("press(1): ", e15.6, "  delp: ", e15.6, "  (p+delp)/p: ", e15.6, "  log((p+delp)/p): ", e15.6, "  Pavg: ", e15.6)') &
!      atm%press_lev(1), delp, dble(atm%press_lev(1) + delp)/dble(atm%press_lev(1)), log(dble(atm%press_lev(1) + delp)/dble(atm%press_lev(1))), Pavg
  endif
!end WSO 8/27/21
  
  rhov_top = 100.*atm%qv_lev(1)*Pavg/(Rd*Tavg*(1.+atm%qv_lev(1)*(Rv/Rd-1.)))
  !print*, rhov_top
  tpw=1000.*(rhov_top)*(atm%hgt_lev(1)-atm%sfc_elev)
  !calculate rhov in subsequent layers
  rhov_bot=rhov_top
  
  do i=2,nlev
    rhov_top = 100.*atm%qv_lev(i)*atm%press_lev(i)/(Rd*atm%temp_lev(i)*(1.+atm%qv_lev(i)*(Rv/Rd-1.)))
    !print*, i,atm%hgt_lev(i), atm%qv_lev(i), atm%temp_lev(i), rhov_top, rhov_bot
    tpw=tpw+1000.*0.5*(rhov_top+rhov_bot)*(atm%hgt_lev(i)-atm%hgt_lev(i-1))
    rhov_bot=rhov_top
  end do
  !print*, tpw

end subroutine calc_tpw

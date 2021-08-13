

subroutine rtm_emission_csu(atm,nlyr,ch_mask,Tb_out)

  use profile_def
  use RSS_RTM

  implicit none
  
  include 'parametersMERRA.inc'
  include 'parametersGMI.inc'
  include 'constants_sjm.inc'
  
  logical, intent(in) :: ch_mask(13)
  type(profile), intent(inout) :: atm
  integer, intent(in) :: nlyr
  real, intent(out) :: Tb_out(13)

  logical, parameter :: diffuse   = .TRUE.
  integer, parameter :: verbose   = 0
  real,    parameter :: salinity  = 35.0
  real,    parameter :: mis_val = -999.9
  
  integer :: ichan, ilyr
  integer :: ipol(13)
  data ipol /1,2,1,2,1,1,2,1,2,1,2,1,1/
  
  real    :: emis
  real    :: ebar
  real    :: refl
  real    :: view_angle
  real    :: azim
  real    :: e0(2), ewind(2), eharm(2,4), land_emis(13)
  real    :: omega(2), trans, optdepth
  
  real    :: tavg, pavg, qavg, clw_avg, atm_ext, kext_clw, salb_clw, asym_clw
  real    :: hgt_lev(0:MERRA_NLEV), temp_lev(0:MERRA_NLEV), press_lev(0:MERRA_NLEV), kext(MERRA_NLEV)
  real    :: tb, tbdown
  real    :: tb_down(13)
  real    :: cfreq(13)
  data cfreq /10.65,10.65,18.70,18.70,23.80,36.64,36.64,89.00,89.00,166.00,166.00,186.31,190.31/
  
  
  press_lev(0) = atm%psfc
  press_lev(1:nlyr) = atm%press_lev(1:nlyr)
  hgt_lev(0) = atm%sfc_elev
  temp_lev(0) = atm%t2m
  hgt_lev(1:nlyr) = atm%hgt_lev(1:nlyr)
  temp_lev(1:nlyr) = atm%temp_lev(1:nlyr)
  !print*, atm%sfc_type
  do ichan = 1, 13
    !print*, ichan, ch_mask(ichan)
    if(ch_mask(ichan)) then
      !print*, ichan, GMI_nf(ichan)
      if(ichan .le. 9) then
        view_angle = atm%eia(1)
      else
        view_angle = atm%eia(2)
      endif
      optdepth = 0.0
      !print*, view_angle
      do ilyr = 1, nlyr   
        if(press_lev(ilyr) .eq. press_lev(ilyr-1)) then
          pavg = press_lev(ilyr)
        else
          pavg = (press_lev(ilyr) - press_lev(ilyr-1)) / log(press_lev(ilyr)/press_lev(ilyr-1))
        endif
        tavg = (temp_lev(ilyr) + temp_lev(ilyr-1)) / 2.0
        if(ilyr .gt. 1) then
          qavg = (atm%qv_lev(ilyr)+atm%qv_lev(ilyr-1)) / 2.0
          clw_avg = 0.5*(atm%cloud_water(ilyr)+atm%cloud_water(ilyr-1))
        else
          qavg = atm%qv_lev(ilyr)
          clw_avg = atm%cloud_water(ilyr)
        endif
        !print*, ilyr, atm%cloud_water(ilyr)
        !call monortm_lut(ifreq(ichan), pavg, tavg, mix_ratio(ilyr), atm_ext)
        !print*, ichan, ilyr, tavg, pavg, qavg
        call intplte_gas(ichan,tavg,pavg,qavg,atm_ext)
        !call mie_clw(cfreq(ichan), tavg, cloud_water(ilyr), kext_clw, salb_clw, asym_clw)
        call intplte_clw(GMI_nf(ichan),Tavg,kext_clw)
        kext_clw = clw_avg*kext_clw
        kext(ilyr) = atm_ext + kext_clw
        optdepth = optdepth + kext(ilyr) * (hgt_lev(ilyr) - hgt_lev(ilyr-1))
        !print*, pavg, tavg, qavg, atm_ext, clw_avg, kext_clw, hgt_lev(ilyr), optdepth
      end do 
      
      trans = exp(-optdepth)
      !print*, trans
      !print*, atm%relAz,atm%w10m
      if((atm%sfc_type .eq. 0)) then! .or. (minval(atm%emis) .le. 0)) then !ocean simulation
        !print*, nlyr
        !print*, temp_lev
        !print*, kext
        call radtranCSU(nlyr, atm%sfc_type, view_angle, atm%tskin, hgt_lev, temp_lev, kext, 0.5, 0.5, tb, tbdown)
        !write(6,*) cfreq(ichan)!, tbdown, tb
        !print*, atm%tskin, atm%w10m, view_angle, atm%relAz, tbdown
        call find_surface_tb(freq=cfreq(ichan), surtep=atm%tskin, ssws=atm%w10m, tht=view_angle, &
                             phir=atm%relAz, sal=salinity, e0=e0, ewind=ewind, eharm=eharm, &
                             tran=trans, tbdw=tbdown, omega=omega)
        emis = e0(ipol(ichan)) + ewind(ipol(ichan)) + eharm(1,ipol(ichan))*cos(atm%relAz*dtor) + eharm(2,ipol(ichan))*cos(2.*atm%relAz*dtor)
        atm%emis(ichan) = emis
        !call intplte_emis(GMI_nf(ichan),ipol(ichan),atm%tskin,atm%w10m,atm%relAz,view_angle,emis,ebar)
        if (diffuse) then
          refl = (1.0 - emis) * (1.0 + omega(ipol(ichan)))
        else
          refl = 1.0 - emis
        endif
        if (refl .lt. (1.0 - emis)) refl = 1.0 - emis ! Reflectivity cannot be less than 1.0 - emis
        if (refl .gt. 1.0) refl = 1.0                 ! Reflectivity cannot be greater than 1.0
      else !land simulation
        emis = atm%emis(ichan)
        refl = 1.0-emis
      endif
      !print*, ichan, nlyr, atm%sfc_type, view_angle, atm%tskin, hgt_lev(0:nlyr), temp_lev(0:nlyr), kext(1:nlyr), emis, refl
      
      call radtranCSU(nlyr, atm%sfc_type, view_angle, atm%tskin, hgt_lev(0:nlyr), temp_lev(0:nlyr), kext(1:nlyr), emis, refl, tb, tbdown)
      !print*, ichan, tb, tbdown
      !stop
      !if ((tb .gt. 50.0) .and. (tb .lt. 350.0)) then
      tb_out(ichan)  = tb
      tb_down(ichan) = tbdown
      !else
      !  write(6,*) ' Tb is outside physical range'
      !  print*, tb, tbdown, emis
      !  stop
      !endif    
    endif !ch_mask
  end do  ! end loop over nchan
  
  !stop
end subroutine rtm_emission_csu


  

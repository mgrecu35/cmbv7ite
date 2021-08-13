subroutine calcPIA(atm,nlev,PIAKu,PIAKa)
  !Calculate the PIA from gases and clouds at Ku and Ka band. This is to correct sigma_zero in analyses.
  
  use profile_def
  !use Type_Kinds
  !use RTV_Define
  !use Emission_Module
  !use CRTM_Parameters
  implicit none
  
  include 'constants_sjm.inc'
  include 'parametersGMI.inc'
  include 'parametersMERRA.inc'
  
  !I/O parameters
  type(profile), intent(inout) :: atm
  real, intent(out) :: PIAKu, PIAKa
  
  integer :: k, nlev, binsfc
  real :: kextatm, kextclw, dz, Pavg, Tavg, qv, pv
  
  !nlev=0
  
  PIAKu = 0.
  PIAKa = 0.
  

  
  do k=1,nlev
    
    kextatm = 0.
    kextclw = 0.
    
    !calculate gas extinction
    !average temperature and pressure in this layer for H2O and O2 extinction
    if(k .eq. 1) then
      Tavg = (atm%temp_lev(k)+atm%t2m)/2.
      if(atm%psfc-atm%press_lev(k) .gt. 1.) then
        Pavg = (atm%psfc-atm%press_lev(k)) / log(atm%psfc/atm%press_lev(k))
      else
        Pavg = 0.5*(atm%psfc+atm%press_lev(k))
      endif
      qv = atm%qv_lev(k)
      dz = (atm%hgt_lev(k)-atm%sfc_elev)
    else
      Tavg = (atm%temp_lev(k)+atm%temp_lev(k-1))/2.
      Pavg = (atm%press_lev(k) - atm%press_lev(k-1))/log(atm%press_lev(k)/atm%press_lev(k-1))
      qv = (atm%qv_lev(k)+atm%qv_lev(k-1))/2.
      dz = (atm%hgt_lev(k)-atm%hgt_lev(k-1))
    endif
    !pv = Qv*(Pavg*Rv/Rd)/(1.-Qv*(1.+Rv/Rd))
    !print*, k, Tavg, Pavg,Qv
    call intplte_gas(14,Tavg,Pavg,Qv,kextatm)
    call intplte_clw(11,Tavg,kextclw)
    kextclw = atm%cloud_water(k)*kextclw
    PIAKu = PIAKu+(kextclw+kextatm)*dz*2.
    
    call intplte_gas(15,Tavg,Pavg,Qv,kextatm)
    call intplte_clw(12,Tavg,kextclw)
    kextclw = atm%cloud_water(k)*kextclw
    PIAKa = PIAKa+(kextclw+kextatm)*dz*2.
  end do
  !convert to dB
  PIAKu = 10.*alog10(exp(PIAKu))
  PIAKa = 10.*alog10(exp(PIAKa))
  !print*, PIAKu, PIAKa
  
  
end subroutine
module profile_def !define required profile properties
  type profile
    real*4 :: sclat, sclon, scalt, scyaw
    real*4 :: lat, lon
    integer*4 :: scan, pos
    real*4 :: gmi_tb(13), emis(13), eia(2), sga
    real*4 :: tskin, t2m, psfc, w10m, relAz, sfc_elev
    !ancillary data
    real*4 :: anc_tskin, anc_t2m, anc_psfc, anc_u10m, anc_v10m, anc_sfc_elev, anc_clwp
    real*4 :: press_lev(27), anc_temp_lev(27), anc_hgt_lev(27), anc_qv_lev(27), anc_clw_lev(27)
    !retrieved atmosphere
    real*4 :: temp_lev(27), hgt_lev(27), qv_lev(27), cloud_water(27), rain_water(27), cloud_ice(27)
    integer*4 :: sfc_type, ice_type, tskin_index
    integer*4 :: year, month, dayOfMonth
    real*4 :: secondOfDay
    !observed data
    real*4 :: isst, iws, iwd, tb89std
  end type profile
  
  type ret_data
    real*4 :: sim_tb0(13), sim_tb(13), emis(13), semis(13), aemis(13)
    real*4 :: w10m_init, w10m_ret, sw10m_init, sw10m_ret, aw10m
    real*4 :: relAz_init, relAz_ret, srelAz_init, srelAz_ret, arelAz
    real*4 :: tpw_init, tpw_ret, stpw_init, stpw_ret
    real*4 :: tskin_init, tskin_ret, stskin_init, stskin_ret, atskin_ret
    real*4 :: clwp_init, clwp_ret, sclwp_init, sclwp_ret, aclwp
    real*4 :: wc_ret(5), swc(5), awc(5)
    real*4 :: rlwp_ret, iwp_ret, prate_sfc
    real*4 :: chi_squared, obs_err, x_err, chisq_ocean, chisq_land, chisq_ice
    real*4 :: PIAKu, PIAKa
    real*4 :: specularity
    real*4 :: err_eof(10), aerr_eof(10), serr_eof(10)
    real*4 :: press_lev(27), temp_lev(27), hgt_lev(27), qv_lev(27), cloud_water(27), t2m_ret
    integer :: niter, sfc_type
  end type ret_data
end module profile_def

! module limit_def
!   type limit
!     real :: theta(2), chi(2), phi(2), caph(2), nu(2), mu(2)
!   end type limit
! end module limit_def

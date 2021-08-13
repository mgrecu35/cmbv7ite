  real :: Rd, Rv !gas constants
  data Rv /461.5/
  data Rd /286.9/
  
  real :: radius_earth !radius of the earth [km]
  data radius_earth /6378./
  
  real :: dtor !degree to radians conversion
  data dtor /0.0174533/
  
  real :: PR_MISSING, PR_CLUTTER
  data PR_MISSING /-99.99/
  data PR_CLUTTER /-88.88/
  
  real :: default_melt_density_strat, default_melt_density_conv
  data default_melt_density_strat /200./
  data default_melt_density_conv /800./
  
  real :: Dm_Z_a_strat, Dm_Z_b_strat, Dm_Z_a_conv, Dm_Z_b_conv
  parameter(Dm_Z_a_strat = 0.5973)
  parameter(Dm_Z_b_strat = 0.1073)
  parameter(Dm_Z_a_conv = 0.4778)
  parameter(Dm_Z_b_conv = 0.1210)

  !Shape parameter power law: sigma_m = aDm^b
  real :: sigmaM_Dm_a
  parameter(sigmaM_Dm_a = 0.29)
  real :: sigmaM_Dm_b
  parameter(sigmaM_Dm_b = 1.5)
  
  !maximum Ku-band PIA allowed (to prevent profiles from "blowing up")
  real :: PIA_limit
  parameter(PIA_limit=50.0)
  
  !maximum dBZ allowed (to prevent profiles from "blowing up")
  real :: dBZ_limit
  parameter(dBZ_limit=65.0)
  
  !maximum 85 GHz kext allowed (to prevent RTM error)
  real :: kext85_limit
  parameter(kext85_limit=100.0)
  
  !minimum reflectivity to perform retrieval
  real :: precip_thres_Ku
  parameter(precip_thres_Ku=18.0)
  
  !Beta for NUBF correction
  real :: beta_nubf
  parameter(beta_nubf=0.77)
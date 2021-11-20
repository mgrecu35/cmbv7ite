subroutine radarRetSub4(nmu2,  nmfreq2,   icL, tbRgrid,               &
      dprrain,ichunk,orbNumb,ialg,idir,nscans_c)
!  SFM  end    12/13/2013
  use local_RD_var
  use globalData
  use f90DataTypes
  use f90Types
  use cldclass
  use ran_mod
  use geophysEns
  use nbinMod
  !use tables2
  use weight
  Use BMCVparameters
  use emissMod
!begin  MG 10/29/15 add gEnv module
  use gEnv
!end    MG 10/29/15
!begin  WSO 9/14/13 incorporate missing flags
  use missingMod
!end    WSO 9/14/13
!begin  WSO 6/5/18 add limits to output variables
  use outputminmax
!end    WSO 6/5/18
  use LUT_def !SJM 7/9/2015
  use oetbs
  implicit none
!  type (geoDataType) :: geoData
!  type (gridEnvDataType) :: gridEnvData
!  type (dPRDataType)     :: dPRData
!  type (dPRRetType)      :: dPRRet
!  type (gmi2GridType)    :: gmi2Grid
!  type (cgMIDataType)     :: gmiData
  integer :: nmu2, nmfreq2
  integer*4 :: ichunk
!  integer :: st_2adpr              ! file open status for 2adpr file
  real :: tbRgrid(14,49,9300), dprrain(49,300)
 ! integer :: tbRgridIn(9,49,9300)

  real                   :: meansfcRain,stddevSfcRain, tbout(14)
  !type (radarRetType)    :: radarRet
  !type (radarDataType)   :: radarData
  type (retParamType)    :: retParam
  real :: xin363(363)
 
  integer :: ii, jj, iGMI, jGMI
  integer :: di(8), dj(8)  
  integer :: i, j, ig, jg, ntpw, nmemb1, itop, irand
  real    :: pia13m, rms1, rms2,  unSortedRR(200), corrcoef, sfcRain2, tpw_ij
  integer :: iy(200), kflag, it
  real    :: sysdNl, pia13mean
  integer :: iLandSea(5,5), i1, j1, igetlandsea, ic2, nobs, iit
  real :: a0(2,8), emiss(2), tb19v, tb19h, tb37v, tb37h, tb22v, tb85v, tb85h
  real :: meanTb85
  real :: stdTb85,kgain,kgain37, ymean(3)
!...Ensemble parameters
  !real, allocatable ::  Yens(:,:), Xens(:,:), Yobs(:), Xup(:)
  integer :: ibatch
  real    :: stddev, srtpiaf
  real    :: FWHMx, FWHMy, tbconv(2), tbconvEns(2,100)
  integer :: dnx,dny, ik
  !real, allocatable :: ndn(:), ndnp(:,:), xscalev(:), logdNwf(:), randemiss(:), dwind(:)
  !real, allocatable  :: rhPCij(:,:), cldwPCij(:,:)
  real :: cldw(nlayer), rh(nlayer), pia13s
  integer :: nx,ny, icount, imin
  real ::  xm1,xs,rmsmin, prob, probtot, rmstot
  real :: piaR(100), fPIA, z13m
  integer :: ntbpix, ntbpix2
  real :: emtbm(9)
  real :: zminsc
  real :: realOut(49)
  !real :: w10(49,300), w10_out_NS(49,300), w10_out_MS(49,300), w10_min, w10_max, emis, relAz
  !real :: w10_rms_NS(49,300), emis_rms_NS(49,300,13), w10_rms_MS(49,300), emis_rms_MS(49,300,13)
  !real :: dZms(49,300) !! MS addition Feb 10, 2017
  !integer :: msFlag(49, 300) !!WSO addition Feb 11, 2017
!begin  WSO 2/8/17 new variables
  integer :: multiscatcalc_NS(49, 300), multiscatcalc_MS(49, 300)
  integer :: algotype_NS(49, 300), algotype_MS(49, 300)
  integer :: profclass_NS(49, 300), profclass_MS(49, 300)
  real :: subfootvariability_NS(49, 300), subfootvariability_MS(49, 300)
  real :: multiscatsurface_NS(49, 300), multiscatsurface_MS(49, 300)
  real :: skintempsigma_NS(49, 300), skintempsigma_MS(49, 300)
  real :: columnvaporsigma_NS(49, 300), columnvaporsigma_MS(49, 300)
  real :: columncloudliqsigma_NS(49, 300), columncloudliqsigma_MS(49, 300)
  real :: errorofdatafit_NS(49, 300), errorofdatafit_MS(49, 300)
  real :: initnw_NS(nbin, 49, 300), initnw_MS(nbin, 49, 300) 
  real :: princomp_NS(5, 49, 300), princomp_MS(5, 49, 300)
  real :: surfprecipbiasratio_NS(49, 300), surfprecipbiasratio_MS(49, 300)
!end    WSO 2/8/17 
  integer :: l, ipias(2)
  character*3 :: ifdpr, iftest
  character*90 :: outfile
  integer :: ink
  DOUBLE PRECISION input(6)
  DOUBLE PRECISION output(2)
  real :: wfract(5,5), wfractm, wfractsd
  real                    emissv(n_chan)
  real                    emissh(n_chan)
  real                    emissv_std(n_chan)
  real                    emissh_std(n_chan)
  integer :: stype!SJM 7/9/2015
  real  :: vLand(18,18), vOcean(10,10)
  real  :: pMLand(18), pMOcean(10)
  real  :: mTbLand(9), mTbOcean(9)
  real  :: stTbLand(9), stTbOcean(9)
  double precision  :: xin(18), xpred(18), yout(9)
  !real, allocatable :: emissoutL(:,:,:), emis_out_NS(:,:,:), emis_out_MS(:,:,:) !sjm 8/10/15
!begin  WSO 8/19/13 change Nw variable name (not dN) and add mu
  real :: cldwprof(88), cldiprof(88), log10NwMean(88), mu_mean_prof(88)
  integer *2 :: env_nodes(10, 49)
  real :: env_levs(10), ray_angle, pi
!end    WSO 8/19/13
  real :: lFract(49,300), sprobs, probs(100), rmsS(100)
  real :: covar, xf
  integer  :: orbNumb
  !begin SJM 7/25/14
  real :: s0Ku, s0Ka, s0stdKu, s0stdKa, s0corr, ds0Ku, ds0Ka
  !real :: sigmaZeroVarKu(49,300), sigmaZeroVarKa(49,300), sigmaZeroCov(49,300)
  !end SJM 7/25/2014
!begin WSO 8/8/13
  real :: gatelength
  real :: depthBB, depthML, depth
  real :: mu_mean(49, 300)
  real :: mu_meanMS(49, 300)
!  real :: scLatPR(49,300),scLonPR(49,300),wfmap(49,300), fpmap(49,300,15), fpmapN(49,300,15)
  !real :: S1eiaPR(49,300), S2eiaPR(49,300)
  real :: mlwc_frac(10, 49, 300)
  real :: mrate_frac(10, 49, 300)
  real :: mlwc_fracMS(10, 49, 300)
  real :: mrate_fracMS(10, 49, 300)
  real :: sfcRainLiqFrac(49, 300)
  real :: sfcRainLiqFracMS(49, 300)
  real :: tbMax1(15), tbMin1(15)

!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
!  SFM  begin  06/22/2014
  real :: wfractPix, windPert(100), windPertU(100), windPertV(100), qvPert(100), dnqv
!  SFM  end    06/22/2014
!  SFM  end    07/29/2014
!end   WSO 8/8/13
  

  integer :: actOb(49,300), iactOb
  integer :: jk, nf
  integer :: dig               ! SFM  04/16/2014  for M.Grecu
  real   :: cl(9,25), xin25(25),dtb(9)
  real   :: ebar, minl

  !integer*4 :: istart, iend
  integer :: iconv, ialg, icL
  real :: nubfc, stdpia35

  integer,parameter :: nscans=300, npixs=49, nlev=88, nchans=13
  integer :: nfreq, idir, nscans_c
  integer :: pType(nscans,npixs)
  real :: sfcTemp(nscans,npixs), cldw3d(nscans,npixs,nlev)
  integer :: clutFree(nscans,npixs)
  real :: pRate(nscans,npixs,nlev), swc3d(nscans,npixs,nlev), tbobsT(nscans,npixs,nchans)
  real :: z13(nscans,npixs,nlev),emiss2d(nscans,npixs,nchans)
  real :: nw3d(nscans,npixs,nlev), press3d(nscans,npixs,nlev), &
       airTemp3d(nscans,npixs,nlev),qv3d(nscans,npixs,nlev)
  integer :: binNodes(nscans,npixs,5)
  integer :: envNode(nscans,npixs,10)
  real    :: pRateOut(nscans,npixs,nlev), swcOut(nscans,npixs,nlev), nwOut(nscans,npixs,nlev)
  integer :: sfcBin(nscans,npixs)
  real    :: tbsim(nscans,npixs,nchans)

end subroutine radarRetSub4

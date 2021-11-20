!  SFM  04/06/2013  Code changes from M.Grecui
!
subroutine radarRetSubEmpty(geoData, dPRData, dPRRet, gridEnvData, gmi2Grid,        &
      gmiData, nmu2, nMemb, nmfreq2, ngates, sysdN, ic, tbRgrid,               &
      dprrain,ichunk,st_2adpr,orbNumb,istart,iend)
!  SFM  end    12/13/2013
  use f90DataTypes
  use f90Types
  use cldclass
  use ran_mod
  use geophysEns
  use nbinMod
  use tables2
  use weight
  Use BMCVparameters
  use emissMod
!begin  MG 10/29/15 add gEnv module
  use gEnv
!end    MG 10/29/15
!begin  WSO 9/14/13 incorporate missing flags
  use missingMod
!end    WSO 9/14/13

  implicit none
  type (geoDataType) :: geoData
  type (gridEnvDataType) :: gridEnvData
  type (dPRDataType)     :: dPRData
  type (dPRRetType)      :: dPRRet
  type (gmi2GridType)    :: gmi2Grid
  type (cgMIDataType)     :: gmiData
  integer :: nmu2, nMemb, nmfreq2, ngates
  integer*4 :: ichunk
  integer :: st_2adpr              ! file open status for 2adpr file
  real :: tbRgrid(14,49,9300), dprrain(49,300)
 ! integer :: tbRgridIn(9,49,9300)

  real                   :: meansfcRain,stddevSfcRain, tbout(14)
  type (radarRetType)    :: radarRet
  type (radarDataType)   :: radarData
  type (stormStructType) :: stormStruct
  type (retParamType)    :: retParam
  real :: xin363(363)
  real :: sfcRain(49,300),sfcRainStd(49,300)
  real :: rRate3D(nbin,49,300),  rRate3Dstd(nbin,49,300)
  real :: pwc3D(nbin,49,300),  pwc3Dstd(nbin,49,300)
  real :: zcKu3D(nbin,49,300), d03D(nbin,49,300), piaOut(49,300)

  real :: sfcRainMS(49,300),sfcRainStdMS(49,300),pia35m(49,300)
  real :: rRate3DMS(nbin,49,300),  rRate3DstdMS(nbin,49,300)
  real :: pwc3DMS(nbin,49,300),  pwc3DstdMS(nbin,49,300)
  real :: zcKu3DMS(nbin,49,300), zcKa3DMS(nbin,49,300), d03DMS(nbin,49,300), &
          piaOutKuMS(49,300), piaOutKaMS(49,300)

  integer :: ic, ii, jj, iGMI, jGMI
  integer :: di(8), dj(8)  
  integer :: i, j, k, ig, jg, ntpw, nmemb1, itop, irand
  real    :: pia13m, rms1, rms2, sysdN, unSortedRR(200), corrcoef, sfcRain2
  integer :: iy(200), kflag, it
  real    :: sysdNl, pia13mean
  integer :: iLandSea(5,5), i1, j1, igetlandsea, ic2, nobs, iit
  real :: a0(2,8), emiss(2), tb19v, tb19h, tb37v, tb37h, tb22v, tb85v, tb85h
  real :: meanTb85
  real :: stdTb85,kgain,kgain37, ymean(3)
!...Ensemble parameters
  real, allocatable ::  Yens(:,:), Xens(:,:), Yobs(:), Xup(:)
  integer :: ibatch
  real    :: stddev, srtpiaf
  real    :: FWHMx, FWHMy, tbconv(2), tbconvEns(2,100)
  integer :: dnx,dny, ik
  integer :: ipol(15), ifreq(15), iobs(15), ifreq1
  real, allocatable :: ndn(:), ndnp(:,:), xscalev(:), logdNwf(:), randemiss(:), dwind(:)
  real, allocatable  :: rhPCij(:,:), cldwPCij(:,:)
  real :: cldw(nlayer), rh(nlayer), pia13s
  integer :: nx,ny, icount, imin
  real ::  xm1,xs,rmsmin, prob, probtot, rmstot
  real :: piaR(100), fPIA, z13m
  integer :: ntbpix, ntbpix2
  real :: emtbm(9)
  real :: zminsc
  real :: realOut(49)
  real :: w10(49,300), w10_out_NS(49,300), w10_out_MS(49,300), w10_min, w10_max
  integer :: l, ipias(2)
  character*3 :: ifdpr, iftest
  character*90 :: outfile
  integer :: ink
  integer :: ifreqG(15), ipolG(15) 
  DOUBLE PRECISION input(6)
  DOUBLE PRECISION output(2)
  real :: wfract(5,5), wfractm, wfractsd
  real                    emissv(n_chan)
  real                    emissh(n_chan)
  real                    emissv_std(n_chan)
  real                    emissh_std(n_chan)
  real  :: vLand(18,18), vOcean(10,10)
  real  :: pMLand(18), pMOcean(10)
  real  :: mTbLand(9), mTbOcean(9)
  real  :: stTbLand(9), stTbOcean(9)
  double precision  :: xin(18), xpred(18), yout(9)
  real, allocatable :: emissout(:,:,:)
!begin  WSO 8/19/13 change Nw variable name (not dN) and add mu
  real :: cldwprof(88), cldiprof(88), log10NwMean(88), mu_mean_prof(88)
  integer *2 :: env_nodes(10, 49)
  real :: env_levs(10), ray_angle, pi
!end    WSO 8/19/13
  real :: lFract(49,300), sprobs, probs(100), rmsS(100)
  real :: covar
  integer  :: orbNumb
  !begin SJM 7/25/14
  real :: s0Ku, s0Ka, s0stdKu, s0stdKa, s0corr, ds0Ku, ds0Ka
  real :: sigmaZeroVarKu(49,300), sigmaZeroVarKa(49,300), sigmaZeroCov(49,300)
  !end SJM 7/25/2014
!begin WSO 8/8/13
  real :: gatelength
  real :: depthBB, depthML, depth
  real :: mu_mean(49, 300)
  real :: mu_meanMS(49, 300)
  real :: scLatPR(49,300),scLonPR(49,300),wfmap(49,300), fpmap(49,300,15), fpmapN(49,300,15)
  real :: mlwc_frac(10, 49, 300)
  real :: mrate_frac(10, 49, 300)
  real :: mlwc_fracMS(10, 49, 300)
  real :: mrate_fracMS(10, 49, 300)
  real :: sfcRainLiqFrac(49, 300)
  real :: sfcRainLiqFracMS(49, 300)
  real :: tbMax1(15), tbMin1(15)

!  SFM  begin  07/29/2014; for M.Grecu  eliminate NANs
!  SFM  begin  06/22/2014
  real :: wfractPix, windPert(100), qvPert(100), dnqv
!  SFM  end    06/22/2014
!  SFM  end    07/29/2014
!end   WSO 8/8/13
  real :: covTb(49,300,15,15), tbMax(49,300,15), tbMin(49,300,15), &
       tbMean(49,300,15)
  real :: invCovTb(49,300,15,15)
  real :: tbout2D(49,300,15), tb(49,300,15), tbNoOcean(49,300,15), &
       tbout2DNoOcean(49,300,15), tbObs(49,300,15), fobj
  real :: dfdtb(49,300,15), rerr(15), tb0(49,300,15), fem(15) , &
       tb0MS(49,300,15), tbNoOceanMS(49,300,15), tbout2DNoOceanMS(49,300,15),&
       tbout2DMS(49,300,15)

  integer :: actOb(49,300), iactOb
  integer :: jk, nf
  integer :: dig               ! SFM  04/16/2014  for M.Grecu
  real   :: cl(9,25), xin25(25),dtb(9)
  real   :: ebar, minl
  real, allocatable :: geoloc(:), hFreqTbs(:,:), PRgeoloc(:), hFreqPRg(:,:,:)
  integer*4 :: istart, iend
  integer :: iconv
  real :: nubfc, stdpia35
!begin  WSO 8/30/13 prescribe levels for environmental parameters

end subroutine radarRetSubEmpty


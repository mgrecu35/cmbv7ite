//  SFM 04/06/2013  Code module added in merge from M.Grecu's code
//  SFM 05/06/2013  Modifications from LW to facilitate using job names
//  SFM 06/27/2013  Parameter name changes from W.Olson; reduce unused code
//  SFM 07/19/2013  Large volume of code added for M.Grecu
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf.h>
#include <mfhdf.h>
#include "TKheaders.h"
#include "TK_2BCMB_hdf5.h"
#ifdef GFOR 
extern int __nbinmod_MOD_imemb;
#define nbins __nbinmod_MOD_imemb
//begin  WSO 9/15/13 
extern float __missingmod_MOD_missing_r4;
#define missing_r4c __missingmod_MOD_missing_r4
extern short __missingmod_MOD_missing_i2;
#define missing_i2c __missingmod_MOD_missing_i2
extern long __missingmod_MOD_missing_i4;
#define missing_i4c __missingmod_MOD_missing_i4
extern int __nbinmod_MOD_ntransition;
#define ntransitions __nbinmod_MOD_ntransition
//end    WSO 9/15/13
#endif

#ifdef IFORT 
extern int nbinmod_mp_nbin_;
//begin  WSO 8/8/13
extern int nbinmod_mp_ntransition_;
//end    WSO 8/8/13
#define nbins nbinmod_mp_nbin_
//begin  WSO 8/8/13
#define ntransitions nbinmod_mp_ntransition_
//end    WSO 8/8/13
//begin  WSO 9/15/13
extern float missingmod_mp_missing_r4_;
#define missing_r4c missingmod_mp_missing_r4_
extern short missingmod_mp_missing_i2_;
#define missing_i2c missingmod_mp_missing_i2_
extern long  missingmod_mp_missing_i4_;
#define missing_i4c missingmod_mp_missing_i4_
//end    WSO 9/15/13
#endif

//begin WSO 04/07/2013
//Note that there were many structure/variable name changes in this
//version to be compatible with TKIO 3.50.8
//All S1 and S2 were changed to NS and MS, respectively
//The variable ending "Out" was removed because a separate Input structure
//was created
//end WSO 04/07/2013

extern TKINFO dprtkfile;
TKINFO ctkfile;
TKINFO ctkfileIn;

L2BCMB_SWATHS swath;
L2ADPR_SWATHS dprswath;
L2ADPRX_SWATHS dprxswath;
L2BCMB_SWATHS swath1;
L2AKu_NS      L2AKuData;
L2AKuX_FS      L2AKuDataX;





void setlatlons1_fs_300_(int *isc,float *lat, float *lon, float *sfcPrecip, 
                  float *sfcPrecipStd, float *piaOut)
{
  int i;
  extern L2BCMBX_SWATHS swathx300[300];

  for(i=0;i<49;i++)
    {
      swathx300[*isc].NS.Latitude[i]=lat[i];
      swathx300[*isc].NS.Longitude[i]=lon[i];
      if(swathx300[*isc].NS.Longitude[i]>180)
	swathx300[*isc].NS.Longitude[i]-=360;
      swathx300[*isc].NS.surfPrecipTotRate[i]=sfcPrecip[i];
      swathx300[*isc].NS.surfPrecipTotRateSigma[i]=sfcPrecipStd[i];
      swathx300[*isc].NS.pia[i]=piaOut[i];
    }
}

void setlatlons2_fs_300_(int *isc,float *lat, float *lon, float *sfcPrecip, 
                  float *sfcPrecipStd, float *piaOutKu, float *piaOutKa)
{
  int i;
  extern L2BCMBX_SWATHS swathx300[300];

  for(i=0;i<49;i++)
    {
      swathx300[*isc].FS.Latitude[i]=lat[i];
      swathx300[*isc].FS.Longitude[i]=lon[i];
      if(swathx300[*isc].FS.Longitude[i]>180)
	swathx300[*isc].FS.Longitude[i]-=360;
      swathx300[*isc].FS.surfPrecipTotRate[i]=sfcPrecip[i];
      swathx300[*isc].FS.surfPrecipTotRateSigma[i]=sfcPrecipStd[i];
      swathx300[*isc].FS.pia[i][0]=piaOutKu[i];
      swathx300[*isc].FS.pia[i][1]=piaOutKa[i];
    }
}

void copyrrates1_fs_300_(int *isc,float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<nbins;k++)
    {
      swathx300[*isc].NS.precipTotRate[*i][k]=rrate[k];
      swathx300[*isc].NS.precipTotRateSigma[*i][k]=rratestd[k];
    }
}

void copysflfract_fs_300_(int *isc,float *lfract, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].NS.surfLiqRateFrac[*i]=*lfract;
  if(*i>=12 && *i<=37)
    swathx300[*isc].FS.surfLiqRateFrac[*i]=*lfract;
  
}


//begin  WSO 8/30/13
void copyenvsfqvs1_fs_300_(int *isc,float *envQv, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.surfaceVaporDensity[*i]=envQv[nbins-1];
}

void copyenvsfqvs2_fs_300_(int *isc,float *envQv, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].FS.surfaceVaporDensity[*i]=envQv[nbins-1];
}

void copyenvqvs1_fs_300_(int *isc,float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  for(k=0;k<10;k++)
      swathx300[*isc].NS.vaporDensity[*i][k]=envQv[envnodes[k]-1];
}

void copyenvqvs2_fs_300_(int *isc,float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  for(k=0;k<10;k++)
      swathx300[*isc].FS.vaporDensity[*i][k]=envQv[envnodes[k]-1];
}

void copyenvpresss1_fs_300_(int *isc,float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  for(k=0;k<10;k++)
    swathx300[*isc].NS.airPressure[*i][k]=envQv[envnodes[k]-1];
}

void copyenvpresss2_fs_300_(int *isc,float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  for(k=0;k<10;k++)
    swathx300[*isc].FS.airPressure[*i][k]=envQv[envnodes[k]-1];
}

void copyenvtemps1_fs_300_(int *isc,float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  for(k=0;k<10;k++)
     {
      swathx300[*isc].NS.envParamNode[*i][k]=envnodes[k]-1;
      swathx300[*isc].NS.airTemperature[*i][k]=envQv[envnodes[k]-1];
     }
}

void copyenvtemps2_fs_300_(int *isc,float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  for(k=0;k<10;k++)
    {
     swathx300[*isc].FS.envParamNode[*i][k]=envnodes[k]-1;
     swathx300[*isc].FS.airTemperature[*i][k]=envQv[envnodes[k]-1];
    }
}

void copyenvsftemps1_fs_300_(int *isc,float *envQv, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.surfaceAirTemperature[*i]=envQv[nbins-1];
}

void copyenvsftemps2_fs_300_(int *isc,float *envQv, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].FS.surfaceAirTemperature[*i]=envQv[nbins-1];
}

//end    WSO 8/30/13

void copypwcs1_fs_300_(int *isc,float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<nbins;k++)
    {
      swathx300[*isc].NS.precipTotWaterCont[*i][k]=rrate[k];
      swathx300[*isc].NS.precipTotWaterContSigma[*i][k]=rratestd[k];
    }
}

//begin  WSO 8/7/13
void copylwcfracs1_fs_300_(int *isc,float *mlwc_frac, float *mrate_frac, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<ntransitions;k++)
    {
      swathx300[*isc].NS.liqMassFracTrans[*i][k]=mlwc_frac[k];
      swathx300[*isc].NS.liqRateFracTrans[*i][k]=mrate_frac[k];
    }
}

void copysfcrainliqfracs1_fs_300_(int *isc,float *sfcrainliq_frac, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.surfLiqRateFrac[*i]=*sfcrainliq_frac;

}
//end    WSO 8/7/13

void copyd0s1_fs_300_(int *isc,float *dm, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<nbins;k++)
    {
      swathx300[*isc].NS.precipTotPSDparamHigh[*i][k]=dm[k];
    }
}

void copyzckus1_fs_300_(int *isc,float *zc, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

   for(k=0;k<nbins;k++)
    {
      if(zc[k] > -90.)
        swathx300[*isc].NS.correctedReflectFactor[*i][k] = zc[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swathx300[*isc].NS.correctedReflectFactor[*i][k] = missing_r4c;
//end    WSO 9/17/13
    }
}

void copynodess1_fs_300_(int *isc,int *node, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  for(k=0;k<5;k++)
    {
      swathx300[*isc].NS.phaseBinNodes[*i][k]=node[k];
    }
}

void copyrrates2_fs_300_(int *isc,float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<nbins;k++)
    {
      swathx300[*isc].FS.precipTotRate[*i][k]=rrate[k];
      swathx300[*isc].FS.precipTotRateSigma[*i][k]=rratestd[k];
    }
}

void copypwcs2_fs_300_(int *isc,float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<nbins;k++)
    {
//begin WSO 4/18/2013
//changed NS to MS
      swathx300[*isc].FS.precipTotWaterCont[*i][k]=rrate[k];
      swathx300[*isc].FS.precipTotWaterContSigma[*i][k]=rratestd[k];
//end  WSO 4/18/2013
    }
}

//begin  WSO 8/7/13
void copylwcfracs2_fs_300_(int *isc,float *mlwc_frac, float *mrate_frac, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<ntransitions;k++)
    {
      swathx300[*isc].FS.liqMassFracTrans[*i][k]=mlwc_frac[k];
      swathx300[*isc].FS.liqRateFracTrans[*i][k]=mrate_frac[k];
    }
}

void copysfcrainliqfracs2_fs_300_(int *isc,float *sfcrainliq_frac, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].FS.surfLiqRateFrac[*i]=*sfcrainliq_frac;

}
//end   WSO 8/7/13


void copyd0s2_fs_300_(int *isc,float *dm, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<nbins;k++)
    {
// SFM 05/06/2013 Changed NS to MS to match M.Grecu code from 04/19/2013
      swathx300[*isc].FS.precipTotPSDparamHigh[*i][k]=dm[k];
    }
}

void copyzckus2_fs_300_(int *isc,float *zku, float *zka, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

   for(k=0;k<nbins;k++)
    {
      if(zku[k] > -90.)
        swathx300[*isc].FS.correctedReflectFactor[*i][k][0] = zku[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swathx300[*isc].FS.correctedReflectFactor[*i][k][0] = missing_r4c;
//end    WSO 9/17/13
      if(zka[k] > -90.)
        swathx300[*isc].FS.correctedReflectFactor[*i][k][1] = zka[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swathx300[*isc].FS.correctedReflectFactor[*i][k][1] = missing_r4c;
//end    WSO 9/17/13
    }
}
void copynodess2_fs_300_(int *isc,int *node, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  for(k=0;k<5;k++)
    {
      swathx300[*isc].FS.phaseBinNodes[*i][k]=node[k];
    }
}

//  SFM  begin  12/13/2013; add flag to call sequence
void frominput_fs_300_(int *isc,long *st_2adpr)
{
//  SFM  begin  12/13/2013
  extern TKINFO       granuleHandle2AKu;
  extern L2AKu_NS        L2AKuData;
  extern L2AKuX_FS        L2AKuDataX;
  extern L2BCMBX_SWATHS swathx300[300];
//begin  WSO 9/1/13
  extern L2ADPR_SWATHS dprswath;
  extern L2ADPRX_SWATHS dprxswath;
//end    WSO 9/1/13
  int j;
  int status, status_dpr ;
//for diagnostic
  float dummyPIA[49];
  int k, printPIA[49];
//end for diagnostic
//

//  SFM  begin  12/13/2013; add conditional to dpr read
  status=TKreadScan(&granuleHandle2AKu,&L2AKuDataX);
  if (*st_2adpr == 0) status_dpr=TKreadScan(&dprtkfile,&dprxswath);
  //  return;

//  SFM  begin  12/13/2013

  for( j=0; j<49; j++)
    {
      //swathx300[*isc].NS.Input.piaEffective[j]=L2AKuData.SRT.pathAtten[j];
      swathx300[*isc].NS.Input.piaEffective[j]=L2AKuDataX.SRT.PIAhybrid[j];  //MG  7/31/18, use hybrid PIA
//begin  WSO 9/5/13 remove flag assignment
//       swathx300[*isc].NS.Input.piaEffectiveSigma[j]=-99;
//end    WSO 9/5/13
   //   swathx300[*isc].NS.Input.piaEffectiveReliabFlag[j]=
	// L2AKuData.SRT.reliabFlag[j];
      swathx300[*isc].NS.Input.piaEffectiveReliabFlag[j]=
	L2AKuDataX.SRT.reliabFlagHY[j];                              //WSO  8/2/18 use hybrid flag
      swathx300[*isc].NS.Input.precipitationType[j]=
	L2AKuDataX.CSF.typePrecip[j];
      swathx300[*isc].NS.Input.precipTypeQualityFlag[j]=
	L2AKuDataX.CSF.qualityTypePrecip[j];
      swathx300[*isc].NS.Input.surfaceElevation[j]=L2AKuDataX.PRE.elevation[j];
      swathx300[*isc].NS.Input.localZenithAngle[j]=L2AKuDataX.PRE.localZenithAngle[j];
      swathx300[*isc].NS.Input.surfaceType[j]=L2AKuDataX.PRE.landSurfaceType[j];
//begin  WSO 9/28/13 use alternate rain flag that includes missing for bad scans
//      swathx300[*isc].NS.Input.precipitationFlag[j]=L2AKuData.PRE.flagPrecip[j];
//end    WSO 9/28/13
      swathx300[*isc].NS.Input.surfaceRangeBin[j]=(L2AKuDataX.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
      swathx300[*isc].NS.Input.stormTopBin[j]=(L2AKuDataX.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
      if(swathx300[*isc].NS.Input.stormTopBin[j]<0)
	swathx300[*isc].NS.Input.stormTopBin[j]=missing_i2c;
      swathx300[*isc].NS.Input.stormTopAltitude[j]=L2AKuDataX.PRE.heightStormTop[j];
//begin  WSO 09/30/15 add one bin to the binClutterFreeBottom to temporarily compensate for
//the subtraction of one bin by the radar team
//      swathx300[*isc].NS.Input.lowestClutterFreeBin[j]=
//	(L2AKuData.PRE.binClutterFreeBottom[j]-1)/2; // MG 04/11/2014
//begin  WSO 10/19/15 subtract one 125 m bin from binClutterFreeBottom, and 
//restore V3 definition of lowestClutterFreeBin as a test
//      swathx300[*isc].NS.Input.lowestClutterFreeBin[j]=
//	(L2AKuData.PRE.binClutterFreeBottom[j])/2;
      swathx300[*isc].NS.Input.lowestClutterFreeBin[j]=
	(L2AKuDataX.PRE.binClutterFreeBottom[j] - 2)/2;
//end    WSO 10/15/15
//end    WSO 09/30/15
//begin  WSO 9/17/13 correction for two bin average location in combined
      swathx300[*isc].NS.Input.ellipsoidBinOffset[j]=
	    L2AKuDataX.PRE.ellipsoidBinOffset[j] + 0.125/2.;
//end    WSO 9/17/13
//begin  WSO 8/19/13
      swathx300[*isc].NS.Input.zeroDegAltitude[j] = L2AKuDataX.VER.heightZeroDeg[j];
      swathx300[*isc].NS.Input.zeroDegBin[j] = (L2AKuDataX.VER.binZeroDeg[j]-1)/2; // MG 04/11/2014
//end    WSO 8/19/13
      if(j>=0 && j<49)
	{
	  swathx300[*isc].FS.Input.surfaceElevation[j]=
	    L2AKuDataX.PRE.elevation[j];
	  swathx300[*isc].FS.Input.localZenithAngle[j]=
	    L2AKuDataX.PRE.localZenithAngle[j];
	  swathx300[*isc].FS.Input.surfaceType[j]=
	    L2AKuDataX.PRE.landSurfaceType[j];
	  swathx300[*isc].FS.Input.surfaceRangeBin[j][0]=
	    (L2AKuDataX.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
	  swathx300[*isc].FS.Input.surfaceRangeBin[j][1]=
	    (L2AKuDataX.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
	  swathx300[*isc].FS.Input.stormTopBin[j][0]=
	    (L2AKuDataX.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
	  swathx300[*isc].FS.Input.stormTopBin[j][1]=  // MG 04/11/2014
	    (L2AKuDataX.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
	  if(swathx300[*isc].FS.Input.stormTopBin[j][0]<0)
	    swathx300[*isc].FS.Input.stormTopBin[j][0]=missing_i2c;
	  if(swathx300[*isc].FS.Input.stormTopBin[j][1]<0)
	    swathx300[*isc].FS.Input.stormTopBin[j][1]=missing_i2c;
	  swathx300[*isc].FS.Input.stormTopAltitude[j][0]=
	    L2AKuDataX.PRE.heightStormTop[j];
	  swathx300[*isc].FS.Input.stormTopAltitude[j][1]=
	    L2AKuDataX.PRE.heightStormTop[j];
	  swathx300[*isc].FS.Input.lowestClutterFreeBin[j][0]=
	    (L2AKuDataX.PRE.binClutterFreeBottom[j] - 2)/2;
	  swathx300[*isc].FS.Input.lowestClutterFreeBin[j][1]=
	    (L2AKuDataX.PRE.binClutterFreeBottom[j] - 2)/2;
//end    WSO 10/19/15
//end    WSO 09/30/15
//begin  WSO 9/17/13 correction for two bin average location in combined
	  swathx300[*isc].FS.Input.ellipsoidBinOffset[j][0]=
	    L2AKuDataX.PRE.ellipsoidBinOffset[j] + 0.125/2.;
	  swathx300[*isc].FS.Input.ellipsoidBinOffset[j][1]=
	    L2AKuDataX.PRE.ellipsoidBinOffset[j] + 0.125/2.;
//end    WSO 9/17/13
//begin  WSO 9/5/13 reset pia's using DPR output
	  swathx300[*isc].FS.Input.piaEffective[j][0]=  
	    dprxswath.FS.SRT.pathAtten[j][0];
	  swathx300[*isc].FS.Input.piaEffective[j][1]=
	    dprxswath.FS.SRT.pathAtten[j][1];
//begin  WSO 9/5/13 remove flag assignment
//	  swathx300[*isc].FS.Input.piaEffectiveSigma[j][0]=-99;
//end    WSO 9/5/13
	  swathx300[*isc].FS.Input.piaEffectiveReliabFlag[j][0]=
	    dprxswath.FS.SRT.reliabFlag[j];
	  swathx300[*isc].FS.Input.piaEffectiveReliabFlag[j][1]=
	    dprxswath.FS.SRT.reliabFlag[j];
//end    WSO 9/5/13
	  swathx300[*isc].FS.Input.precipitationType[j]=
	    L2AKuDataX.CSF.typePrecip[j];
	  swathx300[*isc].FS.Input.precipTypeQualityFlag[j]=
	    L2AKuDataX.CSF.qualityTypePrecip[j];
//begin  WSO 8/19/13 need to update toolkit
          swathx300[*isc].FS.Input.zeroDegAltitude[j] =
            L2AKuDataX.VER.heightZeroDeg[j];
          swathx300[*isc].FS.Input.zeroDegBin[j][0] =
            (L2AKuDataX.VER.binZeroDeg[j]-1)/2;   // MG 04/11/2014
//end    WSO 8/19/13
	}
//diagnostic assignment
      //dummyPIA[j] = dprxswath.FS.SRT.pathAtten[j];
//end diagnostic 
    }

//diagnostic
//       if(L2AKuData.Latitude[24] > 30. &&  L2AKuData.Latitude[24] < 40. && L2AKuData.Longitude[24] > -165. && L2AKuData.Longitude[24] <-155.)
//         {
//           for(k=0;k<49;k++)
//             if(dummyPIA[k] < -99.)
//               {
//                 printPIA[k] = 99;
//               }
//             else
//               printPIA[k] = dummyPIA[k]*10.;
//           printf("lon: %10.2f,  ", L2AKuData.Longitude[24]);
//           for(k=0;k<49;k++)
//             printf("%2i", printPIA[k]);
//           printf("\n");
//         }
//end diagnostic

//begin  WSO 9/1/13 scanStatus variables copied from 2AKu
    swathx300[*isc].NS.scanStatus.FractionalGranuleNumber =    
     L2AKuDataX.scanStatus.FractionalGranuleNumber;
    swathx300[*isc].NS.scanStatus.SCorientation =
     L2AKuDataX.scanStatus.SCorientation;
    swathx300[*isc].NS.scanStatus.acsModeMidScan =
     L2AKuDataX.scanStatus.acsModeMidScan;
    swathx300[*isc].NS.scanStatus.dataQuality =
     L2AKuDataX.scanStatus.dataQuality;
    swathx300[*isc].NS.scanStatus.dataWarning =
     L2AKuDataX.scanStatus.dataWarning;
    swathx300[*isc].NS.scanStatus.geoError =
     L2AKuDataX.scanStatus.geoError;
    swathx300[*isc].NS.scanStatus.geoWarning =
     L2AKuDataX.scanStatus.geoWarning;
    swathx300[*isc].NS.scanStatus.limitErrorFlag =
     L2AKuDataX.scanStatus.limitErrorFlag;
    swathx300[*isc].NS.scanStatus.missing =
     L2AKuDataX.scanStatus.missing;
    swathx300[*isc].NS.scanStatus.modeStatus =
     L2AKuDataX.scanStatus.modeStatus;
    swathx300[*isc].NS.scanStatus.operationalMode =
     L2AKuDataX.scanStatus.operationalMode;
    swathx300[*isc].NS.scanStatus.pointingStatus =
     L2AKuDataX.scanStatus.pointingStatus;
    swathx300[*isc].NS.scanStatus.targetSelectionMidScan =
     L2AKuDataX.scanStatus.targetSelectionMidScan;
//from 2ADPRX
    swathx300[*isc].FS.scanStatus.FractionalGranuleNumber =
     dprxswath.FS.scanStatus.FractionalGranuleNumber;
    swathx300[*isc].FS.scanStatus.SCorientation =
     dprxswath.FS.scanStatus.SCorientation;
    swathx300[*isc].FS.scanStatus.acsModeMidScan =
     dprxswath.FS.scanStatus.acsModeMidScan;
    swathx300[*isc].FS.scanStatus.dataQuality =
     dprxswath.FS.scanStatus.dataQuality[1];
    swathx300[*isc].FS.scanStatus.dataWarning =
     dprxswath.FS.scanStatus.dataWarning[1];
    swathx300[*isc].FS.scanStatus.geoError =
     dprxswath.FS.scanStatus.geoError[1];
    swathx300[*isc].FS.scanStatus.geoWarning =
     dprxswath.FS.scanStatus.geoWarning[1];
    swathx300[*isc].FS.scanStatus.limitErrorFlag =
     dprxswath.FS.scanStatus.limitErrorFlag[1];
    swathx300[*isc].FS.scanStatus.missing =
     dprxswath.FS.scanStatus.missing[1];
    swathx300[*isc].FS.scanStatus.modeStatus =
     dprxswath.FS.scanStatus.modeStatus[1];
    swathx300[*isc].FS.scanStatus.operationalMode =
     dprxswath.FS.scanStatus.operationalMode[1];
    swathx300[*isc].FS.scanStatus.pointingStatus =
     dprxswath.FS.scanStatus.pointingStatus[1];
    swathx300[*isc].FS.scanStatus.targetSelectionMidScan =
     dprxswath.FS.scanStatus.targetSelectionMidScan;
//end    WSO 9/1/13

}

void copyscantime_fs_300_(int *isc,int *i)
{
  extern L2BCMBX_SWATHS swathx300[300];
  extern int DayOfMonth[300], DayOfYear[300], Hour[300], MilliSecond[300],
    Minute[300], Month[300], Second[300], Year[300], SecondOfDay[300];
  extern NAVIGATION navigation[300];

 swathx300[*isc].NS.ScanTime.DayOfMonth=DayOfMonth[*i];
 swathx300[*isc].NS.ScanTime.DayOfYear=DayOfYear[*i];
 swathx300[*isc].NS.ScanTime.Hour=Hour[*i];
 swathx300[*isc].NS.ScanTime.MilliSecond=MilliSecond[*i];
 swathx300[*isc].NS.ScanTime.Minute=Minute[*i];
 swathx300[*isc].NS.ScanTime.Month=Month[*i];
 swathx300[*isc].NS.ScanTime.Second=Second[*i];
 swathx300[*isc].NS.ScanTime.Year=Year[*i];
 swathx300[*isc].NS.ScanTime.SecondOfDay=SecondOfDay[*i];
 memcpy(&swathx300[*isc].NS.navigation, &navigation[*i], sizeof(NAVIGATION));

//begin WSO 04/07/2013
//added MS swath scantimes
 swathx300[*isc].FS.ScanTime.DayOfMonth=DayOfMonth[*i];
 swathx300[*isc].FS.ScanTime.DayOfYear=DayOfYear[*i];
 swathx300[*isc].FS.ScanTime.Hour=Hour[*i];
 swathx300[*isc].FS.ScanTime.MilliSecond=MilliSecond[*i];
 swathx300[*isc].FS.ScanTime.Minute=Minute[*i];
 swathx300[*isc].FS.ScanTime.Month=Month[*i];
 swathx300[*isc].FS.ScanTime.Second=Second[*i];
 swathx300[*isc].FS.ScanTime.Year=Year[*i];
 swathx300[*isc].FS.ScanTime.SecondOfDay=SecondOfDay[*i];
 memcpy(&swathx300[*isc].FS.navigation, &navigation[*i], sizeof(NAVIGATION));
//end WSO 04/07/2013
}

void copypreciptype_fs_300_(int *isc,int *ptype, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  //swath.S1.precipitationType[*i]=*ptype;
}

void copyw10_fs_300_(int *isc,float *w10, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.tenMeterWindSpeed[*i]=*w10;
}

void copyw10sigma_fs_300_(int *isc,float *w10s, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.tenMeterWindSigma[*i]=*w10s;
}

void copyw10small_fs_300_(int *isc,float *w10, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].FS.tenMeterWindSpeed[*i]=*w10;
}

void copyw10smallsigma_fs_300_(int *isc,float *w10s, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].FS.tenMeterWindSigma[*i]=*w10s;
}



//  begin  SFM  12/26/2013
void write_empty_fs_300_(int *isc)

//    brief utility to put empty keyword into output file header
//    when needed
{
  char emptygranuletext[100];

  strcpy(emptygranuletext,"EMPTY") ;
  TKsetMetaString(&ctkfile, "FileHeader", "EmptyGranule", 
                  emptygranuletext);	  
}
//  end    SFM  12/26/2013

//  begin  SFM  11/27/2013
void writescan_fs_300_(int *isc)
{
  int ret;
  char emptygranuletext[100];
  extern L2BCMBX_SWATHS swathx300[300];
  // TKgetMetaString(&ctkfile, "FileHeader", "EmptyGranule", 
  //                emptygranuletext);	  
  // if (strncmp(emptygranuletext,"NOT_EMPTY",9) == 0)
  ret= TKwriteScan(&ctkfile,&swathx300[*isc]);
}
//  end    SFM  11/27/2013

void copysfcairtemps1_fs_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.surfaceAirTemperature[*i]=*sfcVar;
}

void copysfcairtemps2_fs_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].FS.surfaceAirTemperature[*i]=*sfcVar;
}

void copysfcairpresss1_fs_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.surfaceAirPressure[*i]=*sfcVar;
}

void copysfcairpresss2_fs_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].FS.surfaceAirPressure[*i]=*sfcVar;
}

void copyskintemps1_fs_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.skinTemperature[*i]=*sfcVar;
}

void copyskintemps2_fs_300_(int *isc,float *sfcVar, int *i)
{  
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].FS.skinTemperature[*i]=*sfcVar;
}

//write skin temperature estimate uncertainty
void copyskintempsigmas1_fs_300_(int *isc,float *skinsigma, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].NS.skinTempSigma[*i] = *skinsigma;
}

void copyskintempsigmas2_fs_300_(int *isc,float *skinsigma, int *i)
{  
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].FS.skinTempSigma[*i] = *skinsigma;
}

//write column vapor estimate uncertainty
void copycolumnvaporsigmas1_fs_300_(int *isc,float *colvaporsigma, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].NS.columnVaporSigma[*i] = *colvaporsigma;
}

void copycolumnvaporsigmas2_fs_300_(int *isc,float *colvaporsigma, int *i)
{  
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].FS.columnVaporSigma[*i] = *colvaporsigma;
}

//write column cloud liquid estimate uncerainty
void copycolumncloudliqsigmas1_fs_300_(int *isc,float *colcldsigma, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  
  swathx300[*isc].NS.columnCloudLiqSigma[*i] = *colcldsigma;
}

void copycolumncloudliqsigmas2_fs_300_(int *isc,float *colcldsigma, int *i)
{  
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].FS.columnCloudLiqSigma[*i] = *colcldsigma;
}

//write algorithm type flag
void copyalgotypes1_fs_300_(int *isc,int *algotype, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].NS.FLG.algoType[*i] = *algotype;
}

void copyalgotypes2_fs_300_(int *isc,int *algotype, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].FS.FLG.algoType[*i] = *algotype;
}


//write error of non-raining data fit
void copyerrorofdatafits1_fs_300_(int *isc,float *erroroffit, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  
  swathx300[*isc].NS.errorOfDataFit[*i] = *erroroffit;
}

void copyerrorofdatafits2_fs_300_(int *isc,float *erroroffit, int *i)
{  
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].FS.errorOfDataFit[*i] = *erroroffit;
}

void copysfcemissouts1_fs_300_(int *isc,float *tbout, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swathx300[*isc].NS.surfEmissivity[*i][k]=tbout[k];
    else
      swathx300[*isc].NS.surfEmissivity[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swathx300[*isc].NS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts1sigma_fs_300_(int *isc,float *tbout, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swathx300[*isc].NS.surfEmissSigma[*i][k]=tbout[k];
    else
      swathx300[*isc].NS.surfEmissSigma[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swathx300[*isc].NS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts2_fs_300_(int *isc,float *tbout, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swathx300[*isc].FS.surfEmissivity[*i][k]=tbout[k];
    else
      swathx300[*isc].FS.surfEmissivity[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//    for(k=13;k<15;k++)
//      swathx300[*isc].FS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts2sigma_fs_300_(int *isc,float *tbout, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swathx300[*isc].FS.surfEmissSigma[*i][k]=tbout[k];
    else
      swathx300[*isc].FS.surfEmissSigma[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//    for(k=13;k<15;k++)
//      swathx300[*isc].FS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copytbouts1_fs_300_(int *isc,float *tbout, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
//begin  WSO 9/16/13
  for(k=0;k<13;k++)
    if(tbout[k] > -90.)
      swathx300[*isc].NS.simulatedBrightTemp[*i][k]=tbout[k];
    else
      swathx300[*isc].NS.simulatedBrightTemp[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels  
//  for(k=13;k<15;k++)
//    swathx300[*isc].NS.simulatedBrightTemp[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
}

void copytbouts2_fs_300_(int *isc,float *tbout, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
//begin  WSO 9/16/13
  for(k=0;k<13;k++)
    if(tbout[k] > -90.)
      swathx300[*isc].FS.simulatedBrightTemp[*i][k]=tbout[k];
    else
      swathx300[*isc].FS.simulatedBrightTemp[*i][k]=missing_r4c;
  //for(k=0;k<2;k++)
  //  swathx300[*isc].FS.simulatedBrightTemp[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swathx300[*isc].FS.simulatedBrightTemp[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end  WSO 9/16/13
}

void copyrainflags1_fs_300_(int *isc,int *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.Input.precipitationFlag[*i]=*sfcVar;
}

void copyrainflags2_fs_300_(int *isc,int *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].FS.Input.precipitationFlag[*i][0]=*sfcVar;
  swathx300[*isc].FS.Input.precipitationFlag[*i][1]=*sfcVar;
}

//begin  WSO 8/20/14 write new ioquality flags
void copyioqualitys1_fs_300_(int *isc,int *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.FLG.ioQuality[*i]=*sfcVar;
}

void copyioqualitys2_fs_300_(int *isc,int *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].FS.FLG.ioQuality[*i]=*sfcVar;
}
//end    WSO 8/20/14
//
//begin  WSO 3/17/17 write snow ice cover flags
void copysnowices1_fs_300_(int *isc,int *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.Input.snowIceCover[*i]=*sfcVar;
}

void copysnowices2_fs_300_(int *isc,int *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].FS.Input.snowIceCover[*i]=*sfcVar;
}
//end    WSO 3/17/17

void copysfcliqfracts1_fs_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.surfLiqRateFrac[*i]=*sfcVar;
}

void copysfcliqfracts2_fs_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].FS.surfLiqRateFrac[*i]=*sfcVar;
}

void copycldwaters1_fs_300_(int *isc,float *var1d, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<nbins;k++)
    {
        swathx300[*isc].NS.cloudLiqWaterCont[*i][k]=var1d[k];
    }
}

void copycldwaters2_fs_300_(int *isc,float *var1d, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300]; 

  for(k=0;k<nbins;k++)
    {
        swathx300[*isc].FS.cloudLiqWaterCont[*i][k]=var1d[k];
    }
}

void copycldices1_fs_300_(int *isc,float *var1d, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<nbins;k++)
    {
      swathx300[*isc].NS.cloudIceWaterCont[*i][k]=var1d[k];
    }
}

void copycldices2_fs_300_(int *isc,float *var1d, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<nbins;k++)
    {
      swathx300[*isc].FS.cloudIceWaterCont[*i][k]=var1d[k];
    }
}

//begin  WSO 9/5/13 new copy routine for SRT and DSRT pia effective sigma's
void copysigmapias1_fs_300_(int *isc,float *sigmapia, int *i)
{
  extern L2BCMBX_SWATHS swathx300[300];
//diagnostic
//  printf("in writeCMBOut i: %5i,  sigmapia: %10.4f\n", *sigmapia, *i);
//end diagnostic
  swathx300[*isc].NS.Input.piaEffectiveSigma[*i] = *sigmapia;
}
void copysigmapias2_fs_300_(int *isc,float *sigmapiaku, float *sigmapiaka, int *i)
{
  extern L2BCMBX_SWATHS swathx300[300];
//    diagnostic
//  printf("in writeCMBOut i: %5i,  sigmapiaku: %10.4f,  sigmapiaka: %10.4f\n", 
//     *sigmapiaku, *sigmapiaka, *i);
//end diagnostic
    swathx300[*isc].FS.Input.piaEffectiveSigma[*i][0] = *sigmapiaku;
    swathx300[*isc].FS.Input.piaEffectiveSigma[*i][1] = *sigmapiaka;
}
//end    WSO 9/5/13

//write principal components
void copyprincomps1_fs_300_(int *isc,float *princomp, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<5;k++)
    {
      swathx300[*isc].NS.aPriori.prinComp[*i][k] = princomp[k];
    }
}
//
void copyprincomps2_fs_300_(int *isc,float *princomp, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<5;k++)
    {
     swathx300[*isc].FS.aPriori.prinComp[*i][k] = princomp[k];
    }
}

//write profile class
void copyprofclasss1_fs_300_(int *isc,int *profclass, int *i)
{
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].NS.aPriori.profClass[*i] = *profclass;
}

void copyprofclasss2_fs_300_(int *isc,int *profclass, int *i)
{
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].FS.aPriori.profClass[*i] = *profclass;
}

//write surface precip bias ratio
void copysurfprecipbiasratios1_fs_300_(int *isc,float *biasratio, int *i)
{ 
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].NS.aPriori.surfPrecipBiasRatio[*i] = *biasratio;
}

void copysurfprecipbiasratios2_fs_300_(int *isc,float *biasratio, int *i)
{
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].FS.aPriori.surfPrecipBiasRatio[*i] = *biasratio;
} 


//write initial log10 of the PSD intercept
void copyinitnws1_fs_300_(int *isc,float *initlogNw, int *n9, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<9;k++)
    {
     if(n9[k]>1 && n9[k]<=88)
       {
         swathx300[*isc].NS.aPriori.initNw[*i][k] = initlogNw[n9[k]-1];
       }
     else
       {
         swathx300[*isc].NS.aPriori.initNw[*i][k] = missing_r4c;
       }
    }
}

void copyinitnws2_fs_300_(int *isc,float *initlogNw, int *n9, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<9;k++)
    {
     if(n9[k]>1 && n9[k]<=88)
       {
         swathx300[*isc].FS.aPriori.initNw[*i][k] = initlogNw[n9[k]-1];
       }
     else
       { 
         swathx300[*isc].FS.aPriori.initNw[*i][k] = missing_r4c;
       }
    }
}

//write sub-footprint variability parameter
void copysubfootvariabilitys1_fs_300_(int *isc,float *subfoot, int *i)
{ 
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].NS.nubfPIAfactor[*i] = *subfoot;
} 
  
void copysubfootvariabilitys2_fs_300_(int *isc,float *subfoot, int *i)
{
  extern L2BCMBX_SWATHS swathx300[300];
    
  swathx300[*isc].FS.nubfPIAfactor[*i] = *subfoot;
}

//write multiple scattering flag
void copymultiscatcalcs1_fs_300_(int *isc,int *multiscat, int *i)
{ 
  extern L2BCMBX_SWATHS swathx300[300];
  
  swathx300[*isc].NS.FLG.multiScatCalc[*i] = *multiscat;
}

void copymultiscatcalcs2_fs_300_(int *isc,int *multiscat, int *i)
{ 
  extern L2BCMBX_SWATHS swathx300[300];
  
  swathx300[*isc].FS.FLG.multiScatCalc[*i] = *multiscat;
}

//write multiple scattering surface parameter
void copymultiscatsurfaces1_fs_300_(int *isc,float *multisfc, int *i)
{
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].NS.multiScatMaxContrib[*i] = *multisfc;
}

void copymultiscatsurfaces2_fs_300_(int *isc,float *multisfc, int *i)
{
  extern L2BCMBX_SWATHS swathx300[300];

  swathx300[*isc].FS.multiScatMaxContrib[*i] = *multisfc;
}

//
//begin  WSO 2/8/17 copy routine for measured sigma-zeros
void copysigmazeros1_fs_300_(int *isc,float *sigmazeroku, int *i)
{
  extern L2BCMBX_SWATHS swathx300[300];
  swathx300[*isc].NS.Input.sigmaZeroMeasured[*i] = *sigmazeroku;
}
void copysigmazeros2_fs_300_(int *isc,float *sigmazeroku, float *sigmazeroka, int *i)
{
    extern L2BCMBX_SWATHS swathx300[300];
//        swathx300[*isc].FS.Input.sigmaZeroMeasured[*i][0] = *sigmazeroku;
//            swathx300[*isc].FS.Input.sigmaZeroMeasured[*i][1] = *sigmazeroka;
          swathx300[*isc].FS.Input.sigmaZeroMeasured[*i] = *sigmazeroka;
}
//end    WSO 2/8/17

//begin  WSO 8/19/13 modified copy routines to include nodes
void copylognws1_fs_300_(int *isc,float *logNw, int *n9, int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<9;k++)
    {
//begin WSO 9/7/14 added upper limit for n9 in next line
      if(n9[k]>1 && n9[k]<=88)
        {
         swathx300[*isc].NS.PSDparamLowNode[*i][k] = n9[k]-1;
	 swathx300[*isc].NS.precipTotPSDparamLow[*i][k][0]=logNw[n9[k]-1];
        }
      else
        {
         swathx300[*isc].NS.PSDparamLowNode[*i][k] = missing_i2c;
         swathx300[*isc].NS.precipTotPSDparamLow[*i][k][0]= missing_r4c;
        }
    }
}

void copylognws2_fs_300_(int *isc,float *logNw, int *n9,int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];
  
  for(k=0;k<9;k++)
    {
      if(n9[k]>0)
        {
         swathx300[*isc].FS.PSDparamLowNode[*i][k] = n9[k]-1;
	 swathx300[*isc].FS.precipTotPSDparamLow[*i][k][0]=logNw[n9[k]-1];
        }
      else
        {
         swathx300[*isc].FS.PSDparamLowNode[*i][k] = missing_i2c;
	 swathx300[*isc].FS.precipTotPSDparamLow[*i][k][0]= missing_r4c;
        }
    }
}
//end    WSO 8/19/13

//begin  WSO 8/19/13 add mu as second low-resolution parameter
void copymus1_fs_300_(int *isc,float *mu_prof, int *n9,int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<9;k++)
    {
      if(n9[k]>1)
        {
         swathx300[*isc].NS.precipTotPSDparamLow[*i][k][1]=mu_prof[n9[k]-1];
        }
      else
        {
         swathx300[*isc].NS.precipTotPSDparamLow[*i][k][1]=missing_r4c;
        }
    }
}

void copymus2_fs_300_(int *isc,float *mu_prof, int *n9,int *i)
{
  int k;
  extern L2BCMBX_SWATHS swathx300[300];

  for(k=0;k<9;k++)
    {
      if(n9[k]>0)
        {
         swathx300[*isc].FS.precipTotPSDparamLow[*i][k][1]=mu_prof[n9[k]-1];
        }
      else
        {
         swathx300[*isc].FS.precipTotPSDparamLow[*i][k][1]=missing_r4c;
        }
    }
}
//end    WSO 8/19/13



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

//L2BCMB_SWATHS swath;
//L2BCMB_SWATHS swath300[300];
//L2BCMBX_SWATHS swathx300[300];
L2ADPR_SWATHS dprswath;
L2ADPRX_SWATHS dprxswath;
L2BCMB_SWATHS swath1;
L2AKu_NS      L2AKuData;
L2AKuX_FS      L2AKuDataX;


void setlatlons1_300_(int *isc,float *lat, float *lon, float *sfcPrecip, 
                  float *sfcPrecipStd, float *piaOut)
{
  int i;
  extern L2BCMB_SWATHS swath300[300];

  for(i=0;i<49;i++)
    {
      swath300[*isc].NS.Latitude[i]=lat[i];
      swath300[*isc].NS.Longitude[i]=lon[i];
      if(swath300[*isc].NS.Longitude[i]>180)
	swath300[*isc].NS.Longitude[i]-=360;
      swath300[*isc].NS.surfPrecipTotRate[i]=sfcPrecip[i];
      swath300[*isc].NS.surfPrecipTotRateSigma[i]=sfcPrecipStd[i];
      swath300[*isc].NS.pia[i]=piaOut[i];
    }
}

void setlatlons2_300_(int *isc,float *lat, float *lon, float *sfcPrecip, 
                  float *sfcPrecipStd, float *piaOutKu, float *piaOutKa)
{
  int i;
  extern L2BCMB_SWATHS swath300[300];

  for(i=12;i<37;i++)
    {
      swath300[*isc].MS.Latitude[i-12]=lat[i];
      swath300[*isc].MS.Longitude[i-12]=lon[i];
      if(swath300[*isc].MS.Longitude[i-12]>180)
	swath300[*isc].MS.Longitude[i-12]-=360;
      swath300[*isc].MS.surfPrecipTotRate[i-12]=sfcPrecip[i];
      swath300[*isc].MS.surfPrecipTotRateSigma[i-12]=sfcPrecipStd[i];
      swath300[*isc].MS.pia[i-12][0]=piaOutKu[i];
      swath300[*isc].MS.pia[i-12][1]=piaOutKa[i];
    }
}

void copyrrates1_300_(int *isc,float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<nbins;k++)
    {
      swath300[*isc].NS.precipTotRate[*i][k]=rrate[k];
      swath300[*isc].NS.precipTotRateSigma[*i][k]=rratestd[k];
    }
}

void copysflfract_300_(int *isc,float *lfract, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].NS.surfLiqRateFrac[*i]=*lfract;
  if(*i>=12 && *i<=37)
    swath300[*isc].MS.surfLiqRateFrac[*i-12]=*lfract;
  
}

void copyzka_300_(int *isc,float *zka, int *i)
{
  int k;
  extern L2ADPR_SWATHS dprswath;

  for(k=0;k<nbins;k++)
    {
      dprswath.MS.PRE.zFactorMeasured[*i][2*k]=(int)(zka[k]*100);
      dprswath.MS.PRE.zFactorMeasured[*i][2*k+1]=(int)(zka[k]*100);
      //      printf("%g ",zka[k]);
    }
}

void copypiaka_300_(int *isc,float *piaKa, int *i)
{
  int k;
  extern L2ADPR_SWATHS dprswath;
 
  for(k=0;k<nbins;k++)
    {
      dprswath.MS.SRT.pathAtten[*i]=*piaKa;
    }
}

void copytruerrate_300_(int *isc,float *rrate, int *i)
{
  int k;
  extern L2ADPR_SWATHS dprswath;

  for(k=0;k<nbins;k++)
    {
      dprswath.NS.SLV.precipRate[*i][2*k]=(int)(rrate[k]*100);
      dprswath.NS.SLV.precipRate[*i][2*k+1]=(int)(rrate[k]*100);
    }
}

//begin  WSO 8/30/13
void copyenvsfqvs1_300_(int *isc,float *envQv, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.surfaceVaporDensity[*i]=envQv[nbins-1];
}

void copyenvsfqvs2_300_(int *isc,float *envQv, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].MS.surfaceVaporDensity[*i]=envQv[nbins-1];
}

void copyenvqvs1_300_(int *isc,float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  for(k=0;k<10;k++)
      swath300[*isc].NS.vaporDensity[*i][k]=envQv[envnodes[k]-1];
}

void copyenvqvs2_300_(int *isc,float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  for(k=0;k<10;k++)
      swath300[*isc].MS.vaporDensity[*i][k]=envQv[envnodes[k]-1];
}

void copyenvpresss1_300_(int *isc,float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  for(k=0;k<10;k++)
    swath300[*isc].NS.airPressure[*i][k]=envQv[envnodes[k]-1];
}

void copyenvpresss2_300_(int *isc,float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  for(k=0;k<10;k++)
    swath300[*isc].MS.airPressure[*i][k]=envQv[envnodes[k]-1];
}

void copyenvtemps1_300_(int *isc,float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  for(k=0;k<10;k++)
     {
      swath300[*isc].NS.envParamNode[*i][k]=envnodes[k]-1;
      swath300[*isc].NS.airTemperature[*i][k]=envQv[envnodes[k]-1];
     }
}

void copyenvtemps2_300_(int *isc,float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  for(k=0;k<10;k++)
    {
     swath300[*isc].MS.envParamNode[*i][k]=envnodes[k]-1;
     swath300[*isc].MS.airTemperature[*i][k]=envQv[envnodes[k]-1];
    }
}

void copyenvsftemps1_300_(int *isc,float *envQv, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.surfaceAirTemperature[*i]=envQv[nbins-1];
}

void copyenvsftemps2_300_(int *isc,float *envQv, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].MS.surfaceAirTemperature[*i]=envQv[nbins-1];
}

//end    WSO 8/30/13

void copypwcs1_300_(int *isc,float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<nbins;k++)
    {
      swath300[*isc].NS.precipTotWaterCont[*i][k]=rrate[k];
      swath300[*isc].NS.precipTotWaterContSigma[*i][k]=rratestd[k];
    }
}

//begin  WSO 8/7/13
void copylwcfracs1_300_(int *isc,float *mlwc_frac, float *mrate_frac, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<ntransitions;k++)
    {
      swath300[*isc].NS.liqMassFracTrans[*i][k]=mlwc_frac[k];
      swath300[*isc].NS.liqRateFracTrans[*i][k]=mrate_frac[k];
    }
}

void copysfcrainliqfracs1_300_(int *isc,float *sfcrainliq_frac, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.surfLiqRateFrac[*i]=*sfcrainliq_frac;

}
//end    WSO 8/7/13

void copyd0s1_300_(int *isc,float *dm, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<nbins;k++)
    {
      swath300[*isc].NS.precipTotPSDparamHigh[*i][k]=dm[k];
    }
}

void copyzckus1_300_(int *isc,float *zc, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

   for(k=0;k<nbins;k++)
    {
      if(zc[k] > -90.)
        swath300[*isc].NS.correctedReflectFactor[*i][k] = zc[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swath300[*isc].NS.correctedReflectFactor[*i][k] = missing_r4c;
//end    WSO 9/17/13
    }
}

void copynodess1_300_(int *isc,int *node, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  for(k=0;k<5;k++)
    {
      swath300[*isc].NS.phaseBinNodes[*i][k]=node[k];
    }
}

void copyrrates2_300_(int *isc,float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<nbins;k++)
    {
      swath300[*isc].MS.precipTotRate[*i][k]=rrate[k];
      swath300[*isc].MS.precipTotRateSigma[*i][k]=rratestd[k];
    }
}

void copypwcs2_300_(int *isc,float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<nbins;k++)
    {
//begin WSO 4/18/2013
//changed NS to MS
      swath300[*isc].MS.precipTotWaterCont[*i][k]=rrate[k];
      swath300[*isc].MS.precipTotWaterContSigma[*i][k]=rratestd[k];
//end  WSO 4/18/2013
    }
}

//begin  WSO 8/7/13
void copylwcfracs2_300_(int *isc,float *mlwc_frac, float *mrate_frac, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<ntransitions;k++)
    {
      swath300[*isc].MS.liqMassFracTrans[*i][k]=mlwc_frac[k];
      swath300[*isc].MS.liqRateFracTrans[*i][k]=mrate_frac[k];
    }
}

void copysfcrainliqfracs2_300_(int *isc,float *sfcrainliq_frac, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].MS.surfLiqRateFrac[*i]=*sfcrainliq_frac;

}
//end   WSO 8/7/13


void copyd0s2_300_(int *isc,float *dm, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<nbins;k++)
    {
// SFM 05/06/2013 Changed NS to MS to match M.Grecu code from 04/19/2013
      swath300[*isc].MS.precipTotPSDparamHigh[*i][k]=dm[k];
    }
}

void copyzckus2_300_(int *isc,float *zku, float *zka, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

   for(k=0;k<nbins;k++)
    {
      if(zku[k] > -90.)
        swath300[*isc].MS.correctedReflectFactor[*i][k][0] = zku[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swath300[*isc].MS.correctedReflectFactor[*i][k][0] = missing_r4c;
//end    WSO 9/17/13
      if(zka[k] > -90.)
        swath300[*isc].MS.correctedReflectFactor[*i][k][1] = zka[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swath300[*isc].MS.correctedReflectFactor[*i][k][1] = missing_r4c;
//end    WSO 9/17/13
    }
}
void copynodess2_300_(int *isc,int *node, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  for(k=0;k<5;k++)
    {
      swath300[*isc].MS.phaseBinNodes[*i][k]=node[k];
    }
}

void rewind_300_(int *isc,int *ic)
{
  extern TKINFO       granuleHandle2AKu;
  int status = TKseek(&granuleHandle2AKu, *ic, TK_ABS_SCAN_OFF); 
}

void rewindc_300_(int *isc,int *ic)
{
  extern TKINFO       granuleHandle2AKu;
  printf("rewind cmb\n");
  int status = TKseek(&ctkfile, *ic, TK_ABS_SCAN_OFF); 
}

//begin WSO 9/8/13 rewind DPR file
void rewind_dpr_300_(int *isc,int *ic)
{
  extern TKINFO       dprtkfile;
    int status_dpr = TKseek(&dprtkfile, *ic, TK_ABS_SCAN_OFF);
}
//end WSO 9/8/13

//  SFM  begin  12/13/2013; add flag to call sequence
void frominput_300_(int *isc,long *st_2adpr)
{
//  SFM  begin  12/13/2013
  extern TKINFO       granuleHandle2AKu;
  extern L2AKu_NS        L2AKuData;
  extern L2BCMB_SWATHS swath300[300];
//begin  WSO 9/1/13
  extern L2ADPR_SWATHS dprswath;
//end    WSO 9/1/13
  int j;
  int status, status_dpr ;
//for diagnostic
  float dummyPIA[49];
  int k, printPIA[49];
//end for diagnostic
//

//  SFM  begin  12/13/2013; add conditional to dpr read
  status=TKreadScan(&granuleHandle2AKu,&L2AKuData);
  if (*st_2adpr == 0) status_dpr=TKreadScan(&dprtkfile,&dprswath);
//  SFM  begin  12/13/2013

  for( j=0; j<49; j++)
    {
      //swath300[*isc].NS.Input.piaEffective[j]=L2AKuData.SRT.pathAtten[j];
      swath300[*isc].NS.Input.piaEffective[j]=L2AKuData.SRT.PIAhybrid[j];  //MG  7/31/18, use hybrid PIA
//begin  WSO 9/5/13 remove flag assignment
//       swath300[*isc].NS.Input.piaEffectiveSigma[j]=-99;
//end    WSO 9/5/13
   //   swath300[*isc].NS.Input.piaEffectiveReliabFlag[j]=
	// L2AKuData.SRT.reliabFlag[j];
      swath300[*isc].NS.Input.piaEffectiveReliabFlag[j]=
	L2AKuData.SRT.reliabFlagHY[j];                              //WSO  8/2/18 use hybrid flag
      swath300[*isc].NS.Input.precipitationType[j]=
	L2AKuData.CSF.typePrecip[j];
      swath300[*isc].NS.Input.precipTypeQualityFlag[j]=
	L2AKuData.CSF.qualityTypePrecip[j];
      swath300[*isc].NS.Input.surfaceElevation[j]=L2AKuData.PRE.elevation[j];
      swath300[*isc].NS.Input.localZenithAngle[j]=L2AKuData.PRE.localZenithAngle[j];
      swath300[*isc].NS.Input.surfaceType[j]=L2AKuData.PRE.landSurfaceType[j];
//begin  WSO 9/28/13 use alternate rain flag that includes missing for bad scans
//      swath300[*isc].NS.Input.precipitationFlag[j]=L2AKuData.PRE.flagPrecip[j];
//end    WSO 9/28/13
      swath300[*isc].NS.Input.surfaceRangeBin[j]=(L2AKuData.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
      swath300[*isc].NS.Input.stormTopBin[j]=(L2AKuData.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
      if(swath300[*isc].NS.Input.stormTopBin[j]<0)
	swath300[*isc].NS.Input.stormTopBin[j]=missing_i2c;
      swath300[*isc].NS.Input.stormTopAltitude[j]=L2AKuData.PRE.heightStormTop[j];
//begin  WSO 09/30/15 add one bin to the binClutterFreeBottom to temporarily compensate for
//the subtraction of one bin by the radar team
//      swath300[*isc].NS.Input.lowestClutterFreeBin[j]=
//	(L2AKuData.PRE.binClutterFreeBottom[j]-1)/2; // MG 04/11/2014
//begin  WSO 10/19/15 subtract one 125 m bin from binClutterFreeBottom, and 
//restore V3 definition of lowestClutterFreeBin as a test
//      swath300[*isc].NS.Input.lowestClutterFreeBin[j]=
//	(L2AKuData.PRE.binClutterFreeBottom[j])/2;
      swath300[*isc].NS.Input.lowestClutterFreeBin[j]=
	(L2AKuData.PRE.binClutterFreeBottom[j] - 2)/2;
//end    WSO 10/15/15
//end    WSO 09/30/15
//begin  WSO 9/17/13 correction for two bin average location in combined
      swath300[*isc].NS.Input.ellipsoidBinOffset[j]=
	    L2AKuData.PRE.ellipsoidBinOffset[j] + 0.125/2.;
//end    WSO 9/17/13
//begin  WSO 8/19/13
      swath300[*isc].NS.Input.zeroDegAltitude[j] = L2AKuData.VER.heightZeroDeg[j];
      swath300[*isc].NS.Input.zeroDegBin[j] = (L2AKuData.VER.binZeroDeg[j]-1)/2; // MG 04/11/2014
//end    WSO 8/19/13
      if(j>=12 && j<37)
	{
	  swath300[*isc].MS.Input.surfaceElevation[j-12]=
	    L2AKuData.PRE.elevation[j];
	  swath300[*isc].MS.Input.localZenithAngle[j-12]=
	    L2AKuData.PRE.localZenithAngle[j];
	  swath300[*isc].MS.Input.surfaceType[j-12]=
	    L2AKuData.PRE.landSurfaceType[j];
//begin  WSO 9/28/13 use alternate rain flag that includes missing for bad scans
//	  swath300[*isc].MS.Input.precipitationFlag[j-12][0]=
//	    L2AKuData.PRE.flagPrecip[j];
//	  swath300[*isc].MS.Input.precipitationFlag[j-12][1]=
//	    L2AKuData.PRE.flagPrecip[j];
//end    WSO 9/28/13
	  swath300[*isc].MS.Input.surfaceRangeBin[j-12][0]=
	    (L2AKuData.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
	  swath300[*isc].MS.Input.surfaceRangeBin[j-12][1]=
	    (L2AKuData.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
	  swath300[*isc].MS.Input.stormTopBin[j-12][0]=
	    (L2AKuData.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
	  swath300[*isc].MS.Input.stormTopBin[j-12][1]=  // MG 04/11/2014
	    (L2AKuData.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
	  if(swath300[*isc].MS.Input.stormTopBin[j-12][0]<0)
	    swath300[*isc].MS.Input.stormTopBin[j-12][0]=missing_i2c;
	  if(swath300[*isc].MS.Input.stormTopBin[j-12][1]<0)
	    swath300[*isc].MS.Input.stormTopBin[j-12][1]=missing_i2c;
	  swath300[*isc].MS.Input.stormTopAltitude[j-12][0]=
	    L2AKuData.PRE.heightStormTop[j];
	  swath300[*isc].MS.Input.stormTopAltitude[j-12][1]=
	    L2AKuData.PRE.heightStormTop[j];
//begin  WSO 09/30/15 add one bin to the binClutterFreeBottom to temporarily compensate for
//the subtraction of one bin by the radar team
//	  swath300[*isc].MS.Input.lowestClutterFreeBin[j-12][0]=
//    (L2AKuData.PRE.binClutterFreeBottom[j]-1)/2; // MG 04/11/2014
//	  swath300[*isc].MS.Input.lowestClutterFreeBin[j-12][1]=
//	    (L2AKuData.PRE.binClutterFreeBottom[j]-1)/2; // MG 04/11/2014
//begin  WSO 10/19/15 subtract one 125 m bin from binClutterFreeBottom, and 
// restore V3 definition of lowestClutterFreeBin in test
//	  swath300[*isc].MS.Input.lowestClutterFreeBin[j-12][0]=
//	    (L2AKuData.PRE.binClutterFreeBottom[j])/2;
//	  swath300[*isc].MS.Input.lowestClutterFreeBin[j-12][1]=
//	    (L2AKuData.PRE.binClutterFreeBottom[j])/2;
	  swath300[*isc].MS.Input.lowestClutterFreeBin[j-12][0]=
	    (L2AKuData.PRE.binClutterFreeBottom[j] - 2)/2;
	  swath300[*isc].MS.Input.lowestClutterFreeBin[j-12][1]=
	    (L2AKuData.PRE.binClutterFreeBottom[j] - 2)/2;
//end    WSO 10/19/15
//end    WSO 09/30/15
//begin  WSO 9/17/13 correction for two bin average location in combined
	  swath300[*isc].MS.Input.ellipsoidBinOffset[j-12][0]=
	    L2AKuData.PRE.ellipsoidBinOffset[j] + 0.125/2.;
	  swath300[*isc].MS.Input.ellipsoidBinOffset[j-12][1]=
	    L2AKuData.PRE.ellipsoidBinOffset[j] + 0.125/2.;
//end    WSO 9/17/13
//begin  WSO 9/5/13 reset pia's using DPR output
	  swath300[*isc].MS.Input.piaEffective[j-12][0]=  
	    dprswath.NS.SRT.pathAtten[j];
	  swath300[*isc].MS.Input.piaEffective[j-12][1]=
	    dprswath.MS.SRT.pathAtten[j-12];
//begin  WSO 9/5/13 remove flag assignment
//	  swath300[*isc].MS.Input.piaEffectiveSigma[j-12][0]=-99;
//end    WSO 9/5/13
	  swath300[*isc].MS.Input.piaEffectiveReliabFlag[j-12][0]=
	    dprswath.NS.SRT.reliabFlag[j];
	  swath300[*isc].MS.Input.piaEffectiveReliabFlag[j-12][1]=
	    dprswath.MS.SRT.reliabFlag[j-12];
//end    WSO 9/5/13
	  swath300[*isc].MS.Input.precipitationType[j-12]=
	    L2AKuData.CSF.typePrecip[j];
	  swath300[*isc].MS.Input.precipTypeQualityFlag[j-12]=
	    L2AKuData.CSF.qualityTypePrecip[j];
//begin  WSO 8/19/13 need to update toolkit
          swath300[*isc].MS.Input.zeroDegAltitude[j-12] =
            L2AKuData.VER.heightZeroDeg[j];
          swath300[*isc].MS.Input.zeroDegBin[j-12][0] =
            (L2AKuData.VER.binZeroDeg[j]-1)/2;   // MG 04/11/2014
//end    WSO 8/19/13
	}
//diagnostic assignment
          dummyPIA[j] = dprswath.NS.SRT.pathAtten[j];
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
    swath300[*isc].NS.scanStatus.FractionalGranuleNumber =    
     L2AKuData.scanStatus.FractionalGranuleNumber;
    swath300[*isc].NS.scanStatus.SCorientation =
     L2AKuData.scanStatus.SCorientation;
    swath300[*isc].NS.scanStatus.acsModeMidScan =
     L2AKuData.scanStatus.acsModeMidScan;
    swath300[*isc].NS.scanStatus.dataQuality =
     L2AKuData.scanStatus.dataQuality;
    swath300[*isc].NS.scanStatus.dataWarning =
     L2AKuData.scanStatus.dataWarning;
    swath300[*isc].NS.scanStatus.geoError =
     L2AKuData.scanStatus.geoError;
    swath300[*isc].NS.scanStatus.geoWarning =
     L2AKuData.scanStatus.geoWarning;
    swath300[*isc].NS.scanStatus.limitErrorFlag =
     L2AKuData.scanStatus.limitErrorFlag;
    swath300[*isc].NS.scanStatus.missing =
     L2AKuData.scanStatus.missing;
    swath300[*isc].NS.scanStatus.modeStatus =
     L2AKuData.scanStatus.modeStatus;
    swath300[*isc].NS.scanStatus.operationalMode =
     L2AKuData.scanStatus.operationalMode;
    swath300[*isc].NS.scanStatus.pointingStatus =
     L2AKuData.scanStatus.pointingStatus;
    swath300[*isc].NS.scanStatus.targetSelectionMidScan =
     L2AKuData.scanStatus.targetSelectionMidScan;
//from 2ADPR
    swath300[*isc].MS.scanStatus.FractionalGranuleNumber =
     dprswath.MS.scanStatus.FractionalGranuleNumber;
    swath300[*isc].MS.scanStatus.SCorientation =
     dprswath.MS.scanStatus.SCorientation;
    swath300[*isc].MS.scanStatus.acsModeMidScan =
     dprswath.MS.scanStatus.acsModeMidScan;
    swath300[*isc].MS.scanStatus.dataQuality =
     dprswath.MS.scanStatus.dataQuality;
    swath300[*isc].MS.scanStatus.dataWarning =
     dprswath.MS.scanStatus.dataWarning;
    swath300[*isc].MS.scanStatus.geoError =
     dprswath.MS.scanStatus.geoError;
    swath300[*isc].MS.scanStatus.geoWarning =
     dprswath.MS.scanStatus.geoWarning;
    swath300[*isc].MS.scanStatus.limitErrorFlag =
     dprswath.MS.scanStatus.limitErrorFlag;
    swath300[*isc].MS.scanStatus.missing =
     dprswath.MS.scanStatus.missing;
    swath300[*isc].MS.scanStatus.modeStatus =
     dprswath.MS.scanStatus.modeStatus;
    swath300[*isc].MS.scanStatus.operationalMode =
     dprswath.MS.scanStatus.operationalMode;
    swath300[*isc].MS.scanStatus.pointingStatus =
     dprswath.MS.scanStatus.pointingStatus;
    swath300[*isc].MS.scanStatus.targetSelectionMidScan =
     dprswath.MS.scanStatus.targetSelectionMidScan;
//end    WSO 9/1/13

}

void copyscantime_300_(int *isc,int *i)
{
  extern L2BCMB_SWATHS swath300[300];
  extern int DayOfMonth[300], DayOfYear[300], Hour[300], MilliSecond[300],
    Minute[300], Month[300], Second[300], Year[300], SecondOfDay[300];
  extern NAVIGATION navigation[300];

 swath300[*isc].NS.ScanTime.DayOfMonth=DayOfMonth[*i];
 swath300[*isc].NS.ScanTime.DayOfYear=DayOfYear[*i];
 swath300[*isc].NS.ScanTime.Hour=Hour[*i];
 swath300[*isc].NS.ScanTime.MilliSecond=MilliSecond[*i];
 swath300[*isc].NS.ScanTime.Minute=Minute[*i];
 swath300[*isc].NS.ScanTime.Month=Month[*i];
 swath300[*isc].NS.ScanTime.Second=Second[*i];
 swath300[*isc].NS.ScanTime.Year=Year[*i];
 swath300[*isc].NS.ScanTime.SecondOfDay=SecondOfDay[*i];
 memcpy(&swath300[*isc].NS.navigation, &navigation[*i], sizeof(NAVIGATION));

//begin WSO 04/07/2013
//added MS swath300[300] scantimes
 swath300[*isc].MS.ScanTime.DayOfMonth=DayOfMonth[*i];
 swath300[*isc].MS.ScanTime.DayOfYear=DayOfYear[*i];
 swath300[*isc].MS.ScanTime.Hour=Hour[*i];
 swath300[*isc].MS.ScanTime.MilliSecond=MilliSecond[*i];
 swath300[*isc].MS.ScanTime.Minute=Minute[*i];
 swath300[*isc].MS.ScanTime.Month=Month[*i];
 swath300[*isc].MS.ScanTime.Second=Second[*i];
 swath300[*isc].MS.ScanTime.Year=Year[*i];
 swath300[*isc].MS.ScanTime.SecondOfDay=SecondOfDay[*i];
 memcpy(&swath300[*isc].MS.navigation, &navigation[*i], sizeof(NAVIGATION));
//end WSO 04/07/2013
}

void copypreciptype_300_(int *isc,int *ptype, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  //swath300[*isc].S1.precipitationType[*i]=*ptype;
}

void copyw10_300_(int *isc,float *w10, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.tenMeterWindSpeed[*i]=*w10;
}

void copyw10sigma_300_(int *isc,float *w10s, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.tenMeterWindSigma[*i]=*w10s;
}

void copyw10small_300_(int *isc,float *w10, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].MS.tenMeterWindSpeed[*i]=*w10;
}

void copyw10smallsigma_300_(int *isc,float *w10s, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].MS.tenMeterWindSigma[*i]=*w10s;
}

void writedprscan_300_(int *isc)
{
  int ret;
  ret= TKwriteScan(&dprtkfile,&dprswath);
}

//  begin  SFM  12/26/2013
void write_empty_300_(int *isc)

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
void writescan_300_(int *isc)
{
  int ret;
  char emptygranuletext[100];
  extern L2BCMB_SWATHS swath300[300];
  // TKgetMetaString(&ctkfile, "FileHeader", "EmptyGranule", 
  //                emptygranuletext);	  
  // if (strncmp(emptygranuletext,"NOT_EMPTY",9) == 0)
  ret= TKwriteScan(&ctkfile,&swath300);
}
//  end    SFM  11/27/2013

void copysfcairtemps1_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.surfaceAirTemperature[*i]=*sfcVar;
}

void copysfcairtemps2_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].MS.surfaceAirTemperature[*i]=*sfcVar;
}

void copysfcairpresss1_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.surfaceAirPressure[*i]=*sfcVar;
}

void copysfcairpresss2_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].MS.surfaceAirPressure[*i]=*sfcVar;
}

void copyskintemps1_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.skinTemperature[*i]=*sfcVar;
}

void copyskintemps2_300_(int *isc,float *sfcVar, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].MS.skinTemperature[*i]=*sfcVar;
}

//write skin temperature estimate uncertainty
void copyskintempsigmas1_300_(int *isc,float *skinsigma, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].NS.skinTempSigma[*i] = *skinsigma;
}

void copyskintempsigmas2_300_(int *isc,float *skinsigma, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].MS.skinTempSigma[*i] = *skinsigma;
}

//write column vapor estimate uncertainty
void copycolumnvaporsigmas1_300_(int *isc,float *colvaporsigma, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].NS.columnVaporSigma[*i] = *colvaporsigma;
}

void copycolumnvaporsigmas2_300_(int *isc,float *colvaporsigma, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].MS.columnVaporSigma[*i] = *colvaporsigma;
}

//write column cloud liquid estimate uncerainty
void copycolumncloudliqsigmas1_300_(int *isc,float *colcldsigma, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  
  swath300[*isc].NS.columnCloudLiqSigma[*i] = *colcldsigma;
}

void copycolumncloudliqsigmas2_300_(int *isc,float *colcldsigma, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].MS.columnCloudLiqSigma[*i] = *colcldsigma;
}

//write algorithm type flag
void copyalgotypes1_300_(int *isc,int *algotype, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].NS.FLG.algoType[*i] = *algotype;
}

void copyalgotypes2_300_(int *isc,int *algotype, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].MS.FLG.algoType[*i] = *algotype;
}


//write error of non-raining data fit
void copyerrorofdatafits1_300_(int *isc,float *erroroffit, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  
  swath300[*isc].NS.errorOfDataFit[*i] = *erroroffit;
}

void copyerrorofdatafits2_300_(int *isc,float *erroroffit, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].MS.errorOfDataFit[*i] = *erroroffit;
}

void copysfcemissouts1_300_(int *isc,float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swath300[*isc].NS.surfEmissivity[*i][k]=tbout[k];
    else
      swath300[*isc].NS.surfEmissivity[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swath300[*isc].NS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts1sigma_300_(int *isc,float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swath300[*isc].NS.surfEmissSigma[*i][k]=tbout[k];
    else
      swath300[*isc].NS.surfEmissSigma[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swath300[*isc].NS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts2_300_(int *isc,float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swath300[*isc].MS.surfEmissivity[*i][k]=tbout[k];
    else
      swath300[*isc].MS.surfEmissivity[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//    for(k=13;k<15;k++)
//      swath300[*isc].MS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts2sigma_300_(int *isc,float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swath300[*isc].MS.surfEmissSigma[*i][k]=tbout[k];
    else
      swath300[*isc].MS.surfEmissSigma[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//    for(k=13;k<15;k++)
//      swath300[*isc].MS.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copytbouts1_300_(int *isc,float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
//begin  WSO 9/16/13
  for(k=0;k<13;k++)
    if(tbout[k] > -90.)
      swath300[*isc].NS.simulatedBrightTemp[*i][k]=tbout[k];
    else
      swath300[*isc].NS.simulatedBrightTemp[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels  
//  for(k=13;k<15;k++)
//    swath300[*isc].NS.simulatedBrightTemp[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
}

void copytbouts2_300_(int *isc,float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
//begin  WSO 9/16/13
  for(k=0;k<13;k++)
    if(tbout[k] > -90.)
      swath300[*isc].MS.simulatedBrightTemp[*i][k]=tbout[k];
    else
      swath300[*isc].MS.simulatedBrightTemp[*i][k]=missing_r4c;
  //for(k=0;k<2;k++)
  //  swath300[*isc].MS.simulatedBrightTemp[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swath300[*isc].MS.simulatedBrightTemp[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end  WSO 9/16/13
}

void copyrainflags1_300_(int *isc,int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.Input.precipitationFlag[*i]=*sfcVar;
}

void copyrainflags2_300_(int *isc,int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].MS.Input.precipitationFlag[*i][0]=*sfcVar;
  swath300[*isc].MS.Input.precipitationFlag[*i][1]=*sfcVar;
}

//begin  WSO 8/20/14 write new ioquality flags
void copyioqualitys1_300_(int *isc,int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.FLG.ioQuality[*i]=*sfcVar;
}

void copyioqualitys2_300_(int *isc,int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].MS.FLG.ioQuality[*i]=*sfcVar;
}
//end    WSO 8/20/14
//
//begin  WSO 3/17/17 write snow ice cover flags
void copysnowices1_300_(int *isc,int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.Input.snowIceCover[*i]=*sfcVar;
}

void copysnowices2_300_(int *isc,int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].MS.Input.snowIceCover[*i]=*sfcVar;
}
//end    WSO 3/17/17

void copysfcliqfracts1_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.surfLiqRateFrac[*i]=*sfcVar;
}

void copysfcliqfracts2_300_(int *isc,float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].MS.surfLiqRateFrac[*i]=*sfcVar;
}

void copycldwaters1_300_(int *isc,float *var1d, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<nbins;k++)
    {
        swath300[*isc].NS.cloudLiqWaterCont[*i][k]=var1d[k];
    }
}

void copycldwaters2_300_(int *isc,float *var1d, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300]; 

  for(k=0;k<nbins;k++)
    {
        swath300[*isc].MS.cloudLiqWaterCont[*i][k]=var1d[k];
    }
}

void copycldices1_300_(int *isc,float *var1d, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<nbins;k++)
    {
      swath300[*isc].NS.cloudIceWaterCont[*i][k]=var1d[k];
    }
}

void copycldices2_300_(int *isc,float *var1d, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<nbins;k++)
    {
      swath300[*isc].MS.cloudIceWaterCont[*i][k]=var1d[k];
    }
}

//begin  WSO 9/5/13 new copy routine for SRT and DSRT pia effective sigma's
void copysigmapias1_300_(int *isc,float *sigmapia, int *i)
{
  extern L2BCMB_SWATHS swath300[300];
//diagnostic
//  printf("in writeCMBOut i: %5i,  sigmapia: %10.4f\n", *sigmapia, *i);
//end diagnostic
  swath300[*isc].NS.Input.piaEffectiveSigma[*i] = *sigmapia;
}
void copysigmapias2_300_(int *isc,float *sigmapiaku, float *sigmapiaka, int *i)
{
  extern L2BCMB_SWATHS swath300[300];
//    diagnostic
//  printf("in writeCMBOut i: %5i,  sigmapiaku: %10.4f,  sigmapiaka: %10.4f\n", 
//     *sigmapiaku, *sigmapiaka, *i);
//end diagnostic
    swath300[*isc].MS.Input.piaEffectiveSigma[*i][0] = *sigmapiaku;
    swath300[*isc].MS.Input.piaEffectiveSigma[*i][1] = *sigmapiaka;
}
//end    WSO 9/5/13

//write principal components
void copyprincomps1_300_(int *isc,float *princomp, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<5;k++)
    {
      swath300[*isc].NS.aPriori.prinComp[*i][k] = princomp[k];
    }
}
//
void copyprincomps2_300_(int *isc,float *princomp, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<5;k++)
    {
     swath300[*isc].MS.aPriori.prinComp[*i][k] = princomp[k];
    }
}

//write profile class
void copyprofclasss1_300_(int *isc,int *profclass, int *i)
{
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].NS.aPriori.profClass[*i] = *profclass;
}

void copyprofclasss2_300_(int *isc,int *profclass, int *i)
{
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].MS.aPriori.profClass[*i] = *profclass;
}

//write surface precip bias ratio
void copysurfprecipbiasratios1_300_(int *isc,float *biasratio, int *i)
{ 
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].NS.aPriori.surfPrecipBiasRatio[*i] = *biasratio;
}

void copysurfprecipbiasratios2_300_(int *isc,float *biasratio, int *i)
{
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].MS.aPriori.surfPrecipBiasRatio[*i] = *biasratio;
} 


//write initial log10 of the PSD intercept
void copyinitnws1_300_(int *isc,float *initlogNw, int *n9, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<9;k++)
    {
     if(n9[k]>1 && n9[k]<=88)
       {
         swath300[*isc].NS.aPriori.initNw[*i][k] = initlogNw[n9[k]-1];
       }
     else
       {
         swath300[*isc].NS.aPriori.initNw[*i][k] = missing_r4c;
       }
    }
}

void copyinitnws2_300_(int *isc,float *initlogNw, int *n9, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<9;k++)
    {
     if(n9[k]>1 && n9[k]<=88)
       {
         swath300[*isc].MS.aPriori.initNw[*i][k] = initlogNw[n9[k]-1];
       }
     else
       { 
         swath300[*isc].MS.aPriori.initNw[*i][k] = missing_r4c;
       }
    }
}

//write sub-footprint variability parameter
void copysubfootvariabilitys1_300_(int *isc,float *subfoot, int *i)
{ 
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].NS.nubfPIAfactor[*i] = *subfoot;
} 
  
void copysubfootvariabilitys2_300_(int *isc,float *subfoot, int *i)
{
  extern L2BCMB_SWATHS swath300[300];
    
  swath300[*isc].MS.nubfPIAfactor[*i] = *subfoot;
}

//write multiple scattering flag
void copymultiscatcalcs1_300_(int *isc,int *multiscat, int *i)
{ 
  extern L2BCMB_SWATHS swath300[300];
  
  swath300[*isc].NS.FLG.multiScatCalc[*i] = *multiscat;
}

void copymultiscatcalcs2_300_(int *isc,int *multiscat, int *i)
{ 
  extern L2BCMB_SWATHS swath300[300];
  
  swath300[*isc].MS.FLG.multiScatCalc[*i] = *multiscat;
}

//write multiple scattering surface parameter
void copymultiscatsurfaces1_300_(int *isc,float *multisfc, int *i)
{
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].NS.multiScatMaxContrib[*i] = *multisfc;
}

void copymultiscatsurfaces2_300_(int *isc,float *multisfc, int *i)
{
  extern L2BCMB_SWATHS swath300[300];

  swath300[*isc].MS.multiScatMaxContrib[*i] = *multisfc;
}

//
//begin  WSO 2/8/17 copy routine for measured sigma-zeros
void copysigmazeros1_300_(int *isc,float *sigmazeroku, int *i)
{
  extern L2BCMB_SWATHS swath300[300];
  swath300[*isc].NS.Input.sigmaZeroMeasured[*i] = *sigmazeroku;
}
void copysigmazeros2_300_(int *isc,float *sigmazeroku, float *sigmazeroka, int *i)
{
    extern L2BCMB_SWATHS swath300[300];
//        swath300[*isc].MS.Input.sigmaZeroMeasured[*i][0] = *sigmazeroku;
//            swath300[*isc].MS.Input.sigmaZeroMeasured[*i][1] = *sigmazeroka;
          swath300[*isc].MS.Input.sigmaZeroMeasured[*i] = *sigmazeroka;
}
//end    WSO 2/8/17

//begin  WSO 8/19/13 modified copy routines to include nodes
void copylognws1_300_(int *isc,float *logNw, int *n9, int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<9;k++)
    {
//begin WSO 9/7/14 added upper limit for n9 in next line
      if(n9[k]>1 && n9[k]<=88)
        {
         swath300[*isc].NS.PSDparamLowNode[*i][k] = n9[k]-1;
	 swath300[*isc].NS.precipTotPSDparamLow[*i][k][0]=logNw[n9[k]-1];
        }
      else
        {
         swath300[*isc].NS.PSDparamLowNode[*i][k] = missing_i2c;
         swath300[*isc].NS.precipTotPSDparamLow[*i][k][0]= missing_r4c;
        }
    }
}

void copylognws2_300_(int *isc,float *logNw, int *n9,int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];
  
  for(k=0;k<9;k++)
    {
      if(n9[k]>0)
        {
         swath300[*isc].MS.PSDparamLowNode[*i][k] = n9[k]-1;
	 swath300[*isc].MS.precipTotPSDparamLow[*i][k][0]=logNw[n9[k]-1];
        }
      else
        {
         swath300[*isc].MS.PSDparamLowNode[*i][k] = missing_i2c;
	 swath300[*isc].MS.precipTotPSDparamLow[*i][k][0]= missing_r4c;
        }
    }
}
//end    WSO 8/19/13

//begin  WSO 8/19/13 add mu as second low-resolution parameter
void copymus1_300_(int *isc,float *mu_prof, int *n9,int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<9;k++)
    {
      if(n9[k]>1)
        {
         swath300[*isc].NS.precipTotPSDparamLow[*i][k][1]=mu_prof[n9[k]-1];
        }
      else
        {
         swath300[*isc].NS.precipTotPSDparamLow[*i][k][1]=missing_r4c;
        }
    }
}

void copymus2_300_(int *isc,float *mu_prof, int *n9,int *i)
{
  int k;
  extern L2BCMB_SWATHS swath300[300];

  for(k=0;k<9;k++)
    {
      if(n9[k]>0)
        {
         swath300[*isc].MS.precipTotPSDparamLow[*i][k][1]=mu_prof[n9[k]-1];
        }
      else
        {
         swath300[*isc].MS.precipTotPSDparamLow[*i][k][1]=missing_r4c;
        }
    }
}
//end    WSO 8/19/13


void idiot_check_300_(int *isc,int *number, char *ident)
{
printf(" sfm idiot check %i %s \n",number,ident);
}

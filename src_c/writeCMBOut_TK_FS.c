//  SFM 04/06/2013  Code module added in merge from M.Grecu's code
//  SFM 05/06/2013  Modifications from LW to facilitate using job names
//  SFM 06/27/2013  Parameter name changes from W.Olson; reduce unused code
//  SFM 07/19/2013  Large volume of code added for M.Grecu
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TKheaders.h"
#include "TK_2BCMB.h"
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
L2ADPR_SWATHS dprxswath;
L2BCMB_SWATHS swath1;
L2AKu_FS     L2AKuDataX;

void openoutputfile_fs_(char *jobname, char *fname)
{
  char ranstring[80], ranMsg[80];
  int status;
  printf(" Output file = %s \n",fname);

  int ret;
//  SFM  04/06/2013  Changed file type to 2BCMB   WSO 04/07/2013
//  SFM  09/04/2013  Changed jobid to jobname
//  SFM  09/11/2013  Moved metadta settings to centralized location

  ret = TKopen(fname, "2BCMB", TKWRITE, "HDF5", jobname, &ctkfile,1); //WSO 04/07/2013
}


void closeoutputfile_fs_(void)
{
  int ret;
  ret=TKclose(&ctkfile);
}



void setlatlons1_fs_(float *lat, float *lon, float *sfcPrecip, 
                  float *sfcPrecipStd, float *piaOut)
{
  int i;
  extern L2BCMB_SWATHS swathx;

  for(i=0;i<49;i++)
    {
      swathx.KuGMI.Latitude[i]=lat[i];
      swathx.KuGMI.Longitude[i]=lon[i];
      if(swathx.KuGMI.Longitude[i]>180)
	swathx.KuGMI.Longitude[i]-=360;
      swathx.KuGMI.nearSurfPrecipTotRate[i]=sfcPrecip[i];
      swathx.KuGMI.nearSurfPrecipTotRateSigma[i]=sfcPrecipStd[i];
      swathx.KuGMI.pia[i]=piaOut[i];
    }
}

void setlatlons2_fs_(float *lat, float *lon, float *sfcPrecip, 
		     float *sfcPrecipStd, float *piaOutKu, float *piaOutKa,
		     int *scanPatternFlag)
{
  int i;
  extern L2BCMB_SWATHS swathx;

  for(i=0;i<49;i++)
    {
      swathx.KuKaGMI.Latitude[i]=lat[i];
      swathx.KuKaGMI.Longitude[i]=lon[i];
      if(swathx.KuKaGMI.Longitude[i]>180)
	swathx.KuKaGMI.Longitude[i]-=360;
      swathx.KuKaGMI.nearSurfPrecipTotRate[i]=sfcPrecip[i];
      swathx.KuKaGMI.nearSurfPrecipTotRateSigma[i]=sfcPrecipStd[i];
      swathx.KuKaGMI.pia[i][0]=piaOutKu[i];
      swathx.KuKaGMI.pia[i][1]=piaOutKa[i];
      if( *scanPatternFlag==0 && (i<12 || i>=37))
	{
	  swathx.KuKaGMI.Latitude[i]=missing_r4c;
	  swathx.KuKaGMI.Longitude[i]=missing_r4c;
	  swathx.KuKaGMI.nearSurfPrecipTotRate[i]=missing_r4c;
	  swathx.KuKaGMI.nearSurfPrecipTotRateSigma[i]=missing_r4c;
	  swathx.KuKaGMI.pia[i][0]=missing_r4c;
	  swathx.KuKaGMI.pia[i][1]=missing_r4c;
	}
    }
}

void copy_oe_missing_staff_(int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEestimSurfPrecipLiqRate[*i]=missing_r4c;
swathx.KuGMI.OptEst.OEestimSurfPrecipTotRate[*i]=missing_r4c;
swathx.KuKaGMI.OptEst.OEestimSurfPrecipLiqRate[*i]=missing_r4c;
swathx.KuKaGMI.OptEst.OEestimSurfPrecipTotRate[*i]=missing_r4c;
swathx.KuKaGMI.OptEst.OEestimSurfPrecipTotRateSigma[*i]=missing_r4c;
swathx.KuGMI.OptEst.OEestimSurfPrecipTotRateSigma[*i]=missing_r4c;
}


void copy_oecloudicepath_kugmi_(float *OECloudLiqPath, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEcolumnCloudIceWater[*i]=missing_r4c;
}

void copy_oecloudicepath_kukagmi_(float *OECloudLiqPath, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEcolumnCloudIceWater[*i]=missing_r4c;
}

void copyrrates1_fs_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
      swathx.KuGMI.precipTotRate[*i][k]=rrate[k];
      swathx.KuGMI.precipTotRateSigma[*i][k]=rratestd[k];
    }
}

void copysflfract_fs_(float *lfract, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  swathx.KuGMI.nearSurfPrecipLiqRate[*i]=*lfract;
  if(*i>=12 && *i<=37)
    swathx.KuKaGMI.nearSurfPrecipLiqRate[*i]=*lfract;
  
}


//begin  WSO 8/30/13
void copyenvsfqvs1_fs_(float *envQv, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.surfaceVaporDensity[*i]=envQv[nbins-1];
}

void copyenvsfqvs2_fs_(float *envQv, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuKaGMI.surfaceVaporDensity[*i]=envQv[nbins-1];
}

void copyenvqvs1_fs_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  for(k=0;k<10;k++)
      swathx.KuGMI.vaporDensity[*i][k]=envQv[envnodes[k]-1];
}

void copyenvqvs2_fs_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  for(k=0;k<10;k++)
      swathx.KuKaGMI.vaporDensity[*i][k]=envQv[envnodes[k]-1];
}

void copyenvpresss1_fs_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  for(k=0;k<10;k++)
    swathx.KuGMI.airPressure[*i][k]=envQv[envnodes[k]-1];
}

void copyenvpresss2_fs_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  for(k=0;k<10;k++)
    swathx.KuKaGMI.airPressure[*i][k]=envQv[envnodes[k]-1];
}

void copyenvtemps1_fs_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  for(k=0;k<10;k++)
     {
      swathx.KuGMI.envParamNode[*i][k]=envnodes[k]-1;
      swathx.KuGMI.airTemperature[*i][k]=envQv[envnodes[k]-1];
     }
}

void copyenvtemps2_fs_(float *envQv, short *envnodes, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  for(k=0;k<10;k++)
    {
     swathx.KuKaGMI.envParamNode[*i][k]=envnodes[k]-1;
     swathx.KuKaGMI.airTemperature[*i][k]=envQv[envnodes[k]-1];
    }
}

void copyenvsftemps1_fs_(float *envQv, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.surfaceAirTemperature[*i]=envQv[nbins-1];
}

void copyenvsftemps2_fs_(float *envQv, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuKaGMI.surfaceAirTemperature[*i]=envQv[nbins-1];
}

//end    WSO 8/30/13

void copypwcs1_fs_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
      swathx.KuGMI.precipTotWaterCont[*i][k]=rrate[k];
      swathx.KuGMI.precipTotWaterContSigma[*i][k]=rratestd[k];
    }
}

//begin  WSO 8/7/13
void copylwcfracs1_fs_(float *mlwc_frac, float *mrate_frac, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<ntransitions;k++)
    {
      //swathx.KuGMI.liqMassFracTrans[*i][k]=mlwc_frac[k];
      //swathx.KuGMI.liqRateFracTrans[*i][k]=mrate_frac[k];
    }
}

void copysfcrainliqfracs1_fs_(float *sfcrainliq_frac, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.nearSurfPrecipLiqRate[*i]=*sfcrainliq_frac;

}
//end    WSO 8/7/13

void copyd0s1_fs_(float *dm, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
      swathx.KuGMI.precipTotDm[*i][k]=dm[k];
    }
}

void copyd0s1_a_fs_(float *dm, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
      swathx.KuGMI.aPriori.initPrecipTotDm[*i][k]=dm[k];
    }
}

void copyd0s2_a_fs_(float *dm, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
      swathx.KuKaGMI.aPriori.initPrecipTotDm[*i][k]=dm[k];
    }
}


void copynws1_fs_(float *Nw, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
      swathx.KuGMI.precipTotLogNw[*i][k]=Nw[k];
      swathx.KuGMI.precipTotMu[*i][k]=2.0;
    }
}

void copynws1_a_fs_(float *Nw, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
      swathx.KuGMI.aPriori.initPrecipTotLogNw[*i][k]=Nw[k];
    }
}

void copynws2_fs_(float *Nw, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
      swathx.KuKaGMI.precipTotLogNw[*i][k]=Nw[k];
      swathx.KuKaGMI.precipTotMu[*i][k]=2.0;
    }
}


void copynws2_a_fs_(float *Nw, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
       swathx.KuKaGMI.aPriori.initPrecipTotLogNw[*i][k]=Nw[k];
    }
}

void copyzckus1_fs_(float *zc, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

   for(k=0;k<nbins;k++)
    {
      if(zc[k] > -90.)
        swathx.KuGMI.correctedReflectFactor[*i][k] = zc[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swathx.KuGMI.correctedReflectFactor[*i][k] = missing_r4c;
//end    WSO 9/17/13
    }
}

void copynodess1_fs_(int *node, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  for(k=0;k<5;k++)
    {
      swathx.KuGMI.phaseBinNodes[*i][k]=node[k];
    }
}

void copyrrates2_fs_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
      swathx.KuKaGMI.precipTotRate[*i][k]=rrate[k];
      swathx.KuKaGMI.precipTotRateSigma[*i][k]=rratestd[k];
    }
}

void copypwcs2_fs_(float *rrate, float *rratestd, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
//begin WSO 4/18/2013
//changed KuGMI to MS
      swathx.KuKaGMI.precipTotWaterCont[*i][k]=rrate[k];
      swathx.KuKaGMI.precipTotWaterContSigma[*i][k]=rratestd[k];
//end  WSO 4/18/2013
    }
}

//begin  WSO 8/7/13
void copylwcfracs2_fs_(float *mlwc_frac, float *mrate_frac, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<ntransitions;k++)
    {
      //      swathx.KuKaGMI.liqMassFracTrans[*i][k]=mlwc_frac[k];
      //      swathx.KuKaGMI.liqRateFracTrans[*i][k]=mrate_frac[k];
    }
}

void copy_tot_to_liqwatercont_ku_(int *ij, float *precipTot, int *node)
{
  int i;
  extern L2BCMB_SWATHS swathx;
  for(i=node[0];i<=node[1];i++)
    swathx.KuGMI.precipLiqWaterCont[*ij][i]=0.0;
  int n3=node[3];
  if(n3>node[4])
    n3=node[4];
  for(i=node[1]+1;i<=n3;i++)
    if(precipTot[i]>=-1e-9)
      {
	float f=(i-node[1])/(node[3]-node[1]+1e-3);
	swathx.KuGMI.precipLiqWaterCont[*ij][i]=f*precipTot[i];
      }
    else
      swathx.KuGMI.precipLiqWaterCont[*ij][i]=missing_r4c;
  int node4=node[4]+1;
  if(node4>87)
    node4=87;
  for(i=node[3]+1;i<=node4;i++)
    if(precipTot[i]>=-1e-9)
      swathx.KuGMI.precipLiqWaterCont[*ij][i]=precipTot[i];
    else
      swathx.KuGMI.precipLiqWaterCont[*ij][i]=missing_r4c;
  for(i=node[4]+1;i<88;i++)
    swathx.KuGMI.precipLiqWaterCont[*ij][i]=missing_r4c;
}

void copy_tot_to_liqrate_ku_(int *ij, float *precipTot, int *node)
{
  int i;
  extern L2BCMB_SWATHS swathx;
  for(i=node[0];i<=node[1];i++)
    swathx.KuGMI.precipLiqRate[*ij][i]=0.0;
  int n3=node[3];
  if(n3>node[4])
    n3=node[4];
  for(i=node[1]+1;i<=n3;i++)
    if(precipTot[i]>=-1e-9)
      {
	float f=(i-node[1])/(node[3]-node[1]+1e-3);
	swathx.KuGMI.precipLiqRate[*ij][i]=f*precipTot[i];
      }
    else
      swathx.KuGMI.precipLiqRate[*ij][i]=missing_r4c;
  int node4=node[4]+1;
  if(node4>87)
    node4=87;
  for(i=node[3]+1;i<=node4;i++)
    if(precipTot[i]>=-1e-9)
      	swathx.KuGMI.precipLiqRate[*ij][i]=precipTot[i];
    else
      swathx.KuGMI.precipLiqRate[*ij][i]=missing_r4c;
  for(i=node[4]+1;i<88;i++)
    swathx.KuGMI.precipLiqRate[*ij][i]=missing_r4c;
}
void copy_tot_to_liqrate_kuka_(int *ij, float *precipTot, int *node, int *scanPatternFlag)
{
 int i;
 extern L2BCMB_SWATHS swathx;

 if(*scanPatternFlag==0 &&(*ij<12 || *ij>37))
   {
     for(i=0;i<88;i++)
       swathx.KuKaGMI.precipLiqRate[*ij][i]=missing_r4c;
     return;
   }
 for(i=node[0];i<=node[1];i++)
   swathx.KuKaGMI.precipLiqRate[*ij][i]=0.0;
 int n3=node[3];
 if(n3>node[4])
   n3=node[4];
 for(i=node[1]+1;i<=n3;i++)
   if(precipTot[i]>=-1e-9)
     {
       float f=(i-node[1])/(node[3]-node[1]+1e-3);
       swathx.KuKaGMI.precipLiqRate[*ij][i]=f*precipTot[i];
     }
   else
     swathx.KuKaGMI.precipLiqRate[*ij][i]=missing_r4c;
 int node4=node[4]+1;
 if(node4>87)
   node4=87;
 for(i=node[3]+1;i<=node4;i++)
   if(precipTot[i]>=-1e-9)
     swathx.KuKaGMI.precipLiqRate[*ij][i]=precipTot[i];
   else
     swathx.KuKaGMI.precipLiqRate[*ij][i]=missing_r4c;
 for(i=node[4]+1;i<88;i++)
   swathx.KuKaGMI.precipLiqRate[*ij][i]=missing_r4c;
}
			  //'precipLiqWaterCont'

void copy_tot_to_liqwatercont_kuka_(int *ij, float *precipTot, int *node, int *scanPatternFlag)
{
 int i;
 extern L2BCMB_SWATHS swathx;

 if(*scanPatternFlag==0 &&(*ij<12 || *ij>37))
   {
     for(i=0;i<88;i++)
       swathx.KuKaGMI.precipLiqWaterCont[*ij][i]=missing_r4c;
     return;
   }
 for(i=node[0];i<=node[1];i++)
   swathx.KuKaGMI.precipLiqWaterCont[*ij][i]=0.0;
 int n3=node[3];
 if(n3>node[4])
   n3=node[4];
 for(i=node[1]+1;i<=n3;i++)
   if(precipTot[i]>=-1e-9)
     {
       float f=(i-node[1])/(node[3]-node[1]+1e-3);
       swathx.KuKaGMI.precipLiqWaterCont[*ij][i]=f*precipTot[i];
     }
   else
     swathx.KuKaGMI.precipLiqWaterCont[*ij][i]=missing_r4c;
 int node4=node[4]+1;
 if(node4>87)
   node4=87;
 for(i=node[3]+1;i<=node4;i++)
   if(precipTot[i]>=-1e-9)
     swathx.KuKaGMI.precipLiqWaterCont[*ij][i]=precipTot[i];
   else
     swathx.KuKaGMI.precipLiqWaterCont[*ij][i]=missing_r4c;
 for(i=node[4]+1;i<88;i++)
   swathx.KuKaGMI.precipLiqWaterCont[*ij][i]=missing_r4c;
}
			  //'precipLiqWaterCont'

void copysfcrainliqfracs2_fs_(float *sfcrainliq_frac, int *i)
{   
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuKaGMI.nearSurfPrecipLiqRate[*i]=*sfcrainliq_frac;

}
//end   WSO 8/7/13


void copyd0s2_fs_(float *dm, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
// SFM 05/06/2013 Changed KuGMI to MS to match M.Grecu code from 04/19/2013
      swathx.KuKaGMI.precipTotDm[*i][k]=dm[k];
    }
}

void copyzckus2_fs_(float *zku, float *zka, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

   for(k=0;k<nbins;k++)
    {
      if(zku[k] > -90.)
        swathx.KuKaGMI.correctedReflectFactor[*i][k][0] = zku[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swathx.KuKaGMI.correctedReflectFactor[*i][k][0] = missing_r4c;
//end    WSO 9/17/13
      if(zka[k] > -90.)
        swathx.KuKaGMI.correctedReflectFactor[*i][k][1] = zka[k];
      else
//begin  WSO 9/17/13 standardized missing flags
        swathx.KuKaGMI.correctedReflectFactor[*i][k][1] = missing_r4c;
//end    WSO 9/17/13
    }
}
void copynodess2_fs_(int *node, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  for(k=0;k<5;k++)
    {
      swathx.KuKaGMI.phaseBinNodes[*i][k]=node[k];
    }
}

void rewind_fs_(int *ic)
{
  extern TKINFO       granuleHandle2AKu;
  int status = TKseek(&granuleHandle2AKu, *ic, TK_ABS_SCAN_OFF); 
}

//begin WSO 9/8/13 rewind DPR file
void rewind_dpr_fs_(int *ic)
{
  extern TKINFO       dprtkfile;
  int status_dpr = TKseek(&dprtkfile, *ic, TK_ABS_SCAN_OFF);
}
//end WSO 9/8/13

//  SFM  begin  12/13/2013; add flag to call sequence
void frominput_fs_(long *st_2adpr, int *flagScanPattern)
{
//  SFM  begin  12/13/2013
  extern TKINFO       granuleHandle2AKu;
  extern L2AKu_FS        L2AKuDataX;
  extern L2BCMB_SWATHS swathx;
//begin  WSO 9/1/13
  extern L2ADPR_SWATHS dprswath;
  extern L2ADPR_SWATHS dprxswath;
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
  
  int c_flagScanPattern=(int) dprxswath.FS.FLG.flagScanPattern[1];
  //printf("read flagScanPattern %i \n",c_flagScanPattern);
  *flagScanPattern= c_flagScanPattern;
  for( j=0; j<49; j++)
    {
      //swathx.KuGMI.Input.piaEffective[j]=L2AKuData.SRT.pathAtten[j];
      swathx.KuGMI.Input.piaEffective[j]=L2AKuDataX.SRT.PIAhybrid[j];  //MG  7/31/18, use hybrid PIA
//begin  WSO 9/5/13 remove flag assignment
//       swathx.KuGMI.Input.piaEffectiveSigma[j]=-99;
//end    WSO 9/5/13
   //   swathx.KuGMI.Input.piaEffectiveReliabFlag[j]=
	// L2AKuData.SRT.reliabFlag[j];
      swathx.KuGMI.Input.piaEffectiveReliabFlag[j]=
	L2AKuDataX.SRT.reliabFlagHY[j];                              //WSO  8/2/18 use hybrid flag
      swathx.KuGMI.Input.precipitationType[j]=
	L2AKuDataX.CSF.typePrecip[j];
      swathx.KuGMI.Input.precipTypeQualityFlag[j]=
	L2AKuDataX.CSF.qualityTypePrecip[j];
      swathx.KuGMI.Input.surfaceElevation[j]=L2AKuDataX.PRE.elevation[j];
      swathx.KuGMI.Input.localZenithAngle[j]=L2AKuDataX.PRE.localZenithAngle[j];
      swathx.KuGMI.Input.surfaceType[j]=L2AKuDataX.PRE.landSurfaceType[j];
//begin  WSO 9/28/13 use alternate rain flag that includes missing for bad scans
//      swathx.KuGMI.Input.precipitationFlag[j]=L2AKuData.PRE.flagPrecip[j];
//end    WSO 9/28/13
      swathx.KuGMI.Input.surfaceRangeBin[j]=(L2AKuDataX.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
      swathx.KuGMI.Input.stormTopBin[j]=(L2AKuDataX.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
      if(swathx.KuGMI.Input.stormTopBin[j]<0)
	swathx.KuGMI.Input.stormTopBin[j]=missing_i2c;
      swathx.KuGMI.Input.stormTopAltitude[j]=L2AKuDataX.PRE.heightStormTop[j];
//begin  WSO 09/30/15 add one bin to the binClutterFreeBottom to temporarily compensate for
//the subtraction of one bin by the radar team
//      swathx.KuGMI.Input.lowestClutterFreeBin[j]=
//	(L2AKuData.PRE.binClutterFreeBottom[j]-1)/2; // MG 04/11/2014
//begin  WSO 10/19/15 subtract one 125 m bin from binClutterFreeBottom, and 
//restore V3 definition of lowestClutterFreeBin as a test
//      swathx.KuGMI.Input.lowestClutterFreeBin[j]=
//	(L2AKuData.PRE.binClutterFreeBottom[j])/2;
      swathx.KuGMI.Input.lowestClutterFreeBin[j]=
	(L2AKuDataX.PRE.binClutterFreeBottom[j] - 2)/2;
//end    WSO 10/15/15
//end    WSO 09/30/15
//begin  WSO 9/17/13 correction for two bin average location in combined
      swathx.KuGMI.Input.ellipsoidBinOffset[j]=
	    L2AKuDataX.PRE.ellipsoidBinOffset[j] + 0.125/2.;
//end    WSO 9/17/13
//begin  WSO 8/19/13
      swathx.KuGMI.Input.zeroDegAltitude[j] = L2AKuDataX.VER.heightZeroDeg[j];
      swathx.KuGMI.Input.zeroDegBin[j] = (L2AKuDataX.VER.binZeroDeg[j]-1)/2; // MG 04/11/2014
//end    WSO 8/19/13
      if((j>=0 && j<49 &&  *flagScanPattern==1) || (j>=12 && j<37))
	{
	  swathx.KuKaGMI.Input.surfaceElevation[j]=
	    L2AKuDataX.PRE.elevation[j];
	  swathx.KuKaGMI.Input.localZenithAngle[j]=
	    L2AKuDataX.PRE.localZenithAngle[j];
	  swathx.KuKaGMI.Input.surfaceType[j]=
	    L2AKuDataX.PRE.landSurfaceType[j];
	  swathx.KuKaGMI.Input.surfaceRangeBin[j][0]=
	    (L2AKuDataX.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
	  swathx.KuKaGMI.Input.surfaceRangeBin[j][1]=
	    (L2AKuDataX.PRE.binRealSurface[j]-1)/2; // MG 04/11/2014
	  swathx.KuKaGMI.Input.stormTopBin[j][0]=
	    (L2AKuDataX.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
	  swathx.KuKaGMI.Input.stormTopBin[j][1]=  // MG 04/11/2014
	    (L2AKuDataX.PRE.binStormTop[j]-1)/2; // MG 04/11/2014
	  if(swathx.KuKaGMI.Input.stormTopBin[j][0]<0)
	    swathx.KuKaGMI.Input.stormTopBin[j][0]=missing_i2c;
	  if(swathx.KuKaGMI.Input.stormTopBin[j][1]<0)
	    swathx.KuKaGMI.Input.stormTopBin[j][1]=missing_i2c;
	  swathx.KuKaGMI.Input.stormTopAltitude[j][0]=
	    L2AKuDataX.PRE.heightStormTop[j];
	  swathx.KuKaGMI.Input.stormTopAltitude[j][1]=
	    L2AKuDataX.PRE.heightStormTop[j];
	  swathx.KuKaGMI.Input.lowestClutterFreeBin[j][0]=
	    (L2AKuDataX.PRE.binClutterFreeBottom[j] - 2)/2;
	  swathx.KuKaGMI.Input.lowestClutterFreeBin[j][1]=
	    (L2AKuDataX.PRE.binClutterFreeBottom[j] - 2)/2;
//end    WSO 10/19/15
//end    WSO 09/30/15
//begin  WSO 9/17/13 correction for two bin average location in combined
	  swathx.KuKaGMI.Input.ellipsoidBinOffset[j][0]=
	    L2AKuDataX.PRE.ellipsoidBinOffset[j] + 0.125/2.;
	  swathx.KuKaGMI.Input.ellipsoidBinOffset[j][1]=
	    L2AKuDataX.PRE.ellipsoidBinOffset[j] + 0.125/2.;
//end    WSO 9/17/13
//begin  WSO 9/5/13 reset pia's using DPR output
	  swathx.KuKaGMI.Input.piaEffective[j][0]=  
	    dprxswath.FS.SRT.pathAtten[j][0];
	  swathx.KuKaGMI.Input.piaEffective[j][1]=
	    dprxswath.FS.SRT.pathAtten[j][1];
//begin  WSO 9/5/13 remove flag assignment
//	  swathx.KuKaGMI.Input.piaEffectiveSigma[j][0]=-99;
//end    WSO 9/5/13
	  swathx.KuKaGMI.Input.piaEffectiveReliabFlag[j][0]=
	    dprxswath.FS.SRT.reliabFlag[j];
	  swathx.KuKaGMI.Input.piaEffectiveReliabFlag[j][1]=
	    dprxswath.FS.SRT.reliabFlag[j];
//end    WSO 9/5/13
	  swathx.KuKaGMI.Input.precipitationType[j]=
	    L2AKuDataX.CSF.typePrecip[j];
	  swathx.KuKaGMI.Input.precipTypeQualityFlag[j]=
	    L2AKuDataX.CSF.qualityTypePrecip[j];
//begin  WSO 8/19/13 need to update toolkit
          swathx.KuKaGMI.Input.zeroDegAltitude[j] =
            L2AKuDataX.VER.heightZeroDeg[j];
          swathx.KuKaGMI.Input.zeroDegBin[j][0] =
            (L2AKuDataX.VER.binZeroDeg[j]-1)/2;   // MG 04/11/2014]
	  swathx.KuKaGMI.Input.zeroDegBin[j][1] =
            (L2AKuDataX.VER.binZeroDeg[j]-1)/2;   // MG 04/11/2014
//end    WSO 8/19/13
	}
      else
	{
	  swathx.KuKaGMI.Input.surfaceElevation[j]=missing_r4c;
	  swathx.KuKaGMI.Input.localZenithAngle[j]=missing_r4c;
	  swathx.KuKaGMI.Input.surfaceType[j]=missing_i2c;
	  swathx.KuKaGMI.Input.surfaceRangeBin[j][0]=missing_i2c;
	  swathx.KuKaGMI.Input.surfaceRangeBin[j][1]=missing_i2c;
	  swathx.KuKaGMI.Input.stormTopBin[j][0]=missing_i2c;
	  swathx.KuKaGMI.Input.stormTopBin[j][1]=missing_i2c;
	  swathx.KuKaGMI.Input.stormTopAltitude[j][0]=missing_r4c;
	  swathx.KuKaGMI.Input.stormTopAltitude[j][1]=missing_r4c;
	  swathx.KuKaGMI.Input.lowestClutterFreeBin[j][0]=missing_i2c;
	  swathx.KuKaGMI.Input.lowestClutterFreeBin[j][1]=missing_i2c;
	  swathx.KuKaGMI.Input.ellipsoidBinOffset[j][0]=missing_i2c;
	  swathx.KuKaGMI.Input.ellipsoidBinOffset[j][1]=missing_i2c;
	  swathx.KuKaGMI.Input.piaEffective[j][0]=missing_r4c; 
	  swathx.KuKaGMI.Input.piaEffective[j][1]=missing_r4c; 
	  swathx.KuKaGMI.Input.piaEffectiveReliabFlag[j][0]=missing_i2c;
	  swathx.KuKaGMI.Input.piaEffectiveReliabFlag[j][1]=missing_i2c;
	  swathx.KuKaGMI.Input.precipitationType[j]=missing_i4c;
	  swathx.KuKaGMI.Input.precipTypeQualityFlag[j]=missing_i2c;
          swathx.KuKaGMI.Input.zeroDegAltitude[j] =missing_r4c;
          swathx.KuKaGMI.Input.zeroDegBin[j][0] =missing_i2c;
          swathx.KuKaGMI.Input.zeroDegBin[j][1] =missing_i2c;
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
    swathx.KuGMI.scanStatus.FractionalGranuleNumber =    
     L2AKuDataX.scanStatus.FractionalGranuleNumber;
    swathx.KuGMI.scanStatus.SCorientation =
     L2AKuDataX.scanStatus.SCorientation;
    swathx.KuGMI.scanStatus.acsModeMidScan =
     L2AKuDataX.scanStatus.acsModeMidScan;
    swathx.KuGMI.scanStatus.dataQuality =
     L2AKuDataX.scanStatus.dataQuality;
    swathx.KuGMI.scanStatus.dataWarning =
     L2AKuDataX.scanStatus.dataWarning;
    swathx.KuGMI.scanStatus.geoError =
     L2AKuDataX.scanStatus.geoError;
    swathx.KuGMI.scanStatus.geoWarning =
     L2AKuDataX.scanStatus.geoWarning;
    swathx.KuGMI.scanStatus.limitErrorFlag =
     L2AKuDataX.scanStatus.limitErrorFlag;
    swathx.KuGMI.scanStatus.missing =
     L2AKuDataX.scanStatus.missing;
    swathx.KuGMI.scanStatus.modeStatus =
     L2AKuDataX.scanStatus.modeStatus;
    swathx.KuGMI.scanStatus.operationalMode =
     L2AKuDataX.scanStatus.operationalMode;
    swathx.KuGMI.scanStatus.pointingStatus =
     L2AKuDataX.scanStatus.pointingStatus;
    swathx.KuGMI.scanStatus.targetSelectionMidScan =
     L2AKuDataX.scanStatus.targetSelectionMidScan;
//from 2ADPRX
    swathx.KuKaGMI.scanStatus.FractionalGranuleNumber =
     dprxswath.FS.scanStatus.FractionalGranuleNumber;
    swathx.KuKaGMI.scanStatus.SCorientation =
     dprxswath.FS.scanStatus.SCorientation;
    swathx.KuKaGMI.scanStatus.acsModeMidScan =
     dprxswath.FS.scanStatus.acsModeMidScan;
    swathx.KuKaGMI.scanStatus.dataQuality =
     dprxswath.FS.scanStatus.dataQuality[1];
    swathx.KuKaGMI.scanStatus.dataWarning =
     dprxswath.FS.scanStatus.dataWarning[1];
    swathx.KuKaGMI.scanStatus.geoError =
     dprxswath.FS.scanStatus.geoError[1];
    swathx.KuKaGMI.scanStatus.geoWarning =
     dprxswath.FS.scanStatus.geoWarning[1];
    swathx.KuKaGMI.scanStatus.limitErrorFlag =
     dprxswath.FS.scanStatus.limitErrorFlag[1];
    swathx.KuKaGMI.scanStatus.missing =
     dprxswath.FS.scanStatus.missing[1];
    swathx.KuKaGMI.scanStatus.modeStatus =
     dprxswath.FS.scanStatus.modeStatus[1];
    swathx.KuKaGMI.scanStatus.operationalMode =
     dprxswath.FS.scanStatus.operationalMode[1];
    swathx.KuKaGMI.scanStatus.pointingStatus =
     dprxswath.FS.scanStatus.pointingStatus[1];
    swathx.KuKaGMI.scanStatus.targetSelectionMidScan =
     dprxswath.FS.scanStatus.targetSelectionMidScan;
//end    WSO 9/1/13

}

void copyscantime_fs_(int *i)
{
  extern L2BCMB_SWATHS swathx;
  extern int DayOfMonth[300], DayOfYear[300], Hour[300], MilliSecond[300],
    Minute[300], Month[300], Second[300], Year[300], SecondOfDay[300];
  extern NAVIGATION navigation[300];

 swathx.KuGMI.ScanTime.DayOfMonth=DayOfMonth[*i];
 swathx.KuGMI.ScanTime.DayOfYear=DayOfYear[*i];
 swathx.KuGMI.ScanTime.Hour=Hour[*i];
 swathx.KuGMI.ScanTime.MilliSecond=MilliSecond[*i];
 swathx.KuGMI.ScanTime.Minute=Minute[*i];
 swathx.KuGMI.ScanTime.Month=Month[*i];
 swathx.KuGMI.ScanTime.Second=Second[*i];
 swathx.KuGMI.ScanTime.Year=Year[*i];
 swathx.KuGMI.ScanTime.SecondOfDay=SecondOfDay[*i];
 memcpy(&swathx.KuGMI.navigation, &navigation[*i], sizeof(NAVIGATION));

//begin WSO 04/07/2013
//added MS swath scantimes
 swathx.KuKaGMI.ScanTime.DayOfMonth=DayOfMonth[*i];
 swathx.KuKaGMI.ScanTime.DayOfYear=DayOfYear[*i];
 swathx.KuKaGMI.ScanTime.Hour=Hour[*i];
 swathx.KuKaGMI.ScanTime.MilliSecond=MilliSecond[*i];
 swathx.KuKaGMI.ScanTime.Minute=Minute[*i];
 swathx.KuKaGMI.ScanTime.Month=Month[*i];
 swathx.KuKaGMI.ScanTime.Second=Second[*i];
 swathx.KuKaGMI.ScanTime.Year=Year[*i];
 swathx.KuKaGMI.ScanTime.SecondOfDay=SecondOfDay[*i];
 memcpy(&swathx.KuKaGMI.navigation, &navigation[*i], sizeof(NAVIGATION));
//end WSO 04/07/2013
}

void copypreciptype_fs_(int *ptype, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  //swath.S1.precipitationType[*i]=*ptype;
}

void copyw10_fs_(float *w10, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.tenMeterWindSpeed[*i]=*w10;
}

void copyw10sigma_fs_(float *w10s, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.tenMeterWindSigma[*i]=*w10s;
}

void copyw10small_fs_(float *w10, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuKaGMI.tenMeterWindSpeed[*i]=*w10;
}

void copyw10smallsigma_fs_(float *w10s, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuKaGMI.tenMeterWindSigma[*i]=*w10s;
}



//  begin  SFM  12/26/2013
void write_empty_fs_(void)

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
void writescan_fs_(void)
{
  int ret;
  char emptygranuletext[100];
  extern L2BCMB_SWATHS swathx;
  // TKgetMetaString(&ctkfile, "FileHeader", "EmptyGranule", 
  //                emptygranuletext);	  
  // if (strncmp(emptygranuletext,"NOT_EMPTY",9) == 0)
  ret= TKwriteScan(&ctkfile,&swathx);
}
//  end    SFM  11/27/2013

void copysfcairtemps1_fs_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.surfaceAirTemperature[*i]=*sfcVar;
}

void copysfcairtemps2_fs_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuKaGMI.surfaceAirTemperature[*i]=*sfcVar;
}

void copysfcairpresss1_fs_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.surfaceAirPressure[*i]=*sfcVar;
}

void copysfcairpresss2_fs_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuKaGMI.surfaceAirPressure[*i]=*sfcVar;
}

void copyskintemps1_fs_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.skinTemperature[*i]=*sfcVar;
}

void copyskintemps2_fs_(float *sfcVar, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuKaGMI.skinTemperature[*i]=*sfcVar;
}

//write skin temperature estimate uncertainty
void copyskintempsigmas1_fs_(float *skinsigma, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  swathx.KuGMI.skinTempSigma[*i] = *skinsigma;
}

void copyskintempsigmas2_fs_(float *skinsigma, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swathx;

  swathx.KuKaGMI.skinTempSigma[*i] = *skinsigma;
}

//write column vapor estimate uncertainty
void copycolumnvaporsigmas1_fs_(float *colvaporsigma, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  //  swathx.KuGMI.columnVaporSigma[*i] = *colvaporsigma;
}

void copycolumnvaporsigmas2_fs_(float *colvaporsigma, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swathx;

  //  swathx.KuKaGMI.columnVaporSigma[*i] = *colvaporsigma;
}

//write column cloud liquid estimate uncerainty
void copycolumncloudliqsigmas1_fs_(float *colcldsigma, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  
  //  swathx.KuGMI.columnCloudLiqSigma[*i] = *colcldsigma;
}

void copycolumncloudliqsigmas2_fs_(float *colcldsigma, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swathx;

  //  swathx.KuKaGMI.columnCloudLiqSigma[*i] = *colcldsigma;
}

//write algorithm type flag
void copyalgotypes1_fs_(int *algotype, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  //  swathx.KuGMI.FLG.algoType[*i] = *algotype;
}

void copyalgotypes2_fs_(int *algotype, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  //  swathx.KuKaGMI.FLG.algoType[*i] = *algotype;
}


//write error of non-raining data fit
void copyerrorofdatafits1_fs_(float *erroroffit, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  
  //  swathx.KuGMI.errorOfDataFit[*i] = *erroroffit;
}

void copyerrorofdatafits2_fs_(float *erroroffit, int *i)
{  
  int k;
  extern L2BCMB_SWATHS swathx;

  //  swathx.KuKaGMI.errorOfDataFit[*i] = *erroroffit;
}

void copysfcemissouts1_fs_(float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swathx.KuGMI.surfEmissivity[*i][k]=tbout[k];
    else
      swathx.KuGMI.surfEmissivity[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swathx.KuGMI.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts1sigma_fs_(float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
    {
//      printf("%g ",tbout[k]);
      if(tbout[k] > 0.)
	swathx.KuGMI.surfEmissSigma[*i][k]=tbout[k];
      else
	swathx.KuGMI.surfEmissSigma[*i][k]=missing_r4c;
    }
//  printf(" from c\n");
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swathx.KuGMI.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts2_fs_(float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swathx.KuKaGMI.surfEmissivity[*i][k]=tbout[k];
    else
      swathx.KuKaGMI.surfEmissivity[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//    for(k=13;k<15;k++)
//      swathx.KuKaGMI.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copysfcemissouts2sigma_fs_(float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
//begin  WSO 10/13/15 change loops to output HF
//                    emissivities
  for(k=0;k<13;k++)
//begin  WSO 9/16/13
    if(tbout[k] > 0.)
      swathx.KuKaGMI.surfEmissSigma[*i][k]=tbout[k];
    else
      swathx.KuKaGMI.surfEmissSigma[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//    for(k=13;k<15;k++)
//      swathx.KuKaGMI.surfEmissivity[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
//end    WSO 10/13/15
}

void copytbouts1_fs_(float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
//begin  WSO 9/16/13
  for(k=0;k<13;k++)
    if(tbout[k] > -90.)
      swathx.KuGMI.simulatedBrightTemp[*i][k]=tbout[k];
    else
      swathx.KuGMI.simulatedBrightTemp[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels  
//  for(k=13;k<15;k++)
//    swathx.KuGMI.simulatedBrightTemp[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end    WSO 9/16/13
}

void copytbouts2_fs_(float *tbout, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
//begin  WSO 9/16/13
  for(k=0;k<13;k++)
    if(tbout[k] > 1.)
      swathx.KuKaGMI.simulatedBrightTemp[*i][k]=tbout[k];
    else
      swathx.KuKaGMI.simulatedBrightTemp[*i][k]=missing_r4c;
  //for(k=0;k<2;k++)
  //  swathx.KuKaGMI.simulatedBrightTemp[*i][k]=missing_r4c;
//begin  WSO 7/28/16 remove extra channels
//  for(k=13;k<15;k++)
//    swathx.KuKaGMI.simulatedBrightTemp[*i][k]=missing_r4c;
//end    WSO 7/28/16
//end  WSO 9/16/13
}

void copyrainflags1_fs_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.Input.precipitationFlag[*i]=*sfcVar;
}

void copyrainflags2_fs_(int *sfcVar, int *i, int *scanPatternFlag)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuKaGMI.Input.precipitationFlag[*i][0]=*sfcVar;
  swathx.KuKaGMI.Input.precipitationFlag[*i][1]=*sfcVar;
  if(*scanPatternFlag==0 &&(*i<12 || *i>36))
    {
      swathx.KuKaGMI.Input.precipitationFlag[*i][0]=missing_i2c;
      swathx.KuKaGMI.Input.precipitationFlag[*i][1]=missing_i2c;
    }
}

//begin  WSO 8/20/14 write new ioquality flags
void copyioqualitys1_fs_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.FLG.ioQuality[*i]=*sfcVar;
}

void copyioqualitys2_fs_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuKaGMI.FLG.ioQuality[*i]=*sfcVar;
}
//end    WSO 8/20/14
//
//begin  WSO 3/17/17 write snow ice cover flags
void copysnowices1_fs_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.Input.snowIceCover[*i]=*sfcVar;
}

void copysnowices2_fs_(int *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuKaGMI.Input.snowIceCover[*i]=*sfcVar;
}
//end    WSO 3/17/17

void copysfcliqfracts1_fs_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.nearSurfPrecipLiqRate[*i]=*sfcVar;
}

void copysfcliqfracts2_fs_(float *sfcVar, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  swathx.KuKaGMI.nearSurfPrecipLiqRate[*i]=*sfcVar;
}

void copycldwaters1_fs_(float *var1d, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
        swathx.KuGMI.cloudLiqWaterCont[*i][k]=var1d[k];
    }
}

void copycldwaters2_fs_(float *var1d, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx; 

  for(k=0;k<nbins;k++)
    {
        swathx.KuKaGMI.cloudLiqWaterCont[*i][k]=var1d[k];
    }
}

void copycldices1_fs_(float *var1d, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
      swathx.KuGMI.cloudIceWaterCont[*i][k]=var1d[k];
    }
}

void copycldices2_fs_(float *var1d, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<nbins;k++)
    {
      swathx.KuKaGMI.cloudIceWaterCont[*i][k]=var1d[k];
    }
}

//begin  WSO 9/5/13 new copy routine for SRT and DSRT pia effective sigma's
void copysigmapias1_fs_(float *sigmapia, int *i)
{
  extern L2BCMB_SWATHS swathx;
//diagnostic
//  printf("in writeCMBOut i: %5i,  sigmapia: %10.4f\n", *sigmapia, *i);
//end diagnostic
  swathx.KuGMI.Input.piaEffectiveSigma[*i] = *sigmapia;
}
void copysigmapias2_fs_(float *sigmapiaku, float *sigmapiaka, int *i)
{
  extern L2BCMB_SWATHS swathx;
//    diagnostic
//  printf("in writeCMBOut i: %5i,  sigmapiaku: %10.4f,  sigmapiaka: %10.4f\n", 
//     *sigmapiaku, *sigmapiaka, *i);
//end diagnostic
    swathx.KuKaGMI.Input.piaEffectiveSigma[*i][0] = *sigmapiaku;
    swathx.KuKaGMI.Input.piaEffectiveSigma[*i][1] = *sigmapiaka;
}
//end    WSO 9/5/13

//write principal components
void copyprincomps1_fs_(float *princomp, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<5;k++)
    {
      //.. swathx.KuGMI.aPriori.prinComp[*i][k] = princomp[k];
    }
}
//
void copyprincomps2_fs_(float *princomp, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<5;k++)
    {
      //     swathx.KuKaGMI.aPriori.prinComp[*i][k] = princomp[k];
    }
}

//write profile class
void copyprofclasss1_fs_(int *profclass, int *i)
{
  extern L2BCMB_SWATHS swathx;

  swathx.KuGMI.aPriori.profClass[*i] = *profclass;
}

void copyprofclasss2_fs_(int *profclass, int *i)
{
  extern L2BCMB_SWATHS swathx;

  swathx.KuKaGMI.aPriori.profClass[*i] = *profclass;
}

//write surface precip bias ratio
void copysurfprecipbiasratios1_fs_(float *biasratio, int *i)
{ 
  extern L2BCMB_SWATHS swathx;

  //  swathx.KuGMI.aPriori.surfPrecipBiasRatio[*i] = *biasratio;
}

void copysurfprecipbiasratios2_fs_(float *biasratio, int *i)
{
  extern L2BCMB_SWATHS swathx;

  //  swathx.KuKaGMI.aPriori.surfPrecipBiasRatio[*i] = *biasratio;
} 


//write initial log10 of the PSD intercept
void copyinitnws1_fs_(float *initlogNw, int *n9, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<9;k++)
    {
     if(n9[k]>1 && n9[k]<=88)
       {
         //swathx.KuGMI.aPriori.initNw[*i][k] = initlogNw[n9[k]-1];
       }
     else
       {
         //swathx.KuGMI.aPriori.initNw[*i][k] = missing_r4c;
       }
    }
}

void copyinitnws2_fs_(float *initlogNw, int *n9, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<9;k++)
    {
     if(n9[k]>1 && n9[k]<=88)
       {
	 // swathx.KuKaGMI.aPriori.initNw[*i][k] = initlogNw[n9[k]-1];
       }
     else
       { 
	 //         swathx.KuKaGMI.aPriori.initNw[*i][k] = missing_r4c;
       }
    }
}

//write sub-footprint variability parameter
void copysubfootvariabilitys1_fs_(float *subfoot, int *i)
{ 
  extern L2BCMB_SWATHS swathx;

  swathx.KuGMI.nubfPIAfactor[*i] = *subfoot;
} 
  
void copysubfootvariabilitys2_fs_(float *subfoot, int *i)
{
  extern L2BCMB_SWATHS swathx;
    
  swathx.KuKaGMI.nubfPIAfactor[*i] = *subfoot;
}

//write multiple scattering flag
void copymultiscatcalcs1_fs_(int *multiscat, int *i)
{ 
  extern L2BCMB_SWATHS swathx;
  
  swathx.KuGMI.FLG.multiScatCalc[*i] = *multiscat;
}

void copymultiscatcalcs2_fs_(int *multiscat, int *i)
{ 
  extern L2BCMB_SWATHS swathx;
  
  swathx.KuKaGMI.FLG.multiScatCalc[*i] = *multiscat;
}

//write multiple scattering surface parameter
void copymultiscatsurfaces1_fs_(float *multisfc, int *i)
{
  extern L2BCMB_SWATHS swathx;

  swathx.KuGMI.multiScatMaxContrib[*i] = *multisfc;
}

void copymultiscatsurfaces2_fs_(float *multisfc, int *i)
{
  extern L2BCMB_SWATHS swathx;

  swathx.KuKaGMI.multiScatMaxContrib[*i] = *multisfc;
}

//
//begin  WSO 2/8/17 copy routine for measured sigma-zeros
void copysigmazeros1_fs_(float *sigmazeroku, int *i)
{
  extern L2BCMB_SWATHS swathx;
  swathx.KuGMI.Input.sigmaZeroMeasured[*i] = *sigmazeroku;
}
void copysigmazeros2_fs_(float *sigmazeroku, float *sigmazeroka, int *i)
{
    extern L2BCMB_SWATHS swathx;
//        swathx.KuKaGMI.Input.sigmaZeroMeasured[*i][0] = *sigmazeroku;
//            swathx.KuKaGMI.Input.sigmaZeroMeasured[*i][1] = *sigmazeroka;
          swathx.KuKaGMI.Input.sigmaZeroMeasured[*i] = *sigmazeroka;
}
//end    WSO 2/8/17

//begin  WSO 8/19/13 modified copy routines to include nodes
void copylognws1_fs_(float *logNw, int *n9, int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<9;k++)
    {
//begin WSO 9/7/14 added upper limit for n9 in next line
      if(n9[k]>1 && n9[k]<=88)
        {
	  //swathx.KuGMI.PSDparamLowNode[*i][k] = n9[k]-1;
	  //swathx.KuGMI.precipTotPSDparamLow[*i][k][0]=logNw[n9[k]-1];
        }
      else
        {
	  //swathx.KuGMI.PSDparamLowNode[*i][k] = missing_i2c;
	  //swathx.KuGMI.precipTotPSDparamLow[*i][k][0]= missing_r4c;
        }
    }
}

void copylognws2_fs_(float *logNw, int *n9,int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;
  
  for(k=0;k<9;k++)
    {
      if(n9[k]>0)
        {
	  //swathx.KuKaGMI.PSDparamLowNode[*i][k] = n9[k]-1;
	  //swathx.KuKaGMI.precipTotPSDparamLow[*i][k][0]=logNw[n9[k]-1];
        }
      else
        {
	  //swathx.KuKaGMI.PSDparamLowNode[*i][k] = missing_i2c;
	  //swathx.KuKaGMI.precipTotPSDparamLow[*i][k][0]= missing_r4c;
        }
    }
}
//end    WSO 8/19/13

//begin  WSO 8/19/13 add mu as second low-resolution parameter
void copymus1_fs_(float *mu_prof, int *n9,int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<9;k++)
    {
      if(n9[k]>1)
        {
	  //         swathx.KuGMI.precipTotPSDparamLow[*i][k][1]=mu_prof[n9[k]-1];
        }
      else
        {
	  //         swathx.KuGMI.precipTotPSDparamLow[*i][k][1]=missing_r4c;
        }
    }
}

void copymus2_fs_(float *mu_prof, int *n9,int *i)
{
  int k;
  extern L2BCMB_SWATHS swathx;

  for(k=0;k<9;k++)
    {
      if(n9[k]>0)
        {
	  //swathx.KuKaGMI.precipTotPSDparamLow[*i][k][1]=mu_prof[n9[k]-1];
        }
      else
        {
	  //swathx.KuKaGMI.precipTotPSDparamLow[*i][k][1]=missing_r4c;
        }
    }
}
//end    WSO 8/19/13



#include <stdlib.h>
#include <stdio.h>
#include "TK_2BCMB.h"
#include "math.h"

extern float missingmod_mp_missing_r4_;
extern short missingmod_mp_missing_i2_;
extern long missingmod_mp_missing_i4_;

#define missing_r4c missingmod_mp_missing_r4_
#define missing_i2c missingmod_mp_missing_i2_
#define missing_i4c missingmod_mp_missing_i4_;
float pRateCCTables[39][50];
float pRate_bzd_Tables[2][2][71];
FILE *foutD;
void readcluttertables_(void)
{
  FILE *fout;
  int i,j;
  fout=fopen("AncData/cluttCorectTables.txt","r");
  for(i=0; i<39;i++)
    for(j=0;j<50;j++)
      fscanf(fout,"%g",&pRateCCTables[38-i][j]);
  fclose(fout);
  foutD=fopen("cluttCorectDiag.txt","w");
}

void readclutter_bzd_tables_(void)
{
  FILE *fout;
  int i,j,idum;
  fout=fopen("AncData/cluttCorectTables_bzd.txt","r");
  for(i=0; i<71;i++)
    fscanf(fout,"%i %g %g %g %g",&idum,&pRate_bzd_Tables[0][0][i],
	   &pRate_bzd_Tables[0][1][i],&pRate_bzd_Tables[1][0][i],
	   &pRate_bzd_Tables[1][1][i]);
  fclose(fout);
  foutD=fopen("cluttCorectDiag.txt","w");
}

void estimated_sfc_precip1_(int *i, float *pRate1d, float *pRateStd1d, float *sfcRain, 
			    int *sfc, int *bzd, int *cfb, int *ptype, float *liqFract)
{
  extern L2BCMB_SWATHS swathx;
  int j0, itop, ibott, fzClass;
  float est_Surf_Precip=*sfcRain;
  float pRateCS[40];
  swathx.KuGMI.lowestUnclutteredBin[*i]=swathx.KuGMI.phaseBinNodes[*i][4];
  swathx.KuGMI.lowestEstimateBin[*i]=swathx.KuGMI.phaseBinNodes[*i][4];
  float liqFractR;
  int sfc2=2*(*sfc), sfcType,iType;
  if(swathx.KuGMI.Input.surfaceType[*i]==0)
    {
      sfc2=175;
      sfcType=0;
    }
  else
    sfcType=1;
  //printf("%i %i %i %i\n",sfc2,*bzd,*cfb,swathx.KuGMI.Input.surfaceType[*i]);
  int i1,i2;
  double localZenith=swathx.KuGMI.Input.localZenithAngle[*i];
  i1=(int)(0.5+((*cfb)-(*bzd))*cos(localZenith/180*3.1416));
  i2=(int)(0.5+(sfc2-(*bzd))*cos(localZenith/180*3.1416));
  i1+=30;
  i2+=30;
  i2-=4;
  if(i1>70)
    i1=70;
  if(i2>70)
    i2=70;
  if(i1<0)
    i1=0;
  if(i2<0)
    i2=0;
  if(i2<i1)
    i2=i1;
  if(*ptype>0)
    {
      int n2;
      itop=*cfb-0-130;
      ibott=(*sfc+1)*2-1-130;
      //printf("isfc=%i %i\n",*sfc,*cfb);
      fzClass=*bzd-128;
      if (fzClass>49) fzClass=49;

      swathx.KuGMI.lowestEstimateBin[*i]=(int)(*sfc);
      int dRange=(int)(*sfc-(*bzd)/2);
      liqFractR=*liqFract;
      if(dRange>=4)
	liqFractR=1;
      else
	if(dRange>=0 && dRange<4)
	  liqFractR=dRange/4.0;
      
      swathx.KuGMI.FLG.estimPrecipInClutter[*i]=1;
      float rRatio;
      if(*ptype==1)
	iType=0;
      else
	iType=1;
      rRatio=pRate_bzd_Tables[sfcType][iType][i2]/	\
	pRate_bzd_Tables[sfcType][iType][i1];
      if(rRatio>1.5)
	rRatio=1.5;
      swathx.KuGMI.estimSurfPrecipTotRate[*i]=rRatio*
	swathx.KuGMI.nearSurfPrecipTotRate[*i];
      swathx.KuGMI.estimSurfPrecipTotRateSigma[*i]=rRatio*
	pRateStd1d[(int)(*cfb/2)-1];
      swathx.KuGMI.estimSurfPrecipLiqRate[*i]=rRatio*
	swathx.KuGMI.nearSurfPrecipTotRate[*i]*liqFractR;
    }
  else
    {
      //missing_r4c;
      swathx.KuGMI.lowestEstimateBin[*i]=87;
      swathx.KuGMI.estimSurfPrecipTotRate[*i]=0;
      swathx.KuGMI.estimSurfPrecipLiqRate[*i]=0;
      swathx.KuGMI.estimSurfPrecipTotRateSigma[*i]=0;
      swathx.KuGMI.FLG.estimPrecipInClutter[*i]=0;
    }

  if(swathx.KuGMI.estimSurfPrecipTotRate[*i]<0)
    swathx.KuGMI.estimSurfPrecipTotRate[*i]=missing_r4c;
  if(swathx.KuGMI.estimSurfPrecipLiqRate[*i]<0)
    swathx.KuGMI.estimSurfPrecipLiqRate[*i]=missing_r4c;
  if(swathx.KuGMI.estimSurfPrecipTotRateSigma[*i]<0)
    swathx.KuGMI.estimSurfPrecipTotRateSigma[*i]=missing_r4c;
  //printf("%i %i \n",  swathx.KuGMI.lowestEstimateBin[*i],swathx.KuGMI.phaseBinNodes[*i][4]);
  //swathx.KuGMI.estimSurfPrecipTotWaterCont[*i]=missing_r4c;
 //swathx.KuGMI.estimSurfPrecipLiqWaterCont[*i]=missing_r4c;
}

/*
87 88 82 0.201421 0.256668 
87 88 82 0.645152 0.321897 
87 88 83 0 0 
87 88 82 0 0 
87 88 82 0 0 
87 88 83 0.246041 0 
87 88 82 0.381081 0.403427 
87 88 81 0.66163 0.161783 
87 88 83 0.257529 0.158413 
87 88 83 0 0 
87 88 82 0 0 
87 88 83 0.237464 0.294474 
87 88 81 0.323459 0.444358 
87 88 81 0.332284 0.195314 
*/
/*
if itop<39 and pType[i,j]>0:
pRateCS=pRateShape[itop:min(ibott,39),fzClass]
  if pRateCS[0]>0.1:
  cPRate=precipRate[i,j,bcf[i,j]-50]*pRateCS/pRateCS[0]
    nExt=cPRate.shape[0]
    kmax=1
    sfcPrecipRateC[i,j]=precipRate[i,j,bcf[i,j]-50]
    pRL=precipRate[i,j,bcf[i,j]-50]
    for k in range(1,binSfc[i,j]-bcf[i,j]+1):
    if 2*k<nExt:
      precipRate[i,j,bcf[i,j]-50+k]=cPRate[2*k]
	pRL=precipRate[i,j,bcf[i,j]-50+k]
	kmax=k
	s3+=pRL
	sfcPrecipRateC[i,j]=pRL
*/


void estimated_sfc_precip2_(int *i, float *pRate1d, float *pRateStd1d, float *sfcRain, 
			    int *sfc, int *bzd, int *cfb, int *ptype, \
			    int *flagScanPattern, float *liqFract)
{
  extern L2BCMB_SWATHS swathx;
  int j0, itop, ibott, fzClass;
  float est_Surf_Precip=*sfcRain;
  float pRateCS[40];
  float liqFractR;
  //float liqFractR;
  int sfc2=2*(*sfc),sfcType, iType;
  if(swathx.KuGMI.Input.surfaceType[*i]==0)
    {
      sfc2=175;
      sfcType=0;
    }
  else
    sfcType=1;
  //printf("%i %i %i %i\n",sfc2,*bzd,*cfb,swathx.KuGMI.Input.surfaceType[*i]);
  int i1,i2;
  double localZenith=swathx.KuGMI.Input.localZenithAngle[*i];
  i1=(int)(0.5+((*cfb)-(*bzd))*cos(localZenith/180*3.1416));
  i2=(int)(0.5+(sfc2-(*bzd))*cos(localZenith/180*3.1416));
  i1+=30;
  i2+=30;
  i2-=4;
  if(i1>70)
    i1=70;
  if(i2>70)
    i2=70;
  if(i1<0)
    i1=0;
  if(i2<0)
    i2=0;
  if(i2<i1)
    i2=i1;

  //printf("scanPattern=%i %i \n",*flagScanPattern,*i);
  swathx.KuKaGMI.FLG.scanPatternFlag=*flagScanPattern;
  swathx.KuGMI.FLG.scanPatternFlag=0;//*flagScanPattern;
  if(*flagScanPattern==0 &&(*i<11 || *i>37))
    {
      swathx.KuKaGMI.estimSurfPrecipTotRate[*i]=missing_r4c;
      swathx.KuKaGMI.estimSurfPrecipLiqRate[*i]=missing_r4c;
      swathx.KuKaGMI.estimSurfPrecipTotRateSigma[*i]=missing_r4c;
      swathx.KuKaGMI.FLG.estimPrecipInClutter[*i]=missing_i2c;
      return;
    }
  swathx.KuKaGMI.lowestUnclutteredBin[*i]=swathx.KuKaGMI.phaseBinNodes[*i][4];
  swathx.KuKaGMI.lowestEstimateBin[*i]=swathx.KuKaGMI.phaseBinNodes[*i][4];
  if(*ptype>0)
    {
      int n2;
      itop=*cfb-0-130;
      ibott=(*sfc+1)*2-1-130;
      fzClass=*bzd-128;
      if (fzClass>49) fzClass=49;
      swathx.KuKaGMI.lowestEstimateBin[*i]=(int)(*sfc);
      int dRange=(int)(*sfc-(*bzd)/2);
      liqFractR=*liqFract;
      if(dRange>=4)
	liqFractR=1;
      else
	if(dRange>=0 && dRange<4)
	  liqFractR=dRange/4.0;
      
      swathx.KuKaGMI.FLG.estimPrecipInClutter[*i]=1;
      float rRatio;
      if(*ptype==1)
	iType=0;
      else
	iType=1;
      rRatio=pRate_bzd_Tables[sfcType][iType][i2]/	\
	pRate_bzd_Tables[sfcType][iType][i1];
      if(rRatio>1.5)
	rRatio=1.5;
      swathx.KuKaGMI.estimSurfPrecipTotRate[*i]=rRatio*
	swathx.KuKaGMI.nearSurfPrecipTotRate[*i];
      swathx.KuKaGMI.estimSurfPrecipLiqRate[*i]=rRatio*
	swathx.KuKaGMI.nearSurfPrecipTotRate[*i]*liqFractR;
      swathx.KuKaGMI.estimSurfPrecipTotRateSigma[*i]=rRatio*
	pRateStd1d[(int)(*cfb/2)-1];
    }
  else
    {
      //missing_r4c;
      swathx.KuKaGMI.lowestEstimateBin[*i]=87;
      swathx.KuKaGMI.estimSurfPrecipTotRate[*i]=0;
      swathx.KuKaGMI.estimSurfPrecipLiqRate[*i]=0;
      swathx.KuKaGMI.estimSurfPrecipTotRateSigma[*i]=0;
      swathx.KuKaGMI.FLG.estimPrecipInClutter[*i]=0;
    }

  if(swathx.KuKaGMI.estimSurfPrecipTotRate[*i]<0)
    swathx.KuKaGMI.estimSurfPrecipTotRate[*i]=missing_r4c;
  if(swathx.KuKaGMI.estimSurfPrecipLiqRate[*i]<0)
    swathx.KuKaGMI.estimSurfPrecipLiqRate[*i]=missing_r4c;
  if(swathx.KuKaGMI.estimSurfPrecipTotRateSigma[*i]<0)
    swathx.KuKaGMI.estimSurfPrecipTotRateSigma[*i]=missing_r4c;
  //  swathx.KuKaGMI.estimSurfPrecipTotWaterCont[*i]=missing_r4c;
  //swathx.KuKaGMI.estimSurfPrecipLiqWaterCont[*i]=missing_r4c;
}


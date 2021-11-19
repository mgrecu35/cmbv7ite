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
void copy_oeqv_kugmi_(float *OEQv,short *n10,int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<10;k++)
{
swathx.KuGMI.OptEst.OEvaporDensity[*i][k]=OEQv[n10[k]-1];
}
}

void copy_oetemp_kugmi_(float *OETemp,short *n10,int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<10;k++)
{
swathx.KuGMI.OptEst.OEairTemperature[*i][k]=OETemp[n10[k]-1];
}
}

void copy_oecloud_kugmi_(float *OECloud, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<nbins;k++)
{
swathx.KuGMI.OptEst.OEcloudLiqWaterCont[*i][k]=OECloud[k];
}
}

void copy_oecloudliqpath_kugmi_(float *OECloudLiqPath, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEcolumnCloudLiqWater[*i]=*OECloudLiqPath;
}

void copy_oecloudliqsigma_kugmi_(float *OEcloudLiqSigma, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEcolumnCloudLiqSigma[*i]=*OEcloudLiqSigma;
}

void copy_oesfcwind_kugmi_(float *OESfcWind, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEtenMeterWindSpeed[*i]=*OESfcWind;
}

void copy_oesfcwindsigma_kugmi_(float *OESfcWindSigma, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEtenMeterWindSigma[*i]=*OESfcWindSigma;
}

void copy_oeskntemp_kugmi_(float *OESknTemp, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEskinTemperature[*i]=*OESknTemp;
}

void copy_oeskintempsigma_kugmi_(float *OEskinTempSigma, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEskinTempSigma[*i]=*OEskinTempSigma;
}

void copy_oesfctemp_kugmi_(float *OESfcTemp, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEsurfaceAirTemperature[*i]=*OESfcTemp;
}

void copy_oesfcqv_kugmi_(float *OESfcQv, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEsurfaceVaporDensity[*i]=*OESfcQv;
}

void copy_oetpw_kugmi_(float *OEtpw, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEcolumnWaterVapor[*i]=*OEtpw;
}

void copy_oetpwsigma_kugmi_(float *OEtpwSigma, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEcolumnVaporSigma[*i]=*OEtpwSigma;
}

void copy_oechisq_kugmi_(float *OEchiSq, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEerrorOfDataFit[*i]=*OEchiSq;
}

void copy_oepianonrain_kugmi_(float *OEpiaNonRain, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<2;k++)
{
swathx.KuGMI.OptEst.OEpiaNoPrecip[*i]=OEpiaNonRain[0];
}
}

void copy_oesimtbnonrain_kugmi_(float *OEsimTbNonRain, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<13;k++)
{
swathx.KuGMI.OptEst.OEsimulatedBrightTemp[*i][k]=OEsimTbNonRain[k];
}
}

void copy_oeemis_kugmi_(float *OEemis, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<13;k++)
{
swathx.KuGMI.OptEst.OEsurfEmissivity[*i][k]=OEemis[k];
}
}

void copy_oeemissigma_kugmi_(float *OEemisSigma, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<13;k++)
{
swathx.KuGMI.OptEst.OEsurfEmissSigma[*i][k]=OEemisSigma[k];
}
}

void copy_oeemisa_kugmi_(float *OEemisA, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<13;k++)
{
swathx.KuGMI.OptEst.OEemissivityAvgKernel[*i][k]=OEemisA[k];
}
}

void copy_oestype_kugmi_(int *OEstype, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuGMI.OptEst.OEsurfaceClass[*i]=*OEstype;
}

void copy_oeqv_kukagmi_(float *OEQv,short *n10,int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<10;k++)
{
swathx.KuKaGMI.OptEst.OEvaporDensity[*i][k]=OEQv[n10[k]-1];
}
}

void copy_oetemp_kukagmi_(float *OETemp,short *n10,int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<10;k++)
{
swathx.KuKaGMI.OptEst.OEairTemperature[*i][k]=OETemp[n10[k]-1];
}
}

void copy_oecloud_kukagmi_(float *OECloud, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<nbins;k++)
{
swathx.KuKaGMI.OptEst.OEcloudLiqWaterCont[*i][k]=OECloud[k];
}
}

void copy_oecloudliqpath_kukagmi_(float *OECloudLiqPath, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEcolumnCloudLiqWater[*i]=*OECloudLiqPath;
}

void copy_oecloudliqsigma_kukagmi_(float *OEcloudLiqSigma, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEcolumnCloudLiqSigma[*i]=*OEcloudLiqSigma;
}

void copy_oesfcwind_kukagmi_(float *OESfcWind, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEtenMeterWindSpeed[*i]=*OESfcWind;
}

void copy_oesfcwindsigma_kukagmi_(float *OESfcWindSigma, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEtenMeterWindSigma[*i]=*OESfcWindSigma;
}

void copy_oeskntemp_kukagmi_(float *OESknTemp, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEskinTemperature[*i]=*OESknTemp;
}

void copy_oeskintempsigma_kukagmi_(float *OEskinTempSigma, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEskinTempSigma[*i]=*OEskinTempSigma;
}

void copy_oesfctemp_kukagmi_(float *OESfcTemp, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEsurfaceAirTemperature[*i]=*OESfcTemp;
}

void copy_oesfcqv_kukagmi_(float *OESfcQv, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEsurfaceVaporDensity[*i]=*OESfcQv;
}

void copy_oetpw_kukagmi_(float *OEtpw, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEcolumnWaterVapor[*i]=*OEtpw;
}

void copy_oetpwsigma_kukagmi_(float *OEtpwSigma, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEcolumnVaporSigma[*i]=*OEtpwSigma;
}

void copy_oechisq_kukagmi_(float *OEchiSq, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEerrorOfDataFit[*i]=*OEchiSq;
}

void copy_oepianonrain_kukagmi_(float *OEpiaNonRain, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<2;k++)
{
swathx.KuKaGMI.OptEst.OEpiaNoPrecip[*i][k]=OEpiaNonRain[k];
}
}

void copy_oesimtbnonrain_kukagmi_(float *OEsimTbNonRain, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<13;k++)
{
swathx.KuKaGMI.OptEst.OEsimulatedBrightTemp[*i][k]=OEsimTbNonRain[k];
}
}

void copy_oeemis_kukagmi_(float *OEemis, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<13;k++)
{
swathx.KuKaGMI.OptEst.OEsurfEmissivity[*i][k]=OEemis[k];
}
}

void copy_oeemissigma_kukagmi_(float *OEemisSigma, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<13;k++)
{
swathx.KuKaGMI.OptEst.OEsurfEmissSigma[*i][k]=OEemisSigma[k];
}
}

void copy_oeemisa_kukagmi_(float *OEemisA, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
for(k=0;k<13;k++)
{
swathx.KuKaGMI.OptEst.OEemissivityAvgKernel[*i][k]=OEemisA[k];
}
}

void copy_oestype_kukagmi_(int *OEstype, int *i)
{int k;
extern L2BCMB_SWATHS swathx;
swathx.KuKaGMI.OptEst.OEsurfaceClass[*i]=*OEstype;
}


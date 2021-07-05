#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef GFOR 
#define nb 5
#endif
#ifdef IFORT 
#define nb 8
#endif
typedef struct stormStructType
{
  int     *nodes;     // 5 nodes that define the storm structure 
                      // for stratiform, with BB, nodes[0] is the storm top, 
                      // nodes[1] is the BB top,
                      // nodes[2] is the BB peak, nodes[3] is the BB bottom, 
                      // and nodes[4] is the lowest
                      // clutter free gate (0<=nodes[0]<nodes[1]<nodes[2]
                      // <nodes[3]<nodes[4]<ngates)
                      // for convective, nodes[3] is not used
                      // the C convention is used the gates are 
                      // numbered from 0 to ngates-1
  long    inf1[nb];
  int     iSurf;      // surface gate number
  float   freezH;     // freezing height - not used
  int     rainType;   // rainType -1 stratiform
                      //          -2 convective
} stormStructType;

typedef struct radarDataType
{
  int     ngates;     // number of gates
  float   *z13obs;    // observed Ku-reflectivities
  long    inf1[nb];
  float   *z35obs;    // observed Ka-reflectivities
  long    inf2[nb];
  float   xlong;      // longitude -not used
  float   xlat;       // latitude -not used
  float   pia13srt;   // Ku-band SRT PIA 
  float   relPia13srt;   // Ku-band SRT PIA 
  float   pia35srt;   // Ka-band SRT PIA
  float   relPia35srt;   // Ka-band SRT PIA
  float   dr;         // gate size
  float   *hh;        // hh[i] is the height of gate [i]
  long    inf3[nb];
  float   hfreez;
  float   sigmaZeroKu; //SJM 12/3/2014
  float   sigmaZeroKa; //SJM 12/3/2014
} radarDataType;

typedef struct radarRetType
{

  int     ngates;     // number of gates
  int     nMemb;
  int     nmfreq;     // # of simulated passive microwave frequencies
  float   *z13c;      // effective reflectivity (attenuation corrected)
                      // at Ku-band
  long    inf1[nb];
  float   *z35mod0;   // simulated observations at Ka-band
  long    inf2[nb];
//  SFM  begin  06/16/2014; for M.Grecu, multiple scattering
  float   *dz35ms;   // simulated observations at Ka-band
  long    inf2t[nb];
//  SFM  end    06/16/2014
  float   *z35;   // simulated observations at Ka-band
  long    inf2s[nb];
  float   *pwc;       // precipitation water content (g/m3) 
  long    inf3[nb];
  float   *rrate;     // rain rate (mm/h)
  long    inf4[nb];
  float   *d0;        // median diameter (mm)
  long    inf5[nb];
  float   *log10dNw;  // 
  long    inf6[nb];
  float   *tb;        // simulated brightness temperatures
  long    inf7[nb];
  float   *emTb;      // simulated emission brightness temperatures
  long    inf71[nb];
  float   *emis;      // simulated emissivity
  long    inf72[nb];
  int     *imu;       // mu index of the look up table
  long    inf8[nb];
  int     *iwc;       // index of surface wind profile (from 1 to nc)
  long    inf8p[nb];
  int     *icc;       // index of RH profile (from 1 to nc)
                      // nc is the number of possible RH profiles see cloud.f90
  long    inf9[nb];
  int     *jcc;       // index of cloud profile (from 1 to nc)
                      // nc is the number of possible RH profiles see cloud.f90
  long    inf10[nb];
  float  *sfc_wind;   // surface wind speed
  long    inf11[nb];
  float  *sfc_windU;   // surface wind speed U
  long    inf11a[nb];
  float  *sfc_windV;   // surface wind speed V
  long    inf11b[nb];
  float  *pia13;   // surface wind speed
  long    inf12[nb];
  float  *pia35;   // surface wind speed
  long    inf12b[nb];
  float  *simSigmaZeroKu;   // simulated backscatter
  long    inf12c[nb];
  float  *simSigmaZeroKa;   // simulated backscatter
  long    inf13[nb];
  float  *z35mMean;   // surface wind speed
  long    inf14[nb];
  float  *z35mSig;    // surface wind speed
  long    inf15[nb];
  float  pia13mMean, pia35mMean, pia13mSig, pia35mSig;    // surface wind speed
  int    ntpw;
  float  *tpw;
  long   inf16[nb];
  float  *tpwCldMod;
  long   inf17[nb];
  float   *logdNw;
 
} radarRetType;


typedef struct retParamType
{
  /*
    F=wz*SUM(zsim,ka-zobs,ka)**2+w13*(pia13-pia13srt)**2+w35*(pia35,sim-pia35srt)**2
    wz the weight of the reflectivity term
    w13 the weight of the pia squared difference at ku-band
    w35 the weight of the pia squared difference at ka-band
    z13thresh - threshold (dBZ) to determine the Ku-band observations used in the algorithm
    z35thresh - threshold (dBZ) to determine the Ka-band observations used in the algorithm
  */
  
  float wz, w13,  w35,  z13thresh,  z35thresh;
} retParamType;

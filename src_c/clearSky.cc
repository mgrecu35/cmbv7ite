#include <math.h>
#include "typedef.h"
#define WANT_STREAM
#define WANT_MATH                    // include.h will get math fns
                                     // newmatap.h will get include.h

#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif

#ifdef GFOR 
//#define ran1() __ran_mod_MOD_ran1()
#define ran1() __ran_mod__ran1()
#endif

#ifdef IFORT 
#define ran1() ran_mod_mp_ran1_()
#endif

extern "C" double ran_mod_mp_ran1_();  // returns a uniformly distributed 
                                       // random number between
                                       // 0. and 1.0
extern "C" double __ran_mod__ran1();   // returns a uniformly distributed 
                                       // random number between
                                       // 0. and 1.0

extern "C" double __ran_mod_MOD_ran1();// returns a uniformly distributed 
                                       // random number between
                                       // 0. and 1.0


extern "C" void calccstb_(float *sfc_wind,float *umu,int *ic,int *jc,
			  float *hh,int *ngates, int *nmfreqm, float *tb,
			  float *tpw);

/* this is a fortran subroutine (file rterain.f90) that simulates brigthness temperatures
as a function of retrieved precipitation profiles 
sfc_wind         : surface wind -- randomly set
umu              : cosine of the viewing agle
kext, salb, asym : are the extinction coefficient, scattering albedo, and asymmetry
                   factors determined during the radar retrieval process
ngates,hh,
nmfreqm          : as defined above
tb               : simulated brightness temperatures
*/           






extern "C" void clearsky_( radarDataType   *radarData, 
			   stormStructType *stormStruct,
			   retParamType    *retParam,  
			   int *nmu, radarRetType    *radarRet)
  
{
  extern float *cldclass_mp_tpw_;
  int imemb;
  int nMemb=radarRet->nMemb;
  float umu;
  for(int i=0;i<radarRet->ngates;i++)
    {
      radarData->hh[i]=(80-i)*radarData->dr;
     
    }
 
  umu=cos(53/180.*3.14);
  float c0=105.352+3, c1=0.2418, c2=0.3218, c3=-0.0292, c4=-0.0124,
    c5=-1.1742, c6=0.5742;
 
#pragma omp parallel for default (shared) private(imemb)
  for(imemb=0;imemb<radarRet->nMemb;imemb++)
    {
      /*  calccstb_(&radarRet->sfc_wind[imemb],&umu, &radarRet->icc[imemb],
		&radarRet->jcc[imemb],
		radarData->hh,&radarData->ngates,
		&radarRet->nmfreq,
		&radarRet->tb[imemb*2*radarRet->nmfreq],
		&radarRet->tpw[radarRet->icc[imemb]-1]);
      */
    }
 for(imemb=0;imemb<radarRet->nMemb;imemb++)
   {
     float tb10h=radarRet->tb[imemb*2*radarRet->nmfreq];
     float tb19h=radarRet->tb[imemb*2*radarRet->nmfreq+1];
     float tb21h=radarRet->tb[imemb*2*radarRet->nmfreq+2];
     float tb37h=radarRet->tb[imemb*2*radarRet->nmfreq+3];


     float tb10v=radarRet->tb[(imemb*2+1)*radarRet->nmfreq];
     float tb19v=radarRet->tb[(imemb*2+1)*radarRet->nmfreq+1];
     float tb21v=radarRet->tb[(imemb*2+1)*radarRet->nmfreq+2];
     float tb37v=radarRet->tb[(imemb*2+1)*radarRet->nmfreq+3];
     float tb85v=radarRet->tb[(imemb*2+1)*radarRet->nmfreq+4];
     float ws=c0+c1*tb10v+c2*tb10h+c3*tb19v+c4*tb21v+c5*tb37v+c6*tb37h;
     float tpw=10*(23.828-4.7515*(log(280-tb21v)+0.01384*(tb37v-pow(tb85v,0.8f))));
     tpw=10*(22.73-3.969*(log(280-tb21v)-0.003423*(tb37v-pow(tb85v,0.8f))));
     tpw=10*(23.82-4.059*log(280-tb21v)+0.02451*(log(280-tb21v)-tb37v));
     // tpw=10*(16.7636-3.3427*(log(280-tb21v)));
     // ws=147.90 + 1.0969*tb19v-0.4555*tb21v
     //  -1.7600*tb37v+0.7860*tb37h;

     /*     printf("%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f  %6.2f  %6.2f  %6.2f  %6.2f\n",
     	    tb10v,tb10h,tb37v,tb37h,ws,radarRet->sfc_wind[imemb],
     	    cldclass_mp_tpw_[radarRet->icc[imemb]-1]/100.,tpw,tb21v, tb10h, tb10v);
     */
   }
 //  exit(0);
  /*


  calccstb(sfc_wind,umu,ic,jc,&
  hh,ngates,nmfreqm,tb)
  */
  
}

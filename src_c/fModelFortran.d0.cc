#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef GFOR 
extern int nbinmod_mp_nbin_;
#define nbins nbinmod_mp_nbin_
#endif

#ifdef IFORT 
extern int nbinmod_mp_nbin_;
#define nbins nbinmod_mp_nbin_
#endif

extern "C" void fhb11_(float *z13,float *z35,float *z13obs,
		       float *pia13,float *pia35,int *ic,int *jc,
		       float *z35mod,float *pwc,float *n0w,float *dr,int *node,
		       int *isurf, int *imu,int *ngates,int *nmfreqm,
		       float *hh,int *itype,float *kext,float *salb,float *asym,
		       float *rrate, float *d0, float *hfreez, float *pia13srt);

extern "C" void fhb12_(float *z13,float *z35,float *z13obs,
		       float *pia13,float *pia35,int *ic,int *jc,
		       float *z35mod,float *pwc,float *n0w,float *dr,int *node,
		       int *isurf, int *imu,int *ngates,int *nmfreqm,
		       float *hh,int *itype,float *kext,float *salb,float *asym,
		       float *rrate, float *d0, float *hfreez, float *pia13srt);

extern "C" void  fhb11d0_(float *z13,float *z35,float *z13obs,
			  float *pia13,float *pia35,int *ic,int *jc,
			  float *z35mod,float *lwc,float *log10dNw,float *dr,
			  int *node,int *isurf,int *imu,
			  int *ngates,int *nmfreqm,float *hh,int *itype,
			  float *kext,float *salb,float *asym,
			  float *rrate,float *d0,float *hfreez,
			  float *pia13srt);

void fModelFortranD0(float *z13obs, int nodes[5], int isurf, int imu,
		   float *d0n, int nNodes, 
		   float *pia35M, float *pia13M,
		   float *z35mod, float *pwc, float dr, int ic, int jc, float *hh,
		   float delta, int iNode, int nmfreq, 
		   float *salb, float *kext,float *asym, int itype,
		   int ngates,float *rrate,float *d0,float *log10dN,
		   float *z13, float *z35,
		   int *imuv, float *hfreez,float *dz, float *pia13srt, float *relPia13srt)
/**********************************************************************************************

**** on input    ******************************************************************************
float *z13obs    : vector, observed Ku-band reflectivity
float *hh        : vector, height associated with the reflectivity observations
int ngates       : scalar, # of gates in the profile
int itype        : scalar, i
float dr         : scalar, radar gatesize
int nodes[5]     : vector, the 5 nodes defining the storm structure
int nNodes       : scalar, # nodes defining the log10dN profile
int *nodeP       : vector, contains the locations of the nNodes nodes
float *log10dNP  : vector, contains the values of log10dN of the nNodes nodes
int ic           : scalar, index of RH profile (from 1 to nc)
                 : nc is the number of possible RH profiles see cloud.f90
int jc           : scalar, index of cloud profile (from 1 to nc)
                 : nc is the number of possible cloud profiles see cloud.f90
int imu          : scalar, mu index of look up table (from 1 to nmu)
int iz           : integer, generally 0, but if an ensemble of solution
                 : is derived can be used to
                 : differentiate among ensemble members

float delta      : scalar, used to calculate the jacobian
int iNode        : scalar, index of the node where log10dN is perturbed
int nmfreq       : in


**** on output   *****************************************************************************


float *pia35M    : reference to scalar, model pia at Ka-band
float *pia13M    : reference to scalar, model pia at Ku-band
float *z35mod    : vector, model Ka-band reflectivities
float *pwc       : vector, precipitation water content (g/m3)

float *salb      : vector, fortran 2-d array, returns the scattering albedo at nmfreq mw freq
                 : these frequencies are set in mie2 subroutine in file Mie/mie3.f90 and 
                 : read by readTables subroutine in readTables.f90
                 : they can be accessed in variable mfreq in f90 
                 : module microwFreq in readTables.f90
float *kext      : vector, fortran 2-d array, returns the extinction coefficient (1/km)
                 : at nmfreq mw freq
float *asym      : vector, fortran 2-d array, returns the asymmetry fact

float *rrate     : vector, returns retrieved rain rate (mm/h)
float *d0        : vector, returns retrieved d0 (mm)
float *log10dN   : vector, returns retrieved log10(N0/N0ref)
********************************************************************************************/  
{
  int iLev, i, i1, i2, node[5];
  float att1,att351,dndum,z351;
  float fi, fi1, log10dNPi, f, ntot;
  float atm_extKav,cld_extKav,att35,att13,z13obsP[100], d0P[100],rrateP[100];;
 
 
 
  for(i=0;i<5;i++)
    {
      node[i]=nodes[i]+1;
      //printf("%3i ",node[i]);
      if(node[i]<1) node[i]=1;
      if(node[i]>100) node[i]>=100;
    }  
  if(rrate[nodes[4]-1]>0.0001 && d0[nodes[4]-1]<=0.0001)
    {
      printf("at the top exception in fModelFortran.cc d0 \n");
      printf("%g %g %g \n",z13obs[nodes[4]-1],rrate[nodes[4]-1],d0[nodes[4]-1]);
      exit(0);
    }


  for(i=0;i<1;i++)
    {
      for(iLev=nodes[i];iLev<=nodes[i+1];iLev++)
	{
	  f=(iLev-nodes[i]+0.)/(nodes[i+1]-nodes[i]+0.01);
	  if(f>1) f=1;
	  d0P[iLev]=(1-f)*(d0n[i])+f*(d0n[i+1]);
	}
    }
  
  *pia35M=0;
  *pia13M=0;
  //printf("%3i \n",imu);
  for(i=0;i<ngates;i++)
    imuv[i]=imu;
  for(i=0;i<ngates;i++)
    hh[i]=(nbins-i)*dr;

  if(node[4]>nbins) node[4]=nbins;

  for(i=node[0]-1;i<node[1]-1;i++)
    {

      hh[i]=(nbins-i)*dr;
      z13[i]=z13obs[i];
      z13obsP[i]=z13obs[i]+dz[i];
    }
  for(i=node[0];i<node[4];i++)
    {
      
      //printf("%6.2f %6.2f %6.2f %6.2f \n", hh[i], z13obsP[i],rrate[i],d0[i]);
    }


  /*if(*pia13srt<10)
    fhb11_(z13,z35,z13obsP,
	   pia13M,pia35M,&ic,&jc,z35mod,pwc,log10dN, &dr,node,&isurf,imuv,
	   &ngates,&nmfreq,hh,&itype,kext,salb,asym,rrate,d0,hfreez, pia13srt);
  else
    if(*relPia13srt>3)
      fhb12_(z13,z35,z13obsP,
	     pia13M,pia35M,&ic,&jc,z35mod,pwc,log10dN, &dr,node,&isurf,imuv,
	     &ngates,&nmfreq,hh,&itype,kext,salb,asym,rrate,d0,hfreez, pia13srt);
    else
  */
 
 
  // printf("%2i %2i \n",ic,jc);

  /*
    printf("in d0 %6.2f %6.2f %6.2f %6.2f %6.2f \n", 					
    z13obs[nodes[0]],z13obs[nodes[1]],z13obs[nodes[2]],z13obs[nodes[3]],z13obs[nodes[4]]);
    printf("in d0 %2i %2i %2i %2i %2i \n", 					
    nodes[0],nodes[1],nodes[2],nodes[3],nodes[4]);
    printf("in d0 %2i %2i %2i %2i %2i \n", 					
    node[0],node[1],node[2],node[3],node[4]);
  */
 

  fhb11d0_(z13,z35,z13obs,
	   pia13M,pia35M,&ic,&jc,z35mod,pwc,log10dN,&dr,node,&isurf,imuv,
	   &ngates,&nmfreq,hh,&itype,kext,salb,asym,rrateP,d0P,hfreez,pia13srt);
  printf("%6.3f %6.3f \n",d0n[0],d0n[1]);
  for(i=nodes[0];i<nodes[1];i++)
    {
      printf("%i %6.2f %6.3f %6.3f %6.3f %6.3f \n",i, z13obs [i],d0[i],d0P[i],rrate[i],rrateP[i]);
    }
  // printf("in d0 %2i %2i %2i %2i %2i \n", 					
  //	 node[0],node[1],node[2],node[3],node[4]);
  for(i=nodes[3];i<nodes[4];i++)
    {
      
      // printf("%6.2f %6.2f  %6.2f %6.2f %6.2f \n", hh[i], z13obs[i],z35mod[i], rrate[i],d0[i]);
    }
 // exit(0);
  // fhb11d0_(z13,z35,z13obsP,
  //	 pia13M,pia35M,&ic,&jc,z35mod,pwc,log10dN, &dr,node,&isurf,imuv,
  //	 &ngates,&nmfreq,hh,&itype,kext,salb,asym,rrate,d0,hfreez, pia13srt);
 
  // printf("%g %g %g %g\n",z13[node[4]-1],rrate[node[4]-1],d0[node[4]-1], *pia13M);

  if(rrate[nodes[4]-1]>0.0001 && d0[nodes[4]-1]<=0.0001)
    {
      printf("exception in fModelFortran.cc d0 \n");
      printf("%g %g %g \n",z13obs[nodes[4]-1],rrate[nodes[4]-1],d0[nodes[4]-1]);
      exit(0);
    }

  for(i=node[0];i<node[4];i++)
    {
      //      printf("%g %g \n",z13obs[i],rrate[i]);
    }
  /*
    printf("%g %g ",z13obs[nodes[4]],rrate[nodes[4]]);
    fhb1 is the Hitschfeld Bordan procedure in fhb1.f90 file 
  */ 
 

  /* free(z13);
  free(z35);
  free(imuv);
  */
  
}

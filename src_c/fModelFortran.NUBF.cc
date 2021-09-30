// SFM   07/26/2013  Code mods for M.Grecu


#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef GFOR 
extern int  __nbinmod_MOD_nbin;
#define nbins __nbinmod_MOD_nbin
#endif

#ifdef IFORT 
extern int nbinmod_mp_nbin_;
#define nbins nbinmod_mp_nbin_
#endif

//begin  MG 10/29/15 include appA
#include "appA.hpp"
//end    MG 10/29/15

int min1(int i1, int i2)
{
  if(i1<i2)
    return i1;
  else
    return i2;
}
#define min(a,b) (((a) < (b)) ? (a) : (b))

extern "C" void fhb11_(float *z13,float *z35,float *z13obs,
		       float *pia13,float *pia35,int *ic,int *jc,
		       float *z35mod,float *pwc,float *n0w,float *dr,int *node,
		       int *isurf, int *imu,int *ngates,int *nmfreqm,
		       float *hh,int *itype,float *kext,float *salb,float *asym,
		       float *rrate, float *d0, float *hfreez, float *pia13srt, 
		       int *imemb);

extern "C" void fhb12_(float *z13,float *z35,float *z13obs,
		       float *pia13,float *pia35,int *ic,int *jc,
		       float *z35mod,float *pwc,float *n0w,float *dr,int *node,
		       int *isurf, int *imu,int *ngates,int *nmfreqm,
		       float *hh,int *itype,float *kext,float *salb,float *asym,
		       float *rrate, float *d0, float *hfreez, float *pia13srt);

extern "C" float ran_mod_mp_normal2_(float *nm, float *nstd); 

void fModelFortran(float *z13obs, float *z35obs, 
		   int nodes[5], int isurf, int imu,
		   float *log10dNP, int *nodeP, int nNodes, 
		   float *pia35M, float *pia13M,
		   float *z35mod, float *pwc, float dr, int ic, int jc, 
		   float *hh,
		   float delta, int iNode, int nmfreq, 
		   float *salb, float *kext,float *asym, int itype,
		   int ngates,float *rrate,float *d0,float *log10dN,
		   float *z13, float *z35,
		   int *imuv, float *hfreez,float *dz, 
		   float *pia13srt, float *relPia13srt,
		   float *pia35srt, float *relPia35srt,
//  SFM  begin  06/22/2014; for M.Grecu  (unknown justification)
//  SFM  begin  07/01/2014; for M.Grecu  random sequences
		   int *imemb,float *localZAngle, float *wfractPix, 
		   float *xs, long *ichunk, float nstdA)
//  SFM  end    07/01/2014
//  SFM  end    06/22/2014
/**********************************************************************************************

**** on input    ******************************************************************************
float *z13obs    : vector, observed Ku-band reflectivity
float *hh        : vector, height associated with the reflectivity observations
int ngates       : scalar, # of gates in the profile
int itype        : scalar, i
float dr         : scalar, radar gatesize
int nodes[5]     : vector, the 5 nodes defining the storm structure
int nNodes       : scalar, # nodes defining the log10dN profile
int *nodeP      : vector, contains the locations of the nNodes nodes
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
long *ichunk     : for diagnostics only, number of segments into current orbit

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
  float atm_extKav,cld_extKav,att35,att13,z13obsP[88];
  int nSub=7;
//  SFM  begin  07/01/2014; for M.Grecu  random sequences
//  float xs[15];
//  SFM  end    07/01/2014
//  SFM  begin  06/22/2014; for M.Grecu (unknown justification)
  float nm=0., nstd=.125, sxs=0;
//  SFM  end    06/22/2014
//  SFM  begin  08/13/2014; for M.Grecu   "long orbit" correction
  int countSub[88];
  float kextSub[88][8][50],salbSub[88][8][50],asymSub[88][8][50],
    kextAvg[88][8], salbAvg[88][8], asymAvg[88][8], rainAvg[88], pwcAvg[88],
    d0Avg[88], rrateSub[88][50], pwcSub[88][50], d0Sub[88][50],
    z35modSub[88][50], z35Sub[88][50], z13Sub[88][50], pia35Sub[50],
    pia13Sub[50];
//  SFM  end    08/13/2014
//
//  SFM  begin  06/22/2014; for M.Grecu (unknown justification)
  float z13sfc, z13sfcA, z35sfc, z35sfcA;
  float pia13tot=0, pia35tot=0;
//begin  MG 10/29/15
  list<long double> _13,_35;
//end    MG 10/29/15
  nstd=nstdA;
//  SFM  begin  07/01/2014; for M.Grecu  random sequences
  /*if(*wfractPix>10.)
    nstd=0.125;
  */
  //nstd=0;
  nSub=8;
  //  nstd*=1;
//  SFM  end    06/22/2014
  for(i=0;i<nSub;i++)
    {
//  SFM  begin  07/01/2014; for M.Grecu  random sequences
      xs[i]=exp(nstd*xs[i]); 
//  SFM  end    07/01/2014
      sxs+=xs[i];
    }
  sxs/=nSub;
  for(i=0;i<nSub;i++)
    {
      xs[i]/=sxs;
      xs[i]=10.*log10(xs[i]);
      // printf("%g \n",xs[i]);
    }
  //  exit(0);
  for(i=0;i<ngates;i++)
    {
      log10dN[i]=-99.9;
      d0[i]=0;
      rrate[i]=0;
      z13obsP[i]=-99;
    }
  
  for(i=0;i<5;i++)
    {
      if(i<5)
	node[i]=nodes[i]+1;
      else
	node[i]=nodes[i];
      //printf("%3i ",node[i]);
      if(node[i]<1) node[i]=1;
      if(node[i]>100) node[i]>=100;
    } 
  //printf("\n");

  for(i=0;i<nNodes-1;i++)
    {
      fi=0.;
      fi1=0.;
      if(i==iNode)
	fi=delta;
      if(i+1==iNode)
	fi1=delta;
      for(iLev=nodeP[i];iLev<=min(87,nodeP[i+1]);iLev++)
	//	if(iLev<87)
	{
	  f=(iLev-nodeP[i]+0.)/(nodeP[i+1]-nodeP[i]+0.01);
	  if(f>1) f=1;
	  log10dNPi=(1-f)*(log10dNP[i]+fi)+f*(log10dNP[i+1]+fi1);

// **************** crash site
//printf("   A fModelFortran: %i %i %i %12.6f  \n",nbins,ngates,nNodes,delta);
//printf("   B fModelFortran: %i %i %i %i      \n",i,nodeP[i],nodeP[i+1],iLev);
//printf("   C fModelFortran: %i %12.6f %12.6f \n",iNode,f,log10dNPi);
	  log10dN[iLev]=(log10dNPi);
// **************** crash site

	}
    }
  
  *pia35M=0;
  *pia13M=0;

  for(i=0;i<ngates;i++)
    imuv[i]=imu;
  //  for(i=0;i<ngates;i++)
  //  hh[i]=(nbins-i)*dr*cos(*localZAngle/180.*3.1415);
  //printf("nbins=%i \n",nbins);
  if(node[4]>nbins) node[4]=nbins;
  //printf("node0=%i node1=%i \n",node[0],node[4]);


  for(i=node[0]-1;i<node[4];i++)
    {

      //hh[i]=(nbins-i)*dr;
      z13[i]=z13obs[i];
      z13obsP[i]=z13obs[i];
      // printf("%6.2f \n",z13obsP[i]);
    }
  int j,ifreq;
  for(j=0;j<88;j++)
    for(i=0;i<nSub;i++)
      {
	rrateSub[j][i]=-99.;
	d0Sub[j][i]=-99.;
	pwcSub[j][i]=-99;
	z35modSub[j][i]=-99;
	z35Sub[j][i]=-99;
	z13Sub[j][i]=-99;
      }

  for(j=0;j<88;j++)
    for(ifreq=0;ifreq<8;ifreq++)
      {
	for(i=0;i<nSub;i++)
	  {
	    kextSub[j][ifreq][i]=-99.;
	    salbSub[j][ifreq][i]=-99.;
	    asymSub[j][ifreq][i]=-99.;
	  }
	
	kextAvg[j][ifreq]=0.;
	salbAvg[j][ifreq]=0.;
	asymAvg[j][ifreq]=0.;
      }

  
//  SFM  begin  06/22/2014; for M.Grecu  (unknown justification)
  z13sfc=0;
  z13sfcA=0;
  z35sfc=0;
  z35sfcA=0;
//  SFM  end    07/29/2014; for M.Grecu  to eliminate NANs
  int ithresh=0;
//  SFM  end  07/29/2014
//  SFM  end  06/22/2014

  for(i=0;i<nSub;i++)
    {
      for(j=node[0]-1;j<node[4];j++)
	{
//  SFM  begin  06/22/2014; for M.Grecu  (unknown justification)
//  SFM  begin  07/01/2014; for M.Grecu  random sequences
//  SFM  begin  08/13/2014; for M.Grecu   "long orbit" correction
	  z13obsP[j]= z13obs[j]+xs[i];
//  SFM  end    08/13/2014
//  SFM  end    07/01/2014
	  if(z13obsP[j]>50)
	    z13obsP[j]=50;
//  SFM  end    06/22/2014
	}
      fhb11_(z13,z35,z13obsP,
	     pia13M,pia35M,&ic,&jc,z35mod,pwc,log10dN, 
	     &dr,node,&isurf,imuv,
	     &ngates,&nmfreq,hh,&itype,kext,salb,asym,
	     rrate,d0,hfreez, pia13srt,imemb);
//  SFM  begin  07/29/2014; for M.Grecu  to eliminate NANs
      if(isnan(*pia35M))
	{
	  printf("%g \n",*pia35M);
	  printf("fmodel");
	  exit(0);
	}
      if(z13obs[node[4]-1]>12)
	{
	  
	  z13sfc=z13sfc+pow(10,0.1*z13[node[4]-1]);
	  z13sfcA=z13sfcA+pow(10,0.1*z13[node[4]-1]-0.1*(*pia13M));
	  z35sfc=z35sfc+pow(10,0.1*z35[node[4]-1]);
	  z35sfcA=z35sfcA+pow(10,0.1*z35[node[4]-1]-0.1*(*pia35M));
	  if(z35sfcA<1)
	    {
	      // printf("%g %g \n",z35sfcA,*pia35M);
	      // printf("%g %g \n",z35[node[4]-1],0.1*(*pia35M));
	    }
	  ithresh++;
	}
//  SFM  end    07/29/2014; for M.Grecu  to eliminate NANs
      ic=0;
      pia13tot+=*pia13M;
      pia35tot+=*pia35M;
//begin  MG 10/29/15
      _13.push_back(pow(10.,-0.1*(*pia13M)));
      _35.push_back(pow(10.,-0.1*(*pia35M)));
//end    MG 10/29/15

      for(ifreq=0;ifreq<8;ifreq++)
	for(j=0;j<88;j++)
	  {
	    if(isnan(asym[ic]))
	      asym[ic]=0;
	    kextSub[j][ifreq][i]=kext[ic];
	    salbSub[j][ifreq][i]=salb[ic];
	    asymSub[j][ifreq][i]=asym[ic];
	    if(isnan(asym[ic]))
	      {
		printf("%g %g %g %i %i\n",asym[ic],salb[ic],kext[ic],ifreq,j);
		printf("fmodel");
		exit(0);
	      }
	    ic++;
	  }

//  SFM  begin  06/22/2014; for M.Grecu  (unknown justification)
      float tau=0;
      for(j=0;j<88;j++)
	{
	  rrateSub[j][i]=rrate[j];
	  d0Sub[j][i]=d0[j];
	  pwcSub[j][i]=pwc[j];
	  if(salbSub[j][3][i]>0 && 
	     kextSub[j][3][i]>0)
	    tau+=salbSub[j][3][i]* 
	      kextSub[j][3][i]*0.25;
	  z35mod[j]+=tau*7;
	  z35modSub[j][i]=z35mod[j];
	  z35Sub[j][i]=z35[j];
	  z13Sub[j][i]=z13[j];
//  SFM  end  06/22/2014
	}
    }

  for(j=0;j<88;j++)
    {
      countSub[j]=0;
      rrate[j]=0.;
      pwc[j]=0.;
      d0[j]=0.;
      z35[j]=0;
      z35mod[j]=0;
      z13[j]=0;
    }
  int iSub;  
//  SFM  start  06/22/2014; for M.Grecu; unknown justification)
//  SFM  begin  08/13/2014; for M.Grecu   "long orbit" correction
  for(i=node[0]-1;i<node[4];i++)
    if(z13obs[i]>12)
//  SFM  end    08/13/2014
      {
	for(iSub=0;iSub<nSub;iSub++)
	  if(kextSub[i][1][iSub]>0)
	    {
	      z35mod[i]+=pow(10.,0.1*z35modSub[i][iSub]);
	      z13[i]+=pow(10.,0.1*z13Sub[i][iSub]);
	      z35[i]+=pow(10.,0.1*z35Sub[i][iSub]);
	      rrate[i]+=rrateSub[i][iSub];
//  SFM  begin  07/29/2014; for M.Grecu  to eliminate NAN data
	      pwc[i]+=pwcSub[i][iSub];
//  SFM  end    07/29/2014
	      d0[i]+=d0Sub[i][iSub];
	      countSub[i]++;
	    }
	//*pia13M=pia13tot/nSub;
	//*pia35M=pia35tot/nSub;
      }
//  SFM  begin  07/29/2014; for M.Grecu  to eliminate NAN data
  if(ithresh>0)
    {
//  SFM  begin  08/13/2014; for M.Grecu   "long orbit" correction
      if(z13sfcA>0.01 && z13sfc>0.01)
//  SFM  end    08/13/2014
	*pia13M=log10(z13sfc/z13sfcA)*10.;
      else
	*pia13M=pia13tot/nSub;
      if(z35sfcA>0.01 && z35sfc>0.01)
	*pia35M=log10(z35sfc/z35sfcA)*10.;
      else
	*pia35M=pia35tot/nSub;
      if(isnan(*pia35M))
	{
	  //printf("%g %g \n",z35sfc,z35sfcA);
	  //exit(0);
	}
      if(isinf(*pia35M))
	{
	  //printf("%g %g \n",z35sfc,z35sfcA);
	  // printf("%g %g %g %g\n",log10(z35[node[4]-1]/countSub[node[4]-1]),
	  //log10(z35mod[node[4]-1]/countSub[node[4]-1]),
	  //		 pia35tot/nSub,*pia35M);
	  //exit(0);
	}
    }
  else
    {
      *pia13M=pia13tot/nSub;
      *pia35M=pia35tot/nSub;
    }
//  SFM  end    06/22/2014
//  SFM  end    07/29/2014

//  SFM  begin  08/13/2014; for M.Grecu   "long orbit" correction
  for(i=node[0]-1;i<node[4];i++)
//  SFM  end    08/13/2014
    {
      if(countSub[i]>0)
	{
	  z35mod[i]=10.*log10(z35mod[i]/countSub[i]);
	  z13[i]=10.*log10(z13[i]/countSub[i]);
	  z35[i]=10.*log10(z35[i]/countSub[i]);
	  rrate[i]/=countSub[i];
//  SFM  begin  07/29/2014; for M.Grecu  to eliminate NAN data
	  pwc[i]=(pwc[i]/countSub[i]);
//  SFM  end    07/29/2014
	  d0[i]/=countSub[i];
//  SFM  start  06/22/2014; for M.Grecu; unknown justification)
	  if(z13[i]>155)
//  SFM  end    06/22/2014
	    {
	    printf(" fModelFortran-a %6.2f %6.2f \n",z13[i],z35[i],rrate[i]);
	    for(iSub=0;iSub<nSub;iSub++)
	      printf("z13=%6.2f z35=%6.2f ",z13Sub[i][iSub],z35Sub[i][iSub]);
            printf("\n");
	    }
	}
//  SFM  begin  07/29/2014; for M.Grecu  to eliminate NAN data
      else
	{
	  z35mod[i]=-99;
	  z13[i]=-99;
	  z35[i]=-99;
	  rrate[i]=-99;
	  pwc[i]=-99;
//  SFM  end    07/29/2014

	}

    }
     
//  SFM  deletion  08/13/2014; for M.Grecu, "long orbit" correction

  for(i=0;i<nNodes;i++)
    {
      //log10dNP[i]=log10dN[nodeP[i]];
    }

//  SFM  begin  07/29/2014; for M.Grecu  to eliminate NAN data
  int count;
  for(j=0;j<88;j++)
    {
      for(ifreq=0;ifreq<8;ifreq++)
	{
	  count=0;
	  kextAvg[j][ifreq]=0;
	  salbAvg[j][ifreq]=0;
	  asymAvg[j][ifreq]=0;
	  for(iSub=0;iSub<nSub;iSub++)
	    if(kextSub[j][ifreq][iSub]>0)
	      {
		kextAvg[j][ifreq]+=kextSub[j][ifreq][iSub];
		salbAvg[j][ifreq]+=salbSub[j][ifreq][iSub];
		asymAvg[j][ifreq]+=asymSub[j][ifreq][iSub];
		count++;
	      }
	  if(count>0)
	    {
	      kextAvg[j][ifreq]/=count;
	      salbAvg[j][ifreq]/=count;
	      asymAvg[j][ifreq]/=count;
	    }
	  else
	    {
	      kextAvg[j][ifreq]=-99;
	      salbAvg[j][ifreq]=-99;
	      asymAvg[j][ifreq]=-99;
	    }
	}
 
    }
  ic=0;

  for(ifreq=0;ifreq<8;ifreq++)
    for(j=0;j<88;j++)
      {
	kext[ic]=kextAvg[j][ifreq];
	salb[ic]=salbAvg[j][ifreq];
	asym[ic]=asymAvg[j][ifreq];
	ic++;
      }
//begin  MG 10/29/15
  *pia13M=-10.*log10(_(_13)/nSub);
  *pia35M=-10.*log10(_(_35)/nSub);
//end    MG 10/29/15

//  SFM  end    07/29/2014
//  SFM  deletion  08/13/2014; for M.Grecu, "long orbit" correction
if(isinf(*pia35M))
    {
      printf("%g %g %g\n",*pia13M,*pia35M,pia35tot/nSub);
      *pia35M=pia35tot/nSub;
      //      exit(0);
    }
  
  for(j=0;j<88;j++)
    {
      if(log10dN[j]>4.5) log10dN[j]=4.5; 
      if(log10dN[j]<-4.5) log10dN[j]=-4.5; 
    }

  return;    
 
}

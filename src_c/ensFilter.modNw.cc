
 //  SFM  04/06/2013  Minor changes to calls involving "int" for compile errors
 //                     in modules setNodeP and setrandclass
 // SFM   07/26/2013  Code mods for M.Grecu

#include "typedef.h"

#ifdef GFOR 
#define ran1() __ran_mod_MOD_ran1()
//#define ran1() __ran_mod__ran1()
#endif

#ifdef IFORT 
#define ran1() ran_mod_mp_ran1_()
#endif

#ifdef IFORT 

#define iMemb nbinmod_mp_imemb_
#endif

extern int nbinmod_mp_nbin_;
#define nbins nbinmod_mp_nbin_

float Dm[60],NwSt[60],NwCv[60];


void runEns(radarDataType   *radarData, 
	    stormStructType *stormStruct,
	    retParamType    *retParam,  
	    int *nmu, radarRetType    *radarRet, 
	    float *xscalev,
	    float *randemiss, float *localZAngle, 
	    float *wfractPix, long *ichunk, 
	    float **xs, float **logdNw,float **kext, float **asym, 
	    float **salb,float *rsurf, float **dz, int *imuv, int *nodeP, 
	    int nNodes, float nstdA);

extern "C" void readdmnw_(void)
{
  int i;
  FILE *fdm;
  float dN=log10(8.e6);
  fdm=fopen("AncData/NwDm.txt","r");
  for(i=0;i<60;i++)
    {
      fscanf(fdm,"%g %g %g",&Dm[i],&NwCv[i],&NwSt[i]);
      NwCv[i]-=dN;
      NwSt[i]-=dN;
      NwCv[i]+=0.1;
    }
  fclose(fdm);
}

void printrefprof()
{

}

#ifdef GFOR 
#define iMemb __nbinmod_MOD_imemb
//#define ran1() __ran_mod__ran1()
#endif
void filterZ(float *zku, int node5[5])
{
  int i, j;

  float zkuf[188];
  for(i=node5[0];i<node5[4];i++)
    zkuf[i]=zku[i];
  for(i=node5[0]+1;i<node5[4]-1;i++)
    {
      //zkuf[i]=zku[i];
      if(zku[i] <17 || zku[i-1]<17 || zku[i+1]<17 )
	zkuf[i]=-99;
    }
  for(i=node5[0];i<node5[4];i++)
    zku[i]=zkuf[i];
}
extern "C" int nbinmod_mp_imemb_;
extern "C" int __nbinmod_MOD_imemb;
extern "C" float ran_mod_mp_normal2_(float *nm, float *nstd); 
extern "C" float __ran_mod_MOD_normal2(float *nm, float *nstd); 
extern "C" double ran_mod_mp_ran1_();  // returns a uniformly distributed 
                                       // random number between
                                       // 0. and 1.0
extern "C" double __ran_mod__ran1();   // returns a uniformly distributed 
                                       // random number between
                                       // 0. and 1.0x

extern "C" double __ran_mod_MOD_ran1();// returns a uniformly distributed 
                                       // random number between
                                       // 0. and 1.0

void setNodeP(int **nodeP,int *nNodes,stormStructType *stormStruct);

void gaussNewton(float pia13obs, float pia35obs,
		 float *z35obs,
		 float **logdNP, float *pia13ModEns, 
		 float *pia35ModEns, float *z35ModEns, float *d0Ens, float *logNwEns,
		 int *nodes, int nNodes, int ngates, int nMemb, float sigPIA35,  float sigPIA13, 
		 int it, float *rsurf, float *rms1, float relPIA13);

void gaussNewtonDepAgg(float pia13obs, float pia35obs,
		       float *z13obs, float *rrate, float *d0,
		       float **logdNP, float *pia13ModEns, 
		       float *pia35ModEns, float *z35ModEns,
		       int *nodes, int nNodes, int ngates, int nMemb, 
		       float sigPIA35,  float sigPIA13, int it, float *rsurf, float *rms1,
		       float relPia13srt);

void gaussNewtonE(float pia13obs, float pia35obs,
		  float *z35obs,
		  float **logdNP, float *pia13ModEns, 
		  float *pia35ModEns, float *z35ModEns,
		  int *nodes, int nNodes, int ngates, int nMemb, 
		  float sigPIA35,  float sigPIA13,
		  int it, float *rms2);

void enkf(float pia13obs, float pia35obs,
	  float **logdNP, float *pia13Mens, 
	  float *pia35Mens, float *z35sim,
	  int nNodes, int nMemb, float sigPIA13, int it);



void enkf35(float pia13obs, float pia35obs,
	    float **logdNP, float *pia13Mens, 
	    float *pia35Mens, float *z35sim,
	    int nNodes, int nMemb, float sigPIA35, int it);

void printRet(radarDataType *radarData,stormStructType *stormStruct, radarRetType *radarRet);
void fModelFortranD0(float *z13obs, int nodes[5], int isurf, int imu,
		     float *d0n, int nNodes, 
		     float *pia35M, float *pia13M,
		     float *z35mod, float *pwc, float dr, int ic, int jc, float *hh,
		     float delta, int iNode, int nmfreq,
		     float *salb, float *kext,float *asym, int rainType,
		     int ngates,float *rrate,float *d0,float *logN,
		     float *z13, float *z35,
		     int *imuv, float *hfreez, float *dz, float *pia13srt, 
		     float *relPia13srt);


//  SFM  begin  06/16/2014; for M. Grecu, multiple scattering
extern "C" FILE *fout;
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
//  SFM  begin  06/22/2014; for M. Grecu, (unknown justification)
//  SFM  begin  07/01/2014; for M. Grecu, random sequences
		   int *imemb,float *localZAngle, float *wfract,
		   float *xs, long *ichunk, float nstdA);

void fModelFortranCv(float *z13obs, float *z35obs, 
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
		     //  SFM  begin  06/22/2014; for M. Grecu, (unknown justification)
		     //  SFM  begin  07/01/2014; for M. Grecu, random sequences
		     int *imemb,float *localZAngle, float *wfract,
		     float *xs, long *ichunk, int *i0, int *j0);



//  SFM  end    07/01/2014
//  SFM  end    06/22/2014
//  SFM  end    06/16/2014
/**********************************************************************************************

**** on input    ******************************************************************************
float *z13obs    : vector, observed Ku-band reflectivity
float *hh        : vector, height associated with the reflectivity observations
int ngates       : scalar, # of gates in the profile
int rainType     : scalar, 1 stratiform, 2 convective
float dr         : scalar, radar gatesize
int nodes[5]     : vector, the 5 nodes defining the storm structure
int nNodes       : scalar, # nodes defining the logN profile
int *nodeP       : vector, contains the locations of the nNodes nodes
float *logdNP    : vector, contains the values of logdN of the nNodes nodes
int ic           : scalar, index of RH profile (from 1 to nc)
                 : nc is the number of possible RH classes see cloud.f90
int jc           : scalar, index of cloud profile (from 1 to nc)
int imu          : scalar, mu index of look up table (from 1 to nmu)
int iz           : integer, generally 0, but if an ensemble of solution is derived can be used to
                 : differentiate among ensemble members

float delta      : scalar, used to calculate the jacobian
int iNode        : scalar, index of the node where logdN is perturbed
int nmfreq       : integer, the number of microwave frequencies (read in readTables and 
                   set in mie2)


**** on output   *****************************************************************************


float *pia35M    : reference to scalar, model pia at Ka-band
float *pia13M    : reference to scalar, model pia at Ku-band
float *z35mod    : vector, model Ka-band reflectivities  
 
 

float *pwc       : vector, precipitation water content (g/m3)

float *salb      : vector, fortran 2-d array, returns the scattering albedo at nmfreq mw freq
                 : these frequencies are set in mie2 subroutine in file Mie/mie3.f90 and 
                 : read by readTables subroutine in readTables.f90
                 : they can be accessed in variable mfreq in f90 module microwFreq in readTables.f90
float *kext      : vector, fortran 2-d array, returns the extinction coefficient (1/km)
                 : at nmfreq mw freq
float *asym      : vector, fortran 2-d array, returns the asymmetry fact

float *rrate     : vector, returns retrieved rain rate (mm/h)
float *d0        : vector, currently empty 
float *logN      : vector, currently empty
********************************************************************************************/  

void calctbc(float *sfc_wind,float *umu,float *kext,
	     float *salb,float *asym,int *node,int *ic,int *jc,
	     int *ngates,int *nmfreqm,float *hh,float *tb,float *emtb, float *emis,
	     float *tpw, float *hfreez,
	     float *randemiss, int *imemb, int *nmemb, int *i0, int *j0,
	     float *w1rand, float *jrand);

void avgScProp(float *sfc_wind,float *umu,float *kext,
	       float *salb,float *asym, int *node,int *ic,int *jc,
	       int *ngates,int *nmfreqm,float *hh,float *tb, float *tpw, float *hfreez,
	       float *randemiss, int *imemb,int *nmemb);

extern "C" void calctb_(float *sfc_wind,float *umu,float *kext,
			float *salb,float *asym,int *node,int *ic,int *jc,
			int *ngates,int *nmfreqm,float *hh,float *tb, 
			float *emtb,float *emis,
			float *tpw, float *hfreez, 
			float *randemiss, int *imemb,int *nmemb,
			int *i0, int *j0, float *w1ran, float *jran);


extern "C" int multiscatter_(int *nrange, float *extFort, 
			     float *ext2bscatt, float *salbFort, float *gFort,
			     float *bscatFort, int noMS);

extern "C" void avgscprop_(float *sfc_wind,float *umu,float *kext,
			   float *salb,float *asym, int *node,int *ic,int *jc,
			   int *ngates,int *nmfreqm,float *hh,float *tb, float *tpw, float *hfreez,
			   float *randemiss, int *imemb,int *nmemb);

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






void setNodeP(int **nodeP,int *nNodes,stormStructType *stormStruct)
{
  extern int nbinmod_mp_n9_[9];
  int i;
  *nNodes=9;
  (*nodeP)=(int *) malloc(sizeof(int)*(*nNodes));
  (*nodeP)[0]=stormStruct->nodes[0];                                   // set nodeP
  (*nodeP)[1]=int(0.75*stormStruct->nodes[0]+0.25*stormStruct->nodes[1]);
  (*nodeP)[2]=int(0.5*stormStruct->nodes[0]+0.5*stormStruct->nodes[1]);
  (*nodeP)[3]=int(0.25*stormStruct->nodes[0]+0.75*stormStruct->nodes[1]);
  (*nodeP)[4]=stormStruct->nodes[1];
  (*nodeP)[5]=stormStruct->nodes[2];
  (*nodeP)[6]=stormStruct->nodes[3];
  (*nodeP)[*nNodes-1]=stormStruct->nodes[4];
  (*nodeP)[*nNodes-2]=int(0.5*stormStruct->nodes[4]+0.5*stormStruct->nodes[3]);
  for(i=0;i<9;i++)
    nbinmod_mp_n9_[i]=(*nodeP)[i];
}

extern "C" void setrandclass_(radarRetType *radarRet, int *nmu)
{

  int i ;
  double imu, icc, jcc ;

  for(i=0;i<radarRet->nMemb;i++)
    {
      imu=(ran1()*(*nmu-0.5))+1;
      if(imu>*nmu) imu=*nmu;
      icc= ran1()*(49-0.1)+1;
      if(icc>49) icc=49;
      jcc= ran1()*49+1;
      if(jcc>49) jcc=49;

      radarRet->imu[i]=int(imu);
      radarRet->icc[i]=int(icc);
      radarRet->jcc[i]=int(jcc);
    }
}




/*
  subroutine calctb simulates brightness temperatures 
  at frequencies at radarRet.nmfreq
  the subroutine resides in rterain.f90
*/
/* calctb_(&radarRet->sfc_wind[0],&umu,kext[0],
   salb[0],asym[0],node,&radarRet->icc[0],
   &radarRet->jcc[0],&radarData->ngates,&radarRet->nmfreq,
   radarData->hh,radarRet->tb);
*/
extern "C" float getextka_(int *i);

extern "C" void calctb2_(float *sfc_wind,float *umu,int *node,int *ic,int *jc,
			int *ngates,int *nmfreqm,float *hh,float *tb);

void calctbc(float *sfc_wind,float *umu,float *kext,
	     float *salb,float *asym,int *node,int *ic,int *jc,
	     int *ngates,int *nmfreqm,float *hh,float *tb,float *emtb, float *emis,
	     float *tpw, float *hfreez,
	     float *randemiss, int *imemb, int *nmemb, int *i0, int *j0,
	     float *w1rand, float *jrand)
{
  calctb_(sfc_wind,umu,kext,salb,asym,node,ic,jc, ngates,nmfreqm,hh,tb,
	  emtb,emis,tpw, hfreez,
	  randemiss,imemb,nmemb,i0,j0,w1rand,jrand);
 
}

void avgScProp(float *sfc_wind,float *umu,float *kext,
	     float *salb,float *asym,int *node,int *ic,int *jc,
	     int *ngates,int *nmfreqm,float *hh,float *tb, float *tpw, float *hfreez,
	     float *randemiss, int *imemb, int *nmemb)
{
  avgscprop_(sfc_wind,umu,kext,salb,asym,node,ic,jc, ngates,nmfreqm,hh,tb,tpw, hfreez,
	     randemiss,imemb,nmemb);
 
}


extern "C" void ensradretstcvku_( radarDataType   *radarData, 
				  stormStructType *stormStruct,
				  retParamType    *retParam,  
				  int *nmu, radarRetType    *radarRet, 
				  int *iflag, float *rms1,
				  float *rms2, float *dn, int *iit, 
				  float *xscalev,
//  SFM  begin  06/22/2104; for M.Grecu (unknown justification)				  
				  float *randemiss, float *localZAngle, 
				  float *wfractPix, long *ichunk, 
				  int *i0, int *j0, float *dZms, int *msFlag) //!! MS&WSO addition Feb 11, 2017
//  SFM  end    06/22/2104; for M.Grecu (unknown justification				  
//  SFM  general 08/01/2014  ichunk added for diagnostics  
{
 
  
  float **kext,    // pointer to vector, fortran 2-d array, 
                   // contains the extinction coefficient (1/km) profile
    **salb,        // pointer to vector, fortran 2-d array, 
                   // contains the scattering albedo profile  
    **asym;        // pointer to vector, fortran 2-d array, 
                   // contains the asymmetry factor profile   
 
 

  int nNodes;      // integer, number of points that define the parametrized logNw profile
  float **logdNw;  // pointer to vector, this is a vector of vectors 
                   // logdNw[][*] is a nNodes-long vector containing log10(Nw/NwRef) 
                   // currently first dimension is 1 and the first index is always 0,
                   // but the first dimension can be used in ensemble retrievals
  int *nodeP;      // vector, nodeP[i], i<nNodes, is the node associated with logdNw[i]
                   // the logdNw profile is derived by interpolation between these nodes
 
  int iNode;       // local variable


  float *z35mod, pia35M, pia13M;  // simulated observations
  int nl;                         // counts the number of DF observations
  float delta;                    // used to perturb the Nw profile
  int k, icount, iLev,i;          // local variables
  int j, ic;
  float pia35M0,                  // unperturbed simulated PIA at Ka-band
    pia13M0;                      // unperturbed simulated PIA at Ka-band
  int nMemb=radarRet->nMemb;
  float *z13, *z35, **dz, *rsurf;
  int *imuv;
  int imemb;
//  SFM  begin  07/01/2014; for M.Grecu, random sequences
  float umu, **xs;
//  SFM  end    07/01/2014
  int  node[5];
  float logdNw_pert;

  umu=cos(53/180.0*3.141592654);
  for(i=0;i<5;i++) {
    node[i]=stormStruct->nodes[i]+1;  
  }

  rsurf=(float *) malloc(sizeof(float)*nMemb);
  kext=(float**)malloc(sizeof(float*)*nMemb);
  salb=(float**)malloc(sizeof(float*)*nMemb);
  asym=(float**)malloc(sizeof(float*)*nMemb);
//  SFM  begin  07/01/2014; for M.Grecu, random sequences
  xs=(float**)malloc(sizeof(float*)*nMemb);
//  SFM  end    07/01/2014

  z35=(float*) malloc(sizeof(float)*radarData->ngates*nMemb);
  z13=(float*) malloc(sizeof(float)*radarData->ngates*nMemb);
  dz=(float**) malloc(sizeof(float*)*nMemb);
  imuv=(int*) malloc(sizeof(int)*radarData->ngates*nMemb);

  for(i=0;i<nMemb;i++)
    {
      kext[i]=(float*)malloc(sizeof(float)*radarData->ngates*radarRet->nmfreq);
      salb[i]=(float*)malloc(sizeof(float)*radarData->ngates*radarRet->nmfreq);
      asym[i]=(float*)malloc(sizeof(float)*radarData->ngates*radarRet->nmfreq);
      dz[i]=(float*) malloc(sizeof(float)*radarData->ngates);
//  SFM  begin  07/01/2014; for M.Grecu, random sequences
      xs[i]=(float*) malloc(sizeof(float)*50);
      for(j=0;j<50;j++)
	{
	  float nm=0, nstd=1.;
	  xs[i][j]=ran_mod_mp_normal2_(&nm,&nstd);
	}
//  SFM  end    07/01/2014
    }

  z35mod=(float *) malloc(sizeof(float)*radarData->ngates);
  pia13M=0;
  nl=0;

  /*
  for(iLev=0;iLev<5;iLev++)
    fprintf(fout,"%i,",stormStruct->nodes[iLev]);
  fprintf(fout,"node5\n");
  for(iLev=stormStruct->nodes[0];iLev<stormStruct->nodes[4];iLev++)
    fprintf(fout,"%g,",radarData->z13obs[iLev]);
  fprintf(fout,"z13obs\n");
  for(iLev=stormStruct->nodes[0];iLev<stormStruct->nodes[4];iLev++)
    fprintf(fout,"%g,",radarData->z35obs[iLev]);
  fprintf(fout,"z35obs\n");
  fprintf(fout,"%g,%g,pia13srt\n",radarData->pia13srt,radarData->pia35srt);
  fprintf(fout,"%i,itype\n",stormStruct->rainType);
  fprintf(fout,"%g,hfreez\n",radarData->hfreez);*/
  for(iLev=0;iLev<=stormStruct->nodes[4];iLev++)
    {
      radarRet->log10dNw[iLev]=0.;
      for(j=0;j<nMemb;j++)
	dz[j][iLev]=(ran1()-0.5);
      z35mod[iLev]=-99;
      if(radarData->z13obs[iLev]>retParam->z13thresh)
        nl++;  

    }

  // allocates memory
  //stormStruct->nodes[0]=0;
  int iset=0;
  //stormStruct->nodes[0]=0;
  for (i=0;i<-radarData->ngates-1;i++) 
    {
      radarRet->rrate[i]=-99.9;
      radarRet->d0[i]=-99.9;
	if(radarData->z13obs[i]>retParam->z13thresh && 
	   radarData->z13obs[i+1]>retParam->z13thresh && iset==0)
	{
	  stormStruct->nodes[0]=i;
	  iset=1;
	}
	//printf("%3i %6.2f \n",i-8,radarData->z13obs[i]);
    }
  //  printf("PIASRT= %g %g\n", radarData->pia13srt, radarData->relPia13srt);
  

  if(stormStruct->nodes[0]>=stormStruct->nodes[1]) stormStruct->nodes[0]=stormStruct->nodes[1]-1;
  
  setNodeP(&nodeP,&nNodes,stormStruct);
  logdNw=(float**)malloc(sizeof(float*)*nMemb);
  float nm=0., nstd=1;
  float xscale;


  for(i=0;i<radarRet->nMemb;i++)
    {
      logdNw[i]=(float *)malloc(sizeof(float)*nNodes);
      for(j=0;j<nNodes;j++)
	{
	  logdNw[i][j]=radarRet->logdNw[i*nNodes+j];
	}
      for(j=0;j<nNodes;j++)
	{
	  float nm=0, nstd=1.;
	  logdNw[i][j]+=0.25*ran_mod_mp_normal2_(&nm,&nstd);
	  //logdNw[i][j]+=0.75*__ran_mod_MOD_normal2(&nm,&nstd);
	}

      for(j=0;j<nNodes;j++)
	{
	  if(logdNw[i][j]>4.)
	    logdNw[i][j]=4.;
	  if(logdNw[i][j]<-4.)
	    logdNw[i][j]=-4.;
	}
    }

  for(i=0;i<radarRet->nMemb;i++)
    if(stormStruct->nodes[2]>stormStruct->nodes[0] &&
       stormStruct->nodes[2]>stormStruct->nodes[1])
      {
	float dn1=(-stormStruct->nodes[0]+stormStruct->nodes[2])*0.025*0.25;
	float dn2=(-stormStruct->nodes[1]+stormStruct->nodes[2])*0.025*0.25;
	logdNw[i][0]+=dn1;
	logdNw[i][1]+=0.25*dn2+0.75*dn1;
	logdNw[i][2]+=0.5*(dn1+dn2);
	logdNw[i][3]+=0.75*dn2+0.25*dn1;
	logdNw[i][4]+=dn2;
      }
	

  if(stormStruct->rainType!=2)// and radarData->hfreez>=2)
    {
      //printf("%i \n",stormStruct->rainType);
      for(j=0;j<9;j++)
	{
	 float logdnm=0;
	 for(i=0;i<radarRet->nMemb;i++)
	   logdnm+=logdNw[i][j];
	 logdnm/=radarRet->nMemb;
	 //printf("%g \n",logdnm);
	 for(i=0;i<radarRet->nMemb;i++)
	   if(*wfractPix<10)
	     logdNw[i][j]=logdnm+0.75*(logdNw[i][j]-logdnm);
	   else
	     logdNw[i][j]=logdnm+0.75*(logdNw[i][j]-logdnm);
       }
    }

  if(*wfractPix>90)
    {
      for(j=0;j<9;j++)
	for(i=0;i<radarRet->nMemb;i++)
	  logdNw[i][j]+=0.1;
    }
  else
    if(*wfractPix<10 and radarData->hfreez>2.5 and
       stormStruct->rainType==2)
      for(j=0;j<9;j++)
	for(i=0;i<radarRet->nMemb;i++)
	  logdNw[i][j]-=0.1;

  if(stormStruct->rainType==2 and radarData->hfreez>=-2)
    for(i=0;i<radarRet->nMemb;i++)
      for(j=0;j<9;j++)
	logdNw[i][j]=0.85*logdNw[i][j]+0.15*
	  (0.25*ran_mod_mp_normal2_(&nm,&nstd));

   float pia13EnsM=0;
  // printf("Nmemb=%i %i\n",nMemb,stormStruct->iSurf);

//here, modify logdNw to determine impact
  logdNw_pert = 0.2;
  for(i=0;i<radarRet->nMemb;i++)
    {
      for(j=0;j<nNodes;j++)
        {
          logdNw[i][j]+=logdNw_pert;
        }
      for(j=0;j<nNodes;j++)
       {
        if(logdNw[i][j]>4.)
           logdNw[i][j]=4.;
        if(logdNw[i][j]<-4.)
          logdNw[i][j]=-4.;
       }
     }
//end here

  for(i=0;i<radarData->ngates;i++)
    radarData->hh[i]=(nbins-i)*radarData->dr*cos(*localZAngle/180.*3.1415);

  float nstdA=0.125;
  runEns(radarData, stormStruct,retParam, nmu,radarRet, 
	 xscalev, randemiss, localZAngle, wfractPix, ichunk, 
	 xs, logdNw,kext,asym,salb,rsurf,dz, imuv, nodeP, 
	 nNodes,nstdA);

  float dm=0,rsfcMean;
  for( imemb=0;imemb<radarRet->nMemb;imemb++)
    {
      dm+=radarRet->d0[imemb*radarData->ngates+stormStruct->nodes[4]];
      rsfcMean+=rsurf[imemb];
    }
  dm/=radarRet->nMemb;
  rsfcMean/=radarRet->nMemb;
  
  if(stormStruct->rainType==1 && dm>0.5 && *wfractPix<10)
    {
      //printf("%g %g %g %i %i\n",dm,rsfcMean,*cBEst,*i0,*j0);
      if(rsfcMean>5)
	{
	  nstdA=0.125+(rsfcMean-5)/5.*0.125;
	  if(nstdA>0.5)
	    nstdA=0.5;
	}
		     
      for(i=0;i<radarRet->nMemb;i++)
	for(j=0;j<9;j++)
	  {
	    if(stormStruct->rainType==1 && dm>0.5)
	      logdNw[i][j]-=0.3*(dm-1.5);
	    int i0dm=(int)((dm-0.5)/0.04);
	    if(i0dm<0)
	      i0dm=0;
	    if(i0dm>59)
	      i0dm=59;
	    if(stormStruct->rainType==1 && dm>0.5)
	      logdNw[i][j]+=0.6*(NwSt[i0dm]-logdNw[i][j]);
	    //if(stormStruct->rainType==1 && dm>0.5)
	    //	logdNw[i][j]-=0.5;
	  }
      runEns(radarData, stormStruct,retParam, nmu,radarRet, 
	     xscalev, randemiss, localZAngle, wfractPix, ichunk, 
	     xs, logdNw,kext,asym,salb,rsurf,dz, imuv, nodeP, 
	     nNodes,nstdA);
    }

  if(stormStruct->rainType==2 && dm>0.5 && *wfractPix<10)
    {
      float cdN=0.56826200674023877-0.6907993*dm;
      //printf("%g %g %g %i %i\n",dm,rsfcMean,*cBEst,*i0,*j0);
      for(i=0;i<radarRet->nMemb;i++)
	for(j=0;j<9;j++)
	  {
	    int i0dm=(int)((dm-0.5)/0.04);
	    if(i0dm<0)
	      i0dm=0;
	    if(i0dm>59)
	      i0dm=59;
	    if(stormStruct->rainType==2 && dm>0.5)
	      logdNw[i][j]+=0.2*(cdN-logdNw[i][j]);
	    //if(stormStruct->rainType==1 && dm>0.5)
	    //	logdNw[i][j]-=0.5;
	    if(rsfcMean>15)
	      {
		nstdA=0.125+(rsfcMean-15)/10.*0.125;
		if(nstdA>0.35)
		  nstdA=0.35;
	      }

	  }
      runEns(radarData, stormStruct,retParam, nmu,radarRet, 
	     xscalev, randemiss, localZAngle, wfractPix, ichunk, 
	     xs, logdNw,kext,asym,salb,rsurf,dz, imuv, nodeP, 
	     nNodes,nstdA);
    }

  float genv_mp_wvext;
  float lamb=0.00857;
  float kextFort[88],z35med[88], ext2bscattFort[88], salbFort[88], gFort[88]; 
  float pi4=97.409091034;
  float lamb4=5.39415e-9;

  for(i=0;i<88;i++)
    {
      z35med[i]=0;
      kextFort[i]=0;
      salbFort[i]=0;
      gFort[i]=0;
      genv_mp_wvext=getextka_(&i);
      for( imemb=0;imemb<radarRet->nMemb;imemb++)
	{
	  z35med[i]+=pow(10,0.1*radarRet->z35[imemb*radarData->ngates+i]);
	  if(kext[imemb][nbins*3+i]>0)
	    kextFort[i]+=kext[imemb][nbins*3+i];
	  if(salb[imemb][nbins*3+i]>0)
	    salbFort[i]+=salb[imemb][nbins*3+i];
	  if(asym[imemb][nbins*3+i]>-10)
	    gFort[i]+=asym[imemb][nbins*3+i];
	}
      
      if(z35med[i]>0)
	{
	  z35med[i]=10.*log10(z35med[i]/radarRet->nMemb);
	  
	}
      salbFort[i]/=radarRet->nMemb;
      gFort[i]/=radarRet->nMemb;
      float ff=kextFort[i]/(genv_mp_wvext*radarRet->nMemb+kextFort[i]);
      salbFort[i]*=ff;
      kextFort[i]=(genv_mp_wvext+ kextFort[i]/radarRet->nMemb)*1e-3;
      if(z35med[i]>0)
	{
	  float Z=pow(10.,0.1*z35med[i])*pi4/1e18/lamb4/4*0.93;
	  ext2bscattFort[i]=kextFort[i]/Z;
	  /*printf("%i %g %g %g %g %g \n",i,z35med[i],kextFort[i],
		 salbFort[i],gFort[i],ext2bscattFort[i]);
	  */
	}
      else
	ext2bscattFort[i]=1000;
      if(ext2bscattFort[i]>1000)
	ext2bscattFort[i]=1000;
    }

 int nrange=80;
 float bscatFort[88],bscatFortnoMS[88];
 float dZ=log10(pi4/1e18/lamb4/4*0.93)*10;
 int imsEns=0;
 // goto skipMS;
 
 multiscatter_(&nrange, &kextFort[8], 
	       &ext2bscattFort[8], &salbFort[8], &gFort[8],
	       &bscatFort[8],0);
 multiscatter_(&nrange, &kextFort[8], 
	       &ext2bscattFort[8], &salbFort[8], &gFort[8],
	       &bscatFortnoMS[8],1);


 *dZms=0; //!! MS addition Feb 10, 2017
 *msFlag = 0;  //!!WSO addition Feb 11, 2017
 for(i=8;i<88;i++)
   if(z35med[i]>0)
     {
       float z1=log10(bscatFort[i]+1e-18)*10-dZ;
       float z2=log10(bscatFortnoMS[i]+1e-18)*10-dZ;
       if(z1-z2>0.1)
        {
	       imsEns=1;
          *msFlag = 1; //!!WSO addition Feb 11, 2017
        } 
       if(z1-z2>*dZms) //!! MS addition Feb 10, 2017
	 *dZms=z1-z2;  //!! MS addition Feb 10, 2017
       //if(z1-z2>30)
       // printf("%g %g %g %2i %2i %2i\n",z35med[i],z1,z2,*i0,*j0,i);
       
     }

 float kextFortM[50][88], ext2bscattFortM[50][88], salbFortM[50][88], 
   gFortM[50][88]; 
 
 for(i=0;i<88;i++)
   {
     z35med[i]=0;
     genv_mp_wvext=getextka_(&i);
     for( imemb=0;imemb<radarRet->nMemb;imemb++)
       {
	 z35med[i]=radarRet->z35[imemb*radarData->ngates+i];
	 if(kext[imemb][nbins*3+i]>0)
	   kextFortM[imemb][i]=kext[imemb][nbins*3+i];
	 else
	   kextFortM[imemb][i]=0;
	 if(salb[imemb][nbins*3+i]>0)
	   salbFortM[imemb][i]=salb[imemb][nbins*3+i];
	 else
	   salbFortM[imemb][i]=0;
	 if(asym[imemb][nbins*3+i]>-10)
	   gFortM[imemb][i]=asym[imemb][nbins*3+i];
	 else
	   gFortM[imemb][i]=0;
	 float ff=kextFortM[imemb][i]/(genv_mp_wvext+kextFortM[imemb][i]);
	 salbFortM[imemb][i]*=ff;
	 kextFortM[imemb][i]=(genv_mp_wvext+ kextFortM[imemb][i])*1e-3;
	 if(z35med[i]>0)
	   {
	     float Z=pow(10.,0.1*z35med[i])*pi4/1e18/lamb4/4*0.93;
	     ext2bscattFortM[imemb][i]=kextFortM[imemb][i]/Z;
	   }
	 else
	   ext2bscattFortM[imemb][i]=1000;
	 if(ext2bscattFortM[imemb][i]>1000)
	   ext2bscattFortM[imemb][i]=1000;
       }
   }

 float bscatFortM[50][88];
 if(imsEns==1)
   {

#pragma omp parallel for default (shared) private(imemb)
     for(imemb=0;imemb<radarRet->nMemb;imemb++)
       multiscatter_(&nrange, &kextFortM[imemb][8], 
		     &ext2bscattFortM[imemb][8], &salbFortM[imemb][8], 
		     &gFortM[imemb][8],
		     &bscatFortM[imemb][8],0);

       for(i=8;i<88;i++)
	 if(z35med[i]>0)
	   {
	     float z1=log10(bscatFort[i]+1e-18)*10-dZ;
	     float z2=log10(bscatFortnoMS[i]+1e-18)*10-dZ;
	     for(imemb=0;imemb<radarRet->nMemb;imemb++)
	       {
		 float z3=log10(bscatFortM[imemb][i]+1e-18)*10-dZ;
		 radarRet->z35mod0[imemb*radarData->ngates+i]=z3;
	       }
	     
	   }
   }

 imemb=5;
 if(*j0==231 && *i0==-24)
   {
     for(i=8;i<88;i++)
       printf("%6.2f %6.2f\n",radarData->z35obs[i],
	      radarRet->z35mod0[imemb*radarData->ngates+i]);
     //exit(0);
   }
 skipMS:
//  SFM  end  06/16/2014; for M.Grecu, multiple scattering
  float *w1rand, *jrand;

  w1rand=(float*)malloc(sizeof(float)*radarRet->nMemb);
  jrand=(float*)malloc(sizeof(float)*radarRet->nMemb*50);

  for(i=0;i<radarRet->nMemb;i++)
    w1rand[i]=ran1();
  for(i=0;i<radarRet->nMemb*50;i++)
    jrand[i]=ran1();

 

  for(imemb=0;imemb<-radarRet->nMemb;imemb++)
    for(k=0;k<2*radarRet->nmfreq;k++)
      radarRet->tb[imemb*2*radarRet->nmfreq+k]=-99;


 #pragma omp parallel for default (shared) private(imemb)  
  for(imemb=0;imemb<radarRet->nMemb;imemb++)
    {
      //      if( pia13EnsM>0.95)
	{
	  calctbc(&radarRet->sfc_wind[1*imemb],&umu,kext[1*imemb],
		  salb[1*imemb],asym[1*imemb],node,&radarRet->icc[1*imemb],
		  &radarRet->jcc[1*imemb],&radarData->ngates,&radarRet->nmfreq,
		  radarData->hh,&radarRet->tb[imemb*2*radarRet->nmfreq],
		  &radarRet->emTb[imemb*2*radarRet->nmfreq],
		  &radarRet->emis[imemb*2*radarRet->nmfreq],
		  &radarRet->tpw[radarRet->icc[imemb]], &radarData->hfreez, 
		  randemiss, &imemb, 
		  &radarRet->nMemb,i0,j0,w1rand,jrand);
	  if(imemb==-100)
	    {
	      for(k=0;k<2*radarRet->nmfreq;k++)
		printf("%g ",radarRet->tb[imemb*2*radarRet->nmfreq+k]);
	      printf("\n");
	    }
	  /*
	//	exit(0);
	*/
	
      }
    }
  free(w1rand);
  free(jrand);
  for(i=0;i<nMemb;i++)
    {
      free(kext[i]);
      free(salb[i]);
      free(asym[i]);
      free(logdNw[i]);
      free(dz[i]);
//  SFM  begin  07/01/2014; for M.Grecu, random sequences
      free(xs[i]);
//  SFM  end    07/01/2014
    }

  free(dz);
  free(kext);
  free(salb);
  free(asym);
//  SFM  begin  07/01/2014; for M.Grecu, random sequences
  free(xs);
//  SFM  end    07/01/2014
  free(z35);
  free(z13);
  free(z35mod);
  free(logdNw);
  free(nodeP);
  free(imuv);
  free(rsurf);

  return;
  
}


void runEns(radarDataType   *radarData, 
	    stormStructType *stormStruct,
	    retParamType    *retParam,  
	    int *nmu, radarRetType    *radarRet, 
	    float *xscalev,
	    float *randemiss, float *localZAngle, 
	    float *wfractPix, long *ichunk, 
	    float **xs, float **logdNw,float **kext, float **asym, 
	    float **salb,float *rsurf, float **dz, int *imuv, int *nodeP, 
	    int nNodes,float nstdA)
{

  int imemb;
  float pia35M0,pia13M0,delta;
  int i,iNode;


#pragma omp parallel for default (shared) private(imemb,delta,i,pia35M0,pia13M0)
  for( imemb=0;imemb<radarRet->nMemb;imemb++)
    {
      delta=0.;
      for(i=0;i<radarData->ngates*radarRet->nmfreq;i++)
	kext[imemb][i]=salb[imemb][i]=asym[imemb][i]=-99.9;
      /* 
	 subroutine fModelFortran is a C-wraper around the Hitschfeld Bordan retrieval
	 subroutine which is written in Fortran
	 fModelFortran resides in fModelFortran.cc
	 The Fortran implementation of HB solution is resides in file fhb1.f90
      */
      iMemb=imemb;
      iNode=1;

// printf(" before fmodelfortran %i %i %i \n",*ichunk,imemb,radarData);
      if(stormStruct->rainType!=2 or radarData->hfreez<2)
	fModelFortran(radarData->z13obs, 
		      radarData->z35obs,      //  SFM  06/16/2014, for M.Grecu, multiple scattering
		      stormStruct->nodes,  
		      stormStruct->iSurf,
		      radarRet->imu[imemb],                           // this is the "unperturbed"
		      logdNw[imemb], nodeP,  nNodes,                  // logdNw  run
		      &pia35M0, &pia13M0,                             // the "unperturbed" simulated 
		      &radarRet->z35mod0[imemb*radarData->ngates], 
		      &radarRet->pwc[imemb*radarData->ngates], 	    
		      radarData->dr, radarRet->icc[imemb], 
		      radarRet->jcc[imemb],                           
		      radarData->hh, delta, iNode,
		      radarRet->nmfreq,
		      salb[imemb], kext[imemb], asym[imemb], 
		      stormStruct->rainType, 
		      radarData->ngates, 
		      &radarRet->rrate[imemb*radarData->ngates], 
		      &radarRet->d0[imemb*radarData->ngates], 	    
		      &radarRet->log10dNw[imemb*radarData->ngates],		    
		      &radarRet->z13c[imemb*radarData->ngates], 
		      &radarRet->z35[imemb*radarData->ngates], 
		      &imuv[imemb*radarData->ngates],
		      &radarData->hfreez,dz[imemb],&radarData->pia13srt,
		      //  SFM  begin  06/16/2014; for M.Grecu, multiple scattering
		      &radarData->relPia13srt, 
		      &radarData->pia35srt,
		      &radarData->relPia13srt,
		      //  SFM  begin  07/01/2014; for M.Grecu, random sequences
		      &imemb, localZAngle, wfractPix, xs[imemb], ichunk,nstdA);
      else
	fModelFortran(radarData->z13obs, 
			radarData->z35obs,      //  SFM  06/16/2014, for M.Grecu, multiple scattering
			stormStruct->nodes,  
			stormStruct->iSurf,
			radarRet->imu[imemb],                           // this is the "unperturbed"
			logdNw[imemb], nodeP,  nNodes,                  // logdNw  run
			&pia35M0, &pia13M0,                             // the "unperturbed" simulated 
			&radarRet->z35mod0[imemb*radarData->ngates], 
			&radarRet->pwc[imemb*radarData->ngates], 	    
			radarData->dr, radarRet->icc[imemb], 
			radarRet->jcc[imemb],                           
			radarData->hh, delta, iNode,
			radarRet->nmfreq,
			salb[imemb], kext[imemb], asym[imemb], 
			stormStruct->rainType, 
			radarData->ngates, 
			&radarRet->rrate[imemb*radarData->ngates], 
			&radarRet->d0[imemb*radarData->ngates], 	    
			&radarRet->log10dNw[imemb*radarData->ngates],		    
			&radarRet->z13c[imemb*radarData->ngates], 
			&radarRet->z35[imemb*radarData->ngates], 
			&imuv[imemb*radarData->ngates],
			&radarData->hfreez,dz[imemb],&radarData->pia13srt,
			//  SFM  begin  06/16/2014; for M.Grecu, multiple scattering
			&radarData->relPia13srt, 
			&radarData->pia35srt,
			&radarData->relPia13srt,
			//  SFM  begin  07/01/2014; for M.Grecu, random sequences
		      &imemb, localZAngle, wfractPix, xs[imemb], ichunk, nstdA);//, i0, j0);

      //  SFM  end    07/01/2014
      rsurf[imemb]=radarRet->rrate[imemb*radarData->ngates+stormStruct->nodes[4]];
      radarRet->pia13[imemb]=pia13M0;
      radarRet->pia35[imemb]=pia35M0;
      //printf("%6.3f %6.3f ",pia13M0,pia35M0);
//  SFM  end    06/22/2014
    }  
  //  printf("\n");
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "typedef.h"                 // this file contains the definitions of the
                                     // radar data and retrieval structures
#define WANT_STREAM


extern "C" void init_random_seed_();

extern "C" void mie2_(float *mu);     // builts the electromagnetic 
                                      // scattering look up tables 
                                      // as a function of mu




void printRet(radarDataType *radarData,stormStructType *stormStruct, radarRetType *radarRet);

extern "C" void readtables2_(int *nmu,int *nmfreq);  // reads the look up tables
extern "C" void cloud_init_(int *nmfreq);            // computes electromagnetic properties 
                                                     // at passive
                                                     // microwave frequencies for a predefined
                                                     // set of cloude profiles
                                                  
void RetStCv( radarDataType   *radarData, stormStructType *stormStruct,
	      retParamType    *retParam,  int *nmu, radarRetType    *radarRet);
int main(void)
{
  FILE *fout;
 
  stormStructType stormStruct;
  radarDataType   radarData;
  radarRetType    radarRet;
  retParamType    retParam;

  float mu;      // local variable

  int nmu=5,   // # of mu values in the look up tables 
    nmfreq=5;  // # of passive microwave frequencies in the look up tables

 
  int i, j;

  init_random_seed_();  /* initialize random number generator seed */
  
  radarData.ngates  = 101;
  radarData.z35obs  = (float*) malloc(sizeof(float)*radarData.ngates);
  radarData.z13obs  = (float*) malloc(sizeof(float)*radarData.ngates);
  radarData.hh      = (float*) malloc(sizeof(float)*radarData.ngates);
  radarData.dr      =  0.09;

  radarRet.z13c     = (float*) malloc(sizeof(float)*radarData.ngates);
  radarRet.z35mod0  = (float*) malloc(sizeof(float)*radarData.ngates);
  radarRet.pwc      = (float*) malloc(sizeof(float)*radarData.ngates);
  radarRet.rrate    = (float*) malloc(sizeof(float)*radarData.ngates);
  radarRet.d0       = (float*) malloc(sizeof(float)*radarData.ngates);
  radarRet.tb       = (float*) malloc(2*sizeof(float)*nmfreq);
  radarRet.log10dNw = (float*) malloc(sizeof(float)*radarData.ngates);

  for(i=0;i<nmu;i++)
    {
      float mu=i-2.;   
      //     mie2_(&mu);    // builds the look up tables if not already built
                     // the look up files are saved in directory TablesN
                     // subroutine mie2 is contained in file Mie/mie3.f90
                     // comment out mie2_ if the look up tables are already built
      
    } 
 
  readtables2_(&nmu,&nmfreq);  // reads the look up tables
                               // see file readTables.f90
  cloud_init_(&nmfreq);        // computes electromagnetic properties at passive
                               // microwave frequencies for a predefined set
                               // of cloud and relative humidity profiles
                               // see file cloud.f90
                                                  
  radarRet.nmfreq=nmfreq;
  
  stormStruct.nodes=(int*) malloc(sizeof(int)*5);

  stormStruct.nodes[4]=100;
 

  int ist=0;  // counts the number of stratiform "BB" profiles

  retParam.wz=1.;
  retParam.w13=1.;
  retParam.w35=1.;
  retParam.z13thresh=5;
  retParam.z35thresh=5;

  fout=fopen("MAPR/mapr2.dat","r"); //see comments in the data file
  char c1,c2,sline[100];
  c1=getc(fout);
  c2=getc(fout);

  while(1)                      // read comments in the data files
    {                           // the comment section needs to end with "*/"
      c1=c2;                    // and can contain only one "*/"
      c2=getc(fout);
      if(int(c1=='*' && c2=='/')==1) break;
    }
  fscanf(fout,"%s\n",sline);   // read the rest of the line, i.e. "\n"

  for(j=0;j<268;j++)
    {
      stormStruct.nodes[4]=radarData.ngates-1;     
                                     // this is lowest cluter free gate
      fscanf(fout,"%i %i %g %g\n",&stormStruct.nodes[0],&stormStruct.nodes[2], 
	     &radarData.pia13srt, &radarData.pia35srt);  // reads data
      stormStruct.nodes[2]=stormStruct.nodes[2]-1;
                                     // this and the next 2 statements are
                                     // a quick and dirty procedure to determine the
      stormStruct.nodes[1]=stormStruct.nodes[2]-5;
                                     // storm structure based on the BB node,
                                     // i.e. nodes[2] 
      stormStruct.nodes[3]=stormStruct.nodes[2]+4;           
                                     // the nodes are supposed to be provided
                                     // by an external, 2A23-like, 
                                     // procedure 
      for(i=0;i<=stormStruct.nodes[4];i++)   
	{
	  fscanf(fout,"%g %g %g \n", &radarData.hh[i], &radarData.z13obs[i], &
		 radarData.z35obs[i]);   
                                     // reads reflectivity data
	 
	}
    
      stormStruct.iSurf=104;
      stormStruct.nodes[4]-=3;        // substract 3 gates to make sure 
                                      // last node is outside 
                                      // the clutter region
      if(stormStruct.nodes[2]>0)      // stratiform precipitation
	{ 
	  stormStruct.rainType=1;     //  1 if stratiform 
	  RetStCv( &radarData, &stormStruct, &retParam,  &nmu, &radarRet);
	  printRet(&radarData,&stormStruct, &radarRet);
	  ist++;
	}
      else                            // "convective" precipitation, 
                                      // i.e. everything that 
                                      // does not have a brightband
	{
	  stormStruct.nodes[0]=5;
	  stormStruct.nodes[1]=32;   // nodes[1] and nodes[3] are quite 
                                     // arbitrarily set here
                                     // it is assumed that a 2A23 like procedure
                                     // will provide more rigorous values
	  stormStruct.nodes[3]=52;
	  stormStruct.nodes[2]=(stormStruct.nodes[3]+stormStruct.nodes[1])/2;
	  stormStruct.rainType=2;     // 2 if convective
	  if(radarData.pia35srt>00)     
	    {
	      RetStCv( &radarData, &stormStruct, &retParam,  &nmu, &radarRet);
	      printRet(&radarData,&stormStruct, &radarRet);
	
	    }
	}
      
    }

  free(radarData.z35obs);
  free(radarData.z13obs);
  free(radarRet.z35mod0);
  free(radarRet.z13c);
  free(radarRet.pwc);
  free(radarData.hh);
  free(radarRet.tb);
  free(radarRet.log10dNw);
}

extern "C" bool isnan_(float *var)
{
    volatile double d = *var;
    return d != d;
}

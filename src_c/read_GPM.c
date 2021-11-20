//  SFM 04/06/2013  Code module added in merge from M.Grecu's code
//  SFM 04/16/2013  Modified file name lengths from 80, 100, and 256 to 
//                    1000 for all occurrences; including "malloc" useage
//  SFM 05/06/2013  Code from LW to pass jobname into modules
//  SFM 06/19/2013  Undocumented code changes from M.Grecu
//  SFM  6/27/2013  Parameter name changes from W.Olson; reduce unused code
//  SFM  7/30/2013  zFactorMeasured value reduced by factor 100 to accommodate
//                    changes in 2AKu data format
//  LAW  5/02/2014  Received update readpr.GPM.model_05_01_14.c.gz from Bill Olson

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//begin  WSO 9/5/13 include math library for power calculation
#include <math.h>
//end    WSO 9/5/13
#include "TKheaders.h"
//begin  aaa LAW 9/4/13 1GMI to 1CGMI_r1
#include "TK_1CGMI.h"
//end    aaa LAW 9/4/13
#ifdef GFOR 
extern int __nbinmod_MOD_nbin;
#define nbins __nbinmod_MOD_nbin
//begin  WSO 9/15/13 
extern float __missingmod_MOD_missing_r4;
#define missing_r4c __missingmod_MOD_missing_r4
extern short __missingmod_MOD_missing_i2;
#define missing_i2c __missingmod_MOD_missing_i2
extern long __missingmod_MOD_missing_i4;
#define missing_i4c __missingmod_MOD_missing_i4
//end    WSO 9/15/13
#endif

int DayOfMonth[300], DayOfYear[300], Hour[300], MilliSecond[300],
  Minute[300], Month[300], Second[300], Year[300], SecondOfDay[300];
//  begin  LAW  05/08/2014; add navigation data to output file
NAVIGATION navigation[300];
//  end    LAW  05/08/2014

#ifdef IFORT 
extern int nbinmod_mp_nbin_;
#define nbins nbinmod_mp_nbin_
//begin  WSO 9/15/13
extern float missingmod_mp_missing_r4_;
#define missing_r4c missingmod_mp_missing_r4_
extern short missingmod_mp_missing_i2_;
#define missing_i2c missingmod_mp_missing_i2_
extern long missingmod_mp_missing_i4_;
#define missing_i4c missingmod_mp_missing_i4_
//end    WSO 9/15/13
#endif
TKINFO granuleHandle2AKu;
TKINFO tkfileinfo;

TKINFO dprtkfile;
extern L2ADPR_SWATHS dprswath;

//  begin  SFM  09/12/2013
extern TKINFO ctkfile;
int iYearSet=0, iDaySet=0;  // 4/14/MG
//  end  SFM  09/12/2013

//  begin WSO 8/14/2014
unsigned long getbits(unsigned long x, int p, int n);
unsigned long temp_byte;
//  end   WSO 8/14/2014

//  SFM  begin 12/06/2013; to pass out TK status message
int readgmi_(char *jobname, char *f1cgmi,  int *n1b11, float *tmilow,
             float *tmihigh, float *xlon,float *xlat, float *gmilon, 
	     float *gmilat, int *mm, int *yyyy, int *jday, 
	     int *dd, int *ifile, float *secondOfDay,
	     float *SCLon, float *SCLat, float *SCAlt, int *SCYaw,
	     float *xlonS2, float *xlatS2,
	     float *eiaS1, float *eiaS2, float *sgaS1)  //4/14/14 MG
{
  int             status_alpha ;
//  SFM  end   12/06/2013
  char            *granuleID;
  int             dataType;
  char            filemode;
  int             status, orbNumber,i, ind, ipia, indh,j, k, nscans, nscanmax;

//  TKINFO          tkfileinfo;
  L1CGMI_SWATHS   data;
  char emptygranuletext[100];
 
  i=0;

  granuleID=(char*)malloc(sizeof(char)*1000);
  strcpy(&granuleID[0],f1cgmi);
  printf(" readgmi: f1cgmi %s \n",f1cgmi);
  
//  SFM  begin 12/06/2013; to pass out TK status message
//begin  LAW  9/4/13 1CGMI to 1CGMI_r1
//aa   status = TKopen(granuleID,"1CGMI", TKREAD, "HDF5", jobname, &tkfileinfo, 1);
    status_alpha = TKopen(granuleID,"1CGMI", TKREAD, "HDF5", jobname, &tkfileinfo, 1);
    if (status_alpha != 0)
      {
       printf("WARNING: Unable to access 1CGMI %i \n",status_alpha);
       return status_alpha ;
      }
//  SFM  end   12/06/2013   
//end    LAW  9/4/13

//  begin  SFM  09/12/2013
    TKtransferMetaData(&tkfileinfo, &ctkfile);
    printf("  STATUS, meta transfer 1CGMI %i  \n",status_alpha);
//  end  SFM  09/12/2013
 
  TKgetMetaString(&tkfileinfo, "FileHeader", "EmptyGranule", 
                  emptygranuletext);	  
  if (strncmp(emptygranuletext,"EMPTY",5) == 0)
     {
      status_alpha = -3;
      printf(" STATUS, 1CGMI has empty granule \n",status_alpha);
      return status_alpha ;
      }
      
  status = TKgetScanCount ( &tkfileinfo, "S1" );
    
  *n1b11=0;
  
  status = TKseek(&tkfileinfo, 0,
                  TK_ABS_SCAN_OFF);
  nscans=0;
  while (TKendOfFile (&tkfileinfo) != TK_EOF)
    {
      status=TKreadScan(&tkfileinfo,&data);
      nscans++;
    }

  printf("   number scans = %i \n",nscans);
  if(*ifile==1)
    {
      status = TKseek(&tkfileinfo, nscans-200,
		      TK_ABS_SCAN_OFF);
    }
  else
    {
      status = TKseek(&tkfileinfo, 0,
		      TK_ABS_SCAN_OFF);
    }
  if(*ifile==3)
    nscanmax=200;
  else
    nscanmax=3500;
      
  ind=0;
  ipia=0;
  indh=0;

  int dday;  //4/14/14 MG
  while (TKendOfFile (&tkfileinfo) != TK_EOF && *n1b11<nscanmax)
    {
      status=TKreadScan(&tkfileinfo,&data);

      *mm=data.S1.ScanTime.Month;
      *yyyy=data.S1.ScanTime.Year;
      *jday=data.S1.ScanTime.DayOfYear;
      *dd=data.S1.ScanTime.DayOfMonth;
      
//  SFM  begin  04/16/2014; for M.Grecu, revision of nodes processing
      if(iYearSet==0)                                     //4/14/14 MG Begin
	{
	  iYearSet=*yyyy;
	  iDaySet=*jday;
	  dday=0;
	}
      else
	{
	  if(iYearSet/4*4==iYearSet)
	    dday=(iYearSet-*yyyy)*366*24*3600;
	  else
	    dday=(iYearSet-*yyyy)*365*24*3600.;
	  dday+=(iDaySet-*jday)*24*3600.;
	}
      secondOfDay[*n1b11]=dday+data.S1.ScanTime.SecondOfDay;     //4/14/14 MG End
      SCLat[*n1b11]=data.S1.SCstatus.SClatitude;
      SCLon[*n1b11]=data.S1.SCstatus.SClongitude;
      SCAlt[*n1b11]=data.S1.SCstatus.SCaltitude;
      SCYaw[*n1b11]=data.S1.SCstatus.SCorientation;

//  SFM  end  04/16/2014

      (*n1b11)+=1;

      for(j=0; j<221; j++)
	{
	  for(k=0;k<4;k++)
	    {
	      tmihigh[indh]=data.S2.Tc[j][k]/1.;
	      indh++;
	    }
	  for(k=0;k<9;k++)
	    {
	      tmilow[ind]=data.S1.Tc[j][k]/1.;
//  SFM  begin  12/24/2013; diagnostic to insert temperature defaults
//              if (j<110) tmilow[ind] = -999.0 ;
//  SFM  end    12/24/2013
	      ind++;
	    }
	}
      for(j=0;j<221;j++)
	{
	  eiaS1[ipia]=*data.S1.incidenceAngle[j];
	  eiaS2[ipia]=*data.S2.incidenceAngle[j];
	  sgaS1[ipia]=*data.S1.sunGlintAngle[j];
	  xlon[ipia]=data.S1.Longitude[j];
	  xlat[ipia]=data.S1.Latitude[j];
	  xlonS2[ipia]=data.S2.Longitude[j];
	  xlatS2[ipia]=data.S2.Latitude[j];
	  ipia++;
	}

    }

//  SFM  begin 12/06/2013; to pass out TK status message
  status_alpha = TKclose(&tkfileinfo);
  free(granuleID);    
  return status_alpha ;
//  SFM  end   12/06/2013; to pass out TK status message
}

void get_scAngle(float *scAngle);
TKINFO       granuleHandle2AKu;
char fileName[1000];

//  SFM  begin 12/12/2013; to pass out TK status message
//  SFM  begin  04/16/2014; for M.Grecu, revision of nodes processing
//begin  WSO 08/18/14; add ioqualityflags pass to main routine



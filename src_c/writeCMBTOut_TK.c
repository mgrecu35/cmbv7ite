//  SFM 04/06/2013  Code module added in merge from M.Grecu's code
//  SFM 05/06/2013  Modifications from LW to facilitate using job names
//  SFM 06/27/2013  Parameter name changes from W.Olson; reduce unused code
//  SFM 07/19/2013  Large volume of code added for M.Grecu
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TKheaders.h"
#include "TK_2BCMBT.h"
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


void openoutputfilet_(char *jobname, char *fname)
{
  char ranstring[80], ranMsg[80];
  int status;
  printf(" Output file = %s \n",fname);

  int ret;
//  SFM  04/06/2013  Changed file type to 2BCMB   WSO 04/07/2013
//  SFM  09/04/2013  Changed jobid to jobname
//  SFM  09/11/2013  Moved metadta settings to centralized location

  ret = TKopen(fname, "2BCMBT", TKWRITE, "HDF5", jobname, &ctkfile,1); //WSO 04/07/2013
}

void openinputfilet_(char *jobname, char *fname)
{
  int ret;

//  SFM  04/06/2013  Changed file type to 2BCMB   WSO 04/07/2013
  ret = TKopen(fname, "2BCMB", TKREAD, "HDF5", jobname, &ctkfileIn, 1);
  printf("%s %i \n",fname,ret);
}



void closeoutputfilet_(void)
{
  int ret;
  ret=TKclose(&ctkfile);
}

void closedproutputfilet_(void)
{
  int ret;
  ret=TKclose(&dprtkfile);
}


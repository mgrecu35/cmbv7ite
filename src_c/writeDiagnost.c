//  SFM  04/06/2013  Code changes from M.Grecu
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "hdf5.h"

#define RANK  3
#define DIM0_SUB  1                         /* subset dimensions */ 
#define DIM1_SUB  49 
#define DIM2_SUB  9 


#define DIM0     8                          /* size of dataset */       
#define DIM2     9 
#define DIM1     49 

int writediag_(char *filename, int *itime, int *ndpr, int *istrec, int *nrec, float *tbrgrid) 
{
  hid_t       file_id, group_id, group1_id, dataset_id, dataspace_id, memspace_id;  /* identif
									  iers */
  hsize_t     dims[3];
  herr_t      status;
  int         j, dset1_data[9300][49][9], dset2_data[2][10];
  int         sdata[1][49][9];
  

 
 dims[0] = *ndpr;
 dims[1] = 49;
 dims[2] = 9;
 if(*itime==0)
   {
     file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
     dataspace_id = H5Screate_simple(3, dims, NULL);
     
     /* Create a dataset in group "MyGroup". */
     
     group1_id = H5Gcreate2(file_id, "/DiagGroup",  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     dataset_id = H5Dcreate2(file_id, "/DiagGroup/dset1", H5T_STD_I32BE, dataspace_id,
			     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     
     status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		       dset1_data);
     
     
     status = H5Sclose(dataspace_id);
 
 /* Close the second dataset */
     status = H5Dclose(dataset_id);
     
     /* Close the group. */
     status = H5Gclose(group1_id);
     status = H5Fclose(file_id);
   }


 hsize_t     count[3];              /* size of subset in the file */
 hsize_t     offset[3];             /* subset offset in the file */
 hsize_t     stride[3];
 hsize_t     block[3];
 hsize_t     dimsm[3];
 
// printf("%s \n",filename);
 file_id = H5Fopen (filename, H5F_ACC_RDWR, H5P_DEFAULT);
 dataset_id = H5Dopen (file_id,"/DiagGroup/dset1", H5P_DEFAULT);
 
 /* Specify size and shape of subset to write. */
 
 offset[0] = 1;
 offset[1] = 0;
 offset[2] = 0;
 
 count[0]  = DIM0_SUB;  
 count[1]  = DIM1_SUB;
 count[2]  = DIM2_SUB;
 
 stride[0] = 1;
 stride[1] = 1;
 stride[2] = 1;
 
 block[0] = 1;
 block[1] = 1;
 block[2] = 1;
 
 /* Create memory space with size of subset. Get file dataspace 
    and select subset from file dataspace. */
 
 dimsm[0] = DIM0_SUB;
 dimsm[1] = DIM1_SUB;
 dimsm[2] = DIM2_SUB;
 memspace_id = H5Screate_simple (RANK, dimsm, NULL); 
 
 dataspace_id = H5Dget_space (dataset_id);
 int i,ii,jj;
 int ic=0;
 for(i=0;i<*nrec;i++)
   {
     offset[0] = *istrec+i;
     status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset,
			       stride, count, block);
 
     /* Write a subset of data to the dataset, then read the 
	entire dataset back from the file.  */
     
 
       for(ii=0;ii<49;ii++)
	 for(jj=0;jj<9;jj++)
	   {
	     sdata[0][ii][jj]=tbrgrid[ic];
	     ic++;
	   }
     status = H5Dwrite (dataset_id, H5T_NATIVE_INT, memspace_id,
			dataspace_id, H5P_DEFAULT, sdata);
   }
 
 status = H5Sclose (memspace_id);
 status = H5Sclose (dataspace_id);
 status = H5Dclose (dataset_id);
 status = H5Fclose (file_id);
 
}


int writediagnew_(char *filename, int *itime, int *ndpr, int *istrec, int *nrec, float *tbrgrid) 
{
  hid_t       file_id, group_id, group1_id, dataset_id, dataspace_id, memspace_id;  /* identif
									  iers */
  hsize_t     dims[3];
  herr_t      status;
  int         j, dset1_data[9300][49][9], dset2_data[2][10];
  int         sdata[1][49][9];
  

 
 dims[0] = *ndpr;
 dims[1] = 49;
 dims[2] = 9;
 //printf("%s %i\n",filename,*itime);
 if(*itime==0)
   {
     //  printf("%s \n",filename);
     file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
     dataspace_id = H5Screate_simple(3, dims, NULL);
     
     /* Create a dataset in group "MyGroup". */
     
     group1_id = H5Gcreate2(file_id, "/DiagGroup",  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     dataset_id = H5Dcreate2(file_id, "/DiagGroup/dset1", H5T_STD_I32BE, dataspace_id,
			     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     
     
     status = H5Sclose(dataspace_id);
 
 /* Close the second dataset */
     status = H5Dclose(dataset_id);
     
     /* Close the group. */
     status = H5Gclose(group1_id);
     status = H5Fclose(file_id);
   }
 //return;
 
 hsize_t     count[3];              /* size of subset in the file */
 hsize_t     offset[3];             /* subset offset in the file */
 hsize_t     stride[3];
 hsize_t     block[3];
 hsize_t     dimsm[3];
 
 file_id = H5Fopen (filename, H5F_ACC_RDWR, H5P_DEFAULT);
 dataset_id = H5Dopen (file_id,"/DiagGroup/dset1", H5P_DEFAULT);
 
 /* Specify size and shape of subset to write. */
 
 offset[0] = 1;
 offset[1] = 0;
 offset[2] = 0;
 
 count[0]  = DIM0_SUB;  
 count[1]  = DIM1_SUB;
 count[2]  = DIM2_SUB;
 
 stride[0] = 1;
 stride[1] = 1;
 stride[2] = 1;
 
 block[0] = 1;
 block[1] = 1;
 block[2] = 1;
 
 /* Create memory space with size of subset. Get file dataspace 
    and select subset from file dataspace. */
 
 dimsm[0] = DIM0_SUB;
 dimsm[1] = DIM1_SUB;
 dimsm[2] = DIM2_SUB;
 memspace_id = H5Screate_simple (RANK, dimsm, NULL); 
 
 dataspace_id = H5Dget_space (dataset_id);
 int i,ii,jj;
 int ic=0;
 for(i=0;i<*nrec;i++)
   {
     offset[0] = *istrec+i;
     status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset,
			       stride, count, block);
 
     /* Write a subset of data to the dataset, then read the 
	entire dataset back from the file.  */
     
 
       for(ii=0;ii<49;ii++)
	 for(jj=0;jj<9;jj++)
	   {
	     sdata[0][ii][jj]=tbrgrid[ic];
	     ic++;
	   }
     status = H5Dwrite (dataset_id, H5T_NATIVE_INT, memspace_id,
			dataspace_id, H5P_DEFAULT, sdata);
   }
 
 status = H5Sclose (memspace_id);
 status = H5Sclose (dataspace_id);
 status = H5Dclose (dataset_id);
 status = H5Fclose (file_id);
 
}

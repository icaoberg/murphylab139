/*
 * Copyright (C) 2006 Murphy Lab,Carnegie Mellon University
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 * 
 * For additional information visit http://murphylab.web.cmu.edu or
 * send email to murphy@cmu.edu
 */
/*/////////////////////////////////////////////////////////////////////////
//
//
//                            ml_3Dtexture.c 
//
//
//                            Xiang Chen
//                            modified from mb_texture.c
//                            Aug 19, 2002
//
//  Revisions:
//
/////////////////////////////////////////////////////////////////////////*/


#include "mex.h"
#include "matrix.h"
#include "Include/ppgm.h"
#include "Include/CVIPtexture.h"
#include <sys/types.h>
#include <stdlib.h>

#define y 0
#define x 1
#define z 2

#define row 0
#define col 1

#define ind(x, y, z)             (y) + (x) * ny + (z) * ny * nx

extern TEXTURE * Extract_Texture_Features();

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  
  int         distance;                 /*parameter for texture calculations*/
  u_int8_t*   p_img;                    /*The image from Matlab*/
  u_int8_t*** p_gray;                   /*Image converted for texture calcs*/
  int         ny;                       /*Image y*/
  int         nx;                       /*Image x*/
  int         nz;                       /*Image z*/
  TEXTURE_FEATURE_MAP* features_used ;  /*Indicate which features to calc.*/
  TEXTURE*    features ;                 /*Returned struct of features*/
  int         NDims;
  long        Nvoxels;
  const int*  dims;
  int         Dims[3];
  int         i, j, k ;
  int         imgsize[3] ;              
  int         imgindex[3] ;
  long        offset ;                  
  int         outputsize[2] ;           /*Dimensions of TEXTURE struct*/
  int         outputindex[2] ;
  float*      output ;                  /*Features to return*/
  u_int8_t*   outtest ;

  if (nrhs != 1) {
    mexErrMsgTxt("ml_3Dtexture requires a single input argument.\n") ;
  } else if (nlhs != 1) {
    mexErrMsgTxt("ml_3Dtexture returns a single output.\n") ;
  }

  if (!mxIsNumeric(prhs[0])) {
    mexErrMsgTxt("ml_3Dtexture requires a single numeric input.\n") ;
  }

  if (!mxIsUint8(prhs[0])) {
    mexErrMsgTxt("ml_3Dtexture requires a single input of type unsigned 8-bit integer.\n") ;
  }

  NDims = mxGetNumberOfDimensions(prhs[0]);
  
  if (NDims != 3) {
    mexErrMsgTxt("ml_3Dtexture requires a 3D image as the input.\n");
  }

  dims = mxGetDimensions(prhs[0]);
  Dims[0] = dims[0];
  Dims[1] = dims[1];
  Dims[2] = dims[2];
  ny = Dims[0];
  nx = Dims[1];
  nz = Dims[2];

  if(!(nx > 1) || !(ny > 1) || !(nz > 1)) {
    mexErrMsgTxt("ml_3Dtexture requires an input 3D image, not a scalar.\n") ;
  }

  Nvoxels = x * y * z;

  if (!mxIsUint8(prhs[0])) {
    mexErrMsgTxt("Input image should be of class Uint8.\n");
  }

  p_img = (u_int8_t*)mxGetData(prhs[0]) ;

  distance = 1 ;

  features_used = mxCalloc(1, sizeof(TEXTURE_FEATURE_MAP)) ;
  if(!features_used) 
    mexErrMsgTxt("ml_3Dtexture: error allocating features_used.") ;

  features_used->ASM = 1 ;
  features_used->contrast = 1 ;
  features_used->correlation = 1 ;
  features_used->variance = 1 ;
  features_used->IDM = 1 ;
  features_used->sum_avg = 1 ;
  features_used->sum_var = 1 ;
  features_used->sum_entropy = 1 ;
  features_used->entropy = 1 ;
  features_used->diff_var = 1 ;
  features_used->diff_entropy = 1 ;
  features_used->meas_corr1 = 1 ;
  features_used->meas_corr2 = 1 ;
  features_used->max_corr_coef = 0 ;

  imgsize[x] = nx ;
  imgsize[y] = ny ;
  imgsize[z] = nz;
  
  p_gray = mxCalloc(nx, sizeof(u_int8_t**)) ;
  if(p_gray) {
    for(i=0; i<nx ; i++) {
      p_gray[i] = mxCalloc(ny, sizeof(u_int8_t*)) ;
      if (p_gray[i]) {
	for (j = 0; j < ny; j ++) {
	  p_gray[i][j] = mxCalloc(nz, sizeof(u_int8_t));
	  if(!p_gray[i]) mexErrMsgTxt("ml_3Dtexture : error allocating p_gray[i]") ;
	}
      } else mexErrMsgTxt("ml_3Dtexture : error allocating p_gray") ;
    }
  } else mexErrMsgTxt("ml_3Dtexture : error allocating p_gray");

  for(k = 0 ; k < nz ; k++) 
    for(j = 0 ; j < ny ; j++) 
      for (i = 0; i < nx; i++) { 
	offset = ind(i, j, k) ;
	p_gray[i][j][k] = p_img[offset] ;
      }
 
  features=Extract_Texture_Features(distance,p_gray,nx,ny,nz,features_used) ;
  /*
  outputsize[row] = mrows ;
  outputsize[col] = ncols ;
  plhs[0] = mxCreateNumericArray(2, outputsize, mxUINT8_CLASS, mxREAL) ;
  if (!plhs[0]) mexErrMsgTxt("mb_texture: error allocating return variable.") ;
  outtest = (u_int8_t*)mxGetData(plhs[0]) ;
  for(imgindex[row] = 0 ; imgindex[row] < imgsize[row] ; imgindex[row]++) 
    for(imgindex[col] = 0 ; imgindex[col] < imgsize[col] ; imgindex[col]++ ) {
      offset = mxCalcSingleSubscript(prhs[0], 2, imgindex) ;
      outtest[offset] = p_gray[imgindex[row]][imgindex[col]] ;
    }
  */


  outputsize[row] = 14 ;
  outputsize[col] = 15 ;
 
  plhs[0] = mxCreateNumericArray(2, outputsize, mxSINGLE_CLASS, mxREAL) ;
  if (!plhs[0]) mexErrMsgTxt("ml_3Dtexture: error allocating return variable.") ;

  output = (float*)mxGetData(plhs[0]) ;

  /* Copy the features into the return variable */

  for(outputindex[col]=0 ; outputindex[col] < outputsize[col] ; 
      outputindex[col]++) {
    outputindex[row]=0 ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->ASM[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->contrast[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->correlation[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->variance[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->IDM[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->sum_avg[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->sum_var[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->sum_entropy[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->entropy[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->diff_var[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->diff_entropy[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->meas_corr1[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->meas_corr2[outputindex[col]] ;

    outputindex[row]++ ;
    offset =  mxCalcSingleSubscript(plhs[0], 2, outputindex) ;
    output[offset] = features->max_corr_coef[outputindex[col]] ;
    
  }

 

  /*
    Memory clean-up.
  */
  for(i=0; i<nx ; i++) {
    for( j = 0; j < ny; j++) {
      mxFree(p_gray[i][j]) ;
    }
    mxFree(p_gray[i]);
  }
  mxFree(p_gray) ;  
  mxFree(features_used) ;
  /* val is calloc'd inside of Extract_Texture_Features */
  free(features) ;
  
}

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
#ifdef _MEX_
#include "mex.h"
#include "matrix.h"
#endif /* #ifdef _MEX_ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <memory.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/types.h>

#ifdef _MEX_
void mexFunction(
		 int nlhs,
		 mxArray *plhs[], 
		 int nrhs,
		 const mxArray *prhs[])
{
  unsigned char* output ;     /* Return value to Matlab */

    int x, y, z, i;
    long n_voxels; /* The total number of voxels in the image */

    unsigned char *Image;
    long NDims;
    const int *Dims;
    int Width, Height, Depth;

    unsigned char MostCommonPixelValue;
    int newvalue;
    unsigned char *threshp;
    unsigned char thresh;
    

    /* Argument checking */
    if (nrhs != 2) {
      mexErrMsgTxt("ml_binarize() requires two input arguments.") ;
    } else if (nlhs != 1) {
      mexErrMsgTxt("ml_binarize() returns a single output.") ;
    }

    /* Get Image size info etc */
    NDims = mxGetNumberOfDimensions( prhs[0]);
    Dims = mxGetDimensions( prhs[0]);
    n_voxels = 1;
    for( i = 0; i < NDims; i++) {
      n_voxels = n_voxels * Dims[i]; 
    }

    /* Get hold of the input image */
    Image = (unsigned char*) mxGetData( prhs[0]);

    /* Get threshold value */
    if( !mxIsUint8( prhs[1])) 
      mexErrMsgTxt("THRESH argument must be of class uint8") ;
    threshp = (unsigned char*) mxGetData( prhs[1]);
    thresh = *threshp;
    /*printf("Thresh = %i",thresh); */

    /* Create output image */
    plhs[0] = mxCreateNumericArray(NDims, Dims, mxUINT8_CLASS, mxREAL) ;
    output = (unsigned char*) mxGetData(plhs[0]) ;

    /*////////////////////////////////////////////////////////////////*/
    for( i = 0; i < n_voxels; i++) {
      if( Image[i] >= thresh) {
	output[i] = 1;
      } else {
	output[i] = 0;
      }
    }
    /*////////////////////////////////////////////////////////////////*/

    return ;
}

#endif /* #ifdef _MEX_ */

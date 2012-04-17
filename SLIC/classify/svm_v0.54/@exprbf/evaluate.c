/******************************************************************************

File        : @exprbf/evaluate.c

Author      : Kai Huang

Description : Mex file reimplemetation of the kernel evaluation method of a
              MATLAB class implementing an exponential radial basis kernel.

******************************************************************************/

#include <math.h>

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   double *x1, *x2, *y, gamma;

   int row, i, j, k, m, n, o, N;

   const mxArray *net;

   /* check number of input and output arguments */

   if (nrhs != 3)
   {
      mexErrMsgTxt("Wrong number input arguments.");
   }
   else if (nlhs > 1)
   {
      mexErrMsgTxt("Too many output arguments.");
   }

   /* get kernel structure */

   gamma = mxGetScalar(mxGetField(prhs[0],0,"gamma"));

   /* get input arguments */

   if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
   {
      mexErrMsgTxt("x1 must be a double matrix.");
   }

   m  = mxGetM(prhs[1]);
   x1 = mxGetPr(prhs[1]);

   if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]))
   {
      mexErrMsgTxt("x2 must be a double matrix.");
   }

   n  = mxGetM(prhs[2]);
   o  = mxGetN(prhs[2]);
   x2 = mxGetPr(prhs[2]);

   /* allocate and initialise output matrix */

   plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);

   y = mxGetPr(plhs[0]);

   /* compute kernel matrix */

   for (i = 0; i < m; i++)
   {
      for (j = 0; j < n; j++)
      {
         for (k = 0; k < o; k++)
         {
            y[i+j*m] += fabs(x1[i+k*m] - x2[j+k*n]);
         }
      }
   }

   for (i = 0; i < m*n; i++)
   {
      y[i] = exp(-gamma*y[i]);
   }

   /* bye bye... */
}

/**************************** That's all Folks!  *****************************/


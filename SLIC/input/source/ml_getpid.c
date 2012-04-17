#include "mex.h"

/*Get Process ID for current matlab*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double* p_pid;
  int outputsize[2];
  outputsize[0] = 1;
  outputsize[1] = 1;
  plhs[0] = mxCreateNumericArray(1,outputsize,mxDOUBLE_CLASS,mxREAL);
  if (!plhs[0]) 
    mexErrMsgTxt("ml_getpid: error allocating return variable.");
  p_pid = (double *)mxGetData(plhs[0]);
  p_pid[0] = getpid();
}

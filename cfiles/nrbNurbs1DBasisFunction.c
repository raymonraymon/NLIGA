#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "nrbNurbs.h"
#include "mex.h"

#define tol    1e-15

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*
	* The function will return all related 1D NURBS basis functions
	* at point xi to matlab.
	* 
	* It can be recalled in form of:
	* R = nrbNurbs1DBasisFunction(xi, p, uknot, weights)
	* Input:
	*   xi      - point, xi, where we want to interpolate
	*   p       - degree of basis
	*   uknot   - knot vectors
	*   weights - vector of weights
	* Output:
	*   R       - all non-zero basis functions
	*/
	
	// get the inputs
	
	double *xi = mxGetPr(prhs[0]);	
	double *pd = mxGetPr(prhs[1]);
    
	int p = (int) *pd;
	
	double *uknot = mxGetPr(prhs[2]);    
	int uknot_num = mxGetN(prhs[2]);
    
	int u_ctrlPtsNum = uknot_num - 1 - p - 1;
    int noFuncs      = p+1; 
	
	double *weight = mxGetPr(prhs[3]);
    int weight_num = mxGetN(prhs[3]);
	
	double *Nu = (double *)malloc(sizeof(double)*(p+1));
	
	int uspan = nrbFindSpan(u_ctrlPtsNum, p, xi[0], uknot); 
	
	int i, c;   
    double fac, wgt, w = 0.0;
	
	double *R;
	
	
	if(nrhs != 4) mexErrMsgTxt("Please pass in 4 arguments to the function."
		"We expect it to be in the form R = nrbNurbs1DBasisFunction(xi, p uknot, weights)\n");
	
	if(fabs(xi[0]-uknot[u_ctrlPtsNum-1]) < tol) 
		xi[0] = uknot[u_ctrlPtsNum-1] - tol;

    // Calculate the non-zero univariate B-spline basis functions
	
	nrbAllBsplineBasisFuns(uspan, xi[0], p, uknot, Nu);

	// create NURBS approximation 
		
	for(i = 0; i <= p; i++)
	{      
		c = uspan - p + i;	           
		wgt = weight[c];
		
		w += Nu[i] * wgt;
	}
	
	// create output 
    
    plhs[0] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
	
	R      = mxGetPr(plhs[0]);

	for(i = 0; i <= p; i++)
	{   
		c = uspan - p + i;        
		
		R[i] = Nu[i]*weight[c]/w;
		
		fac = weight[c]/(w*w);
	}
	
	free(Nu);
	
}
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
	* The function will return all related 1D NURBS basis functions and 
	* their first derivatives at point xi to matlab.
	* 
	* It can be recalled in form of:
	* [R dRdxi] - nrbNurbs1DBasisDerivs(xi, p, uknot, weights)
	* Input:
	*   xi      - point in form of [xi] for calculation
	*   p       - degree of basis
	*   uknot   - knot vectors
	*   weights - vector of weights
	* Output:
	*   R       - all non-zero basis functions
	*   dRdxi   - derivatives
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
	double **Nu_ders = init2DArray(u_ctrlPtsNum+1, p+1);
	
	int uspan = nrbFindSpan(u_ctrlPtsNum, p, xi[0], uknot); 
	
	int i, c;   
    double fac, wgt, w = 0.0, dwdxi = 0.0;
	
	double *R, *dRdxi;
	
	
	if(nrhs != 4) mexErrMsgTxt("Please pass in 4 arguments to the function."
		"We expect it to be in the form [R dRdxi] = nrbNurbs1DBasisDerivs(xi, p uknot, weights)\n");
	
	if(fabs(xi[0]-uknot[u_ctrlPtsNum-1]) < tol) 
		xi[0] = uknot[u_ctrlPtsNum-1] - tol;

    // Calculate the non-zero univariate B-spline basis functions and first derivatives 
	
	nrbAllBsplineBasisFuns(uspan, xi[0], p, uknot, Nu);
	
	nrbAllBsplineBasisFunsDers(uspan, xi[0], p, uknot, u_ctrlPtsNum, Nu_ders);

	// create NURBS approximation
		
	for(i = 0; i <= p; i++)
	{      
		c = uspan - p + i;	           
		wgt = weight[c];
		
		w += Nu[i] * wgt;
		dwdxi  += Nu_ders[1][i] * wgt;
	}
	
	// create output
    
    plhs[0] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
	
	R      = mxGetPr(plhs[0]);
    dRdxi  = mxGetPr(plhs[1]);

	for(i = 0; i <= p; i++)
	{   
		c = uspan - p + i;        
		
		R[i] = Nu[i]*weight[c]/w;
		
		fac = weight[c]/(w*w);
		dRdxi[i] = (Nu_ders[1][i]*w - Nu[i]*dwdxi) * fac;
	}
	
	free(Nu);
    free2Darray(Nu_ders, (u_ctrlPtsNum+1));
	
}
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
	* The function will return all related 2D NURBS basis functions
    * at point [xi, eta]to matlab.
	* 
	* It can be recalled in form of:
	* [R] = nrbNurbs2DBasisFunction(xi, p, q, uknot, vknot, weights)
	* Input:
	*   xi            - point in form of [xi, eta] for calculation
	*   p,q           - degrees of basis
	*   uknot,vknot   - knot vectors
	*   weights       - vector of weights
	* Output:
	*   R             - all non-zero basis functions
	*/
	
	// get the inputs 
	
	double *xi = mxGetPr(prhs[0]);	
	double *pd = mxGetPr(prhs[1]);
	double *qd = mxGetPr(prhs[2]);
    
	int p = (int) *pd;
	int q = (int) *qd;
	
	double *uknot = mxGetPr(prhs[3]);
    double *vknot = mxGetPr(prhs[4]);
    
	int uknot_num = mxGetN(prhs[3]);
    int vknot_num = mxGetN(prhs[4]);
    
	int u_ctrlPtsNum = uknot_num - 1 - p - 1;
    int v_ctrlPtsNum = vknot_num - 1 - q - 1;
    int noFuncs      = (p+1)*(q+1); 
	
	double *weight = mxGetPr(prhs[5]);
    int weight_num = mxGetN(prhs[5]);
	
	double *Nu = (double *)malloc(sizeof(double)*(p+1));
    double *Nv = (double *)malloc(sizeof(double)*(q+1));
	
	int uspan = nrbFindSpan(u_ctrlPtsNum, p, xi[0], uknot); 
    int vspan = nrbFindSpan(v_ctrlPtsNum, q, xi[1], vknot); 
	
	int i, j, c, k = 0;
	int uind, vind;   
    double fac, wgt, w = 0.0;
	
	double *R;
	
	if(nrhs != 6) mexErrMsgTxt("Please pass in 6 arguments to the function."
		"We expect it to be in the form [R dRdxi dRdeta] = nrbNurbs2DBasisDerivs(xi, p, q, uknot, vknot, weights)\n");
	
	if(fabs(xi[0]-uknot[u_ctrlPtsNum-1]) < tol) 
		xi[0] = uknot[u_ctrlPtsNum-1] - tol;
    
	if(fabs(xi[1]-vknot[v_ctrlPtsNum-1]) < tol) 
		xi[1] = vknot[v_ctrlPtsNum-1] - tol;
	
    // Calculate the non-zero univariate B-spline basis functions and first derivatives 
	
	nrbAllBsplineBasisFuns(uspan, xi[0], p, uknot, Nu);
	nrbAllBsplineBasisFuns(vspan, xi[1], q, vknot, Nv);
	
	// create NURBS approximation 
	
	for(j = 0; j <= q; j++)
    {
		vind = vspan - q + j;
        for(i = 0; i <= p; i++)
        {      
			uind = uspan - p + i; 	
            c = vind*(u_ctrlPtsNum+1)+uind;            
            wgt = weight[c];
            
            w += Nu[i] * Nv[j] * wgt;
        }
    }
	
	// create output 
    
    plhs[0] = mxCreateDoubleMatrix(1,noFuncs,mxREAL);
	R = mxGetPr(plhs[0]);
	
	for(j = 0; j <= q; j++)
    {
		vind = vspan - q + j;
        for(i = 0; i <= p; i++)
        {      
			uind = uspan - p + i; 	
            c = vind*(u_ctrlPtsNum+1)+uind;            
            fac = weight[c]/w;
            
            R[k] = Nu[i]*Nv[j]*fac;
            
            k += 1;
        }
    }
	
	free(Nu);
    free(Nv);
	
}
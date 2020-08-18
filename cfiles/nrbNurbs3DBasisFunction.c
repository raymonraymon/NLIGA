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
	* The function will return all non-zero 3D NURBS basis functions 
    *	at point [xi,eta,zeta] to matlab.
	* 
	* It can be recalled in form of:
	* [R] = nrbNurbs3DBasisDerivs(xi, p, q, r, uknot, vknot, wknot, weights)
	* Input:
	*   xi                 - point in form of [xi,eta,zeta] for calculation
	*   p,q,r              - degree of basis
	*   uknot,vknot,wknot  - knot vectors
	*   weights            - vector of weights
	* Output:
	*   R                  - all non-zero basis functions
	*/
	
	// get the inputs
	
	double *xi = mxGetPr(prhs[0]);	
	double *pd = mxGetPr(prhs[1]);
	double *qd = mxGetPr(prhs[2]);
	double *rd = mxGetPr(prhs[3]);
    
	int p = (int) *pd;
	int q = (int) *qd;
	int r = (int) *rd;
	
	double *uknot = mxGetPr(prhs[4]);
    double *vknot = mxGetPr(prhs[5]);
	double *wknot = mxGetPr(prhs[6]);
    
	int uknot_num = mxGetN(prhs[4]);
    int vknot_num = mxGetN(prhs[5]);
	int wknot_num = mxGetN(prhs[6]);
    
	int u_ctrlPtsNum = uknot_num - 1 - p - 1;
    int v_ctrlPtsNum = vknot_num - 1 - q - 1;
	int w_ctrlPtsNum = wknot_num - 1 - r - 1;
    int noFuncs      = (p+1)*(q+1)*(r+1); 
	
	double *weight = mxGetPr(prhs[7]);
    int weight_num = mxGetN(prhs[7]);
	
	double *Nu = (double *)malloc(sizeof(double)*(p+1));
    double *Nv = (double *)malloc(sizeof(double)*(q+1));
	double *Nw = (double *)malloc(sizeof(double)*(r+1));
	
	int uspan = nrbFindSpan(u_ctrlPtsNum, p, xi[0], uknot); 
    int vspan = nrbFindSpan(v_ctrlPtsNum, q, xi[1], vknot); 
	int wspan = nrbFindSpan(w_ctrlPtsNum, r, xi[2], wknot);
	
	int i, j, k, c, l = 0;
	int uind, vind, wind;   
    double fac, wgt, w = 0.0;
	
	double *R;
	
	if(nrhs != 8) mexErrMsgTxt("Please pass in 8 arguments to the function."
		"We expect it to be in the form [R dRdxi dRdeta dRdzeta] = nrbNurbs3DBasisDerivs(xi, p, q, r, uknot, vknot, wknot, weights)\n");
	
	if(fabs(xi[0]-uknot[u_ctrlPtsNum-1]) < tol) 
		xi[0] = uknot[u_ctrlPtsNum-1] - tol;
    
	if(fabs(xi[1]-vknot[v_ctrlPtsNum-1]) < tol) 
		xi[1] = vknot[v_ctrlPtsNum-1] - tol;
	
	if(fabs(xi[2]-wknot[w_ctrlPtsNum-1]) < tol) 
		xi[2] = wknot[w_ctrlPtsNum-1] - tol;
	
    // Calculate the non-zero univariate B-spline basis functions and first derivatives 
	
	nrbAllBsplineBasisFuns(uspan, xi[0], p, uknot, Nu);
	nrbAllBsplineBasisFuns(vspan, xi[1], q, vknot, Nv);
	nrbAllBsplineBasisFuns(wspan, xi[2], r, wknot, Nw);
	
	// create NURBS approximation
	
	for(k = 0; k<=r; k++)
	{
		wind = wspan - r + k;
		for(j = 0; j <= q; j++)
		{
			vind = vspan - q + j;
			for(i = 0; i <= p; i++)
			{
				uind = uspan - p + i;			
				c = (wind*(v_ctrlPtsNum+1)+vind)*(u_ctrlPtsNum+1)+uind;
				wgt = weight[c];
				
				w += Nu[i] * Nv[j] * Nw[k] * wgt;
			}
		}
	}
	
	// create output 
    
    plhs[0] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
	R = mxGetPr(plhs[0]);
	
	for(k = 0; k<=r; k++)
	{
		wind = wspan - r + k;
		for(j = 0; j <= q; j++)
		{
			vind = vspan - q + j;
			for(i = 0; i <= p; i++)
			{
				uind = uspan - p + i;			
				c = (wind*(v_ctrlPtsNum+1)+vind)*(u_ctrlPtsNum+1)+uind;
				fac = weight[c]/w;
				
				R[l] = Nu[i]*Nv[j]*Nw[k]*fac;
				
				l += 1;
			}
		}
	}
	
	free(Nu);
    free(Nv);
	free(Nw);
}
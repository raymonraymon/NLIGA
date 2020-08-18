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
	* The function will return all related 2D NURBS basis functions and 
	* their first derivatives at point [xi,eta] to matlab.
	* 
	* It can be recalled in form of:
	* [R, dRdxi, dRdeta] = nrbNurbs2DBasisDerivs(xi, p, q, uknot, vknot, weights)
	* Input:
	*   xi            - point in form of [xi eta] for calculation
	*   p,q           - degrees of basis
	*   uknot,vknot   - knot vectors
	*   weights       - vector of weights
	* Output:
	*   R             - all non-zero basis functions
	*   dRdxi,dRdeta  - derivatives
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
	double **Nu_ders = init2DArray(u_ctrlPtsNum+1, p+1);
    double **Nv_ders = init2DArray(v_ctrlPtsNum+1, q+1);
	
	int uspan = nrbFindSpan(u_ctrlPtsNum, p, xi[0], uknot); 
    int vspan = nrbFindSpan(v_ctrlPtsNum, q, xi[1], vknot);
	
	int i, j, c, k = 0;
	int uind, vind;   
    double fac, wgt, w = 0.0, dwdxi = 0.0, dwdeta = 0.0;
	
	double *R, *dRdxi, *dRdeta;
	
	
	if(nrhs != 6) mexErrMsgTxt("Please pass in 6 arguments to the function."
		"We expect it to be in the form [R dRdxi dRdeta] = nrbNurbs2DBasisDerivs(xi, p, q, uknot, vknot, weights)\n");
	
	if(fabs(xi[0]-uknot[u_ctrlPtsNum-1]) < tol) 
		xi[0] = uknot[u_ctrlPtsNum-1] - tol;
    
	if(fabs(xi[1]-vknot[v_ctrlPtsNum-1]) < tol) 
		xi[1] = vknot[v_ctrlPtsNum-1] - tol;

    // Calculate the non-zero univariate B-spline basis functions and first derivatives 
	
	nrbAllBsplineBasisFuns(uspan, xi[0], p, uknot, Nu);
	nrbAllBsplineBasisFuns(vspan, xi[1], q, vknot, Nv);
	
	nrbAllBsplineBasisFunsDers(uspan, xi[0], p, uknot, u_ctrlPtsNum, Nu_ders);
	nrbAllBsplineBasisFunsDers(vspan, xi[1], q, vknot, v_ctrlPtsNum, Nv_ders);

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
            dwdxi  += Nu_ders[1][i] * Nv[j] * wgt;
            dwdeta += Nu[i] * Nv_ders[1][j] * wgt;
        }
    }
	
	// create output
    
    plhs[0] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
    plhs[2] = mxCreateDoubleMatrix(1,noFuncs,mxREAL); 
	
	R      = mxGetPr(plhs[0]);
    dRdxi  = mxGetPr(plhs[1]);
    dRdeta = mxGetPr(plhs[2]);

	for(j = 0; j <= q; j++)
    {
		vind = vspan - q + j;
        for(i = 0; i <= p; i++)
        {      
			uind = uspan - p + i; 	
            c = vind*(u_ctrlPtsNum+1)+uind;            
            fac = weight[c]/(w*w);
            
            R[k] = Nu[i]*Nv[j]*fac*w;
            dRdxi[k] = (Nu_ders[1][i]*Nv[j]*w - Nu[i]*Nv[j]*dwdxi) * fac;
            dRdeta[k] = (Nu[i]*Nv_ders[1][j]*w - Nu[i]*Nv[j]*dwdeta) * fac;
            
            k += 1;
        }
    }
	
	free(Nu);
    free(Nv);
	free2Darray(Nu_ders, (u_ctrlPtsNum+1));
    free2Darray(Nv_ders, (v_ctrlPtsNum+1));
	
}
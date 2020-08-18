#include "nrbNurbs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

int nrbFindSpan(int n, int p, double u, double knot[])
{
	/* 
	* This function finds the knot span.
	* u{i} <= u < (not equal) u{i+1} 
	*/
	
	int low = p;
	int high = n + 1;
	int middle = (low + high) / 2;
	
	if (u>=knot[n+1])
		return n;
	if (u<knot[p])
		return p;   
	
	while ( u<knot[middle] || u>=knot[middle+1] )  
	{
		if ( u<knot[middle] )
			high = middle;
		else
			low = middle;
		middle = (low + high) / 2;
	}
	return middle;
}

void nrbAllBsplineBasisFuns(int i, double u, int p, double knot[], double *N)
{
	/*
	* Calculate all non zero basis functions at point u
	*/
	 
	int j, r;
	double saved, temp;  
	double *lef = (double *)malloc(sizeof(double)*(p+1));
	double *rig = (double *)malloc(sizeof(double)*(p+1)); 

	N[0] = 1.0;
	for (j = 1; j <= p; j++)
	{
		lef[j] = u - knot[i + 1 - j];
		rig[j] = knot[i + j] - u;
		saved = 0.0;
		for (r = 0; r < j; r++)
		{
			temp = N[r] / (rig[r + 1] + lef[j - r]);
			N[r] = saved + rig[r + 1] * temp;
			saved = lef[j - r] * temp;
		}
		N[j] = saved;
	}

	free(lef);
	free(rig);
}

void nrbAllBsplineBasisFunsDers(int i, double u, int p, double knot[], int nd, double **ders)
{
	/*
	* Calculate all non-zero derivatives of the b-spline functions
	*/
	
	int j, r, k, s1, s2, rk, pk, j1, j2;
	double saved, temp, d;	
	double *lef  = (double *)malloc(sizeof(double)*(p+1));
    double *rig  = (double *)malloc(sizeof(double)*(p+1));
	double **ndu  = init2DArray(p+1, p+1);
    double **a    = init2DArray(p+1, p+1);

	ndu[0][0] = 1.0;
	for (j = 1; j <= p; j++)
	{
		lef[j] = u - knot[i + 1 - j];
		rig[j] = knot[i + j] - u;
		saved = 0.0;
		for (r = 0; r < j; r++)
		{
			ndu[j][r] = rig[r + 1] + lef[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];
			ndu[r][j] = saved + rig[r + 1] * temp;
			saved = lef[j - r] * temp;
		}
		ndu[j][j] = saved;
	}

	for (j = 0; j <= p; j++)
		ders[0][j] = ndu[j][p];

	for (r = 0; r <= p; r++)
	{
		s1 = 0; s2 = 1;
		a[0][0] = 1.0;

		for (k = 1; k <= nd; k++)
		{
			d = 0.0;
			rk = r - k;
			pk = p - k;
			if (r >= k)
			{
				a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}
			if (rk >= -1)  j1 = 1;
			else          j1 = -rk;
			if (r - 1 <= pk) j2 = k - 1;
			else          j2 = p - r;

			for (j = j1; j <= j2; j++)
			{
				a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
				d += a[s2][j] * ndu[rk + j][pk];
			}

			if (r <= pk)
			{
				a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
				d += a[s2][k] * ndu[r][pk];
			}
			ders[k][r] = d;
			j = s1; s1 = s2; s2 = j;
		}
	}

	r = p;
	for (k = 1; k <= nd; k++)
	{
		for (j = 0; j <= p; j++)
			ders[k][j] *= r;
		r *= (p - k);
	}

	free(lef); 
    free(rig);
   
    free2Darray(ndu, p+1);
    free2Darray(a, p+1);
}

double** init2DArray(int x, int y)
{
	/*
	* Allocate space for two-dimensional array
	*/
	 
 	double **array = (double **)malloc(x * sizeof(double *));
 	
 	int c;
 	for(c = 0; c < x; c++)
 	{
 		array[c] = (double*)malloc(y * sizeof(double));
 	}
 	
 	return array;
}

void free2Darray(double **array, int x)
{
	/*
	* Free space for two-dimensional array
	*/
	
	int c;
	for(c = 0; c < x; c++)
	{
		free(array[c]);
	}
	free(array);
}



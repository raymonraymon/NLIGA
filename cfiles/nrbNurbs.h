// The following functions could be found in the book 'The NURBS book' L Piegl, W Tiller
/**
* @brief Determine the index of the node interval the parameter u located
* @param n    number of ctrlpts-1
* @param p    degree of Bspline Basis
* @param u    parameter
* @param knot node vector
*
* @return the index 
*/
int nrbFindSpan(int n, int p, double u, double *knot);

/**
* @brief Evaluate all non-zero basis functions
* @param i    the index
* @param u    parameter
* @param p    degree of Bspline Basis
* @param knot node vector
* @param N    all non-zero basis function values,*N = new double[p+1].
*/
void nrbAllBsplineBasisFuns(int i, double u, int p, double *knot, double *N);

/**
* @brief Evaluate all non-zero basis functions and 1-nd derivatives
* @param i    the index
* @param u    parameter
* @param p    degree of Bspline Basis
* @param knot node vector
* @param nd   order of derivatives
*/
void nrbAllBsplineBasisFunsDers(int i, double u, int p, double *knot, int nd, double **ders);

/**
* @brief Allocate space for two-dimensional array
*/
double** init2DArray(int x, int y);

/**
* @brief Free space for two-dimensional array
*/
void free2Darray(double **array, int x);





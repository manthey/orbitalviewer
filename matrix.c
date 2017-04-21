/* This file contains all the routines necessary for manipulating matrices.
 *  The header MATRIX.H should be added to any program which calls these
 *  routines.  Any matrix with a total size of less than 16383 elements can
 *  be manipulated unless the routine states otherwise.        8/16/94-DWM */

#include <math.h>
#include "matrix.h"

void abgtomat(float *r, float alpha, float beta, float gamma)
/* Convert a set of alpha, beta, gammera angles to a 3x3 orthonormal rotation
 *  matrix.
 * Enter: float *r: array of nine values for matrix.
 *        float alpha, beta, gamma: angles of rotation.         1/7/98-DWM */
{
  double ca=cos(alpha), cb=cos(beta), cg=cos(gamma);
  double sa=sin(alpha), sb=sin(beta), sg=sin(gamma);

  r[0] =  cg*ca-sg*cb*sa;  r[1] =  sg*ca+cg*cb*sa;  r[2] = sb*sa;
  r[3] = -cg*sa-sg*cb*ca;  r[4] = -sg*sa+cg*cb*ca;  r[5] = sb*ca;
  r[6] =  sg*sb;           r[7] = -cg*sb;           r[8] = cb;
}

void angtomat(float *r, float ta, float pa, float sa)
/* Convert a set of three angles to a 3x3 orthonormal rotation matrix.
 * Enter: float *r: array of nine values for matrix.
 *        float ta, pa, sa: theta, phi, and psi -- angles of rotation about
 *                          the z, x, and y axes, respectively.11/1/94-DWM */
{
  r[0] = sin(pa)*sin(sa)*sin(ta)+cos(sa)*cos(ta);
  r[1] = sin(pa)*sin(sa)*cos(ta)-cos(sa)*sin(ta);
  r[2] = cos(pa)*sin(sa);  r[3] = cos(pa)*sin(ta);
  r[4] = cos(pa)*cos(ta);  r[5] = -sin(pa);
  r[6] = sin(pa)*cos(sa)*sin(ta)-sin(sa)*cos(ta);
  r[7] = sin(pa)*cos(sa)*cos(ta)+sin(sa)*sin(ta);
  r[8] = cos(pa)*cos(sa);
}

float *cross_product(float *v1, float *v2, float *r)
/* Compute the cross product of two vectors.  The result may replace one of
 *  the source vectors.
 * Enter: float *v1: pointer to array of float values of the first vector.
 *        float *v2: pointer to second vector.
 *        float *r: pointer to location for result.
 * Exit:  float *r: pointer to location of result (same as entered value).
 *                                                             8/20/91-DWM */
{
  float t[3];            /* Uses extra vector to allow r to equal v1 or v2 */

  t[0] = v1[1]*v2[2] - v1[2]*v2[1];
  t[1] = v1[2]*v2[0] - v1[0]*v2[2];
  t[2] = v1[0]*v2[1] - v1[1]*v2[0];
  r[0] = t[0];  r[1] = t[1];  r[2] = t[2];
  return(r);
}

long cubic(float *coef, float *roots)
/* Compute the roots of a cubic equation.  Only real roots are returned.
 *  This works properly with lesser degree equations if the higher order
 *  coefficients are set to zero.
 * Enter: float *coef: coefficients of the cubix equation.  These are of the
 *                     form coef[0]*x^3+coef[1]*x^2+coef[2]*x+coef[3]=0.
 *        float *roots: location to store up to three roots.
 * Exit:  long numroots: number of real roots (0,1,2,or 3).    1/18/97-DWM */
{
  double a, b, c, d, n, o, s, sc, t;
  long neg=0;

  a = coef[0];  b = coef[1];  c = coef[2];  d = coef[3];
  if (!a) {
    if (b) {
      roots[0] = -c+sqrt(fabs(c*c-4*b*d))/(2*b);
      roots[1] = -c-sqrt(fabs(c*c-4*b*d))/(2*b);
      return(2); }
    if (c) {
      roots[0] = -d/c;  return(1); }
    return(0); }
  b /= -a;  c /= -a;  d /= -a;
  n = b*b+3*c;  o = -(2*b*b*b+9*b*c+27*d);
  if (n<0) {
    n = fabs(n);  sc = 0.5*o*pow(n,-1.5);  t = 2*sqrt(n)/3;
    s = log(fabs(sc)+sqrt(sc*sc+1))*((sc>0)-(sc<0))/3;
    roots[0] = b/3-t*0.5*(exp(s)-exp(-s));  return(1); }
  else if (n>0) {
    sc = 0.5*o*pow(n,-1.5);  t = 2*sqrt(n)/3;
    if (sc>=-1 && sc<=1) {
      s = asin(sc)/3;
      roots[0] = t*sin(s-2*PI/3)+b/3;
      roots[1] = t*sin(   s    )+b/3;
      roots[2] = t*sin(s+2*PI/3)+b/3;  return(3); }
    else {
      s = log(fabs(sc)+sqrt(sc*sc-1))/3;
      roots[0] = t*0.5*(exp(s)+exp(-s))+b/3;  return(1); } }
  return(0);
}

float dot_product(float *v1, float *v2)
/* Compute the dot product of two vectors.
 * Enter: float *v1: array of float values representing the first vector.
 *        float *v2: the second vector.
 * Exit:  float dot: dot product answer.                       8/20/91-DWM */
{
  return(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}

void eigen(float *m, float *l)
/* Compute the eigenvalues of a 3x3 matrix.  A pointer to the matrix's array
 *  is passed, not an actual matrix.  The three eigenvalues are stored with
 *  the greatest one first, then in decreasing order.  It is possible that
 *  there is only one real-valued eigenvalue.  In this case, it is computed,
 *  and the other two eigenvalues are stored as zeroes.
 * Enter: float *m: 3x3 matrix to compute eigenvalues for.
 *        float *l: 3x1 array to store resultant values.      10/12/94-DWM */
{
  double b, c, d, n, o, s, sc, t;
  long neg=0;

  b = m[0]+m[4]+m[8];
  c = m[1]*m[3]+m[2]*m[6]+m[5]*m[7]-m[0]*m[4]-m[4]*m[8]-m[0]*m[8];
  d = m[0]*(m[4]*m[8]-m[5]*m[7]) + m[1]*(m[5]*m[6]-m[3]*m[8]) +
      m[2]*(m[3]*m[7]-m[4]*m[6]);
  n = b*b+3*c;  o = -(2*b*b*b+9*b*c+27*d);
  if (n<0) {
    n = fabs(n);  sc = 0.5*o*pow(n,-1.5);  t = 2*sqrt(n)/3;
    s = log(fabs(sc)+sqrt(sc*sc+1))*((sc>0)-(sc<0))/3;
    l[0] = b/3-t*0.5*(exp(s)-exp(-s));
    l[1] = l[2] = 0; }
  else if (!n) {
    l[0]=l[1]=l[2]=0; return; }
  else {
    sc = 0.5*o*pow(n,-1.5);  t = 2*sqrt(n)/3;
    if (sc>=-1 && sc<=1) {
      s = asin(sc)/3;
      l[0] = t*sin(s-2*PI/3)+b/3;  if (l[0]<0)  l[0] = 0;
      l[1] = t*sin(   s    )+b/3;  if (l[1]<0)  l[1] = 0;
      l[2] = t*sin(s+2*PI/3)+b/3;  if (l[2]<0)  l[2] = 0; }
    else {
      s = log(fabs(sc)+sqrt(sc*sc-1))/3;
      l[0] = t*0.5*(exp(s)+exp(-s))+b/3;
      l[1] = l[2] = 0; } }
  if (l[0]<l[1]) { d = l[1];  l[1] = l[0];  l[0] = d; }
  if (l[1]<l[2]) { d = l[1];  l[1] = l[2];  l[2] = d; }
  if (l[0]<l[1]) { d = l[1];  l[1] = l[0];  l[0] = d; }
}

long eigenvect(float *m, float *l, float *ev, long column)
/* Calculate the three unit eigenvectors for a 3x3 matrix.  The eigenvalues
 *  must have already been calculated.  The three vectors are stored in a 3x3
 *  matrix, either as row or column vectors.  Note that the matrices are
 *  arrays, not actual matrix structures.
 * Enter: float *m: pointer to initial 3x3 matrix.
 *        float *l: eigenvalues.
 *        float *ev: location to store eigenvectors.
 *        long column: 0 to store eigenvectors as row vectors (the first
 *                    vector is in ev[0],ev[1],ev[2]).  1 to store as column
 *                    vectors (the first vector is in ev[0],ev[3],ev[6]).
 * Exit:  long failed: 0 if successful, 1 if failed.          10/13/94-DWM */
{
  double emc[10], avg;
  lmatrix *em=lmat(emc);
  long i;

  for (i=0; i<3; i++) {
    em->w = 3;  em->h = 2;
    em->m[0]=m[0]-l[i];  em->m[1]=m[1];       em->m[2]=m[2];
    em->m[3]=m[3];       em->m[4]=m[4]-l[i];  em->m[5]=m[5];
    if (lrowreduce(em)) {
      em->h = 3;
      em->m[0]=m[0]-l[i];  em->m[1]=m[1];       em->m[2]=m[2];
      em->m[3]=m[3];       em->m[4]=m[4]-l[i];  em->m[5]=m[5];
      em->m[6]=m[6];       em->m[7]=m[7];       em->m[8]=m[8]-l[i];
      lrowreduce(em); }
    avg = sqrt(em->m[2]*em->m[2]+em->m[5]*em->m[5]+1);
    if (!avg)
      return(1);
    avg = 1 / avg;
    ev[i*(3-2*column)]              = -em->m[2]*avg;
    ev[i*(3-2*column)+(1+2*column)] = -em->m[5]*avg;
    ev[i*(3-2*column)+(2+4*column)] =           avg; }
  return(0);
}

void labgtomat(double *r, double alpha, double beta, double gamma)
/* Convert a set of alpha, beta, gammera angles to a 3x3 orthonormal rotation
 *  matrix.
 * Enter: double *r: array of nine values for matrix.
 *        double alpha, beta, gamma: angles of rotation.        1/7/98-DWM */
{
  double ca=cos(alpha), cb=cos(beta), cg=cos(gamma);
  double sa=sin(alpha), sb=sin(beta), sg=sin(gamma);

  r[0] =  cg*ca-sg*cb*sa;  r[1] =  sg*ca+cg*cb*sa;  r[2] = sb*sa;
  r[3] = -cg*sa-sg*cb*ca;  r[4] = -sg*sa+cg*cb*ca;  r[5] = sb*ca;
  r[6] =  sg*sb;           r[7] = -cg*sb;           r[8] = cb;
}

void langtomat(double *r, double ta, double pa, double sa)
/* Convert a set of three angles to a 3x3 orthonormal rotation matrix.
 * Enter: double *r: array of nine values for matrix.
 *        double ta, pa, sa: theta, phi, and psi -- angles of rotation about
 *                          the z, x, and y axes, respectively.11/1/94-DWM */
{
  r[0] = sin(pa)*sin(sa)*sin(ta)+cos(sa)*cos(ta);
  r[1] = sin(pa)*sin(sa)*cos(ta)-cos(sa)*sin(ta);
  r[2] = cos(pa)*sin(sa);  r[3] = cos(pa)*sin(ta);
  r[4] = cos(pa)*cos(ta);  r[5] = -sin(pa);
  r[6] = sin(pa)*cos(sa)*sin(ta)-sin(sa)*cos(ta);
  r[7] = sin(pa)*cos(sa)*cos(ta)+sin(sa)*sin(ta);
  r[8] = cos(pa)*cos(sa);
}

double *lcross_product(double *v1, double *v2, double *r)
/* Compute the cross product of two vectors.  The result may replace one of
 *  the source vectors.
 * Enter: double *v1: pointer to array of double values of the first vector.
 *        double *v2: pointer to second vector.
 *        double *r: pointer to location for result.
 * Exit:  double *r: pointer to location of result (same as entered value).
 *                                                             8/20/91-DWM */
{
  double t[3];          /* Uses extra vector to allow r to equal v1 or v2 */

  t[0] = v1[1]*v2[2] - v1[2]*v2[1];
  t[1] = v1[2]*v2[0] - v1[0]*v2[2];
  t[2] = v1[0]*v2[1] - v1[1]*v2[0];
  r[0] = t[0];  r[1] = t[1];  r[2] = t[2];
  return(r);
}

long lcubic(double *coef, double *roots)
/* Compute the roots of a cubic equation.  Only real roots are returned.
 *  This works properly with lesser degree equations if the higher order
 *  coefficients are set to zero.
 * Enter: double *coef: coefficients of the cubix equation.  These are of the
 *                      form coef[0]*x^3+coef[1]*x^2+coef[2]*x+coef[3]=0.
 *        double *roots: location to store up to three roots.
 * Exit:  long numroots: number of real roots (0,1,2,or 3).    1/18/97-DWM */
{
  double a, b, c, d, n, o, s, sc, t;
  long neg=0;

  a = coef[0];  b = coef[1];  c = coef[2];  d = coef[3];
  if (!a) {
    if (b) {
      roots[0] = -c+sqrt(fabs(c*c-4*b*d))/(2*b);
      roots[1] = -c-sqrt(fabs(c*c-4*b*d))/(2*b);
      return(2); }
    if (c) {
      roots[0] = -d/c;  return(1); }
    return(0); }
  b /= -a;  c /= -a;  d /= -a;
  n = b*b+3*c;  o = -(2*b*b*b+9*b*c+27*d);
  if (n<0) {
    n = fabs(n);  sc = 0.5*o*pow(n,-1.5);  t = 2*sqrt(n)/3;
    s = log(fabs(sc)+sqrt(sc*sc+1))*((sc>0)-(sc<0))/3;
    roots[0] = b/3-t*0.5*(exp(s)-exp(-s));  return(1); }
  else if (n>0) {
    sc = 0.5*o*pow(n,-1.5);  t = 2*sqrt(n)/3;
    if (sc>=-1 && sc<=1) {
      s = asin(sc)/3;
      roots[0] = t*sin(s-2*PI/3)+b/3;
      roots[1] = t*sin(   s    )+b/3;
      roots[2] = t*sin(s+2*PI/3)+b/3;  return(3); }
    else {
      s = log(fabs(sc)+sqrt(sc*sc-1))/3;
      roots[0] = t*0.5*(exp(s)+exp(-s))+b/3;  return(1); } }
  return(0);
}

double ldot_product(double *v1, double *v2)
/* Compute the dot product of two vectors.
 * Enter: double *v1: array of double values representing the first vector.
 *        double *v2: the second vector.
 * Exit:  double dot: dot product answer.                      8/20/91-DWM */
{
  return(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}

long least(matrix *a, float *coef, float weight)
/* Adds a set of coefficients to a matrix for use in a least squares process.
 *  This only uses coefficients in the "upper" triangle of the matrix.  The
 *  matsym function must be used to copy this upper triangle into the lower
 *  triangle.  There must be the same number of coefficients as there are
 *  column in the matrix, and there must be one more column than the number
 *  of rows.  The coefficeints are the coefficients of the equation to solve.
 *  For example, if the equations is of the form Ax+By+Cz=D, to solve for x,
 *  y, and z, the coefficients are {A,B,C,D}.  To perform a complete least
 *  squares solution: (1) set a matrix to the appropriate size.  In the
 *  example, width=4, height=3.  (2) call matzero.  (3) call least with the
 *  appropriate coefficient array and weight for each data point.  (4) call
 *  matsym.  (5) call rowreduce.  The answers wil be located in the last
 *  column of the matrix.
 * Enter: matrix *a: pointer to least squares matrix.
 *        float *coef: pointer to array of coefficients to use.  This array
 *                      is not modified.
 *        float weight: weighting factor for these coefficients.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, j;

  if (a->w<=a->h)
    return(1);
  for (j=0; j<a->h; j++)
    for (i=j; i<a->w; i++)
      a->m[j*a->w+i] += coef[i]*coef[j]*weight;
  return(0);
}

void leigen(double *m, double *l)
/* Compute the eigenvalues of a 3x3 matrix.  A pointer to the matrix's array
 *  is passed, not an actual matrix.  The three eigenvalues are stored with
 *  the greatest one first, then in decreasing order.  It is possible that
 *  there is only one real-valued eigenvalue.  In this case, it is computed,
 *  and the other two eigenvalues are stored as zeroes.
 * Enter: double *m: 3x3 matrix to compute eigenvalues for.
 *        double *l: 3x1 array to store resultant values.     10/12/94-DWM */
{
  double b, c, d, n, o, s, sc, t;
  long neg=0;

  b = m[0]+m[4]+m[8];
  c = m[1]*m[3]+m[2]*m[6]+m[5]*m[7]-m[0]*m[4]-m[4]*m[8]-m[0]*m[8];
  d = m[0]*(m[4]*m[8]-m[5]*m[7]) + m[1]*(m[5]*m[6]-m[3]*m[8]) +
      m[2]*(m[3]*m[7]-m[4]*m[6]);
  n = b*b+3*c;  o = -(2*b*b*b+9*b*c+27*d);
  if (n<0) {
    n = fabs(n);  sc = 0.5*o*pow(n,-1.5);  t = 2*sqrt(n)/3;
    s = log(fabs(sc)+sqrt(sc*sc+1))*((sc>0)-(sc<0))/3;
    l[0] = b/3-t*0.5*(exp(s)-exp(-s));
    l[1] = l[2] = 0; }
  else if (!n) {
    l[0]=l[1]=l[2]=0; return; }
  else {
    sc = 0.5*o*pow(n,-1.5);  t = 2*sqrt(n)/3;
    if (sc>=-1 && sc<=1) {
      s = asin(sc)/3;
      l[0] = t*sin(s-2*PI/3)+b/3;  if (l[0]<0)  l[0] = 0;
      l[1] = t*sin(   s    )+b/3;  if (l[1]<0)  l[1] = 0;
      l[2] = t*sin(s+2*PI/3)+b/3;  if (l[2]<0)  l[2] = 0; }
    else {
      s = log(fabs(sc)+sqrt(sc*sc-1))/3;
      l[0] = t*0.5*(exp(s)+exp(-s))+b/3;
      l[1] = l[2] = 0; } }
  if (l[0]<l[1]) { d = l[1];  l[1] = l[0];  l[0] = d; }
  if (l[1]<l[2]) { d = l[1];  l[1] = l[2];  l[2] = d; }
  if (l[0]<l[1]) { d = l[1];  l[1] = l[0];  l[0] = d; }
}

long leigenvect(double *m, double *l, double *ev, long column)
/* Calculate the three unit eigenvectors for a 3x3 matrix.  The eigenvalues
 *  must have already been calculated.  The three vectors are stored in a 3x3
 *  matrix, either as row or column vectors.  Note that the matrices are
 *  arrays, not actual matrix structures.
 * Enter: double *m: pointer to initial 3x3 matrix.
 *        double *l: eigenvalues.
 *        double *ev: location to store eigenvectors.
 *        long column: 0 to store eigenvectors as row vectors (the first
 *                    vector is in ev[0],ev[1],ev[2]).  1 to store as column
 *                    vectors (the first vector is in ev[0],ev[3],ev[6]).
 * Exit:  long failed: 0 if successful, 1 if failed.          10/13/94-DWM */
{
  double emc[10], avg;
  lmatrix *em=lmat(emc);
  long i;

  for (i=0; i<3; i++) {
    em->w = 3;  em->h = 2;
    em->m[0]=m[0]-l[i];  em->m[1]=m[1];       em->m[2]=m[2];
    em->m[3]=m[3];       em->m[4]=m[4]-l[i];  em->m[5]=m[5];
    if (lrowreduce(em)) {
      em->h = 3;
      em->m[0]=m[0]-l[i];  em->m[1]=m[1];       em->m[2]=m[2];
      em->m[3]=m[3];       em->m[4]=m[4]-l[i];  em->m[5]=m[5];
      em->m[6]=m[6];       em->m[7]=m[7];       em->m[8]=m[8]-l[i];
      lrowreduce(em); }
    avg = sqrt(em->m[2]*em->m[2]+em->m[5]*em->m[5]+1);
    if (!avg)
      return(1);
    avg = 1 / avg;
    ev[i*(3-2*column)]              = -em->m[2]*avg;
    ev[i*(3-2*column)+(1+2*column)] = -em->m[5]*avg;
    ev[i*(3-2*column)+(2+4*column)] =           avg; }
  return(0);
}

long lleast(lmatrix *a, double *coef, double weight)
/* Adds a set of coefficients to the upper triangle of a matrix for use in a
 *  least squares process.  See the least function for a complete
 *  description.
 * Enter: lmatrix *a: pointer to least squares matrix.
 *        double *coef: pointer to array of coefficients to use.  This array
 *                      is not modified.
 *        double weight: weighting factor for these coefficients.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, j;

  if (a->w<=a->h)
    return(1);
  for (j=0; j<a->h; j++)
    for (i=j; i<a->w; i++)
      a->m[j*a->w+i] += coef[i]*coef[j]*weight;
  return(0);
}

lmatrix *lmat(double *a)
/* Set a pointer to an lmatrix to point to an array of double values.  The
 *  array of doubles must be at least two larger than the array within the
 *  matrix.  Width and height are arbitrarily set to 3.
 * Enter: double *a: pointer to double array for matrix.
 * Exit:  lmatrix *c: pointer to lmatrix structure.           10/14/94-DWM */
{
  lmatrix *c=(lmatrix *)a;

  c->m = a+2;
  c->w = c->h = 3;
  return(c);
}

long lmatadd(double c, lmatrix *a, double d, lmatrix *b, lmatrix *caplusdb)
/* Computes the matrix function c*a+d*b, where a and b are both n x m
 *  matrices, and c and d are both scalar values.
 * Enter: double c: value to scale matrix a by.
 *        lmatrix *a: pointer to matrix to add.
 *        double d: value to scale matrix b by.
 *        lmatrix *b: pointer to matrix to add.
 *        lmatrix *caplusdb: pointer to resultant matrix.  This may be one of
 *                           the calling matrices.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, len=a->w*a->h;

  if (b->w!=a->w || a->h!=b->h)
    return(1);
  for (i=0; i<len; i++)
    caplusdb->m[i] = c*a->m[i]+d*b->m[i];
  caplusdb->w = a->w;
  caplusdb->h = a->h;
  return(0);
}

double lmatdet(lmatrix *a)
/* Find the determinant of a square matrix.  There is no error return for
 *  non-square matrices.  Maximum size is 5x5.  To change the maximum size,
 *  the size of the subm matrix must be changed.
 * Enter: lmatrix *a: pointer to matrix to calculate the determinant for.
 * Exit:  double det: determinant of the supplied matrix.      8/16/94|DWM */
{
  long i, j, k, s=a->w-1, size=a->w;
  double d=0;
  char array[16*8+8];                                   /* (n-1)*(n-1)*8+8 */
  lmatrix *subm=lmat((double *)array);

  if (size==2)
    return((double)a->m[0]*a->m[3]-(double)a->m[1]*a->m[2]);
  for (k=0; k<size; k++) {
    subm->w = subm->h = s;
    for (i=0; i<s; i++)
      for (j=0; j<s; j++)
        subm->m[j*s+i] = a->m[j*size+i+(i>=s-k)];
    d = d + (1-(k%2)*2) * lmatdet(subm) * a->m[s*size+(s-k)]; }
  return(d);
}

void lmatdup(lmatrix *a, lmatrix *adup)
/* Duplicate the matrix a.
 * Enter: lmatrix *a: original matrix.
 *        lmatrix *adup: location to place copy of a.          8/16/94-DWM */
{
  adup->w = a->w;  adup->h = a->h;
  memcpy(adup->m, a->m, sizeof(double)*a->w*a->h);
}

void lmatident(lmatrix *i, long size)
/* Create a size x size identity matrix.
 * Enter: lmatrix *i: pointer to matrix.
 *        long size: size of matrix to create.                  8/16/94-DWM */
{
  long j;

  i->w = i->h = size;
  for (j=1; j<size*size-1; j++)
    i->m[j] = 0;
  for (j=0; j<size; j++)
    i->m[j*(size+1)] = 1;
}

long lmatinv(lmatrix *a, lmatrix *ainv)
/* Invert a square matrix.  The input matrix is changed to the identity
 *  matrix.
 * Enter: lmatrix *a: point to matrix to invert.
 *        lmatrix *ainv: location to store inverse matrix.  This can not be
 *                        the calling matrix.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  if (a->w!=a->h)
    return(1);
  lmatident(ainv, a->w);
  return(lrowreduce2(a, ainv));
}

long lmatmul(lmatrix *a, lmatrix *b, lmatrix *amulb)
/* Computes the matrix function a*b, where a is an n x m matrix and b is an
 *   p x n matrix.  The result is a p x m matrix.
 * Enter: lmatrix *a, *b: pointer to matrices to multiply.
 *        lmatrix *amulb: pointer to resultant matrix.  This can not be one
 *                        of the calling matrices.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, j, k, w, h, s=a->w;

  if (a->w!=b->h)
    return(1);
  w = amulb->w = b->w;
  h = amulb->h = a->h;
  for (j=0; j<h; j++)
    for (i=0; i<w; i++)
      for (k=amulb->m[j*w+i]=0; k<s; k++)
        amulb->m[j*w+i] += a->m[j*s+k]*b->m[k*w+i];
  return(0);
}

void lmatscale(lmatrix *a, double scalar, lmatrix *result)
/* Multiply a matrix by a scalar value.
 * Enter: lmatrix *a: pointer to matrix.
 *        double scalar: scalar multiplier.
 *        lmatrix *result: location to store result, may be a. 8/16/94-DWM */
{
  long i, len=a->w*a->h;

  for (i=0; i<len; i++)
    result->m[i] = scalar*a->m[i];
  result->w = a->w;
  result->h = a->h;
}

long lmatsub(lmatrix *a, lmatrix *b, lmatrix *aminusb)
/* Computes the matrix function a-b, where a and b are both n x m matrices.
 * Enter: lmatrix *a, *b: pointer to matrices to subtract.
 *        lmatrix *aminusb: pointer to resultant matrix.  This may be one of
 *                        the calling matrices.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, len=a->w*a->h;

  if (b->w*b->h!=len)
    return(1);
  for (i=0; i<len; i++)
    aminusb->m[i] = a->m[i]-b->m[i];
  aminusb->w = a->w;
  aminusb->h = a->h;
  return(0);
}

long lmatsym(lmatrix *a)
/* Make the specified matrix symmetrical.  This is done by copying the upper
 *  triangle to the lower triangle.  The matrix must be at least as wide as
 *  it is tall.
 * Enter: lmatrix *a: pointer to matrix to make symmetric.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, j;

  if (a->w<a->h)
    return(1);
  for (j=1; j<a->h; j++)
    for (i=0; i<j; i++)
      a->m[j*a->w+i] = a->m[i*a->w+j];
  return(0);
}

void lmattoabg(double *r, double *abg)
/* Convert a 3x3 orthonormal matrix to a set of angles.  The angles are in
 *  the order alpha, beta, gammera
 * Enter: double *r: pointer to array of nine values.
 *        double *abg: pointer to array to store 3 angles.      1/8/98-DWM */
{
  abg[1] = acos(r[8]);
  abg[2] = atan2(r[6], -r[7]);
  abg[0] = atan2(r[2], r[5]);
  if (fabs(sin(abg[2])*cos(abg[0])+cos(abg[2])*cos(abg[1])*sin(abg[0])-r[1])>
      1e-6) {
    abg[1] = -acos(r[8]);
    abg[2] = atan2(-r[6], r[7]);
    abg[0] = atan2(-r[2], -r[5]); }
}

void lmattoang(double *r, double *a)
/* Convert a 3x3 orthonormal matrix to a set of angles.  See the Derive file
 *  MTOSIN.MTH for the exact form of these angles.  The angles are in the
 *  order theta, phi, psi.
 * Enter: double *r: pointer to array of nine values.
 *        double *a: pointer to array to store 3 angles.       11/1/94-DWM */
{
  double temp;

  temp = -r[5];  if (fabs(temp)>1)  temp = (temp>0)-(temp<0);
  a[1] = asin(temp);
  if (fabs(cos(a[1]))>1e-30) {
    temp = r[3]/cos(a[1]);  if (fabs(temp)>1)  temp = (temp>0)-(temp<0);
    a[0] = asin(temp);
    if (cos(a[0])*cos(a[1])*r[4]<0)
      a[0] = PI-a[0];
    temp = r[2]/cos(a[1]);  if (fabs(temp)>1)  temp = (temp>0)-(temp<0);
    a[2] = asin(temp);
    if (cos(a[2])*cos(a[1])*r[8]<0)
      a[2] = PI-a[2];
    if (fabs(sin(a[1])*sin(a[2])*sin(a[0])+cos(a[2])*cos(a[0])-r[0])>1e-30 ||
        fabs(sin(a[1])*sin(a[2])*cos(a[0])-cos(a[2])*sin(a[0])-r[1])>1e-30) {
      a[0] = PI+a[0];  a[1] = PI-a[1];  a[2] = PI+a[2]; } }
  else {
    a[2] = 0;
    temp = r[0];  if (fabs(temp)>1)  temp = (temp>0)-(temp<0);
    a[0] = acos(temp);
    if (fabs(-sin(a[0])-r[1])>1e-30)
      a[0] *= -1;  }
}

void lmattoypr(double *r, double *a)
/* Convert a 3x3 orthonormal matrix to a set of yaw, pitch, and roll angles.
 * Enter: double *r: pointer to array of nine values.
 *        double *a: pointer to array to store 3 angles.       9/10/95-DWM */
{
  if (r[5])  a[2] = atan2(r[2], r[5]);
  else       a[2] = PI*0.5*((r[2]>0)-(r[2]<0));
  if (r[7])  a[0] = atan2(r[6], -r[7]);
  else       a[0] = PI*0.5*((r[6]>0)-(r[7]<0));
  a[1] = asin(r[8]);
  if (cos(a[1])*sin(a[0])*r[6]<0)  a[1] = PI-a[1];
}

void lmattrans(lmatrix *a, lmatrix *atrans)
/* Compute the transpose of a matrix.
 * Enter: lmatrix *a: pointer to matrix.
 *        lmatrix *atrans: location to store result, can not be a.
 *                                                             8/16/94|DWM */
{
  long i, j;

  atrans->w = a->h;
  atrans->h = a->w;
  for (j=0; j<a->w; j++)
    for (i=0; i<a->h; i++)
      atrans->m[j*atrans->w+i] = a->m[i*a->w+j];
}

void lmatzero(lmatrix *a)
/* Zero a matrix.
 * Enter: lmatrix *a: pointer to matrix to zero.               10/5/94-DWM */
{
  long i, s=a->w*a->h;

  for (i=0; i<s; i++)
    a->m[i] = 0;
}

double *lqconj(double *a, double *astar)
/* Conjugate a quaternion.  The quarternions are 3 dimensional.  The format
 *  is (s, v1, v2, v3), where s is the scalar part of the quarternion, and v1
 *  through v3 is the vector part.
 * Enter: double *a: array holding initial quarternion.
 *        double *astar: array to place result.
 * Exit:  double *astar: same as input.                       10/24/94-DWM */
{
  astar[0] = a[0];  astar[1] = -a[1];  astar[2] = -a[2];  astar[3] = -a[3];
  return(astar);
}

double lqdot(double *a, double *b)
/* Take the dot product of two quaternions.  The quarternions are 3
 *  dimensional.  The format is (s, v1, v2, v3), where s is the scalar part
 *  of the quarternion, and v1 through v3 is the vector part.
 * Enter: double *a, *b: array holding the two quarternions.
 * Exit:  double adotb: dot product.                          10/24/94-DWM */
{
  return(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]);
}

double *lqmul(double *a, double *b, double *amulb)
/* Multiply two quaternions.  The quarternions are 3 dimensional.  The format
 *  is (s, v1, v2, v3), where s is the scalar part of the quarternion, and v1
 *  through v3 is the vector part.
 * Enter: double *a, *b: array holding the two multipicands.
 *        double *amulb: array to place result.  This may be one of the
 *                       calling quarternions.
 * Exit:  double *amulb: same as input.                       10/24/94-DWM */
{
  double ab[4];

  ab[0] = a[0]*b[0] - (a[1]*b[1]+a[2]*b[2]+a[3]*b[3]);
  ab[1] = a[0]*b[1] + b[0]*a[1] + a[2]*b[3]-a[3]*b[2];
  ab[2] = a[0]*b[2] + b[0]*a[2] + a[3]*b[1]-a[1]*b[3];
  ab[3] = a[0]*b[3] + b[0]*a[3] + a[1]*b[2]-a[2]*b[1];
  amulb[0] = ab[0];  amulb[1] = ab[1];  amulb[2] = ab[2];  amulb[3] = ab[3];
  return(amulb);
}

double lqnorm(double *a)
/* Take the normal of a quaternion -- ºaº.  The quarternion is 3 dimensional.
 *  The format is (s, v1, v2, v3), where s is the scalar part of the
 *  quarternion, and v1 through v3 is the vector part.
 * Enter: double *a: quarternion to take the normal of.
 * Exit:  double anorm: norm.                                 10/24/94-DWM */
{
  return(sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]+a[3]*a[3]));
}

void lqtomat(double *q, double *rot)
/* Convert a quaternion to a rotation matrix.  The quarternion is 3
 *  dimensional.  The format is (s, v1, v2, v3), where s is the scalar part
 *  of the quarternion, and v1 through v3 is the vector part.
 * Enter: double *q: quarternion to convert.
 *        double *rot: array of nine doubles to store the resulting matrix.
 *                                                            10/18/97-DWM */
{
  rot[0] = q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  rot[1] = 2*(q[1]*q[2]-q[0]*q[3]);
  rot[2] = 2*(q[1]*q[3]+q[0]*q[2]);
  rot[3] = 2*(q[1]*q[2]+q[0]*q[3]);
  rot[4] = q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  rot[5] = 2*(q[2]*q[3]-q[0]*q[1]);
  rot[6] = 2*(q[1]*q[3]-q[0]*q[2]);
  rot[7] = 2*(q[2]*q[3]+q[0]*q[1]);
  rot[8] = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
}

double *lqunit(double *a, double *aunit)
/* Make a unit quaternion.  The quarternion is 3 dimensional.  The format is
 *  (s, v1, v2, v3), where s is the scalar part of the quarternion, and v1
 *  through v3 is the vector part.
 * Enter: double *a: quarternion to unitize.
 *        double *aunit: array to place result.  This may be a.
 * Exit:  double *aunit: same as input.                       10/24/94-DWM */
{
  double norm=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]+a[3]*a[3]);

  if (!norm)
    aunit[0] = aunit[1] = aunit[2] = aunit[3] = 0;
  else {
    aunit[0] = a[0]/norm;  aunit[1] = a[1]/norm;
    aunit[2] = a[2]/norm;  aunit[3] = a[3]/norm; }
  return(aunit);
}

void lrotate_fit(double *XYZ, double *xyz, long fitstart, long fitnum,
           long modstart, long modnum, long scale, double *rot, double *tran)
/* See rotate_fit for details.  The actual code is identical; only the
 *  function definition is different.
 * Enter: double *XYZ: array containing reference coordinates to fit to.
 *        double *xyz: array containing data to be fitted.  If the x
 *                     coordinate of either XYZ or xyz is greater than 1e20,
 *                     that data point is unused.
 *        long fitstart: starting node number for fit (usually 0).
 *        long fitnum: number of nodes to use in fit.
 *        long modstart: starting node number for modifying data.
 *        long modnum: number of nodes to fit and modify.
 *        long scale: 0 to only rotate, 1 to scale and rotate, +2 to only
 *                    compute rotation and translation values.
 *        double *rot: location to store rotation matrix (9 values), or 0 for
 *                     don't store.  Scale is also in this matrix
 *        double *tran: location to store traslation (3 values), or 0 for
 *                      don't store.                          10/12/94-DWM */
{
  double x[3]={0,0,0}, X[3]={0,0,0}, y[3], Y[3], l[3], err, best=1e300;
  long i, j, total=0, besti;
  double cc[132], t[3], sx2=0, c, tr;
  lmatrix *h=lmat(cc), *ht=lmat(cc+11), *hth=lmat(cc+22), *hht=lmat(cc+33);
  lmatrix *r=lmat(cc+44), *ut=lmat(cc+55), *v=lmat(cc+66), *s=lmat(cc+77);
  lmatrix *t1=lmat(cc+88), *t2=lmat(cc+99), *d=lmat(cc+110), *u=lmat(cc+121);
  int mul[]={1,1,1,1,1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1};

  for (i=fitstart; i<fitnum; i++)
    if (xyz[i*3]<1e20 && XYZ[i*3]<1e20) {
      x[0] += xyz[i*3];  x[1] += xyz[i*3+1];  x[2] += xyz[i*3+2];
      X[0] += XYZ[i*3];  X[1] += XYZ[i*3+1];  X[2] += XYZ[i*3+2];
      total++; }
  x[0] /= total;  x[1] /= total;  x[2] /= total;
  X[0] /= total;  X[1] /= total;  X[2] /= total;
  for (i=fitstart; i<fitnum; i++)
    if (xyz[i*3]<1e20 && XYZ[i*3]<1e20)
      for (j=0; j<3; j++)
        sx2 += (xyz[i*3+j]-x[j])*(xyz[i*3+j]-x[j]);
  lmatzero(h);
  for (i=fitstart; i<fitnum; i++)
    if (xyz[i*3]<1e20 && XYZ[i*3]<1e20) {
      y[0] = xyz[i*3]  -x[0];  Y[0] = XYZ[i*3]  -X[0];
      y[1] = xyz[i*3+1]-x[1];  Y[1] = XYZ[i*3+1]-X[1];
      y[2] = xyz[i*3+2]-x[2];  Y[2] = XYZ[i*3+2]-X[2];
      h->m[0]+=y[0]*Y[0];  h->m[1]+=y[0]*Y[1];  h->m[2]+=y[0]*Y[2];
      h->m[3]+=y[1]*Y[0];  h->m[4]+=y[1]*Y[1];  h->m[5]+=y[1]*Y[2];
      h->m[6]+=y[2]*Y[0];  h->m[7]+=y[2]*Y[1];  h->m[8]+=y[2]*Y[2]; }
  lmattrans(h, ht);
  lmatmul(ht, h, hth);
  leigen(hth->m, l);
  leigenvect(hth->m, l, v->m, 1);                   /* Column eigenvectors */
  lmatmul(h, ht, hht);
  leigen(hht->m, l);
  leigenvect(hht->m, l, ut->m, 0);                     /* Row eigenvectors */
  lmattrans(ut, u);
  lmatzero(d);  lmatzero(s);
  d->m[0] = sqrt(l[0]);  d->m[4] = sqrt(l[1]);  d->m[8] = sqrt(l[2]);
  for (i=0; i<8; i++) {
    s->m[0]=mul[i*3];  s->m[4]=mul[i*3+1];  s->m[8]=mul[i*3+2];
    lmattrans(v, t1);             /* t1 = v`         */
    lmatmul(s, t1, t2);           /* t2 = s.v`       */
    lmatmul(d, t2, t1);           /* t1 = d.s.v`     */
    lmatmul(u, t1, t2);           /* t2 = u.d.s.v`   */
    lmatadd(1, t2, -1, h, t2);    /* t2 = u.d.s.v`-h */
    for (j=0, err=0; j<9; j++)
      err += t2->m[j]*t2->m[j];
    if (err<best) {
      best = err;
      besti = i; } }
  s->m[0]=mul[besti*3];  s->m[4]=mul[besti*3+1];  s->m[8]=mul[besti*3+2];
  lmatmul(s, ut, t1);
  lmatmul(v, t1, r);
  tr = d->m[0]+d->m[4]+d->m[8];
  if (lmatdet(r)<0) {
    tr -= d->m[8]*2;
    s->m[8] *= -1;
    lmatmul(s, ut, t1);
    lmatmul(v, t1, r); }
  c = tr/sx2;
  if (scale&1)
    for (i=0; i<9; i++)
      r->m[i] *= c;
  t[0] = X[0]-x[0]*r->m[0]-x[1]*r->m[1]-x[2]*r->m[2];
  t[1] = X[1]-x[0]*r->m[3]-x[1]*r->m[4]-x[2]*r->m[5];
  t[2] = X[2]-x[0]*r->m[6]-x[1]*r->m[7]-x[2]*r->m[8];
  if (rot)   memcpy(rot, r->m, 9*sizeof(double));
  if (tran)  memcpy(tran, t, 3*sizeof(double));
  if (scale>=2)  return;
  for (i=modstart; i<modnum; i++)
    if (xyz[i*3]<1e20) {
      y[0] = xyz[i*3];  y[1] = xyz[i*3+1];  y[2] = xyz[i*3+2];
      xyz[i*3]   = y[0]*r->m[0]+y[1]*r->m[1]+y[2]*r->m[2]+t[0];
      xyz[i*3+1] = y[0]*r->m[3]+y[1]*r->m[4]+y[2]*r->m[5]+t[1];
      xyz[i*3+2] = y[0]*r->m[6]+y[1]*r->m[7]+y[2]*r->m[8]+t[2]; }
}

long lrowreduce(lmatrix *a)
/* Row reduce a matrix.  The matrix must be at least as wide as it is high.
 *  The main diagonal is changed to 1, and all other elements in the main
 *  square are changed to zero.  If the process fails due to an undetermined
 *  value, the matrix is only partially row reduced.
 * Enter: lmatrix *a: matrix to row reduce.  The matrix is changed.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, j, k, w=a->w, h=a->h;
  double temp;

  if (w<h)
    return(1);
  for (i=0; i<h; i++) {
    for (j=i; j<h; j++)
      if (a->m[j*w+i])
        break;
    if (j==h)
      return(1);
    if (j!=i)
      for (k=i; k<w; k++) {
        temp = a->m[i*w+k];
        a->m[i*w+k] = a->m[j*w+k];
        a->m[j*w+k] = temp;  }
    for (j=w-1; j>=i; j--)
      a->m[i*w+j] /= a->m[i*w+i];
    for (j=0; j<h; j++)
      if (j!=i) if (temp=a->m[j*w+i])
        for (k=i; k<w; k++)
          a->m[j*w+k] -= a->m[i*w+k]*temp; }
  return(0);
}

long lrowreduce2(lmatrix *a, lmatrix *b)
/* Row reduce matrix a, modifying matrix b in the same manner.  Matrix a must
 *  be at least as wide as it is high, and matrix b must be the same height
 *  as matrix a.  The main diagonal of a is changed to 1, and all other
 *  elements in the main square are changed to zero.  If the process fails
 *  due to an undetermined value, the matrices are only partially row
 *  reduced.
 * Enter: lmatrix *a: matrix to row reduce.  The matrix is changed.
 *        lmatrix *b: matrix to modify in the same manner as matrix a.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, j, k, w=a->w, h=a->h, w2=b->w, iw, iw2;
  double temp;

  if (w<h || b->h!=h)
    return(1);
  for (i=iw=iw2=0; i<h; i++, iw+=w, iw2+=w2) {
    for (j=i; j<h; j++)
      if (a->m[j*w+i])
        break;
    if (j==h)
      return(1);
    if (j!=i) {
      for (k=i; k<w; k++) {
        temp = a->m[iw+k];
        a->m[iw+k] = a->m[j*w+k];
        a->m[j*w+k] = temp;  }
      for (k=0; k<w2; k++) {
        temp = b->m[iw2+k];
        b->m[iw2+k] = b->m[j*w2+k];
        b->m[j*w2+k] = temp;  } }
    temp = 1/a->m[iw+i];
    a->m[iw+i] = 1;
    for (j=i+1; j<w; j++)
      a->m[iw+j] *= temp;
    for (j=0; j<w2; j++)
      b->m[iw2+j] *= temp;
    for (j=0; j<h; j++)
      if (j!=i) if (temp=a->m[j*w+i]) {
        a->m[j*w+i] = 0;
        for (k=i+1; k<w; k++)
          a->m[j*w+k] -= a->m[iw+k]*temp;
        for (k=0; k<w2; k++)
          b->m[j*w2+k] -= b->m[iw2+k]*temp; } }
  return(0);
}

double *lunit(double *v, double *r)
/* Compute the unit vector of the given vector.  If the vector is the zero
 *  vector, the zero vector is returned.
 * Enter: double *v: vector to compute unit vector from.
 *        double *r: location to place resultant vector.
 * Exit:  double *r: location of resultant vector (same as entered value).
 *                                                             8/20/91|DWM */
{
  double len=0;
  short i;

  for (i=0; i<3; i++)
    len += v[i] * v[i];
  len = sqrt(len);
  if (len)
    for (i=0; i<3; i++)
      r[i] = v[i] / len;
  else
    r[0] = r[1] = r[2] = 0;
  return(r);
}

void lyprtomat(double *rot, double y, double p, double r)
/* Convert a set of roll, pitch, and yaw angles to a 3x3 orthonormal rotation
 *  matrix.  This is actually the matrix
 *   [  cos(r) sin(r) 0 ] [ 1   0       0    ] [ cos(y)  sin(y) 0 ]
 *   [ -sin(r) cos(r) 0 ] [ 0 cos(p) -sin(p) ] [   0       0    1 ]
 *   [    0      0    1 ] [ 0 sin(p)  cos(p) ] [ sin(y) -cos(y) 0 ]
 * Enter: double *rot: array of nine values for matrix.
 *        double r: roll about the x axis.
 *        double p: pitch about the y axis.
 *        double y: yaw about the z axis.                      9/10/95-DWM */
{
  rot[0] = cos(r)*cos(y)-sin(p)*sin(r)*sin(y);
  rot[1] = sin(p)*sin(r)*cos(y)+cos(r)*sin(y);  rot[2] = cos(p)*sin(r);
  rot[3] = -sin(r)*cos(y)-sin(p)*cos(r)*sin(y);
  rot[4] = sin(p)*cos(r)*cos(y)-sin(r)*sin(y);  rot[5] = cos(p)*cos(r);
  rot[6] = cos(p)*sin(y);  rot[7] = -cos(p)*cos(y);  rot[8] = sin(p);
}

matrix *mat(float *a)
/* Set a pointer to a matrix to point to an array of float values.  The array
 *  of floats must be at least four larger than the array within the matrix.
 *  Width and height are arbitrarily set to 3.
 * Enter: float *a: pointer to float array for matrix.
 * Exit:  matrix *c: pointer to matrix structure.             10/14/94-DWM */
{
  matrix *c=(matrix *)a;

  c->m = a+4;
  c->w = c->h = 3;
  return(c);
}

long matadd(float c, matrix *a, float d, matrix *b, matrix *caplusdb)
/* Computes the matrix function c*a+d*b, where a and b are both n x m
 *  matrices, and c and d are both scalar values.
 * Enter: float c: value to scale matrix a by.
 *        matrix *a: pointer to matrix to add.
 *        float d: value to scale matrix b by.
 *        matrix *b: pointer to matrix to add.
 *        matrix *caplusdb: pointer to resultant matrix.  This may be one of
 *                          the calling matrices.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, len=a->w*a->h;

  if (b->w!=a->w || a->h!=b->h)
    return(1);
  for (i=0; i<len; i++)
    caplusdb->m[i] = c*a->m[i]+d*b->m[i];
  caplusdb->w = a->w;
  caplusdb->h = a->h;
  return(0);
}

double matdet(matrix *a)
/* Find the determinant of a square matrix.  There is no error return for
 *  non-square matrices.  Maximum size is 5x5.  To change the maximum size,
 *  the size of the subm matrix must be changed.
 * Enter: matrix *a: pointer to matrix to calculate the determinant for.
 * Exit:  double det: determinant of the supplied matrix.      8/16/94|DWM */
{
  long i, j, k, s=a->w-1, size=a->w;
  double d=0;
  char array[16*4+8];                                   /* (n-1)*(n-1)*4+8 */
  matrix *subm=mat((float *)array);

  if (size==2)
    return((double)a->m[0]*a->m[3]-(double)a->m[1]*a->m[2]);
  for (k=0; k<size; k++) {
    subm->w = subm->h = s;
    for (i=0; i<s; i++)
      for (j=0; j<s; j++)
        subm->m[j*s+i] = a->m[j*size+i+(i>=s-k)];
    d = d + (1-(k%2)*2) * matdet(subm) * a->m[s*size+(s-k)]; }
  return(d);
}

void matdup(matrix *a, matrix *adup)
/* Duplicate the matrix a.
 * Enter: matrix *a: original matrix.
 *        matrix *adup: location to place copy of a.           8/16/94-DWM */
{
  adup->w = a->w;  adup->h = a->h;
  memcpy(adup->m, a->m, sizeof(float)*a->w*a->h);
}

void matident(matrix *i, long size)
/* Create a size x size identity matrix.
 * Enter: matrix *i: pointer to matrix.
 *        long size: size of matrix to create.                  8/16/94-DWM */
{
  long j;

  i->w = i->h = size;
  for (j=1; j<size*size-1; j++)
    i->m[j] = 0;
  for (j=0; j<size; j++)
    i->m[j*(size+1)] = 1;
}

long matinv(matrix *a, matrix *ainv)
/* Invert a square matrix.  The input matrix is changed to the identity
 *  matrix.
 * Enter: matrix *a: point to matrix to invert.
 *        matrix *ainv: location to store inverse matrix.  This can not be
 *                        the calling matrix.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  if (a->w!=a->h)
    return(1);
  matident(ainv, a->w);
  return(rowreduce2(a, ainv));
}

long matmul(matrix *a, matrix *b, matrix *amulb)
/* Computes the matrix function a*b, where a is an n x m matrix and b is an
 *   p x n matrix.  The result is a p x m matrix.
 * Enter: matrix *a, *b: pointer to matrices to multiply.
 *        matrix *amulb: pointer to resultant matrix.  This can not be one of
 *                       the calling matrices.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, j, k, w, h, s=a->w;

  if (a->w!=b->h)
    return(1);
  w = amulb->w = b->w;
  h = amulb->h = a->h;
  for (j=0; j<h; j++)
    for (i=0; i<w; i++)
      for (k=amulb->m[j*w+i]=0; k<s; k++)
        amulb->m[j*w+i] += a->m[j*s+k]*b->m[k*w+i];
  return(0);
}

void matscale(matrix *a, float scalar, matrix *result)
/* Multiply a matrix by a scalar value.
 * Enter: matrix *a: pointer to matrix.
 *        float scalar: scalar multiplier.
 *        matrix *result: location to store result, may be a.  8/16/94-DWM */
{
  long i, len=a->w*a->h;

  for (i=0; i<len; i++)
    result->m[i] = scalar*a->m[i];
  result->w = a->w;
  result->h = a->h;
}

long matsub(matrix *a, matrix *b, matrix *aminusb)
/* Computes the matrix function a-b, where a and b are both n x m matrices.
 * Enter: matrix *a, *b: pointer to matrices to subtract.
 *        matrix *aminusb: pointer to resultant matrix.  This may be one of
 *                        the calling matrices.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, len=a->w*a->h;

  if (b->w*b->h!=len)
    return(1);
  for (i=0; i<len; i++)
    aminusb->m[i] = a->m[i]-b->m[i];
  aminusb->w = a->w;
  aminusb->h = a->h;
  return(0);
}

long matsym(matrix *a)
/* Make the specified matrix symmetrical.  This is done by copying the upper
 *  triangle to the lower triangle.  The matrix must be at least as wide as
 *  it is tall.
 * Enter: matrix *a: pointer to matrix to make symmetric.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, j;

  if (a->w<a->h)
    return(1);
  for (j=1; j<a->h; j++)
    for (i=0; i<j; i++)
      a->m[j*a->w+i] = a->m[i*a->w+j];
  return(0);
}

void mattoabg(float *r, float *abg)
/* Convert a 3x3 orthonormal matrix to a set of angles.  The angles are in
 *  the order alpha, beta, gammera
 * Enter: float *r: pointer to array of nine values.
 *        float *abg: pointer to array to store 3 angles.       1/8/98-DWM */
{
  abg[1] = acos(r[8]);
  abg[2] = atan2(r[6], -r[7]);
  abg[0] = atan2(r[2], r[5]);
  if (fabs(sin(abg[2])*cos(abg[0])+cos(abg[2])*cos(abg[1])*sin(abg[0])-r[1])>
      1e-6) {
    abg[1] = -acos(r[8]);
    abg[2] = atan2(-r[6], r[7]);
    abg[0] = atan2(-r[2], -r[5]); }
}

void mattoang(float *r, float *a)
/* Convert a 3x3 orthonormal matrix to a set of angles.  See the Derive file
 *  MTOSIN.MTH for the exact form of these angles.  The angles are in the
 *  order theta, phi, psi.
 * Enter: float *r: pointer to array of nine values.
 *        float *a: pointer to array to store 3 angles.        11/1/94-DWM */
{
  double temp;

  temp = -r[5];  if (fabs(temp)>1)  temp = (temp>0)-(temp<0);
  a[1] = asin(temp);
  if (fabs(cos(a[1]))>1e-30) {
    temp = r[3]/cos(a[1]);  if (fabs(temp)>1)  temp = (temp>0)-(temp<0);
    a[0] = asin(temp);
    if (cos(a[0])*cos(a[1])*r[4]<0)
      a[0] = PI-a[0];
    temp = r[2]/cos(a[1]);  if (fabs(temp)>1)  temp = (temp>0)-(temp<0);
    a[2] = asin(temp);
    if (cos(a[2])*cos(a[1])*r[8]<0)
      a[2] = PI-a[2];
    if (fabs(sin(a[1])*sin(a[2])*sin(a[0])+cos(a[2])*cos(a[0])-r[0])>1e-30 ||
        fabs(sin(a[1])*sin(a[2])*cos(a[0])-cos(a[2])*sin(a[0])-r[1])>1e-30) {
          a[0] = PI+a[0];  a[1] = PI-a[1];  a[2] = PI+a[2]; } }
  else {
    a[2] = 0;
    temp = r[0];  if (fabs(temp)>1)  temp = (temp>0)-(temp<0);
    a[0] = acos(temp);
    if (fabs(-sin(a[0])-r[1])>1e-30)
      a[0] *= -1;  }
}

void mattoypr(float *r, float *a)
/* Convert a 3x3 orthonormal matrix to a set of yaw, pitch, and roll angles.
 * Enter: float *r: pointer to array of nine values.
 *        float *a: pointer to array to store 3 angles.        9/10/95-DWM */
{
  if (r[5])  a[2] = atan2(r[2], r[5]);
  else       a[2] = PI*0.5*((r[2]>0)-(r[2]<0));
  if (r[7])  a[0] = atan2(r[6], -r[7]);
  else       a[0] = PI*0.5*((r[6]>0)-(r[7]<0));
  a[1] = asin(r[8]);
  if (cos(a[1])*sin(a[0])*r[6]<0)  a[1] = PI-a[1];
}

void mattrans(matrix *a, matrix *atrans)
/* Compute the transpose of a matrix.
 * Enter: matrix *a: pointer to matrix.
 *        matrix *atrans: location to store result, can not be a.
 *                                                             8/16/94|DWM */
{
  long i, j;

  atrans->w = a->h;
  atrans->h = a->w;
  for (j=0; j<a->w; j++)
    for (i=0; i<a->h; i++)
      atrans->m[j*atrans->w+i] = a->m[i*a->w+j];
}

void matzero(matrix *a)
/* Zero a matrix.
 * Enter: matrix *a: pointer to matrix to zero.                10/5/94-DWM */
{
  long i, s=a->w*a->h;

  for (i=0; i<s; i++)
    a->m[i] = 0;
}

float *qconj(float *a, float *astar)
/* Conjugate a quaternion.  The quarternions are 3 dimensional.  The format
 *  is (s, v1, v2, v3), where s is the scalar part of the quarternion, and v1
 *  through v3 is the vector part.
 * Enter: float *a: array holding initial quarternion.
 *        float *astar: array to place result.
 * Exit:  float *astar: same as input.                        10/24/94-DWM */
{
  astar[0] = a[0];  astar[1] = -a[1];  astar[2] = -a[2];  astar[3] = -a[3];
  return(astar);
}

float qdot(float *a, float *b)
/* Take the dot product of two quaternions.  The quarternions are 3
 *  dimensional.  The format is (s, v1, v2, v3), where s is the scalar part
 *  of the quarternion, and v1 through v3 is the vector part.
 * Enter: float *a, *b: array holding the two quarternions.
 * Exit:  float adotb: dot product.                           10/24/94-DWM */
{
  return(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]);
}

float *qmul(float *a, float *b, float *amulb)
/* Multiply two quaternions.  The quarternions are 3 dimensional.  The format
 *  is (s, v1, v2, v3), where s is the scalar part of the quarternion, and v1
 *  through v3 is the vector part.
 * Enter: float *a, *b: array holding the two multipicands.
 *        float *amulb: array to place result.  This may be one of the
 *                       calling quarternions.
 * Exit:  float *amulb: same as input.                        10/24/94-DWM */
{
  float ab[4];

  ab[0] = a[0]*b[0] - (a[1]*b[1]+a[2]*b[2]+a[3]*b[3]);
  ab[1] = a[0]*b[1] + b[0]*a[1] + a[2]*b[3]-a[3]*b[2];
  ab[2] = a[0]*b[2] + b[0]*a[2] + a[3]*b[1]-a[1]*b[3];
  ab[3] = a[0]*b[3] + b[0]*a[3] + a[1]*b[2]-a[2]*b[1];
  amulb[0] = ab[0];  amulb[1] = ab[1];  amulb[2] = ab[2];  amulb[3] = ab[3];
  return(amulb);
}

float qnorm(float *a)
/* Take the normal of a quaternion ºaº.  The quarternion is 3 dimensional.
 *  The format is (s, v1, v2, v3), where s is the scalar part of the
 *  quarternion, and v1 through v3 is the vector part.
 * Enter: float *a: quarternion to take the normal of.
 * Exit:  float anorm: norm.                                  10/24/94-DWM */
{
  return(sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]+a[3]*a[3]));
}

void qtomat(float *q, float *rot)
/* Convert a quaternion to a rotation matrix.  The quarternion is 3
 *  dimensional.  The format is (s, v1, v2, v3), where s is the scalar part
 *  of the quarternion, and v1 through v3 is the vector part.
 * Enter: float *q: quarternion to convert.
 *        float *rot: array of nine floats to store the resulting matrix.
 *                                                            10/18/97-DWM */
{
  rot[0] = q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  rot[1] = 2*(q[1]*q[2]-q[0]*q[3]);
  rot[2] = 2*(q[1]*q[3]+q[0]*q[2]);
  rot[3] = 2*(q[1]*q[2]+q[0]*q[3]);
  rot[4] = q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  rot[5] = 2*(q[2]*q[3]-q[0]*q[1]);
  rot[6] = 2*(q[1]*q[3]-q[0]*q[2]);
  rot[7] = 2*(q[2]*q[3]+q[0]*q[1]);
  rot[8] = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
}

float *qunit(float *a, float *aunit)
/* Make a unit quaternion.  The quarternion is 3 dimensional.  The format is
 *  (s, v1, v2, v3), where s is the scalar part of the quarternion, and v1
 *  through v3 is the vector part.
 * Enter: float *a: quarternion to unitize.
 *        float *aunit: array to place result.  This may be a.
 * Exit:  float *aunit: same as input.                        10/24/94-DWM */
{
  float norm=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]+a[3]*a[3]);

  if (!norm)
    aunit[0] = aunit[1] = aunit[2] = aunit[3] = 0;
  else {
    aunit[0] = a[0]/norm;  aunit[1] = a[1]/norm;
    aunit[2] = a[2]/norm;  aunit[3] = a[3]/norm; }
  return(aunit);
}

void rotate_fit(float *XYZ, float *xyz, long fitstart, long fitnum,
             long modstart, long modnum, long scale, float *rot, float *tran)
/* Rotate and translate the calculated geometry based on a least-squares fit
 *  to the known geometry.  The method used is taken from two papers: K. S.
 *  Arun, T. S. Huang, and S. D. Blostein, "Least-Squares Fitting of Two 3-D
 *  Point Sets", IEEE Transactions on Pattern Analysis and Machine
 *  Intelligence, Vol. 9, No. 5, Sept. 1987, pp. 698-700; and Shinji
 *  Umeyama, "Least-Squares Estimation of Transform Parameters Between Two
 *  Point Patterns", same magazine, Vol. 13, No. 4, April 1991, pp. 376-380.
 *  The singular value decomposition method (SVD) used to find u.d.v`=h was
 *  taken from Bert W. Rust and Walter R. Burrus, "Mathematical Programming
 *  and the Numerical Solution of Linear Equations", American Elsevier
 *  Publishing Company, Inc., New York 1972.
 * Enter: float *XYZ: array containing reference coordinates to fit to.
 *        float *xyz: array containing data to be fitted.  If the x
 *                    coordinate of either XYZ or xyz is greater than 1e20,
 *                    that data point is unused.
 *        long fitstart: starting node number for fit (usually 0).
 *        long fitnum: number of nodes to use in fit.
 *        long modstart: starting node number for modifying data.
 *        long modnum: number of nodes to fit and modify.
 *        long scale: 0 to only rotate, 1 to scale and rotate, +2 to only
 *                    compute rotation and translation values.
 *        float *rot: location to store rotation matrix (9 values), or 0 for
 *                    don't store.  Scale is also in this matrix
 *        float *tran: location to store traslation (3 values), or 0 for
 *                     don't store.                           10/12/94-DWM */
{
  double x[3]={0,0,0}, X[3]={0,0,0}, y[3], Y[3], l[3], err, best=1e300;
  long i, j, total=0, besti;
  double cc[132], t[3], sx2=0, c, tr;
  lmatrix *h=lmat(cc), *ht=lmat(cc+11), *hth=lmat(cc+22), *hht=lmat(cc+33);
  lmatrix *r=lmat(cc+44), *ut=lmat(cc+55), *v=lmat(cc+66), *s=lmat(cc+77);
  lmatrix *t1=lmat(cc+88), *t2=lmat(cc+99), *d=lmat(cc+110), *u=lmat(cc+121);
  int mul[]={1,1,1,1,1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1};

  for (i=fitstart; i<fitnum; i++)
    if (xyz[i*3]<1e20 && XYZ[i*3]<1e20) {
      x[0] += xyz[i*3];  x[1] += xyz[i*3+1];  x[2] += xyz[i*3+2];
      X[0] += XYZ[i*3];  X[1] += XYZ[i*3+1];  X[2] += XYZ[i*3+2];
      total++; }
  x[0] /= total;  x[1] /= total;  x[2] /= total;
  X[0] /= total;  X[1] /= total;  X[2] /= total;
  for (i=fitstart; i<fitnum; i++)
    if (xyz[i*3]<1e20 && XYZ[i*3]<1e20)
      for (j=0; j<3; j++)
        sx2 += (xyz[i*3+j]-x[j])*(xyz[i*3+j]-x[j]);
  lmatzero(h);
  for (i=fitstart; i<fitnum; i++)
    if (xyz[i*3]<1e20 && XYZ[i*3]<1e20) {
      y[0] = xyz[i*3]  -x[0];  Y[0] = XYZ[i*3]  -X[0];
      y[1] = xyz[i*3+1]-x[1];  Y[1] = XYZ[i*3+1]-X[1];
      y[2] = xyz[i*3+2]-x[2];  Y[2] = XYZ[i*3+2]-X[2];
      h->m[0]+=y[0]*Y[0];  h->m[1]+=y[0]*Y[1];  h->m[2]+=y[0]*Y[2];
      h->m[3]+=y[1]*Y[0];  h->m[4]+=y[1]*Y[1];  h->m[5]+=y[1]*Y[2];
      h->m[6]+=y[2]*Y[0];  h->m[7]+=y[2]*Y[1];  h->m[8]+=y[2]*Y[2]; }
  lmattrans(h, ht);
  lmatmul(ht, h, hth);
  leigen(hth->m, l);
  leigenvect(hth->m, l, v->m, 1);                   /* Column eigenvectors */
  lmatmul(h, ht, hht);
  leigen(hht->m, l);
  leigenvect(hht->m, l, ut->m, 0);                     /* Row eigenvectors */
  lmattrans(ut, u);
  lmatzero(d);  lmatzero(s);
  d->m[0] = sqrt(l[0]);  d->m[4] = sqrt(l[1]);  d->m[8] = sqrt(l[2]);
  for (i=0; i<8; i++) {
    s->m[0]=mul[i*3];  s->m[4]=mul[i*3+1];  s->m[8]=mul[i*3+2];
    lmattrans(v, t1);             /* t1 = v`         */
    lmatmul(s, t1, t2);           /* t2 = s.v`       */
    lmatmul(d, t2, t1);           /* t1 = d.s.v`     */
    lmatmul(u, t1, t2);           /* t2 = u.d.s.v`   */
    lmatadd(1, t2, -1, h, t2);    /* t2 = u.d.s.v`-h */
    for (j=0, err=0; j<9; j++)
      err += t2->m[j]*t2->m[j];
    if (err<best) {
      best = err;
      besti = i; } }
  s->m[0]=mul[besti*3];  s->m[4]=mul[besti*3+1];  s->m[8]=mul[besti*3+2];
  lmatmul(s, ut, t1);
  lmatmul(v, t1, r);
  tr = d->m[0]+d->m[4]+d->m[8];
  if (lmatdet(r)<0) {
    tr -= d->m[8]*2;
    s->m[8] *= -1;
    lmatmul(s, ut, t1);
    lmatmul(v, t1, r); }
  c = tr/sx2;
  if (scale&1)
    for (i=0; i<9; i++)
      r->m[i] *= c;
  t[0] = X[0]-x[0]*r->m[0]-x[1]*r->m[1]-x[2]*r->m[2];
  t[1] = X[1]-x[0]*r->m[3]-x[1]*r->m[4]-x[2]*r->m[5];
  t[2] = X[2]-x[0]*r->m[6]-x[1]*r->m[7]-x[2]*r->m[8];
  if (rot)   memcpy(rot, r->m, 9*sizeof(float));
  if (tran)  memcpy(tran, t, 3*sizeof(float));
  if (scale>=2)  return;
  for (i=modstart; i<modnum; i++)
    if (xyz[i*3]<1e20 && XYZ[i*3]<1e20) {
      y[0] = xyz[i*3];  y[1] = xyz[i*3+1];  y[2] = xyz[i*3+2];
      xyz[i*3]   = y[0]*r->m[0]+y[1]*r->m[1]+y[2]*r->m[2]+t[0];
      xyz[i*3+1] = y[0]*r->m[3]+y[1]*r->m[4]+y[2]*r->m[5]+t[1];
      xyz[i*3+2] = y[0]*r->m[6]+y[1]*r->m[7]+y[2]*r->m[8]+t[2]; }
}

long rowreduce(matrix *a)
/* Row reduce a matrix.  The matrix must be at least as wide as it is high.
 *  The main diagonal is changed to 1, and all other elements in the main
 *  square are changed to zero.  If the process fails due to an undetermined
 *  value, the matrix is only partially row reduced.
 * Enter: matrix *a: matrix to row reduce.  The matrix is changed.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, j, k, w=a->w, h=a->h;
  double temp;

  if (w<h)
    return(1);
  for (i=0; i<h; i++) {
    for (j=i; j<h; j++)
      if (a->m[j*w+i])
        break;
    if (j==h)
      return(1);
    if (j!=i)
      for (k=i; k<w; k++) {
        temp = a->m[i*w+k];
        a->m[i*w+k] = a->m[j*w+k];
        a->m[j*w+k] = temp;  }
    for (j=w-1; j>=i; j--)
      a->m[i*w+j] /= a->m[i*w+i];
    for (j=0; j<h; j++)
      if (j!=i) if (temp=a->m[j*w+i])
        for (k=i; k<w; k++)
          a->m[j*w+k] -= a->m[i*w+k]*temp; }
  return(0);
}

long rowreduce2(matrix *a, matrix *b)
/* Row reduce matrix a, modifying matrix b in the same manner.  Matrix a must
 *  be at least as wide as it is high, and matrix b must be the same height
 *  as matrix a.  The main diagonal of a is changed to 1, and all other
 *  elements in the main square are changed to zero.  If the process fails
 *  due to an undetermined value, the matrices are only partially row
 *  reduced.
 * Enter: matrix *a: matrix to row reduce.  The matrix is changed.
 *        matrix *b: matrix to modify in the same manner as matrix a.
 * Exit:  long failed: 0 if successful, 1 if failed.           8/16/94-DWM */
{
  long i, j, k, w=a->w, h=a->h, w2=b->w, iw, iw2;
  double temp;

  if (w<h || b->h!=h)
    return(1);
  for (i=iw=iw2=0; i<h; i++, iw+=w, iw2+=w2) {
    for (j=i; j<h; j++)
      if (a->m[j*w+i])
        break;
    if (j==h)
      return(1);
    if (j!=i) {
      for (k=i; k<w; k++) {
        temp = a->m[iw+k];
        a->m[iw+k] = a->m[j*w+k];
        a->m[j*w+k] = temp;  }
      for (k=0; k<w2; k++) {
        temp = b->m[iw2+k];
        b->m[iw2+k] = b->m[j*w2+k];
        b->m[j*w2+k] = temp;  } }
    temp = 1/a->m[iw+i];
    a->m[iw+i] = 1;
    for (j=i+1; j<w; j++)
      a->m[iw+j] *= temp;
    for (j=0; j<w2; j++)
      b->m[iw2+j] *= temp;
    for (j=0; j<h; j++)
      if (j!=i) if (temp=a->m[j*w+i]) {
        a->m[j*w+i] = 0;
        for (k=i+1; k<w; k++)
          a->m[j*w+k] -= a->m[iw+k]*temp;
        for (k=0; k<w2; k++)
          b->m[j*w2+k] -= b->m[iw2+k]*temp; } }
  return(0);
}

float *unit(float *v, float *r)
/* Compute the unit vector of the given vector.  If the vector is the zero
 *  vector, the zero vector is returned.
 * Enter: float *v: vector to compute unit vector from.
 *        float *r: location to place resultant vector.
 * Exit:  float *r: location of resultant vector (same as entered value).
 *                                                             8/20/91|DWM */
{
  float len=0;
  short i;

  for (i=0; i<3; i++)
    len += v[i] * v[i];
  len = sqrt(len);
  if (len)
    for (i=0; i<3; i++)
      r[i] = v[i] / len;
  else
    r[0] = r[1] = r[2] = 0;
  return(r);
}

void yprtomat(float *rot, float y, float p, float r)
/* Convert a set of roll, pitch, and yaw angles to a 3x3 orthonormal rotation
 *  matrix.  This is actually the matrix
 *   [  cos(r) sin(r) 0 ] [ 1   0       0    ] [ cos(y)  sin(y) 0 ]
 *   [ -sin(r) cos(r) 0 ] [ 0 cos(p) -sin(p) ] [   0       0    1 ]
 *   [    0      0    1 ] [ 0 sin(p)  cos(p) ] [ sin(y) -cos(y) 0 ]
 * Enter: float *rot: array of nine values for matrix.
 *        float r: roll about the x axis.
 *        float p: pitch about the y axis.
 *        float y: yaw about the z axis.                       9/10/95-DWM */
{
  rot[0] = cos(r)*cos(y)-sin(p)*sin(r)*sin(y);
  rot[1] = sin(p)*sin(r)*cos(y)+cos(r)*sin(y);  rot[2] = cos(p)*sin(r);
  rot[3] = -sin(r)*cos(y)-sin(p)*cos(r)*sin(y);
  rot[4] = sin(p)*cos(r)*cos(y)-sin(r)*sin(y);  rot[5] = cos(p)*cos(r);
  rot[6] = cos(p)*sin(y);  rot[7] = -cos(p)*cos(y);  rot[8] = sin(p);
}

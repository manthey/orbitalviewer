#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "dlt.h"

#define sq(x)  ((x)*(x))

long ldlt_to_physical(double *dlt, double *phys, lmatrix *ang, long ten)
/* Convert a set of 11 dlt parameters to 10 physical parameters.  Distortion
 *  is ignored.  This routine can fail since there are singularities in the
 *  equations.
 * Enter: double *dlt: array with 11 dlt parameters, L1 to L11.
 *        double *phys: array to store physical parameters in the order
 *                      X0, Y0, Z0, Cx, Cy, x0, y0, theta, phi, psi.  Null
 *                      for only matrix output.
 *        lmatrix *ang: 3x3 matrix to store rotation in.  Pass a zero pointer
 *                      if no output is desired.
 *        long ten: 0 to use 11 independent parameters.  1 to use 10.
 * Exit:  long okay: 0 for failed, 1 for okay.                 4/13/95-DWM */
{
  double L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11;
  double X0, Y0, Z0, Cx, Cy, x0, y0, ta, pa, sa;
  double m11, m12, m13, m21, m22, m23, m31, m32, m33, den;
  long ch;

  if (ang) {  ang->w = ang->h = 3; }
  L1=dlt[0]; L2=dlt[1]; L3=dlt[2]; L4=dlt[3]; L5=dlt[4]; L6=dlt[5];
  L7=dlt[6]; L8=dlt[7]; L9=dlt[8]; L10=dlt[9]; L11=dlt[10];
  if (ten)
	L1 = -(L10*L10*L3*L7-L10*(L11*(L2*L7+L3*L6)+L2*L5*L9)+L11*L11*L2*L6-
		   L11*L3*L5*L9+L9*L9*(L2*L6+L3*L7))/
		  (L10*L10*L5-L10*L6*L9+L11*(L11*L5-L7*L9));
  ta = atan2((L10*L3-L11*L2),(L1*L11-L3*L9));
  pa = atan2((L10*cos(ta)+L9*sin(ta)),L11);
  sa = -atan2((L9*cos(pa)-L11*sin(pa)*sin(ta)),(L11*cos(ta)));
  m11= sin(pa)*sin(sa)*sin(ta)+cos(sa)*cos(ta);
  m12= sin(pa)*sin(sa)*cos(ta)-cos(sa)*sin(ta);
  m13= cos(pa)*sin(sa);
  m21= cos(pa)*sin(ta);
  m22= cos(pa)*cos(ta);
  m23= -sin(pa);
  m31= sin(pa)*cos(sa)*sin(ta)-sin(sa)*cos(ta);
  m32= sin(pa)*cos(sa)*cos(ta)+sin(sa)*sin(ta);
  m33= cos(pa)*cos(sa);
  if (ang) {
	ang->m[0] = m11;  ang->m[1] = m12;  ang->m[2] = m13;
	ang->m[3] = m21;  ang->m[4] = m22;  ang->m[5] = m23;
	ang->m[6] = m31;  ang->m[7] = m32;  ang->m[8] = m33; }
  if (!L9 || !(m22*m31-m21*m32) || !(m12*m31-m11*m32) || !m21||!m11||!m31)
	return(0);
  x0 = (L1*L9+L2*L10+L3*L11)/(L9*L9+L10*L10+L11*L11);
  y0 = (L5*L9+L6*L10+L7*L11)/(L9*L9+L10*L10+L11*L11);
  ch = 0;  if (L9*m11<L10*m12)  ch = 1;
  if (L9*m11<L11*m13 && L10*m12<L11*m13 )  ch = 2;
  switch (ch) {
	case 0: Cx = m31*(L9 *x0-L1)/(L9 *m11); break;
	case 1: Cx = m32*(L10*x0-L2)/(L10*m12); break;
	case 2: Cx = m33*(L11*x0-L3)/(L11*m13); }
  ch = 0;  if (L9*m21<L10*m22)  ch = 1;
  if (L9*m21<L11*m23 && L10*m22<L11*m23 )  ch = 2;
  switch (ch) {
	case 0: Cy = m31*(L9 *y0-L5)/(L9 *m21); break;
	case 1: Cy = m32*(L10*y0-L6)/(L10*m22); break;
	case 2: Cy = m33*(L11*y0-L7)/(L11*m23); }
  if (!(den=(Cx*Cy*L9*(m11*(m22*m33-m23*m32)+m12*(m23*m31-m21*m33)+
	  m13*(m21*m32-m22*m31)))))
	return(0);
  Z0 = -m31*(Cx*(Cy*(m11*m22-m12*m21)+(L8-y0)*(m11*m32-m12*m31))+Cy*(L4-x0)*
	   (m22*m31-m21*m32))/den;
  Y0 = (Cy*(L9*Z0*(m21*m33-m23*m31)+m21*m31)+m31*m31*(L8-y0))/
	   (Cy*L9*(m22*m31-m21*m32));
  X0 = -(L9*(Y0*m32+Z0*m33)+m31)/(L9*m31);
  if (fabs(X0*m11+Y0*m12+Z0*m13)<1e-20)  Cx = fabs(Cy)*(1-2*(Cx<0));
  if (fabs(X0*m21+Y0*m22+Z0*m23)<1e-20)  Cy = fabs(Cx)*(1-2*(Cy<0));
  if (Cx<0 || Cy<0) {
	if (Cx<0) {
	  Cx *= -1;  sa *= -1;  pa += PI;
	  if (pa>2*PI)  pa -= 2*PI; }
	if (Cy<0) {
	  Cy *= -1;  sa += PI;
	  if (sa>2*PI)  sa -= 2*PI; }
	if (ang)  langtomat(ang->m, ta, pa, sa); }
  if (!phys)
	return(1);
  phys[0]=X0; phys[1]=Y0; phys[2]=Z0; phys[3]=Cx; phys[4]=Cy;
  phys[5]=x0; phys[6]=y0; phys[7]=ta; phys[8]=pa; phys[9]=sa;
  return(1);
}

long lphysical_to_dlt(double *phys, double *dlt, long focal)
/* Convert a set of physical parameters to DLT parameters.
 * Enter: double *phys: The physical parameters, either in the order
 *                      X0, Y0, Z0, Cx, Cy, x0, y0, theta, phi, psi  or
 *                      X0, Y0, Z0, lx, ly, x0, y0, theta, phi, psi, f.
 *        double *dlt: location to store 18 dlt/distortion parameters.
 *                     Distortion is set to zero.
 *        long focal: 0 for Cx, Cy; 1 for lx, ly, and f.
 * Exit:  long okay: 0 for failed, 1 for okay.                 2/12/96-DWM */
{
  double P, x0, y0, Cx, Cy, X0, Y0, Z0, ta, pa, sa;
  double m11, m12, m13, m21, m22, m23, m31, m32, m33;
  X0 = phys[0];  Y0 = phys[1];  Z0 = phys[2];
  if (!focal || !phys[10]) { Cx = phys[3];           Cy = phys[4]; }
  else                     { Cx = phys[3]/phys[10];  Cy = phys[4]/phys[10]; }
  x0 = phys[5];  y0 = phys[6];
  ta = phys[7];  pa = phys[8];  sa = phys[9];
  while (ta<0)  ta += 2*PI;  while (ta>=2*PI)  ta -= 2*PI;
  while (pa<0)  pa += 2*PI;  while (pa>=2*PI)  pa -= 2*PI;
  while (sa<0)  sa += 2*PI;  while (sa>=2*PI)  sa -= 2*PI;
  if (fabs(ta)<1e-7)  ta = 1e-7;
  if (fabs(pa)<1e-7)  pa = 1e-7;
  if (fabs(sa)<1e-7)  sa = 1e-7;
  if (fabs(ta-PI/2)<1e-7)  ta = PI/2+1e-7;
  if (fabs(pa-PI/2)<1e-7)  pa = PI/2+1e-7;
  if (fabs(sa-PI/2)<1e-7)  sa = PI/2+1e-7;
  if (fabs(ta-PI)<1e-7)  ta = PI+1e-7;
  if (fabs(pa-PI)<1e-7)  pa = PI+1e-7;
  if (fabs(sa-PI)<1e-7)  sa = PI+1e-7;
  if (fabs(ta-3*PI/2)<1e-7)  ta = 3*PI/2+1e-7;
  if (fabs(pa-3*PI/2)<1e-7)  pa = 3*PI/2+1e-7;
  if (fabs(sa-3*PI/2)<1e-7)  sa = 3*PI/2+1e-7;
  m11= sin(pa)*sin(sa)*sin(ta)+cos(sa)*cos(ta);
  m12= sin(pa)*sin(sa)*cos(ta)-cos(sa)*sin(ta);
  m13= cos(pa)*sin(sa);
  m21= cos(pa)*sin(ta);
  m22= cos(pa)*cos(ta);
  m23= -sin(pa);
  m31= sin(pa)*cos(sa)*sin(ta)-sin(sa)*cos(ta);
  m32= sin(pa)*cos(sa)*cos(ta)+sin(sa)*sin(ta);
  m33= cos(pa)*cos(sa);
  P = -m31*X0-m32*Y0-m33*Z0;
  if (fabs(P)<sqrt(fabs(X0*X0+Y0*Y0+Z0*Z0))*1e-4) {
	Z0 = fabs(X0*X0+Y0*Y0+Z0*Z0);
	X0 = Z0*m31;  Y0 = Z0*m32;  Z0 *= m33;
	P = -m31*X0-m32*Y0-m33*Z0; }
  if (!P)  return(0);
  dlt[0] = (x0*m31-Cx*m11)/P;
  dlt[1] = (x0*m32-Cx*m12)/P;
  dlt[2] = (x0*m33-Cx*m13)/P;
  dlt[3] = ((Cx*m11-x0*m31)*X0+(Cx*m12-x0*m32)*Y0+(Cx*m13-x0*m33)*Z0)/P;
  dlt[4] = (y0*m31-Cy*m21)/P;
  dlt[5] = (y0*m32-Cy*m22)/P;
  dlt[6] = (y0*m33-Cy*m23)/P;
  dlt[7] = ((Cy*m21-y0*m31)*X0+(Cy*m22-y0*m32)*Y0+(Cy*m23-y0*m33)*Z0)/P;
  dlt[8] = m31/P;  dlt[9] = m32/P;  dlt[10] = m33/P;
  dlt[11] = x0;  dlt[12] = y0;
  dlt[13] = dlt[14] = dlt[15] = dlt[16] = dlt[17] = 0;
  return(1);
}

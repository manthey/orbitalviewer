#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "dlt.h"

#define sq(x)  ((x)*(x))

#ifdef used
long dlt_to_physical(float *dlt, float *phys, lmatrix *ang, long ten)
/* Convert a set of 11 dlt parameters to 10 physical parameters.  Distortion
 *  is ignored.  This routine can fail since there are singularities in the
 *  equations.
 * Enter: float *dlt: array with 11 dlt parameters, L1 to L11.
 *        float *phys: array to store physical parameters in the order
 *                     X0, Y0, Z0, Cx, Cy, x0, y0, theta, phi, psi.  Null for
 *                     only matrix output.
 *        lmatrix *ang: 3x3 matrix to store rotation in.  Pass a zero pointer
 *                      if no output is desired.
 *        long ten: 0 to use 11 independent parameters.  1 to use 10.
 * Exit:  long okay: 0 for failed, 1 for okay.                 4/13/95-DWM */
{
  double L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11;
  double X0, Y0, Z0, Cx, Cy, x0, y0, ta, pa, sa;
  double m11, m12, m13, m21, m22, m23, m31, m32, m33, den;

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
  Cy = m31*(L9*y0-L5)/(L9*m21);
  Cx = m31*(L9*x0-L1)/(L9*m11);
  if (!(den=(Cx*Cy*L9*(m11*(m22*m33-m23*m32)+m12*(m23*m31-m21*m33)+
	  m13*(m21*m32-m22*m31)))))
	return(0);
  Z0 = -m31*(Cx*(Cy*(m11*m22-m12*m21)+(L8-y0)*(m11*m32-m12*m31))+Cy*(L4-x0)*
	   (m22*m31-m21*m32))/den;
  Y0 = (Cy*(L9*Z0*(m21*m33-m23*m31)+m21*m31)+m31*m31*(L8-y0))/
	   (Cy*L9*(m22*m31-m21*m32));
  X0 = -(L9*(Y0*m32+Z0*m33)+m31)/(L9*m31);
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

long epipolar_line(float *e12, lmatrix *P2P1i, float *xy, float *abc)
/* Calculate the epipolar line in image 2 based on an already located point
 *  in image 1.
 * Enter: float *e12: epipolar point of camera 1 in image 2.  Array of 2
 *                    floats.
 *        lmatrix *P2P1i: P2 P1^-1.  This is a 3x3 matrix as calculated by
 *                        the function epipole().
 *        float *xy: point located in image 1.  Array of 2 floats.
 *        float *abc: location to store resulting line.  This is an array of
 *                    3 floats, in terms of ax+by+c=0.
 * Exit:  long okay: 0 for failed, 1 for okay.                 2/13/96-DWM */
{
  double mc[5], pc[5], p2[2];
  lmatrix *m=lmat(mc), *p=lmat(pc);

  m->w = 1;  m->h = 3;
  m->m[0] = xy[0];  m->m[1] = xy[1];  m->m[2] = 1;
  if (lmatmul(P2P1i, m, p))  return(0);
  if (!p->m[2])              return(0);
  p2[0] = p->m[0]/p->m[2];  p2[1] = p->m[1]/p->m[2];
  abc[0] = p2[1] - e12[1];
  abc[1] = e12[0] - p2[0];
  abc[2] = (p2[0]-e12[0])*e12[1] + (e12[1]-p2[1])*e12[0];
  return(1);
}

long epipolar_segment(float *dlt1, float *dlt2, float *e12, lmatrix *P2P1i,
                      float *xy, float *bounds, float *xyxy)
/* Calculate the epipolar line in image 2 based on an already located point
 *  in image 1.  Only a segment of this line is returned, clipped from the
 *  entire line based on the destination image and some maximum distance from
 *  each camera.  The bounds array contains the screen coordinates (left,
 *  top, right, bottom) and the distance limits (photo 1 maximum, photo 2
 *  maximum).
 * Enter: float *dlt1: 18 dlt parameters/distortion for camera 1.
 *        float *dlt2: 18 dlt parameters/distortion for camera 2.
 *        float *e12: epipolar point of camera 1 in image 2.  Array of 2
 *                    floats.  If null, epipolar_segment will calculate the
 *                    epipole on its own.
 *        lmatrix *P2P1i: P2 P1^-1.  This is a 3x3 matrix as calculated by
 *                        the function epipole().  If null, epipolar_segment
 *                        will calculate the epipole on its own.
 *        float *xy: point located in image 1.  Array of 2 floats.
 *        float *bounds: boundary conditions for the line cropping.  This is
 *                       an array of 6 values: x1,y1,x2,y2,max1,max2.  If
 *                       the maximum values are <=0, no distance cropping is
 *                       performed.
 *        float *xyxy: location to store resulting line segment.  This is an
 *                     array of 4 floats, in terms of x1,y1,x2,y2.  x1 will
 *                     always less than or equal to x2.  If equal to, y1 will
 *                     be less than y2.
 * Exit:  long okay: 0 for failed, 1 for okay.                 2/13/96-DWM */
{
  double P2P1c[11], mc[5], pc[5], dist, dx[3], r, oc[11];
  float abc[3], ae12[2], xyc[8], xyz[6], dltc[36], phys[22], p2[2];
  lmatrix *aP2P1i=lmat(P2P1c), *m=lmat(mc), *p=lmat(pc), *o=lmat(oc);
  long i, j, k, side=1, use;
  short np[]={0,0,0,1,1,0,1,1};

  if (e12 && P2P1i) {
    memcpy(ae12, e12, 2*sizeof(float));
    lmatdup(P2P1i, aP2P1i); }
  else {
    if (!epipole(dlt1, dlt2, ae12, aP2P1i))  return(0); }
  m->w = 1;  m->h = 3;
  m->m[0] = xy[0];  m->m[1] = xy[1];  m->m[2] = 1;
  if (lmatmul(aP2P1i, m, p))  return(0);
  if (!p->m[2])              return(0);
  p2[0] = p->m[0]/p->m[2];  p2[1] = p->m[1]/p->m[2];
  abc[0] = p2[1] - ae12[1];
  abc[1] = ae12[0] - p2[0];
  abc[2] = (p2[0]-ae12[0])*ae12[1] + (ae12[1]-p2[1])*ae12[0];
  if (!abc[0] && !abc[1])  return(0);
  dlt_to_physical(dlt1, phys, 0, 0);    /* Select other side of proj. inf. */
  dlt_to_physical(dlt2, phys+11, o, 0);
  if (o->m[2]*(phys[0]-phys[11])+o->m[5]*(phys[1]-phys[12])+
	  o->m[8]*(phys[2]-phys[13])>0) {
    ae12[0] = p2[0] + (p2[0]-ae12[0])*1000;
    ae12[1] = p2[1] + (p2[1]-ae12[1])*1000; }                  /* End swap */
  if (fabs(p2[0]-ae12[0])>fabs(p2[1]-ae12[1])) { /* Back off a pixel since */
    r = 1-2*(p2[0]>ae12[0]);                     /*  inf. is far away.     */
    p2[0] += r;  p2[1] += r*(p2[1]-ae12[1])/(p2[0]-ae12[0]); }
  else {
    r = 1-2*(p2[1]>ae12[1]);
    p2[1] += r;  p2[0] += r*(p2[0]-ae12[0])/(p2[1]-ae12[1]); }
  xyxy[0] = ae12[0];  xyxy[1] = ae12[1];
  xyxy[2] = p2[0];    xyxy[3] = p2[1];
  if (xyxy[0]>xyxy[2] || (xyxy[0]==xyxy[2] && xyxy[1]>xyxy[3])) {
    xyxy[0] = p2[0];    xyxy[1] = p2[1];
    xyxy[2] = ae12[0];  xyxy[3] = ae12[1];
    side = 0; }
  if (xyxy[0]>bounds[2] || xyxy[2]<bounds[0] || (xyxy[0]<bounds[1] &&
      xyxy[2]<bounds[1]) || (xyxy[0]>bounds[3] && xyxy[2]>bounds[3]))
    return(0);
  if (abc[1]) {
    if (xyxy[0]<bounds[0]) {
      xyxy[0] = bounds[0];  xyxy[1] = (-abc[2]-abc[0]*bounds[0])/abc[1]; }
    if (xyxy[2]>bounds[2]) {
      xyxy[2] = bounds[2];  xyxy[3] = (-abc[2]-abc[0]*bounds[2])/abc[1]; }
    if ((xyxy[1]<bounds[1] && xyxy[3]<bounds[1]) ||
        (xyxy[1]>bounds[3] && xyxy[3]>bounds[3]))  return(0);
    if (xyxy[1]<bounds[1]) {
      xyxy[1] = bounds[1];  xyxy[0] = (-abc[2]-abc[1]*bounds[1])/abc[0]; }
    if (xyxy[3]<bounds[1]) {
      xyxy[3] = bounds[1];  xyxy[2] = (-abc[2]-abc[1]*bounds[1])/abc[0]; }
    if (xyxy[1]>bounds[3]) {
      xyxy[1] = bounds[3];  xyxy[0] = (-abc[2]-abc[1]*bounds[3])/abc[0]; }
    if (xyxy[3]>bounds[3]) {
      xyxy[3] = bounds[3];  xyxy[2] = (-abc[2]-abc[1]*bounds[3])/abc[0]; } }
  else {
    if (xyxy[1]<bounds[1]) {
      xyxy[1] = bounds[1];  xyxy[0] = (-abc[2]-abc[1]*bounds[1])/abc[0]; }
    if (xyxy[3]>bounds[3]); {
      xyxy[3] = bounds[3];  xyxy[2] = (-abc[2]-abc[1]*bounds[3])/abc[0]; }
    if ((xyxy[0]<bounds[0] && xyxy[2]<bounds[0]) ||
        (xyxy[0]>bounds[2] && xyxy[2]>bounds[2]))  return(0); }
  for (k=0; k<2; k++)  if (bounds[k+4]) {
    memcpy(xyc,   xy, 2*sizeof(float));
    memcpy(xyc+2, xyxy, 2*sizeof(float));
    memcpy(xyc+4, xy, 2*sizeof(float));
    memcpy(xyc+6, xyxy+2, 2*sizeof(float));
    memcpy(dltc,    dlt1, 18*sizeof(float));
    memcpy(dltc+18, dlt2, 18*sizeof(float));
    image_to_geo(4, 2, 2, np, xyc, dltc, 0, 0, 2, xyz);
    dlt_to_physical(dltc+k*18, phys, 0, 0);
    dx[0] = sqrt(sq(phys[0]-xyz[0])+sq(phys[1]-xyz[1])+sq(phys[2]-xyz[2]));
    dx[1] = sqrt(sq(phys[0]-(xyz[0]+xyz[3])/2)+sq(phys[1]-(xyz[1]+xyz[4])/2)+
                 sq(phys[2]-(xyz[2]+xyz[5])/2));
    dx[2] = sqrt(sq(phys[0]-xyz[3])+sq(phys[1]-xyz[4])+sq(phys[2]-xyz[5]));
    if (k) {
      abc[0] = 2*dx[0]-4*dx[1]+2*dx[2];
      abc[1] = -3*dx[0]+4*dx[1]-dx[2];
      abc[2] = dx[0];
      dx[0] = 4*abc[0]*bounds[k+4]-4*abc[0]*abc[2]+abc[1]*abc[1];
      if (dx[0]<0)  return(0);
      dx[0] = sqrt(dx[0]);
      dx[1] = (dx[0]-abc[1])/(2*abc[0]);
      dx[0] = (-dx[0]-abc[1])/(2*abc[0]);
      if ((dx[0]<0 && dx[1]<0) || (dx[0]>1 && dx[1]>1))  return(0); }
    else {
      abc[0] = dx[2]-dx[0];  abc[1] = dx[0];
      if (!abc[0])  continue;
      dx[0] = 0;  dx[1] = 1;
      dx[abc[0]>0] = (bounds[k+4]-abc[1])/abc[0];
      if (dx[0]>1 || dx[1]<0)  return(0); }
    for (j=0; j<2; j++)
      if (dx[j]>0 && dx[j]<1) {
        for (i=0; i<3; i++)
          abc[i] = xyz[i]*(1-dx[j])+xyz[3+i]*dx[j];
        geo_to_image(1, abc, dlt2, p2);
        if (p2[0]<1e20)
          memcpy(xyxy+j*2, p2, 2*sizeof(float)); } }
  return(1);
}

long epipole(float *dlt1, float *dlt2, float *e12, lmatrix *P2P1i)
/* Calculate the epipole of camera 1 in image 2.  This is used along with the
 *  function epipolar_line when a point is known in image 1 and is being
 *  searched for in image 2.
 * Enter: float *dlt1: 11 dlt parameters for camera 1.
 *        float *dlt2: 11 dlt parameters for camera 2.
 *        float *e12: epipole of camera 1 in image 2.  Array of 2 floats.
 *        lmatrix *P2P1i: P2 P1^-1, used for epipolar line calculation.  This
 *                        is a 3x3 matrix.
 * Exit:  long okay: 0 for failed, 1 for okay.                 2/13/96-DWM */
{
  double P1c[11], P2c[11], P1ic[11], p1c[5], p2c[5], tempc[5];
  lmatrix *P1=lmat(P1c), *P2=lmat(P2c), *P1i=lmat(P1ic);
  lmatrix *p1=lmat(p1c), *p2=lmat(p2c), *temp=lmat(tempc);

  P1->w = P1->h = P2->w = P2->h = 3;
  P1->m[0] = dlt1[0];  P1->m[1] = dlt1[1];  P1->m[2] = dlt1[2];
  P1->m[3] = dlt1[4];  P1->m[4] = dlt1[5];  P1->m[5] = dlt1[6];
  P1->m[6] = dlt1[8];  P1->m[7] = dlt1[9];  P1->m[8] = dlt1[10];
  P2->m[0] = dlt2[0];  P2->m[1] = dlt2[1];  P2->m[2] = dlt2[2];
  P2->m[3] = dlt2[4];  P2->m[4] = dlt2[5];  P2->m[5] = dlt2[6];
  P2->m[6] = dlt2[8];  P2->m[7] = dlt2[9];  P2->m[8] = dlt2[10];
  p1->w = p2->w = 1;  p1->h = p2->h = 3;
  p1->m[0] = dlt1[3];  p1->m[1] = dlt1[7];  p1->m[2] = 1;
  p2->m[0] = dlt2[3];  p2->m[1] = dlt2[7];  p2->m[2] = 1;
  if (lmatinv(P1, P1i))               return(0);
  if (lmatmul(P2, P1i, P2P1i))        return(0);
  if (lmatmul(P2P1i, p1, temp))       return(0);
  if (lmatadd(-1, temp, 1, p2, temp)) return(0);
  if (!temp->m[2])                    return(0);
  e12[0] = temp->m[0]/temp->m[2];  e12[1] = temp->m[1]/temp->m[2];
  return(1);
}

long foe(long n, float *xy1, float *xy2, long err, float *foexy)
/* Calculate the focus of expansion based on a set of corresponding points
 *  between two images.  The located focus of expansion is returned.  If
 *  requested, the maximal or statistical error is given.
 * Enter: long n: number of points in each image to use.  Must be at least 2.
 *        float *xy1: points in first image in order x1,y1,x2,y2,x3,...
 *        float *xy2: points in second image in order x1,y1,x2,y2,x3,...
 *                    These must have a direct one-to-one correspondence with
 *                    the points from image 1.
 *        long err: 0 for error value equal to zero, 1 for maximal error, 2
 *                  for statistical (3 sigma) error.
 *        float *foexy: location to store resultant foe in image coordinates.
 *                      This is an array of 3 floats: x,y,error.  If error is
 *                      not calculated, the third value is set to zero.
 * Exit:  long okay: 0 for failed, 1 for okay.                 2/21/96-DWM */
{
  double mc[8], l[3], A, B, C, e, se2=0, s=0;
  lmatrix *m=lmat(mc);
  long i;

  if (n<2)  return(0);
  m->w = 3;  m->h = 2;
  lmatzero(m);
  for (i=0; i<n; i++) {
    l[0] = xy1[i*2+1] - xy2[i*2+1];
    l[1] = xy2[i*2]   - xy1[i*2];
    l[2] = -( (xy1[i*2]-xy2[i*2])*xy2[i*2+1] +
              (xy2[i*2+1]-xy1[i*2+1])*xy2[i*2] );
    lleast(m, l, 1); }
  lmatsym(m);
  if (lrowreduce(m))
    return(0);
  foexy[0] = m->m[2];  foexy[1] = m->m[5];  foexy[2] = 0;
  if (!err)  return(1);
  for (i=0; i<n; i++) {
    A = xy1[i*2+1] - xy2[i*2+1];
    B = xy2[i*2]   - xy1[i*2];
    C = (xy1[i*2]-xy2[i*2])*xy2[i*2+1] +
        (xy2[i*2+1]-xy1[i*2+1])*xy2[i*2];
    if (A || B) {
      e = fabs(A*foexy[0]+B*foexy[1]+C)/sqrt(A*A+B*B);
      se2 += e*e;  s ++;
      if (e>foexy[2])  foexy[2] = e; } }
  if (err==2 && s)
    foexy[2] = 3*sqrt(se2/s);
  return(1);
}

long geo_to_dlt(long n, float *xyz, float *xy, float *prel, float *dlt,
                long distort)
/* Convert a set of nodes known both in 2D and 3D to a set of 11 dlt
 *  parameters, two centering values, and up to 5 distortion parameters.
 *  Note that the centering values are actually included in the dlt values,
 *  but are stored separately for convenience.
 * Enter: long n: number of nodes.
 *        float *xyz: three-dimensional coordinates.  Stored x1,y1,z1,x2,y2,
 *                    z2,x3,...
 *        float *xy: two-dimensional pixel coordinates.  Stored x1,y1,x2,...
 *        float *prel: reliability of each x and y coordinate.  Must
 *                     correspond with xy array.  Stored wx1,wy1,wx2,wy2,wx3,
 *                     ...  Pass uniform values (such as 1) if reliability is
 *                     unknown.  This pointer may be set to zero to not use
 *                     point reliability.
 *        float *dlt: location for dlt/distortion.  Must hold 18 floats.
 *        long distort: Distortion parameters to solve for.  Note that
 *                      solving for distortion is iterative and slow.  Bits
 *                      0-1: Number of Kappa (radial) values to solve for
 *                      (0-3).  Bit 2: Solve for P (asymmetric) values.
 * Exit:  long okay: 0 for failed, 1 for okay.                 2/12/96-DWM */
{
  double mc[274], l[17], x, y, X, Y, Z, A, x0, y0, r2, xp, yp, d[5], fac;
  lmatrix *m=lmat(mc);
  long i, j, k, dist, darr[5]={0,0,0,0,0};

  m->w = 12;  m->h = 11;               /* Determine initial DLT parameters */
  lmatzero(m);
  for (i=0; i<n; i++) {
    x = xy[i*2];  y = xy[i*2+1];
    X = xyz[i*3];  Y = xyz[i*3+1];  Z = xyz[i*3+2];
    if (x>1e29 || y>1e29 || X>1e29 || Y>1e29 || Z>1e29)  continue;
    l[0] = X;  l[1] = Y;  l[2] = Z;
    l[3]=1; l[4]=l[5]=l[6]=l[7]=0; l[8]=-x*l[0]; l[9]=-x*l[1];
    l[10]=-x*l[2]; l[11]=x;
    fac = 1;
    if (prel)  if (prel[i*2])  fac = prel[i*2];
    lleast(m, l, fac);
    l[4]=l[0]; l[5]=l[1]; l[6]=l[2]; l[7]=1; l[0]=l[1]=l[2]=l[3]=0;
    l[8]=-y*l[4]; l[9]=-y*l[5]; l[10]=-y*l[6]; l[11]=y;
    fac = 1;
    if (prel)  if (prel[i*2+1])  fac = prel[i*2+1];
    lleast(m, l, fac); }
  lmatsym(m);
  if (lrowreduce(m))
    return(0);
  for (i=0; i<11; i++)
    dlt[i] = m->m[i*12+11];
  memset(dlt+11, 0, 7*sizeof(float));
  A = sq(dlt[8])+sq(dlt[9])+sq(dlt[10]);
  x0 = (dlt[0]*dlt[8]+dlt[1]*dlt[9]+dlt[2]*dlt[10])/A;
  y0 = (dlt[4]*dlt[8]+dlt[5]*dlt[9]+dlt[6]*dlt[10])/A;
  for (i=0; i<(distort&3); i++)  darr[i] = 1;
  if (distort&4)       darr[3] = darr[4] = 1;
  for (i=dist=0; i<5; i++)  dist += darr[i];
  memset(dlt+13, 0, 5*sizeof(float));
  dlt[11] = x0;  dlt[12] = y0;
  for (j=0; j<100 && dist; j++) {
    m->w = 12+dist;  m->h = 11+dist;
    lmatzero(m);
    for (i=0; i<n; i++) {
      x = xy[i*2];  y = xy[i*2+1];
      X = xyz[i*3];  Y = xyz[i*3+1];  Z = xyz[i*3+2];
      if (x>1e29 || y>1e29 || X>1e29 || Y>1e29 || Z>1e29)  continue;
      xp = x-dlt[11];  yp = y-dlt[12];
      r2 = xp*xp+yp*yp;
      A = dlt[8]*X+dlt[9]*Y+dlt[10]*Z+1;
      l[0] = X;  l[1] = Y;  l[2] = Z;
      l[3]=1; l[4]=l[5]=l[6]=l[7]=0; l[8]=-x*l[0]; l[9]=-x*l[1];
      l[10]=-x*l[2];
      d[0]=-xp*r2*A; d[1]=-xp*r2*r2*A; d[2]=-xp*r2*r2*r2*A;
      d[3]=-(r2+2*xp*xp)*A; d[4]=-2*xp*yp*A;
      for (k=dist=0; k<5; k++) if (darr[k]) {
        l[11+dist] = d[k];  dist++; }
      l[11+dist] = x;
      lleast(m, l, 1);
      l[4]=l[0]; l[5]=l[1]; l[6]=l[2]; l[7]=1; l[0]=l[1]=l[2]=l[3]=0;
      l[8]=-y*l[4]; l[9]=-y*l[5]; l[10]=-y*l[6];
      d[0]=-yp*r2*A; d[1]=-yp*r2*r2*A; d[2]=-yp*r2*r2*r2*A;
      d[3]=-2*xp*yp*A; d[4]=-(r2+2*yp*yp)*A;
      for (k=dist=0; k<5; k++) if (darr[k]) {
        l[11+dist] = d[k];  dist++; }
      l[11+dist] = y;
      lleast(m, l, 1); }
    lmatsym(m);
    if (lrowreduce(m)) {
      j = 200; continue; }
    for (i=0; i<11; i++)
      dlt[i] = m->m[i*(12+dist)+11+dist];
    for (k=i=0; k<5; k++) if (darr[k]) {
      dlt[13+i] = m->m[(i+11)*(12+dist)+11+dist];  i++; }
    A = sq(dlt[8])+sq(dlt[9])+sq(dlt[10]);
    x = (dlt[0]*dlt[8]+dlt[1]*dlt[9]+dlt[2]*dlt[10])/A;
    y = (dlt[4]*dlt[8]+dlt[5]*dlt[9]+dlt[6]*dlt[10])/A;
    if (j>10 && fabs(x-dlt[11])+fabs(y-dlt[12])<1e-5)
      return(1);
    dlt[11] = x;  dlt[12] = y; }
  return(1);
}

void geo_to_image(long n, float *xyz, float *dlt, float *xy)
/* Compute the image (2D) coordinates for a point in 3D for a particular
 *  camera.  Currently distortion is not accounted for.
 * Enter: long n: number of points to compute.
 *        float *xyz: 3D coordinates, stored X1,Y1,Z1,X2,Y2,Z2,X3,...
 *        float *dlt: 18 dlt/distortion parameters for camera.
 *        float *xy: array to store calculated 2D data, stored x1,y1,x2,y2,
 *                   x3,...  Erroneous values are stored as 1e30.
 *                                                             2/12/96-DWM */
{
  long i;
  double X, Y, Z, den;

  for (i=0; i<n; i++) {
    X = xyz[i*3];  Y = xyz[i*3+1];  Z = xyz[i*3+2];
    den = dlt[8]*X + dlt[9]*Y + dlt[10]*Z + 1;
    if (!den) { xy[i*2] = xy[i*2+1] = 1e30;  continue; }
    xy[i*2]   = (dlt[0]*X + dlt[1]*Y + dlt[2]*Z + dlt[3])/den;
    xy[i*2+1] = (dlt[4]*X + dlt[5]*Y + dlt[6]*Z + dlt[7])/den; }
}

void image_to_geo(long n2, long p, long minp, short *np, float *xy,
                  float *dlt, float *rel, float *prel, long n, float *xyz)
/* Calculate the three dimensional points from a set of two dimensional
 *  points.  The dlt array contains L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,x0,y0,
 *  K1,K2,K3,P1,P2 for each photograph.
 * Enter: long n2: number of 2D points.  This includes all photos (i.e., if
 *                 each of two photos has 6 points, n is 12).
 *        long p: number of photographs.
 *        long minp: minimum required number of photographs to calculate a
 *                   3D point.  If minp<2 then minp=2; if minp>p, minp=p.
 *        short *np: node and photograph number of a point.  Stored n1,p1,n2,
 *                   p2,n3,...  The node number corresponds with the final
 *                   3D node number of the part.  The photograph number must
 *                   match that of the dlt array.  Any point with a photo
 *                   number outside of the valid range is ignored.
 *        float *xy: 2D coordinate array.  Must correspond with np array.
 *                   Stored x1,y1,x2,y2,x3,...
 *        float *dlt: dlt/distortion parameter array for all photographs.
 *                    Stored L11, L21, ..., P21, L12, L22, ..., P22, L13, ...
 *        float *rel: reliability of each photograph.  Pass uniform values
 *                    (such as 1) if reliability is unknown.  This pointer
 *                    may be set to zero to not use reliability.
 *        float *prel: reliability of each x and y coordinate.  Must
 *                     correspond with np array.  Stored wx1,wy1,wx2,wy2,wx3,
 *                     ...  Pass uniform values (such as 1) if reliability is
 *                     unknown.  This pointer may be set to zero to not use
 *                     point reliability.
 *        long n: number of 3D nodes to calculate.
 *        float *xyz: location to store 3D coordinates.  Stored x1,y1,z1,x2,
 *                    y2,z2,x3,...  Values which were not calculated are set
 *                    with x=1e30.                             2/12/96-DWM */
{
  double mc[14], l[4], xp, yp, fac, fac2, cx, cy, dx, dy, r2;
  lmatrix *m=lmat(mc);
  long i, k, t, pp;

  if (minp<2)  minp = 2;  if (minp>p)  minp = p;
  m->w = 4;  m->h = 3;
  for (i=0; i<n; i++) {
    lmatzero(m);
    for (k=t=0; k<n2; k++)
      if (np[k*2]==i) {
        t++;
        xp = xy[k*2];  yp = xy[k*2+1];  pp = np[k*2+1];
        if (!rel)  fac = 1;
        else       fac = rel[pp];
        if (!fac)  fac = 1;
        cx = xp+dlt[pp*18+11];  cy = yp+dlt[pp*18+12];
        r2 = cx*cx+cy*cy;
        dx = cx*(dlt[pp*18+13]*r2+dlt[pp*18+14]*r2*r2+dlt[pp*18+15]*r2*r2*r2)
             + dlt[pp*18+16]*(r2+2*cx*cx) + 2*dlt[pp*18+17]*cx*cy;
        dy = cy*(dlt[pp*18+13]*r2+dlt[pp*18+14]*r2*r2+dlt[pp*18+15]*r2*r2*r2)
             + 2*dlt[pp*18+16]*cx*cy + dlt[pp*18+17]*(r2+2*cy*cy);
        xp += dx;  yp += dy;
        l[0] = dlt[pp*18  ]-xp*dlt[pp*18+8];
        l[1] = dlt[pp*18+1]-xp*dlt[pp*18+9];
        l[2] = dlt[pp*18+2]-xp*dlt[pp*18+10];
        l[3] = xp-dlt[pp*18+3]-dx;
        fac2 = fac;
        if (prel)  if (prel[k*2])  fac2 = fac * prel[k*2];
        lleast(m, l, fac2);
        l[0] = dlt[pp*18+4]-yp*dlt[pp*18+8];
        l[1] = dlt[pp*18+5]-yp*dlt[pp*18+9];
        l[2] = dlt[pp*18+6]-yp*dlt[pp*18+10];
        l[3] = yp-dlt[pp*18+7]-dy;
        fac2 = fac;
        if (prel)  if (prel[k*2+1])  fac2 = fac * prel[k*2+1];
        lleast(m, l, fac2); }
    if (t<minp) {
      xyz[i*3] = xyz[i*3+1] = xyz[i*3+2] = 1e30;
      continue; }
    lmatsym(m);
    if (!lrowreduce(m)) {
      xyz[i*3]   = m->m[3];
      xyz[i*3+1] = m->m[7];
      xyz[i*3+2] = m->m[11]; }
    else
      xyz[i*3] = xyz[i*3+1] = xyz[i*3+2] = 1e30; }
}

void image2_to_geo(long n, float *xy1, float *xy2, float *dlt1, float *dlt2,
				   float *xyz)
/* Calculate the three dimensional points from two sets of two dimensional
 *  points.  The points from each photograph must correspond exactly.
 * Enter: long n: number of 3D points to calculate.
 *        float *xy1, *xy2: 2D coordinate array for first, second
 *                           photographs.  Stored x1,y1,x2,y2,x3,...
 *        float *dlt1, *dlt2: dlt/distortion parameter array for first,
 *                             second photographs.  Stored L11, L21, ...,
 *                             P21, L12, L22, ..., P22, L13, ...
 *        float *xyz: location to store 3D coordinates.  Stored x1,y1,z1,x2,
 *                     y2,z2,x3,...  Values which were not calculated are set
 *                     with x=1e30.                            2/12/96-DWM */
{
  double mc[14], l[4], xp, yp, cx, cy, dx, dy, r2;
  float *xy, *dlt;
  lmatrix *m=lmat(mc);
  long i, k;

  m->w = 4;  m->h = 3;
  for (i=0; i<n; i++) {
	lmatzero(m);
	for (k=0; k<2; k++) {
	  if (k==0) { xy = xy1;  dlt = dlt1; }
	  else      { xy = xy2;  dlt = dlt2; }
	  xp = xy[i*2];  yp = xy[i*2+1];
	  cx = xp+dlt[11];  cy = yp+dlt[12];  r2 = cx*cx+cy*cy;
	  dx = cx*(dlt[13]*r2+dlt[14]*r2*r2+dlt[15]*r2*r2*r2) +
		   dlt[16]*(r2+2*cx*cx) + 2*dlt[17]*cx*cy;
	  dy = cy*(dlt[13]*r2+dlt[14]*r2*r2+dlt[15]*r2*r2*r2) + 2*dlt[16]*cx*cy +
		   dlt[17]*(r2+2*cy*cy);
	  xp += dx;  yp += dy;
	  l[0] = dlt[0]-xp*dlt[8];  l[1] = dlt[1]-xp*dlt[9];
	  l[2] = dlt[2]-xp*dlt[10]; l[3] = xp-dlt[3]-dx;
	  lleast(m, l, 1);
	  l[0] = dlt[4]-yp*dlt[8];  l[1] = dlt[5]-yp*dlt[9];
	  l[2] = dlt[6]-yp*dlt[10]; l[3] = yp-dlt[7]-dy;
	  lleast(m, l, 1); }
	lmatsym(m);
	if (!lrowreduce(m)) {
	  xyz[i*3]   = m->m[3];
	  xyz[i*3+1] = m->m[7];
	  xyz[i*3+2] = m->m[11]; }
	else
	  xyz[i*3] = xyz[i*3+1] = xyz[i*3+2] = 1e30; }
}
#endif

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

#ifdef used
long lepipolar_line(double *e12, lmatrix *P2P1i, double *xy, double *abc)
/* Calculate the epipolar line in image 2 based on an already located point
 *  in image 1.
 * Enter: double *e12: epipolar point of camera 1 in image 2.  Array of 2
 *                     doubles.
 *        lmatrix *P2P1i: P2 P1^-1.  This is a 3x3 matrix as calculated by
 *                        the function epipole().
 *        double *xy: point located in image 1.  Array of 2 doubles.
 *        double *abc: location to store resulting line.  This is an array of
 *                     3 doubles, in terms of ax+by+c=0.
 * Exit:  long okay: 0 for failed, 1 for okay.                 2/13/96-DWM */
{
  double mc[5], pc[5], p2[2];
  lmatrix *m=lmat(mc), *p=lmat(pc);

  m->w = 1;  m->h = 3;
  m->m[0] = xy[0];  m->m[1] = xy[1];  m->m[2] = 1;
  if (lmatmul(P2P1i, m, p))  return(0);
  if (!p->m[2])              return(0);
  p2[0] = p->m[0]/p->m[2];  p2[1] = p->m[1]/p->m[2];
  abc[0] = p2[1] - e12[1];
  abc[1] = e12[0] - p2[0];
  abc[2] = (p2[0]-e12[0])*e12[1] + (e12[1]-p2[1])*e12[0];
  return(1);
}

long lepipolar_segment(double *dlt1, double *dlt2, double *e12,
                    lmatrix *P2P1i, double *xy, double *bounds, double *xyxy)
/* Calculate the epipolar line in image 2 based on an already located point
 *  in image 1.  Only a segment of this line is returned, clipped from the
 *  entire line based on the destination image and some maximum distance from
 *  each camera.  The bounds array contains the screen coordinates (left,
 *  top, right, bottom) and the distance limits (photo 1 maximum, photo 2
 *  maximum).
 * Enter: double *dlt1: 18 dlt parameters/distortion for camera 1.
 *        double *dlt2: 18 dlt parameters/distortion for camera 2.
 *        double *e12: epipolar point of camera 1 in image 2.  Array of 2
 *                    doubles.  If null, epipolar_segment will calculate the
 *                    epipole on its own.
 *        lmatrix *P2P1i: P2 P1^-1.  This is a 3x3 matrix as calculated by
 *                        the function epipole().  If null, epipolar_segment
 *                        will calculate the epipole on its own.
 *        double *xy: point located in image 1.  Array of 2 doubles.
 *        double *bounds: boundary conditions for the line cropping.  This is
 *                       an array of 6 values: x1,y1,x2,y2,max1,max2.  If
 *                       the maximum values are <=0, no distance cropping is
 *                       performed.
 *        double *xyxy: location to store resulting line segment.  This is an
 *                     array of 4 doubles, in terms of x1,y1,x2,y2.  x1 will
 *                     always less than or equal to x2.  If equal to, y1 will
 *                     be less than y2.
 * Exit:  long okay: 0 for failed, 1 for okay.                 2/13/96-DWM */
{
  double P2P1c[11], mc[5], pc[5], dist, dx[3], r, oc[11];
  double abc[3], ae12[2], xyc[8], xyz[6], dltc[36], phys[22], p2[2];
  lmatrix *aP2P1i=lmat(P2P1c), *m=lmat(mc), *p=lmat(pc), *o=lmat(oc);
  long i, j, k, side=1, use;
  short np[]={0,0,0,1,1,0,1,1};

  if (e12 && P2P1i) {
    memcpy(ae12, e12, 2*sizeof(double));
    lmatdup(P2P1i, aP2P1i); }
  else {
    if (!lepipole(dlt1, dlt2, ae12, aP2P1i))  return(0); }
  m->w = 1;  m->h = 3;
  m->m[0] = xy[0];  m->m[1] = xy[1];  m->m[2] = 1;
  if (lmatmul(aP2P1i, m, p))  return(0);
  if (!p->m[2])               return(0);
  p2[0] = p->m[0]/p->m[2];  p2[1] = p->m[1]/p->m[2];
  abc[0] = p2[1] - ae12[1];
  abc[1] = ae12[0] - p2[0];
  abc[2] = (p2[0]-ae12[0])*ae12[1] + (ae12[1]-p2[1])*ae12[0];
  if (!abc[0] && !abc[1])  return(0);
  ldlt_to_physical(dlt1, phys, 0, 0);   /* Select other side of proj. inf. */
  ldlt_to_physical(dlt2, phys+11, o, 0);
  if (o->m[2]*(phys[0]-phys[11])+o->m[5]*(phys[1]-phys[12])+
	  o->m[8]*(phys[2]-phys[13])>0) {
    ae12[0] = p2[0] + (p2[0]-ae12[0])*1000;
    ae12[1] = p2[1] + (p2[1]-ae12[1])*1000; }                  /* End swap */
  if (fabs(p2[0]-ae12[0])>fabs(p2[1]-ae12[1])) { /* Back off a pixel since */
    r = 1-2*(p2[0]>ae12[0]);                     /*  inf. is far away.     */
    p2[0] += r;  p2[1] += r*(p2[1]-ae12[1])/(p2[0]-ae12[0]); }
  else {
    r = 1-2*(p2[1]>ae12[1]);
    p2[1] += r;  p2[0] += r*(p2[0]-ae12[0])/(p2[1]-ae12[1]); }
  xyxy[0] = ae12[0];  xyxy[1] = ae12[1];
  xyxy[2] = p2[0];    xyxy[3] = p2[1];
  if (xyxy[0]>xyxy[2] || (xyxy[0]==xyxy[2] && xyxy[1]>xyxy[3])) {
    xyxy[0] = p2[0];    xyxy[1] = p2[1];
    xyxy[2] = ae12[0];  xyxy[3] = ae12[1];
    side = 0; }
  if (xyxy[0]>bounds[2] || xyxy[2]<bounds[0] || (xyxy[0]<bounds[1] &&
      xyxy[2]<bounds[1]) || (xyxy[0]>bounds[3] && xyxy[2]>bounds[3]))
    return(0);
  if (abc[1]) {
    if (xyxy[0]<bounds[0]) {
      xyxy[0] = bounds[0];  xyxy[1] = (-abc[2]-abc[0]*bounds[0])/abc[1]; }
    if (xyxy[2]>bounds[2]) {
      xyxy[2] = bounds[2];  xyxy[3] = (-abc[2]-abc[0]*bounds[2])/abc[1]; }
    if ((xyxy[1]<bounds[1] && xyxy[3]<bounds[1]) ||
        (xyxy[1]>bounds[3] && xyxy[3]>bounds[3]))  return(0);
    if (xyxy[1]<bounds[1]) {
      xyxy[1] = bounds[1];  xyxy[0] = (-abc[2]-abc[1]*bounds[1])/abc[0]; }
    if (xyxy[3]<bounds[1]) {
      xyxy[3] = bounds[1];  xyxy[2] = (-abc[2]-abc[1]*bounds[1])/abc[0]; }
    if (xyxy[1]>bounds[3]) {
      xyxy[1] = bounds[3];  xyxy[0] = (-abc[2]-abc[1]*bounds[3])/abc[0]; }
    if (xyxy[3]>bounds[3]) {
      xyxy[3] = bounds[3];  xyxy[2] = (-abc[2]-abc[1]*bounds[3])/abc[0]; } }
  else {
    if (xyxy[1]<bounds[1]) {
      xyxy[1] = bounds[1];  xyxy[0] = (-abc[2]-abc[1]*bounds[1])/abc[0]; }
    if (xyxy[3]>bounds[3]); {
      xyxy[3] = bounds[3];  xyxy[2] = (-abc[2]-abc[1]*bounds[3])/abc[0]; }
    if ((xyxy[0]<bounds[0] && xyxy[2]<bounds[0]) ||
        (xyxy[0]>bounds[2] && xyxy[2]>bounds[2]))  return(0); }
  for (k=0; k<2; k++)  if (bounds[k+4]) {
    memcpy(xyc,   xy, 2*sizeof(double));
    memcpy(xyc+2, xyxy, 2*sizeof(double));
    memcpy(xyc+4, xy, 2*sizeof(double));
    memcpy(xyc+6, xyxy+2, 2*sizeof(double));
    memcpy(dltc,    dlt1, 18*sizeof(double));
    memcpy(dltc+18, dlt2, 18*sizeof(double));
    limage_to_geo(4, 2, 2, np, xyc, dltc, 0, 0, 2, xyz);
    ldlt_to_physical(dltc+k*18, phys, 0, 0);
    dx[0] = sqrt(sq(phys[0]-xyz[0])+sq(phys[1]-xyz[1])+sq(phys[2]-xyz[2]));
    dx[1] = sqrt(sq(phys[0]-(xyz[0]+xyz[3])/2)+sq(phys[1]-(xyz[1]+xyz[4])/2)+
                 sq(phys[2]-(xyz[2]+xyz[5])/2));
    dx[2] = sqrt(sq(phys[0]-xyz[3])+sq(phys[1]-xyz[4])+sq(phys[2]-xyz[5]));
    if (k) {
      abc[0] = 2*dx[0]-4*dx[1]+2*dx[2];
      abc[1] = -3*dx[0]+4*dx[1]-dx[2];
      abc[2] = dx[0];
      dx[0] = 4*abc[0]*bounds[k+4]-4*abc[0]*abc[2]+abc[1]*abc[1];
      if (dx[0]<0)  return(0);
      dx[0] = sqrt(dx[0]);
      dx[1] = (dx[0]-abc[1])/(2*abc[0]);
      dx[0] = (-dx[0]-abc[1])/(2*abc[0]);
      if ((dx[0]<0 && dx[1]<0) || (dx[0]>1 && dx[1]>1))  return(0); }
    else {
      abc[0] = dx[2]-dx[0];  abc[1] = dx[0];
      if (!abc[0])  continue;
      dx[0] = 0;  dx[1] = 1;
      dx[abc[0]>0] = (bounds[k+4]-abc[1])/abc[0];
      if (dx[0]>1 || dx[1]<0)  return(0); }
    for (j=0; j<2; j++)
      if (dx[j]>0 && dx[j]<1) {
        for (i=0; i<3; i++)
          abc[i] = xyz[i]*(1-dx[j])+xyz[3+i]*dx[j];
        lgeo_to_image(1, abc, dlt2, p2);
        if (p2[0]<1e20)
          memcpy(xyxy+j*2, p2, 2*sizeof(double)); } }
  return(1);
}

long lepipole(double *dlt1, double *dlt2, double *e12, lmatrix *P2P1i)
/* Calculate the epipole of camera 1 in image 2.  This is used along with the
 *  function epipolar_line when a point is known in image 1 and is being
 *  searched for in image 2.
 * Enter: double *dlt1: 11 dlt parameters for camera 1.
 *        double *dlt2: 11 dlt parameters for camera 2.
 *        double *e12: epipole of camera 1 in image 2.  Array of 2 doubles.
 *        lmatrix *P2P1i: P2 P1^-1, used for epipolar line calculation.  This
 *                        is a 3x3 matrix.
 * Exit:  long okay: 0 for failed, 1 for okay.                 2/13/96-DWM */
{
  double P1c[11], P2c[11], P1ic[11], p1c[5], p2c[5], tempc[5];
  lmatrix *P1=lmat(P1c), *P2=lmat(P2c), *P1i=lmat(P1ic);
  lmatrix *p1=lmat(p1c), *p2=lmat(p2c), *temp=lmat(tempc);

  P1->w = P1->h = P2->w = P2->h = 3;
  P1->m[0] = dlt1[0];  P1->m[1] = dlt1[1];  P1->m[2] = dlt1[2];
  P1->m[3] = dlt1[4];  P1->m[4] = dlt1[5];  P1->m[5] = dlt1[6];
  P1->m[6] = dlt1[8];  P1->m[7] = dlt1[9];  P1->m[8] = dlt1[10];
  P2->m[0] = dlt2[0];  P2->m[1] = dlt2[1];  P2->m[2] = dlt2[2];
  P2->m[3] = dlt2[4];  P2->m[4] = dlt2[5];  P2->m[5] = dlt2[6];
  P2->m[6] = dlt2[8];  P2->m[7] = dlt2[9];  P2->m[8] = dlt2[10];
  p1->w = p2->w = 1;  p1->h = p2->h = 3;
  p1->m[0] = dlt1[3];  p1->m[1] = dlt1[7];  p1->m[2] = 1;
  p2->m[0] = dlt2[3];  p2->m[1] = dlt2[7];  p2->m[2] = 1;
  if (lmatinv(P1, P1i))               return(0);
  if (lmatmul(P2, P1i, P2P1i))        return(0);
  if (lmatmul(P2P1i, p1, temp))       return(0);
  if (lmatadd(-1, temp, 1, p2, temp)) return(0);
  if (!temp->m[2])                    return(0);
  e12[0] = temp->m[0]/temp->m[2];  e12[1] = temp->m[1]/temp->m[2];
  return(1);
}

long lfoe(long n, double *xy1, double *xy2, long err, double *foexy)
/* Calculate the focus of expansion based on a set of corresponding points
 *  between two images.  The located focus of expansion is returned.  If
 *  requested, the maximal or statistical error is given.
 * Enter: long n: number of points in each image to use.  Must be at least 2.
 *        double *xy1: points in first image in order x1,y1,x2,y2,x3,...
 *        double *xy2: points in second image in order x1,y1,x2,y2,x3,...
 *                    These must have a direct one-to-one correspondence with
 *                    the points from image 1.
 *        long err: 0 for error value equal to zero, 1 for maximal error, 2
 *                  for statistical (3 sigma) error.
 *        double *foexy: location to store resultant foe in image coordinates.
 *                      This is an array of 3 doubles: x,y,error.  If error is
 *                      not calculated, the third value is set to zero.
 * Exit:  long okay: 0 for failed, 1 for okay.                 2/21/96-DWM */
{
  double mc[8], l[3], A, B, C, e, se2=0, s=0;
  lmatrix *m=lmat(mc);
  long i;

  if (n<2)  return(0);
  m->w = 3;  m->h = 2;
  lmatzero(m);
  for (i=0; i<n; i++) {
    l[0] = xy1[i*2+1] - xy2[i*2+1];
    l[1] = xy2[i*2]   - xy1[i*2];
    l[2] = -( (xy1[i*2]-xy2[i*2])*xy2[i*2+1] +
              (xy2[i*2+1]-xy1[i*2+1])*xy2[i*2] );
    lleast(m, l, 1); }
  lmatsym(m);
  if (lrowreduce(m))
    return(0);
  foexy[0] = m->m[2];  foexy[1] = m->m[5];  foexy[2] = 0;
  if (!err)  return(1);
  for (i=0; i<n; i++) {
    A = xy1[i*2+1] - xy2[i*2+1];
    B = xy2[i*2]   - xy1[i*2];
    C = (xy1[i*2]-xy2[i*2])*xy2[i*2+1] +
        (xy2[i*2+1]-xy1[i*2+1])*xy2[i*2];
    if (A || B) {
      e = fabs(A*foexy[0]+B*foexy[1]+C)/sqrt(A*A+B*B);
      se2 += e*e;  s ++;
      if (e>foexy[2])  foexy[2] = e; } }
  if (err==2 && s)
    foexy[2] = 3*sqrt(se2/s);
  return(1);
}

long lgeo_to_dlt(long n, double *xyz, double *xy, double *prel, double *dlt,
                 long distort)
/* Convert a set of nodes known both in 2D and 3D to a set of 11 dlt
 *  parameters, two centering values, and up to 5 distortion parameters.
 *  Note that the centering values are actually included in the dlt values,
 *  but are stored separately for convenience.
 * Enter: long n: number of nodes.
 *        double *xyz: three-dimensional coordinates.  Stored x1,y1,z1,x2,y2,
 *                     z2,x3,...
 *        double *xy: two-dimensional pixel coordinates.  Stored x1,y1,x2,...
 *        double *prel: reliability of each x and y coordinate.  Must
 *                      correspond with xy array.  Stored wx1,wy1,wx2,wy2,
 *                      wx3,...  Pass uniform values (such as 1) if
 *                      reliability is unknown.  This pointer may be set to
 *                      zero to not use point reliability.
 *        double *dlt: location for dlt/distortion.  Must hold 18 doubles.
 *        long distort: Distortion parameters to solve for.  Note that
 *                      solving for distortion is iterative and slow.  Bits
 *                      0-1: Number of Kappa (radial) values to solve for
 *                      (0-3).  Bit 2: Solve for P (asymmetric) values.
 * Exit:  long okay: 0 for failed, 1 for okay.                 2/12/96-DWM */
{
  double mc[274], l[17], x, y, X, Y, Z, A, x0, y0, r2, xp, yp, d[5], fac;
  lmatrix *m=lmat(mc);
  long i, j, k, dist, darr[5]={0,0,0,0,0};

  m->w = 12;  m->h = 11;               /* Determine initial DLT parameters */
  lmatzero(m);
  for (i=0; i<n; i++) {
    x = xy[i*2];  y = xy[i*2+1];
    X = xyz[i*3];  Y = xyz[i*3+1];  Z = xyz[i*3+2];
    if (x>1e29 || y>1e29 || X>1e29 || Y>1e29 || Z>1e29)  continue;
    l[0] = X;  l[1] = Y;  l[2] = Z;
    l[3]=1; l[4]=l[5]=l[6]=l[7]=0; l[8]=-x*l[0]; l[9]=-x*l[1];
    l[10]=-x*l[2]; l[11]=x;
    fac = 1;
    if (prel)  if (prel[i*2])  fac = prel[i*2];
    lleast(m, l, fac);
    l[4]=l[0]; l[5]=l[1]; l[6]=l[2]; l[7]=1; l[0]=l[1]=l[2]=l[3]=0;
    l[8]=-y*l[4]; l[9]=-y*l[5]; l[10]=-y*l[6]; l[11]=y;
    fac = 1;
    if (prel)  if (prel[i*2+1])  fac = prel[i*2+1];
    lleast(m, l, fac); }
  lmatsym(m);
  if (lrowreduce(m))
    return(0);
  for (i=0; i<11; i++)
    dlt[i] = m->m[i*12+11];
  memset(dlt+11, 0, 7*sizeof(double));
  A = sq(dlt[8])+sq(dlt[9])+sq(dlt[10]);
  x0 = (dlt[0]*dlt[8]+dlt[1]*dlt[9]+dlt[2]*dlt[10])/A;
  y0 = (dlt[4]*dlt[8]+dlt[5]*dlt[9]+dlt[6]*dlt[10])/A;
  for (i=0; i<(distort&3); i++)  darr[i] = 1;
  if (distort&4)       darr[3] = darr[4] = 1;
  for (i=dist=0; i<5; i++)  dist += darr[i];
  memset(dlt+13, 0, 5*sizeof(double));
  dlt[11] = x0;  dlt[12] = y0;
  for (j=0; j<100 && dist; j++) {
    m->w = 12+dist;  m->h = 11+dist;
    lmatzero(m);
    for (i=0; i<n; i++) {
      x = xy[i*2];  y = xy[i*2+1];
      X = xyz[i*3];  Y = xyz[i*3+1];  Z = xyz[i*3+2];
      if (x>1e29 || y>1e29 || X>1e29 || Y>1e29 || Z>1e29)  continue;
      xp = x-dlt[11];  yp = y-dlt[12];
      r2 = xp*xp+yp*yp;
      A = dlt[8]*X+dlt[9]*Y+dlt[10]*Z+1;
      l[0] = X;  l[1] = Y;  l[2] = Z;
      l[3]=1; l[4]=l[5]=l[6]=l[7]=0; l[8]=-x*l[0]; l[9]=-x*l[1];
      l[10]=-x*l[2];
      d[0]=-xp*r2*A; d[1]=-xp*r2*r2*A; d[2]=-xp*r2*r2*r2*A;
      d[3]=-(r2+2*xp*xp)*A; d[4]=-2*xp*yp*A;
      for (k=dist=0; k<5; k++) if (darr[k]) {
        l[11+dist] = d[k];  dist++; }
      l[11+dist] = x;
      lleast(m, l, 1);
      l[4]=l[0]; l[5]=l[1]; l[6]=l[2]; l[7]=1; l[0]=l[1]=l[2]=l[3]=0;
      l[8]=-y*l[4]; l[9]=-y*l[5]; l[10]=-y*l[6];
      d[0]=-yp*r2*A; d[1]=-yp*r2*r2*A; d[2]=-yp*r2*r2*r2*A;
      d[3]=-2*xp*yp*A; d[4]=-(r2+2*yp*yp)*A;
      for (k=dist=0; k<5; k++) if (darr[k]) {
        l[11+dist] = d[k];  dist++; }
      l[11+dist] = y;
      lleast(m, l, 1); }
    lmatsym(m);
    if (lrowreduce(m)) {
      j = 200; continue; }
    for (i=0; i<11; i++)
      dlt[i] = m->m[i*(12+dist)+11+dist];
    for (k=i=0; k<5; k++) if (darr[k]) {
      dlt[13+i] = m->m[(i+11)*(12+dist)+11+dist];  i++; }
    A = sq(dlt[8])+sq(dlt[9])+sq(dlt[10]);
    x = (dlt[0]*dlt[8]+dlt[1]*dlt[9]+dlt[2]*dlt[10])/A;
    y = (dlt[4]*dlt[8]+dlt[5]*dlt[9]+dlt[6]*dlt[10])/A;
    if (j>10 && fabs(x-dlt[11])+fabs(y-dlt[12])<1e-5)
      return(1);
    dlt[11] = x;  dlt[12] = y; }
  return(1);
}

void lgeo_to_image(long n, double *xyz, double *dlt, double *xy)
/* Compute the image (2D) coordinates for a point in 3D for a particular
 *  camera.  Currently distortion is not accounted for.
 * Enter: long n: number of points to compute.
 *        double *xyz: 3D coordinates, stored X1,Y1,Z1,X2,Y2,Z2,X3,...
 *        double *dlt: 18 dlt/distortion parameters for camera.
 *        double *xy: array to store calculated 2D data, stored x1,y1,x2,y2,
 *                    x3,...  Erroneous values are stored as 1e30.
 *                                                             2/12/96-DWM */
{
  long i;
  double X, Y, Z, den;

  for (i=0; i<n; i++) {
    X = xyz[i*3];  Y = xyz[i*3+1];  Z = xyz[i*3+2];
    den = dlt[8]*X + dlt[9]*Y + dlt[10]*Z + 1;
    if (!den) { xy[i*2] = xy[i*2+1] = 1e30;  continue; }
    xy[i*2]   = (dlt[0]*X + dlt[1]*Y + dlt[2]*Z + dlt[3])/den;
    xy[i*2+1] = (dlt[4]*X + dlt[5]*Y + dlt[6]*Z + dlt[7])/den; }
}

void limage_to_geo(long n2, long p, long minp, short *np, double *xy,
                 double *dlt, double *rel, double *prel, long n, double *xyz)
/* Calculate the three dimensional points from a set of two dimensional
 *  points.  The dlt array contains L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,x0,y0,
 *  K1,K2,K3,P1,P2 for each photograph.
 * Enter: long n2: number of 2D points.  This includes all photos (i.e., if
 *                 each of two photos has 6 points, n is 12).
 *        long p: number of photographs.
 *        long minp: minimum required number of photographs to calculate a
 *                   3D point.  If minp<2 then minp=2; if minp>p, minp=p.
 *        short *np: node and photograph number of a point.  Stored n1,p1,n2,
 *                   p2,n3,...  The node number corresponds with the final
 *                   3D node number of the part.  The photograph number must
 *                   match that of the dlt array.  Any point with a photo
 *                   number outside of the valid range is ignored.
 *        double *xy: 2D coordinate array.  Must correspond with np array.
 *                    Stored x1,y1,x2,y2,x3,...
 *        double *dlt: dlt/distortion parameter array for all photographs.
 *                     Stored L11, L21, ..., P21, L12, L22, ..., P22, L13,
 *                     ...
 *        double *rel: reliability of each photograph.  Pass uniform values
 *                     (such as 1) if reliability is unknown.  This pointer
 *                     may be set to zero to not use reliability.
 *        double *prel: reliability of each x and y coordinate.  Must
 *                      correspond with np array.  Stored wx1,wy1,wx2,wy2,
 *                      wx3,...  Pass uniform values (such as 1) if
 *                      reliability is unknown.  This pointer may be set to
 *                      zero to not use point reliability.
 *        long n: number of 3D nodes to calculate.
 *        double *xyz: location to store 3D coordinates.  Stored x1,y1,z1,x2,
 *                     y2,z2,x3,...  Values which were not calculated are set
 *                     with x=1e30.                            2/12/96-DWM */
{
  double mc[14], l[4], xp, yp, fac, fac2, cx, cy, dx, dy, r2;
  lmatrix *m=lmat(mc);
  long i, k, t, pp;

  if (minp<2)  minp = 2;  if (minp>p)  minp = p;
  m->w = 4;  m->h = 3;
  for (i=0; i<n; i++) {
    lmatzero(m);
    for (k=t=0; k<n2; k++)
      if (np[k*2]==i) {
        t++;
        xp = xy[k*2];  yp = xy[k*2+1];  pp = np[k*2+1];
        if (!rel)  fac = 1;
        else       fac = rel[pp];
        if (!fac)  fac = 1;
        cx = xp+dlt[pp*18+11];  cy = yp+dlt[pp*18+12];
        r2 = cx*cx+cy*cy;
        dx = cx*(dlt[pp*18+13]*r2+dlt[pp*18+14]*r2*r2+dlt[pp*18+15]*r2*r2*r2)
             + dlt[pp*18+16]*(r2+2*cx*cx) + 2*dlt[pp*18+17]*cx*cy;
        dy = cy*(dlt[pp*18+13]*r2+dlt[pp*18+14]*r2*r2+dlt[pp*18+15]*r2*r2*r2)
             + 2*dlt[pp*18+16]*cx*cy + dlt[pp*18+17]*(r2+2*cy*cy);
        xp += dx;  yp += dy;
        l[0] = dlt[pp*18  ]-xp*dlt[pp*18+8];
        l[1] = dlt[pp*18+1]-xp*dlt[pp*18+9];
        l[2] = dlt[pp*18+2]-xp*dlt[pp*18+10];
        l[3] = xp-dlt[pp*18+3]-dx;
        fac2 = fac;
        if (prel)  if (prel[k*2])  fac2 = fac * prel[k*2];
        lleast(m, l, fac2);
        l[0] = dlt[pp*18+4]-yp*dlt[pp*18+8];
        l[1] = dlt[pp*18+5]-yp*dlt[pp*18+9];
        l[2] = dlt[pp*18+6]-yp*dlt[pp*18+10];
        l[3] = yp-dlt[pp*18+7]-dy;
        fac2 = fac;
        if (prel)  if (prel[k*2+1])  fac2 = fac * prel[k*2+1];
        lleast(m, l, fac2); }
    if (t<minp) {
      xyz[i*3] = xyz[i*3+1] = xyz[i*3+2] = 1e30;
      continue; }
    lmatsym(m);
    if (!lrowreduce(m)) {
      xyz[i*3]   = m->m[3];
      xyz[i*3+1] = m->m[7];
      xyz[i*3+2] = m->m[11]; }
    else
      xyz[i*3] = xyz[i*3+1] = xyz[i*3+2] = 1e30; }
}

void limage2_to_geo(long n, double *xy1, double *xy2, double *dlt1,
					double *dlt2, double *xyz)
/* Calculate the three dimensional points from two sets of two dimensional
 *  points.  The points from each photograph must correspond exactly.
 * Enter: long n: number of 3D points to calculate.
 *        double *xy1, *xy2: 2D coordinate array for first, second
 *                           photographs.  Stored x1,y1,x2,y2,x3,...
 *        double *dlt1, *dlt2: dlt/distortion parameter array for first,
 *                             second photographs.  Stored L11, L21, ...,
 *                             P21, L12, L22, ..., P22, L13, ...
 *        double *xyz: location to store 3D coordinates.  Stored x1,y1,z1,x2,
 *                     y2,z2,x3,...  Values which were not calculated are set
 *                     with x=1e30.                            2/12/96-DWM */
{
  double mc[14], l[4], xp, yp, cx, cy, dx, dy, r2, *xy, *dlt;
  lmatrix *m=lmat(mc);
  long i, k;

  m->w = 4;  m->h = 3;
  for (i=0; i<n; i++) {
	lmatzero(m);
	for (k=0; k<2; k++) {
	  if (k==0) { xy = xy1;  dlt = dlt1; }
	  else      { xy = xy2;  dlt = dlt2; }
	  xp = xy[i*2];  yp = xy[i*2+1];
	  cx = xp+dlt[11];  cy = yp+dlt[12];  r2 = cx*cx+cy*cy;
	  dx = cx*(dlt[13]*r2+dlt[14]*r2*r2+dlt[15]*r2*r2*r2) +
		   dlt[16]*(r2+2*cx*cx) + 2*dlt[17]*cx*cy;
	  dy = cy*(dlt[13]*r2+dlt[14]*r2*r2+dlt[15]*r2*r2*r2) + 2*dlt[16]*cx*cy +
		   dlt[17]*(r2+2*cy*cy);
	  xp += dx;  yp += dy;
	  l[0] = dlt[0]-xp*dlt[8];  l[1] = dlt[1]-xp*dlt[9];
	  l[2] = dlt[2]-xp*dlt[10]; l[3] = xp-dlt[3]-dx;
	  lleast(m, l, 1);
	  l[0] = dlt[4]-yp*dlt[8];  l[1] = dlt[5]-yp*dlt[9];
	  l[2] = dlt[6]-yp*dlt[10]; l[3] = yp-dlt[7]-dy;
	  lleast(m, l, 1); }
	lmatsym(m);
	if (!lrowreduce(m)) {
	  xyz[i*3]   = m->m[3];
	  xyz[i*3+1] = m->m[7];
	  xyz[i*3+2] = m->m[11]; }
	else
	  xyz[i*3] = xyz[i*3+1] = xyz[i*3+2] = 1e30; }
}

double lphoto_reliability(long n, double *xy, double *xyz, double *dlt)
/* Calculate a normalized reliability parameter for a photograph.  The
 *  normalized reliability is the recipricol of the average pixel error in
 *  the image compared to where it 'should' be.
 * Enter: long n: number of nodes.
 *        double *xy: 2D nodal coordinates.  Stored x1,y1,x2,y2,x3,...
 *        double *xyz: 3D nodal coordinates.  Stored x1,y1,z1,x2,y2,z3,x3,...
 *        double *dlt: 18 dlt/distortion parameters for this photograph.
 * Exit:  double rel: reliability of photo.  Higher is better. 2/12/96-DWM */
{
  double x, y, X, Y, Z, calcx, calcy, cx, cy, r2, dx, dy, total=0, den;
  long i, count=0;

  for (i=0; i<n; i++) {
    x = xy[i*2];   y = xy[i*2+1];
    X = xyz[i*3];  Y = xyz[i*3+1];  Z = xyz[i*3+2];
    if (x>1e29 || y>1e29 || X>1e29 || Y>1e29 || Z>1e29)  continue;
    den = dlt[8]*X + dlt[9]*Y + dlt[10]*Z + 1;
    if (den) {
      calcx = (dlt[0]*X + dlt[1]*Y + dlt[2]*Z + dlt[3])/den;
      calcy = (dlt[4]*X + dlt[5]*Y + dlt[6]*Z + dlt[7])/den;
      cx = calcx+dlt[11];  cy = calcy+dlt[12];       /** Sign could be off */
      r2 = cx*cx+cy*cy;
      dx = dlt[16]*(r2+2*cx*cx) + 2*dlt[17]*cx*cy +
           cx*(dlt[13]*r2+dlt[14]*r2*r2+dlt[15]*r2*r2*r2);
      dy = dlt[17]*(r2+2*cy*cy) + 2*dlt[16]*cx*cy +
           cy*(dlt[13]*r2+dlt[14]*r2*r2+dlt[15]*r2*r2*r2);
      calcx -= dx;  calcy -= dy;
      total += sqrt(sq(calcx-x) + sq(calcy-y));
      count++; } }
  if (count)  total /= count;
  else        total = 1e6;
  if (total)  total = 1./total;
  else        total = 1e6;
  return(total);
}
#endif

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

#ifdef used
void lrotate_dlt(double *dlt, double *org, double *axis, double angle,
				 double *findlt)
/* Rotate a camera position about a point and axis by a given angle.  This
 *  takes a set of DLT parameters and produces a new set of DLT parameters.
 * Enter: double *dlt: set of 18 dlt parameters.
 *        double *org: origin of rotation.  This is in the original camera
 *                     coordinate system.  This is an array of 3 values.
 *                     This may be null if the origin is at 0,0,0.
 *        double *axis: normalized axis of rotation.  This may be null if the
 *                      axis is 0,0,1.
 *        double angle: angle of rotation in radians.
 *        double *findlt: location to output 18 dlt parameters.  This may be
 *                        the same as the input.                6/5/97-DWM */
{
  double phys[10], mc[55], v1[3], v2[3], x, y, z;
  lmatrix *m=lmat(mc), *r=lmat(mc+11), *a=lmat(mc+22), *ainv=lmat(mc+33);
  lmatrix *ar=lmat(mc+44);

  memcpy(findlt, dlt, 18*sizeof(double));
  dlt_to_physical(dlt, phys, m, 0);
  if (org) {
	phys[0] -= org[0];  phys[1] -= org[1];  phys[2] -= org[2]; }
  a->w = a->h = ainv->w = ainv->h = r->w = r->h = 3;
  if (axis) {
	v1[0] = v1[1] = v1[2] = 1;
	if (axis[0]==axis[1] && axis[0]==axis[2])
	  v1[0] = -1;
	lcross_product(axis, v1, v2);  lunit(v2, v2);
	lcross_product(v2, axis, v1);  lunit(v1, v1);
	a->m[0] = v1[0];  a->m[1] = v2[0];  a->m[2] = axis[0];
	a->m[3] = v1[1];  a->m[4] = v2[1];  a->m[5] = axis[1];
	a->m[6] = v1[2];  a->m[7] = v2[2];  a->m[8] = axis[2];
	lmatdup(a, r);  lmatinv(r, ainv); }
  else {
	lmatident(a, 3);  lmatident(ainv, 3); }
  r->m[0] = cos(angle);  r->m[1] = -sin(angle);  r->m[2] = 0;
  r->m[3] = sin(angle);  r->m[4] =  cos(angle);  r->m[5] = 0;
  r->m[6] = 0;           r->m[7] =  0;           r->m[8] = 1;
  lmatmul(a, r, ar);  lmatmul(ar, ainv, r);
  x = phys[0]*r->m[0]+phys[1]*r->m[3]+phys[2]*r->m[6];
  y = phys[0]*r->m[1]+phys[1]*r->m[4]+phys[2]*r->m[7];
  z = phys[0]*r->m[2]+phys[1]*r->m[5]+phys[2]*r->m[8];
  phys[0] = x;  phys[1] = y;  phys[2] = z;
  lmatmul(m, r, a);
  lmattoang(a->m, v1);
  phys[7] = v1[0];  phys[8] = v1[1];  phys[9] = v1[2];
  if (org) {
	phys[0] += org[0];  phys[1] += org[1];  phys[2] += org[2]; }
  physical_to_dlt(phys, findlt, 0);
}

double photo_reliability(long n, float *xy, float *xyz, float *dlt)
/* Calculate a normalized reliability parameter for a photograph.  The
 *  normalized reliability is the recipricol of the average pixel error in
 *  the image compared to where it 'should' be.
 * Enter: long n: number of nodes.
 *        float *xy: 2D nodal coordinates.  Stored x1,y1,x2,y2,x3,...
 *        float *xyz: 3D nodal coordinates.  Stored x1,y1,z1,x2,y2,z3,x3,...
 *        float *dlt: 18 dlt/distortion parameters for this photograph.
 * Exit:  double rel: reliability of photo.  Higher is better. 2/12/96-DWM */
{
  double x, y, X, Y, Z, calcx, calcy, cx, cy, r2, dx, dy, total=0, den;
  long i, count=0;

  for (i=0; i<n; i++) {
    x = xy[i*2];   y = xy[i*2+1];
    X = xyz[i*3];  Y = xyz[i*3+1];  Z = xyz[i*3+2];
    if (x>1e29 || y>1e29 || X>1e29 || Y>1e29 || Z>1e29)  continue;
    den = dlt[8]*X + dlt[9]*Y + dlt[10]*Z + 1;
    if (den) {
      calcx = (dlt[0]*X + dlt[1]*Y + dlt[2]*Z + dlt[3])/den;
      calcy = (dlt[4]*X + dlt[5]*Y + dlt[6]*Z + dlt[7])/den;
      cx = calcx+dlt[11];  cy = calcy+dlt[12];       /** Sign could be off */
      r2 = cx*cx+cy*cy;
      dx = dlt[16]*(r2+2*cx*cx) + 2*dlt[17]*cx*cy +
           cx*(dlt[13]*r2+dlt[14]*r2*r2+dlt[15]*r2*r2*r2);
      dy = dlt[17]*(r2+2*cy*cy) + 2*dlt[16]*cx*cy +
           cy*(dlt[13]*r2+dlt[14]*r2*r2+dlt[15]*r2*r2*r2);
      calcx -= dx;  calcy -= dy;
      total += sqrt(sq(calcx-x) + sq(calcy-y));
      count++; } }
  if (count)  total /= count;
  else        total = 1e6;
  if (total)  total = 1./total;
  else        total = 1e6;
  return(total);
}

long physical_to_dlt(float *phys, float *dlt, long focal)
/* Convert a set of physical parameters to DLT parameters.
 * Enter: float *phys: The physical parameters, either in the order
 *                     X0, Y0, Z0, Cx, Cy, x0, y0, theta, phi, psi  or
 *                     X0, Y0, Z0, lx, ly, x0, y0, theta, phi, psi, f.
 *        float *dlt: location to store 18 dlt/distortion parameters.
 *                    Distortion is set to zero.
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

void rotate_dlt(float *dlt, float *org, float *axis, float angle,
				float *findlt)
/* Rotate a camera position about a point and axis by a given angle.  This
 *  takes a set of DLT parameters and produces a new set of DLT parameters.
 * Enter: float *dlt: set of 18 dlt parameters.
 *        float *org: origin of rotation.  This is in the original camera
 *                    coordinate system.  This is an array of 3 values.  This
 *                    may be null if the origin is at 0,0,0.
 *        float *axis: normalized axis of rotation.  This may be null if the
 *                     axis is 0,0,1.
 *        float angle: angle of rotation in radians.
 *        float *findlt: location to output 18 dlt parameters.  This may be
 *                       the same as the input.                 6/5/97-DWM */
{
  double phys[10], mc[55], v1[3], v2[3], x, y, z;
  lmatrix *m=lmat(mc), *r=lmat(mc+11), *a=lmat(mc+22), *ainv=lmat(mc+33);
  lmatrix *ar=lmat(mc+44);

  memcpy(findlt, dlt, 18*sizeof(float));
  dlt_to_physical(dlt, phys, m, 0);
  if (org) {
	phys[0] -= org[0];  phys[1] -= org[1];  phys[2] -= org[2]; }
  a->w = a->h = ainv->w = ainv->h = r->w = r->h = 3;
  if (axis) {
	v1[0] = v1[1] = v1[2] = 1;
	if (axis[0]==axis[1] && axis[0]==axis[2])
	  v1[0] = -1;
	lcross_product(axis, v1, v2);  lunit(v2, v2);
	lcross_product(v2, axis, v1);  lunit(v1, v1);
	a->m[0] = v1[0];  a->m[1] = v2[0];  a->m[2] = axis[0];
	a->m[3] = v1[1];  a->m[4] = v2[1];  a->m[5] = axis[1];
	a->m[6] = v1[2];  a->m[7] = v2[2];  a->m[8] = axis[2];
	lmatdup(a, r);  lmatinv(r, ainv); }
  else {
	lmatident(a, 3);  lmatident(ainv, 3); }
  r->m[0] = cos(angle);  r->m[1] = -sin(angle);  r->m[2] = 0;
  r->m[3] = sin(angle);  r->m[4] =  cos(angle);  r->m[5] = 0;
  r->m[6] = 0;           r->m[7] =  0;           r->m[8] = 1;
  lmatmul(a, r, ar);  lmatmul(ar, ainv, r);
  x = phys[0]*r->m[0]+phys[1]*r->m[3]+phys[2]*r->m[6];
  y = phys[0]*r->m[1]+phys[1]*r->m[4]+phys[2]*r->m[7];
  z = phys[0]*r->m[2]+phys[1]*r->m[5]+phys[2]*r->m[8];
  phys[0] = x;  phys[1] = y;  phys[2] = z;
  lmatmul(m, r, a);
  lmattoang(a->m, v1);
  phys[7] = v1[0];  phys[8] = v1[1];  phys[9] = v1[2];
  if (org) {
	phys[0] += org[0];  phys[1] += org[1];  phys[2] += org[2]; }
  physical_to_dlt(phys, findlt, 0);
}
#endif

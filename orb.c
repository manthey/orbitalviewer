/* This contains the orbital specific routines for the Oribital Viewer
 *  program.  These routines are not operating system or compiler dependant,
 *  and should be fully ANSI C compatible.                     6/11/97-DWM */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "orb.h"

real CutGrad[3], CutProximity, SurfaceXYZ[3], ToJump;
long CutSurf;

real calc_grad(real *r, OATOM *p, long mag)
/* Calculate the cartesian gradient at the last point that probability was
 *  calculated for a particular atom.  Note that the gradient calculated is
 *  for the psi^2 function, not the psi function.
 * Enter: real *r: pointer to array of three reals to store gradient in.
 *        OATOM *p: atom to calculate gradient for.
 *        long mag: 0 to compute directional gradient, 1 to compute magnitude
 *                 of gradient only.
 * Exit:  real magnitude: magnitude of gradient, or 0 if mag is 0.
 *                                                             12/5/93-DWM */
{
  real gr, gt, gp, sum, st=p->lts, ct=p->ltc, ct2=p->ltc2, cons, sp, cp;
  long j;

  if (!p->lr) {
    r[0] = r[1] = r[2] = 0;
    return(0); }
  cons = 2*p->lpsi*p->cc;
  if (p->nml1)  sum = p->crg[0];
  else          sum = 0;
  for (j=1; j<=p->nml1-1; j++)
    sum = sum*p->lrzna+p->crg[j];
  gr = cons*p->ltf*p->lpf*(p->lrf*(p->l/p->lr-p->zna)+p->zna*p->lrle*sum);
  if (st) {
    if (ct) {
      sum = p->ctg[0];
      for (j=1; j<=p->lmam12; j++)
        sum = sum*ct2+p->ctg[j];
      if (!(p->lmam&1))
        sum *= ct; }
    else
      sum = 0;
    if (p->m)
      gt = cons*p->lrf*p->lpf*(p->am*ct*p->ltf/st-st*p->lstm*sum)/p->lr;
    else
      gt = -cons*p->lrf*p->lpf*st*sum/p->lr;
    if (p->m>0)  gp =  cos(p->m*p->lp);
    else         gp = -sin(p->m*p->lp);
    gp = cons*p->lrf*p->ltf*gp*sqrt2*p->m/(p->lr*st); }
  else
    gt = gp = 0;
  if (mag)
    return(sqrt(gt*gt+gp*gp+gr*gr));
  cp = cos(p->lp);  sp = sin(p->lp);
  r[0] = gr*st*cp - gp*sp + gt*ct*cp;
  r[1] = gr*st*sp + gp*cp + gt*ct*sp;
  r[2] = gr*ct            - gt*st;
  if (!p->identity)  transform2(r, r, p->matp);
  return(0);
}

real calc_grad_mag(MOLECULE *mol)
/* Calculate the magnitude of the total gradient at the last point whose
 *  probability was calculated.
 * Enter: MOLECULE *mol: pointer to molecular record.
 * Exit:  real mag: magnitude of gradient.                     1/22/95-DWM */
{
  real r[3];

  if (mol->nump==1)
    return(calc_grad(r, mol->orb, 1));
  calc_grad_total(r, mol);
  return(sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
}

void calc_grad_total(real *r, MOLECULE *mol)
/* Calculate the total gradient at the last point whose probability was
 *  calculated.
 * Enter: real *r: location to store result.
 *        MOLECULE *mol: pointer to molecule to calculate.
 * Exit:  real *r: location of result.                        12/12/93-DWM */
{
  long i;
  real p[3];

  if (mol->nump==1) {
    calc_grad(r, mol->orb, 0);
    return; }
  r[0] = r[1] = r[2] = 0;
  for (i=0; i<mol->nump; i++)
    if (mol->orb[i].lpsi) {
      calc_grad(p, mol->orb+i, 0);
      r[0] += p[0]/mol->orb[i].lpsi;
      r[1] += p[1]/mol->orb[i].lpsi;
      r[2] += p[2]/mol->orb[i].lpsi; }
  r[0] *= mol->lastpsi;  r[1] *= mol->lastpsi;  r[2] *= mol->lastpsi;
}

real calc_prob(real *x, OATOM *p)
/* Calculate the probability at the select point in space.  Actually
 *  calculates psi.  This should be squared for the actual probability.
 * Enter: real *x: coordinates.  These are not modified.
 *        OATOM *p: atom to calculate probability for.         12/5/93-DWM */
{
  real psi, r, ph, tx[3], sum, ct, st, ct2, lrzna;
  long j;

  tx[0] = x[0]-p->x[0];  tx[1] = x[1]-p->x[1];  tx[2] = x[2]-p->x[2];
  if (!p->identity)  transform2(tx, tx, p->matm);
  p->lr = r = sqrt((sum=(tx[0]*tx[0]+tx[1]*tx[1]))+tx[2]*tx[2]);
  if (!r) return(0);
  p->lts = st = sqrt(sum)/r;
  p->ltc = ct = tx[2]/r;
  if (tx[0])  p->lp = ph = atan2(tx[1], tx[0]);
  else        p->lp = ph = PI2 + PI*(tx[1]<0);
  p->lrzna = lrzna = r*p->zna;
  sum = p->crs[0];
  for (j=1; j<=p->nml1; j++)
    sum = sum*lrzna+p->crs[j];
  p->lrf = (p->lrle=powr(r, p->l)*exp(-lrzna))*sum;
  if (p->lmam2)  p->ltc2 = ct2 = ct*ct;
  sum = p->cts[0];
  for (j=1; j<=p->lmam2; j++)
    sum = sum*ct2+p->cts[j];
  if (p->lmam&1)
    sum *= ct;
  if (!ct) { if (p->lmam<2)  sum = 1;  else sum = 0; }
  if (st)          p->ltf = (p->lstm=powr(st, p->am))*sum;
  else             p->ltf = (!p->m)*sum;
  if (p->m>0)      p->lpf = sqrt2*sin(p->m*ph);
  else if (!p->m)  p->lpf = 1;
  else             p->lpf = sqrt2*cos(p->m*ph);
  p->lpsi = psi = p->cc*p->lrf*p->ltf*p->lpf;
  return(psi);
}

real calc_prob_total(real *x, MOLECULE *mol, long *phase)
/* Calculate the total probability (psi^2) based on all orbitals at a given
 *  point.  The sign of the returned number is always positive.
 * Enter: real *x: location for which to calculate probability.
 *        MOLECULE *mol: location of molecular record.
 *        long *phase: location to store phase, either 1 or -1.
 * Exit:  real psi2: probability.                             12/12/93-DWM */
{
  long i;

  for (i=mol->lastpsi=0; i<mol->nump; i++)
    mol->lastpsi += calc_prob(x, mol->orb+i);
  if (mol->lastpsi>=0)
    phase[0] = 1;
  else
    phase[0] = -1;
  return(mol->lastpsi*mol->lastpsi);
}

real *calc_unit(real *v)
/* Calculate a unit vector.
 * Enter: real *v: vector to normalize.
 * Exit:  real *v: location of result.                        12/12/93-DWM */
{
  real total;

  total = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (total) {
    v[0] /= total;
    v[1] /= total;
    v[2] /= total; }
  return(v);
}

real *cross(real *v1, real *v2, real *r)
/* Compute the cross product of two vectors.  The result may replace one of
 *  the source vectors.
 * Enter: real *v1: pointer to array of real values of the first vector.
 *        real *v2: pointer to second vector.
 *        real *r: pointer to location for result.
 * Exit:  real *r: pointer to location of result (same as entered value).
 *                                                             8/20/91-DWM */
{
  real t[3];            /* Uses extra vector to allow r to equal v1 or v2 */

  t[0] = v1[1]*v2[2] - v1[2]*v2[1];
  t[1] = v1[2]*v2[0] - v1[0]*v2[2];
  t[2] = v1[0]*v2[1] - v1[1]*v2[0];
  r[0] = t[0];  r[1] = t[1];  r[2] = t[2];
  return(r);
}

long cut_intercept(real *x, real *v, CUTAWAY *cut)
/* Find the first surface intersection between the cutaway planes and a
 *  particular vector, starting at a given starting point.  The intersection
 *  is prohibited from being at the starting point.
 * Enter: real *x: starting point and location of result.
 *        real *v: vector to travel in.  Modified.
 *        CUTAWAY *cut: pointer to cutaway record, or null for none.
 * Exit:  long hit: 0 for no intersection, 1 for intersection.11/23/97-DWM */
{
  real x2[3], x3[3]={0,0,0}, *mat, inter[3], den, d, best=1e30, dist;
  long i, j, k, l, numcut, last=-1, inv;

  if (!cut)  return(0);
  calc_unit(v);
  memcpy(x2, x, 3*sizeof(real));
  cut_zone(x3, cut, -1);              /* Just to initalize cutaway matrix */
  for (i=0; cut[i].notlast; i++);  numcut = i+1;
  for (k=0; k<numcut; k++)  if (k!=last) {
    mat = cut[k].mat;  inv = 1-2*(cut[k].invert!=0);
    for (i=2; i>2-cut[k].type; i--) {
      den = v[0]*mat[i*3]+v[1]*mat[i*3+1]+v[2]*mat[i*3+2];
      if (!den)  continue;
      d = (mat[i*3]*(cut[k].xyz[0]-x2[0])+mat[i*3+1]*(cut[k].xyz[1]-x2[1])+
           mat[i*3+2]*(cut[k].xyz[2]-x2[2]))/den;
      if (d<0)  continue;
      inter[0] = x2[0]+v[0]*d;
      inter[1] = x2[1]+v[1]*d;
      inter[2] = x2[2]+v[2]*d;
      if (cut[k].type>=2) {
        for (j=0; j<3; j++)
          x3[j] = mat[j*3]*(inter[0]-cut[k].xyz[0])+mat[j*3+1]*(inter[1]-
                  cut[k].xyz[1])+mat[j*3+2]*(inter[2]-cut[k].xyz[2]);
        if (cut[k].type==2) {
          if (i==2 && x3[1]<0)  continue;
          if (i==1 && x3[2]<0)  continue; }
        if (cut[k].type==3) {
          if (i==2 && (x3[0]<0 || x3[1]<0))  continue;
          if (i==1 && (x3[0]<0 || x3[2]<0))  continue;
          if (i==0 && (x3[1]<0 || x3[2]<0))  continue; } }
      dist = sq(inter[0]-x2[0])+sq(inter[1]-x2[1])+sq(inter[2]-x2[2]);
      if (dist<best) {
        best = dist;
        CutGrad[0] = -mat[i]  *inv;
        CutGrad[1] = -mat[i+3]*inv;
        CutGrad[2] = -mat[i+6]*inv;
        CutSurf = i;
        if (!inter[0])  inter[0] = 1e-100;
        if (!inter[1])  inter[1] = 1e-100;
        if (!inter[2])  inter[2] = 1e-100;
        memcpy(x, inter, 3*sizeof(real)); } }
    if (best!=1e30)
      for (l=0; l<numcut; l++)  if (l!=k)
        if (cut_zone(x, cut, 0)) {
          best = 1e30;  last = k;  k = -1; } }
  return(best!=1e30);
}

long cut_zone(real *x, CUTAWAY *cut, long zero)
/* Check if a point is within the cutaway zone.  If the point is within the
 *  cutaway zone, CutProximity contains the distance to the cutaway.  Note
 *  that for multiple cutaway zones, this distance is the distance to the
 *  first zone that the point was within.
 * Enter: real *x: three coordinates of the point in meters.
 *        CUTAWAY *cut: cutaway zone.  Null for none.
 *        long zero: 0-actual cross point is cut, 1-cross point is kept, -1
 *                   to just initialize matrix.
 * Exit:  long cut: 0 for outside of the cut zone (point is kept), 1 for
 *                  within it (point is cut).                   7/6/97-DWM */
{
  real u1[3], u2[3], u3[3], x2[3], x3[3], *mat;
  long first=1;

  if (!cut)  return(0);
  do {
    if (first)  first = 0;
    else        cut = cut+1;
    if (!cut->type)  continue;
    mat = cut->mat;
    if (!mat[0] && !mat[1] && !mat[2]) {
      u1[1] = u1[2] = u2[0] = u2[2] = u3[0] = u3[1] = 0;
      u1[0] = u2[1] = u3[2] = 1;
      transform(u1, u1, cut->ang, 1);
      transform(u2, u2, cut->ang, 1);
      transform(u3, u3, cut->ang, 1);
      mat[0]=u1[0];  mat[3]=u2[0];  mat[6]=u3[0];
      mat[1]=u1[1];  mat[4]=u2[1];  mat[7]=u3[1];
      mat[2]=u1[2];  mat[5]=u2[2];  mat[8]=u3[2]; }
    if (zero==-1)  continue;
    x2[0] = x[0]-cut->xyz[0];
    x2[1] = x[1]-cut->xyz[1];
    x2[2] = x[2]-cut->xyz[2];
    x3[0] = x2[0]*mat[0]+x2[1]*mat[1]+x2[2]*mat[2];
    x3[1] = x2[0]*mat[3]+x2[1]*mat[4]+x2[2]*mat[5];
    x3[2] = x2[0]*mat[6]+x2[1]*mat[7]+x2[2]*mat[8];
    CutProximity = x3[2];
    if (cut->type>=2 && x3[1]<CutProximity)  CutProximity = x3[1];
    if (cut->type>=3 && x3[0]<CutProximity)  CutProximity = x3[0];
    if (!cut->invert) {
      if (!zero) { if (CutProximity<0)          continue; }
      else       { if (CutProximity<a0*0.0002)  continue; } }
    else {
      if (!zero) { if (CutProximity>0)          continue; }
      else       { if (CutProximity>-a0*0.0002) continue; } }
    return(1); }
  while (cut->notlast);
  return(0);
}

long intercept(real *x, real *v, MOLECULE *mol, CUTAWAY *cut)
/* Find the first surface intersection with a particular vector, starting at
 *  a given starting point.
 * Enter: real *x: starting point and location of result.
 *        real *v: vector to travel in.  Modified.
 *        MOLECULE *mol: pointer to molecular record.
 *        CUTAWAY *cut: pointer to cutaway record, or null for none.
 * Exit:  long found: 0 for none, 1 for intersection with positive phase, -1
 *                    for intersection with negative phase.  The phase number
 *                    is double (2 or -2) if the intersection occurs on the
 *                    shear plane of a cutaway.               12/12/93-DWM */
{
  long phase;
  real mul=0, psi1=0, psi2, total, v2[3], x2[3], psi3, mul2;

  if (!to_sphere(x, v, mol->maxcheck, mol->zsize, mol->syscntr)) {
    mol->lasthit = 0;
    return(0); }
  total = to_jump()/mol->zsize;
  while (1) {
    psi2 = psi1;
    if ((psi1=calc_prob_total(x, mol, &phase))>mol->Psi) {
      while (mul>mol->minres) {
        mul2 = (psi1-mol->Psi)/(psi1-psi2)*mul;
        if (mul2<mol->minres)           mul2 = mol->minres;
        else if (mul-mul2<mol->minres)  mul2 = mul-mol->minres;
        x2[0]=x[0]-v[0]*mul2;  x2[1]=x[1]-v[1]*mul2;  x2[2]=x[2]-v[2]*mul2;
        psi3 = calc_prob_total(x2, mol, &phase);
        if (psi3<mol->Psi) {
          psi2 = psi3;  mul = mul2; }
        else if (psi3>mol->Psi) {
          psi1 = psi3;  mul -= mul2;  total -= mul2;
          x[0] = x2[0];  x[1] = x2[1];  x[2] = x2[2]; }
        else {
          psi1 = psi3;  mul = 0;  total -= mul2;
          x[0] = x2[0];  x[1] = x2[1];  x[2] = x2[2]; } }
      if (cut)  if (cut_zone(x, cut, 0)) {
        memcpy(x2, x, 3*sizeof(real));  memcpy(v2, v, 3*sizeof(real));
        if (!cut_intercept(x, v2, cut)) {
          mol->lasthit = 0;
          return(0); }
        total += sqrt(sq(x[0]-x2[0])+sq(x[1]-x2[1])+sq(x[2]-x2[2]))/
                 mol->zsize;
        if ((psi1=calc_prob_total(x, mol, &phase))<=mol->Psi)
          continue;
        mol->lasthit = total;
        return(phase*2); }
      mol->lasthit = total;
      return(phase); }
    if (mol->nump==1) { if (mol->orb[0].lr>=mol->maxcheck)  break; }
    else if (sq(x[0]-mol->syscntr[0])+sq(x[1]-mol->syscntr[1])+sq(x[2]-
             mol->syscntr[2])>=sq(mol->maxcheck))  break;
    mul = calc_grad_mag(mol)*mol->zsize;
    if (mul) {
      if (psi1>psi2)
        mul = (mol->Psi-psi1)/mul;
      else
        mul = (mol->Psi+psi1)/mul; }
    if (mul<mol->minres+total*mol->back)  mul = mol->minres+total*mol->back;
    else if (mul>mol->maxjump)  mul = mol->maxjump;
    total += mul;
    x[0] += v[0]*mul;  x[1] += v[1]*mul;  x[2] += v[2]*mul; }
  mol->lasthit = 0;
  return(0);
}

void matrix_inverse(real *mat, real *inv)
/* Compute the inverse of a 3x3 matrix.  If the determinant of the original
 *  matrix is zero, the transpose is stored instead.  The destination can not
 *  be the same as the original matrix.
 * Enter: real *mat: matrix in such a form that [x y z][0 1 2] = [x' y' z'].
 *        real *inv: location to store inverse.        [3 4 5]
 *                                                      [6 7 8]1/11/94-DWM */
{
  real det=mat[0]*mat[4]*mat[8]+mat[3]*mat[7]*mat[2]+mat[6]*mat[1]*mat[5]-
          (mat[0]*mat[7]*mat[5]+mat[3]*mat[1]*mat[8]+mat[6]*mat[4]*mat[2]);

  if (det) {
    inv[0] = (mat[4]*mat[8]-mat[5]*mat[7])/det;
    inv[3] = (mat[5]*mat[6]-mat[3]*mat[8])/det;
    inv[6] = (mat[3]*mat[7]-mat[4]*mat[6])/det;
    inv[1] = (mat[2]*mat[7]-mat[1]*mat[8])/det;
    inv[4] = (mat[0]*mat[8]-mat[2]*mat[6])/det;
    inv[7] = (mat[1]*mat[6]-mat[0]*mat[7])/det;
    inv[2] = (mat[1]*mat[5]-mat[2]*mat[4])/det;
    inv[5] = (mat[2]*mat[3]-mat[0]*mat[5])/det;
    inv[8] = (mat[0]*mat[4]-mat[1]*mat[3])/det; }
  else {
    inv[0] = mat[0];  inv[3] = mat[1];  inv[6] = mat[2];  inv[1] = mat[3];
    inv[4] = mat[4];  inv[7] = mat[5];  inv[2] = mat[6];  inv[5] = mat[7];
    inv[8] = mat[8]; }
}

void oangtomat(double *r, double ta, double pa, double sa)
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

long orb_asymptote(ASYMPTOTE *as, CUTAWAY *cut, MOLECULE *mol, float time,
                   long process)
/* Calculate polygons for an asymptote plot.  This takes into account cutaway
 *  information, and can be run as a partial process.  If this is the first
 *  time this routine is called, the ASYMPTOTE structure should be zeroed
 *  except for opacity and density.  The maximum number of triangles is
 *  dependant on the orbital and the density.  Also, this should be called
 *  with process==4 before the program exits to free allocated memory.
 * Internally, the asymptote structure has a process flag.  This tells what
 *  stage the calculation is as follows: 0-nothing done, 1-finished prepping
 *  atoms, 2-finished prepping check, 3-finished computing fast lookup
 *  tables, 4-processing, 5-finished.
 * Enter: ASYMPTOTE *as: structure to store results.  Also specifies opacity
 *                       and density.
 *        CUTAWAY *cut: cutaway to use.  Note that if mat is all zero, it is
 *                      calculated from ang, otherwise, ang is ignored.  Null
 *                      for no cutaway (same as type==0).
 *        MOLECULE *mol: molecule to calculate.  This structure may be
 *                       modified in the initialization process.
 *        float time: amount of time in seconds that this routine is allowed
 *                    to use.  0 requests that the routine runs until
 *                    finished.  The time request is merely a suggestion; the
 *                    actual time used may differ significantly.
 *        long process: 0 to continue processing.
 *                      1 for start completely new molecule.
 *                      2 for change cutaway.
 *                      4 to free all allocated resources.
 * Exit:  long done: 0 to signal that more processing is needed, 1 to
 *                   indicate that processing is finished, -1 for memory
 *                   error.                                    11/9/97-DWM */
{
  clock_t cl=clock();
  long i, j, k, d, phase, cas, *tempelem, a, b, pos[4];
  long *tetra;
  long tetra1[]={0,1,2,4, 1,4,5,7, 1,3,2,7, 2,7,6,4, 2,7,1,4};
  long tetra2[]={3,5,6,7, 0,3,6,2, 0,4,6,5, 0,1,3,5, 0,3,5,6};
  long poly[]={-1,0,0,0,0,0,0,0,  0,3,1,3,2,3,-1,0,  0,2,1,2,3,2,-1,0,
               0,2,2,1,1,3,3,0,   0,1,2,1,3,1,-1,0,  1,0,0,3,3,2,2,1,
               1,0,0,2,2,3,3,1,   1,0,2,0,3,0,-1,0};
  OATOM *p;
  real r, ph, th, ct, st, sum, di;
  real xyz[24], psi2[8], psi;
  real pa, pb, pc, dist, xa[3], xb[3], xc[3];
  float *tempxyz;

  if (!process && as->process==5)  return(1);
  if (as->process>5 && as->opacity && !as->xyz)  as->process = 0;
  if (!as->opacity || as->process>5) { as->process = 5;  return(1); }
  if (cut)  if (!cut->type)  cut = 0;
  switch (process) {
    case 1: as->process = 0; break;
    case 2: as->process = 3; break;
    case 4: if (as->xyz) { free2(as->xyz); as->xyz = 0; }
      if (as->elem)   { free2(as->elem);   as->elem = 0; }
      if (as->scrxyz) { free2(as->scrxyz); as->scrxyz = 0; }
      if (as->norm)   { free2(as->norm);   as->norm = 0; }
      if (as->enorm)  { free2(as->enorm);  as->enorm = 0; }
      if (as->rfac)   { free2(as->rfac);   as->rfac = 0; }
      if (as->tfac)   { free2(as->tfac);   as->tfac = 0; }
      if (as->pfac)   { free2(as->pfac);   as->pfac = 0; }
      as->n = as->e = as->process = 0;  return(1); }
  if (as->process>=2)  if (as->checkpsi!=mol->checkpsi) {
    psi = mol->Psi;  mol->Psi = as->checkpsi;
    prep_check(mol);
    mol->Psi = psi; }
  switch (as->process) {
    case 0: for (i=0; i<mol->nump; i++)
        prep_atom(mol->orb+i, mol->orb[i].factor);
      as->process = 1;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 1: psi = mol->Psi;
      mol->Psi = 0.01;
      prep_check(mol);
      while (mol->maxpsi<1e-30) {
        mol->Psi /= 10;
        prep_check(mol); }
      mol->Psi = mol->maxpsi/1000;
      prep_check(mol);
      as->checkpsi = mol->checkpsi;
      mol->Psi = psi;
      as->process = 2;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 2: if (as->rfac) { free2(as->rfac);  as->rfac = 0; }
      if (as->tfac) { free2(as->tfac);  as->tfac = 0; }
      if (as->pfac) { free2(as->pfac);  as->pfac = 0; }
      if (as->density<6)  as->density = 6;
      if (as->density&1)  as->density += 1;
      if (mol->nump==1 && !mol->orb[0].x[0] && !mol->orb[0].x[1] &&
          !mol->orb[0].x[2] && !mol->orb[0].ang[0] && !mol->orb[0].ang[1] &&
          !mol->orb[0].ang[2]) {
        as->rfac = malloc2((size_t)(as->density+1)*sizeof(real));
        as->tfac = malloc2((size_t)(as->density+1)*sizeof(real));
        as->pfac = malloc2((size_t)(as->density*2+1)*sizeof(real)); }
      if ((!as->tfac || !as->pfac) && as->rfac) {
        free2(as->rfac);  as->rfac = 0; }
      if (as->rfac) {
        p = mol->orb;
        for (i=0; i<=as->density; i++) {
          if (i) {
            r = mol->maxcheck*i/as->density;
            p->lrzna = r*p->zna;
            for (j=sum=0; j<=p->nml1; j++)
              sum += p->crs[j]*powr(p->lrzna, p->nml1-j);
            as->rfac[i] = p->cc*powr(r, p->l)*exp(-p->lrzna)*sum; }
          else
            as->rfac[i] = 0;
          th = PI*i/as->density;
          if (i==as->density)  th = PI;
          st = sin(th);  ct = cos(th);
          for (j=sum=0; j<=p->lmam/2 && ct; j++)
            sum += p->cts[j]*powr(ct, p->lmam-j-j);
          if (!ct && p->lmam<2)  sum = 1;
          if (st)  as->tfac[i] = powr(st, p->am)*sum;
          else     as->tfac[i] = (!p->m)*sum; }
        for (i=0; i<=as->density*2; i++) {
          ph = 2*PI*i/(as->density*2);
          if (i==as->density*2)  ph = 0;
          if (p->m>0)      as->pfac[i] = sqrt2*sin(p->m*ph);
          else if (!p->m)  as->pfac[i] = 1;
          else             as->pfac[i] = sqrt2*cos(p->m*ph); } }
      as->process = 3;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 3: if (!as->xyz)  as->xyz = malloc2(1);
      if (!as->elem)  as->elem = malloc2(1);
      if (!as->xyz || !as->elem)  return(-1);
      as->n = as->e = 0;
      as->index[0] = 1;  as->index[1] = as->index[2] = 0;
      as->process = 4;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 4: for (; as->index[2]<as->density*2; as->index[1]=0,as->index[2]++)
      for (; as->index[1]<as->density; as->index[0]=0, as->index[1]++)
      for (; as->index[0]<as->density; as->index[0]++) {
        if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
        for (k=d=0; k<=1; k++)  for (j=0; j<=1; j++)
          for (i=0; i<=1; i++, d++) {
            r = mol->maxcheck*(as->index[0]+i)/as->density;
            th = PI*(as->index[1]+j)/as->density;
            ph = 2*PI*(as->index[2]+k)/(as->density*2);
            if (as->index[1]+j==as->density)   th = PI;
            if (as->index[2]+k==as->density*2) ph = 0;
            xyz[d*3]   = r*cos(ph)*sin(th)+mol->syscntr[0];
            xyz[d*3+1] = r*sin(ph)*sin(th)+mol->syscntr[1];
            xyz[d*3+2] = r*        cos(th)+mol->syscntr[2];
            if (as->rfac) {
              psi = as->rfac[as->index[0]+i]*as->tfac[as->index[1]+j]*
                    as->pfac[as->index[2]+k];
              psi2[d] = psi*fabs(psi); }
            else {
              psi2[d] = calc_prob_total(xyz+d*3, mol, &phase);
              if (phase<0)  psi2[d] = -psi2[d]; } }
        for (i=0; i<8; i++) {
          if (cut_zone(xyz+i*3, cut, 0))  break;
          if (!psi2[i])  psi2[i] = 1e-30; }
        if (i!=8)  continue;
        tetra = tetra1;
        if ((as->index[0]+as->index[1]+as->index[2])&1)
          tetra = tetra2;
        for (i=0; i<5; i++) {
          cas = (psi2[tetra[i*4  ]]>0)*8+(psi2[tetra[i*4+1]]>0)*4+
                (psi2[tetra[i*4+2]]>0)*2+(psi2[tetra[i*4+3]]>0);
          if (cas>=8)  cas = 15-cas;
          if (!cas)  continue;
          if (!(tempxyz=realloc2(as->xyz,(size_t)(as->n+4)*3*sizeof(float))))
            return(-1);
          as->xyz = tempxyz;
          if (!(tempelem=realloc2(as->elem,(size_t)(as->e+2)*3*sizeof(long))))
            return(-1);
          as->elem = tempelem;
          if (as->norm)   { free2(as->norm);   as->norm = 0; }
          if (as->enorm)  { free2(as->enorm);  as->enorm = 0; }
          for (j=0; j<4; j++) {
            if (poly[cas*8+j*2]<0) break;
            a = tetra[i*4+poly[cas*8+j*2]]; b = tetra[i*4+poly[cas*8+j*2+1]];
            if (psi2[a]<psi2[b]) {
              pa = psi2[a];  pb = psi2[b];
              memcpy(xa, xyz+a*3, 3*sizeof(real));
              memcpy(xb, xyz+b*3, 3*sizeof(real)); }
            else {
              pa = psi2[b];  pb = psi2[a];
              memcpy(xa, xyz+b*3, 3*sizeof(real));
              memcpy(xb, xyz+a*3, 3*sizeof(real)); }
            dist = sqrt(sq(xa[0]-xb[0])+sq(xa[1]-xb[1])+sq(xa[2]-xb[2]));
            do {
              xc[0] = (xa[0]+xb[0])*0.5;
              xc[1] = (xa[1]+xb[1])*0.5;
              xc[2] = (xa[2]+xb[2])*0.5;
              pc = calc_prob_total(xc, mol, &phase);
              pc *= phase;
              if (pc>0) { pb = pc;  memcpy(xb, xc, 3*sizeof(real)); }
              else {      pa = pc;  memcpy(xa, xc, 3*sizeof(real)); }
              dist /= 2; }
            while (dist>mol->minres*mol->zsize);
            for (k=0; k<3; k++)  as->xyz[as->n*3+k] = xc[k];
            pos[j] = as->n;  as->n++;
            cut_zone(xc, cut, 0);
            for (k=0; k<as->n-1; k++)
              if ((di=fabs(sq(as->xyz[k*3]-as->xyz[pos[j]*3])+
                           sq(as->xyz[k*3+1]-as->xyz[pos[j]*3+1])+
                           sq(as->xyz[k*3+2]-as->xyz[pos[j]*3+2])))<
                  sq(mol->maxcheck/as->density/4)) {
                if (cut)
                  if (CutProximity<mol->maxcheck*0.001 &&
                      sqrt(di)>mol->minres*mol->zsize*2)
                    continue;
                pos[j] = k;  as->n--;  break; }
            if (j==2)
              if (pos[0]!=pos[1] && pos[1]!=pos[2] && pos[2]!=pos[0]) {
                as->elem[as->e*3] = pos[1];
                as->elem[as->e*3+1] = pos[2];
                as->elem[as->e*3+2] = pos[0];
                if ((as->xyz[pos[0]*3]-as->xyz[pos[1]*3])*(as->xyz[pos[2]*3
                    +1]-as->xyz[pos[1]*3+1])-(as->xyz[pos[2]*3]-as->xyz[
                    pos[1]*3])*(as->xyz[pos[0]*3+1]-as->xyz[pos[1]*3+1])>0) {
                  as->elem[as->e*3+1] = pos[0];
                  as->elem[as->e*3+2] = pos[2]; }
                as->e++; }
            if (j==3)
              if (pos[0]!=pos[2] && pos[2]!=pos[3] && pos[3]!=pos[0]) {
                as->elem[as->e*3] = pos[3];
                as->elem[as->e*3+1] = pos[0];
                as->elem[as->e*3+2] = pos[2];
                if ((as->xyz[pos[2]*3]-as->xyz[pos[3]*3])*(as->xyz[pos[0]*3
                    +1]-as->xyz[pos[3]*3+1])-(as->xyz[pos[0]*3]-as->xyz[
                    pos[3]*3])*(as->xyz[pos[2]*3+1]-as->xyz[pos[3]*3+1])>0) {
                  as->elem[as->e*3+1] = pos[2];
                  as->elem[as->e*3+2] = pos[0]; }
                as->e++; } } } }
      as->process = 5; }
  return(1);
}

long orb_points(OPOINTS *pt, CUTAWAY *cut, MOLECULE *mol, float time,
                long process)
/* Calculate points for a dot plot.  This takes into account cutaway
 *  information, and can be run as a partial process.  If this is the first
 *  time this routine is called, the OPOINTS structure should be zeroed
 *  except for maxn, which should be set to the maximum number of points to
 *  calculate.  Also, this should be called with process==4 before the
 *  program exits to free allocated memory.
 * Internally, the point structure has a process flag.  This tells what stage
 *  the calculation is as follows: 0-nothing done, 1-finished prepping atoms,
 *  2-finished prepping check, 3-finished clearing old points, 4-finished
 *  allocating space and now processing, 5-finished.
 * Enter: OPOINTS *pt: structure to store results.  Also specifies maxn.
 *        CUTAWAY *cut: cutaway to use.  Note that if mat is all zero, it is
 *                      calculated from ang, otherwise, ang is ignored.  Null
 *                      for no cutaway (same as type==0).
 *        MOLECULE *mol: molecule to calculate.  This structure may be
 *                       modified in the initialization process.
 *        float time: amount of time in seconds that this routine is allowed
 *                    to use.  0 requests that the routine runs until
 *                    finished.  The time request is merely a suggestion; the
 *                    actual time used may differ significantly.
 *        long process: 0 to continue processing.
 *                      1 for start completely new molecule.
 *                      2 for change cutaway.
 *                      3 for change maxn.
 *                      4 to free all allocated resources.
 * Exit:  long done: 0 to signal that more processing is needed, 1 to
 *                   indicate that processing is finished, -1 for memory
 *                   error.                                    11/9/97-DWM */
{
  clock_t cl=clock();
  long i;
  float *tempxyz;
  char *tempphase;
  real x[3], val, psi;

  if (!process && pt->process==5)  return(1);
  if (pt->process>5) { pt->process = 5;  return(1); }
  if (cut)  if (!cut->type)  cut = 0;
  switch (process) {
    case 1: pt->process = 0; break;
    case 2: pt->process = 2; break;
    case 3: pt->process = 3; break;
    case 4: if (pt->xyz) { free2(pt->xyz);  pt->xyz = 0; }
      if (pt->scrxyz) { free2(pt->scrxyz);  pt->scrxyz = 0; }
      if (pt->phase)  { free2(pt->phase);   pt->phase = 0; }
      pt->n = pt->process = 0;  return(1); }
  if (pt->process>=2)  if (pt->checkpsi!=mol->checkpsi) {
    psi = mol->Psi;  mol->Psi = pt->checkpsi;
    prep_check(mol);
    mol->Psi = psi; }
  switch (pt->process) {
    case 0: for (i=0; i<mol->nump; i++)
        prep_atom(mol->orb+i, mol->orb[i].factor);
      pt->process = 1;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 1: psi = mol->Psi;
      mol->Psi = 0.01;
      prep_check(mol);
      while (mol->maxpsi<1e-30) {
        mol->Psi /= 10;
        prep_check(mol); }
      mol->Psi = mol->maxpsi/1000;
      prep_check(mol);
      pt->checkpsi = mol->checkpsi;
      mol->Psi = psi;
      pt->process = 2;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 2: if (pt->n)  pt->n = 0;
      pt->process = 3;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 3: if (pt->n>pt->maxn) {
        pt->n = pt->maxn; }
      if (!pt->xyz || !pt->phase)  pt->n = 0;
      if (!pt->xyz)  tempxyz = malloc2((size_t)pt->maxn*3*sizeof(float));
      else           tempxyz = realloc2(pt->xyz,(size_t)pt->maxn*3*sizeof(float));
      if (!tempxyz)  return(-1);
      pt->xyz = tempxyz;
      if (!pt->phase)  tempphase=malloc2((size_t)pt->maxn*sizeof(char));
      else             tempphase=realloc2(pt->phase,(size_t)pt->maxn*sizeof(char));
      if (!tempphase)  return(-1);
      pt->phase = tempphase;
      pt->process = 4;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 4: while (pt->n<pt->maxn) {
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
      x[0] = rnd(-mol->maxcheck, mol->maxcheck, 3);
      x[1] = rnd(-mol->maxcheck, mol->maxcheck, 3);
      x[2] = rnd(-mol->maxcheck, mol->maxcheck, 3);
      if (sq(x[0])+sq(x[1])+sq(x[2])>mol->maxcheck2)  continue;
      x[0]+=mol->syscntr[0];  x[1]+=mol->syscntr[1];  x[2]+=mol->syscntr[2];
      if (cut_zone(x, cut, 0))  continue;
      val = calc_prob_total(x, mol, &i);
      if (val<rnd(0, mol->maxpsi, 3))  continue;
      pt->xyz[pt->n*3]   = x[0];
      pt->xyz[pt->n*3+1] = x[1];
      pt->xyz[pt->n*3+2] = x[2];
      pt->phase[pt->n] = i;
      pt->n++; } }
  pt->process = 5;
  return(1);
}

long orb_polygons(POLYGON *pl, CUTAWAY *cut, MOLECULE *mol, float time,
                  long process)
/* Calculate polygons for an orbital plot.  This takes into account cutaway
 *  information, and can be run as a partial process.  If this is the first
 *  time this routine is called, the POLYGON structure should be zeroed
 *  except for density.  The maximum number of triangles is dependant on the
 *  orbital and the density.  Also, this should be called with process==4
 *  before the program exits to free allocated memory.
 * Internally, the polygon structure has a process flag.  This tells what
 *  stage the calculation is as follows: 0-nothing done, 1-finished prepping
 *  atoms, 2-finished prepping check, 3-finished computing fast lookup tables,
 *  4-processing, 5-refining triangles, 6-finished.
 * If process==5, then the processmeter variable contains a number indicating
 *  how refined the triangles are.  This is a value which ranges from around
 *  sq(mol->maxcheck/pl->density) to around
 *  sq(mol->maxcheck/pl->density/pl->refine).
 * Enter: POLYGON *pl: structure to store results.  Also specifies density.
 *        CUTAWAY *cut: cutaway to use.  Note that if mat is all zero, it is
 *                      calculated from ang, otherwise, ang is ignored.  Null
 *                      for no cutaway (same as type==0).
 *        MOLECULE *mol: molecule to calculate.  This structure may be
 *                       modified in the initialization process.
 *        float time: amount of time in seconds that this routine is allowed
 *                    to use.  0 requests that the routine runs until
 *                    finished.  The time request is merely a suggestion; the
 *                    actual time used may differ significantly.
 *        long process: 0 to continue processing.
 *                      1 for start completely new molecule.
 *                      2 for change cutaway.
 *                      4 to free all allocated resources.
 *                      5 for change psi^2.
 *                      6 for refine triangles.
 * Exit:  long done: 0 to signal that more processing is needed, 1 to
 *                   indicate that processing is finished, -1 for memory
 *                   error.                                    11/9/97-DWM */
{
  clock_t cl=clock();
  long i, j, k, l, d, phase, cas, *tempelem, a, b, pos[4], e1, e2;
  long *tetra, *tempptph;
  long tetra1[]={0,1,2,4, 1,4,5,7, 1,3,2,7, 2,7,6,4, 2,7,1,4};
  long tetra2[]={3,5,6,7, 0,3,6,2, 0,4,6,5, 0,1,3,5, 0,3,5,6};
  long poly[]={-1,0,0,0,0,0,0,0,  0,3,1,3,2,3,-1,0,  0,2,1,2,3,2,-1,0,
               0,2,2,1,1,3,3,0,   0,1,2,1,3,1,-1,0,  1,0,0,3,3,2,2,1,
               1,0,0,2,2,3,3,1,   1,0,2,0,3,0,-1,0};
  OATOM *p;
  real r, ph, th, ct, st, sum, best, di;
  real xyz[24], psi2[8], psi;
  real pa, pb, pc, dist, xa[3], xb[3], xc[3];
  float *tempxyz;
  char *tempph;

  if (!process && pl->process==6)  return(1);
  if (pl->process>6) { pl->process = 6;  return(1); }
  if (cut)  if (!cut->type)  cut = 0;
  switch (process) {
    case 1: pl->process = 0; break;
    case 2: pl->process = 3; break;
    case 4: if (pl->xyz) { free2(pl->xyz);   pl->xyz = 0; }
      if (pl->elem)    { free2(pl->elem);    pl->elem = 0; }
      if (pl->scrxyz)  { free2(pl->scrxyz);  pl->scrxyz = 0; }
      if (pl->norm)    { free2(pl->norm);    pl->norm = 0; }
      if (pl->enorm)   { free2(pl->enorm);   pl->enorm = 0; }
      if (pl->phase)   { free2(pl->phase);   pl->phase = 0; }
      if (pl->ptphase) { free2(pl->ptphase); pl->ptphase = 0; }
      if (pl->rfac)    { free2(pl->rfac);    pl->rfac = 0; }
      if (pl->tfac)    { free2(pl->tfac);    pl->tfac = 0; }
      if (pl->pfac)    { free2(pl->pfac);    pl->pfac = 0; }
      pl->n = pl->e = pl->process = 0;  return(1);
    case 5: pl->process = 1; break;
    case 6: pl->process = 5; }
  if (pl->process>=2)  if (pl->checkpsi!=mol->checkpsi) {
    psi = mol->Psi;  mol->Psi = pl->checkpsi;
    prep_check(mol);
    if (!mol->maxcheck)  mol->maxcheck = a0;
    mol->Psi = psi; }
  switch (pl->process) {
    case 0: for (i=0; i<mol->nump; i++)
        prep_atom(mol->orb+i, mol->orb[i].factor);
      pl->process = 1;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 1: prep_check(mol);
      pl->checkpsi = mol->checkpsi;
      if (!mol->maxcheck)  mol->maxcheck = a0;
      pl->process = 2;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 2: if (pl->rfac) { free2(pl->rfac);  pl->rfac = 0; }
      if (pl->tfac) { free2(pl->tfac);  pl->tfac = 0; }
      if (pl->pfac) { free2(pl->pfac);  pl->pfac = 0; }
      if (pl->density<6)  pl->density = 6;
      if (pl->density&1)  pl->density += 1;
      if (mol->nump==1 && !mol->orb[0].x[0] && !mol->orb[0].x[1] &&
          !mol->orb[0].x[2] && !mol->orb[0].ang[0] && !mol->orb[0].ang[1] &&
          !mol->orb[0].ang[2]) {
        pl->rfac = malloc2((size_t)(pl->density+1)*sizeof(real));
        pl->tfac = malloc2((size_t)(pl->density+1)*sizeof(real));
        pl->pfac = malloc2((size_t)(pl->density*2+1)*sizeof(real)); }
      if ((!pl->tfac || !pl->pfac) && pl->rfac) {
        free2(pl->rfac);  pl->rfac = 0; }
      if (pl->rfac) {
        p = mol->orb;
        for (i=0; i<=pl->density; i++) {
          if (i) {
            r = mol->maxcheck*i/pl->density;
            p->lrzna = r*p->zna;
            for (j=sum=0; j<=p->nml1; j++)
              sum += p->crs[j]*powr(p->lrzna, p->nml1-j);
            pl->rfac[i] = p->cc*powr(r, p->l)*exp(-p->lrzna)*sum; }
          else
            pl->rfac[i] = 0;
          th = PI*i/pl->density;
          if (i==pl->density)  th = PI;
          st = sin(th);  ct = cos(th);
          for (j=sum=0; j<=p->lmam/2 && ct; j++)
            sum += p->cts[j]*powr(ct, p->lmam-j-j);
          if (!ct && p->lmam<2)  sum = 1;
          if (st)  pl->tfac[i] = powr(st, p->am)*sum;
          else     pl->tfac[i] = (!p->m)*sum; }
        for (i=0; i<=pl->density*2; i++) {
          ph = 2*PI*i/(pl->density*2);
          if (i==pl->density*2)  ph = 0;
          if (p->m>0)      pl->pfac[i] = sqrt2*sin(p->m*ph);
          else if (!p->m)  pl->pfac[i] = 1;
          else             pl->pfac[i] = sqrt2*cos(p->m*ph); } }
      pl->process = 3;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 3: if (!pl->xyz)  pl->xyz = malloc2(1);
      if (!pl->elem)  pl->elem = malloc2(1);
      if (!pl->phase)  pl->phase = malloc2(1);
      if (!pl->ptphase)  pl->ptphase = malloc2(1);
      if (!pl->xyz || !pl->elem || !pl->phase || !pl->ptphase)  return(-1);
      pl->n = pl->e = 0;
      pl->index[0] = 1;  pl->index[1] = pl->index[2] = 0;
      pl->process = 4;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 4: for (; pl->index[2]<pl->density*2; pl->index[1]=0, pl->index[2]++)
      for (; pl->index[1]<pl->density; pl->index[0]=0, pl->index[1]++)
      for (; pl->index[0]<pl->density; pl->index[0]++) {
        if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
        for (k=d=0; k<=1; k++)  for (j=0; j<=1; j++)
          for (i=0; i<=1; i++, d++) {
            r = mol->maxcheck*(pl->index[0]+i)/pl->density;
            th = PI*(pl->index[1]+j)/pl->density;
            ph = 2*PI*(pl->index[2]+k)/(pl->density*2);
            if (pl->index[1]+j==pl->density)   th = PI;
            if (pl->index[2]+k==pl->density*2) ph = 0;
            xyz[d*3]   = r*cos(ph)*sin(th)+mol->syscntr[0];
            xyz[d*3+1] = r*sin(ph)*sin(th)+mol->syscntr[1];
            xyz[d*3+2] = r*        cos(th)+mol->syscntr[2];
            if (pl->rfac) {
              psi = pl->rfac[pl->index[0]+i]*pl->tfac[pl->index[1]+j]*
                    pl->pfac[pl->index[2]+k];
              psi2[d] = psi*fabs(psi); }
            else {
              psi2[d] = calc_prob_total(xyz+d*3, mol, &phase);
              if (phase<0)  psi2[d] = -psi2[d]; } }
        for (i=0; i<8; i++) {
          if (cut_zone(xyz+i*3, cut, 1))  break;
          if (!psi2[i])  psi2[i] = 1e-30; }
        if (i!=8)  continue;
        tetra = tetra1;
        if ((pl->index[0]+pl->index[1]+pl->index[2])&1)
          tetra = tetra2;
        for (i=0; i<5; i++)
          for (l=-1; l<=1; l+=2) {
            cas = (psi2[tetra[i*4  ]]>mol->Psi*l)*8+(psi2[tetra[i*4+1]]>
                   mol->Psi*l)*4+(psi2[tetra[i*4+2]]>mol->Psi*l)*2+
                   (psi2[tetra[i*4+3]]>mol->Psi*l);
            if (cas>=8)  cas = 15-cas;
            if (!cas)  continue;
            if (!(tempxyz=realloc2(pl->xyz,(size_t)(pl->n+4)*3*
                sizeof(float))))
              return(-1);
            pl->xyz = tempxyz;
            if (!(tempelem=realloc2(pl->elem,(size_t)(pl->e+2)*3*
                sizeof(long))))
              return(-1);
            pl->elem = tempelem;
            if (pl->norm)   { free2(pl->norm);   pl->norm = 0; }
            if (pl->enorm)  { free2(pl->enorm);  pl->enorm = 0; }
            if (!(tempph=realloc2(pl->phase,(size_t)pl->e+2)))
              return(-1);
            pl->phase = tempph;
            if (!(tempptph=realloc2(pl->ptphase,(size_t)(pl->n+4)*
                sizeof(long))))
              return(-1);
            pl->ptphase = tempptph;
            for (j=0; j<4; j++) {
              if (poly[cas*8+j*2]<0) break;
              a = tetra[i*4+poly[cas*8+j*2]];
              b = tetra[i*4+poly[cas*8+j*2+1]];
              if (psi2[a]<psi2[b]) {
                pa = psi2[a];  pb = psi2[b];
                memcpy(xa, xyz+a*3, 3*sizeof(real));
                memcpy(xb, xyz+b*3, 3*sizeof(real)); }
              else {
                pa = psi2[b];  pb = psi2[a];
                memcpy(xa, xyz+b*3, 3*sizeof(real));
                memcpy(xb, xyz+a*3, 3*sizeof(real)); }
              dist = sqrt(sq(xa[0]-xb[0])+sq(xa[1]-xb[1])+sq(xa[2]-xb[2]));
              do {
                xc[0] = (xa[0]+xb[0])*0.5;
                xc[1] = (xa[1]+xb[1])*0.5;
                xc[2] = (xa[2]+xb[2])*0.5;
                pc = calc_prob_total(xc, mol, &phase);
                pc *= phase;
                if (pc>mol->Psi*l) {
                  pb = pc;  memcpy(xb, xc, 3*sizeof(real)); }
                else {
                  pa = pc;  memcpy(xa, xc, 3*sizeof(real)); }
                dist /= 2; }
              while (dist>mol->minres*mol->zsize);
              for (k=0; k<3; k++)  pl->xyz[pl->n*3+k] = xc[k];
              pl->ptphase[pl->n] = l;
              pos[j] = pl->n;  pl->n++;
              cut_zone(xc, cut, 2);
              for (k=0; k<pl->n-1; k++)  if (pl->ptphase[k]*l>0)
                if ((di=fabs(sq(pl->xyz[k*3]-pl->xyz[pos[j]*3])+
                             sq(pl->xyz[k*3+1]-pl->xyz[pos[j]*3+1])+
                             sq(pl->xyz[k*3+2]-pl->xyz[pos[j]*3+2])))<
                    sq(mol->maxcheck/pl->density/2)) {
                  if (cut)
                    if (CutProximity<mol->maxcheck*0.001 &&
                        sqrt(di)>mol->minres*mol->zsize*2)
                      continue;
                  pos[j] = k;  pl->n--;  break; }
              if (j==2)
                if (pos[0]!=pos[1] && pos[1]!=pos[2] && pos[2]!=pos[0]) {
                  pl->elem[pl->e*3] = pos[1];
                  pl->elem[pl->e*3+1] = pos[2];
                  pl->elem[pl->e*3+2] = pos[0];
                  if (orb_polygons_normal(pl->elem+pl->e*3, pl->xyz, mol)) {
                    pl->elem[pl->e*3+1] = pos[0];
                    pl->elem[pl->e*3+2] = pos[2]; }
                  pl->phase[pl->e] = l;  pl->e++; }
              if (j==3)
                if (pos[0]!=pos[2] && pos[2]!=pos[3] && pos[3]!=pos[0]) {
                  pl->elem[pl->e*3] = pos[3];
                  pl->elem[pl->e*3+1] = pos[0];
                  pl->elem[pl->e*3+2] = pos[2];
                  if (orb_polygons_normal(pl->elem+pl->e*3, pl->xyz, mol)) {
                    pl->elem[pl->e*3+1] = pos[2];
                    pl->elem[pl->e*3+2] = pos[0]; }
                  pl->phase[pl->e] = l;  pl->e++; } } } }
      pl->process = 5;  pl->processmeter = 1e30;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 5: if (pl->refine) while (1) {
        if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
        for (i=0, best=0, e1=0; i<pl->e; i++)
          for (j=0; j<3; j++) {
            dist=sq(pl->xyz[pl->elem[i*3+j]*3  ]-pl->xyz[pl->elem[i*3+((j+1)%3)]*3  ])+
                 sq(pl->xyz[pl->elem[i*3+j]*3+1]-pl->xyz[pl->elem[i*3+((j+1)%3)]*3+1])+
                 sq(pl->xyz[pl->elem[i*3+j]*3+2]-pl->xyz[pl->elem[i*3+((j+1)%3)]*3+2]);
            if (dist>best) {
              best = dist;  e1 = i*3+j;  e2 = -1; }
            else if (dist==best)
              if (pl->elem[i*3+((j+1)%3)]==pl->elem[e1] &&
                  pl->elem[i*3+j]==pl->elem[(e1/3)*3+(e1+1)%3])
                e2 = i*3+j; }
        if (best>pl->processmeter)  break;
        pl->processmeter = best;
        if (best<sq(mol->maxcheck/pl->density/pl->refine))  break;
        if (!(tempxyz=realloc2(pl->xyz,(size_t)(pl->n+1)*3*sizeof(float))))
          return(-1);
        pl->xyz = tempxyz;
        if (!(tempelem=realloc2(pl->elem,(size_t)(pl->e+2)*3*sizeof(long))))
          return(-1);
        pl->elem = tempelem;
        if (pl->norm)   { free2(pl->norm);   pl->norm = 0; }
        if (pl->enorm)  { free2(pl->enorm);  pl->enorm = 0; }
        if (!(tempph=realloc2(pl->phase,(size_t)pl->e+2)))
          return(-1);
        pl->phase = tempph;
        if (!(tempptph=realloc2(pl->ptphase,(size_t)(pl->n+1)*sizeof(long))))
          return(-1);
        pl->ptphase = tempptph;
        pl->ptphase[pl->n] = pl->ptphase[pl->elem[e1]];
        for (i=0; i<3; i++)
          pl->xyz[pl->n*3+i] = (pl->xyz[pl->elem[e1]*3+i]+
                               pl->xyz[pl->elem[(e1/3)*3+(e1+1)%3]*3+i])*0.5;
        pl->n++;
        snap_to_surface(pl->xyz+(pl->n-1)*3, mol, mol->Psi*pl->phase[e1/3],
                        sqrt(best)/2.2);
        pl->elem[pl->e*3] = pl->elem[e1];
        pl->elem[pl->e*3+1] = pl->n-1;
        pl->elem[pl->e*3+2] = pl->elem[(e1/3)*3+(e1+2)%3];
        pl->phase[pl->e] = pl->phase[e1/3];
        pl->e++;
        pl->elem[e1] = pl->n-1;
        if (e2>=0) {
          pl->elem[pl->e*3] = pl->elem[e2];
          pl->elem[pl->e*3+1] = pl->n-1;
          pl->elem[pl->e*3+2] = pl->elem[(e2/3)*3+(e2+2)%3];
          pl->phase[pl->e] = pl->phase[e2/3];
          pl->e++;
          pl->elem[e2] = pl->n-1; } }
      pl->process = 6; }
  return(1);
}

long orb_polygons_normal(long *e, float *xyz, MOLECULE *mol)
/* Determine if the element is in the correct orientation based on its
 *  normal and the local gradient.  Note that the element is guaranteed to
 *  have non-zero length sides.
 * Enter: long *e: array of three node numbers.
 *        float *xyz: nodal coordinates.
 *        MOLECULE *mol: molecule to check.
 * Exit:  long reversed: 0 for the element is correct, 1 for the element is
 *                       reversed.                            11/14/97-DWM */
{
  real norm[3], x[3], grad[3];
  long phase;

  norm[0] = (xyz[e[0]*3+1]-xyz[e[1]*3+1])*(xyz[e[2]*3+2]-xyz[e[1]*3+2])-
            (xyz[e[0]*3+2]-xyz[e[1]*3+2])*(xyz[e[2]*3+1]-xyz[e[1]*3+1]);
  norm[1] = (xyz[e[0]*3+2]-xyz[e[1]*3+2])*(xyz[e[2]*3  ]-xyz[e[1]*3  ])-
            (xyz[e[0]*3  ]-xyz[e[1]*3  ])*(xyz[e[2]*3+2]-xyz[e[1]*3+2]);
  norm[2] = (xyz[e[0]*3  ]-xyz[e[1]*3  ])*(xyz[e[2]*3+1]-xyz[e[1]*3+1])-
            (xyz[e[0]*3+1]-xyz[e[1]*3+1])*(xyz[e[2]*3  ]-xyz[e[1]*3  ]);
  x[0] = xyz[e[1]*3];  x[1] = xyz[e[1]*3+1];  x[2] = xyz[e[1]*3+2];
  calc_prob_total(x, mol, &phase);
  calc_grad_total(grad, mol);
  return(norm[0]*grad[0]+norm[1]*grad[1]+norm[2]*grad[2]<0);
}

long orb_render(RENDER *re, STEREO *st, CUTAWAY *cut, MOLECULE *mol,
                float time, long process)
/* Produce a rendered image.  This takes into account cutaway information,
 *  asymptotes, opacity, and refraction, and can be run as a partial process.
 *  If this is the first time this routine is called, the RENDER structure
 *  should be zeroed except for type, antialias, color, and camera.  Note
 *  that the angles used with camera are in the "standard" math system, not
 *  in the "standard" physics system used in the rest of this code.  Also,
 *  this should be called with process==4 before the program exits to free
 *  allocated memory.  Note that the routine is vastly faster for an iso-
 *  surface with no asymptotes.
 * Internally, the render structure has a process flag.  This tells what
 *  stage the calculation is as follows: 0-nothing done, 1-finished prepping
 *  atoms, 2-finished prepping check, 3-finished allocating image space,
 *  4-finished computing 1/16th image, 5-finished computing 1/4 image,
 *  6-finished non-antialiased image, 7-finished.
 * Enter: RENDER *re: structure to store results.  Also specifies type,
 *                    antialias, color, and camera.
 *        STEREO *st: stereo mode to render in.  Null for monoscopic (same as
 *                    mode==0);
 *        CUTAWAY *cut: cutaway to use.  Note that if mat is all zero, it is
 *                      calculated from ang, otherwise, ang is ignored.  Null
 *                      for no cutaway (same as type==0).
 *        MOLECULE *mol: molecule to calculate.  This structure may be
 *                       modified in the initialization process.
 *        float time: amount of time in seconds that this routine is allowed
 *                    to use.  0 requests that the routine runs until
 *                    finished.  The time request is merely a suggestion; the
 *                    actual time used may differ significantly.
 *        long process: 0 to continue processing.
 *                      1 for start completely new molecule.
 *                      2 for change cutaway.
 *                      4 to free all allocated resources.
 *                      6 for start completely new molecule, but stop after
 *                        prep_check.  This allows the camera to be adjusted
 *                        for size prior to drawing.  Set process to zero to
 *                        resume.
 * Exit:  long done: 0 to signal that more processing is needed, 1 to
 *                   indicate that processing is finished, -1 for memory
 *                   error.                                   11/16/97-DWM */
{
  clock_t cl=clock();
  long i, j, k, stop=0, coarse=0, scan, dir, t, precise=0, full, old=0;
  uchar *img;
  real clrin[12], clrout[3], total[3], oldpsi, psi;

  if (!process && re->process==7)  return(1);
  if (cut)  if (!cut->type)  cut = 0;
  if (st)   if (!st->mode)   st  = 0;
  if (re->opacity[0] || (re->opacity[1]!=1 && re->opacity[1]) ||
      (re->opacity[2] && re->opacity[1]!=1) || re->opacity[3] ||
      re->opacity[4]!=re->opacity[1] || (re->opacity[5] && re->opacity[4]!=1)
      || re->opacity[6])
    precise = 1;
  if (cut)  if (cut->nosurface)  precise = 1;
  if (re->refrac[0]!=re->refrac[1] || re->refrac[0]!=re->refrac[2] ||
      (re->refrac[0]!=0 && re->refrac[0]!=1))
    precise = 1;
  if (cut)  for (i=0; i<mol->numl; i++)
    if (mol->ls[i].a!=1)  precise = 1;
  for (i=0; i<12; i++)
    clrin[i] = (real)re->color[i]/255;
  switch (process) {
    case 6: stop = 1;
    case 1: re->process = 0; break;
    case 2: re->process = 2; break;
    case 4: if (re->buf) { free2(re->buf);  re->buf = 0; }
      if (re->zbuf) { free2(re->zbuf);  re->zbuf = 0; }
      if (re->phase) { free2(re->phase);  re->phase = 0; }
      re->process = 0;  return(1); }
  if (re->process>=2)  if (re->checkpsi!=mol->checkpsi) {
    psi = mol->Psi;  mol->Psi = re->checkpsi;
    prep_check(mol);
    mol->Psi = psi; }
  if ((re->lastw!=re->w || re->lasth!=re->h) && re->process>2)
    re->process = 2;
  switch (re->process) {
    case 0: for (i=0; i<mol->nump; i++)
        prep_atom(mol->orb+i, mol->orb[i].factor);
      re->process = 1;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 1: oldpsi = mol->Psi;
      if (precise) {
        mol->Psi = 0.01;
        prep_check(mol);
        while (mol->maxpsi<1e-30) {
          mol->Psi /= 10;
          prep_check(mol); }
        mol->Psi = mol->maxpsi/1000; }
      mol->Psi = min(mol->Psi, oldpsi);
      prep_check(mol);
      re->checkpsi = mol->checkpsi;
      mol->Psi = oldpsi;
      re->process = 2;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
      if (stop)  return(1);
    case 2: if (re->buf) {
        if (re->lastw==re->w && re->lasth==re->h && re->lasttype==re->type)
          old = 1;
        else {
          free2(re->buf);  re->buf = 0; } }
      if (re->zbuf) {
        free2(re->zbuf);  re->zbuf = 0; }
      if (re->type) {
        re->scan = (re->w*3+3)&0x7FFFFFFC;
        if (!old) {
          if (!(re->buf=malloc2(re->scan*re->h+40)))  return(-1);
          memset(re->buf, 0, 40);
          ((long *)re->buf)[0] = 40;
          ((long *)re->buf)[1] = re->w;  ((long *)re->buf)[2] = re->h;
          ((long *)re->buf)[3] = 24*0x10000+1;
          ((long *)re->buf)[5] = re->scan*re->h; }
        re->image = re->buf+40+re->scan*(re->h-1);
        re->scan *= -1; }
      else {
        re->scan = re->w*3;
        if (!old)
          if (!(re->buf=malloc2(re->scan*re->h)))  return(-1);
        re->image = re->buf; }
      if (!re->phase || !old) {
        if (re->phase)  free2(re->phase);
        if (!(re->phase=malloc2(re->h*re->w*sizeof(long))))  return(-1); }
      re->x = re->y = 0;
      if (!old)
        for (j=0; j<re->h; j++)
          for (i=0; i<re->w; i++) {
            re->image[re->scan*j+i*3]   = re->color[6+re->type];
            re->image[re->scan*j+i*3+1] = re->color[7];
            re->image[re->scan*j+i*3+2] = re->color[8-re->type]; }
      if (!re->brightness)  re->brightness = 1;
      if (re->autobright) {
        re->autobright = 1;  re->brightness = 1;
        re->antialias &= 0xFFFFFFF3; }
      oangtomat(re->cammat, re->camera[7], re->camera[8], re->camera[9]);
      for (i=0; i<mol->numl; i++)
        if (!mol->ls[i].local) {
          mol->ls[i].lx[0] = -mol->ls[i].x[0]*re->cammat[0]+mol->ls[i].x[1]*
                             re->cammat[3]-mol->ls[i].x[2]*re->cammat[6];
          mol->ls[i].lx[1] = -mol->ls[i].x[0]*re->cammat[1]+mol->ls[i].x[1]*
                             re->cammat[4]-mol->ls[i].x[2]*re->cammat[7];
          mol->ls[i].lx[2] = -mol->ls[i].x[0]*re->cammat[2]+mol->ls[i].x[1]*
                             re->cammat[5]-mol->ls[i].x[2]*re->cammat[8]; }
        else
          memcpy(mol->ls[i].lx, mol->ls[i].x, 3*sizeof(real));
      re->lastw = re->w;  re->lasth = re->h;  re->lasttype = re->type;
      if (st)  if (st->mode==StereoSTEREOGRAM)
        re->zbuf = malloc2(re->w*re->h*sizeof(ushort));
      re->process = 3;
      if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
    case 3: scan = re->scan;  img = re->image;
      if (((re->antialias>>2)&3)>0)  re->process = 4;
      else do {
        for (; re->y<re->h; re->x=0, re->y+=16)
          for (; re->x<re->w; re->x+=16) {
            if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
            re->phase[re->x+re->y*re->w] = orb_render_ray(re, st, cut, mol,
                                       re->x, re->y, precise, clrin, clrout);
            if (re->brightness==1) {
              img[re->y*scan+re->x*3]  =clrout[re->type]  *255;
              img[re->y*scan+re->x*3+1]=clrout[1]         *255;
              img[re->y*scan+re->x*3+2]=clrout[2-re->type]*255; }
            else {
              img[re->y*scan+re->x*3]  =min(clrout[re->type]  *255*re->brightness,255);
              img[re->y*scan+re->x*3+1]=min(clrout[1]         *255*re->brightness,255);
              img[re->y*scan+re->x*3+2]=min(clrout[2-re->type]*255*re->brightness,255); }
            for (j=re->y; j<re->h && j<re->y+16; j++)
              for (i=re->x; i<re->w && i<re->x+16; i++)
                if (i!=re->x || j!=re->y)
                  for (k=0; k<3; k++)
                    img[j*scan+i*3+k] = img[re->y*scan+re->x*3+k]; }
        if (re->autobright==1) {
          re->autobright = 2;  re->brightness = 0;
          for (re->y=0; re->y<re->h; re->y+=16)
            for (re->x=0; re->x<re->w; re->x+=16)
              re->brightness=max(max(re->brightness,img[re->y*scan+re->x*3]),
                  max(img[re->y*scan+re->x*3+1], img[re->y*scan+re->x*3+2]));
          re->brightness = max(max(re->brightness,clrin[6]*255),
                               max(clrin[7]*255, clrin[8]*255));
          if (!re->brightness)  re->brightness = 224;
          re->brightness = 224./re->brightness; }
        else if (re->autobright==2)  re->autobright = 3;
        re->x = re->y = 0; }
      while (re->autobright==2);
      re->process = 4;
    case 4: scan = re->scan;  img = re->image;
      full = (((re->antialias>>2)&3)>1);
      if (full && !(re->antialias&1))  re->process = 5;
      else for (; re->y<re->h; re->x=0, re->y+=4)
        for (; re->x<re->w; re->x+=4)
          if (re->x%16 || re->y%16 || ((re->antialias>>2)&3)>0) {
            if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
            re->phase[re->x+re->y*re->w] = orb_render_ray(re, st, cut, mol,
                                       re->x, re->y, precise, clrin, clrout);
            if (re->brightness==1) {
              img[re->y*scan+re->x*3]  =clrout[re->type]  *255;
              img[re->y*scan+re->x*3+1]=clrout[1]         *255;
              img[re->y*scan+re->x*3+2]=clrout[2-re->type]*255; }
            else {
              img[re->y*scan+re->x*3]  =min(clrout[re->type]  *255*re->brightness,255);
              img[re->y*scan+re->x*3+1]=min(clrout[1]         *255*re->brightness,255);
              img[re->y*scan+re->x*3+2]=min(clrout[2-re->type]*255*re->brightness,255); }
            if (!full)
              for (j=re->y; j<re->h && j<re->y+4; j++)
                for (i=re->x; i<re->w && i<re->x+4; i++)
                  if (i!=re->x || j!=re->y)
                    for (k=0; k<3; k++)
                      img[j*scan+i*3+k] = img[re->y*scan+re->x*3+k]; }
      if (full && re->process==4)  re->antialias -= (1<<2);
      re->x = re->y = 0;
      re->process = 5;
    case 5: scan = re->scan;  img = re->image;
      full = (((re->antialias>>2)&3)>1);
      for (; re->y<re->h; re->x=0, re->y++)
        for (; re->x<re->w; re->x++)
          if (re->x%4 || re->y%4 || full) {
            if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
            coarse = 0;
            if ((re->antialias&1) && ((re->antialias>>2)&3)<=1) {
              coarse = 1;
              if (re->w-re->x<=4 || re->h-re->y<=4)  coarse = 0; }
            if (st)  if (st->mode==StereoSTEREOGRAM)  coarse = 0;
            if (coarse) {
              i = re->x&0xFFFFFFFC;  j = re->y&0xFFFFFFFC;
              if (re->phase[j*re->w+i]!=re->phase[(j+4)*re->w+i] ||
                  re->phase[j*re->w+i]!=re->phase[(j+4)*re->w+i+4] ||
                  re->phase[j*re->w+i]!=re->phase[j*re->w+i+4])
                coarse = 0; }
            for (k=0; coarse && k<3; k++) {
              if (abs(img[j*scan+i*3+k]-img[j*scan+(i+4)*3+k])>15 ||
                  abs(img[j*scan+i*3+k]-img[(j+4)*scan+(i+4)*3+k])>15 ||
                  abs(img[j*scan+i*3+k]-img[(j+4)*scan+i*3+k])>15)
                coarse = 0;
              if (abs(img[j*scan+(i+4)*3+k]-img[(j+4)*scan+i*3+k])>15 ||
                  abs(img[j*scan+(i+4)*3+k]-img[(j+4)*scan+(i+4)*3+k])>15||
                  abs(img[(j+4)*scan+i*3+k]-img[(j+4)*scan+(i+4)*3+k])>15)
                coarse = 0; }
            if (coarse)  if (st)
              if (st->mode==StereoINTERLACED)
                coarse = 0;
            if (coarse) {
              for (k=0; k<3; k++)
                img[re->y*scan+re->x*3+k] =
                    (img[ j   *scan+ i   *3+k]*(float)(4-(re->x-i))/4 +
                     img[ j   *scan+(i+4)*3+k]*(float)   (re->x-i) /4)*
                                                     (float)(4-(re->y-j))/4 +
                    (img[(j+4)*scan+ i   *3+k]*(float)(4-(re->x-i))/4 +
                     img[(j+4)*scan+(i+4)*3+k]*(float)   (re->x-i) /4)*
                                                     (float)   (re->y-j) /4;
              re->phase[re->x+re->y*re->w] = re->phase[i+j*re->w];
              continue; }
            re->phase[re->x+re->y*re->w] = orb_render_ray(re, st, cut, mol,
                                       re->x, re->y, precise, clrin, clrout);
            if (re->brightness==1) {
              img[re->y*scan+re->x*3]  =clrout[re->type]  *255;
              img[re->y*scan+re->x*3+1]=clrout[1]         *255;
              img[re->y*scan+re->x*3+2]=clrout[2-re->type]*255; }
            else {
              img[re->y*scan+re->x*3]  =min(clrout[re->type]  *255*re->brightness,255);
              img[re->y*scan+re->x*3+1]=min(clrout[1]         *255*re->brightness,255);
              img[re->y*scan+re->x*3+2]=min(clrout[2-re->type]*255*re->brightness,255); } }
      re->x = re->y = 0;
      re->process = 6;
      if (!(re->antialias&2)) {
        re->process = 7;  return(1); }
    case 6: scan = re->scan;  img = re->image;
      #define ANTI_THRESHOLD 2
      #define ANTI_STEPS 4
      for (; re->y<re->h-1; re->x=0, re->y++)
        for (; re->x<re->w-1; re->x++) {
          for (k=dir=0; k<3; k++) {
            if (abs(img[re->y*scan+re->x*3+k]-img[re->y*scan+(re->x+1)*3+k])>=
                ANTI_THRESHOLD)  dir |= 1;
            if (abs(img[re->y*scan+re->x*3+k]-img[(re->y+1)*scan+re->x*3+k])>=
                ANTI_THRESHOLD)  dir |= 2; }
          if (!dir)  continue;
          if (time && clock()-cl>time*CLOCKS_PER_SEC)  return(0);
          t = 1;
          for (k=0; k<3; k++)
            total[k] = img[re->y*scan+re->x*3+k];
          for (j=0; j<1+(ANTI_STEPS-1)*(dir>>1); j++)
            for (i=0; i<1+(ANTI_STEPS-1)*(dir&1); i++)  if (i || j) {
              orb_render_ray(re, st, cut, mol, re->x+i*1.0/ANTI_STEPS,
                             re->y+j*1.0/ANTI_STEPS, precise, clrin, clrout);
              if (re->brightness==1) {
                total[0] += clrout[re->type]  *255;
                total[1] += clrout[1]         *255;
                total[2] += clrout[2-re->type]*255; }
              else {
                total[0] += min(clrout[re->type]  *255*re->brightness,255);
                total[1] += min(clrout[1]         *255*re->brightness,255);
                total[2] += min(clrout[2-re->type]*255*re->brightness,255); }
              t++; }
          for (k=0; k<3; k++)
            img[re->y*scan+re->x*3+k] = total[k]/t; }
      re->process = 7; }
  return(1);
}

long orb_render_ray(RENDER *re, STEREO *st, CUTAWAY *cut, MOLECULE *mol,
                   real xx, real yy, long precise, real *clrin, real *clrout)
/* Compute the color along a ray for a specified pixel.  This takes into
 *  account all stereo modes.
 * Enter: RENDER *re: structure to store results.  Also specifies type,
 *                    antialias, color, and camera.
 *        STEREO *st: stereo mode to render in.  Null for monoscopic (same as
 *                    mode==0);
 *        CUTAWAY *cut: cutaway to use.  Note that if mat is all zero, it is
 *                      calculated from ang, otherwise, ang is ignored.  Null
 *                      for no cutaway (same as type==0).
 *        MOLECULE *mol: molecule to calculate.  This structure may be
 *                       modified in the initialization process.
 *        real xx, yy: pixel location to calculate, in screen coordinates.
 *        long precise: 0 to use surface raytracing, 1 to use step-through
 *                      raytracing.
 *        real *clrin: array of four RGB real triplets containing the colors
 *                     of the positive phase, the negative phase, the
 *                     background, and the asymptote.  1.0==full intensity.
 *        real *clrout: location to store RGB color.  Note that the values
 *                      can potentially exceed the [0,1] range.
 * Exit:  long phase: this is either the return value of ray or ray_precise.
 *                                                              1/3/98-DWM */
{
  real x[3], v[3], c2[3], c3[3];
  long camera=0, second=0, p1, p2, i;
  double z, H, V;

  if (st) switch (st->mode) {
    case StereoSTEREOSCOPE: if (xx<re->w/2)  camera = 1;
                            else             camera = 2; break;
    case StereoINTERLACED:
      if ((((long)(yy))^st->flags)&1) camera = 1;
      else                            camera = 2; break;
    case StereoREDBLUE: case StereoOVERLAY: second = 1;  camera = 1; break;
    default: ; }
  switch (camera) {
    case 1:
      orb_render_vector(st->lrcamera, st->lrcammat, xx,yy, mol->zsize, x, v);
      break;
    case 2:
      orb_render_vector(st->lrcamera+10, st->lrcammat+9, xx, yy, mol->zsize,
                        x, v);  break;
    default:
      orb_render_vector(re->camera, re->cammat, xx, yy, mol->zsize, x, v); }
  if (!precise)
    p1 = ray(x, v, clrin, mol, cut, clrout);
  else
    p1 = ray_precise(x, v, clrin, re, mol, cut, clrout);
  if (st) {
    if (second) {
      orb_render_vector(st->lrcamera+10, st->lrcammat+9, xx, yy, mol->zsize,
                        x, v);
      if (!precise)  p2 = ray(x, v, clrin, mol, cut, c2);
      else           p2 = ray_precise(x, v, clrin, re, mol, cut, c2);
      switch (st->mode) {
        case StereoREDBLUE:
          c3[0] = max(max(clrout[0], clrout[1]), clrout[2]);
          c3[1] = max(max(c2[0], c2[1]), c2[2]);
          if (!p1 && !p2) {
            clrout[0] = clrin[6];
            clrout[1] = clrin[7];
            clrout[2] = clrin[8]; break; }
          clrout[0] = clrout[1] = clrout[2] = 0;
          if (p1)  clrout[st->flags&2] = c3[0];
          if (p2)  clrout[2-(st->flags&2)] = c3[1]; break;
        case StereoOVERLAY: if (p2) {
            if (!p1) {
              clrout[0] = c2[0]*0.5;
              clrout[1] = c2[1]*0.5;
              clrout[2] = c2[2]*0.5;  p1 = p2; }
            else {
              clrout[0] = (clrout[0]+c2[0])*0.5;
              clrout[1] = (clrout[1]+c2[1])*0.5;
              clrout[2] = (clrout[2]+c2[2])*0.5; } }
          else
            if (p1) {
              clrout[0] *= 0.5;
              clrout[1] *= 0.5;
              clrout[2] *= 0.5; } break;
        default: ; } }
    if (st->mode==StereoSTEREOGRAM && re->zbuf) {
      if (p1) {
        z = re->cammat[6]*SurfaceXYZ[0]+re->cammat[7]*SurfaceXYZ[1]+
            re->cammat[8]*SurfaceXYZ[2];
        z = (z+mol->maxcheck)/(mol->maxcheck*2)*65535;
        if (z<0)  z = 0;  if (z>65535)  z = 65535; }
      else
        z = 65535;
      re->zbuf[((long)yy)*re->w+((long)xx)] = z; }
    if (st->mode==StereoCHROMADEPTH && p1) {
      z = re->cammat[6]*SurfaceXYZ[0]+re->cammat[7]*SurfaceXYZ[1]+
          re->cammat[8]*SurfaceXYZ[2];
      z = (z+mol->maxcheck)/(mol->maxcheck*2);
      if (z<0)  z = 0;  if (z>1)  z = 1;
      H = z*5;
      i = floor(H);  H = 1-(H-i);
      V = max(max(clrout[0], clrout[1]), clrout[2]);
      switch (i) {
        case 0: clrout[0]=V;       clrout[1]=V*(1-H); clrout[2]=0; break;
        case 1: clrout[0]=V*H;     clrout[1]=V;       clrout[2]=0; break;
        case 2: clrout[0]=0;       clrout[1]=V;       clrout[2]=V*(1-H); break;
        case 3: clrout[0]=0;       clrout[1]=V*H;     clrout[2]=V; break;
        case 4: clrout[0]=V*(1-H); clrout[1]=0;       clrout[2]=V; break;
        case 5: clrout[0]=V;       clrout[1]=0;       clrout[2]=V*H; } } }
  return(p1);
}

void orb_render_vector(double *c, double *cammat, real x, real y, real zsize,
                       real *xyz, real *v)
/* Determine the camera location and the vector for a specific pixel used for
 *  rendering.
 * Enter: double *c: pointer to 10 camera elements.  These are identical to
 *                   the physical parameters converted from a set of 11 dlt
 *                   values.
 *        double *cammat: array of 9 camera matrix values.  This is the
 *                        matrix produced by angtomat(cam+7).
 *        real x, y: pixel coordinate to determine a vector for.
 *        real zsize: scaling factor for vector.  This is usually mol->zsize.
 *        real *xyz: location to store the camera location.
 *        real *v: location to store the vector.  This points from the camera
 *                 location forward through the pixel.  It is of length
 *                 zsize.                                     11/24/97-DWM */
{
  real ct=cos(c[7]), st=sin(c[7]), cp=cos(c[8]), sp=sin(c[8]);
  real cs=cos(c[9]), ss=sin(c[9]), den;

  den = (c[4]*cp*(c[3]*cs+(c[5]-x)*ss)+c[3]*(y-c[6])*sp);
  if (!den)  den = 1;
  v[0] = (cp*(c[3]*c[4]*c[0]*cs+c[4]*c[0]*(c[5]-x)*ss+c[3]*(c[2]-1)*(y-c[6])*
         st)-sp*(c[3]*c[4]*(c[2]-1)*cs*st+c[4]*(1-c[2])*(x-c[5])*ss*st-c[3]*
         c[0]*(y-c[6]))+c[4]*(c[2]-1)*ct*((x-c[5])*cs+c[3]*ss))/den;
  v[1] = (cp*(c[3]*c[4]*c[1]*cs+c[4]*c[1]*(c[5]-x)*ss+c[3]*(c[2]-1)*(y-c[6])*
         ct)-sp*(c[3]*c[4]*(c[2]-1)*cs*ct+c[4]*(1-c[2])*(x-c[5])*ss*ct-c[3]*
         c[1]*(y-c[6]))+c[4]*(1-c[2])*st*((x-c[5])*cs+c[3]*ss))/den;
  v[2] = zsize/sqrt(v[0]*v[0]+v[1]*v[1]+1);
  if (v[0]*cammat[6]+v[1]*cammat[7]+cammat[8]<0)
    v[2] *= -1;
  v[0] *= v[2];  v[1] *= v[2];
  xyz[0] = c[0];
  xyz[1] = c[1];
  xyz[2] = c[2];
}

real pow1(long minusone, long pow)
/* Raise -1 to an integer power.
 * Enter: long minusone: dummy to allow compatibility with pow().
 *        long pow: power to raise -1 to.
 * Exit:  real result: -1^pow.                                 7/20/97-DWM */
{
  if (pow&1)
    return(-1);
  return(1);
}

real powr(real val, long pow)
/* Raise a number to a non-negative integer power.
 * Enter: real val: number to raise to a power.
 *        long pow: power.  n^0 = 1.  0^0 = 0.
 * Exit:  real valpow: val^pow.                                7/20/97-DWM */
{
  #ifdef ANSIC
  real ret;

  if (!pow)  return(val!=0);
  ret = 1;
  while (pow) {
    if (pow&1)  ret *= val;
    val *= val;  pow >>= 1; }
  return(ret);
  #else
  _asm {
          mov eax, pow
          sub eax, 01
          je     ipower6
          jc     ipower3
          fld val
          fld st(0)
ipower1:  shr eax, 1
          jnc    ipower2
          fmul st(1), st
          je     ipower4
ipower2:  fmul st, st(0)
          jmp    ipower1
ipower3:  fld1
          jmp    ipower5
ipower4:  fstp st(0)
ipower5:  fstp val
ipower6:
    }
  return(val);
  #endif
}

void prep_atom(OATOM *p, real factor)
/* Set up all the constant used for swift calculation involved with a
 *  particular atom.
 * Enter: OATOM *p: atom record.  Must contain n,l,m,Z,mass,x,y,z,theta,phi.
 *                  Mass may be N+Z.  x,y,z may be in meters or a0.  theta
 *                  and phi may be in degrees, but only if >2pi or <2pi;
 *                  otherwise they must be in radians.
 *        real factor: multiplicative factor for orbitals.  Use +1/-1 for
 *                     standard/negative phases.               12/5/93-DWM */
{
  long i, j;
  real a, u1[3], u2[3], u3[3];
  static long ever=0;
  static real fac[MAXN*2]={1}, snfac;

  if (!ever) {
    for (i=1; i<MAXN*2; i++)
      fac[i] = i*fac[i-1];
    snfac = sqrt(nfac);
    ever = 1; }
  if (p->n<1)                p->n = 1;
  if (p->l>=p->n || p->l<0)  p->l = 0;
  if (abs(p->m)>p->l)        p->m = 0;
  if (p->Z<1)                p->Z = 1;
  if (!p->mass)              p->mass = p->Z+p->N;
  if (p->mass<ma)            p->mass = ma;
  if (p->mass>0.5) {                             /* must be in terms of au */
    if (p->mass=((long)p->mass))
      p->mass = p->Z*mz + (p->mass-p->Z)*mn;
    else
      p->mass *= ma; }
  if (fabs(p->x[0])>1e-4)  p->x[0] *= a0;
  if (fabs(p->x[1])>1e-4)  p->x[1] *= a0;
  if (fabs(p->x[2])>1e-4)  p->x[2] *= a0;
  if (fabs(p->ang[0])>2*PI)  p->ang[0] *= deg;
  if (fabs(p->ang[1])>2*PI)  p->ang[1] *= deg;
  if (fabs(p->ang[2])>2*PI)  p->ang[2] *= deg;
  p->mu = me*p->mass/(me+p->mass);
  a = a0*me/p->mu;
  p->nml1 = p->n - p->l - 1;
  p->am = abs(p->m);
  p->lmam = p->l - p->am;  p->lmam2 = p->lmam/2;  p->lmam12 = (p->lmam-1)/2;
  p->cc = factor * pow1(-1, (p->m+p->am+2)/2) * powr(p->Z, p->l+1) /
                                   (powr(p->n, p->l+2)*powr(a, p->l+1)) *
          sqrt(p->Z*fac[p->nml1]*fac[p->n+p->l]*fac[p->lmam]*(2*p->l+1) /
                (PI*a*fac[p->l+p->am])) * snfac;
  p->zna = p->Z/(p->n*a);
  for (j=0; j<=p->nml1; j++) {
    p->crs[j] = powr(-2, p->nml1-j)/(fac[p->n+p->l-j]*fac[j]*fac[p->nml1-j]);
    p->crg[j] = p->crs[j]*(p->nml1-j); }
  for (j=0; j<=p->lmam/2; j++) {
    p->cts[j] = pow1(-1, j)*fac[2*(p->l-j)]/
                                       (fac[j]*fac[p->l-j]*fac[p->lmam-2*j]);
    p->ctg[j] = p->cts[j]*(p->lmam-2*j); }
  p->maxr = 1.e7*a0*p->n/((real)p->Z*p->Z*p->Z);
  u1[1] = u1[2] = u2[0] = u2[2] = u3[0] = u3[1] = 0;
  u1[0] = u2[1] = u3[2] = 1;
  transform(u1, u1, p->ang, 1);
  transform(u2, u2, p->ang, 1);
  transform(u3, u3, p->ang, 1);
  p->matp[0] = u1[0];  p->matp[3] = u2[0];  p->matp[6] = u3[0];
  p->matp[1] = u1[1];  p->matp[4] = u2[1];  p->matp[7] = u3[1];
  p->matp[2] = u1[2];  p->matp[5] = u2[2];  p->matp[8] = u3[2];
  matrix_inverse(p->matp, p->matm);
  p->identity = 0;
  if (fabs(p->matp[0]+p->matp[4]+p->matp[8]-3)<1e-15)
    p->identity = 1;
}

void prep_check(MOLECULE *mol)
/* Calculate the maximum distance away from an orbital where there might be a
 *  surface.
 * Enter: MOLECULE *mol: pointer to molecular record.         12/18/93-DWM */
{
  long i, j, bad;
  real max, z, psi, low=mol->Psi*1e-30, th, ph, ct, st, sum, r, sep;
  OATOM *orb=mol->orb;

  mol->maxcheck = 0;  mol->maxpsi = 0;
  mol->syscntr[0] = mol->syscntr[1] = mol->syscntr[2] = 0;
  mol->maxjump = 1e7;
  for (i=0; i<mol->nump; i++)
    for (j=0; j<3; j++)
      mol->syscntr[j] += orb[i].x[j]/mol->nump;
  for (i=0; i<mol->nump; i++) {
    sep = sqrt(sq(mol->syscntr[0]-orb[i].x[0])+
               sq(mol->syscntr[1]-orb[i].x[1])+
               sq(mol->syscntr[2]-orb[i].x[2]));
    ph = sqrt2*orb[i].cc;
    if (!orb[i].m)  ph = orb[i].cc;
    for (th=psi=0; th<PI2; th+=0.0001) {
      ct = cos(th);
      for (j=sum=0; j<=orb[i].lmam/2 && ct; j++)
        sum += orb[i].cts[j]*powr(ct, orb[i].lmam-j-j);
      if (!ct && orb[i].lmam<2)  sum = 1;
      if (st=sin(th))  orb[i].ltf = powr(st, orb[i].am)*sum;
      else              orb[i].ltf = (!orb[i].m)*sum;
      if (fabs(orb[i].ltf)>psi)   psi = fabs(orb[i].ltf); }
    ph = fabs(ph*psi);
    z = bad = max = 0;
    while (bad<10*orb[i].n) {
      z += 0.01*orb[i].n;
      r = z*a0;
      for (j=sum=0; j<=orb[i].nml1; j++) {
        sum += orb[i].crs[j]*powr(r*orb[i].zna, orb[i].nml1-j);
        if (r>orb[i].maxr)
          bad = 10000; }
      psi = powr(r, orb[i].l)*exp(-r*orb[i].zna)*sum*ph;
      if (fabs(psi)<(LDBL_MIN*1e6)) {
        psi = 0;  bad += 2*orb[i].n; }
      psi *= psi;
      if (psi>mol->Psi/sq(mol->nump) && bad<9999) {
        if (psi*sq(mol->nump)>mol->maxpsi)
          mol->maxpsi = psi*sq(mol->nump);
        max = z + 0.04*orb[i].n*orb[i].n; }
      if (psi<low)  bad ++;
      else          bad = 0;
      if ((psi<mol->Psi*1e-10 || psi>mol->Psi*1e10) && z>4*orb[i].n*orb[i].n)
        z += 1*orb[i].n;
      if ((psi<mol->Psi*1e-20 || psi>mol->Psi*1e20) && z>4*orb[i].n*orb[i].n)
        z += 100*orb[i].n; }
    if (sep+max*a0>mol->maxcheck)
      mol->maxcheck = sep+max*a0;
    if (sqrt(max)*0.250<mol->maxjump)
      mol->maxjump = sqrt(max)*0.250; }
  mol->maxjump = mol->maxjump*mol->EffScale/(mol->nump*mol->maxjumpadj)/2;
  mol->maxcheck2 = sq(mol->maxcheck);
  mol->zsize = 10*a0 * 11.25/17.5/mol->EffScale/2;
  mol->checkpsi = mol->Psi;
}

long quickint(real *x, real *v, MOLECULE *mol, CUTAWAY *cut, real *y)
/* Check if there is an obstruction between a given point and the last
 *  surface point found.  Before tracking, the routine ensures that it is off
 *  of the last surface.
 * Enter: real *x: starting point, not modified.
 *        real *v: vector last travelled in, not modified.
 *        MOLECULE *mol: pointer to molecular record.
 *        CUTAWAY *cut: pointer to cutaway record, or null for none.
 *        real *y: location to head toward.
 * Exit: long found: 0 for none, 1 for hit.                     1/22/95-DWM */
{
  real res, point[3], newv[3], dist, psi1, psi2, newvl, mul, v2[3], last[3];
  long phase;

  mul = res = mol->minres + mol->lasthit*mol->back;
  point[0] = x[0]-res*v[0];
  point[1] = x[1]-res*v[1];
  point[2] = x[2]-res*v[2];
  newv[0] = y[0]-point[0];
  newv[1] = y[1]-point[1];
  newv[2] = y[2]-point[2];
  newvl = sqrt(newv[0]*newv[0]+newv[1]*newv[1]+newv[2]*newv[2]);
  dist = newvl/mol->zsize;
  newv[0] /= dist;  newv[1] /= dist;  newv[2] /= dist;
  psi1 = 0;
  while (dist-mul>0) {
    dist -= mul;
    point[0] += newv[0]*mul;
    point[1] += newv[1]*mul;
    point[2] += newv[2]*mul;
    psi2 = psi1;
    if ((psi1=calc_prob_total(point, mol, &phase))>mol->Psi) {
      memcpy(last, point, 3*sizeof(real));
      if (cut) if (cut_zone(point, cut, 0)) {
        memcpy(v2, newv, 3*sizeof(real));
        if (!cut_intercept(point, v2, cut))
          return(0);
        if ((psi1=calc_prob_total(point, mol, &phase))<=mol->Psi) {
          mul = res;
          dist -= sqrt(sq(last[0]-y[0])+sq(last[1]-y[1])+sq(last[2]-y[2]));
          continue; } }
      return(1); }
    if (mol->nump==1) {
      if (mol->orb[0].lr>=mol->maxcheck)
        return(0); }
    else {
      if (sq(point[0]-mol->syscntr[0])+sq(point[1]-mol->syscntr[1])+
          sq(point[2]-mol->syscntr[2])>=mol->maxcheck2)
        return(0); }
    mul = mol->zsize*calc_grad_mag(mol);
    if (!mul)  mul = res;
    if (psi1>psi2)
      mul = (mol->Psi-psi1)/mul;
    else
      mul = (mol->Psi+psi1)/mul;
    if (mul<res)
      mul = res;
    else if (mul>mol->maxjump)
      mul = mol->maxjump;
    res += mul*mol->back; }
  return(0);
}

long ray(real *x, real *v, real *colorin, MOLECULE *mol, CUTAWAY *cut,
         real *colorout)
/* Track a ray starting at x[] and heading in v[].  Determine what color it
 *  should be.
 * Enter: real *x: ray starting position.
 *        real *v: ray direction.
 *        real *colorin: array of three RGB real triplets containing the
 *                       colors of the positive phase, the negative phase,
 *                       and the background.  1.0==full intensity.
 *        MOLECULE *mol: pointer to molecular record.
 *        CUTAWAY *cut: pointer to cutaway record, or null for none.
 *        real *colorout: location to store RGB color.  Note that the values
 *                        can potentially exceed the [0,1] range.
 * Exit:  long hyperphase: this is a number which is intended to be unique
 *                         for each lobe and shadow condition.  Smooth
 *                         shading should only be applied between points
 *                         with identical hyperphases.  The actual location
 *                         of the surafce is stored in SurfaceXYZ.
 *                                                            12/12/93-DWM */
{
  long i, lastphase;
  real mag, g[3], magt=0, spc=0, l[3];
  LIGHT *ls=mol->ls;

  lastphase = intercept(x, v, mol, cut);
  switch (lastphase) {
    case -2: case -1: colorout[0] = colorin[3];
                      colorout[1] = colorin[4];
                      colorout[2] = colorin[5];  break;
    case  0:          colorout[0] = colorin[6];
                      colorout[1] = colorin[7];
                      colorout[2] = colorin[8];  return(0);
    case  1: case  2: colorout[0] = colorin[0];
                      colorout[1] = colorin[1];
                      colorout[2] = colorin[2]; }
  memcpy(SurfaceXYZ, x, 3*sizeof(real));
  if (lastphase!=-2 && lastphase!=2) {
    calc_grad_total(g, mol);
    lastphase *= 4; }
  else {
    memcpy(g, CutGrad, 3*sizeof(real));
    lastphase = (lastphase/2)*(1+CutSurf); }
  calc_unit(g);
  for (i=0; i<mol->numl; i++) {
    l[0] = ls[i].lx[0]-x[0];
    l[1] = ls[i].lx[1]-x[1];
    l[2] = ls[i].lx[2]-x[2];
    calc_unit(l);
    mag = (0.5-0.5*(g[0]*l[0]+g[1]*l[1]+g[2]*l[2]))*ls[i].i;
    if (ls[i].a!=1)
      if (quickint(x, v, mol, cut, ls[i].lx)) {
        lastphase *= 2;
        mag *= ls[i].a; }
    magt += mag; }
  if (!i)  magt = 1;
  if (magt<0)  magt = 0;   if (magt>1)  magt = 1;
  lastphase = ((long)magt)*128 + lastphase*256;
  for (i=0; i<mol->nump; i++)
    if (mol->orb[i].m)
      lastphase += 2048*(long)((mol->orb[i].lp+PI)*mol->orb[i].am/PI);
  colorout[0] *= magt;  colorout[1] *= magt;  colorout[2] *= magt;
  return(lastphase);
}

real ray_opacity(real *x, real *v, long steps, RENDER *re, MOLECULE *mol,
                 CUTAWAY *cut)
/* Track a ray starting at x[] and heading in v[] for a specific number of
 *  steps.  Determine what percentage of the light makes it through this
 *  area.  This does a slow tracking, taking into account opacity but not
 *  refraction.
 * Enter: real *x: ray starting position, modified.
 *        real *v: ray direction.  The magnitude must be the step size.  Not
 *                 modified.
 *        long steps: maximum number of steps to take.
 *        RENDER *re: render specification, including opacity.
 *        MOLECULE *mol: pointer to molecular record.
 *        CUTAWAY *cut: pointer to cutaway record, or null for none.
 * Exit:  real opacity: 0 for completely clear (all light transmitted), 1 for
 *                      completely opaque.                     12/1/97-DWM */
{
  long i=0, phase1, phase2;
  real opac=1, o, o2, o3, oa;
  real psi1, psi2, v2[3], gr[3], gra[3];
  real v3[3];

  memcpy(v3, v, 3*sizeof(real));
  if (!to_sphere(x, v3, mol->maxcheck, mol->zsize, mol->syscntr))
    return(0);
  if (cut_zone(x, cut, 0)) {
    memcpy(v2, v, 3*sizeof(real));
    if (!cut_intercept(x, v2, cut))  return(0);
    else {
      x[0] -= v[0]*2;  x[1] -= v[1]*2;  x[2] -= v[2]*2;
      to_sphere(x, v3, mol->maxcheck, mol->zsize, mol->syscntr); } }
  psi2 = calc_prob_total(x, mol, &phase2);
  while (sq(x[0]-mol->syscntr[0])+sq(x[1]-mol->syscntr[1])+
         sq(x[2]-mol->syscntr[2])<mol->maxcheck2 && opac>0.001 && i<steps) {
    x[0] += v[0];  x[1] += v[1];  x[2] += v[2];  i++;
    psi1 = psi2;  phase1 = phase2;
    if (cut_zone(x, cut, 0)) {
      psi2 = 0;  phase2 = 2*(phase1>0)-2*(phase1<0); }
    else
      psi2 = calc_prob_total(x, mol, &phase2);
    if (abs(phase2)>=2 && abs(phase1)>=2)  continue;
    o = o2 = o3 = oa = 0;
    if (re->opacity[3*(phase2<0)])                          /* Probability */
      o = (1-o)*re->opacity[3*(phase2<0)]*(psi2/mol->maxpsi);
    if ((psi1<mol->Psi && psi2>=mol->Psi) || (psi1>=mol->Psi &&
        psi2<mol->Psi))                                         /* Surface */
      o2 = re->opacity[1+3*(phase2<0)];
    if (psi2>mol->Psi && re->opacity[2+3*(phase2<0)])          /* Interior */
      o3 = (1-o-o2)*re->opacity[2+3*(phase2<0)];
    if (psi1 && psi2 && phase1*phase2<0 && re->opacity[6])    /* Asymptote */
      oa = 1;
    if (o <0.0001)  o  = 0;
    if (o2<0.0001)  o2 = 0;
    if (o3<0.0001)  o3 = 0;
    if (o || o2 || o3 || oa) {
      if (oa) {
        calc_grad_total(gr, mol);
        memcpy(v2, v, 3*sizeof(real));
        memcpy(gra, gr, 3*sizeof(real));
        calc_unit(v2);  calc_unit(gra);
        oa = fabs(v2[0]*gra[0]+v2[1]*gra[1]+v2[2]*gra[2]);
        if (oa)  oa = 1-(pow(1-re->opacity[6], 1./oa));
        else     oa = 1;
        oa = (1-o-o2-o3)*oa; }
      opac = (1-(o+o2+o3+oa))*opac; } }
  return(1-opac);
}

long ray_precise(real *x, real *v, real *colorin, RENDER *re, MOLECULE *mol,
                 CUTAWAY *cut, real *colorout)
/* Track a ray starting at x[] and heading in v[].  Determine what color it
 *  should be.  This does a slow tracking, taking into account opacity and
 *  refraction.
 * Enter: real *x: ray starting position.
 *        real *v: ray direction.
 *        real *colorin: array of four RGB real triplets containing the
 *                       colors of the positive phase, the negative phase,
 *                       the background, and the asymptote.  1.0==full
 *                       intensity.
 *        RENDER *re: render specification, including opacity, refraction,
 *                    and asymptote "thickness".
 *        MOLECULE *mol: pointer to molecular record.
 *        CUTAWAY *cut: pointer to cutaway record, or null for none.
 *        real *colorout: location to store RGB color.  Note that the values
 *                        can potentially exceed the [0,1] range.
 * Exit:  long hit: zero if the total opacity was less than 50%, 1 if >=50%.
 *                  The actual xyz coordinates where the opacity crossed 50%
 *                  is stored in SurfaceXYZ[3].               11/28/97-DWM */
{
  LIGHT *ls=mol->ls;
  long i, phase1, phase2, refrac=0;
  real oo=1, o, o2, o3, oa, n1, n2, shadow, old;
  real psi1, psi2, v2[3], mag, magt, maga, magat, l[3], gr[3], gra[3];
  real x2[3], v3[3], st2, ct2;

  if (!to_sphere(x, v, mol->maxcheck, mol->zsize, mol->syscntr)) {
    memcpy(colorout, colorin+6, 3*sizeof(real));
    mol->lasthit = 0;  return(0); }
  colorout[0] = colorout[1] = colorout[2] = 0;
  if (!re->steps)  re->steps = 1000;
  calc_unit(v);
  v[0] = v2[0] = v[0]*mol->maxcheck/re->steps;
  v[1] = v2[1] = v[1]*mol->maxcheck/re->steps;
  v[2] = v2[2] = v[2]*mol->maxcheck/re->steps;
  if (cut) {
    if (cut_zone(x, cut, 0)) {
      if (!cut_intercept(x, v2, cut)) {
        memcpy(colorout, colorin+6, 3*sizeof(real));
        mol->lasthit = 0;  return(0); }
      else {
        x[0] -= v[0]*2;  x[1] -= v[1]*2;  x[2] -= v[2]*2;
        if (!to_sphere(x, v2, mol->maxcheck, mol->zsize, mol->syscntr)) {
          memcpy(colorout, colorin+6, 3*sizeof(real));
          mol->lasthit = 0;  return(0); } } }
    else {
      memcpy(x2, x, 3*sizeof(real));
      if (!cut->notlast)  if (!cut_intercept(x2, v2, cut))  cut = 0; } }
  psi2 = calc_prob_total(x, mol, &phase2);
  for (i=0; i<3; i++)
    refrac |= (re->refrac[i]!=0 && re->refrac[i]!=1);
  if (refrac) {
    n2 = re->refrac[phase2<0]*(psi2>mol->Psi &&
         (re->opacity[1+3*(phase2<0)] || re->opacity[2+3*(phase2<0)]));
    if (re->opacity[3*(phase2<0)] && re->refrac[phase2<0])
      n2 = (re->refrac[phase2<0]-1)*(psi2/mol->maxpsi)+1;
    if (!n2)  n2 = 1; }
  while (sq(x[0]-mol->syscntr[0])+sq(x[1]-mol->syscntr[1])+
         sq(x[2]-mol->syscntr[2])<mol->maxcheck2 && oo>0.001) {
    x[0] += v[0];  x[1] += v[1];  x[2] += v[2];
    psi1 = psi2;  phase1 = phase2;  n1 = n2;
    if (cut_zone(x, cut, 0)) {
      psi2 = 0;  phase2 = 2*(phase1>0)-2*(phase1<0);  n2 = 1; }
    else
      psi2 = calc_prob_total(x, mol, &phase2);
    if (abs(phase2)>=2 && abs(phase1)>=2)  continue;
    if (cut) if (cut->nosurface)
      if (abs(phase2)>=2 || abs(phase1)>=2)  continue;
    o = o2 = o3 = oa = 0;
    if (re->opacity[3*(phase2<0)])                          /* Probability */
      o = (1-o)*re->opacity[3*(phase2<0)]*(psi2/mol->maxpsi);
    if ((psi1<mol->Psi && psi2>=mol->Psi) || (psi1>=mol->Psi &&
        psi2<mol->Psi))                                         /* Surface */
      o2 = re->opacity[1+3*(phase2<0)];
    if (psi2>mol->Psi && re->opacity[2+3*(phase2<0)])          /* Interior */
      o3 = (1-o-o2)*re->opacity[2+3*(phase2<0)];
    if (psi1 && psi2 && phase1*phase2<0 && re->opacity[6])    /* Asymptote */
      oa = 1;
    if (o <0.0001)  o  = 0;
    if (o2<0.0001)  o2 = 0;
    if (o3<0.0001)  o3 = 0;
    if (refrac) {
      n2 = re->refrac[phase2<0]*(psi2>mol->Psi &&
           (re->opacity[1+3*(phase2<0)] || re->opacity[2+3*(phase2<0)]));
      if (re->opacity[3*(phase2<0)] && re->refrac[phase2<0])
        n2 = (re->refrac[phase2<0]-1)*(psi2/mol->maxpsi)+1;
      if (!n2 && oa)
        n2 = re->refrac[2];
      if (!n2)  n2 = 1; }
    if (o || o2 || o3 || oa) {
      if (((o || o2) && mol->numl) || oa || n1!=n2) {
        if ((abs(phase1)<2 && abs(phase2)<2) || oa)
          calc_grad_total(gr, mol);
        if (oa) {
          memcpy(v2, v, 3*sizeof(real));
          memcpy(gra, gr, 3*sizeof(real));
          calc_unit(v2);  calc_unit(gra);
          oa = fabs(v2[0]*gra[0]+v2[1]*gra[1]+v2[2]*gra[2]);
          if (oa)  oa = 1-(pow(1-re->opacity[6], 1./oa));
          else     oa = 1;
          oa = (1-o-o2-o3)*oa; }
        if (abs(phase1)>=2 || abs(phase2)>=2) {
          memcpy(x2, x, 3*sizeof(real));
          memcpy(v2, v, 3*sizeof(real));
          if (abs(phase2)<2) {
            x2[0] -= v2[0];  x2[1] -= v2[1];  x2[2] -= v2[2]; }
          cut_intercept(x2, v2, cut);
          memcpy(gr, CutGrad, 3*sizeof(real)); }
        calc_unit(gr);
        for (i=0, magt=magat=0; i<mol->numl; i++) {
          l[0]=ls[i].lx[0]-x[0];
          l[1]=ls[i].lx[1]-x[1];
          l[2]=ls[i].lx[2]-x[2];
          calc_unit(l);
          mag = (0.5-0.5*(gr[0]*l[0]+gr[1]*l[1]+gr[2]*l[2]))*ls[i].i;
          if (oa)
            maga = (0.5-0.5*(gra[0]*l[0]+gra[1]*l[1]+gra[2]*l[2]))*ls[i].i;
          if (ls[i].a!=1) {
            memcpy(x2, x, 3*sizeof(real));
            if (psi2>=mol->Psi && o2) {
              x2[0] -= v[0];  x2[1] -= v[1];  x2[2] -= v[2]; }
            v2[0] = ls[i].lx[0]-x2[0];
            v2[1] = ls[i].lx[1]-x2[1];
            v2[2] = ls[i].lx[2]-x2[2];
            shadow = sqrt(sq(v2[0])+sq(v2[1])+sq(v2[2]))/
                     (mol->maxcheck/re->steps);
            if (shadow>re->steps*2)  shadow = re->steps*2;
            calc_unit(v2);
            v2[0] *= mol->maxcheck/re->steps;
            v2[1] *= mol->maxcheck/re->steps;
            v2[2] *= mol->maxcheck/re->steps;
            shadow = ray_opacity(x2, v2, (long)shadow, re, mol, cut);
            shadow = 1-shadow+ls[i].a*shadow;
            mag *= shadow;
            if (oa) maga *= shadow; }
          magt += mag;
          if (oa)  magat += maga; }
        if (!mol->numl)  magt = magat = 1;
        if (magt<0)  magt = 0;   if (magt>1)  magt = 1;
        if (refrac) {
          if (n1!=n2) {
            calc_unit(cross(gr, v, v3));              /* Basis: v2, v3, gr */
            calc_unit(cross(v3, gr, v2));
            memcpy(v3, v, 3*sizeof(real));  calc_unit(v3);
            st2 = n1/n2 * (v3[0]*v2[0]+v3[1]*v2[1]+v3[2]*v2[2]);
            if (st2<-1 || st2>1) {                  /* Internal reflection */
              st2 = v3[0]*v2[0]+v3[1]*v2[1]+v3[2]*v2[2];
              ct2 = -(v3[0]*gr[0]+v3[1]*gr[1]+v3[2]*gr[2]); }
            else {
              ct2 = cos(asin(st2));
              if (ct2*(v3[0]*gr[0]+v3[1]*gr[1]+v3[2]*gr[2])<0)  ct2 *= -1; }
            v[0] = (st2*v2[0]+ct2*gr[0])*mol->maxcheck/re->steps;
            v[1] = (st2*v2[1]+ct2*gr[1])*mol->maxcheck/re->steps;
            v[2] = (st2*v2[2]+ct2*gr[2])*mol->maxcheck/re->steps; } } }
      else
        magt = 1;
      if (oa) {
        if (magat<0)  magat = 0;   if (magat>1)  magat = 1;
        maga = oa*magat;
        colorout[0] += colorin[ 9]*oo*maga;
        colorout[1] += colorin[10]*oo*maga;
        colorout[2] += colorin[11]*oo*maga; }
      magt = (o+o2)*magt+o3;
      colorout[0] += colorin[  3*(phase2<0)]*oo*magt;
      colorout[1] += colorin[1+3*(phase2<0)]*oo*magt;
      colorout[2] += colorin[2+3*(phase2<0)]*oo*magt;
      old = oo;
      o = o+o2+o3+oa;
      oo = (1-o)*oo;
      if (old>0.5 && oo<=0.5)
        memcpy(SurfaceXYZ, x, 3*sizeof(real)); } }
  colorout[0] += oo*colorin[6];
  colorout[1] += oo*colorin[7];
  colorout[2] += oo*colorin[8];
  if (colorout[0]*re->brightness>1)  colorout[0] = 1./re->brightness;
  if (colorout[1]*re->brightness>1)  colorout[1] = 1./re->brightness;
  if (colorout[2]*re->brightness>1)  colorout[2] = 1./re->brightness;
  return(oo<=0.5);
}

real rnd(real low, real high, long inc)
/* Pick a random number between two values.
 * Enter: real low, high: bounds for the random number.  Low must be less
 *                        than high.
 *        long inc: inclusion: 0-(low,high), 1-[), 2-(], 3-[].
 * Exit:  real rand: random number.                             7/6/97-DWM */
{
  real r, s;
  static long ever=0;
  long i;

  if (high<=low)  return(low);
  if (!ever) {
    srand((unsigned int)((time(0)<<16)+(time(0)>>16)));
    ever = 1; }
  do {
    r = 0;
    s = (real)RAND_MAX;
    for (i = 0; i < 4; i++) {
      r += (real)rand() / s;
      s *= ((real)RAND_MAX + 1);
    }
  } while (r > 1 || (r == 0 && (inc == 0 || inc == 2)) ||
           (r >= 1 && inc <= 1));
  return r * (high - low) + low;
}

void snap_to_surface(float *xyz, MOLECULE *mol, real psi2, real uncertain)
/* Adjust a point so that it lies on the surface of the part.
 * Enter: float *xyz: point to modify.  This must be "close" to the
 *                    surface.
 *        MOLECULE *mol: molecule to analyze.
 *        real psi2: probability surface to locate.
 *        real uncertain: uncertainty of initial surface guess in meters.
 *                                                            11/11/97-DWM */
{
  real x[3], v[3], mul=uncertain, psi1=0, total=0;
  long phase;

  x[0] = xyz[0];  x[1] = xyz[1];  x[2] = xyz[2];
  psi1 = calc_prob_total(x, mol, &phase);
  psi1 *= phase;
  calc_grad_total(v, mol);
  calc_unit(v);
  if (phase<0) {
    v[0] *= -1;  v[1] *= -1;  v[2] *= -1; }
  while (mul>mol->minres*mol->zsize) {
    mul /= 2;
    if (psi1<psi2) {
      total += mul;
      x[0] += v[0]*mul;  x[1] += v[1]*mul;  x[2] += v[2]*mul; }
    else {
      total -= mul;
      x[0] -= v[0]*mul;  x[1] -= v[1]*mul;  x[2] -= v[2]*mul; }
    psi1 = calc_prob_total(x, mol, &phase);
    psi1 *= phase; }
  xyz[0] = x[0];  xyz[1] = x[1];  xyz[2] = x[2];
}

real to_jump(void)
/* Pass back the distance moved in the to_sphere() routine.
 * Exit:  real dist: actual distance moved.                    1/22/95-DWM */
{
  return(ToJump);
}

long to_sphere(real *x, real *v, real radius, real siz, real *cntr)
/* Based on a ray starting a point x and extending in direction v, calculate
 *  the point where it first reaches the critical radius away from any
 *  orbital.  Use to_jump() to retreive the actual distance moved.
 * Enter: real *x: starting coordinates.  Final coordinates are stored here.
 *        real *v: initial vector.  Normalized then multiplied by siz.
 *        real radius: critical radius.
 *        real siz: distance to pass critical radius and scale factor for v.
 *        real *cntr: sphere center.
 * Exit:  long hit: 0 for missed, 1 for within critical radius.  1/5/94-DWM */
{
  real p[3], min, cost, curr;

  calc_unit(v);
  p[0] = cntr[0]-x[0];  p[1] = cntr[1]-x[1];  p[2] = cntr[2]-x[2];
  min = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
  cost = p[0]*v[0]+p[1]*v[1]+p[2]*v[2];
  curr = radius*radius+(cost*cost-min);
  if (curr<0)
    return(0);
  curr = cost - sqrt(curr);
  if (curr<0) {
    ToJump = 0;
    v[0] *= siz;        v[1] *= siz;          v[2] *= siz;
    return(1); }
  curr += siz;
  x[0] += v[0]*curr;    x[1] += v[1]*curr;    x[2] += v[2]*curr;
  v[0] *= siz;          v[1] *= siz;          v[2] *= siz;
  ToJump = curr;
  return(1);
}

real *transform(real *x, real *t, real *a, real mag)
/* Transform a set of coordinates by a set of angles.
 * Enter: real *x: initial coordinates.
 *        real *t: transformed coordinates.
 *        real *a: set of three angles.
 *        real mag: magnitude of transformation (1 for positive, -1 for
 *                  negative).
 * Exit:  real *t: transformed coordinates.                   12/11/93-DWM */
{
  real ca=cos(a[0]*mag), cb=cos(a[1]*mag), cg=cos(a[2]*mag);
  real sa=sin(a[0]*mag), sb=sin(a[1]*mag), sg=sin(a[2]*mag);
  real tx, ty, tz;

  if (mag<0) {
    tx = sa;  sa = sg;  sg = tx;
    tx = ca;  ca = cg;  cg = tx; }
  tx = x[0]*(cg*ca-sg*cb*sa)+x[1]*(-cg*sa-sg*cb*ca)+x[2]*  sg*sb;
  ty = x[0]*(sg*ca+cg*cb*sa)+x[1]*(-sg*sa+cg*cb*ca)+x[2]*(-cg*sb);
  tz = x[0]* sb*sa          +x[1]*  sb*ca          +x[2]*  cb;
  t[0] = tx;  t[1] = ty;  t[2] = tz;
  return(t);
}

real *transform2(real *x, real *t, real *m)
/* Transform a set of coordinates by a matrix.
 * Enter: real *x: initial coordinates.
 *        real *t: transformed coordinates.
 *        real *m: transformation matrix.
 * Exit:  real *t: transformed coordinates.                   12/18/93-DWM */
{
  #ifdef ANSIC
  real tx, ty, tz;

  tx = x[0]*m[0] + x[1]*m[3] + x[2]*m[6];
  ty = x[0]*m[1] + x[1]*m[4] + x[2]*m[7];
  tz = x[0]*m[2] + x[1]*m[5] + x[2]*m[8];
  t[0] = tx;  t[1] = ty;  t[2] = tz;
  #else
  _asm {
          mov esi, x
          mov edi, m
          fld QWORD PTR [esi]
          fld QWORD PTR [edi]
          fmul
          fld QWORD PTR [esi+8]
          fld QWORD PTR [edi+0x18]
          fmul
          fadd
          fld QWORD PTR [esi+0x10]
          fld QWORD PTR [edi+0x30]
          fmul
          fadd
          fld QWORD PTR [esi]
          fld QWORD PTR [edi+0x08]
          fmul
          fld QWORD PTR [esi+8]
          fld QWORD PTR [edi+0x20]
          fmul
          fadd
          fld QWORD PTR [esi+0x10]
          fld QWORD PTR [edi+0x38]
          fmul
          fadd
          fld QWORD PTR [esi]
          fld QWORD PTR [edi+0x10]
          fmul
          fld QWORD PTR [esi+8]
          fld QWORD PTR [edi+0x28]
          fmul
          fadd
          fld QWORD PTR [esi+0x10]
          fld QWORD PTR [edi+0x40]
          fmul
          fadd
          mov edi, t
          fstp QWORD PTR [edi+0x10]
          fstp QWORD PTR [edi+8]
          fstp QWORD PTR [edi]
    }
  #endif
  return(t);
}

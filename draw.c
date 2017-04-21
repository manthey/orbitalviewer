/* This file contains all routines for updating the window.  It also handles
 *  mouse movement, both in the main window and in preview dialog.
 *                                                             6/11/97-DWM */

#include <windows.h>
#include <commctrl.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "matrix.h"
#include "dlt.h"
#include "preview.h"
#include "ovrc.h"
#include "ov.h"

long CutZ, Dither=1, dithernum=-1, dither16num=-1, dithertbl2[900], Down=0,
     DownX, DownY, DownBuf[7], LastMove=0, opacitytbl[OPACDITHER], OpacPos=0;
long CameraCoor[]={76,4,100,100}, CutawayCoor[]={76,4,100,100},
     LightCoor[]={96,34,100,100};
char *DispInfo[]={"Use mouse to rotate/pan/zoom."};
char *DispMode[]={"Render Mode"};
uchar *ColorTable, dither16tbl[2048], MasterPal[768];
uchar imask[] = {0x80, 0x40, 0x20, 0x10, 8, 4, 2, 1};
uchar dithertbl1[2048], dithertbl3[]={0,0,1,2,3,4,5,6,7,8,9,10,12,14,16,18};

bgr_to_rgb(uchar *image, long siz)
/* Convert a packed array of BGR values to RGB values.
 * Enter: uchar *image: pointer to array to convert.
 *        long size: number of values to convert.  Each value is 3 bytes.
 *                                                             2/26/96-DWM */
{
  _asm {
          mov esi, image
          mov ecx, siz
bgrrgb1:  mov al, [esi]
          xchg al, [esi+2]
          mov [esi], al
          add esi, 03
          loop   bgrrgb1
    }
}

camera_figure(HWND hdlg, HDC hdc, DATA *data, double *phys)
/* Create the figure showing the current state of the camera.
 * Enter: HWND hdlg: handle to dialog to work with.
 *        HDC hdc: device context handle to draw with.
 *        DATA *data: pointer to owner window.
 *        double *phys: array containing modified camera parameters.
 *                                                              1/5/98-DWM */
{
  long w, h, backclr, *lpic, s, i, j, x, y;
  RECT rect;
  HBITMAP wpic, wpic2;
  uchar *pic2;
  LOGPALETTE *lpal;
  long tim=clock();
  static long init=0;
  double fac;

  if (!init) {
    init = 1;
    for (i=0; i<26; i++) {
      for (j=0; j<4; j++) {
        CamElem[i*15+j+26*15] = CamElem[i*15+j]+32;
        CamElem[i*15+j+52*15] = CamElem[i*15+j]+64; }
      for (j=4; j<12; j+=2) {
        x = CamElem[i*15+j];  y = CamElem[i*15+j+1];
        if (y<30)  x += 90;
        else       x += 30;
        CamElem[i*15+j+26*15] = x;  CamElem[i*15+j+26*15+1] = y;
        x = CamElem[i*15+j];  y = CamElem[i*15+j+1];
        if (y<30)                           y += 30;
        else if (x<30 || (x>=90 && x<120))  x += 60;
        else                                x -= 30;
        CamElem[i*15+j+52*15] = x;  CamElem[i*15+j+52*15+1] = y; }
      for (j=0; j<3; j++) {
        CamElem[i*15+26*15+12+((j+1)%3)] = CamElem[i*15+12+j];
        CamElem[i*15+52*15+12+((j+2)%3)] = CamElem[i*15+12+j]; } } }
  if (data->mol.maxcheck) {
    fac = data->mol.maxcheck/-CamXYZ[0];
    for (i=0; i<96*3; i++)
      CamXYZ[i] *= fac; }
  camera_rect(&rect, hdlg);
  w = rect.right-rect.left;  h = rect.bottom-rect.top;
  s = ((w*3+3)/4)*4;
  wpic = GlobalAlloc(GMEM_MOVEABLE, s*h+40);
  if (!wpic)  return(0);
  pic2 = GlobalLock(wpic);  lpic = (long *)pic2;
  lpic[0] = 40;  lpic[1] = w;  lpic[2] = h;
  ((short *)pic2)[6] = 1;  ((short *)pic2)[7] = 24;
  lpic[4] = lpic[6] = lpic[7] = 0;  lpic[5] = s*h;
  lpic[8] = lpic[9] = 0;
  backclr = GetSysColor(COLOR_BTNFACE);
  for (i=0; i<h; i++)
    fill_zone24(pic2+40+s*i, backclr, w);
  camera_figure_draw(w, h, -s, pic2+40+(h-1)*s, data, phys);
  GlobalUnlock(wpic);
  wpic = reduce_color_space(wpic, hdc);
  draw_bmp(hdc, rect.left, rect.top, wpic);
  GlobalFree(wpic);
  tim = 1000*(tim-clock())/CLOCKS_PER_SEC;
  if (tim<100)       PreviewDelay = 100;
  else if (tim<200)  PreviewDelay = 500;
  else if (tim>1000) PreviewDelay = 2000;
  else               PreviewDelay = 1000;
}

camera_figure_draw(long w, long h, long scan, uchar *buf, DATA *data,
                   double *phys)
/* Draw the camera figure to a buffer.
 * Enter: long w, h: size of figure to draw.
 *        long scan: bytes per scan line.
 *        uchar *buf: pointer to first scan line.
 *        DATA *data: pointer to owner window.
 *        double *phys: array containing modified camera parameters.
 *                                                              1/6/98-DWM */
{
  ushort *zbuf;
  long sx[96], sy[96], sz[96], i, j, k, coor[15];
  long numnode=96, numsq=26*3, tbl[]={0,1,3,2,3,1};
  float *fz;
  double val[10], DLT[18], phys2[10], bright;
  double L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, temp[10], tr;
  double mc[9], den, minz=1e30, maxz=-1e30, cz, x, y, z;
  lmatrix m;

  if (!(zbuf=malloc2(w*h*sizeof(short))))  return(0);
  memset(zbuf, 0xFFFF, w*h*sizeof(short));
  fz = (float *)sz;
  memcpy(temp, phys, 10*sizeof(double));
  memcpy(val, data->renderval, 10*sizeof(double));
  langtomat(mc, temp[7]-PI, temp[8]-PI, temp[9]);
  temp[2] -= fabs(val[0])*val[1];
  val[4] = -(temp[0]*mc[0]+temp[1]*mc[3]+temp[2]*mc[6]);
  val[5] = -(temp[0]*mc[1]+temp[1]*mc[4]+temp[2]*mc[7]);
  val[6] = -(temp[0]*mc[2]+temp[1]*mc[5]+temp[2]*mc[8]);
  render_physical(data->renderval, data->renderdlt, phys2, 0);
  phys2[7] = temp[7]-PI;  phys2[8] = temp[8]-PI;  phys2[9] = temp[9];
  render_move(0,0,0, -w,h, val, phys2, DLT);
  m.m = mc;  ldlt_to_physical(DLT, temp, &m, 0);
  tr = mc[6]*temp[0]+mc[7]*temp[1]+mc[8]*temp[2];
  L1 = DLT[0];  L2 = DLT[1];  L3 = DLT[2];  L4 = DLT[3];  L5 = DLT[4];
  L6 = DLT[5];  L7 = DLT[6];  L8 = DLT[7];  L9 = DLT[8];  L10 = DLT[9];
  L11 = DLT[10];
  for (i=j=0; i<numnode; i++,j+=3) {
    den = L9*CamXYZ[j]+L10*CamXYZ[j+1]+L11*CamXYZ[j+2]+1;
    if (den) {
      fz[i] = mc[6]*CamXYZ[j]+mc[7]*CamXYZ[j+1]+mc[8]*CamXYZ[j+2]-tr;
      if (fz[i]<minz)  minz = fz[i];  if (fz[i]>maxz)  maxz=fz[i];
      den = 1./den;
      sx[i] = (L1*CamXYZ[j]+L2*CamXYZ[j+1]+L3*CamXYZ[j+2]+L4)*den;
      sy[i] = (L5*CamXYZ[j]+L6*CamXYZ[j+1]+L7*CamXYZ[j+2]+L8)*den; } }
  if (maxz==minz)  maxz = minz+1;
  maxz = 65335./(maxz-minz);
  for (i=0; i<numnode; i++)
    sz[i] = (fz[i]-minz)*maxz+100;
  cz = (-minz*maxz+100+1)*256;
  if (cz<LONG_MIN)       CutZ = LONG_MIN;
  else if (cz>LONG_MAX)  CutZ = LONG_MAX;
  else                   CutZ = cz;
  for (i=0; i<numsq; i++)
    for (j=0; j<6; j+=3) {
      for (k=0; k<3; k++) {
        coor[k*5]   = sx[CamElem[i*15+tbl[k+j]]];
        coor[k*5+1] = sy[CamElem[i*15+tbl[k+j]]];
        coor[k*5+2] = sz[CamElem[i*15+tbl[k+j]]];
        coor[k*5+3] = CamElem[i*15+tbl[k+j]*2+4];
        coor[k*5+4] = CamElem[i*15+tbl[k+j]*2+5]; }
      x=CamElem[i*15+12]*mc[0]+CamElem[i*15+13]*mc[1]+CamElem[i*15+14]*mc[2];
      y=CamElem[i*15+12]*mc[3]+CamElem[i*15+13]*mc[4]+CamElem[i*15+14]*mc[5];
      z=CamElem[i*15+12]*mc[6]+CamElem[i*15+13]*mc[7]+CamElem[i*15+14]*mc[8];
      bright = (x+y-z)*0.28867+0.5;
      render_triangle_texture(w, h, scan, buf, w, zbuf, coor, PREVIEWWIDTH,
                              90, 0, PreviewGraphic+3*PREVIEWWIDTH*
                              PREVIEWHEIGHT*PREVIEWNUM, 1, 1, bright); }
  free2(zbuf);
}

camera_move(HWND hdlg, long delx, long dely, long w, long h, long down,
            DATA *data, long ptr)
/* Actually move the camera image, modifying the phys array.
 * Enter: HWND hdlg: handle of current dialog.
 *        long delx, dely: amount of mouse movement.
 *        long w, h: preview image size.
 *        long down: 1-rotate, 2-zoom, 3-shift.
 *        DATA *data: pointer to current window.
 *        long ptr: actually double *phys.  This is the modified camera
 *                  parameters.
 * Exit:  long moved: 0 for no change in object position, 1 for change.
 *                                                              1/6/98-DWM */
{
  double temp[10], val[10], phys2[10], dlt[18];
  float dx[3]={0,0,0}, dang[3]={0,0,0};
  double mc[9];
  long mod=0;

  memcpy(temp, (double *)ptr, 10*sizeof(double));
  memcpy(val, data->renderval, 10*sizeof(double));
  langtomat(mc, temp[7]-PI, temp[8]-PI, temp[9]);
  temp[2] -= fabs(val[0])*val[1];
  val[4] = -(temp[0]*mc[0]+temp[1]*mc[3]+temp[2]*mc[6]);
  val[5] = -(temp[0]*mc[1]+temp[1]*mc[4]+temp[2]*mc[7]);
  val[6] = -(temp[0]*mc[2]+temp[1]*mc[5]+temp[2]*mc[8]);
  render_physical(data->renderval, data->renderdlt, phys2, 0);
  phys2[7] = temp[7]-PI;  phys2[8] = temp[8]-PI;  phys2[9] = temp[9];
  render_move(0,0,0, -w,h, val, phys2, dlt);
  switch (down) {
    case 1: dang[0] = -PI*dely/w;  dang[1] = -PI*delx/h; break;
    case 2: dx[2] = 10.*dely/h;  break;
    case 3: dx[0] = (double)delx/w;  dx[1] = (double)dely/h; }
  if (dx[0] || dx[1] || dx[2] || dang[0] || dang[1])  mod = 1;
  if (mod) {
    render_move(dx, dang, 0, 0,0, val, dlt, dlt);
    camera_rotate(val, dlt, temp, 0);
    memcpy((double *)ptr, temp, 3*sizeof(double));
    memcpy(((double *)ptr)+7, temp+7, 3*sizeof(double)); }
  return(mod);
}

camera_rect(RECT *rect, HWND hdlg)
/* Locate the rectangle used in a camera preview.
 * Enter: RECT *rect: pointer to store result.  If null, the rectangle is
 *                    invalidated (refreshed).
 *        HWND hdlg: dialog the rectangle is located in.       10/7/96-DWM */
{
  long inv=0;
  RECT rect2;

  if (!rect) { rect = &rect2;  inv = 1; }
  rect->left = CameraCoor[0];
  rect->top = CameraCoor[1];
  rect->right = CameraCoor[0]+CameraCoor[2];
  rect->bottom = CameraCoor[1]+CameraCoor[3];
  MapDialogRect(hdlg, rect);
  if (inv)
    InvalidateRect(hdlg, rect, 0);
}

cutaway_figure(HWND hdlg, HDC hdc, real *pos, long type)
/* Create the figure showing the current state of the cutaway.
 * Enter: HWND hdlg: handle to dialog to work with.
 *        HDC hdc: device context handle to draw with.
 *        real *pos: x,y,z,alpha,beta,gamma,radius,nosurface,invert
 *                   specifying cutaway.
 *        long type: 0-none, 1-plane, 2-corner, 3-wedge.       7/26/97-DWM */
{
  long w, h, backclr, *lpic, s, i, lastpart=8;
  RECT rect;
  HBITMAP wpic, wpic2;
  uchar *pic2;
  LOGPALETTE *lpal;
  long tim=clock();

  cutaway_rect(&rect, hdlg);
  w = rect.right-rect.left;  h = rect.bottom-rect.top;
  if (h>w)  h = w;  if (w>h)  w = h;
  backclr = GetSysColor(COLOR_BTNFACE);
  s = ((w*3+3)/4)*4;
  wpic = GlobalAlloc(GMEM_MOVEABLE, s*h+40);
  if (!wpic)  return(0);
  pic2 = GlobalLock(wpic);  lpic = (long *)pic2;
  lpic[0] = 40;  lpic[1] = w;  lpic[2] = h;
  ((short *)pic2)[6] = 1;  ((short *)pic2)[7] = 24;
  lpic[4] = lpic[6] = lpic[7] = 0;  lpic[5] = s*h;
  lpic[8] = lpic[9] = 0;
  cutaway_figure_draw(w, h, s, pic2+40, pos, type, backclr);
  GlobalUnlock(wpic);
  wpic = reduce_color_space(wpic, hdc);
  draw_bmp(hdc, rect.left, rect.top, wpic);
  GlobalFree(wpic);
  tim = 1000*(tim-clock())/CLOCKS_PER_SEC;
  if (tim<100)       PreviewDelay = 100;
  else if (tim<200)  PreviewDelay = 500;
  else if (tim>1000) PreviewDelay = 2000;
  else               PreviewDelay = 1000;
}

cutaway_figure_draw(long w, long h, long scan, uchar *buf, real *pos,
                    long type, long color)
/* Draw the cutaway figure to a buffer.
 * Enter: long w, h: size of figure to draw.
 *        long scan: bytes per scan line.
 *        uchar *buf: pointer to first scan line.
 *        real *pos: x,y,z,alpha,beta,gamma,radius,nosurface,invert
 *                   specifying cutaway.
 *        long type: 0-none, 1-plane, 2-corner, 3-wedge.
 *        long color: background color.                        7/29/97-DWM */
{
  long i, j, k, l, d, delx, dely, delxy, total[3];
  uchar buf2[3];
  real x[3], v[3];
  DATA data;

  memset(&data, 0, sizeof(DATA));
  memcpy(&data.pref, &DispData->pref, sizeof(PREF));
  memcpy(data.cut.xyz, pos, 3*sizeof(real));
  memcpy(data.cut.ang, pos+3, 3*sizeof(real));
  data.cut.type = type;
  data.cut.invert = pos[8];
  for (j=0; j<h; j++) {
    d = j*scan;
    for (i=0; i<w; i++, d+=3) {
      x[0] = x[1] = 0;  x[2] = pos[6]*5;
      v[0] = i-0.5*w;  v[1] = j-0.5*h;  v[2] = -2.4*w;
      cutaway_figure_point(x, v, pos, color, &data, buf+d);
      if (i && j && BitsPixel>8) {
        delx = abs(buf[d]-buf[d-3])+abs(buf[d+1]-buf[d-2])+
               abs(buf[d+2]-buf[d-1]);
        dely = abs(buf[d]-buf[d-scan])+abs(buf[d+1]-buf[d+1-scan])+
               abs(buf[d+2]-buf[d+2-scan]);
        delxy = abs(buf[d]-buf[d-scan-3])+abs(buf[d+1]-buf[d-scan-2])+
               abs(buf[d+2]-buf[d-scan-1]);
        if ((delx>8 && dely>8) || delxy>8) {
          total[0] = buf[d];  total[1] = buf[d+1];  total[2] = buf[d+2];
          for (k=0; k<4; k++)  for (l=0; l<4; l++)  if (k || l) {
            x[0] = x[1] = 0;  x[2] = pos[6]*5;
            v[0] = i-0.5*w-k*0.25;  v[1] = j-0.5*h-l*0.25;  v[2] = -2.4*w;
            cutaway_figure_point(x, v, pos, color, &data, buf2);
            total[0]+=buf2[0];  total[1]+=buf2[1];  total[2]+=buf2[2];  }
          buf[d]   = (total[0]+8)/16;
          buf[d+1] = (total[1]+8)/16;
          buf[d+2] = (total[2]+8)/16; }
        else if (delx>8) {
          total[0] = buf[d];  total[1] = buf[d+1];  total[2] = buf[d+2];
          for (k=0; k<4; k++)  if (k) {
            x[0] = x[1] = 0;  x[2] = pos[6]*5;
            v[0] = i-0.5*w-k*0.25;  v[1] = j-0.5*h;  v[2] = -2.4*w;
            cutaway_figure_point(x, v, pos, color, &data, buf2);
            total[0]+=buf2[0];  total[1]+=buf2[1];  total[2]+=buf2[2];  }
          buf[d]   = (total[0]+2)/4;
          buf[d+1] = (total[1]+2)/4;
          buf[d+2] = (total[2]+2)/4; }
        else if (dely>8) {
          total[0] = buf[d];  total[1] = buf[d+1];  total[2] = buf[d+2];
          for (l=0; l<4; l++)  if (l) {
            x[0] = x[1] = 0;  x[2] = pos[6]*5;
            v[0] = i-0.5*w;  v[1] = j-0.5*h-l*0.25;  v[2] = -2.4*w;
            cutaway_figure_point(x, v, pos, color, &data, buf2);
            total[0]+=buf2[0];  total[1]+=buf2[1];  total[2]+=buf2[2];  }
          buf[d]   = (total[0]+2)/4;
          buf[d+1] = (total[1]+2)/4;
          buf[d+2] = (total[2]+2)/4; } } } }
}

cutaway_figure_point(real *x, real *v, real *pos, long color, DATA *data,
                     uchar *dest)
/* Calculate the color of the cutaway along a specific ray.
 * Enter: real *x: starting point of ray.
 *        real *v: direction vector of ray.
 *        real *pos: x,y,z,alpha,beta,gamma,radius,nosurface,invert
 *                   specifying cutaway.
 *        long color: background color of cutaway figure.
 *        DATA *data: memory area containing cutaway info.
 *        uchar *dest: location to store color.                8/22/97-DWM */
{
  long i, j, k, l, d, num, num2, clr;
  real inter[30], x2[3], v2[3], inter2[30], cost, acost, mul;
  real r, g, b, r2, g2, b2, o;
  float oo=0.7, rr, gg, bb, amb=0.6;  /* opacity, color, ambient */

  if ((num=cutaway_inter(x, v, pos, data, inter))<1) {
    dest[0] = color&0xFF;  ((ushort *)(dest+1))[0] = color>>8;  return(0); }
  clr = get_color(data, PREVIEWCOLOR);
  rr = (float)(clr>>16)/255;
  gg = (float)((clr>>8)&0xFF)/255;
  bb = (float)(clr&0xFF)/255;
  for (k=0, r=g=b=0, o=1; k<num; k++) {
    memcpy(x2, inter+k*6, 3*sizeof(real));
    v2[0] = -1;  v2[1] = 1;  v2[2] = 1;
    num2 = cutaway_inter(x2, v2, pos, data, inter2);
    r2 = g2 = b2 = 1;
    for (l=0; l<num2; l++) {
      r2 *= rr;  g2 *= gg;  b2 *= bb;
      cost = v2[0]*inter2[l*6+3]+v2[1]*inter2[l*6+4]+v2[2]*inter2[l*6+5];
      acost = (fabs(cost)+1)*0.5;
      r2 *= (1-oo)*acost;  g2 *= (1-oo)*acost;  b2 *= (1-oo)*acost; }
    cost = v2[0]*inter[k*6+3]+v2[1]*inter[k*6+4]+v2[2]*inter[k*6+5];
    acost = (fabs(cost)+1)*0.5;
    cost = (cost+1)*0.5;
    r2 = (r2*(1-amb)*acost+amb*cost)*oo*rr;
    g2 = (g2*(1-amb)*acost+amb*cost)*oo*gg;
    b2 = (b2*(1-amb)*acost+amb*cost)*oo*bb;
    r += r2*o;  g += g2*o;  b += b2*o;  o *= (1-oo); }
  dest[0] = 255*b;  dest[1] = 255*g;  dest[2] = 255*r;
}

cutaway_inter(real *x, real *v, real *pos, DATA *data, real *inter)
/* Calculate the intercepts anf surface normals along a sphere with a cutout.
 *  The intersection may not be at the starting point.
 * Enter: real *x: starting point of ray to use for intersection.
 *        real *v: vector to use for intersection.  This is normalized.
 *        real *pos: x,y,z,alpha,beta,gamma,radius,nosurface specifying
 *                   cutaway point.
 *        DATA *data: pointer to cutaway type and matrix.
 *        real *inter: location to store up to five intersections.  Each
 *                     intersection contains an x,y,z coordinate and a dx,dy,
 *                     dz unit surface normal.  This must be able to hold 30
 *                     values.
 * Exit:  long numinter: number of intersections.              7/29/97-DWM */
{
  real si[6], x2[3], x3[2], dist[5], den, d, *mat=data->cut.mat;
  real *xyz=data->cut.xyz;
  long i, j, num=0, num2, inv=1-2*(data->cut.invert!=0);
  #define cm0 i*3
  #define cm1 i*3+1
  #define cm2 i*3+2

  calc_unit(v);
  x2[0] = x[0]+v[0]*0.001*pos[6];
  x2[1] = x[1]+v[1]*0.001*pos[6];
  x2[2] = x[2]+v[2]*0.001*pos[6];
  if (!(num=sphere_inter(x2, v, pos[6], si)))  return(0);
  for (i=0; i<num; i++)
    if (cut_zone(si+3*i, &data->cut, 0)) {
      if (!i)  memcpy(si, si+3, 3*sizeof(real));
      i--;  num--; }
  if (num>0) {
    memcpy(inter,   si,   3*sizeof(real));
    memcpy(inter+3, si,   3*sizeof(real)); }
  if (num>1) {
    memcpy(inter+6, si+3, 3*sizeof(real));
    memcpy(inter+9, si+3, 3*sizeof(real)); }
  num2 = num;
  if (!pos[7]) for (i=2; i>2-data->cut.type; i--) {
    den = v[0]*mat[cm0]+v[1]*mat[cm1]+v[2]*mat[cm2];
    if (!den)  continue;
    d = (mat[cm0]*(xyz[0]-x2[0])+mat[cm1]*(xyz[1]-x2[1])+
         mat[cm2]*(xyz[2]-x2[2]))/den;
    if (d<0)  continue;
    inter[num*6]   = x2[0]+v[0]*d;
    inter[num*6+1] = x2[1]+v[1]*d;
    inter[num*6+2] = x2[2]+v[2]*d;
    if (data->cut.type>=2) {
      for (j=0; j<3; j++)
        x3[j] = mat[j*3]*(inter[num*6]-xyz[0])+mat[j*3+1]*(inter[num*6+1]-
                xyz[1])+mat[j*3+2]*(inter[num*6+2]-xyz[2]);
      if (data->cut.type==2) {
        if (i==2 && x3[1]<0)  continue;
        if (i==1 && x3[2]<0)  continue; }
      if (data->cut.type==3) {
        if (i==2 && (x3[0]<0 || x3[1]<0))  continue;
        if (i==1 && (x3[0]<0 || x3[2]<0))  continue;
        if (i==0 && (x3[1]<0 || x3[2]<0))  continue; } }
    if (sq(inter[num*6])+sq(inter[num*6+1])+sq(inter[num*6+2])>sq(pos[6]))
      continue;
    inter[num*6+3] = mat[i]  *inv;
    inter[num*6+4] = mat[i+3]*inv;
    inter[num*6+5] = mat[i+6]*inv; num++; }
  if (num!=num2 && num>=2) {
    for (i=0; i<num; i++)
      dist[i] = sq(inter[i*6]-x2[0])+sq(inter[i*6+1]-x2[1])+
                sq(inter[i*6+2]-x2[2]);
    for (i=num-2; i<num-1; i++)
      for (j=i+1; j<num; j++)
        if (dist[i]>dist[j]) {
          memcpy(si, inter+j*6, 6*sizeof(real));
          memcpy(inter+j*6, inter+i*6, 6*sizeof(real));
          memcpy(inter+i*6, si, 6*sizeof(real));
          d = dist[i];  dist[i] = dist[j];  dist[j] = d; } }
  for (i=0; i<num; i++)
    calc_unit(inter+i*6+3);
  return(num);
}

cutaway_move(HWND hdlg, long delx, long dely, long w, long h, long down,
             DATA *data, real *pos)
/* Actually move the cutaway, updating the dialog elements..
 * Enter: HWND hdlg: handle of current dialog.
 *        long delx, dely: amount of mouse movement.
 *        long w, h: preview image size.
 *        long down: 1-rotate, 2-zoom, 3-shift.
 *        DATA *data: pointer to current window.
 *        real *pos: x,y,z,alpha,beta,gamma,radius,nosurface specifying
 *                   cutaway point.
 * Exit:  long moved: 0 for no change in object position, 1 for change.
 *                                                              1/7/98-DWM */
{
  long i, ru, du;
  double m[27], ang, ang1[3], pos1[3], pos2[3], x;
  lmatrix r, dr, rp;
  char mv[]={4,8,5,7, 8,0,6,2};
  char text[80], text2[80];

  if (down==1) {
    if (!delx && !dely)  return(0);
    if (sizeof(real)!=sizeof(double))  return(0);
    pos1[0] = pos[0];  pos1[1] = pos[1];  pos1[2] = pos[2];
    r.m = m;  dr.m = m+9;  rp.m = m+18;  r.w = r.h = 3;
    labgtomat(r.m, pos[3], pos[4], pos[5]);
    for (i=0; i<2; i++) {
      if (i)  ang = PI*delx/w;
      else    ang = PI*dely/h;
      if (!ang)  continue;
      lmatident(&dr, 3);
      dr.m[mv[i*4]] = dr.m[mv[i*4+1]] = cosl(ang);
      dr.m[mv[i*4+2]] = sinl(ang);
      dr.m[mv[i*4+3]] = -sinl(ang);
      pos2[0] = pos1[0]*dr.m[0]+pos1[1]*dr.m[3]+pos1[2]*dr.m[6];
      pos2[1] = pos1[0]*dr.m[1]+pos1[1]*dr.m[4]+pos1[2]*dr.m[7];
      pos2[2] = pos1[0]*dr.m[2]+pos1[1]*dr.m[5]+pos1[2]*dr.m[8];
      memcpy(pos1, pos2, 3*sizeof(double));
      lmatmul(&r, &dr, &rp);
      lmatdup(&rp, &r); }
    lmattoabg(r.m, ang1);
    ru = SendDlgItemMessage(hdlg, CutRadUnit, CB_GETCURSEL, 0, 0);
    for (i=0; i<3; i++) {
      ang = ang1[i]/RadVal[ru];
      if (fabs(ang)<1e-6)  ang = 0;  sprintf(text, "%g", ang);
      sprintf(text2, "%g", ang+1e-6-2e-6*(ang<0));
      if (strlen(text2)<strlen(text)) strcpy(text, text2);
      SetDlgItemText(hdlg, CutAng0+i, text); }
    du = SendDlgItemMessage(hdlg, CutPosUnit, CB_GETCURSEL, 0, 0);
    for (i=0; i<3; i++) {
      sprintf(text, "%g", pos1[i]/DistVal[du]);
      SetDlgItemText(hdlg, CutPosX+i, text); }
    return(1); }
  if (down==2) {
    if (!dely)  return(0);
    pos1[2] = pos[2]+pos[6]*2*dely/h;
    du = SendDlgItemMessage(hdlg, CutPosUnit, CB_GETCURSEL, 0, 0);
    sprintf(text, "%g", pos1[2]/DistVal[du]);
    SetDlgItemText(hdlg, CutPosZ, text);
    return(1); }
  if (down==3) {
    if (!delx && !dely)  return(0);
    pos1[0] = pos[0]+pos[6]*2*delx/w;
    pos1[1] = pos[1]-pos[6]*2*dely/h;
    du = SendDlgItemMessage(hdlg, CutPosUnit, CB_GETCURSEL, 0, 0);
    sprintf(text, "%g", pos1[0]/DistVal[du]);
    SetDlgItemText(hdlg, CutPosX, text);
    sprintf(text, "%g", pos1[1]/DistVal[du]);
    SetDlgItemText(hdlg, CutPosY, text);
    return(1); }
  return(0);
}

cutaway_rect(RECT *rect, HWND hdlg)
/* Locate the rectangle used in a cutaway preview.
 * Enter: RECT *rect: pointer to store result.  If null, the rectangle is
 *                    invalidated (refreshed).
 *        HWND hdlg: dialog the rectangle is located in.       10/7/96-DWM */
{
  long inv=0;
  RECT rect2;

  if (!rect) { rect = &rect2;  inv = 1; }
  rect->left = CutawayCoor[0];
  rect->top = CutawayCoor[1];
  rect->right = CutawayCoor[0]+CutawayCoor[2];
  rect->bottom = CutawayCoor[1]+CutawayCoor[3];
  MapDialogRect(hdlg, rect);
  if (inv)
    InvalidateRect(hdlg, rect, 0);
}

draw_bmp(HDC hdc, long x, long y, HANDLE bmp)
/* Draw a BMP format graphic.
 * Enter: HDC hdc: Windows style pointer to a window.
 *        long x, y: location to draw BMP.
 *        HANDLE bmp: handle to BMP format graphic.             5/8/96-DWM */
{
  LPBITMAPINFOHEADER header;
  LPSTR buf;
  long w, h, pal;

  header = GlobalLock(bmp);
  w = header->biWidth;  h = header->biHeight;
  pal = header->biBitCount;
  if (pal<=8)  pal = (1<<pal);  else pal = 0;
  buf = ((LPSTR)header) + header->biSize + 4*pal;
  SetDIBitsToDevice(hdc, x, y, w, h, 0, 0, 0, h, buf, (LPBITMAPINFO)header,
                    DIB_RGB_COLORS);
  GlobalUnlock(bmp);
}

draw_bmp2(HDC hdc, long x, long y, uchar *bmp)
/* Draw a BMP format graphic.
 * Enter: HDC hdc: Windows style pointer to a window.
 *        long x, y: location to draw BMP.
 *        uchar *bmp: pointer to BMP format graphic.            5/8/96-DWM */
{
  LPBITMAPINFOHEADER header;
  LPSTR buf;
  long w, h, pal;

  header = bmp;
  w = header->biWidth;  h = header->biHeight;
  pal = header->biBitCount;
  if (pal<=8)  pal = (1<<pal);  else pal = 0;
  buf = ((LPSTR)header) + header->biSize + 4*pal;
  SetDIBitsToDevice(hdc, x, y, w, h, 0, 0, 0, h, buf, (LPBITMAPINFO)header,
                    DIB_RGB_COLORS);
}

fill_norm(DATA *data)
/* Fill the arrays which contain the surface normals.  These are actual
 *  surface normals, not screen surface normals.
 * Enter: DATA *data: pointer to window's data area.          12/29/97-DWM */
{
  float *node, **normptr, **enormptr, *norm, *enorm;
  double x, y, z, V1[3], V2[3], V3[3];
  long i, j, *elem, num, n, e, mode;

  mode = (data->dflag>>7)&3;
  if (mode==2 || (!mode && !data->asym.opacity))  return(0);
  for (i=0; i<2; i++) {
    if (!i) {
      if (!mode)  continue;
      normptr = &data->poly.norm;  enormptr = &data->poly.enorm;
      n = data->poly.n;  e = data->poly.e;
      node = data->poly.xyz;  elem = data->poly.elem; }
    else {
      if (!data->asym.opacity)  continue;
      normptr = &data->asym.norm;  enormptr = &data->asym.enorm;
      n = data->asym.n;  e = data->asym.e;
      node = data->asym.xyz;  elem = data->asym.elem; }
    if (!e)  continue;
    if (normptr[0] && enormptr[0])  continue;
    if (normptr[0]) {   free2(normptr[0]);   normptr[0] = 0; }
    if (enormptr[0]) {  free2(enormptr[0]);  enormptr[0] = 0; }
    normptr[0]  = malloc2(n*3*sizeof(float));
    enormptr[0] = malloc2(e*3*sizeof(float));
    if (!normptr[0] || !enormptr[0]) {
      if (normptr[0])  free2(normptr[0]);   normptr[0]  = 0;
      if (enormptr[0]) free2(enormptr[0]);  enormptr[0] = 0;  continue; }
    norm = normptr[0];  enorm = enormptr[0];
    memset(norm, 0, n*3*sizeof(float));
    for (i=0; i<e; i++) {
      V1[0] = node[elem[i*3+1]*3]   - node[elem[i*3]*3];
      V1[1] = node[elem[i*3+1]*3+1] - node[elem[i*3]*3+1];
      V1[2] = node[elem[i*3+1]*3+2] - node[elem[i*3]*3+2];
      V2[0] = node[elem[i*3+2]*3]   - node[elem[i*3]*3];
      V2[1] = node[elem[i*3+2]*3+1] - node[elem[i*3]*3+1];
      V2[2] = node[elem[i*3+2]*3+2] - node[elem[i*3]*3+2];
      lunit(lcross_product(V1, V2, V3), V3);
      enorm[i*3] = V3[0];  enorm[i*3+1] = V3[1];  enorm[i*3+2] = V3[2];
      if (!V3[0] && !V3[1] && !V3[2])  continue;
      for (j=0; j<3; j++) {
        norm[elem[i*3+j]*3]   += enorm[i*3]+100;
        norm[elem[i*3+j]*3+1] += enorm[i*3+1];
        norm[elem[i*3+j]*3+2] += enorm[i*3+2]; } }
    for (i=0; i<n; i++) {
      num = floor((norm[i*3]+50)/100);
      if (num<=0)  continue;
      norm[i*3] = (norm[i*3]-100*num)/num;
      norm[i*3+1] /= num;
      norm[i*3+2] /= num; } }
}

fill_zone24(uchar *dest, long bclr, long count)
/* Fill a memory area with a 24-bit value.  This is useful for "erasing" a
 *  memory array to a non-grey color.
 * Enter: char *dest: location to store results.
 *        long bclr: 24-bit value to fill memory with.
 *        long count: number of 24-bit values to store.  Note that three
 *                    times this many bytes are written.        6/9/96-DWM */
{
  _asm {
          mov edi, dest
          mov eax, bclr
          mov bx, ax
          shr eax, 0x10
          mov ecx, count
fill1:    mov es:[edi], bx
          mov es:[edi+2], al
          add edi, 0x03
          loop    fill1
    }
}

focal_geo(long axis, float jump)
/* Shift the camera's internal parameters for the geometry associated with
 *  the top window.
 * Enter: long axis: 0-x0, 1-y0, 2-f.
 *        float jump: amount to change axis.                   3/31/97-DWM */
{
  DATA *data;
  float dx0[4]={0,0,0,0};
  long i;
  HWND hwnd;

  hwnd = SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0);
  data = lock_window(hwnd);
  if (!data)  return(0);
  for (i=0; i<18; i++)  if (data->renderdlt[i]) break;
  if (i==18) {
    unlock_window(data);  return(0); }
  dx0[axis] = jump;
  render_move(0, 0, dx0, 0, 0, data->renderval, data->renderdlt,
              data->renderdlt);
  update_process(data);
  unlock_window(data);
  InvalidateRect(hwnd, 0, 0);
}

light_figure(HWND hdlg, HDC hdc, LIGHT *ls)
/* Create the figure showing the current state of the light.
 * Enter: HWND hdlg: handle to dialog to work with.
 *        HDC hdc: device context handle to draw with.
 *        LIGHT *ls: light source to draw, or null to blank screen.
 *                                                            12/27/97-DWM */
{
  long w, h, backclr, *lpic, s, i, lastpart=8;
  RECT rect;
  HBITMAP wpic, wpic2;
  uchar *pic2;
  LOGPALETTE *lpal;
  long tim=clock();

  light_rect(&rect, hdlg);
  w = rect.right-rect.left;  h = rect.bottom-rect.top;
  if (h>w)  h = w;  if (w>h)  w = h;
  backclr = GetSysColor(COLOR_BTNFACE);
  s = ((w*3+3)/4)*4;
  wpic = GlobalAlloc(GMEM_MOVEABLE, s*h+40);
  if (!wpic)  return(0);
  pic2 = GlobalLock(wpic);  lpic = (long *)pic2;
  lpic[0] = 40;  lpic[1] = w;  lpic[2] = h;
  ((short *)pic2)[6] = 1;  ((short *)pic2)[7] = 24;
  lpic[4] = lpic[6] = lpic[7] = 0;  lpic[5] = s*h;
  lpic[8] = lpic[9] = 0;
  light_figure_draw(w, h, s, pic2+40, ls, backclr);
  GlobalUnlock(wpic);
  wpic = reduce_color_space(wpic, hdc);
  draw_bmp(hdc, rect.left, rect.top, wpic);
  GlobalFree(wpic);
  tim = 1000*(tim-clock())/CLOCKS_PER_SEC;
  if (tim<100)       PreviewDelay = 100;
  else if (tim<200)  PreviewDelay = 500;
  else if (tim>1000) PreviewDelay = 2000;
  else               PreviewDelay = 1000;
}

light_figure_draw(long w, long h, long scan, uchar *buf, LIGHT *ls,
                  long color)
/* Draw the light figure to a buffer.
 * Enter: long w, h: size of figure to draw.
 *        long scan: bytes per scan line.
 *        uchar *buf: pointer to first scan line.
 *        LIGHT *ls: pointer to light source to draw, or null to just blank
 *                   screen.
 *        long type: 0-none, 1-plane, 2-corner, 3-wedge.
 *        long color: background color.                        7/29/97-DWM */
{
  long clr, br, bg, bb, fr, fg, fb;
  long i, j, k, l, d, total[3];
  float r, r2, w2=w*0.5, h2=h*0.5;
  real lv[3], norm[3], b;

  clr = get_color(0, PREVIEWCOLOR);
  br = color&0xFF;  bg = (color>>8)&0xFF;  bb = color>>16;
  fr = clr  &0xFF;  fg = (clr  >>8)&0xFF;  fb = clr  >>16;
  if (ls) {
    lv[0] = ls->x[0];  lv[1] = ls->x[1];  lv[2] = ls->x[2];
    calc_unit(lv); }
  r = min(w*0.5, h*0.5);  r2 = sq(r);
  for (j=0; j<h; j++) {
    d = j*scan;
    for (i=0; i<w; i++, d+=3) {
      total[0] = total[1] = total[2] = 0;
      for (k=0; k<2; k++)  for (l=0; l<2; l++)
        if (sq(i+k*0.5-w2)+sq(j+l*0.5-h2)>r2 || !ls) {
          total[0] += br;  total[1] += bg;  total[2] += bb; }
        else {
          norm[0] = (i+k*0.5-w2)/r;  norm[1] = (j+l*0.5-h2)/r;
          norm[2] = sqrtl(fabs(1-sq(norm[0])-sq(norm[1])));
          b = (norm[0]*lv[0]+norm[1]*lv[1]+norm[2]*lv[2])*0.5+0.5;
          if (b>=0.5)  b *= ls->i;
          else         b *= ls->i*ls->a;
          if (b<0)  b = 0;  if (b>1)  b = 1;
          total[0] += fr*b;  total[1] += fg*b;  total[2] += fb*b; }
      buf[d]   = total[0]/4;
      buf[d+1] = total[1]/4;
      buf[d+2] = total[2]/4; } }
}

light_move(HWND hdlg, long delx, long dely, long w, long h, long down)
/* Actually move the light source, modifying the dialog text fields.
 * Enter: HWND hdlg: handle of current dialog.
 *        long delx, dely: amount of mouse movement.
 *        long w, h: preview image size.
 *        long down: 1-rotate, 2-zoom, 3-shift.
 * Exit:  long moved: 0 for no change in object position, 1 for change.
 *                                                              1/7/98-DWM */
{
  real val[3], fac, x, y, z, st, ct;
  char text[80];
  long i;

  if ((!dely && down==2) || (!dely && !delx) || down==3)  return(0);
  for (i=0; i<3; i++) {
    val[i] = 0;
    GetDlgItemText(hdlg, LightEditX+i, text, 79);
    sscanf(text, "%lf", val+i); }
  if (down==2 && dely) {
    if (dely>0)  fac = pow(100, (float)dely/h);
    else         fac = 1./pow(100, (float)fabs(dely)/h);
    val[0] *= fac;  val[1] *= fac;  val[2] *= fac; }
  if (down==1 && delx) {
    fac = PI*delx/w;  st = sinl(fac);  ct = cosl(fac);
    z = val[2]*ct-val[0]*st;
    x = val[2]*st+val[0]*ct;
    val[0] = x;  val[2] = z; }
  if (down==1 && dely) {
    fac = PI*dely/h;  st = sinl(fac);  ct = cosl(fac);
    y = val[1]*ct-val[2]*st;
    z = val[1]*st+val[2]*ct;
    val[1] = y;  val[2] = z; }
  for (i=0; i<3; i++) {
    sprintf(text, "%g", val[i]);
    SetDlgItemText(hdlg, LightEditX+i, text); }
  return(1);
}

light_rect(RECT *rect, HWND hdlg)
/* Locate the rectangle used in a light preview.
 * Enter: RECT *rect: pointer to store result.  If null, the rectangle is
 *                    invalidated (refreshed).
 *        HWND hdlg: dialog the rectangle is located in.       10/7/96-DWM */
{
  long inv=0;
  RECT rect2;

  if (!rect) { rect = &rect2;  inv = 1; }
  rect->left = LightCoor[0];
  rect->top = LightCoor[1];
  rect->right = LightCoor[0]+LightCoor[2];
  rect->bottom = LightCoor[1]+LightCoor[3];
  MapDialogRect(hdlg, rect);
  if (inv)
    InvalidateRect(hdlg, rect, 0);
}

make_palette()
/* Set up a palette for use in preview mode if this is an 8-bit display.
 *  This also initializes the dither tables.                   3/11/97-DWM */
{
  LOGPALETTE *lpal;
  uchar *pal2, *pal;
  long i, j;

  pal = smooth_palette(0, SMOOTHLOW);
  use_palette(pal);
  lpal = LocalAlloc(LPTR, sizeof(LOGPALETTE) + 256*sizeof(PALETTEENTRY));
  if (!lpal) { free2(pal);  return(0); }
  lpal->palVersion = 0x300;  lpal->palNumEntries = 256;
  for (i=0; i<256; i++) {
    lpal->palPalEntry[i].peRed   = pal[i*3];
    lpal->palPalEntry[i].peGreen = pal[i*3+1];
    lpal->palPalEntry[i].peBlue  = pal[i*3+2];
    lpal->palPalEntry[i].peFlags = 0; }
  if (Hpal) { DeleteObject(Hpal);  Hpal = 0; }
  Hpal = CreatePalette(lpal);
  LocalFree(lpal);
  free2(pal);
  if (dithernum<0 && Dither) {
    for (i=0; i<2048; i++)  dithertbl1[i] = dithertbl3[rand()%16]+22;
    for (i=0; i<256+22*2; i++) {
      j = i-22;  if (j<0) j = 0;  if (j>255) j =255;
      dithertbl2[i]     = (j&0xF8)<<8;
      dithertbl2[i+300] = (j&0xFC)<<3;
      dithertbl2[i+600] =  j      >>3; }
    dithernum = 0; }
}

HPALETTE make_windows_palette(uchar *pal, long bgr)
/* Make a standard palette into a windows palette.
 * Enter: uchar *pal: 768 bytes containg a palette.
 *        long bgr: 0 for RGB, 1 for bgr.
 * Exit:  HPALETTE wpal: globally alloced windows palette.  Null for error.
 *                                                             7/26/97-DWM */
{
  LOGPALETTE *lpal;
  long b=2, r=0, i;
  HPALETTE wpal;

  if (bgr) { b = 0;  r = 2; }
  lpal = LocalAlloc(LPTR, sizeof(LOGPALETTE) + 256*sizeof(PALETTEENTRY));
  if (!lpal)  return(0);
  lpal->palVersion = 0x300;  lpal->palNumEntries = 256;
  for (i=0; i<256; i++) {
    lpal->palPalEntry[i].peRed   = pal[i*3+r];
    lpal->palPalEntry[i].peGreen = pal[i*3+1];
    lpal->palPalEntry[i].peBlue  = pal[i*3+b];
    lpal->palPalEntry[i].peFlags = 0; }
  wpal = CreatePalette(lpal);
  LocalFree(lpal);
  return(wpal);
}

mouse(HWND hwnd, long but, long up, long flag, long x, long y)
/* Handle mouse button pushes and free mouse movement.
 * Enter: HWND hwnd: handle of window this is for.
 *        long but: button number.  -1 for movement, 0 for left, 1 for right,
 *                  2 for middle.
 *        long up: 0 for down, 1 for up (undefined for movement).
 *        long flag: MK_CONTROL, MK_LBUTTON, MK_MBUTTON, MK_RBUTTON, and
 *                   MK_SHIFT flags.
 *        long x, y: position of the mouse within the window.   4/2/97-DWM */
{
  static long cur=-1;
  long newcur=0, i;
  RECT rect;
  DATA *data;

  if (!hwnd)  return(0);
  if (x>32768)  x -= 65536;  if (y>32768)  y -= 65536;
  if (but>=0 || Down)  LastMove = 0;
  if (Down && (but<0 || up)) {
    mouse_move(hwnd, up, flag, x, y);  return(0); }
  if (but>=0 && up) return(0);
  data = lock_window(hwnd);
  if (!data)  return(0);
  DownX = x;  DownY = y;
  zone_geo(hwnd, data, x, y, DownBuf);
  if ((!but && (flag&MK_RBUTTON)) || (but==1 && (flag&MK_LBUTTON)))
    but = 2;
  if (!but)        { Down = 1;  newcur = 3;  SetCapture(hwnd); }
  else if (but==1) { Down = 2;  newcur = 4;  SetCapture(hwnd); }
  else if (but==2) { Down = 3;  newcur = 2;  SetCapture(hwnd); }
  unlock_window(data);
  if (newcur!=cur)
    cursor(cur=newcur);
}

mouse_move(HWND hwnd, long up, long flag, long x, long y)
/* Handle mouse button releases and constrained mouse movement.  The current
 *  mouse process is defined by the Down value.  This is:
 *   1 - rotate geo.
 *   2 - shift geo.
 *   3 - zoom geo.
 *  See zone_geo for a description of the values in DownBuf[].
 * Enter: HWND hwnd: handle of window this is for.
 *        long up: 0 for down or movement, 1 for up.
 *        long flag: MK_CONTROL, MK_LBUTTON, MK_MBUTTON, MK_RBUTTON, and
 *                   MK_SHIFT flags.
 *        long x, y: position of the mouse within the window.   4/2/97-DWM */
{
  switch (Down) {
    case 1: if (x!=DownX && DownBuf[5])
        rotate_geo(1, 180.*(DownX-x)/DownBuf[5]);
      if (y!=DownY && DownBuf[6])
        rotate_geo(0, 180.*(DownY-y)/DownBuf[6]);
      DownX = x;  DownY = y; break;
    case 2: if (y!=DownY && DownBuf[6])
      shift_geo(2, 20.*(y-DownY)/DownBuf[6]);
      DownY = y; break;
    case 3: if (x!=DownX && DownBuf[5])
        shift_geo(0, (float)(x-DownX)/DownBuf[5]);
      if (y!=DownY && DownBuf[6])
        shift_geo(1, (float)(y-DownY)/DownBuf[6]);
      DownX = x;  DownY = y; break;
    /** Additional mouse movements go here **/
    default: ; }
  if (up) {
    Down = 0;
    ReleaseCapture(); }
}

opacity_dither()
/* Initialize the opacity dithering table.                    12/28/97-DWM */
{
  long i;

  for (i=0; i<OPACDITHER; i++)
    opacitytbl[i] = rnd(0, 16776959, 3);
}

preview_mouse(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp, DATA *data,
              long num, long ptr)
/* Handle mouse movement in the Camera, Cutaway, or Light Options dialog.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.
 *        DATA *data: pointer to owner window.
 *        long num: 0-camera, 1-cutaway, 2-light.
 *        long ptr: this can be a pointer to a relevant data array, or any
 *                  other long value to get passed on to the actual
 *                  movement routines.
 * Exit:  long moved: 0 for no change in object position, 1 for change.
 *                                                              1/6/98-DWM */
{
  static long down=0, dx, dy;
  long *box, x, y, out[4], but, up, cur=-1, newcur=0, ret=0;

  switch (num) {
    case 0: box = CameraCoor; break;
    case 1: box = CutawayCoor; break;
    case 2: box = LightCoor; }
  x = lp&0xFFFF;  y = lp>>16;
  if (x>32768)  x -= 65536;  if (y>32768)  y -= 65536;
  map_dialog_rect(hdlg, box[0], box[1], box[2], box[3], out);
  if (!down)
    if (x<out[0] || x>=out[0]+out[2] || y<out[1] || y>=out[1]+out[3])
      return(0);
  up = (msg==WM_LBUTTONUP || msg==WM_MBUTTONUP || msg==WM_RBUTTONUP);
  but = -1;
  if (msg==WM_LBUTTONUP || msg==WM_LBUTTONDOWN)  but = 0;
  if (msg==WM_RBUTTONUP || msg==WM_RBUTTONDOWN)  but = 1;
  if (msg==WM_MBUTTONUP || msg==WM_MBUTTONDOWN)  but = 2;
  if (down && (but<0 || up)) {
    if (dx!=x || dy!=y)
      switch (num) {
        case 0: ret = camera_move(hdlg, x-dx,y-dy, out[2],out[3], down, data,
                                  ptr);  break;
        case 1: ret = cutaway_move(hdlg, x-dx,y-dy, out[2],out[3], down,
                                   data, (real *)ptr);  break;
        case 2: ret = light_move(hdlg, x-dx,y-dy, out[2],out[3], down); }
    if (up) {
      down = 0;  Busy &= 0xFFFFFFEF;
      ReleaseCapture();
      cursor(cur=0); }
    dx = x;  dy = y;  return(ret); }
  if (but>=0 && up) return(down=0);
  dx = x;  dy = y;
  if ((!but && (wp&MK_RBUTTON)) || (but==1 && (wp&MK_LBUTTON)))
    but = 2;
  if (num==2 && but>=2)  return(0);
  if (!but)        { down = 1;  newcur = 3;  SetCapture(hdlg); }
  else if (but==1) { down = 2;  newcur = 4;  SetCapture(hdlg); }
  else if (but==2) { down = 3;  newcur = 2;  SetCapture(hdlg); }
  if (down)  Busy |= 16;
  if (newcur!=cur)
    cursor(cur=newcur);
  return(0);
}

HBITMAP reduce_color_space(HBITMAP wpic, HDC hdc)
/* Convert a 24-bit per pixel image into a palettize image, if appropriate.
 *  Also dithers images for 16-bit display modes.
 * Enter: HBITMAP wpic: handle to bitmap to reduce.
 *        HDC hdc: device context for palette.
 * Exit:  HBITMAP wpic: handle of reduced bitmap.              1/24/98-DWM */
{
  HBITMAP wpic2;
  LOGPALETTE *lpal;
  uchar *pic;
  long w, h, i;

  if (BitsPixel<=8) {
    wpic2 = PalettizeBMP(wpic, 1, 0);
    if (wpic2) {
      GlobalFree(wpic);  wpic = wpic2;
      lpal = BMPPalette(wpic);
      if (lpal) {
        if (Hpal) { DeleteObject(Hpal);  Hpal = 0; }
        Hpal = CreatePalette(lpal);
        LocalFree(lpal); }
      if (Hpal) {
        SelectPalette(hdc, Hpal, 0);
        RealizePalette(hdc); } } }
  if (BitsPixel>8 && BitsPixel<=16) {
    pic = lock2(wpic);
    if (dither16num<0) {
      for (i=0; i<2048; i++)
        dither16tbl[i] = rand()&7;
      dither16num = 0; }
    w = ((long *)pic)[1];  h = ((long *)pic)[2];
    for (i=40; i<w*h*3+40; i+=3) {
      pic[i] = min(255,pic[i]+dither16tbl[dither16num]);
      pic[i+1] = min(255,pic[i+1]+(dither16tbl[dither16num+1]&3));
      pic[i+2] = min(255,pic[i+2]+dither16tbl[dither16num+2]);
      dither16num = (dither16num+3)%2046; }
    unlock2(pic); }
  return(wpic);
}

reframe_geo()
/* Return the geometry to the default position.                 4/1/97-DWM */
{
  DATA *data;
  long i;
  HWND hwnd;
  float dx0[4]={0,0,0,0};

  hwnd = SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0);
  data = lock_window(hwnd);
  if (!data)  return(0);
  if (!data->mol.maxcheck) {
    unlock_window(data);
    error("Atom has not been computed yet.");  return(0); }
  DefSize = data->mol.maxcheck;
  dx0[3] = 2*data->mol.maxcheck;
  render_move(0, 0, dx0, 0, 0, data->renderval, data->renderdlt,
              data->renderdlt);
  update_process(data);
  unlock_window(data);
  InvalidateRect(hwnd, 0, 0);
}

render_asymptote(long w, long h, long scan, uchar *scrbuf, ushort *zbuf,
                 DATA *data)
/* Render the asymptote on an orbital.  This can either be a clear wireframe
 *  or a solid/translucent lighted surface.
 * Enter: long w, h: size of output screen in pixels.
 *        long scan: number of bytes per scan line on output screen.  Must be
 *                   at least 3*w.  This is negative for bgr organized data.
 *        uchar *scrbuf: buffer pointing to the top scanline on the screen.
 *        ushort *zbuf: zbuffer.  Lower is 'on top'.
 *        DATA *data: pointer to window data to draw.  This includes screen
 *                    coordinates of the nodes and viewing flags.
 *                                                            12/28/97-DWM */
{
  long i, j, k, color, *elem, *node, e, n, n2, clr, numl, coor[15];
  float *norm, val;
  LIGHT *ls;

  if (!data->asym.opacity || !data->asym.e || !data->asym.scrxyz)  return(0);
  node = data->asym.scrxyz;  elem = data->asym.elem;
  e = data->asym.e*3;  n = data->asym.n;  n2 = n+n;
  color = get_color(data, ASYMCOLOR);
  if (data->asym.wire==1) {
    for (i=0; i<e; i+=3) {
      render_line(w, h, scan, scrbuf, w, zbuf, node[elem[i]],
                  node[elem[i]+n], node[elem[i]+n2], node[elem[i+1]],
                  node[elem[i+1]+n], node[elem[i+1]+n2], color);
      render_line(w, h, scan, scrbuf, w, zbuf, node[elem[i]],
                  node[elem[i]+n], node[elem[i]+n2], node[elem[i+2]],
                  node[elem[i+2]+n], node[elem[i+2]+n2], color);
      render_line(w, h, scan, scrbuf, w, zbuf, node[elem[i+1]],
                  node[elem[i+1]+n], node[elem[i+1]+n2], node[elem[i+2]],
                  node[elem[i+2]+n], node[elem[i+2]+n2], color); }
    return(0); }
  if (data->asym.wire==2) {
    for (i=0; i<e; i++)
      render_point(w, h, scan, scrbuf, zbuf, node[elem[i]], node[elem[i]+n],
                   node[elem[i]+n2], Pref.pointsize, color, 1);
    return(0); }
  norm = data->asym.norm;
  ls = data->mol.ls;  numl = data->mol.numl;
  for (i=0; i<e && norm; i+=3) {
    for (j=0; j<3; j++) {
      coor[j*5]   = node[elem[i+j]];
      coor[j*5+1] = node[elem[i+j]+n];
      coor[j*5+2] = node[elem[i+j]+n2];
      if (!numl)
        coor[j*5+3] = 65535;
      else
        for (k=0, coor[j*5+3]=0; k<numl; k++)
          coor[j*5+3] += (0.5+0.5*fabs(norm[elem[i+j]*3]*ls[k].ux[0]+
                         norm[elem[i+j]*3+1]*ls[k].ux[1]+norm[elem[i+j]*3+2]*
                         ls[k].ux[2]))*ls[k].i*65535;
      coor[j*5+4] = 65535*data->asym.opacity; }
    render_triangle_opac(w, h, scan, scrbuf, w, zbuf, coor, &color); }
}

render_coor(DATA *data, double *dlt, long **scrxyz, float *node,
            long numnode, long **scrxyz2, float *node2, long numnode2)
/* Convert a set of coordinates into screen coordinates, allocating or
 *  reallocating the scrxyz array as necessary.  This also determines the
 *  actual position of all light sources.
 * Enter: DATA *data: pointer to winodw data.  This is used to add for
 *                    light source adjustment.
 *        double *dlt: set of 11 dlt parameters to use.
 *        long **scrxyz: pointer to location of scrxyz pointer.
 *        float *node: pointer to an array of xyz nodes.
 *        long numnode: number of nodes.
 *        long **scrxyz2: pointer to location of second scrxyz pointer, null
 *                        for only one set of nodes.
 *        float *node2: pointer to an array of second xyz nodes.
 *        long numnode2: number of second nodes.                7/7/97-DWM */
{
  float *fz, *fz2, minz=1e30, maxz=-1e30, den, cz;
  double L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, temp[10], tr;
  long *newx, *newy, *newz, i, j, *newx2, *newy2, *newz2;
  double mc[9];
  lmatrix m;

  if (!scrxyz || !node)  numnode = 0;
  if (numnode) {
    if (!scrxyz[0])
      newx = malloc2(numnode*sizeof(long)*3);
    else
      newx = realloc2(scrxyz[0], numnode*sizeof(long)*3);
    if (!newx)  return(0);
    scrxyz[0] = newx;  newy = newx+numnode;  newz = newy+numnode; }
  if (!scrxyz2 || !node2)  numnode2 = 0;
  if (numnode2) {
    if (!scrxyz2[0])
      newx2 = malloc2(numnode2*sizeof(long)*3);
    else
      newx2 = realloc2(scrxyz2[0], numnode2*sizeof(long)*3);
    if (!newx2)  numnode2 = 0;
    scrxyz2[0] = newx2;  newy2 = newx2+numnode2;  newz2 = newy2+numnode2; }
  if (numnode)   fz  = (float *)newz;
  if (numnode2)  fz2 = (float *)newz2;
  m.m = mc;  ldlt_to_physical(dlt, temp, &m, 0);
  tr = mc[6]*temp[0]+mc[7]*temp[1]+mc[8]*temp[2];
  L1 = dlt[0];  L2 = dlt[1];  L3 = dlt[2];  L4 = dlt[3];  L5 = dlt[4];
  L6 = dlt[5];  L7 = dlt[6];  L8 = dlt[7];  L9 = dlt[8];  L10 = dlt[9];
  L11 = dlt[10];
  for (i=j=0; i<numnode; i++,j+=3) {
    den = L9*node[j]+L10*node[j+1]+L11*node[j+2]+1;
    if (den) {
      fz[i] = mc[6]*node[j]+mc[7]*node[j+1]+mc[8]*node[j+2]-tr;
      if (fz[i]<minz)  minz = fz[i];  if (fz[i]>maxz)  maxz = fz[i];
      den = 1./den;
      newx[i] = (L1*node[j]+L2*node[j+1]+L3*node[j+2]+L4)*den;
      newy[i] = (L5*node[j]+L6*node[j+1]+L7*node[j+2]+L8)*den; } }
  for (i=j=0; i<numnode2; i++,j+=3) {
    den = L9*node2[j]+L10*node2[j+1]+L11*node2[j+2]+1;
    if (den) {
      fz2[i] = mc[6]*node2[j]+mc[7]*node2[j+1]+mc[8]*node2[j+2]-tr;
      if (fz2[i]<minz)  minz = fz2[i];  if (fz2[i]>maxz)  maxz = fz2[i];
      den = 1./den;
      newx2[i] = (L1*node2[j]+L2*node2[j+1]+L3*node2[j+2]+L4)*den;
      newy2[i] = (L5*node2[j]+L6*node2[j+1]+L7*node2[j+2]+L8)*den; } }
  if (maxz==minz)  maxz = minz+1;
  maxz = 65335./(maxz-minz);
  for (i=0; i<numnode; i++)
    newz[i] = (fz[i]-minz)*maxz+100;
  for (i=0; i<numnode2; i++)
    newz2[i] = (fz2[i]-minz)*maxz+100;
  cz = (-minz*maxz+100+1)*256;
  if (cz<LONG_MIN)       CutZ = LONG_MIN;
  else if (cz>LONG_MAX)  CutZ = LONG_MAX;
  else                   CutZ = cz;
  for (i=0; i<data->mol.numl; i++) {
    if (!data->mol.ls[i].local) {
      data->mol.ls[i].ux[0] = data->mol.ls[i].x[0]*mc[0]-data->mol.ls[i].x[1]
                              *mc[3]+data->mol.ls[i].x[2]*mc[6];
      data->mol.ls[i].ux[1] = data->mol.ls[i].x[0]*mc[1]-data->mol.ls[i].x[1]
                              *mc[4]+data->mol.ls[i].x[2]*mc[7];
      data->mol.ls[i].ux[2] = data->mol.ls[i].x[0]*mc[2]-data->mol.ls[i].x[1]
                              *mc[5]+data->mol.ls[i].x[2]*mc[8]; }
    else
      memcpy(data->mol.ls[i].ux, data->mol.ls[i].x, 3*sizeof(real));
    calc_unit(data->mol.ls[i].ux); }
}

render_line(long w, long h, long scan, uchar *scrbuf, long zscan,
            ushort *zbuf, long x0, long y0, long z0, long x1, long y1,
            long z1, long color)
/* Draw a line in 3D.  This line is drawn at the "topmost" point, such that
 *  an equivalent surface will not overwrite it.  Note that coordinates
 *  greater than +/-32767 in x or y will have bizarre effects.
 * Enter: long w, h: size of output screen in pixels.  These values, along
 *                   with scrbuf, can be set to make a clipping rectangle.
 *        long scan: number of bytes per scan line on output screen.  Must be
 *                   at least 3*w.  This is negative for bgr organized data.
 *        uchar *scrbuf: buffer pointing to the top scanline on the screen.
 *        long zscan: number of shorts per scan line of zbuffer.  Must be at
 *                    least w.
 *        ushort *zbuf: zbuffer.  Lower is 'on top'.
 *        long x0, y0, z0: coordinates of first end point on line, in screen
 *                         coordinates.  The z value is in z buffer units,
 *                         which range from 0 to CutZ.
 *        long x1, y1, z1: coordinates of second end point on line.
 *        long color: 0xBBGGRR or 0xRRGGBB color for the line.3/18/97-DWM */
{
  long temp, dx, dy, d, zd, h16, w16, dz;
  char clrb;
  short clrrg;

  dx = x1-x0;  dy = y1-y0;  dz = z1-z0;  h16 = h<<16;  w16 = w<<16;
  clrb = (color&0xFF);  clrrg = color>>8;
  if (abs(dx)>10000 || abs(dy)>10000)  return(0);
  if (abs(dx)>abs(dy)) {
    if (x0>x1) {
      temp = x0;  x0 = x1;  x1 = temp;
      temp = y0;  y0 = y1;  y1 = temp;
      temp = z0;  z0 = z1;  z1 = temp; }
    y0 = (y0<<16)+0x8000;  z0 = (z0<<8);
    y1 = (y1<<16)+0x8000;  z1 = (z1<<8);
    dy = (dy<<16)/dx;  dz = (dz<<8)/dx;
    if (x0<0) {
      y0 += -x0*dy;  z0 += -x0*dz;  x0 = 0; }
    if (x1>=w)  x1 = w-1;
    z0 -= abs(dz)-1;
    for (; x0<=x1; x0++) {
      if (y0>=0 && y0<h16) {
        zd = (y0>>16)*zscan+x0;
        if ((z0>>8)<=zbuf[zd] && z0>CutZ) {
          zbuf[zd] = z0>>8;
          d = (y0>>16)*scan+x0*3;
          scrbuf[d] = clrb;
          ((short *)(scrbuf+1+d))[0] = clrrg; } }
      y0 += dy;  z0 += dz; } }
  else {
    if (y0>y1) {
      temp = x0;  x0 = x1;  x1 = temp;
      temp = y0;  y0 = y1;  y1 = temp;
      temp = z0;  z0 = z1;  z1 = temp; }
    if (y0==y1) {
      dy = 1;  z0 = z1 = min(z0,z1); }
    x0 = (x0<<16)+0x8000;  z0 = (z0<<8);
    x1 = (x1<<16)+0x8000;  z1 = (z1<<8);
    dx = (dx<<16)/dy;  dz = (dz<<8)/dy;
    if (y0<0) {
      x0 += -y0*dx;  z0 += -y0*dz;  y0 = 0; }
    if (y1>=h)  y1 = h-1;
    z0 -= abs(dz)-1;
    for (; y0<=y1; y0++) {
      if (x0>=0 && x0<w16) {
        zd = y0*zscan+(x0>>16);
        if ((z0>>8)<=zbuf[zd] && z0>CutZ) {
          zbuf[zd] = z0>>8;
          d = y0*scan+(x0>>16)*3;
          scrbuf[d] = clrb;
          ((short *)(scrbuf+1+d))[0] = clrrg; } }
      x0 += dx;  z0 += dz; } }
}

render_point(long w, long h, long scan, uchar *scrbuf, ushort *zbuf,
             long x, long y, long z, long size, long clr, long bgr)
/* Draw a single point.
 * Enter: long w, h: size of output screen in pixels.
 *        long scan: number of bytes per scan line on output screen.  Must be
 *                   at least 3*w.  This is negative for bgr organized data.
 *        uchar *scrbuf: buffer pointing to the top scanline on the screen.
 *        ushort *zbuf: zbuffer.  Lower is 'on top'.
 *        long x, y, z: location in screen coordinates.
 *        long size: 0-1x1, 1-3x3, 2-5x5, etc.
 *        long clr: color of point.
 *        long bgr: 0 for rgb, 1 for bgr.                       1/7/98-DWM */
{
  long r, g, b, dx, dy, zd, d, rr=2*bgr, bb=2-2*bgr;

  if (x<0 || x>=w || y<0 || y>=h || z<=CutZ)  return(0);
  r = clr>>16;  g = (clr>>8)&0xFF;  b = clr&0xFF;
  for (dx=x-size; dx<=x+size; dx++)
    for (dy=y-size; dy<=y+size; dy++) {
      if (size)
        if (dx<0 || dx>=w || dy<0 || dy>=h)  continue;
      zd = dy*w+dx;
      if (z<=zbuf[zd]) {
        zbuf[zd] = z;
        d = dy*scan+dx*3;
        scrbuf[d+rr] = r;  scrbuf[d+1]  = g;  scrbuf[d+bb] = b; } }
}

render_points(long w, long h, long scan, uchar *scrbuf, ushort *zbuf,
              DATA *data, long bgr)
/* Render the nodal points on a geometry.  The color is taken from the phase.
 * Enter: long w, h: size of output screen in pixels.
 *        long scan: number of bytes per scan line on output screen.  Must be
 *                   at least 3*w.  This is negative for bgr organized data.
 *        uchar *scrbuf: buffer pointing to the top scanline on the screen.
 *        ushort *zbuf: zbuffer.  Lower is 'on top'.
 *        DATA *data: pointer to window data to draw.  This includes screen
 *                    coordinates of the nodes and viewing flags.
 *        long bgr: 0 for rgb, 1 for bgr.                      3/25/97-DWM */
{
  long pclr, nclr, *sx, *sy, *sz, clr, i;

  if (!data->points.n || !data->points.scrxyz)  return(0);
  pclr = get_color(data, POSCOLOR);  nclr = get_color(data, NEGCOLOR);
  sx = data->points.scrxyz;  sy = sx+data->points.n;  sz = sy+data->points.n;
  for (i=0; i<data->points.n; i++) {
    if (data->points.phase[i]>0)  clr = pclr;
    else                          clr = nclr;
    render_point(w, h, scan, scrbuf, zbuf, sx[i], sy[i], sz[i],
                 Pref.pointsize, clr, bgr); }
}

render_polygon(long w, long h, long scan, uchar *scrbuf, ushort *zbuf,
               DATA *data)
/* Render the polygon orbital.  This can either be a clear wireframe or a
 *  solid/translucent lighted surface.
 * Enter: long w, h: size of output screen in pixels.
 *        long scan: number of bytes per scan line on output screen.  Must be
 *                   at least 3*w.  This is negative for bgr organized data.
 *        uchar *scrbuf: buffer pointing to the top scanline on the screen.
 *        ushort *zbuf: zbuffer.  Lower is 'on top'.
 *        DATA *data: pointer to window data to draw.  This includes screen
 *                    coordinates of the nodes and viewing flags.
 *                                                            12/29/97-DWM */
{
  long i, j, k, colorp, colorn, *elem, *node, e, n, n2, clr, numl, coor[15];
  long color;
  float *norm, val;
  real opac;
  LIGHT *ls;

  if (!data->poly.e || !data->poly.scrxyz)  return(0);
  node = data->poly.scrxyz;  elem = data->poly.elem;
  e = data->poly.e*3;  n = data->poly.n;  n2 = n+n;
  colorp = get_color(data, POSCOLOR);
  colorn = get_color(data, NEGCOLOR);
  if (data->poly.wire==1) {
    for (i=0; i<e; i+=3) {
      if (data->poly.phase[i/3]>0)  color = colorp;
      else                          color = colorn;
      render_line(w, h, scan, scrbuf, w, zbuf, node[elem[i]],
                  node[elem[i]+n], node[elem[i]+n2], node[elem[i+1]],
                  node[elem[i+1]+n], node[elem[i+1]+n2], color);
      render_line(w, h, scan, scrbuf, w, zbuf, node[elem[i]],
                  node[elem[i]+n], node[elem[i]+n2], node[elem[i+2]],
                  node[elem[i+2]+n], node[elem[i+2]+n2], color);
      render_line(w, h, scan, scrbuf, w, zbuf, node[elem[i+1]],
                  node[elem[i+1]+n], node[elem[i+1]+n2], node[elem[i+2]],
                  node[elem[i+2]+n], node[elem[i+2]+n2], color); }
    return(0); }
  if (data->poly.wire==2) {
    for (i=0; i<e; i++) {
      if (data->poly.phase[i/3]>0)  color = colorp;
      else                          color = colorn;
      render_point(w, h, scan, scrbuf, zbuf, node[elem[i]], node[elem[i]+n],
                   node[elem[i]+n2], Pref.pointsize, color, 1); }
    return(0); }
  norm = data->poly.norm;
  ls = data->mol.ls;  numl = data->mol.numl;
  for (i=0; i<e; i+=3) {
    if (data->poly.phase[i/3]>0) {
      color = colorp;  opac = data->poly.opacity[0]; }
    else {
      color = colorn;  opac = data->poly.opacity[1]; }
    for (j=0; j<3; j++) {
      coor[j*5]   = node[elem[i+j]];
      coor[j*5+1] = node[elem[i+j]+n];
      coor[j*5+2] = node[elem[i+j]+n2];
      if (!numl)
        coor[j*5+3] = 65535;
      else
        for (k=0, coor[j*5+3]=0; k<numl; k++)
          coor[j*5+3] += (0.5-0.5*(norm[elem[i+j]*3]*ls[k].ux[0]+
                         norm[elem[i+j]*3+1]*ls[k].ux[1]+norm[elem[i+j]*3+2]*
                         ls[k].ux[2]))*ls[k].i*65535;
      coor[j*5+4] = 65535*opac; }
    render_triangle_opac(w, h, scan, scrbuf, w, zbuf, coor, &color); }
}

render_sort_coor(long *coor, long *x, long s, long d)
/* Sort the three coordinates so the lowest y value is first.
 * Enter: long *coor: pointer to coordinate array.
 *        long *x: pointer to output storage array.
 *        long s: number of components per source coordinate.
 *        long d: number of components per destination coordinate.
 *                                                             5/28/97-DWM */
{
  long y2=s+1, y3=s+s+1, s2=s+s, d2=d+d;

  switch (((coor[1]<coor[y2])<<2)+((coor[y2]<coor[y3])<<1)+
          (coor[1]<coor[y3])) {
    case 0: memcpy(x+d2, coor,    s*sizeof(long));
            memcpy(x+d,  coor+s,  s*sizeof(long));
            memcpy(x,    coor+s2, s*sizeof(long)); break;
    case 2: memcpy(x+d2, coor,    s*sizeof(long));
            memcpy(x,    coor+s,  s*sizeof(long));
            memcpy(x+d,  coor+s2, s*sizeof(long)); break;
    case 4: memcpy(x+d,  coor,    s*sizeof(long));
            memcpy(x+d2, coor+s,  s*sizeof(long));
            memcpy(x,    coor+s2, s*sizeof(long)); break;
    case 3: memcpy(x+d,  coor,    s*sizeof(long));
            memcpy(x,    coor+s,  s*sizeof(long));
            memcpy(x+d2, coor+s2, s*sizeof(long)); break;
    case 5: memcpy(x,    coor,    s*sizeof(long));
            memcpy(x+d2, coor+s,  s*sizeof(long));
            memcpy(x+d,  coor+s2, s*sizeof(long)); break;
    default:memcpy(x,    coor,    s*sizeof(long));
            memcpy(x+d,  coor+s,  s*sizeof(long));
            memcpy(x+d2, coor+s2, s*sizeof(long)); }
}

render_triangle_opac(long w, long h, long scan, uchar *scrbuf, long zscan,
                     ushort *zbuf, long coor[15], uchar *color)
/* Draw a triangle shaded based on a scalar value at each point.  Note that
 *  coordinates greater than +/-32767 in x or y will have bizarre effects.
 * Enter: long w, h: size of output screen in pixels.  These values, along
 *                   with scrbuf, can be set to make a clipping rectangle.
 *        long scan: number of bytes per scan line on output screen.  Must be
 *                   at least 3*w.  This is negative for bgr organized data.
 *        uchar *scrbuf: buffer pointing to the top scanline on the screen.
 *        long zscan: number of shorts per scan line of zbuffer.  Must be at
 *                    least w.
 *        ushort *zbuf: zbuffer.  Lower is 'on top'.
 *        long coor[15]: coordinates of triangle points, itensities, and
 *                       opacities.  Points are in screen coordinates, with
 *                       z values in z buffer units, which range from 0 to
 *                       CutZ.  Intensity and opacity is on a scale of 0 to
 *                       65535.  This array contains x0,y0,z0,i0,o0,x1,y1,z1,
 *                       i1,o1,x2,y2,z2,i2,o2.
 *        uchar *color: array containg either RGB or BGR maximum intensity
 *                      color.  RGB is the same sense as the destination
 *                      screen.                               12/28/97-DWM */
{
  long x[20], dx[15], i, j, k, l, y, xi, src, dest, za;
  long x1, z1, I1, o1, x2, z2, I2, o2, xi1, xi2, dz, ds, dt, z, I, o, zdest;
  double den, inten;

  render_sort_coor(coor, x, 5, 5);
  if (x[11]<0 || x[1]>=h)  return(0);
  for (j=k=0; j<3; j++,k+=5) {
    x[0+k] = (x[0+k]<<8)+0x80;  x[2+k] = (x[2+k]<<8)+0x80;
    if (x[3+k]<0) x[3+k] = 0;  if (x[3+k]>=65535)  x[3+k] = 65535;
    if (x[4+k]<0) x[4+k] = 0;  if (x[4+k]>=65535)  x[4+k] = 65535;
    x[3+k] <<= 8;  x[4+k] <<= 8; }
  memcpy(x+15, x, 5*sizeof(long));
  for (j=k=0; j<3; j++,k+=5) {
    den = x[k+6]-x[k+1];  if (!den)  den = 1;  else  den = 1./den;
    for (l=0; l<5; l++)  if (l!=1)
      dx[l+k] = (x[k+5+l]-x[k+l])*den; }
  x1 = x2 = x[0];  z1 = z2 = x[2];  I1 = I2 = x[3];  o1 = o2 = x[4];
  if (z1<CutZ && x[7]<CutZ && x[12]<CutZ)  return(0);
  if (x[11]>h-1)  x[11] = h-1;
  if (x[11]-x[1]>10000)  return(0);
  za = abs(dx[2])+abs(dx[12])-1;  if (za>25499)  za = 25499;
  for (y=x[1]; y<=x[11]; y++) {
    if (y==x[6]) {
      x2 = x[5];  z2 = x[7];  I2 = x[8];  o2 = x[9];
      memcpy(dx, dx+5, 5*sizeof(long)); }
    if (y>=0) {
      if (x1<x2) {
        xi1 = x1>>8;  xi2 = x2>>8;
        if (j=(xi2-xi1)) {
          den = 1./j;
          dz = (z2-z1)*den;  ds = (I2-I1)*den;  dt = (o2-o1)*den; }
        z = z1+abs(dz)+za;  I = I1;  o = o1; }
      else {
        xi1 = x2>>8;  xi2 = x1>>8;
        if (j=(xi2-xi1)) {
          den = 1./j;
          dz = (z1-z2)*den;  ds = (I1-I2)*den;  dt = (o1-o2)*den; }
        z = z2+abs(dz)+za;  I = I2;  o = o2; }
      if (z1>0xFFFEFF || z2>0xFFFEFF) {
        if (z>0xFFFEFF && z1<=0xFFFEFF)  z = 0xFFFEFF;
        if (z2>z1 && z+z2-z1>0xFFFEFF && z2<=0xFFFEFF)  z = 0xFFFEFF+z1-z2; }
      dest = y*scan+xi1*3;  zdest = y*zscan+xi1;
      if (xi2>w-1)  xi2 = w-1;
      for (xi=xi1; xi<=xi2; xi++, dest+=3, zdest++) {
        if (xi>=0)
          if ((z>>8)<zbuf[zdest] && z>CutZ) {
            if (opacitytbl[OpacPos]<o) {
              zbuf[zdest] = z>>8;
              inten = (float)I/16776960;  if (inten>1)  inten = 1;
              scrbuf[dest]   = color[0]*inten;
              scrbuf[dest+1] = color[1]*inten;
              scrbuf[dest+2] = color[2]*inten; }
            OpacPos = (OpacPos+1)%OPACDITHER; }
        z += dz;  I += ds;  o += dt; } }
    x1 += dx[10];  x2 += dx[0];
    z1 += dx[12];  z2 += dx[2];
    I1 += dx[13];  I2 += dx[3];
    o1 += dx[14];  o2 += dx[4]; }
}

render_triangle_texture(long w, long h, long scan, uchar *scrbuf, long zscan,
                    ushort *zbuf, long coor[15], long pw, long ph, long ppal,
                    uchar *photo, long bgr, long scale, float bright)
/* Draw a triangle textured from a photograph in 3D.  Note that coordinates
 *  greater than +/-32767 in x or y will have bizarre effects.  If photo
 *  coordinates are not within the photograph, the photograph will be
 *  stretched to cover the entire triangle.
 * Enter: long w, h: size of output screen in pixels.  These values, along
 *                   with scrbuf, can be set to make a clipping rectangle.
 *        long scan: number of bytes per scan line on output screen.  Must be
 *                   at least 3*w.  This is negative for bgr organized data.
 *        uchar *scrbuf: buffer pointing to the top scanline on the screen.
 *        long zscan: number of shorts per scan line of zbuffer.  Must be at
 *                    least w.
 *        ushort *zbuf: zbuffer.  Lower is 'on top'.
 *        long coor[15]: coordinates of triangle points and photo points, in
 *                       screen and photo coordinates, respectively.  The z
 *                       values are in z buffer units, which range from 0 to
 *                       CutZ.  This array contains x0,y0,z0,s0,t0,x1,y1,z1,
 *                       s1,t1,x2,y2,z2,s2,t2.
 *        long pw, ph: the size of the photograph in pixels.
 *        long ppal: 0 for a 24-bit photograph, 1 for an 8-bit photograph.
 *        uchar *photo: pointer to the photograph.  8-bit photos have 768
 *                      bytes of palette followed by the actual photo data.
 *        long bgr: 0 if the photograph is the same sense as the output
 *                  screen, 1 to reverse sense.
 *        long scale: interpolation parameter (0-none, 1-bicubic).
 *        float bright: brightness to adjust picture by.       3/20/97-DWM */
{
  long x[20], dx[15], i, j, k, l, rr=2*bgr, bb=2-2*bgr, y, xi, src, dest;
  long x1, z1, s1, t1, x2, z2, s2, t2, xi1, xi2, dz, ds, dt, z, s, t, zdest;
  long src2, src3, src4, m1, m2, m3, m4, za, bl;
  double den;

  render_sort_coor(coor, x, 5, 5);
  if (x[11]<0 || x[1]>=h)  return(0);
  for (j=k=0; j<3; j++,k+=5) {
    x[0+k] = (x[0+k]<<8)+0x80;  x[2+k] = (x[2+k]<<8)+0x80;
    if (x[3+k]<0) x[3+k] = 0;  if (x[3+k]>=pw-scale)  x[3+k] = pw-1-scale;
    if (x[4+k]<0) x[4+k] = 0;  if (x[4+k]>=ph-scale)  x[4+k] = ph-1-scale;
    x[3+k] = (x[3+k]<<8)+0x80*(!scale); x[4+k] = (x[4+k]<<8)+0x80*(!scale); }
  memcpy(x+15, x, 5*sizeof(long));
  if (!ppal)  pw *= 3;
  for (j=k=0; j<3; j++,k+=5) {
    den = x[k+6]-x[k+1];  if (!den)  den = 1;  else  den = 1./den;
    for (l=0; l<5; l++)  if (l!=1)
      dx[l+k] = (x[k+5+l]-x[k+l])*den; }
  x1 = x2 = x[0];  z1 = z2 = x[2];  s1 = s2 = x[3];  t1 = t2 = x[4];
  if (z1<CutZ && x[7]<CutZ && x[12]<CutZ)  return(0);
  if (x[11]>h-1)  x[11] = h-1;
  if (x[11]-x[1]>10000)  return(0);
  za = abs(dx[2])+abs(dx[12])-1;  if (za>25499)  za = 25499;
  for (y=x[1]; y<=x[11]; y++) {
    if (y==x[6]) {
      x2 = x[5];  z2 = x[7];  s2 = x[8];  t2 = x[9];
      memcpy(dx, dx+5, 5*sizeof(long)); }
    if (y>=0) {
      if (x1<x2) {
        xi1 = x1>>8;  xi2 = x2>>8;
        if (j=(xi2-xi1)) {
          den = 1./j;
          dz = (z2-z1)*den;  ds = (s2-s1)*den;  dt = (t2-t1)*den; }
        z = z1+abs(dz)+za;  s = s1;  t = t1; }
      else {
        xi1 = x2>>8;  xi2 = x1>>8;
        if (j=(xi2-xi1)) {
          den = 1./j;
          dz = (z1-z2)*den;  ds = (s1-s2)*den;  dt = (t1-t2)*den; }
        z = z2+abs(dz)+za;  s = s2;  t = t2; }
      if (z1>0xFFFEFF || z2>0xFFFEFF) {
        if (z>0xFFFEFF && z1<=0xFFFEFF)  z = 0xFFFEFF;
        if (z2>z1 && z+z2-z1>0xFFFEFF && z2<=0xFFFEFF)  z = 0xFFFEFF+z1-z2; }
      dest = y*scan+xi1*3;  zdest = y*zscan+xi1;
      if (xi2>w-1)  xi2 = w-1;
      for (xi=xi1; xi<=xi2; xi++, dest+=3, zdest++) {
        if (xi>=0)
          if ((z>>8)<zbuf[zdest] && z>CutZ) {
            zbuf[zdest] = z>>8;
            if (!ppal)  src = (t>>8)*pw+(s>>8)*3;
            else        src = photo[bl=((t>>8)*pw+(s>>8)+768)]*3;
            if (!scale) {
              scrbuf[dest+rr] = photo[src]*bright;
              scrbuf[dest+1]  = photo[src+1]*bright;
              scrbuf[dest+bb] = photo[src+2]*bright; }
            else {
              if (!ppal) {
                src2 = src+3;  src3 = src+pw;  src4 = src3+3; }
              else {
                src2 = photo[bl+1]*3;
                src3 = photo[bl+pw]*3;
                src4 = photo[bl+pw+1]*3; }
              m1 = (256-(s&0xFF))*(256-(t&0xFF));
              m2 =      (s&0xFF) *(256-(t&0xFF));
              m3 = (256-(s&0xFF))*     (t&0xFF);
              m4 =      (s&0xFF) *     (t&0xFF);
              scrbuf[dest+rr] = ((long)((photo[src]*m1+photo[src2]*m2+
                         photo[src3]*m3+photo[src4]*m4)*bright))>>16;
              scrbuf[dest+1]  = ((long)((photo[src+1]*m1+photo[src2+1]*m2+
                         photo[src3+1]*m3+photo[src4+1]*m4)*bright))>>16;
              scrbuf[dest+bb] = ((long)((photo[src+2]*m1+photo[src2+2]*m2+
                         photo[src3+2]*m3+photo[src4+2]*m4)*bright))>>16; } }
        z += dz;  s += ds;  t += dt; } }
    x1 += dx[10];  x2 += dx[0];
    z1 += dx[12];  z2 += dx[2];
    s1 += dx[13];  s2 += dx[3];
    t1 += dx[14];  t2 += dx[4]; }
}

reset_geo()
/* Return the geometry to the default position.                 4/1/97-DWM */
{
  DATA *data;
  long i;
  HWND hwnd;

  hwnd = SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0);
  data = lock_window(hwnd);
  if (!data)  return(0);
  DefSize = data->mol.maxcheck;
  render_dlt_new(0, 0, 100, 100, Pref.perspec, data->renderval,
                 data->renderdlt);
  update_process(data);
  unlock_window(data);
  InvalidateRect(hwnd, 0, 0);
}

rotate_geo(long axis, float ang)
/* Rotate the geometry associated with the top window.
 * Enter: long axis: 0-x, 1-y, 2-z.
 *        float ang: angle in degrees.                         3/31/97-DWM */
{
  DATA *data;
  float dang[3]={0,0,0};
  long i;
  HWND hwnd;

  hwnd = SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0);
  data = lock_window(hwnd);
  if (!data)  return(0);
  for (i=0; i<18; i++)  if (data->renderdlt[i]) break;
  if (i==18) {
    unlock_window(data);  return(0); }
  dang[axis] = ang*deg;
  render_move(0, dang, 0, 0, 0, data->renderval, data->renderdlt,
              data->renderdlt);
  update_process(data);
  unlock_window(data);
  InvalidateRect(hwnd, 0, 0);
}

scale_pic24(uchar *dest, uchar *source, short srcw, short destw, short simgw,
            short dimgw, short simgh, short dimgh)
/* Scale a 24-bit color image.  Interpolation is always performed.  If the
 *  source image is not the same shape as the destination image, the edges
 *  are letterboxed.  The image's aspect ratio is maintained.
 * Enter: long dest, source: flat address values for images.
 *        int srcw, destw: width of arrays to hold images (bytes).
 *        int simgw, dimgw: width of actual image to scale (pixels).
 *        int simgh, dimgh: height of actual image to scale.    1/5/98-DWM */
{
  long *horz, *vert, i, j, d, w, h, inter=1;
  long w2, w3, h2, h3;
  float step;

  w = dimgw;  h = simgh*w/simgw;
  if (h>dimgh || w==simgw*dimgh/simgh) {
    h = dimgh;  w = simgw*h/simgh;  if (w>dimgw) w = dimgw; }
  if (w==simgw && h==simgh)  inter = 0;
  if (!(horz=malloc2((w+h)*sizeof(long)*4)))  return(0);
  vert = horz+w*4;
  step = min((float)simgw/w, (float)simgh/h);
  for (i=0; i<w; i++)
    horz[i] = floor(step*i+(simgw-(w-1)*step)/2)*3;
  for (i=0; i<h; i++)
    vert[i] = (long)(floor(step*i+(simgh-(h-1)*step)/2))*srcw;
  if (inter) {
    w2 = w+w;  w3 = w2+w;  h2 = h+h;  h3 = h2+h;
    for (i=0; i<w; i++) {
      horz[i+w] = min(ceil(step*i+(simgw-(w-1)*step)/2),simgw-1)*3;
      horz[i+w3] = 256*(step*i+(simgw-(w-1)*step)/2-horz[i]/3);
      horz[i+w2] = 256-horz[i+w3]; }
    for (i=0; i<h; i++) {
      vert[i+h]=(long)(min(ceil(step*i+(simgh-(h-1)*step)/2),simgh-1))*srcw;
      vert[i+h3] = 256*(step*i+(simgh-(h-1)*step)/2-vert[i]/srcw);
      vert[i+h2] = 256-vert[i+h3]; } }
  dest += ((dimgw-w)/2)*3+((dimgh-h)/2)*destw;
  if (!inter)
    for (j=d=0; j<h; j++, d+=destw-w*3)
      for (i=0; i<w; i++, d+=3) {
        ((short *)(dest+d))[0] = ((short *)(source+vert[j]+horz[i]))[0];
        dest[d+2] = source[vert[j]+horz[i]+2]; }
  else
    for (j=d=0; j<h; j++, d+=destw-w*3)
      for (i=0; i<w; i++, d+=3) {
        dest[d] = ((source[vert[j]  +horz[i]]  *horz[i+w2]+
                    source[vert[j]  +horz[i+w]]*horz[i+w3])*vert[j+h2]+
                   (source[vert[j+h]+horz[i]]  *horz[i+w2]+
                    source[vert[j+h]+horz[i+w]]*horz[i+w3])*vert[j+h3])>>16;
        dest[d+1]=((source[vert[j]+horz[i]+1]*horz[i+w2]+
                   source[vert[j]+horz[i+w]+1]*horz[i+w3])*vert[j+h2]+
                  (source[vert[j+h]+horz[i]+1]*horz[i+w2]+
                   source[vert[j+h]+horz[i+w]+1]*horz[i+w3])*vert[j+h3])>>16;
        dest[d+2]=((source[vert[j]+horz[i]+2]*horz[i+w2]+
                   source[vert[j]+horz[i+w]+2]*horz[i+w3])*vert[j+h2]+
                  (source[vert[j+h]+horz[i]+2]*horz[i+w2]+
                 source[vert[j+h]+horz[i+w]+2]*horz[i+w3])*vert[j+h3])>>16; }
  free2(horz);
}

scale_zbuf(ushort *dest, ushort *source, short simgw, short dimgw,
           short simgh, short dimgh)
/* Scale a zbuffer image.  When scaled down, rows are dropped; when scaled
 *  up, rows are duplicated.  If the source image is not the same shape as
 *  the destination image, the edges are letterboxed.  The image's aspect
 *  ratio is maintained.  The letterboxed area is not modified.
 * Enter: ushort *dest, *source: address values for zbuffers.
 *        int simgw, dimgw: width of image to scale.
 *        int simgh, dimgh: height of image to scale.           1/5/98-DWM */
{
  long *horz, *vert, i, j, d, w, h;
  float step;

  w = dimgw;  h = simgh*w/simgw;
  if (h>dimgh || w==simgw*dimgh/simgh) {
    h = dimgh;  w = simgw*h/simgh;  if (w>dimgw) w = dimgw; }
  if (!(horz=malloc2((w+h)*sizeof(long))))  return(0);
  vert = horz+w;
  step = min((float)simgw/w, (float)simgh/h);
  for (i=0; i<w; i++)
    horz[i] = floor(step*i+(simgw-(w-1)*step)/2);
  for (i=0; i<h; i++)
    vert[i] = (long)(floor(step*i+(simgh-(h-1)*step)/2))*simgw;
  dest += (dimgw-w)/2+(dimgh-h)/2*dimgw;
  for (j=d=0; j<h; j++, d+=dimgw-w)
    for (i=0; i<w; i++, d++)
      dest[d] = source[vert[j]+horz[i]];
  free2(horz);
}

shift_geo(long axis, float jump)
/* Shift the geometry associated with the top window.
 * Enter: long axis: 0-x, 1-y, 2-z.
 *        float jump: portion of the screen to jump.           3/31/97-DWM */
{
  DATA *data;
  float dx[3]={0,0,0};
  long i;
  HWND hwnd;

  hwnd = SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0);
  data = lock_window(hwnd);
  if (!data)  return(0);
  for (i=0; i<18; i++)  if (data->renderdlt[i]) break;
  if (i==18) {
    unlock_window(data);  return(0); }
  dx[axis] = jump;
  render_move(dx, 0, 0, 0, 0, data->renderval, data->renderdlt,
              data->renderdlt);
  update_process(data);
  unlock_window(data);
  InvalidateRect(hwnd, 0, 0);
}

uchar *smooth_palette(DATA *data, long numcolors)
/* Create a smoothly changing palette based on spectrum data.
 * Enter: DATA *data: either a pointer to a structure which contains the
 *                    colors to use, or null to use the default colors.
 *        long numcolors: the number of colors to use for the smooth
 *                        palette.  This must be a multiple of 3.  The
 *                        palette will be 3 times this length (in bytes).  An
 *                        extra byte is allocated at the end of the palette
 *                        to allow access to the palette entries using (long)
 *                        variables.
 * Exit:  uchar *pal: allocated palette, or null for an error. 3/21/97-DWM */
{
  long i, j, d, s, clr[4];
  PREF *pref;
  float v;
  uchar *pal, *cclr;
  uchar winpal[]={0,0,0, 128,0,0, 0,128,0, 128,128,0, 0,0,128, 128,0,128, 0,128,128, 151,151,151, 128,128,128, 255,0,0, 192,192,192, 255,255,0, 0,0,255, 128,128,255, 0,255,255, 255,255,255};

  clr[0] = get_color(data, BACKCOLOR);
  clr[1] = get_color(data, POSCOLOR);
  clr[2] = get_color(data, NEGCOLOR);
  clr[3] = get_color(data, ASYMCOLOR);
  if (!(pal=malloc2((numcolors+16)*3)))  return(0);
  cclr = (uchar *)clr;
  numcolors -= numcolors%3;
  for (i=d=0; i<numcolors; i++, d+=3) {
    j = i;
    s = j*3/numcolors;
    v = ((float)j*3/numcolors-s)*(numcolors/3)/((numcolors/3)-1);
    pal[d]   = cclr[2]*(1-v)+cclr[s*4+4+2]*v;
    pal[d+1] = cclr[1]*(1-v)+cclr[s*4+4+1]*v;
    pal[d+2] = cclr[0]*(1-v)+cclr[s*4+4]*v; }
  memcpy(pal+d, winpal, 16*3);
  return(pal);
}

long sphere_inter(real *x, real *v, real radius, real *inter)
/* Based on a ray starting a point x and extending in direction v, calculate
 *  the point(s) where it cross a sphere centered at 0,0,0 and of the
 *  specified radius.  Up to two points can be returned.  If a single point
 *  is returned, the starting coordinate is within the sphere.
 * Enter: real *x: starting coordinates.
 *        real *v: initial vector.
 *        real radius: sphere radius.
 *        real *inter: location to store up to 2 points.
 * Exit:  long hit: 0 for missed, 1 for within radius, 2 for two intercepts.
 *                                                             7/29/97-DWM */
{
  real v1[3], temp[3], min, cost, curr, curr1, curr2;

  memcpy(v1, v, 3*sizeof(real));
  calc_unit(v1);
  min = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
  cost = -(x[0]*v1[0]+x[1]*v1[1]+x[2]*v1[2]);
  curr = radius*radius+(cost*cost-min);
  if (curr<0)
    return(0);
  curr1 = cost - sqrtl(curr);
  curr2 = cost + sqrtl(curr);
  if (curr2<0)  return(0);
  if (curr1<0) {
    inter[0] = x[0]+v1[0]*curr2;
    inter[1] = x[1]+v1[1]*curr2;
    inter[2] = x[2]+v1[2]*curr2;
    return(1); }
  inter[0] = x[0]+v1[0]*curr2;  inter[3] = x[0]+v1[0]*curr1;
  inter[1] = x[1]+v1[1]*curr2;  inter[4] = x[1]+v1[1]*curr1;
  inter[2] = x[2]+v1[2]*curr2;  inter[5] = x[2]+v1[2]*curr1;
  if (sq(inter[0]-x[0])+sq(inter[1]-x[1])+sq(inter[2]-x[2])>
      sq(inter[3]-x[0])+sq(inter[4]-x[1])+sq(inter[5]-x[2])) {
    inter[0] = x[0]+v1[0]*curr1;  inter[3] = x[0]+v1[0]*curr2;
    inter[1] = x[1]+v1[1]*curr1;  inter[4] = x[1]+v1[1]*curr2;
    inter[2] = x[2]+v1[2]*curr1;  inter[5] = x[2]+v1[2]*curr2; }
  return(2);
}

stereo_figure(HWND hdlg, long mode, long swap, STEREO *st)
/* Draw a preview picture based on the current stereo mode.
 * Enter: HWND hdlg: pointer to owner dialog.
 *        long mode: mode number to draw (0-6).
 *        long swap: 0 for no swap, 1 for swap.  Other values are invalid.
 *        STEREO *st: pointer to stereo record containing image to use with
 *                    stereogram.                               1/1/98-DWM */
{
  RECT rect;
  HDC hdc;
  HANDLE dest=0, bmp=0;
  uchar *src;
  long box[]={94, 4, 94+136, 4+96}, out[4];
  long w, h, pal=0, ow=PREVIEWWIDTH, oh=PREVIEWHEIGHT;
  long i, j, k, od, rr, gg, bb;
  long tim=clock();

  if (!PreviewGraphic)  return(0);
  hdc = GetDC(hdlg);
  map_dialog_rect(hdlg, box[0], box[1], box[2], box[3], out);
  w = out[2]-out[0];  h = out[3]-out[1];
  if (w<ow || h<oh)  return(0);
  out[0] += (w-ow)/2;  out[1] += (h-oh)/2;  w = ow;  h = oh;
  src = malloc2(ow*oh*3);
  if (src) {
    switch (mode) {
      case StereoSTEREOGRAM: for (j=0; j<oh; j++)
          for (i=0; i<ow; i++) {
            if (!st->image || !(st->flags&8)) {
              od = (j*ow+i)%OPACDITHER;
              rr = opacitytbl[od]&0xFF;
              gg = (opacitytbl[od]>>8)&0xFF;
              bb = (opacitytbl[od]>>16)&0xFF; }
            else {
              if (st->pal)  od = st->image[768+(j%st->h)*st->w+(i%st->w)]*3;
              else          od = ((j%st->h)*st->w+(i%st->w))*3;
              rr = st->image[od];
              gg = st->image[od+1];
              bb = st->image[od+2]; }
            src[j*ow*3+i*3] = rr;
            src[j*ow*3+i*3+1] = gg;
            src[j*ow*3+i*3+2] = bb; } break;
      /** Additional preview modes go here **/
      default:
        k = mode+swap;  if (mode==StereoINTERLACED && swap)  k = 7;
        memcpy(src, PreviewGraphic+k*ow*oh*3, ow*oh*3); }
    dest = unlock2(src); }
  if (dest && BitsPixel<=8) {
    src = PalettizeGraphic(w, h, 1, dest, 1, 0);
    free2(lock2(dest));
    pal = 1;
    dest = src; }
  if (dest) {
    bmp = GraphicToBMP(w, h, pal, dest);
    free2(lock2(dest)); }
  if (bmp) {
    draw_bmp(hdc, out[0], out[1], bmp);
    free2(lock2(bmp)); }
  ReleaseDC(hdlg, hdc);
  tim = 1000*(tim-clock())/CLOCKS_PER_SEC;
  if (tim<100)       PreviewDelay = 100;
  else if (tim<200)  PreviewDelay = 500;
  else if (tim>1000) PreviewDelay = 2000;
  else               PreviewDelay = 1000;
}

stereogram(long w, long h, long scan, uchar *scrbuf, ushort *zbuf,
           DATA *data)
/* Convert information in the zbuffer into a single image image-based or
 *  random dot stereogram.
 * Enter: long w, h: size of output screen in pixels.
 *        long scan: number of bytes per scan line on output screen.  Must be
 *                   at least 3*w.  This is negative for bgr organized data.
 *        uchar *scrbuf: buffer pointing to the top scanline on the screen.
 *        ushort *zbuf: zbuffer.  Lower is 'on top'.
 *        DATA *data: pointer to window data to draw, used for interocular
 *                    distance and the stereogram source image. 1/3/98-DWM */
{
  long i, j, d, zd, r, od;
  float *buf, l, s, z;
  real inter=data->stereo.interocular[0];
  uchar rr, gg, bb;

  if (!(buf=malloc2(w*sizeof(float))))  return(0);
  for (j=0; j<h; j++) {
    zd = j*w;  d = j*scan;
    memset(buf, -1, w*sizeof(float));
    for (i=0; i<w; i++) {
      z = 1-(float)zbuf[zd+i]/65535;
      s = inter/2-(inter/6)*z;
      r = (long)(i+s/2);  l = r-s;
      if (l<0 || r<0 || l>w-1 || r>w-1)  continue;
      buf[r] = l; }
    for (i=0; i<w; i++) {
      if (buf[i]<0) {
        if (!data->stereo.image || !(data->stereo.flags&8)) {
          od = (j*w+i)%OPACDITHER;
          rr = opacitytbl[od]&0xFF;
          gg = (opacitytbl[od]>>8)&0xFF;
          bb = (opacitytbl[od]>>16)&0xFF; }
        else {
          if (data->stereo.pal)
            od = data->stereo.image[768+(j%data->stereo.h)*data->stereo.w+
                                    (i%data->stereo.w)]*3;
          else
            od = ((j%data->stereo.h)*data->stereo.w+(i%data->stereo.w))*3;
          rr = data->stereo.image[od];
          gg = data->stereo.image[od+1];
          bb = data->stereo.image[od+2]; }
        scrbuf[d+i*3] = bb;
        scrbuf[d+i*3+1] = gg;
        scrbuf[d+i*3+2] = rr; }
      else {
        od = floor(buf[i]);  s = buf[i]-od;  od *= 3;
        scrbuf[d+i*3]   = scrbuf[d+od]  *(1-s)+scrbuf[d+od+3]*s;
        scrbuf[d+i*3+1] = scrbuf[d+od+1]*(1-s)+scrbuf[d+od+4]*s;
        scrbuf[d+i*3+2] = scrbuf[d+od+2]*(1-s)+scrbuf[d+od+5]*s; } } }
  free2(buf);
}

uchar *update_bmp(DATA *data, long w, long h, long *freeimg, long bits)
/* Draw the appropriate picture for the current data, filling a space of the
 *  specified width, height, and bit depth.
 * Enter: DATA *data: pointer to window data.
 *        long w, h: size of desitination screen.
 *        long *freeimg: pointer to location to store if the returned image
 *                       was allocated by this routine (1) or not (0).
 *                       Set this initially to zero if the calling program
 *                       will not modify (including unlock) the return data
 *                       if it was not allocated.  Set to non-zero if the
 *                       returned data will be manipulated.
 *        long bits: bits per pixel (either 8 or 24).
 * Exit:  uchar *bmp: pointer to store windows-style bmp, with 40 byte
 *                    header.                                 12/26/97-DWM */
{
  long i, j, d, s, ster, mode, scan, scan8, refresh=0, twice=1, r, b, I;
  long color;
  uchar *pic, *pic2, *new;
  ushort *zbuf=0, *zbuf2;
  double *dlt, V, H;
  LPBITMAPINFOHEADER bmp;

  OpacPos = 0;
  mode = (data->dflag>>5)&3;
  data->dflag = (data->dflag&0xFFFFFFE7F)|(mode<<7);
  if (mode==1 || (!mode && data->asym.opacity))  fill_norm(data);
  if (data->stereo.mode==StereoSTEREOGRAM && mode==2)  mode = 3;
  if (mode==2 && bits==24 && data->render.buf && !freeimg[0] &&
      w==data->render.lastw && h==data->render.lasth)
    return(data->render.buf);
  if (mode!=2 && (data->stereo.mode==StereoREDBLUE ||
      data->stereo.mode==StereoOVERLAY))
    twice = 2;
  scan = (w*3+3)&0x7FFFFFFC;
  if (!(pic=malloc2(scan*h*twice+40+1024*(bits<=8)+3)))  return(0);
  freeimg[0] = 1;
  pic2 = pic+40+1024*(bits<=8);
  if (mode!=2 || !data->render.buf || data->render.lastw!=w ||
      data->render.lasth!=h) {
    color = get_color(data, BACKCOLOR);
    if (data->stereo.mode!=StereoSTEREOGRAM)
      for (j=0; j<h; j++)
        fill_zone24(pic2+scan*j, color, w); }
  if (mode!=2) {
    if (!(zbuf=malloc2(w*h*twice*sizeof(short)))) { free2(pic); return(0); }
    memset(zbuf, 0xFF, w*h*twice*sizeof(short)); }
  for (ster=0; ster<2-(mode>=2); ster++) {
    if (ster==1 && data->stereo.mode!=StereoSTEREOSCOPE && data->stereo.mode
        !=StereoINTERLACED && data->stereo.mode!=StereoREDBLUE &&
        data->stereo.mode!=StereoOVERLAY)
      continue;
    dlt = data->renderdlt;
    if (mode<2)
      switch (data->stereo.mode) {
        case StereoSTEREOSCOPE: if (zbuf) for (i=0; i<h; i++) {
            memset(zbuf+i*w,  0xFF*(!ster), (w/2)*sizeof(short));
            memset(zbuf+i*w+w/2, 0xFF*ster, (w/2)*sizeof(short)); }
          dlt = data->renderdlt+18*(ster+1); break;
        case StereoINTERLACED: if (zbuf) for (i=0; i<h; i++)
            memset(zbuf+i*w,  0xFF*(ster!=((i^data->stereo.flags)&1)),
                   w*sizeof(short));
          dlt = data->renderdlt+18*(ster+1); break;
        case StereoREDBLUE:  case StereoOVERLAY:
          if (ster) { pic2 += scan*h;  zbuf2 = zbuf;  zbuf += w*h; }
          dlt = data->renderdlt+18*(ster+1); break;
        default: ; }
    switch (mode) {
      case 0: if (!data->asym.opacity)
          render_coor(data, dlt, &data->points.scrxyz, data->points.xyz,
                      data->points.n, 0, 0, 0);
        else
          render_coor(data, dlt, &data->points.scrxyz, data->points.xyz,
                      data->points.n, &data->asym.scrxyz, data->asym.xyz,
                      data->asym.n);
        render_points(w, h, -scan, pic2+scan*(h-1), zbuf, data, 1);
        OpacPos = 0;
        if (data->asym.opacity)
          render_asymptote(w, h, -scan, pic2+scan*(h-1), zbuf, data); break;
      case 1: if (!data->asym.opacity)
          render_coor(data, dlt, &data->poly.scrxyz, data->poly.xyz,
                      data->poly.n, 0, 0, 0);
        else
          render_coor(data, dlt, &data->poly.scrxyz, data->poly.xyz,
                      data->poly.n, &data->asym.scrxyz, data->asym.xyz,
                      data->asym.n);
        render_polygon(w, h, -scan, pic2+scan*(h-1), zbuf, data);
        OpacPos = 0;
        if (data->asym.opacity)
          render_asymptote(w, h, -scan, pic2+scan*(h-1), zbuf, data); break;
      case 2: if (!data->render.buf)  break;
        if (data->render.lastw==w && data->render.lasth==h)
          memcpy(pic2, data->render.buf+40, fabs(data->render.scan)*
                 data->render.lasth);
        else
          scale_pic24(pic2, data->render.buf+40, -data->render.scan, scan,
                      data->render.lastw, w, data->render.lasth, h); break;
      case 3: if (!data->render.zbuf) break;
        if (data->render.lastw==w && data->render.lasth==h)
          memcpy(zbuf, data->render.zbuf, w*h*sizeof(ushort));
        else
          scale_zbuf(zbuf, data->render.zbuf, data->render.lastw, w,
                     data->render.lasth, h); } }
  if (mode!=2)  switch (data->stereo.mode) {
    case StereoREDBLUE: zbuf = zbuf2;  pic2 -= scan*h;
      for (j=0; j<h; j++)  for (i=0; i<w; i++) {
        r = b = -1;
        if (zbuf[j*w+i]!=0xFFFF)
          r = max(max(pic2[scan*(h-1)+j*-scan+i*3], pic2[scan*(h-1)+j*-scan+
                  i*3+1]), pic2[scan*(h-1)+j*-scan+i*3+2]);
        if (zbuf[j*w+i+w*h]!=0xFFFF)
          b = max(max(pic2[scan*(2*h-1)+j*-scan+i*3], pic2[scan*(2*h-1)+j*
                  -scan+i*3+1]), pic2[scan*(2*h-1)+j*-scan+i*3+2]);
        if (r<0 && b<0)  continue;
        if (!(data->stereo.flags&2)) { d = b;  b = r;  r = d; }
        if (r<0)  pic2[scan*(h-1)+j*-scan+i*3] = 0;
        else      pic2[scan*(h-1)+j*-scan+i*3] = r;
        pic2[scan*(h-1)+j*-scan+i*3+1] = 0;
        if (b<0)  pic2[scan*(h-1)+j*-scan+i*3+2] = 0;
        else      pic2[scan*(h-1)+j*-scan+i*3+2] = b; } break;
    case StereoSTEREOGRAM:
      stereogram(w, h, -scan, pic2+scan*(h-1), zbuf, data);  break;
    case StereoOVERLAY: zbuf = zbuf2;  pic2 -= scan*h;
      for (j=0; j<h; j++)  for (i=0; i<w; i++) {
        d = scan*(h-1)+j*-scan+i*3;
        for (s=0; s<3; s++)
          if (zbuf[j*w+i+w*h]!=0xFFFF) {
            if (zbuf[j*w+i]==0xFFFF)
              pic2[d+s] = pic2[d+scan*h+s]*0.5;
            else
              pic2[d+s] = (pic2[d+s]+pic2[d+scan*h+s])*0.5; }
          else if (zbuf[j*w+i]!=0xFFFF)
            pic2[d+s] *= 0.5; } break;
    case StereoCHROMADEPTH:
      for (j=0; j<h; j++)  for (i=0; i<w; i++) if (zbuf[j*w+i]!=0xFFFF) {
        d = scan*(h-1)+j*-scan+i*3;
        H = 5.*zbuf[j*w+i]/65535;
        I = floor(H);  H = 1-(H-I);
        V = max(max(pic2[d], pic2[d+1]), pic2[d+2]);
        switch (I) {
          case 0: pic2[d+2]=V;       pic2[d+1]=V*(1-H); pic2[d]=0; break;
          case 1: pic2[d+2]=V*H;     pic2[d+1]=V;       pic2[d]=0; break;
          case 2: pic2[d+2]=0;       pic2[d+1]=V;     pic2[d]=V*(1-H); break;
          case 3: pic2[d+2]=0;       pic2[d+1]=V*H;     pic2[d]=V; break;
          case 4: pic2[d+2]=V*(1-H); pic2[d+1]=0;       pic2[d]=V; break;
          case 5: pic2[d+2]=V;       pic2[d+1]=0;       pic2[d]=V*H; } }
      break;
    default: ; }
  if (zbuf)  free2(zbuf);
  bmp = pic;
  if ((bits<=8 && ColorTable && Hpal)) {
    scan8 = (w+3)&0x7FFFFFFC;
    s = d = 40+1024;
    for (j=0; j<h; j++, s+=scan-w*3, d+=scan8-w)
      for (i=0; i<w; i++, s+=3, d++)
        if (((short *)(pic+s))[0]==1)
          pic[d] = pic[s+2];
        else if (!Dither)
          pic[d] = ColorTable[((pic[s+2]&0xF8)<<8)+((pic[s+1]&0xFC)<<3)+
                              (pic[s]>>3)];              /* bgr to palette */
        else {
          pic[d] = ColorTable[dithertbl2[pic[s+2]+dithertbl1[dithernum]]+
                   dithertbl2[pic[s+1]+dithertbl1[dithernum+1]+300]+
                   dithertbl2[pic[s]+dithertbl1[dithernum+2]+600]];
          dithernum+=3;  if (dithernum>2040) dithernum = rand()%100; }
    new = realloc2(pic, 40+1024+scan8*h);
    if (new)  pic = new;  scan = scan8;
    bmp = pic;
    bmp->biBitCount = 8;
    for (i=0; i<256; i++) {
      pic[40+i*4]   = MasterPal[i*3+2];
      pic[40+i*4+1] = MasterPal[i*3+1];
      pic[40+i*4+2] = MasterPal[i*3];
      pic[40+i*4+3] = 0; } }
  else
    bmp->biBitCount = 24;
  if (bits>8 && bits<=16) {
    if (dither16num<0) {
      for (i=0; i<2048; i++)
        dither16tbl[i] = rand()&7;
      dither16num = 0; }
    for (i=40; i<w*h*3+40; i+=3) {
      pic[i] = min(255,pic[i]+dither16tbl[dither16num]);
      pic[i+1] = min(255,pic[i+1]+(dither16tbl[dither16num+1]&3));
      pic[i+2] = min(255,pic[i+2]+dither16tbl[dither16num+2]);
      dither16num = (dither16num+3)%2046; } }
  bmp->biSize = 40;  bmp->biWidth = w;  bmp->biHeight = h;
  bmp->biPlanes = 1;  bmp->biCompression = 0;
  bmp->biSizeImage = scan*h;
  bmp->biClrUsed = bmp->biClrImportant = 0;
  return(pic);
}

update_geo(HWND hwnd, DATA *data, HDC dc)
/* Redraw a window containing a geometry, or just ensure that the DLT
 *  parameters are correct.
 * Enter: HWND hwnd: handle of window.
 *        DATA *data: locked down data area for this window.
 *        HDC dc: device context for drawing, or 0 to just check DLTs.
 *                                                             3/12/97-DWM */
{
  RECT rect;
  long w, h, w2, h2, i, freeimg=0;
  LPBITMAPINFOHEADER bmp;
  static lastw=-1, lasth=-1;

  if (IsIconic(hwnd))  return(0);
  if (hwnd==SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0))
    if (HwndStat)
      SendMessage(HwndStat, SB_SETTEXT, 3, (LPARAM)DispInfo[0]);
  GetClientRect(hwnd, &rect);
  w = rect.right;  h = rect.bottom;
  if (w<2 || h<2) { w = lastw;  h = lasth; }
  if (w<1) w = 1;  if (h<1) h = 1;
  w2 = w;  h2 = h;
  if (((data->dflag>>5)&3)==2)
    if (((data->dflag>>9)&1) && data->w>0 && data->h>0) {
      w2 = data->w;  h2 = data->h; }
  for (i=0; i<18; i++)  if (data->renderdlt[i]) break;
  if (i==18) {
    DefSize = data->mol.maxcheck;
    render_dlt_new(0, 0, w2, h2, Pref.perspec, data->renderval,
                   data->renderdlt);
    update_process(data); }
  else {
    render_move(0, 0, 0, w2, h2, data->renderval, data->renderdlt,
                data->renderdlt);
    render_physical(data->renderval,data->renderdlt, data->render.camera, 0);
    if (data->render.w!=w2 || data->render.h!=h2) {
      data->render.w = w2;  data->render.h = h2;
      if (data->render.process>=3)  data->render.process = 2; }
    prep_stereo(data); }
  if (!dc)  return(0);
  bmp = update_bmp(data, lastw=w, lasth=h, &freeimg, BitsPixel);
  if (bmp)  draw_bmp2(dc, 0, 0, bmp);
  if (((data->dflag>>7)&3)==1 && data->render.process==2 &&
      data->render.lastw==w2 && data->render.lasth==h2 && data->render.buf
      && BitsPixel>16) {
    free2(data->render.buf);  freeimg = 0;  data->render.buf = bmp;
    data->render.antialias = (data->render.antialias&0xFFFFFFF3)|(2<<2); }
  if (freeimg)  free2(bmp);
}

update_window(HWND hwnd)
/* Main paint routine for all windows.  This locks down the window, starts
 *  the paint process, then calls the appropriate sub-update routine.
 * Enter: HWND hwnd: handle to window to redraw.  The window is guaranteed
 *                   not to be iconic.
 * Exit:  long okay: 0 for failed, 1 for okay.                 3/12/97-DWM */
{
  DATA *data;
  HDC hdc;
  PAINTSTRUCT paint;
  HBRUSH new, old;
  static HWND lasthwnd=-1;

  data = lock_window(hwnd);
  hdc = BeginPaint(hwnd, &paint);
  if (Hpal) {
    SelectPalette(hdc, Hpal, 0);
    RealizePalette(hdc); }
  LastMove = 0;
  if (data)
    update_geo(hwnd, data, hdc);
  EndPaint(hwnd, &paint);
  unlock_window(data);
  if (hwnd==SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0))
    recheck(hwnd, 0);
  return(1);
}

use_palette(uchar *pal)
/* Use a palette for palettization purposes.  This stores the palette in
 *  MasterPal and sets up the references in ColorTable.
 * Enter: uchar *pal: palette to use.
 * Exit:  long exact: 0 for insufficient memory, 2 for okay.   3/14/97-DWM */
{
  uchar xr[256], xg[256], xb[256], nr[256], ng[256], nb[256];
  long i, r, g, b, c, changed=1;

  if (!ColorTable)
    if (!(ColorTable=malloc2(65536))) return(0);
  memcpy(MasterPal, pal, 768);
  memset(ColorTable, 0, 65536);
  for (i=0; i<256; i++) {
    xr[i] = nr[i] = (pal[i*3]>>3);
    xg[i] = ng[i] = (pal[i*3+1]>>2);
    xb[i] = nb[i] = (pal[i*3+2]>>3); }
  while (changed) {
    changed = 0;
    for (i=1; i<256; i++)
      if (nr[i]!=0xFF) {
        for (c=0, r=nr[i], g=ng[i]; g<=xg[i]; g++)
          for (b=nb[i]; b<=xb[i]; b++)
            if (!ColorTable[(r<<11)+(g<<5)+b])
              ColorTable[(r<<11)+(g<<5)+b] = c = changed = i;
        for (r++; r<xr[i]; r++) {
          for (g=ng[i], b=nb[i]; b<=xb[i]; b++)
            if (!ColorTable[(r<<11)+(g<<5)+b])
              ColorTable[(r<<11)+(g<<5)+b] = c = changed = i;
          for (g++; g<xg[i]; g++) {
            b = nb[i];  if (!ColorTable[(r<<11)+(g<<5)+b])
                          ColorTable[(r<<11)+(g<<5)+b] = c = changed = i;
            b = xb[i];  if (!ColorTable[(r<<11)+(g<<5)+b])
                          ColorTable[(r<<11)+(g<<5)+b] = c = changed = i; }
          for (b=nb[i]; b<=xb[i] && g<=xg[i]; b++)
            if (!ColorTable[(r<<11)+(g<<5)+b])
              ColorTable[(r<<11)+(g<<5)+b] = c = changed = i; }
        for (g=ng[i]; g<=xg[i] && r<=xr[i]; g++)
          for (b=nb[i]; b<=xb[i]; b++)
            if (!ColorTable[(r<<11)+(g<<5)+b])
              ColorTable[(r<<11)+(g<<5)+b] = c = changed = i;
        if (!c) { nr[i] = 0xFF; }
        else {
          nr[i] = max(0, nr[i]-1);  xr[i] = min(31, xr[i]+1);
          ng[i] = max(0, ng[i]-1);  xg[i] = min(63, xg[i]+1);
          nb[i] = max(0, nb[i]-1);  xb[i] = min(31, xb[i]+1); } } }
  r = nr[0];  g = ng[0];  b = nb[0];
  ColorTable[(r<<11)+(g<<5)+b] = 0;
  return(2);
}

zone_geo(HWND hwnd, DATA *data, long xc, long yc, long *buf)
/* Determine which zone of a geo window a point is in, and what the relative
 *  size and position within that zone is.
 * Enter: HWND hwnd: handle of window.
 *        DATA *data: locked down data area for this window.
 *        long xc, yc: absolute coordinates within window.
 *        long *buf: array for 7 return values.  These are: 0-null, 1-x
 *                   within zone; 2-y within zone, 3-x of zone's left edge,
 *                   4-y of zone's top edge, 5-width of zone in pixels,
 *                   6-height of zone in pixels.                4/2/97-DWM */
{
  RECT rect;

  GetClientRect(hwnd, &rect);
  buf[0] = 0;  buf[1] = xc;  buf[2] = yc;
  buf[3] = 0;  buf[4] = 0;  buf[5] = rect.right;  buf[6] = rect.bottom;
}

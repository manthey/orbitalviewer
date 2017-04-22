/* This contains fully ANSI-compatible command line parameter access to the
 *  orbital routines.                                          11/9/97-DWM */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "matrix.h"
#include "dlt.h"
#include "ansi.h"

FILE *Outfptr;

#include "common.c"

void free2(void *memblock)
/* Call free.
 * Enter: void *memblock: pointer to memory to be freed.       11/9/97-DWM */
{
  free(memblock);
}

int main(int argc, char *argv[])
{
  DATA data;
  FILE *fptr;
  long seq=0, first=1, ppm=1, ovw=0, ovh=0;
  char name[ONENAMELEN], *ext[2]={".WRL",".PPM"};
  long tim=clock();

  if (argc<1 || argc>4 || argv[argc-1][1]=='?') {
    printf("Command line orbital calculator.\n\n"
           "Syntax:  ANSIORB (orb file) (output file) [screen size]\n"
           "         ANSIORB (orb file) > output file\n"
           "         ANSIORB < orb file > output file\n"
           "         ANSIORB /?\n\n"
           "The screen size is given as widthxheight (e.g, 640x480).\n"
           "See OV.PDF for orb file format.\n");
    exit(0); }
  Endian = 1;
  if (((uchar *)&Endian)[0]==1)  Endian = 0;
  memset(&data, 0, sizeof(DATA));
  if (!new_window_mid(&data)) {
    fprintf(stderr, "Error initializing memory.\n");  exit(0); }
  if (argc>=4) {
    sscanf(argv[3], "%d%*c%d", &ovw, &ovh);
    if (ovw>0 && ovh<=0)  ovh = ovw; }
  if (argc>=2) {
    if (!(fptr=fopen(argv[1], "rt"))) {
      fprintf(stderr, "Can't open input file.\n");  exit(0); } }
  else
    fptr = stdin;
  if (open_orb(fptr, &data, 0)) {
    if (argc>=2)  fclose(fptr);
    fprintf(stderr, "Error reading input file.\n");  exit(0); }
  if (ovw>0) {
    data.w = ovw;  data.h = ovh;
    data.dflag |= (1<<9); }
  if (argc>=2)  fclose(fptr);
  if (data.seq[0] && data.seq[1] && data.frame!=data.lastframe)
    seq = 1;
  do {
    if (!first && seq) {
      if (data.frame<data.lastframe)  data.frame++;
      else                            data.frame--;
      play_frame(&data, data.frame); }
    switch ((data.dflag>>3)&3) {
      case 0: if (orb_points(&data.points,&data.cut,&data.mol,0,0)!=1) {
          fprintf(stderr, "Failed to create points.\n");  exit(0); }
        if (data.asym.opacity)
          if (orb_asymptote(&data.asym,&data.cut, &data.mol,0,0)!=1) {
            fprintf(stderr, "Failed to create asymptotes.\n");  exit(0); }
        ppm = 0; break;
      case 1: if (orb_polygons(&data.poly,&data.cut,&data.mol,0,0)!=1) {
          fprintf(stderr, "Failed to create polygons.\n");  exit(0); }
        if (data.asym.opacity)
          if (orb_asymptote(&data.asym,&data.cut, &data.mol,0,0)!=1) {
            fprintf(stderr, "Failed to create asymptotes.\n");  exit(0); }
        ppm = 0; break;
      case 2: if (orb_render(&data.render,&data.stereo,&data.cut,&data.mol,
                             0,6)!=1) {
          fprintf(stderr, "Failed to render orbital.\n");  exit(0); }
        update_geo(&data);
        if (orb_render(&data.render,&data.stereo,&data.cut,&data.mol,
                       0,0)!=1) {
          fprintf(stderr, "Failed to render orbital.\n");  exit(0); } }
    if (seq && data.basename[0]) {
      sprintf(name, "%s%d.%s", data.basename, data.frame, ext[ppm]);
      if (!(fptr=fopen(name, "wt"))) {
        fprintf(stderr, "Failed to open sequence output file %s.\n", name);
        exit(0); } }
    else if (argc>=3) {
      if (!(fptr=fopen(argv[2], "wt"))) {
        fprintf(stderr, "Failed to open output file.\n");  exit(0); } }
    else
      fptr = stdout;
    if (!ppm) {
      Outfptr = fptr;
      if (save_vrml(&data)) {
        fprintf(stderr, "Failed to write VRML file.\n");
        if ((seq && data.basename[0]) || argc>=3)  fclose(fptr);
        exit(0); } }
    else
      if (save_ppm(&data, fptr)) {
        fprintf(stderr, "Failed to write PPM file.\n");
        if ((seq && data.basename[0]) || argc>=3)  fclose(fptr);
        exit(0); }
    if ((seq && data.basename[0]) || argc>=3)  fclose(fptr);
    orb_points(&data.points, 0, &data.mol, 0, 4);
    orb_polygons(&data.poly, 0, &data.mol, 0, 4);
    orb_asymptote(&data.asym, 0, &data.mol, 0, 4);
    orb_render(&data.render, 0, 0, &data.mol, 0, 4); }
  while (seq && data.frame!=data.lastframe);
  if (argc>=3)
    printf("Elapsed time: %4.2f seconds.\n", (float)(clock()-tim)/CLOCKS_PER_SEC);
}


void *malloc2(size_t size)
/* Call malloc.
 * Enter: size_t size: amount of memory to allocate.
 * Exit:  void *memblock: pointer to allocated memory or null for failed.
 *                                                             11/9/97-DWM */
{
  return(malloc(size));
}

void *realloc2(void *memblock, size_t size)
/* Call realloc.
 * Enter: void *memblock: current pointer to memory.
 *        size_t size: amount of memory to actually have allocated.
 * Exit:  void *memblock: pointer to realloced array or 0 for failed.
 *                                                             11/9/97-DWM */
{
  return(realloc(memblock, size));
}

int save_ppm(DATA *data, FILE *fptr)
/* Write an ascii P3 style PPM file.
 * Enter: DATA *data: pointer to rendered data.
 *        FILE *fptr: pointer to open file.  The file is not closed.
 * Exit:  int failed: 0 for okay, 1 for failed.                1/18/97-DWM */
{
  if (data->render.w<=0 || data->render.h<=0 || !data->render.buf)
    return(1);
  fprintf(fptr, "P6 %d %d 255\n", data->render.w, data->render.h);
  fwrite(data->render.buf, 1, data->render.w*data->render.h*3, fptr);
  return(0);
}

void update_geo(DATA *data)
/* Update the parameters necessary for raytracing.  This should only be
 *  called after orb_render has been called with process==6.
 * Enter: DATA *data: locked down data area for this window.   3/12/97-DWM */
{
  long w, h, w2, h2, i;

  w2 = w = 640;  h2 = h = 480;
  if (((data->dflag>>9)&1) && data->w>0 && data->h>0) {
    w2 = data->w;  h2 = data->h; }
  for (i=0; i<18; i++)  if (data->renderdlt[i]) break;
  if (i==18) {
    DefSize = data->mol.maxcheck;
    render_dlt_new(0, 0, w2, h2, data->pref.perspec, data->renderval,
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
}

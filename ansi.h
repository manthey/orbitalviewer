#define real   double
#define uchar  unsigned char
#define ulong  unsigned long
#define ushort unsigned short

#define ANSIC           1
#define NAMELEN      1024
#define ONENAMELEN    256
#define SUBVERSION      0
#define VERSION         1

#define BACKCOLOR       0
#define POSCOLOR        1
#define NEGCOLOR        2
#define ASYMCOLOR       3
#define PREVIEWCOLOR    4
#define NUMCOLORS       5

#define sq(x)  ((x)*(x))

#define DEFSIZE     27*a0
#define DIGISCALE    1e10
#define VRMLSCALE    1e10

#include "orb.h"

typedef struct PREF {
  long flags;
  float perspec;
  long colors[NUMCOLORS]; } PREF;
typedef struct DATA {
  long dflag;                         /* see comments after structure list */
  char name[ONENAMELEN];
  PREF pref;
  double renderdlt[54], renderval[10];          /* DLT center, left, right */
  long w, h;

  MOLECULE mol;
  CUTAWAY cut;
  ASYMPTOTE asym;
  POLYGON poly;
  OPOINTS points;
  RENDER render;
  STEREO stereo;

  long frame, lastframe, bezier, incr, seqtype, seqfps;
  char basename[ONENAMELEN];
  struct DATA *seq[4]; } DATA;

/* flags contains the following bit fields:
 *  Bit 0: 0-default to 4f0 orbital, 1-default to DEFAULT.ORB
 *      1: 0-open non-maximized window, 1-open maximized (not used by SG)
 *      2: 0-use diffuse color in vrml, 1-don't use diffuse color
 *      3: 0-use ambient color in vrml, 1-don't use ambient color
 *      4: 0-use emissive point color in vrml, 1-don't use emissive color
 *      5: --not used--
 *      6: 0-change only current window preferences, 1-change all windows
 *      7: 0-change global preferences, 1-do not change global preferences.
 *    8-9: --not used--
 *     10: 0-show splash screen, 1-do not show splash screen.
 * dflag contains the following bit fields:
 * Bit 0: 0-don't quick rendering method, 1-use quick rendering
 *   1-2: Quick render method: 0-point, 1-polygon, 2-raytrace
 *   3-4: Precise render method: 0-point, 1-polygon, 2-raytrace
 *   5-6: Current requested drawing mode: 0-point, 1-polygon, 2-raytrace
 *   7-8: Current actual drawing mode: 0-point, 1-polygon, 2-raytrace
 *     9: 0-use screen size, 1-use specified size.
 *    10: 0-not playing a sequence, 1-playing a sequence.                  */

void   camera_rotate  (double *renderval, double *renderdlt, double *phys,
                       long dir);
long   get_color      (DATA *data, long num);
double interpolate    (double *val, long *time);
int    main           (int argc, char *argv[]);
long   new_window_mid (DATA *data);
long   open_orb       (FILE *fptr, DATA *data, long noheader);
void   play_frame     (DATA *data, long frame);
void   prep_stereo    (DATA *data);
void   render_dlt     (long numnode, float *node, long w, long h,
                       float perspec, double *renderval, double *initdlt,
                       double *findlt);
void   render_dlt_new (long numnode, float *node, long w, long h,
                       float perspec, double *renderval, double *findlt);
void   render_move    (float *dx, float *dang, float *dx0, long w, long h,
                       double *renderval, double *initdlt, double *findlt);
void   render_physical(double *val, double *dlt, double *phys, lmatrix *ang);
long   save_ppm       (DATA *data, FILE *fptr);
long   save_vrml      (DATA *data);
void   save_vrml_color(FILE *fptr, char *text, long clr);
void   update_geo     (DATA *data);
void   update_process (DATA *data);

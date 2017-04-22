#ifndef OVHEADER                                 /* Prevent double inclusion */
#define OVHEADER 1

#define real   double
#define uchar  unsigned char
#define ulong  unsigned long
#define ushort unsigned short

#define DIRLEN         256
#define MAXLOCK        256
#define MAXWINDOWS     256
#define NAMELEN       1024
#define NUMMALLOC      500
#define ONENAMELEN     256
#define OPACDITHER    2137
#define PROGRESSDELAY    3
#define SMOOTHLOW      240
#define SUBVERSION       3
#define VERSION          1

#define BACKCOLOR        0
#define POSCOLOR         1
#define NEGCOLOR         2
#define ASYMCOLOR        3
#define PREVIEWCOLOR     4
#define NUMPREFCOLORS    5

#define sq(x)    ((x)*(x))

#define DEFSIZE      27*a0
#define DIGISCALE     1e10
#define VRMLSCALE     1e10

#define PREVIEWWIDTH   180
#define PREVIEWHEIGHT  130
#define PREVIEWNUM       8

#define FileOrb          0
#define FileOV           1
#define FilePPM          2
#define FileTIF          3
#define FileBMP          4
#define FileVRML         5
#define FileDigistar     6

#include "matrix.h"
#include "orb.h"

typedef struct PREF {
  long flags;
  long error, status, toolbar, tooltip, warn;
  long kgunit, munit, radunit;      /* 0-amu,1-kg,2-g,3-yg; 0-angstrom,1-nm,2-m,3-a; 0-rad,1-deg */
  float radstep, panstep, zoomstep;
  long pointsize;                           /* 0-1 pixel, 1-3 pixels, etc. */
  float perspec;
  long colors[NUMPREFCOLORS]; } PREF;
typedef struct DATA {
  long dflag;                         /* see comments after structure list */
  long changed, filetype;
  char name[ONENAMELEN];
  WINDOWPLACEMENT place;
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

  long frame, lastframe, incr, bezier, seqtype, seqfps;
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

#define debug(x) MessageBox(Hwnd, (x), "--", MB_OK)

BOOL CALLBACK asymptote_dialog   (HWND hdlg,ulong msg,WPARAM wp,LPARAM lp);
LRESULT CALLBACK main_loop       (HWND hdlg,ulong msg,WPARAM wp,LPARAM lp);
BOOL CALLBACK orbital_dialog     (HWND hdlg,ulong msg,WPARAM wp,LPARAM lp);
BOOL CALLBACK point_dialog       (HWND hdlg,ulong msg,WPARAM wp,LPARAM lp);
BOOL CALLBACK polygon_dialog     (HWND hdlg,ulong msg,WPARAM wp,LPARAM lp);
BOOL CALLBACK progress_dialog    (HWND hdlg,ulong msg,WPARAM wp,LPARAM lp);
BOOL CALLBACK render_dialog      (HWND hdlg,ulong msg,WPARAM wp,LPARAM lp);
BOOL CALLBACK render_opt_dialog  (HWND hdlg,ulong msg,WPARAM wp,LPARAM lp);
BOOL CALLBACK sequence_dialog    (HWND hdlg,ulong msg,WPARAM wp,LPARAM lp);
BOOL CALLBACK stereo_dialog      (HWND hdlg,ulong msg,WPARAM wp,LPARAM lp);
VOID CALLBACK timer              (HWND hwnd,ulong msg,ulong id,long time);

void     camera_rotate  (double *renderval, double *renderdlt, double *phys,
                         long dir);
void     free2          (void *mem);
double   interpolate    (double *val, long *time);
DATA    *lock_window    (HANDLE hnd);
void    *lock2          (HANDLE hmem);
void     play_frame     (DATA *data, long frame);
         progress       (HWND hwnd, char *title, char *message, long current,
                         long maximum);
real     rnd            (real low, real high, long inc);
void    *realloc2       (void *current, long size);
void     render_dlt     (long numnode, float *node, long w, long h,
                         float perspec, double *renderval, double *initdlt,
                         double *findlt);
void     render_dlt_new (long numnode, float *node, long w, long h,
                         float perspec, double *renderval, double *findlt);
void     render_move    (float *dx, float *dang, float *dx0, long w, long h,
                         double *renderval, double *initdlt, double *findlt);
void     render_physical(double *val, double *dlt, double *phys,
                         lmatrix *ang);
void     save_vrml_color(FILE *fptr, char *text, long clr);
HANDLE   unlock_window  (DATA *data);
HANDLE   unlock2        (void *mem);
void     update_process (DATA *data);
HWND     window_handle  (DATA *data);

extern char *DistText[], HelpFile[], lastview[], *MassText[], OrbLet[],
       Program[], *RadText[], Untitled[], WinName2[];
extern DATA *DispData;
extern float DefSize;
extern HPALETTE Hpal, HpalSplash;
extern HINSTANCE havi, hcom, hinst;
extern HMENU Hmenu;
extern HWND Hwnd, HwndC, HwndList[], HwndSplash, HwndStat, HwndTool;
extern long AVIHeader[], AVILoad, BitsPixel, Busy, CloseMode, EffScale,
       LockList[], NumLocked, NumWindows, Perspect, PreviewDelay,
       PreviewPic[], Speed, Splash[], *ToolCustom, ToolInit[];
extern PREF Pref;
extern real DistVal[], RadVal[];
extern TBBUTTON ToolList[], ToolSep;
extern WINDOWPLACEMENT WinPlace;
extern uchar *PreviewGraphic;

#endif                                    /* End of prevent double inclusion */

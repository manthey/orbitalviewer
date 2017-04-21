#ifndef ORBITALHEADER
#define ORBITALHEADER

#define real   double
#define uchar  unsigned char
#define ulong  unsigned long
#define ushort unsigned short

#define ONENAMELEN    256                         /* File system dependant */

#define sq(x)     ((x)*(x))
#define abs(x)    ((((x)>0)?(x):-(x)))
#ifndef max
  #define max(a,b)  (((a)>(b))?(a):(b))
#endif
#ifndef min
  #define min(a,b)  (((a)<(b))?(a):(b))
#endif

#define MAXN  30                    /* Maximum of 30 with double precision */
#define MAXP 100                    /* Max of    877 with long double      */
#define MAXL 100

#ifndef PI
  #define PI  3.14159265358979323846264338327950L
#endif
#define PI2   (PI/2)
#define deg   0.0174532925199
#define me    9.109389427e-31                      /* mass of electron, kg */
#define mz    1.6726230e-27                          /* mass of proton, kg */
#define mn    1.6749286e-27                         /* mass of neutron, kg */
#define ma    1.6605402e-27                        /* atomic unit mass, kg */
#define a0    5.29177249e-11              /* radius of first Bohr orbit, m */
#define nfac  2*PI*a0*a0*a0
#define sqrt2 1.41421356237309504876L

#ifdef ANSIC
  #define atan2l atan2
  #define cosl   cos
  #define expl   exp
  #define sinl   sin
  #define sqrtl  sqrt
#else
  double atan2l(double y, double x);
  double cosl  (double x);
  double expl  (double x);
  double sinl  (double x);
  double sqrtl (double x);
#endif
#define fabsl  fabs
#define powl   pow

#define StereoMONOSCOPIC    0
#define StereoSTEREOSCOPE   1
#define StereoINTERLACED    2
#define StereoREDBLUE       3
#define StereoSTEREOGRAM    4
#define StereoOVERLAY       5
#define StereoCHROMADEPTH   6
#define StereoMODES         7                       /* Max number of modes */

typedef struct OATOM {
  real mass;        /* Mass in kg (can be in ma, but will be changed by prep_atom) */
  long n, l, m, N, Z;   /* nlm=quantum numbers, N=num neutrons, Z=num protons */
  real x[3], ang[3];                           /* Location and orientation */
  real factor;          /* Multiplcation factor for atom.  Usually 1 or -1 */
  real mu, matp[9], matm[9];                  /* These are used internally */
  real cc, zna, crs[MAXN], cts[MAXN/2], crg[MAXN], ctg[MAXN/2], maxr;
  long nml1, am, lmam, lmam2, lmam12, identity;
  real lr, lts, ltc, lp, lrf, ltf, lpf, lpsi, lrzna, ltc2, lrle;
  real lstm; } OATOM;
typedef struct LIGHT {
  real x[3], lx[3], i, a;     /* vector, local vector, intensity, ambiance */
  real ux[3];                       /* unit vector, used only in rendering */
  long local; } LIGHT;
typedef struct MOLECULE {
  OATOM *orb;
  LIGHT *ls;
  long nump, numl;             /* number of atoms, number of light sources */
  real Psi;                        /* Psi^2 value to calculate surface for */
  long EffScale;          /* Approximate pixel width of destination screen */
  real maxjumpadj;            /* Quality number   1 is good, 1.1 is better */
  real minres;                /* Quality number  0.1 if no shadows, 0.004 if shadows.  Reduce to improve. */
  real back;                  /* 11.25/17.5/EffScale/2.4*(1 if no shadows or 0.04 if shadows) */
  real maxpsi;  /* Max psi^2 value in orbital.  Set for caller's reference */
  real checkpsi;            /* Psi value for which prep_check was last run */
  real maxcheck, maxcheck2, maxjump, syscntr[3], zsize, /* Used internally */
       lastpsi, lasthit; } MOLECULE;

typedef struct ASYMPTOTE {
  real opacity;                                 /* 0=no asymptote, 1=solid */
  long n, e, process, index[3];
  float *xyz;
  long wire;       /* rendering only: 0 to draw solid, 1 to draw wireframe */
  long *elem;
  long *scrxyz;        /* Used for rendering, not for orbital calculations */
  float *norm, *enorm; /* Used for rendering, not for orbital calculations */
  long density;         /* density=1 is 1 point across radius, 2 is denser */
  real checkpsi;
  real *rfac, *tfac, *pfac; } ASYMPTOTE;
typedef struct CUTAWAY {
  long type;                        /* 0-none, 1-planar, 2-corner, 3-wedge */
  long nosurface;            /* 0-cutaway produces a surface, 1-no surface */
  long invert; /* 0-cutaway is all positive octants, 1-any negative octant */
  long notlast;    /* 0-last cutaway in list, 1-additional cutaways follow */
  real xyz[3], ang[3], mat[9]; } CUTAWAY;
typedef struct OPOINTS {
  long n, maxn, process;   /* n==num points, xyz==coordinates for each point, phase==-1 or 1 for each point, maxn==number of points which will be generated, process is used internally */
  real checkpsi;
  float *xyz;
  long *scrxyz;        /* Used for rendering, not for orbital calculations */
  char *phase; } OPOINTS;
typedef struct POLYGON {
  long n, e, process, index[3];
  float *xyz;
  long *elem;
  long *scrxyz;        /* Used for rendering, not for orbital calculations */
  float *norm, *enorm; /* Used for rendering, not for orbital calculations */
  char *phase;
  long *ptphase;        /* This may be modified, provided the sign is kept */
  long density;         /* density=1 is 1 point across radius, 2 is denser */
  real refine;                          /* refine=0 for none, 4 is typical */
  real processmeter; /* Amount of refinement done. Purely for information. */
  real opacity[2];                   /* rendering only: positive, negative */
  long wire;       /* rendering only: 0 to draw solid, 1 to draw wireframe */
  real checkpsi;
  real *rfac, *tfac, *pfac; } POLYGON;
typedef struct RENDER {
  double camera[10], cammat[9];         /* X,Y,Z,Cx,Cy,x0,y0,theta,phi,psi */
  uchar color[12];                   /* +phase,-phase,background,asymptote */
  real opacity[7]; /* positive:overall, surface, interior; negative:o,s,i; asymptote */
  real refrac[3];                         /* positive, negative, asymptote */
  long steps;                                   /* # of steps per maxcheck */
  long autobright; /* 0 for none (use brightness), 1 to set, 2 to use last */
  real brightness;                       /* Multiplicand factor for colors */
  long antialias; /* bit 0:0-normal,1-coarse, 1:0-no anti,1-anti, 2-3:0-16 blocks,1-4 blocks,2-pixel blocks */
  real checkpsi;
  long process, x, y;
  long w, h, scan, type;                 /* type==0 for Graphic, 2 for BMP */
  long lastw, lasth, lasttype;               /* Used for memory efficiency */
  uchar *buf;  /* Pointer to image.  This can either be a Graphic or a BMP */
  uchar *image;                                /* Start of first scan line */
  ushort *zbuf;                            /* Filled for STEREOGRAM images */
  long *phase; } RENDER;                                    /* hyperphases */
typedef struct STEREO {
  long mode;                                     /* See StereoXX constants */
  long flags;                                                 /* See below */
  real interocular[2];             /* 0-all but stereoscope, 1-stereoscope */
  real separation;                             /* eye separation in meters */
  long w, h, pal;      /* width, height, and palette state of source image */
  uchar *image;                                 /* STEREOGRAM source image */
  char name[ONENAMELEN];                  /* path and name of source image */
  double lrcamera[20], lrcammat[18]; } STEREO;         /* Used with render */

/* STEREO flags:
/*  Bit 0: 0-top scan line is left, 1-top scan line is right (INTERLACED)
/*      1: 0-red is left, 1-red is right (REDBLUE)
/*      2: auto-separation flag, used by calling routines
/*      3: 0-random dot, 1-seed image (STEREOGRAM)                         */

real  calc_grad            (real *r, OATOM *p, long mag);
real  calc_grad_mag        (MOLECULE *mol);
void  calc_grad_total      (real *r, MOLECULE *mol);
real  calc_prob            (real *x, OATOM *p);
real  calc_prob_total      (real *x, MOLECULE *mol, long *phase);
real *calc_unit            (real *v);
real *cross                (real *v1, real *v2, real *r);
long  cut_intercept        (real *x, real *v, CUTAWAY *cut);
long  cut_zone             (real *x, CUTAWAY *cut, long zero);
long  intercept            (real *x, real *v, MOLECULE *mol, CUTAWAY *cut);
void  matrix_inverse       (real *mat, real *inv);
void  oangtomat            (double *r, double ta, double pa, double sa);
long  orb_asymptote        (ASYMPTOTE *as, CUTAWAY *cut, MOLECULE *mol,
                            float time, long process);
long  orb_points           (OPOINTS *pt, CUTAWAY *cut, MOLECULE *mol,
                            float time, long process);
long  orb_polygons         (POLYGON *pl, CUTAWAY *cut, MOLECULE *mol,
                            float time, long process);
long  orb_polygons_normal  (long *e, float *xyz, MOLECULE *mol);
long  orb_render           (RENDER *re, STEREO *st, CUTAWAY *cut,
                            MOLECULE *mol, float time, long process);
long  orb_render_ray       (RENDER *re, STEREO *st, CUTAWAY *cut,
                            MOLECULE *mol, real xx, real yy, long precise,
                            real *clrin, real *clrout);
void  orb_render_stereogram(RENDER *re, STEREO *st);
void  orb_render_vector    (double *c, double *cammat, real x, real y,
                            real zsize, real *xyz, real *v);
real  pow1                 (long minusone, long pow);
real  powr                 (real val, long pow);
void  prep_atom            (OATOM *p, real factor);
void  prep_check           (MOLECULE *mol);
long  quickint             (real *x, real *v, MOLECULE *mol, CUTAWAY *cut,
                            real *y);
long  ray                  (real *x, real *v, real *colorin, MOLECULE *mol,
                            CUTAWAY *cut, real *colorout);
real  ray_opacity          (real *x, real *v, long steps, RENDER *re,
                            MOLECULE *mol, CUTAWAY *cut);
long  ray_precise          (real *x, real *v, real *colorin, RENDER *re,
                            MOLECULE *mol, CUTAWAY *cut, real *colorout);
real  rnd                  (real low, real high, long inc);
void  snap_to_surface      (float *xyz, MOLECULE *mol, real psi2,
                            real uncertain);
real  to_jump              (void);
long  to_sphere            (real *x, real *v, real radius, real siz,
                            real *cntr);
real *transform            (real *x, real *t, real *a, real mag);
real *transform2           (real *x, real *t, real *m);

void  free2                (void *memblock);
void *malloc2              (size_t size);
void *realloc2             (void *memblock, size_t size);

#endif

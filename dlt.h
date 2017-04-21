/* This header file should accompany programs using the DLT photogrammetry
 *  equations.  There are two version of each routine.  The version preceeded
 *  by an 'l' use double values, the others use floats.        2/12/95-DWM */

long   dlt_to_physical   (float *dlt, float *phys, lmatrix *ang, long ten);
long   epipolar_line     (float *e12, lmatrix *P2P1i, float *xy,
                          float *abc);
long   epipolar_segment  (float *dlt1, float *dlt2, float *e12,
                          lmatrix *P2P1i, float *xy, float *bounds,
                          float *xyxy);
long   epipole           (float *dlt1, float *dlt2, float *e12,
                          lmatrix *P2P1i);
long   foe               (long n, float *xy1, float *xy2, long err,
                          float *foexy);
long   geo_to_dlt        (long n, float *xyz, float *xy, float *prel,
                          float *dlt, long distort);
void   geo_to_image      (long n, float *xyz, float *dlt, float *xy);
void   image_to_geo      (long n2, long p, long minp, short *np, float *xy,
                          float *dlt, float *rel, float *prel, long n,
                          float *xyz);
void   image2_to_geo     (long n, float *xy1, float *xy2, float *dlt1,
                          float *dlt2, float *xyz);
long   ldlt_to_physical  (double *dlt, double *phys, lmatrix *ang, long ten);
long   lepipolar_line    (double *e12, lmatrix *P2P1i, double *xy,
                          double *abc);
long   lepipolar_segment (double *dlt1, double *dlt2, double *e12,
                          lmatrix *P2P1i, double *xy, double *bounds,
                          double *xyxy);
long   lepipole          (double *dlt1, double *dlt2, double *e12,
                          lmatrix *P2P1i);
long   lfoe              (long n, double *xy1, double *xy2, long err,
                          double *foexy);
long   lgeo_to_dlt       (long n, double *xyz, double *xy, double *prel,
                          double *dlt, long distort);
void   lgeo_to_image     (long n, double *xyz, double *dlt, double *xy);
void   limage_to_geo     (long n2, long p, long minp, short *np, double *xy,
                          double *dlt, double *rel, double *prel, long n,
                          double *xyz);
void   limage2_to_geo    (long n, double *xy1, double *xy2, double *dlt1,
                          double *dlt2, double *xyz);
double lphoto_reliability(long n, double *xy, double *xyz, double *dlt);
long   lphysical_to_dlt  (double *phys, double *dlt, long focal);
void   lrotate_dlt       (double *dlt, double *org, double *axis,
                          double angle, double *findlt);
double photo_reliability (long n, float *xy, float *xyz, float *dlt);
long   physical_to_dlt   (float *phys, float *dlt, long focal);
void   rotate_dlt        (float *dlt, float *org, float *axis, float angle,
                          float *findlt);

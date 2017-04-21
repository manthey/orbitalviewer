/* This header file contains the matrix and lmatrix structures.  All commands
 *  working with double precision values are prefixed by a lower case L.
 *                                                             8/17/94-DWM */

typedef struct matrix {
  long w, h;
  float *m; } matrix;

typedef struct lmatrix {
  long w, h;
  double *m; } lmatrix;

#ifndef PI
  #define PI 3.14159265358979323846264338327950
#endif

void     abgtomat      (float *r, float alpha, float beta, float gamma);
void     angtomat      (float *r, float ta, float pa, float sa);
float   *cross_product (float *v1, float *v2, float *r);
long     cubic         (float *coef, float *roots);
float    dot_product   (float *v1, float *v2);
void     eigen         (float *m, float *l);
long     eigenvect     (float *m, float *l, float *ev, long column);
void     labgtomat     (double *r, double alpha, double beta, double gamma);
void     langtomat     (double *r, double ta, double pa, double sa);
double  *lcross_product(double *v1, double *v2, double *r);
long     lcubic        (double *coef, double *roots);
double   ldot_product  (double *v1, double *v2);
long     least         (matrix *a, float *coef, float weight);
void     leigen        (double *m, double *l);
long     leigenvect    (double *m, double *l, double *ev, long column);
long     lleast        (lmatrix *a, double *coef, double weight);
lmatrix *lmat          (double *a);
long     lmatadd       (double c, lmatrix *a, double d, lmatrix *b,
                        lmatrix *caplusdb);
double   lmatdet       (lmatrix *a);
long     lmatdup       (lmatrix *a, lmatrix *adup);
void     lmatident     (lmatrix *i, long size);
long     lmatinv       (lmatrix *a, lmatrix *ainv);
long     lmatmul       (lmatrix *a, lmatrix *b, lmatrix *amulb);
void     lmatscale     (lmatrix *a, double scalar, lmatrix *result);
long     lmatsub       (lmatrix *a, lmatrix *b, lmatrix *aminusb);
long     lmatsym       (lmatrix *a);
void     lmattoabg     (double *r, double *abg);
void     lmattoang     (double *r, double *a);
void     lmattoypr     (double *r, double *a);
void     lmattrans     (lmatrix *a, lmatrix *atrans);
void     lmatzero      (lmatrix *a);
double  *lqconj        (double *a, double *astar);
double   lqdot         (double *a, double *b);
double  *lqmul         (double *a, double *b, double *amulb);
double   lqnorm        (double *a);
void     lqtomat       (double *q, double *rot);
double  *lqunit        (double *a, double *aunit);
void     lrotate_fit   (double *XYZ, double *xyz, long fitstart, long fitnum,
                        long modstart, long modnum, long scale, double *rot,
                        double *tran);
long     lrowreduce    (lmatrix *a);
long     lrowreduce2   (lmatrix *a, lmatrix *b);
double  *lunit         (double *v, double *r);
void     lyprtomat     (double *rot, double y, double p, double r);
matrix  *mat           (float *a);
long     matadd        (float c, matrix *a, float d, matrix *b,
                        matrix *caplusdb);
double   matdet        (matrix *a);
long     matdup        (matrix *a, matrix *adup);
void     matident      (matrix *i, long size);
long     matinv        (matrix *a, matrix *ainv);
long     matmul        (matrix *a, matrix *b, matrix *amulb);
void     matscale      (matrix *a, float scalar, matrix *result);
long     matsub        (matrix *a, matrix *b, matrix *aminusb);
long     matsym        (matrix *a);
void     mattoabg      (float *r, float *abg);
void     mattoang      (float *r, float *a);
void     mattoypr      (float *r, float *a);
void     mattrans      (matrix *a, matrix *atrans);
void     matzero       (matrix *a);
float   *qconj         (float *a, float *astar);
float    qdot          (float *a, float *b);
float   *qmul          (float *a, float *b, float *amulb);
float    qnorm         (float *a);
void     qtomat        (float *q, float *rot);
float   *qunit         (float *a, float *aunit);
void     rotate_fit    (float *XYZ, float *xyz, long fitstart, long fitnum,
                        long modstart, long modnum, long scale, float *rot,
                        float *tran);
long     rowreduce     (matrix *a);
long     rowreduce2    (matrix *a, matrix *b);
float   *unit          (float *v, float *r);
void     yprtomat      (float *rot, float y, float p, float r);

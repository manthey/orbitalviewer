#ifndef DRAWHEADER                               /* Prevent double inclusion */
#define DRAWHEADER 1

#include <windows.h>
#include "ov.h"

void   bgr_to_rgb    (uchar *image, long siz);
void   camera_figure (HWND hdlg, HDC hdc, DATA *data, double *phys);
void   camera_figure_draw(long w, long h, long scan, uchar *buf, DATA *data,
                      double *phys);
long   camera_move   (HWND hdlg, long delx, long dely, long w, long h,
                      long down, DATA *data, long ptr);
void   camera_rect   (RECT *rect, HWND hdlg);
void   cutaway_figure(HWND hdlg, HDC hdc, real *pos, long type);
void   cutaway_figure_draw(long w, long h, long scan, uchar *buf, real *pos,
                      long type, long color);
void   cutaway_figure_point(real *x, real *v, real *pos, long color,
                      DATA *data, uchar *dest);
long   cutaway_inter (real *x, real *v, real *pos, DATA *data, real *inter);
long   cutaway_move  (HWND hdlg, long delx, long dely, long w, long h,
                      long down, DATA *data, real *pos);
void   cutaway_rect  (RECT *rect, HWND hdlg);
void   draw_bmp      (HDC hdc, long x, long y, HANDLE bmp);
void   draw_bmp2     (HDC hdc, long x, long y, uchar *bmp);
void   fill_norm     (DATA *data);
void   fill_zone24   (uchar *dest, long bclr, long count);
void   focal_geo     (long axis, float jump);
void   light_figure  (HWND hdlg, HDC hdc, LIGHT *ls);
void   light_figure_draw(long w, long h, long scan, uchar *buf, LIGHT *ls,
                      long color);
long   light_move    (HWND hdlg, long delx, long dely, long w, long h,
                      long down);
void   light_rect    (RECT *rect, HWND hdlg);
void   make_palette  (void);
HPALETTE make_windows_palette(uchar *pal, long bgr);
void   mouse         (HWND hwnd, long but, long up, long flag, long x, long y);
void   mouse_move    (HWND hwnd, long up, long flag, long x, long y);
void   opacity_dither(void);
long   preview_mouse (HWND hdlg, ulong msg, WPARAM wp, LPARAM lp, DATA *data,
                      long num, long ptr);
HBITMAP reduce_color_space(HBITMAP wpic, HDC hdc);
void   reframe_geo   (void);
void   render_asymptote(long w, long h, long scan, uchar *scrbuf, ushort *zbuf,
                      DATA *data);
void   render_coor   (DATA *data, double *dlt, long **scrxyz, float *node,
                      long numnode, long **scrxyz2, float *node2,
                      long numnode2);
void   render_line   (long w, long h, long scan, uchar *scrbuf, long zscan,
                      ushort *zbuf, long x0, long y0, long z0, long x1,
                      long y1, long z1, long color);
void   render_point  (long w, long h, long scan, uchar *scrbuf, ushort *zbuf,
                      long x, long y, long z, long size, long clr, long bgr);
void   render_points (long w, long h, long scan, uchar *scrbuf, ushort *zbuf,
                      DATA *data, long bgr);
void   render_polygon(long w, long h, long scan, uchar *scrbuf, ushort *zbuf,
                      DATA *data);
void   render_sort_coor(long *coor, long *x, long s, long d);
void   render_triangle_opac(long w, long h, long scan, uchar *scrbuf,
                      long zscan, ushort *zbuf, long coor[15], uchar *color);
void   render_triangle_texture(long w, long h, long scan, uchar *scrbuf,
                      long zscan, ushort *zbuf, long coor[15], long pw,
                      long ph, long ppal, uchar *photo, long bgr, long scale,
                      float bright);
void   reset_geo     (void);
void   rotate_geo    (long axis, float ang);
void   scale_pic24   (uchar *dest, uchar *source, short srcw, short destw,
                      short simgw, short dimgw, short simgh, short dimgh);
void   scale_zbuf    (ushort *dest, ushort *source, short simgw, short dimgw,
                      short simgh, short dimgh);
void   shift_geo     (long axis, float jump);
uchar *smooth_palette(DATA *data, long numcolors);
long   sphere_inter  (real *x, real *v, real radius, real *inter);
void   stereo_figure (HWND hdlg, long mode, long swap, STEREO *st);
void   stereogram    (long w, long h, long scan, uchar *scrbuf, ushort *zbuf,
                      DATA *data);
uchar *update_bmp    (DATA *data, long w, long h, long *freeimg, long bits);
void   update_geo    (HWND hwnd, DATA *data, HDC dc);
long   update_window (HWND hwnd);
long   use_palette   (uchar *pal);

extern char *DispInfo[];
extern long LastMove;
extern uchar *ColorTable;

#endif                                    /* End of prevent double inclusion */

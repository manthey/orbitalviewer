#ifndef FILEHEADER                               /* Prevent double inclusion */
#define FILEHEADER 1

#include <windows.h>
#include <vfw.h>

#include "ov.h"

BOOL CALLBACK check_dialog      (HWND hdlg, ulong msg, WPARAM wp, LPARAM lp);
BOOL CALLBACK color_dialog      (HWND hdlg, ulong msg, WPARAM wp, LPARAM lp);
VOID CALLBACK hide_splash       (HWND hwnd, ulong msg, ulong id, long time);
BOOL CALLBACK preferences_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp);
BOOL CALLBACK save_dialog       (HWND hdlg, ulong msg, WPARAM wp, LPARAM lp);
BOOL CALLBACK splash_dialog     (HWND hdlg, ulong msg, WPARAM wp, LPARAM lp);

void   add_extensions   (HWND hwnd);
long   close_window     (HWND hwnd);
void   compress_avi     (void);
void   compress_avi_stream(HWND hwnd, PAVISTREAM stream, AVISTREAMINFO *info,
                         PAVIFILE cavi);
void   copy             (void);
void   drop_file        (HANDLE drop);
char  *find_space       (char *text);
void   join_names       (char *text);
void   new_window       (void);
long   open_ov          (FILE *fptr, DATA *data);
void   open_window      (void);
void   open_window_mid  (void);
void   read_ini         (void);
void   read_ini_defaults(void);
void   read_ini_toolbar (long num);
void   recheck          (HWND hwnd, long refresh);
void   remove_quotes    (char *text);
long   save_avi         (DATA *data);
long   save_graphic     (DATA *data, long type);
long   save_orb         (DATA *data);
void   save_orb_write   (DATA *data, FILE *fptr, long mode);
long   save_ov          (DATA *data);
long   save_window      (HWND hwnd, long saveas);
void   sequence_browse  (HWND hdlg);
DATA  *sequence_file    (HWND hdlg, DATA *seq);
void   set_default      (void);
void   setup            (void);
void   show_splash      (void);
void   status           (char *text);
void   stereo_image     (HWND hdlg, STEREO *stereo, long *replace);
void   strpath          (char *dest, char *src);
long   windows_file     (char *name);
void   write_ini        (void);
void   write_ini_toolbar(long numt);

void   camera_rotate  (double *renderval, double *renderdlt, double *phys,
                       long dir);
long   get_color      (DATA *data, long num);
double interpolate    (double *val, long *time);
long   new_window_mid (DATA *data);
long   open_orb       (FILE *fptr, DATA *data, long noheader);
void   play_frame     (DATA *data, long frame);
void   prep_stereo    (DATA *data);
void   render_dlt     (long numnode, float *node, long w, long h,
                       float perspec, double *renderval, double *initdlt,
                       double *findlt);
void   render_dlt_new (long numnode, float *node, long w, long h,
                       double perspec, double *renderval, double *findlt);
void   render_move    (float *dx, float *dang, float *dx0, long w, long h,
                       double *renderval, double *initdlt, double *findlt);
void   render_physical(double *val, double *dlt, double *phys, lmatrix *ang);
long   save_digistar  (DATA *data);
long   save_vrml      (DATA *data);
void   save_vrml_color(FILE *fptr, char *text, long clr);
void   update_process (DATA *data);

extern char OpenName[];

#endif                                    /* End of prevent double inclusion */

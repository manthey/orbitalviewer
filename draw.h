#ifndef DRAWHEADER                               /* Prevent double inclusion */
#define DRAWHEADER 1

#include <windows.h>
#include "ov.h"

void   stereo_figure(HWND hdlg, long mode, long swap, STEREO *st);
void   stereogram   (long w, long h, long scan, uchar *scrbuf, ushort *zbuf,
                     DATA *data);
uchar *update_bmp   (DATA *data, long w, long h, long *freeimg, long bits);
void   update_geo   (HWND hwnd, DATA *data, HDC dc);
long   update_window(HWND hwnd);
long   use_palette  (uchar *pal);

#endif                                    /* End of prevent double inclusion */

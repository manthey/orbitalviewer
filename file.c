/* This contains all file access functions, preference functions, and
/*  clipboard functions for the Orbital Viewer program.         1/3/98-DWM */

#include <windows.h>
#include <commctrl.h>
#include <vfw.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "matrix.h"
#include "dlt.h"
#include "ovrc.h"
#include "ov.h"
#include "\p\lib\gvlib.h"

#include "common.c"

char OpenDir[DIRLEN*7], OpenName[NAMELEN], opfilter[NAMELEN];
HBITMAP SplashBMP=0;
long CustColor[16], FilterType[7]={1,1,1,1,1,1,1};
OPENFILENAME openfile;
char *extlist[]={"ORB","OV",0};
char DefFile[262]="DEFAULT.ORB", inifile[256]="OV.INI",
     *inihead[]={"Orbital Viewer","Toolbar"};
char No[]="No", Yes[]="Yes", SZero[]="";
char *key[]={"WindowPosition",SZero,0,"ErrorMessages",Yes,No,"WarningMessages",Yes,No,"StatusBar",Yes,No,"ToolBar",Yes,No,"ToolTip",Yes,No,"CustomTools",SZero,0,"OpenDirectory",SZero,0,"OpenSaveFilter",SZero,0,
     "MainFlags",SZero,0,"DefaultPerspective",SZero,0,"Colors",SZero,0,"LastView1",SZero,0,"LastView2",SZero,0,"LastView3",SZero,0,"LastView4",SZero,0,"SaveDirectory",SZero,0,"StereogramDirectory",SZero,0,"RotationStep",SZero,0,"PanStep",SZero,0,
     "ZoomStep",SZero,0,"PointSize",SZero,0,"AngleUnits",SZero,0,"DistanceUnits",SZero,0,"MassUnits",SZero,0,"SequenceOpenDirectory",SZero,0,"SequenceSaveDirectory",SZero,0,"OpenAVIDirectory",SZero,0,"SaveAVIDirectory",SZero,0, 0};
char OpenFilter[]="All files|*.*|Orbital specification files (*.ORB)|*.ORB|"
     "Orbital viewer files (*.OV)|*.OV|";
char OpenFilter2[]="All files|*.*|Graphics files|*.BMP;*.CUR;*.DIB;*.GRB;"
     "*.GRO;*.ICO;*.JIF;*.JFF;*.JPE;*.JPG;*.PBM;*.PCC;*.PCX;*.PGM;*.PNM;"
     "*.PPM;*.TGA;*.TIF|";
char OpenFilter3[]="All files|*.*|AVI files (*.AVI)|*.AVI|";
char OrbCheck[]="OrbitalFileV1.0", OrbCheck2[]="OrbitalViewerFileV1.0";
char SaveFilter[]="Orbital specifcation file (*.ORB)|*.ORB|"
     "Orbital viewer file (*.OV)|*.OV|";
char SaveFilterA[]="Portable pixel map file (*.PPM)|*.PPM|TIFF file (*.TIF)|"
     "*.TIF|Windows bitmap file (*.BMP)|*.BMP|";
char SaveFilterB[]="VRML 1.0 file (*.WRL)|*.WRL|";
char SaveFilterC[]="AVI file (*.AVI)|*.AVI|";
char SaveFilterD[]="Digistar file (*.TXT)|*.TXT|";

compress_avi_stream(HWND hwnd, PAVISTREAM stream, AVISTREAMINFO *info,
                    PAVIFILE cavi);

add_extensions(hwnd)
/* Enter: HWND hwnd: handle of window.                         5/27/96-DWM */
{
  char ext[]="Extensions", base[256], out[256], *cmd, ekey[5];
  long i, j;
  HKEY hkey;

  if (MsgBox(hwnd, "This function will associate all appropriate file\n"
                   "extensions with the Orbital Viewer program.  When an\n"
                   "ORB file is run, Orbital Viewer will automatically\n"
                   "be used to load it.", "Link Extensions",
      MB_OKCANCEL)!=IDOK)
    return(0);
  cursor(1);
  cmd = GetCommandLine();
  if (strlen(cmd)<245)  strcpy(base, cmd);  else  strcpy(base, Program);
  if (find_space(base))  find_space(base)[0] = 0;
  for (i=0; extlist[i]; i++) {
    sprintf(out, "%s ^.%s", base, extlist[i]);
    WriteProfileString(ext, extlist[i], out); }
  sprintf(out, "%s \"%%1\"", base);
  RegCreateKeyEx(HKEY_CLASSES_ROOT, "OV\\shell\\open\\command", 0, "",
                 REG_OPTION_NON_VOLATILE, KEY_ALL_ACCESS, 0, &hkey, &j);
  RegSetValueEx(hkey, "", 0, REG_SZ, out, strlen(out)+1);
  RegCloseKey(hkey);
  RegCreateKeyEx(HKEY_CLASSES_ROOT, "OV\\DefaultIcon", 0, "",
                 REG_OPTION_NON_VOLATILE, KEY_ALL_ACCESS, 0, &hkey, &j);
  strcpy(out, base);  remove_quotes(out);
  RegSetValueEx(hkey, "", 0, REG_SZ, out, strlen(out)+1);
  RegCloseKey(hkey);
  RegCreateKeyEx(HKEY_CLASSES_ROOT, Program, 0, "",
                 REG_OPTION_NON_VOLATILE, KEY_ALL_ACCESS, 0, &hkey, &j);
  RegSetValueEx(hkey, "", 0, REG_SZ, "Orbital Viewer File", 14);
  RegCloseKey(hkey);
  for (i=0; extlist[i]; i++) {
    sprintf(out, ".%s", extlist[i]);
    RegCreateKeyEx(HKEY_CLASSES_ROOT, out, 0, "", REG_OPTION_NON_VOLATILE,
                   KEY_ALL_ACCESS, 0, &hkey, &j);
    RegSetValueEx(hkey, "", 0, REG_SZ, "OV", 3);
    RegCloseKey(hkey); }
  cursor(0);
  MsgBox(hwnd, "Extensions Linked", "Link Extensions", MB_OK);
}

BOOL CALLBACK check_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the file has changed check dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.        5/26/96-DWM */
{
  HWND ctrl;
  char text[ONENAMELEN+80], name[ONENAMELEN];

  switch (msg) {
    case WM_COMMAND: EndDialog(hdlg, wp&0xFFFF); break;
    case WM_INITDIALOG: if (CloseMode!=1) {
        ctrl = GetDlgItem(hdlg, CheckAll);   DestroyWindow(ctrl);
        ctrl = GetDlgItem(hdlg, CheckNone);  DestroyWindow(ctrl); }
      GetWindowText(SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0), name,
                    ONENAMELEN-1);
      if (!strcmp(text, Untitled))
        strcpy(text, "This file has not been saved.  Save it?");
      else
        sprintf(text, "The file %s has changed.  Save these changes?", name);
      SetDlgItemText(hdlg, CheckText, text);
      SetFocus(hdlg);
      return(1); }
  return(0);
}

long close_window(HWND hwnd)
/* Close a window.  If the data in the window has changed, ofter to save it
 *  before quiting.  Handle cancelling as well as saving.
 * Enter: HWND hwnd: handle of window to close.
 * Exit:  long abort: 0 for close the window, 1 for cancel.    2/13/97-DWM */
{
  DATA *data;
  long res;
  char text[ONENAMELEN+80];

  if (hwnd==SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0)) {
    if (IsZoomed(hwnd))        Pref.flags |= 2;
    else if (!IsIconic(hwnd))  Pref.flags &= 0xFFFFFFFD; }
  data = lock_window(hwnd);
  GetWindowPlacement(Hwnd, &data->place);
  if (!data->changed) {
    unlock_window(data);
    free_window(hwnd);
    return(0); }
  unlock_window(data);
  switch (CloseMode) {
    case 0: case 1: Busy++;
      res = DialogBox(hinst, CheckDLG, hwnd, check_dialog);
      Busy--;
      if (res==CheckAll) {  CloseMode = 3;  res = IDYES; }
      if (res==CheckNone) { CloseMode = 2;  res = IDNO; } break;
    case 2: res = IDNO; break;
    case 3: res = IDYES; }
  if (res==IDCANCEL)  return(1);
  if (res==IDYES)
    if (!save_window(hwnd, 0))
      return(1);
  free_window(hwnd);
  return(0);
}

BOOL CALLBACK color_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the Colors dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.         1/8/98-DWM */
{
  HWND hwnd;
  DATA *data=0, *data3;
  static DATA data2;
  char text[80];
  long *clr, box[]={80,6,64,18, 24}, out[4], i, j, k, c, c2, old[NUMPREFCOLORS];
  RECT rect;
  PAINTSTRUCT paint;
  HDC hdc;
  HBRUSH oldb, new;
  CHOOSECOLOR cc;

  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case ColorClr1: case ColorClr2: case ColorClr3: case ColorClr4:
      case ColorClr5: cc.lStructSize = sizeof(CHOOSECOLOR);
        cc.hwndOwner = hdlg;  cc.hInstance = hinst;
        cc.rgbResult = get_color(&data2, (wp&0xFFFF)-ColorClr1);
        bgr_to_rgb(&cc.rgbResult, 1);
        cc.lpCustColors = CustColor;
        cc.Flags = CC_FULLOPEN|CC_RGBINIT;  cc.lCustData = 0;
        if (ChooseColor(&cc)) {
          i = cc.rgbResult;  bgr_to_rgb(&i, 1);
          data2.pref.colors[(wp&0xFFFF)-ColorClr1] = i;
          i = (wp&0xFFFF)-ColorClr1;
          map_dialog_rect(hdlg, box[0],box[1]+box[4]*i, box[2],box[3], out);
          rect.left = out[0]-2;   rect.top = out[1]-2;
          rect.right = out[2]+out[0]+2;  rect.bottom = out[3]+out[1]+2;
          InvalidateRect(hdlg, &rect, 1); } break;
      case ColorReset: memset(data2.pref.colors, -1, NUMPREFCOLORS*sizeof(long));
        map_dialog_rect(hdlg, box[0], box[1], box[0]+box[2],
                        box[1]+box[3]+(box[4]*NUMPREFCOLORS-1), out);
        rect.left = out[0]-2;   rect.top = out[1]-2;
        rect.right = out[2]+2;  rect.bottom = out[3]+2;
        InvalidateRect(hdlg, &rect, 1); break;
      case HelpHelp: WinHelp(Hwnd,HelpFile,HELP_CONTEXT,HelpColor); break;
      case IDOK: clr = Pref.colors;
        if (hwnd=SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0))
          if (data=lock_window(hwnd))
            clr = data->pref.colors;
        Pref.flags &= 0xFFFFFF3F;
        if (SendDlgItemMessage(hdlg,ColorAll,BM_GETCHECK,0,0))
          Pref.flags |= (1<<6);
        if (SendDlgItemMessage(hdlg,ColorPref,BM_GETCHECK,0,0))
          Pref.flags |= (1<<7);
        memcpy(old, clr, NUMPREFCOLORS*sizeof(long));
        if (memcmp(old, data2.pref.colors, NUMPREFCOLORS*sizeof(long))) {
          if ((Pref.flags>>7)&1)
            for (i=0; i<NUMPREFCOLORS; i++)
              if (old[i]!=data2.pref.colors[i])
                Pref.colors[i] = data2.pref.colors[i];
          if ((Pref.flags>>6)&1)
            for (i=0; i<NumWindows; i++) {
              data3 = lock_window(HwndList[i]);
              for (j=0; j<NUMPREFCOLORS; j++)
                if (old[j]!=data2.pref.colors[j])
                  data3->pref.colors[j] = data2.pref.colors[j];
              out[0] = get_color(data3, POSCOLOR);
              out[1] = get_color(data3, NEGCOLOR);
              out[2] = get_color(data3, BACKCOLOR);
              out[3] = get_color(data3, ASYMCOLOR);
              for (k=0; k<4; k++)  for (j=0; j<3; j++)
                data3->render.color[k*3+j] = ((uchar *)(out+k))[2-j];
              update_process(data3);
              unlock_window(data3);
              InvalidateRect(HwndList[i], 0, 0); }
          memcpy(clr, data2.pref.colors, NUMPREFCOLORS*sizeof(long));
          out[0] = get_color(data, POSCOLOR);
          out[1] = get_color(data, NEGCOLOR);
          out[2] = get_color(data, BACKCOLOR);
          out[3] = get_color(data, ASYMCOLOR);
          for (i=0; i<4; i++)  for (j=0; j<3; j++)
            data->render.color[i*3+j] = ((uchar *)(out+i))[2-j];
          update_process(data); }
        unlock_window(data);
      case IDCANCEL: EndDialog(hdlg, 1);  return(1);
      default: return(0); } break;
    case WM_INITDIALOG: clr = Pref.colors;
      if (hwnd=SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0))
        if (data=lock_window(hwnd))
          clr = data->pref.colors;
      memset(&data2, 0, sizeof(DATA));
      memcpy(data2.pref.colors, clr, NUMPREFCOLORS*sizeof(long));
      SendDlgItemMessage(hdlg, ColorAll, BM_SETCHECK, (Pref.flags>>6)&1, 0);
      SendDlgItemMessage(hdlg, ColorPref, BM_SETCHECK, (Pref.flags>>7)&1, 0);
      map_dialog_rect(hdlg, box[0], box[1], box[0]+box[2],
                      box[1]+box[3]+(box[4]*NUMPREFCOLORS-1), out);
      rect.left = out[0]-2;   rect.top = out[1]-2;
      rect.right = out[2]+2;  rect.bottom = out[3]+2;
      InvalidateRect(hdlg, &rect, 1);
      unlock_window(data);
      SetFocus(hdlg);
      return(1);
    case WM_PAINT: hdc = BeginPaint(hdlg, &paint);
      for (i=0; i<NUMPREFCOLORS; i++) {
        map_dialog_rect(hdlg, box[0], box[1]+box[4]*i, box[2], box[3], out);
        c2 = get_color(&data2, i);
        c = (c2>>16)+(c2&0xFF00)+((c2&0xFF)<<16);
        if ((c>>16)>=0x80 || ((c>>8)&0xFF)>=0x80 || (c&0xFF)>=0x80)
          SelectObject(hdc, GetStockObject(BLACK_PEN));
        else
          SelectObject(hdc, GetStockObject(WHITE_PEN));
        new = CreateSolidBrush(c);
        oldb = SelectObject(hdc, new);
        Rectangle(hdc, out[0], out[1], out[2]+out[0], out[3]+out[1]);
        SelectObject(hdc, oldb);
        DeleteObject(new); }
      EndPaint(hdlg, &paint); }
  return(0);
}

compress_avi()
/* Present a dialog to select an AVI, then ask for compression options, then
 *  for a save name, then compress it.                         1/21/98-DWM */
{
  long i, len;
  FILE *fptr;
  char buf[12], oname[ONENAMELEN], cname[ONENAMELEN];
  PAVIFILE oavi=0, cavi=0;
  PAVISTREAM stream=0, sarray[1];
  AVICOMPRESSOPTIONS opt;
  LPAVICOMPRESSOPTIONS oarray[1];
  AVISTREAMINFO info;

  start_avi();
  if (!havi) {
    recheck(0, 0);
    error(Hwnd, "Can't access Windows AVI function calls.");  return(0); }
  strcpy(opfilter, OpenFilter3);
  len = strlen(opfilter);
  for (i=0; i<len; i++)  if (opfilter[i]=='|')  opfilter[i] = 0;
  openfile.lStructSize = sizeof(OPENFILENAME);
  openfile.hwndOwner = Hwnd;  openfile.hInstance = hinst;
  openfile.lpstrFilter = opfilter;
  openfile.nFilterIndex = FilterType[5];
  openfile.lpstrFile = OpenName;
  openfile.nMaxFile = NAMELEN-1;
  openfile.lpstrFileTitle = openfile.lpstrCustomFilter = 0;
  openfile.lpstrInitialDir = OpenDir+DIRLEN*5;
  openfile.Flags = OFN_FILEMUSTEXIST|OFN_HIDEREADONLY|OFN_PATHMUSTEXIST|
                   OFN_EXPLORER;
  memset(OpenName, 0xFF, NAMELEN);  OpenName[0] = 0;
  if (!GetOpenFileName(&openfile))  return(0);
  FilterType[5] = openfile.nFilterIndex;
  if (!(fptr=fopen(OpenName, "rb"))) {
    error(Hwnd, "Can't read file.");  return(0); }
  if (strchr(OpenName, '\\'))
    strpath(OpenDir+DIRLEN*5, OpenName);
  fread(buf, 1, 12, fptr);
  fclose(fptr);
  if (memcmp(buf, "RIFF", 4) || memcmp(buf+8, "AVI ", 4)) {
    error(Hwnd, "File is not a recognized AVI file.");  return(0); }
  strcpy(oname, OpenName);
  strcpy(opfilter, SaveFilterC);
  len = strlen(opfilter);
  for (i=0; i<len; i++)  if (opfilter[i]=='|')  opfilter[i] = 0;
  openfile.lpstrFilter = opfilter;
  openfile.nFilterIndex = FilterType[6];
  openfile.lpstrFileTitle = openfile.lpstrCustomFilter = 0;
  openfile.lpstrInitialDir = OpenDir+DIRLEN*6;
  openfile.lpstrTitle = 0;
  openfile.lpfnHook = save_dialog;
  openfile.Flags = OFN_ENABLEHOOK|OFN_HIDEREADONLY|OFN_PATHMUSTEXIST|
                   OFN_EXPLORER|(OFN_OVERWRITEPROMPT*Pref.warn);
  memset(OpenName, 0xFF, NAMELEN);  OpenName[0] = 0;
  if (!GetSaveFileName(&openfile))  return(0);
  if (!stricmp(OpenName, oname)) {
    error(Hwnd, "Destination file can not be the same as source file.");
    return(0); }
  FilterType[6] = openfile.nFilterIndex;
  if (strchr(OpenName, '\\'))
    strpath(OpenDir+DIRLEN*6, OpenName);
  strcpy(cname, OpenName);
  end_ctl3d();              /* If CTL3D is enabled, AVISaveOptions crashes */
  if ((*GetProcAddress(havi,"AVIFileOpen"))(&oavi, oname, OF_READ, 0)) {
    error(Hwnd, "Couldn't open AVI file using Windows AVI function.");
    start_ctl3d();  return(0); }
  remove(cname);                                /* delete destination file */
  if ((*GetProcAddress(havi,"AVIFileOpen"))(&cavi, cname, OF_CREATE|OF_WRITE, 0)) {
    error(Hwnd, "Couldn't create destination file.");
    (*GetProcAddress(havi,"AVIFileRelease"))(oavi);  start_ctl3d();
    return(0); }
  if ((*GetProcAddress(havi,"AVIStreamOpenFromFile"))(&stream, oname,
      streamtypeVIDEO, 0, OF_READ,0)) {
    error(Hwnd, "Couldn't find source video.");
    (*GetProcAddress(havi,"AVIFileRelease"))(oavi);
    (*GetProcAddress(havi,"AVIFileRelease"))(cavi);
    start_ctl3d();  return(0); }
  sarray[0] = stream;
  oarray[0] = &opt;
  memset(oarray[0], 0, sizeof(AVICOMPRESSOPTIONS));
  if (!((*GetProcAddress(havi,"AVISaveOptions"))(Hwnd,
      ICMF_CHOOSE_KEYFRAME|ICMF_CHOOSE_DATARATE|ICMF_CHOOSE_PREVIEW, 1,
      sarray, oarray))) {
    (*GetProcAddress(havi,"AVIFileRelease"))(oavi);
    (*GetProcAddress(havi,"AVIFileRelease"))(cavi);
    start_ctl3d();  return(0); }
  stream = 0;
  if ((*GetProcAddress(havi,"AVIMakeCompressedStream"))(&stream, sarray[0],
      oarray[0], 0)==AVIERR_OK) {
    (*GetProcAddress(havi,"AVIStreamInfo"))(stream, &info,
                     sizeof(AVISTREAMINFO));
    compress_avi_stream(Hwnd, stream, &info, cavi);
    if (stream)  (*GetProcAddress(havi,"AVIStreamRelease"))(stream); }
  else {
    (*GetProcAddress(havi,"AVIStreamInfo"))(sarray[0], &info,
                     sizeof(AVISTREAMINFO));
    compress_avi_stream(Hwnd, sarray[0], &info, cavi); }
  (*GetProcAddress(havi,"AVISaveOptionsFree"))(1, oarray);
  (*GetProcAddress(havi,"AVIFileRelease"))(oavi);
  (*GetProcAddress(havi,"AVIFileRelease"))(cavi);
  start_ctl3d();
}

compress_avi_stream(HWND hwnd, PAVISTREAM stream, AVISTREAMINFO *info,
                    PAVIFILE cavi)
/* Actually copy the original file to the compressed file.
 * Enter: HWND hwnd: handle of calling window.
 *        PAVISTREAM: pointer to input stream.
 *        AVISTREAMINFO *info: stream information block.
 *        PAVIFILE cavi: essentially *fptr for output file.    1/20/98-DWM */
{
  long fsize, ssize, stsize, pos, read, *buf, bsize, *buf2, end, start;
  PAVISTREAM new;
  char text[80];

  (*GetProcAddress(havi,"AVIStreamReadFormat"))(stream, 0, 0, &fsize);
  bsize = fsize;
  if (!(buf=malloc2(bsize)))  return(0);
  memset(buf, 0, bsize);
  (*GetProcAddress(havi,"AVIStreamReadFormat"))(stream, 0, buf, &fsize);
  if ((*GetProcAddress(havi,"AVIFileCreateStream"))(cavi, &new, info)) {
    free2(buf);  return(0); }
  if ((*GetProcAddress(havi,"AVIStreamSetFormat"))(new, 0, buf, fsize)) {
    free2(buf);  return(0); }
  bsize = fsize;  stsize = 0;
  (*GetProcAddress(havi,"AVIStreamRead"))(stream,
                   (*GetProcAddress(havi,"AVIStreamStart"))(stream),
                   AVISTREAMREAD_CONVENIENT, 0, 0, &stsize, 0);
  if (fsize==40)
    if (((((long *)buf)[1]*3+3)&0xFFFFFFFC)*((long *)buf)[2]>stsize)
      stsize = ((((long *)buf)[1]*3+3)&0xFFFFFFFC)*((long *)buf)[2];
  if (stsize>bsize) {
    if (!(buf2=realloc2(buf, stsize))) { free2(buf);  return(0); }
    buf = buf2;  bsize = stsize;  memset(buf, 0, bsize); }
  start = (*GetProcAddress(havi,"AVIStreamStart"))(stream);
  end = (*GetProcAddress(havi,"AVIStreamStart"))(stream) +
        (*GetProcAddress(havi,"AVIStreamLength"))(stream);
  progress(hwnd, "Compressing AVI", "Compressing AVI", 0, 1);
  for (pos=start; pos<end; ) {
    sprintf(text, "Compressing frame %d of %d.", pos-start+1, end-start);
    if (!progress(0,0, text, pos-start, end-start))  break;
    (*GetProcAddress(havi,"AVIStreamRead"))(stream, pos, 1,0,0, &ssize, 0);
    (*GetProcAddress(havi,"AVIStreamRead"))(stream, pos, bsize/ssize, buf,
                     bsize, &stsize, &read);
    (*GetProcAddress(havi,"AVIStreamWrite"))(new, pos, read, buf, stsize,
                     AVIIF_KEYFRAME, &read, &stsize);
    pos += read; }
  free2(buf);
  progress(0,0,0,0,0);
  (*GetProcAddress(havi,"AVIStreamRelease"))(new);
}

copy()
/* Copy the current image to the clipboard.                     1/8/98-DWM */
{
  HGLOBAL gmem;
  HANDLE wpic;
  char *dest, *src;
  HWND hwnd;
  DATA *data;
  RECT rect;
  long w, h, ow, oh, old=0, freeimg=1;

  if (!(hwnd=SendMessage(HwndC,WM_MDIGETACTIVE,0,0)))  return(0);
  if (!(data=lock_window(hwnd)))  return(0);
  cursor(1);
  if (OpenClipboard(hwnd)) {
    GetClientRect(hwnd, &rect);
    w = rect.right;  h = rect.bottom;
    if (w<1) w = 1;  if (h<1) h = 1;
    if (((data->dflag>>9)&1) && data->w>0 && data->h>0) {
      w = data->w;  h = data->h; }
    if (w!=data->renderval[2] || h!=data->renderval[3]) {
      old = 1;  ow = data->renderval[2];  oh = data->renderval[3];
      render_move(0, 0, 0, w, h, data->renderval, data->renderdlt,
                  data->renderdlt); }
    wpic = unlock2(update_bmp(data, w, h, &freeimg, 24));
    if (old)
      render_move(0, 0, 0, ow, oh, data->renderval, data->renderdlt,
                  data->renderdlt);
    gmem = GlobalAlloc(GHND|GMEM_DDESHARE, GlobalSize(wpic));
    if (!gmem) {
      error(hwnd, "Low memory: Clipboard Allocation");
      free2(lock2(wpic));  unlock_window(data);
      cursor(0);
      return(0); }
    dest = GlobalLock(gmem);
    src = lock2(wpic);
    memmove(dest, src, GlobalSize(wpic));
    free2(src);
    GlobalUnlock(gmem);
    EmptyClipboard();
    SetClipboardData(CF_DIB, gmem);
    CloseClipboard(); }
  cursor(0);
  unlock_window(data);
}

drop_file(HANDLE drop)
/* Handle dropped-in files just as if they were opened.        10/1/96-DWM */
{
  long n, i, len=0;
  char name[NAMELEN], name2[NAMELEN];

  n = DragQueryFile(drop, -1, 0, 0);
  if (n==0) {
    DragFinish(drop);  return(0); }
  for (i=0; i<n && len<NAMELEN-1; i++) {
    name[len] = 0;
    DragQueryFile(drop, i, name2, NAMELEN-len-2);
    if (name2[0]!='"')  sprintf(name+strlen(name), "\"%s\"", name2);
    else                strcat(name, name2);
    len = strlen(name);
    if (i!=n-1 && len<NAMELEN-1) {
      name[len] = ' ';  len++;  name[len] = 0; } }
  DragFinish(drop);
  strcpy(OpenName, name);
  open_window_mid();
}

long find_file(char *name, long type, long recurse, WIN32_FIND_DATA *find,
               char *path)
/* Locate a file, recursing subdirectories if requested.  Note that the
 *  maximum recursion is 20 subdirectories.  This requires 1k of data space.
 *  To repeat a call to find_file to locate the next file of a group, pass a
 *  null to the name, type, and recurse variables.  To cancel a find_file
 *  call, pass a null to all variables.
 * Enter: char *name: name including wildcards and possibly a path.  If the
 *                    first character of this is a null, then the second
 *                    character specifies how many extensions to search for.
 *                    Starting with the third character is the path, then
 *                    the extensions, which may be a maximum of three
 *                    characters and do not include the period.  The path and
 *                    extensions are separated by nulls.  A maximum of 63
 *                    extensions may be specified.
 *        long type: low byte contains required field values, high byte
 *                   contains field masks.  The fields are: bit 0-read only,
 *                   1-hidden, 2-system, 3-volume ID, 4-subdirectory,
 *                   5-archive.  0x1E00 locates only 'normal' files.
 *        long recurse: 0 for do not recurse subdirectories, 1 for recurse
 *                      all subdirectories, 2 for recurse only subdirectories
 *                      which would have been found by the name
 *                      specification.
 *        WIN32_FIND_DATA *find: structure to return information in.  0 for
 *                               no return.
 *        char *path: string to store path in.  If find is 0, the entire
 *                    file name, including path, is stored here.  This must
 *                    be at least 64 bytes, including the terminating null.
 *                    0 for no return.  There is never a terminating \ on the
 *                    path string.
 * Exit:  long none: 0 for file found, 1 for no more files.      6/7/95-DWM */
{
  #define maxrecursedepth 20
  static WIN32_FIND_DATA f[maxrecursedepth], *last=0;
  static HANDLE h[maxrecursedepth];
  static long depth, new, multi;
  static char root[256], spec[256], xtype, atype, rec;
  char fname[256], *temp;
  long i, found=0;

  if (!name && !path && !find) {
    for (i=0; i<=depth; i++)
      FindClose(h[depth]);
    return(0); }
  if (name) {
    if (!name[0]) {
      multi = name[1];  if (multi<0) multi=0;  if (multi>63) multi=63;
      strcpy(root, name+2);
      strcpy(spec, "*.*");
      name += 2+strlen(name+2)+1;
      for (i=0; i<multi; i++) {
        strcpy(spec+4+i*4, name);
        name += strlen(name)+1; }
      name = 0;
      depth = 0;  new = 1;
      rec = recurse;  atype = type>>8;  xtype = type&0x3F; } }
  if (name) {
    multi = 0;
    getcwd(root, 256-1);
    if (strchr(name, '\\')) {
      strcpy(spec, strrchr(name, '\\')+1);
      strrchr(name, '\\')[0] = 0;
      if (name[1]!=':' && name[0]!='\\') {
        if (strlen(name)+strlen(root)>256-2)
          name[256-1-strlen(root)-1] = 0;
        sprintf(root+strlen(root), "\\%s", name); }
      else
        strcpy(root, name); }
    else
      strcpy(spec, name);
    if (root[strlen(root)-1]=='\\')
      root[strlen(root)-1] = 0;
    depth = 0;  new = 1;
    rec = recurse;  atype = type>>8;  xtype = type&0x3F; }

  while (1) {
    if (new) {                                                /* find file */
      new = 0;
      strcpy(fname, root);
      for (i=0; i<depth && strlen(fname)<256; i++)
        sprintf(fname+strlen(fname), "\\%s", f[i].cFileName);
      if (path)
        strcpy(path, fname);
      sprintf(fname+strlen(fname), "\\%s", spec);
      h[depth] = INVALID_HANDLE_VALUE;
      f[depth].dwFileAttributes = 0x3F;
      if ((h[depth]=FindFirstFile(fname, &f[depth]))!=INVALID_HANDLE_VALUE)
          do {
        if (!((f[depth].dwFileAttributes&atype)^xtype) &&
            strcmp(f[depth].cFileName, ".")) {
          if (multi) {
            temp = strrchr(f[depth].cFileName, '.');
            if (temp)  temp++;
            for (i=0; i<multi && temp; i++)
              if (!stricmp(temp, spec+i*4+4))
                found = 1; }
          else
            found = 1; }
        if (found) {
          if (find)
            memcpy(find, &f[depth], sizeof(WIN32_FIND_DATA));
          if (path && !find)
            sprintf(path+strlen(path), "\\%s", f[depth].cFileName);
          return(0); } }
      while (FindNextFile(h[depth], &f[depth]));
      if (h[depth]!=INVALID_HANDLE_VALUE)
        FindClose(h[depth]); }
    else {
      while (FindNextFile(h[depth], &f[depth])) {
        if (!((f[depth].dwFileAttributes&atype)^xtype) &&
            strcmp(f[depth].cFileName, ".")) {
          if (multi) {
            temp = strrchr(f[depth].cFileName, '.');
            if (temp)  temp++;
            for (i=0; i<multi && temp; i++)
              if (!stricmp(temp, spec+i*4+4))
                found = 1; }
          else
            found = 1; }
        if (found) {
          if (find)
            memcpy(find, &f[depth], sizeof(WIN32_FIND_DATA));
          if (path) {
            sprintf(path, "%s", root);
            for (i=0; i<depth && strlen(fname)<256; i++)
              sprintf(path+strlen(path), "\\%s", f[i].cFileName); }
          if (path && !find)
            sprintf(path+strlen(path), "\\%s", f[depth].cFileName);
          return(0); } }
      if (h[depth]!=INVALID_HANDLE_VALUE)
        FindClose(h[depth]); }
    if (!rec)  return(1);                              /* find directories */
    sprintf(fname, "%s\\", root);
    for (i=0; i<depth && strlen(fname)<256; i++)
      sprintf(fname+strlen(fname), "%s\\", f[i].cFileName);
    if (rec==2)  strcat(fname, spec);
    else         strcat(fname, "*.*");
    if (depth<maxrecursedepth) {
      h[depth] = INVALID_HANDLE_VALUE;
      f[depth].dwFileAttributes = 0x10;
      if ((h[depth]=FindFirstFile(fname, &f[depth]))!=INVALID_HANDLE_VALUE)
          do
        if ((f[depth].dwFileAttributes&0x10) && strcmp(f[depth].cFileName,
            ".") && strcmp(f[depth].cFileName, "..")) {
          depth++;  new = 1;  break; }
      while (FindNextFile(h[depth], &f[depth]));
      if (h[depth]!=INVALID_HANDLE_VALUE && !new)
        FindClose(h[depth]); }
    while (!new) {
      if (!depth) return(1);
      depth--;
      while (FindNextFile(h[depth], &f[depth]))
        if ((f[depth].dwFileAttributes&0x10) && strcmp(f[depth].cFileName,
            ".") && strcmp(f[depth].cFileName, "..")) {
          depth++;  new = 1;  break; }
      if (h[depth]!=INVALID_HANDLE_VALUE && !new)
        FindClose(h[depth]); } }
}

char *find_space(char *text)
/* Locate the first space in a list of files.  Note that this does not count
 *  spaces in long file names.  Any long file name must be surrounded by
 *  quotes.
 * Enter: char *text: text to search.
 * Exit:  char *space: location of first space, or zero for none.
 *                                                             11/9/96-DWM */
{
  long first=1;

  while (text[first-1]=='"')  if (text[first]=='"') first++; else break;
  if (text[0]=='"')
    if (strchr(text+first, '"'))
      text = strchr(text+first, '"');
  return(strchr(text, ' '));
}

VOID CALLBACK hide_splash(HWND hwnd, ulong msg, ulong id, long time)
/* Destroy the splash screen window.
 * Enter: HWND hwnd, ulong msg, id, time: ignored.             11/5/97-DWM */
{
  HDC hdc;

  KillTimer(HwndSplash, 5);
  if (HwndSplash) {
    free2(lock2(SplashBMP));
    if (HpalSplash)  DeleteObject(HpalSplash);
    HpalSplash = 0;
    DestroyWindow(HwndSplash);
    HwndSplash = 0; }
}

join_names(char *text)
/* Combine a directory and a file name which are separated by a space.
 *  Either portion may be surrounded by quotes, in which case the final part
 *  will have quotes.
 * Enter: char *text: text to join.  There may be additional text after the
 *                    file name.                               11/9/96-DWM */
{
  char *loc, *loc2;

  if (loc=find_space(text)) {
    if (text[0]=='"' && loc[1]=='"') {
      loc[-1] = '\\';  memmove(loc, loc+2, strlen(loc+1)); }
    else if (text[0]=='"') {
      if (!(loc2=find_space(loc+1)))  loc2 = strlen(loc+1);
      loc[-1] = '\\';  memmove(loc, loc+1, loc2-loc-1);  loc2[-1] = '"'; }
    else if (loc[1]=='"') {
      memmove(text+1, text, loc-text);
      text[0] = '"';  loc[1] = '\\'; }
    else
      loc[0] = '\\'; }
}

new_window()
/* Open a new window.  The actual work is done by open_window_mid() and
 *  new_window_mid().                                          6/18/97-DWM */
{
  OpenName[0] = 0;
  open_window_mid();
}

open_ov(FILE *fptr, DATA *data)
/* Read in an OV file.  This is an OV header, followed by a complete ORB
 *  file, including header, followed by OV specific data.  This routine does
 *  not close the file.
 * Enter: FILE *fptr: pointer to the open file to read.
 *        DATA *data: pointer to the window's primary data array.
 * Exit:  long error: 0 for okay, 1 for memory error, 2 for read error.
 *                                                              1/9/98-DWM */
{
  long err, i=256, header[2]={0,0}, val, nump, numl;
  uchar buf[257];
  DATA *seq[4];
  OATOM *orb;
  LIGHT *ls;

  if (err=open_orb(fptr, data, 0))  return(err);
  while (fread(buf, 1, 256, fptr)==256) {
    for (i=0; i<256; i++)  if (buf[i]==0xFF)  break;
    if (i!=256) break; }
  if (i>=256 || buf[i]!=0xFF)  return(0);
  fseek(fptr, -(256-i)+4, SEEK_CUR);
  fread(header, 4, 2, fptr);
  if (header[0]!=sizeof(DATA))  return(0);
  for (i=0; i<4; i++)  seq[i] = data->seq[i];
  orb = data->mol.orb;  nump = data->mol.nump;
  ls = data->mol.ls;    numl = data->mol.numl;
  fread(data, sizeof(DATA), 1, fptr);
  val = sizeof(DATA) +
      data->mol.nump*sizeof(OATOM) + data->mol.numl*sizeof(LIGHT) +
      data->asym.n*3*sizeof(float) + data->asym.e*3*sizeof(long) +
      (data->asym.rfac!=0)*(data->asym.density+1)*sizeof(real) +
      (data->asym.tfac!=0)*(data->asym.density+1)*sizeof(real) +
      (data->asym.pfac!=0)*(data->asym.density*2+1)*sizeof(real) +
      data->poly.n*3*sizeof(float) + data->poly.e*3*sizeof(long) +
      data->poly.n*sizeof(long) + data->poly.e +
      (data->poly.rfac!=0)*(data->poly.density+1)*sizeof(real) +
      (data->poly.tfac!=0)*(data->poly.density+1)*sizeof(real) +
      (data->poly.pfac!=0)*(data->poly.density*2+1)*sizeof(real) +
      data->points.n*3*sizeof(float) + data->points.n +
      (data->render.buf!=0)*(-data->render.scan*data->render.h+40) +
      (data->render.zbuf!=0)*data->render.w*data->render.h*sizeof(ushort) +
      (data->render.phase!=0)*data->render.h*data->render.w*sizeof(long) +
      (data->stereo.image!=0)*(data->stereo.w*data->stereo.h*
                               (3-2*data->stereo.pal)+768*data->stereo.pal);
  if (val!=header[1] || data->mol.nump!=nump || data->mol.numl!=numl) {
    memset(data, 0, sizeof(DATA));
    for (i=0; i<4; i++)  if (seq[i]) {
      free2(seq[i]->mol.orb);
      free2(seq[i]->mol.ls);
      free2(seq[i]); }
    data->mol.orb = orb;  data->mol.ls = ls;
    rewind(fptr);
    new_window_mid(data);
    return(open_orb(fptr, data, 0)); }
  for (i=0; i<4; i++)
    data->seq[i] = seq[i];
  data->mol.orb = orb;  data->mol.ls = ls;
  fread(data->mol.orb, sizeof(OATOM), data->mol.nump, fptr);
  if (data->mol.numl)     fread(data->mol.ls, sizeof(LIGHT), data->mol.numl, fptr);
  data->asym.xyz = data->asym.elem = data->asym.scrxyz = data->asym.norm = 0;
  data->asym.enorm = 0;
  if (data->asym.n)       data->asym.xyz = malloc2(sizeof(float)*3*data->asym.n);
  if (data->asym.xyz)     fread(data->asym.xyz, sizeof(float)*3, data->asym.n, fptr);
  if (data->asym.e)       data->asym.elem = malloc2(sizeof(long)*3*data->asym.e);
  if (data->asym.elem)    fread(data->asym.elem, sizeof(long)*3, data->asym.e, fptr);
  if (data->asym.rfac)    data->asym.rfac = malloc2(sizeof(real)*(data->asym.density+1));
  if (data->asym.rfac)    fread(data->asym.rfac, sizeof(real), data->asym.density+1, fptr);
  if (data->asym.tfac)    data->asym.tfac = malloc2(sizeof(real)*(data->asym.density+1));
  if (data->asym.tfac)    fread(data->asym.tfac, sizeof(real), data->asym.density+1, fptr);
  if (data->asym.pfac)    data->asym.pfac = malloc2(sizeof(real)*(data->asym.density*2+1));
  if (data->asym.pfac)    fread(data->asym.pfac, sizeof(real), data->asym.density*2+1, fptr);
  if (!data->asym.xyz || !data->asym.elem) {
    free2(data->asym.xyz);  free2(data->asym.elem);
    data->asym.n = data->asym.e = 0; }
  data->poly.xyz = data->poly.elem = data->poly.scrxyz = data->poly.norm = 0;
  data->poly.enorm = data->poly.phase = data->poly.ptphase = 0;
  if (data->poly.n)       data->poly.xyz = malloc2(sizeof(float)*3*data->poly.n);
  if (data->poly.xyz)     fread(data->poly.xyz, sizeof(float)*3, data->poly.n, fptr);
  if (data->poly.e)       data->poly.elem = malloc2(sizeof(long)*3*data->poly.e);
  if (data->poly.elem)    fread(data->poly.elem, sizeof(long)*3, data->poly.e, fptr);
  if (data->poly.n)       data->poly.ptphase = malloc2(sizeof(long)*data->poly.n);
  if (data->poly.ptphase) fread(data->poly.ptphase, sizeof(long), data->poly.n, fptr);
  if (data->poly.e)       data->poly.phase = malloc2(data->poly.e);
  if (data->poly.phase)   fread(data->poly.phase, 1, data->poly.e, fptr);
  if (data->poly.rfac)    data->poly.rfac = malloc2(sizeof(real)*(data->poly.density+1));
  if (data->poly.rfac)    fread(data->poly.rfac, sizeof(real), data->poly.density+1, fptr);
  if (data->poly.tfac)    data->poly.tfac = malloc2(sizeof(real)*(data->poly.density+1));
  if (data->poly.tfac)    fread(data->poly.tfac, sizeof(real), data->poly.density+1, fptr);
  if (data->poly.pfac)    data->poly.pfac = malloc2(sizeof(real)*(data->poly.density*2+1));
  if (data->poly.pfac)    fread(data->poly.pfac, sizeof(real), data->poly.density*2+1, fptr);
  if (!data->poly.xyz || !data->poly.elem || !data->poly.ptphase ||
      !data->poly.phase) {
    free2(data->poly.xyz);      free2(data->poly.elem);
    free2(data->poly.ptphase);  free2(data->poly.phase);
    data->poly.n = data->poly.e = 0; }
  data->points.scrxyz = data->points.xyz = data->points.phase = 0;
  if (data->points.n)     data->points.xyz = malloc2(sizeof(float)*3*data->points.n);
  if (data->points.xyz)   fread(data->points.xyz, sizeof(float)*3, data->points.n, fptr);
  if (data->points.n)     data->points.phase = malloc2(data->points.n);
  if (data->points.phase) fread(data->points.phase, 1, data->points.n, fptr);
  if (!data->points.xyz || !data->points.phase) {
    free2(data->points.xyz);  free2(data->points.phase);
    data->points.n = 0; }
  if (data->render.buf)   data->render.buf = malloc2(-data->render.scan*data->render.h+40);
  if (data->render.buf)   fread(data->render.buf, 1, -data->render.scan*data->render.h+40, fptr);
  if (data->render.zbuf)  data->render.zbuf = malloc2(sizeof(ushort)*data->render.w*data->render.h);
  if (data->render.zbuf)  fread(data->render.zbuf, sizeof(ushort), data->render.w*data->render.h, fptr);
  if (data->render.phase) data->render.phase = malloc2(sizeof(long)*data->render.w*data->render.h);
  if (data->render.phase) fread(data->render.phase, sizeof(long), data->render.w*data->render.h, fptr);
  if (data->stereo.image) data->stereo.image = malloc2(data->stereo.w*data->stereo.h*(3-2*data->stereo.pal)+768*data->stereo.pal);
  if (data->stereo.image) fread(data->stereo.image, 1, data->stereo.w*data->stereo.h*(3-2*data->stereo.pal)+768*data->stereo.pal, fptr);
  return(0);
}

open_window()
/* Open a new file.                                             3/23/97-DWM */
{
  long len, i, j;
  char name[NAMELEN*2];

  strcpy(opfilter, OpenFilter);
  len = strlen(opfilter);
  for (i=0; i<len; i++)  if (opfilter[i]=='|')  opfilter[i] = 0;
  openfile.lStructSize = sizeof(OPENFILENAME);
  openfile.hwndOwner = Hwnd;  openfile.hInstance = hinst;
  openfile.lpstrFilter = opfilter;
  openfile.nFilterIndex = FilterType[0];
  openfile.lpstrFile = OpenName;
  openfile.nMaxFile = NAMELEN-1;
  openfile.lpstrFileTitle = openfile.lpstrCustomFilter = 0;
  openfile.lpstrInitialDir = OpenDir;
  openfile.Flags = OFN_FILEMUSTEXIST|OFN_HIDEREADONLY|OFN_PATHMUSTEXIST|
                   OFN_ALLOWMULTISELECT|OFN_EXPLORER;
  memset(OpenName, 0xFF, NAMELEN);
  OpenName[0] = 0;
  if (!GetOpenFileName(&openfile))
    return(0);
  FilterType[0] = openfile.nFilterIndex;
  for (i=len=j=0; i<NAMELEN; i++)
    if (OpenName[i]) { if (len>j)  j = len;  len = 0; }
    else               len++;
  if (j==2) {
    name[0] = i = 0;
    do {
      sprintf(name+strlen(name), "\"%s\" ", OpenName+i);
      i += strlen(OpenName+i)+1; }
    while (OpenName[i]);
    name[strlen(name)-1] = 0;
    strncpy(OpenName, name, NAMELEN-1);
    OpenName[NAMELEN-1] = 0; }
  join_names(OpenName);
  open_window_mid();
}

open_window_mid()
/* Open the file(s) specified in the OpenName array.  If the OpenName array
 *  is the null string, open a new window.                     2/12/97-DWM */
{
  FILE *fptr;
  char name[NAMELEN], text[NAMELEN], buf[51];
  long type;
  long i, err=0, winmode=SW_SHOW;
  DATA *data;
  HANDLE hnd;
  MDICREATESTRUCT m;
  HWND hwnd;
  static untitled=1;

  if (OpenName[0]) {
    if (strchr(OpenName, '\\'))
      strpath(OpenDir, OpenName);
    strcpy(name, OpenName);
    remove_quotes(name);
    if (!(fptr=fopen(name, "rb"))) {
      sprintf(text, "Can't open file %s.", name);
      error(Hwnd, text);
      return(0); }
    memset(buf, 0, 51);
    fread(buf, 1, 50, fptr);
    rewind(fptr);
    type = -1;
    if (!memcmp(buf, OrbCheck, strlen(OrbCheck)))    type = FileOrb;
    if (!memcmp(buf, OrbCheck2, strlen(OrbCheck2)))  type = FileOV;
    /** Check for additional types here **/
    if (type==-1) {
      error(Hwnd, "Can't identify file format.");
      fclose(fptr);  return(0); } }
  else {
    sprintf(name, "Untitled %d", untitled);
    untitled++; }
  if (hwnd=SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0)) {
    if (IsZoomed(hwnd))        Pref.flags |= 2;
    else if (!IsIconic(hwnd))  Pref.flags &= 0xFFFFFFFD; }
  if (!(data=malloc2(sizeof(DATA)))) {
    error(Hwnd, "Low memory.\n\nCan't open window.");
    if (OpenName[0])  fclose(fptr);  return(0); }
  memset(data, 0, sizeof(DATA));
  for (i=0; i<strlen(name); i++)
    if (name[i]==toupper(name[i])) break;
  if (i==strlen(name))
    for (i=0; i<strlen(name); i++)
      name[i] = toupper(name[i]);
  strcpy(data->name, name);
  strcpy(text, name);
  if (strchr(name, '\\'))  strcpy(text, strrchr(name, '\\')+1);
  m.szClass = WinName2;  m.szTitle = text;
  m.hOwner = hinst;
  m.x = m.y = m.cx = m.cy = CW_USEDEFAULT;
  m.style = (WS_MAXIMIZE*((Pref.flags&2)!=0));
  m.lParam = 0;
  hwnd = SendMessage(HwndC, WM_MDICREATE, 0, (long)(&m));
  hnd = unlock2(data);
  SetWindowLong(hwnd, GWL_USERDATA, (long)hnd);
  ShowWindow(hwnd, SW_SHOW);
  status("Loading file.");
  data = lock_window(hwnd);
  new_window_mid(data);
  if (OpenName[0])
    switch (type) {
      case FileOrb: err = open_orb(fptr, data, 0); break;
      case FileOV: err = open_ov(fptr, data); break;
      default: err = 3; }
  else {
    if (Pref.flags&1)
      if (fptr=fopen(DefFile, "rb")) {
        open_orb(fptr, data, 0);
        sprintf(data->name, "Untitled %d", untitled-1);
        fclose(fptr); }
    if (!data->mol.nump)  data->mol.nump = 1; }
  if (data->place.length)
    SetWindowPlacement(hwnd, &data->place);
  data->filetype = type;
  if (data->mol.nump<1)  err = 2;
  update_process(data);
  unlock_window(data);
  if (OpenName[0]) {
    strcpy(data->name, name);
    fclose(fptr); }
  if (err) {
    switch (err) {
      case 1: error(Hwnd, "Insufficient memory."); break;
      case 2: error(Hwnd, "Couldn't read file."); break;
      case 3: error(Hwnd, "Function not implemented."); }
    SendMessage(HwndC, WM_MDIDESTROY, (long)(hwnd), 0);
    return(0); }
  if (OpenName[0]) {
    if (!strchr(name, '\\')) {
      getcwd(text, NAMELEN-1);
      if (strlen(text))
        if (text[strlen(text)-1]!='\\')
          strcat(text, "\\");
      strcat(text, name);
      strcpy(name, text); }
    for (i=0; i<3; i++)
      if (!stricmp(name, lastview+i*256))
        break;
    memmove(lastview+256, lastview, 256*i);
    strcpy(lastview, name); }
  status(0);
  SetTimer(hwnd, 1, 0, timer);
  if (OpenName[0]) {
    if (find_space(OpenName)) {
      strcpy(name, find_space(OpenName)+1);
      strcpy(OpenName, name);
      open_window_mid(); } }
}

BOOL CALLBACK preferences_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the Preferences dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.         1/8/98-DWM */
{
  HWND hwnd;
  DATA *data;
  char text[80], text2[80];
  double x2;
  long ru, i, j, init=0;

  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case HelpHelp: WinHelp(Hwnd,HelpFile,HELP_CONTEXT,HelpPreferences);
        break;
      case IDOK: i = Pref.toolbar*2+Pref.status;
        Pref.error = SendDlgItemMessage(hdlg, PrefError, BM_GETCHECK, 0, 0);
        Pref.warn  = SendDlgItemMessage(hdlg, PrefWarn,  BM_GETCHECK, 0, 0);
        Pref.toolbar=SendDlgItemMessage(hdlg, PrefTool,  BM_GETCHECK, 0, 0);
        Pref.tooltip=SendDlgItemMessage(hdlg, PrefTips,  BM_GETCHECK, 0, 0);
        Pref.status= SendDlgItemMessage(hdlg, PrefStatus,BM_GETCHECK, 0, 0);
        if (Pref.toolbar*2+Pref.status!=i)  setup();
        Pref.flags &= 0xFFFFFBE3;
        Pref.flags|=(1<<10)*!SendDlgItemMessage(hdlg,PrefSplash,BM_GETCHECK,0,0);
        Pref.flags|=(1<<2)*!SendDlgItemMessage(hdlg,PrefDiffuse,BM_GETCHECK,0,0);
        Pref.flags|=(1<<3)*!SendDlgItemMessage(hdlg,PrefAmbient,BM_GETCHECK,0,0);
        Pref.flags|=(1<<4)*!SendDlgItemMessage(hdlg,PrefEmissive,BM_GETCHECK,0,0);
        i = Pref.pointsize;
        Pref.pointsize = SendDlgItemMessage(hdlg,PrefPoint,CB_GETCURSEL,0,0);
        GetDlgItemText(hdlg, PrefRad, text, 79);
        x2 = Pref.radstep;
        sscanf(text, "%lf", &x2);  x2 *= RadVal[Pref.radunit];
        Pref.radstep = x2;
        GetDlgItemText(hdlg, PrefPan, text, 79);
        x2 = Pref.panstep*100;  sscanf(text, "%lf", &x2);
        Pref.panstep = x2*0.01;
        GetDlgItemText(hdlg, PrefZoom, text, 79);
        x2 = Pref.zoomstep*100;  sscanf(text, "%lf", &x2);
        Pref.zoomstep = x2*0.01;
        if (i!=Pref.pointsize)
          for (i=0; i<NumWindows; i++) {
            data = lock_window(HwndList[i]);
            if (data->points.process>=5)  data->points.process = 6;
            if (data->asym.process>=5)    data->asym.process = 6;
            if (data->poly.process>=6)    data->poly.process = 7;
            unlock_window(data); }
      case IDCANCEL: EndDialog(hdlg, 1);  return(1);
      case PrefRadUnit: if ((wp>>16)!=CBN_SELCHANGE) break;
        ru = Pref.radunit;
        GetDlgItemText(hdlg, PrefRad, text, 79);
        x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= RadVal[ru];
        ru = SendDlgItemMessage(hdlg, PrefRadUnit, CB_GETCURSEL, 0, 0);
        Pref.radunit = ru;  x2 /= RadVal[ru];
        if (fabs(x2)<1e-6)  x2 = 0;  sprintf(text, "%g", x2);
        sprintf(text2, "%g", x2+1e-6-2e-6*(x2<0));
        if (strlen(text2)<strlen(text)) strcpy(text, text2);
        SetDlgItemText(hdlg, PrefRad, text); break;
      case PrefReset:
        if (!warning(hdlg, "Reset all preferences and colors?"))  break;
        read_ini_defaults();  init = 1; break;
      default: return(0); } break;
    case WM_INITDIALOG: init = 1;
      for (i=0; i<3; i++) {
        sprintf(text, "%d", i*2+1);
        SendDlgItemMessage(hdlg,PrefPoint,CB_ADDSTRING, 0,(long)text); }
      for (i=0; i<2; i++)
        SendDlgItemMessage(hdlg,PrefRadUnit,CB_ADDSTRING, 0,(long)RadText[i]);
      SetFocus(hdlg); }
  if (init) {
    SendDlgItemMessage(hdlg, PrefError, BM_SETCHECK, Pref.error, 0);
    SendDlgItemMessage(hdlg, PrefWarn,  BM_SETCHECK, Pref.warn, 0);
    SendDlgItemMessage(hdlg, PrefTool,  BM_SETCHECK, Pref.toolbar, 0);
    SendDlgItemMessage(hdlg, PrefTips,  BM_SETCHECK, Pref.tooltip, 0);
    SendDlgItemMessage(hdlg, PrefStatus,BM_SETCHECK, Pref.status, 0);
    SendDlgItemMessage(hdlg,PrefSplash,BM_SETCHECK,!((Pref.flags>>10)&1),0);
    SendDlgItemMessage(hdlg,PrefDiffuse,BM_SETCHECK,!((Pref.flags>>2)&1),0);
    SendDlgItemMessage(hdlg,PrefAmbient,BM_SETCHECK,!((Pref.flags>>3)&1),0);
    SendDlgItemMessage(hdlg,PrefEmissive,BM_SETCHECK,!((Pref.flags>>4)&1),0);
    SendDlgItemMessage(hdlg,PrefPoint, CB_SETCURSEL,Pref.pointsize,0);
    SendDlgItemMessage(hdlg,PrefRadUnit, CB_SETCURSEL,Pref.radunit,0);
    ru = Pref.radunit;  x2 = Pref.radstep/RadVal[ru];
    if (fabs(x2)<1e-6)  x2 = 0;  sprintf(text, "%g", x2);
    sprintf(text2, "%g", x2+1e-6-2e-6*(x2<0));
    if (strlen(text2)<strlen(text)) strcpy(text, text2);
    SetDlgItemText(hdlg, PrefRad, text);
    sprintf(text, "%g", Pref.panstep*100);
    SetDlgItemText(hdlg, PrefPan, text);
    sprintf(text, "%g", Pref.zoomstep*100);
    SetDlgItemText(hdlg, PrefZoom, text); }
  return(0);
}

read_ini()
/* Locate and read in values from the INI file.                5/26/96-DWM */
{
  char *cmd, text[256], text2[32], *t2;
  long i, j, k, def, def2, numt=0, buf[5];

  WinPlace.length = 0;
  read_ini_defaults();
  OpenDir[0] = OpenDir[DIRLEN] = OpenDir[DIRLEN*2] = OpenDir[DIRLEN*3] = 0;
  OpenDir[DIRLEN*4] = OpenDir[DIRLEN*5] = OpenDir[DIRLEN*6] = 0;
  cmd = GetCommandLine();
  if (strlen(cmd)<256) {
    strcpy(inifile, cmd);
    remove_quotes(inifile);
    if (strchr(inifile, '\\'))  strrchr(inifile, '\\')[0] = 0;
    sprintf(DefFile, "%s\\DEFAULT.ORB", inifile);
    strcat(inifile, "\\OV.INI"); }
  for (i=0; key[i*3]; i++) {
    GetPrivateProfileString(inihead[0],key[i*3],key[i*3+1],text,256,inifile);
    def = 1-(!strcmpi(text, key[i*3+1]));
    def2 = def;
    if (key[i*3+2]) if (def && strcmpi(text, key[i*3+2]))  def2++;
    switch (i) {
      case 0: for (j=0; j<sizeof(WINDOWPLACEMENT) && j*2<strlen(text); j++) {
          sscanf(text+j*2, "%2X", &k);
          ((char *)(&WinPlace.length))[j] = k; } break;
      case 1: Pref.error = def^1; break;
      case 2: Pref.warn = 1-def; break;
      case 3: Pref.status = 1-def; break;
      case 4: Pref.toolbar = 1-def; break;
      case 5: Pref.tooltip = 1-def; break;
      case 6: sscanf(text, "%d", &numt); break;
      case 7: if (strlen(text)) strcpy(OpenDir, text); break;
      case 8: sscanf(text, "%d%*c%d%*c%d%*c%d%*c%d%*c%d%*c%d", FilterType,
              FilterType+1, FilterType+2, FilterType+3, FilterType+4,
              FilterType+5, FilterType+6); break;
      case 9: sscanf(text, "%d", &Pref.flags); break;
      case 10: sscanf(text, "%g", &Pref.perspec); break;
      case 11: if (!strlen(text)) break;
        t2 = text;
        for (j=0; j<16 && j<NUMPREFCOLORS && t2; j++) {
          sscanf(t2, "%X", Pref.colors+j);
          if (strchr(t2, ','))  t2 = strchr(t2, ',')+1;
          else                  t2 = 0; } break;
      case 12: if (strlen(text)) strcpy(lastview,     text); break;
      case 13: if (strlen(text)) strcpy(lastview+256, text); break;
      case 14: if (strlen(text)) strcpy(lastview+512, text); break;
      case 15: if (strlen(text)) strcpy(lastview+768, text); break;
      case 16: if (strlen(text)) strcpy(OpenDir+DIRLEN, text); break;
      case 17: if (strlen(text)) strcpy(OpenDir+DIRLEN*2, text); break;
      case 18: sscanf(text, "%g", &Pref.radstep); break;
      case 19: sscanf(text, "%g", &Pref.panstep); break;
      case 20: sscanf(text, "%g", &Pref.zoomstep); break;
      case 21: if (strlen(text)) {
        sscanf(text, "%d", &Pref.pointsize);
        Pref.pointsize /= 2; } break;
      case 22: for (j=0; j<2; j++)
          if (!stricmp(text, RadText[j])) Pref.radunit = j; break;
      case 23: for (j=0; j<4; j++)
          if (!stricmp(text, DistText[j])) Pref.munit = j; break;
      case 24: for (j=0; j<4; j++)
          if (!stricmp(text, MassText[j])) Pref.kgunit = j; break;
      case 25: if (strlen(text)) strcpy(OpenDir+DIRLEN*3, text); break;
      case 26: if (strlen(text)) strcpy(OpenDir+DIRLEN*4, text); break;
      case 27: if (strlen(text)) strcpy(OpenDir+DIRLEN*5, text); break;
      case 28: if (strlen(text)) strcpy(OpenDir+DIRLEN*6, text); break;
      default: ; } }
  read_ini_toolbar(numt-(!numt));
}

read_ini_defaults()
/* Setup the global Pref structure to the default values.       1/9/98-DWM */
{
  long i;

  memset(&Pref, 0, sizeof(PREF));
  Pref.error = Pref.warn = Pref.status = Pref.toolbar = Pref.tooltip = 1;
  Pref.perspec = 25;
  Pref.pointsize = 0;
  Pref.radstep = 30*deg;  Pref.panstep = 0.15;  Pref.zoomstep = 0.5;
  for (i=0; i<NUMPREFCOLORS; i++)
    Pref.colors[i] = -1;
}

read_ini_toolbar(long num)
/* Read in a customized toolbar.
 * Enter: long num: number of toolbar entries to setup function (no reading
 *                  is done), or zero to actually read ini file.
 *                                                            12/10/97-DWM */
{
  char text[256], text2[32], t2[256];
  long i, j, k, l;
  static long numt=0;

  if (num!=0) { numt = num;  return(0); }
  if (numt<=0)  return(0);
  if (ToolCustom)  free2(ToolCustom);
  if (!(ToolCustom=malloc2((numt+1)*sizeof(long))))  return(0);
  memset(ToolCustom, 0, (numt+1)*sizeof(long));
  for (i=l=0; l<numt; l++) {
    sprintf(text2, "ToolBar%d", l);
    GetPrivateProfileString(inihead[1], text2, "", text, 256, inifile);
    if (!strcmp(text, "---")) {
      ToolCustom[i] = -1;  i++;  continue; }
    for (k=0; ToolList[k].idCommand; k++) {
      GetMenuString(Hmenu, ToolList[k].idCommand, t2, 255, MF_BYCOMMAND);
      if (!strncmp(t2, "Redo", 4) || !strncmp(t2, "Undo", 4))
        strcpy(t2, "Undo");
      if (strchr(t2, '\t'))  strchr(t2, '\t')[0] = 0;
      for (j=0; j<strlen(t2); j++)
        if (t2[j]=='&' || t2[j]=='.' || t2[j]==' ') {
          memmove(t2+j, t2+j+1, strlen(t2+j));  j--; }
      if (!strcmp(t2, text)) break; }
    if (!ToolList[k].idCommand)  continue;
    ToolCustom[i] = ToolList[k].idCommand;
    i++; }
  if (!i)  { if (ToolCustom)  free2(ToolCustom);  ToolCustom = 0; }
}

recheck(HWND hwnd, long refresh)
/* Refresh all menu check marks to ensure that they are correct.
 * Enter: HWND hwnd: handle of window.  0 for no window.
 *        long refresh: 1 to redraw window afterwards.         5/26/96-DWM */
{
  DATA *data;
  long i, en;
  HMENU sub;
  char text[300], *text2;
  long first=0;
  long enlist[]={MenuClose,MenuSave,MenuSaveAs,MenuCopy,MenuRenderOpt,
       MenuPointOpt,MenuPolyOpt,MenuRayOpt,MenuAsymOpt,MenuOrb,MenuLight,
       MenuStereo,MenuCamera,MenuLeft,MenuRight,MenuUp,MenuDown,MenuRotXM,
       MenuRotXP,MenuRotYM,MenuRotYP,MenuRotZM,MenuRotZP,MenuZoomIn,
       MenuZoomOut,MenuFocalIn,MenuFocalOut,MenuResetPos,MenuReframe,
       MenuCutAway,MenuDefault,MenuSequence,MenuStop,-1};

  en = (!AVILoad || havi);
  EnableMenuItem(Hmenu, MenuCompAVI, MF_BYCOMMAND|(MF_GRAYED*(!en)));
  SendMessage(HwndTool, TB_ENABLEBUTTON, MenuCompAVI, en);
  data = lock_window(hwnd);
  en = (NumWindows!=0);
  for (i=0; enlist[i]>=0; i++) {
    EnableMenuItem(Hmenu, enlist[i], MF_BYCOMMAND|(MF_GRAYED*(!en)));
    SendMessage(HwndTool, TB_ENABLEBUTTON, enlist[i], en); }
  if (data) {
    en = (en && ((data->dflag>>10)&1));
    EnableMenuItem(Hmenu, MenuStop, MF_BYCOMMAND|(MF_GRAYED*(!en)));
    SendMessage(HwndTool, TB_ENABLEBUTTON, MenuStop, en);
    /** Additional check marks go here **/
    unlock_window(data); }
  sub = GetSubMenu(Hmenu, MenuNew);
  if (GetMenuItemID(sub, 0)==1)  first = 1;
  sub = GetSubMenu(Hmenu, MenuNew-1+first);
  if (refresh)
    InvalidateRect(hwnd, 0, 0);
  for (i=0; i<4; i++)
    RemoveMenu(sub, MenuLastView+i, MF_BYCOMMAND);
  for (i=0; i<4; i++)
    if (lastview[i*256]) {
      text2 = lastview+i*256;
      if (strrchr(text2, '\\'))  text2 = strrchr(text2, '\\')+1;
      sprintf(text, "&%d: %s", i+1, text2);
      AppendMenu(sub, MF_STRING, MenuLastView+i, text); }
}

remove_quotes(char *text)
/* If a text string is surrounded by quotes, they are removed.  The string
 *  is first limited to the text before the first space outside of quotes.
 * Enter: char *text: text string to modify.                   11/9/96-DWM */
{
  long i;

  for (i=0; i<strlen(text)-1; i++)
    if (text[i]=='"' && text[i+1]=='"') {
      memmove(text+i, text+1+i, strlen(text)-i);
      i--; }
  if (find_space(text))   find_space(text)[0] = 0;
  if (text[0]=='"' && text[strlen(text)-1]=='"' && strlen(text)>=2) {
    memmove(text, text+1, strlen(text));
    text[strlen(text)-1] = 0; }
}

save_avi(DATA *data)
/* Save the rendered picture by appending it to an AVI file.  The filename is
 *  in OpenName.
 * Enter: DATA *data: pointer to window data.
 * Exit:  long error: 0 for okay, 1 for insufficient memory, 2 for write
 *                    error.                                    1/4/98-DWM */
{
  uchar *bmp, head[4096];
  long freeimg=1, w, h, old=0, ow, oh, sc, numfr, i;
  RECT rect;
  FILE *fptr;

  GetClientRect(window_handle(data), &rect);
  w = rect.right;  h = rect.bottom;
  if (w<1) w = 1;  if (h<1) h = 1;
  if (((data->dflag>>9)&1) && data->w>0 && data->h>0) {
    w = data->w;  h = data->h; }
  if (w!=data->renderval[2] || h!=data->renderval[3]) {
    old = 1;  ow = data->renderval[2];  oh = data->renderval[3];
    render_move(0, 0, 0, w, h, data->renderval, data->renderdlt,
                data->renderdlt); }
  if (!(bmp=update_bmp(data, w, h, &freeimg, 24)))  return(1);
  if (old)
    render_move(0, 0, 0, ow, oh, data->renderval, data->renderdlt,
                data->renderdlt);
  sc = (w*3+3)&0xFFFFFFFC;
  if (!(fptr=fopen(OpenName, "r+b"))) {
    if (!(fptr=fopen(OpenName, "w+b")))  return(2);
    if (data->seqfps<=0)  data->seqfps = 30;
    memcpy(head, AVIHeader, 4096);
    ((long *)(head+0x4))[0] = 4096-8;
    ((long *)(head+0x20))[0] = 1e6/data->seqfps;
    ((long *)(head+0x24))[0] = sc*h*data->seqfps;
    ((long *)(head+0x30))[0] = 0;
    ((long *)(head+0x3C))[0] = ((long *)(head+0x90))[0] = sc*h;
    ((long *)(head+0x40))[0] = ((short *)(head+0xA0))[0] = w;
    ((long *)(head+0x44))[0] = ((short *)(head+0xA2))[0] = h;
    ((long *)(head+0x84))[0] = data->seqfps;
    memcpy(head+0xAC, bmp, 40); }
  else {
    fread(head, 1, 4096, fptr);
    if (((long *)(head+0x40))[0]!=w || ((long *)(head+0x44))[0]!=h ||
        ((long *)(head+0x24))[0]!= sc*h*data->seqfps) {
      fclose(fptr);  return(2); } }
  numfr = ((long *)(head+0x30))[0];
  ((long *)(head+0x4))[0] = ((sc*h+8)+0x10)*(numfr+1)+4096;
  ((long *)(head+0x30))[0] = ((long *)(head+0x8C))[0] = numfr+1;
  ((long *)(head+0xFF8))[0] = (sc*h+8)*(numfr+1)+4;
  rewind(fptr);
  fwrite(head, 1, 4096, fptr);
  fseek(fptr, 4096+(sc*h+8)*numfr, SEEK_SET);
  fwrite("00db", 1, 4, fptr);
  i = sc*h;
  fwrite(&i, 1, 4, fptr);
  fwrite(bmp+40, i, 1, fptr);
  fwrite("idx1", 1, 4, fptr);
  i = (numfr+1)*0x10;
  fwrite(&i, 1, 4, fptr);
  strcpy(head, "00db");
  ((long *)(head+0x4))[0] = 0x10;
  ((long *)(head+0xC))[0] = sc*h;
  for (i=0; i<numfr+1; i++) {
    ((long *)(head+0x8))[0] = 4+(sc*h+8)*i;
    fwrite(head, 1, 16, fptr); }
  fclose(fptr);
  free2(bmp);
  return(0);
}

BOOL CALLBACK save_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Adjust the extensions in the save dialog.
 * Enter: HWND hwnd: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.        6/11/96-DWM */
{
  char *ext, name[NAMELEN];
  long pos, i;

  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case cmb1:
        if ((wp>>16)!=CBN_SELCHANGE) break;
        GetDlgItemText(hdlg, edt1, name, NAMELEN-1);
        if (strchr(name, '*') || strchr(name, '?')) break;
        ext = name;
        if (strchr(ext, '\\'))  ext = strrchr(ext, '\\')+1;
        if (!strchr(ext, '.'))  ext += strlen(ext);
        else                    ext = strchr(ext, '.');
        ext[0] = 0;
        ext = opfilter;
        pos = SendDlgItemMessage(hdlg, cmb1, CB_GETCURSEL, 0, 0);
        for (i=0; i<pos*2+1; i++)
          ext += strlen(ext)+1;
        sprintf(name+strlen(name), ext+1);
        SetDlgItemText(hdlg, edt1, name); break;
      case IDOK:
        GetDlgItemText(hdlg, edt1, name, NAMELEN-1);
        ext = name;
        if (strchr(ext, '\\'))  ext = strrchr(ext, '\\')+1;
        if (!strchr(ext, '.')) {
          ext = opfilter;
          pos = SendDlgItemMessage(hdlg, cmb1, CB_GETCURSEL, 0, 0);
          for (i=0; i<pos*2+1; i++)
            ext += strlen(ext)+1;
          sprintf(name+strlen(name), ext+1);
          SetDlgItemText(hdlg, edt1, name);
          return(1); } break; } break;
    case WM_INITDIALOG: SetFocus(hdlg); break;
    case WM_NOTIFY: switch (((LPOFNOTIFY)lp)->hdr.code) {
      case CDN_INITDONE:
        ext = opfilter;
        pos = SendDlgItemMessage(((LPOFNOTIFY)lp)->hdr.hwndFrom, cmb1,
                                 CB_GETCURSEL, 0, 0);
        for (i=0; i<pos*2+1; i++)
          ext += strlen(ext)+1;
        SendMessage(((LPOFNOTIFY)lp)->hdr.hwndFrom, CDM_SETDEFEXT, 0,
                    (long)(ext+2)); break;
      case CDN_TYPECHANGE:
        if (!SendMessage(((LPOFNOTIFY)lp)->hdr.hwndFrom, CDM_GETSPEC,
            NAMELEN-1, (long)name)) break;
        if (strchr(name, '*') || strchr(name, '?')) break;
        ext = name;
        if (strchr(ext, '\\'))  ext = strrchr(ext, '\\')+1;
        if (!strchr(ext, '.'))  ext += strlen(ext);
        else                    ext = strchr(ext, '.');
        ext[0] = 0;
        ext = opfilter;
        pos = SendDlgItemMessage(((LPOFNOTIFY)lp)->hdr.hwndFrom, cmb1,
                                 CB_GETCURSEL, 0, 0);
        for (i=0; i<pos*2+1; i++)
          ext += strlen(ext)+1;
        sprintf(name+strlen(name), ext+1);
        SendMessage(((LPOFNOTIFY)lp)->hdr.hwndFrom, CDM_SETDEFEXT, 0,
                    (long)(ext+2));
        SendMessage(((LPOFNOTIFY)lp)->hdr.hwndFrom, CDM_SETCONTROLTEXT, edt1,
                    (long)name); break; } }
  return(0);
}

save_graphic(DATA *data, long type)
/* Save the rendered picture as a graphics file.  The filename is in
 *  OpenName.
 * Enter: DATA *data: pointer to window data.
 *        long type: Either GV_PPM, GV_TIF, or GV_BMP.
 * Exit:  long error: 0 for okay, 1 for insufficient memory, 2 for write
 *                    error.                                    1/4/98-DWM */
{
  GV gv;
  HANDLE bmp;
  long freeimg=1, err, w, h, old=0, ow, oh;
  RECT rect;

  memset(&gv, 0, sizeof(GV));
  gv.name = OpenName;
  gv.format = type;
  switch (type) {
    case GV_PPM: gv.fs.ppm.ppmtype = PPM_BINARY; break;
    case GV_TIF: gv.fs.tif.compression = TIF_RLE; }
  GetClientRect(window_handle(data), &rect);
  w = rect.right;  h = rect.bottom;
  if (w<1) w = 1;  if (h<1) h = 1;
  if (((data->dflag>>9)&1) && data->w>0 && data->h>0) {
    w = data->w;  h = data->h; }
  if (w!=data->renderval[2] || h!=data->renderval[3]) {
    old = 1;  ow = data->renderval[2];  oh = data->renderval[3];
    render_move(0, 0, 0, w, h, data->renderval, data->renderdlt,
                data->renderdlt); }
  if (!(bmp=unlock2(update_bmp(data, w, h, &freeimg, 24))))
    return(1);
  if (old)
    render_move(0, 0, 0, ow, oh, data->renderval, data->renderdlt,
                data->renderdlt);
  err = !SaveBMP(&gv, bmp);
  free2(lock2(bmp));
  if (err)  return(2);
  return(0);
}

save_orb(DATA *data)
/* Save the current orbital as a specification file.  The file name is in
 *  OpenName.
 * The file format starts with the header "OrbitalFileV1.0".  This is
 *  followed by a series of keywords each followed by either a number, a
 *  string, or { } containing additional keywords and predicates.
 * Enter: DATA *data: pointer to window data.
 * Exit:  long error: 0 for okay, 1 for insufficient memory, 2 for write
 *                    error.                                    1/9/98-DWM */
{
  char *text, *dif;
  FILE *fptr;

  if (!(fptr=fopen(OpenName, "wt")))  return(2);
  save_orb_write(data, fptr, 0);
  fclose(fptr);
  return(0);
}

save_orb_write(DATA *data, FILE *fptr, long mode)
/* Write an orbital specification to a file.
 * The type variable contains: 0-none, 1-text, 2-lval with decimal, 3-rval,
 *  4-lval with hex, 5-lval with phrase, 6-lval with char.
 * Enter: DATA *data: pointer to window data.
 *        FILE *fptr: pointer to file to write to.
 *        long mode: bit 0: 0-header, 1-no header; bit 1: 0-sub orbitals,
 *                   1-no sub orbitals.                         1/9/98-DWM */
{
  long i, j=0, k=-1, type, lval, seq=0;
  double rval;
  char text[ONENAMELEN+2], *key, *sp;
  double phys[10];
  char space[]="    ";

  if (!(mode&1))  fprintf(fptr, "%s\n", OrbCheck);
  camera_rotate(data->renderval, data->renderdlt, phys, 0);
  for (i=0; OrbKey[i]; i++) {
    key = OrbKey[i];  sp = space+4-(2*((mode&2)!=0));  type = 0;
    switch (i) {
      case 0: type = 3;  rval = data->pref.perspec; break;
      case 1: case 2: case 3: case 4: case 5:
        if (data->pref.colors[i-1]>=0) {
          type = 4;  lval = data->pref.colors[i-1]; } break;
      case 6: type = 5;  lval = (data->dflag&1); break;
      case 7: if (!(data->dflag&1)) break;
        type = 5;  lval = ((data->dflag>>1)&3)+2; break;
      case 8: type = 5;  lval = ((data->dflag>>3)&3)+2; break;
      case 9: type = 5;  lval = ((data->dflag>>9)&1); break;
      case 10: if ((data->dflag>>9)&1)  type = 2; lval = data->w; break;
      case 11: if ((data->dflag>>9)&1)  type = 2; lval = data->h; break;
      case 12: type = 1;  sprintf(text, "\"%s\"", data->name); break;
      case 13: type = 3;  rval = fabs(data->renderval[0]); break;
      case 14: type = 3;  rval = data->renderval[1]; break;
      case 15: type = 3;  rval = data->renderval[2]; break;
      case 16: type = 3;  rval = data->renderval[3]; break;
      case 17: case 18: case 19: type = 3;  rval = phys[i-17]; break;
      case 20: case 21: case 22: type = 3;  rval = phys[i+7-20]; break;
      case 23: type = 3;  rval = phys[3]; break;
      case 24: for (; j<data->mol.nump; j++) {
        if (k<0)  fprintf(fptr, "%s%s {     # Atom %d\n", sp, key, j+1);
        k++;
        if (!AtomKey[k]) {
          k = -1;  fprintf(fptr, "%s  }\n", sp);  continue; }
        key = AtomKey[k];  sp -= 2;  i--;
        switch (k) {
          case 0: type = 2;  lval = data->mol.orb[j].n; break;
          case 1: type = 2;  lval = data->mol.orb[j].l;
            if (lval<strlen(OrbLet)) {
              type = 6; lval = OrbLet[lval]; } break;
          case 2: type = 2;  lval = data->mol.orb[j].m; break;
          case 3: type = 2;  lval = data->mol.orb[j].N; break;
          case 4: type = 2;  lval = data->mol.orb[j].Z; break;
          case 5: type = 3;  rval = data->mol.orb[j].mass; break;
          case 6: case 7: case 8:
            type = 3;  rval = data->mol.orb[j].x[k-6];
            if (!rval)  type = 0;  break;
          case 9: case 10: case 11:
            type = 3;  rval = data->mol.orb[j].ang[k-9];
            if (!rval)  type = 0;  break;
          case 12: type = 3;  rval = data->mol.orb[j].factor;
            if (rval==1)  type = 0;  break;
          default: type = 0; } break; } break;
      case 25: j = 0;  k = -1; break;          /* Used to reset for lights */
      case 26: if (!data->mol.numl) { type = 5;  lval = 5;  break; }
        for (; j<data->mol.numl; j++) {
          if (k<0)  fprintf(fptr, "%s%s {     # Light %d\n", sp, key, j+1);
          k++;
          if (!LightKey[k]) {
            k = -1;  fprintf(fptr, "%s  }\n", sp);  continue; }
          key = LightKey[k];  sp -= 2;  i--;
          switch (k) {
            case 0: case 1: case 2:
              type = 3;  rval = data->mol.ls[j].x[k]; break;
            case 3: type = 3;  rval = data->mol.ls[j].i; break;
            case 4: type = 3;  rval = data->mol.ls[j].a; break;
            case 5: type = 5;  lval = data->mol.ls[j].local; break;
            default: type = 0; } break; } break;
      case 27: j = 0;  k = -1; break;                     /* Used to reset */
      case 28: type = 3; rval = log10(data->mol.Psi); break;
      case 29: if (!1) { type = 5;  lval = 5;  break; }        /** numcut **/
        for (; j<1; j++) {                                     /** numcut **/
          if (k<0)  fprintf(fptr, "%s%s {     # Cutaway %d\n", sp, key, j+1);
          k++;
          if (!CutKey[k]) {
            k = -1;  fprintf(fptr, "%s  }\n", sp);  continue; }
          key = CutKey[k];  sp -= 2;  i--;
          switch (k) {
            case 0: type = 5;  lval = data->cut.type+5; break;
            case 1: type = 5;  lval = !data->cut.nosurface; break;
            case 2: type = 5;  lval = data->cut.invert; break;
            case 3: case 4: case 5:  type = 3;  rval = data->cut.xyz[k-3];
              if (!rval)  type = 0;  break;
            case 6: case 7: case 8:  type = 3;  rval = data->cut.ang[k-6];
              if (!rval)  type = 0;  break;
            default: type = 0; } break; } break;
      case 30: j = 0;  k = -1; break;                     /* Used to reset */
      case 31: for (; j<1; j++) {
        if (k<0)  fprintf(fptr, "%s%s {\n", sp, key);
        k++;
        if (!StereoKey[k]) {
          k = -1;  fprintf(fptr, "%s  }\n", sp);  continue; }
        key = StereoKey[k];  sp -= 2;  i--;
        switch (k) {
          case 0: type = 5;  lval = 9+data->stereo.mode; break;
          case 1: type = 5;  lval = 16+(data->stereo.flags&1);  break;
          case 2: type = 5;  lval = 16+((data->stereo.flags&2)!=0);  break;
          case 3: type = 5;  lval = ((data->stereo.flags&4)!=0);  break;
          case 4: type = 5;  lval = ((data->stereo.flags&8)!=0);  break;
          case 5: case 6:
            type = 3;  rval = data->stereo.interocular[k-5]; break;
          case 7: type = 3;  rval = data->stereo.separation; break;
          case 8: if (!data->stereo.name[0]) type = 0;
            else {
              type = 1;  sprintf(text, "\"%s\"", data->stereo.name); }
            break;
          default: type = 0; } break; } break;
      case 32: j = 0;  k = -1; break;
      case 33: for (; j<1; j++) {
        if (k<0)  fprintf(fptr, "%s%s {\n", sp, key);
        k++;
        if (!PointsKey[k]) {
          k = -1;  fprintf(fptr, "%s  }\n", sp);  continue; }
        key = PointsKey[k];  sp -= 2;  i--;
        switch (k) {
          case 0: type = 2;  lval = data->points.maxn; break;
          default: type = 0; } break; } break;
      case 34: j = 0;  k = -1; break;
      case 35: for (; j<1; j++) {
        if (k<0)  fprintf(fptr, "%s%s {\n", sp, key);
        k++;
        if (!PolyKey[k]) {
          k = -1;  fprintf(fptr, "%s  }\n", sp);  continue; }
        key = PolyKey[k];  sp -= 2;  i--;
        switch (k) {
          case 0: type = 2;  lval = data->poly.density; break;
          case 1: type = 3;  rval = data->poly.refine; break;
          case 2: case 3: type = 3;  rval = data->poly.opacity[k-2]; break;
          case 4: type = 5;  lval = 18+data->poly.wire; break;
          default: type = 0; } break; } break;
      case 36: j = 0;  k = -1; break;
      case 37: for (; j<1; j++) {
        if (k<0)  fprintf(fptr, "%s%s {\n", sp, key);
        k++;
        if (!AsymKey[k]) {
          k = -1;  fprintf(fptr, "%s  }\n", sp);  continue; }
        key = AsymKey[k];  sp -= 2;  i--;
        switch (k) {
          case 0: type = 2;  lval = data->asym.density; break;
          case 1: type = 3;  rval = data->asym.opacity; break;
          case 2: type = 5;  lval = 18+data->asym.wire; break;
          default: type = 0; } break; } break;
      case 38: j = 0;  k = -1; break;
      case 39: for (; j<1; j++) {
        if (k<0)  fprintf(fptr, "%s%s {\n", sp, key);
        k++;
        if (!RenderKey[k]) {
          k = -1;  fprintf(fptr, "%s  }\n", sp);  continue; }
        key = RenderKey[k];  sp -= 2;  i--;
        switch (k) {
          case 0: case 1: case 2: case 3: case 4: case 5: case 6:
            type = 3;  rval = data->render.opacity[k]; break;
          case 7: case 8: case 9:
            type = 3;  rval = data->render.refrac[k-7]; break;
          case 10: type = 2;  lval = data->render.steps; break;
          case 11: type = 5;  lval = (data->render.autobright!=0); break;
          case 12: type = 5;  lval = data->render.antialias&1; break;
          case 13: type = 5;  lval = ((data->render.antialias&2)!=0); break;
          default: type = 0; } break; } break;
      case 40: j = 0;  k = -1; break;
      case 41: if ((data->seq[0] && data->seq[1]) || (mode&2))  type = 2;
        lval = data->frame; break;
      case 42: if (data->seq[0] && data->seq[1] && !(mode&2))  type = 2;
        lval = data->lastframe; break;
      case 43: if (seq<4 && !(mode&2))
        if (data->seq[seq]) {
          fprintf(fptr, "%s%s {     # Sequence %d\n", sp, key, seq+1);
          save_orb_write(data->seq[seq], fptr, 3);
          fprintf(fptr, "%s  }\n", sp, key);
          i--;  seq++; } break;
      case 44: if (data->seq[0] && data->seq[1] && !(mode&2) &&
                   data->basename[0])  type = 1;
        sprintf(text, "\"%s\"", data->basename); break;
      case 45: if (data->seq[0] && data->seq[1] && !(mode&2))  type = 5;
        lval = data->incr; break;
      case 46: if (data->seq[0] && data->seq[1] && !(mode&2))  type = 5;
        lval = data->seqtype+21; break;
      case 48: if (data->seq[0] && data->seq[1] && !(mode&2))  type = 5;
        lval = data->bezier; break;
      case 49: if ((data->seq[0] && data->seq[1]) || (mode&2))  type = 2;
        lval = data->seqfps; break;
      default: ; }
    switch (type) {
      case 1: fprintf(fptr, "%s%-20s %s\n", sp, key, text); break;
      case 2: fprintf(fptr, "%s%-20s %d\n", sp, key, lval); break;
      case 3: fprintf(fptr, "%s%-20s %20.12lg\n", sp, key, rval); break;
      case 4: fprintf(fptr, "%s%-20s 0x%X\n", sp, key, lval); break;
      case 5: fprintf(fptr, "%s%-20s %s\n", sp, key, OrbPhrase[lval]); break;
      case 6: fprintf(fptr, "%s%-20s %c\n", sp, key, lval); } }
}

save_ov(DATA *data)
/* Save the current orbital as a complete orbital viewer file.  The file name
 *  is in OpenName.
 * This is the header "OrbitalViewerFileV1.0", followed by the exact same
 *  format as .ORB files (without the header), followed by four bytes of
 *  0xFF.  After this is as follows:
 *   Byte 0: sizeof(DATA)
 *        4: total length of stored file.
 *        8: additional file information (see code).
 *  Note that the sequence specification is solely stored in the textual part
 *  of the file.
 * Enter: DATA *data: pointer to window data.
 * Exit:  long error: 0 for okay, 1 for insufficient memory, 2 for write
 *                    error.                                    1/9/98-DWM */
{
  FILE *fptr;
  long m1=-1, header[2], hptr;

  if (!(fptr=fopen(OpenName, "wt")))  return(2);
  fprintf(fptr, "%s\n", OrbCheck2);
  save_orb_write(data, fptr, 1);
  fprintf(fptr, "%-20s %s\n", OrbKey[47], OrbPhrase[24]);
  fclose(fptr);
  if (!(fptr=fopen(OpenName, "r+b")))  return(2);
  fseek(fptr, 0, SEEK_END);
  fwrite(&m1, 4, 1, fptr);
  header[0] = sizeof(DATA);
  header[1] = sizeof(DATA) +
      data->mol.nump*sizeof(OATOM) + data->mol.numl*sizeof(LIGHT) +
      data->asym.n*3*sizeof(float) + data->asym.e*3*sizeof(long) +
      (data->asym.rfac!=0)*(data->asym.density+1)*sizeof(real) +
      (data->asym.tfac!=0)*(data->asym.density+1)*sizeof(real) +
      (data->asym.pfac!=0)*(data->asym.density*2+1)*sizeof(real) +
      data->poly.n*3*sizeof(float) + data->poly.e*3*sizeof(long) +
      data->poly.n*sizeof(long) + data->poly.e +
      (data->poly.rfac!=0)*(data->poly.density+1)*sizeof(real) +
      (data->poly.tfac!=0)*(data->poly.density+1)*sizeof(real) +
      (data->poly.pfac!=0)*(data->poly.density*2+1)*sizeof(real) +
      data->points.n*3*sizeof(float) + data->points.n +
      (data->render.buf!=0)*(-data->render.scan*data->render.h+40) +
      (data->render.zbuf!=0)*data->render.w*data->render.h*sizeof(ushort) +
      (data->render.phase!=0)*data->render.h*data->render.w*sizeof(long) +
      (data->stereo.image!=0)*(data->stereo.w*data->stereo.h*
                               (3-2*data->stereo.pal)+768*data->stereo.pal);
  fwrite(header, 4, 2, fptr);
  fwrite(data, sizeof(DATA), 1, fptr);
  fwrite(data->mol.orb, sizeof(OATOM), data->mol.nump, fptr);
  if (data->mol.numl)     fwrite(data->mol.ls, sizeof(LIGHT), data->mol.numl, fptr);
  if (data->asym.n)       fwrite(data->asym.xyz, sizeof(float)*3, data->asym.n, fptr);
  if (data->asym.e)       fwrite(data->asym.elem, sizeof(long)*3, data->asym.e, fptr);
  if (data->asym.rfac)    fwrite(data->asym.rfac, sizeof(real), data->asym.density+1, fptr);
  if (data->asym.tfac)    fwrite(data->asym.tfac, sizeof(real), data->asym.density+1, fptr);
  if (data->asym.pfac)    fwrite(data->asym.pfac, sizeof(real), data->asym.density*2+1, fptr);
  if (data->poly.n)       fwrite(data->poly.xyz, sizeof(float)*3, data->poly.n, fptr);
  if (data->poly.e)       fwrite(data->poly.elem, sizeof(long)*3, data->poly.e, fptr);
  if (data->poly.n)       fwrite(data->poly.ptphase, sizeof(long), data->poly.n, fptr);
  if (data->poly.e)       fwrite(data->poly.phase, 1, data->poly.e, fptr);
  if (data->poly.rfac)    fwrite(data->poly.rfac, sizeof(real), data->poly.density+1, fptr);
  if (data->poly.tfac)    fwrite(data->poly.tfac, sizeof(real), data->poly.density+1, fptr);
  if (data->poly.pfac)    fwrite(data->poly.pfac, sizeof(real), data->poly.density*2+1, fptr);
  if (data->points.n)     fwrite(data->points.xyz, sizeof(float)*3, data->points.n, fptr);
  if (data->points.n)     fwrite(data->points.phase, 1, data->points.n, fptr);
  if (data->render.buf)   fwrite(data->render.buf, 1, -data->render.scan*data->render.h+40, fptr);
  if (data->render.zbuf)  fwrite(data->render.zbuf, sizeof(ushort), data->render.w*data->render.h, fptr);
  if (data->render.phase) fwrite(data->render.phase, sizeof(long), data->render.w*data->render.h, fptr);
  if (data->stereo.image) fwrite(data->stereo.image, 1, data->stereo.w*data->stereo.h*(3-2*data->stereo.pal)+768*data->stereo.pal, fptr);
  fclose(fptr);
  return(0);
}

save_window(HWND hwnd, long saveas)
/* Save the data in the specified window.  If there is currently no filename,
 *  get one.
 * Enter: HWND hwnd: handle of window to save.  Null to save the topmost
 *                   window.
 *        long saveas: 0 for save, 1 for saveas.
 * Exit:  long saved: 1 for done, 0 for canceled or error.     2/13/97-DWM */
{
  DATA *data;
  long i, len, err=3, type;
  char dir[DIRLEN], name[NAMELEN], orig[NAMELEN], text[100];
  long out[]={FileOrb, FileOV, FilePPM, FileTIF, FileBMP, FileVRML,
              FileDigistar, -1};

  if (!hwnd)
    hwnd = SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0);
  if (!hwnd)  return(0);
  Busy+=10;
  data = lock_window(hwnd);
  GetWindowPlacement(hwnd, &data->place);
  strcpy(opfilter, SaveFilter);  strcat(opfilter, SaveFilterA);
  if (data->points.n || data->poly.n || data->asym.n)
    strcat(opfilter, SaveFilterB);
  if (data->points.n || data->poly.n || data->asym.n)
    strcat(opfilter, SaveFilterD);
  len = strlen(opfilter);
  for (i=0; i<len; i++)  if (opfilter[i]=='|')  opfilter[i] = 0;
  openfile.lStructSize = sizeof(OPENFILENAME);
  openfile.hwndOwner = Hwnd;  openfile.hInstance = hinst;
  openfile.lpstrFilter = opfilter;
  openfile.nFilterIndex = FilterType[1];
  openfile.lpstrFile = OpenName;
  openfile.nMaxFile = NAMELEN-1;
  openfile.lpstrFileTitle = openfile.lpstrCustomFilter = 0;
  openfile.lpstrInitialDir = OpenDir+DIRLEN;
  openfile.lpstrTitle = 0;
  openfile.lpfnHook = save_dialog;
  openfile.Flags = OFN_ENABLEHOOK|OFN_HIDEREADONLY|OFN_PATHMUSTEXIST|
                   OFN_EXPLORER|(OFN_OVERWRITEPROMPT*Pref.warn);
  memset(OpenName, 0xFF, NAMELEN);
  OpenName[0] = 0;
  strcpy(orig, data->name);
  if (strnicmp(orig, "Untitled", 8)) {
    if (strchr(orig, '\\')) {
      strcpy(dir, orig);
      strrchr(dir, '\\')[1] = 0;
      if (dir[strlen(dir)-2]!=':')  dir[strlen(dir)-1] = 0;
      if (saveas!=1)
        openfile.lpstrInitialDir = dir;
      strcpy(OpenName, strrchr(orig, '\\')+1); }
    else
      strcpy(OpenName, orig); }
  else
    saveas = 1;
  if (OpenName[0])
    for (i=0; out[i]>=0; i++)
      if (data->filetype==out[i])
        openfile.nFilterIndex = i+1;
  if (saveas) {
    if (!GetSaveFileName(&openfile)) {
      Busy-=10;  return(0); } }
  else
    strcpy(OpenName, data->name);
  FilterType[1] = openfile.nFilterIndex;
  if (strchr(OpenName, '\\'))
    strpath(OpenDir+DIRLEN, OpenName);
  if (type==FileOrb || type==FileOV)  strcpy(data->name, OpenName);
  type = out[FilterType[1]-1];
  cursor(1);
  switch (type) {
    case FileOrb: err = save_orb(data); break;
    case FileOV: err = save_ov(data); break;
    case FilePPM: err = save_graphic(data, GV_PPM); break;
    case FileTIF: err = save_graphic(data, GV_TIF); break;
    case FileBMP: err = save_graphic(data, GV_BMP); break;
    case FileVRML: err = save_vrml(data); break;
    case FileDigistar: err = save_digistar(data); break;
    /** Add additional file formats here **/
    default: ; }
  cursor(0);
  if (err)
    switch (err) {
      case 1: error(hwnd, "Insufficient memory."); break;
      case 2: error(hwnd, "Couldn't write file."); break;
      case 3: error(hwnd, "Function not implemented."); break;
      case 4: error(hwnd, "Window does not contain a forming limit curve.");
        break; }
  data->filetype = type;
  strcpy(name, OpenName);
  remove_quotes(name);
  getcwd(dir, NAMELEN-1);
  if (!strchr(name, '\\') && strlen(dir)+strlen(name)+2<NAMELEN) {
    sprintf(dir+strlen(dir), "\\%s", name);
    strcpy(name, dir); }
  for (i=0; i<strlen(name); i++)
    if (name[i]==toupper(name[i])) break;
  if (i!=strlen(name))
    for (i=0; i<strlen(name); i++)
      name[i] = toupper(name[i]);
  strcpy(data->name, name);
  set_window_name(hwnd, data);
  data->changed = 0;
  unlock_window(data);
  recheck(hwnd, 1);
  Busy-=10;
  return(1);
}

sequence_browse(HWND hdlg)
/* Obtain the base file name for sequence output.  This is stored in the
 *  appropriate spot in the dialog.
 * Enter: HWND hdlg: handle of dialog.                         1/13/98-DWM */
{
  long i, len;
  char *name;

  Busy += 10;
  strcpy(opfilter, SaveFilterA);
  strcat(opfilter, SaveFilterC);
  len = strlen(opfilter);
  for (i=0; i<len; i++)  if (opfilter[i]=='|')  opfilter[i] = 0;
  openfile.lStructSize = sizeof(OPENFILENAME);
  openfile.hwndOwner = Hwnd;  openfile.hInstance = hinst;
  openfile.lpstrFilter = opfilter;
  openfile.nFilterIndex = FilterType[4];
  openfile.lpstrFile = OpenName;
  openfile.nMaxFile = NAMELEN-1;
  openfile.lpstrFileTitle = openfile.lpstrCustomFilter = 0;
  openfile.lpstrInitialDir = OpenDir+DIRLEN*4;
  openfile.lpstrTitle = 0;  openfile.lpfnHook = 0;
  openfile.Flags = OFN_HIDEREADONLY|OFN_PATHMUSTEXIST|
                   OFN_EXPLORER|(OFN_OVERWRITEPROMPT*Pref.warn);
  memset(OpenName, 0xFF, NAMELEN);
  OpenName[0] = 0;
  GetDlgItemText(hdlg, SeqBase, OpenName, NAMELEN-1);
  if (!GetSaveFileName(&openfile)) {
    Busy -=10;  return(0); }
  name = OpenName;
  if (strchr(OpenName, '\\'))  name = strrchr(OpenName, '\\');
  if (strchr(name, '.'))       strchr(name, '.')[0] = 0;
  SetDlgItemText(hdlg, SeqBase, OpenName);
  FilterType[4] = openfile.nFilterIndex;
  if (strchr(OpenName, '\\'))
    strpath(OpenDir+DIRLEN*4, OpenName);
  SendDlgItemMessage(hdlg, SeqType, CB_SETCURSEL, FilterType[4]-1, 0);
  Busy -= 10;
}

DATA *sequence_file(HWND hdlg, DATA *seq)
/* Obtain an orbital specification for use with a sequence.
 * Enter: HWND hdlg: handle of dialog.
 *        DATA *seq: the current orbital specification.  This is freed if a
 *                   new one is loaded.
 * Exit:  DATA *seq: the new orbital specification, if changed.  Otherwise,
 *                   the original pointer.                     1/13/98-DWM */
{
  long i, len;
  DATA *new;
  FILE *fptr;

  strcpy(opfilter, OpenFilter);
  len = strlen(opfilter);
  for (i=0; i<len; i++)  if (opfilter[i]=='|')  opfilter[i] = 0;
  openfile.lStructSize = sizeof(OPENFILENAME);
  openfile.hwndOwner = hdlg;  openfile.hInstance = hinst;
  openfile.lpstrFilter = opfilter;
  openfile.nFilterIndex = FilterType[3];
  openfile.lpstrFile = OpenName;
  openfile.nMaxFile = NAMELEN-1;
  openfile.lpstrFileTitle = openfile.lpstrCustomFilter = 0;
  openfile.lpstrInitialDir = OpenDir+DIRLEN*3;
  openfile.Flags = OFN_FILEMUSTEXIST|OFN_HIDEREADONLY|OFN_PATHMUSTEXIST|
                   OFN_EXPLORER;
  memset(OpenName, 0xFF, NAMELEN);  OpenName[0] = 0;
  if (!GetOpenFileName(&openfile))
    return(seq);
  FilterType[3] = openfile.nFilterIndex;
  if (strchr(OpenName, '\\'))
    strpath(OpenDir+DIRLEN*3, OpenName);
  if (!(new=malloc2(sizeof(DATA)))) {
    error(hdlg, "Couldn't allocate memory.");  return(seq); }
  memset(new, 0, sizeof(DATA));
  if (!new_window_mid(new)) {
    free2(new);
    error(hdlg, "Couldn't allocate memory.");  return(seq); }
  if (!(fptr=fopen(OpenName, "rb"))) {
    free2(new->mol.orb);  free2(new->mol.ls);  free2(new);
    error(hdlg, "Couldn't read file.");  return(seq); }
  open_orb(fptr, new, 2);
  fclose(fptr);
  strcpy(new->name, OpenName);
  if (seq) {
    free2(seq->mol.orb);  free2(seq->mol.ls);  free2(seq); }
  return(new);
}

set_default()
/* Set the current orbital as the default orbital.             1/11/98-DWM */
{
  HWND hwnd;
  DATA *data;

  if (!(hwnd=SendMessage(HwndC,WM_MDIGETACTIVE,0,0))) return(0);
  if (!(data=lock_window(hwnd)))  return(0);
  strcpy(OpenName, DefFile);
  if (save_orb(data))
    error(hwnd, "Can't create default file.");
  else
    Pref.flags |= 1;
  unlock_window(data);
  if (!Pref.warn)  return(0);
  Busy += 10;
  MsgBox(hwnd, "The default orbital has been set.", "Set Default", MB_OK);
  Busy -= 10;
}

setup()
/* Configure the window to match the current user settings, which may have
 *  changed.  This adjusts the status bar and the tool bar.     2/5/97-DWM */
{
  RECT rect, rect2, rect3;
  FARPROC create;
  TBBUTTON *but;
  TBSAVEPARAMS tbs;
  long i, j, numbut, numtool, *toolinit=ToolInit;
  long partlist[]={1, 130, 200, -1};

  if (Pref.status && !HwndStat) {
    create = GetProcAddress(hcom, "CreateStatusWindow");
    HwndStat = (*create)(WS_CHILD|WS_VISIBLE|CCS_BOTTOM|SBARS_SIZEGRIP, "",
                         Hwnd, StatWindow);
    SendMessage(HwndStat, SB_SETPARTS, 4, (LPARAM)partlist); }
  if (!Pref.status && HwndStat) {
    DestroyWindow(HwndStat);  HwndStat = 0; }
  if (Pref.toolbar && !HwndTool) {
    if (ToolCustom)  toolinit = ToolCustom;
    for (i=numbut=0; ToolList[i].idCommand; i++)
      numbut = max(numbut, ToolList[i].iBitmap+1);
    for (numtool=0; toolinit[numtool]; numtool++);
    if (!(but=malloc2(numtool*sizeof(TBBUTTON)))) {
      error(Hwnd, "Low memory.\n\nCan't initialize toolbar.");
      Pref.toolbar = 0; } }
  if (Pref.toolbar && !HwndTool) {
    for (i=0; i<numtool; i++)
      if (toolinit[i]==-1)
        memcpy(but+i, &ToolSep, sizeof(TBBUTTON));
      else {
        for (j=0; j<numbut-1; j++)
          if (toolinit[i]==ToolList[j].idCommand)
                break;
        memcpy(but+i, ToolList+j, sizeof(TBBUTTON)); }
    create = GetProcAddress(hcom, "CreateToolbarEx");
    HwndTool = (*create)(Hwnd, WS_CHILD|CCS_ADJUSTABLE|TBSTYLE_ALTDRAG|
            TBSTYLE_WRAPABLE|TBSTYLE_TOOLTIPS|WS_VISIBLE, ToolWindow, numbut,
            hinst, TOOLBMP, but, numtool, 0, 0, 16, 16, sizeof(TBBUTTON));
    free2(but); }
  if (!Pref.toolbar && HwndTool) {
    write_ini();
    DestroyWindow(HwndTool);  HwndTool = 0; }
  compute_client(Hwnd, 0, 0);
  recheck(Hwnd, 0);
}

show_splash()
/* Open a small window in the center of the screen and show the splash
 *  message.                                                   11/5/97-DWM */
{
  HWND hdlg;
  long fw, x, y, w, h;
  GV gv;
  char name[9];
  HBITMAP hbmp, hbmp2;
  HDC hdc;
  LPLOGPALETTE lp;

  if (!((Pref.flags>>10)&1)) {
    fw = GetSystemMetrics(SM_CYDLGFRAME)*2;
    x = GetSystemMetrics(SM_CXSCREEN);
    y = GetSystemMetrics(SM_CYSCREEN);
    memset(&gv, 0, sizeof(GV));
    name[0] = 0;
    ((long *)(name+1))[0]=(long)(Splash+1);  ((long *)(name+1))[1]=Splash[0];
    gv.name = name;
    hbmp = LoadBMP(&gv);
    if (hbmp) {
      w = gv.width;  h = gv.height;
      hdlg = CreateDialog(hinst, SplashDLG, HWND_DESKTOP, splash_dialog);
      SetWindowPos(hdlg, HWND_TOPMOST, (x-w-fw)/2, (y-h-fw)/2, w+fw, h+fw,
                   SWP_NOACTIVATE);
      hdc = GetDC(hdlg);
      BitsPixel = GetDeviceCaps(hdc, BITSPIXEL);
      ReleaseDC(hdlg, hdc);
      if (BitsPixel<=8) {
        hbmp2 = PalettizeBMP(hbmp, 1, 0);
        if (hbmp2) {
          free2(lock2(hbmp));  hbmp = hbmp2; }
        lp = BMPPalette(hbmp);
        if (lp) {
          HpalSplash = CreatePalette(lp);  LocalFree(lp); } }
      ShowWindow(hdlg, SW_SHOWNOACTIVATE);
      HwndSplash = hdlg;
      SplashBMP = hbmp;
      UpdateWindow(hdlg);
      SetTimer(HwndSplash, 5, 5000, hide_splash); } }
  memset(&gv, 0, sizeof(GV));                /* Decompress preview picture */
  name[0] = 0;
  ((long *)(name+1))[0]=(long)(PreviewPic+1);
  ((long *)(name+1))[1]=PreviewPic[0];
  gv.name = name;
  PreviewGraphic = lock2(LoadGraphic(&gv));
}

BOOL CALLBACK splash_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Repaint the splash screen.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.
 * Exit:  long okay: 0 for cancelled, 1 for okay.              11/5/97-DWM */
{
  HDC hdc;
  PAINTSTRUCT paint;
  long out[4];

  switch (msg) {
    case WM_LBUTTONDOWN:  case WM_MBUTTONDOWN:  case WM_RBUTTONDOWN:
      hide_splash(0,0,0,0); break;
    case WM_PAINT: hdc = BeginPaint(hdlg, &paint);
      if (HpalSplash) {
        SelectPalette(hdc, HpalSplash, 0);
        RealizePalette(hdc); }
      map_dialog_rect(hdlg, 0, 0, 320, 240, out);
      draw_bmp(hdc, out[0], out[1], SplashBMP);
      EndPaint(hdlg, &paint); }
  return(0);
}

status(char *text)
/* Set the text of the status bar.
 * Enter: char *text: text to set.  Null to clear bar.         2/12/97-DWM */
{
  if (HwndStat) {
    if (text)  SendMessage(HwndStat, SB_SETTEXT, 3, (LPARAM)text);
    else       SendMessage(HwndStat, SB_SETTEXT, 3, (LPARAM)""); }
}

stereo_image(HWND hdlg, STEREO *stereo, long *replace)
/* Obtain an image for use with stereograms.
 * Enter: HWND hdlg: handle of dialog.
 *        STEREO *stereo: pointer to structure to store name, picture.
 *        long *replace: flag to indicate if this is the original picture in
 *                       the dialog (0) or not (1).             1/1/98-DWM */
{
  long i, len;
  uchar *img;
  GV gv;

  strcpy(opfilter, OpenFilter2);
  len = strlen(opfilter);
  for (i=0; i<len; i++)  if (opfilter[i]=='|')  opfilter[i] = 0;
  openfile.lStructSize = sizeof(OPENFILENAME);
  openfile.hwndOwner = hdlg;  openfile.hInstance = hinst;
  openfile.lpstrFilter = opfilter;
  openfile.nFilterIndex = FilterType[2];
  openfile.lpstrFile = OpenName;
  openfile.nMaxFile = NAMELEN-1;
  openfile.lpstrFileTitle = openfile.lpstrCustomFilter = 0;
  openfile.lpstrInitialDir = OpenDir+DIRLEN*2;
  openfile.Flags = OFN_FILEMUSTEXIST|OFN_HIDEREADONLY|OFN_PATHMUSTEXIST|
                   OFN_EXPLORER;
  memset(OpenName, 0xFF, NAMELEN);  OpenName[0] = 0;
  if (!GetOpenFileName(&openfile)) {
    if (!stereo->name[0])
      SendDlgItemMessage(hdlg, SterImage, BM_SETCHECK, 0, 0);
    return(0); }
  FilterType[2] = openfile.nFilterIndex;
  if (strchr(OpenName, '\\'))
    strpath(OpenDir+DIRLEN*2, OpenName);
  memset(&gv, 0, sizeof(GV));
  gv.name = OpenName;
  if (!(img=lock2(LoadGraphic(&gv)))) {
    error(hdlg, "The file is not a recognized graphics format.");
    if (!stereo->name[0])
      SendDlgItemMessage(hdlg, SterImage, BM_SETCHECK, 0, 0);
    return(0); }
  if (replace[0]) {
    free2(stereo->image);
    stereo->image = 0; }
  replace[0] = 1;
  strcpy(stereo->name, OpenName);
  stereo->w = gv.width;
  stereo->h = gv.height;
  stereo->pal = gv.palette;
  stereo->image = img;
  SetDlgItemText(hdlg, SterImageName, stereo->name);
}

strpath(char *dest, char *src)
/* If the source is different from the destination, copy the source to the
 *  destination.  Then, remove the trailing filename from the destination.
 *  This works properly for the c:\ type condition.
 * Enter: char *dest: location to store result.
 *        char *src: location of initial string.  Can be the same as dest.
 *                                                            11/15/96-DWM */
{
  if (src!=dest)
    strcpy(dest, src);
  remove_quotes(dest);
  if (strchr(dest, '\\'))
    strrchr(dest, '\\')[0] = 0;
  if (dest[strlen(dest)-1]==':')
    strcat(dest, "\\");
}

long windows_file(char *name)
/* Determine if a specified file is present in the WINDOWS/SYSTEM or
 *  SYSTEM32 path.  This is usually used with DLL files.
 * Enter: char *name: name of file to search for.  Wildcards are allowed.
 * Exit:  long present: 0 for not found, 1 for present.        10/2/96-DWM */
{
  char path[256], *cmd;

  return(SearchPath(0, name, 0, 255, path, &cmd));
}

write_ini()
/* Write values to the INI file.                               5/26/96-DWM */
{
  long i, j, numt;
  char *text, text2[32], t2[256], *alt;
  RECT rect;
  TBBUTTON tb;

  if (HwndTool) {
    numt = SendMessage(HwndTool, TB_BUTTONCOUNT, 0, 0);
    if (ToolCustom)  free2(ToolCustom);  ToolCustom = 0;
    if (ToolCustom=malloc2((numt+1)*sizeof(long))) {
      memset(ToolCustom, 0, (numt+1)*sizeof(long));
      for (i=0; i<numt; i++) {
        SendMessage(HwndTool, TB_GETBUTTON, i, (LPARAM)&tb);
        if (tb.idCommand>0)  ToolCustom[i] = tb.idCommand;
        else                 ToolCustom[i] = -1; } } }
  for (i=0; key[i*3]; i++) {
    text = key[i*3+1];  alt = key[i*3+2];
    switch (i) {
      case 0: text = t2;  t2[0] = 0;
        GetWindowPlacement(Hwnd, &WinPlace);
        for (j=0; j<sizeof(WINDOWPLACEMENT); j++)
          sprintf(t2+strlen(t2), "%02X", ((uchar *)(&WinPlace.length))[j]);
        break;
      case 1: if (!Pref.error)    text = alt;
        if (Pref.error==2) { sprintf(text2, "Status");  text = t2; } break;
      case 2: if (!Pref.warn)     text = alt;  break;
      case 3: if (!Pref.status)   text = alt;  break;
      case 4: if (!Pref.toolbar)  text = alt;  break;
      case 5: if (!Pref.tooltip)  text = alt;  break;
      case 6: if (!HwndTool || !ToolCustom)  continue;
        sprintf(text2, "%d", numt);  text = text2; break;
      case 7: text = OpenDir; break;
      case 8: text = t2;
        sprintf(t2, "%d,%d,%d,%d,%d,%d,%d", FilterType[0], FilterType[1],
                FilterType[2], FilterType[3], FilterType[4], FilterType[5],
                FilterType[6]); break;
      case 9: text = t2;  sprintf(t2, "%d", Pref.flags); break;
      case 10: text = t2;  sprintf(t2, "%g", Pref.perspec); break;
      case 11: text = t2;  t2[0] = 0;
        for (j=0; j<16 && j<NUMPREFCOLORS; j++)
          sprintf(t2+strlen(t2), "%06X,", Pref.colors[j]);
        t2[strlen(t2)-1] = 0; break;
      case 12: text = lastview;     break;
      case 13: text = lastview+256; break;
      case 14: text = lastview+512; break;
      case 15: text = lastview+768; break;
      case 16: text = OpenDir+DIRLEN; break;
      case 17: text = OpenDir+DIRLEN*2; break;
      case 18: text = t2;  sprintf(t2, "%g", Pref.radstep); break;
      case 19: text = t2;  sprintf(t2, "%g", Pref.panstep); break;
      case 20: text = t2;  sprintf(t2, "%g", Pref.zoomstep); break;
      case 21: text = t2;  sprintf(t2, "%d", Pref.pointsize*2+1); break;
      case 22: text = RadText[Pref.radunit]; break;
      case 23: text = DistText[Pref.munit]; break;
      case 24: text = MassText[Pref.kgunit]; break;
      case 25: text = OpenDir+DIRLEN*3; break;
      case 26: text = OpenDir+DIRLEN*4; break;
      case 27: text = OpenDir+DIRLEN*5; break;
      case 28: text = OpenDir+DIRLEN*6; break;
      default: continue; }
    WritePrivateProfileString(inihead[0], key[i*3], text, inifile); }
  write_ini_toolbar(numt);
}

write_ini_toolbar(long numt)
/* Write the current toolbar configuration to the ini file.
 * Enter: long numt: number of items in toolbar, including separators.
 *                                                            12/10/97-DWM */
{
  long i, j;
  char text2[32], t2[256];

  if (!HwndTool || !ToolCustom)  return(0);
  for (i=0; i<numt; i++) {
    sprintf(text2, "ToolBar%d", i);
    if (ToolCustom[i]<0)  strcpy(t2, "---");
    else {
      GetMenuString(Hmenu, ToolCustom[i], t2, 255, MF_BYCOMMAND);
      if (!strncmp(t2, "Redo", 4) || !strncmp(t2, "Undo", 4))
        strcpy(t2, "Undo");
      if (strchr(t2, '\t'))  strchr(t2, '\t')[0] = 0;
      for (j=0; j<strlen(t2); j++)
        if (t2[j]=='&' || t2[j]=='.' || t2[j]==' ') {
          memmove(t2+j, t2+j+1, strlen(t2+j));  j--; } }
    WritePrivateProfileString(inihead[1], text2, t2, inifile); }
}

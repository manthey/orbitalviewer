/* This contains the main windows functions and all file access functions for
 *  the Orbital Viewer program.                                 6/3/97-DWM */

#include <windows.h>
#include <commctrl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "matrix.h"
#include "ovrc.h"
#include "ov.h"
#include "data.h"
#include "\p\lib\gvlib.h"

char HelpFile[]="OV.HLP", lastview[1024], Program[]="Orbital Viewer",
     Untitled[]="Untitled", WinName[]="ORBwin", WinName2[]="ORBw2";
DATA *DispData;
HINSTANCE havi=0, hctl3d=0, hinst, hcom=0;
HMENU Hmenu;
HPALETTE Hpal=0, HpalSplash=0;
HWND DispHwnd, Hwnd, HwndC, HwndList[MAXWINDOWS], HwndSplash=0, HwndStat=0,
     HwndTool=0;
long AVILoad=0, BitsPixel, Busy=0, CloseMode=0, LockList[MAXLOCK*3],
     malloctbl[NUMMALLOC*4], mallocnum=0, NumLocked=0, NumWindows=0,
     PalChange=0, PreviewDelay=1500, ProgressCancel, ProgressPos,
     *ToolCustom=0;
PREF Pref;
WINDOWPLACEMENT WinPlace;
uchar *PreviewGraphic=0;
char IMSP[]="Insufficient memory to start program.\n";
real DistVal[]={1e-10, 1, 1e-9, a0};
char *DistText[]={"angstroms","meters","nanometers","Bohr atomic radii"};
real MassVal[]={ma, 1, 1e-3, 1e-27};
char *MassText[]={"amu","kilograms","grams","yoctograms"};
real RadVal[]={1,deg};
char *RadText[]={"radians","degrees"};
TBBUTTON ToolList[]={
    { 2,MenuNew,      TBSTATE_ENABLED,TBSTYLE_BUTTON},
    { 0,MenuOpen,     TBSTATE_ENABLED,TBSTYLE_BUTTON},
    { 1,MenuSave,     TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {21,MenuSaveAs,   TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {22,MenuClose,    TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {34,MenuCustom,   TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {42,MenuPref,     TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {33,MenuColors,   TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {29,MenuExt,      TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {30,MenuExit,     TBSTATE_ENABLED,TBSTYLE_BUTTON},
    { 4,MenuCopy,     TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {27,MenuRenderOpt,TBSTATE_ENABLED,TBSTYLE_BUTTON},
    { 5,MenuPointOpt, TBSTATE_ENABLED,TBSTYLE_BUTTON},
    { 6,MenuPolyOpt,  TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {26,MenuRayOpt,   TBSTATE_ENABLED,TBSTYLE_BUTTON},
    { 3,MenuAsymOpt,  TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {38,MenuOrb,      TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {23,MenuLight,    TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {41,MenuStereo,   TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {24,MenuCamera,   TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {40,MenuCutAway,  TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {44,MenuSequence, TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {45,MenuStop,     TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {46,MenuCompAVI,  TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {11,MenuLeft,     TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {12,MenuUp,       TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {13,MenuRight,    TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {14,MenuDown,     TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {15,MenuRotXM,    TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {16,MenuRotXP,    TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {17,MenuRotYM,    TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {18,MenuRotYP,    TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {19,MenuRotZP,    TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {20,MenuRotZM,    TBSTATE_ENABLED,TBSTYLE_BUTTON},
    { 9,MenuZoomIn,   TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {10,MenuZoomOut,  TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {36,MenuFocalIn,  TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {35,MenuFocalOut, TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {37,MenuResetPos, TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {39,MenuReframe,  TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {43,MenuDefault,  TBSTATE_ENABLED,TBSTYLE_BUTTON},
    { 8,MenuCascade,  TBSTATE_ENABLED,TBSTYLE_BUTTON},
    { 7,MenuTile,     TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {31,MenuArrange,  TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {25,MenuHelp,     TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {32,MenuHelpS,    TBSTATE_ENABLED,TBSTYLE_BUTTON},
    {28,MenuAbout,    TBSTATE_ENABLED,TBSTYLE_BUTTON},
    { 0,0}}; /* Next is 47 */
TBBUTTON ToolSep={0,0,TBSTATE_ENABLED,TBSTYLE_SEP,0,0,0};
long ToolInit[]={-1,MenuNew,MenuOpen,MenuSave,MenuSaveAs,MenuClose,-1,
     MenuCopy,-1,MenuRenderOpt,MenuPointOpt,MenuPolyOpt,MenuRayOpt,
     MenuAsymOpt,-1,MenuOrb,MenuLight,MenuStereo,MenuCamera,MenuCutAway,-1,
     MenuColors,-1,MenuResetPos,MenuReframe,-1,MenuTile,MenuCascade,-1,
     MenuHelp,0};

BOOL CALLBACK about(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the About dialog box.
 * Enter: HWND hwnd: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.        5/26/96-DWM */
{
  char text[80];

  switch (msg) {
    case WM_COMMAND: EndDialog(hdlg, 1); break;
    case WM_INITDIALOG:
      sprintf(text, "Version %d.%02d - %s", VERSION, SUBVERSION, __DATE__);
      SetDlgItemText(hdlg, AbDate, text);
      SetFocus(hdlg);
      return(1); }
  return(0);
}

advance_orbital()
/* Change the orbital to the next one in a logical sequence.  This also
 *  modifies the Psi^2 setting.                                1/31/98-DWM */
{
  HWND hwnd;
  DATA *data;
  double tbl[]={1.9,2.7,3.5,4.3,4.9,5.3,5.7,6.1,6.4,6.7,
                6.9,7.1,7.3,7.5,7.7,7.8,8.0,8.1,8.3,8.4, 8.5};

  if (!(hwnd=SendMessage(HwndC,WM_MDIGETACTIVE,0,0)))  return(0);
  if (!(data=lock_window(hwnd)))  return(0);
  if (data->mol.nump!=1) {      unlock_window(data);  return(0); }
  if (data->mol.orb[0].n>=21) { unlock_window(data);  return(0); }
  data->mol.orb[0].m++;
  if (data->mol.orb[0].m>data->mol.orb[0].l) {
    data->mol.orb[0].l++;  data->mol.orb[0].m = 0; }
  if (data->mol.orb[0].l>=data->mol.orb[0].n) {
    data->mol.orb[0].n++;  data->mol.orb[0].l = 0;  data->mol.orb[0].m = 0; }
  data->mol.Psi = pow(10, -tbl[data->mol.orb[0].n-1]);
  sprintf(data->name, "%d%c%d.ORB", data->mol.orb[0].n,
          OrbLet[data->mol.orb[0].l], data->mol.orb[0].m);
  data->asym.process = data->points.process = 0;
  data->poly.process = data->render.process = 0;
  orb_points(&data->points, 0, &data->mol, 0, 4);
  orb_polygons(&data->poly, 0, &data->mol, 0, 4);
  orb_asymptote(&data->asym, 0, &data->mol, 0, 4);
  orb_render(&data->render, 0, 0, &data->mol, 0, 4);
  data->changed = 1;
  unlock_window(data);
  recheck(hwnd, 0);
}

BOOL CALLBACK asymptote_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the Asymptote Options dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.         7/1/97-DWM */
{
  HWND hwnd;
  DATA *data;
  long x, y, i, j, scr;
  real opac, val;
  char text[80];

  hwnd = SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0);
  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case AsymAuto: if (!(data=lock_window(hwnd)))  break;
        for (i=j=0; i<data->mol.nump; i++) {
          if (data->mol.orb[i].n-data->mol.orb[i].l-1>j)
            j = data->mol.orb[i].n-data->mol.orb[i].l-1;
          if (data->mol.orb[i].l-fabs(data->mol.orb[i].m)-1>j)
            j = data->mol.orb[i].l-fabs(data->mol.orb[i].m)-1;
          if (fabs(data->mol.orb[i].m)>j)
            j = fabs(data->mol.orb[i].m); }
        j = j*2*data->mol.nump;  if (j<6)  j = 6;
        SetScrollPos(GetDlgItem(hdlg, AsymDScr), SB_CTL, (long)j/2, 1);
        sprintf(text, "%d", j);  SetDlgItemText(hdlg, AsymDens, text);
        unlock_window(data);  return(0);
      case AsymDens: GetDlgItemText(hdlg, AsymDens, text, 79);
        i = GetScrollPos(GetDlgItem(hdlg, AsymDScr), SB_CTL)*2;
        sscanf(text, "%d", &i);
        SetScrollPos(GetDlgItem(hdlg, AsymDScr), SB_CTL, i/2, 1);
        return(0);
      case AsymEdit: GetDlgItemText(hdlg, AsymEdit, text, 79);
        i = GetScrollPos(GetDlgItem(hdlg, AsymScr), SB_CTL);
        sscanf(text, "%d", &i);
        SetScrollPos(GetDlgItem(hdlg, AsymScr), SB_CTL, i, 1);
        return(0);
      case AsymWire: if (SendDlgItemMessage(hdlg,AsymWire,BM_GETCHECK,0,0))
          SendDlgItemMessage(hdlg,AsymWire2,BM_SETCHECK,0,0); break;
      case AsymWire2: if (SendDlgItemMessage(hdlg,AsymWire2,BM_GETCHECK,0,0))
          SendDlgItemMessage(hdlg,AsymWire,BM_SETCHECK,0,0); break;
      case HelpHelp: WinHelp(Hwnd, HelpFile, HELP_CONTEXT, HelpAsymptote);
        break;
      case IDOK: if (data=lock_window(hwnd)) {
        GetDlgItemText(hdlg, AsymEdit, text, 79);
        val = GetScrollPos(GetDlgItem(hdlg, AsymScr), SB_CTL);
        sscanf(text, "%lf", &val);
        GetDlgItemText(hdlg, AsymDens, text, 79);
        j = GetScrollPos(GetDlgItem(hdlg, AsymDScr), SB_CTL)*2;
        sscanf(text, "%d", &j);
        opac = val*0.01;
        if (opac!=data->asym.opacity || j!=data->asym.density) {
          data->changed = 1;
          if (data->asym.process>1 && j!=data->asym.density)
            data->asym.process = 1;
          if (data->points.process>=5)  data->points.process = 6;
          if (data->asym.process>=5)    data->asym.process = 6;
          if (data->poly.process>=6)    data->poly.process = 7; }
        if (opac!=data->render.opacity[6]) {
          data->changed = 1;
          orb_render(&data->render, 0, 0, &data->mol, 0, 4); }
        data->asym.opacity = data->render.opacity[6] = opac;
        data->asym.density = j;
        i = data->asym.wire;
        data->asym.wire = SendDlgItemMessage(hdlg,AsymWire,BM_GETCHECK,0,0)+
                        2*SendDlgItemMessage(hdlg,AsymWire2,BM_GETCHECK,0,0);
        if (i!=data->asym.wire) {
          if (data->asym.process>=5)    data->asym.process = 6; }
        unlock_window(data); }
      case IDCANCEL: EndDialog(hdlg, 1);  return(1);
      default: return(0); } break;
    case WM_HSCROLL: scr = GetDlgCtrlID(lp);
      GetScrollRange(lp, SB_CTL, &x, &y);
      i = GetScrollPos(lp, SB_CTL);
      switch (wp&0xFFFF) {
        case SB_BOTTOM: i = x; break;
        case SB_LINELEFT: if (i>x) i--; break;
        case SB_LINERIGHT: if (i<y) i++; break;
        case SB_PAGELEFT: i -=2*(1+4*(scr==AsymScr));  if (i<x) i = x; break;
        case SB_PAGERIGHT: i+=2*(1+4*(scr==AsymScr));  if (i>y) i = y; break;
        case SB_THUMBPOSITION: case SB_THUMBTRACK: i = (wp>>16); break;
        case SB_TOP: i = y; }
      SetScrollPos(lp, SB_CTL, i, 1);
      switch (scr) {
        case AsymScr: sprintf(text, "%d", i);
          SetDlgItemText(hdlg, AsymEdit, text); break;
        case AsymDScr: sprintf(text, "%d", i*2);
          SetDlgItemText(hdlg, AsymDens, text); } break;
    case WM_INITDIALOG: if (!(data=lock_window(hwnd)))
        return(1);
      if (data->asym.density<6)  data->asym.density = 6;
      SetScrollRange(GetDlgItem(hdlg, AsymScr), SB_CTL, 0, 100, 0);
      SetScrollPos(GetDlgItem(hdlg, AsymScr), SB_CTL,
                   (long)data->asym.opacity*100, 1);
      sprintf(text, "%g", data->asym.opacity*100);
      SetDlgItemText(hdlg, AsymEdit, text);
      SetScrollRange(GetDlgItem(hdlg, AsymDScr), SB_CTL, 3, 30, 0);
      SetScrollPos(GetDlgItem(hdlg, AsymDScr), SB_CTL,
                   (long)data->asym.density/2, 1);
      sprintf(text, "%d", data->asym.density);
      SetDlgItemText(hdlg, AsymDens, text);
      SendDlgItemMessage(hdlg, AsymWire, BM_SETCHECK, data->asym.wire==1, 0);
      SendDlgItemMessage(hdlg, AsymWire2,BM_SETCHECK, data->asym.wire==2, 0);
      SetFocus(hdlg);
      unlock_window(data);
      return(1); }
  return(0);
}

BOOL CALLBACK camera_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the Camera Options dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.        7/26/97-DWM */
{
  HWND hwnd=DispHwnd;
  DATA *data=DispData;
  char text[80], text2[80];
  static long du, ru, update, init;
  real val[3], x2;
  static double phys[10], lastphys[10];
  static real maxc;
  long i;
  HDC hdc;

  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case CamAng0: case CamAng1: case CamAng2: case CamPosX: case CamPosY:
      case CamPosZ:  if (init) break;
        for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, CamAng0+i, text, 79);
          x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= RadVal[ru];
          if (x2<1e10)  phys[i+7] = x2; }
        for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, CamPosX+i, text, 79);
          x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= DistVal[du];
          if (x2<1e10)  phys[i] = x2; }
        KillTimer(hdlg, 2);
        SetTimer(hdlg, 2, PreviewDelay, 0);
        update |= 1; break;
      case CamPosUnit: if ((wp>>16)!=CBN_SELCHANGE) break;
        du = SendDlgItemMessage(hdlg, CamPosUnit, CB_GETCURSEL, 0, 0);
        Pref.munit = du;  init = 1;
        for (i=0; i<3; i++) {
          sprintf(text, "%g", phys[i]/DistVal[du]);
          SetDlgItemText(hdlg, CamPosX+i, text); }
        sprintf(text, "%g %s", data->mol.maxcheck/DistVal[du], DistText[du]);
        SetDlgItemText(hdlg, CamRadius, text);  init = 0;  break;
      case CamRadUnit: if ((wp>>16)!=CBN_SELCHANGE) break;
        ru = SendDlgItemMessage(hdlg, CamRadUnit, CB_GETCURSEL, 0, 0);
        Pref.radunit = ru;  init = 1;
        for (i=0; i<3; i++) {
          x2 = phys[i+7]/RadVal[ru];
          if (fabs(x2)<1e-6)  x2 = 0;  sprintf(text, "%g", x2);
          sprintf(text2, "%g", x2+1e-6-2e-6*(x2<0));
          if (strlen(text2)<strlen(text)) strcpy(text, text2);
          SetDlgItemText(hdlg, CamAng0+i, text); }  init = 0;  break;
      case HelpHelp: WinHelp(Hwnd,HelpFile,HELP_CONTEXT,HelpCamera); break;
      case IDOK: for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, CamAng0+i, text, 79);
          x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= RadVal[ru];
          if (x2<1e10)  phys[i+7] = x2; }
        for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, CamPosX+i, text, 79);
          x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= DistVal[du];
          if (x2<1e10)  phys[i] = x2; }
        data->dflag = (data->dflag&0xFFFFFDFF)|(1<<9)*
                      SendDlgItemMessage(hdlg,CamFixed,BM_GETCHECK,0,0);
        GetDlgItemText(hdlg, CamWidth,  text, 79);
        i = data->w;  sscanf(text, "%d", &i);  if (i>0)  data->w = i;
        GetDlgItemText(hdlg, CamHeight, text, 79);
        i = data->h;  sscanf(text, "%d", &i);  if (i>0)  data->h = i;
        camera_rotate(data->renderval, data->renderdlt, phys, 1);
        update_process(data);
      case IDCANCEL:  KillTimer(hdlg, 2);
        EndDialog(hdlg, 1);  make_palette();  return(0);
      default: return(0); } break;
    case WM_INITDIALOG: update = 2;
      if (PreviewDelay<0 || PreviewDelay>3000)  PreviewDelay = 1000;
      KillTimer(hdlg, 2);
      camera_rotate(data->renderval, data->renderdlt, phys, 0);
      for (i=0; i<4; i++) {
        if (i<2) SendDlgItemMessage(hdlg, CamRadUnit, CB_ADDSTRING, 0,
                                    (long)RadText[i]);
        SendDlgItemMessage(hdlg, CamPosUnit, CB_ADDSTRING, 0,
                           (long)DistText[i]); }
      SendDlgItemMessage(hdlg,CamRadUnit, CB_SETCURSEL,Pref.radunit,0);
      SendDlgItemMessage(hdlg,CamPosUnit, CB_SETCURSEL,Pref.munit,  0);
      ru = Pref.radunit;  du = Pref.munit;
      init = 1;
      for (i=0; i<3; i++) {
        x2 = phys[i+7]/RadVal[ru];
        if (fabs(x2)<1e-6)  x2 = 0;  sprintf(text, "%g", x2);
        sprintf(text2, "%g", x2+1e-6-2e-6*(x2<0));
        if (strlen(text2)<strlen(text)) strcpy(text, text2);
        SetDlgItemText(hdlg, CamAng0+i, text); }
      for (i=0; i<3; i++) {
        sprintf(text, "%g", phys[i]/DistVal[du]);
        SetDlgItemText(hdlg, CamPosX+i, text); }
      init = 0;
      maxc = 0;
      SendDlgItemMessage(hdlg, CamFixed, BM_SETCHECK, (data->dflag>>9)&1, 0);
      sprintf(text, "%d", data->w);  SetDlgItemText(hdlg, CamWidth, text);
      sprintf(text, "%d", data->h);  SetDlgItemText(hdlg, CamHeight, text);
      camera_rect(0, hdlg);  memset(lastphys, 0, 10*sizeof(double));
      SetTimer(hdlg, 2, 0, 0);
      SetFocus(hdlg); break;
    case WM_LBUTTONDOWN: case WM_LBUTTONUP: case WM_MBUTTONDOWN:
    case WM_MBUTTONUP: case WM_MOUSEMOVE: case WM_RBUTTONDOWN:
    case WM_RBUTTONUP:
      if (preview_mouse(hdlg, msg, wp, lp, data, 0, (long)phys)) {
        init = 1;
        for (i=0; i<3; i++) {
          x2 = phys[i+7]/RadVal[ru];
          if (fabs(x2)<1e-6)  x2 = 0;  sprintf(text, "%g", x2);
          sprintf(text2, "%g", x2+1e-6-2e-6*(x2<0));
          if (strlen(text2)<strlen(text)) strcpy(text, text2);
          SetDlgItemText(hdlg, CamAng0+i, text); }
        for (i=0; i<3; i++) {
          sprintf(text, "%g", phys[i]/DistVal[du]);
          SetDlgItemText(hdlg, CamPosX+i, text); }
        init = 0;  camera_rect(0, hdlg); } break;
    case WM_PAINT: update = 2;
    case WM_TIMER: KillTimer(hdlg, 2);
      if (!update) return(0);
      if (!memcmp(phys, lastphys, 10*sizeof(double)) && update<2)  return(0);
      if (maxc!=data->mol.maxcheck) {
        sprintf(text, "%g %s", data->mol.maxcheck/DistVal[du], DistText[du]);
        SetDlgItemText(hdlg, CamRadius, text);
        maxc = data->mol.maxcheck; }
      memcpy(lastphys, phys, 10*sizeof(double));
      update = 0;
      hdc = GetDC(hdlg);
      camera_figure(hdlg, hdc, data, phys);
      ReleaseDC(hdlg, hdc);
      return(0); }
  return(0);
}

compute_client(HWND hwnd, long w, long h)
/* Adjust the MDI client window to account for the tool and status bars if
 *  they exist.
 * Enter: HWND hwnd: frame window handle.
 *        long w, h: size of the MDI frame window's client area.  If 0, this
 *                   is obtained using a GetClientRect call.    2/7/97-DWM */
{
  RECT rect, rect2, rect3;

  if (!w || !h) {
    GetClientRect(hwnd, &rect);
    w = rect.right;  h = rect.bottom; }
  if (HwndStat) {
    GetClientRect(Hwnd, &rect);
    GetWindowRect(HwndStat, &rect2);
    MoveWindow(HwndStat, 0, rect.bottom-rect2.bottom, w, h, 1);
    GetWindowRect(HwndStat, &rect2); }
  else
    rect2.top = rect2.bottom = 0;
  if (HwndTool) {
    GetWindowRect(HwndTool, &rect3);
    MoveWindow(HwndTool, 0, 0, w, rect3.bottom-rect3.top, 1);
    GetWindowRect(HwndTool, &rect3);
    InvalidateRect(HwndTool, 0, 1); }
  else
    rect3.top = rect3.bottom = 0;
  MoveWindow(HwndC, 0, rect3.bottom-rect3.top, w,
             h-(rect2.bottom-rect2.top)-(rect3.bottom-rect3.top), 1);
  InvalidateRect(hwnd, 0, 1);
}

cursor(long shape)
/* Set the cursor to a specified shape.
 * Enter: long shape: 0 for arrow, 1 for wait, 2 for move, 3 for rotate, 4
 *                    for zoom.                                5/26/96-DWM */
{
  HANDLE cur=0;
  HWND hwnd, hwnd2;

  switch (shape) {
    case 1:  SetCursor(cur=LoadCursor(0, IDC_WAIT)); break;
    case 2:  SetCursor(cur=LoadCursor(0, IDC_SIZE)); break;
    case 3:  SetCursor(cur=LoadCursor(hinst, CURSOR1)); break;
    case 4:  SetCursor(cur=LoadCursor(hinst, CURSOR2)); break;
    default: SetCursor(cur=LoadCursor(0, IDC_ARROW)); }
  if (cur) {
    hwnd = GetForegroundWindow();
    do {
      hwnd2 = GetTopWindow(hwnd);
      if (hwnd2)  hwnd = hwnd2; }
    while (hwnd2);
    SetClassLong(hwnd, GCL_HCURSOR, (long)cur); }
}

BOOL CALLBACK cutaway_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the Cutaway Options dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.        7/26/97-DWM */
{
  HWND hwnd=DispHwnd;
  DATA *data=DispData;
  char text[80], text2[80];
  static long du, ru, type, update, lasttype;
  static real pos[9], lastpos[9];
  real val[3], x2;
  long i;
  HDC hdc;

  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case CutAng0: case CutAng1: case CutAng2: case CutInvert: case CutPosX:
      case CutPosY: case CutPosZ: case CutSurface: case CutType0:
      case CutType1: case CutType2: case CutType3:
        for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, CutAng0+i, text, 79);
          x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= RadVal[ru];
          if (x2<1e10)  pos[i+3] = x2; }
        for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, CutPosX+i, text, 79);
          x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= DistVal[du];
          if (x2<1e10)  pos[i] = x2; }
        for (i=0; i<4; i++)
          if (SendDlgItemMessage(hdlg,CutType0+i,BM_GETCHECK,0,0))
            type = i;
        KillTimer(hdlg, 3);
        SetTimer(hdlg, 3, PreviewDelay, 0);
        update = 1; break;
      case CutPosUnit: for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, CutPosX+i, text, 79);
          x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= DistVal[du];
          if (x2<1e10)  val[i] = x2;
          else          val[i] = data->cut.xyz[i]; }
        du = SendDlgItemMessage(hdlg, CutPosUnit, CB_GETCURSEL, 0, 0);
        Pref.munit = du;
        for (i=0; i<3; i++) {
          sprintf(text, "%g", val[i]/DistVal[du=Pref.munit]);
          SetDlgItemText(hdlg, CutPosX+i, text); } break;
      case CutRadUnit: for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, CutAng0+i, text, 79);
          x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= RadVal[ru];
          if (x2<1e10)  val[i] = x2;
          else          val[i] = data->cut.ang[i]; }
        ru = SendDlgItemMessage(hdlg, CutRadUnit, CB_GETCURSEL, 0, 0);
        Pref.radunit = ru;
        for (i=0; i<3; i++) {
          x2 = val[i]/RadVal[ru];
          if (fabs(x2)<1e-6)  x2 = 0;  sprintf(text, "%g", x2);
          sprintf(text2, "%g", x2+1e-6-2e-6*(x2<0));
          if (strlen(text2)<strlen(text)) strcpy(text, text2);
          SetDlgItemText(hdlg, CutAng0+i, text); } break;
      case HelpHelp: WinHelp(Hwnd,HelpFile,HELP_CONTEXT,HelpCutaway); break;
      case IDOK: for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, CutAng0+i, text, 79);
          x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= RadVal[ru];
          if (x2<1e10)  data->cut.ang[i] = x2; }
        for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, CutPosX+i, text, 79);
          x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= DistVal[du];
          if (x2<1e10)  data->cut.xyz[i] = x2; }
        for (i=0; i<4; i++)
          if (SendDlgItemMessage(hdlg,CutType0+i,BM_GETCHECK,0,0))
            data->cut.type = i;
        data->cut.nosurface = !SendDlgItemMessage(hdlg, CutSurface,
                                                  BM_GETCHECK, 0, 0);
        data->cut.invert=SendDlgItemMessage(hdlg,CutInvert,BM_GETCHECK,0,0);
        memset(data->cut.mat, 0, 9*sizeof(real));
        if (data->asym.process>2)    data->asym.process = 2;
        if (data->points.process>2)  data->points.process = 2;
        if (data->poly.process>2)    data->poly.process = 2;
        if (data->render.process>=3) {
          data->render.x = data->render.y = 0;
          data->render.antialias = (data->render.antialias&0xFFFFFFF3)|
                                   (min(data->render.process-3,2)<<2);
          data->render.process = 3; }
      case IDCANCEL:  KillTimer(hdlg, 3);
        EndDialog(hdlg, 1);  make_palette();  return(0);
      default: return(0); } break;
    case WM_INITDIALOG: update = 1;
      if (PreviewDelay<0 || PreviewDelay>3000)  PreviewDelay = 1000;
      KillTimer(hdlg, 3);
      lasttype = -1;
      pos[6] = data->mol.maxcheck;
      if (pos[6]<=0)  pos[6] = a0;
      for (i=0; i<4; i++) {
        if (i<2) SendDlgItemMessage(hdlg, CutRadUnit, CB_ADDSTRING, 0,
                                    (long)RadText[i]);
        SendDlgItemMessage(hdlg, CutPosUnit, CB_ADDSTRING, 0,
                           (long)DistText[i]); }
      SendDlgItemMessage(hdlg,CutRadUnit, CB_SETCURSEL,Pref.radunit,0);
      SendDlgItemMessage(hdlg,CutPosUnit, CB_SETCURSEL,Pref.munit,  0);
      for (i=0; i<3; i++) {
        x2 = data->cut.ang[i]/RadVal[ru=Pref.radunit];
        if (fabs(x2)<1e-6)  x2 = 0;  sprintf(text, "%g", x2);
        sprintf(text2, "%g", x2+1e-6-2e-6*(x2<0));
        if (strlen(text2)<strlen(text)) strcpy(text, text2);
        SetDlgItemText(hdlg, CutAng0+i, text); }
      for (i=0; i<3; i++) {
        sprintf(text, "%g", data->cut.xyz[i]/DistVal[du=Pref.munit]);
        SetDlgItemText(hdlg, CutPosX+i, text); }
      SendDlgItemMessage(hdlg, CutType0+data->cut.type, BM_SETCHECK, 1, 0);
      memcpy(pos, data->cut.xyz, 3*sizeof(real));
      memcpy(pos+3, data->cut.ang, 3*sizeof(real));
      type = data->cut.type;
      SendDlgItemMessage(hdlg,CutSurface,BM_SETCHECK,!data->cut.nosurface,0);
      SendDlgItemMessage(hdlg,CutInvert,BM_SETCHECK,data->cut.invert,0);
      cutaway_rect(0, hdlg);  update = 1;
      SetTimer(hdlg, 3, 0, 0);
      SetFocus(hdlg); return(0);
    case WM_LBUTTONDOWN: case WM_LBUTTONUP: case WM_MBUTTONDOWN:
    case WM_MBUTTONUP: case WM_MOUSEMOVE: case WM_RBUTTONDOWN:
    case WM_RBUTTONUP:
      preview_mouse(hdlg, msg, wp, lp, data, 1, (long)pos);  return(0);
    case WM_PAINT: KillTimer(hdlg, 3);
      SetTimer(hdlg, 3, 0, 0);  update = 1;  lasttype = -1;  return(0);
    case WM_TIMER: KillTimer(hdlg, 3);
      if (!update) return(0);
      pos[7] = !SendDlgItemMessage(hdlg, CutSurface, BM_GETCHECK, 0, 0);
      pos[8] = SendDlgItemMessage(hdlg, CutInvert, BM_GETCHECK, 0, 0);
      if (!memcmp(pos, lastpos, 9*sizeof(real)) && lasttype==type)
        return(0);
      memcpy(lastpos, pos, 9*sizeof(real));
      lasttype = type;
      update = 0;
      hdc = GetDC(hdlg);
      cutaway_figure(hdlg, hdc, pos, type);
      ReleaseDC(hdlg, hdc);
      return(0); }
  return(0);
}

end_avi()
/* Terminate the AVI control functions, if started.            1/22/98-DWM */
{
  if (!havi)  return(0);
  FreeLibrary(havi);
  havi = 0;
}

end_common()
/* Terminate the common control functions, if started.         1/22/98-DWM */
{
  if (!hcom)  return(0);
  FreeLibrary(hcom);
  havi = 0;
}

end_ctl3d()
/* Terminate the faux-3D control functions, if started.        10/4/96-DWM */
{
  FARPROC reg;

  if (!hctl3d)  return(0);
  reg = GetProcAddress(hctl3d, "Ctl3dUnregister");
  (*reg)(hinst);
  FreeLibrary(hctl3d);
  hctl3d = 0;
}

error(HWND hwnd2, char *text)
/* Report an error, if errors are currently viewed.
 * Enter: HWND hwnd: handle of the window or null to use best-guess top
 *                   window.
 *        char *text: error text.                              2/10/97-DWM */
{
  char t2[256];
  static char *last;
  HWND hwnd;

  if (!hwnd2 && NumWindows)
    hwnd = GetTopWindow(0);
  if (!hwnd2 && !NumWindows)  hwnd = Hwnd;
  if (hwnd2)  hwnd = hwnd2;
  status(0);
  if (!Pref.error)  return(0);
  cursor(0);
  Busy += 10;
  if (Pref.error==1 || !HwndStat)
    MsgBox(hwnd, text, "Error", MB_OK);
  else {
    sprintf(t2, "Error: %s", text);
    status(t2); }
  Busy -= 10;
}

free_window(HWND hwnd)
/* Free all memory associated with a given window.  This is done just before
 *  closing.
 * Enter: HWND hwnd: handle of window to free.                 2/13/97-DWM */
{
  DATA *data;
  long i;

  if (!hwnd)  return(0);
  if (!(data=lock_window(hwnd))) return(0);
  for (i=0; i<NumLocked; i++)
    if (LockList[i*3+1]==data)
      break;
  if (i<NumLocked) {
    memmove(LockList+i*3, LockList+i*3+3, (NumLocked-i-1)*3*sizeof(long));
    NumLocked--; }
  orb_points(&data->points, 0, &data->mol, 0, 4);
  orb_polygons(&data->poly, 0, &data->mol, 0, 4);
  orb_asymptote(&data->asym, 0, &data->mol, 0, 4);
  orb_render(&data->render, 0, 0, &data->mol, 0, 4);
  free2(data->mol.orb);
  free2(data->mol.ls);
  free2(data->stereo.image);
  for (i=0; i<4; i++)  if (data->seq[i]) {
    free2(data->seq[i]->mol.orb);
    free2(data->seq[i]->mol.ls);
    free2(data->seq[i]); }
  /** Additional data items get freed here **/
  free2(data);
}

void free2(void *mem)
/* Actually call free.  Provided for ease and convenience.
 * Enter: void *mem: pointer to memory.                         5/4/96-DWM */
{
  long i;
  HGLOBAL hmem;

  if (!mem)  return;
  for (i=0; i<mallocnum; i++)
    if ((long)mem==malloctbl[i*2+1]) {
      hmem = (HGLOBAL)malloctbl[i*2];
      GlobalUnlock(hmem);
      GlobalFree(hmem);
      memmove(malloctbl+i*2, malloctbl+i*2+2,
              (mallocnum-i-1)*2*sizeof(long));
      mallocnum--;
      return; }
}

void free2_all(void)
/* Free all pointers allocated using malloc2.                  2/11/97-DWM */
{
  long i;

  for (i=mallocnum-1; i>=0; i--)
    free2(malloctbl[i*2+1]);
}

BOOL CALLBACK light_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the light source dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.         7/8/97-DWM */
{
  HWND hwnd=DispHwnd;
  DATA *data=DispData;
  long scr, init=0, x, y, i;
  char text[80];
  LIGHT *ls2;
  MOLECULE *mol;
  static LIGHT *ls;
  static long numl, cur, listed, du, update, last=0;
  static count=0;
  real val;
  HDC hdc;

  mol = &data->mol;
  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case HelpHelp: WinHelp(Hwnd,HelpFile,HELP_CONTEXT,HelpLight); break;
      case IDOK: KillTimer(hdlg, 4);
        if (listed>=0) { light_read(hdlg, ls+listed);  last = listed; }
        if (numl!=mol->numl || memcmp(mol->ls, ls, sizeof(LIGHT)*min(numl,
            mol->numl))) {
          if (data->asym.process>=5)    data->asym.process = 6;
          if (data->poly.process>=6)    data->poly.process = 7;
          orb_render(&data->render, 0, 0, &data->mol, 0, 4);
          for (i=0; i<numl; i++)  if (ls[i].a!=1) break;
          data->mol.minres = (0.1/data->mol.maxjumpadj)*((i==numl)+
                             0.04*(i!=numl));
          data->mol.back = (11.25/17.5/data->mol.EffScale/2/
                           data->mol.maxjumpadj)*((i==numl)+0.04*(i!=numl));
          data->changed = 1; }
        free2(mol->ls);
        mol->ls = ls;  mol->numl = numl;
        EndDialog(hdlg, 1);  return(0);
      case IDCANCEL: KillTimer(hdlg, 4);
        EndDialog(hdlg, 1);
        free2(ls);  return(1);
      case LightAdd: if (numl==MAXL) {
          error(hdlg, "The maximum number of light sources is already present.");
          return(0); }
        if (!(ls2=realloc2(ls, sizeof(LIGHT)*(numl+1)))) {
          error(hdlg, "Insufficient memory.");  return(0); }
        ls = ls2;
        if (listed>=0) light_read(hdlg, ls+listed);  listed = -1;
        if (cur>=0)
          memcpy(ls+numl, ls+cur, sizeof(LIGHT));
        else {
          ls[numl].x[0] = -80;  ls[numl].x[1] = ls[numl].x[2] = 80;
          ls[numl].i = ls[numl].a = 1;  ls[numl].local = 0; }
        numl++;
        cur = numl-1;
        init = 20; break;
      case LightDel: if (numl==0) {
          error(hdlg, "There are no light sources to delete.");  return(0); }
        if (!warning(hdlg, "Delete the current light source?"))  return(0);
        if (cur!=numl-1)
          memmove(ls+cur, ls+cur+1, (numl-cur-1)*sizeof(LIGHT));
        numl--;  listed = -1;
        if (cur==numl)  cur--;
        init = 20; break;
      case LightEditA: GetDlgItemText(hdlg, LightEditA, text, 79);
        val = GetScrollPos(GetDlgItem(hdlg, LightScrA), SB_CTL)*0.01;
        sscanf(text, "%lf", &val);
        SetScrollPos(GetDlgItem(hdlg, LightScrA), SB_CTL, val*100+0.5, 1);
        if (!update) {
          KillTimer(hdlg, 4);
          SetTimer(hdlg, 4, PreviewDelay, 0);  update = 1; }
        return(0);
      case LightEditI: GetDlgItemText(hdlg, LightEditI, text, 79);
        val = GetScrollPos(GetDlgItem(hdlg, LightScrI), SB_CTL)*0.01;
        sscanf(text, "%lf", &val);
        SetScrollPos(GetDlgItem(hdlg, LightScrI), SB_CTL, val*100+0.5, 1);
        if (!update) {
          KillTimer(hdlg, 4);
          SetTimer(hdlg, 4, PreviewDelay, 0);  update = 1; }
        return(0);
      case LightEditX: case LightEditY: case LightEditZ: KillTimer(hdlg, 4);
        SetTimer(hdlg, 4, PreviewDelay, 0);  update = 1;
        return(0);
      case LightList: if (listed>=0)  light_read(hdlg, ls+listed);
        else return(0);
        i = SendDlgItemMessage(hdlg, LightList, CB_GETCURSEL, 0, 0);
        if (i!=cur) {
          cur = i;
          init = 10; } break;
      case LightUnits: for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, LightEditX+i, text, 79);
          val = 1e30;  sscanf(text, "%lf", &val);  val *= DistVal[du];
          if (val<1e10)  ls[cur].x[i] = val; }
        du = SendDlgItemMessage(hdlg, LightUnits, CB_GETCURSEL, 0, 0);
        Pref.munit = du;
        for (i=0; i<3; i++) {
          sprintf(text, "%g", ls[cur].x[i]/DistVal[du=Pref.munit]);
          SetDlgItemText(hdlg, LightEditX+i, text); } break;
      default: return(0); } break;
    case WM_HSCROLL: scr = GetDlgCtrlID(lp);
      GetScrollRange(lp, SB_CTL, &x, &y);
      i = GetScrollPos(lp, SB_CTL);
      switch (wp&0xFFFF) {
        case SB_BOTTOM: i = x; break;
        case SB_LINELEFT: if (i>x) i--; break;
        case SB_LINERIGHT: if (i<y) i++; break;
        case SB_PAGELEFT: i -= 10;  if (i<x) i = x; break;
        case SB_PAGERIGHT: i += 10; if (i>y) i = y; break;
        case SB_THUMBPOSITION: case SB_THUMBTRACK: i = (wp>>16); break;
        case SB_TOP: i = y; }
      SetScrollPos(lp, SB_CTL, i, 1);
      sprintf(text, "%3.2f", i*0.01);
      SetDlgItemText(hdlg, scr-1, text); break;
    case WM_INITDIALOG: numl = mol->numl;  cur = last;  listed = -1;
      if (PreviewDelay<0 || PreviewDelay>3000)  PreviewDelay = 1000;
      KillTimer(hdlg, 4);
      if (!(ls=malloc2(sizeof(LIGHT)*(numl+1)))) {
        error(hdlg, "Insufficient memory.");
        EndDialog(hdlg, 1);  return(0); }
      memcpy(ls, mol->ls, sizeof(LIGHT)*mol->numl);
      SetScrollRange(GetDlgItem(hdlg, LightScrI), SB_CTL, 0, 100, 0);
      SetScrollRange(GetDlgItem(hdlg, LightScrA), SB_CTL, 0, 100, 0);
      for (i=0; i<4; i++)
        SendDlgItemMessage(hdlg, LightUnits, CB_ADDSTRING, 0,
                           (long)DistText[i]);
      SendDlgItemMessage(hdlg,LightUnits, CB_SETCURSEL, Pref.munit, 0);
      init = 20;
      SetFocus(hdlg); break;
    case WM_LBUTTONDOWN: case WM_LBUTTONUP: case WM_MBUTTONDOWN:
    case WM_MBUTTONUP: case WM_MOUSEMOVE: case WM_RBUTTONDOWN:
    case WM_RBUTTONUP:
      if (preview_mouse(hdlg, msg, wp, lp, data, 2, 0)) {
        count ^= 1;
        if (!count)  light_rect(0, hdlg); } return(0);
    case WM_PAINT: update = 1;
    case WM_TIMER: if (!update) return(0);
      KillTimer(hdlg, 4);
      update = 0;
      hdc = GetDC(hdlg);
      if (listed>=0 && listed<numl) {
        light_read(hdlg, ls+listed);
        light_figure(hdlg, hdc, ls+listed); }
      else
        light_figure(hdlg, hdc, 0);
      ReleaseDC(hdlg, hdc);
      return(0); }
  if (init>=20) {                                     /* Update light list */
    listed = -1;
    SendDlgItemMessage(hdlg, LightList, CB_RESETCONTENT, 0, 0);
    if (cur>=numl)  cur = numl-1;
    for (i=0; i<numl; i++) {
      sprintf(text, "Light %d", i+1);
      SendDlgItemMessage(hdlg, LightList, CB_ADDSTRING, 0, (long)text); }
    SendDlgItemMessage(hdlg, LightList, CB_SETCURSEL, cur, 0); }
  if (init>=10) {                           /* Update light-specific items */
    listed = -1;
    if (cur>=0 && numl>0) {
      SetScrollPos(GetDlgItem(hdlg, LightScrA), SB_CTL, ls[cur].a*100, 1);
      sprintf(text,"%3.2f",ls[cur].a);  SetDlgItemText(hdlg,LightEditA,text);
      SetScrollPos(GetDlgItem(hdlg, LightScrI), SB_CTL, ls[cur].i*100, 1);
      sprintf(text,"%3.2f",ls[cur].i);  SetDlgItemText(hdlg,LightEditI,text);
      SendDlgItemMessage(hdlg, LightLocal, BM_SETCHECK, ls[cur].local, 0);
      for (i=0; i<3; i++) {
        sprintf(text, "%g", ls[cur].x[i]/DistVal[du=Pref.munit]);
        SetDlgItemText(hdlg, LightEditX+i, text); } }
    KillTimer(hdlg, 4);
    SetTimer(hdlg, 4, PreviewDelay, 0);  update = 1;
    listed = cur; }
  return(0);
}

light_read(HWND hdlg, LIGHT *ls)
/* Read in the current settings of the dialog.  This is used when the user
 *  selects done or switches which light or units they are looking at.
 * Enter: HWND hdlg: handle of the dialog.
 *        LIGHT *ls: pointer to the current light source.     12/27/97-DWM */
{
  long i;
  char text[80];
  double val;

  GetDlgItemText(hdlg, LightEditA, text, 79);
  val = GetScrollPos(GetDlgItem(hdlg, LightScrA), SB_CTL)*0.01;
  sscanf(text, "%lf", &val);  ls->a = val;
  GetDlgItemText(hdlg, LightEditI, text, 79);
  val = GetScrollPos(GetDlgItem(hdlg, LightScrI), SB_CTL)*0.01;
  sscanf(text, "%lf", &val);  ls->i = val;
  i = SendDlgItemMessage(hdlg, LightUnits, CB_GETCURSEL, 0, 0);
  if (i>=0 && i<=3)  Pref.munit = i;
  ls->local = SendDlgItemMessage(hdlg, LightLocal, BM_GETCHECK, 0, 0);
  for (i=0; i<3; i++) {
    val = 1e30;
    GetDlgItemText(hdlg, LightEditX+i, text, 79);
    sscanf(text, "%lf", &val);  val *= DistVal[Pref.munit];
    if (val<1e20)  ls->x[i] = val; }
}

DATA *lock_window(HWND hwnd)
/* Lock down all of the memory areas of a window.
 * Enter: HWND hwnd: handle to window.
 * Exit:  DATA *data: pointer to main data array.              3/11/97-DWM */
{
  DATA *data;
  long i;
  HANDLE hnd;

  if (!hwnd)  return(0);
  for (i=0; i<NumLocked; i++)
    if (LockList[i*3]==hwnd) {
      LockList[i*3+2] ++;
      return(LockList[i*3+1]); }
  if (NumLocked==MAXLOCK)  return(0);
  if (!(hnd=GetWindowLong(hwnd, GWL_USERDATA)))  return(0);
  if (!(data=lock2(hnd)))  return(0);
  LockList[NumLocked*3]   = (long)(hwnd);
  LockList[NumLocked*3+1] = (long)(data);
  LockList[NumLocked*3+2] = 1;
  NumLocked++;
  data->mol.orb      = lock2(data->mol.orb);
  data->mol.ls       = lock2(data->mol.ls);
  data->asym.xyz     = lock2(data->asym.xyz);
  data->asym.elem    = lock2(data->asym.elem);
  data->asym.scrxyz  = lock2(data->asym.scrxyz);
  data->asym.norm    = lock2(data->asym.norm);
  data->asym.enorm   = lock2(data->asym.enorm);
  data->asym.rfac    = lock2(data->asym.rfac);
  data->asym.tfac    = lock2(data->asym.tfac);
  data->asym.pfac    = lock2(data->asym.pfac);
  data->points.xyz   = lock2(data->points.xyz);
  data->points.scrxyz= lock2(data->points.scrxyz);
  data->points.phase = lock2(data->points.phase);
  data->poly.xyz     = lock2(data->poly.xyz);
  data->poly.elem    = lock2(data->poly.elem);
  data->poly.scrxyz  = lock2(data->poly.scrxyz);
  data->poly.norm    = lock2(data->poly.norm);
  data->poly.enorm   = lock2(data->poly.enorm);
  data->poly.phase   = lock2(data->poly.phase);
  data->poly.ptphase = lock2(data->poly.ptphase);
  data->poly.rfac    = lock2(data->poly.rfac);
  data->poly.tfac    = lock2(data->poly.tfac);
  data->poly.pfac    = lock2(data->poly.pfac);
  data->render.buf   = lock2(data->render.buf);
  data->render.phase = lock2(data->render.phase);
  data->stereo.image = lock2(data->stereo.image);
  for (i=0; i<4; i++)  if (data->seq[i]) {
    data->seq[i]         = lock2(data->seq[i]);
    data->seq[i]->mol.orb = lock2(data->seq[i]->mol.orb);
    data->seq[i]->mol.ls  = lock2(data->seq[i]->mol.ls); }
  /** Additional data items get locked here **/
  return(data);
}

void *lock2(HANDLE hmem)
/* Lock down a windows handle so that it can be accessed via pointers and
 *  manipulated with the free2 and realloc2 commands.  This can only fail if
 *  there are too many pointers allocated.
 * Enter: HANDLE hmem: handle to windows memory.
 * Exit:  void *mem: pointer to memory.                        2/11/97-DWM */
{
  char *mem;
  long i;

  if (!hmem)  return(0);
  for (i=0; i<mallocnum; i++)
    if ((long)hmem==malloctbl[i*2])
      return(malloctbl[i*2+1]);
  if (mallocnum>=NUMMALLOC*2)  return(0);
  mem = GlobalLock(hmem);
  malloctbl[mallocnum*2]   = (long)hmem;
  malloctbl[mallocnum*2+1] = (long)mem;
  mallocnum++;
  return(mem);
}

LRESULT CALLBACK main_loop(HWND hwnd, ulong msg, WPARAM wp, LPARAM lp)
/* Main processing loop.
 * Enter: HWND hwnd: handle of current window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.         5/6/96-DWM */
{
  long i, j;
  static char text[256];
  char *t2;
  HWND hwnd2;

  switch (msg) {
    case WM_CLOSE: case WM_QUERYENDSESSION: if (NumWindows>1) CloseMode = 1;
      while (hwnd2=SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0)) {
        SendMessage(hwnd2, WM_CLOSE, 0, 0);
        if (hwnd2==SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0)) {
          CloseMode = 0;
          return(0); } }
      CloseMode = 0;
      write_ini();
      return(DefWindowProc(hwnd, msg, wp, lp));
    case WM_COMMAND: if (HwndStat) {
        LoadString(hinst, wp&0xFFFF, text, 255);
        status(text); }
      if (HwndSplash) hide_splash(0,0,0,0);
      switch (wp&0xFFFF) {
        case MenuAbout: Busy++; DialogBox(hinst,AbDLG,hwnd,about); Busy--;
          return(0);
        case MenuAdvance: advance_orbital(); return(0);
        case MenuArrange: SendMessage(HwndC,WM_MDIICONARRANGE,0,0);return(0);
        case MenuAsymOpt: Busy++;
          DialogBox(hinst, AsymDLG, hwnd, asymptote_dialog); Busy--;
          return(0);
        case MenuCamera:  DispHwnd = SendMessage(HwndC,WM_MDIGETACTIVE,0,0);
          if (!(DispData=lock_window(DispHwnd)))  return(0);
          Busy++;
          DialogBox(hinst, CamDLG, hwnd, camera_dialog); Busy--;
          unlock_window(DispData); return(0);
        case MenuCascade: SendMessage(HwndC, WM_MDICASCADE, 0, 0); return(0);
        case MenuClose: SendMessage(SendMessage(HwndC, WM_MDIGETACTIVE, 0,0),
                                    WM_CLOSE, 0, 0); return(0);
        case MenuColors: Busy++;
          DialogBox(hinst, ColorDLG, hwnd, color_dialog); Busy--; return(0);
        case MenuCompAVI: Busy+=10;  compress_avi();  Busy-=10;  return(0);
        case MenuCopy: copy();  return(0);
        case MenuCustom: if (HwndTool)
            SendMessage(HwndTool, TB_CUSTOMIZE, 0, 0);
          else
            error(hwnd, "Toolbar must be visible to allow customization.");
          return(0);
        case MenuCutAway:  DispHwnd = SendMessage(HwndC,WM_MDIGETACTIVE,0,0);
          if (!(DispData=lock_window(DispHwnd)))  return(0);
          Busy++;
          DialogBox(hinst, CutDLG, hwnd, cutaway_dialog); Busy--;
          unlock_window(DispData); return(0);
        case MenuDefault: set_default(); return(0);
        case MenuDown: shift_geo(1, Pref.panstep); return(0);
        case MenuExit: PostMessage(hwnd,WM_SYSCOMMAND,SC_CLOSE,0); return(0);
        case MenuExt: add_extensions(hwnd); return(0);
        case MenuFocalIn: focal_geo(2, 0.5); return(0);
        case MenuFocalOut: focal_geo(2, -0.5); return(0);
        case MenuHelp: WinHelp(hwnd, HelpFile, HELP_CONTENTS, 0); return(0);
        case MenuHelpS: WinHelp(hwnd, HelpFile, HELP_FINDER, 0); return(0);
        case MenuLastView: case MenuLastView+1: case MenuLastView+2:
        case MenuLastView+3:
          sprintf(OpenName,"\"%s\"",lastview+256*((wp&0xFFFF)-MenuLastView));
          open_window_mid(); break;
        case MenuLeft: shift_geo(0, -Pref.panstep); return(0);
        case MenuLight: DispHwnd = SendMessage(HwndC,WM_MDIGETACTIVE,0,0);
          if (!(DispData=lock_window(DispHwnd)))  return(0);
          Busy++;
          DialogBox(hinst, LightDLG, hwnd, light_dialog); Busy--;
          unlock_window(DispData); return(0);
        case MenuNew: new_window(); return(0);
        case MenuOpen: open_window(); return(0);
        case MenuOrb: DispHwnd = SendMessage(HwndC,WM_MDIGETACTIVE,0,0);
          if (!(DispData=lock_window(DispHwnd)))  return(0);
          Busy++;
          DialogBox(hinst, OrbDLG, hwnd, orbital_dialog); Busy--;
          unlock_window(DispData); return(0);
        case MenuPointOpt: Busy++;
          DialogBox(hinst, PointDLG, hwnd, point_dialog); Busy--; return(0);
        case MenuPolyOpt: Busy++;
          DialogBox(hinst, PolyDLG, hwnd, polygon_dialog); Busy--; return(0);
        case MenuPref: Busy++;
          DialogBox(hinst, PrefDLG, hwnd, preferences_dialog); Busy--;
          return(0);
        case MenuRayOpt: Busy++;
          DialogBox(hinst, RenderDLG, hwnd, render_opt_dialog); Busy--;
          return(0);
        case MenuRedirect: process_arguments();  return(0);
        case MenuReframe: reframe_geo();  return(0);
        case MenuRenderOpt: Busy++;
          DialogBox(hinst, RendDLG, hwnd, render_dialog); Busy--; return(0);
        case MenuResetPos: reset_geo(); return(0);
        case MenuRight: shift_geo(0, Pref.panstep); return(0);
        case MenuRotXM: rotate_geo(0, -Pref.radstep/deg); return(0);
        case MenuRotXP: rotate_geo(0,  Pref.radstep/deg); return(0);
        case MenuRotYM: rotate_geo(1, -Pref.radstep/deg); return(0);
        case MenuRotYP: rotate_geo(1,  Pref.radstep/deg); return(0);
        case MenuRotZM: rotate_geo(2, -Pref.radstep/deg); return(0);
        case MenuRotZP: rotate_geo(2,  Pref.radstep/deg); return(0);
        case MenuSave: save_window(0, 0); return(0);
        case MenuSaveAs: save_window(0, 1); return(0);
        case MenuSequence: DispHwnd = SendMessage(HwndC,WM_MDIGETACTIVE,0,0);
          if (!(DispData=lock_window(DispHwnd)))  return(0);
          Busy++;
          DialogBox(hinst, SeqDLG, hwnd, sequence_dialog); Busy--;
          unlock_window(DispData); return(0);
        case MenuStereo:  DispHwnd = SendMessage(HwndC,WM_MDIGETACTIVE,0,0);
          if (!(DispData=lock_window(DispHwnd)))  return(0);
          Busy++;
          DialogBox(hinst, SterDLG, hwnd, stereo_dialog); Busy--;
          unlock_window(DispData); return(0);
        case MenuStop: sequence_stop(); return(0);
        case MenuSwitch: SendMessage(HwndC, WM_MDINEXT, SendMessage(HwndC,
                                     WM_MDIGETACTIVE, 0, 0), 0); return(0);
        case MenuSwitchBack: SendMessage(HwndC, WM_MDINEXT, SendMessage(
                                HwndC, WM_MDIGETACTIVE, 0, 0), 1); return(0);
        case MenuTile: SendMessage(HwndC, WM_MDITILE, MDITILE_HORIZONTAL|
                                   MDITILE_VERTICAL, 0); return(0);
        case MenuUp: shift_geo(1, -Pref.panstep); return(0);
        case MenuZoom: if (IsZoomed(SendMessage(HwndC,WM_MDIGETACTIVE,0,0)))
            SendMessage(HwndC, WM_MDIRESTORE, SendMessage(HwndC,
                        WM_MDIGETACTIVE, 0, 0), 0);
          else
            SendMessage(HwndC, WM_MDIMAXIMIZE, SendMessage(HwndC,
                        WM_MDIGETACTIVE, 0, 0), 0);
          return(0);
        case MenuZoomIn: shift_geo(2, Pref.zoomstep); return(0);
        case MenuZoomOut: shift_geo(2, -Pref.zoomstep); return(0);
        default: return(DefFrameProc(hwnd, HwndC, msg, wp, lp)); } break;
    case WM_DESTROY: quit(); return(0);
    case WM_DROPFILES: drop_file(wp); return(0);
    case WM_INITMENU: if (HwndSplash) hide_splash(0,0,0,0); break;
    case WM_MDIACTIVATE: case WM_MDICASCADE: case WM_MDICREATE:
    case WM_MDIDESTROY: case WM_MDIGETACTIVE: case WM_MDIICONARRANGE:
    case WM_MDIMAXIMIZE: case WM_MDINEXT: case WM_MDIREFRESHMENU:
    case WM_MDIRESTORE: case WM_MDISETMENU: case WM_MDITILE:
      return(DefFrameProc(hwnd, HwndC, msg, wp, lp));
    case WM_MENUSELECT: if (!HwndStat) break;
      if ((wp>>16)&(MF_POPUP|MF_SYSMENU))
        text[0] = 0;
      else
        LoadString(hinst, wp&0xFFFF, text, 255);
      LastMove = 0;
      status(text); break;
    case WM_NOTIFY: if (((LPNMHDR)lp)->idFrom==ToolWindow)
      switch (((LPTBNOTIFY)lp)->hdr.code) {
        case TBN_CUSTHELP: WinHelp(hwnd, HelpFile, HELP_CONTENTS, 0); break;
        case TBN_ENDADJUST: compute_client(Hwnd, 0, 0);  break;
        case TBN_GETBUTTONINFO: for (i=0; ToolList[i].idCommand; i++);
          if (((LPTBNOTIFY)lp)->iItem>=i)  return(0);
          ((LPTBNOTIFY)lp)->tbButton = ToolList[((LPTBNOTIFY)lp)->iItem];
          GetMenuString(Hmenu, ToolList[((LPTBNOTIFY)lp)->iItem].idCommand,
                        text, 255, MF_BYCOMMAND);
          if (strchr(text, '\t'))  strchr(text, '\t')[0] = 0;
          if (t2=strchr(text, '&'))  memmove(t2, t2+1, strlen(t2));
          if (((LPTBNOTIFY)lp)->cchText>0) {
            if (strlen(text)+1>((LPTBNOTIFY)lp)->cchText)
              text[((LPTBNOTIFY)lp)->cchText-1] = 0;
            strcpy(((LPTBNOTIFY)lp)->pszText, text);
            ((LPTBNOTIFY)lp)->cchText = strlen(((LPTBNOTIFY)lp)->pszText); }
          return(1);
        case TBN_QUERYDELETE: return(1);
        case TBN_QUERYINSERT: return(1);
        case TBN_RESET: j = SendMessage(HwndTool, TB_BUTTONCOUNT, 0, 0);
          for (i=0; i<j; i++)
            SendMessage(HwndTool, TB_DELETEBUTTON, 0, 0);
          for (i=0; ToolInit[i]; i++)
            if (ToolInit[i]==-1)
              SendMessage(HwndTool, TB_ADDBUTTONS, 1, (LPARAM)&ToolSep);
            else {
              for (j=0; ToolList[j].idCommand &&
                   ToolInit[i]!=ToolList[j].idCommand; j++);
              SendMessage(HwndTool, TB_ADDBUTTONS, 1,(LPARAM)&ToolList[j]); }
          break; }
      if (HwndTool)
        switch (((LPTOOLTIPTEXT)lp)->hdr.code) {
          case TTN_NEEDTEXT: if (HwndStat) {
              LoadString(hinst, ((LPTOOLTIPTEXT)lp)->hdr.idFrom, text, 255);
              status(text); }
            if (!Pref.tooltip) break;
            GetMenuString(Hmenu, ((LPTOOLTIPTEXT)lp)->hdr.idFrom, text, 255,
                          MF_BYCOMMAND);
            if (strchr(text, '\t'))  strchr(text, '\t')[0] = 0;
            ((LPTOOLTIPTEXT)lp)->lpszText = text;
            return(1);
          default: ; } break;
    case WM_PALETTECHANGED: PalChange = 1; break;
    case WM_QUERYNEWPALETTE: if (!Hpal)  break;
      if (PalChange)  InvalidateRect(hwnd, 0, 0);
      PalChange = 0; break;
    case WM_SIZE: compute_client(hwnd, lp&0xFFFF, lp>>16); break;
    default: return(DefWindowProc(hwnd, msg, wp, lp)); }
  return(DefWindowProc(hwnd, msg, wp, lp));
}

void *malloc2(long size)
/* Actually call malloc.  Provided for source code compatibility.
 * Enter: long size: number of bytes to allocate.
 * Exit:  void *mem: pointer to memory.                         5/4/96-DWM */
{
  HGLOBAL hmem;
  char *mem;

  if (mallocnum>=NUMMALLOC || !size)  return(0);
  hmem = GlobalAlloc(GMEM_MOVEABLE, size);
  if (!hmem)  return(0);
  mem = GlobalLock(hmem);
  malloctbl[mallocnum*2]   = (long)hmem;
  malloctbl[mallocnum*2+1] = (long)mem;
  mallocnum++;
  return(mem);
}

map_dialog_rect(HWND hdlg, long x1, long y1, long x2, long y2, long *out)
/* Convert from the goofy dialog units to pixels.
 * Enter: HWND hdlg: handle of dialog with coordinates.
 *        long x1, y1, x2, y2: coordinates to transform.
 *        long *out: array of 4 longs for output.              1/12/97-DWM */
{
  RECT rect;

  rect.left = x1;   rect.top = y1;
  rect.right = x2;  rect.bottom = y2;
  MapDialogRect(hdlg, &rect);
  out[0] = rect.left;   out[1] = rect.top;
  out[2] = rect.right;  out[3] = rect.bottom;
}

long MsgBox(HWND hwnd, char *text, char *title, ulong style)
/* Perform the standard MessageBox procedure, but disable automatic functions
 *  while this is happening.
 * Enter: HWND hwnd: owner window.
 *        char *text: text within box.
 *        char *title: title on top of box.
 *        ulong style: see windows help for the meaning of this.
 * Exit:  long selected: value of selected button.            11/13/96-DWM */
{
  long sel;

  Busy++;
  BringWindowToTop(Hwnd);
  sel = MessageBox(hwnd, text, title, style);
  Busy--;
  return(sel);
}

next_frame(DATA *data)
/* Save the current finished frame, and advance to the next frame in a
 *  sequence.
 * Enter: DATA *data: pointer to window data.                  1/14/98-DWM */
{
  long typetbl[]={GV_PPM, GV_TIF, GV_BMP};
  char *ext[]={"PPM","TIF","BMP","AVI"};

  if (data->basename[0] && data->seqtype!=3) {
    sprintf(OpenName, "%s%d.%s", data->basename, data->frame,
            ext[data->seqtype]);
    save_graphic(data, typetbl[data->seqtype]); }
  if (data->basename[0] && data->seqtype==3) {
    sprintf(OpenName, "%s.%s", data->basename, ext[data->seqtype]);
    save_avi(data); }
  data->dflag &= 0xFFFFFBFF;
  if (data->frame<data->lastframe)       play_frame(data, data->frame+1);
  else if (data->frame>data->lastframe)  play_frame(data, data->frame-1);
}

BOOL CALLBACK orbital_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the Orbital dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.         7/8/97-DWM */
{
  HWND hwnd=DispHwnd;
  DATA *data=DispData;
  long x, y, i, j, scr, init=0;
  char text[80], text2[80];
  real x2, r[9], ang[6];
  OATOM *orb2;
  MOLECULE *mol;
  static OATOM *orb;
  static long nump, cur, listed, mu, du, ru, last=0;

  mol = &data->mol;
  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case HelpHelp: WinHelp(Hwnd,HelpFile,HELP_CONTEXT,HelpOrb); break;
      case IDOK:
        if (listed>=0) { orbital_read(hdlg, orb+listed);  last = listed; }
        if (nump!=mol->nump || memcmp(mol->orb, orb, sizeof(OATOM)*min(nump,
            mol->nump))) {
          data->asym.process = data->points.process = 0;
          data->poly.process = data->render.process = 0;
          orb_points(&data->points, 0, &data->mol, 0, 4);
          orb_polygons(&data->poly, 0, &data->mol, 0, 4);
          orb_asymptote(&data->asym, 0, &data->mol, 0, 4);
          orb_render(&data->render, 0, 0, &data->mol, 0, 4);
          data->changed = 1; }
        free2(mol->orb);
        mol->orb = orb;  mol->nump = nump;
        EndDialog(hdlg, 1);  return(0);
      case IDCANCEL: free2(orb); EndDialog(hdlg, 1);  return(1);
      case OrbAdd: if (nump==MAXP) {
          error(hdlg, "The maximum number of atoms is already present.");
          return(0); }
        if (!(orb2=realloc2(orb, sizeof(OATOM)*(nump+1)))) {
          error(hdlg, "Insufficient memory.");  return(0); }
        orb = orb2;
        if (listed>=0) orbital_read(hdlg, orb+listed);  listed = -1;
        memcpy(orb+nump, orb+cur, sizeof(OATOM));
        nump++;
        cur = nump-1;
        init = 20; break;
      case OrbDel: if (nump==1) {
          error(hdlg, "Can't delete the last atom.");  return(0); }
        if (!warning(hdlg, "Delete the current atom?"))  return(0);
        if (cur!=nump-1)
          memmove(orb+cur, orb+cur+1, (nump-cur-1)*sizeof(OATOM));
        nump--;  listed = -1;
        if (cur==nump)  cur--;
        init = 20; break;
      case OrbEditl: GetDlgItemText(hdlg, OrbEditl, text, 79);
        i = GetScrollPos(GetDlgItem(hdlg, OrbScrl), SB_CTL);
        sscanf(text, "%d", &i);
        for (j=0; j<strlen(OrbLet); j++)
          if (tolower(text[0])==OrbLet[j])  i = j;
        SetScrollPos(GetDlgItem(hdlg, OrbScrl), SB_CTL, i, 1); return(0);
      case OrbEditm: GetDlgItemText(hdlg, OrbEditm, text, 79);
        i = GetScrollPos(GetDlgItem(hdlg, OrbScrm), SB_CTL)-MAXN;
        sscanf(text, "%d", &i);
        SetScrollPos(GetDlgItem(hdlg, OrbScrm), SB_CTL, i+MAXN,1); return(0);
      case OrbEditn: GetDlgItemText(hdlg, OrbEditn, text, 79);
        i = GetScrollPos(GetDlgItem(hdlg, OrbScrn), SB_CTL);
        sscanf(text, "%d", &i);
        SetScrollPos(GetDlgItem(hdlg, OrbScrn), SB_CTL, i, 1); return(0);
      case OrbList: if (listed>=0) orbital_read(hdlg, orb+listed);
        else return(0);
        i = SendDlgItemMessage(hdlg, OrbList, CB_GETCURSEL, 0, 0);
        if (i!=cur) {
          cur = i;
          init = 10; } break;
      case OrbMassUnit: GetDlgItemText(hdlg, OrbEditMass, text, 79);
        x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= MassVal[mu];
        if (x2>0 && x2<1e-10)  orb[cur].mass = x2;
        mu = SendDlgItemMessage(hdlg, OrbMassUnit, CB_GETCURSEL, 0, 0);
        Pref.kgunit = mu;
        sprintf(text, "%g", orb[cur].mass/MassVal[mu]);
        SetDlgItemText(hdlg, OrbEditMass, text); break;
      case OrbRadUnit: for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, OrbEditt+i, text, 79);
          x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= RadVal[ru];
          if (x2<1e10)  orb[cur].ang[i] = x2; }
        ru = SendDlgItemMessage(hdlg, OrbRadUnit, CB_GETCURSEL, 0, 0);
        Pref.radunit = ru;
        for (i=0; i<3; i++) {
          x2 = orb[cur].ang[i]/RadVal[ru=Pref.radunit];
          if (fabs(x2)<1e-6)  x2 = 0;  sprintf(text, "%g", x2);
          sprintf(text2, "%g", x2+1e-6-2e-6*(x2<0));
          if (strlen(text2)<strlen(text)) strcpy(text, text2);
          SetDlgItemText(hdlg, OrbEditt+i, text); } break;
      case OrbUnits: for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, OrbEditX+i, text, 79);
          x2 = 1e30;  sscanf(text, "%lf", &x2);  x2 *= DistVal[du];
          if (x2<1e10)  orb[cur].x[i] = x2; }
        du = SendDlgItemMessage(hdlg, OrbUnits, CB_GETCURSEL, 0, 0);
        Pref.munit = du;
        for (i=0; i<3; i++) {
          sprintf(text, "%g", orb[cur].x[i]/DistVal[du=Pref.munit]);
          SetDlgItemText(hdlg, OrbEditX+i, text); } break;
      default: return(0); } break;
    case WM_HSCROLL: scr = GetDlgCtrlID(lp);
      GetScrollRange(lp, SB_CTL, &x, &y);
      i = GetScrollPos(lp, SB_CTL);
      switch (wp&0xFFFF) {
        case SB_BOTTOM: i = x; break;
        case SB_LINELEFT: if (i>x) i--; break;
        case SB_LINERIGHT: if (i<y) i++; break;
        case SB_PAGELEFT: i -= 1;  if (i<x) i = x; break;
        case SB_PAGERIGHT: i += 1; if (i>y) i = y; break;
        case SB_THUMBPOSITION: case SB_THUMBTRACK: i = (wp>>16); break;
        case SB_TOP: i = y; }
      SetScrollPos(lp, SB_CTL, i, 1);
      switch (scr) {
        case OrbScrl: sprintf(text, "%d", i);
          if (scr==OrbScrl && i>=0 && i<strlen(OrbLet))
            sprintf(text, "%c", OrbLet[i]); break;
        case OrbScrm: sprintf(text, "%d", i-MAXN); break;
        case OrbScrn: sprintf(text, "%d", i); }
      SetDlgItemText(hdlg, scr-1, text); break;
    case WM_INITDIALOG: nump = mol->nump;  cur = last;  listed = -1;
      if (!(orb=malloc2(sizeof(OATOM)*nump))) {
        error(hdlg, "Insufficient memory.");
        EndDialog(hdlg, 1);  return(0); }
      memcpy(orb, mol->orb, sizeof(OATOM)*mol->nump);
      SetScrollRange(GetDlgItem(hdlg, OrbScrn), SB_CTL, 1, MAXN, 0);
      SetScrollRange(GetDlgItem(hdlg, OrbScrl), SB_CTL, 0, MAXN-1, 0);
      SetScrollRange(GetDlgItem(hdlg, OrbScrm), SB_CTL, 1, 2*MAXN-1, 0);
      for (i=0; i<4; i++) {
        SendDlgItemMessage(hdlg, OrbMassUnit, CB_ADDSTRING, 0,
                           (long)MassText[i]);
        if (i<2) SendDlgItemMessage(hdlg, OrbRadUnit, CB_ADDSTRING, 0,
                                    (long)RadText[i]);
        SendDlgItemMessage(hdlg, OrbUnits, CB_ADDSTRING, 0,
                           (long)DistText[i]); }
      SendDlgItemMessage(hdlg,OrbMassUnit,CB_SETCURSEL,Pref.kgunit, 0);
      SendDlgItemMessage(hdlg,OrbRadUnit, CB_SETCURSEL,Pref.radunit,0);
      SendDlgItemMessage(hdlg,OrbUnits,   CB_SETCURSEL,Pref.munit,  0);
      init = 20;
      SetFocus(hdlg); }
  if (init>=20) {                                      /* Update atom list */
    listed = -1;
    if (cur>=nump)  cur = nump-1;
    SendDlgItemMessage(hdlg, OrbList, CB_RESETCONTENT, 0, 0);
    for (i=0; i<nump; i++) {
      sprintf(text, "Atom %d", i+1);
      SendDlgItemMessage(hdlg, OrbList, CB_ADDSTRING, 0, (long)text); }
    SendDlgItemMessage(hdlg, OrbList, CB_SETCURSEL, cur, 0); }
  if (init>=10) {                            /* Update atom-specific items */
    listed = -1;
    SetScrollPos(GetDlgItem(hdlg, OrbScrn), SB_CTL, orb[cur].n, 1);
    sprintf(text, "%d", orb[cur].n);  SetDlgItemText(hdlg, OrbEditn, text);
    SetScrollPos(GetDlgItem(hdlg, OrbScrl), SB_CTL, orb[cur].l, 1);
    sprintf(text, "%d", orb[cur].l);
    if (orb[cur].l<strlen(OrbLet))  sprintf(text, "%c", OrbLet[orb[cur].l]);
    SetDlgItemText(hdlg, OrbEditl, text);
    SetScrollPos(GetDlgItem(hdlg, OrbScrm), SB_CTL,orb[cur].m+orb[cur].l, 1);
    sprintf(text, "%d", orb[cur].m);  SetDlgItemText(hdlg, OrbEditm, text);
    sprintf(text, "%d", orb[cur].Z);  SetDlgItemText(hdlg, OrbEditP, text);
    sprintf(text, "%g", orb[cur].mass/MassVal[mu=Pref.kgunit]);
    SetDlgItemText(hdlg, OrbEditMass, text);
    sprintf(text,"%g",orb[cur].factor);  SetDlgItemText(hdlg,OrbEditF,text);
    for (i=0; i<3; i++) {
      x2 = orb[cur].ang[i]/RadVal[ru=Pref.radunit];
      if (fabs(x2)<1e-6)  x2 = 0;  sprintf(text, "%g", x2);
      sprintf(text2, "%g", x2+1e-6-2e-6*(x2<0));
      if (strlen(text2)<strlen(text)) strcpy(text, text2);
      SetDlgItemText(hdlg, OrbEditt+i, text); }
    for (i=0; i<3; i++) {
      sprintf(text, "%g", orb[cur].x[i]/DistVal[du=Pref.munit]);
      SetDlgItemText(hdlg, OrbEditX+i, text); }
    listed = cur; }
  return(0);
}

orbital_read(HWND hdlg, OATOM *orb)
/* Read in the current settings of the dialog.  This is used when the user
 *  selects done or switches which atom or units they are looking at.
 * Enter: HWND hdlg: handle of the dialog.
 *        OATOM *orb: pointer to the current orbital.           7/10/97-DWM */
{
  long i;
  char text[80];
  double x2;

  GetDlgItemText(hdlg, OrbEditm, text, 79);
  i = GetScrollPos(GetDlgItem(hdlg, OrbScrm), SB_CTL)-orb->l;
  sscanf(text, "%d", &i);  orb->m = i;
  GetDlgItemText(hdlg, OrbEditl, text, 79);
  i = GetScrollPos(GetDlgItem(hdlg, OrbScrl), SB_CTL);
  sscanf(text, "%d", &i);  orb->l = i;
  GetDlgItemText(hdlg, OrbEditn, text, 79);
  i = GetScrollPos(GetDlgItem(hdlg, OrbScrn), SB_CTL);
  sscanf(text, "%d", &i);  orb->n = i;
  if (orb->n<=0)       orb->n = 0;
  if (orb->n>MAXN)     orb->n = MAXN;
  if (orb->l<0)        orb->l = 0;
  if (orb->l>=orb->n)  orb->l = orb->n-1;
  if (orb->m<-orb->l)  orb->m = -orb->l;
  if (orb->m>orb->l)   orb->m = orb->l;
  i = 0;
  GetDlgItemText(hdlg, OrbEditP, text, 79);
  sscanf(text, "%d", &i);  if (i>=1)  orb->Z = i;
  i = SendDlgItemMessage(hdlg, OrbMassUnit, CB_GETCURSEL, 0, 0);
  if (i>=0 && i<=3)  Pref.kgunit = i;
  x2 = 0;
  GetDlgItemText(hdlg, OrbEditMass, text, 79);
  sscanf(text, "%lf", &x2);  x2 *= MassVal[Pref.kgunit];
  if (x2>0 && x2<1e-10)  orb->mass = x2;
  x2 = 0;
  GetDlgItemText(hdlg, OrbEditF, text, 79);
  sscanf(text, "%lf", &x2);  if (x2)  orb->factor = x2;
  i = SendDlgItemMessage(hdlg, OrbRadUnit, CB_GETCURSEL, 0, 0);
  if (i>=0 && i<=1)  Pref.radunit = i;
  for (i=0; i<3; i++) {
    x2 = 1e30;
    GetDlgItemText(hdlg, OrbEditt+i, text, 79);
    sscanf(text, "%lf", &x2);  x2 *= RadVal[Pref.radunit];
    if (x2<1e20)  orb->ang[i] = x2; }
  i = SendDlgItemMessage(hdlg, OrbUnits, CB_GETCURSEL, 0, 0);
  if (i>=0 && i<=3)  Pref.munit = i;
  for (i=0; i<3; i++) {
    x2 = 1e30;
    GetDlgItemText(hdlg, OrbEditX+i, text, 79);
    sscanf(text, "%lf", &x2);  x2 *= DistVal[Pref.munit];
    if (x2<1e20)  orb->x[i] = x2; }
}

void pass_arguments(HWND hwnd, char *arg)
/* Pass the arguments on the command line to an already running version of
 *  the program.  This is done by placing the arguments on the clipboard,
 *  then sending a message to the already running version to get them.
 * Enter: HWND hwnd: handle to already running version of the program.
 *        char *arg: argument string to pass.  Can be null.    1/30/98-DWM */
{
  HANDLE gmem;
  char *mem;

  if (!arg)  return;
  if (!strlen(arg))  return;
  if (!(gmem=CreateFileMapping(0xFFFFFFFF, 0, PAGE_READWRITE, 0,
        strlen(arg)+1, "PassArguments")))  return;
  mem = MapViewOfFile(gmem, FILE_MAP_WRITE, 0, 0, 0);
  if (mem) {
    strcpy(mem, arg);
    SendMessage(hwnd, WM_COMMAND, MenuRedirect, 0);
    UnmapViewOfFile(mem); }
  CloseHandle(gmem);
}

BOOL CALLBACK point_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the Point Options dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.        9/26/97-DWM */
{
  HWND hwnd;
  DATA *data;
  long x, y, i;
  char text[80], *newp;
  real temp2;
  float *newn;

  hwnd = SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0);
  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case HelpHelp: WinHelp(Hwnd,HelpFile,HELP_CONTEXT,HelpPointOpt); break;
      case IDOK: if (data=lock_window(hwnd)) {
        GetDlgItemText(hdlg, PointPoint, text, 79);
        temp2 = GetScrollPos(GetDlgItem(hdlg, PointScr), SB_CTL)*1000;
        sscanf(text, "%lf", &temp2);
        if (temp2>0 && temp2!=data->points.maxn) {
          data->changed = 1;
          data->points.maxn = temp2;
          if (data->points.process>3)  data->points.process = 3; }
        unlock_window(data); }
      case IDCANCEL: EndDialog(hdlg, 1);  return(1);
      case PointOpt: DialogBox(hinst, AsymDLG, hdlg, asymptote_dialog);
        return(0);
      case PointPoint: GetDlgItemText(hdlg, PointPoint, text, 79);
        temp2 = GetScrollPos(GetDlgItem(hdlg, PointScr), SB_CTL)*1000;
        sscanf(text, "%lf", &temp2);
        SetScrollPos(GetDlgItem(hdlg, PointScr), SB_CTL,
                     (long)(temp2*0.001), 1);  return(0);
      default: return(0); } break;
    case WM_HSCROLL: GetScrollRange(lp, SB_CTL, &x, &y);
      i = GetScrollPos(lp, SB_CTL);
      switch (wp&0xFFFF) {
        case SB_BOTTOM: i = x; break;
        case SB_LINELEFT: if (i>x) i--; break;
        case SB_LINERIGHT: if (i<y) i++; break;
        case SB_PAGELEFT: i -= 10;  if (i<x) i = x; break;
        case SB_PAGERIGHT: i += 10; if (i>y) i = y; break;
        case SB_THUMBPOSITION: case SB_THUMBTRACK: i = (wp>>16); break;
        case SB_TOP: i = y; }
      SetScrollPos(lp, SB_CTL, i, 1);
      sprintf(text, "%d", i*1000);
      SetDlgItemText(hdlg, PointPoint, text);
      break;
    case WM_INITDIALOG: if (!(data=lock_window(hwnd)))
        return(1);
      SetScrollRange(GetDlgItem(hdlg, PointScr), SB_CTL, 1, 100,0);
      SetScrollPos(GetDlgItem(hdlg, PointScr), SB_CTL,
                   (long)(data->points.maxn*0.001), 1);
      sprintf(text, "%ld", data->points.maxn);
      SetDlgItemText(hdlg, PointPoint, text);
      SetFocus(hdlg);
      unlock_window(data);
      return(1); }
  return(0);
}

BOOL CALLBACK polygon_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the Polygon Options dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.       12/29/97-DWM */
{
  HWND hwnd;
  DATA *data;
  long x, y, i, j, scr;
  real opacp, opacn, refine, psi, psi2, val;
  char text[80];

  hwnd = SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0);
  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case HelpHelp: WinHelp(Hwnd, HelpFile, HELP_CONTEXT, HelpPolygon);
        break;
      case IDOK: if (data=lock_window(hwnd)) {
        GetDlgItemText(hdlg, PolyPsi, text, 79);
        psi = GetScrollPos(GetDlgItem(hdlg, PolyScrPsi), SB_CTL)*-0.1;
        sscanf(text, "%lf", &psi);  psi2 = pow(10, psi);
        GetDlgItemText(hdlg, PolyEditP, text, 79);
        val = GetScrollPos(GetDlgItem(hdlg, PolyScrP), SB_CTL);
        sscanf(text, "%lf", &val);  opacp = val*0.01;
        GetDlgItemText(hdlg, PolyEditN, text, 79);
        val = GetScrollPos(GetDlgItem(hdlg, PolyScrN), SB_CTL);
        sscanf(text, "%lf", &val);  opacn = val*0.01;
        GetDlgItemText(hdlg, PolyDens, text, 79);
        j = GetScrollPos(GetDlgItem(hdlg, PolyScrD), SB_CTL)*2;
        sscanf(text, "%d", &j);
        GetDlgItemText(hdlg, PolyRefine, text, 79);
        refine = GetScrollPos(GetDlgItem(hdlg, PolyScrR), SB_CTL)*0.01;
        sscanf(text, "%lf", &refine);
        if (opacp!=data->poly.opacity[0] || opacn!=data->poly.opacity[1] ||
            j!=data->poly.density || refine!=data->poly.refine ||
            psi2!=data->mol.Psi) {
          data->changed = 1;
          if (data->poly.process>1 && (j!=data->poly.density ||
              psi2!=data->mol.Psi))
            data->poly.process = 1;
          if (psi2!=data->mol.Psi)
            orb_render(&data->render, 0, 0, &data->mol, 0, 4);
          if (data->poly.process>=6)    data->poly.process = 5; }
        data->poly.opacity[0] = opacp;
        data->poly.opacity[1] = opacn;
        data->poly.density = j;
        data->poly.refine = refine;
        data->mol.Psi = psi2;
        i = data->poly.wire;
        data->poly.wire = SendDlgItemMessage(hdlg,PolyWire,BM_GETCHECK,0,0)+
                        2*SendDlgItemMessage(hdlg,PolyWire2,BM_GETCHECK,0,0);
        if (i!=data->poly.wire) {
          if (data->poly.process>=6)    data->poly.process = 7; }
        unlock_window(data); }
      case IDCANCEL: EndDialog(hdlg, 1);  return(1);
      case PolyAsym: DialogBox(hinst, AsymDLG, hdlg, asymptote_dialog);
        return(0);
      case PolyAuto: if (!(data=lock_window(hwnd)))  break;
        if (!data->mol.maxpsi) {
          Busy += 10;
          psi = data->mol.Psi;
          for (i=0; i<data->mol.nump; i++)
            prep_atom(data->mol.orb+i, data->mol.orb[i].factor);
          data->mol.Psi = 0.01;
          prep_check(&data->mol);
          while (data->mol.maxpsi<1e-30) {
            data->mol.Psi /= 10;
            prep_check(&data->mol); }
          data->mol.Psi = psi;
          Busy -= 10; }
        psi = log10(data->mol.maxpsi)-1-0.25*(data->mol.orb[0].n-
              data->mol.orb[0].l)-0.4*(data->mol.nump-1);
        SetScrollPos(GetDlgItem(hdlg, PolyScrPsi), SB_CTL, (long)psi*-10, 1);
        sprintf(text, "%3.2f", psi);
        SetDlgItemText(hdlg, PolyPsi, text);
        for (i=j=0; i<data->mol.nump; i++) {
          if (data->mol.orb[i].n-data->mol.orb[i].l>j)
            j = data->mol.orb[i].n-data->mol.orb[i].l;
          if (data->mol.orb[i].l-fabs(data->mol.orb[i].m)>j)
            j = data->mol.orb[i].l-fabs(data->mol.orb[i].m);
          if (fabs(data->mol.orb[i].m)>j)
            j = fabs(data->mol.orb[i].m); }
        j = j*2*data->mol.nump;  if (j<6)  j = 6;
        SetScrollPos(GetDlgItem(hdlg, PolyScrD), SB_CTL, (long)j/2, 1);
        sprintf(text, "%d", j);  SetDlgItemText(hdlg, PolyDens, text);
        unlock_window(data);  return(0);
      case PolyDens: GetDlgItemText(hdlg, PolyDens, text, 79);
        i = GetScrollPos(GetDlgItem(hdlg, PolyScrD), SB_CTL)*2;
        sscanf(text, "%d", &i);
        SetScrollPos(GetDlgItem(hdlg, PolyScrD), SB_CTL, i/2, 1);
        return(0);
      case PolyEditN: GetDlgItemText(hdlg, PolyEditN, text, 79);
        i = GetScrollPos(GetDlgItem(hdlg, PolyScrN), SB_CTL);
        sscanf(text, "%d", &i);
        SetScrollPos(GetDlgItem(hdlg, PolyScrN), SB_CTL, i, 1);
        return(0);
      case PolyEditP: GetDlgItemText(hdlg, PolyEditP, text, 79);
        i = GetScrollPos(GetDlgItem(hdlg, PolyScrP), SB_CTL);
        sscanf(text, "%d", &i);
        SetScrollPos(GetDlgItem(hdlg, PolyScrP), SB_CTL, i, 1);
        return(0);
      case PolyPsi: GetDlgItemText(hdlg, PolyPsi, text, 79);
        psi = GetScrollPos(GetDlgItem(hdlg, PolyScrPsi), SB_CTL)*-0.1;
        sscanf(text, "%lf", &psi);
        SetScrollPos(GetDlgItem(hdlg, PolyScrPsi), SB_CTL, -psi*10+0.5, 1);
        return(0);
      case PolyRefine: GetDlgItemText(hdlg, PolyRefine, text, 79);
        refine = GetScrollPos(GetDlgItem(hdlg, PolyScrR), SB_CTL)*0.1;
        sscanf(text, "%lf", &refine);
        SetScrollPos(GetDlgItem(hdlg, PolyScrR), SB_CTL, refine*10, 1);
        return(0);
      case PolyWire: if (SendDlgItemMessage(hdlg,PolyWire,BM_GETCHECK,0,0))
          SendDlgItemMessage(hdlg,PolyWire2,BM_SETCHECK,0,0); break;
      case PolyWire2: if (SendDlgItemMessage(hdlg,PolyWire2,BM_GETCHECK,0,0))
          SendDlgItemMessage(hdlg,PolyWire,BM_SETCHECK,0,0); break;
      default: return(0); } break;
    case WM_HSCROLL: scr = GetDlgCtrlID(lp);
      GetScrollRange(lp, SB_CTL, &x, &y);
      i = GetScrollPos(lp, SB_CTL);
      switch (wp&0xFFFF) {
        case SB_BOTTOM: i = x; break;
        case SB_LINELEFT: if (i>x) i--; break;
        case SB_LINERIGHT: if (i<y) i++; break;
        case SB_PAGELEFT: i -=2*(1+4*(scr!=PolyScrD)); if (i<x) i = x; break;
        case SB_PAGERIGHT: i+=2*(1+4*(scr!=PolyScrD)); if (i>y) i = y; break;
        case SB_THUMBPOSITION: case SB_THUMBTRACK: i = (wp>>16); break;
        case SB_TOP: i = y; }
      SetScrollPos(lp, SB_CTL, i, 1);
      switch (scr) {
        case PolyScrD: sprintf(text, "%d", i*2);
          SetDlgItemText(hdlg, PolyDens, text); break;
        case PolyScrN: sprintf(text, "%d", i);
          SetDlgItemText(hdlg, PolyEditN, text); break;
        case PolyScrP: sprintf(text, "%d", i);
          SetDlgItemText(hdlg, PolyEditP, text); break;
        case PolyScrPsi: sprintf(text, "%3.2f", -i*0.1);
          SetDlgItemText(hdlg, PolyPsi, text); break;
        case PolyScrR: sprintf(text, "%3.1f", i*0.1);
          SetDlgItemText(hdlg, PolyRefine, text); } break;
    case WM_INITDIALOG: if (!(data=lock_window(hwnd)))
        return(1);
      if (data->poly.density<6)  data->poly.density = 6;
      SetScrollRange(GetDlgItem(hdlg, PolyScrPsi), SB_CTL, 1, 120, 0);
      psi = log10(data->mol.Psi);
      SetScrollPos(GetDlgItem(hdlg, PolyScrPsi), SB_CTL, (long)psi*-10, 1);
      sprintf(text, "%3.2f", psi);
      SetDlgItemText(hdlg, PolyPsi, text);
      SetScrollRange(GetDlgItem(hdlg, PolyScrP), SB_CTL, 0, 100, 0);
      SetScrollPos(GetDlgItem(hdlg, PolyScrP), SB_CTL,
                   (long)data->poly.opacity[0]*100, 1);
      sprintf(text, "%g", data->poly.opacity[0]*100);
      SetDlgItemText(hdlg, PolyEditP, text);
      SetScrollRange(GetDlgItem(hdlg, PolyScrN), SB_CTL, 0, 100, 0);
      SetScrollPos(GetDlgItem(hdlg, PolyScrN), SB_CTL,
                   (long)data->poly.opacity[1]*100, 1);
      sprintf(text, "%g", data->poly.opacity[1]*100);
      SetDlgItemText(hdlg, PolyEditN, text);
      SetScrollRange(GetDlgItem(hdlg, PolyScrD), SB_CTL, 3, 30, 0);
      SetScrollPos(GetDlgItem(hdlg, PolyScrD), SB_CTL,
                   (long)data->poly.density/2, 1);
      sprintf(text, "%d", data->poly.density);
      SetDlgItemText(hdlg, PolyDens, text);
      SetScrollRange(GetDlgItem(hdlg, PolyScrR), SB_CTL, 10, 40, 0);
      SetScrollPos(GetDlgItem(hdlg, PolyScrR), SB_CTL,
                   (long)data->poly.refine*10, 1);
      sprintf(text, "%3.1f", data->poly.refine);
      SetDlgItemText(hdlg, PolyRefine, text);
      SendDlgItemMessage(hdlg, PolyWire, BM_SETCHECK, data->poly.wire==1, 0);
      SendDlgItemMessage(hdlg, PolyWire2,BM_SETCHECK, data->poly.wire==2, 0);
      SetFocus(hdlg);
      unlock_window(data);
      return(1); }
  return(0);
}

process_arguments()
/* Process a set of command line arguments from another copy of the program.
 *                                                             3/14/01-DWM */
{
  char *temp=0, *mem;
  HANDLE gmem;

  if (!(gmem=CreateFileMapping(0xFFFFFFFF, 0, PAGE_READWRITE, 0,
        0, "PassArguments")))  return;
  mem = MapViewOfFile(gmem, FILE_MAP_WRITE, 0, 0, 0);
  if (mem) {
    if (strlen(mem))
      if (temp=malloc2(strlen(mem)+1))
        strcpy(temp, mem);
    UnmapViewOfFile(mem); }
  CloseHandle(gmem);
  if (temp) {
    strncpy(OpenName, temp, NAMELEN-1);
    OpenName[NAMELEN-1] = 0;
    open_window_mid();
    free2(temp); }
}

progress(HWND hwnd, char *title, char *message, long current, long maximum)
/* Show a progress dialog if a process takes a significant amount of time, or
 *  update the progress dialog, or dispose of the progress dialog.  This also
 *  checks if the user has canceled the dialog.  Note that the progress will
 *  never decline.  To reduce the value in the Progress meter, change the
 *  global value of ProgressPos, which has a range from 0 (0%) to 1000
 *  (100%), then call this routine.
 * Enter: HWND hwnd: owner of the dialog.  If this is non-zero, the dialog
 *                   will be created (assuming the process takes long
 *                   enough).
 *        char *title: title of window.  If this is non-zero, the title bar
 *                     of the progress dialog is changed.  If this is -1, the
 *                     dialog is changed to a topmost window.
 *        char *message: message for progress dialog.  If this is non-zero,
 *                       the message is changed.
 *        long current: current number of units processed (0 to max,
 *                      inclusive).  If this is -1, the old number is
 *                      maintained.
 *        long maximum: maximum number of units to process.  If this is zero,
 *                      the dialog is closed.  If this is negative, the old
 *                      number is maintained.  If this is less than -1, the
 *                      dialog is displayed even if it has not been
 *                      PROGRESSDELAY seconds since process was first invoked.
 * Exit:  long proceed: 0 for user has canceled, 1 for proceed.9/19/98-DWM */
{
  static HWND owner=0, hdlg=0;
  static long cur, max, topmost;
  static char tit[256], mes[1024];
  static clock_t time1, time0;
  clock_t time2;
  MSG msg;
  long i, h, m, s, out[4], cbox[4]={8,28,144,10}, first=0;
  float pos;
  char text[80];
  RECT rect;

  if (!maximum) {
    cursor(0);
    if (hdlg)  DestroyWindow(hdlg);  hdlg = 0;  owner = 0;
    return(1); }
  if (hwnd && !hdlg) {
    owner = hwnd;  cursor(1);
    time0 = time1 = clock(); }
  if (!owner)  return(1);
  cursor(1);
  if ((long)title==-1) {
    topmost = 1;
    if (hdlg)
      SetWindowPos(hdlg, HWND_TOPMOST, 0,0,0,0, SWP_NOACTIVATE|SWP_NOMOVE|
                   SWP_NOSIZE); }
  else if (title)   {
    strncpy(tit, title, 255);  tit[255] = 0;  topmost = 0;
    if (hdlg)  SetWindowText(hdlg, tit); }
  if (message) {
    strncpy(mes, message, 1023);  mes[1023] = 0;
    if (hdlg)  SetDlgItemText(hdlg, ProgMessage, mes); }
  if (maximum>0)   max = maximum;
  if (current>=0)  cur = current;
  if (cur>max)  cur = max;
  if (maximum<-1)  time0 = 0;
  time2 = clock();
  if (!hdlg) {
    if (time2-time0<CLOCKS_PER_SEC*PROGRESSDELAY)  return(1);
    if (max==cur)  return(1);
    ProgressCancel = 0;  ProgressPos = 0;  first = 1;
    hdlg = CreateDialog(hinst, ProgDLG, owner, progress_dialog);
    if (topmost)
      SetWindowPos(hdlg, HWND_TOPMOST, 0,0,0,0, SWP_NOACTIVATE|SWP_NOMOVE|
                   SWP_NOSIZE);
    SetWindowText(hdlg, tit);
    SetDlgItemText(hdlg, ProgMessage, mes); }
  for (i=0; i<1+5*(!first); i++)
    if (PeekMessage(&msg, hdlg, 0, 0, PM_NOREMOVE))
      if (GetMessage(&msg, hdlg, 0, 0)) {
        TranslateMessage(&msg);
        if (msg.message==WM_KEYDOWN && msg.wParam==VK_ESCAPE)
          ProgressCancel = 1;
        DispatchMessage(&msg); }
  if (ProgressCancel) {
    DestroyWindow(hdlg);  hdlg = 0;  owner = 0;  return(0); }
  if (cur==max && ProgressPos!=1000)  time1 = 0;
  if ((time2-time1)*3<CLOCKS_PER_SEC && !first)  return(1);
  pos = 1000.*cur/max;
  if (pos<ProgressPos)  pos = ProgressPos;
  if (((long)pos)!=ProgressPos) {
    map_dialog_rect(hdlg, cbox[0], cbox[1], cbox[2], cbox[3], out);
    rect.left = out[0]+1;  rect.right = out[0]+out[2];
    rect.top = out[1]+1;   rect.bottom = out[1]+out[3];
    ProgressPos = pos;
    InvalidateRect(hdlg, &rect, 1);
    UpdateWindow(hdlg); }
  if (pos) {
    s = ((float)time2-time0)/CLOCKS_PER_SEC / pos * (1000-pos) + 0.999;
    h = s / 3600;  s -= h*3600;
    m = s / 60;    s -= m*60;
    if (h)  sprintf(text, "%d:%02d", h, m);
    else    sprintf(text, "%d", m);
    sprintf(text+strlen(text), ":%02d estimated time remaining", s);
    SetDlgItemText(hdlg, ProgLeft, text); }
  time1 = time2;
  return(1);
}

BOOL CALLBACK progress_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the progress dialog box.
 * Enter: HWND hwnd: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.        9/19/98-DWM */
{
  char text[80];
  long out[4], cbox[4]={8,28,144,10}, w;
  RECT rect;
  PAINTSTRUCT paint;
  HDC hdc;
  HBRUSH old, new;

  switch (msg) {
    case WM_CLOSE: case WM_COMMAND: ProgressCancel = 1; break;
    case WM_INITDIALOG: SetFocus(hdlg); return(1);
    case WM_PAINT: hdc = BeginPaint(hdlg, &paint);
      map_dialog_rect(hdlg, cbox[0], cbox[1], cbox[2], cbox[3], out);
      SelectObject(hdc, GetStockObject(NULL_BRUSH));
      Rectangle(hdc, out[0], out[1], out[0]+out[2], out[1]+out[3]);
      sprintf(text, "%2d%%", ProgressPos/10);
      rect.left = 0;  rect.right = out[2]-2;
      rect.top = 0;   rect.bottom = out[3]-2;
      DrawText(hdc, text, strlen(text), &rect, DT_LEFT|DT_NOPREFIX|
               DT_SINGLELINE|DT_CALCRECT);
      w = rect.right;
      rect.left = out[0]+1+(out[2]-w)/2;  rect.right = out[0]+out[2]-1;
      rect.top = out[1]+1;   rect.bottom = out[1]+out[3]-1;
      SetTextColor(hdc, GetSysColor(COLOR_WINDOWTEXT));
      SetBkColor(hdc, GetSysColor(COLOR_BTNFACE));
      DrawText(hdc, text, strlen(text), &rect, DT_LEFT|DT_NOPREFIX|
               DT_SINGLELINE|DT_VCENTER);
      SelectObject(hdc, GetStockObject(NULL_PEN));
      new = CreateSolidBrush(GetSysColor(COLOR_ACTIVECAPTION));
      old = SelectObject(hdc, new);
      Rectangle(hdc, out[0]+1, out[1]+1, w=(out[0]+(out[2]*ProgressPos/1000)),
                out[1]+out[3]);
      SelectObject(hdc, old);
      DeleteObject(new);
      rect.right = w;
      SetTextColor(hdc, GetSysColor(COLOR_BTNFACE));
      SetBkColor(hdc, GetSysColor(COLOR_ACTIVECAPTION));
      DrawText(hdc, text, strlen(text), &rect, DT_LEFT|DT_NOPREFIX|
               DT_SINGLELINE|DT_VCENTER);
      EndPaint(hdlg, &paint);  return(1); }
  return(0);
}

quit()
/* Exit the program.                                            5/8/96-DWM */
{
  PostQuitMessage(0);
}

void *realloc2(void *current, long size)
/* Actually call realloc.  Provided for source code compatibility.
 * Enter: void *current: pointer to current data.
 *        long size: number of bytes to allocate for the new block.
 * Exit:  void *mem: pointer to memory.                         5/4/96-DWM */
{
  long i;
  HGLOBAL hmem, hmem2;

  for (i=0; i<mallocnum; i++)
    if ((long)current==malloctbl[i*2+1]) {
      hmem = (HGLOBAL)malloctbl[i*2];
      GlobalUnlock(hmem);
      hmem2 = GlobalReAlloc(hmem, size, GMEM_MOVEABLE);
      if (!hmem2) { GlobalLock(hmem);  return(0); }
      current = GlobalLock(hmem2);
      malloctbl[i*2] = (long)hmem2;
      malloctbl[i*2+1] = (long)current;
      return(current); }
  return(0);
}

BOOL CALLBACK render_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the Render settings dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.        5/26/96-DWM */
{
  HWND hwnd;
  DATA *data;
  long old;

  hwnd = SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0);
  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case HelpHelp: WinHelp(Hwnd,HelpFile,HELP_CONTEXT,HelpRender); break;
      case IDOK: if (data=lock_window(hwnd)) {
        old = data->dflag;
        data->dflag &= (~0x1F);
        if (SendDlgItemMessage(hdlg,RendQuick,BM_GETCHECK,0,0))
          data->dflag |= 1;
        if (SendDlgItemMessage(hdlg,RendA2,BM_GETCHECK,0,0)) data->dflag|=2;
        if (SendDlgItemMessage(hdlg,RendA3,BM_GETCHECK,0,0)) data->dflag|=4;
        if (SendDlgItemMessage(hdlg,RendB2,BM_GETCHECK,0,0)) data->dflag|=8;
        if (SendDlgItemMessage(hdlg,RendB3,BM_GETCHECK,0,0)) data->dflag|=16;
        if (data->dflag!=old)  data->changed = 1;
        if ((old^data->dflag)&0x18)  InvalidateRect(hwnd, 0, 1);
        unlock_window(data); }
      case IDCANCEL: EndDialog(hdlg, 1);  return(1);
      default: return(0); } break;
    case WM_INITDIALOG: if (!(data=lock_window(hwnd)))
        return(1);
      SendDlgItemMessage(hdlg, RendQuick, BM_SETCHECK, data->dflag&1, 0);
      SendDlgItemMessage(hdlg,RendA1,BM_SETCHECK,((data->dflag>>1)&3)==0,0);
      SendDlgItemMessage(hdlg,RendA2,BM_SETCHECK,((data->dflag>>1)&3)==1,0);
      SendDlgItemMessage(hdlg,RendA3,BM_SETCHECK,((data->dflag>>1)&3)==2,0);
      SendDlgItemMessage(hdlg,RendB1,BM_SETCHECK,((data->dflag>>3)&3)==0,0);
      SendDlgItemMessage(hdlg,RendB2,BM_SETCHECK,((data->dflag>>3)&3)==1,0);
      SendDlgItemMessage(hdlg,RendB3,BM_SETCHECK,((data->dflag>>3)&3)==2,0);
      SetFocus(hdlg);
      unlock_window(data);
      return(1); }
  return(0);
}

BOOL CALLBACK render_opt_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the Raytrace Options dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.       12/29/97-DWM */
{
  HWND hwnd;
  DATA *data;
  long x, y, i, j, scr, steps, anti, aut, mod=0;
  real psi, psi2, val, opac[7], refrac[3];
  char text[80], text2[80];
  static long link=1;
  long opactbl[]={RendEditP0,RendEditP1,RendEditP2,RendEditN0,RendEditN1,
                  RendEditN2,RendEditA1,RendEditPR,RendEditNR,RendEditAR};

  hwnd = SendMessage(HwndC, WM_MDIGETACTIVE, 0, 0);
  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case HelpHelp: WinHelp(Hwnd, HelpFile, HELP_CONTEXT, HelpRenderOpt);
        break;
      case IDOK: if (data=lock_window(hwnd)) {
        GetDlgItemText(hdlg, RendPsi, text, 79);
        psi = GetScrollPos(GetDlgItem(hdlg, RendScrPsi), SB_CTL)*-0.1;
        sscanf(text, "%lf", &psi);  psi2 = pow(10, psi);
        for (i=0; i<7; i++) {
          GetDlgItemText(hdlg, opactbl[i], text, 79);
          opac[i] = GetScrollPos(GetDlgItem(hdlg, opactbl[i]+1), SB_CTL);
          sscanf(text, "%lf", opac+i);  opac[i] *= 0.01; }
        for (i=0; i<3; i++) {
          GetDlgItemText(hdlg, opactbl[i+7], text, 79);
          refrac[i] = GetScrollPos(GetDlgItem(hdlg, opactbl[i+7]+1),
                                   SB_CTL)*0.01;
          sscanf(text, "%lf", refrac+i); }
        GetDlgItemText(hdlg, RendStep, text, 79);
        steps = GetScrollPos(GetDlgItem(hdlg, RendScrStep), SB_CTL)*10;
        sscanf(text, "%d", &steps);
        link = SendDlgItemMessage(hdlg, RendLink, BM_GETCHECK, 0, 0);
        aut = SendDlgItemMessage(hdlg, RendAutoB, BM_GETCHECK, 0, 0);
        if (aut && data->render.autobright)  aut = data->render.autobright;
        anti = SendDlgItemMessage(hdlg, RendCoarse, BM_GETCHECK, 0, 0)+
               2*SendDlgItemMessage(hdlg, RendAnti, BM_GETCHECK, 0, 0);
        for (i=0; i<7; i++) {
          if (opac[i]<0)  opac[i] = 0;  if (opac[i]>1)  opac[i] = 1; }
        for (i=0; i<3; i++)  if (refrac[i]<=0)  refrac[i] = 1;
        if (steps<10)  steps = 10;
        opac[2] *= 0.01;  opac[5] *= 0.01;
        if (steps!=data->render.steps)  mod |= 0xF;
        if (memcmp(opac, data->render.opacity, 7*sizeof(real)))  mod |= 1;
        if (memcmp(refrac, data->render.refrac, 3*sizeof(real)))  mod |= 2;
        if (psi2!=data->mol.Psi)  mod |= 4;
        if (aut!=data->render.autobright)  mod |= 8;
        if (anti!=(data->render.antialias&3))  mod |= 16;
        if (mod) {
          data->changed = 1;
          if (mod&0xF)
            orb_render(&data->render, 0, 0, &data->mol, 0, 4);
          if (((anti^data->render.antialias)&1) && !(anti&1) &&
              data->render.process>=5) {
            data->render.x = data->render.y = 0;
            data->render.process = 5; }
          if (((anti^data->render.antialias)&2) && (anti&2) &&
              data->render.process>=6) {
            data->render.x = data->render.y = 0;
            data->render.process = 6; }
            }
        if (opac[6]!=data->asym.opacity) {
          data->asym.opacity = opac[6];
          if (data->points.process>=5)  data->points.process = 6;
          if (data->asym.process>=5)    data->asym.process = 6;
          if (data->poly.process>=6)    data->poly.process = 7; }
        if (mod&4)
          if (data->poly.process>1)  data->poly.process = 1;
        memcpy(data->render.opacity, opac, 7*sizeof(real));
        memcpy(data->render.refrac, refrac, 3*sizeof(real));
        data->render.steps = steps;
        data->render.autobright = aut;
        if (!aut)  data->render.brightness = 1;
        data->render.antialias = (data->render.antialias&0xFFFFFFFC)|anti;
        data->mol.Psi = psi2;
        unlock_window(data); }
      case IDCANCEL: EndDialog(hdlg, 1);  return(1);
      case RendAuto: if (!(data=lock_window(hwnd)))  break;
        if (!data->mol.maxpsi) {
          Busy += 10;
          psi = data->mol.Psi;
          for (i=0; i<data->mol.nump; i++)
            prep_atom(data->mol.orb+i, data->mol.orb[i].factor);
          data->mol.Psi = 0.01;
          prep_check(&data->mol);
          while (data->mol.maxpsi<1e-30) {
            data->mol.Psi /= 10;
            prep_check(&data->mol); }
          data->mol.Psi = psi;
          Busy -= 10; }
        psi = log10(data->mol.maxpsi)-1-0.25*(data->mol.orb[0].n-
              data->mol.orb[0].l)-0.4*(data->mol.nump-1);
        SetScrollPos(GetDlgItem(hdlg, RendScrPsi), SB_CTL, (long)psi*-10, 1);
        sprintf(text, "%3.2f", psi);
        SetDlgItemText(hdlg, RendPsi, text);
        unlock_window(data);  return(0);
      case RendEditA1: case RendEditN0: case RendEditN1: case RendEditN2:
      case RendEditP0: case RendEditP1: case RendEditP2:
        GetDlgItemText(hdlg, wp&0xFFFF, text, 79);
        i = GetScrollPos(GetDlgItem(hdlg, (wp&0xFFFF)+1), SB_CTL);
        sscanf(text, "%d", &i);
        SetScrollPos(GetDlgItem(hdlg, (wp&0xFFFF)+1), SB_CTL, i, 1);
        if (SendDlgItemMessage(hdlg,RendLink,BM_GETCHECK,0,0) &&
            (wp&0xFFFF)<=RendScrNR) {
          GetDlgItemText(hdlg, (wp&0xFFFF)^2, text2, 79);
          if (strcmp(text, text2))
            SetDlgItemText(hdlg, (wp&0xFFFF)^2, text); }
        return(0);
      case RendEditAR: case RendEditNR: case RendEditPR:
        GetDlgItemText(hdlg, wp&0xFFFF, text, 79);
        val = GetScrollPos(GetDlgItem(hdlg, (wp&0xFFFF)+1), SB_CTL)*0.01;
        sscanf(text, "%lf", &val);
        SetScrollPos(GetDlgItem(hdlg, (wp&0xFFFF)+1), SB_CTL, val*100+0.1,1);
        if (SendDlgItemMessage(hdlg,RendLink,BM_GETCHECK,0,0) &&
            (wp&0xFFFF)<=RendScrNR) {
          GetDlgItemText(hdlg, (wp&0xFFFF)^2, text2, 79);
          if (strcmp(text, text2))
            SetDlgItemText(hdlg, (wp&0xFFFF)^2, text); }
        return(0);
      case RendPsi: GetDlgItemText(hdlg, RendPsi, text, 79);
        psi = GetScrollPos(GetDlgItem(hdlg, RendScrPsi), SB_CTL)*-0.1;
        sscanf(text, "%lf", &psi);
        SetScrollPos(GetDlgItem(hdlg, RendScrPsi), SB_CTL, -psi*10+0.5, 1);
        return(0);
      case RendStep: GetDlgItemText(hdlg, RendStep, text, 79);
        i = GetScrollPos(GetDlgItem(hdlg, RendScrStep), SB_CTL)*10;
        sscanf(text, "%d", &i);
        SetScrollPos(GetDlgItem(hdlg, RendScrStep), SB_CTL, i*0.1, 1);
        return(0);
      default: return(0); } break;
    case WM_HSCROLL: scr = GetDlgCtrlID(lp);
      GetScrollRange(lp, SB_CTL, &x, &y);
      i = GetScrollPos(lp, SB_CTL);
      switch (wp&0xFFFF) {
        case SB_BOTTOM: i = x; break;
        case SB_LINELEFT: if (i>x) i--; break;
        case SB_LINERIGHT: if (i<y) i++; break;
        case SB_PAGELEFT: i -= 10; if (i<x) i = x; break;
        case SB_PAGERIGHT: i+= 10; if (i>y) i = y; break;
        case SB_THUMBPOSITION: case SB_THUMBTRACK: i = (wp>>16); break;
        case SB_TOP: i = y; }
      SetScrollPos(lp, SB_CTL, i, 1);
      switch (scr) {
        case RendScrAR: case RendScrNR: case RendScrPR:
          sprintf(text, "%3.2f", i*0.01);
          SetDlgItemText(hdlg, scr-1, text); break;
        case RendScrPsi: sprintf(text, "%3.2f", -i*0.1);
          SetDlgItemText(hdlg, RendPsi, text); break;
        case RendScrStep: sprintf(text, "%d", i*10);
          SetDlgItemText(hdlg, RendStep, text); break;
        default: sprintf(text, "%d", i);
          SetDlgItemText(hdlg, scr-1, text); } break;
    case WM_INITDIALOG: if (!(data=lock_window(hwnd)))
        return(1);
      SetScrollRange(GetDlgItem(hdlg, RendScrPsi), SB_CTL, 1, 120, 0);
      psi = log10(data->mol.Psi);
      SetScrollPos(GetDlgItem(hdlg, RendScrPsi), SB_CTL, (long)psi*-10, 1);
      sprintf(text, "%3.2f", psi);
      SetDlgItemText(hdlg, RendPsi, text);
      for (i=0; i<7; i++) {
        SetScrollRange(GetDlgItem(hdlg, opactbl[i]+1), SB_CTL, 0, 100, 0);
        SetScrollPos(GetDlgItem(hdlg, opactbl[i]+1), SB_CTL,
                   (long)data->render.opacity[i]*100*(1+99*(i==2||i==5)), 1);
        sprintf(text, "%g", data->render.opacity[i]*100*(1+99*(i==2||i==5)));
        SetDlgItemText(hdlg, opactbl[i], text); }
      for (i=0; i<3; i++) {
        SetScrollRange(GetDlgItem(hdlg, opactbl[i+7]+1), SB_CTL, 90, 150, 0);
        SetScrollPos(GetDlgItem(hdlg, opactbl[i+7]+1), SB_CTL,
                     (long)data->render.refrac[i]*100, 1);
        sprintf(text, "%3.2f", data->render.refrac[i]);
        SetDlgItemText(hdlg, opactbl[i+7], text); }
      SetScrollRange(GetDlgItem(hdlg, RendScrStep), SB_CTL, 10, 200, 0);
      SetScrollPos(GetDlgItem(hdlg, RendScrStep), SB_CTL,
                   data->render.steps*0.01, 1);
      sprintf(text, "%d", data->render.steps);
      SetDlgItemText(hdlg, RendStep, text);
      SendDlgItemMessage(hdlg, RendLink, BM_SETCHECK, link, 0);
      SendDlgItemMessage(hdlg, RendAutoB, BM_SETCHECK,
                         data->render.autobright, 0);
      SendDlgItemMessage(hdlg, RendAnti, BM_SETCHECK,
                         data->render.antialias&2, 0);
      SendDlgItemMessage(hdlg, RendCoarse, BM_SETCHECK,
                         data->render.antialias&1, 0);
      SetFocus(hdlg);
      unlock_window(data);
      return(1); }
  return(0);
}

LRESULT CALLBACK second_loop(HWND hwnd, ulong msg, WPARAM wp, LPARAM lp)
/* Secondary processing loop.
 * Enter: HWND hwnd: handle of current window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.         5/6/96-DWM */
{
  long i;

  switch (msg) {
    case WM_CLOSE: if (close_window(hwnd)) return(0);
      for (i=0; i<NumWindows; i++)
        if (hwnd==HwndList[i])
          break;
      if (i<NumWindows)
        memmove(HwndList+i, HwndList+i+1, (NumWindows-i-1)*sizeof(HWND));
      NumWindows--;
      if (!NumWindows)
        recheck(0, 0); break;
    case WM_CREATE: NumWindows++;  HwndList[NumWindows-1] = hwnd;  break;
    case WM_LBUTTONDOWN: mouse(hwnd, 0, 0, wp, lp&0xFFFF, lp>>16); break;
    case WM_LBUTTONUP: mouse(hwnd, 0, 1, wp, lp&0xFFFF, lp>>16); break;
    case WM_MBUTTONDOWN: mouse(hwnd, 2, 0, wp, lp&0xFFFF, lp>>16); break;
    case WM_MBUTTONUP: mouse(hwnd, 2, 1, wp, lp&0xFFFF, lp>>16); break;
    case WM_MOUSEMOVE: mouse(hwnd,-1, 0, wp, lp&0xFFFF, lp>>16); break;
    case WM_PAINT: if (IsIconic(hwnd)) break;
      if (update_window(hwnd))  return(0);
      break;
    case WM_RBUTTONDOWN: mouse(hwnd, 1, 0, wp, lp&0xFFFF, lp>>16); break;
    case WM_RBUTTONUP: mouse(hwnd, 1, 1, wp, lp&0xFFFF, lp>>16); break;
    case WM_SIZE: InvalidateRect(hwnd, 0, 1); break;
    default: ; }
  return(DefMDIChildProc(hwnd, msg, wp, lp));
}

BOOL CALLBACK sequence_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the Play Sequence dialog box.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.        1/11/98-DWM */
{
  DATA *data=DispData;
  static DATA *seq[4];
  long i, j;
  char text[80];
  char *form[]={"PPM","TIF","BMP","AVI"};

  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case HelpHelp: WinHelp(Hwnd,HelpFile,HELP_CONTEXT,HelpSequence); break;
      case IDCANCEL: for (i=0; i<4; i++) if (seq[i]) {
          free2(seq[i]->mol.orb);  free2(seq[i]->mol.ls);
          free2(seq[i]); }
        EndDialog(hdlg, 1);  return(1);
      case SeqBrowse: sequence_browse(hdlg); break;
      case SeqCheck1: case SeqCheck2: case SeqCheck3: case SeqCheck4:
        if (SendDlgItemMessage(hdlg,wp&0xFFFF,BM_GETCHECK,0,0)) {
          i = ((wp&0xFFFF)-SeqCheck1)/3;
          seq[i] = sequence_file(hdlg, seq[i]);
          SetDlgItemText(hdlg, SeqName1+i*3, seq[i]->name); } break;
      case SeqPlay: for (i=j=0; i<4; i++)
          if (seq[i]&&SendDlgItemMessage(hdlg,SeqCheck1+i*3,BM_GETCHECK,0,0))
            j++;
        if (j<2) {
          error(hdlg, "At least two orbitals must be specified to play a sequence.");
          break; }
      case IDOK:
        for (i=j=0; i<4; i++) {
          if (data->seq[i]) {
            free2(data->seq[i]->mol.orb);  free2(data->seq[i]->mol.ls);
            free2(data->seq[i]);  data->seq[i] = 0; }
          if (!SendDlgItemMessage(hdlg,SeqCheck1+i*3,BM_GETCHECK,0,0) &&
              seq[i]) {
            free2(seq[i]->mol.orb);  free2(seq[i]->mol.ls);
            free2(seq[i]);  seq[i] = 0; }
          if (seq[i]) {
            seq[j]->frame = 0;
            GetDlgItemText(hdlg, SeqFrame1+i*3, text, 79);
            sscanf(text, "%d", &seq[j]->frame);
            data->seq[j] = seq[i];  j++; } }
        data->incr   = SendDlgItemMessage(hdlg,SeqIncr,  BM_GETCHECK,0,0);
        data->bezier = SendDlgItemMessage(hdlg,SeqBezier,BM_GETCHECK,0,0);
        GetDlgItemText(hdlg, SeqBase, data->basename, ONENAMELEN-1);
        if (strchr(data->basename, '.'))  strchr(data->basename, '.')[0] = 0;
        data->seqtype = SendDlgItemMessage(hdlg,SeqType,CB_GETCURSEL,0,0);
        GetDlgItemText(hdlg, SeqFPS, text, 79);
        sscanf(text, "%d", &data->seqfps);
        GetDlgItemText(hdlg, SeqStart, text, 79);
        sscanf(text, "%d", &data->frame);
        GetDlgItemText(hdlg, SeqEnd, text, 79);
        sscanf(text, "%d", &data->lastframe);
        EndDialog(hdlg, 1);
        if ((wp&0xFFFF)==SeqPlay)
          play_frame(data, data->frame);
        return(1);
      default: return(0); } break;
    case WM_INITDIALOG: seq[0] = seq[1] = seq[2] = seq[3] = 0;
      for (i=0; i<4; i++)
        if (data->seq[i])
          if (seq[i]=malloc2(sizeof(DATA))) {
            memcpy(seq[i], data->seq[i], sizeof(DATA));
            seq[i]->mol.orb = 0;  seq[i]->mol.ls = 0;
            if (seq[i]->mol.nump) {
              if (!(seq[i]->mol.orb=malloc2(seq[i]->mol.nump*
                  sizeof(OATOM)))) {
                free2(seq[i]);  seq[i] = 0;  continue; }
              memcpy(seq[i]->mol.orb, data->seq[i]->mol.orb,
                     seq[i]->mol.nump*sizeof(OATOM)); }
            if (seq[i]->mol.numl) {
              if (!(seq[i]->mol.ls=malloc2(seq[i]->mol.numl*
                  sizeof(LIGHT)))) {
                free2(seq[i]->mol.orb);
                free2(seq[i]);  seq[i] = 0;  continue; }
              memcpy(seq[i]->mol.ls, data->seq[i]->mol.ls,
                     seq[i]->mol.numl*sizeof(LIGHT)); } }
      for (i=0; i<4; i++)  if (seq[i]) {
        SendDlgItemMessage(hdlg, SeqCheck1+i*3, BM_SETCHECK, 1, 0);
        sprintf(text, "%d", seq[i]->frame);
        SetDlgItemText(hdlg, SeqFrame1+i*3, text);
        SetDlgItemText(hdlg, SeqName1+i*3, seq[i]->name); }
      SendDlgItemMessage(hdlg, SeqIncr,   BM_SETCHECK, data->incr, 0);
      SendDlgItemMessage(hdlg, SeqBezier, BM_SETCHECK, data->bezier, 0);
      SetDlgItemText(hdlg, SeqBase, data->basename);
      for (i=0; i<4; i++)
        SendDlgItemMessage(hdlg, SeqType, CB_ADDSTRING, 0, (long)form[i]);
      SendDlgItemMessage(hdlg, SeqType, CB_SETCURSEL, data->seqtype, 0);
      if (!data->seqfps)  data->seqfps = 30;
      sprintf(text, "%d", data->seqfps);
      SetDlgItemText(hdlg, SeqFPS, text);
      sprintf(text, "%d", data->frame);
      SetDlgItemText(hdlg, SeqStart, text);
      sprintf(text, "%d", data->lastframe);
      SetDlgItemText(hdlg, SeqEnd, text);
      SetFocus(hdlg);
      return(1); }
  return(0);
}

sequence_stop()
/* Stop a currently playing sequence.                          1/11/98-DWM */
{
  HWND hwnd;
  DATA *data;

  if (!(hwnd=SendMessage(HwndC,WM_MDIGETACTIVE,0,0)))  return(0);
  if (!(data=lock_window(hwnd)))  return(0);
  data->dflag &= (~(1<<10));
  unlock_window(data);
  recheck(hwnd, 0);
}

set_window_name(HWND hwnd, DATA *data)
/* Reset the current name of the window to match the name of the actual file.
 * Enter: HWND hwnd: window handle.
 *        DATA *data: pointer to window.                      10/24/97-DWM */
{
  char text[NAMELEN];

  strcpy(text, data->name);
  if (strchr(text, ':'))   strcpy(text, strchr(data->name, ':')+1);
  if (strchr(text, '\\'))  strcpy(text, strrchr(data->name, '\\')+1);
  SetWindowText(hwnd, text);
}

long size2(void *mem)
/* Return the size of the specified memory block.  If the pointer does not
 *  point to the base address of a memory block, a value of zero is returned.
 * Enter: void *mem: pointer to a locked memory block.
 * Exit:  long size: size of the memory block in bytes.        2/27/01-DWM */
{
  HANDLE hmem;

  if (!mem)  return(0);
  if (!(hmem=GlobalHandle(mem)))  return(0);
  return(GlobalSize(hmem));
}

start_avi()
/* Start up the AVI control functions, if available.           1/22/98-DWM */
{
  FARPROC proc;
  char name[]="AVIFIL32.DLL";

  AVILoad = 1;
  if (havi)  return(0);
  if (!windows_file(name))  return(0);
  if (!(havi=LoadLibrary(name)))
    return(0);
}

start_common()
/* Start up the common controls.  If unavailable, print a snide message, then
 *  bail.                                                       2/5/97-DWM */
{
  FARPROC init;
  char name[]="COMCTL32.DLL";

  if (!windows_file(name)) {
    MsgBox(Hwnd, "Program aborted.\n\n"
           "Can't use common controls.  The file COMCTL32.DLL\n"
           "is missing from the system directory.", "Terminal Error", MB_OK);
    return(0); }
  if (!(hcom=LoadLibrary(name)))
        return(0);
  init = GetProcAddress(hcom, "InitCommonControls");
  (*init)();
}

start_ctl3d()
/* Start up the faux-3D control functions, if available.       10/2/96-DWM */
{
  FARPROC reg, sub;
  char name[]="CTL3D32.DLL";

  if (!windows_file(name))  return(0);
  if (!(hctl3d=LoadLibrary(name)))
    return(0);
  reg = GetProcAddress(hctl3d, "Ctl3dRegister");
  sub = GetProcAddress(hctl3d, "Ctl3dAutoSubclass");
  if ((*reg)(hinst))
    (*sub)(hinst);
}

BOOL CALLBACK stereo_dialog(HWND hdlg, ulong msg, WPARAM wp, LPARAM lp)
/* Handle the controls in the Stereo Options dialog box.  This is primarily
 *  3D effects.
 * Enter: HWND hdlg: handle of current dialog window.
 *        long msg: message to process.
 *        WPARAM wp, LPARAM lp: parameters for message.       12/29/97-DWM */
{
  HWND ctrl;
  DATA *data=DispData;
  char text[80];
  long i, grey, mode=-1, len;
  real val;
  static long du, lastmode, update, lastswap, replace;
  static real ocular[2];
  static STEREO stereo;
  float dx0[4]={0,0,0,0};

  switch (msg) {
    case WM_COMMAND: switch (wp&0xFFFF) {
      case SterImage: if (SendDlgItemMessage(hdlg,SterImage,BM_GETCHECK,0,0))
          stereo_image(hdlg, &stereo, &replace);
        mode = lastmode;  break;
      case SterMode1: case SterMode2: case SterMode3: case SterMode4:
      case SterMode5: case SterMode6: case SterMode7:
        if ((wp&0xFFFF)-SterMode1==lastmode)  break;
        mode = (wp&0xFFFF)-SterMode1;
        GetDlgItemText(hdlg, SterOcular, text, 79);
        if (mode==StereoSTEREOSCOPE)          sscanf(text, "%lg", ocular);
        else if (lastmode==StereoSTEREOSCOPE) sscanf(text, "%lg", ocular+1);
        break;
      case SterSwap3: case SterSwap4: mode = lastmode; break;
      case SterUnits: GetDlgItemText(hdlg, SterSep, text, 79);
        val = 1e30;  sscanf(text, "%lf", &val);  val *= DistVal[du];
        du = SendDlgItemMessage(hdlg, SterUnits, CB_GETCURSEL, 0, 0);
        Pref.munit = du;
        if (val<1e10) {
          sprintf(text, "%g", val/DistVal[du]);
          SetDlgItemText(hdlg, SterSep, text); } break;
      case HelpHelp: WinHelp(Hwnd, HelpFile, HELP_CONTEXT, HelpStereo);
        break;
      case IDOK: KillTimer(hdlg, 6);
        for (i=0; i<StereoMODES; i++)
          if (SendDlgItemMessage(hdlg,SterMode1+i,BM_GETCHECK,0,0))
            stereo.mode = i;
        stereo.flags &= 0xFFFFFFF0;
        if (SendDlgItemMessage(hdlg,SterSwap3,BM_GETCHECK,0,0))
          stereo.flags |= 1;
        if (SendDlgItemMessage(hdlg,SterSwap4,BM_GETCHECK,0,0))
          stereo.flags |= 2;
        if (ocular[0]<1e20 && ocular[0]>0)
          stereo.interocular[0] = ocular[0];
        if (ocular[1]<1e20 && ocular[1]>0)
          stereo.interocular[1] = ocular[1];
        GetDlgItemText(hdlg, SterOcular, text, 79);
        val = 1e30;  sscanf(text, "%lf", &val);
        if (val<1e20 && val>0)
          stereo.interocular[stereo.mode==StereoSTEREOSCOPE] = val;
        GetDlgItemText(hdlg, SterSep, text, 79);
        val = 1e30;  sscanf(text, "%lf", &val);  val *= DistVal[du];
        if (val<1e20 && val>0)
          stereo.separation = val;
        if (SendDlgItemMessage(hdlg,SterSepAuto,BM_GETCHECK,0,0))
          stereo.flags |= 4;
        if (SendDlgItemMessage(hdlg,SterImage,BM_GETCHECK,0,0))
          stereo.flags |= 8;
        GetDlgItemText(hdlg, SterPerspect, text, 79);
        val = 1e30;  sscanf(text, "%lf", &val);
        if (val<0 || val>1e10)  val = data->pref.perspec;
        if (!replace)  stereo.image = data->stereo.image;
        if (memcmp(&data->stereo, &stereo, sizeof(STEREO)) ||
            val!=data->pref.perspec) {
          if (replace)  free2(data->stereo.image);
          memcpy(&data->stereo, &stereo, sizeof(STEREO));
          if (val!=data->pref.perspec) {
            dx0[2] = val-data->renderval[1];
            render_move(0, 0, dx0, 0, 0, data->renderval, data->renderdlt,
                        data->renderdlt); }
          data->pref.perspec = val;
          data->changed = 1;
          if (data->points.process>=5)  data->points.process = 6;
          if (data->asym.process>=5)    data->asym.process = 6;
          if (data->poly.process>=6)    data->poly.process = 7;
          orb_render(&data->render, 0, 0, &data->mol, 0, 4);
          update_process(data); }
        EndDialog(hdlg, 1);  return(1);
      case IDCANCEL: KillTimer(hdlg, 6);
        if (replace)  free2(stereo.image);
        EndDialog(hdlg, 1);  return(1);
      default: return(0); } break;
    case WM_INITDIALOG: memcpy(&stereo, &data->stereo, sizeof(STEREO));
      if (PreviewDelay<0 || PreviewDelay>3000)  PreviewDelay = 1000;
      KillTimer(hdlg, 6);
      replace = 0;
      mode = stereo.mode;
      SendDlgItemMessage(hdlg, SterMode1+mode, BM_SETCHECK, 1,0);
      SendDlgItemMessage(hdlg, SterSwap3, BM_SETCHECK, stereo.flags&1, 0);
      SendDlgItemMessage(hdlg, SterSwap4, BM_SETCHECK, stereo.flags&2, 0);
      ocular[0] = stereo.interocular[0];
      ocular[1] = stereo.interocular[1];
      sprintf(text, "%g", stereo.separation/DistVal[du=Pref.munit]);
      SetDlgItemText(hdlg, SterSep, text);
      for (i=0; i<4; i++)
        SendDlgItemMessage(hdlg, SterUnits, CB_ADDSTRING, 0,
                           (long)DistText[i]);
      SendDlgItemMessage(hdlg, SterUnits, CB_SETCURSEL, Pref.munit, 0);
      SendDlgItemMessage(hdlg, SterSepAuto, BM_SETCHECK, stereo.flags&4, 0);
      sprintf(text, "%g", data->pref.perspec);
      SetDlgItemText(hdlg, SterPerspect, text);
      if (stereo.name[0]) {
        SendDlgItemMessage(hdlg, SterImage, BM_SETCHECK, stereo.flags&8, 0);
        SetDlgItemText(hdlg, SterImageName, stereo.name); }
      else
        SetDlgItemText(hdlg, SterImageName, "None");
      SetFocus(hdlg); break;
    case WM_PAINT: KillTimer(hdlg, 6);
      SetTimer(hdlg, 6, 0, 0);  update = 1;  return(0);
    case WM_TIMER: KillTimer(hdlg, 6);
      if (lastmode<0 || !update)  return(0);
      update = i = 0;
      if (lastmode==2)  i=SendDlgItemMessage(hdlg,SterSwap3,BM_GETCHECK,0,0);
      if (lastmode==3)  i=SendDlgItemMessage(hdlg,SterSwap4,BM_GETCHECK,0,0);
      stereo.flags = (stereo.flags&0xFFFFFFF7)|(8*
                     SendDlgItemMessage(hdlg,SterImage,BM_GETCHECK,0,0));
      stereo_figure(hdlg, lastmode, i, &stereo);
      return(0); }
  if (mode>=0) {
    sprintf(text, "%g", ocular[mode==StereoSTEREOSCOPE]);
    SetDlgItemText(hdlg, SterOcular, text);
    grey = (mode==StereoSTEREOSCOPE || mode==StereoINTERLACED ||
            mode==StereoREDBLUE || mode==StereoSTEREOGRAM ||
            mode==StereoOVERLAY);
    for (i=SterOcularText; i<=SterOcularPix; i++) {
      ctrl = GetDlgItem(hdlg, i);  EnableWindow(ctrl, grey); }
    grey = (mode==StereoSTEREOSCOPE || mode==StereoINTERLACED ||
            mode==StereoREDBLUE || mode==StereoOVERLAY);
    for (i=SterSepText; i<=SterSepAuto; i++) {
      ctrl = GetDlgItem(hdlg, i);  EnableWindow(ctrl, grey); }
    grey = (mode==StereoSTEREOGRAM);
    for (i=SterImage; i<=SterImageName; i++) {
      ctrl = GetDlgItem(hdlg, i);  EnableWindow(ctrl, grey);
    lastmode = mode; } }
  if (mode!=-1) {
    KillTimer(hdlg, 6);  SetTimer(hdlg, 6, PreviewDelay, 0);  update = 1; }
  return(0);
}

test(char *data)
/* Standard test message entry point.  This routine allows diagnostic
 *  messages to be printed by non-Windows source routines.
 * Enter: char *data: string to print in a dialog.             10/1/96-DWM */
{
  FILE *fptr;

/**/  fptr = fopen("C:\\TEMP\\TEST.TXT", "at");
  fprintf(fptr, "%s\n", data);
  fclose(fptr);
  return(0); /**/
/**  SetWindowText(Hwnd, data); /**/
 **  debug(data);  /**/
}

VOID CALLBACK timer(HWND hwnd, ulong msg, ulong id, long time)
/* Enter: HWND hwnd: handle of calling window.
 *        ulong msg: ignored.
 *        ulong id: id of calling timer.
 *        long time: time in some windows format.               7/6/97-DWM */
{
  static long inloop=0, wnd=0;
  long quick, mode, update=0, i, status=0, lastm;
  char stext[80];
  float speed;
  real e, s, y;
  HDC dc;
  DATA *data;

  if (Busy<0)  Busy = 0;
  if (inloop || Busy>=10 || !NumWindows)  return;
  if ((wnd%3)==2) { wnd++;  return; }
  hwnd = HwndList[(wnd/3)%NumWindows];
  if (!(wnd%3))  hwnd = SendMessage(HwndC, WM_MDIGETACTIVE,0,0);
  wnd++;  if (wnd>2000)  wnd = 0;
  if (!(data=lock_window(hwnd)))  return;
  inloop = 1;
  if (Busy)  LastMove = 0;
  LastMove ++;
  if (LastMove<20)       speed = 0.1;
  else if (LastMove<40)  speed = 1;
  else                   speed = 5;
  quick = (data->dflag&1);
  lastm = (data->dflag>>5)&3;
  for (i=0; i<2; i++) {
    if ((!i && !quick) || (i && update))  continue;
    mode = (data->dflag>>(1+2*i))&3;
    switch (mode) {
      case 0: if (data->points.process!=5 || (lastm && i)) {
          update = 1;  data->dflag = (data->dflag&0xFFFFFFF9F);
          if (orb_points(&data->points,&data->cut,&data->mol,speed,0)==1)
            update = 2;
          status = 1; }
        if (update!=1 && data->asym.process!=5 && data->asym.opacity) {
          update = 1;  status = 2;  data->dflag = (data->dflag&0xFFFFFFF9F);
          orb_asymptote(&data->asym, &data->cut, &data->mol, speed, 0); }
        break;
      case 1: if (data->poly.process!=6 || (lastm!=1 && i)) {
          update = 1;  status = 3;
          data->dflag = (data->dflag&0xFFFFFFF9F)|(1<<5);
          if (orb_polygons(&data->poly,&data->cut,&data->mol,speed,0)==1)
            update = 2; }
        if (update!=1 && data->asym.process!=5 && data->asym.opacity) {
          update = 1;  status = 2;
          data->dflag = (data->dflag&0xFFFFFFF9F)|(1<<5);
          orb_asymptote(&data->asym, &data->cut, &data->mol, speed, 0); }
        break;
      case 2: if (data->render.process!=7 || (lastm!=2 && i)) {
        update_geo(hwnd, data, 0);
        update = 1;  status = 4;
        data->dflag = (data->dflag&0xFFFFFFF9F)|(2<<5);
        if (data->stereo.mode && data->render.process<=1) {
          orb_render(&data->render, &data->stereo, &data->cut, &data->mol,
                     0, 6);
          prep_stereo(data); }
        orb_render(&data->render, &data->stereo, &data->cut, &data->mol,
                   speed, 0); } } }
  if (update) {
    dc = GetDC(hwnd);
    if (Hpal) {
      SelectPalette(dc, Hpal, 0);
      RealizePalette(dc); }
    update_geo(hwnd, data, dc);
    ReleaseDC(hwnd, dc); }
  else if ((data->dflag>>10)&1)
    next_frame(data);
  if (hwnd==SendMessage(HwndC, WM_MDIGETACTIVE,0,0) && HwndStat) {
    stext[0] = 0;
    switch (status) {
      case 1: if (data->points.process<4) {
          strcpy(stext, "Points: initializing");  break; }
        if (data->points.maxn!=data->points.n && data->points.maxn)
          sprintf(stext, "Points: %3.1f%% done",
                  100.*data->points.n/data->points.maxn); break;
      case 2: if (data->asym.process<4) {
          strcpy(stext, "Asymptote: initializing");  break; }
        if (data->asym.density)
          sprintf(stext, "Asymptote: %3.1f%% done",
                  100.*((float)data->asym.index[2]*sq(data->asym.density)+
                  (float)data->asym.index[1]*data->asym.density+
                  (float)data->asym.index[0])/(data->asym.density*
                  data->asym.density*data->asym.density*2)); break;
      case 3: if (data->poly.process<4) {
          strcpy(stext, "Polygons: initializing");  break; }
        if (data->poly.process==4) {
          if (data->poly.density)
            sprintf(stext, "Polygons: %3.1f%% done",
                    100.*((float)data->poly.index[2]*sq(data->poly.density)+
                    (float)data->poly.index[1]*data->poly.density+
                    (float)data->poly.index[0])/(data->poly.density*
                    data->poly.density*data->poly.density*4)); break; }
        if (!data->poly.density || !data->poly.refine) break;
        s = sq(1.75*data->mol.maxcheck/data->poly.density);
        e = sq(data->mol.maxcheck/data->poly.density/data->poly.refine);
        if (e==s) break;
        y = e*s/(data->poly.processmeter*2*(s-e))+0.5*s/(e-s)+1;
        if (y<0.5)  y = 0.5;  if (y>1)  y = 1;
        sprintf(stext, "Polygons: %3.1f%% done", y*100); break;
      case 4: if (data->render.process<3) {
          strcpy(stext, "Raytrace: initializing");  break; }
        if (!data->render.w || !data->render.h) break;
        y = (float)(data->render.y*data->render.w+data->render.x)/
            (data->render.w*data->render.h);
        if (data->render.process==3)  y /= 256;
        if (data->render.process==4)  y = 1./256+255./256*(y/16);
        if (data->render.process==5)  y = 1./16+15./16*y;
        if (data->render.process==6)  y = 0.5+0.5*y;
        if ((data->render.antialias&2) && data->render.process<6)
          y *= 0.5;
        sprintf(stext, "Raytrace: %3.1f%% done", y*100); break; }
    SendMessage(HwndStat, SB_SETTEXT, 1, (LPARAM)stext);
    if (data->mol.nump>1)
      strcpy(stext, "Molecule");
    else if (data->mol.orb[0].l<21)
      sprintf(stext, "%d%c%d", data->mol.orb[0].n,
              OrbLet[data->mol.orb[0].l], data->mol.orb[0].m);
    else
      sprintf(stext, "%d,%d,%d", data->mol.orb[0].n, data->mol.orb[0].l,
              data->mol.orb[0].m);
    if ((data->dflag>>10)&1)
      sprintf(stext, "Frame %d", data->frame);
    SendMessage(HwndStat, SB_SETTEXT, 2, (LPARAM)stext); }
  inloop = 0;
  unlock_window(data);
}

HANDLE unlock_window(DATA *data)
/* Lock down all of the memory areas of a window.
 * Enter: HANDLE hnd: handle to main data array for window.
 * Exit:  DATA *data: pointer to main data array.              3/11/97-DWM */
{
  HANDLE hnd;
  long i;

  if (!data)  return(0);
  for (i=0; i<NumLocked; i++)
    if (LockList[i*3+1]==data) {
      LockList[i*3+2]--;
      if (LockList[i*3+2])
        return(GetWindowLong(LockList[i*3], GWL_USERDATA));
      break; }
  if (i<NumLocked) {
    memmove(LockList+i*3, LockList+i*3+3, (NumLocked-i-1)*3*sizeof(long));
    NumLocked--; }
  data->mol.orb      = unlock2(data->mol.orb);
  data->mol.ls       = unlock2(data->mol.ls);
  data->asym.xyz     = unlock2(data->asym.xyz);
  data->asym.elem    = unlock2(data->asym.elem);
  data->asym.scrxyz  = unlock2(data->asym.scrxyz);
  data->asym.norm    = unlock2(data->asym.norm);
  data->asym.enorm   = unlock2(data->asym.enorm);
  data->asym.rfac    = unlock2(data->asym.rfac);
  data->asym.tfac    = unlock2(data->asym.tfac);
  data->asym.pfac    = unlock2(data->asym.pfac);
  data->points.xyz   = unlock2(data->points.xyz);
  data->points.scrxyz= unlock2(data->points.scrxyz);
  data->points.phase = unlock2(data->points.phase);
  data->poly.xyz     = unlock2(data->poly.xyz);
  data->poly.elem    = unlock2(data->poly.elem);
  data->poly.scrxyz  = unlock2(data->poly.scrxyz);
  data->poly.norm    = unlock2(data->poly.norm);
  data->poly.enorm   = unlock2(data->poly.enorm);
  data->poly.phase   = unlock2(data->poly.phase);
  data->poly.ptphase = unlock2(data->poly.ptphase);
  data->poly.rfac    = unlock2(data->poly.rfac);
  data->poly.tfac    = unlock2(data->poly.tfac);
  data->poly.pfac    = unlock2(data->poly.pfac);
  data->render.buf   = unlock2(data->render.buf);
  data->render.phase = unlock2(data->render.phase);
  data->stereo.image = unlock2(data->stereo.image);
  for (i=0; i<4; i++)  if (data->seq[i]) {
    data->seq[i]->mol.orb = unlock2(data->seq[i]->mol.orb);
    data->seq[i]->mol.ls  = unlock2(data->seq[i]->mol.ls);
    data->seq[i]          = unlock2(data->seq[i]); }
  /** Additional data items get unlocked here **/
  hnd = unlock2(data);
  return(hnd);
}

HANDLE unlock2(void *mem)
/* Unlock a memory pointer which was allocated using malloc2, returning the
 *  Windows handle associated with it.  This piece of memory is no longer
 *  allocated using free2 and malloc2, but instead requires the dumb Windows
 *  routines.
 * Enter: void *mem: pointer to memory.
 *        HANDLE handle: handle to the memory.                 2/12/97-DWM */
{
  long i;
  HGLOBAL hmem;

  if (!mem)  return(0);
  for (i=0; i<mallocnum; i++)
    if ((long)mem==malloctbl[i*2+1]) {
      hmem = (HGLOBAL)malloctbl[i*2];
      GlobalUnlock(hmem);
      memmove(malloctbl+i*2, malloctbl+i*2+2,
                      (mallocnum-i-1)*2*sizeof(long));
      mallocnum--;
      return(hmem); }
  return(0);
}

warning(HWND hwnd2, char *text)
/* Report a warning, if warnings are currently viewed.
 * Enter: HWND hwnd: handle of the window or null to use best-guess top
 *                   window.
 *        char *text: warning text.
 * Exit:  long okay: 0 for canceled, 1 for okay.               4/24/97-DWM */
{
  long i;
  HWND hwnd;

  if (!hwnd2 && NumWindows)
    hwnd = GetTopWindow(0);
  if (!hwnd2 && !NumWindows)  hwnd = Hwnd;
  if (hwnd2)  hwnd = hwnd2;
  if (!Pref.warn)  return(0);
  Busy += 10;
  cursor(0);
  i = MsgBox(hwnd, text, "Warning", MB_OKCANCEL);
  Busy -= 10;
  if (i==IDCANCEL)
    return(0);
  else
    return(1);
}

HWND window_handle(DATA *data)
/* Return the window handle associated with a particular locked data array.
 * Enter: DATA *data: pointer to window data.
 * Exit:  HWND hwnd: handle of owner window.                    7/3/97-DWM */
{
  long i;

  for (i=0; i<NumLocked; i++)
    if (LockList[i*3+1]==data)
      return(LockList[i*3]);
  return(0);
}

APIENTRY WinMain(HINSTANCE inst, HINSTANCE prev, LPSTR argv, long winmode)
/* Enter: HINSTANCE inst: number identifying this program.
 *        HINSTANCE prev: 0.
 *        LPSTR argv: single string containing command line parameters.
 *        long winmode: how window is displayed at start.       5/4/96-DWM */
{
  CLIENTCREATESTRUCT ccs;
  MSG msg;
  WNDCLASS wcl;
  HACCEL acc;
  HDC hdc;
  RECT rect;
  long x, y;
  HWND hwnd;

  memset(lastview, 0, 1024);
  hinst = wcl.hInstance = inst;
  if (hwnd=FindWindow(WinName, 0)) {
    SetForegroundWindow(hwnd);
    if (IsIconic(hwnd))  ShowWindow(hwnd, SW_RESTORE);
    pass_arguments(hwnd, argv);
    return(0); }
  read_ini();
  show_splash();
  wcl.lpszClassName = WinName;
  wcl.lpfnWndProc = main_loop;
  wcl.style = 0;
  wcl.hIcon = LoadIcon(inst, ICON1);
  wcl.hCursor = LoadCursor(0, IDC_ARROW);
  wcl.lpszMenuName = "MainMenu";
  wcl.cbClsExtra = wcl.cbWndExtra = 0;
  wcl.hbrBackground = (HBRUSH)(COLOR_APPWORKSPACE+1);
  if (!RegisterClass(&wcl))  return(0);
  Hwnd = CreateWindow(WinName, Program, WS_OVERLAPPEDWINDOW|WS_CLIPCHILDREN,
                      0, 0, 640, 480, HWND_DESKTOP, 0, inst, 0);
  Hmenu = GetMenu(Hwnd);
  read_ini_toolbar(0);
  if (WinPlace.length) {
    SetWindowPlacement(Hwnd, &WinPlace);
    winmode = WinPlace.showCmd; }
  if (!start_common())  return(0);
  ShowWindow(Hwnd, winmode);
  DragAcceptFiles(Hwnd, 1);
  hdc = GetDC(Hwnd);
  BitsPixel = GetDeviceCaps(hdc, BITSPIXEL);
  ReleaseDC(Hwnd, hdc);
  opacity_dither();
  make_palette();
  start_ctl3d();
  acc = LoadAccelerators(inst, "MainMenu");
  ccs.hWindowMenu = GetSubMenu(Hmenu, MenuCascade-1);
  ccs.idFirstChild = 1000;
  HwndC = CreateWindow("MDICLIENT", 0, WS_VISIBLE|WS_CHILD|WS_CLIPCHILDREN|
                    WS_VSCROLL|WS_HSCROLL, 0, 0, 0, 0, Hwnd, 0, hinst, &ccs);
  GetWindowPlacement(HwndC, &WinPlace);
  GetClientRect(Hwnd, &WinPlace.rcNormalPosition);
  SetWindowPlacement(HwndC, &WinPlace);
  ShowWindow(HwndC, SW_SHOW);
  UpdateWindow(HwndC);
  InvalidateRect(HwndC, 0, 0);
  wcl.lpszClassName = WinName2;
  wcl.lpfnWndProc = second_loop;
  wcl.lpszMenuName = wcl.cbClsExtra = wcl.cbWndExtra = 0;
  wcl.hbrBackground = 0;
  if (!RegisterClass(&wcl))
    return(0);
  setup();
  if (argv) if (strlen(argv)) {
    strncpy(OpenName, argv, NAMELEN-1);
    OpenName[NAMELEN-1] = 0;
    open_window_mid(); }
  while (GetMessage(&msg, 0, 0, 0)) {
    if (!TranslateMDISysAccel(HwndC, &msg) &&
        !TranslateAccelerator(Hwnd, acc, &msg)) {
      TranslateMessage(&msg);
      DispatchMessage(&msg); } }
  WinHelp(Hwnd, HelpFile, HELP_QUIT, 0);
  if (Hpal)  DeleteObject(Hpal);
  if (HpalSplash)  DeleteObject(HpalSplash);
  DragAcceptFiles(Hwnd, 0);
  end_ctl3d();
  end_avi();
  end_common();
  free2_all();
  return(msg.wParam);
}

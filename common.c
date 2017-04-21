/* This contains all routines not contained in ORB.C, MATRIX.C, or DLT.C that
 *  are identical or very similar between the ANSI C and Windows versions.
 *  The ANSI version has the constant ANSIC defined.           1/18/98-DWM */

/* Uncomment the below two lines for crummy c compilers.                   */
/* #define stricmp strcasecmp */
/* #define strnicmp strncasecmp */

long Endian=0;                      /* 0-little endian (IBM), 1-big endian */
float DefSize=0;
char OrbLet[]="spdfghiklmnoqrtuvwxyz";
char VRMLCheck1[]="#VRML V1.0 ascii";

char *AsymKey[]={"Density","Opacity","Style",0};
char *AtomKey[]={
     "n",             "l",              "m",     "Neutrons(N)","Protons(Z)","Mass(kg)","CenterX(m)","CenterY(m)","CenterZ(m)","AngleAlpha(rad)",
     "AngleBeta(rad)","AngleGamma(rad)","Factor", 0};
char *CutKey[]={"Type","Surface","Invert","PositionX(m)","PositionY(m)","PositionZ(m)","AngleAlpha(rad)","AngleBeta(rad)","AngleGamma(rad)",0};
char *LightKey[]={"PositionX(m)","PositionY(m)","PositionZ(m)","Intensity","Ambience","Local",0};
char *OrbKey[]={
     "DefaultPerspective","BackgroundColor","PositiveColor", "NegativeColor","AsymptoteColor","PreviewColor",     "UseQuickRender",  "QuickRenderMode", "RenderMode",      "FixedSize",
     "FixedWidth",        "FixedHeight",    "FileName",      "Scale(m)",     "Perspective",   "LastWidth",        "LastHeight",      "CameraCenterX(m)","CameraCenterY(m)","CameraCenterZ(m)",
     "CameraTheta(rad)",  "CameraPhi(rad)", "CameraPsi(rad)","CameraCx",     "Atom",          "EndAtom",          "Light",           "EndLight",        "Psi^2(log10)",    "Cutaway",
     "EndCutaway",        "Stereo",         "EndStereo",     "Points",       "EndPoints",     "Polygons",         "EndPolygons",     "Asymptotes",      "EndAsymptotes",   "Raytrace",
     "EndRaytrace",       "Frame",          "LastFrame",     "Sequence",     "SequenceBase",  "SequenceIncrement","SequenceFileType","EndOfFile",       "SequenceBezier",  "FramesPerSecond",
     0};
char *OrbPhrase[]={
     "No",         "Yes",       "Points", "Polygons",  "Raytrace","None",       "Plane","Corner","Wedge","Monoscopic",
     "Stereoscope","Interlaced","RedBlue","Stereogram","Overlay", "Chromadepth","Left", "Right", "Solid","Wireframe",
     "Points",     "PPM",       "TIF",    "BMP",       "AVI",     "End",     0};
char *PointsKey[]={"Points",0};
char *PolyKey[]={"Density","Refinement","PositiveOpacity","NegativeOpacity","Style",0};
char *RenderKey[]={
     "PosProbOpacity","PosSurfaceOpacity","PosInteriorOpacity","NegProbOpacity","NegSurfaceOpacity","NegInteriorOpacity","AsymptoteOpacity","PositiveRefraction","NegativeRefraction","AsymptoteRefraction",
     "Steps",         "AutoBrightness","Coarse",         "Antialias",     0};
char *StereoKey[]={"Mode","TopScanline","RedSide","AutoSeparation","UseImage","Interocular","StereoscopeInter","Separation(m)","ImageName",0};

void camera_rotate(double *renderval, double *renderdlt, double *phys,
                   long dir)
/* Rotate the camera coordinates from the orbital coordinate system to the
 *  camera coordinate system, or vice versa.
 * Enter: double *renderval: pointer to initial renderval.
 *        double *renderdlt: pointer to initial dlt parameters.
 *        double *phys: array of 10 physical parameters, a la DLT.  The first
 *                      three and last three values are modified.
 *        long dir: 0 to rotate to camera coordinates, 1 to rotate to
 *                  orbital coordinates.                        1/4/98-DWM */
{
  double m[10], x, y, z, maxx, phys2[10];
  lmatrix r;
  long i;

  r.w = r.h = 3;  r.m = m;
  if (!dir) {
    render_physical(renderval, renderdlt, phys, &r);
    x = -renderval[4];
    y = -renderval[5];
    z = -renderval[6];
    phys[0] = x*m[0]+y*m[1]+z*m[2];
    phys[1] = x*m[3]+y*m[4]+z*m[5];
    phys[2]=x*m[6]+y*m[7]+z*m[8]+fabs(renderval[0])*renderval[1];
    phys[4] = phys[3];
    maxx = max(max(fabs(phys[0]), fabs(phys[1])), fabs(phys[2]));
    if (!maxx)  maxx = 1;
    phys[7] += PI;  phys[8] += PI;
    for (i=0; i<3; i++) {
      if (fabs(phys[i])/maxx<1e-6)  phys[i] = 0;
      while (phys[7+i]<-PI)  phys[7+i] += 2*PI;
      while (phys[7+i]>PI)   phys[7+i] -= 2*PI;
      if (fabs(phys[7+i])<1e-5)  phys[7+i] = 0; } }
  else {
    langtomat(m, phys[7]-PI, phys[8]-PI, phys[9]);
    phys[2] -= fabs(renderval[0])*renderval[1];
    renderval[4] = -(phys[0]*m[0]+phys[1]*m[3]+phys[2]*m[6]);
    renderval[5] = -(phys[0]*m[1]+phys[1]*m[4]+phys[2]*m[7]);
    renderval[6] = -(phys[0]*m[2]+phys[1]*m[5]+phys[2]*m[8]);
    render_physical(renderval, renderdlt, phys2, 0);
    phys2[7] = phys[7]-PI;  phys2[8] = phys[8]-PI;  phys2[9] = phys[9];
    phys2[4] = phys2[3];
    render_move(0,0,0,-1,0,renderval,phys2,renderdlt); }
}

long get_color(DATA *data, long num)
/* Get a color value.  If the color has not been defined, select an
 *  appropriate color based on the background color.
 * Enter: DATA *data: pointer to window structure, or null to use global
 *                    preferences
 *        long num: number of the color to obtain.  These are defined in
 *                  OV.H.
 * Exit:  long color: 0xBBGGRR color to use.                   3/18/97-DWM */
{
  long clr;
  long def[]={0x000000,0x0000FF,0xFF8000,0x808080,0x8080FF};
  /** Add more colors here **/

  if (!data)
    #ifndef ANSIC
      clr = Pref.colors[num];
    #else
      clr = -1;
    #endif
  else
    clr = data->pref.colors[num];
  if (clr!=-1)  return(clr);
  return(def[num]);
}

double interpolate(double *val, long *time)
/* Interpolate or extrapolate between 2, 3, or 4 values using either a
 *  linear, quadratic, cubic function, or bezier cubic.  The interpolation
 *  functions are:
 *   First order:  v = C t + D
 *   Second order: v = B t^2 + C t + D
 *   Third order:  v = A t^3 + B t^2 + C t + D
 *   Bezier:       see program
 * Enter: double *val: array of 2, 3, or 4 values to interpolate between.
 *        long *time: time[0] is the order (1, 2, or 3) plus 10 for bezier
 *                    (11, 12, or 13).  time[1] is the current time (t) to
 *                    interpolate the value for.  time[2] through time[5]
 *                    contain 2, 3, or 4 time values specifying when the
 *                    values occured.  The times must not be equal to each
 *                    other.                                   1/13/98-DWM */
{
  double A=0, B=0, C, D, value, rc[20], l[5], coef[4], roots[3], t;
  lmatrix r;
  long i, order=time[0], rn;

  if (order==11)  order = 1;
  r.m = rc;
  C = (val[0]-val[1])/(time[2]-time[3]);
  D = (time[2]*val[1]-time[3]*val[0])/(time[2]-time[3]);
  t = time[1];
  switch (order) {
    case 1: if (val[0]==val[1])  return(val[0]); break;
    case 2: if (val[0]==val[1] && val[0]==val[2])  return(val[0]);
      r.w = 4;  r.h = 3;
      lmatzero(&r);
      for (i=0; i<3; i++) {
        l[0] = sq(time[i+2]);  l[1] = time[i+2];  l[2] = 1;  l[3] = val[i];
        lleast(&r, l, 1); }
      lmatsym(&r);
      lrowreduce(&r);
      B = rc[3];  C = rc[7];  D = rc[11]; break;
    case 3:
      if (val[0]==val[1] && val[0]==val[2] && val[0]==val[3]) return(val[0]);
      r.w = 5;  r.h = 4;
      lmatzero(&r);
      for (i=0; i<4; i++) {
        l[0] = time[i+2]*time[i+2]*time[i+2];  l[1] = sq(time[i+2]);
        l[2] = time[i+2];  l[3] = 1;  l[4] = val[i];
        lleast(&r, l, 1); }
      lmatsym(&r);
      lrowreduce(&r);
      A = rc[4];  B = rc[9];  C = rc[14];  D = rc[19]; break;
    case 12: if (val[0]==val[1] && val[0]==val[2])  return(val[0]);
      if (time[1]==time[2])  return(val[0]);
      if (time[1]==time[4])  return(val[2]);
      coef[0] = (double)time[4]-(double)time[2];
      coef[1] = -3.*time[3]+3.*time[2];
      coef[2] = 3.*time[3]-3.*time[2];
      coef[3] = (double)time[2]-(double)time[1];
      if (!(rn=lcubic(coef, roots))) return(val[0]);
      t = roots[0];
      for (i=0; i<rn; i++)
        if (roots[i]>=0 && roots[i]<=1)  t = roots[i];
      A = val[2]-val[0];      B = -3*val[1]+3*val[0];
      C = 3*val[1]-3*val[0];  D = val[0]; break;
    case 13:
      if (val[0]==val[1] && val[0]==val[2] && val[0]==val[3]) return(val[0]);
      if (time[1]==time[2])  return(val[0]);
      if (time[1]==time[5])  return(val[3]);
      coef[0] = (double)time[5]-3.*time[4]+3.*time[3]-(double)time[2];
      coef[1] = 3.*time[4]-6.*time[3]+3.*time[2];
      coef[2] = 3.*time[3]-3.*time[2];
      coef[3] = (double)time[2]-(double)time[1];
      if (!(rn=lcubic(coef, roots))) return(val[0]);
      t = roots[0];
      for (i=0; i<rn; i++)
        if (roots[i]>=0 && roots[i]<=1)  t = roots[i];
      A = val[3]-3*val[2]+3*val[1]-val[0];  B = 3*val[2]-6*val[1]+3*val[0];
      C = 3*val[1]-3*val[0];                D = val[0]; }
  value = A*t*t*t+B*t*t+C*t+D;
  return(value);
}

long new_window_mid(DATA *data)
/* Set up a data array for a new window.  This is prior to loading or to a
 *  new window.  The data array has already been zeroed and certain values
 *  (such as file name) set.
 * Enter: DATA *data: pointer to window data.
 * Exit:  long okay: 0 for failed, 1 for okay.                 6/18/97-DWM */
{
  OATOM *orb;
  LIGHT *ls;
  long clr[4], i, j;

  #ifndef ANSIC
    memcpy(&data->pref, &Pref, sizeof(PREF));
  #else
    memset(&data->pref, 0, sizeof(PREF));
    data->pref.perspec = 25;
    for (i=0; i<NUMCOLORS; i++)
      data->pref.colors[i] = -1;
  #endif
  data->mol.Psi = 1e-4;  data->mol.EffScale = 1024;
  data->mol.maxjumpadj = 1.2;
  data->mol.minres = 0.1/data->mol.maxjumpadj;
  data->mol.back = (11.25/17.5/data->mol.EffScale/2/data->mol.maxjumpadj);
  if (data->mol.orb)  free2(data->mol.orb);
  if (data->mol.ls)   free2(data->mol.ls);
  data->mol.orb = malloc2(5*sizeof(OATOM));
  data->mol.ls = malloc2(5*sizeof(LIGHT));
  if (!data->mol.orb || !data->mol.ls) {
    if (data->mol.orb)  free2(data->mol.orb);  data->mol.orb = 0;
    if (data->mol.ls)   free2(data->mol.ls);   data->mol.ls  = 0;
    return(0); }
  memset(data->mol.orb, 0, 5*sizeof(OATOM));
  orb = data->mol.orb;  ls = data->mol.ls;
  orb->mass = ma;  orb->n = 4;  orb->l = 3;  orb->m = 0;  orb->Z = 1;
  orb->factor = 1;  orb->matm[0] = orb->matm[4] = orb->matm[8] = 1;
  orb->matp[0] = orb->matp[4] = orb->matp[8] = 1;
  data->mol.nump = 1;
  data->points.maxn = 10000;
  memset(data->mol.ls, 0, 5*sizeof(LIGHT));
  ls->x[0] = -80;
  ls->x[1] = ls->x[2] = 80;
  ls->i = ls->a = 1;
  data->mol.numl = 1;
  data->dflag = 9;
  data->poly.opacity[0] = data->poly.opacity[1] = 1;
  data->poly.density = 8;  data->poly.refine = 1.3;
  #ifndef ANSIC
    data->render.type = 2;
  #else
    data->render.type = 0;
  #endif
  data->render.w = data->render.h = 100;
  data->render.antialias = 1;  data->render.steps = 1000;
  data->render.brightness = 1;
  clr[0] = get_color(data, POSCOLOR);
  clr[1] = get_color(data, NEGCOLOR);
  clr[2] = get_color(data, BACKCOLOR);
  clr[3] = get_color(data, ASYMCOLOR);
  for (i=0; i<4; i++)  for (j=0; j<3; j++)
    if (!Endian)  data->render.color[i*3+j] = ((uchar *)(clr+i))[2-j];
    else          data->render.color[i*3+j] = ((uchar *)(clr+i))[j];
  data->render.opacity[1] = data->render.opacity[4] = 1;
  data->render.refrac[0]=data->render.refrac[1]=data->render.refrac[2]=1;
  data->stereo.interocular[0] = 250;  data->stereo.interocular[1] = 800;
  data->stereo.separation = 10*a0;  data->stereo.flags = 4;
  return(1);
}

long open_orb(FILE *fptr, DATA *data, long noheader)
/* Read in an ORB file.  This routine does not close the file.
 * Enter: FILE *fptr: pointer to the open file to read.
 *        DATA *data: pointer to the window's primary data array.
 *        long noheader: 0 to remove header, 1 to assume file position is
 *                       past header, +2 for no sub data arrays.
 * Exit:  long error: 0 for okay, 1 for memory error, 2 for read error.
 *                                                              1/9/98-DWM */
{
  char key[80], text[NAMELEN+2], **keylist;
  long mode=0, lval, pval, i, kval, nump, numl, seq=0;
  double rval, phys[10];
  OATOM *newp;
  LIGHT *newl;
  uchar *img;
  #ifndef ANSIC
    GV gv;
  #endif

  data->mol.nump = data->mol.numl = nump = numl = 0;
  memset(phys, 0, 10*sizeof(double));
  if (!(noheader&1))
    fscanf(fptr, "%*s");                                /* Burn off header */
  while (1) {
    if (fscanf(fptr, "%s", key)!=1)  break;
    if (key[0]=='#') {
      if (strchr(key, '\r') || strchr(key, '\n'))  continue;
      fscanf(fptr, "%*[^\r\n]");  continue; }
    if (key[0]=='}') {                            /* Nesting isn't allowed */
      if (!mode)  break;
      mode = 0;  continue; }
    fscanf(fptr, "%s", text);
    while (text[0]=='#') {
      if (!strchr(text, '\r') && !strchr(text, '\n'))
        fscanf(fptr, "%*[^\r\n]");
      text[0] = 0;
      fscanf(fptr, "%s", text); }
    if (strlen(text))
      if (text[strlen(text)-1]=='\r' || text[strlen(text)-1]=='\n')
        text[strlen(text)-1] = 0;
    if (strlen(text))
      if (text[strlen(text)-1]=='\r' || text[strlen(text)-1]=='\n')
        text[strlen(text)-1] = 0;
    if (text[0]=='{') {                           /* Nesting isn't allowed */
      if (mode) {
        fscanf(fptr, "%*[^}]%*c");  continue; }
      mode = -1; }
    switch (mode) {
      case 24: keylist = AtomKey; break;
      case 26: keylist = LightKey; break;
      case 29: keylist = CutKey; break;
      case 31: keylist = StereoKey; break;
      case 33: keylist = PointsKey; break;
      case 35: keylist = PolyKey; break;
      case 37: keylist = AsymKey; break;
      case 39: keylist = RenderKey; break;
      default: keylist = OrbKey; }
    for (i=0; keylist[i]; i++)
      if (!stricmp(key, keylist[i]))
        break;
    if (!keylist[i]) {
      if (mode==-1) {  fscanf(fptr, "%*[^}]%*c");  mode = 0; }
      continue; }
    kval = i;
    if (kval==47) break;                                            /* EOF */
    if (mode==-1) {
      mode = i;  kval = 99; }
    if (text[0]=='"') {
      if (!strchr(text+1, '"'))
        fscanf(fptr, "%[^\"]%*c", text+strlen(text));
      else
        strchr(text+1, '"')[0] = 0;
      memmove(text, text+1, strlen(text)); }
    lval = rval = pval = 0;
    sscanf(text, "%d", &lval);
    if (text[0]=='0' && toupper(text[1])=='X')
      sscanf(text+2, "%x", &lval);
    sscanf(text, "%lf", &rval);
    for (i=0; OrbPhrase[i]; i++)
      if (!stricmp(text, OrbPhrase[i])) {
        pval = i;  break; }
    switch (mode*100+kval) {
      case 0: if (rval>0)  data->pref.perspec = rval; break;
      case 1: case 2: case 3: case 4: case 5:
        if (lval>=0)  data->pref.colors[kval-1] = lval; break;
      case 6: if (pval>=0 && pval<=1)
          data->dflag = (data->dflag&0xFFFFFFFE)|pval; break;
      case 7: if (pval>=2 && pval<=4)
          data->dflag = (data->dflag&0xFFFFFFF9)|((pval-2)<<1); break;
      case 8: if (pval>=2 && pval<=4)
          data->dflag = (data->dflag&0xFFFFFFE7)|((pval-2)<<3); break;
      case 9: if (pval>=0 && pval<=1)
          data->dflag = (data->dflag&0xFFFFFDFF)|(pval<<9); break;
      case 10: if (lval>0)  data->w = lval; break;
      case 11: if (lval>0)  data->h = lval; break;
      case 12: strcpy(data->name, text); break;
      case 13: case 14: case 15: case 16:
        if (rval>0)  data->renderval[kval-13] = rval; break;
      case 17: case 18: case 19: phys[kval-17] = rval; break;
      case 20: case 21: case 22: phys[kval-20+7] = rval; break;
      case 23: phys[3] = phys[4] = rval; break;
      case 28: data->mol.Psi = pow(10, rval); break;
      case 41: data->frame = lval; break;
      case 42: data->lastframe = lval; break;
      case 44: strcpy(data->basename, text); break;
      case 45: data->incr = (pval==1); break;
      case 46: if (pval>=21 && pval<=24)  data->seqtype = pval-21; break;
      case 48: data->bezier = (pval==1); break;
      case 49: if (lval>0)  data->seqfps = lval; break;
      case 2400: if (lval>0 && lval<=MAXN)
          data->mol.orb[nump-1].n = lval; break;
      case 2401: for (i=0; i<strlen(OrbLet); i++)
          if (toupper(text[0])==toupper(OrbLet[i]))  lval = i;
        if (lval>=0 && lval<MAXN)
          data->mol.orb[nump-1].l = lval; break;
      case 2402: if (lval>-MAXN && lval<MAXN)
          data->mol.orb[nump-1].m = lval; break;
      case 2403: if (lval>=0)  data->mol.orb[nump-1].N = lval; break;
      case 2404: if (lval>0)   data->mol.orb[nump-1].Z = lval; break;
      case 2405: if (rval>0)  data->mol.orb[nump-1].mass = rval; break;
      case 2406: case 2407: case 2408:
        data->mol.orb[nump-1].x[kval-6] = rval; break;
      case 2409: case 2410: case 2411:
        data->mol.orb[nump-1].ang[kval-9] = rval; break;
      case 2412: data->mol.orb[nump-1].factor = rval; break;
      case 2499: if (data->mol.orb)
          newp = realloc2(data->mol.orb, (size_t)((nump+1)*sizeof(OATOM)));
        else
          newp = malloc2(sizeof(OATOM));
        if (!newp) {
          fscanf(fptr, "%*[^}]%*c");  mode = 0;  break; }
        data->mol.orb = newp;
        memset(data->mol.orb+nump, 0, sizeof(OATOM));
        data->mol.orb[nump].mass = ma;  data->mol.orb[nump].n = 4;
        data->mol.orb[nump].l = 3;  data->mol.orb[nump].m = 0;
        data->mol.orb[nump].Z = 1;  data->mol.orb[nump].factor = 1;
        data->mol.orb[nump].matm[0] = data->mol.orb[nump].matm[4] =
          data->mol.orb[nump].matm[8] = data->mol.orb[nump].matp[0] =
          data->mol.orb[nump].matp[4] = data->mol.orb[nump].matp[8] = 1;
        data->mol.nump++;  nump++; break;
      case 2600: case 2601: case 2602:
        data->mol.ls[numl-1].x[kval] = rval; break;
      case 2603:
        data->mol.ls[numl-1].i = rval; break;
      case 2604: data->mol.ls[numl-1].a = rval; break;
      case 2605: if (pval>=0 && pval<=1)
        data->mol.ls[numl-1].local = pval; break;
      case 2699: if (data->mol.ls)
          newl = realloc2(data->mol.ls, (size_t)((numl+1)*sizeof(LIGHT)));
        else
          newl = malloc2(sizeof(LIGHT));
        if (!newl) {
          fscanf(fptr, "%*[^}]%*c");  mode = 0;  break; }
        data->mol.ls = newl;
        memset(data->mol.ls+numl, 0, sizeof(LIGHT));
        data->mol.ls[numl].x[0] = -80;
        data->mol.ls[numl].x[1] = data->mol.ls[numl].x[2] = 80;
        data->mol.ls[numl].i = data->mol.ls[numl].a = 1;
        data->mol.numl++;  numl++; break;
      case 2900: if (pval>=5 && pval<=8)  data->cut.type = pval-5; break;
      case 2901: if (pval>=0 && pval<=1)  data->cut.nosurface = !pval; break;
      case 2902: if (pval>=0 && pval<=1)  data->cut.invert = pval; break;
      case 2903: case 2904: case 2905:
        data->cut.xyz[kval-3] = rval; break;
      case 2906: case 2907: case 2908:
        data->cut.ang[kval-6] = rval; break;
      case 3100: if (pval>=9 && pval<=15)  data->stereo.mode = pval-9; break;
      case 3101: if (pval>=16 && pval<=17)
          data->stereo.flags = (data->stereo.flags&0xFFFFFFFE)|(pval-16);
        break;
      case 3102: if (pval>=16 && pval<=17)
          data->stereo.flags=(data->stereo.flags&0xFFFFFFFD)|((pval-16)<<1);
        break;
      case 3103: if (pval>=0 && pval<=1)
          data->stereo.flags=(data->stereo.flags&0xFFFFFFFB)|(pval<<2);
        break;
      case 3104: if (pval>=0 && pval<=1)
          data->stereo.flags=(data->stereo.flags&0xFFFFFFF7)|(pval<<3);
        break;
      case 3105: case 3106: if (rval>=0)
          data->stereo.interocular[kval-5] = rval; break;
      case 3107: if (rval>0)
        data->stereo.separation = rval; break;
      case 3108: strcpy(data->stereo.name, text); break;
      case 3300: if (lval>0) data->points.maxn = lval; break;
      case 3500: if (lval>=6)  data->poly.density = lval&0xFFFFFFFE; break;
      case 3501: if (rval>=0)  data->poly.refine = rval; break;
      case 3502: case 3503:
        if (rval>=0 && rval<=1)  data->poly.opacity[kval-2] = rval; break;
      case 3504: if ((pval>=18 && pval<=20) || pval==2)
        data->poly.wire = pval-18;  if (pval==2)  data->poly.wire = 2; break;
      case 3700: if (lval>=6)  data->asym.density = lval&0xFFFFFFFE; break;
      case 3701: if (rval>=0 && rval<=1)  data->asym.opacity = rval; break;
      case 3702: if ((pval>=18 && pval<=20) || pval==2)
        data->asym.wire = pval-18;  if (pval==2)  data->asym.wire = 2; break;
      case 3900: case 3901: case 3902: case 3903: case 3904: case 3905:
      case 3906: if (rval>=0 && rval<=1)
          data->render.opacity[kval] = rval; break;
      case 3907: case 3908: case 3909:
        data->render.refrac[kval-7] = rval; break;
      case 3910: if (lval>0)  data->render.steps = lval; break;
      case 3911: if (pval>=0 && pval<=1)
          data->render.autobright = pval; break;
      case 3912: if (pval>=0 && pval<=1)
          data->render.antialias = (data->render.antialias&0xFFFFFFFE)|pval;
        break;
      case 3913: if (pval>=0 && pval<=1)
          data->render.antialias = (data->render.antialias&0xFFFFFFFD)|
                                   (pval<<1); break;
      case 4399: if (seq>=4 || (noheader&2)) break;
        if (!data->seq[seq])  data->seq[seq] = malloc2(sizeof(DATA));
        if (!data->seq[seq]) {
          fscanf(fptr, "%*[^}]%*c");  mode = 0;  break; }
        memset(data->seq[seq], 0, sizeof(DATA));
        new_window_mid(data->seq[seq]);
        open_orb(fptr, data->seq[seq], 3);
        mode = 0;
        if (!data->seq[seq]->mol.nump) {
          free2(data->seq[seq]->mol.orb);  free2(data->seq[seq]->mol.ls);
          free2(data->seq[seq]);
          data->seq[seq] = 0; }
        else
          seq++; break;
      default: ; } }
  DefSize = data->renderval[0]*0.5;
  if (!data->renderval[2])  data->renderval[2] = 100;
  if (!data->renderval[3])  data->renderval[3] = 100;
  render_dlt_new(10, 0, data->renderval[2], data->renderval[3],
                 data->renderval[1], data->renderval, data->renderdlt);
  camera_rotate(data->renderval, data->renderdlt, phys, 1);
  #ifndef ANSIC
    if (data->stereo.name[0] && !(noheader&2)) {
      memset(&gv, 0, sizeof(GV));
      gv.name = data->stereo.name;
      if (img=lock2(LoadGraphic(&gv))) {
        data->stereo.w = gv.width;
        data->stereo.h = gv.height;
        data->stereo.pal = gv.palette;
        data->stereo.image = img; } }
  #endif
  if (!strnicmp(data->name, "Untitled", 8) || !data->name[0]) {
    if (data->mol.nump>1)
      strcpy(data->name, "Molecule");
    else if (data->mol.orb[0].l<21)
      sprintf(data->name, "%d%c%d", data->mol.orb[0].n,
              OrbLet[data->mol.orb[0].l], data->mol.orb[0].m);
    else
      sprintf(data->name, "%d,%d,%d", data->mol.orb[0].n, data->mol.orb[0].l,
              data->mol.orb[0].m); }
  update_process(data);
  return(0);
}

void play_frame(DATA *data, long frame)
/* Set the orbital to a specified frame and begin rendering, setting the flag
 *  that indicates a sequence is playing.
 * Enter: DATA *data: pointer to window.
 *        long frame: frame number.                            1/13/98-DWM */
{
  long order=1, i, j, k, color[4], time[6], refresh=0;
  DATA *seq[4], old;
  double val[4], phys[50];
  OATOM *newp, oldp;
  LIGHT *newl, oldl;

  if (!data->seq[0] || !data->seq[1])  return;
  seq[0] = seq[2] = seq[3] = data->seq[0];  seq[1] = data->seq[1];
  if (data->seq[2]) {
    order = 2;  seq[2] = data->seq[2];
    if (data->seq[3]) {
      order = 3;  seq[3] = data->seq[3]; } }
  if (seq[1]->frame==seq[0]->frame)  seq[1]->frame++;
  if (order==2)
    if (seq[2]->frame==seq[0]->frame || seq[2]->frame==seq[1]->frame)
      seq[2]->frame = max(seq[0]->frame,seq[1]->frame)+1;
  if (order==3)
    if (seq[3]->frame==seq[0]->frame || seq[3]->frame==seq[1]->frame ||
        seq[3]->frame==seq[2]->frame)
      seq[3]->frame = max(seq[2]->frame,max(seq[0]->frame, seq[1]->frame))+1;
  time[0] = order+10*data->bezier;  time[1] = frame;
  for (j=0; j<4; j++)
    time[j+2] = seq[j]->frame;
  memcpy(&old, data, sizeof(DATA));
  memcpy(&oldp, data->mol.orb, sizeof(OATOM));
  memcpy(&oldl, data->mol.ls, sizeof(LIGHT));
  for (i=0; i<NUMCOLORS; i++) {                                  /* Colors */
    for (j=0; j<4; j++)
      if (seq[j]->pref.colors[i]!=-1) break;
    if (j==4)
      data->pref.colors[i] = -1;
    else {
      for (j=0; j<4; j++)
        color[j] = get_color(seq[j], i);
      data->pref.colors[i] = 0;
      for (j=0; j<3; j++) {
        for (k=0; k<4; k++)  val[k] = ((uchar *)(color+k))[j];
        ((uchar *)(data->pref.colors+i))[j] = interpolate(val, time); } } }
  for (i=0; i<4; i++)                                   /* Camera position */
    camera_rotate(seq[i]->renderval, seq[i]->renderdlt, phys+10*i+10, 0);
  camera_rotate(data->renderval, data->renderdlt, phys, 0);
  for (i=0; i<10; i++) {
    if (!data->incr || (i>=3 && i<7)) {
      for (j=0; j<4; j++)  val[j] = phys[j*10+10+i];
      phys[i] = interpolate(val, time); }
    else
      phys[i] = phys[10+i]+(phys[20+i]-phys[10+i])*(time[1]-time[2]);
    for (j=0; j<4; j++)  val[j] = seq[j]->renderval[i];
    data->renderval[i] = interpolate(val, time); }
  render_dlt_new(10, 0, data->renderval[2], data->renderval[3],
                 data->renderval[1], data->renderval, data->renderdlt);
  camera_rotate(data->renderval, data->renderdlt, phys, 1);
  for (i=0; i<4; i++)  val[i] = seq[i]->w;             /* Fixed image size */
  data->w = interpolate(val, time);
  for (i=0; i<4; i++)  val[i] = seq[i]->h;
  data->h = interpolate(val, time);
  for (i=1,j=seq[0]->mol.nump; i<4; i++)                          /* Atoms */
    if (seq[i]->mol.nump<j)  j = seq[i]->mol.nump;
  if (j<1)  j = 1;
  if (!(newp=realloc2(data->mol.orb, (size_t)(j*sizeof(OATOM)))))
    return;
  data->mol.orb = newp;  data->mol.nump = j;
  for (j=0; j<data->mol.nump; j++) {
    for (i=0; i<4; i++)  val[i] = seq[i]->mol.orb[j].mass;
    data->mol.orb[j].mass = interpolate(val, time);
    for (i=0; i<4; i++)  val[i] = seq[i]->mol.orb[j].n;
    data->mol.orb[j].n = interpolate(val, time);
    for (i=0; i<4; i++)  val[i] = seq[i]->mol.orb[j].l;
    data->mol.orb[j].l = interpolate(val, time);
    for (i=0; i<4; i++)  val[i] = seq[i]->mol.orb[j].m;
    data->mol.orb[j].m = interpolate(val, time);
    for (i=0; i<4; i++)  val[i] = seq[i]->mol.orb[j].N;
    data->mol.orb[j].N = interpolate(val, time);
    for (i=0; i<4; i++)  val[i] = seq[i]->mol.orb[j].Z;
    data->mol.orb[j].Z = interpolate(val, time);
    for (i=0; i<4; i++)  val[i] = seq[i]->mol.orb[j].factor;
    data->mol.orb[j].factor = interpolate(val, time);
    for (k=0; k<3; k++) {
      for (i=0; i<4; i++)  val[i] = seq[i]->mol.orb[j].x[k];
      data->mol.orb[j].x[k] = interpolate(val, time);
      for (i=0; i<4; i++)  val[i] = seq[i]->mol.orb[j].ang[k];
      data->mol.orb[j].ang[k] = interpolate(val, time); } }
  for (i=1,j=seq[0]->mol.numl; i<4; i++)                         /* Lights */
    if (seq[i]->mol.numl<j)  j = seq[i]->mol.numl;
  if (!(newl=realloc2(data->mol.ls, (size_t)(max(j,1)*sizeof(LIGHT)))))
    return;
  data->mol.ls = newl;  data->mol.numl = j;
  if (j)
    memcpy(data->mol.ls, seq[0]->mol.ls, (size_t)(max(j,1)*sizeof(LIGHT)));
  for (j=0; j<data->mol.numl; j++) {
    for (i=0; i<4; i++)  val[i] = seq[i]->mol.ls[j].i;
    data->mol.ls[j].i = interpolate(val, time);
    for (i=0; i<4; i++)  val[i] = seq[i]->mol.ls[j].a;
    data->mol.ls[j].a = interpolate(val, time);
    for (k=0; k<3; k++) {
      for (i=0; i<4; i++)  val[i] = seq[i]->mol.ls[j].x[k];
      data->mol.ls[j].x[k] = interpolate(val, time); }
    data->mol.ls[j].local = seq[0]->mol.ls[j].local; }
  for (i=0; i<4; i++)  val[i] = log10(seq[i]->mol.Psi);           /* Psi^2 */
  data->mol.Psi = pow(10, interpolate(val, time));
  for (i=0; i<4; i++)  val[i] = seq[i]->asym.opacity;         /* Asymptote */
  data->asym.opacity = interpolate(val, time);
  for (i=0; i<4; i++)  val[i] = seq[i]->asym.density;
  data->asym.density = interpolate(val, time);
  for (k=0; k<3; k++) {                                         /* Cutaway */
    for (i=0; i<4; i++)  val[i] = seq[i]->cut.xyz[k];
    data->cut.xyz[k] = interpolate(val, time);
    for (i=0; i<4; i++)  val[i] = seq[i]->cut.ang[k];
    data->cut.ang[k] = interpolate(val, time); }
  data->cut.type = seq[0]->cut.type;
  data->cut.nosurface = seq[0]->cut.nosurface;
  data->cut.invert = seq[0]->cut.invert;
  for (i=0; i<4; i++)  val[i] = seq[i]->points.maxn;             /* Points */
  data->points.maxn = interpolate(val, time);
  for (i=0; i<4; i++)  val[i] = seq[i]->poly.density;          /* Polygons */
  data->poly.density = interpolate(val, time);
  for (i=0; i<4; i++)  val[i] = seq[i]->poly.refine;
  data->poly.refine = interpolate(val, time);
  for (k=0; k<2; k++) {
    for (i=0; i<4; i++)  val[i] = seq[i]->poly.opacity[k];
    data->poly.opacity[k] = interpolate(val, time); }
  for (k=0; k<7; k++) {                                        /* Raytrace */
    for (i=0; i<4; i++)  val[i] = seq[i]->render.opacity[k];
    data->render.opacity[k] = interpolate(val, time); }
  for (k=0; k<3; k++) {
    for (i=0; i<4; i++)  val[i] = seq[i]->render.refrac[k];
    data->render.refrac[k] = interpolate(val, time); }
  for (i=0; i<4; i++)  val[i] = seq[i]->render.steps;
  data->render.steps = interpolate(val, time);
  data->render.autobright = seq[0]->render.autobright;
  data->render.antialias  = seq[0]->render.antialias;
  data->stereo.mode  = seq[0]->stereo.mode;                      /* Stereo */
  data->stereo.flags = seq[0]->stereo.flags;
  for (k=0; k<2; k++) {
    for (i=0; i<4; i++)  val[i] = seq[i]->stereo.interocular[k];
    data->stereo.interocular[k] = interpolate(val, time); }
  for (i=0; i<4; i++)  val[i] = seq[i]->stereo.separation;
  data->stereo.separation = interpolate(val, time);
  data->frame = frame;                                           /* Update */
  data->dflag |= (1<<10);
  if (data->mol.nump>1)  refresh = 0xF;
  if (memcmp(old.mol.orb, data->mol.orb, sizeof(OATOM)))  refresh = 0xF;
  if (data->mol.numl>1)  refresh |= 8;
  if (memcmp(old.mol.ls, data->mol.ls, sizeof(LIGHT)))  refresh |= 8;
  if (data->mol.Psi!=old.mol.Psi)  refresh |= 10;
  if (old.asym.density!=data->asym.density)  refresh |= 4;
  if (memcmp(&old.cut, &data->cut, sizeof(CUTAWAY)))  refresh = 0xF;
  if (old.points.maxn!=data->points.maxn)  refresh |= 1;
  if (old.poly.density!=data->poly.density ||
      old.poly.refine!=data->poly.refine)  refresh |= 2;
  if (memcmp(&old.render, &data->render, sizeof(RENDER)))  refresh |= 8;
  if (memcmp(&old.stereo, &data->stereo, sizeof(STEREO)))  refresh |= 8;
  if (refresh!=0xF)  update_process(data);
  if (refresh&1) {
    data->points.process = 0;
    orb_points(&data->points, 0, &data->mol, 0, 4); }
  if (refresh&2) {
    data->poly.process = 0;
    orb_polygons(&data->poly, 0, &data->mol, 0, 4); }
  if (refresh&4) {
    data->asym.process = 0;
    orb_asymptote(&data->asym, 0, &data->mol, 0, 4); }
  if (refresh&8) {
    data->render.process = 0;
    orb_render(&data->render, 0, 0, &data->mol, 0, 4); }
  #ifndef ANSIC
    data->changed = 1;
    recheck(window_handle(data), 0);
  #endif
}

void prep_stereo(DATA *data)
/* Compute the left and right camera DLTs, physical parameters, and camera
 *  matrices.
 * Enter: DATA *data: pointer to window to compute.             1/3/98-DWM */
{
  double val[10];
  float dxl[3], dxr[3], dx0l[4], dx0r[4];
  real oc, sep;

  if (data->stereo.mode!=StereoSTEREOSCOPE &&
      data->stereo.mode!=StereoINTERLACED &&
      data->stereo.mode!=StereoREDBLUE &&
      data->stereo.mode!=StereoOVERLAY)
    return;
  oc = data->stereo.interocular[data->stereo.mode==StereoSTEREOSCOPE]/
       (2+6*(data->stereo.mode!=StereoSTEREOSCOPE));
  memcpy(val, data->renderval, 10*sizeof(double));
  sep = data->stereo.separation;
  if (data->stereo.flags&4)  sep = data->mol.maxcheck*0.25;
  dxl[0] = -0.5*sep/fabs(val[0]);
  dxl[1] = dxl[2] = dxr[1] = dxr[2] = 0;  dxr[0] = -dxl[0];
  dx0l[0] = -oc*0.5;  dx0r[0] = -dx0l[0];
  dx0l[1] = dx0l[2] = dx0l[3] = dx0r[1] = dx0r[2] = dx0r[3] = 0;
  render_move(dxl, 0, dx0l, 0, 0, val, data->renderdlt, data->renderdlt+18);
  ldlt_to_physical(data->renderdlt+18, data->stereo.lrcamera, 0, 0);
  data->stereo.lrcamera[4] = data->stereo.lrcamera[3];
  oangtomat(data->stereo.lrcammat, data->stereo.lrcamera[7],
            data->stereo.lrcamera[8], data->stereo.lrcamera[9]);
  memcpy(val, data->renderval, 10*sizeof(double));
  render_move(dxr, 0, dx0r, 0, 0, val, data->renderdlt, data->renderdlt+36);
  ldlt_to_physical(data->renderdlt+36, data->stereo.lrcamera+10, 0, 0);
  data->stereo.lrcamera[14] = data->stereo.lrcamera[13];
  oangtomat(data->stereo.lrcammat+9, data->stereo.lrcamera[17],
            data->stereo.lrcamera[18], data->stereo.lrcamera[19]);
}

void render_dlt(long numnode, float *node, long w, long h, float perspec,
                double *renderval, double *initdlt, double *findlt)
/* Compute a set of dlt parameters based on an original dlt set, a geometry,
 *  and the screen size.  The initial DLT values are only used for
 *  orientation.
 * The values stored in renderval are:
 *  0 - scale: diameter of object.  Negative to indicate that the normals
 *             need updating.
 *  1 - perspec: perspective distance.
 *  2,3 - width,height of the screen.
 *  4,5,6 - view center x,y,z.
 *  7,8,9 - object center x,y,z (center of rotation).
 * Enter: long numnode: number of 3D nodes.  These will fill the screen as
 *                      best possible.  If this number is negative, the
 *                      positive value minus one is used instead, and initdlt
 *                      is used for just an initial rotation matrix.
 *        float *node: pointer to xyzxyzxyz array of coordinates.
 *        long w, h: size of destination screen.
 *        float perspec: perspective distance.  A value of 5 works well.
 *                       Smaller values produce more "bulge".
 *        double *renderval: location to store 10 double values containing
 *                           the values needed for render_move.
 *        double *initdlt: if numnode is positive, this is the initial DLT
 *                         parameters for orientation.  If numnode is
 *                         negative, this is a 3x3 rotation matrix (double
 *                         *r).
 *        double *findlt: location to store result.  May be the same as
 *                        initdlt.                            11/11/96-DWM */
{
  double size=0.9;                                        /* Size on screen */
  double minx[6], *center=renderval+7, phys[10], z;
  double scale, back=perspec;
  long i, j;
  double rc[11], ang[3];
  lmatrix *r=lmat(rc);

  if (numnode>=0)
    ldlt_to_physical(initdlt, phys, r, 0);
  else {
    numnode = -numnode-1;
    for (i=0; i<9; i++)
      r->m[i] = initdlt[i]; }
  if (!DefSize)
    DefSize = DEFSIZE;
  for (j=0; j<3; j++) {
    minx[j] = -DefSize;  minx[j+3] = DefSize; }
  for (j=0, scale=0; j<3; j++) {
    center[j] = (minx[j+3]+minx[j])/2;
    if (minx[j+3]==minx[j])  minx[j+3] += 1;
    if (minx[j+3]-minx[j]>scale)
      scale = minx[j+3]-minx[j]; }
  for (j=0; j<3; j++) { r->m[j] *= -1;  r->m[j+6] *= -1; }
  phys[0] = r->m[6]*scale*back+center[0];
  phys[1] = r->m[7]*scale*back+center[1];
  phys[2] = r->m[8]*scale*back+center[2];
  phys[3] = phys[4] = back*size*(double)min(w, h);
  phys[4] *= -1;
  phys[5] = w/2;  phys[6] = h/2;
  lmattoang(r->m, ang);
  phys[7] = ang[0];  phys[8] = ang[1];  phys[9] = ang[2];
  lphysical_to_dlt(phys, findlt, 0);
  renderval[0] = scale;  renderval[1] = back;
  renderval[2] = w;      renderval[3] = h;
  renderval[4] = renderval[5] = renderval[6] = 0;
}

void render_dlt_new(long numnode, float *node, long w, long h, float perspec,
                    double *renderval, double *findlt)
/* Compute a set of dlt parameters based on a geometry and the screen size.
 *  The default orientation is viewing the geometry 'straight on'.
 * Enter: long numnode: number of 3D nodes.  These will fill the screen as
 *                      best possible.  Zero for use default size (DefSize).
 *        float *node: pointer to xyzxyzxyz array of coordinates.
 *        long w, h: size of destination screen.
 *        float perspec: perspective distance.  A value of 25 works well.
 *                       Smaller values produce more "bulge".
 *        double *renderval: location to store 10 double values containing
 *                           the values needed for render_move.
 *        double *findlt: location to store result.  May be the same as
 *                        initdlt.                             3/13/97-DWM */
{
  double r[9]={1,0,0, 0,1,0, 0,0,-1};

  render_dlt(-numnode-1, node, w, h, perspec, renderval, r, findlt);
}

void render_move(float *dx, float *dang, float *dx0, long w, long h,
                 double *renderval, double *initdlt, double *findlt)
/* Adjust a set of dlt parameters based on specific movement, rotation, and
 *  screen size changes.
 * Enter: float *dx: array of three values.  These are translations in the
 *                   screen coordinate system, and are proportional to the
 *                   size of the object.  A value of +/- 1 would move the
 *                   object by its diameter.  Values are dx,dy,dz.  Pass a
 *                   null for no translation.
 *        float *dang: array of three angles in radians.  These are the
 *                    rotations about the x, y, and z screen axes.  Pass a
 *                    null for no rotation.
 *        float *dx0: four values.  The first two translate the screen
 *                    center, and are in screen coordinaes.  The third value
 *                    is the amount to change the perspective.  The fourth
 *                    value is a new scale parameter.  Values are dx0,dy0,
 *                    dperspec,scale.  Pass a null for no screen translation/
 *                    perspective change.
 *        long w, h: size of destination screen.  Send nulls for no screen
 *                   size change.  Width is negative if initdlt actually
 *                   contains physical values.  -1 is treated as 0 with
 *                   physical values.
 *        double *renderval: 10 double values returned by render_dlt.
 *        double *initdlt: initial DLT parameters for orientation.  If w<0,
 *                         these are the physical parameters which would
 *                         be returned by ldlt_to_physical.
 *        double *findlt: location to store result.  May be the same as
 *                        initdlt.                            11/11/96-DWM */
{
  double phys[10], *center=renderval+7, x0, y0;
  double scale=renderval[0], back=renderval[1], *view=renderval+4, v[3];
  double rc[44], ang[3];
  lmatrix *r=lmat(rc), *dr=lmat(rc+11), *rp=lmat(rc+22), *rt=lmat(rc+33);
  long ow=renderval[2], oh=renderval[3], j;
  char mv[]={4,8,5,7, 8,0,6,2, 0,4,1,3};

  if (!renderval[0])  return;
  scale = fabs(scale);
  if (w>=0)  ldlt_to_physical(initdlt, phys, r, 0);
  else {
    memcpy(phys, initdlt, 10*sizeof(double));
    langtomat(r->m, phys[7], phys[8], phys[9]);
    if (w==-1)  w = 0;   w *= -1; }
  if (!w)  w = ow;  if (!h)  h = oh;
  if (phys[4]>0) {
    phys[4] = -phys[4];  phys[9] += PI;
    langtomat(r->m, phys[7], phys[8], phys[9]); }
  phys[4] = -phys[3];                        /* Prevents non-square pixels */
  x0 = (phys[5]-ow/2)/min(ow, oh);  y0 = (phys[6]-oh/2)/min(ow, oh);
  if (dx) {
    view[0] += scale*(-dx[0]*r->m[0]+dx[1]*r->m[3]-dx[2]*r->m[6]);
    view[1] += scale*(-dx[0]*r->m[1]+dx[1]*r->m[4]-dx[2]*r->m[7]);
    view[2] += scale*(-dx[0]*r->m[2]+dx[1]*r->m[5]-dx[2]*r->m[8]); }
  if (dang) {
    renderval[0] = -scale;
    v[0] = view[0]*r->m[0]+view[1]*r->m[1]+view[2]*r->m[2];
    v[1] = view[0]*r->m[3]+view[1]*r->m[4]+view[2]*r->m[5];
    v[2] = view[0]*r->m[6]+view[1]*r->m[7]+view[2]*r->m[8];
    lmattrans(r, rt);
    for (j=0; j<3; j++)
      if (dang[j]) {
        lmatident(dr, 3);
        dr->m[mv[j*4]] = dr->m[mv[j*4+1]] = cos(dang[j]);
        dr->m[mv[j*4+2]] = sin(dang[j]);
        dr->m[mv[j*4+3]] = -sin(dang[j]);
        lmatmul(dr, r, rp);
        lmatmul(rt, rp, dr);
        lmatmul(r, dr, rp);
        lmatdup(rp, r); }
    lmattoang(r->m, ang);
    phys[7] = ang[0];  phys[8] = ang[1];  phys[9] = ang[2];
    view[0] = v[0]*r->m[0]+v[1]*r->m[3]+v[2]*r->m[6];
    view[1] = v[0]*r->m[1]+v[1]*r->m[4]+v[2]*r->m[7];
    view[2] = v[0]*r->m[2]+v[1]*r->m[5]+v[2]*r->m[8]; }
  if (dx0) {
    x0 += dx0[0]/min(w, h);  y0 += dx0[1]/min(w, h);
    if (dx0[2] && fabs(dx0[2]+back)>0.01) {
      phys[3] /= back;  phys[4] /= back;
      back += dx0[2];
      phys[3] *= back;  phys[4] *= back; }
    if (dx0[3])
      renderval[0] = scale = dx0[3]; }
  phys[0] = r->m[6]*scale*back+center[0]+view[0];
  phys[1] = r->m[7]*scale*back+center[1]+view[1];
  phys[2] = r->m[8]*scale*back+center[2]+view[2];
  phys[3] *= (double)min(w, h)/min(ow, oh);
  phys[4] *= (double)min(w, h)/min(ow, oh);
  phys[5] = w/2+x0*min(w,h);  phys[6] = h/2+y0*min(w,h);
  lphysical_to_dlt(phys, findlt, 0);
  renderval[1] = back;
  renderval[2] = w;  renderval[3] = h;
}

void render_physical(double *val, double *dlt, double *phys, lmatrix *ang)
/* Convert a set of 11 dlt parameters and 10 renderval parameters to 10
 *  physical parameters.  Distortion is ignored.  This routine can fail since
 *  there are singularities in the equations.
 * Enter: double *val: array of 10 renderval parameters.
 *        double *dlt: array with 11 dlt parameters, L1 to L11.
 *        double *phys: array to store physical parameters in the order
 *                      X0, Y0, Z0, Cx, Cy, x0, y0, theta, phi, psi.  Null
 *                      for only matrix output.
 *        lmatrix *ang: 3x3 matrix to store rotation in.  Pass a zero pointer
 *                      if no output is desired.               1/28/98-DWM */
{
  double size=0.9, rc[11];
  lmatrix *r=lmat(rc);

  if (!ang)  ang = r;
  ldlt_to_physical(dlt, phys, ang, 0);
  if (!phys)  return;
  phys[3] = phys[4] = size*val[1]*(double)min(val[2], val[3]);
  phys[5] = val[2]*0.5;  phys[6] = val[3]*0.5;
  return;
}

long save_digistar(DATA *data)
/* Save the polygon and point data to a Digistar file.  The file name is in
 *  OpenName.  If there is both point and polygon data, the point data is
 *  saved.
 * Enter: DATA *data: pointer to window data.
 * Exit:  long error: 0 for okay, 1 for insufficient memory, 2 for write
 *                    error.                                    1/4/98-DWM */
{
  FILE *fptr;
  long clr, i, j, k;

  #ifndef ANSIC
    if (!(fptr=fopen(OpenName, "wt")))  return(2);
  #else
    fptr = Outfptr;
  #endif
  if (data->points.n) {
    for (i=0; i<data->points.n; i++)
      fprintf(fptr, "D %g %g %g 1.0\n", data->points.xyz[i*3]*DIGISCALE,
             data->points.xyz[i*3+1]*DIGISCALE,
             data->points.xyz[i*3+2]*DIGISCALE); }
  if (!data->points.n && data->asym.e) {
    for (i=0; i<data->asym.n; i++)
      fprintf(fptr, "D %g %g %g 1.0\n", data->asym.xyz[i*3]*DIGISCALE,
              data->asym.xyz[i*3+1]*DIGISCALE,
              data->asym.xyz[i*3+2]*DIGISCALE); }
  if (!data->points.n && data->poly.e) {
    for (i=0; i<data->poly.n; i++)
      fprintf(fptr, "D %g %g %g 1.0\n", data->poly.xyz[i*3]*DIGISCALE,
              data->poly.xyz[i*3+1]*DIGISCALE,
              data->poly.xyz[i*3+2]*DIGISCALE); }
  #ifndef ANSIC
    fclose(fptr);
  #endif
  return(0);
}

long save_vrml(DATA *data)
/* Save the polygon and point data to a VRML version 1.0 file.  The file name
 *  is in OpenName.
 * Enter: DATA *data: pointer to window data.
 * Exit:  long error: 0 for okay, 1 for insufficient memory, 2 for write
 *                    error.                                    1/4/98-DWM */
{
  FILE *fptr;
  long clr, i, j, k;

  #ifndef ANSIC
    if (!(fptr=fopen(OpenName, "wt")))  return(2);
  #else
    fptr = Outfptr;
  #endif
  fprintf(fptr, "%s\n", VRMLCheck1);
  if (data->asym.e) {
    clr = get_color(data, ASYMCOLOR);
    fprintf(fptr, "  Separator {\nMaterial {\n");
    if (!((data->pref.flags>>3)&1))
      save_vrml_color(fptr, "ambientColor [ %g %g %g ]\n", clr);
    if (!((data->pref.flags>>2)&1))
      save_vrml_color(fptr, "diffuseColor [ %g %g %g ]\n", clr);
    fprintf(fptr, "transparency %lg }\n"
           "Coordinate3 { point [\n", 1-data->asym.opacity);
    for (i=0; i<data->asym.n; i++)
      fprintf(fptr, "%g %g %g,\n", data->asym.xyz[i*3]*VRMLSCALE,
           data->asym.xyz[i*3+1]*VRMLSCALE, data->asym.xyz[i*3+2]*VRMLSCALE);
    fprintf(fptr, "] } IndexedFaceSet { coordIndex [\n");
    for (i=0; i<data->asym.e; i++)
      fprintf(fptr, "%d, %d, %d, -1,\n", data->asym.elem[i*3],
             data->asym.elem[i*3+1], data->asym.elem[i*3+2]);
    fprintf(fptr, "] } }\n"); }
  if (data->poly.e) {
    for (i=0,j=k=1; i<data->poly.n; i++)
      if (data->poly.ptphase[i]>0) { data->poly.ptphase[i] = j;  j++; }
      else {                         data->poly.ptphase[i] =-k;  k++; }
    j--;  k--;
    if (j) {
      clr = get_color(data, POSCOLOR);
      fprintf(fptr, "  Separator {\nMaterial {\n");
      if (!((data->pref.flags>>3)&1))
        save_vrml_color(fptr, "ambientColor [ %g %g %g ]\n", clr);
      if (!((data->pref.flags>>2)&1))
        save_vrml_color(fptr, "diffuseColor [ %g %g %g ]\n", clr);
      fprintf(fptr, "transparency %lg }\n"
              "Coordinate3 { point [\n", 1-data->poly.opacity[0]);
      for (i=0; i<data->poly.n; i++)  if (data->poly.ptphase[i]>0)
        fprintf(fptr, "%g %g %g,\n", data->poly.xyz[i*3]*VRMLSCALE,
           data->poly.xyz[i*3+1]*VRMLSCALE, data->poly.xyz[i*3+2]*VRMLSCALE);
      fprintf(fptr, "] } IndexedFaceSet { coordIndex [\n");
      for (i=0; i<data->poly.e; i++)  if (data->poly.phase[i]>0)
        fprintf(fptr, "%d, %d, %d, -1,\n",
               data->poly.ptphase[data->poly.elem[i*3]]-1,
               data->poly.ptphase[data->poly.elem[i*3+1]]-1,
               data->poly.ptphase[data->poly.elem[i*3+2]]-1);
      fprintf(fptr, "] } }\n"); }
    if (k) {
      clr = get_color(data, NEGCOLOR);
      fprintf(fptr, "  Separator {\nMaterial {\n");
      if (!((data->pref.flags>>3)&1))
        save_vrml_color(fptr, "ambientColor [ %g %g %g ]\n", clr);
      if (!((data->pref.flags>>2)&1))
        save_vrml_color(fptr, "diffuseColor [ %g %g %g ]\n", clr);
      fprintf(fptr, "transparency %lg }\n"
              "Coordinate3 { point [\n", 1-data->poly.opacity[1]);
      for (i=0; i<data->poly.n; i++)  if (data->poly.ptphase[i]<0)
        fprintf(fptr, "%g %g %g,\n", data->poly.xyz[i*3]*VRMLSCALE,
           data->poly.xyz[i*3+1]*VRMLSCALE, data->poly.xyz[i*3+2]*VRMLSCALE);
      fprintf(fptr, "] } IndexedFaceSet { coordIndex [\n");
      for (i=0; i<data->poly.e; i++)  if (data->poly.phase[i]<0)
        fprintf(fptr, "%d, %d, %d, -1,\n",
               -data->poly.ptphase[data->poly.elem[i*3]]-1,
               -data->poly.ptphase[data->poly.elem[i*3+1]]-1,
               -data->poly.ptphase[data->poly.elem[i*3+2]]-1);
      fprintf(fptr, "] } }\n"); } }
  if (data->points.n && !data->poly.e) {
    clr = get_color(data, POSCOLOR);
    fprintf(fptr, "  Separator {\nMaterial {\n");
    if (!((data->pref.flags>>3)&1))
      save_vrml_color(fptr, "ambientColor [ %g %g %g ]\n", clr);
    if (!((data->pref.flags>>2)&1))
      save_vrml_color(fptr, "diffuseColor [ %g %g %g ]\n", clr);
    if (!((data->pref.flags>>4)&1))
      save_vrml_color(fptr, "emissiveColor [ %g %g %g ]\n", clr);
    fprintf(fptr, "}\nCoordinate3 { point [\n");
    for (i=0; i<data->points.n; i++)  if (data->points.phase[i]>0)
      fprintf(fptr, "%g %g %g,\n", data->points.xyz[i*3]*VRMLSCALE,
             data->points.xyz[i*3+1]*VRMLSCALE,
             data->points.xyz[i*3+2]*VRMLSCALE);
    fprintf(fptr, "] }\nPointSet { startIndex 0 numPoints -1 } }\n");
    clr = get_color(data, NEGCOLOR);
    fprintf(fptr, "  Separator {\nMaterial {\n");
    if (!((data->pref.flags>>3)&1))
      save_vrml_color(fptr, "ambientColor [ %g %g %g ]\n", clr);
    if (!((data->pref.flags>>2)&1))
      save_vrml_color(fptr, "diffuseColor [ %g %g %g ]\n", clr);
    if (!((data->pref.flags>>4)&1))
      save_vrml_color(fptr, "emissiveColor [ %g %g %g ]\n", clr);
    fprintf(fptr, "}\nCoordinate3 { point [\n");
    for (i=0; i<data->points.n; i++)  if (data->points.phase[i]<0)
      fprintf(fptr, "%g %g %g,\n", data->points.xyz[i*3]*VRMLSCALE,
             data->points.xyz[i*3+1]*VRMLSCALE,
             data->points.xyz[i*3+2]*VRMLSCALE);
    fprintf(fptr, "] }\nPointSet { startIndex 0 numPoints -1 } }\n"); }
  #ifndef ANSIC
    fclose(fptr);
  #endif
  return(0);
}

void save_vrml_color(FILE *fptr, char *text, long clr)
/* Write a color to a VRML file.
 * Enter: FILE *fptr: pointer to file to write to.
 *        char *text: text including three %g's for color output.
 *        long color: color in whatever endian we're in.       1/18/98-DWM */
{
  if (!Endian)
    fprintf(fptr, text, (float)(clr>>16)/255, (float)((clr>>8)&0xFF)/255,
            (float)(clr&0xFF)/255);
  else
    fprintf(fptr, text, (float)(clr&0xFF)/255, (float)((clr>>8)&0xFF)/255,
            (float)(clr>>16)/255);
}

void update_process(DATA *data)
/* Update the processing of all of the different rendering techniques, such
 *  that the appropriate pictures will be drawn.
 * Enter: DATA *data: pointer to window data.                 12/29/97-DWM */
{
  long m, i, j, clr[4];

  clr[0] = get_color(data, POSCOLOR);
  clr[1] = get_color(data, NEGCOLOR);
  clr[2] = get_color(data, BACKCOLOR);
  clr[3] = get_color(data, ASYMCOLOR);
  for (i=0; i<4; i++)  for (j=0; j<3; j++)
    if (!Endian)  data->render.color[i*3+j] = ((uchar *)(clr+i))[2-j];
    else          data->render.color[i*3+j] = ((uchar *)(clr+i))[j];
  render_physical(data->renderval, data->renderdlt, data->render.camera, 0);
  data->render.camera[4] = data->render.camera[3];
  if (!((data->dflag>>9)&1) || data->w<=0 || data->h<=0) {
    data->render.w = data->renderval[2];
    data->render.h = data->renderval[3]; }
  else {
    data->render.w = data->w;
    data->render.h = data->h; }
  if (data->points.process>=5)  data->points.process = 6;
  if (data->asym.process>=5)    data->asym.process = 6;
  if (data->poly.process>=6)    data->poly.process = 7;
  if (data->render.process>=3)  data->render.process = 2;
  if (data->dflag&1) {
    m = (data->dflag>>1)&3;
    data->dflag = (data->dflag&0xFFFFFF9F)|(m<<5); }
  data->renderval[0] = -fabs(data->renderval[0]);
  data->render.antialias &= 0xFFFFFFF3;
  prep_stereo(data);
}

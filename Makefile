CL_OPTS=/O2 /GL /W2 /w24101 /w34244 /w34305 /nologo /I C:\P\lib
LINK_OPTS=/LTCG /libpath:C:\P\lib /nologo
RC_OPTS=/n /nologo

all: ansiorb.exe ov.exe


ansiorb.exe: ansi.obj dlt.obj matrix.obj orbansi.obj
	link $(LINK_OPTS) /OUT:ansiorb.exe ansi.obj dlt.obj matrix.obj orbansi.obj 
ov.exe: dlt.obj draw.obj file.obj matrix.obj orb.obj ov.obj ov.res
	link $(LINK_OPTS) /OUT:ov.exe dlt.obj draw.obj file.obj matrix.obj orb.obj ov.obj ov.res  gvlib.lib mem.lib  advapi32.lib comctl32.lib comdlg32.lib gdi32.lib kernel32.lib shell32.lib user32.lib


ansi.obj: ansi.c common.c ansi.h orb.h
	cl /c /Za $(CL_OPTS) ansi.c
data.h: u\splash.jpg u\preview.jpg u\header.avi
	data2c u\splash.jpg   data.h Splash     104 0 0 254
	data2c u\preview.jpg +data.h PreviewPic 104 0 0 254
	data2c u\header.avi  +data.h AVIHeader    4 0 0 254
dlt.obj: dlt.c dlt.h matrix.h
	cl /c /Za $(CL_OPTS) dlt.c
draw.obj: draw.c draw.h orb.h ov.h ovrc.h preview.h
	cl /c $(CL_OPTS) draw.c
file.obj: file.c common.c file.h draw.h orb.h ov.h ovrc.h
	cl /c $(CL_OPTS) file.c
matrix.obj: matrix.c matrix.h
	cl /c /Za $(CL_OPTS) matrix.c
orbansi.obj: orb.c orb.h
	cl /c /Za /DANSIC=true $(CL_OPTS) /Foorbansi.obj orb.c
orb.obj: orb.c orb.h
	cl /c $(CL_OPTS) orb.c
ov.obj: ov.c file.h orb.h ov.h ovrc.h data.h draw.h
	cl /c $(CL_OPTS) ov.c
ov.res: ov.rc ovrc.h u\tools.bmp
	rc $(RC_OPTS) ov.rc


clean:
	del *.obj
	del *.res
	del data.h


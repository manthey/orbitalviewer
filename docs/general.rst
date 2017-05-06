************
Introduction
************

Orbital Viewer can display almost any atomic or molecular electron probability function.  For atoms, this plotted function is the hydrogenic solution of Schrodinger's Equation.  For molecules, the linear combination of atomic orbitals (LCAO) method is used.

There are three broad categories of what can be displayed: surfaces of constant probability, probability distributions, and asymptotes (or locations of zero probability).  Each of these show something different.  The probability distribution is the closest to what might be measured by actual experiments.

There are also three general types of display: point plots (only for probability distributions and asymptotes), polygon plots (only for surfaces of constant probability and asymptotes), and raytracing.  Generally, points are the fastest to manipulate (rotate, zoom, etc.), polygons are the fastest to coarsely compute, and raytracing is the best looking.

Display options include:

- *Atoms*: any energy state with *n* <= 30.  Atomic weight and number can be specified.

- *Molecules*: any number of atoms

- *Cutaways*: remove part of the orbital to better see the internal structure

- *Light sources*: any number of lights with varying ambiance and intensity

- *Opacity*: make positive phases, negative phase, and asymptotes translucent

- *Refraction*: change the index of refraction

- *Stereo*: see orbitals in 3D.  Many different stereo algorithms are available

- *Sequences*: record a series of images as the orbital rotates or otherwise changes

- *Output*: VRML and standard graphics files

*******************
System Requirements
*******************

Orbital Viewer requires:

- Windows 3.x or newer.  Windows 3.x requires WIN32s.

- A 386 or better processor

- 2 Mb of RAM

- 1 Mb of disk space

The following is recommended:

- The fastest processor available

- 24-bit graphics

********
Versions
********

Although this version of Orbital Viewer has version numbers beginning at 1.00, the program dates back many years.  The following is a rough history:

- 1986-87:  A primitive point probability program was written for the Apple II series of computers.  This included a fixed animation (no user control) spinning in real-time on the Apple II plus (admittedly, it was only about 80 pixels square with a few hundred points).  This version could only handle with orbitals of *n* <= 4.  This was used as a demonstration for a high-school chemistry coarse I was taking at the time.  Programming was done in Basic and assembly language.

- 1988-89:  A version which allowed any orbital with *m* = 0 was written.  This produced a point surface plot (like the polygon mode with points showing).  A series of fixed animations were created for the Apple IIe computer (this time at the full resolution of 280x192).  Programming was exclusively in assembly language.

- 1990-92:  A preliminary polygon version was written for the Apple IIgs using a 3D graphics library the author also created.  This was much, much too slow.  At the same time, a computationally inefficient raytracing version was being run on an IBM mainframe.  Programming was in assembly language, Pascal, and C.

- 1993-96:  A very cryptic, but powerful, command line version of the raytracing method was programmed for Intel machines.  This was used extensively to produce MPEG animations.  Programming was in C with the actual computational aspects in very finely tuned assembly language.  An ANSI-C compatible version was also programmed.

- 1997-98:  Finally a user interface for viewing orbitals.  This program is written in C with some assembly language.  An ANSI-C command line version was made.

- 1999-02:  Only maintainence on the program has been done.  Numerous features have been suggested, and will probably be implemented some day.

See changes.txt for a version history since 1998.

**************
Mouse Commands
**************

Whenever an orbital window is open, the mouse can perform the following actions:

- *Left Button*: mouse movement rotates the orbital.

- *Right Button*: vertical mouse movement zooms the orbital.

- *Middle Button* or *Both Buttons*: mouse movement pans (shifts) the orbital.

These actions also work in the previews in the Cutaway Dialog, the Lighting Dialog, and the Camera Dialog.



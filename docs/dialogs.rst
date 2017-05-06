*******
Dialogs
*******

Asymptote Options Dialog
========================

The asymptotes are the points at which the probability function goes to zero.  This occurs between the positive and negative phases.  For raytracing, only the opacity has relevance.

- **Opacity**: If set to less than 100%, the surface of the asymptote will be translucent.  When using raytracing, this has a uniform opacity.  When in polugon or point mode, this is not a true translucence, but is done by drawing only some of the pixels to produce quick results.  If 0%, no asymptotes are drawn.

- **Wireframe**: If selected, then only the edges of the asymptote polygons are drawn.  The opacity settings are ignored.

- **Points**: If selected, then only the verticies of the asymptote polygons are drawn.  The opacity settings are ignored.

- **Polygon Density**:- When an asymptote is drawn in either the point mode or the polygon mode, it is approximated by a series of polygons (much the way the polygon mode works).  The density is the maximum number of polygon edges across the asymptote radius.  This number must be even and must be at least 6.

  Generally, the density is required to be at least the maximum of :math:`2 (n - l)`, :math:`2 (l - |m|)`, and :math:`2 |m|`, or features will be missing.  If the asymptote is not complete, try increasing the density.

- **Auto**: Automatically sets the polygon density to "good" values.  The density is a best guess, and will almost never be too high.  It may, however, be too low.

Camera Options Dialog
=====================

When the orbital is rotated using the mouse or menu commands, it's coordinate system actually remains the same.  Rather, the camera is moved.

The axes figure shows the orientation, position, and size of the current camera view.  It shows the three cardinal axes.  It can be rotated using the left mouse button, zoomed using the right mouse button, and panned using the middle mouse button (or both right and left mouse buttons together).

- **Fixed Image Size**: If this option is not selected, the image on the screen is stored internally at the same resolution as the screen.  This means that the screen resolution is used for saving, copying, and raytracing.  If the window is resized, the size of the image changes.  When using raytracing, the computer must restart the computations.

  When fixed size is turned on and a positive **width** and **height** are specified, the orbital is internally drawn at the specified resolution.  This is the resolution saved to a graphics file or copied to the clipboard.  Raytracing draws the orbital at the specified resolution, then scales the result to fit the visible window.  This allows images larger than the screen to be produced (limited only by available memory), and prevents resizing a window from affecting processing.

- **Current Radius**: This is current length of the axes drawn to the right.  It is also the size of the orbital that is considered when Reset Position or Reframe Orbital are selected.  If asymptotes or probability density plots are used, this radius is located at point where the probability is 1/1000 of the maximum probability present in the orbital.  If only probability surfaces are used, this radius is at the point which encloses this surface.  Note that the radius will be artificially large for molecules, as the computations to determine surface extents and maximum probability are too time-consuming to do otherwise.

- **theta**, **phi**, and **psi** angles: This is the orientation of the camera.  At the 0, 0, 0 orientation, the coordinate system has positive *z* out of the screen, positive *x* to the right, and positive *y* pointing upward.  A rotation in *theta* rotates about the *z* axis.  A rotation in *phi* rotates about the transformed *x* axis (transformed by *theta*).  A rotation in *psi* rotates about the transformed *y* axis (transformed by *theta* and *phi*).

- **x**, **y**, and **z** position: This is the camera's position with respect to the center of the orbital, in the camera's rotation system.  Zero *x* and *y* values and a positive *z* value always points directly at the center of the orbital, regardless of the camera orientation.

Colors Dialog
=============

The specified colors are used for all drawing modes (points, polygons, and raytracing).  On 8-bit displays, the actual specified colors may not be available, and an approximation is used.  In stereo display modes, the colors may be ignored.

- **Background**: This is the color used in areas where no orbital is visible.  It is not used on stereograms.

- **Positive Phase**: This is the color used for the positive phase of the orbital when the light intensity is of value 1.  The positive phase is the part of the orbital where the value of *Psi* (:math:`\Psi`) is greater than zero (:math:`\Psi^2` is plotted).  It is not used in red-blue, stereograms, or chromadepth stereo displays.

- **Negative Phase**: This is the color used for the negative phase of the orbital when the light intensity is of value 1.  The negative phase is the part of the orbital where the value of Psi is less than zero (:math:`\Psi^2` is plotted).  It is not used in red-blue, stereograms, or chromadepth stereo displays.

- **Asymptote Color**: This is the color used for the orbital asymptote when the light intensity is of value 1. It is not used in red-blue, stereograms, or chromadepth stereo displays.

- **Preview**: This is the color of the figures shown in the Lighting Dialog and in the Cutaway Dialog.

- **Apply to all windows**: If this is selected, changes in color will affect all open windows.  Otherwise, changes will only affect the active window.

- **Change preferences**: If this is selected, all new orbitals will have the new colors; otherwise this change is a one-time event.

- **Reset**: Sets the colors back to the default values.

Customize Toolbar Dialog
========================

This is a standard Windows dialog; see Windows documentation for details.

Cutaway Options Dialog
======================

A cutaway allows the interior of an orbital to be seen.  For polygons, this results in a very coarse approximation unless the cutaway is exactly aligned with the coordinate system axes.

The figure shows what part of the orbital will be removed.  It is drawn in the default camera orientation (not the current camera orientation).  The sphere has the same radius as the current orbital (this is displayed in the Camera Dialog).

The cutaway can be rotated using the left mouse button, moved front and back using the right mouse button, and panned using the middle mouse button (or both right and left mouse buttons together).

- **Cutaway Type**: The cutaway type determines how much of the orbital is removed.

  - **None** turns off cutaways.

  - **Plane** removes everything with a positive *z* value.

  - **Corner** removes everything that has a positive *z* and a positive *y* value.

  - **Wedge** removes everything which has positive *x*, *y*, and *z* values.

- **Surface**: (Applies only to raytracing)  This determines if a part of the orbital which is cut by the cutaway will have a surface.  Without the surface, the "interior" of a surface plot can be seen.  With a surface, the orbital will appear as a solid object with a piece removed.

- **Invert**: This reverses what is kept versus what is removed.  The coordinates for this are determined by the cutaway orientation and position.  This has no effect if the None option is selected.

- **alpha**, **beta**, and **gamma** angles: This is the orientation of the cutaway.  At the 0, 0, 0 orientation, the coordinate system has positive *z* out of the screen, positive *x* to the right, and positive *y* pointing upward.  A rotation in *alpha* rotates about the *z* axis.  A rotation in *beta* rotates about the transformed *x* axis (transformed by *alpha*).  A rotation in *gamma* rotates about the transformed *z* axis (transformed by *alpha* and *beta*).

- **x**, **y**, and **z** position: This is the cutaway's position with respect to the center of the orbital, in the orbital's rotation system.

Light Source Dialog
===================

Lighting determines how a polygon or raytraced orbital looks.  The figure shows a preview of the currently selected light's position, intensity, and ambience.  It is always shaded with a planar light source, not a point source.  The light's position can be rotated using the left mouse button, and can be moved closer or further from the origin with the right mouse button.

- **Light**: This is the currently displayed light source.  Any number of light sources may be specified, limited only to available memory.  Additional light sources will slow down orbital computation for all but point displays.

- **Add**: Adds an additional light source, exactly duplicating the currently displayed light source.  In no light sources are currently defined, the light source is the default upper-left light source.

- **Delete**: Discards the current light source.

- **x**, **y**, and **z** position: This is the location of the current light source in the orbital coordinate system  If the **Rotate with viewpoint** box is not checked, then the light source will be always appear to be from the same location, regardless of how the camera is rotated.  The light source is never moved aside from this rotation.

  For polygon displays, all light source are planar light sources.  Their position only determines the direction.  As such, a light source at 0, 0, 0 will not properly illuminate a polygon set.

  For raytraced displays, all light sources are point light sources.  If the light source is sufficiently far away, it is effectively a planar light source.

- **Rotate with viewpoint**: If unselected, the light source will remain in a visually fixed location while the orbital is rotated.  If selected, then the light source remains in the orbital coordinate system, which rotates the same way the orbital rotates.

- **Intensity**: This is the amount of light put out by the current light source.  An intensity of 1 which is normal to a surface will produce exactly the color specified in the Colors Dialog.  Intensities may exceed 1.  Typically, the sum of the intensities of all light sources should be from 1 to 1.1.  If an orbital is too dark, the intensity can be increased to brighten it.

- **Ambiance**: This determines how sharp the shadows are.  It only applies to raytracing.  An ambiance of 1 produces no shadows.  An ambiance of 0 produces completely black shadows.

Orbital Dialog
==============

This dialog determines exactly what atom or molecule is drawn.  Different elements can be specified by explicitly giving the number of protons and the mass of the atom.  Molecules are produced using a linear combination of atomic orbitals (LCAO).  To reverse the positive and negative phases of an atom, change the factor to -1.

- **Atom**: This is the currently displayed atom.  Any positive number of atoms may be specified, limited only to available memory.  Calculation time is directly affected by the number of atoms used.

- **Add**: Adds an additional atom, exactly duplicating the currently displayed atom.

- **Delete**: Discards the current atom.

- **n**: Principal quantum number.  This is the most significant quantum number.  It is a positive integer.  Due to precision limitations of 64-bit floating point numbers, *n* is restricted to values from 1 to 30.

- **l**: Orbital quantum number.  This is a non-negative integer that is less than *n*.  It is often referred to by a letter (for historical reasons), as follows:

  0 - s;  1 - p;  2 - d;  3 - f;  4 - g;  5 - h;  6 - i;  7 - k;  8 - l;  9 - m;  10 - n;  11 - o;  12 - q;  13 - r;  14 - t;  15 - u;  16 - v;  17 - w;  18 - x;  19 - y;  20 - z

- **m**: Angular momentum quantum number.  This is an integer from *-l* to *+l*.

- **Protons**: Number of protons in the atom.  This affects the scale of the orbital.

- **Atomic Mass**: The mass of the atom affects the radius of the orbital.  It is dependent on the number of protons and neutrons.  As an approximation, it is (number of protons) + (number of neutrons) amu.

- **Factor**: This is a multiplicative factor for the orbital probability, Psi (:math:`\Psi`).  It is primarily of use in constructing molecules.  For example, the orbital for H2O contains two hydrogen atoms, one which has a positive phase and one which has a negative.  For the first, the factor is 1, and for the second it is -1.

- **alpha**, **beta**, and **gamma** angles: This is the orientation of the atom.  At the :math:`0, 0, 0` orientation, the coordinate system has positive *z* out of the screen, positive *x* to the right, and positive *y* pointing upward.  A rotation in *alpha* rotates about the *z* axis.  A rotation in *beta* rotates about the transformed *x* axis (transformed by *alpha*).  A rotation in *gamma* rotates about the transformed *z* axis (transformed by *alpha* and *beta*).

- **x**, **y**, and **z** position: This is the location of the center of the atom.  It is primarily useful in specifying molecules, since each atom must be at a different location.

Play Sequence Dialog
====================

A sequence gradually transforms the orbital.  This can include rotation and position, lighting, probability values, opacity, color, atomic values, cutaway location, and more.

- **Frame**: This is the frame at which point the sequence is exactly like the orbital file specified to the right.  These points are used to interpolate all of the other frames.

  If Bezier interpolation is used, only the first and last frames will exactly match the selected orbital files.  The other frames are construction frames.
  
  The check box shows which of the four sequence specification orbitals are in use.  At least two orbitals must be used to define a sequence.  If two are used, all values are scaled linearly between the two sequences.  If three are used, values are interpolated using either a bezier or quadratic function, and if four are specified, values are interpolated using a bezier or cubic function.

  Clicking on an empty box will allow the orbital to be selected from a file.  To keep an existing orbital, cancel the file dialog.

  Files can be in either .ORB or .OV format.

  The frame number must be an integer.

- **Incremental positions**:  If this box is turned off, all values, including camera position and orientation are determined by either a linear, quadratic, cubic, or Bezier fit depending on the number of orbitals used to define the sequence.  Since angles are specified on a [0, 2pi) range, using interpolation does not always produce the desired results.

  When the box is checked, the first orbital is used as a starting position, and the difference between the first and second orbital is used as an increment for the camera angle and position on a per frame basis.  Even if the second orbital is at an earlier frame, the difference is the increment per positive frame.

- **Bezier interpolation**: This box only applies when three or four orbital specification files are used to define the sequence.

  When this box is not checked, each frame is determined by a quadratic or cubic interpolation.  The interpolated values are guaranteed to match the specified values at the specified times.  For example, the cubic function uses the formula
  
  :math:`v = A t^3 + B t^2 + C t + D`
  
  where :math:`v` is the computed value, :math:`A`, :math:`B`, :math:`C`, and :math:`D` are the coefficients of the cubic function, and :math:`t` is the current time.

  When the Bezier box is checked, a Bezier spline is used to determine the interpolated values.  The interpolated values are guaranteed to match the specified values at the specified times for the first and last specifications only.  The slope of the interpolation function at the first and last points is the same as the slope between the first and second points and the second to last and last points, respectively.  The Bezier interpolation satisfies the equations

  :math:`t = A_t p^3 + B_t p^2 + C_t p + D_t`

  :math:`v = A_v p^3 + B_v p^2 + C_v P + D_v`

  where :math:`v` is the computed value, :math:`p` is a parametric value, :math:`t` is the current time, :math:`A_t`, :math:`B_t`, :math:`C_t`, and :math:`D_t` are coefficients based on the times of the specification files, and :math:`A_v`, :math:`B_v`, :math:`C_v`, and :math:`D_v` are coefficients based on the values of the specification files.


- **Base save file name**: After a sequence frame has been computed, it can optionally be saved.  If the base is set, the frame is saved as a graphic file with the name ``base number.type``, where ``base`` is the specified base name, ``number`` is the frame number, and ``type`` is the file type.  The base name should be chosen such that it produces valid filenames in the current system using the specified frame numbers.

  For .AVI files, the entire sequence is saved to a file with the name ``base.AVI``.

- **Browse**: This shows a file dialog so that the base file name can be graphically selected.

- **File format**: Sequence files can be saved as any of these graphics formats:

  - .PPM - Portable Pixel Map Files - save a 24-bit per pixel uncompressed graphic.  This is a format popular on Unix machines.

  - .TIF - Tagged Image File Format - save a 24-bit per pixel run-length-encoded compressed graphic.  This is a very widespread format that can be read on almost any platform.

  - .BMP - Windows Bitmap Files - save a 24-bit per pixel uncompressed graphic.  This is compatible with most Windows programs.

  - .AVI - Audio Video Interleaved Files - save an AVI animation file.  This is an uncompressed file (use Compress AVI to compress it) contain every frame.  Frames are appended to the file, never destroying anything already in the file.

- **Frames per second**: This is the rate at which AVI files will be played back.  It only applies to AVI files.

- **Start frame** and **End frame**: The sequence will begin with the specified start frame, and advance or retreat one frame at a time until the end frame is reached.  The orbital for each frame is determined by interpolating between the specification files.  If the end frame is less than the start frame, then the sequence is run "backwards", counting down.

- **Done**: Closes the dialog, saving changes but not starting the sequence.

- **Play**: Closes the dialog, storing the changes made, and starts playing the sequence.  Each frame is generated and saved, then the window advances to the next frame.

Point Options Dialog
====================

The point display is a probability density plot.  This means that there are more points in areas with greater probabilities.

- **Number of Points**: This is the number of points which will be generated.  Depending on the exact orbital, this may be very quick or very slow.  All points are guaranteed to be within the radius where the orbital's probability is 1/1000th of the maximum probability.

- **Asymptote Options**: Show the Asymptote Options Dialog

Polygon Options Dialog
======================

Polygons are an approximation of a surface of constant probability.  They are fast to draw, but are not as good looking as raytracing.

- **Psi^2 probability**: This is the value for which the constant probability is drawn, expressed at the base-10 exponent.  A more negative number results in a probability surface which has a larger radius.  If the probability is too high, no surfaces of the orbital may be visible.

- **Auto**: Automatically sets the probability and polygon density to "good" values.  The density is a best guess, and will almost never be too high.  It may, however, be too low.

- **Asymptote Options**: Show the Asymptote Options Dialog

- **Positive phase opacity**:  If set to less than 100%, the surface of the polygon will be functionally translucent by not drawing all pixels.  This is not a true translucence like raytracing, but is done to still produce quick results.  If 0%, the positive phase surface will not be drawn.

- **Negative phase opacity**:  If set to less than 100%, the surface of the polygon will be functionally translucent by not drawing all pixels.  This is not a true translucence like raytracing, but is done to still produce quick results.

- **Density**: A probability surface is approximated by a series of polygons.  The density is the maximum number of initial polygon edges across the orbital radius.  This number must be even and must be at least 6.  After the initial polygons are computed, they are divided into smaller polygons based on the Refine setting.

  Generally, the density is required to be at least the maximum of :math:`2 (n - l)`, :math:`2 (l - |m|)`, and :math:`2 |m|`, or features will be missing.  If the orbital is not complete, try increasing the density.

- **Refine**: The probability surface is initially approximated by a series of polygons, the quantity of which is determined by the density.  After the initial polygons are computed, they are divided into smaller polygons based on the Refine setting.  A setting of 0 will prevent any refinement from being performed.  A refinement of :math:`x` will ensure that no edge is longer than :math:`1/x` of the original density specification.  The larger the refinement number, the more polygons will be generated.

- **Wireframe**: If selected, then only the edges of the orbital polygons are drawn.  The opacity settings are ignored.

- **Points**: If selected, then only the verticies of the orbital polygons are drawn.  The opacity settings are ignored.

Preferences Dialog
==================

- **Show Error Messages**: These are messages informing of problems such as insufficient memory or a problem reading or writing a file.

- **Show Warning Messages**: These are dialogs which confirm that a certain action should be taken, such as deleting an atom or a light source, overwriting an existing file, or reseting the preferences.

- **Show Toolbar**: This is the line of controls immediately below the menu.

- **Show Tool Tips**: These are the small rectangles of text which identify the controls on the toolbar.

- **Show Status Line**: This is the bottom line on the main window showing the percent completion of an orbital, the frame or atom, and informational text about menu items.

- **Show Start Screen**: This is the window which is displayed for five seconds when the program starts.  When it is shown, it can be removed sooner by clicking on it or using a menu.

- **VRML File Options**: Some VRML require a particular point or polygon color style or the points and polygons will not appear.  The default is to use ambient and diffuse color for both point and polygons, and to use emissive color for the points (emissive color is never used for polygons).  This works in many, but not all VRML viewers.

- **Point Size**: In the point drawing style, and for asymptotes and polygons drawn as points, this is the size of the point on the screen.  The point will be this many pixels square.  Larger sizes improve visibility, but produce a chunky look.

- **Menu / Toolbar Step Sizes**: When the menu items are used to rotate, pan, or zoom the orbital, these controls specify by how much.

- **Reset All**: This restores all colors and preferences to the default setting.  It also disregards the DEFAULT.ORB file, making new windows with the original default orbital.

Raytrace Options Dialog
=======================

Raytracing is the highest quality method of drawing an orbital, but it is also the slowest.  It is the only method which supports shadows, index of refraction, and true opacity.

- **Psi^2 probability**: This is the value for which the constant probability is drawn, expressed at the base-10 exponent.  A more negative number results in a probability surface which has a larger radius.  If the probability is too high, no surfaces of the orbital may be visible.  This does not affect drawings based on probability opacity.

- **Auto**: Automatically select surface probability.  This is the computer's best guess at what will make an interesting figure.

- **Probability opacity per step**: This determines the opacity of a particular pixel based on the probability of the orbital along the ray from the camera point through that pixel.  The higher the probability, the more opaque the image.  This is computed by taking a number of steps through the orbital (specified in the Number of Steps parameter), and evaluating the probability at each point.  A point with the maximum probability in the orbital will have the specified opacity within that step.  Typically small values work well, such as between 5 % and 20 %.

- **Surface opacity**: This is the opacity of a surface of constant probability.  For positive and negative phases, the probability is specified above, whereas for asymptotes the probability is zero.  For positive and negative phases, the opacity is constant regardless of view angle.  For asymptotes the opacity increases as the view becomes oblique.

- **Interior opacity per step**: If a point is within the surface of constant probability, than it has this opacity.  This is computed by taking a number of steps through the orbital (specified in the Number of Steps parameter), and evaluating the probability at each point.  Any point above the specified probability will have this opacity.  Numbers are specified in percent of percent (1%% equals 0.0001).

- **Index of refraction**: A non-unity index of refraction creates a "lensing" effect, as if the orbital were made of glass.  The index of refraction is relative to that of free space.  This means that to model a glass orbital in air, an index of refraction of 1.54 (crown glass) might be specified, whereas to model an air orbital in glass, the index of refraction would be :math:`1/1.54=0.65`.  Note that shadow computations are not technically perfect when a non-unity index of refraction is specified.

- **Link positive and negative phases**: When this box is checked, any change to the positive phase will result in an equal change to the negative phase, and vice versa.

- **Number of steps**: For any plots other than 100% opacity surface plots, each pixel is computed by taking a series of small steps through its interior and determining the cumulative effect.  This determines the number of steps through the radius of the orbital.  Twice this many will be taken at the center, whereas almost no steps are taken at the edge.

- **Auto Brightness**: When turned on, this will draw the orbital in a very coarse manner (16x16 pixel blocks), and, based on this picture will adjust the brightness to a "good" value.  This necessitates redrawing the coarse pixels.

  To adjust brightness manually, change the intensity of the light source.

- **Antialias**: An orbital is raytraced using a multipass technique.  First, 16x16 pixel blocks are drawn, then 4x4 pixel blocks, then individual pixels.  If antialias is turned on, one more pass is performed, where pixels at color transitions are analyzed by oversampling (measuring the pixels at non-integer locations).  The average of the oversampling is then drawn.  This has the effect of smoothing the edges of the orbital.

- **Coarse Render**: An orbital is raytraced using a multipass technique.  First, 16x16 pixel blocks are drawn, then 4x4 pixel blocks, then individual pixels.  If coarse rendering is turned on, then the individual pixels may be approximated from the 4x4 blocks.  This is only done in regions where there is relatively smooth color transition between the 4x4 blocks.  Coarse rendering can substantially speed up drawing, but may blur or miss small features.

Rendering Method Dialog
=======================

Point mode draws a point probability density plot.

Polygon mode draws a probablity surface approximated by polygons.

Raytracing produces the best drawings, but is generally the slowest.  It has more options than point or polygon modes, including index of refraction and interior opacity.

- **Use Quick Rendering**:- When the orbital is drawn, one of the three rendering techniques can be used.  When quick rendering is turned on, any time the orbital is moved or rotated, the quick method will be used to draw the orbital.  When the quick method has finished, the precise method will then be used to redraw it.

  This allows, for example, a polygon drawing method be used for spinning the orbital in real time, while drawing a higher precision raytracing when the orbital s not moving.

Stereo Options Dialog
=====================

A stereo display allows the orbital to be seen three-dimensionally.  This is always at the expense of either resolution or color.  

The accompanying graphic shows an example of the stereo mode.  With the exception of stereograms, these are all pregenerated images in the actual stereo mode.  Stereograms only show random dots (no image within it) or the selected image file.

- **Monoscopic**: The standard one image, non-stereo method of drawing.

- **Stereoscope**: A stereoscope is a device with mirrors to allows two images from different vantage points to be seen easily.  Stereoscopes can be cheaply constructed from easily obtainable materials (see the elsewhere inthe documentation for instructions), or can be bought with high-quality optics.  The two views are present side by side.

- **Interlaced**: LCD shutter glasses often work by showing alternate scan lines of the monitor to each eye.  The **Swap** option switches which line is drawn on the left.

- **Red-blue** (Anaglyph): The old classic 3D glasses use one red lens and one blue lens.  This allows either the right or the left to contain either color.  Color information is lost in this mode.

- **Stereogram**: The Single Image Random Dot Stereogram (SIRDS) has the benefit that it does not require special optics to view.  This can either use random pixels or can use a source image.  The preview shown to the right is NOT an actual stereo image.  Color and intensity information is lost in this mode.

- **Overlay**: For those people blessed with excellent stereo vision, this method is a more attractive alternative to stereograms.

- **Chromadepth**: Chromadepth glasses are similar to anaglyph (red-blue) glasses, but they use diffraction to create depth.  The depth is therefore based on the color, with red being closest.  Original color information is lost in this mode.

- **Interocular distance**: This is the separation between the eyes in pixels.  For red-blue, overlay, stereogram, and interlaced modes, this is the actual physical distance between the eyes (which is dependent both on the viewer and on the monitor used).  For stereoscopes, this is the distance between the eyes, plus the mirror separation distance.

- **Actual separation**: For stereoscope, red-blue, overlay, and interlaced modes, this is the "physical" separation between the eyes in the orbital coordinate system.  Generally, this is a very small number since the orbitals are very small.    More separation will give a greater perception of depth.  If **Auto** is turned on, the separation is selected automatically to a "good" value.

- **Perspective factor**: This number determines the effective focal length of the camera used in viewing the orbital.  A larger number is a larger focal length, which results in less perspective distortion.  Typically, values between 5 and 25 work well.  Larger numbers have little noticeable effect, while smaller numbers have pronounced keystone effects.

- **Stereogram Image**: Stereograms can either be based on random dots or on a source image.  If this box is not checked, it is based on random dots.  Checking the box will present a file dialog along any graphic file in the following formats to be selected: .BMP (Windows bitmap), .CUR (Windows cursor), .GROB (HP calculator), .ICO (Windows icon), .JPG (jpeg), .PCX (Paintbrush), .PPM (portable pixel map), .TGA (targa), .TIF (tagged image file format).  The image is displayed in the preview area.

  Some additional formats may be usable if Python with the Pillow module are installed such that Python can be run from the system path.


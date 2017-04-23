# Orbital Viewer

[![Build status](https://ci.appveyor.com/api/projects/status/sxe3531fs0qjh5g5/branch/master?svg=true)](https://ci.appveyor.com/project/manthey/orbital-viewer/branch/master)

A program for drawing orbitals.

## Atomic Orbitals

Electron orbitals are the probability distribution of an electron in a atom or molecule.

The electron orbitals presented here represent a volume of space within which an electron would have a certain probability of being based on particular energy states and atoms.  For example, in a simple lowest-energy state hydrogen atom, the electrons are most likely to be found within a sphere around the nucleus of an atom.  In a higher energy state, the shapes become lobes and rings, due to the interaction of the quantum effects between the different atomic particles.  In addition to technical merits, they make pretty pictures.

The shape of the orbital depends on many factors.  The most important are the quantum numbers associated with the particular energy state.  These are _n_, the principal quantum number, _l_, the orbital quantum number, and _m_, the angular momentum quantum number.

## Historical Note

This program was started in either 1986 or 1987, initially as a simple black-and-white dot-density plot on an Apple \]\[ computer that only rendered orbitals where `m = 0`.  A few years later there was an Apple \]\[gs version that could render orbitals where `n < 6`.  This version first came out in 1998.

The coding style is not what I would chose today.

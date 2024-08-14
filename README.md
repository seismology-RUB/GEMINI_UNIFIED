# GEMINI_UNIFIED
GEMINI (Green's function of the Earth by MInors INtegration) is a code that computes Green functions for a spherically symmetric, elastic and isotropic earth model. The method is described in Friederich and Dalkolmo (1995).

GEMINI_UNIFIED is new version of GEMINI that offers capabilities for both global and shallow regional computations. In the former case, it performs a spherical harmonics expansion of the displacement field, in the latter case, an approximation of the spherical harmonics by Bessel function is applied. With this version it is also possible to compute synthetic seismograms at receivers inside the earth. Moreover, by using the reciprocity theorem, Green functions at the surface for many sources at different depths can be obtained in a single run.

GEMINI works in a two-phase procedure: first, radially varying spherical harmonics expansion coefficients of the Fourier transform of the displacement field are computed which depend on the earth model and the source depth. Based on these coefficients, synthetic seismograms for any moment tensor or single force, epicentral distance and receiver depth can be rapidly calculated by spherical harmionics summation and inverse Fourier transform.

GEMINI_UNIFIED runs reliably for global computations up to 0.2 Hz and has been tested for global body wave computations up to frequencies of 0.5 Hz. It is not unconditionally stable and may fail for certain combinations of frequencies and harmonic degrees. Much higher frequencies can be obtained for shallow, short-distance regional computations. However, this aspect of the code still needs some thorough testing.

GEMINI_UNIFIED is written in modern Fortran, fully modularized, fully parallelized using MPI and produces convenient HDF output files. Most executables are controlled by one parameter file (``parfile_gemini``). In addition, users need to provide a JSON earth model file (see folder ``earthmodels``) and an event file (see the templates folder for examples). To compute synthetic seismograms, a station file is required (see ``ASKI_stations_S`` in ``templates``).

For computation of the radially dependent spherical harmonics expansion coefficients in the frequency domain (these are called here Green frequency-wavenumber (FK) spectra), you may want to run the program ``computeGreenFKSpectraForSynthetics`` which computes Green FK spectra for a fixed receiver at some depth and sources at specified radii (EXTERNAL_RADIAL_NODES). Seismograms are then obtained by first running ``computeStationFilter`` (in case of a new station file) and then ``computeSyntheticSeismograms``.

You need a Fortran compiler that supports at least the Fortran 2003 standard, an MPI installation and HDF5 with enabled parallel option. A Makefile is offered for compilation. There, you should modify the paths to the diverse libraries to values appropriate for your system.

------------
Contributors
------------
Wolfgang Friederich, Jörg Dalkolmo, Thomas Möller, Kasper D. Fischer

--------------
Latest release
--------------

[![DOI](https://zenodo.org/badge/841969358.svg)](https://zenodo.org/doi/10.5281/zenodo.13321290)

----------
References
----------
Friederich, W. and Dalkolmo, J. (1995), Complete synthetic seismograms for a spherically symmetric earth by a numerical computation of Green's function in the frequency domain, Geophys. J. Int., 537--550, 122, https://doi.org/10.1111/j.1365-246X.1995.tb07012.x

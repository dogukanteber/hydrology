# Hydrology solvers for OpenFOAM

_The OpenFOAM user guide can be found on https://www.openfoam.com/documentation/user-guide/_

- _Author_:  Laurent Orgogozo
- _Contact_: Laurent.Orgogozo (@) get (.) omp (.) eu


## Introduction

### Source code

Licence: GPL-3.0-or-later


## Installation

0. You need an installation of OpenFOAM on your computer and
   compiler/make tools.
   - Tested with [OpenFOAM-v1912](https://www.openfoam.com),
     but should also work with newer version.

1. Clone the repository or download a tar-archive to a user-writeable
   location

2. Run the top-level `./Allwmake` script, which will compile the 
   hydrology library and the solvers.


The hydrology library primarily supplies boundary
conditions such as _rainFlux_ and _noRainFlux_ as well as some
additional functs for seasonal time variation of boundary conditions:
(seasonal, thawing, cosine1, sine1, square1).


## Testing

To test the solver permaFoam (see Orgogozo et al., Permafrost and
Periglacial Processes 2019 for details), there are a few demonstration
cases.
You can copy any of these cases to the "run/" directory of your
OpenFOAM workspace (create it if it does not already exist), or simply
run inplace. The usual `./Allrun` and `./Allclean` commands are used to
run the case (create mesh, run solver, postProcess).

The post-processing output files will contain averaged pressure heads,
temperature, total water content, liquid water content, and ice
content (in postProcessing/internal/0/volFieldValue.dat), in terms of
water fluxes at the top boundary (in
postProcessing/top/0/surfaceFieldValue.dat), and, of course, in terms
of all the output fields at each time of writing of the results, which
will be contained in the time directories (named along the associated
time coordinate).

_Note_:
The directory `expected_postProcessing` contains the expected results
of the postprocessing operations (permaFoam -postProcess), so that
they can be compared to the ones that the execution of the commandsRun
procedure provides, for validating your installation.


The full fields can be visualized in paraview, with the regular
ParaView/VTK OpenFOAM reader and/or with the so-called _paraFoam_
reader module.


### tutorials/permaFoam/demoCase

This case is a simple setup.


### tutorials/permaFoam/demoCase_sinusoidalClimatForcings

This case illustrates the handling of seasonal variabilities in the
boundary conditions.


## Further References

The system of equations, the input fields and the post-processing
information are further described in the [peraFoam wiki][wiki-permaFoam]
and the [RichardsFoam wiki][wiki-richardsFoam].


---

[wiki-permaFoam]:     https://develop.openfoam.com/Community/hydrology/-/wikis/permafoam
[wiki-richardsFoam]:  https://develop.openfoam.com/Community/hydrology/-/wikis/RichardsFoam

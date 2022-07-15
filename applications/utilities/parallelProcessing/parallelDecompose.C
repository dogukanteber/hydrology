/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    parallelDecompose

Group
    grpParallelUtilities

Description
    Automatically decomposes a mesh and fields of a case in parallel
    execution of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char* argv[]) {

    argList::addNote
    (
        "Decompose a mesh and fields of a case in parallel execution"
    );
    argList::noParallel();
    argList::addOption
    (
        "decomposeParDict",
        "file",
        "Use specified file for decomposePar dictionary"
    );
    
    #include "setRootCase.H"

    #include "createTime.H"

    // get custom decomposeParDict location
    fileName decompDictFile(args.get<fileName>("decomposeParDict", ""));
    if (!decompDictFile.empty() && !decompDictFile.isAbsolute())
    {
        decompDictFile = runTime.globalPath()/decompDictFile;
    }

    return 0;
}

// ************************************************************************* //
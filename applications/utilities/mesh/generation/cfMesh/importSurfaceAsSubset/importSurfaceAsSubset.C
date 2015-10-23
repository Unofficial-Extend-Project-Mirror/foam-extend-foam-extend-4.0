/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Finds feature edges and corners of a triangulated surface

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fileName.H"
#include "triSurf.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "demandDrivenData.H"
#include "triSurfaceImportSurfaceAsSubset.H"
#include <cstdlib>
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("master surface file");
    argList::validArgs.append("import surface file");

    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName importFileName(args.args()[2]);

    triSurf originalSurface(inFileName);

    triSurf importedSurface(importFileName);

    triSurfaceImportSurfaceAsSubset importSurf(originalSurface);

    importSurf.addSurfaceAsSubset(importedSurface, importFileName.lessExt());

    if( inFileName.ext() == "fms" )
    {
        originalSurface.writeSurface(inFileName);
    }
    else
    {
        fileName newName = inFileName.lessExt();
        newName.append(".fms");
        Warning << "Writting surface as " << newName
            << " to preserve the subset!!" << endl;

        originalSurface.writeSurface(newName);
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

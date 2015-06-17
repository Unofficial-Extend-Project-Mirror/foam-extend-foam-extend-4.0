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
    Reads the surface mesh, remove the selected facets
    and writes the modified mesh into a new file

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "triSurf.H"
#include "triSurfaceRemoveFacets.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");
    argList::validArgs.append("Subset/patch name");

    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);
    word patchName(args.args()[3]);

    //- read the input surface
    triSurf origSurf(inFileName);

    //- remove the selected facets
    triSurfaceRemoveFacets rsf(origSurf);

    rsf.selectFacetsInPatch(patchName);

    rsf.removeFacets();

    //- write the modified surface mesh
    origSurf.writeSurface(outFileName);

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //

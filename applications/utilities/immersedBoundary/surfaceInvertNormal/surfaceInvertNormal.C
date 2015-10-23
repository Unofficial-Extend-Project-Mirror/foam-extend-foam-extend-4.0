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
    Write surface and normal vector

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurface.H"
#include "IFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");
    argList args(argc, argv);

    fileName surfFileName(args.additionalArgs()[0]);

    Info<< "Reading surf from " << surfFileName << " ..." << endl;

    fileName outFileName(args.additionalArgs()[1]);

    Info<< "Writing surf to " << outFileName << " ..." << endl;

    triSurface ts(surfFileName);

    triFaceList invertedFaces(ts.size());

    forAll (ts, tsI)
    {
        invertedFaces[tsI] = ts[tsI].reverseFace();
    }

    triSurface invertedTs
    (
        invertedFaces,
        ts.points()
    );

    invertedTs.write(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

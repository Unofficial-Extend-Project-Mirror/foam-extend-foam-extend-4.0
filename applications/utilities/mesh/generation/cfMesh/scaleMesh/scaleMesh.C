/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Scales the mesh into other units.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMeshGen.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::validArgs.append("scalingFactor");

#   include "setRootCase.H"
#   include "createTime.H"

    const scalar scalingFactor(help::textToScalar(args.args()[1]));

    Info << "Scaling mesh vertices by a factor " << scalingFactor << endl;

    //- read the mesh from disk
    polyMeshGen pmg(runTime);

    Info << "Reading mesh" << endl;
    pmg.read();

    //- scale the points
    pointFieldPMG& pts = pmg.points();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(pts, pointI)
        pts[pointI] *= scalingFactor;

    //- write the mesh back on disk
    Info << "Writting scaled mesh" << endl;
    pmg.write();

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

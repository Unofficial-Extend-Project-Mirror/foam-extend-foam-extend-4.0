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

Application
    deformedTetFemGeom

Description
    Deforms a polyMesh using a displacement field U from the tetFem method
    and a scaling factor supplied as an argument.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "pointFields.H"
#include "tetFemMatrix.H"
#include "tetPointFields.H"
#include "IStringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("scaling factor");

#   include "setRootCase.H"

    scalar scaleFactor(readScalar(IStringStream(args.additionalArgs()[0])()));

#   include "createTime.H"
#   include "createPolyMesh.H"

    tetPolyMesh tetMesh(mesh);

    // Get times list
    instantList Times = runTime.times();

    pointField zeroPoints(mesh.points());

    runTime.setTime(Times[0], 0);

    for (int i = 1; i<Times.size(); i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check U exists
        if (Uheader.headerOk())
        {
            Info<< "    Reading U" << endl;
            tetPointVectorField U(Uheader, tetMesh);

            vectorField Uslice
            (
                vectorField::subField
                (
                    U.internalField(),
                    zeroPoints.size()
                )
            );

            pointField newPoints = zeroPoints + scaleFactor*Uslice;

            mesh.movePoints(newPoints);

            mesh.write();
        }
        else
        {
            Info<< "    No U" << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

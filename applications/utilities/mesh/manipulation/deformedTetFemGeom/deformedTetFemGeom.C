/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
#include "timeSelector.H"
#include "pointFields.H"
#include "tetFemMatrix.H"
#include "tetPointFields.H"
#include "IStringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "addRegionOption.H"

    argList::validArgs.append("scaling factor");
    argList::validOptions.insert("UName", "name");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    scalar scaleFactor(readScalar(IStringStream(args.additionalArgs()[0])()));

    word UName("U");
    args.optionReadIfPresent("UName", UName);
    
    instantList timeDirs = timeSelector::select0(runTime, args);

    tetPolyMesh tetMesh(mesh);

    pointField zeroPoints(mesh.points());

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        polyMesh::readUpdateState state = mesh.readUpdate();

        IOobject Uheader
        (
            UName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check U exists
        if (Uheader.headerOk())
        {
            Info<< "    Reading " << UName << endl;
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
            Info<< "    No " << UName << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

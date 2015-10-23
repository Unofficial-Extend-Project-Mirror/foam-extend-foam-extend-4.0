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
    calcLevelSet

Description
    Calculates and writes the scalar magnitude of the levelSet variable

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("normal", "vector");
#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"

    // Get normal
    vector normal(0, 0, 1);

    if (args.options().found("normal"))
    {
        normal = vector(IStringStream(args.options()["normal"])());

        if (mag(normal) > SMALL)
        {
            normal /= mag(normal);
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Invalid normal given: " << normal
                << abort(FatalError);
        }

    }

    Info << "Normal vector for height field: " << normal << endl;

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        Info<< "    Calculating height" << endl;
        volScalarField height
        (
            IOobject
            (
                "height",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh.C() & normal
        );

        height.write();
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

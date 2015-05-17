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
    stressDifferences

Description
    Calculates and writes the first (N1) and second (N2) stress difference
    for each time.

Author
    Jovani L. Favero.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"

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

        IOobject tauHeader
        (
            "tau",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check tau exists
        if (tauHeader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading tau" << endl;
            volSymmTensorField tau(tauHeader, mesh);


            Info<< "    Calculating N1"<< endl;
            volScalarField N1
            (
                IOobject
                (
                    "N1",
                     runTime.timeName(),
                     mesh,
                     IOobject::NO_READ
                ),
                tau.component(symmTensor::XX) - tau.component(symmTensor::YY)
            );
            N1.write();

            Info<< "    Calculating N2"<< endl;
            volScalarField N2
            (
                IOobject
                (
                    "N2",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                tau.component(symmTensor::YY) - tau.component(symmTensor::ZZ)
            );
            N2.write();
        }
        else
        {
            Info<< "    No tau" << endl;
        }

        Info<< endl;

    } //for time

    Info<< endl;
    Info<< "    End"<< nl << endl;

    return(0);
}


// ************************************************************************* //

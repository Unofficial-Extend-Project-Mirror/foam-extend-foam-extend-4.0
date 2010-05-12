/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    stressSymmComponents

Description
    Calculates and writes the scalar fields of the six components of the
    stress tensor for each time.

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

            for (direction i=0; i<symmTensor::nComponents; i++)
            {
                Info<< "    Calculating tau"
                    << symmTensor::componentNames[i] << endl;

                volScalarField taui
                (
                    IOobject
                    (
                        "tau" + word(symmTensor::componentNames[i]),
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    tau.component(i)
                );
                taui.write();
            }
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

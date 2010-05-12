/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010 H. Jasak. All rights reserved
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
    calcPTot

Description
    Calculate total pressure

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"

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

        IOobject pHeader
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check all fields exists
        if
        (
            pHeader.headerOk()
         && Uheader.headerOk()
        )
        {
            mesh.readUpdate();

            Info<< "    Reading p" << endl;
            volScalarField p(pHeader, mesh);

            Info<< "    Reading U" << endl;
            volVectorField U(Uheader, mesh);

            if (p.dimensions() == dimPressure)
            {
                // Dynamic pressure, read rho
                IOobject rhoHeader
                (
                    "rho",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ
                );

                if (rhoHeader.headerOk())
                {
                    Info<< "    Reading rho" << endl;
                    volScalarField rho(rhoHeader, mesh);

                    volScalarField pTot
                    (
                        IOobject
                        (
                            "pTot",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        p + 0.5*rho*magSqr(U)
                    );

                    pTot.write();
                }
                else
                {
                    Info<< "Not all fields are present.  "
                        << "p is dynamic pressure and rho cannot be found"
                        << endl;
                }
            }
            else
            {
                // Kinematic pressure
                volScalarField pTot
                (
                    IOobject
                    (
                        "pTot",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    p + 0.5*magSqr(U)
                );

                pTot.write();
            }
        }
        else
        {
            Info << "Not all fields are present.  " << endl;

            if (!pHeader.headerOk())
            {
                Info << "pd ";
            }

            if (!Uheader.headerOk())
            {
                Info << "U ";
            }

            Info << "missing." << endl;
        }
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

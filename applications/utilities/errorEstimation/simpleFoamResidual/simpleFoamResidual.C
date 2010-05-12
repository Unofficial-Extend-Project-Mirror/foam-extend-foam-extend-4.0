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
    simpleFoamResidual

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Description
    Calculate and write residual for a simpleFoam momentum equation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RASModel/RASModel.H"

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

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject pHeader
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check U exists
        if (Uheader.headerOk() && pHeader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading U" << endl;
            volVectorField U(Uheader, mesh);

            Info<< "    Reading p" << endl;
            volScalarField p(pHeader, mesh);

#           include "createPhi.H"

            singlePhaseTransportModel laminarTransport(U, phi);

            autoPtr<incompressible::RASModel> turbulence
            (
                incompressible::RASModel::New(U, phi, laminarTransport)
            );

            Info<< "    Calculating uResidual" << endl;
            volScalarField uResidual
            (
                IOobject
                (
                    "uResidual",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                mag(fvc::div(phi, U) + fvc::div(turbulence->R()) + fvc::grad(p))
            );

            Info<< "uResidual max: " << max(uResidual.internalField())
                << " mean: "
                << sum(uResidual.internalField()*mesh.V())/sum(mesh.V()).value()
                << endl;

            uResidual.write();
        }
        else
        {
            Info<< "    No U or p" << endl;
        }

        Info<< endl;
            
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

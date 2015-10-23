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
    simpleFoamResidual

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Description
    Calculate and write residual for a simpleFoam momentum equation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

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
                mag
                (
                    fvc::div(phi, U)
                  + fvc::div(turbulence->R())
                  + fvc::grad(p)
                )
            );

            Info<< "uResidual max: " << max(uResidual.internalField())
                << " mean: "
                << sum(uResidual.internalField()*mesh.V())/
                       sum(mesh.V()).value()
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

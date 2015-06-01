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
    shapeDerivative

Description
    Calculate continuity error on immersed boundary mesh

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    IOobject phiHeader
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject faceIbMaskHeader
    (
        "faceIbMask",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (phiHeader.headerOk() && faceIbMaskHeader.headerOk())
    {
        Info<< "    Reading phi" << endl;
        surfaceScalarField phi(phiHeader, mesh);

        Info<< "    Reading faceIbMask" << endl;
        surfaceScalarField faceIbMask(faceIbMaskHeader, mesh);

        volScalarField contErr = fvc::div(faceIbMask*phi);

        scalar sumLocalContErr = runTime.deltaT().value()*
            mag(contErr)().weightedAverage(mesh.V()).value();

        scalar globalContErr = runTime.deltaT().value()*
            contErr.weightedAverage(mesh.V()).value();

        Info<< "IB time step continuity errors : sum local = "
            << sumLocalContErr
            << ", global = " << globalContErr
            << endl;

        volScalarField magContErr
        (
            "magContErr",
            mag(contErr)
        );
        magContErr.write();

    }
    else
    {
        Info<< "    No phi or faceIbMask" << endl;
    }

}


// ************************************************************************* //

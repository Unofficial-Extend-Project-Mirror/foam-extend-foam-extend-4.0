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
    checkSurfaceCurvature

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Description
    Check surface curvature using Finite Area tools

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

#   include "createFaMesh.H"

    volSurfaceMapping vsm(aMesh);

    // Curvature calculated from point data

    volScalarField curvature
    (
        IOobject
        (
            "curvature",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    );

    const areaScalarField& c = aMesh.faceCurvatures();

    Info<< "Curvature: min = " << Foam::min(c).value()
        << " max = " << Foam::max(c).value() << endl;

    vsm.mapToVolume(c, curvature.boundaryField());

    curvature.write();

    return(0);
}

// ************************************************************************* //

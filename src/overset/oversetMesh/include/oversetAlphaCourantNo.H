/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 Wikki Ltd
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Global
    oversetAlphaCourantNo

Description
    Calculates and outputs the mean and maximum Courant numbers for interface
    cells with overset mesh correction.

\*---------------------------------------------------------------------------*/

scalar maxAlphaCo
(
    readScalar(runTime.controlDict().lookup("maxAlphaCo"))
);

scalar alphaCoNum = 0.0;
scalar meanAlphaCoNum = 0.0;

surfaceScalarField alpha1f =
    fvc::interpolate(min(max(alpha1, scalar(0)), scalar(1)));

const dimensionedScalar alphaOffset("alphaOffset", dimless, dAlpha);

if (mesh.nInternalFaces())
{
    surfaceScalarField magAlphaPhi
    (
        pos(alpha1f - alphaOffset)*
        pos(scalar(1) - alphaOffset - alpha1f)*
        mag(faceOversetMask*phi)
    );

    surfaceScalarField SfUfbyDelta =
        mesh.surfaceInterpolation::deltaCoeffs()*magAlphaPhi;

    const scalar deltaT = runTime.deltaT().value();

    alphaCoNum = max(SfUfbyDelta/mesh.magSf()).value()*deltaT;

    meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf())).value()*deltaT;
}

Info<< "Interface Courant Number mean: " << meanAlphaCoNum
    << " max: " << alphaCoNum << endl;

// ************************************************************************* //

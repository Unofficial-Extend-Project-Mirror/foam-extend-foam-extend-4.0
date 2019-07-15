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

const surfaceScalarField alpha1f =
    fvc::interpolate(min(max(alpha1, scalar(0)), scalar(1)));

const dimensionedScalar alphaOffset
(
    "alphaOffset",
    dimless,
    runTime.controlDict().lookupOrDefault("dAlpha", 0.01)
);

if (mesh.nInternalFaces())
{
    const oversetMesh& om = oversetMesh::New(mesh);

    const surfaceScalarField magAlphaPhi
    (
        pos(alpha1f - alphaOffset)*
        pos(scalar(1) - alphaOffset - alpha1f)*
        mag(om.sGamma()*phi)
    );

    const surfaceScalarField SfUfbyDelta =
        mesh.surfaceInterpolation::deltaCoeffs()*magAlphaPhi;

    const scalar deltaT = runTime.deltaT().value();

    alphaCoNum = max(SfUfbyDelta/mesh.magSf()).value()*deltaT;

    meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf())).value()*deltaT;
}

Info<< "Interface Courant Number mean: " << meanAlphaCoNum
    << " max: " << alphaCoNum << endl;

// ************************************************************************* //

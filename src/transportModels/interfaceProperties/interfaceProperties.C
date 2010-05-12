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
    interfaceProperties

Description
    Properties to aid interFoam :
    1. Correct the gamma boundary condition for dynamic contact angle.
    2. Calculate interface curvature.

\*---------------------------------------------------------------------------*/

#include "interfaceProperties.H"
#include "gammaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const scalar interfaceProperties::convertToRad =
    mathematicalConstant::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void interfaceProperties::correctContactAngle
(
    surfaceVectorField::GeometricBoundaryField& nHatb
) const
{
    const fvMesh& mesh = gamma_.mesh();
    const volScalarField::GeometricBoundaryField& gbf = gamma_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<gammaContactAngleFvPatchScalarField>(gbf[patchi]))
        {
            const gammaContactAngleFvPatchScalarField& gcap = 
                refCast<const gammaContactAngleFvPatchScalarField>
                (gbf[patchi]);

            fvsPatchVectorField& nHatp = nHatb[patchi];
            scalarField theta = 
                convertToRad*gcap.theta(U_.boundaryField()[patchi], nHatp);

            vectorField nf = boundary[patchi].nf();

            // Reset nHatp to correspond to the contact angle

            scalarField a12 = nHatp & nf;

            scalarField b1 = cos(theta);

            scalarField b2(nHatp.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det = 1.0 - a12*a12;

            scalarField a = (b1 - a12*b2)/det;
            scalarField b = (b2 - a12*b1)/det;

            nHatp = a*nf + b*nHatp;

            nHatp /= (mag(nHatp) + deltaN_.value());
        }
    }
}


void interfaceProperties::calculateK()
{
    const fvMesh& mesh = gamma_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of gamma
    volVectorField gradGamma = fvc::grad(gamma_);

    // Interpolated face-gradient of gamma
    surfaceVectorField gradGammaf = fvc::interpolate(gradGamma);
    //gradGammaf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(gamma_) - (mesh.Sf() & gradGammaf)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nHatfv = gradGammaf/(mag(gradGammaf) + deltaN_);
    correctContactAngle(nHatfv.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    /*
    volVectorField nHat = gradGamma/(mag(gradGamma) + deltaN_);
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryField()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interfaceProperties::interfaceProperties
(
    const volScalarField& gamma,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cGamma_
    (
        readScalar
        (
            gamma.mesh().solutionDict().subDict("PISO").lookup("cGamma")
        )
    ),
    sigma_(dict.lookup("sigma")),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(gamma.mesh().V()), 1.0/3.0)
    ),

    gamma_(gamma),
    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            gamma_.time().timeName(),
            gamma_.mesh()
        ),
        gamma_.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),

    K_
    (
        IOobject
        (
            "K",
            gamma_.time().timeName(),
            gamma_.mesh()
        ),
        gamma_.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    )
{
    calculateK();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

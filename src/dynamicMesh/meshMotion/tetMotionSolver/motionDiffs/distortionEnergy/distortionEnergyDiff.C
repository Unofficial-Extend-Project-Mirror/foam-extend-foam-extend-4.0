/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Description
    Non-orthogonality-based motion diffusivity.

\*---------------------------------------------------------------------------*/

#include "distortionEnergyDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "elementFields.H"
#include "cellQuality.H"
#include "tetFec.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distortionEnergyDiff, 0);
    addToRunTimeSelectionTable(motionDiff, distortionEnergyDiff, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::distortionEnergyDiff::distortionEnergyDiff
(
    const tetMotionSolver& mSolver
)
:
    motionDiff(mSolver),
    motionGamma_
    (
        IOobject
        (
            "motionGamma",
            tetMesh().time().timeName(),
            tetMesh()(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tetMesh(),
        dimensionedScalar("1.0", dimless, 1.0)
    )
{
    exponent_ = readInt(mSolver.lookup("diffusivityExponent"));

    Info << "Value of exponent for distortion energy based motion diffusivity: "
        << exponent_ << endl;


    mSolver.storeTotDisplacement();

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distortionEnergyDiff::~distortionEnergyDiff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distortionEnergyDiff::correct()
{
    motionGamma_.internalField() =
        pow(mSolver().totDistortionEnergy()().internalField(), exponent_)
      + SMALL;

    motionGamma_.internalField() /= max(motionGamma_.internalField());
}

// ************************************************************************* //

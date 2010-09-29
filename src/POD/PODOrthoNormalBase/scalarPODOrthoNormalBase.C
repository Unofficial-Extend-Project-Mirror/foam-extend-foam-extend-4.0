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

Class
    PODOrthoNormalBase

Description
    Establish POD ortho-normal base and interpolation coefficients give a list
    of fields. Size of ortho-normal base is calculated from the desired
    accuracy, e.g. 0.99-0.99999 (in energy terms)

\*---------------------------------------------------------------------------*/

#include "scalarPODOrthoNormalBase.H"
#include "PODEigenBase.H"
#include "IOmanip.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<>
void Foam::PODOrthoNormalBase<Foam::scalar>::calcOrthoBase
(
    const PtrList<volScalarField>& snapshots,
    const scalar accuracy
)
{
    PODEigenBase eigenBase(snapshots);

    label baseSize = 0;

    const scalarField& cumEigenValues = eigenBase.cumulativeEigenValues();

    forAll (cumEigenValues, i)
    {
        baseSize++;

        if (cumEigenValues[i] > accuracy)
        {
            break;
        }
    }

    Info<< "Cumulative eigen-values: "
        << setprecision(14) << cumEigenValues << nl
        << "Base size: " << baseSize << endl;

    // Establish orthonormal base
    orthoFields_.setSize(baseSize);

    for (label baseI = 0; baseI < baseSize; baseI++)
    {
        const scalarField& eigenVector = eigenBase.eigenVectors()[baseI];

        volScalarField* onBasePtr
        (
            new volScalarField
            (
                IOobject
                (
                    snapshots[0].name() + "POD" + name(baseI),
                    snapshots[0].time().timeName(),
                    snapshots[0].mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                snapshots[0].mesh(),
                dimensionedScalar("zero", snapshots[0].dimensions(), 0)
            )
        );
        volScalarField& onBase = *onBasePtr;

        forAll (eigenVector, eigenI)
        {
            onBase += eigenVector[eigenI]*snapshots[eigenI];
        }

        // Re-normalise ortho-normal vector
        scalar magSumSquare = Foam::sqrt(sumSqr(onBase));
        if (magSumSquare > SMALL)
        {
            onBase /= magSumSquare;
            onBase.correctBoundaryConditions();
        }

        orthoFields_.set(baseI, onBasePtr);
    }

    // Calculate interpolation coefficients
    interpolationCoeffsPtr_ =
        new RectangularMatrix<scalar>(snapshots.size(), orthoFields_.size());
    RectangularMatrix<scalar>& coeffs = *interpolationCoeffsPtr_;

    forAll (snapshots, snapshotI)
    {
        forAll (orthoFields_, baseI)
        {
            coeffs[snapshotI][baseI] =
                POD::projection
                (
                    snapshots[snapshotI],
                    orthoFields_[baseI]
                );
        }
    }
}


// ************************************************************************* //

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

#include "PODOrthoNormalBase.H"
#include "POD.H"
#include "PODEigenBase.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::PODOrthoNormalBase<Type>::calcOrthoBase
(
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& snapshots,
    const scalar accuracy
)
{
    // Calculate ortho-normal base for each component
    PtrList<PODEigenBase> eigenBase(pTraits<Type>::nComponents);

    const label nSnapshots = snapshots.size();

    typename
    powProduct<Vector<label>, pTraits<Type>::rank>::type validComponents
    (
        pow
        (
            snapshots[0].mesh().solutionD(),
            pTraits
            <
                typename powProduct<Vector<label>,
                pTraits<Type>::rank
            >::type>::zero
        )
    );

    label nValidCmpts = 0;

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        // Component not valid, skipping
        if (validComponents[cmpt] == -1) continue;

        // Create a list of snapshots
        PtrList<volScalarField> sf(nSnapshots);

        forAll (snapshots, i)
        {
            sf.set(i, new volScalarField(snapshots[i].component(cmpt)));
        }

        // Create eigen base
        eigenBase.set(cmpt, new PODEigenBase(sf));

        Info<< "Cumulative eigen-values for component " << cmpt
            << ": " << setprecision(14)
            << eigenBase[nValidCmpts].cumulativeEigenValues() << endl;

        nValidCmpts++;
    }

    eigenBase.setSize(nValidCmpts);

    Info << "Number of valid eigen components: " << nValidCmpts << endl;

    label baseSize = 0;
    for (label snapI = 0; snapI < nSnapshots; snapI++)
    {
        baseSize++;

        // Get minimum cumulative eigen value for all valid components
        scalar minCumEigen = 1.0;

        nValidCmpts = 0;

        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
        {
            // Skip invalid components
            if (validComponents[cmpt] != -1)
            {
                minCumEigen =
                    Foam::min
                    (
                        minCumEigen,
                        eigenBase[nValidCmpts].cumulativeEigenValues()[snapI]
                    );

                nValidCmpts++;
            }
        }

        if (minCumEigen > accuracy)
        {
            break;
        }
    }

    Info << "Base size: " << baseSize << endl;

    // Establish orthonormal base
    orthoFields_.setSize(baseSize);

    for (label baseI = 0; baseI < baseSize; baseI++)
    {
        GeometricField<Type, fvPatchField, volMesh>* onBasePtr
        (
            new GeometricField<Type, fvPatchField, volMesh>
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
                dimensioned<Type>
                (
                    "zero",
                    snapshots[0].dimensions(),
                    pTraits<Type>::zero
                )
            )
        );
        GeometricField<Type, fvPatchField, volMesh>& onBase = *onBasePtr;

        nValidCmpts = 0;

        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
        {
            if (validComponents[cmpt] != -1)
            {
                // Valid component, grab eigen-factors

                const scalarField& eigenVector =
                    eigenBase[nValidCmpts].eigenVectors()[baseI];
                nValidCmpts++;

                volScalarField onBaseCmpt = onBase.component(cmpt);

                forAll (eigenVector, eigenI)
                {
                    onBaseCmpt +=
                        eigenVector[eigenI]*snapshots[eigenI].component(cmpt);
                }

                // Re-normalise ortho-normal vector
                onBaseCmpt /= Foam::sqrt(sumSqr(onBaseCmpt));

                onBase.replace(cmpt, onBaseCmpt);
            }
            else
            {
                // Component invalid.  Grab first snapshot.
                onBase.replace
                (
                    cmpt,
                    snapshots[0].component(cmpt)
                );
            }
        }

        orthoFields_.set(baseI, onBasePtr);
    }

    // Calculate interpolation coefficients
    interpolationCoeffsPtr_ =
        new RectangularMatrix<Type>(snapshots.size(), orthoFields_.size());
    RectangularMatrix<Type>& coeffs = *interpolationCoeffsPtr_;

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// given list of snapshots and accuracy
template<class Type>
Foam::PODOrthoNormalBase<Type>::PODOrthoNormalBase
(
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& snapshots,
    const scalar accuracy
)
:
    orthoFields_(),
    interpolationCoeffsPtr_(NULL)
{
    calcOrthoBase(snapshots, accuracy);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::PODOrthoNormalBase<Type>::~PODOrthoNormalBase()
{
    deleteDemandDrivenData(interpolationCoeffsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //

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
    PODEigenBase

Description

\*---------------------------------------------------------------------------*/

#include "PODEigenBase.H"
#include "volFields.H"
#include "POD.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PODEigenBase::calcEigenBase(const scalarSquareMatrix& orthMatrix)
{
    // Calculate eigen-values

    EigenSolver<scalar> eigenSolver(orthMatrix);

    // Sort and assemble

    SortableList<scalar> sortedList(orthMatrix.n());

    forAll (sortedList, i)
    {
        sortedList[i] = eigenSolver.eigenValue(i);
    }

    // Sort will sort the values in descending order and insert values
    sortedList.sort();

    label n = 0;
    forAllReverse(sortedList, i)
    {
        eigenValues_[n] = sortedList[i];
        eigenVectors_.set
        (
            n,
            new scalarField(eigenSolver.eigenVector(sortedList.indices()[i]))
        );

        n++;
    }

    // Assemble cumulative relative eigen-values
    cumEigenValues_[0] = eigenValues_[0];

    // Assemble accumulated normalised eigenvalues
    for (label i = 1; i < cumEigenValues_.size(); i++)
    {
        cumEigenValues_[i] = cumEigenValues_[i - 1] + eigenValues_[i];
    }

    // Renormalise
    cumEigenValues_ /= sum(eigenValues_);

//     // Check products
//     for (label i = 0; i < orthMatrix.m(); i++)
//     {
//         const scalarField& eVector = eigenVectors_[i];

//         Info<< "Scalar product: "
//             << eigenValues_[i]*eVector
//             << endl;

//         scalarField vp(orthMatrix.m(), 0);

//         forAll (vp, i)
//         {
//             forAll (vp, j)
//             {
//                 vp[i] += orthMatrix[i][j]*eVector[j];
//             }
//         }

//         Info << "Vector product: " << vp << nl << endl;
//     }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given a list of fields
Foam::PODEigenBase::PODEigenBase(const PtrList<volScalarField>& snapshots)
:
    eigenValues_(snapshots.size()),
    cumEigenValues_(snapshots.size()),
    eigenVectors_(snapshots.size())
{
    // Calculate the snapshot of the field with all available fields
    label nSnapshots = snapshots.size();

    scalarSquareMatrix orthMatrix(nSnapshots);

    for (label snapI = 0; snapI < nSnapshots; snapI++)
    {
        for (label snapJ = 0; snapJ <= snapI; snapJ++)
        {
            // Calculate the inner product and insert it into the matrix
            orthMatrix[snapI][snapJ] =
                POD::projection
                (
                    snapshots[snapI],
                    snapshots[snapJ]
                );

            if (snapI != snapJ)
            {
                orthMatrix[snapJ][snapI] = orthMatrix[snapI][snapJ];

//                 Info << "Product: " << orthMatrix[snapI][snapJ]
//                     << " relative: "
//                     <<
//                     orthMatrix[snapI][snapJ]/
//                     (
//                         Foam::sqrt(sumSqr(snapshots[snapI]))*
//                         Foam::sqrt(sumSqr(snapshots[snapJ]))
//                       + SMALL
//                     ) << endl;
            }
        }
    }

    calcEigenBase(orthMatrix);
}


// ************************************************************************* //

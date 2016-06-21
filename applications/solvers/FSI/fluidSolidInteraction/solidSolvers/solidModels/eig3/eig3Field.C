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
    Function to calculate eigen values and eigen vectors of volSymmTensorField
    Using the main procedure/code from here:
    http://barnesc.blogspot.ie/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
    and code here:
    http://www.connellybarnes.com/code/c/eig3-1.0.0.zip
    Note: built-in OpenFOAM functions mess-up on a number of different tensors.

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "eig3Field.H"
#include "eig3.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void eig3Field
(
    const volTensorField& A, volTensorField& V, volVectorField& d
)
{
    const tensorField& AI = A.internalField();
    tensorField& VI = V.internalField();
    vectorField& dI = d.internalField();

    forAll(AI, cellI)
    {
        eigen_decomposition(AI[cellI], VI[cellI], dI[cellI]);
    }

    forAll(A.boundaryField(), patchI)
    {
        if
        (
           !A.boundaryField()[patchI].coupled()
         && A.boundaryField()[patchI].type() != "empty"
        )
        {
            const tensorField& AB = A.boundaryField()[patchI];
            tensorField& VB = V.boundaryField()[patchI];
            vectorField& dB = d.boundaryField()[patchI];

            forAll(AB, faceI)
            {
                eigen_decomposition(AB[faceI], VB[faceI], dB[faceI]);
            }
        }
    }

//     V.correctBoundaryConditions();
//     d.correctBoundaryConditions();
}


void eig3Field
(
    const volSymmTensorField& A, volTensorField& V, volVectorField& d
)
{
    const symmTensorField& AI = A.internalField();
    tensorField& VI = V.internalField();
    vectorField& dI = d.internalField();

    forAll(AI, cellI)
    {
        eigen_decomposition(AI[cellI], VI[cellI], dI[cellI]);
    }

    forAll(A.boundaryField(), patchI)
    {
        if
        (
            !A.boundaryField()[patchI].coupled()
            && A.boundaryField()[patchI].type() != "empty"
        )
        {
            const symmTensorField& AB = A.boundaryField()[patchI];
            tensorField& VB = V.boundaryField()[patchI];
            vectorField& dB = d.boundaryField()[patchI];

            forAll(AB, faceI)
            {
                eigen_decomposition(AB[faceI], VB[faceI], dB[faceI]);
            }
        }
    }

//     V.correctBoundaryConditions();
//     d.correctBoundaryConditions();
}


void eig3Field
(
    const volSymmTensorField& A, volTensorField& V, volDiagTensorField& d
)
{
    const symmTensorField& AI = A.internalField();
    tensorField& VI = V.internalField();
    diagTensorField& dI = d.internalField();

    forAll(AI, cellI)
    {
        eigen_decomposition(AI[cellI], VI[cellI], dI[cellI]);
    }

    forAll(A.boundaryField(), patchI)
    {
        if
        (
            !A.boundaryField()[patchI].coupled()
            && A.boundaryField()[patchI].type() != "empty"
        )
        {
            const symmTensorField& AB = A.boundaryField()[patchI];
            tensorField& VB = V.boundaryField()[patchI];
            diagTensorField& dB = d.boundaryField()[patchI];

            forAll(AB, faceI)
            {
                eigen_decomposition(AB[faceI], VB[faceI], dB[faceI]);
            }
        }
    }

//     V.correctBoundaryConditions();
//     d.correctBoundaryConditions();
}


void eig3Field
(
    const pointSymmTensorField& A,
    pointTensorField& V,
    pointVectorField& d
)
{
    const symmTensorField& AI = A.internalField();
    tensorField& VI = V.internalField();
    vectorField& dI = d.internalField();

    forAll(AI, pointI)
    {
        eigen_decomposition(AI[pointI], VI[pointI], dI[pointI]);
    }

//     V.correctBoundaryConditions();
//     d.correctBoundaryConditions();
}

tmp<volSymmTensorField> log(const volSymmTensorField& vf)
{
    tmp<volSymmTensorField> tresult
    (
        new volSymmTensorField
        (
            IOobject
            (
                "log("+vf.name()+")",
                vf.time().timeName(),
                vf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vf
        )
    );

    volSymmTensorField& result = tresult();

    // Calculate eigen values and eigen vectors
    // The OpenFOAM eigenValues/eigenVectors sometimes give wrong results when
    // eigenValues are repeated or zero, so I will use my own implementation.
    // The efficiency of the implementation may need to be revisited, however,
    // it is fine for creation of post processing fields e.g calculate true
    // strain

    volDiagTensorField eigenVal
    (
        IOobject
        (
            "eigenVal("+vf.name()+")",
            vf.time().timeName(),
            vf.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vf.mesh(),
        dimensionedDiagTensor("zero", vf.dimensions(), diagTensor::zero)
    );

    volTensorField eigenVec
    (
        IOobject
        (
            "eigenVec("+vf.name()+")",
            vf.time().timeName(),
            vf.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vf.mesh(),
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    // Calculate eigen values and eigen vectors of vf
    eig3Field(vf, eigenVec, eigenVal);

    // Now we will calculate the log og the eigenValues and then rotate the
    // tensor back to the physcial configuration

    const diagTensorField& eigenValI = eigenVal.internalField();
    const tensorField& eigenVecI = eigenVec.internalField();
    symmTensor logEigenVal = symmTensor::zero;
    symmTensorField& resultI = result.internalField();

    forAll(eigenValI, cellI)
    {
        // We cannot have log of zeros
        if
        (
            eigenValI[cellI][diagTensor::XX] < SMALL
         || eigenValI[cellI][diagTensor::YY] < SMALL
         || eigenValI[cellI][diagTensor::ZZ] < SMALL
        )
        {
            FatalError
                << "log of zero is not allowed"
                << abort(FatalError);
        }

        // Calculate log
        logEigenVal[symmTensor::XX] =
            Foam::log(eigenValI[cellI][diagTensor::XX]);
        logEigenVal[symmTensor::YY] =
            Foam::log(eigenValI[cellI][diagTensor::YY]);
        logEigenVal[symmTensor::ZZ] =
            Foam::log(eigenValI[cellI][diagTensor::ZZ]);

        // Rotate back
        resultI[cellI] =
            symm(eigenVecI[cellI].T() & logEigenVal & eigenVecI[cellI]);
//         resultI[cellI] = transform(eigenVecI[cellI].T(), logEigenVal);
    }

    forAll(eigenVal.boundaryField(), patchI)
    {
        if
        (
           !vf.boundaryField()[patchI].coupled()
         &&
            vf.boundaryField()[patchI].type()
         != emptyFvPatchField<symmTensor>::typeName
        )
        {
            const diagTensorField& eigenValB =
                eigenVal.boundaryField()[patchI];
            const tensorField& eigenVecB = eigenVec.boundaryField()[patchI];
            symmTensorField& resultB = result.boundaryField()[patchI];

            forAll(eigenValB, faceI)
            {
                // We cannot have log of zeros
                if
                (
                    eigenValB[faceI][diagTensor::XX] < SMALL
                 || eigenValB[faceI][diagTensor::YY] < SMALL
                 || eigenValB[faceI][diagTensor::ZZ] < SMALL
                )
                {
                    FatalError
                        << "log of zero is not allowed"
                        << abort(FatalError);
                }

                // Calculate log
                logEigenVal[symmTensor::XX] =
                    Foam::log(eigenValB[faceI][diagTensor::XX]);
                logEigenVal[symmTensor::YY] =
                    Foam::log(eigenValB[faceI][diagTensor::YY]);
                logEigenVal[symmTensor::ZZ] =
                    Foam::log(eigenValB[faceI][diagTensor::ZZ]);

                // Rotate back
                resultB[faceI] =
                    symm(eigenVecB[faceI].T() & logEigenVal & eigenVecB[faceI]);
//                 resultB[faceI] = transform(eigenVecB[faceI].T(), logEigenVal);
            }
        }
    }

    result.correctBoundaryConditions();

    return tresult;
}

symmTensor log(const symmTensor& t)
{
    symmTensor result = symmTensor::zero;

    diagTensor eigenVal = diagTensor::zero;
    tensor eigenVec = tensor::zero;

    // Calculate eigen values and eigen vectors of t
    eigen_decomposition(t, eigenVec, eigenVal);

    symmTensor logEigenVal = symmTensor::zero;

    if
    (
        eigenVal[diagTensor::XX] < SMALL
     || eigenVal[diagTensor::YY] < SMALL
     || eigenVal[diagTensor::ZZ] < SMALL
    )
    {
        FatalError
            << "log of zero is not allowed"
                << abort(FatalError);
    }

    // Calculate log
    logEigenVal[symmTensor::XX] =
        Foam::log(eigenVal[diagTensor::XX]);
    logEigenVal[symmTensor::YY] =
        Foam::log(eigenVal[diagTensor::YY]);
    logEigenVal[symmTensor::ZZ] =
        Foam::log(eigenVal[diagTensor::ZZ]);

    // Rotate back
    result = symm(eigenVec.T() & logEigenVal & eigenVec);

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

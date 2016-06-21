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

\*---------------------------------------------------------------------------*/

#include "extendedLduMatrix.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedLduMatrix, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::extendedLduMatrix::clearOut()
{
    deleteDemandDrivenData(extendedLowerPtr_);
    deleteDemandDrivenData(extendedUpperPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedLduMatrix::extendedLduMatrix
(
    const lduMatrix& ldum,
    const extendedLduAddressing& extLduAddr
)
:
    basicLduMatrix_(ldum),
    extLduAddr_(extLduAddr),
    extendedLowerPtr_(NULL),
    extendedUpperPtr_(NULL)
{
    if (debug)
    {
        InfoIn("extendedLduMatrix(lduMatrix&, label, polyMesh&)")
            << "Constructing extendedLduMatrix."
            << endl;
    }

    if (ldum.diagonal())
    {
        WarningIn("extendedLduMatrix(lduMatrix&, label, polyMesh&)")
            << "Attempted to create extended lower/upper coeffs for matrix "
            << "that is diagonal."
            << nl << endl;
    }
    else if (ldum.symmetric())
    {
        // Get reference to faceMap in extended addressing
        const unallocLabelList& faceMap = extLduAddr_.faceMap();

        // Matrix is considered symmetric if the upper is allocated and lower
        // is not allocated. Allocating extended upper only.
        extendedUpperPtr_ = new scalarField
        (
            extLduAddr_.extendedUpperAddr().size(),
            0.0
        );
        scalarField& extUpper = *extendedUpperPtr_;

        // Get upper coeffs from underlying lduMatrix
        const scalarField& upper = ldum.upper();

        // Copy non-zero coeffs from basic lduMatrix into corresponding
        // positions
        forAll (upper, faceI)
        {
            extUpper[faceMap[faceI]] = upper[faceI];
        }
    }
    else
    {
        // Get reference to faceMap in extended addressing
        const unallocLabelList& faceMap = extLduAddr_.faceMap();

        // Get number of extended faces
        const label nExtFaces = extLduAddr_.extendedUpperAddr().size();

        // Allocate extended upper and lower
        extendedUpperPtr_ = new scalarField(nExtFaces, scalar(0));
        scalarField& extUpper = *extendedUpperPtr_;

        extendedLowerPtr_ = new scalarField(nExtFaces, scalar(0));
        scalarField& extLower = *extendedLowerPtr_;

        // Get upper and lower coeffs from underlying lduMatrix
        const scalarField& upper = ldum.upper();
        const scalarField& lower = ldum.lower();

        // Copy non-zero coeffs from basic lduMatrix into corresponding
        // positions
        forAll (upper, faceI)
        {
            extUpper[faceMap[faceI]] = upper[faceI];
            extLower[faceMap[faceI]] = lower[faceI];
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::extendedLduMatrix::~extendedLduMatrix()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::extendedLduMatrix::extendedLower()
{
    if (!extendedLowerPtr_)
    {
        if (extendedUpperPtr_)
        {
            extendedLowerPtr_ = new scalarField(*extendedUpperPtr_);
        }
        else
        {
            extendedLowerPtr_ = new scalarField
            (
                extLduAddr_.extendedLowerAddr().size(),
                scalar(0)
            );
        }
    }

    return *extendedLowerPtr_;
}


Foam::scalarField& Foam::extendedLduMatrix::extendedUpper()
{
    if (!extendedUpperPtr_)
    {
        if (extendedLowerPtr_)
        {
            extendedUpperPtr_ = new scalarField(*extendedLowerPtr_);
        }
        else
        {
            extendedUpperPtr_ = new scalarField
            (
                extLduAddr_.extendedUpperAddr().size(),
                scalar(0)
            );
        }
    }

    return *extendedUpperPtr_;
}


const Foam::scalarField& Foam::extendedLduMatrix::extendedLower() const
{
    if (!extendedLowerPtr_ && !extendedUpperPtr_)
    {
        FatalErrorIn
        (
            "const scalarfield& extendedLduMatrix::extendedLower() const"
        )   << "extendedLowerPtr_ or extendedUpperPtr_ unallocated"
            << abort(FatalError);
    }

    if (extendedLowerPtr_)
    {
        return *extendedLowerPtr_;
    }
    else
    {
        return *extendedUpperPtr_;
    }
}


const Foam::scalarField& Foam::extendedLduMatrix::extendedUpper() const
{
    if (!extendedLowerPtr_ && !extendedUpperPtr_)
    {
        FatalErrorIn
        (
            "const scalarfield& extendedLduMatrix::extendedUpper() const"
        )   << "extendedLowerPtr_ or extendedUpperPtr_ unallocated"
            << abort(FatalError);
    }

    if (extendedUpperPtr_)
    {
        return *extendedUpperPtr_;
    }
    else
    {
        return *extendedLowerPtr_;
    }
}


// ************************************************************************* //

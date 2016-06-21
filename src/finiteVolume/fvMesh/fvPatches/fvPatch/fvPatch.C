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

#include "fvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvBoundaryMesh.H"
#include "fvMesh.H"
#include "primitiveMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fvPatch, 0);
defineRunTimeSelectionTable(fvPatch, polyPatch);
addToRunTimeSelectionTable(fvPatch, fvPatch, polyPatch);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fvPatch::fvPatch(const polyPatch& p, const fvBoundaryMesh& bm)
:
    polyPatch_(p),
    boundaryMesh_(bm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fvPatch::~fvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool fvPatch::constraintType(const word& pt)
{
    return fvPatchField<scalar>::patchConstructorTablePtr_->found(pt);
}


wordList fvPatch::constraintTypes()
{
    wordList cTypes(polyPatchConstructorTablePtr_->size());

    label i = 0;

    for
    (
        polyPatchConstructorTable::iterator cstrIter =
            polyPatchConstructorTablePtr_->begin();
        cstrIter != polyPatchConstructorTablePtr_->end();
        ++cstrIter
    )
    {
        if (constraintType(cstrIter.key()))
        {
            cTypes[i++] = cstrIter.key();
        }
    }

    cTypes.setSize(i);

    return cTypes;
}


const unallocLabelList& fvPatch::faceCells() const
{
    return polyPatch_.faceCells();
}


const vectorField& fvPatch::Cf() const
{
    return boundaryMesh().mesh().Cf().boundaryField()[index()];
}


tmp<vectorField> fvPatch::Cn() const
{
    // Since empty patch has a different size that a polyPatch, cutting
    // needs to be performed again.  HJ, 28/Dec/2006

    tmp<vectorField> tcc(new vectorField(size()));
    vectorField& cc = tcc();

    const unallocLabelList& faceCells = this->faceCells();

    // Get reference to global cell centres
    const vectorField& gcc = boundaryMesh().mesh().cellCentres();

    forAll (faceCells, faceI)
    {
        cc[faceI] = gcc[faceCells[faceI]];
    }

    return tcc;
}


tmp<vectorField> fvPatch::nf() const
{
    return Sf()/magSf();
}


const vectorField& fvPatch::Sf() const
{
    return boundaryMesh().mesh().Sf().boundaryField()[index()];
}


const scalarField& fvPatch::magSf() const
{
    return boundaryMesh().mesh().magSf().boundaryField()[index()];
}


tmp<vectorField> fvPatch::delta() const
{
    return Cf() - Cn();
}


void fvPatch::makeWeights(scalarField& w) const
{
    w = 1.0;
}


void fvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    dc = 1.0/(nf() & delta());
}


void fvPatch::makeCorrVecs(vectorField& cv) const
{
    cv = vector::zero;
}


void fvPatch::initMovePoints()
{}


void fvPatch::movePoints()
{}


const scalarField& fvPatch::deltaCoeffs() const
{
    return boundaryMesh().mesh().deltaCoeffs().boundaryField()[index()];
}


const scalarField& fvPatch::weights() const
{
    return boundaryMesh().mesh().weights().boundaryField()[index()];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

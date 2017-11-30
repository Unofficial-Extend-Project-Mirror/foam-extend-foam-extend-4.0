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

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "slicedVolFields.H"
#include "slicedSurfaceFields.H"
#include "immersedBoundaryFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryFvPatch, 0);

    addToRunTimeSelectionTable(fvPatch, immersedBoundaryFvPatch, polyPatch);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::immersedBoundaryFvPatch::makeCf(slicedSurfaceVectorField& Cf) const
{
    // Correct face centres by re-cutting and inserting the immersed patch
    Info<< "Re-slicing Cf" << endl;
    Cf.reset(ibPolyPatch_.correctedFaceCentres());

    // Insert the patch data for the immersed boundary
    // Note: use the face centres from the stand-alone patch within the IB
    // HJ, 30/Nov/2017
    Cf.boundaryField()[index()].UList::operator=
    (
        vectorField::subField(ibPolyPatch_.ibPatch().faceCentres(), size())
    );
}


void Foam::immersedBoundaryFvPatch::makeSf(slicedSurfaceVectorField& Sf) const
{
    // Correct face centres by re-cutting and inserting the immersed patch
    Info<< "Re-slicing Sf" << endl;
    Sf.reset(ibPolyPatch_.correctedFaceAreas());

    // Insert the patch data for the immersed boundary
    // Note: use the face centres from the stand-alone patch within the IB
    // HJ, 30/Nov/2017
    Sf.boundaryField()[index()].UList::operator=
    (
        vectorField::subField(ibPolyPatch_.ibPatch().areas(), size())
    );
}


void Foam::immersedBoundaryFvPatch::makeC(slicedVolVectorField& C) const
{
    // Correct face centres by re-cutting and inserting the immersed patch
    Info<< "Re-slicing C" << endl;
    C.reset
    (
        ibPolyPatch_.correctedCellCentres(),
        ibPolyPatch_.correctedFaceCentres()
    );

    // Insert the patch data for the immersed boundary
    // Note: use the face centres from the stand-alone patch within the IB
    // HJ, 30/Nov/2017
    C.boundaryField()[index()].UList::operator=
    (
        vectorField::subField(ibPolyPatch_.ibPatch().faceCentres(), size())
    );
}


void Foam::immersedBoundaryFvPatch::makeV(scalarField& V) const
{
    // Correct face centres by re-cutting and inserting the immersed patch
    Info<< "Re-slicing V" << endl;
    V = ibPolyPatch_.correctedCellVolumes();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryFvPatch::immersedBoundaryFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    fvPatch(patch, bm),
    ibPolyPatch_(refCast<const immersedBoundaryPolyPatch>(patch)),
    mesh_(bm.mesh())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::immersedBoundaryFvPatch::size() const
{
    // Immersed boundary patch size equals to the number of intersected cells
    // HJ, 28/Nov/2017
    return ibPolyPatch_.ibCells().size();
}


const Foam::unallocLabelList&
Foam::immersedBoundaryFvPatch::faceCells() const
{
    return ibPolyPatch_.ibCells();
}


// ************************************************************************* //

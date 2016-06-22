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

#include "mgMeshLevel.H"
#include "coarseMgMeshLevel.H"

extern "C"
{
#   include "mgridgen.h"
#   ifdef darwin
#       undef FALSE
#       undef TRUE
#   endif

#undef sign
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::mgMeshLevel, 0);


const Foam::debug::optimisationSwitch Foam::mgMeshLevel::mgMinClusterSize_
(
    "mgMinClusterSize",
    2
);

const Foam::debug::optimisationSwitch Foam::mgMeshLevel::mgMaxClusterSize_
(
    "mgMaxClusterSize",
    8
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mgMeshLevel::makeCompactAddressingAndWeights
(
    labelList& cellCells,
    labelList& cellCellOffsets,
    scalarField& faceWeights,
    scalarField& boundaryAreas
) const
{
    const unallocLabelList& upperAddr = owner();
    const unallocLabelList& lowerAddr = neighbour();

    const label nFineCells = nCells();
    const label nFineFaces = nInternalFaces();

    // Number of neighbours for each cell
    labelList nNbrs(nFineCells, 0);

    forAll (upperAddr, facei)
    {
        nNbrs[upperAddr[facei]]++;
    }

    forAll (lowerAddr, facei)
    {
        nNbrs[lowerAddr[facei]]++;
    }

    // Set the sizes of the addressing and faceWeights arrays
    cellCellOffsets.setSize(nFineCells + 1);
    cellCells.setSize(2*nFineFaces);
    faceWeights.setSize(2*nFineFaces);


    // Assemble cell offsets
    cellCellOffsets[0] = 0;
    forAll (nNbrs, celli)
    {
        cellCellOffsets[celli+1] = cellCellOffsets[celli] + nNbrs[celli];
    }

    // Reset the whole list to use as counter
    nNbrs = 0;

    const vectorField& Si = faceAreas();

    forAll (upperAddr, facei)
    {
        label own = upperAddr[facei];
        label nei = lowerAddr[facei];

        label l1 = cellCellOffsets[own] + nNbrs[own]++;
        label l2 = cellCellOffsets[nei] + nNbrs[nei]++;

        cellCells[l1] = nei;
        cellCells[l2] = own;

        faceWeights[l1] = mag(Si[facei]);
        faceWeights[l2] = mag(Si[facei]);
    }

    // Sum up boundary areas
    boundaryAreas.setSize(nFineCells, 0);

    for (label patchI = 0; patchI < nPatches(); patchI++)
    {
        const labelList& fc = faceCells(patchI);
        const scalarField& patchAreas = magPatchFaceAreas(patchI);

        forAll (fc, faceI)
        {
            boundaryAreas[fc[faceI]] += patchAreas[faceI];
        }
    }
}


void Foam::mgMeshLevel::makeChild() const
{
    // Compact addressing for cellCells
    labelList cellCells;
    labelList cellCellOffsets;

    // Face weights = face areas of the internal faces (double size)
    scalarField faceWeights;

    // Area sum of boundary faces for each cell
    scalarField boundaryAreas;

    makeCompactAddressingAndWeights
    (
        cellCells,
        cellCellOffsets,
        faceWeights,
        boundaryAreas
    );

    // MGridGen agglomeration options.
    List<int> options(4, 0);
    options[0] = 4;                   // globular agglom
    options[1] = 6;                   // objective F3 and F2
    options[2] = 128;                 // debugging output level
    options[3] = nGeometricD();       // Dimensionality of the grid

    // Output: cell to coarse clusted addressing
    int nCoarseCells = 0;
    child_.setSize(nCells());
    int nMoves = -1;

#   ifdef WM_DP

    MGridGen
    (
        nCells(),
        cellCellOffsets.begin(),
        const_cast<scalar*>(cellVolumes().begin()),
        const_cast<scalar*>(boundaryAreas.begin()),
        cellCells.begin(),
        faceWeights.begin(),
        mgMinClusterSize_(),
        mgMaxClusterSize_(),
        options.begin(),
        &nMoves,
        &nCoarseCells,
        child_.begin()
    );

#   else

    // Conversion of type for MGridGen interface
    const scalarField& vols = cellVolumes();

    Field<double> dblVols(vols.size());
    forAll (dblVols, i)
    {
        dblVols[i] = vols[i];
    }

    Field<double> dblAreas(boundaryAreas.size());
    forAll (dblAreas, i)
    {
        dblAreas[i] = boundaryAreas[i];
    }

    Field<double> dblFaceWeights(faceWeights.size());
    forAll (dblFaceWeights, i)
    {
        dblFaceWeights[i] = faceWeights[i];
    }

    MGridGen
    (
        nCells(),
        cellCellOffsets.begin(),
        dblVols.begin(),
        dblAreas.begin(),
        cellCells.begin(),
        dblFaceWeights.begin(),
        mgMinClusterSize_(),
        mgMaxClusterSize_(),
        options.begin(),
        &nMoves,
        &nCoarseCells_,
        child_.begin()
    );

#   endif

    Info<< "Number of coarse cells = " << nCoarseCells_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::mgMeshLevel&  Foam::mgMeshLevel::finestLevel() const
{
    if (finest())
    {
        return *this;
    }
    else
    {
        return fineLevel();
    }
}

Foam::label Foam::mgMeshLevel::nCoarseCells() const
{
    if (nCoarseCells_ < 0)
    {
        makeChild();
    }

    return nCoarseCells_;
}


const Foam::labelList& Foam::mgMeshLevel::child() const
{
    if (child_.empty())
    {
        makeChild();
    }

    return child_;
}


Foam::autoPtr<Foam::mgMeshLevel> Foam::mgMeshLevel::makeNextLevel() const
{
    return autoPtr<mgMeshLevel>(new coarseMgMeshLevel(*this));
}


// ************************************************************************* //

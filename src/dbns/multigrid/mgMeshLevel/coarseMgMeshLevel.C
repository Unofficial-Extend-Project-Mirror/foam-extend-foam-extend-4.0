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

#include "coarseMgMeshLevel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::coarseMgMeshLevel, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coarseMgMeshLevel::calcLevel()
{
    // Check agglomeration
    const labelList& child = fineLevel_.child();

    if (max(child) + 1 != fineLevel_.nCoarseCells())
    {
        FatalErrorIn("void coarseMgMeshLevel::calcLevel()")
            << "Error in child array assembly: "
            << max(child) + 1
            << " and " << fineLevel_.nCoarseCells()
            << abort(FatalError);
    }

    // Record number of coarse cells
    nCells_ = fineLevel_.nCoarseCells();

    // Construct the coarse mesh and ldu addressing for the next level
    // Algorithm:
    // 1) Loop through all fine faces. If the child labels on two sides are
    //    different, this creates a coarse face. Define owner and neighbour
    //    for this face based on cluster IDs.
    // 2) Check if the face has been seen before. If yes, add the fine face
    //    to the appropriate face cluster.  If no, create
    //    a new face with neighbour ID and add the fine face to it
    // 3) Once all the faces have been created, loop through all clusters and
    //    insert the faces in the upper order. At the same time, collect the
    //    owner and neighbour addressing.

    // Get addressing
    const unallocLabelList& upperAddr = fineLevel_.owner();
    const unallocLabelList& lowerAddr = fineLevel_.neighbour();

    const label nFineFaces = fineLevel_.nInternalFaces();

    // Storage for cell neighbours and addressing

    // Guess initial maximum number of neighbours in block
    label maxNnbrs = 10;

    // Number of neighbours per block
    labelList cellNnbrs(nCells_, 0);

    // Setup initial packed storage for neighbours
    labelList cellNbrsData(maxNnbrs*nCells_);

    // Set face-restriction addressing
    faceRestrictAddr_.setSize(nFineFaces);

    // Initial neighbour array (not in upper-triangle order)
    labelList initCoarseNeighb(nFineFaces);

    // Counter for coarse faces
    label nCoarseFaces = 0;

    // Loop through all fine faces
    forAll (upperAddr, fineFaceI)
    {
        label rmUpperAddr = child[upperAddr[fineFaceI]];
        label rmLowerAddr = child[lowerAddr[fineFaceI]];

        if (rmUpperAddr == rmLowerAddr)
        {
            // For each fine face inside of a coarse cluster keep the address
            // of the cluster corresponding to the face in the
            // faceRestrictAddr_ as a negative index
            faceRestrictAddr_[fineFaceI] = -(rmUpperAddr + 1);
        }
        else
        {
            // This face is a part of a coarse face

            label cOwn = rmUpperAddr;
            label cNei = rmLowerAddr;

            // Get coarse owner and neighbour
            if (rmUpperAddr > rmLowerAddr)
            {
                cOwn = rmLowerAddr;
                cNei = rmUpperAddr;
            }

            // Check the neighbour to see if this face has already been found
            bool nbrFound = false;
            label& ccnFaces = cellNnbrs[cOwn];

            for (int i = 0; i < ccnFaces; i++)
            {
                if (initCoarseNeighb[cellNbrsData[maxNnbrs*cOwn + i]] == cNei)
                {
                    nbrFound = true;
                    faceRestrictAddr_[fineFaceI] =
                        cellNbrsData[maxNnbrs*cOwn + i];
                    break;
                }
            }

            if (!nbrFound)
            {
                if (ccnFaces >= maxNnbrs)
                {
                    // Double the size of list and copy data
                    label oldMaxNnbrs = maxNnbrs;
                    maxNnbrs *= 2;

                    // Resize and copy list
                    const labelList oldCellNbrsData = cellNbrsData;
                    cellNbrsData.setSize(maxNnbrs*nCells_);

                    forAll (cellNnbrs, i)
                    {
                        for (int j = 0; j < cellNnbrs[i]; j++)
                        {
                            cellNbrsData[maxNnbrs*i + j] =
                                oldCellNbrsData[oldMaxNnbrs*i + j];
                        }
                    }
                }

                cellNbrsData[maxNnbrs*cOwn + ccnFaces] = nCoarseFaces;
                initCoarseNeighb[nCoarseFaces] = cNei;
                faceRestrictAddr_[fineFaceI] = nCoarseFaces;
                ccnFaces++;

                // New coarse face created
                nCoarseFaces++;
            }
        }
    } // End for all fine faces


    // Renumber into upper-triangular order

    // All coarse owner-neighbour storage
    nInternalFaces_ = nCoarseFaces;
    owner_.setSize(nCoarseFaces);
    neighbour_.setSize(nCoarseFaces);
    labelList coarseFaceMap(nCoarseFaces);

    label coarseFaceI = 0;

    forAll (cellNnbrs, cci)
    {
        label* cFaces = &cellNbrsData[maxNnbrs*cci];
        label ccnFaces = cellNnbrs[cci];

        for (int i = 0; i < ccnFaces; i++)
        {
            owner_[coarseFaceI] = cci;
            neighbour_[coarseFaceI] = initCoarseNeighb[cFaces[i]];
            coarseFaceMap[cFaces[i]] = coarseFaceI;
            coarseFaceI++;
        }
    }

    forAll(faceRestrictAddr_, fineFaceI)
    {
        if (faceRestrictAddr_[fineFaceI] >= 0)
        {
            faceRestrictAddr_[fineFaceI] =
                coarseFaceMap[faceRestrictAddr_[fineFaceI]];
        }
    }

    // Assemble patch faceCells.  Note: no coarsening on boundary
    faceCells_.setSize(nPatches());

    labelList finestChild = fineLevel().prolongToFinest(child);

    for (label patchI = 0; patchI < nPatches(); patchI++)
    {
        const labelList& finestPatchFaceCells =
            finestLevel().faceCells(patchI);

        labelList& curPatchFaceCells = faceCells_[patchI];

        curPatchFaceCells.setSize(finestPatchFaceCells.size());

        forAll (curPatchFaceCells, faceI)
        {
            curPatchFaceCells[faceI] =
                finestChild[finestPatchFaceCells[faceI]];
        }
    }

    // Calculate geometry. Beware: order is important
    cellCentres_.setSize(nCells_);
    faceCentres_.setSize(nInternalFaces_);
    cellVolumes_.setSize(nCells_);
    faceAreas_.setSize(nInternalFaces_);

    // Cell volumes, face areas
    fineLevel_.restrict(cellVolumes_, fineLevel_.cellVolumes());
    restrictFace(faceAreas_, fineLevel_.faceAreas());

    // Magnitude of face area calculated from coarse area vector
    // Do not use agglomeration
    magFaceAreas_ = mag(faceAreas_);

    // Calculate aglomerated cell centres
    vectorField agglomFineCellCentres =
        fineLevel_.cellCentres()*fineLevel_.cellVolumes();

    fineLevel_.restrict(cellCentres_, agglomFineCellCentres);
    cellCentres_ /= cellVolumes_;

    // Calculate agglomerated face centres
    vectorField agglomFineFaceCentres =
        fineLevel_.faceCentres()*fineLevel_.magFaceAreas();

    restrictFace(faceCentres_, agglomFineFaceCentres);

    scalarField agglomFineFaceAreas(nInternalFaces_);
    restrictFace(agglomFineFaceAreas, fineLevel_.magFaceAreas());

    faceCentres_ /= agglomFineFaceAreas;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coarseMgMeshLevel::coarseMgMeshLevel(const mgMeshLevel& fineLevel)
:
    fineLevel_(fineLevel),
    // Sizes
    nCells_(0),
    nInternalFaces_(0),
    // Addressing
    owner_(),
    neighbour_(),
    // Geometrical data
    cellCentres_(),
    faceCentres_(),
    cellVolumes_(),
    faceAreas_(),
    magFaceAreas_(),
    faceRestrictAddr_()
{
    calcLevel();
}


// ************************************************************************* //

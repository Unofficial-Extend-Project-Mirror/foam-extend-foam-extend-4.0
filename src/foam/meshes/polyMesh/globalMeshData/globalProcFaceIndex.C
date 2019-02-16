/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "globalProcFaceIndex.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::globalProcFaceIndex::calcFaceIndex()
{
    // Count number of unique faces on this processor
    nUniqueFaces_[Pstream::myProcNo()] = mesh_.nInternalFaces();

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Assign unique face label to all master processor faces

    // Count faces and processor faces per processor
    forAll (patches, patchI)
    {
        // Only skip slave processor patches
        if (isA<processorPolyPatch>(patches[patchI]))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);

            if (procPatch.master())
            {
                // Found unique faces
                nUniqueFaces_[Pstream::myProcNo()] += procPatch.size();
            }
            // else Slave processor patch.  Skip it
        }
        else
        {
            // Regular patch.  Add faces to count
            nUniqueFaces_[Pstream::myProcNo()] += patches[patchI].size();
        }
    }

    // Gather data to master processor
    Pstream::gatherList(nUniqueFaces_);
    Pstream::scatterList(nUniqueFaces_);

    // Adjust all lists to calculate offsets on the master processor only
    if (Pstream::master())
    {
        for (label procI = 1; procI < procFaceOffset_.size(); procI++)
        {
            // Number of unique faces for this processor is equal to the number
            // of faces on previous processor + number of local unique faces
            procFaceOffset_[procI] =
                procFaceOffset_[procI - 1] + nUniqueFaces_[procI - 1];
        }
    }

    // Scatter offset data to all processors
    Pstream::scatter(procFaceOffset_);

    // Assemble global label list for mesh faces
    label globalFaceIndex = procFaceOffset_[Pstream::myProcNo()];

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        globalLabel_[faceI] = globalFaceIndex;
        globalFaceIndex++;
    }

    // Assign and send from master
    forAll (patches, patchI)
    {
        const label patchSize = patches[patchI].size();
        const label patchStart = patches[patchI].start();

        if (isA<processorPolyPatch>(patches[patchI]))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);

            if (procPatch.master())
            {
                // Master processor patch: assign face index
                for
                (
                    label patchFaceI = patchStart;
                    patchFaceI < patchStart + patchSize;
                    patchFaceI++
                )
                {
                    globalLabel_[patchFaceI] = globalFaceIndex;
                    globalFaceIndex++;
                }

                // Make a slice and send it
                labelList curFaceLabels(procPatch.patchSlice(globalLabel_));

                OPstream toProc
                (
                    // HR 12.12.18: nonBlocking fails on PLB of Aachen bomb
                    Pstream::blocking,
                    procPatch.neighbProcNo()
                );
                toProc<< curFaceLabels;
            }
        }
        else
        {
            // Regular patch: assign face index
            for
            (
                label patchFaceI = patchStart;
                patchFaceI < patchStart + patchSize;
                patchFaceI++
            )
            {
                globalLabel_[patchFaceI] = globalFaceIndex;
                globalFaceIndex++;
            }
        }
    }

    // Receive on slave
    forAll (patches, patchI)
    {
        if (isA<processorPolyPatch>(patches[patchI]))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);

            if (!procPatch.master())
            {
                // Slave processor patch
                // Receive the data from master and insert into the list
                IPstream fromProc
                (
                    // HR 12.12.18: nonBlocking fails on PLB of Aachen bomb
                    Pstream::blocking,
                    procPatch.neighbProcNo()
                );

                labelList masterFaceLabels(fromProc);

                // Insert the data into the list
                const label patchStart = patches[patchI].start();

                forAll (masterFaceLabels, patchFaceI)
                {
                    globalLabel_[patchStart + patchFaceI] =
                        masterFaceLabels[patchFaceI];
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalProcFaceIndex::globalProcFaceIndex(const polyMesh& mesh)
:
    mesh_(mesh),
    nUniqueFaces_(Pstream::nProcs(), 0),
    procFaceOffset_(Pstream::nProcs(), 0),
    globalLabel_(mesh_.nFaces())
{
    calcFaceIndex();
}


// ************************************************************************* //

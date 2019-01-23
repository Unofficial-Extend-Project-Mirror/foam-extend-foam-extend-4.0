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

#include "globalProcPointIndex.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::globalProcPointIndex::calcPointIndex()
{
    // Count and mark points in the following order
    // 1. global points
    // 2. processor patches
    // 3. regular patches
    // 4. anything else ("internal points")
    //
    // The marking is a follows
    // -4 : Not visited yet
    // -3 : Found in a regular patch
    // -2 : Found in slave processor patch
    // -1 : Found in a master processor patch
    //  0 - nGlobalPoints-1 : Global point ID

    // 1. Mark the global points
    const labelList& spl = mesh_.globalData().sharedPointLabels();
    const labelList& spa = mesh_.globalData().sharedPointAddr();

    forAll(spl, spI)
    {
        globalLabel_[spl[spI]] = spa[spI];
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelList procPointCounter(Pstream::nProcs(), 0);
    labelList patchPointCounter(patches.size(), 0);
    label& nProcPoints = procPointCounter[Pstream::myProcNo()];

    // Master processor patches - first pass
    forAll (patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const labelList& patchMeshPoints = pp.meshPoints();

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            label& nPatchPoints = patchPointCounter[patchI];

            if (procPatch.master())
            {
                forAll (patchMeshPoints, mpI)
                {
                    label& pl = globalLabel_[patchMeshPoints[mpI]];

                    if (pl < 0)
                    {
                        nProcPoints++;
                        nPatchPoints++;
                        pl = -1;
                    }
                }
            }
        }
    }

    // Slave processor patches - first pass
    forAll (patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const labelList& patchMeshPoints = pp.meshPoints();

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            if (!procPatch.master())
            {
                forAll (patchMeshPoints, mpI)
                {
                    label& pl = globalLabel_[patchMeshPoints[mpI]];

                    if (pl == -4)
                    {
                        pl = -2;
                    }
                }
            }
        }
    }

    // Regular patches - first pass
    forAll (patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const labelList& patchMeshPoints = pp.meshPoints();

        if (!isA<processorPolyPatch>(pp))
        {
            label& nPatchPoints = patchPointCounter[patchI];

            forAll (patchMeshPoints, mpI)
            {
                label& pl = globalLabel_[patchMeshPoints[mpI]];

                if (pl == -4)
                {
                    nProcPoints++;
                    nPatchPoints++;
                    pl = -3;
                }
            }
        }
    }

    // "Internal" points - first pass
    label nInternalPoints = 0;
    forAll (globalLabel_, pI)
    {
        if (globalLabel_[pI] == -4)
        {
            nInternalPoints++;
            nProcPoints++;
        }
    }

    // Gather-Scatter counter data
    Pstream::gatherList(procPointCounter);
    Pstream::scatterList(procPointCounter);

    // Convert counters to offsets
    procPointOffset_[0] = mesh_.globalData().nGlobalPoints();
    for (label procI = 1; procI < procPointOffset_.size(); procI++)
    {
        procPointOffset_[procI] =
            procPointCounter[procI-1] + procPointOffset_[procI-1];
    }

    patchPointOffset_[0] =
        procPointOffset_[Pstream::myProcNo()] + nInternalPoints;
    for (label patchI = 1; patchI < patchPointOffset_.size(); patchI++)
    {
        patchPointOffset_[patchI] =
            patchPointCounter[patchI-1] + patchPointOffset_[patchI-1];
    }

    const pointField& points = mesh_.points();

    // Master processor patches - second pass
    // Send patch offset to slave and assign global labels to points on
    // master processor patches and regular patches
    forAll (patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const labelList& patchMeshPoints = pp.meshPoints();
        const label patchPO = patchPointOffset_[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            if (procPatch.master())
            {
                // Master processor patch
                // Send the offset to the slave and mark the points through a
                // face loop to establish an order that can be untangled on
                // slave side

                OPstream toProc
                (
                    Pstream::nonBlocking,
                    procPatch.neighbProcNo()
                );
                toProc<< patchPO;

                pointField pointLocs(patchMeshPoints.size());
                label glPointLabelsI = 0;
                forAll (pp, faceI)
                {
                    const face& curFace = pp[faceI];

                    forAll (curFace, fpI)
                    {
                        const label pointI = curFace[fpI];
                        label& pl = globalLabel_[pointI];

                        if (pl == -1)
                        {
                            pl = patchPO + glPointLabelsI;
                            pointLocs[glPointLabelsI] = points[pointI];
                            glPointLabelsI++;
                        }
                    }
                }
            }
        }
    }

    // Slave processor and regular patches - second pass
    // Receive data on slave and assign global labels to points on
    // master processor patches and regular patches
    forAll (patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            if (!procPatch.master())
            {
                // Slave processor patch - second pass
                // Receive the offset from master and mark the points by taking
                // into account that the points in the faces are in reverse
                // order

                IPstream fromProc
                (
                    Pstream::nonBlocking,
                    procPatch.neighbProcNo()
                );

                label masterPatchPO(readLabel(fromProc));

                label glPointLabelsI = 0;
                forAll (pp, faceI)
                {
                    const face curFace = pp[faceI].reverseFace();

                    forAll (curFace, fpI)
                    {
                        const label pointI = curFace[fpI];
                        label& pl = globalLabel_[pointI];

                        if (pl == -2)
                        {
                            pl = masterPatchPO + glPointLabelsI;
                            glPointLabelsI++;
                        }
                    }
                }
            }
        }
        else
        {
            // Regular patch - second pass
            const labelList& patchMeshPoints = pp.meshPoints();
            const label patchPO = patchPointOffset_[patchI];
            label glPointLabelsI = 0;
            forAll (patchMeshPoints, mpI)
            {
                label& pl = globalLabel_[patchMeshPoints[mpI]];

                if (pl == -3)
                {
                    pl = patchPO + glPointLabelsI;
                    glPointLabelsI++;
                }
            }
        }
    }

    // "Internal" points - second pass
    nInternalPoints = 0;
    forAll (globalLabel_, pI)
    {
        label& pl = globalLabel_[pI];

        if (pl == -4)
        {
            pl = procPointOffset_[Pstream::myProcNo()] + nInternalPoints;
            nInternalPoints++;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalProcPointIndex::globalProcPointIndex(const polyMesh& mesh)
:
    mesh_(mesh),
    procPointOffset_(Pstream::nProcs()+1),
    patchPointOffset_(mesh_.boundaryMesh().size()+1),
    globalLabel_(mesh_.nPoints(), -4)
{
    calcPointIndex();
}


// ************************************************************************* //

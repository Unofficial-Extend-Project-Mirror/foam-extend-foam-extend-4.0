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

#include "sharedPoints.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelListList Foam::sharedPoints::procPatchPairs() const
{
    labelListList patchPairs(meshes_.size());

    // Initialise patch pair indices to -1
    forAll (meshes_, meshI)
    {
        if (meshes_.set(meshI))
        {
            patchPairs[meshI].setSize(meshes_[meshI].boundaryMesh().size(), -1);
        }
    }

    forAll (meshes_, meshI)
    {
        if (meshes_.set(meshI))
        {
            const polyMesh& curMesh = meshes_[meshI];
            const polyBoundaryMesh& curPatches = curMesh.boundaryMesh();

            forAll (curPatches, patchI)
            {
                if (isA<processorPolyPatch>(curPatches[patchI]))
                {
                    // Found processor patch
                    if (patchPairs[meshI][patchI] == -1)
                    {
                        // Neighbour not found.  Find one
                        const processorPolyPatch& myProcPatch =
                            refCast<const processorPolyPatch>
                            (
                                curPatches[patchI]
                            );

                        const int myProcID = meshI;
                        const int nbrProcID = myProcPatch.neighbProcNo();

                        // Get the other mesh
                        if (!meshes_.set(nbrProcID))
                        {
                            FatalErrorIn
                            (
                                "labelListList sharedPoints::procPatchPairs()"
                            )   << "Neighbour mesh does not exist for proc "
                                << meshI << " patch " << patchI
                                << " and neighbour " << nbrProcID
                                << abort(FatalError);
                        }

                        const polyMesh& nbrMesh = meshes_[nbrProcID];
                        const polyBoundaryMesh& nbrPatches =
                            nbrMesh.boundaryMesh();

                        // Check all neighbour processor patches until a match
                        // is found
                        bool found = false;

                        forAll (nbrPatches, nbrPatchI)
                        {
                            if (isA<processorPolyPatch>(nbrPatches[nbrPatchI]))
                            {
                                const processorPolyPatch& nbrProcPatch =
                                    refCast<const processorPolyPatch>
                                    (
                                        nbrPatches[nbrPatchI]
                                    );

                                if (nbrProcPatch.neighbProcNo() == myProcID)
                                {
                                    // Pair found.  Record it twice
                                    patchPairs[myProcID][patchI] = nbrPatchI;

                                    patchPairs[nbrProcID][nbrPatchI] = patchI;

                                    found = true;
                                }
                            }

                            if (found)
                            {
                                break;
                            }
                        }

                        if (!found)
                        {
                            FatalErrorIn
                            (
                                "labelListList sharedPoints::procPatchPairs()"
                            )   << "Neighbour patch does not exist for proc "
                                << meshI << " patch " << patchI
                                << " and neighbour " << nbrProcID
                                << abort(FatalError);
                        }
                    }
                }
            }
        }
    }

    return patchPairs;
}


void Foam::sharedPoints::syncMark
(
    labelListList& markedPoints,
    const labelListList& patchPairs,
    const label fromMesh
) const
{
    label nSynced;

    do
    {
        nSynced = 0;

        // Sync mark across processor boundaries
        for (label meshI = fromMesh; meshI < meshes_.size(); meshI++)
        {
            if (meshes_.set(meshI))
            {
                const polyMesh& curMesh = meshes_[meshI];
                const polyBoundaryMesh& patches = curMesh.boundaryMesh();

                labelList& curMarkedPoints = markedPoints[meshI];

                // Collect the points that have been addressed multiple times
                forAll (patches, patchI)
                {
                    const polyPatch& pp = patches[patchI];

                    if (isA<processorPolyPatch>(pp))
                    {
                        // My processor patch
                        const processorPolyPatch& myProcPatch =
                            refCast<const processorPolyPatch>(pp);

                        // Get local mesh points
                        const labelList& patchMeshPoints = pp.meshPoints();

                        // Get my local faces
                        const faceList& patchLocalFaces = pp.localFaces();

                        // Get neighbour processor ID
                        const int nbrProcID = myProcPatch.neighbProcNo();

                        // Neighbour patch
                        const polyPatch& nbrPatch =
                            meshes_[nbrProcID].boundaryMesh()
                            [patchPairs[meshI][patchI]];

                        // Get neighbour mesh points
                        const labelList& nbrMeshPoints = nbrPatch.meshPoints();

                        // Get neighbour local faces
                        const faceList& nbrLocalFaces = nbrPatch.localFaces();

                        labelList& nbrMarkedPoints = markedPoints[nbrProcID];

                        // Note:
                        // Cannot loop over mesh points because they are sorted
                        // Use face loops as they are synchronised
                        // HJ, 21/May/2018
                        forAll (patchLocalFaces, faceI)
                        {
                            const face& curFace = patchLocalFaces[faceI];

                            // Reverse neighbour face (copy)
                            const face nbrFace =
                                nbrLocalFaces[faceI].reverseFace();

                            forAll (curFace, fpI)
                            {
                                const label patchMpI =
                                    patchMeshPoints[curFace[fpI]];

                                const label nbrMpI =
                                    nbrMeshPoints[nbrFace[fpI]];

                                if
                                (
                                    curMarkedPoints[patchMpI]
                                 != nbrMarkedPoints[nbrMpI]
                                )
                                {
                                    const label maxMark =
                                        Foam::max
                                        (
                                            curMarkedPoints[patchMpI],
                                            nbrMarkedPoints[nbrMpI]
                                        );

                                    nbrMarkedPoints[nbrMpI] = maxMark;
                                    curMarkedPoints[patchMpI] = maxMark;

                                    nSynced++;
                                }
                            }
                        }
                    }
                }
            }
        }
    } while (nSynced > 0);
}


void Foam::sharedPoints::calcSharedPoints()
{
    // Algorithm
    // Go through all processor patches and mark local points that are used
    // by more than one processor patch and mark them as globally shared
    // Pass the data to other processors.  Mark the locally multiply shared
    // points and pass on the data
    // Once all the data is passed forwards and back, check all points on
    // all processors.  Record globally shared point, its local label and its
    // slot in the globally shared point list

    // Mark-up:
    // 0 = point does not touch a processor boundary
    // 1 = point on only one processor boundary: not locally shared
    // 2 = locally detected global point

    // Mark-up array: procI, procJ,
    labelListList markedPoints(meshes_.size());

    forAll (meshes_, meshI)
    {
        if (meshes_.set(meshI))
        {
            markedPoints[meshI].setSize(meshes_[meshI].nPoints(), 0);
        }
    }

    // Mark up points for the first time
    forAll (meshes_, meshI)
    {
        if (meshes_.set(meshI))
        {
            const polyMesh& curMesh = meshes_[meshI];
            const polyBoundaryMesh& patches = curMesh.boundaryMesh();

            // Mark points belonging to processor patches.  If the point
            // is marked more than once, it may be a globally shared point
            labelList& curMarkedPoints = markedPoints[meshI];

            forAll (patches, patchI)
            {
                const polyPatch& pp = patches[patchI];

                if (isA<processorPolyPatch>(pp))
                {
                    // Found processor patch
                    const labelList& patchMeshPoints = pp.meshPoints();

                    forAll (patchMeshPoints, mpI)
                    {
                        // Mark the point
                        curMarkedPoints[patchMeshPoints[mpI]]++;
                    }
                }
            }
        }
    }

    // Get processor patch to neighbour processor patch addressing
    labelListList patchPairs = procPatchPairs();

    // Communicate and count global points

    // Identify, count and communicate points across processor boundaries
    // Repeat until the number of points per processor stabilises,
    // ie. no further points are found through communication

    // syncMark across all processors
    syncMark(markedPoints, patchPairs, 0);

    // Grab marked points
    forAll (meshes_, meshI)
    {
        if (meshes_.set(meshI))
        {
            // Count the points and then collect
            label nShared = 0;

            const labelList& curMarkedPoints = markedPoints[meshI];

            forAll (curMarkedPoints, pointI)
            {
                if (curMarkedPoints[pointI] > 1)
                {
                    nShared++;
                }
            }

            labelList& curShared = sharedPointLabels_[meshI];
            curShared.setSize(nShared);

            // Re-use the counter for insertion
            nShared = 0;

            forAll (curMarkedPoints, pointI)
            {
                if (curMarkedPoints[pointI] > 1)
                {
                    curShared[nShared] = pointI;
                    nShared++;
                }
            }
        }
    }

    // Clear markup list.  It will be used for the global processor
    // point index
    forAll (markedPoints, meshI)
    {
        markedPoints[meshI] = -1;
    }

    // Provide global mark for all processors and communicate it across
    // processor boundaries
    forAll (meshes_, meshI)
    {
        if (meshes_.set(meshI))
        {
            // Get shared points and assign global shared point index
            const labelList& curSharedPoints = sharedPointLabels_[meshI];

            // Prepare addressing into the global shared point list
            labelList& curSharedAddr = sharedPointAddr_[meshI];
            curSharedAddr.setSize(curSharedPoints.size());

            labelList& curMarkedPoints = markedPoints[meshI];

            // Collect shared point addressing
            forAll (curSharedPoints, spI)
            {
                if (curMarkedPoints[curSharedPoints[spI]] == -1)
                {
                    // Found new point.  Mark it and collect addressing
                    curMarkedPoints[curSharedPoints[spI]] = nGlobalPoints_;

                    // Collect addressing
                    curSharedAddr[spI] = nGlobalPoints_;

                    nGlobalPoints_++;
                }
                else
                {
                    // Point already marked.  Collect addressing
                    curSharedAddr[spI] = curMarkedPoints[curSharedPoints[spI]];
                }
            }

            // Communicate labels accross the boundary using processor patches
            // from this processor to all higher
            syncMark(markedPoints, patchPairs, meshI);
        }
    }

    // debug
    {
        pointField gp(nGlobalPoints_);
        boolList gpSet(nGlobalPoints_, false);

        forAll (meshes_, meshI)
        {
            if (meshes_.set(meshI))
            {
                // Get points
                const pointField& P = meshes_[meshI].points();

                // Get list of local point labels that are globally shared
                const labelList& curShared = sharedPointLabels_[meshI];

                // Get index in global point list
                const labelList& curAddr = sharedPointAddr_[meshI];

                // Loop through all local points
                forAll (curShared, i)
                {
                    if (!gpSet[curAddr[i]])
                    {
                        // Point not set: set it
                        gp[curAddr[i]] = P[curShared[i]];
                        gpSet[curAddr[i]] = true;
                    }
                    else
                    {
                        // Point already set: check location
                        if (mag(gp[curAddr[i]] - P[curShared[i]]) > SMALL)
                        {
                            Info<< "MERGE MISMATCH: mesh" << meshI
                                << " point: " << curShared[i]
                                << " dist: " << gp[curAddr[i]] << " "
                                << P[curShared[i]]
                                << endl;
                        }
                    }
                }
            }
        }

        // Grab marked points
        OFstream ppp("points.vtk");
        ppp << "# vtk DataFile Version 2.0" << nl
            << "points.vtk" << nl
            << "ASCII" << nl
            << "DATASET POLYDATA" << nl
            << "POINTS " << nGlobalPoints_ << " float" << nl;

        forAll (gp, i)
        {
            ppp << float(gp[i].x()) << ' '
                << float(gp[i].y()) << ' '
                << float(gp[i].z())
                << nl;
        }
    }
}


void Foam::sharedPoints::calcSharedPoints
(
    const labelListList& globalPointIndex
)
{
    typedef HashTable<Pair<label>, label, Hash<label> > markTable;
    markTable marker;

    List<labelHashSet> procSharedPoints(meshes_.size());

    // Mark up points for the first time
    forAll (meshes_, meshI)
    {
        if (meshes_.set(meshI))
        {
            const labelList& curGlobalPointIndex = globalPointIndex[meshI];

            labelHashSet& sharedPoints = procSharedPoints[meshI];

            forAll (curGlobalPointIndex, pI)
            {
                const label gpi = curGlobalPointIndex[pI];

                markTable::iterator iter = marker.find(gpi);
                if (iter == Map<label>::end())
                {
                    // Record meshI and pI
                    marker.insert(gpi, Pair<label>(-meshI - 1, pI));
                }
                else
                {
                    if (iter().first() < 0)
                    {
                        if(iter().first() != -meshI - 1)
                        {
                            // Shared point detected. Add for both meshes

                            // Add for first mesh
                            const label firstMesh = -iter().first() - 1;
                            const label firstPointI = iter().second();
                            procSharedPoints[firstMesh].insert(firstPointI);

                            // Add for current mesh
                            sharedPoints.insert(pI);

                            // Count shared points and update bookkeeping
                            iter().first() = nGlobalPoints_;
                            nGlobalPoints_++;

                        }
                    }
                    else
                    {
                        // Existing shared point. Add for current mesh
                        sharedPoints.insert(pI);
                    }
                }
            }
        }
    }

    forAll (meshes_, meshI)
    {
        if (meshes_.set(meshI))
        {
            const labelList& curGlobalPointIndex = globalPointIndex[meshI];

            labelList& curSharedPoints = sharedPointLabels_[meshI];
            curSharedPoints = procSharedPoints[meshI].toc();

            labelList& curSharedAddr = sharedPointAddr_[meshI];
            curSharedAddr.setSize(curSharedPoints.size());

            forAll(curSharedPoints, i)
            {
                const label gpi = curSharedPoints[i];
                curSharedAddr[i] = marker.find(curGlobalPointIndex[gpi])().first();
            }
        }
    }

    // debug
    {
        pointField gp(nGlobalPoints_);
        boolList gpSet(nGlobalPoints_, false);

        forAll (meshes_, meshI)
        {
            if (meshes_.set(meshI))
            {
                // Get points
                const pointField& P = meshes_[meshI].points();

                // Get list of local point labels that are globally shared
                const labelList& curShared = sharedPointLabels_[meshI];

                // Get index in global point list
                const labelList& curAddr = sharedPointAddr_[meshI];

                // Loop through all local points
                forAll (curShared, i)
                {
                    if (!gpSet[curAddr[i]])
                    {
                        // Point not set: set it
                        gp[curAddr[i]] = P[curShared[i]];
                        gpSet[curAddr[i]] = true;
                    }
                    else
                    {
                        // Point already set: check location
                        if (mag(gp[curAddr[i]] - P[curShared[i]]) > SMALL)
                        {
                            Info<< "MERGE MISMATCH: mesh" << meshI
                                << " point: " << curShared[i]
                                << " dist: " << gp[curAddr[i]] << " "
                                << P[curShared[i]]
                                << endl;
                        }
                    }
                }
            }
        }

    /*
            // Grab marked points
            OFstream ppp("points.vtk");
            ppp << "# vtk DataFile Version 2.0" << nl
                << "points.vtk" << nl
                << "ASCII" << nl
                << "DATASET POLYDATA" << nl
                << "POINTS " << nGlobalPoints_ << " float" << nl;

            Pout << "Global marking " << nGlobalPoints_ << endl;
            forAll (gp, i)
            {
                ppp << float(gp[i].x()) << ' '
                    << float(gp[i].y()) << ' '
                    << float(gp[i].z())
                    << nl;

                Pout << gp[i] << endl;
            }
            */
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sharedPoints::sharedPoints
(
    const PtrList<fvMesh>& meshes
)
:
    meshes_(meshes),
    sharedPointAddr_(meshes_.size()),
    sharedPointLabels_(meshes_.size()),
    nGlobalPoints_(0)
{
    calcSharedPoints();
}

Foam::sharedPoints::sharedPoints
(
    const PtrList<fvMesh>& meshes,
    const labelListList& globalPointIndex
)
:
    meshes_(meshes),
    sharedPointAddr_(meshes_.size()),
    sharedPointLabels_(meshes_.size()),
    nGlobalPoints_(0)
{
    calcSharedPoints(globalPointIndex);
}


// ************************************************************************* //

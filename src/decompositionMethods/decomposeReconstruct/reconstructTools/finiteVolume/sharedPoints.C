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
    Info<< "patch pairs: " << patchPairs << endl;
    return patchPairs;
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
    labelList nGlobalPointsPerProc(meshes_.size(), 0);

    // Identify, count and communicate points across processor boundaries
    // Repeat until the number of points per processor stabilises,
    // ie. no further points are found through communication
    label oldNTotalPoints, newNTotalPoints;
    do
    {
        oldNTotalPoints = sum(nGlobalPointsPerProc);

        // Reset the list
        nGlobalPointsPerProc = 0;

        forAll (meshes_, meshI)
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
                        // Found processor patch
                        const labelList& patchMeshPoints = pp.meshPoints();

                        // My processor patch
                        const processorPolyPatch& myProcPatch =
                            refCast<const processorPolyPatch>(pp);

                        // Get neighbour processor ID
                        const int nbrProcID = myProcPatch.neighbProcNo();

                        // Neighbour patch
                        const polyPatch& nbrPatch =
                            meshes_[nbrProcID].boundaryMesh()
                            [patchPairs[meshI][patchI]];

                        const labelList& nbrMeshPoints = nbrPatch.meshPoints();

                        forAll (patchMeshPoints, mpI)
                        {
                            if (curMarkedPoints[patchMeshPoints[mpI]] > 1)
                            {
                                // Mark the point on the other processor/side
                                markedPoints[nbrProcID][nbrMeshPoints[mpI]] =
                                    curMarkedPoints[patchMeshPoints[mpI]];
                            }
                        }
                    }
                }

                // Count number of shared points per processor
                forAll (curMarkedPoints, cpI)
                {
                    if (curMarkedPoints[cpI] > 1)
                    {
                        nGlobalPointsPerProc[meshI]++;
                    }
                }
            }
        }

        newNTotalPoints = sum(nGlobalPointsPerProc);

        Info<< "Proc merge pass: " << oldNTotalPoints << " "
            << newNTotalPoints << endl;
    } while (oldNTotalPoints != newNTotalPoints);

    Info<< "Number of shared points per processor: " << nGlobalPointsPerProc
        << endl;

    // Collect points for every processor, in order to re-use the markedPoints
    // list.  Note: the list of global labels of shared points
    // will be collected later
    forAll (meshes_, meshI)
    {
        if (meshes_.set(meshI))
        {
            labelList& curSharedPoints = sharedPointLabels_[meshI];

            curSharedPoints.setSize(nGlobalPointsPerProc[meshI]);

            // Count inserted points
            label nShared = 0;

            // Get point marking
            const labelList& curMarkedPoints = markedPoints[meshI];

            forAll (curMarkedPoints, pointI)
            {
                if (curMarkedPoints[pointI] > 1)
                {
                    curSharedPoints[nShared] = pointI;
                    nShared++;
                }
            }
        }
    }

    // Clear markup list.  It will be used for the global processor point
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
            const polyMesh& curMesh = meshes_[meshI];
            const polyBoundaryMesh& patches = curMesh.boundaryMesh();

            // Get shared points and assign global shared point index
            const labelList& curSharedPoints = sharedPointLabels_[meshI];

            // Prepare addressing into the global shared point list
            labelList& curSharedAddr = sharedPointAddr_[meshI];
            curSharedAddr.setSize(curSharedPoints.size());

            labelList& curMarkedPoints = markedPoints[meshI];

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
            forAll (patches, patchI)
            {
                const polyMesh& curMesh = meshes_[meshI];
                const polyBoundaryMesh& patches = curMesh.boundaryMesh();

                // Get point marking
                const labelList& curMarkedPoints = markedPoints[meshI];

                const polyPatch& pp = patches[patchI];

                if (isA<processorPolyPatch>(pp))
                {
                    // Found processor patch

                    // My processor patch
                    const processorPolyPatch& myProcPatch =
                        refCast<const processorPolyPatch>(pp);

                    const labelList& patchMeshPoints = pp.meshPoints();

                    // Get neighbour processor ID
                    const int nbrProcID = myProcPatch.neighbProcNo();

                    // Neighbour patch
                    const polyPatch& nbrPatch =
                        meshes_[nbrProcID].boundaryMesh()
                        [patchPairs[meshI][patchI]];

                    const labelList& nbrMeshPoints = nbrPatch.meshPoints();

                    forAll (patchMeshPoints, mpI)
                    {
                        if (curMarkedPoints[patchMeshPoints[mpI]] > -1)
                        {
                            // Mark opposite side
                            markedPoints[nbrProcID][nbrMeshPoints[mpI]] =
                                curMarkedPoints[patchMeshPoints[mpI]];
                        }
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sharedPoints::sharedPoints(const PtrList<fvMesh>& meshes)
:
    meshes_(meshes),
    sharedPointAddr_(meshes_.size()),
    sharedPointLabels_(meshes_.size()),
    nGlobalPoints_(0)
{
    calcSharedPoints();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Foam::sharedPoints::~sharedPoints()
// {}


// ************************************************************************* //

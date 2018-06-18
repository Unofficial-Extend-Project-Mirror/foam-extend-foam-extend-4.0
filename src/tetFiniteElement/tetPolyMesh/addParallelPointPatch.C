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

Description

\*---------------------------------------------------------------------------*/

#include "tetPolyMesh.H"
#include "foamTime.H"
#include "PstreamCombineReduceOps.H"
#include "processorTetPolyPatch.H"
#include "globalTetPolyPatch.H"
#include "Pstream.H"
#include "labelIOList.H"
#include "globalMeshData.H"

using namespace Foam;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Global patch table combination class
template<class Type>
class globalPatchTableCombine
{

public:

    typedef HashTable<label, Type, Hash<Type> > HashTableType;

    void operator()
    (
        HashTableType& globalObject,
        const HashTableType& myObject
    ) const
    {
        label mySize = globalObject.size();

        forAllConstIter (typename HashTableType, myObject, iter)
        {
            const Type& myCurType = iter.key();

            if (!globalObject.found(myCurType))
            {
                globalObject.insert(myCurType, mySize++);
            }
        }
    }
};


void Foam::tetPolyMesh::addParallelPointPatch()
{
    // Note:
    // The processor point patch will be added if processor boundaries
    // exist in the case.  If the mesh with processor boundaries is
    // not created during a parallel run (e.g. decomposePar), the
    // addressing will be dummy.  HJ, date deleted

    if (mesh_.globalData().parallel())
    {
        // Get global point numbering
        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                mesh_.time().findInstance(operator()().meshDir(), "faces"),
                mesh_.meshSubDir,
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                mesh_.time().findInstance(operator()().meshDir(), "faces"),
                mesh_.meshSubDir,
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Adjust the faceProcAddressing for the offset from decomposition
        forAll (faceProcAddressing, faceI)
        {
            faceProcAddressing[faceI] = mag(faceProcAddressing[faceI]) - 1;
        }

        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                mesh_.time().findInstance(operator()().meshDir(), "faces"),
                mesh_.meshSubDir,
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        const label globalFaceOffset = mesh_.globalData().nTotalPoints();
        const label globalCellOffset =
            mesh_.globalData().nTotalFaces() + globalFaceOffset;

        // Add the cell centres to the lookup list

        // Bug fix: only live points/faces/cells included.  HJ, 24/Aug/2012
        const label oldPointSize = mesh_.nPoints();
        const label oldFaceSize = mesh_.nFaces();
        const label oldCellSize = mesh_.nCells();

        pointProcAddressing.setSize
        (
            oldPointSize + oldFaceSize + oldCellSize
        );

        label oldSize = oldPointSize;

        for (label faceI = 0; faceI < oldFaceSize; faceI++)
        {
            pointProcAddressing[oldSize] =
                faceProcAddressing[faceI] + globalFaceOffset;
            oldSize++;
        }

        forAll (cellProcAddressing, cellI)
        {
            pointProcAddressing[oldSize] =
                cellProcAddressing[cellI] + globalCellOffset;
            oldSize++;
        }

        // Do parallel points

        // Get the list of local parallel processor points
        const labelList& localParPoints =
            mesh_.globalData().sharedPointLabels();

        // Do parallel edges

        // Get the list of local parallel processor edges
        const edgeList& localParEdges = parallelEdges();

        // Label to store global parallel edges count
        label nGlobalParEdges = 0;

        // Create the local-to-global edge mapping
        labelList localEdgeMapping(localParEdges.size(), -1);

        if (Pstream::parRun())
        {
            // Global numbers for each of the parallel edges
            HashTable<label, edge, Hash<edge> > globalParEdgesTable;

            // Insert local parallel edges
            forAll (localParEdges, edgeI)
            {
                globalParEdgesTable.insert
                (
                    edge
                    (
                        pointProcAddressing[localParEdges[edgeI].start()],
                        pointProcAddressing[localParEdges[edgeI].end()]
                    ),
                    edgeI
                );
            }

            // Create a global set of edges by reduction from all processors
            combineReduce
            (
                globalParEdgesTable,
                globalPatchTableCombine<edge>()
            );

            // Find out which of the parallel edges are local
            forAll (localParEdges, edgeI)
            {
                edge curGlobalEdge =
                    edge
                    (
                        pointProcAddressing[localParEdges[edgeI].start()],
                        pointProcAddressing[localParEdges[edgeI].end()]
                    );

                localEdgeMapping[edgeI] = globalParEdgesTable[curGlobalEdge];
            }

            // Store global parallel edges count
            nGlobalParEdges = globalParEdgesTable.size();

            // Clear global parallel edge table
            globalParEdgesTable.clear();
        }

        // Debug check
        if (debug)
        {
            if (localEdgeMapping.size() > 0)
            {
                if (min(localEdgeMapping) < 0)
                {
                    FatalErrorIn
                    (
                        "void Foam::tetPolyMesh::"
                        "addParallelPointPatch()"
                    )   << "Error in parallel points patch: edges"
                        << abort(FatalError);
                }
            }
        }

        // Do parallel cut edges

        // Algorithm
        // Parallel cut edges come in two flavours: the ones local to
        // the processor and the ones shared between several
        // processors.  The first lot cause no trouble, but for the
        // second we need to make sure that only one processor
        // multiplies out the global sum part of the matrix.
        // At the same time, the local part of the product needs to be
        // done on every processor seprately.  This is catered for by
        // two mechanisms:
        // 1) Local edges will be calculated here for every patch.
        //    The list needs to be complete.
        // 2) As the lists are calculated, a global list of cut edges
        // is assembled.  The first time the edge is added into the
        // global list, it is accepted for global sum multiplication.
        // If the same global edge is found, the contribution to the
        // global sum is blocked.

        // Count the maximum number of edges coming from the patch
        label maxEdgesOnPatch = 0;

        // For every local point get a list of edges
        forAll (localParPoints, pointI)
        {
            maxEdgesOnPatch += nEdgesForPoint(localParPoints[pointI]);
        }

        edgeList localCutEdges(maxEdgesOnPatch);
        label nCutEdges = 0;

        // Go through all the local points and get all the edges coming
        // from that point.  Check if the edge is local, if not, it is cut

        const unallocLabelList& L = lduAddr().lowerAddr();
        const unallocLabelList& U = lduAddr().upperAddr();

        labelHashSet localParPointsSet(localParPoints);

        forAll (localParPoints, pointI)
        {
            const label curPoint = localParPoints[pointI];

            labelList curEdges = edgesForPoint(curPoint);

            forAll (curEdges, edgeI)
            {
                edge e(L[curEdges[edgeI]], U[curEdges[edgeI]]);
                label otherEnd = -1;

                if (e.start() == curPoint)
                {
                    otherEnd = e.end();
                }
                else
                {
                    otherEnd = e.start();
                }

                // If the other end is not a local point,
                // this is a cut edge
                if (!localParPointsSet.found(otherEnd))
                {
                    // This is a cut edge. Add it to the list
                    localCutEdges[nCutEdges] = e;
                    nCutEdges++;
                }
            }
        }

        // Clear parallel points set
        localParPointsSet.clear();

        // Reset the size of the local cut edge list
        localCutEdges.setSize(nCutEdges);

        // Create a masking array.  For the edges that do not contribute
        // to the global product, the mask will be set to zero
        scalarField localCutEdgeMask(localCutEdges.size(), 1);

        if (Pstream::parRun())
        {
            // Creating the mask
            // Every processor goes through all of its edges and creates the
            // global version of the edge. If such an edge does not
            // already exist, it is added to the list, if it does, its
            // mask is set to zero (somebody is already doing the
            // multiplication)
            if (Pstream::master())
            {
                // Master's mask is always one. Add all edges to the list
                HashSet<edge, Hash<edge> > globalCutEdgesSet;

                forAll (localCutEdges, edgeI)
                {
                    globalCutEdgesSet.insert
                    (
                        edge
                        (
                            pointProcAddressing[localCutEdges[edgeI].start()],
                            pointProcAddressing[localCutEdges[edgeI].end()]
                        )
                    );
                }

                // Send the list to the first slave
                OPstream toFirstSlave
                (
                    Pstream::blocking,
                    Pstream::firstSlave()
                );

                toFirstSlave << globalCutEdgesSet;
            }
            else
            {
                // Slave.  Read the list from the previous slave, do
                // local processing and send it to the next slave
                for
                (
                    int slave = Pstream::firstSlave();
                    slave <= Pstream::lastSlave();
                    slave++
                )
                {
                    if (Pstream::myProcNo() == slave)
                    {
                        // My turn.
                        int receiveFrom = slave - 1;
                        int sendTo = slave + 1;

                        if (Pstream::myProcNo() == Pstream::firstSlave())
                        {
                            // Receive from master
                            receiveFrom = Pstream::masterNo();
                        }

                        IPstream rf(Pstream::blocking, receiveFrom);

                        HashSet<edge, Hash<edge> > globalCutEdgesSet(rf);

                        // Check local cut edges against the list
                        forAll (localCutEdges, edgeI)
                        {
                            edge e
                            (
                                pointProcAddressing
                                    [localCutEdges[edgeI].start()],
                                pointProcAddressing[localCutEdges[edgeI].end()]
                            );

                            if (globalCutEdgesSet.found(e))
                            {
                                // Edge already exists.  Set mask to zero
                                localCutEdgeMask[edgeI] = 0;
                            }
                            else
                            {
                                // Edge not found. Add it to the list
                                globalCutEdgesSet.insert(e);
                            }
                        }

                        // Re-transmit the list to the next processor
                        if (slave < Pstream::lastSlave())
                        {
                            OPstream st(Pstream::blocking, sendTo);
                            st << globalCutEdgesSet;
                        }
                    }
                }
            }
        }

        // Add the processor point patch
        boundary_.setSize(boundary_.size() + 1);

        boundary_.set
        (
            boundary_.size() - 1,
            new globalTetPolyPatch
            (
                mesh_.globalData().nGlobalPoints(),
                localParPoints,
                mesh_.globalData().sharedPointAddr(),
                nGlobalParEdges,
                localParEdges,
                localEdgeMapping,
                localCutEdges,
                localCutEdgeMask,
                boundary_,
                boundary_.size() - 1
            )
        );

        // LDU addressing is required to create the boundary, but on addition
        // of a global patch it is incorrect.  Deleting it will
        // force recalculation.  HJ, 4/Nov/2007
        deleteDemandDrivenData(lduPtr_);
    }
}


// ************************************************************************* //

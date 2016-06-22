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

Class
    dynamicTopoFvMesh

Description
    Functions specific to coupled connectivity

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "dynamicTopoFvMesh.H"

#include "foamTime.H"
#include "eMesh.H"
#include "Random.H"
#include "triFace.H"
#include "volFields.H"
#include "changeMap.H"
#include "topoMapper.H"
#include "coupledInfo.H"
#include "matchPoints.H"
#include "SortableList.H"
#include "motionSolver.H"
#include "surfaceFields.H"
#include "globalMeshData.H"
#include "cyclicPolyPatch.H"
#include "fvMeshDistribute.H"
#include "faceSetAlgorithm.H"
#include "cellSetAlgorithm.H"
#include "coordinateSystem.H"
#include "processorPolyPatch.H"
#include "decompositionMethod.H"
#include "mapDistributePolyMesh.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupleMap, 0);

//! \cond fileScope
// Geometric relative match tolerance
static const Foam::debug::tolerancesSwitch geomMatchTol_
(
    "geomMatchTol",
    1e-4
);

// Priority scheme enumerants
enum priorityScheme
{
    LINEAR,
    RANDOM,
    LIST,
    ROTATE,
    MAX_PRIORITY_SCHEMES
};

// Priority scheme names
static const char* prioritySchemeNames_[MAX_PRIORITY_SCHEMES + 1] =
{
    "Linear",
    "Random",
    "List",
    "Rotate",
    "Invalid"
};

//! \endcond

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set coupled modification
void dynamicTopoFvMesh::setCoupledModification() const
{
    coupledModification_ = true;
}


// Unset coupled modification
void dynamicTopoFvMesh::unsetCoupledModification() const
{
    coupledModification_ = false;
}


// Initialize the coupled stack
void dynamicTopoFvMesh::initCoupledStack
(
    const labelHashSet& entities,
    bool useEntities
)
{
    // Clear existing lists/stacks.
    stack(0).clear();

    if (useEntities)
    {
        bool emptyEntity;

        // Initialize the stack with entries
        // in the labelHashSet and return
        forAllConstIter(labelHashSet, entities, eIter)
        {
            // Add only valid entities
            emptyEntity =
            (
                is2D() ?
                faces_[eIter.key()].empty() :
                edgeFaces_[eIter.key()].empty()
            );

            if (emptyEntity)
            {
                continue;
            }

            stack(0).insert(eIter.key());
        }

        if (debug > 3 && Pstream::parRun())
        {
            Pout<< nl << "Entity stack size: " << stack(0).size() << endl;

            if (debug > 4)
            {
                // Write out stack entities
                labelList stackElements(stack(0).size(), -1);

                forAll(stackElements, elemI)
                {
                    stackElements[elemI] = stack(0)[elemI];
                }

                label elemType = is2D() ? 2 : 1;

                writeVTK
                (
                    "entityStack_"
                  + Foam::name(Pstream::myProcNo()),
                    stackElements,
                    elemType
                );
            }
        }

        return;
    }

    // Loop though boundary faces and check whether
    // they belong to master/slave coupled patches.
    for (label faceI = nOldInternalFaces_; faceI < faces_.size(); faceI++)
    {
        // Add only valid faces
        if (faces_[faceI].empty())
        {
            continue;
        }

        label pIndex = whichPatch(faceI);

        if (pIndex == -1)
        {
            continue;
        }

        // Check if this is a locally coupled master face.
        if (patchCoupling_.size())
        {
            if (patchCoupling_(pIndex))
            {
                // Is patch cyclic?
                bool addIndex = true;

                if (isA<cyclicPolyPatch>(boundaryMesh()[pIndex]))
                {
                    // Check if this is a cyclic master face
                    const coupleMap& cMap = patchCoupling_[pIndex].map();
                    const Map<label>& fMap = cMap.entityMap(coupleMap::FACE);

                    if (!fMap.found(faceI))
                    {
                        addIndex = false;
                    }
                }

                // Add this to the coupled modification stack.
                if (addIndex)
                {
                    if (is2D())
                    {
                        stack(0).push(faceI);
                    }
                    else
                    {
                        const labelList& mfEdges = faceEdges_[faceI];

                        forAll(mfEdges, edgeI)
                        {
                            // Add this to the coupled modification stack.
                            stack(0).push(mfEdges[edgeI]);
                        }
                    }
                }
            }
        }

        // Check if this is a processor patch.
        label myProcID = Pstream::myProcNo();
        label neiProcID = getNeighbourProcessor(pIndex);

        // Check if this is a master processor patch.
        if
        (
            neiProcID > -1
         && priority(neiProcID, greaterOp<label>(), myProcID)
        )
        {
            // Add this to the coupled modification stack.
            if (is2D())
            {
                stack(0).push(faceI);
            }
            else
            {
                const label edgeEnum = coupleMap::EDGE;
                const label pointEnum = coupleMap::POINT;

                const labelList& mfEdges = faceEdges_[faceI];

                forAll(mfEdges, edgeI)
                {
                    label eIndex = mfEdges[edgeI];

                    const edge& checkEdge = edges_[eIndex];

                    // Need to avoid this edge if it is
                    // talking to higher-priority processors.
                    bool permissible = true;

                    forAll(procIndices_, pI)
                    {
                        label nProc = procIndices_[pI];

                        // List is specified in ascending order
                        // of neighbour processors, so break out.
                        if (priority(nProc, greaterOp<label>(), myProcID))
                        {
                            break;
                        }

                        // Fetch reference to subMeshes
                        const coupledMesh& recvMesh = recvMeshes_[pI];
                        const coupleMap& rcMap = recvMesh.map();

                        if (rcMap.findSlave(edgeEnum, eIndex) > -1)
                        {
                            permissible = false;
                            break;
                        }

                        // Check points as well, since there might be
                        // situations where both points lie on a lower
                        // ranked processor, but not the edge itself.
                        if
                        (
                            (rcMap.findSlave(pointEnum, checkEdge[0]) > -1)
                         && (rcMap.findSlave(pointEnum, checkEdge[1]) > -1)
                        )
                        {
                            permissible = false;
                            break;
                        }
                    }

                    if (permissible)
                    {
                        stack(0).push(eIndex);
                    }
                }
            }
        }
    }

    if (debug > 3 && Pstream::parRun())
    {
        Pout<< nl << "Coupled stack size: " << stack(0).size() << endl;

        if (debug > 4)
        {
            // Write out stack entities
            labelList stackElements(stack(0).size(), -1);

            forAll(stackElements, elemI)
            {
                stackElements[elemI] = stack(0)[elemI];
            }

            label elemType = is2D() ? 2 : 1;

            writeVTK
            (
                "coupledStack_"
              + Foam::name(Pstream::myProcNo()),
                stackElements,
                elemType
            );
        }
    }
}


// Execute load balancing operations on the mesh
void dynamicTopoFvMesh::executeLoadBalancing
(
    dictionary& balanceDict
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    bool enabled = readBool(balanceDict.lookup("enabled"));

    if (!enabled)
    {
        return;
    }

    // Read interval
    label interval = readLabel(balanceDict.lookup("interval"));

    // Are we at the interval?
    if ((interval < 0) || (time().timeIndex() % interval != 0))
    {
        return;
    }

    if (debug)
    {
        Info<< " void dynamicTopoFvMesh::executeLoadBalancing() :"
            << " Executing load-balancing operations on the mesh."
            << endl;

        // Write out patch boundaries before re-distribution
        for (label pI = 0; pI < nPatches_; pI++)
        {
            label neiProcNo = getNeighbourProcessor(pI);

            if (neiProcNo == -1)
            {
                continue;
            }

            writeVTK
            (
                "procPatchFacesNew_"
              + Foam::name(Pstream::myProcNo())
              + '_'
              + Foam::name(neiProcNo),
                identity(patchSizes_[pI]) + patchStarts_[pI],
                2
            );

            writeVTK
            (
                "procPatchFacesOld_"
              + Foam::name(Pstream::myProcNo())
              + '_'
              + Foam::name(neiProcNo),
                identity(patchSizes_[pI]) + patchStarts_[pI],
                2, false, true
            );
        }
    }

    // Ensure that the number of processors is consistent,
    // and silently modify dictionary entry if not
    balanceDict.set("numberOfSubdomains", Pstream::nProcs());

    // Read decomposition options and initialize
    autoPtr<decompositionMethod> decomposerPtr
    (
        decompositionMethod::New
        (
            balanceDict,
            (*this)
        )
    );

    // Alias for convenience...
    decompositionMethod& decomposer = decomposerPtr();

    if (!decomposer.parallelAware())
    {
        FatalErrorIn("void dynamicTopoFvMesh::executeLoadBalancing()")
            << "You have selected decomposition method "
            << decomposer.typeName
            << " which is not parallel aware."
            << exit(FatalError);
    }

    // Calculate a merge-distance
    const boundBox& box = polyMesh::bounds();
    scalar mergeTol = readScalar(balanceDict.lookup("mergeTol"));

    // Compute mergeDist
    scalar mergeDist = mergeTol * box.mag();

    // Compute write tolerance
    scalar writeTol =
    (
        std::pow
        (
            scalar(10.0),
           -scalar(IOstream::defaultPrecision())
        )
    );

    Info<< nl
        << " Overall mesh bounding box  : " << box << nl
        << " Relative tolerance         : " << mergeTol << nl
        << " Absolute matching distance : " << mergeDist << nl
        << endl;

    // Check for compatibility
    if (time().writeFormat() == IOstream::ASCII && mergeTol < writeTol)
    {
        FatalErrorIn("void dynamicTopoFvMesh::executeLoadBalancing()")
            << " Your current settings specify ASCII writing with "
            << IOstream::defaultPrecision() << " digits precision." << nl
            << " Your mergeTol (" << mergeTol << ") is finer than this." << nl
            << " Remedies: " << nl
            << "  - Change writeFormat to binary" << nl
            << "  - Increase writePrecision" << nl
            << "  - Adjust the merge tolerance (mergeTol)" << endl
            << exit(FatalError);
    }

    // Decompose the mesh and obtain distribution
    labelList cellDistribution
    (
        decomposer.decompose
        (
            primitiveMesh::cellCentres(),
            scalarField(nCells_, 1)
        )
    );

    // Set the coupledModification switch
    setCoupledModification();

    // Initialize the mesh distribution engine
    fvMeshDistribute distributor(*this, mergeDist);

    // Re-distribute mesh according to new decomposition
    distributor.distribute(cellDistribution);

    // Reset the coupledModification switch
    unsetCoupledModification();

    // Set old / new sizes
    nPoints_ = nOldPoints_ = primitiveMesh::nPoints();
    nEdges_ = nOldEdges_ = 0;
    nFaces_ = nOldFaces_ = primitiveMesh::nFaces();
    nCells_ = nOldCells_ = primitiveMesh::nCells();
    nInternalFaces_ = nOldInternalFaces_ = primitiveMesh::nInternalFaces();
    nInternalEdges_ = nOldInternalEdges_ = 0;

    // Now re-initialize all connectivity structures
    oldPoints_ = polyMesh::points();
    points_ = polyMesh::points();
    faces_ = polyMesh::faces();
    owner_ = polyMesh::faceOwner();
    neighbour_ = polyMesh::faceNeighbour();
    cells_ = primitiveMesh::cells();

    // Check the size of owner / neighbour
    if (owner_.size() != neighbour_.size())
    {
        // Size up to number of faces
        neighbour_.setSize(nFaces_, -1);
    }

    const polyBoundaryMesh& boundary = polyMesh::boundaryMesh();

    // Initialize patch-size information
    nPatches_ = boundary.size();

    oldPatchSizes_.setSize(nPatches_, 0);
    oldPatchStarts_.setSize(nPatches_, -1);
    oldEdgePatchSizes_.setSize(nPatches_, 0);
    oldEdgePatchStarts_.setSize(nPatches_, -1);
    oldPatchNMeshPoints_.setSize(nPatches_, -1);

    patchSizes_.setSize(nPatches_, 0);
    patchStarts_.setSize(nPatches_, -1);
    edgePatchSizes_.setSize(nPatches_, 0);
    edgePatchStarts_.setSize(nPatches_, -1);
    patchNMeshPoints_.setSize(nPatches_, -1);

    for (label i = 0; i < nPatches_; i++)
    {
        patchNMeshPoints_[i] = boundary[i].meshPoints().size();
        oldPatchSizes_[i] = patchSizes_[i] = boundary[i].size();
        oldPatchStarts_[i] = patchStarts_[i] = boundary[i].start();
        oldPatchNMeshPoints_[i] = patchNMeshPoints_[i];
    }

    // Clear pointers and re-initialize
    mapper_.clear();
    eMeshPtr_.clear();
    motionSolver_.clear();
    lengthEstimator_.clear();

    // Initialize edge-related connectivity structures
    initEdges();

    // Load the mesh-motion solver
    loadMotionSolver();

    // Load the field-mapper
    loadFieldMapper();

    // Load the length-scale estimator,
    // and read refinement options
    loadLengthScaleEstimator();

    // Clear parallel structures
    procIndices_.clear();
    procPriority_.clear();
    sendMeshes_.clear();
    recvMeshes_.clear();
}


// Set processor rank priority
void dynamicTopoFvMesh::initProcessorPriority()
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (Pstream::master())
    {
        // Default to linear
        int type = LINEAR;

        const dictionary& meshSubDict = dict_.subDict("dynamicTopoFvMesh");

        if (meshSubDict.found("priorityScheme") || mandatory_)
        {
            word schemeType(meshSubDict.lookup("priorityScheme"));

            type = MAX_PRIORITY_SCHEMES;

            for (label i = 0; i < MAX_PRIORITY_SCHEMES; i++)
            {
                if (prioritySchemeNames_[i] == schemeType)
                {
                    type = i;
                    break;
                }
            }
        }

        switch (type)
        {
            case LINEAR:
            {
                // Initialize to identity map
                procPriority_ = identity(Pstream::nProcs());
                break;
            }

            case RANDOM:
            {
                Random randomizer(std::time(NULL));

                // Initialize to identity map
                procPriority_ = identity(Pstream::nProcs());

                for (label i = 0; i < Pstream::nProcs(); i++)
                {
                    // Specify a random index to shuffle
                    label j = randomizer.integer(0, Pstream::nProcs() - 1);

                    Foam::Swap(procPriority_[i], procPriority_[j]);
                }

                break;
            }

            case ROTATE:
            {
                const label listSize = Pstream::nProcs();

                // Initialize to identity map
                procPriority_ = identity(listSize);

                // In-place reverse macro
#               define inPlaceReverse(list)                                    \
                {                                                              \
                    const label listSize = list.size();                        \
                    const label lastIndex = listSize - 1;                      \
                    const label nIterations = listSize >> 1;                   \
                                                                               \
                    label elemI = 0;                                           \
                    while (elemI < nIterations)                                \
                    {                                                          \
                        Swap(list[elemI], list[lastIndex - elemI]);            \
                        elemI++;                                               \
                    }                                                          \
                }

                // Rotate by time index
                label n = time().timeIndex() % listSize;

                if (n < 0)
                {
                    n += listSize;
                }

                // Rotate list in place
                SubList<label> firstHalf(procPriority_, n, 0);
                SubList<label> secondHalf(procPriority_, listSize - n, n);

                inPlaceReverse(firstHalf);
                inPlaceReverse(secondHalf);

                inPlaceReverse(procPriority_);

                break;
            }

            case LIST:
            {
                procPriority_ = labelList(meshSubDict.lookup("priorityList"));

                if (procPriority_.size() != Pstream::nProcs())
                {
                    FatalErrorIn
                    (
                        "void dynamicTopoFvMesh::initProcessorPriority()"
                    )
                        << " Priority scheme size mismatch." << nl
                        << "   List size: " << procPriority_.size() << nl
                        << "   Processor count: " << Pstream::nProcs() << nl
                        << abort(FatalError);
                }

                break;
            }

            default:
            {
                Info<< "Available schemes: " << endl;

                for (label i = 0; i < MAX_PRIORITY_SCHEMES; i++)
                {
                    Info<< ' ' << prioritySchemeNames_[i] << endl;
                }

                FatalErrorIn
                (
                    "void dynamicTopoFvMesh::initProcessorPriority()"
                )
                    << " Unknown priority scheme." << nl
                    << abort(FatalError);
            }
        }

        // Send priority list to others
        for (label proc = 1; proc < Pstream::nProcs(); proc++)
        {
            meshOps::pWrite(proc, procPriority_);
        }
    }
    else
    {
        // Receive priority from master
        procPriority_.setSize(Pstream::nProcs());
        meshOps::pRead(Pstream::masterNo(), procPriority_);
    }

    // Wait for transfers to complete
    meshOps::waitForBuffers();

    if (debug)
    {
        labelList ranks(Pstream::nProcs(), -1);

        forAll(ranks, procI)
        {
            ranks[procPriority_[procI]] = procI;
        }

        Info<< " Processor priority: " << procPriority_ << nl
            << " Processor rankings: " << ranks
            << endl;
    }
}


// Identify coupled patches.
//  - Also builds global shared point information.
//  - Returns true if no coupled patches were found.
bool dynamicTopoFvMesh::identifyCoupledPatches()
{
    bool coupledPatchesAbsent = true;

    // Check if patches are explicitly coupled
    if (patchCoupling_.size())
    {
        coupledPatchesAbsent = false;
    }

    // Maintain a separate list of processor IDs in procIndices.
    // This is done because this sub-domain may talk to processors
    // that share only edges/points.
    if (!Pstream::parRun() || procIndices_.size())
    {
        return coupledPatchesAbsent;
    }

    // Prepare a list of points for sub-mesh creation.
    //  - Obtain global shared-points information, if necessary.
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Set flag as default for parallel runs
    coupledPatchesAbsent = false;

    // Fetch the list of global points from polyMesh.
    //  - This should be the first evaluation of globalData,
    //    so this involves global communication
    const globalMeshData& gData = polyMesh::globalData();

    const labelList& spA = gData.sharedPointAddr();
    const labelList& spL = gData.sharedPointLabels();

    labelListList spB(Pstream::nProcs(), labelList(0));

    if (gData.nGlobalPoints())
    {
        if (debug)
        {
            Info<< " dynamicTopoFvMesh::identifyCoupledPatches :"
                << " Found " << gData.nGlobalPoints()
                << " global points." << endl;
        }

        // Send others my addressing.
        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            if (proc != Pstream::myProcNo())
            {
                // Send number of entities first.
                meshOps::pWrite(proc, spA.size());

                // Send the buffer.
                if (spA.size())
                {
                    meshOps::pWrite(proc, spA);
                }
            }
        }

        // Receive addressing from others
        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            if (proc != Pstream::myProcNo())
            {
                label procInfoSize = -1;

                // How many entities am I going to be receiving?
                meshOps::pRead(proc, procInfoSize);

                if (procInfoSize)
                {
                    // Size the receive buffer.
                    spB[proc].setSize(procInfoSize, -1);

                    // Schedule for receipt.
                    meshOps::pRead(proc, spB[proc]);
                }
            }
        }
    }
    else
    if (debug)
    {
        Info<< " dynamicTopoFvMesh::identifyCoupledPatches :"
            << " Did not find any global points." << endl;
    }

    Map<label> immNeighbours;
    labelListList procPoints(Pstream::nProcs());

    // Track globally shared points
    List<DynamicList<labelPair> > globalProcPoints
    (
        Pstream::nProcs(),
        DynamicList<labelPair>(5)
    );

    // Insert my immediate neighbours into the list.
    forAll(boundary, pI)
    {
        if (isA<processorPolyPatch>(boundary[pI]))
        {
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pI])
            );

            label neiProcNo = pp.neighbProcNo();

            // Insert all boundary points.
            procPoints[neiProcNo] = pp.meshPoints();

            // Keep track of immediate neighbours.
            immNeighbours.insert(neiProcNo, pI);
        }
    }

    if (gData.nGlobalPoints())
    {
        // Wait for all transfers to complete.
        meshOps::waitForBuffers();

        // Now loop through all processor addressing, and check if
        // any labels coincide with my global shared points.
        // If this is true, we need to be talking to that neighbour
        // as well (if not already).
        for (label proc = 0; proc < Pstream::nProcs(); proc++)
        {
            if (proc == Pstream::myProcNo())
            {
                continue;
            }

            bool foundGlobalMatch = false;

            // Fetch reference to buffer
            const labelList& procBuffer = spB[proc];

            forAll(procBuffer, pointI)
            {
                forAll(spA, pointJ)
                {
                    if (spA[pointJ] == procBuffer[pointI])
                    {
                        // Make an entry, if one wasn't made already
                        if (findIndex(procPoints[proc], spL[pointJ]) == -1)
                        {
                            globalProcPoints[proc].append
                            (
                                labelPair(spL[pointJ], spA[pointJ])
                            );
                        }

                        foundGlobalMatch = true;

                        break;
                    }
                }
            }

            if (!immNeighbours.found(proc))
            {
                if (foundGlobalMatch && debug)
                {
                    Pout<< " dynamicTopoFvMesh::identifyCoupledPatches :"
                        << " Additionally talking to processor: "
                        << proc << endl;
                }
            }
        }
    }

    // Estimate an initial size
    procIndices_.setSize(Pstream::nProcs());

    label nTotalProcs = 0;

    forAll(procPoints, procI)
    {
        if (procPoints[procI].size())
        {
            procIndices_[nTotalProcs++] = procI;
        }

        // Check for point / edge coupling
        if (globalProcPoints[procI].size() && !procPoints[procI].size())
        {
            procIndices_[nTotalProcs++] = procI;
        }
    }

    // Shrink to actual size
    procIndices_.setSize(nTotalProcs);

    // Initialize the processor priority list
    initProcessorPriority();

    // Sort by order of neighbouring
    // processor priority - highest to lowest
    Foam::sort
    (
        procIndices_,
        UList<label>::less(procPriority_)
    );

    // Size the PtrLists.
    sendMeshes_.setSize(nTotalProcs);
    recvMeshes_.setSize(nTotalProcs);

    // Create send/recv patch meshes, and copy
    // the list of points for each processor.
    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI], master = -1, slave = -1;

        // For processors that have only point / edge coupling
        // specify an invalid patch ID for now
        label patchID = immNeighbours.found(proc) ? immNeighbours[proc] : -1;

        // Check if this processor is of higher priority
        if (priority(proc, lessOp<label>(), Pstream::myProcNo()))
        {
            master = proc;
            slave = Pstream::myProcNo();
        }
        else
        {
            master = Pstream::myProcNo();
            slave = proc;
        }

        sendMeshes_.set
        (
            pI,
            new coupledMesh
            (
                *this,               // Reference to this mesh
                is2D(),              // 2D or 3D
                false,               // Not local
                true,                // Sent to neighbour
                patchID,             // Patch index
                master,              // Master index
                slave                // Slave index
            )
        );

        sendMeshes_[pI].map().subMeshPoints() =
        (
            procPoints[procIndices_[pI]]
        );

        sendMeshes_[pI].map().globalProcPoints() =
        (
            globalProcPoints[procIndices_[pI]]
        );

        recvMeshes_.set
        (
            pI,
            new coupledMesh
            (
                *this,               // Reference to this mesh
                is2D(),              // 2D or 3D
                false,               // Not local
                false,               // Not sent to neighbour
                patchID,             // Patch index
                master,              // Master index
                slave                // Slave index
            )
        );

        recvMeshes_[pI].map().subMeshPoints() =
        (
            procPoints[procIndices_[pI]]
        );

        recvMeshes_[pI].map().globalProcPoints() =
        (
            globalProcPoints[procIndices_[pI]]
        );
    }

    if (debug > 3)
    {
        Pout<< " identifyCoupledPatches :"
            << " Talking to processors: "
            << procIndices_ << endl;

        forAll(procIndices_, pI)
        {
            label proc = procIndices_[pI];

            // Write out subMeshPoints as a VTK
            const coupleMap& cMap = sendMeshes_[pI].map();

            writeVTK
            (
                "subMeshPoints_" +
                Foam::name(Pstream::myProcNo()) + "to" +
                Foam::name(proc),
                cMap.subMeshPoints(),
                0, false, true
            );

            // Write out globalProcPoints as a VTK
            const List<labelPair>& gpp = cMap.globalProcPoints();

            if (gpp.size())
            {
                labelList gpPoints(gpp.size(), -1);

                forAll(gpp, pointI)
                {
                    gpPoints[pointI] = gpp[pointI].first();
                }

                writeVTK
                (
                    "globalProcPoints_" +
                    Foam::name(Pstream::myProcNo()) + "to" +
                    Foam::name(proc),
                    gpPoints,
                    0, false, true
                );
            }
        }
    }

    return coupledPatchesAbsent;
}


// Read coupled patch information from dictionary.
void dynamicTopoFvMesh::readCoupledPatches()
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Clear list
    patchCoupling_.clear();

    if (dict_.found("coupledPatches") || mandatory_)
    {
        // Size the initial list
        patchCoupling_.setSize(boundary.size());

        const dictionary& coupledPatches = dict_.subDict("coupledPatches");

        // Determine master and slave patches
        forAllConstIter(dictionary, coupledPatches, dIter)
        {
            const dictionary& dictI = dIter().dict();

            // Lookup the master / slave patches
            word masterPatch = dictI.lookup("master");
            word slavePatch = dictI.lookup("slave");

            // Determine patch indices
            label mPatch = boundary.findPatchID(masterPatch);
            label sPatch = boundary.findPatchID(slavePatch);

            if (debug)
            {
                Info<< " Adding master: " << masterPatch << " : " << mPatch
                    << " with slave: " << slavePatch << " : " << sPatch
                    << endl;
            }

            // Add to the list if entries are legitimate
            if
            (
                mPatch != sPatch &&
                boundary[mPatch].size() == boundary[sPatch].size()
            )
            {
                // Check whether patches are associated with zones.
                Switch specifyZones
                (
                    dictI.lookup("specifyZones")
                );

                label mZone = -1, sZone = -1;

                if (specifyZones)
                {
                    const faceZoneMesh& faceZones = polyMesh::faceZones();

                    mZone = faceZones.findZoneID
                    (
                        dictI.lookup("masterZone")
                    );

                    sZone = faceZones.findZoneID
                    (
                        dictI.lookup("slaveZone")
                    );
                }

                // Configure with regIOobject for check-in
                coupleMap cMap
                (
                    IOobject
                    (
                        "coupleMap_"
                      + Foam::name(mPatch)
                      + "_To_"
                      + Foam::name(sPatch)
                      + "_Local",
                        time().timeName(),
                        *this,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE,
                        true
                    ),
                    is2D(),
                    true,
                    false,
                    mPatch,
                    mPatch,
                    sPatch
                );

                // Set the pointer
                patchCoupling_.set
                (
                    mPatch,
                    new coupledMesh
                    (
                        *this,
                        cMap,
                        mZone,
                        sZone
                    )
                );
            }
            else
            {
                FatalErrorIn("void dynamicTopoFvMesh::readCoupledPatches()")
                    << " Coupled patches are either wrongly specified,"
                    << " or the sizes don't match." << nl
                    << " Master: " << mPatch << ":" << masterPatch
                    << " Size: " << boundary[mPatch].size() << nl
                    << " Slave: " << sPatch << ":" << slavePatch
                    << " Size: " << boundary[sPatch].size() << nl
                    << abort(FatalError);
            }
        }
    }

    // Loop through boundaries and add any cyclic patches
    forAll(boundary, patchI)
    {
        if (!isA<cyclicPolyPatch>(boundary[patchI]))
        {
            continue;
        }

        // Size up patchCoupling, if necessary
        if (patchCoupling_.empty())
        {
            patchCoupling_.setSize(boundary.size());
        }

        if (debug)
        {
            Info<< " Adding cyclic: " << patchI
                << " : " << boundary[patchI].name()
                << endl;
        }

        // Configure with regIOobject for check-in
        coupleMap cMap
        (
            IOobject
            (
                "coupleMap_"
              + Foam::name(patchI)
              + "_To_"
              + Foam::name(patchI)
              + "_Local",
                time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                true
            ),
            is2D(),
            true,
            false,
            patchI,
            patchI,
            patchI
        );

        // Set the pointer
        patchCoupling_.set
        (
            patchI,
            new coupledMesh(*this, cMap, -1, -1)
        );
    }
}


// Initialize coupled patch connectivity for topology modifications.
//  - Send and receive sub meshes for processor patches.
//  - Made static because this function may be run in a separate thread.
void dynamicTopoFvMesh::initCoupledConnectivity
(
    void *argument
)
{
    // Recast the argument
    dynamicTopoFvMesh *mesh = static_cast<dynamicTopoFvMesh*>(argument);

    // Identify coupled patches.
    if (mesh->identifyCoupledPatches())
    {
        return;
    }

    // Build and send patch sub-meshes.
    mesh->buildProcessorPatchMeshes();

    // Build inital maps for locally coupled patches.
    mesh->buildLocalCoupledMaps();

    // Build maps for coupled processor patches.
    mesh->buildProcessorCoupledMaps();
}


// Move coupled subMesh points
void dynamicTopoFvMesh::moveCoupledSubMeshes()
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (debug)
    {
        Info<< " void dynamicTopoFvMesh::moveCoupledSubMeshes() :"
            << " Moving points for coupled subMeshes."
            << endl;
    }

    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        const coupledMesh& sPM = sendMeshes_[pI];
        const coupledMesh& rPM = recvMeshes_[pI];

        // Fetch the coupleMap
        const coupleMap& scMap = sPM.map();
        const coupleMap& rcMap = rPM.map();

        Map<label>& pointMap = scMap.entityMap(coupleMap::POINT);

        // Fill point buffers
        pointField& pBuffer = scMap.pointBuffer();
        pointField& opBuffer = scMap.oldPointBuffer();

        forAllConstIter(Map<label>, pointMap, pIter)
        {
            pBuffer[pIter.key()] = points_[pIter()];
            opBuffer[pIter.key()] = oldPoints_[pIter()];
        }

        // Buffers have already been allocated
        // to the right size, so just transfer points

        // Send point buffers to neighbour
        meshOps::pWrite(proc, scMap.pointBuffer());
        meshOps::pWrite(proc, scMap.oldPointBuffer());

        // Receive point buffers from neighbour
        meshOps::pRead(proc, rcMap.pointBuffer());
        meshOps::pRead(proc, rcMap.oldPointBuffer());
    }

    // Wait for transfers to complete before moving on
    meshOps::waitForBuffers();

    // Set points in mesh
    forAll(procIndices_, pI)
    {
        // Fetch non-const reference to patchSubMesh
        coupledMesh& rPM = recvMeshes_[pI];

        // Fetch the coupleMap
        const coupleMap& rcMap = rPM.map();

        dynamicTopoFvMesh& mesh = rPM.subMesh();
        const polyBoundaryMesh& boundary = mesh.boundaryMesh();

        // Set points / oldPoints
        mesh.points_ = rcMap.pointBuffer();
        mesh.oldPoints_ = rcMap.oldPointBuffer();

        labelList patchSizes(boundary.size(), -1);
        labelList patchStarts(boundary.size(), -1);

        forAll(boundary, patchI)
        {
            patchSizes[patchI] = boundary[patchI].size();
            patchStarts[patchI] = boundary[patchI].start();
        }

        // Reset underlying mesh.
        //  - Use null lists for addressing to avoid over-writes
        //  - Specify non-valid boundary to avoid globalData creation
        mesh.resetPrimitives
        (
            xferCopy(rcMap.oldPointBuffer()),
            (*reinterpret_cast<Xfer<faceList>*>(0)),
            (*reinterpret_cast<Xfer<labelList>*>(0)),
            (*reinterpret_cast<Xfer<labelList>*>(0)),
            patchSizes,
            patchStarts,
            false
        );
    }
}


// Insert the cells around the coupled master entity to the mesh
// - Returns a changeMap with a type specifying:
//     1: Insertion was successful
//    -1: Insertion failed
//    -2: Failed because entity was being handled elsewhere
// - The changeMap index specifies the converted mIndex.
const changeMap dynamicTopoFvMesh::insertCells(const label mIndex)
{
    // Prepare the changeMaps
    changeMap map;

    // Maintain face counts for each inserted cell
    Map<label> nCellFaces;

    // Track inserted entities
    List<Map<label> > pointsToInsert(procIndices_.size());
    List<Map<label> > edgesToInsert(procIndices_.size());
    List<Map<label> > facesToInsert(procIndices_.size());
    List<Map<label> > cellsToInsert(procIndices_.size());

    // Track converted entities
    List<Map<label> > edgesToConvert(procIndices_.size());
    List<Map<label> > facesToConvert(procIndices_.size());

    // Track edge / face patch information
    List<Map<label> > masterConvertEdgePatch(procIndices_.size());
    List<Map<label> > masterConvertFacePatch(procIndices_.size());

    Map<label> createPatch;
    labelList slaveConvertPatch(procIndices_.size(), -1);
    List<Map<Pair<point> > > convertPatchPoints(procIndices_.size());

    // First check to ensure that this case can be handled
    forAll(procIndices_, pI)
    {
        const label cplEnum = is2D() ? coupleMap::FACE : coupleMap::EDGE;

        // Fetch reference to maps
        const coupleMap& cMap = recvMeshes_[pI].map();

        // Does a coupling exist?
        if (cMap.findSlave(cplEnum, mIndex) == -1)
        {
            continue;
        }

        // If this entity is being handled elsewhere, bail out
        if (priority(procIndices_[pI], lessOp<label>(), Pstream::myProcNo()))
        {
            map.type() = -2;

            return map;
        }
    }

    const polyBoundaryMesh& boundary = boundaryMesh();

    // Agglomerate cells from surrounding subMeshes
    forAll(procIndices_, pI)
    {
        const label cplEnum = is2D() ? coupleMap::FACE : coupleMap::EDGE;

        // Fetch reference to maps
        const coupleMap& cMap = recvMeshes_[pI].map();
        const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

        label sIndex = -1;

        if ((sIndex = cMap.findSlave(cplEnum, mIndex)) == -1)
        {
            continue;
        }

        // Fetch references
        Map<label>& procCellMap = cellsToInsert[pI];

        if (is2D())
        {
            // Insert the owner cell
            procCellMap.insert(mesh.owner_[sIndex], -1);
        }
        else
        {
            // Insert all cells connected to this edge
            const labelList& eFaces = mesh.edgeFaces_[sIndex];

            forAll(eFaces, faceI)
            {
                const label own = mesh.owner_[eFaces[faceI]];
                const label nei = mesh.neighbour_[eFaces[faceI]];

                if (!procCellMap.found(own))
                {
                    procCellMap.insert(own, -1);
                }

                if (!procCellMap.found(nei) && (nei != -1))
                {
                    procCellMap.insert(nei, -1);
                }
            }
        }

        // Find a patch that talks to this processor
        label neiProcPatch = -1;

        forAll(boundary, patchI)
        {
            if (getNeighbourProcessor(patchI) == procIndices_[pI])
            {
                neiProcPatch = patchI;
                break;
            }
        }

        // Prepare information about inserted / converted faces and edges
        const polyBoundaryMesh& slaveBoundary = mesh.boundaryMesh();

        // Set the slave conversion patch
        forAll(slaveBoundary, patchI)
        {
            if (mesh.getNeighbourProcessor(patchI) == Pstream::myProcNo())
            {
                slaveConvertPatch[pI] = patchI;
                break;
            }
        }

        forAllIter(Map<label>, procCellMap, cIter)
        {
            const cell& checkCell = mesh.cells_[cIter.key()];

            forAll(checkCell, fI)
            {
                label mfIndex = -1, sfIndex = checkCell[fI];
                label sfPatch = mesh.whichPatch(sfIndex);

                // Check whether a face mapping exists for this face
                if ((mfIndex = cMap.findMaster(coupleMap::FACE, sfIndex)) > -1)
                {
                    // Mark as a candidate for conversion
                    facesToInsert[pI].insert(sfIndex, mfIndex);
                    facesToConvert[pI].insert(sfIndex, mfIndex);
                }
                else
                {
                    // Mark for insertion. Some faces may be duplicated
                    // across processors, but this will be dealt
                    // with at a later stage.
                    facesToInsert[pI].insert(sfIndex, -1);
                }

                // Check face points for coupling
                const face& sFace = mesh.faces_[sfIndex];

                // Configure points which need to be inserted
                forAll(sFace, pointI)
                {
                    label sFacePoint = sFace[pointI];

                    // Skip if added already
                    if (pointsToInsert[pI].found(sFacePoint))
                    {
                        continue;
                    }

                    // Check for coupled points
                    pointsToInsert[pI].insert
                    (
                        sFacePoint,
                        cMap.findMaster
                        (
                            coupleMap::POINT,
                            sFacePoint
                        )
                    );
                }

                // Loop through edges and check whether edge-mapping exists
                const labelList& sfEdges = mesh.faceEdges_[sfIndex];

                forAll(sfEdges, edgeI)
                {
                    const label seIndex = sfEdges[edgeI];
                    const edge& sEdge = mesh.edges_[seIndex];

                    // Meshes in 2D don't have edge-mapping, so check
                    // point maps instead. If either point doesn't exist,
                    // this is an edge that needs to be inserted.
                    label cMs = cMap.findMaster(coupleMap::POINT, sEdge[0]);
                    label cMe = cMap.findMaster(coupleMap::POINT, sEdge[1]);

                    if (!edgesToInsert[pI].found(seIndex))
                    {
                        edgesToInsert[pI].insert(seIndex, -1);
                    }

                    // If both points have maps, this is a conversion edge
                    if (cMs > -1 && cMe > -1)
                    {
                        if (!edgesToConvert[pI].found(seIndex))
                        {
                            edgesToConvert[pI].insert(seIndex, -1);
                        }
                    }

                    // Add a patch entry as well
                    if (masterConvertEdgePatch[pI].found(seIndex))
                    {
                        continue;
                    }

                    label edgePatch = -1;
                    label sePatch = mesh.whichEdgePatch(seIndex);
                    label neiProcNo = mesh.getNeighbourProcessor(sePatch);

                    // Determine patch
                    if (sePatch == -1)
                    {
                        // Slave edge was an interior one
                        if (neiProcPatch > -1)
                        {
                            edgePatch = neiProcPatch;
                        }
                        else
                        {
                            neiProcNo = procIndices_[pI];

                            if (debug > 3)
                            {
                                Pout<< nl << nl
                                    << " No face contact with"
                                    << " processor: " << neiProcNo
                                    << endl;
                            }

                            label pIdx = sendMeshes_[pI].map().patchIndex();

                            // Add a patch creation order, if necessary
                            if (pIdx == -1 && !createPatch.found(neiProcNo))
                            {
                                createPatch.insert(neiProcNo, -1);
                            }

                            // Specify a value for edgePatch:
                            //  - Specify a value that we can use
                            //    to back out the patch after creation
                            //  - Also needs to bypass failure check
                            edgePatch = (-2 - neiProcNo);
                        }
                    }
                    else
                    if (sePatch == (slaveBoundary.size() - 1))
                    {
                        // The 'defaultPatch'
                        edgePatch = neiProcPatch;
                    }
                    else
                    if (neiProcNo > -1)
                    {
                        if (neiProcNo == Pstream::myProcNo())
                        {
                            // Set the master edge patch
                            edgePatch = neiProcPatch;
                        }
                        else
                        {
                            // If this other processor is
                            // lesser-ranked, bail out.
                            if
                            (
                                priority
                                (
                                    neiProcNo,
                                    lessOp<label>(),
                                    Pstream::myProcNo()
                                )
                            )
                            {
                                map.type() = -2;

                                return map;
                            }
                            else
                            if (debug > 3)
                            {
                                Pout<< nl << nl
                                    << " Edge: " << seIndex
                                    << " :: " << mesh.edges_[seIndex]
                                    << " is talking to processor: "
                                    << neiProcNo << endl;
                            }

                            // Find an appropriate boundary patch
                            forAll(boundary, patchI)
                            {
                                label eP = getNeighbourProcessor(patchI);

                                if (eP == neiProcNo)
                                {
                                    edgePatch = patchI;
                                    break;
                                }
                            }

                            if (edgePatch == -1)
                            {
                                if (debug > 3)
                                {
                                    Pout<< nl << nl
                                        << " No direct contact with"
                                        << " processor: " << neiProcNo
                                        << " so adding to: " << neiProcPatch
                                        << endl;
                                }

                                // Add this to neiProcPatch instead.
                                edgePatch = neiProcPatch;
                            }
                        }
                    }
                    else
                    {
                        // Physical type
                        edgePatch = sePatch;
                    }

                    if (edgePatch == -1)
                    {
                        // Write out the edge
                        mesh.writeVTK("seEdge", seIndex, 1);

                        Pout<< " Could not find correct patch info: " << nl
                            << " sEdge: " << sEdge << nl
                            << " seIndex: " << seIndex << nl
                            << " sePatch: " << sePatch << nl
                            << " neiProcPatch: " << neiProcPatch << nl
                            << " neiProcNo: " << neiProcNo << nl
                            << " proc: " << procIndices_[pI] << nl
                            << " cMs: " << cMs << " cMe: " << cMe << nl
                            << " Patch Name: " <<
                            (
                                sePatch > -1 ?
                                slaveBoundary[sePatch].name() :
                                "Internal"
                            )
                            << abort(FatalError);
                    }

                    // Add the patch entry
                    masterConvertEdgePatch[pI].insert
                    (
                        seIndex,
                        edgePatch
                    );
                }

                // Determine patch
                if (masterConvertFacePatch[pI].found(sfIndex))
                {
                    continue;
                }

                label facePatch = -1;
                label neiProcNo = mesh.getNeighbourProcessor(sfPatch);

                if (sfPatch == -1)
                {
                    // Slave face was an interior one
                    if (neiProcPatch > -1)
                    {
                        facePatch = neiProcPatch;
                    }
                    else
                    {
                        neiProcNo = procIndices_[pI];

                        if (debug > 3)
                        {
                            Pout<< nl << nl
                                << " No face contact with"
                                << " processor: " << neiProcNo
                                << endl;
                        }

                        label pIdx = sendMeshes_[pI].map().patchIndex();

                        // Add a patch creation order, if necessary
                        if (pIdx == -1 && !createPatch.found(neiProcNo))
                        {
                            createPatch.insert(neiProcNo, -1);
                        }

                        // Specify a value for facePatch:
                        //  - Specify a value that we can use
                        //    to back out the patch after creation
                        //  - Also needs to bypass failure check
                        facePatch = (-2 - neiProcNo);
                    }
                }
                else
                if (sfPatch == (slaveBoundary.size() - 1))
                {
                    // The 'defaultPatch'
                    facePatch = neiProcPatch;
                }
                else
                if (neiProcNo > -1)
                {
                    if (neiProcNo == Pstream::myProcNo())
                    {
                        // Check to see if this face was converted before
                        if
                        (
                            findIndex
                            (
                                cMap.entityOperations(),
                                coupleMap::CONVERT_PATCH
                            ) > -1
                        )
                        {
                            // Loop through points and check for match
                            const pointField& pts = cMap.moveNewPoints();
                            const face& fCheck = mesh.faces_[sfIndex];
                            const point fC = fCheck.centre(mesh.points_);

                            scalar tol = mag(fC - mesh.points_[fCheck[0]]);

                            forAll(pts, ptI)
                            {
                                scalar dist = mag(fC - pts[ptI]);

                                if (dist < (geomMatchTol_() * tol))
                                {
                                    // Face was converted before
                                    if (debug > 3)
                                    {
                                        Pout<< nl << nl
                                            << " Face: " << sfIndex
                                            << " :: " << mesh.faces_[sfIndex]
                                            << " was converted before."
                                            << endl;
                                    }

                                    map.type() = -2;

                                    return map;
                                }
                            }
                        }

                        // Set the master face patch
                        facePatch = neiProcPatch;
                    }
                    else
                    {
                        // If this other processor is
                        // lesser-ranked, bail out.
                        if
                        (
                            priority
                            (
                                neiProcNo,
                                lessOp<label>(),
                                Pstream::myProcNo()
                            )
                        )
                        {
                            map.type() = -2;

                            return map;
                        }
                        else
                        if (debug > 3)
                        {
                            Pout<< nl << nl
                                << " Face: " << sfIndex
                                << " :: " << mesh.faces_[sfIndex]
                                << " is talking to processor: "
                                << neiProcNo << endl;
                        }

                        // Find an appropriate boundary patch
                        forAll(boundary, patchI)
                        {
                            label fP = getNeighbourProcessor(patchI);

                            if (fP == neiProcNo)
                            {
                                facePatch = patchI;
                                break;
                            }
                        }

                        // Specify a convert-patch operation
                        // for later, if this operation succeeds
                        const face& sfCheck = mesh.faces_[sfIndex];

                        point nfC = sfCheck.centre(mesh.points_);
                        point ofC = sfCheck.centre(mesh.oldPoints_);

                        // Are we talking to this processor already?
                        label prI = findIndex(procIndices_, neiProcNo);

                        // Check for face-contact
                        if (facePatch == -1)
                        {
                            if (debug > 3)
                            {
                                Pout<< nl << nl
                                    << " No face contact with"
                                    << " processor: " << neiProcNo
                                    << endl;
                            }

                            if (prI == -1)
                            {
                                Pout<< " No contact with: " << neiProcNo
                                    << abort(FatalError);
                            }

                            label pIdx = sendMeshes_[prI].map().patchIndex();

                            // Add a patch creation order, if necessary
                            if (pIdx == -1 && !createPatch.found(neiProcNo))
                            {
                                createPatch.insert(neiProcNo, -1);
                            }

                            // Specify a value for facePatch:
                            //  - Specify a value that we can use
                            //    to back out the patch after creation
                            //  - Also needs to bypass failure check
                            facePatch = (-2 - neiProcNo);
                        }

                        // Add to the list
                        convertPatchPoints[prI].insert
                        (
                            sfIndex,
                            Pair<point>(nfC, ofC)
                        );
                    }
                }
                else
                {
                    // Physical type
                    facePatch = sfPatch;
                }

                if (facePatch == -1)
                {
                    // Write out the face
                    mesh.writeVTK("sfFace", sfIndex, 2);

                    Pout<< " Could not find correct patch info: " << nl
                        << " sfIndex: " << sfIndex << nl
                        << " sfPatch: " << sfPatch << nl
                        << " neiProcPatch: " << neiProcPatch << nl
                        << " slavePatch: " <<
                        (
                            sfPatch > -1 ?
                            slaveBoundary[sfPatch].name() :
                            "Internal"
                        )
                        << abort(FatalError);
                }

                // Add the patch entry
                masterConvertFacePatch[pI].insert
                (
                    sfIndex,
                    facePatch
                );
            }
        }
    }

    // Specify a merge tolerance for insertion points
    scalar mergeTol =
    (
        is2D() ? 0.0 : geomMatchTol_() * magSqr(edges_[mIndex].vec(points_))
    );

    // Check for point / edge processor connections
    if (is3D())
    {
        const label myProcNo = Pstream::myProcNo();

        forAll(procIndices_, pI)
        {
            const Map<label>& cEdges = edgesToInsert[pI];
            const coupleMap& cMap = recvMeshes_[pI].map();
            const Map<labelList>& pEdgeMap = cMap.subMeshEdgeMap();
            const Map<labelList>& pPointMap = cMap.subMeshPointMap();
            const dynamicTopoFvMesh& sMesh = recvMeshes_[pI].subMesh();

            // Some points may have been inserted by a prior
            // operation involving disconnected entities.
            // Check if points need to be truly inserted.
            Map<label>& cPoints = pointsToInsert[pI];

            forAllIter(Map<label>, cPoints, pIter)
            {
                if (pIter() > -1)
                {
                    continue;
                }

                label cSp = pIter.key();

                Map<labelList>::const_iterator pIt = pPointMap.find(cSp);

                if (pIt == pPointMap.end())
                {
                    continue;
                }

                const point& pointI = sMesh.points_[pIt.key()];

                bool foundMaster = false;
                const labelList& procs = pIt();

                forAll(procs, procJ)
                {
                    label pJ = -1, nProcNo = procs[procJ];

                    if
                    (
                        (nProcNo == procIndices_[pI]) ||
                        ((pJ = findIndex(procIndices_, nProcNo)) == -1)
                    )
                    {
                        continue;
                    }

                    // Fetch map from the other processor,
                    // and check if a mapping already exists
                    const coupleMap& cMapJ = recvMeshes_[pJ].map();
                    const Map<labelList>& pMapJ = cMapJ.subMeshPointMap();
                    const dynamicTopoFvMesh& sMeshJ = recvMeshes_[pJ].subMesh();

                    forAllConstIter(Map<labelList>, pMapJ, pjIter)
                    {
                        label pointMaster =
                        (
                            cMapJ.findMaster
                            (
                                coupleMap::POINT,
                                pjIter.key()
                            )
                        );

                        if (pointMaster == -1)
                        {
                            continue;
                        }

                        const point& pointJ = sMeshJ.points_[pjIter.key()];

                        if (magSqr(pointI - pointJ) < mergeTol)
                        {
                            if (debug > 2)
                            {
                                Pout<< " Using mapped point: " << pjIter.key()
                                    << " :: " << pointJ
                                    << " procs: " << pjIter()
                                    << " master: " << pointMaster
                                    << endl;
                            }

                            // Assign to existing master point
                            pIter() = pointMaster;

                            // Update maps for the point
                            cMap.mapSlave
                            (
                                coupleMap::POINT,
                                pIter(), pIter.key()
                            );

                            cMap.mapMaster
                            (
                                coupleMap::POINT,
                                pIter.key(), pIter()
                            );

                            foundMaster = true;
                            break;
                        }
                    }

                    if (foundMaster)
                    {
                        break;
                    }
                }
            }

            forAllConstIter(Map<label>, cEdges, eIter)
            {
                label cSe = eIter.key();

                Map<labelList>::const_iterator eIt = pEdgeMap.find(cSe);

                if (eIt != pEdgeMap.end())
                {
                    const labelList& procs = eIt();

                    forAll(procs, procI)
                    {
                        label nProcNo = procs[procI];

                        if (priority(nProcNo, greaterEqOp<label>(), myProcNo))
                        {
                            continue;
                        }

                        if (debug > 3)
                        {
                            Pout<< nl << nl
                                << " Edge: " << cSe
                                << " :: " << sMesh.edges_[cSe]
                                << " is talking to processors: " << procs
                                << " so bailing out."
                                << endl;
                        }

                        map.type() = -2;

                        return map;
                    }
                }

                const edge& checkEdge = sMesh.edges_[cSe];

                forAll(checkEdge, pointI)
                {
                    label cSp = checkEdge[pointI];

                    Map<labelList>::const_iterator pIt = pPointMap.find(cSp);

                    if (pIt == pPointMap.end())
                    {
                        continue;
                    }

                    const labelList& procs = pIt();

                    forAll(procs, procI)
                    {
                        label nProcNo = procs[procI];

                        if (priority(nProcNo, greaterEqOp<label>(), myProcNo))
                        {
                            continue;
                        }

                        if (debug > 3)
                        {
                            Pout<< nl << nl
                                << " Points: " << nl
                                << "  Edge[0]: " << checkEdge[0] << nl
                                << "  Edge[1]: " << checkEdge[1] << nl
                                << " Processor conflict: " << nl
                                << "  Point: " << cSp << nl
                                << "  Processors: " << procs
                                << endl;
                        }

                        map.type() = -2;

                        return map;
                    }
                }
            }
        }
    }

    // Perform a few debug calls before modifications
    if (debug > 1)
    {
        Pout<< nl << nl
            << "Inserting cell(s) around coupled "
            << (is2D() ? "face: " : "edge: ") << mIndex
            << endl;
    }

    // Check to see if any new processor
    // patches need to be created
    forAllIter(Map<label>, createPatch, procIter)
    {
        label neiProcNo = procIter.key();

        // Create a new processor patch
        procIter() = createProcessorPatch(neiProcNo);

        // Renumber edge conversion patches
        forAll(masterConvertEdgePatch, pI)
        {
            Map<label>& convertMap = masterConvertEdgePatch[pI];

            forAllIter(Map<label>, convertMap, mIter)
            {
                if (mIter() < 0)
                {
                    // Back out the neighbouring processor ID
                    label proc = Foam::mag(mIter() + 2);

                    if (proc == neiProcNo)
                    {
                        // Set the index
                        mIter() = procIter();
                    }
                }
            }
        }

        // Renumber face conversion patches
        forAll(masterConvertFacePatch, pI)
        {
            Map<label>& convertMap = masterConvertFacePatch[pI];

            forAllIter(Map<label>, convertMap, mIter)
            {
                if (mIter() < 0)
                {
                    // Back out the neighbouring processor ID
                    label proc = Foam::mag(mIter() + 2);

                    if (proc == neiProcNo)
                    {
                        // Set the index
                        mIter() = procIter();
                    }
                }
            }
        }

        // Find index in processor list
        label pI = findIndex(procIndices_, neiProcNo);

        if (pI == -1)
        {
            Pout<< " Could not find index for processor: " << neiProcNo
                << " in indices: " << procIndices_
                << abort(FatalError);
        }

        // Create patch on subMesh
        dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

        slaveConvertPatch[pI] = mesh.createProcessorPatch(Pstream::myProcNo());
    }

    // Build a list of mapping entities on this processor
    dynamicLabelList mapCells(10);

    if (is2D())
    {
        label own = owner_[mIndex];

        mapCells.append(own);
    }
    else
    {
        const labelList& eFaces = edgeFaces_[mIndex];

        forAll(eFaces, faceI)
        {
            label own = owner_[eFaces[faceI]];
            label nei = neighbour_[eFaces[faceI]];

            if (findIndex(mapCells, own) == -1)
            {
                mapCells.append(own);
            }

            if (nei > -1)
            {
                if (findIndex(mapCells, nei) == -1)
                {
                    mapCells.append(nei);
                }
            }
        }
    }

    // Loop through insertion cells and
    // create an equivalent on this mesh
    forAll(procIndices_, pI)
    {
        const label cplEnum = is2D() ? coupleMap::FACE : coupleMap::EDGE;

        // Fetch reference to maps
        const coupleMap& cMap = recvMeshes_[pI].map();
        const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

        label sIndex = -1;

        if ((sIndex = cMap.findSlave(cplEnum, mIndex)) == -1)
        {
            continue;
        }

        // Fetch references
        Map<label>& procCellMap = cellsToInsert[pI];

        forAllIter(Map<label>, procCellMap, cIter)
        {
            const cell& checkCell = mesh.cells_[cIter.key()];

            scalar sLengthScale = -1.0;

            if (edgeRefinement_)
            {
                sLengthScale = mesh.lengthScale_[cIter.key()];
            }

            // Add an empty cell for now, and update
            // with face information at a later stage.
            cIter() = insertCell(cell(checkCell.size()), sLengthScale);

            // Initialize a face counter
            nCellFaces.insert(cIter(), 0);

            if (debug > 3)
            {
                Pout<< " Map cell: " << cIter()
                    << " for cell: " << cIter.key()
                    << " on proc: " << procIndices_[pI]
                    << endl;
            }

            // Add this cell to the map.
            map.addCell(cIter(), mapCells);

            // Set basic mapping for this cell
            setCellMapping(cIter(), mapCells);
        }

        // Push operation for the slave into coupleMap
        cMap.pushOperation
        (
            sIndex,
            coupleMap::REMOVE_CELL
        );
    }

    // Build a list of edges that need to be converted to interior.
    // - Do this by looking at edges of master face conversion candidates.
    // - Some edges may not need conversion, but deal with this later.
    forAll(facesToConvert, pI)
    {
        // Fetch reference to subMesh
        const coupleMap& cMap = recvMeshes_[pI].map();
        const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

        const Map<label>& procFaceMap = facesToConvert[pI];

        forAllConstIter(Map<label>, procFaceMap, fIter)
        {
            const labelList& mfEdges = faceEdges_[fIter()];
            const labelList& sfEdges = mesh.faceEdges_[fIter.key()];

            forAll(sfEdges, edgeI)
            {
                if (edgesToInsert[pI][sfEdges[edgeI]] != -1)
                {
                    // Already mapped this edge. Move on.
                    continue;
                }

                // Configure the comparison edge
                const edge& sEdge = mesh.edges_[sfEdges[edgeI]];

                label cMs = cMap.findMaster(coupleMap::POINT, sEdge[0]);
                label cMe = cMap.findMaster(coupleMap::POINT, sEdge[1]);

                edge cEdge(cMs, cMe);

                forAll(mfEdges, edgeJ)
                {
                    const edge& mEdge = edges_[mfEdges[edgeJ]];

                    if (mEdge == cEdge)
                    {
                        edgesToInsert[pI][sfEdges[edgeI]] = mfEdges[edgeJ];
                        edgesToConvert[pI][sfEdges[edgeI]] = mfEdges[edgeJ];
                        break;
                    }
                }
            }
        }
    }

    // Insert all points, with merging if necessary
    forAll(pointsToInsert, pI)
    {
        Map<label>& procPointMapI = pointsToInsert[pI];

        // Fetch reference to subMesh
        const coupleMap& cMapI = recvMeshes_[pI].map();
        const dynamicTopoFvMesh& meshI = recvMeshes_[pI].subMesh();

        forAllIter(Map<label>, procPointMapI, pItI)
        {
            if (pItI() != -1)
            {
                continue;
            }

            const point& pointI = meshI.points_[pItI.key()];

            bool foundMerge = false;

            forAll(pointsToInsert, pJ)
            {
                if (pJ == pI)
                {
                    continue;
                }

                Map<label>& procPointMapJ = pointsToInsert[pJ];

                // Fetch reference to subMesh
                const coupleMap& cMapJ = recvMeshes_[pJ].map();
                const dynamicTopoFvMesh& meshJ = recvMeshes_[pJ].subMesh();

                // Compare points with this processor
                forAllIter(Map<label>, procPointMapJ, pItJ)
                {
                    const point& pointJ = meshJ.points_[pItJ.key()];

                    if (magSqr(pointI - pointJ) < mergeTol)
                    {
                        if (pItJ() == -1)
                        {
                            // Make a merge entry
                            label mergePointIndex =
                            (
                                insertPoint
                                (
                                    pointJ,
                                    meshJ.oldPoints_[pItJ.key()],
                                    labelList(1, -1)
                                )
                            );

                            pItI() = mergePointIndex;
                            pItJ() = mergePointIndex;

                            // Update maps for the new point
                            cMapI.mapSlave
                            (
                                coupleMap::POINT,
                                pItI(), pItI.key()
                            );

                            cMapI.mapMaster
                            (
                                coupleMap::POINT,
                                pItI.key(), pItI()
                            );

                            cMapJ.mapSlave
                            (
                                coupleMap::POINT,
                                pItJ(), pItJ.key()
                            );

                            cMapJ.mapMaster
                            (
                                coupleMap::POINT,
                                pItJ.key(), pItJ()
                            );

                            if (debug > 2)
                            {
                                Pout<< " Map point: " << mergePointIndex
                                    << " pointI: " << pItI.key()
                                    << " procI: " << procIndices_[pI]
                                    << " pointJ: " << pItJ.key()
                                    << " procJ: " << procIndices_[pJ]
                                    << endl;
                            }

                            // Add this point to the map.
                            map.addPoint(mergePointIndex);
                        }
                        else
                        {
                            // Point appears to have been inserted
                            // by a previous operation. Map to it.
                            if (debug > 2)
                            {
                                Pout<< " Inserted point: " << pItJ()
                                    << " pointI: " << pItI.key()
                                    << " procI: " << procIndices_[pI]
                                    << " pointJ: " << pItJ.key()
                                    << " procJ: " << procIndices_[pJ]
                                    << endl;
                            }

                            // Set the entry
                            pItI() = pItJ();

                            // Update maps for the point
                            cMapI.mapSlave
                            (
                                coupleMap::POINT,
                                pItI(), pItI.key()
                            );

                            cMapI.mapMaster
                            (
                                coupleMap::POINT,
                                pItI.key(), pItI()
                            );
                        }

                        foundMerge = true;
                        break;
                    }
                }
            }

            if (foundMerge)
            {
                continue;
            }

            // Add a new unique point
            label newPointIndex =
            (
                insertPoint
                (
                    pointI,
                    meshI.oldPoints_[pItI.key()],
                    labelList(1, -1)
                )
            );

            // Set the entry
            pItI() = newPointIndex;

            // Update maps for the new point
            cMapI.mapSlave(coupleMap::POINT, pItI(), pItI.key());
            cMapI.mapMaster(coupleMap::POINT, pItI.key(), pItI());

            if (debug > 2)
            {
                Pout<< " Map point: " << newPointIndex
                    << " for point: " << pItI.key()
                    << " on proc: " << procIndices_[pI]
                    << endl;
            }

            // Add this point to the map.
            map.addPoint(newPointIndex);
        }
    }

    // Occassionally, inserted edges may already be present.
    // Ensure that edges are not added twice.
    if (is3D())
    {
        forAll(facesToInsert, pI)
        {
            // Fetch reference to subMesh
            const coupleMap& cMap = recvMeshes_[pI].map();
            const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

            const Map<label>& procFaceMap = facesToInsert[pI];

            forAllConstIter(Map<label>, procFaceMap, fIter)
            {
                const labelList& sfEdges = mesh.faceEdges_[fIter.key()];

                forAll(sfEdges, edgeI)
                {
                    label seIndex = sfEdges[edgeI];

                    if (edgesToInsert[pI][seIndex] != -1)
                    {
                        // Already mapped this edge. Move on.
                        continue;
                    }

                    label cMe = cMap.findMaster(coupleMap::EDGE, seIndex);

                    if (cMe > -1)
                    {
                        if (debug > 2)
                        {
                            Pout<< " Found master edge: " << cMe
                                << " :: " << edges_[cMe]
                                << " for slave edge: " << seIndex
                                << " :: " << mesh.edges_[seIndex]
                                << " on proc: " << procIndices_[pI]
                                << endl;
                        }

                        edgesToInsert[pI][seIndex] = cMe;
                        edgesToConvert[pI][seIndex] = cMe;

                        continue;
                    }

                    // Check for point-only coupling
                    const edge& sEdge = mesh.edges_[seIndex];

                    edge cEdge
                    (
                        cMap.findMaster(coupleMap::POINT, sEdge[0]),
                        cMap.findMaster(coupleMap::POINT, sEdge[1])
                    );

                    if (cEdge[0] > -1 && cEdge[1] > -1)
                    {
                        label meIndex = -1;

                        // Look at pointEdges info for a boundary edge
                        const labelList& pEdges = pointEdges_[cEdge[0]];

                        forAll(pEdges, edgeJ)
                        {
                            if (edges_[pEdges[edgeJ]] == cEdge)
                            {
                                if (whichEdgePatch(pEdges[edgeJ]) > -1)
                                {
                                    meIndex = pEdges[edgeJ];
                                    break;
                                }
                            }
                        }

                        if (meIndex > -1)
                        {
                            if (debug > 2)
                            {
                                Pout<< " Found master edge: " << meIndex
                                    << " :: " << edges_[meIndex]
                                    << " for slave edge: " << seIndex
                                    << " :: " << mesh.edges_[seIndex]
                                    << " on proc: " << procIndices_[pI]
                                    << endl;
                            }

                            edgesToInsert[pI][seIndex] = meIndex;

                            if (edgesToConvert[pI].found(seIndex))
                            {
                                edgesToConvert[pI][seIndex] = meIndex;
                            }
                            else
                            {
                                edgesToConvert[pI].insert(seIndex, meIndex);
                            }
                        }
                    }
                }
            }
        }
    }

    // Write out some debug info
    if (debug > 2)
    {
        forAll(cellsToInsert, pI)
        {
            label procIndex = procIndices_[pI];

            // Fetch reference to subMesh
            const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

            // Define a convenience output macro
#           define writeOutputVTK(fName, eMap, eType, insert)                  \
            {                                                                  \
                labelList entities(eMap.size());                               \
                                                                               \
                label nEnt = 0;                                                \
                                                                               \
                forAllConstIter(Map<label>, eMap, eIter)                       \
                {                                                              \
                    if (eIter() == -1 && insert)                               \
                    {                                                          \
                        entities[nEnt++] = eIter.key();                        \
                    }                                                          \
                    else                                                       \
                    if (eIter() > -1 && !insert)                               \
                    {                                                          \
                        entities[nEnt++] = eIter.key();                        \
                    }                                                          \
                }                                                              \
                                                                               \
                entities.setSize(nEnt);                                        \
                                                                               \
                if (entities.size())                                           \
                {                                                              \
                    mesh.writeVTK                                              \
                    (                                                          \
                        word(fName) + '_'                                      \
                      + Foam::name(mIndex) + '_'                               \
                      + Foam::name(procIndex),                                 \
                        entities,                                              \
                        eType, false, true                                     \
                    );                                                         \
                }                                                              \
            }

            writeOutputVTK("edgesToConvert", edgesToConvert[pI], 1, false)
            writeOutputVTK("facesToConvert", facesToConvert[pI], 2, false)
            writeOutputVTK("pointsToConvert", pointsToInsert[pI], 0, false)
            writeOutputVTK("pointsToInsert", pointsToInsert[pI], 0, true)
            writeOutputVTK("edgesToInsert", edgesToInsert[pI], 1, true)
            writeOutputVTK("facesToInsert", facesToInsert[pI], 2, true)
            writeOutputVTK("cellsToInsert", cellsToInsert[pI], 3, false)
        }
    }

    // Add edges to the mesh, noting that all
    // required points have already been added.
    forAll(edgesToInsert, pI)
    {
        // Fetch reference to subMesh
        const coupleMap& cMap = recvMeshes_[pI].map();
        const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

        Map<label>& procEdgeMap = edgesToInsert[pI];

        forAllIter(Map<label>, procEdgeMap, eIter)
        {
            // Skip if the edge is a conversion candidate
            if (eIter() != -1)
            {
                continue;
            }

            const edge& sEdge = mesh.edges_[eIter.key()];

            edge newEdge
            (
                cMap.findMaster(coupleMap::POINT, sEdge[0]),
                cMap.findMaster(coupleMap::POINT, sEdge[1])
            );

            // Ensure that this edge hasn't been added before
            label neIndex = -1;
            bool foundDuplicate = false;

            const List<objectMap>& addedEdges = map.addedEdgeList();

            forAll(addedEdges, edgeI)
            {
                neIndex = addedEdges[edgeI].index();

                if (edges_[neIndex] == newEdge)
                {
                    foundDuplicate = true;
                    break;
                }
            }

            if (foundDuplicate)
            {
                // Note the duplicate index for later
                eIter() = neIndex;

                continue;
            }

            // Insert edge with null edgeFaces for now.
            // This can be corrected later.
            eIter() =
            (
                insertEdge
                (
                    masterConvertEdgePatch[pI][eIter.key()],
                    newEdge,
                    labelList(0)
                )
            );

            if (debug > 2)
            {
                Pout<< " Map edge: " << eIter() << "::" << newEdge
                    << " for edge: " << eIter.key() << "::" << sEdge
                    << endl;
            }

            // Add this edge to the map.
            map.addEdge(eIter());
        }
    }

    // Add faces to the mesh, noting that all required
    // points and edges have already been added.
    forAll(facesToInsert, pI)
    {
        // Fetch reference to subMesh and coupleMap
        const coupleMap& cMapI = recvMeshes_[pI].map();
        const dynamicTopoFvMesh& meshI = recvMeshes_[pI].subMesh();

        Map<label>& procFaceMap = facesToInsert[pI];

        forAllIter(Map<label>, procFaceMap, fI)
        {
            // Skip if the face is a conversion candidate
            if (fI() != -1)
            {
                continue;
            }

            const face& sFaceI = meshI.faces_[fI.key()];
            const labelList& sfEdges = meshI.faceEdges_[fI.key()];

            // Configure new / comparison faces
            face nF(sFaceI.size()), cF(sFaceI.size());
            labelList nFaceEdges(sfEdges.size());

            // Configure points from map
            forAll(nF, pointI)
            {
                nF[pointI] =
                (
                    cMapI.findMaster
                    (
                        coupleMap::POINT,
                        sFaceI[pointI]
                    )
                );
            }

            // This may be a face shared by two other processors.
            // Ensure that duplicates are added only once.
            label dupeOwnCell = -1;
            bool foundDuplicate = false, foundInsertedDuplicate = false;

            forAll(facesToInsert, pJ)
            {
                if (pJ == pI)
                {
                    continue;
                }

                // Fetch reference to subMesh and coupleMap
                const coupleMap& cMapJ = recvMeshes_[pJ].map();
                const dynamicTopoFvMesh& meshJ = recvMeshes_[pJ].subMesh();

                const Map<label>& procFaceMapJ = facesToInsert[pJ];

                forAllConstIter(Map<label>, procFaceMapJ, fJ)
                {
                    const face& sFaceJ = meshJ.faces_[fJ.key()];

                    // Discard dissimilar face sizes
                    if (sFaceJ.size() != cF.size())
                    {
                        continue;
                    }

                    // Configure points from map
                    forAll(cF, pointI)
                    {
                        cF[pointI] =
                        (
                            cMapJ.findMaster
                            (
                                coupleMap::POINT,
                                sFaceJ[pointI]
                            )
                        );
                    }

                    if (cF.size() == 3)
                    {
                        // Optimized triangular face comparison
                        if
                        (
                            triFace::compare
                            (
                                triFace(nF[0], nF[1], nF[2]),
                                triFace(cF[0], cF[1], cF[2])
                            )
                        )
                        {
                            foundDuplicate = true;
                        }
                    }
                    else
                    {
                        // Regular face compare
                        if (face::compare(nF, cF))
                        {
                            foundDuplicate = true;
                        }
                    }

                    if (foundDuplicate)
                    {
                        // Record the owner for posterity
                        dupeOwnCell =
                        (
                            cellsToInsert[pJ][meshJ.owner_[fJ.key()]]
                        );

                        // Was the duplicate face inserted before this?
                        if (fJ() != -1)
                        {
                            // Note the duplicate index for later
                            fI() = fJ();

                            // If patch conversion entries were made,
                            // remove them as well
                            if
                            (
                                convertPatchPoints[pI].found(fJ.key())
                             && convertPatchPoints[pJ].found(fI.key())
                            )
                            {
                                convertPatchPoints[pI].erase(fJ.key());
                                convertPatchPoints[pJ].erase(fI.key());
                            }

                            foundInsertedDuplicate = true;
                        }

                        break;
                    }
                }

                if (foundDuplicate)
                {
                    break;
                }
            }

            // Face was inserted before, so don't insert again
            if (foundInsertedDuplicate)
            {
                continue;
            }

            // Configure edges from edgesToInsert
            forAll(sfEdges, edgeI)
            {
                label meIndex = -1;

                // Configure with the appropriate edge
                if (edgesToInsert[pI].found(sfEdges[edgeI]))
                {
                    meIndex = edgesToInsert[pI][sfEdges[edgeI]];
                }

                if (meIndex == -1)
                {
                    // Something is wrong here.
                    Pout<< "  Could not find correspondence for slave edge: "
                        << sfEdges[edgeI]
                        << ":: " << meshI.edges_[sfEdges[edgeI]]
                        << nl << " mIndex: " << mIndex
                        << abort(FatalError);
                }

                nFaceEdges[edgeI] = meIndex;
            }

            // Determine patch, owner and neighbour for this face
            label nPatch = -1, nOwner = -1, nNeighbour = -1;

            const polyBoundaryMesh& slaveBoundary = meshI.boundaryMesh();

            label sfPatch = meshI.whichPatch(fI.key());
            label sFaceOwn = meshI.owner_[fI.key()];
            label sFaceNei = meshI.neighbour_[fI.key()];

            label mFaceOwn =
            (
                cellsToInsert[pI].found(sFaceOwn) ?
                cellsToInsert[pI][sFaceOwn] : -1
            );

            label mFaceNei =
            (
                cellsToInsert[pI].found(sFaceNei) ?
                cellsToInsert[pI][sFaceNei] : -1
            );

            // If a duplicate face was found, over-ride neighbour.
            // Face-flipping will be taken care of automatically.
            if (foundDuplicate)
            {
                mFaceNei = dupeOwnCell;
            }

            if (mFaceOwn != -1 && mFaceNei == -1)
            {
                // Boundary face already has correct orientation
                nOwner = mFaceOwn;
                nNeighbour = -1;

                // Determine patch
                nPatch = masterConvertFacePatch[pI][fI.key()];
            }
            else
            if (mFaceOwn == -1 && mFaceNei != -1)
            {
                // Boundary face is inverted. Flip it
                nF = nF.reverseFace();
                nOwner = mFaceNei;
                nNeighbour = -1;

                // Determine patch
                nPatch = masterConvertFacePatch[pI][fI.key()];
            }
            else
            if (mFaceOwn != -1 && mFaceNei != -1)
            {
                // Interior face. Check if a flip is necessary.
                if (mFaceNei < mFaceOwn)
                {
                    nF = nF.reverseFace();
                    nOwner = mFaceNei;
                    nNeighbour = mFaceOwn;
                }
                else
                {
                    nOwner = mFaceOwn;
                    nNeighbour = mFaceNei;
                }

                nPatch = -1;
            }
            else
            if (mFaceOwn == -1 && mFaceNei == -1)
            {
                // Something is wrong here.
                Pout<< "Could not find correct owner / neighbour info: " << nl
                    << " Face: " << nF << nl
                    << " Owner: " << mFaceOwn << nl
                    << " Neighbour: " << mFaceNei << nl
                    << " - Slave Face: " << sFaceI << nl
                    << " - Slave Patch: " << slaveBoundary[sfPatch].name() << nl
                    << " - Slave Owner: " << sFaceOwn << nl
                    << " - Slave Neighbour: " << sFaceNei << nl
                    << abort(FatalError);
            }

            // Insert the new face
            fI() =
            (
                insertFace
                (
                    nPatch,
                    nF,
                    nOwner,
                    nNeighbour,
                    nFaceEdges
                )
            );

            // Size up edgeFaces for each edge
            forAll(nFaceEdges, edgeI)
            {
                meshOps::sizeUpList
                (
                    fI(),
                    edgeFaces_[nFaceEdges[edgeI]]
                );
            }

            if (debug > 3)
            {
                Pout<< " Map face: " << fI() << "::" << nF
                    << " Own: " << nOwner << " Nei: " << nNeighbour
                    << " fE: " << nFaceEdges << nl
                    << " for face: " << fI.key() << "::" << sFaceI
                    << " Own: " << sFaceOwn << " Nei: " << sFaceNei
                    << " fE: " << sfEdges
                    << endl;
            }

            // Add this face to the map.
            map.addFace(fI());

            if (nPatch > -1 && getNeighbourProcessor(nPatch) == -1)
            {
                // Physical patch on subMesh.
                //  - Set an invalid number so that
                //    an entry is made in facesFromFaces,
                //    while faceParents is empty.
                setFaceMapping(fI(), labelList(1, -1));
            }
            else
            {
                // Interior / processor face
                setFaceMapping(fI());
            }

            // Update cells
            cells_[nOwner][nCellFaces[nOwner]++] = fI();

            if (nNeighbour > -1)
            {
                cells_[nNeighbour][nCellFaces[nNeighbour]++] = fI();
            }
        }
    }

    // Loop through conversion faces
    forAll(facesToConvert, pI)
    {
        // Fetch reference to subMesh
        const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

        const Map<label>& procFaceMap = facesToConvert[pI];

        forAllConstIter(Map<label>, procFaceMap, fIter)
        {
            label fIndex = fIter();

            // Fetch the owner information
            label mFaceOwn = owner_[fIndex];
            label sFaceOwn = mesh.owner_[fIter.key()];
            label newNeighbour = cellsToInsert[pI][sFaceOwn];

            // Insert a new internal face.
            //  - Orientation should be correct, because this is a boundary
            //    face converted to an interior, and adjacent to an added cell.
            //  - Add the new faceEdges from the existing face. This may contain
            //    edges that need to be converted, but that will be done later.
            label newFaceIndex =
            (
                insertFace
                (
                    -1,
                    face(faces_[fIndex]),
                    mFaceOwn,
                    newNeighbour,
                    labelList(faceEdges_[fIndex])
                )
            );

            // Update map
            map.addFace(newFaceIndex, labelList(1, fIndex));

            // Set mapping information
            setFaceMapping(newFaceIndex);

            // Update the owner cell
            meshOps::replaceLabel
            (
                fIndex,
                newFaceIndex,
                cells_[mFaceOwn]
            );

            // Update the neighbour cell
            cells_[newNeighbour][nCellFaces[newNeighbour]++] = newFaceIndex;

            // Replace edgeFaces with the new face index
            const labelList& fEdges = faceEdges_[newFaceIndex];

            forAll(fEdges, edgeI)
            {
                meshOps::replaceLabel
                (
                    fIndex,
                    newFaceIndex,
                    edgeFaces_[fEdges[edgeI]]
                );
            }

            // Remove the old boundary face
            removeFace(fIndex);

            // Update map
            map.removeFace(fIndex);

            // For 2D meshes, the boundary face gets converted
            // to an interior one. Note the index for further operations.
            if ((mIndex == fIndex) && is2D())
            {
                map.index() = newFaceIndex;
            }
        }
    }

    // Sequentially add any convert-patch operations
    forAll(convertPatchPoints, pI)
    {
        // Fetch reference to maps / mesh
        const coupleMap& cMap = recvMeshes_[pI].map();
        dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

        if (convertPatchPoints[pI].empty())
        {
            continue;
        }

        // Find the appropriate processor patch
        // for patch conversion operations
        label procPatch = -1;

        for (label patchI = 0; patchI < mesh.nPatches_; patchI++)
        {
            label neiProc = mesh.getNeighbourProcessor(patchI);

            if (neiProc == Pstream::myProcNo())
            {
                procPatch = patchI;
                break;
            }
        }

        if (procPatch == -1)
        {
            Pout<< " * * * insertCells() * * * " << nl
                << " Could not find patch on slave processor: "
                << procIndices_[pI] << nl
                << " convertPatchPoints: " << convertPatchPoints[pI]
                << abort(FatalError);
        }

        forAllConstIter(Map<Pair<point> >, convertPatchPoints[pI], fIter)
        {
            // Search slave mesh for faces
            const point& newCentre = fIter().first();
            const point& oldCentre = fIter().second();

            label replaceFace = -1;

            // Accumulate stats in case of failure
            scalar minDist = GREAT;
            vector minPoint = vector::zero;
            dynamicLabelList checkedFaces;

            if (debug > 2)
            {
                // Reserve for append
                checkedFaces.setCapacity(50);
            }

            // Loop through all boundary faces,
            // and compute / compare face centres
            label sTot = mesh.faces_.size(), sInt = mesh.nOldInternalFaces_;

            for (label faceI = sInt; faceI < sTot; faceI++)
            {
                const face& fCheck = mesh.faces_[faceI];

                if (fCheck.empty())
                {
                    continue;
                }

                label pIndex = mesh.whichPatch(faceI);

                if (mesh.getNeighbourProcessor(pIndex) == -1)
                {
                    continue;
                }

                // Compute face-centre
                vector fC = fCheck.centre(mesh.points_);

                // Compute tolerance
                scalar tol = mag(mesh.points_[fCheck[0]] - fC);
                scalar dist = mag(fC - newCentre);

                if (dist < (geomMatchTol_() * tol))
                {
                    replaceFace = faceI;
                    break;
                }
                else
                if (dist < minDist)
                {
                    minPoint = fC;
                    minDist = dist;

                    if (debug > 2)
                    {
                        checkedFaces.append(faceI);
                    }
                }
            }

            // Ensure that the face was found
            if (replaceFace == -1)
            {
                mesh.writeVTK
                (
                    "checkedFaces_"
                  + Foam::name(fIter.key()),
                    checkedFaces,
                    2, false, true
                );

                Pout<< " * * * insertCells() * * * " << nl
                    << " Convert patch Op failed." << nl
                    << " Face: " << fIter.key() << nl
                    << " minPoint: " << minPoint << nl
                    << " minDistance: " << minDist << nl
                    << " newCentre: " << newCentre << nl
                    << " oldCentre: " << oldCentre << nl
                    << abort(FatalError);
            }

            // Obtain a copy before adding the new face,
            // since the reference might become
            // invalid during list resizing.
            // Edges don't have to change, since they're
            // all on the boundary anyway.
            face newFace = mesh.faces_[replaceFace];
            label newOwn = mesh.owner_[replaceFace];
            labelList newFaceEdges = mesh.faceEdges_[replaceFace];

            label newFaceIndex =
            (
                mesh.insertFace
                (
                    procPatch,
                    newFace,
                    newOwn,
                    -1,
                    newFaceEdges
                )
            );

            // changeMap update for the new face is not necessary,
            // since it is on the slave subMesh

            // Set mapping information
            mesh.setFaceMapping(newFaceIndex);

            meshOps::replaceLabel
            (
                replaceFace,
                newFaceIndex,
                mesh.cells_[newOwn]
            );

            // Correct edgeFaces with the new face label.
            forAll(newFaceEdges, edgeI)
            {
                meshOps::replaceLabel
                (
                    replaceFace,
                    newFaceIndex,
                    mesh.edgeFaces_[newFaceEdges[edgeI]]
                );
            }

            if (debug > 3)
            {
                Pout<< " Pushing CONVERT_PATCH for face: " << replaceFace
                    << " :: " << mesh.faces_[replaceFace] << nl
                    << " with new point: " << fIter().first()
                    << " and old point: " << fIter().second()
                    << " on proc: " << procIndices_[pI]
                    << endl;
            }

            // Finally remove the face
            mesh.removeFace(replaceFace);

            cMap.pushOperation
            (
                replaceFace,
                coupleMap::CONVERT_PATCH,
                fIter().first(),
                fIter().second()
            );
        }
    }

    // Loop through conversion edges
    forAll(edgesToConvert, pI)
    {
        // Fetch references
        const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

        const Map<label>& procEdgeMap = edgesToConvert[pI];

        // Loop through conversion edges
        forAllConstIter(Map<label>, procEdgeMap, eIter)
        {
            label eIndex = eIter(), physPatch = -1;
            bool allInterior = true, keepEdge = true;

            if (eIndex < 0)
            {
                // Not a conversion edge. Search added edge list for index
                eIndex = edgesToInsert[pI][eIter.key()];

                if (eIndex < 0)
                {
                    // Something's wrong
                    Pout<< " Edge: " << eIter.key()
                        << " :: " << mesh.edges_[eIter.key()]
                        << " was never inserted."
                        << abort(FatalError);
                }
            }

            const labelList& eFaces = edgeFaces_[eIndex];

            if (eFaces.size())
            {
                forAll(eFaces, faceI)
                {
                    if (whichPatch(eFaces[faceI]) > -1)
                    {
                        allInterior = false;
                        break;
                    }
                }
            }
            else
            {
                // Edge was deleted before, so skip it
                allInterior = false;
            }

            if (allInterior)
            {
                physPatch = -1;
                keepEdge = false;
            }
            else
            if (eFaces.size())
            {
                // Check if patches need to be switched
                label edgePatch = whichEdgePatch(eIndex);

                // Is this edge on a processor patch?
                bool edgeOnProc = (getNeighbourProcessor(edgePatch) > -1);

                if (edgeOnProc)
                {
                    bool foundProc = false;

                    forAll(eFaces, faceI)
                    {
                        label fPatch = whichPatch(eFaces[faceI]);

                        if (fPatch == -1)
                        {
                            continue;
                        }

                        if (getNeighbourProcessor(fPatch) == -1)
                        {
                            // Note physical patch for later
                            physPatch = fPatch;
                        }
                        else
                        {
                            // Still on a processor patch.
                            foundProc = true;
                            break;
                        }
                    }

                    if (foundProc)
                    {
                        // Reset physPatch
                        physPatch = -1;
                    }
                }

                if (edgeOnProc && physPatch > -1)
                {
                    // Edge needs to be moved to another patch
                    keepEdge = false;
                }
                else
                if ((mIndex == eIndex) && is3D() && map.index() == -1)
                {
                    // Keep the edge, and note index for later
                    keepEdge = true;
                    map.index() = eIndex;
                }
            }

            if (keepEdge)
            {
                continue;
            }

            // This edge needs to be converted to interior / other patch
            label newEdgeIndex =
            (
                insertEdge
                (
                    physPatch,
                    edge(edges_[eIndex]),
                    labelList(edgeFaces_[eIndex])
                )
            );

            // Update map
            map.addEdge(newEdgeIndex, labelList(1, eIndex));

            // Update faceEdges information for all connected faces
            const labelList& neFaces = edgeFaces_[newEdgeIndex];

            forAll(neFaces, faceI)
            {
                meshOps::replaceLabel
                (
                    eIndex,
                    newEdgeIndex,
                    faceEdges_[neFaces[faceI]]
                );
            }

            // Remove the old boundary edge
            removeEdge(eIndex);

            // Update map
            map.removeEdge(eIndex);

            // For 3D meshes, the boundary edge gets converted
            // to an interior one. Note the index for further operations.
            if ((mIndex == eIndex) && is3D())
            {
                map.index() = newEdgeIndex;
            }

            // Replace the entry in edgesToInsert
            edgesToInsert[pI][eIter.key()] = newEdgeIndex;
        }
    }

    List<changeMap> slaveMaps(procIndices_.size());

    // Loop through all processors, and remove cells
    const label emptyMap = -7;

    forAll(cellsToInsert, pI)
    {
        const Map<label>& procCellMap = cellsToInsert[pI];

        if (procCellMap.empty())
        {
            // Set type to something recognizable
            slaveMaps[pI].type() = emptyMap;

            continue;
        }

        // Prepare a list of cells
        const labelList cellsToRemove = procCellMap.toc();

        // Fetch non-const reference to subMesh
        dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

        // Now remove cells for this processor
        slaveMaps[pI] =
        (
            mesh.removeCells
            (
                cellsToRemove,
                slaveConvertPatch[pI],
                "rc_" + Foam::name(mIndex)
            )
        );
    }

    // Now map entities from the removeCells operation
    forAll(slaveMaps, pI)
    {
        const changeMap& slaveMap = slaveMaps[pI];

        // Skip empty entities
        if (slaveMap.type() == emptyMap)
        {
            continue;
        }

        // Fetch references
        const coupleMap& cMap = recvMeshes_[pI].map();
        const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

        // Map faces from patches.
        const Map<label>& procFaceMap = facesToInsert[pI];

        // Check removed faces
        const labelList& rsfList = slaveMap.removedFaceList();

        forAll(rsfList, indexI)
        {
            label sfIndex = rsfList[indexI], mfIndex = -1;

            // Does an entry exist?
            mfIndex = cMap.findMaster(coupleMap::FACE, sfIndex);

            if (mfIndex > -1)
            {
                if (debug > 2)
                {
                    Pout<< " Removing map for master face: "
                        << mfIndex << " :: " << faces_[mfIndex]
                        << " Master map: "
                        << cMap.findSlave(coupleMap::FACE, mfIndex)
                        << " and slave face: " << sfIndex
                        << " on proc: " << procIndices_[pI]
                        << endl;
                }

                cMap.removeSlave(coupleMap::FACE, sfIndex);
                cMap.removeMaster(coupleMap::FACE, mfIndex);
            }
        }

        // Check added faces
        const List<objectMap>& asfList = slaveMap.addedFaceList();

        forAll(asfList, indexI)
        {
            label sfIndex = asfList[indexI].index(), mfIndex = -1;

            if (mesh.whichPatch(sfIndex) == slaveConvertPatch[pI])
            {
                // Configure a comparison face
                face cFace(mesh.faces_[sfIndex]);

                forAll(cFace, pointI)
                {
                    cFace[pointI] =
                    (
                        cMap.findMaster
                        (
                            coupleMap::POINT,
                            cFace[pointI]
                        )
                    );
                }

                bool foundMatch = false;

                forAllConstIter(Map<label>, procFaceMap, fIter)
                {
                    const face& mFace = faces_[fIter()];

                    // Discard dissimilar face sizes
                    if (mFace.size() != cFace.size())
                    {
                        continue;
                    }

                    if (cFace.size() == 3)
                    {
                        // Optimized triangular face comparison
                        if
                        (
                            triFace::compare
                            (
                                triFace(mFace[0], mFace[1], mFace[2]),
                                triFace(cFace[0], cFace[1], cFace[2])
                            )
                        )
                        {
                            mfIndex = fIter();
                            foundMatch = true;
                            break;
                        }
                    }
                    else
                    {
                        // Regular face compare
                        if (face::compare(mFace, cFace))
                        {
                            mfIndex = fIter();
                            foundMatch = true;
                            break;
                        }
                    }
                }

                if (foundMatch)
                {
                    if (debug > 2)
                    {
                        Pout<< " Found master face: " << mfIndex
                            << " :: " << faces_[mfIndex]
                            << " for slave face: " << sfIndex
                            << " :: " << mesh.faces_[sfIndex]
                            << " on proc: " << procIndices_[pI]
                            << endl;
                    }

                    // Update maps for the new face
                    cMap.mapSlave(coupleMap::FACE, mfIndex, sfIndex);
                    cMap.mapMaster(coupleMap::FACE, sfIndex, mfIndex);
                }
                else
                {
                    // Something is wrong here.
                    Pout<< " Could not find master face for: " << nl
                        << "  Slave face: " << sfIndex
                        << "  :: " << mesh.faces_[sfIndex] << nl
                        << "  cFace: " << cFace << nl
                        << "  on proc: " << procIndices_[pI] << nl
                        << abort(FatalError);
                }
            }
        }

        // Map edges for 3D meshes
        if (is3D())
        {
            // Map edges from patches.
            const Map<label>& procEdgeMap = edgesToInsert[pI];

            // Check removed edges
            const labelList& rseList = slaveMap.removedEdgeList();

            forAll(rseList, indexI)
            {
                label seIndex = rseList[indexI], meIndex = -1;

                // Does an entry exist?
                meIndex = cMap.findMaster(coupleMap::EDGE, seIndex);

                if (meIndex > -1)
                {
                    if (debug > 2)
                    {
                        Pout<< " Removing map for master edge: "
                            << meIndex << " :: " << edges_[meIndex]
                            << " Master map: "
                            << cMap.findSlave(coupleMap::EDGE, meIndex)
                            << " and slave edge: " << seIndex
                            << " on proc: " << procIndices_[pI]
                            << endl;
                    }

                    cMap.removeSlave(coupleMap::EDGE, seIndex);
                    cMap.removeMaster(coupleMap::EDGE, meIndex);
                }
            }

            // Check added edges
            const List<objectMap>& aseList = slaveMap.addedEdgeList();

            forAll(aseList, indexI)
            {
                label seIndex = aseList[indexI].index(), meIndex = -1;

                if (mesh.whichEdgePatch(seIndex) == slaveConvertPatch[pI])
                {
                    // Configure a comparison edge
                    edge cEdge(mesh.edges_[seIndex]);

                    cEdge[0] = cMap.findMaster(coupleMap::POINT, cEdge[0]);
                    cEdge[1] = cMap.findMaster(coupleMap::POINT, cEdge[1]);

                    bool foundMatch = false;

                    forAllConstIter(Map<label>, procEdgeMap, eIter)
                    {
                        const edge& mEdge = edges_[eIter()];

                        if (mEdge == cEdge)
                        {
                            meIndex = eIter();
                            foundMatch = true;
                            break;
                        }
                    }

                    if (foundMatch)
                    {
                        if (debug > 2)
                        {
                            Pout<< " Found master edge: " << meIndex
                                << " :: " << edges_[meIndex]
                                << " for slave edge: " << seIndex
                                << " :: " << mesh.edges_[seIndex]
                                << " on proc: " << procIndices_[pI]
                                << endl;
                        }

                        // Update maps for the new edge
                        cMap.mapSlave(coupleMap::EDGE, meIndex, seIndex);
                        cMap.mapMaster(coupleMap::EDGE, seIndex, meIndex);
                    }
                    else
                    {
                        // Something is wrong here.
                        Pout<< " Could not find master edge for: " << nl
                            << "  Slave edge: " << seIndex
                            << "  :: " << mesh.edges_[seIndex] << nl
                            << "  cEdge: " << cEdge << nl
                            << " on proc: " << procIndices_[pI] << nl
                            << abort(FatalError);
                    }
                }
            }
        }

        // Check removed points
        const labelList& rspList = slaveMap.removedPointList();

        forAll(rspList, indexI)
        {
            label spIndex = rspList[indexI], mpIndex = -1;

            // Does an entry exist?
            mpIndex = cMap.findMaster(coupleMap::POINT, spIndex);

            if (mpIndex > -1)
            {
                if (debug > 2)
                {
                    Pout<< " Removing map for master point: " << mpIndex
                        << " :: new: " << points_[mpIndex]
                        << " :: old: " << oldPoints_[mpIndex] << nl
                        << " Master point map: "
                        << cMap.findSlave(coupleMap::POINT, mpIndex)
                        << " and slave point: " << spIndex
                        << " on proc: " << procIndices_[pI]
                        << endl;
                }

                cMap.removeSlave(coupleMap::POINT, spIndex);
                cMap.removeMaster(coupleMap::POINT, mpIndex);
            }
        }
    }

    // Write out cells for debug, if necessary
    if (debug > 2 || map.index() < 0)
    {
        const List<objectMap>& acList = map.addedCellList();

        labelList addedCells(acList.size(), -1);

        forAll(acList, indexI)
        {
            addedCells[indexI] = acList[indexI].index();
        }

        // Write it out
        writeVTK
        (
            "insertCells(" + Foam::name(mIndex) + ')',
            addedCells, 3, false, true
        );
    }

    // Set mapping information for any added faces
    const List<objectMap>& afList = map.addedFaceList();

    forAll(afList, indexI)
    {
        label mfIndex = afList[indexI].index();
        label patch = whichPatch(mfIndex);
        label neiProc = getNeighbourProcessor(patch);

        // Disregard non-processor faces for coupled mapping
        if (neiProc == -1)
        {
            continue;
        }

        // Update couple maps, if necessary
        const face& mFace = faces_[mfIndex];

        forAll(procIndices_, pI)
        {
            // Skip face if not on this processor
            if (neiProc != procIndices_[pI])
            {
                continue;
            }

            const coupleMap& cMap = recvMeshes_[pI].map();
            const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

            // Configure a comparison face
            face cFace(mFace);

            forAll(cFace, pointI)
            {
                cFace[pointI] =
                (
                    cMap.findSlave
                    (
                        coupleMap::POINT,
                        cFace[pointI]
                    )
                );
            }

            // Ensure all points are found
            if (findIndex(cFace, -1) > -1)
            {
                continue;
            }

            // Loop through boundary faces on slave and look for a match
            label sTot = mesh.faces_.size(), sInt = mesh.nOldInternalFaces_;

            label sfIndex = -1;

            for (label sfI = sInt; sfI < sTot; sfI++)
            {
                const face& sFace = mesh.faces_[sfI];

                if (sFace.empty())
                {
                    continue;
                }

                label pIndex = mesh.whichPatch(sfI);

                if (mesh.getNeighbourProcessor(pIndex) == -1)
                {
                    continue;
                }

                if (is2D())
                {
                    if (face::compare(cFace, sFace))
                    {
                        sfIndex = sfI;
                    }
                }
                else
                {
                    if
                    (
                        triFace::compare
                        (
                            triFace(cFace[0], cFace[1], cFace[2]),
                            triFace(sFace[0], sFace[1], sFace[2])
                        )
                    )
                    {
                        sfIndex = sfI;
                    }
                }

                // Break out if we're done
                if (sfIndex > -1)
                {
                    if (debug > 2)
                    {
                        Pout<< " Matched master face: " << mfIndex
                            << " :: " << mFace
                            << " with slave: " << sfIndex
                            << " :: " << sFace
                            << " on proc: " << procIndices_[pI]
                            << endl;
                    }

                    // Found the slave. Add a map entry
                    cMap.mapSlave
                    (
                        coupleMap::FACE,
                        mfIndex,
                        sfIndex
                    );

                    cMap.mapMaster
                    (
                        coupleMap::FACE,
                        sfIndex,
                        mfIndex
                    );

                    break;
                }
            }

            if (sfIndex == -1)
            {
                Pout<< " Failed match for master face: " << mfIndex
                    << " :: " << mFace
                    << " on proc: " << procIndices_[pI]
                    << abort(FatalError);
            }

            if (is2D())
            {
                continue;
            }

            // Map edges on this face as well
            const labelList& mfEdges = faceEdges_[mfIndex];
            const labelList& sfEdges = mesh.faceEdges_[sfIndex];

            forAll(mfEdges, edgeI)
            {
                label meIndex = mfEdges[edgeI];
                const edge& mEdge = edges_[meIndex];

                // Configure a comparison edge
                edge cEdge
                (
                    cMap.findSlave(coupleMap::POINT, mEdge[0]),
                    cMap.findSlave(coupleMap::POINT, mEdge[1])
                );

                forAll(sfEdges, edgeJ)
                {
                    label seIndex = sfEdges[edgeJ];
                    const edge& sEdge = mesh.edges_[seIndex];

                    if (cEdge == sEdge)
                    {
                        if (debug > 2)
                        {
                            Pout<< " Matched master edge: " << meIndex
                                << " :: " << mEdge
                                << " with slave: " << seIndex
                                << " :: " << sEdge
                                << " on proc: " << procIndices_[pI]
                                << endl;
                        }

                        // Found the slave. Add a map entry
                        cMap.mapSlave
                        (
                            coupleMap::EDGE,
                            meIndex,
                            seIndex
                        );

                        cMap.mapMaster
                        (
                            coupleMap::EDGE,
                            seIndex,
                            meIndex
                        );

                        break;
                    }
                }
            }
        }
    }

    // Specify that the operation was successful
    map.type() = 1;

    // Return the changeMap
    return map;
}


// Remove the specified cells from the mesh,
// and add internal faces/edges to the specified patch
// - Returns a changeMap with a type specifying:
//     1: Operation was successful
//    -1: Operation failed
// - checkOnly performs a feasibility check and returns without modifications.
const changeMap dynamicTopoFvMesh::removeCells
(
    const labelList& cList,
    const label patch,
    const word& rcN,
    bool checkOnly
)
{
    if (cList.empty() || patch < 0)
    {
        if (cList.size())
        {
            writeVTK("removeCellsFailed_" + rcN, cList, 3, false, true);
        }

        FatalErrorIn
        (
            "const changeMap dynamicTopoFvMesh::removeCells"
            "(const labelList&, const label, const word&, bool)"
        )
            << " Wrong arguments. " << nl
            << " cList: " << cList << nl
            << " patch: " << patch << nl
            << abort(FatalError);
    }

    changeMap map;

    labelHashSet pointsToRemove, edgesToRemove, facesToRemove;
    Map<label> facesToConvert, edgesToConvert;

    // First loop through all cells and accumulate
    // a set of faces to be removed/converted.
    forAll(cList, cellI)
    {
        const cell& cellToCheck = cells_[cList[cellI]];

        forAll(cellToCheck, faceI)
        {
            label own = owner_[cellToCheck[faceI]];
            label nei = neighbour_[cellToCheck[faceI]];

            if (nei == -1)
            {
                if (!facesToRemove.found(cellToCheck[faceI]))
                {
                    facesToRemove.insert(cellToCheck[faceI]);
                }
            }
            else
            if
            (
                (findIndex(cList, own) != -1) &&
                (findIndex(cList, nei) != -1)
            )
            {
                if (!facesToRemove.found(cellToCheck[faceI]))
                {
                    facesToRemove.insert(cellToCheck[faceI]);
                }
            }
            else
            {
                facesToConvert.set(cellToCheck[faceI], -1);
            }
        }
    }

    // Add all edges as candidates for conversion.
    // Some of these will be removed altogether.
    forAllConstIter(labelHashSet, facesToRemove, fIter)
    {
        const labelList& fEdges = faceEdges_[fIter.key()];

        forAll(fEdges, edgeI)
        {
            if (whichEdgePatch(fEdges[edgeI]) == patch)
            {
                // Make an identical map
                edgesToConvert.set(fEdges[edgeI], fEdges[edgeI]);
            }
            else
            {
                edgesToConvert.set(fEdges[edgeI], -1);
            }

            if (is2D())
            {
                // Add all points as candidates for removal.
                // Some (or all) of these will be weeded out.
                const edge& eCheck = edges_[fEdges[edgeI]];

                pointsToRemove.set(eCheck[0]);
                pointsToRemove.set(eCheck[1]);
            }
        }
    }

    forAllConstIter(Map<label>, facesToConvert, fIter)
    {
        const labelList& fEdges = faceEdges_[fIter.key()];

        forAll(fEdges, edgeI)
        {
            if (whichEdgePatch(fEdges[edgeI]) == patch)
            {
                // Make an identical map
                edgesToConvert.set(fEdges[edgeI], fEdges[edgeI]);
            }
            else
            {
                edgesToConvert.set(fEdges[edgeI], -1);
            }
        }
    }

    // Build a list of edges to be removed.
    forAllConstIter(Map<label>, edgesToConvert, eIter)
    {
        const labelList& eFaces = edgeFaces_[eIter.key()];

        bool allRemove = true;

        forAll(eFaces, faceI)
        {
            if (!facesToRemove.found(eFaces[faceI]))
            {
                allRemove = false;
                break;
            }
        }

        if (allRemove)
        {
            if (!edgesToRemove.found(eIter.key()))
            {
                edgesToRemove.insert(eIter.key());
            }
        }
    }

    // Weed-out the conversion list.
    forAllConstIter(labelHashSet, edgesToRemove, eIter)
    {
        edgesToConvert.erase(eIter.key());
    }

    // Build a set of points to be removed.
    if (is2D())
    {
        forAllConstIter(Map<label>, edgesToConvert, eIter)
        {
            const edge& eCheck = edges_[eIter.key()];

            if (pointsToRemove.found(eCheck[0]))
            {
                pointsToRemove.erase(eCheck[0]);
            }

            if (pointsToRemove.found(eCheck[1]))
            {
                pointsToRemove.erase(eCheck[1]);
            }
        }
    }
    else
    {
        forAllConstIter(labelHashSet, edgesToRemove, eIter)
        {
            const edge& edgeToCheck = edges_[eIter.key()];

            forAll(edgeToCheck, pointI)
            {
                const labelList& pEdges = pointEdges_[edgeToCheck[pointI]];

                bool allRemove = true;

                forAll(pEdges, edgeI)
                {
                    if (!edgesToRemove.found(pEdges[edgeI]))
                    {
                        allRemove = false;
                        break;
                    }
                }

                if (allRemove)
                {
                    if (!pointsToRemove.found(edgeToCheck[pointI]))
                    {
                        pointsToRemove.insert(edgeToCheck[pointI]);
                    }
                }
            }
        }
    }

    forAllIter(Map<label>, edgesToConvert, eIter)
    {
        const labelList& eFaces = edgeFaces_[eIter.key()];

        label nConvFaces = 0;

        forAll(eFaces, faceI)
        {
            if (facesToConvert.found(eFaces[faceI]))
            {
                nConvFaces++;
            }
        }

        if (nConvFaces > 2)
        {
            Pout<< "Invalid conversion. Bailing out." << endl;
            return map;
        }
    }

    // Are we only performing checks?
    if (checkOnly)
    {
        // Make necessary map entries
        forAllConstIter(Map<label>, facesToConvert, fIter)
        {
            map.removeFace(fIter.key());
        }

        forAllConstIter(Map<label>, edgesToConvert, eIter)
        {
            map.removeEdge(eIter.key());
        }

        forAllConstIter(labelHashSet, facesToRemove, fIter)
        {
            map.removeFace(fIter.key());
        }

        forAllConstIter(labelHashSet, edgesToRemove, eIter)
        {
            map.removeEdge(eIter.key());
        }

        forAllConstIter(labelHashSet, pointsToRemove, pIter)
        {
            map.removePoint(pIter.key());
        }

        forAll(cList, cellI)
        {
            map.removeCell(cList[cellI]);
        }

        // Specify that the operation was successful
        map.type() = 1;

        return map;
    }

    // Perform a few debug calls before modifications
    if (debug > 1)
    {
        Pout<< nl << nl
            << " Removing cell(s): " << cList
            << " and adding internal faces / edges to patch: "
            << boundaryMesh()[patch].name()
            << endl;
    }

    // Write out candidates for post-processing
    if (debug > 2)
    {
        writeVTK("pointsToRemove_" + rcN, pointsToRemove.toc(), 0, false, true);
        writeVTK("edgesToRemove_" + rcN, edgesToRemove.toc(), 1, false, true);
        writeVTK("facesToRemove_" + rcN, facesToRemove.toc(), 2, false, true);
        writeVTK("cellsToRemove_" + rcN, cList, 3, false, true);
        writeVTK("edgesToConvert_" + rcN, edgesToConvert.toc(), 1, false, true);
        writeVTK("facesToConvert_" + rcN, facesToConvert.toc(), 2, false, true);
    }

    // Loop through all faces for conversion, check orientation
    // and create new faces in their place.
    forAllIter(Map<label>, facesToConvert, fIter)
    {
        // Check if this internal face is oriented properly.
        face newFace;
        label newOwner = -1;
        labelList fEdges = faceEdges_[fIter.key()];

        if (findIndex(cList, neighbour_[fIter.key()]) != -1)
        {
            // Orientation is correct
            newFace = faces_[fIter.key()];
            newOwner = owner_[fIter.key()];
        }
        else
        if (findIndex(cList, owner_[fIter.key()]) != -1)
        {
            // Face is to be reversed.
            newFace = faces_[fIter.key()].reverseFace();
            newOwner = neighbour_[fIter.key()];

            setFlip(fIter.key());
        }
        else
        {
            // Something's terribly wrong
            FatalErrorIn
            (
                "\n"
                "const changeMap dynamicTopoFvMesh::removeCells\n"
                "(\n"
                "    const labelList& cList,\n"
                "    const label patch,\n"
                "    const word& rcN,\n"
                "    bool checkOnly\n"
                ")\n"
            )
                << nl << " Invalid mesh. "
                << abort(FatalError);
        }

        // Insert the reconfigured face at the boundary.
        // faceEdges will be corrected later.
        fIter() =
        (
            insertFace
            (
                patch,
                newFace,
                newOwner,
                -1,
                fEdges
            )
        );

        // Add this face to the map.
        map.addFace(fIter());

        // Set mapping information
        //  - But where to map from?
        setFaceMapping(fIter());

        // Replace cell with the new face label
        meshOps::replaceLabel
        (
            fIter.key(),
            fIter(),
            cells_[newOwner]
        );

        // Remove the internal face.
        removeFace(fIter.key());

        // Update map
        map.removeFace(fIter.key());
    }

    // Create a new edge for each converted edge
    forAllIter(Map<label>, edgesToConvert, eIter)
    {
        if (eIter() == -1)
        {
            // Create copies before appending.
            edge newEdge = edges_[eIter.key()];
            labelList eFaces = edgeFaces_[eIter.key()];

            eIter() =
            (
                insertEdge
                (
                    patch,
                    newEdge,
                    eFaces
                )
            );

            // Add this edge to the map.
            map.addEdge(eIter());

            // Remove the edge
            removeEdge(eIter.key());

            // Update map
            map.removeEdge(eIter.key());
        }
    }

    // Loop through all faces for conversion, and replace edgeFaces.
    forAllConstIter(Map<label>, facesToConvert, fIter)
    {
        // Make a copy, because this list is going to
        // be modified within this loop.
        labelList fEdges = faceEdges_[fIter()];

        forAll(fEdges, edgeI)
        {
            if (edgesToConvert.found(fEdges[edgeI]))
            {
                meshOps::replaceLabel
                (
                    fIter.key(),
                    fIter(),
                    edgeFaces_[edgesToConvert[fEdges[edgeI]]]
                );

                meshOps::replaceLabel
                (
                    fEdges[edgeI],
                    edgesToConvert[fEdges[edgeI]],
                    faceEdges_[fIter()]
                );
            }
        }
    }

    // Loop through all edges for conversion, and size-down edgeFaces.
    forAllConstIter(Map<label>, edgesToConvert, eIter)
    {
        // Make a copy, because this list is going to
        // be modified within this loop.
        labelList eFaces = edgeFaces_[eIter()];

        forAll(eFaces, faceI)
        {
            if (facesToRemove.found(eFaces[faceI]))
            {
                meshOps::sizeDownList
                (
                    eFaces[faceI],
                    edgeFaces_[eIter()]
                );
            }

            // Replace old edges with new ones.
            labelList& fEdges = faceEdges_[eFaces[faceI]];

            forAll(fEdges, edgeI)
            {
                if (edgesToConvert.found(fEdges[edgeI]))
                {
                    fEdges[edgeI] = edgesToConvert[fEdges[edgeI]];
                }
            }
        }
    }

    // Remove unwanted faces
    forAllConstIter(labelHashSet, facesToRemove, fIter)
    {
        removeFace(fIter.key());

        // Update map
        map.removeFace(fIter.key());
    }

    // Remove unwanted edges
    forAllConstIter(labelHashSet, edgesToRemove, eIter)
    {
        removeEdge(eIter.key());

        // Update map
        map.removeEdge(eIter.key());
    }

    // Remove unwanted points
    forAllConstIter(labelHashSet, pointsToRemove, pIter)
    {
        removePoint(pIter.key());

        // Update map
        map.removePoint(pIter.key());
    }

    // Remove all cells
    forAll(cList, cellI)
    {
        removeCell(cList[cellI]);

        // Update map
        map.removeCell(cList[cellI]);
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Specify that the operation was successful
    map.type() = 1;

    return map;
}


// Handle topology changes for coupled patches
void dynamicTopoFvMesh::handleCoupledPatches
(
    labelHashSet& entities
)
{
    // Initialize coupled patch connectivity for topology modifications.
    initCoupledConnectivity(this);

    if (patchCoupling_.empty() && procIndices_.empty())
    {
        return;
    }

    // Move coupled subMeshes
    moveCoupledSubMeshes();

    // Exchange length-scale buffers across processors.
    exchangeLengthBuffers();

    if (debug)
    {
        // Check coupled-patch sizes first.
        forAll(patchCoupling_, patchI)
        {
            if (patchCoupling_(patchI))
            {
                const coupleMap& cMap = patchCoupling_[patchI].map();

                label mSize = patchSizes_[cMap.masterIndex()];
                label sSize = patchSizes_[cMap.slaveIndex()];

                if (mSize != sSize)
                {
                    Pout<< "Coupled patch-count is inconsistent." << nl
                        << " Master Patch: " << cMap.masterIndex()
                        << " Count: " << mSize << nl
                        << " Slave Patch: " << cMap.slaveIndex()
                        << " Count: " << sSize
                        << endl;

                    FatalErrorIn
                    (
                        "void dynamicTopoFvMesh::handleCoupledPatches()"
                    )
                        << " Failures were found in connectivity"
                        << " prior to coupled topo-changes."
                        << abort(FatalError);
                }
            }
        }

        Info<< "Handling coupled patches...";
    }

    if (Pstream::parRun())
    {
        // Calculate mapping offsets
        const polyBoundaryMesh& boundary = boundaryMesh();

        // Determine number of physical (non-processor) patches
        label nPhysical = 0;

        forAll(boundary, patchI)
        {
            if (isA<processorPolyPatch>(boundary[patchI]))
            {
                continue;
            }

            nPhysical++;
        }

        // Fetch number of processors
        label nProcs = procIndices_.size();

        // Allocate cell, face and patch offsets
        labelList cellSizes(nProcs, 0), cellStarts(nProcs, 0);
        labelList faceSizes(nProcs, 0), faceStarts(nProcs, 0);
        labelListList patchSizes(nProcs, labelList(nPhysical, 0));
        labelListList patchStarts(nProcs, labelList(nPhysical, 0));

        label nTotalCells = nOldCells_, nTotalIntFaces = nOldInternalFaces_;
        labelList nTotalPatchFaces(SubList<label>(oldPatchSizes_, nPhysical));

        forAll(procIndices_, pI)
        {
            const coupleMap& cMap = recvMeshes_[pI].map();

            // Fetch size from subMesh
            label nCells = cMap.nEntities(coupleMap::CELL);
            label nIntFaces = cMap.nEntities(coupleMap::INTERNAL_FACE);

            // Set size / offset for this processor
            cellSizes[pI] = nCells;
            cellStarts[pI] = nTotalCells;

            faceSizes[pI] = nIntFaces;
            faceStarts[pI] = nTotalIntFaces;

            // Update count
            nTotalCells += nCells;
            nTotalIntFaces += nIntFaces;

            // Fetch patch sizes from subMesh
            const labelList& nPatchFaces =
            (
                cMap.entityBuffer(coupleMap::FACE_SIZES)
            );

            // Loop over physical patches
            forAll(nTotalPatchFaces, patchI)
            {
                // Fetch patch size from subMesh
                label nFaces = nPatchFaces[patchI];

                // Set patch size / offset for this processor
                patchSizes[pI][patchI] = nFaces;
                patchStarts[pI][patchI] = nTotalPatchFaces[patchI];

                // Update patch count
                nTotalPatchFaces[patchI] += nFaces;
            }
        }

        // Set sizes / starts in mapper
        mapper_->setOffsets
        (
            cellSizes,
            cellStarts,
            faceSizes,
            faceStarts,
            patchSizes,
            patchStarts
        );

        if (debug > 3)
        {
            SubList<label> physicalPatches(oldPatchSizes_, nPhysical);

            Pout<< " procIndices: " << procIndices_ << nl
                << " nCells: " << nOldCells_ << nl
                << " proc cellSizes: " << cellSizes << nl
                << " cellStarts: " << cellStarts << nl
                << " proc faceSizes: " << faceSizes << nl
                << " faceStarts: " << faceStarts << nl
                << " patchSizes: " << physicalPatches << nl
                << " proc patchSizes: " << patchSizes << nl
                << " patchStarts: " << patchStarts << endl;
        }
    }

    // Optionally switch off topo-changes
    // on processor patches at run-time
    bool procTopoChanges = true;

    const dictionary& meshSubDict = dict_.subDict("dynamicTopoFvMesh");

    if (meshSubDict.found("procTopoChanges") || mandatory_)
    {
        procTopoChanges = readBool(meshSubDict.lookup("procTopoChanges"));
    }

    // Set coupled modifications.
    setCoupledModification();

    // Loop through the coupled stack and perform changes.
    if (edgeRefinement_)
    {
        // Initialize the coupled stack
        if (procTopoChanges)
        {
            initCoupledStack(entities, false);
        }

        edgeRefinementEngine(&(handlerPtr_[0]));
    }

    // Re-Initialize the stack
    if (procTopoChanges)
    {
        initCoupledStack(entities, false);
    }

    if (is2D())
    {
        swap2DEdges(&(handlerPtr_[0]));
    }
    else
    {
        swap3DEdges(&(handlerPtr_[0]));
    }

    // Build a list of entities that need to be avoided
    // by regular topo-changes.
    buildEntitiesToAvoid(entities, true);

    // Reset coupled modifications.
    unsetCoupledModification();

    if (debug)
    {
        Info<< "Done." << endl;

        // Check coupled-patch sizes after changes.
        forAll(patchCoupling_, patchI)
        {
            if (patchCoupling_(patchI))
            {
                const coupleMap& cMap = patchCoupling_[patchI].map();

                label mSize = patchSizes_[cMap.masterIndex()];
                label sSize = patchSizes_[cMap.slaveIndex()];

                if (mSize != sSize)
                {
                    Pout<< "Coupled patch-count is inconsistent." << nl
                        << " Master Patch: " << cMap.masterIndex()
                        << " Count: " << mSize << nl
                        << " Slave Patch: " << cMap.slaveIndex()
                        << " Count: " << sSize
                        << endl;

                    FatalErrorIn("dynamicTopoFvMesh::handleCoupledPatches()")
                        << " Failures were found in connectivity"
                        << " after coupled topo-changes."
                        << abort(FatalError);
                }
            }
        }
    }

    // Schedule transfer of topology operations across processors
    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        coupledMesh& sPM = sendMeshes_[pI];
        coupledMesh& rPM = recvMeshes_[pI];

        if (priority(proc, lessOp<label>(), Pstream::myProcNo()))
        {
            const coupleMap& cMap = sPM.map();

            // How many entities am I receiving..
            FixedList<label, 3> nOps(-1);

            meshOps::pRead(proc, nOps);

            label nEntities = nOps[0];
            label nMovePoints = nOps[1];
            label nPatchConvert = nOps[2];

            if (debug > 3)
            {
                Pout<< " Op transfer:"
                    << " Receiving from [" << proc << "]::"
                    << " nEntities: " << nEntities
                    << " nMovePoints: " << nMovePoints
                    << " nPatchConvert: " << nPatchConvert
                    << endl;
            }

            if (nEntities)
            {
                // Size up the receipt buffers
                labelList& indices = cMap.entityIndices();
                List<coupleMap::opType>& operations = cMap.entityOperations();

                indices.setSize(nEntities);
                operations.setSize(nEntities);

                // Schedule indices and operations for receipt
                meshOps::pRead(proc, indices);
                meshOps::pRead(proc, operations);
            }

            if (nMovePoints)
            {
                // Size up the receipt buffers
                pointField& newPoints = cMap.moveNewPoints();
                pointField& oldPoints = cMap.moveOldPoints();

                newPoints.setSize(nMovePoints);
                oldPoints.setSize(nMovePoints);

                // Schedule old / new points for receipt
                meshOps::pRead(proc, newPoints);
                meshOps::pRead(proc, oldPoints);
            }

            if (nPatchConvert)
            {
                // Size up the receipt buffers
                labelList& patches = cMap.patchIndices();

                patches.setSize(nPatchConvert);

                // Schedule patch indices for receipt
                meshOps::pRead(proc, patches);
            }
        }
        else
        {
            const coupleMap& cMap = rPM.map();

            FixedList<label, 3> nOps(-1);

            label nEntities = cMap.entityIndices().size();
            label nMovePoints = cMap.moveNewPoints().size();
            label nPatchConvert = cMap.patchIndices().size();

            if (debug > 3)
            {
                Pout<< " Op transfer:"
                    << " Sending to [" << proc << "]:: "
                    << " nEntities: " << nEntities << nl
                    << "  entityIndices: " << cMap.entityIndices() << nl
                    << "  entityOperations: " << cMap.entityOperations() << nl
                    << " nMovePoints: " << nMovePoints << nl
                    << "  moveNewPoints: " << cMap.moveNewPoints() << nl
                    << "  moveOldPoints: " << cMap.moveOldPoints() << nl
                    << " nPatchConvert: " << nPatchConvert << nl
                    << "  patchIndices: " << cMap.patchIndices() << nl
                    << endl;
            }

            nOps[0] = nEntities;
            nOps[1] = nMovePoints;
            nOps[2] = nPatchConvert;

            meshOps::pWrite(proc, nOps);

            if (nEntities)
            {
                // Schedule transfer to processor
                meshOps::pWrite(proc, cMap.entityIndices());
                meshOps::pWrite(proc, cMap.entityOperations());
            }

            if (nMovePoints)
            {
                // Schedule transfer to processor
                meshOps::pWrite(proc, cMap.moveNewPoints());
                meshOps::pWrite(proc, cMap.moveOldPoints());
            }

            if (nPatchConvert)
            {
                // Schedule transfer to processor
                meshOps::pWrite(proc, cMap.patchIndices());
            }
        }
    }

    // We won't wait for transfers to complete for the moment,
    // and will deal with operations once the internal mesh
    // has been dealt with.
}


// Synchronize topology operations across processors
void dynamicTopoFvMesh::syncCoupledPatches(labelHashSet& entities)
{
    if (!Pstream::parRun())
    {
        return;
    }

    // Temporarily reset maxModifications to
    // ensure that synchronization succeeds
    label maxModSave = maxModifications_;

    maxModifications_ = -1;

    const polyBoundaryMesh& boundary = boundaryMesh();

    labelList createPatchOrders(procIndices_.size(), 0);

    forAll(procIndices_, pI)
    {
        const coupleMap& cMap = recvMeshes_[pI].map();

        if (cMap.patchIndex() >= boundary.size())
        {
            // This processor needs a new patch talking to me
            createPatchOrders[pI] = 1;
        }
    }

    // Send / recv create patch orders, if any
    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        if (priority(proc, greaterOp<label>(), Pstream::myProcNo()))
        {
            // Non-blocking send to processor
            meshOps::pWrite(proc, createPatchOrders[pI], false);
        }
        else
        {
            // Non-blocking receive from processor
            meshOps::pRead(proc, createPatchOrders[pI], false);
        }
    }

    // Wait for all transfers to complete.
    meshOps::waitForBuffers();

    Map<label> addedProcPatches;

    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        if (createPatchOrders[pI])
        {
            // Create a new processor patch
            addedProcPatches.insert
            (
                proc,
                createProcessorPatch(proc)
            );
        }
    }

    // Buffer for cell-removal
    dynamicLabelList rCellList(10);

    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        const coupledMesh& sPM = sendMeshes_[pI];

        if (priority(proc, lessOp<label>(), Pstream::myProcNo()))
        {
            const coupleMap& cMap = sPM.map();
            const labelList& indices = cMap.entityIndices();
            const List<coupleMap::opType>& operations = cMap.entityOperations();

            const labelList& patches = cMap.patchIndices();
            const pointField& newPoints = cMap.moveNewPoints();
            const pointField& oldPoints = cMap.moveOldPoints();

            // Find the appropriate processor patch
            // for cell-removal operations
            label procPatch = -1;

            forAll(boundary, pI)
            {
                if (!isA<processorPolyPatch>(boundary[pI]))
                {
                    continue;
                }

                const processorPolyPatch& pp =
                (
                    refCast<const processorPolyPatch>(boundary[pI])
                );

                if (pp.neighbProcNo() == proc)
                {
                    procPatch = pI;
                    break;
                }
            }

            if (procPatch == -1)
            {
                // Check if this was a newly added patch
                if (addedProcPatches.found(proc))
                {
                    procPatch = addedProcPatches[proc];
                }
                else
                if (operations.size())
                {
                    bool required = false;

                    forAll(operations, opI)
                    {
                        const coupleMap::opType op = operations[opI];

                        if
                        (
                            op == coupleMap::REMOVE_CELL ||
                            op == coupleMap::CONVERT_PATCH
                        )
                        {
                            required = true;
                            break;
                        }
                    }

                    // If we actually have operations from
                    // this processor that require patch conversion,
                    // this is a problem.
                    if (required)
                    {
                        Pout<< " * * * Sync Operations * * * " << nl
                            << " Could not find patch for proc: " << proc << nl
                            << " procIndices: " << procIndices_
                            << abort(FatalError);
                    }
                }
                else
                {
                    // Move on to the next processor
                    continue;
                }
            }

            label pointCounter = 0, patchCounter = 0;

            // Sequentially execute operations
            forAll(indices, indexI)
            {
                const label index = indices[indexI];
                const coupleMap::opType op = operations[indexI];

                if (debug > 3)
                {
                    Pout<< " Recv Op index: " << index << endl;
                }

                // Determine the appropriate local index
                label localIndex = -1;

                if (op == coupleMap::MOVE_POINT)
                {
                    localIndex = cMap.entityMap(coupleMap::POINT)[index];
                }
                else
                if
                (
                    op == coupleMap::CONVERT_PATCH ||
                    op == coupleMap::CONVERT_PHYSICAL
                )
                {
                    localIndex =
                    (
                        cMap.entityMap(coupleMap::FACE).found(index) ?
                        cMap.entityMap(coupleMap::FACE)[index] : -1
                    );

                    if (debug > 3)
                    {
                        if (op == coupleMap::CONVERT_PATCH)
                        {
                            const point& newCentre = newPoints[pointCounter];
                            const point& oldCentre = oldPoints[pointCounter];

                            Pout<< " Convert patch: " << index << nl
                                << "  localIndex: " << localIndex << nl
                                << "  newCentre: " << newCentre << nl
                                << "  oldCentre: " << oldCentre << nl
                                << endl;
                        }

                        if (op == coupleMap::CONVERT_PHYSICAL)
                        {
                            Pout<< " Convert patch: " << index << nl
                                << "  localIndex: " << localIndex << nl
                                << "  patch: " << patches[patchCounter] << nl
                                << endl;
                        }
                    }
                }
                else
                {
                    // Pick localIndex based on entity type
                    if (is2D())
                    {
                        const Map<label>& fM = cMap.entityMap(coupleMap::FACE);

                        Map<label>::const_iterator it = fM.find(index);

                        if (it != fM.end())
                        {
                            localIndex = it();
                        }
                        else
                        {
                            Pout<< " * * * Sync Operations * * * " << nl
                                << " Could not find index: " << index
                                << " in faceMap for proc: " << proc << nl
                                << " nSentFaces: "
                                << cMap.nEntities(coupleMap::FACE)
                                << " opType: " << coupleMap::asText(op)
                                << abort(FatalError);
                        }
                    }
                    else
                    {
                        const Map<label>& eM = cMap.entityMap(coupleMap::EDGE);

                        Map<label>::const_iterator it = eM.find(index);

                        if (it != eM.end())
                        {
                            localIndex = it();
                        }
                        else
                        {
                            Pout<< " * * * Sync Operations * * * " << nl
                                << " Could not find index: " << index
                                << " in edgeMap for proc: " << proc << nl
                                << " nSentEdges: "
                                << cMap.nEntities(coupleMap::EDGE)
                                << " opType: " << coupleMap::asText(op)
                                << abort(FatalError);
                        }
                    }
                }

                changeMap opMap;

                switch (op)
                {
                    case coupleMap::BISECTION:
                    {
                        opMap = bisectEdge(localIndex);
                        break;
                    }

                    case coupleMap::COLLAPSE_FIRST:
                    {
                        opMap = collapseEdge(localIndex, 1);
                        break;
                    }

                    case coupleMap::COLLAPSE_SECOND:
                    {
                        opMap = collapseEdge(localIndex, 2);
                        break;
                    }

                    case coupleMap::COLLAPSE_MIDPOINT:
                    {
                        opMap = collapseEdge(localIndex, 3);
                        break;
                    }

                    case coupleMap::REMOVE_CELL:
                    {
                        // Clear existing list
                        rCellList.clear();

                        if (is2D())
                        {
                            // Insert the owner cell
                            rCellList.append(owner_[localIndex]);
                        }
                        else
                        {
                            // Insert all cells connected to this edge
                            const labelList& eFaces = edgeFaces_[localIndex];

                            forAll(eFaces, faceI)
                            {
                                const label own = owner_[eFaces[faceI]];
                                const label nei = neighbour_[eFaces[faceI]];

                                if (findIndex(rCellList, own) == -1)
                                {
                                    rCellList.append(own);
                                }

                                if (nei != -1)
                                {
                                    if (findIndex(rCellList, nei) == -1)
                                    {
                                        rCellList.append(nei);
                                    }
                                }
                            }
                        }

                        opMap =
                        (
                            removeCells
                            (
                                rCellList,
                                procPatch,
                                "rcs_" + Foam::name(localIndex)
                            )
                        );

                        break;
                    }

                    case coupleMap::MOVE_POINT:
                    {
                        const point& newPoint = newPoints[pointCounter];
                        const point& oldPoint = oldPoints[pointCounter];

                        // Move old / new points
                        points_[localIndex] = newPoint;
                        oldPoints_[localIndex] = oldPoint;

                        // Clear the existing map
                        opMap.clear();

                        // Force a successful operation
                        opMap.type() = 1;

                        pointCounter++;

                        break;
                    }

                    case coupleMap::CONVERT_PATCH:
                    {
                        const point& newCentre = newPoints[pointCounter];
                        const point& oldCentre = oldPoints[pointCounter];

                        if (localIndex == -1)
                        {
                            // Accumulate stats in case of failure
                            scalar minDist = GREAT;
                            vector minPoint = vector::zero;
                            dynamicLabelList checkedFaces;

                            if (debug > 2)
                            {
                                // Reserve for append
                                checkedFaces.setCapacity(50);
                            }

                            // New face. Check all boundary faces
                            // and match up centre
                            label sTot = faces_.size();
                            label sInt = nOldInternalFaces_;

                            for (label faceI = sInt; faceI < sTot; faceI++)
                            {
                                const face& fCheck = faces_[faceI];

                                if (fCheck.empty())
                                {
                                    continue;
                                }

                                label pIndex = whichPatch(faceI);

                                if (getNeighbourProcessor(pIndex) == -1)
                                {
                                    continue;
                                }

                                // Compute face-centre
                                vector fC = fCheck.centre(points_);

                                // Compute tolerance
                                scalar tol = mag(points_[fCheck[0]] - fC);
                                scalar dist = mag(fC - newCentre);

                                if (dist < (geomMatchTol_()*tol))
                                {
                                    localIndex = faceI;
                                    break;
                                }
                                else
                                if (dist < minDist)
                                {
                                    minPoint = fC;
                                    minDist = dist;

                                    if (debug > 2)
                                    {
                                        checkedFaces.append(faceI);
                                    }
                                }
                            }

                            // Ensure that the face was found
                            if (localIndex == -1)
                            {
                                writeVTK
                                (
                                    "checkedFaces_"
                                  + Foam::name(index),
                                    checkedFaces,
                                    2, false, true
                                );

                                Pout<< " * * * syncCoupledPatches * * * " << nl
                                    << " Convert patch Op failed." << nl
                                    << " Face: " << index << nl
                                    << " minPoint: " << minPoint << nl
                                    << " minDistance: " << minDist << nl
                                    << " newCentre: " << newCentre << nl
                                    << " oldCentre: " << oldCentre << nl
                                    << abort(FatalError);
                            }
                        }

                        // Fetch reference to face
                        const face& fCheck = faces_[localIndex];

                        point fC = fCheck.centre(points_);

                        // Specify a tolerance
                        scalar tol = mag(points_[fCheck[0]] - fC);
                        scalar dist = mag(fC - newCentre);

                        // Ensure a face-match
                        if (dist > (geomMatchTol_() * tol))
                        {
                            Pout<< " * * * Sync Operations * * * " << nl
                                << " Convert patch Op failed." << nl
                                << " Index: " << index << nl
                                << " localIndex: " << localIndex << nl
                                << " face: " << fCheck << nl
                                << " faceCentre: " << fC << nl
                                << " Master processor: " << proc << nl
                                << " procPatch: " << procPatch << nl
                                << " tolerance: " << (geomMatchTol_() * tol) << nl
                                << " distance: " << dist << nl
                                << " pointCounter: " << pointCounter << nl
                                << " newCentre: " << newCentre << nl
                                << " oldCentre: " << oldCentre << nl
                                << abort(FatalError);
                        }

                        // Obtain a copy before adding the new face,
                        // since the reference might become
                        // invalid during list resizing.
                        // Edges don't have to change, since they're
                        // all on the boundary anyway.
                        face newFace = faces_[localIndex];
                        label newOwn = owner_[localIndex];
                        labelList newFaceEdges = faceEdges_[localIndex];

                        label newFaceIndex =
                        (
                            insertFace
                            (
                                procPatch,
                                newFace,
                                newOwn,
                                -1,
                                newFaceEdges
                            )
                        );

                        meshOps::replaceLabel
                        (
                            localIndex,
                            newFaceIndex,
                            cells_[newOwn]
                        );

                        // Correct edgeFaces with the new face label.
                        forAll(newFaceEdges, edgeI)
                        {
                            meshOps::replaceLabel
                            (
                                localIndex,
                                newFaceIndex,
                                edgeFaces_[newFaceEdges[edgeI]]
                            );
                        }

                        // Finally remove the face
                        removeFace(localIndex);

                        // Clear the existing map
                        opMap.clear();

                        // Update map
                        opMap.addFace
                        (
                            newFaceIndex,
                            labelList(1, localIndex)
                        );

                        // Force a successful operation
                        opMap.type() = 1;

                        pointCounter++;

                        break;
                    }

                    case coupleMap::CONVERT_PHYSICAL:
                    {
                        const label patchIndex = patches[patchCounter];

                        // Obtain a copy before adding the new face,
                        // since the reference might become
                        // invalid during list resizing.
                        // Edges don't have to change, since they're
                        // all on the boundary anyway.
                        face newFace = faces_[localIndex];
                        label newOwn = owner_[localIndex];
                        labelList newFaceEdges = faceEdges_[localIndex];

                        label newFaceIndex =
                        (
                            insertFace
                            (
                                patchIndex,
                                newFace,
                                newOwn,
                                -1,
                                newFaceEdges
                            )
                        );

                        meshOps::replaceLabel
                        (
                            localIndex,
                            newFaceIndex,
                            cells_[newOwn]
                        );

                        // Correct edgeFaces with the new face label.
                        forAll(newFaceEdges, edgeI)
                        {
                            meshOps::replaceLabel
                            (
                                localIndex,
                                newFaceIndex,
                                edgeFaces_[newFaceEdges[edgeI]]
                            );
                        }

                        // Finally remove the face
                        removeFace(localIndex);

                        // Clear the existing map
                        opMap.clear();

                        // Update map
                        opMap.addFace
                        (
                            newFaceIndex,
                            labelList(1, localIndex)
                        );

                        // Force a successful operation
                        opMap.type() = 1;

                        patchCounter++;

                        break;
                    }

                    case coupleMap::INVALID:
                    {
                        Pout<< " * * * Sync Operations * * * " << nl
                            << " Invalid operation." << nl
                            << " Index: " << index << nl
                            << " localIndex: " << localIndex << nl
                            << " operation: " << coupleMap::asText(op) << nl
                            << abort(FatalError);

                        break;
                    }
                }

                if (opMap.type() <= 0)
                {
                    Pout<< " * * * Sync Operations * * * " << nl
                        << " Operation failed." << nl
                        << " Index: " << index << nl
                        << " localIndex: " << localIndex << nl
                        << " operation: " << coupleMap::asText(op) << nl
                        << " opMap.type: " << opMap.type() << nl
                        << abort(FatalError);
                }
            }
        }
    }

    // Reinstate max modifications
    maxModifications_ = maxModSave;

    // Re-Initialize the stack with avoided entities from subMeshes
    // and leave those on processor patches untouched
    labelHashSet procEntities;

    buildEntitiesToAvoid(procEntities, false);

    // First remove processor entries
    forAllConstIter(labelHashSet, procEntities, pIter)
    {
        if (entities.found(pIter.key()))
        {
            entities.erase(pIter.key());
        }
    }

    // Initialize the coupled stack, using supplied entities
    initCoupledStack(entities, true);

    // Loop through the coupled stack and perform changes.
    if (edgeRefinement_)
    {
        edgeRefinementEngine(&(handlerPtr_[0]));
    }

    // Re-Initialize the stack, using supplied entities
    initCoupledStack(entities, true);

    if (is2D())
    {
        swap2DEdges(&(handlerPtr_[0]));
    }
    else
    {
        swap3DEdges(&(handlerPtr_[0]));
    }
}


// Check the state of coupled boundaries
//  - Return true if errors are present
bool dynamicTopoFvMesh::checkCoupledBoundaries(bool report) const
{
    const polyBoundaryMesh& boundary = boundaryMesh();

    bool sizeError = false, misMatchError = false;

    // Maintain a list of master / neighbour anchors
    List<vectorField> mAnchors(boundary.size());
    List<vectorField> nAnchors(boundary.size());

    forAll(boundary, pI)
    {
        if (isA<processorPolyPatch>(boundary[pI]))
        {
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pI])
            );

            label neiProcNo = pp.neighbProcNo();

            // Send point / face sizes
            meshOps::pWrite(neiProcNo, boundary[pI].nPoints());
            meshOps::pWrite(neiProcNo, boundary[pI].size());

            // Send points / centres / areas
            meshOps::pWrite(neiProcNo, boundary[pI].faceAreas());
            meshOps::pWrite(neiProcNo, boundary[pI].faceCentres());
            meshOps::pWrite(neiProcNo, boundary[pI].localPoints());

            // Prepare and send anchor points
            mAnchors[pI].setSize(boundary[pI].size());

            const faceList& lFaces = boundary[pI].localFaces();
            const pointField& lPoints = boundary[pI].localPoints();

            forAll(lFaces, faceI)
            {
                mAnchors[pI][faceI] = lPoints[lFaces[faceI][0]];
            }

            meshOps::pWrite(neiProcNo, mAnchors[pI]);
        }

        if (isA<cyclicPolyPatch>(boundary[pI]))
        {
            const cyclicPolyPatch& cp =
            (
                refCast<const cyclicPolyPatch>(boundary[pI])
            );

            if (cp.size() & 1)
            {
                sizeError = true;

                Pout<< " Incorrect size for cyclic patch: "
                    << cp.size() << endl;
            }

            if (!sizeError)
            {
                // Compute halfSize
                label halfSize = cp.size() / 2;

                vectorField half0Areas
                (
                    SubField<vector>
                    (
                        cp.faceAreas(),
                        halfSize
                    )
                );

                vectorField half0Centres
                (
                    SubField<vector>
                    (
                        cp.faceCentres(),
                        halfSize
                    )
                );

                vectorField half1Areas
                (
                    SubField<vector>
                    (
                        cp.faceAreas(),
                        halfSize,
                        halfSize
                    )
                );

                vectorField half1Centres
                (
                    SubField<vector>
                    (
                        cp.faceCentres(),
                        halfSize,
                        halfSize
                    )
                );

                const faceList& lF = boundary[pI].localFaces();
                const pointField& lP = boundary[pI].localPoints();
                const vectorField& lC = boundary[pI].faceCentres();

                // Prepare anchor points
                mAnchors[pI].setSize(boundary[pI].size());

                forAll(lF, faceI)
                {
                    mAnchors[pI][faceI] = lP[lF[faceI][0]];
                }

                // Transform first-half points and check
                if (cp.transform() == cyclicPolyPatch::ROTATIONAL)
                {
                    forAll(half0Centres, faceI)
                    {
                        half0Centres[faceI] =
                        (
                            Foam::transform
                            (
                                cp.forwardT()[0],
                                half0Centres[faceI]
                            )
                        );

                        mAnchors[pI][faceI] =
                        (
                            Foam::transform
                            (
                                cp.forwardT()[0],
                                mAnchors[pI][faceI]
                            )
                        );
                    }
                }
                else
                if (cp.transform() == cyclicPolyPatch::TRANSLATIONAL)
                {
                    forAll(half0Centres, faceI)
                    {
                        half0Centres[faceI] += cp.separationVector();
                        mAnchors[pI][faceI] += cp.separationVector();
                    }
                }
                else
                {
                    Pout<< " Cyclic check: Unknown transform."
                        << abort(FatalError);
                }

                // Calculate a point-match tolerance per patch
                scalar pTol = -GREAT;

                // Check areas / compute tolerance
                forAll(half0Areas, faceI)
                {
                    scalar fMagSf = mag(half0Areas[faceI]);
                    scalar rMagSf = mag(half1Areas[faceI]);
                    scalar avSf = 0.5 * (fMagSf + rMagSf);

                    if (mag(fMagSf - rMagSf)/avSf > geomMatchTol_())
                    {
                        misMatchError = true;

                        Pout<< " Face: " << faceI
                            << " area does not match neighbour by: "
                            << 100 * mag(fMagSf - rMagSf)/avSf
                            << "% - possible patch ordering problem. "
                            << " Front area:" << fMagSf
                            << " Rear area: " << rMagSf
                            << endl;
                    }

                    pTol =
                    (
                        Foam::max
                        (
                            pTol,
                            geomMatchTol_() * mag(lP[lF[faceI][0]] - lC[faceI])
                        )
                    );
                }

                // Check centres / anchor points
                forAll(half0Centres, faceI)
                {
                    scalar distA =
                    (
                        mag
                        (
                            mAnchors[pI][faceI]
                          - mAnchors[pI][faceI + halfSize]
                        )
                    );

                    scalar distC =
                    (
                        mag
                        (
                            half0Centres[faceI]
                          - half1Centres[faceI]
                        )
                    );

                    if (distA > pTol || distC > pTol)
                    {
                        misMatchError = true;

                        UIndirectList<point> f1(lP, lF[faceI]);
                        UIndirectList<point> f2(lP, lF[faceI + halfSize]);

                        Pout<< " Face: " << faceI << nl
                            << " Points: " << nl << f1 << nl << f2 << nl
                            << " Anchors ::" << nl
                            << mAnchors[pI][faceI] << nl
                            << mAnchors[pI][faceI + halfSize] << nl
                            << " Centres ::" << nl
                            << half0Centres[faceI] << nl
                            << half1Centres[faceI] << nl
                            << " Tolerance: " << pTol << nl
                            << " Measured Anchor distance: " << distA << nl
                            << " Measured Centre distance: " << distC << nl
                            << endl;
                    }
                }

                if (misMatchError)
                {
                    // Write out to disk
                    meshOps::writeVTK
                    (
                        (*this),
                        cp.name(),
                        identity(cp.size()),
                        2,
                        cp.localPoints(),
                        List<edge>(0),
                        cp.localFaces()
                    );
                }
            }
        }
    }

    // Maintain a list of points / areas / centres
    List<vectorField> fAreas(boundary.size());
    List<vectorField> fCentres(boundary.size());
    List<vectorField> neiPoints(boundary.size());

    forAll(boundary, pI)
    {
        if (isA<processorPolyPatch>(boundary[pI]))
        {
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[pI])
            );

            label neiProcNo = pp.neighbProcNo();

            // Receive point / face sizes
            label neiPSize = -1, neiSize = -1;

            meshOps::pRead(neiProcNo, neiPSize);
            meshOps::pRead(neiProcNo, neiSize);

            // Size up lists
            fAreas[pI].setSize(neiSize);
            fCentres[pI].setSize(neiSize);
            neiPoints[pI].setSize(neiPSize);
            nAnchors[pI].setSize(neiSize);

            // Receive points / centres / areas / anchors
            meshOps::pRead(neiProcNo, fAreas[pI]);
            meshOps::pRead(neiProcNo, fCentres[pI]);
            meshOps::pRead(neiProcNo, neiPoints[pI]);
            meshOps::pRead(neiProcNo, nAnchors[pI]);

            if
            (
                neiSize != boundary[pI].size() ||
                neiPSize != boundary[pI].nPoints()
            )
            {
                sizeError = true;

                Pout<< "Incorrect send / recv sizes: " << nl
                    << " Patch: " << boundary[pI].name() << nl
                    << " My nPoints: " << boundary[pI].nPoints() << nl
                    << " Nei nPoints: " << neiPSize << nl
                    << " My nFaces: " << boundary[pI].size() << nl
                    << " Nei nFaces: " << neiSize << nl
                    << endl;
            }
        }
    }

    // Wait for transfers to complete
    meshOps::waitForBuffers();

    // Check centres / areas
    forAll(boundary, pI)
    {
        if (!isA<processorPolyPatch>(boundary[pI]))
        {
            continue;
        }

        const vectorField& myAreas = boundary[pI].faceAreas();
        const vectorField& myCentres = boundary[pI].faceCentres();

        // Fetch local connectivity
        const faceList& myFaces = boundary[pI].localFaces();
        const vectorField& myPoints = boundary[pI].localPoints();

        // Calculate a point-match tolerance per patch
        scalar pTol = -GREAT;

        forAll(myAreas, faceI)
        {
            scalar magSf = mag(myAreas[faceI]);
            scalar nbrMagSf = mag(fAreas[pI][faceI]);
            scalar avSf = 0.5 * (magSf + nbrMagSf);

            if (mag(magSf - nbrMagSf)/avSf > geomMatchTol_())
            {
                misMatchError = true;

                Pout<< " Face: " << faceI
                    << " area does not match neighbour by: "
                    << 100 * mag(magSf - nbrMagSf)/avSf
                    << "% - possible patch ordering problem. "
                    << " My area:" << magSf
                    << " Neighbour area: " << nbrMagSf
                    << " My centre: " << myCentres[faceI]
                    << endl;
            }

            pTol =
            (
                Foam::max
                (
                    pTol,
                    geomMatchTol_() *
                    mag
                    (
                        myPoints[myFaces[faceI][0]]
                      - myCentres[faceI]
                    )
                )
            );
        }

        const processorPolyPatch& pp =
        (
            refCast<const processorPolyPatch>(boundary[pI])
        );

        // Check anchor points
        forAll(mAnchors[pI], faceI)
        {
            scalar dist = mag(mAnchors[pI][faceI] - nAnchors[pI][faceI]);

            if (dist > pTol)
            {
                misMatchError = true;

                Pout<< " Face: " << faceI << "::" << mAnchors[pI][faceI] << nl
                    << " Mismatch for processor: " << pp.neighbProcNo() << nl
                    << " Neighbour Anchor: " << nAnchors[pI][faceI] << nl
                    << " Tolerance: " << pTol << nl
                    << " Measured distance: " << dist << nl
                    << endl;
            }
        }

        // Check neighbPoints addressing
        const labelList& nPts = pp.neighbPoints();

        forAll(myPoints, pointI)
        {
            if (nPts[pointI] < 0)
            {
                Pout<< " Point: " <<  pointI << "::" << myPoints[pointI]
                    << " reported to be multiply connected." << nl
                    << " neighbPoint: " << nPts[pointI] << nl
                    << endl;

                misMatchError = true;

                continue;
            }

            scalar dist = mag(myPoints[pointI] - neiPoints[pI][nPts[pointI]]);

            if (dist > pTol)
            {
                misMatchError = true;

                Pout<< " Point: " << pointI << "::" << myPoints[pointI] << nl
                    << " Mismatch for processor: " << pp.neighbProcNo() << nl
                    << " Neighbour Pt: " << neiPoints[pI][nPts[pointI]] << nl
                    << " Tolerance: " << pTol << nl
                    << " Measured distance: " << dist << nl
                    << " Neighbour Index:" << nPts[pointI] << nl
                    << endl;
            }
        }
    }

    if (!sizeError && !misMatchError && report)
    {
        Pout<< " Coupled boundary check: No problems." << endl;
    }

    return (sizeError || misMatchError);
}


// Build patch sub-meshes for processor patches
void dynamicTopoFvMesh::buildProcessorPatchMeshes()
{
    if (procIndices_.empty())
    {
        return;
    }

    // Maintain a list of cells common to multiple processors.
    Map<labelList> commonCells;

    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        coupledMesh& sPM = sendMeshes_[pI];

        // Build the subMesh.
        buildProcessorPatchMesh(sPM, commonCells);

        const coupleMap& scMap = sPM.map();

        // Send my sub-mesh to the neighbour.
        meshOps::pWrite(proc, scMap.nEntities());

        if (debug > 3)
        {
            Pout<< "Sending to [" << proc << "]:: nEntities: "
                << scMap.nEntities() << endl;
        }

        // Send the pointBuffers
        meshOps::pWrite(proc, scMap.pointBuffer());
        meshOps::pWrite(proc, scMap.oldPointBuffer());

        // Send connectivity (points, edges, faces, cells, etc)
        forAll(scMap.entityBuffer(), bufferI)
        {
            if (scMap.entityBuffer(bufferI).size())
            {
                meshOps::pWrite(proc, scMap.entityBuffer(bufferI));
            }
        }

        // Obtain references
        const coupledMesh& rPM = recvMeshes_[pI];
        const coupleMap& rcMap = rPM.map();

        // First read entity sizes.
        meshOps::pRead(proc, rcMap.nEntities());

        if (debug > 3)
        {
            Pout<< "Receiving from [" << proc << "]:: nEntities: "
                << rcMap.nEntities() << endl;
        }

        // Size the buffers.
        rcMap.allocateBuffers();

        // Receive the pointBuffers
        meshOps::pRead(proc, rcMap.pointBuffer());
        meshOps::pRead(proc, rcMap.oldPointBuffer());

        // Receive connectivity (points, edges, faces, cells, etc)
        forAll(rcMap.entityBuffer(), bufferI)
        {
            if (rcMap.entityBuffer(bufferI).size())
            {
                meshOps::pRead(proc, rcMap.entityBuffer(bufferI));
            }
        }
    }

    // We won't wait for all transfers to complete for the moment.
    // Meanwhile, do some other useful work, if possible.
}


// Build patch sub-mesh for a specified processor
// - At this point, procIndices is available as a sorted list
//   of neighbouring processors.
void dynamicTopoFvMesh::buildProcessorPatchMesh
(
    coupledMesh& subMesh,
    Map<labelList>& commonCells
)
{
    label nP = 0, nE = 0, nF = 0, nC = 0;

    // Obtain references
    const coupleMap& cMap = subMesh.map();
    const labelList& subMeshPoints = cMap.subMeshPoints();
    const List<labelPair>& globalProcPoints = cMap.globalProcPoints();

    Map<label>& rPointMap = cMap.reverseEntityMap(coupleMap::POINT);
    Map<label>& rEdgeMap = cMap.reverseEntityMap(coupleMap::EDGE);
    Map<label>& rFaceMap = cMap.reverseEntityMap(coupleMap::FACE);
    Map<label>& rCellMap = cMap.reverseEntityMap(coupleMap::CELL);

    Map<label>& pointMap = cMap.entityMap(coupleMap::POINT);
    Map<label>& edgeMap = cMap.entityMap(coupleMap::EDGE);
    Map<label>& faceMap = cMap.entityMap(coupleMap::FACE);
    Map<label>& cellMap = cMap.entityMap(coupleMap::CELL);

    // Add all cells connected to points on the subMeshPoints list
    label proc = -1;

    if (cMap.masterIndex() == Pstream::myProcNo())
    {
        proc = cMap.slaveIndex();
    }
    else
    {
        proc = cMap.masterIndex();
    }

    // Check if this is a direct neighbour
    const polyBoundaryMesh& boundary = boundaryMesh();

    // Add sub-mesh points first.
    // Additional halo points will be added later.
    forAll(subMeshPoints, pointI)
    {
        pointMap.insert(nP, subMeshPoints[pointI]);
        rPointMap.insert(subMeshPoints[pointI], nP);
        nP++;
    }

    // Set the number of points (shared) at this point.
    cMap.nEntities(coupleMap::SHARED_POINT) = nP;

    // Size up the point-buffer with global index information.
    // Direct neighbours do not require any addressing.
    labelList& gpBuffer = cMap.entityBuffer(coupleMap::POINT);
    gpBuffer.setSize(globalProcPoints.size(), -1);

    forAll(globalProcPoints, pointI)
    {
        pointMap.insert(nP, globalProcPoints[pointI].first());
        rPointMap.insert(globalProcPoints[pointI].first(), nP);
        nP++;

        // Fill in buffer with global point index
        gpBuffer[pointI] = globalProcPoints[pointI].second();
    }

    // Set the number of points (shared + global) at this point.
    cMap.nEntities(coupleMap::GLOBAL_POINT) = nP;

    labelHashSet localCommonCells;

    // Detect all cells surrounding shared points.
    if (is2D())
    {
        // No pointEdges structure, so loop through all cells
        // to check for those connected to subMeshPoints
        forAll(cells_, cellI)
        {
            const cell& cellToCheck = cells_[cellI];

            if (cellToCheck.empty())
            {
                continue;
            }

            const labelList cellPoints = cellToCheck.labels(faces_);

            forAll(cellPoints, pointI)
            {
                if (rPointMap.found(cellPoints[pointI]))
                {
                    if (commonCells.found(cellI))
                    {
                        // Add locally common cells at the end.
                        if (!localCommonCells.found(cellI))
                        {
                            localCommonCells.insert(cellI);
                        }

                        // Check if the processor exists on the list
                        // and if not, add it.
                        labelList& procList = commonCells[cellI];

                        if (findIndex(procList, proc) == -1)
                        {
                            meshOps::sizeUpList(proc, procList);
                        }
                    }
                    else
                    {
                        cellMap.insert(nC, cellI);
                        rCellMap.insert(cellI, nC);
                        nC++;

                        commonCells.insert(cellI, labelList(1, proc));
                    }

                    break;
                }
            }
        }
    }
    else
    {
        forAllConstIter(Map<label>, rPointMap, pIter)
        {
            // Loop through pointEdges for this point.
            const labelList& pEdges = pointEdges_[pIter.key()];

            forAll(pEdges, edgeI)
            {
                const labelList& eFaces = edgeFaces_[pEdges[edgeI]];

                forAll(eFaces, faceI)
                {
                    label own = owner_[eFaces[faceI]];
                    label nei = neighbour_[eFaces[faceI]];

                    // Check owner cell
                    if (!rCellMap.found(own))
                    {
                        if (commonCells.found(own))
                        {
                            // Add locally common cells at the end.
                            if (!localCommonCells.found(own))
                            {
                                localCommonCells.insert(own);
                            }

                            // Check if the processor exists on the list
                            // and if not, add it.
                            labelList& procList = commonCells[own];

                            if (findIndex(procList, proc) == -1)
                            {
                                meshOps::sizeUpList(proc, procList);
                            }
                        }
                        else
                        {
                            cellMap.insert(nC, own);
                            rCellMap.insert(own, nC);
                            nC++;

                            commonCells.insert(own, labelList(1, proc));
                        }
                    }

                    // Check neighbour cell
                    if (!rCellMap.found(nei) && nei != -1)
                    {
                        if (commonCells.found(nei))
                        {
                            // Add locally common cells at the end.
                            if (!localCommonCells.found(nei))
                            {
                                localCommonCells.insert(nei);
                            }

                            // Check if the processor exists on the list
                            // and if not, add it.
                            labelList& procList = commonCells[nei];

                            if (findIndex(procList, proc) == -1)
                            {
                                meshOps::sizeUpList(proc, procList);
                            }
                        }
                        else
                        {
                            cellMap.insert(nC, nei);
                            rCellMap.insert(nei, nC);
                            nC++;

                            commonCells.insert(nei, labelList(1, proc));
                        }
                    }
                }
            }
        }
    }

    // Now add locally common cells.
    forAllConstIter(labelHashSet, localCommonCells, cIter)
    {
        cellMap.insert(nC, cIter.key());
        rCellMap.insert(cIter.key(), nC);
        nC++;
    }

    // Keep track of inserted boundary face indices
    // - Add exposed internal faces to a 'default' patch
    //   at the end of the list.
    labelList bdyFaceSizes(boundary.size() + 1, 0);
    labelList bdyFaceStarts(boundary.size() + 1, 0);
    List<Map<label> > bdyFaceIndices(boundary.size() + 1);

    // Allocate the faceMap. Since the number of internal faces is unknown,
    // detect internal ones first and update boundaries later.
    label sumNFE = 0;

    forAllConstIter(Map<label>, rCellMap, cIter)
    {
        const cell& thisCell = cells_[cIter.key()];

        forAll(thisCell, faceI)
        {
            label fIndex = thisCell[faceI];

            if (!rFaceMap.found(fIndex))
            {
                // Determine the patch index
                label patchID = whichPatch(fIndex);

                if (patchID == -1)
                {
                    // Internal face. Check if this needs to
                    // be added to the 'default' patch.
                    label own = owner_[fIndex];
                    label nei = neighbour_[fIndex];

                    if (rCellMap.found(own) && rCellMap.found(nei))
                    {
                        faceMap.insert(nF, fIndex);
                        rFaceMap.insert(fIndex, nF);
                        nF++;
                    }
                    else
                    {
                        // Update boundary maps.
                        label bfI = bdyFaceSizes[boundary.size()]++;

                        // Skip faceMap and update the reverseMap for now.
                        // faceMap will be updated once all
                        // internal faces have been detected.
                        rFaceMap.insert(fIndex, bfI);
                        bdyFaceIndices[boundary.size()].insert(bfI, fIndex);
                    }
                }
                else
                {
                    // Update boundary maps.
                    label bfI = bdyFaceSizes[patchID]++;

                    // Skip faceMap and update the reverseMap for now.
                    // faceMap will be updated once all
                    // internal faces have been detected.
                    rFaceMap.insert(fIndex, bfI);
                    bdyFaceIndices[patchID].insert(bfI, fIndex);
                }

                // Accumulate face sizes
                sumNFE += faces_[fIndex].size();
            }
        }
    }

    // Set the number of internal faces at this point
    cMap.nEntities(coupleMap::INTERNAL_FACE) = nF;

    // Set patch starts
    bdyFaceStarts[0] = nF;

    for (label i = 1; i < bdyFaceStarts.size(); i++)
    {
        bdyFaceStarts[i] = bdyFaceStarts[i-1] + bdyFaceSizes[i-1];
    }

    // Update faceMap and reverseFaceMap for boundaries
    forAll(bdyFaceIndices, patchI)
    {
        label pStart = bdyFaceStarts[patchI];

        forAllConstIter(Map<label>, bdyFaceIndices[patchI], fIter)
        {
            faceMap.insert(fIter.key() + pStart, fIter());
            rFaceMap[fIter()] = fIter.key() + pStart;
        }

        // Update face-count
        nF += bdyFaceSizes[patchI];
    }

    // Keep track of inserted boundary edge indices
    // - Add exposed internal edges to a 'default' patch
    //   at the end of the list.
    labelList bdyEdgeSizes(boundary.size() + 1, 0);
    labelList bdyEdgeStarts(boundary.size() + 1, 0);
    List<Map<label> > bdyEdgeIndices(boundary.size() + 1);

    // Allocate the edgeMap. Since the number of internal edges is unknown,
    // detect internal ones first and update boundaries later.
    forAllConstIter(Map<label>, rFaceMap, fIter)
    {
        const labelList& fEdges = faceEdges_[fIter.key()];

        forAll(fEdges, edgeI)
        {
            label eIndex = fEdges[edgeI];

            if (!rEdgeMap.found(eIndex))
            {
                // Determine the patch index
                label patchID = whichEdgePatch(eIndex);

                if (patchID == -1)
                {
                    bool boundaryEdge = false;

                    // Check if any cells touching edgeFaces
                    // do not belong to the cellMap.
                    const labelList& eFaces = edgeFaces_[eIndex];

                    forAll(eFaces, faceI)
                    {
                        label fIndex = eFaces[faceI];

                        label own = owner_[fIndex];
                        label nei = neighbour_[fIndex];

                        if (!rCellMap.found(own) || !rCellMap.found(nei))
                        {
                            boundaryEdge = true;
                            break;
                        }
                    }

                    if (boundaryEdge)
                    {
                        // Update boundary maps.
                        label beI = bdyEdgeSizes[boundary.size()]++;

                        // Skip edgeMap and update the reverseMap for now.
                        // edgeMap will be updated once all
                        // internal edges have been detected.
                        rEdgeMap.insert(eIndex, beI);
                        bdyEdgeIndices[boundary.size()].insert(beI, eIndex);
                    }
                    else
                    {
                        // Internal edge
                        edgeMap.insert(nE, eIndex);
                        rEdgeMap.insert(eIndex, nE);
                        nE++;
                    }
                }
                else
                {
                    // Update boundary maps.
                    label beI = bdyEdgeSizes[patchID]++;

                    // Skip edgeMap and update the reverseMap for now.
                    // edgeMap will be updated once all
                    // internal edges have been detected.
                    rEdgeMap.insert(eIndex, beI);
                    bdyEdgeIndices[patchID].insert(beI, eIndex);
                }
            }
        }
    }

    // Set the number of internal edges at this point
    cMap.nEntities(coupleMap::INTERNAL_EDGE) = nE;

    // Set patch starts
    bdyEdgeStarts[0] = nE;

    for (label i = 1; i < bdyEdgeStarts.size(); i++)
    {
        bdyEdgeStarts[i] = bdyEdgeStarts[i-1] + bdyEdgeSizes[i-1];
    }

    // Update edgeMap and reverseEdgeMap for boundaries
    forAll(bdyEdgeIndices, patchI)
    {
        label pStart = bdyEdgeStarts[patchI];

        forAllConstIter(Map<label>, bdyEdgeIndices[patchI], eIter)
        {
            edgeMap.insert(eIter.key() + pStart, eIter());
            rEdgeMap[eIter()] = eIter.key() + pStart;
        }

        // Update edge-count
        nE += bdyEdgeSizes[patchI];
    }

    // Set additional points in the pointMap
    forAllConstIter(Map<label>, rEdgeMap, eIter)
    {
        const edge& thisEdge = edges_[eIter.key()];

        if (!rPointMap.found(thisEdge[0]))
        {
            pointMap.insert(nP, thisEdge[0]);
            rPointMap.insert(thisEdge[0], nP);
            nP++;
        }

        if (!rPointMap.found(thisEdge[1]))
        {
            pointMap.insert(nP, thisEdge[1]);
            rPointMap.insert(thisEdge[1], nP);
            nP++;
        }
    }

    // Assign sizes to the mesh
    cMap.nEntities(coupleMap::POINT) = nP;
    cMap.nEntities(coupleMap::EDGE) = nE;
    cMap.nEntities(coupleMap::FACE) = nF;
    cMap.nEntities(coupleMap::CELL) = nC;
    cMap.nEntities(coupleMap::NFE_SIZE) = sumNFE;
    cMap.nEntities(coupleMap::NBDY) = boundary.size() + 1;

    // Size up buffers and fill them
    cMap.allocateBuffers();

    pointField& pBuffer = cMap.pointBuffer();
    pointField& opBuffer = cMap.oldPointBuffer();

    forAllConstIter(Map<label>, pointMap, pIter)
    {
        pBuffer[pIter.key()] = points_[pIter()];
        opBuffer[pIter.key()] = oldPoints_[pIter()];
    }

    label index = 0;

    // Edge buffer size: 2 points for every edge
    labelList& eBuffer = cMap.entityBuffer(coupleMap::EDGE);

    for (label i = 0; i < nE; i++)
    {
        label eIndex = edgeMap[i];
        const edge& edgeToCheck = edges_[eIndex];

        eBuffer[index++] = rPointMap[edgeToCheck[0]];
        eBuffer[index++] = rPointMap[edgeToCheck[1]];
    }

    index = 0;
    face thisFace;
    labelList& fBuffer = cMap.entityBuffer(coupleMap::FACE);
    labelList& fOwner = cMap.entityBuffer(coupleMap::OWNER);
    labelList& fNeighbour = cMap.entityBuffer(coupleMap::NEIGHBOUR);
    labelList& feBuffer = cMap.entityBuffer(coupleMap::FACE_EDGE);

    for (label i = 0; i < nF; i++)
    {
        label fIndex = faceMap[i];
        label own = owner_[fIndex];
        label nei = neighbour_[fIndex];

        // Fetch mapped addressing
        label fOwn = rCellMap.found(own) ? rCellMap[own] : -1;
        label fNei = rCellMap.found(nei) ? rCellMap[nei] : -1;

        if (fOwn > -1)
        {
            // Check if this face is pointed the right way
            if ((fNei > -1) && (fNei < fOwn))
            {
                thisFace = faces_[fIndex].reverseFace();
                fOwner[i] = fNei;
                fNeighbour[i] = fOwn;
            }
            else
            {
                thisFace = faces_[fIndex];
                fOwner[i] = fOwn;

                if (fNei > -1)
                {
                    fNeighbour[i] = fNei;
                }
            }
        }
        else
        {
            // This face is pointed the wrong way.
            thisFace = faces_[fIndex].reverseFace();
            fOwner[i] = fNei;
        }

        const labelList& fEdges = faceEdges_[fIndex];

        forAll(fEdges, indexI)
        {
            fBuffer[index] = rPointMap[thisFace[indexI]];
            feBuffer[index] = rEdgeMap[fEdges[indexI]];

            index++;
        }
    }

    // Fill variable size face-list sizes
    labelList& nfeBuffer = cMap.entityBuffer(coupleMap::NFE_BUFFER);

    forAllConstIter(Map<label>, faceMap, fIter)
    {
        nfeBuffer[fIter.key()] = faces_[fIter()].size();
    }

    // Fill in boundary information
    cMap.entityBuffer(coupleMap::FACE_STARTS) = bdyFaceStarts;
    cMap.entityBuffer(coupleMap::FACE_SIZES) = bdyFaceSizes;
    cMap.entityBuffer(coupleMap::EDGE_STARTS) = bdyEdgeStarts;
    cMap.entityBuffer(coupleMap::EDGE_SIZES) = bdyEdgeSizes;

    labelList& ptBuffer = cMap.entityBuffer(coupleMap::PATCH_ID);

    // Fill types for all but the last one (which is default).
    forAll(boundary, patchI)
    {
        if (isA<processorPolyPatch>(boundary[patchI]))
        {
            // Fetch the neighbouring processor ID
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[patchI])
            );

            // Set processor patches to a special type
            ptBuffer[patchI] = -2 - pp.neighbProcNo();
        }
        else
        {
            // Conventional physical patch. Make an identical
            // map, since physical boundaries are present on
            // all processors.
            ptBuffer[patchI] = patchI;
        }
    }

    // Fill the default patch as 'internal'
    ptBuffer[boundary.size()] = -1;

    // Make a temporary dictionary for patch construction
    dictionary patchDict;

    // Specify the list of patch names and types
    wordList patchNames(ptBuffer.size());

    forAll(patchNames, patchI)
    {
        // Add a new subDictionary
        dictionary patchSubDict;

        if (patchI == patchNames.size() - 1)
        {
            // Artificially set the last patch

            // Set name
            patchNames[patchI] = "defaultPatch";

            // Add type
            patchSubDict.add("type", "empty");

            // Add start / size
            patchSubDict.add("startFace", bdyFaceStarts[patchI]);
            patchSubDict.add("nFaces", bdyFaceSizes[patchI]);
        }
        else
        if (ptBuffer[patchI] <= -2)
        {
            // Back out the neighbouring processor ID
            label neiProcNo = Foam::mag(ptBuffer[patchI] + 2);

            // Set name
            patchNames[patchI] =
            (
                "procBoundary"
              + Foam::name(Pstream::myProcNo())
              + "to"
              + Foam::name(neiProcNo)
            );

            // Add type
            patchSubDict.add("type", "subMeshProcessor");

            // Add start / size
            patchSubDict.add("startFace", bdyFaceStarts[patchI]);
            patchSubDict.add("nFaces", bdyFaceSizes[patchI]);

            // Add processor-specific info
            patchSubDict.add("myProcNo", Pstream::myProcNo());
            patchSubDict.add("neighbProcNo", neiProcNo);
        }
        else
        {
            // Set name
            patchNames[patchI] = boundary[ptBuffer[patchI]].name();

            // Add type
            patchSubDict.add("type", boundary[ptBuffer[patchI]].type());

            // Add start / size
            patchSubDict.add("startFace", bdyFaceStarts[patchI]);
            patchSubDict.add("nFaces", bdyFaceSizes[patchI]);
        }

        // Add subdictionary
        patchDict.add(patchNames[patchI], patchSubDict);
    }

    // Set the autoPtr.
    subMesh.setMesh
    (
        proc,
        new dynamicTopoFvMesh
        (
            (*this),
            IOobject
            (
                fvMesh::defaultRegion,
                time().timeName(),
                time()
            ),
            xferCopy(cMap.pointBuffer()),
            xferCopy(cMap.oldPointBuffer()),
            xferCopy(cMap.edges()),
            xferCopy(cMap.faces()),
            xferCopy(cMap.faceEdges()),
            xferCopy(cMap.owner()),
            xferCopy(cMap.neighbour()),
            bdyFaceStarts,
            bdyFaceSizes,
            bdyEdgeStarts,
            bdyEdgeSizes,
            patchNames,
            patchDict
        )
    );

    // Set maps as built.
    subMesh.setBuiltMaps();

    // For debugging purposes...
    if (debug > 3)
    {
        Pout<< "Writing out sent subMesh for processor: "
            << proc << endl;

        writeVTK
        (
            "psMesh_"
          + Foam::name(Pstream::myProcNo())
          + "to"
          + Foam::name(proc),
            rCellMap.toc()
        );

        // Write out patch information
        Pout<< " faceStarts: " << bdyFaceStarts << nl
            << " faceSizes: " << bdyFaceSizes << nl
            << " edgeStarts: " << bdyEdgeStarts << nl
            << " edgeSizes: " << bdyEdgeSizes << nl
            << " patchTypes: " << ptBuffer << endl;
    }
}


// Build coupled maps for locally coupled patches.
//   - Performs a geometric match initially, since the mesh provides
//     no explicit information for topological coupling.
//   - For each subsequent topology change, maps are updated during
//     the re-ordering stage.
void dynamicTopoFvMesh::buildLocalCoupledMaps()
{
    if (patchCoupling_.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "Building local coupled maps...";
    }

    // Check if a geometric tolerance has been specified.
    const boundBox& box = polyMesh::bounds();

    // Compute tolerance
    scalar tol = geomMatchTol_()*box.mag();

    const polyBoundaryMesh& boundary = boundaryMesh();

    forAll(patchCoupling_, patchI)
    {
        if (!patchCoupling_(patchI))
        {
            continue;
        }

        // Build maps only for the first time.
        // Information is subsequently mapped for each topo-change.
        if (patchCoupling_[patchI].builtMaps())
        {
            continue;
        }

        // Deal with cyclics later
        if (isA<cyclicPolyPatch>(boundary[patchI]))
        {
            continue;
        }

        const coupleMap& cMap = patchCoupling_[patchI].map();

        // Map points on each patch.
        const labelList& mP = boundary[cMap.masterIndex()].meshPoints();
        const labelList& sP = boundary[cMap.slaveIndex()].meshPoints();

        // Sanity check: Do patches have equal number of entities?
        if (mP.size() != sP.size())
        {
            FatalErrorIn("void dynamicTopoFvMesh::buildLocalCoupledMaps()")
                << "Patch point sizes are not consistent."
                << abort(FatalError);
        }

        // Check if maps were read from disk.
        // If so, geometric checking is unnecessary.
        if (cMap.entityMap(coupleMap::POINT).size() == 0)
        {
            // Match points geometrically
            labelList pointMap(mP.size(), -1);

            bool matchedAll =
            (
                matchPoints
                (
                    boundary[cMap.masterIndex()].localPoints(),
                    boundary[cMap.slaveIndex()].localPoints(),
                    scalarField(mP.size(), tol),
                    false,
                    pointMap
                )
            );

            if (!matchedAll)
            {
                FatalErrorIn("void dynamicTopoFvMesh::buildLocalCoupledMaps()")
                    << " Failed to match all points"
                    << " within a tolerance of: " << tol << nl
                    << " matchTol: " << geomMatchTol_() << nl
                    << abort(FatalError);
            }

            // Now build maps
            forAll(mP, indexI)
            {
                label masterIndex = mP[indexI];
                label slaveIndex = sP[pointMap[indexI]];

                cMap.mapSlave
                (
                    coupleMap::POINT,
                    masterIndex,
                    slaveIndex
                );

                cMap.mapMaster
                (
                    coupleMap::POINT,
                    slaveIndex,
                    masterIndex
                );
            }
        }

        const labelListList& mpF = boundary[cMap.masterIndex()].pointFaces();
        const labelListList& spF = boundary[cMap.slaveIndex()].pointFaces();

        // Fetch reference to maps
        const Map<label>& pMap = cMap.entityMap(coupleMap::POINT);
        const Map<label>& eMap = cMap.entityMap(coupleMap::EDGE);
        const Map<label>& fMap = cMap.entityMap(coupleMap::FACE);

        // Now that all points are matched,
        // perform topological matching for higher entities.
        forAll(mP, indexI)
        {
            // Fetch global point indices
            label mp = mP[indexI];
            label sp = pMap[mp];

            // Fetch local point indices
            label lmp = boundary[cMap.masterIndex()].whichPoint(mp);
            label lsp = boundary[cMap.slaveIndex()].whichPoint(sp);

            // Fetch patch starts
            label mStart = boundary[cMap.masterIndex()].start();
            label sStart = boundary[cMap.slaveIndex()].start();

            // Match faces for both 2D and 3D.
            const labelList& mpFaces = mpF[lmp];
            const labelList& spFaces = spF[lsp];

            forAll(mpFaces, faceI)
            {
                // Fetch the global face index
                label mfIndex = (mStart + mpFaces[faceI]);

                if (fMap.found(mfIndex))
                {
                    continue;
                }

                const face& mFace = faces_[mfIndex];

                if (is2D())
                {
                    // Set up a comparison face.
                    face cFace(4);

                    // Configure the face for comparison.
                    forAll(mFace, pointI)
                    {
                        cFace[pointI] = pMap[mFace[pointI]];
                    }

                    bool matched = false;

                    forAll(spFaces, faceJ)
                    {
                        // Fetch the global face index
                        label sfIndex = (sStart + spFaces[faceJ]);

                        const face& sFace = faces_[sfIndex];

                        if (face::compare(cFace, sFace))
                        {
                            // Found the slave. Add a map entry
                            cMap.mapSlave
                            (
                                coupleMap::FACE,
                                mfIndex,
                                sfIndex
                            );

                            cMap.mapMaster
                            (
                                coupleMap::FACE,
                                sfIndex,
                                mfIndex
                            );

                            matched = true;

                            break;
                        }
                    }

                    if (!matched)
                    {
                        FatalErrorIn
                        (
                            "void dynamicTopoFvMesh::buildLocalCoupledMaps()"
                        )
                            << " Failed to match face: "
                            << mfIndex << ": " << mFace << nl
                            << " Comparison face: " << cFace
                            << abort(FatalError);
                    }
                }
                else
                {
                    label slaveFaceIndex = -1;

                    // Set up a comparison face.
                    triFace cFace
                    (
                        pMap[mFace[0]],
                        pMap[mFace[1]],
                        pMap[mFace[2]]
                    );

                    forAll(spFaces, faceJ)
                    {
                        // Fetch the global face index
                        label sfIndex = (sStart + spFaces[faceJ]);

                        const face& sFace = faces_[sfIndex];

                        if (triFace::compare(cFace, triFace(sFace)))
                        {
                            // Found the slave. Add a map entry
                            cMap.mapSlave
                            (
                                coupleMap::FACE,
                                mfIndex,
                                sfIndex
                            );

                            cMap.mapMaster
                            (
                                coupleMap::FACE,
                                sfIndex,
                                mfIndex
                            );

                            slaveFaceIndex = sfIndex;

                            break;
                        }
                    }

                    if (slaveFaceIndex == -1)
                    {
                        FatalErrorIn
                        (
                            "void dynamicTopoFvMesh::buildLocalCoupledMaps()"
                        )
                            << " Failed to match face: "
                            << mfIndex << ": " << mFace << nl
                            << " Comparison face: " << cFace
                            << abort(FatalError);
                    }

                    // Match all edges on this face as well.
                    const labelList& mfEdges = faceEdges_[mfIndex];
                    const labelList& sfEdges = faceEdges_[slaveFaceIndex];

                    forAll(mfEdges, edgeI)
                    {
                        if (eMap.found(mfEdges[edgeI]))
                        {
                            continue;
                        }

                        const edge& mEdge = edges_[mfEdges[edgeI]];

                        // Configure a comparison edge.
                        edge cEdge(pMap[mEdge[0]], pMap[mEdge[1]]);

                        bool matchedEdge = false;

                        forAll(sfEdges, edgeJ)
                        {
                            const edge& sEdge = edges_[sfEdges[edgeJ]];

                            if (cEdge == sEdge)
                            {
                                // Found the slave. Add a map entry
                                cMap.mapSlave
                                (
                                    coupleMap::EDGE,
                                    mfEdges[edgeI],
                                    sfEdges[edgeJ]
                                );

                                cMap.mapMaster
                                (
                                    coupleMap::EDGE,
                                    sfEdges[edgeJ],
                                    mfEdges[edgeI]
                                );

                                matchedEdge = true;

                                break;
                            }
                        }

                        if (!matchedEdge)
                        {
                            FatalErrorIn
                            (
                                "void dynamicTopoFvMesh::"
                                "buildLocalCoupledMaps()"
                            )
                                << " Failed to match edge: "
                                << mfEdges[edgeI] << ": "
                                << mEdge << nl
                                << " cEdge: " << cEdge
                                << abort(FatalError);
                        }
                    }
                }
            }
        }

        // Set maps as built
        patchCoupling_[patchI].setBuiltMaps();
    }

    // Build maps for cyclics
    forAll(patchCoupling_, patchI)
    {
        if (!patchCoupling_(patchI))
        {
            continue;
        }

        // Build maps only for the first time.
        // Information is subsequently mapped for each topo-change.
        if (patchCoupling_[patchI].builtMaps())
        {
            continue;
        }

        // Skip if not cyclic
        if (!isA<cyclicPolyPatch>(boundary[patchI]))
        {
            continue;
        }

        const coupleMap& cMap = patchCoupling_[patchI].map();

        // Match faces and points in one go
        label halfSize = boundary[patchI].size() / 2;
        label patchStart = boundary[patchI].start();

        // Fetch reference to map
        const Map<label>& pointMap = cMap.entityMap(coupleMap::POINT);

        for (label i = 0; i < halfSize; i++)
        {
            label half0Index = (i + patchStart);
            label half1Index = (i + halfSize + patchStart);

            // Map the face
            cMap.mapSlave
            (
                coupleMap::FACE,
                half0Index,
                half1Index
            );

            cMap.mapMaster
            (
                coupleMap::FACE,
                half1Index,
                half0Index
            );

            const face& half0Face = faces_[half0Index];
            const face& half1Face = faces_[half1Index];

            label fS = half0Face.size();

            forAll(half0Face, pointI)
            {
                if (pointMap.found(half0Face[pointI]))
                {
                    continue;
                }

                label masterIndex = half0Face[pointI];
                label slaveIndex = half1Face[(fS - pointI) % fS];

                // Map the point
                cMap.mapSlave
                (
                    coupleMap::POINT,
                    masterIndex,
                    slaveIndex
                );

                cMap.mapMaster
                (
                    coupleMap::POINT,
                    slaveIndex,
                    masterIndex
                );
            }

            if (is2D())
            {
                continue;
            }

            // Fetch reference to map
            const Map<label>& edgeMap = cMap.entityMap(coupleMap::EDGE);

            const labelList& f0Edges = faceEdges_[half0Index];
            const labelList& f1Edges = faceEdges_[half1Index];

            forAll(f0Edges, edgeI)
            {
                if (edgeMap.found(f0Edges[edgeI]))
                {
                    continue;
                }

                const edge& mEdge = edges_[f0Edges[edgeI]];

                // Build a comparison edge
                edge cEdge(pointMap[mEdge[0]], pointMap[mEdge[1]]);

                forAll(f1Edges, edgeJ)
                {
                    const edge& sEdge = edges_[f1Edges[edgeJ]];

                    if (sEdge == cEdge)
                    {
                        // Map the edge
                        cMap.mapSlave
                        (
                            coupleMap::EDGE,
                            f0Edges[edgeI],
                            f1Edges[edgeJ]
                        );

                        cMap.mapMaster
                        (
                            coupleMap::EDGE,
                            f1Edges[edgeJ],
                            f0Edges[edgeI]
                        );

                        break;
                    }
                }
            }
        }
    }

    if (debug)
    {
        Info<< "Done." << endl;
    }
}


// Build coupled maps for coupled processor patches
void dynamicTopoFvMesh::buildProcessorCoupledMaps()
{
    if (procIndices_.empty())
    {
        return;
    }

    Map<label> nPrc;

    const polyBoundaryMesh& boundary = boundaryMesh();

    // Build a list of nearest neighbours.
    forAll(boundary, patchI)
    {
        if (isA<processorPolyPatch>(boundary[patchI]))
        {
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[patchI])
            );

            nPrc.insert(pp.neighbProcNo(), patchI);
        }
    }

    // Wait for all transfers to complete.
    meshOps::waitForBuffers();

    // Put un-matched faces in a list.
    labelHashSet unMatchedFaces;

    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        coupledMesh& rPM = recvMeshes_[pI];
        const coupleMap& cMap = rPM.map();
        const labelList& ptBuffer = cMap.entityBuffer(coupleMap::PATCH_ID);

        // Fetch reference to patch starts / sizes
        const labelList& faceStarts = cMap.entityBuffer(coupleMap::FACE_STARTS);
        const labelList& faceSizes = cMap.entityBuffer(coupleMap::FACE_SIZES);
        const labelList& edgeStarts = cMap.entityBuffer(coupleMap::EDGE_STARTS);
        const labelList& edgeSizes = cMap.entityBuffer(coupleMap::EDGE_SIZES);

        // Make a temporary dictionary for patch construction
        dictionary patchDict;

        // Specify the list of patch names and types
        wordList patchNames(ptBuffer.size());

        forAll(patchNames, patchI)
        {
            // Add a new subDictionary
            dictionary patchSubDict;

            if (patchI == patchNames.size() - 1)
            {
                // Artificially set the last patch

                // Set name
                patchNames[patchI] = "defaultPatch";

                // Add type
                patchSubDict.add("type", "empty");

                // Add start / size
                patchSubDict.add("startFace", faceStarts[patchI]);
                patchSubDict.add("nFaces", faceSizes[patchI]);
            }
            else
            if (ptBuffer[patchI] <= -2)
            {
                // Back out the neighbouring processor ID
                label neiProcNo = Foam::mag(ptBuffer[patchI] + 2);

                // Set name
                patchNames[patchI] =
                (
                    "procBoundary"
                  + Foam::name(proc)
                  + "to"
                  + Foam::name(neiProcNo)
                );

                // Add type
                patchSubDict.add("type", "subMeshProcessor");

                // Add start / size
                patchSubDict.add("startFace", faceStarts[patchI]);
                patchSubDict.add("nFaces", faceSizes[patchI]);

                // Add processor-specific info
                patchSubDict.add("myProcNo", proc);
                patchSubDict.add("neighbProcNo", neiProcNo);
            }
            else
            {
                // Set name
                patchNames[patchI] = boundary[ptBuffer[patchI]].name();

                // Add type
                patchSubDict.add("type", boundary[ptBuffer[patchI]].type());

                // Add start / size
                patchSubDict.add("startFace", faceStarts[patchI]);
                patchSubDict.add("nFaces", faceSizes[patchI]);
            }

            // Add subdictionary
            patchDict.add(patchNames[patchI], patchSubDict);
        }

        // Set the autoPtr.
        rPM.setMesh
        (
            proc,
            new dynamicTopoFvMesh
            (
                (*this),
                IOobject
                (
                    fvMesh::defaultRegion,
                    time().timeName(),
                    time()
                ),
                xferCopy(cMap.pointBuffer()),
                xferCopy(cMap.oldPointBuffer()),
                xferCopy(cMap.edges()),
                xferCopy(cMap.faces()),
                xferCopy(cMap.faceEdges()),
                xferCopy(cMap.owner()),
                xferCopy(cMap.neighbour()),
                faceStarts,
                faceSizes,
                edgeStarts,
                edgeSizes,
                patchNames,
                patchDict
            )
        );

        if (debug > 3)
        {
            Pout<< "Writing out received subMesh for processor: "
                << proc << endl;

            rPM.subMesh().writeVTK
            (
                "rPatchMesh_"
              + Foam::name(Pstream::myProcNo())
              + "to"
              + Foam::name(proc),
                identity(cMap.nEntities(coupleMap::CELL)),
                3, false, true
            );
        }

        // First, topologically map points based on subMeshPoints.
        const labelList& mP = cMap.subMeshPoints();

        label nShared = cMap.nEntities(coupleMap::SHARED_POINT);
        label nGlobal = cMap.nEntities(coupleMap::GLOBAL_POINT);

        // Sanity check: Do sub-mesh point sizes match?
        if (mP.size() != nShared)
        {
            FatalErrorIn("void dynamicTopoFvMesh::buildProcessorCoupledMaps()")
                << " Sub-mesh point sizes don't match." << nl
                << " My procID: " << Pstream::myProcNo() << nl
                << " Neighbour processor: " << proc << nl
                << " mP size: " << mP.size() << nl
                << " nShared: " << nShared << nl
                << abort(FatalError);
        }

        if (nPrc.found(proc))
        {
            // This is an immediate neighbour.
            // All patch points must be matched.

            // Fetch addressing for this patch.
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[nPrc[proc]])
            );

            const labelList& neiPoints = pp.neighbPoints();

            if (debug)
            {
                // Find all occurrences of multiply connected points
                labelList mIdx = findIndices(neiPoints, -1);

                if (mIdx.size())
                {
                    // Write out for post-processing
                    UIndirectList<label> myIdx(mP, mIdx);
                    writeVTK("mcPoints", labelList(myIdx), 0);

                    FatalErrorIn
                    (
                        "void dynamicTopoFvMesh::buildProcessorCoupledMaps()"
                    )
                        << " Multiply connected point." << nl
                        << " My procID: " << Pstream::myProcNo() << nl
                        << " Neighbour processor: " << proc << nl
                        << " Neighbour Indices: " << mIdx << nl
                        << " My indices: " << myIdx << nl
                        << abort(FatalError);
                }
            }

            forAll(mP, pointI)
            {
                // Add a map entry
                cMap.mapSlave
                (
                    coupleMap::POINT,
                    mP[pointI],
                    neiPoints[pointI]
                );

                cMap.mapMaster
                (
                    coupleMap::POINT,
                    neiPoints[pointI],
                    mP[pointI]
                );
            }
        }

        // Prepare point maps for globally shared points
        const labelList& gpB = cMap.entityBuffer(coupleMap::POINT);

        // Sanity check: Do global point buffer sizes match?
        if ((mP.size() + gpB.size()) != nGlobal)
        {
            FatalErrorIn("void dynamicTopoFvMesh::buildProcessorCoupledMaps()")
                << " Global point sizes don't match"
                << " for processor: " << proc << nl
                << " mP size: " << mP.size() << nl
                << " gpB size: " << gpB.size() << nl
                << " Total: " << (mP.size() + gpB.size()) << nl
                << " nGlobal: " << nGlobal << nl
                << abort(FatalError);
        }

        // Look at globalPoint addressing for its information.
        const globalMeshData& gD = polyMesh::globalData();
        const labelList& spA = gD.sharedPointAddr();
        const labelList& spL = gD.sharedPointLabels();

        forAll(gpB, pointI)
        {
            // Find my equivalent global point
            label maIndex = findIndex(spA, gpB[pointI]);

            // Add a map entry
            cMap.mapSlave
            (
                coupleMap::POINT,
                spL[maIndex],
                (pointI + nShared)
            );

            cMap.mapMaster
            (
                coupleMap::POINT,
                (pointI + nShared),
                spL[maIndex]
            );
        }

        if (debug > 1)
        {
            const boundBox& box = polyMesh::bounds();
            const pointField& sPoints = rPM.subMesh().points_;
            const Map<label>& pMap = cMap.entityMap(coupleMap::POINT);

            // Compute tolerance
            scalar tol = geomMatchTol_()*box.mag();

            forAllConstIter(Map<label>, pMap, pIter)
            {
                scalar dist = mag(points_[pIter.key()] - sPoints[pIter()]);

                if (dist > tol)
                {
                    FatalErrorIn
                    (
                        "void dynamicTopoFvMesh::buildProcessorCoupledMaps()"
                    )
                        << " Failed to match point: " << pIter.key()
                        << ": " << points_[pIter.key()]
                        << " with point: " << pIter()
                        << ": " << sPoints[pIter()] << nl
                        << " Missed by: " << dist << nl
                        << " Tolerance: " << tol << nl
                        << abort(FatalError);
                }
            }
        }

        // Set up a comparison face.
        face cFace(4);

        // Now match all faces connected to master points.
        if (nPrc.found(proc))
        {
            // This is an immediate neighbour.
            label mStart = boundary[nPrc[proc]].start();
            label mSize  = boundary[nPrc[proc]].size();

            // Abbreviate for convenience
            const dynamicTopoFvMesh& sMesh = rPM.subMesh();
            const Map<label>& pMap = cMap.entityMap(coupleMap::POINT);
            const Map<label>& eMap = cMap.entityMap(coupleMap::EDGE);
            const Map<label>& fMap = cMap.entityMap(coupleMap::FACE);

            // Fetch global pointFaces for the slave.
            const faceList& slaveFaces = sMesh.faces();
            const labelListList& spF = sMesh.pointFaces();

            // Match patch faces for both 2D and 3D.
            for (label i = 0; i < mSize; i++)
            {
                // Fetch the global face index
                label mfIndex = (mStart + i);

                if (fMap.found(mfIndex))
                {
                    continue;
                }

                const face& mFace = faces_[mfIndex];

                // Configure the face for comparison.
                forAll(mFace, pointI)
                {
                    cFace[pointI] = pMap[mFace[pointI]];
                }

                // Fetch pointFaces for the zeroth point.
                const labelList& spFaces = spF[cFace[0]];

                label sFaceIndex = -1, compareVal = -1;

                forAll(spFaces, faceJ)
                {
                    // Fetch the global face index
                    label sfIndex = spFaces[faceJ];

                    const face& sFace = slaveFaces[sfIndex];

                    // Check for triangle face optimization
                    if (mFace.size() == 3)
                    {
                        compareVal =
                        (
                            triFace::compare
                            (
                                triFace(cFace[0], cFace[1], cFace[2]),
                                triFace(sFace[0], sFace[1], sFace[2])
                            )
                        );
                    }
                    else
                    {
                        compareVal = face::compare(cFace, sFace);
                    }

                    if (compareVal)
                    {
                        // Found the slave. Add a map entry
                        cMap.mapSlave
                        (
                            coupleMap::FACE,
                            mfIndex,
                            sfIndex
                        );

                        cMap.mapMaster
                        (
                            coupleMap::FACE,
                            sfIndex,
                            mfIndex
                        );

                        sFaceIndex = sfIndex;

                        break;
                    }
                }

                if (sFaceIndex == -1)
                {
                    unMatchedFaces.insert(mfIndex);

                    continue;
                }

                // No need to match edges in 2D
                if (is2D())
                {
                    continue;
                }

                // Match all edges on this face as well.
                const labelList& mfEdges = faceEdges_[mfIndex];
                const labelList& sfEdges = sMesh.faceEdges_[sFaceIndex];

                forAll(mfEdges, edgeI)
                {
                    if (eMap.found(mfEdges[edgeI]))
                    {
                        continue;
                    }

                    const edge& mEdge = edges_[mfEdges[edgeI]];

                    // Configure a comparison edge.
                    edge cEdge(pMap[mEdge[0]], pMap[mEdge[1]]);

                    bool matchedEdge = false;

                    forAll(sfEdges, edgeJ)
                    {
                        const edge& sEdge = sMesh.edges_[sfEdges[edgeJ]];

                        if (cEdge == sEdge)
                        {
                            // Found the slave. Add a map entry
                            cMap.mapSlave
                            (
                                coupleMap::EDGE,
                                mfEdges[edgeI],
                                sfEdges[edgeJ]
                            );

                            cMap.mapMaster
                            (
                                coupleMap::EDGE,
                                sfEdges[edgeJ],
                                mfEdges[edgeI]
                            );

                            matchedEdge = true;

                            break;
                        }
                    }

                    if (!matchedEdge)
                    {
                        // Write out the edge
                        writeVTK("mEdge", mfEdges[edgeI], 1);

                        FatalErrorIn
                        (
                            "void dynamicTopoFvMesh::"
                            "buildProcessorCoupledMaps()"
                        )
                            << " Failed to match edge: "
                            << mfEdges[edgeI] << ": "
                            << mEdge << nl
                            << " cEdge: " << cEdge
                            << " for processor: " << proc
                            << abort(FatalError);
                    }
                }
            }
        }

        // Attempt to match other edges in 3D, provided any common ones exist.
        //  - This will handle point-only coupling situations for edges.
        if (is3D())
        {
            // Fetch Maps for this processor
            const Map<label>& edgeMap = cMap.entityMap(coupleMap::EDGE);
            const Map<label>& pointMap = cMap.entityMap(coupleMap::POINT);

            forAllConstIter(Map<label>, pointMap, pIter)
            {
                const labelList& pEdges = pointEdges_[pIter.key()];

                forAll(pEdges, edgeI)
                {
                    const label eIndex = pEdges[edgeI];

                    // Only pick boundary edges
                    if (whichEdgePatch(eIndex) == -1)
                    {
                        continue;
                    }

                    // Skip mapped edges
                    if (edgeMap.found(eIndex))
                    {
                        continue;
                    }

                    const edge& eCheck = edges_[eIndex];

                    // Check if a map exists for the other point
                    label sIndex =
                    (
                        cMap.findSlave
                        (
                            coupleMap::POINT,
                            eCheck.otherVertex(pIter.key())
                        )
                    );

                    if (sIndex > -1)
                    {
                        const dynamicTopoFvMesh& mesh = rPM.subMesh();

                        // Configure a check edge
                        edge cEdge(sIndex, pIter());

                        const labelList& spEdges = mesh.pointEdges_[sIndex];

                        forAll(spEdges, edgeJ)
                        {
                            const edge& sEdge = mesh.edges_[spEdges[edgeJ]];

                            if (sEdge == cEdge)
                            {
                                // Found the slave. Add a map entry
                                cMap.mapSlave
                                (
                                    coupleMap::EDGE,
                                    eIndex,
                                    spEdges[edgeJ]
                                );

                                cMap.mapMaster
                                (
                                    coupleMap::EDGE,
                                    spEdges[edgeJ],
                                    eIndex
                                );

                                break;
                            }
                        }
                    }
                }
            }

            // Prepare processor point and edge maps
            //  - For each point / edge on the subMesh, create
            //    a list of processors (other than this one)
            //    that are connected to it
            Map<labelList>& pEdgeMap = cMap.subMeshEdgeMap();
            Map<labelList>& pPointMap = cMap.subMeshPointMap();

            const dynamicTopoFvMesh& sMesh = rPM.subMesh();
            const polyBoundaryMesh& bdy = sMesh.boundaryMesh();

            forAll(bdy, pI)
            {
                if (!isA<processorPolyPatch>(bdy[pI]))
                {
                    continue;
                }

                const processorPolyPatch& pp =
                (
                    refCast<const processorPolyPatch>(bdy[pI])
                );

                label neiProcNo = pp.neighbProcNo();
                label myProcNo = Pstream::myProcNo();

                if (neiProcNo == myProcNo)
                {
                    continue;
                }

                label mStart = pp.start();
                label mSize  = pp.size();

                for (label i = 0; i < mSize; i++)
                {
                    const face& f = sMesh.faces_[i + mStart];
                    const labelList& fe = sMesh.faceEdges_[i + mStart];

                    forAll(f, j)
                    {
                        const label pIndex = f[j];
                        const label eIndex = fe[j];

                        Map<labelList>::iterator pIt = pPointMap.find(pIndex);

                        if (pIt == pPointMap.end())
                        {
                            pPointMap.insert(pIndex, labelList(1, neiProcNo));
                        }
                        else
                        {
                            if (findIndex(pIt(), neiProcNo) == -1)
                            {
                                meshOps::sizeUpList(neiProcNo, pIt());
                            }
                        }

                        Map<labelList>::iterator eIt = pEdgeMap.find(eIndex);

                        if (eIt == pEdgeMap.end())
                        {
                            pEdgeMap.insert(eIndex, labelList(1, neiProcNo));
                        }
                        else
                        {
                            if (findIndex(eIt(), neiProcNo) == -1)
                            {
                                meshOps::sizeUpList(neiProcNo, eIt());
                            }
                        }
                    }
                }
            }
        }

        if (unMatchedFaces.size())
        {
            // Write out the face
            writeVTK("mFaces", unMatchedFaces.toc(), 2);

            FatalErrorIn
            (
                "void dynamicTopoFvMesh::buildProcessorCoupledMaps()"
            )
                << " Unmatched faces were found for processor: " << proc
                << abort(FatalError);
        }
    }

    // Use subMesh and global points to identify additional
    // processors, and update maps accordingly
    if (is3D())
    {
        forAll(procIndices_, pI)
        {
            const label neiProcNo = procIndices_[pI];
            const coupleMap& cMi = recvMeshes_[pI].map();

            const labelList& smPts = cMi.subMeshPoints();
            const List<labelPair>& gpp = cMi.globalProcPoints();

            forAll(procIndices_, pJ)
            {
                if (pJ == pI)
                {
                    continue;
                }

                const coupleMap& cMj = recvMeshes_[pJ].map();

                Map<labelList>& pPointMap = cMj.subMeshPointMap();

                forAll(smPts, pointI)
                {
                    label mP = smPts[pointI];
                    label sP = cMj.findSlave(coupleMap::POINT, mP);

                    if (sP == -1)
                    {
                        continue;
                    }

                    Map<labelList>::iterator pIt = pPointMap.find(sP);

                    if (pIt == pPointMap.end())
                    {
                        pPointMap.insert(sP, labelList(1, neiProcNo));
                    }
                    else
                    {
                        if (findIndex(pIt(), neiProcNo) == -1)
                        {
                            meshOps::sizeUpList(neiProcNo, pIt());
                        }
                    }
                }

                forAll(gpp, pointI)
                {
                    label mP = gpp[pointI].first();
                    label sP = cMj.findSlave(coupleMap::POINT, mP);

                    if (sP == -1)
                    {
                        continue;
                    }

                    Map<labelList>::iterator pIt = pPointMap.find(sP);

                    if (pIt == pPointMap.end())
                    {
                        pPointMap.insert(sP, labelList(1, neiProcNo));
                    }
                    else
                    {
                        if (findIndex(pIt(), neiProcNo) == -1)
                        {
                            meshOps::sizeUpList(neiProcNo, pIt());
                        }
                    }
                }
            }
        }
    }
}


// Introduce a new processor patch to the mesh
label dynamicTopoFvMesh::createProcessorPatch(const label proc)
{
    if (!Pstream::parRun())
    {
        return -2;
    }

    // Find index in list of processors
    label pI = findIndex(procIndices_, proc);

    // Check if an entry already exists
    if (pI > -1)
    {
        label patchIndex = sendMeshes_[pI].map().patchIndex();

        if (patchIndex > -1)
        {
            return patchIndex;
        }
    }

    // Get the new patch index,
    // and increment the number of patches
    label patchID = nPatches_++;

    // Set the new patch index in patchMaps
    if (isSubMesh_)
    {
        if (pI == -1)
        {
            pI = procIndices_.size();

            procIndices_.setSize(pI + 1);
            sendMeshes_.setSize(pI + 1);
        }

        // Create a basic entry on the subMesh
        procIndices_[pI] = proc;

        sendMeshes_.set
        (
            pI,
            new coupledMesh
            (
                *this,               // Reference to this mesh
                is2D(),              // 2D or 3D
                false,               // Not local
                true,                // Sent to neighbour
                patchID,             // Patch index
                proc,                // Master index
                -1                   // Slave index
            )
        );
    }
    else
    {
        if (pI == -1)
        {
            Pout<< " Could not find index for processor: " << proc
                << " in indices: " << procIndices_
                << abort(FatalError);
        }

        sendMeshes_[pI].map().patchIndex() = patchID;
        recvMeshes_[pI].map().patchIndex() = patchID;
    }

    // Size up patches, and copy old information
    label prevPatchID = patchID - 1;

    patchSizes_.setSize(nPatches_, 0);
    oldPatchSizes_.setSize(nPatches_, 0);

    patchStarts_.setSize
    (
        nPatches_,
        patchStarts_[prevPatchID] + patchSizes_[prevPatchID]
    );
    oldPatchStarts_.setSize
    (
        nPatches_,
        oldPatchStarts_[prevPatchID] + oldPatchSizes_[prevPatchID]
    );

    edgePatchSizes_.setSize(nPatches_, 0);
    oldEdgePatchSizes_.setSize(nPatches_, 0);

    edgePatchStarts_.setSize
    (
        nPatches_,
        edgePatchStarts_[prevPatchID] + edgePatchSizes_[prevPatchID]
    );
    oldEdgePatchStarts_.setSize
    (
        nPatches_,
        oldEdgePatchStarts_[prevPatchID] + oldEdgePatchSizes_[prevPatchID]
    );

    patchNMeshPoints_.setSize(nPatches_, 0);
    oldPatchNMeshPoints_.setSize(nPatches_, 0);

    if (debug)
    {
        Pout<< " dynamicTopoFvMesh::createProcessorPatch :"
            << " Created new patch for processor: " << proc << nl
            << " On subMesh: " << isSubMesh_ << nl
            << " pI: " << pI << nl
            << " patchID: " << patchID << nl
            << " oldPatchStarts: " << oldPatchStarts_ << nl
            << " oldPatchSizes: " << oldPatchSizes_ << nl
            << " patchStarts: " << patchStarts_ << nl
            << " patchSizes: " << patchSizes_ << nl
            << endl;
    }

    // Return the new patch index
    return patchID;
}


// If a patch is of processor type, get the neighbouring processor ID
label dynamicTopoFvMesh::getNeighbourProcessor(const label patch) const
{
    if (!Pstream::parRun())
    {
        return -1;
    }

    label neiProcNo = -1;

    const polyBoundaryMesh& boundary = boundaryMesh();

    if (patch == -1)
    {
        return -1;
    }
    else
    if (patch < boundary.size())
    {
        if (isA<processorPolyPatch>(boundary[patch]))
        {
            const processorPolyPatch& pp =
            (
                refCast<const processorPolyPatch>(boundary[patch])
            );

            // Set the neighbour processor ID
            neiProcNo = pp.neighbProcNo();
        }
    }
    else
    {
        // New processor patch

        // Find the neighbour processor ID
        forAll(procIndices_, pI)
        {
            label patchID = sendMeshes_[pI].map().patchIndex();

            if (patch == patchID)
            {
                neiProcNo = procIndices_[pI];
                break;
            }
        }

        if (neiProcNo == -1)
        {
            // An index should have been defined by now
            Pout<< " Patch: " << patch
                << " was not defined on patch subMeshes."
                << abort(FatalError);
        }
    }

    return neiProcNo;
}


// If the number of patches have changed during run-time,
// reset boundaries with new processor patches
void dynamicTopoFvMesh::resetBoundaries()
{
    // Prepare a new list of patches
    List<polyPatch*> patches(nPatches_);

    // Fetch reference to existing boundary
    // - The removeBoundary member merely resets
    //   boundary size, so this reference is safe
    const polyBoundaryMesh& polyBoundary = boundaryMesh();

    // Copy all existing patches first
    label nOldPatches = polyBoundary.size();

    for (label patchI = 0; patchI < nOldPatches; patchI++)
    {
        // Clone the patch
        patches[patchI] = polyBoundary[patchI].clone(polyBoundary).ptr();
    }

    // Create new processor patches
    for (label patchI = nOldPatches; patchI < nPatches_; patchI++)
    {
        // Make a temporary dictionary for patch construction
        dictionary patchDict;

        // Back out the neighbour processor ID
        label neiProcNo = getNeighbourProcessor(patchI);

        // Add relevant info
        patchDict.add("type", "processor");
        patchDict.add("startFace", oldPatchStarts_[patchI]);
        patchDict.add("nFaces", oldPatchSizes_[patchI]);
        patchDict.add("myProcNo", Pstream::myProcNo());
        patchDict.add("neighbProcNo", neiProcNo);

        // Set the pointer
        patches[patchI] =
        (
            polyPatch::New
            (
                "procBoundary"
              + Foam::name(Pstream::myProcNo())
              + "to"
              + Foam::name(neiProcNo),
                patchDict,
                patchI,
                polyBoundary
            ).ptr()
        );
    }

    // Remove the old boundary
    fvMesh::removeFvBoundary();

    // Add patches, but don't calculate geometry, etc
    fvMesh::addFvPatches(patches, false);

    // Since all fvPatches in fvBoundaryMesh have been reset,
    // fvPatch references in all fvPatchFields are no longer
    // valid, and must therefore be remapped.
    const fvBoundaryMesh& bdy = fvMesh::boundary();

    coupledMesh::resizeBoundaries<volScalarField>(*this, bdy);
    coupledMesh::resizeBoundaries<volVectorField>(*this, bdy);
    coupledMesh::resizeBoundaries<volSphericalTensorField>(*this, bdy);
    coupledMesh::resizeBoundaries<volSymmTensorField>(*this, bdy);
    coupledMesh::resizeBoundaries<volTensorField>(*this, bdy);

    coupledMesh::resizeBoundaries<surfaceScalarField>(*this, bdy);
    coupledMesh::resizeBoundaries<surfaceVectorField>(*this, bdy);
    coupledMesh::resizeBoundaries<surfaceSphericalTensorField>(*this, bdy);
    coupledMesh::resizeBoundaries<surfaceSymmTensorField>(*this, bdy);
    coupledMesh::resizeBoundaries<surfaceTensorField>(*this, bdy);
}


// Convenience macro for field sub-setting
#define sendFieldsOfType(type, info, names, types, offset, map, stream)        \
{                                                                              \
    {                                                                          \
        scalar zeroValue = pTraits<scalar>::zero;                              \
        info.send<type##ScalarField>                                           \
        (names[0 + offset], types[0 + offset], zeroValue, map, stream);        \
    }                                                                          \
    {                                                                          \
        vector zeroValue = pTraits<vector>::zero;                              \
        info.send<type##VectorField>                                           \
        (names[1 + offset], types[1 + offset], zeroValue, map, stream);        \
    }                                                                          \
    {                                                                          \
        sphericalTensor zeroValue = pTraits<sphericalTensor>::zero;            \
        info.send<type##SphericalTensorField>                                  \
        (names[2 + offset], types[2 + offset], zeroValue, map, stream);        \
    }                                                                          \
    {                                                                          \
        symmTensor zeroValue = pTraits<symmTensor>::zero;                      \
        info.send<type##SymmTensorField>                                       \
        (names[3 + offset], types[3 + offset], zeroValue, map, stream);        \
    }                                                                          \
    {                                                                          \
        tensor zeroValue = pTraits<tensor>::zero;                              \
        info.send<type##TensorField>                                           \
        (names[4 + offset], types[4 + offset], zeroValue, map, stream);        \
    }                                                                          \
}


// Initialize subMesh field transfers for mapping
void dynamicTopoFvMesh::initFieldTransfers
(
    wordList& types,
    List<wordList>& names,
    List<List<char> >& sendBuffer,
    List<List<char> >& recvBuffer
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (debug)
    {
        Info<< " void dynamicTopoFvMesh::initFieldTransfers() :"
            << " Initializing subMesh field transfers."
            << endl;
    }

    // Clear out mesh geometry, since those
    // are to be wiped out after topo-changes anyway.
    fvMesh::clearOut();
    polyMesh::resetMotion();

    // Size up wordLists
    //  - Five templated volFields
    //  - Five templated surfaceFields
    names.setSize(10);
    types.setSize(10);

    // Fill in field-types
    types[0] = volScalarField::typeName;
    types[1] = volVectorField::typeName;
    types[2] = volSphericalTensorField::typeName;
    types[3] = volSymmTensorField::typeName;
    types[4] = volTensorField::typeName;

    types[5] = surfaceScalarField::typeName;
    types[6] = surfaceVectorField::typeName;
    types[7] = surfaceSphericalTensorField::typeName;
    types[8] = surfaceSymmTensorField::typeName;
    types[9] = surfaceTensorField::typeName;

    // Send / recv buffers for field names
    List<char> fieldNameSendBuffer, fieldNameRecvBuffer;

    if (Pstream::master())
    {
        PtrList<OStringStream> fieldNameStream(1);

        // Set the stream
        fieldNameStream.set(0, new OStringStream(IOstream::BINARY));

        // Alias for convenience
        OStringStream& fNStream = fieldNameStream[0];

        // Fetch field-names by type
        forAll(types, typeI)
        {
            // Get all fields of type
            names[typeI] = objectRegistry::names(types[typeI]);
        }

        // Send field names to Ostream
        fNStream << names;

        // Size up buffers and fill contents
        string contents = fNStream.str();
        const char* ptr = contents.data();

        fieldNameSendBuffer.setSize(contents.size());

        forAll(fieldNameSendBuffer, i)
        {
            fieldNameSendBuffer[i] = *ptr++;
        }

        // Clear the stream
        fieldNameStream.set(0, NULL);

        if (debug > 4)
        {
            Info<< " Registered fields: " << names << endl;
        }

        // Send names to slaves
        for (label proc = 1; proc < Pstream::nProcs(); proc++)
        {
            // Send buffer to processor
            meshOps::pWrite(proc, fieldNameSendBuffer.size());
            meshOps::pWrite(proc, fieldNameSendBuffer);
        }
    }
    else
    {
        // Receive names from master
        label recvBufferSize = -1;
        meshOps::pRead(Pstream::masterNo(), recvBufferSize);

        // Size up buffer and schedule receive
        fieldNameRecvBuffer.setSize(recvBufferSize);
        meshOps::pRead(Pstream::masterNo(), fieldNameRecvBuffer);
    }

    // Wait for transfers to complete
    meshOps::waitForBuffers();

    // Convert names from stringStream
    if (!Pstream::master())
    {
        string contents
        (
            fieldNameRecvBuffer.begin(),
            fieldNameRecvBuffer.size()
        );

        // Set the stream
        IStringStream fieldNameStream(contents, IOstream::BINARY);

        // Get names from stream
        fieldNameStream >> names;
    }

    label nProcs = procIndices_.size();

    // Size up buffers
    sendBuffer.setSize(nProcs);
    recvBuffer.setSize(nProcs);

    // Size up the send stringStream
    PtrList<OStringStream> stream(nProcs);

    // Now fill in subMesh fields
    forAll(procIndices_, pI)
    {
        label proc = procIndices_[pI];

        const coupledMesh& cInfo = sendMeshes_[pI];

        // Initialize stream
        stream.set(pI, new OStringStream(IOstream::BINARY));

        // Subset and send volFields to stream
        const labelList& cMap = cInfo.map().cellMap();

        sendFieldsOfType(vol, cInfo, names, types, 0, cMap, stream[pI]);

        // Subset and send surfaceFields to stream
        const labelList& fMap = cInfo.map().internalFaceMap();

        sendFieldsOfType(surface, cInfo, names, types, 5, fMap, stream[pI]);

        // Size up buffers and fill contents
        string contents = stream[pI].str();
        const char* ptr = contents.data();

        sendBuffer[pI].setSize(contents.size());

        forAll(sendBuffer[pI], i)
        {
            sendBuffer[pI][i] = *ptr++;
        }

        // Clear the stream
        stream.set(pI, NULL);

        // Send buffer to processor
        meshOps::pWrite(proc, sendBuffer[pI].size());
        meshOps::pWrite(proc, sendBuffer[pI]);

        // Receive buffer from processor
        label recvBufferSize = -1;
        meshOps::pRead(proc, recvBufferSize);

        // Size up buffer and schedule receive
        recvBuffer[pI].setSize(recvBufferSize);
        meshOps::pRead(proc, recvBuffer[pI]);
    }

    // We won't wait for all transfers to complete for the moment.
    // Meanwhile, do some other useful work, if possible.
}


// Synchronize field transfers for mapping
void dynamicTopoFvMesh::syncFieldTransfers
(
    wordList& types,
    List<wordList>& names,
    List<List<char> >& recvBuffer
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (debug)
    {
        Info<< " void dynamicTopoFvMesh::syncFieldTransfers() :"
            << " Synchronizing subMesh field transfers."
            << endl;
    }

    // Wait for all transfers to complete.
    meshOps::waitForBuffers();

    label nProcs = procIndices_.size();

    // Size up stringStream
    PtrList<IStringStream> stream(nProcs);

    // Size up field PtrLists

    // Volume Fields
    List<PtrList<volScalarField> > vsF(nProcs);
    List<PtrList<volVectorField> > vvF(nProcs);
    List<PtrList<volSphericalTensorField> > vsptF(nProcs);
    List<PtrList<volSymmTensorField> > vsytF(nProcs);
    List<PtrList<volTensorField> > vtF(nProcs);

    // Surface fields
    List<PtrList<surfaceScalarField> > ssF(nProcs);
    List<PtrList<surfaceVectorField> > svF(nProcs);
    List<PtrList<surfaceSphericalTensorField> > ssptF(nProcs);
    List<PtrList<surfaceSymmTensorField> > ssytF(nProcs);
    List<PtrList<surfaceTensorField> > stF(nProcs);

    const polyBoundaryMesh& polyBoundary = boundaryMesh();

    // Determine number of physical (non-processor) patches
    label nPhysical = 0;

    forAll(polyBoundary, patchI)
    {
        if (isA<processorPolyPatch>(polyBoundary[patchI]))
        {
            continue;
        }

        nPhysical++;
    }

    // Keep track of extra / total entities
    label nTotalCells = nOldCells_, nTotalIntFaces = nOldInternalFaces_;
    labelList nTotalPatchFaces(SubList<label>(oldPatchSizes_, nPhysical));

    // Allocate reverse maps
    List<labelList> irvMaps(nProcs), irsMaps(nProcs);
    List<labelListList> brMaps(nProcs, labelListList(nPhysical));

    forAll(procIndices_, pI)
    {
        const coupledMesh& cInfo = recvMeshes_[pI];

        // Convert buffer to string
        string contents(recvBuffer[pI].begin(), recvBuffer[pI].size());

        // Initialize stream
        stream.set(pI, new IStringStream(contents, IOstream::BINARY));

        // Construct dictionary from stream
        dictionary dict(stream[pI]);

        // Count the number of additional entities
        const coupleMap& cMap = cInfo.map();

        // Fetch size from subMesh
        label nCells = cMap.nEntities(coupleMap::CELL);
        label nIntFaces = cMap.nEntities(coupleMap::INTERNAL_FACE);

        // Set field pointers
        cInfo.setField(names[0], dict.subDict(types[0]), nCells, vsF[pI]);
        cInfo.setField(names[1], dict.subDict(types[1]), nCells, vvF[pI]);
        cInfo.setField(names[2], dict.subDict(types[2]), nCells, vsptF[pI]);
        cInfo.setField(names[3], dict.subDict(types[3]), nCells, vsytF[pI]);
        cInfo.setField(names[4], dict.subDict(types[4]), nCells, vtF[pI]);

        cInfo.setField(names[5], dict.subDict(types[5]), nIntFaces, ssF[pI]);
        cInfo.setField(names[6], dict.subDict(types[6]), nIntFaces, svF[pI]);
        cInfo.setField(names[7], dict.subDict(types[7]), nIntFaces, ssptF[pI]);
        cInfo.setField(names[8], dict.subDict(types[8]), nIntFaces, ssytF[pI]);
        cInfo.setField(names[9], dict.subDict(types[9]), nIntFaces, stF[pI]);

        // Set rmap for this processor
        irvMaps[pI] = (labelField(identity(nCells)) + nTotalCells);
        irsMaps[pI] = (labelField(identity(nIntFaces)) + nTotalIntFaces);

        // Update count
        nTotalCells += nCells;
        nTotalIntFaces += nIntFaces;

        // Fetch patch sizes from subMesh
        const labelList& nPatchFaces =
        (
            cMap.entityBuffer(coupleMap::FACE_SIZES)
        );

        // Loop over physical patches
        forAll(nTotalPatchFaces, patchI)
        {
            // Fetch patch size from subMesh
            label nFaces = nPatchFaces[patchI];

            // Set rmap for this patch
            brMaps[pI][patchI] =
            (
                labelField(identity(nFaces))
              + nTotalPatchFaces[patchI]
            );

            // Update patch-face count
            nTotalPatchFaces[patchI] += nFaces;
        }
    }

    // Prepare internal mappers
    labelList cellAddressing(nTotalCells, 0);
    labelList faceAddressing(nTotalIntFaces, 0);

    // Set identity map for first nCells,
    // and map from cell[0] for the rest
    for (label i = 0; i < nOldCells_; i++)
    {
        cellAddressing[i] = i;
    }

    // Set identity map for first nIntFaces,
    // and map from face[0] for the rest
    for (label i = 0; i < nOldInternalFaces_; i++)
    {
        faceAddressing[i] = i;
    }

    coupledMesh::subMeshMapper vMap
    (
        nOldCells_,
        cellAddressing
    );

    coupledMesh::subMeshMapper sMap
    (
        nOldInternalFaces_,
        faceAddressing
    );

    // Prepare boundary mappers
    labelListList patchAddressing(nPhysical);
    PtrList<coupledMesh::subMeshMapper> bMap(nPhysical);

    forAll(bMap, patchI)
    {
        // Prepare patch mappers
        patchAddressing[patchI].setSize(nTotalPatchFaces[patchI], 0);

        // Set identity map for first nPatchFaces,
        // and map from patch-face[0] for the rest
        for (label i = 0; i < oldPatchSizes_[patchI]; i++)
        {
            patchAddressing[patchI][i] = i;
        }

        // Set the boundary mapper pointer
        bMap.set
        (
            patchI,
            new coupledMesh::subMeshMapper
            (
                oldPatchSizes_[patchI],
                patchAddressing[patchI]
            )
        );
    }

    // Loop through all volFields and re-size
    // to accomodate additional cells / faces
    coupledMesh::resizeMap(names[0], *this, vMap, irvMaps, bMap, brMaps, vsF);
    coupledMesh::resizeMap(names[1], *this, vMap, irvMaps, bMap, brMaps, vvF);
    coupledMesh::resizeMap(names[2], *this, vMap, irvMaps, bMap, brMaps, vsptF);
    coupledMesh::resizeMap(names[3], *this, vMap, irvMaps, bMap, brMaps, vsytF);
    coupledMesh::resizeMap(names[4], *this, vMap, irvMaps, bMap, brMaps, vtF);

    coupledMesh::resizeMap(names[5], *this, sMap, irsMaps, bMap, brMaps, ssF);
    coupledMesh::resizeMap(names[6], *this, sMap, irsMaps, bMap, brMaps, svF);
    coupledMesh::resizeMap(names[7], *this, sMap, irsMaps, bMap, brMaps, ssptF);
    coupledMesh::resizeMap(names[8], *this, sMap, irsMaps, bMap, brMaps, ssytF);
    coupledMesh::resizeMap(names[9], *this, sMap, irsMaps, bMap, brMaps, stF);
}


// Initialize coupled boundary ordering
// - Assumes that faces_ and points_ are consistent
// - Assumes that patchStarts_ and patchSizes_ are consistent
void dynamicTopoFvMesh::initCoupledBoundaryOrdering
(
    List<pointField>& centres,
    List<pointField>& anchors
) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    for (label pI = 0; pI < nPatches_; pI++)
    {
        label neiProcNo = getNeighbourProcessor(pI);

        if (neiProcNo == -1)
        {
            continue;
        }

        // Check if this is a master processor patch.
        label start = patchStarts_[pI];
        label size = patchSizes_[pI];

        // Prepare centres and anchors
        centres[pI].setSize(size, vector::zero);
        anchors[pI].setSize(size, vector::zero);

        if (Pstream::myProcNo() < neiProcNo)
        {
            forAll(centres[pI], fI)
            {
                centres[pI][fI] = faces_[fI + start].centre(points_);
                anchors[pI][fI] = points_[faces_[fI + start][0]];
            }

            {
                if (debug > 3)
                {
                    Pout<< "Sending to: " << neiProcNo
                        << " nCentres: " << size << endl;

                    // Write out my centres to disk
                    meshOps::writeVTK
                    (
                        (*this),
                        "centres_"
                      + Foam::name(Pstream::myProcNo())
                      + "to"
                      + Foam::name(neiProcNo),
                        size, size, size,
                        centres[pI]
                    );
                }

                // Ensure that we're sending the right size
                meshOps::pWrite(neiProcNo, size);
            }

            // Send information to neighbour
            meshOps::pWrite(neiProcNo, centres[pI]);
            meshOps::pWrite(neiProcNo, anchors[pI]);
        }
        else
        {
            {
                label nEntities = -1;

                // Ensure that we're receiving the right size
                meshOps::pRead(neiProcNo, nEntities);

                if (debug > 3)
                {
                    Pout<< " Recving from: " << neiProcNo
                        << " nCentres: " << size << endl;
                }

                if (nEntities != size)
                {
                    // Set the size to what we're receiving
                    centres[pI].setSize(nEntities, vector::zero);
                    anchors[pI].setSize(nEntities, vector::zero);

                    // Issue a warning now, and let
                    // syncCoupledBoundaryOrdering fail later
                    WarningIn
                    (
                        "void dynamicTopoFvMesh::"
                        "initCoupledBoundaryOrdering() const"
                    )
                        << "Incorrect send / recv sizes: " << nl
                        << " nEntities: " << nEntities << nl
                        << " size: " << size << nl
                        << endl;
                }
            }

            // Schedule receive from neighbour
            meshOps::pRead(neiProcNo, centres[pI]);
            meshOps::pRead(neiProcNo, anchors[pI]);
        }
    }
}


// Synchronize coupled boundary ordering
bool dynamicTopoFvMesh::syncCoupledBoundaryOrdering
(
    List<pointField>& centres,
    List<pointField>& anchors,
    labelListList& patchMaps,
    labelListList& rotations
) const
{
    bool anyChange = false, failedPatchMatch = false;

    // Calculate centres and tolerances for any slave patches
    List<scalarField> slaveTols(nPatches_);
    List<pointField> slaveCentres(nPatches_);

    // Handle locally coupled patches
    forAll(patchCoupling_, pI)
    {
        if (patchCoupling_(pI))
        {
            const coupleMap& cMap = patchCoupling_[pI].map();

            label masterPatch = cMap.masterIndex();
            label slavePatch = cMap.slaveIndex();

            label mSize = patchSizes_[masterPatch];
            label sSize = patchSizes_[slavePatch];

            label mStart = patchStarts_[masterPatch];
            label sStart = patchStarts_[slavePatch];

            if (mSize != sSize)
            {
                Pout<< " Patch size mismatch: " << nl
                    << "  mSize: " << mSize << nl
                    << "  sSize: " << sSize << nl
                    << "  mStart: " << mStart << nl
                    << "  sStart: " << sStart << nl
                    << "  master: " << masterPatch << nl
                    << "  slave: " << slavePatch << nl
                    << abort(FatalError);
            }

            // Calculate centres / tolerances
            anchors[masterPatch].setSize(mSize, vector::zero);
            centres[masterPatch].setSize(mSize, vector::zero);

            slaveTols[slavePatch].setSize(sSize, 0.0);
            slaveCentres[slavePatch].setSize(sSize, vector::zero);

            for (label fI = 0; fI < mSize; fI++)
            {
                point& mfa = anchors[masterPatch][fI];
                point& mfc = centres[masterPatch][fI];
                point& sfc = slaveCentres[slavePatch][fI];

                const face& mFace = faces_[fI + mStart];
                const face& sFace = faces_[fI + sStart];

                // Calculate anchor / centres
                mfa = points_[mFace[0]];
                mfc = mFace.centre(points_);
                sfc = sFace.centre(points_);

                scalar maxLen = -GREAT;

                forAll(sFace, fpI)
                {
                    maxLen = max(maxLen, mag(points_[sFace[fpI]] - sfc));
                }

                slaveTols[slavePatch][fI] = geomMatchTol_()*maxLen;
            }

            // For cyclics, additionally test for halves,
            // and transform points as necessary.
            labelList halfMap;

            if (masterPatch == slavePatch)
            {
                scalar featureCos = 0.9;
                label halfSize = (mSize / 2);
                label half0Index = 0, half1Index = 0;

                // Size up halfMap and slaveAnchors
                halfMap.setSize(mSize, -1);

                scalarField half1Tols(halfSize, 0.0);
                vectorField half1Anchors(halfSize, vector::zero);

                // Fetch face-normal for the first face
                vector fhN = boundaryMesh()[masterPatch].faceAreas()[0];

                // Normalize
                fhN /= mag(fhN) + VSMALL;

                // Fetch transform type
                const cyclicPolyPatch& cyclicPatch =
                (
                    refCast<const cyclicPolyPatch>
                    (
                        boundaryMesh()[masterPatch]
                    )
                );

                bool translate =
                (
                    cyclicPatch.transform()
                 == cyclicPolyPatch::TRANSLATIONAL
                );

                for (label fI = 0; fI < mSize; fI++)
                {
                    const face& mFace = faces_[fI + mStart];

                    // Compute normal
                    vector fN = mFace.normal(points_);

                    // Normalize
                    fN /= mag(fN) + VSMALL;

                    // Check which half this face belongs to...
                    bool isFirstHalf = ((fhN & fN) > featureCos);

                    if (isFirstHalf)
                    {
                        vector mA = anchors[masterPatch][fI];
                        vector mC = centres[masterPatch][fI];

                        // Set master anchor
                        anchors[masterPatch][half0Index] = mA;

                        // Transform master points
                        if (translate)
                        {
                            mA += cyclicPatch.separationVector();
                            mC += cyclicPatch.separationVector();
                        }
                        else
                        {
                            // Assume constant transform tensor
                            mA = Foam::transform(cyclicPatch.forwardT()[0], mA);
                            mC = Foam::transform(cyclicPatch.forwardT()[0], mC);
                        }

                        // Set the transformed anchor / centre
                        half1Anchors[half0Index] = mA;
                        centres[masterPatch][half0Index] = mC;
                        half1Tols[half0Index] = slaveTols[masterPatch][fI];

                        // Set halfMap
                        halfMap[fI] = half0Index;

                        half0Index++;
                    }
                    else
                    {
                        // Set the slave centre
                        const scalar& sT = slaveTols[slavePatch][fI];
                        const point& sC = slaveCentres[slavePatch][fI];

                        // Duplicate slaveTols for first / second half
                        slaveTols[slavePatch][half1Index] = sT;
                        slaveCentres[slavePatch][half1Index] = sC;

                        // Set halfMap
                        halfMap[fI] = half1Index + halfSize;

                        half1Index++;
                    }
                }

                // Sanity check for halfSize
                if (half0Index != halfSize || half1Index != halfSize)
                {
                    Pout<< " Master: " << masterPatch << nl
                        << " Slave: " << slavePatch << nl
                        << " half0Index: " << half0Index << nl
                        << " half1Index: " << half1Index << nl
                        << " halfSize: " << halfSize << nl
                        << " failed to divide into halves."
                        << abort(FatalError);
                }

                // Resize lists
                centres[masterPatch].setSize(halfSize);
                slaveCentres[slavePatch].setSize(halfSize);

                // Set slave tols / anchors
                forAll(half1Anchors, fI)
                {
                    slaveTols[masterPatch][fI + halfSize] = half1Tols[fI];
                    anchors[masterPatch][fI + halfSize] = half1Anchors[fI];
                }
            }

            // Fetch reference to slave map
            labelList& patchMap = patchMaps[slavePatch];

            // Try zero separation automatic matching
            bool matchedAll =
            (
                (mSize == sSize)
             &&
                matchPoints
                (
                    slaveCentres[slavePatch],
                    centres[masterPatch],
                    slaveTols[slavePatch],
                    false,
                    patchMap
                )
            );

            // Write out master centres to disk
            if (debug > 3 || !matchedAll)
            {
                meshOps::writeVTK
                (
                    (*this),
                    "localMasterCentres_"
                  + Foam::name(masterPatch)
                  + '_'
                  + Foam::name(slavePatch),
                    mSize, mSize, mSize,
                    centres[masterPatch]
                );

                meshOps::writeVTK
                (
                    (*this),
                    "localSlaveCentres_"
                  + Foam::name(masterPatch)
                  + '_'
                  + Foam::name(slavePatch),
                    sSize, sSize, sSize,
                    slaveCentres[slavePatch]
                );

                writeVTK
                (
                    "localMasterFaces_"
                  + Foam::name(masterPatch)
                  + '_'
                  + Foam::name(slavePatch),
                    identity(mSize) + mStart,
                    2, false, true
                );

                writeVTK
                (
                    "localSlaveFaces_"
                  + Foam::name(masterPatch)
                  + '_'
                  + Foam::name(slavePatch),
                    identity(sSize) + sStart,
                    2, false, true
                );

                if (!matchedAll)
                {
                    Pout<< " Master: " << masterPatch
                        << " Slave: " << slavePatch
                        << " mSize: " << mSize << " sSize: " << sSize
                        << " failed on match for face centres."
                        << endl;

                    // Failed on match-for-all, so patchMaps
                    // will be invalid. Bail out for now.
                    failedPatchMatch = true;

                    continue;
                }
            }

            // Renumber patchMap for cyclics
            if (masterPatch == slavePatch)
            {
                label halfSize = (mSize / 2);

                forAll(halfMap, fI)
                {
                    if (halfMap[fI] >= halfSize)
                    {
                        halfMap[fI] =
                        (
                            patchMap[halfMap[fI] - halfSize]
                          + halfSize
                        );
                    }
                }

                // Transfer map
                patchMap.transfer(halfMap);
            }

            // Fetch reference to rotation map
            labelList& rotation = rotations[slavePatch];

            // Initialize rotation to zero
            rotation.setSize(sSize, 0);

            // Set rotation.
            forAll(patchMap, oldFaceI)
            {
                label newFaceI = patchMap[oldFaceI];

                const point& anchor = anchors[slavePatch][newFaceI];
                const scalar& faceTol = slaveTols[slavePatch][oldFaceI];
                const face& checkFace = faces_[sStart + oldFaceI];

                label anchorFp = -1;
                scalar minDSqr = GREAT;

                forAll(checkFace, fpI)
                {
                    scalar dSqr = magSqr(anchor - points_[checkFace[fpI]]);

                    if (dSqr < minDSqr)
                    {
                        minDSqr = dSqr;
                        anchorFp = fpI;
                    }
                }

                if (anchorFp == -1 || mag(minDSqr) > faceTol)
                {
                    FatalErrorIn
                    (
                        "\n"
                        "void dynamicTopoFvMesh::"
                        "syncCoupledBoundaryOrdering\n"
                        "(\n"
                        "    List<pointField>& centres,\n"
                        "    List<pointField>& anchors,\n"
                        "    labelListList& patchMaps,\n"
                        "    labelListList& rotations\n"
                        ") const\n"
                    )
                        << "Cannot find anchor: " << anchor << nl
                        << " Face: " << checkFace << nl
                        << " Vertices: "
                        << UIndirectList<point>(points_, checkFace) << nl
                        << " on patch: " << slavePatch << nl
                        << " faceTol: " << faceTol << nl
                        << " newFaceI: " << newFaceI << nl
                        << " oldFaceI: " << oldFaceI << nl
                        << " mSize: " << mSize << nl
                        << " sSize: " << sSize << nl
                        << abort(FatalError);
                }
                else
                {
                    // Positive rotation
                    //  - Set for old face. Will be rotated later
                    //    during the shuffling stage
                    rotation[oldFaceI] =
                    (
                        (checkFace.size() - anchorFp) % checkFace.size()
                    );
                }
            }

            // Set the flag
            anyChange = true;
        }
    }

    if (failedPatchMatch)
    {
        FatalErrorIn
        (
            "\n"
            "void dynamicTopoFvMesh::"
            "syncCoupledBoundaryOrdering\n"
            "(\n"
            "    List<pointField>& centres,\n"
            "    List<pointField>& anchors,\n"
            "    labelListList& patchMaps,\n"
            "    labelListList& rotations\n"
            ") const\n"
        )
            << " Matching for local patches failed. " << nl
            << abort(FatalError);
    }

    if (!Pstream::parRun())
    {
        return anyChange;
    }

    for (label pI = 0; pI < nPatches_; pI++)
    {
        label neiProcNo = getNeighbourProcessor(pI);

        if (neiProcNo == -1)
        {
            continue;
        }

        // Check if this is a slave processor patch.
        label start = patchStarts_[pI];
        label size = patchSizes_[pI];

        if (Pstream::myProcNo() > neiProcNo)
        {
            slaveTols[pI].setSize(size, 0.0);
            slaveCentres[pI].setSize(size, vector::zero);

            forAll(slaveCentres[pI], fI)
            {
                point& fc = slaveCentres[pI][fI];

                const face& checkFace = faces_[fI + start];

                // Calculate centre
                fc = checkFace.centre(points_);

                scalar maxLen = -GREAT;

                forAll(checkFace, fpI)
                {
                    maxLen = max(maxLen, mag(points_[checkFace[fpI]] - fc));
                }

                slaveTols[pI][fI] = geomMatchTol_()*maxLen;
            }

            // Write out my centres to disk
            if (debug > 3)
            {
                meshOps::writeVTK
                (
                    (*this),
                    "slaveCentres_"
                  + Foam::name(Pstream::myProcNo())
                  + '_'
                  + Foam::name(neiProcNo),
                    size, size, size,
                    slaveCentres[pI]
                );
            }
        }
    }

    // Wait for transfers before continuing.
    meshOps::waitForBuffers();

    for (label pI = 0; pI < nPatches_; pI++)
    {
        label neiProcNo = getNeighbourProcessor(pI);

        if (neiProcNo == -1)
        {
            continue;
        }

        // Check if this is a master processor patch.
        labelList& patchMap = patchMaps[pI];
        labelList& rotation = rotations[pI];

        // Initialize map and rotation
        patchMap.setSize(patchSizes_[pI], -1);
        rotation.setSize(patchSizes_[pI], 0);

        if (Pstream::myProcNo() < neiProcNo)
        {
            // Do nothing (i.e. identical mapping, zero rotation).
            forAll(patchMap, pfI)
            {
                patchMap[pfI] = pfI;
            }
        }
        else
        {
            // Try zero separation automatic matching
            label mSize = centres[pI].size();
            label sSize = slaveCentres[pI].size();

            bool matchedAll =
            (
                (mSize == sSize)
             &&
                matchPoints
                (
                    slaveCentres[pI],
                    centres[pI],
                    slaveTols[pI],
                    false,
                    patchMap
                )
            );

            // Write out centres to disk
            if (debug > 3 || !matchedAll)
            {
                meshOps::writeVTK
                (
                    (*this),
                    "masterCentres_"
                  + Foam::name(Pstream::myProcNo())
                  + '_'
                  + Foam::name(neiProcNo),
                    mSize, mSize, mSize,
                    centres[pI]
                );

                meshOps::writeVTK
                (
                    (*this),
                    "slaveCentres_"
                  + Foam::name(Pstream::myProcNo())
                  + '_'
                  + Foam::name(neiProcNo),
                    sSize, sSize, sSize,
                    slaveCentres[pI]
                );

                writeVTK
                (
                    "procPatchFaces_"
                  + Foam::name(Pstream::myProcNo())
                  + '_'
                  + Foam::name(neiProcNo),
                    identity(patchSizes_[pI]) + patchStarts_[pI],
                    2, false, true
                );

                if (!matchedAll)
                {
                    Pout<< " Patch: " << pI
                        << " Processor: " << neiProcNo
                        << " mSize: " << mSize << " sSize: " << sSize
                        << " failed on match for face centres."
                        << endl;

                    // Failed on match-for-all, so patchMaps
                    // will be invalid. Bail out for now.
                    failedPatchMatch = true;

                    continue;
                }
            }

            label start = patchStarts_[pI];

            // Set rotation.
            forAll(patchMap, oldFaceI)
            {
                label newFaceI = patchMap[oldFaceI];

                const point& anchor = anchors[pI][newFaceI];
                const scalar& faceTol = slaveTols[pI][oldFaceI];
                const face& checkFace = faces_[start + oldFaceI];

                label anchorFp = -1;
                scalar minDSqr = GREAT;

                forAll(checkFace, fpI)
                {
                    scalar dSqr = magSqr(anchor - points_[checkFace[fpI]]);

                    if (dSqr < minDSqr)
                    {
                        minDSqr = dSqr;
                        anchorFp = fpI;
                    }
                }

                if (anchorFp == -1 || mag(minDSqr) > faceTol)
                {
                    FatalErrorIn
                    (
                        "\n"
                        "void dynamicTopoFvMesh::"
                        "syncCoupledBoundaryOrdering\n"
                        "(\n"
                        "    List<pointField>& centres,\n"
                        "    List<pointField>& anchors,\n"
                        "    labelListList& patchMaps,\n"
                        "    labelListList& rotations\n"
                        ") const\n"
                    )
                        << "Cannot find anchor: " << anchor << nl
                        << " Face: " << checkFace << nl
                        << " Vertices: "
                        << UIndirectList<point>(points_, checkFace) << nl
                        << " on patch: " << pI
                        << " for processor: " << neiProcNo
                        << abort(FatalError);
                }
                else
                {
                    // Positive rotation
                    //  - Set for old face. Will be rotated later
                    //    during the shuffling stage
                    rotation[oldFaceI] =
                    (
                        (checkFace.size() - anchorFp) % checkFace.size()
                    );
                }
            }

            // Set the flag
            anyChange = true;
        }
    }

    // Reduce failures across processors
    reduce(failedPatchMatch, orOp<bool>());

    // Write out processor patch faces if a failure was encountered
    if (failedPatchMatch)
    {
        for (label pI = 0; pI < nPatches_; pI++)
        {
            label neiProcNo = getNeighbourProcessor(pI);

            if (neiProcNo == -1)
            {
                continue;
            }

            writeVTK
            (
                "patchFaces_"
              + Foam::name(Pstream::myProcNo())
              + '_'
              + Foam::name(neiProcNo),
                identity(patchSizes_[pI]) + patchStarts_[pI],
                2, false, true
            );
        }

        FatalErrorIn
        (
            "\n"
            "void dynamicTopoFvMesh::"
            "syncCoupledBoundaryOrdering\n"
            "(\n"
            "    List<pointField>& centres,\n"
            "    List<pointField>& anchors,\n"
            "    labelListList& patchMaps,\n"
            "    labelListList& rotations\n"
            ") const\n"
        )
            << " Matching for processor patches failed. " << nl
            << " Patch faces written out to disk." << nl
            << abort(FatalError);
    }

    return anyChange;
}


// Fill buffers with length-scale info
// and exchange across processors.
void dynamicTopoFvMesh::exchangeLengthBuffers()
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (!edgeRefinement_)
    {
        return;
    }

    forAll(procIndices_, pI)
    {
        coupledMesh& sPM = sendMeshes_[pI];
        coupledMesh& rPM = recvMeshes_[pI];

        const coupleMap& scMap = sPM.map();
        const coupleMap& rcMap = rPM.map();

        if (scMap.slaveIndex() == Pstream::myProcNo())
        {
            // Clear existing buffer
            sPM.subMesh().lengthScale_.clear();

            // Fill in buffers to send.
            sPM.subMesh().lengthScale_.setSize
            (
                scMap.nEntities(coupleMap::CELL),
                -1.0
            );

            Map<label>& cellMap = scMap.entityMap(coupleMap::CELL);

            forAllIter(Map<label>, cellMap, cIter)
            {
                sPM.subMesh().lengthScale_[cIter.key()] = lengthScale_[cIter()];
            }

            meshOps::pWrite
            (
                scMap.masterIndex(),
                sPM.subMesh().lengthScale_
            );

            if (debug > 4)
            {
                Pout<< "Sending to: "
                    << scMap.masterIndex()
                    << " nCells: "
                    << scMap.nEntities(coupleMap::CELL)
                    << endl;
            }
        }

        if (rcMap.masterIndex() == Pstream::myProcNo())
        {
            // Clear existing buffer
            rPM.subMesh().lengthScale_.clear();

            // Schedule receipt from neighbour
            rPM.subMesh().lengthScale_.setSize
            (
                rcMap.nEntities(coupleMap::CELL),
                -1.0
            );

            // Schedule for receipt
            meshOps::pRead
            (
                rcMap.slaveIndex(),
                rPM.subMesh().lengthScale_
            );

            if (debug > 4)
            {
                Pout<< "Receiving from: "
                    << rcMap.slaveIndex()
                    << " nCells: "
                    << rcMap.nEntities(coupleMap::CELL)
                    << endl;
            }
        }
    }

    // Wait for transfers before continuing.
    meshOps::waitForBuffers();

    if (debug > 4)
    {
        forAll(procIndices_, pI)
        {
            const coupledMesh& rPM = recvMeshes_[pI];
            const coupleMap& rcMap = rPM.map();

            if (rcMap.masterIndex() == Pstream::myProcNo())
            {
                rPM.subMesh().writeVTK
                (
                    "lengthScale_" + Foam::name(rcMap.slaveIndex()),
                    identity(rPM.subMesh().nCells()),
                    3,
                    false,
                    false,
                    rPM.subMesh().lengthScale_
                );
            }
        }
    }
}


// Implementing the fillTables operation for coupled edges
bool dynamicTopoFvMesh::coupledFillTables
(
    const label eIndex,
    scalar& minQuality,
    labelList& m,
    PtrList<scalarListList>& Q,
    PtrList<labelListList>& K,
    PtrList<labelListList>& triangulations
) const
{
    bool success = false;

    if (locallyCoupledEntity(eIndex))
    {
        // Fill tables for the slave edge.
        label sI = -1;

        // Determine the slave index.
        forAll(patchCoupling_, patchI)
        {
            if (patchCoupling_(patchI))
            {
                const label edgeEnum  = coupleMap::EDGE;
                const coupleMap& cMap = patchCoupling_[patchI].map();

                if ((sI = cMap.findSlave(edgeEnum, eIndex)) > -1)
                {
                    break;
                }
            }
        }

        if (sI == -1)
        {
            FatalErrorIn
            (
                "\n"
                "bool dynamicTopoFvMesh::coupledFillTables\n"
                "(\n"
                "    const label eIndex,\n"
                "    const scalar minQuality,\n"
                "    labelList& m,\n"
                "    PtrList<scalarListList>& Q,\n"
                "    PtrList<labelListList>& K,\n"
                "    PtrList<labelListList>& triangulations\n"
                ") const\n"
            )
                << "Coupled maps were improperly specified." << nl
                << " Slave index not found for: " << nl
                << " Edge: " << eIndex << nl
                << abort(FatalError);
        }

        // Turn off switch temporarily.
        unsetCoupledModification();

        labelList hullV(10);

        // Call fillTables for the slave edge.
        success = fillTables(sI, minQuality, m, hullV, Q, K, triangulations, 1);

        // Turn it back on.
        setCoupledModification();
    }
    else
    if (processorCoupledEntity(eIndex))
    {
        const edge& checkEdge = edges_[eIndex];
        const labelList& eFaces = edgeFaces_[eIndex];

        // Reset minQuality
        minQuality = GREAT;

        // Need to build alternate addressing / point-list
        // for swaps across processors.
        DynamicList<point> parPts(10);
        dynamicLabelList parVtx(10);

        bool closed = true;
        label nPoints = 0, nProcs = 0;
        label otherPoint = -1, nextPoint = -1;

        // Define a line for this edge
        linePointRef lpr
        (
            points_[checkEdge.start()],
            points_[checkEdge.end()]
        );

        // Define tangent-to-edge / edge-centre
        vector te = -lpr.vec(), xe = lpr.centre();

        // Normalize tangent
        te /= mag(te) + VSMALL;

        // Fill-in vertices for this processor
        forAll(eFaces, faceI)
        {
            label fPatch = whichPatch(eFaces[faceI]);
            label neiProc = getNeighbourProcessor(fPatch);

            if
            (
                neiProc > -1
             && priority(neiProc, lessOp<label>(), Pstream::myProcNo())
            )
            {
                // This edge should not be here
                Pout<< " Edge: " << eIndex
                    << " is talking to processor: " << neiProc
                    << abort(FatalError);
            }

            // This face is either internal or on a physical boundary
            const face& checkFace = faces_[eFaces[faceI]];

            // Find the isolated point
            meshOps::findIsolatedPoint
            (
                checkFace,
                checkEdge,
                otherPoint,
                nextPoint
            );

            // Insert point and index
            parPts.append(points_[otherPoint]);
            parVtx.append(nPoints);

            // Physical patch: Is this an appropriate start face?
            //  - If yes, swap with first index
            if (fPatch > -1 && neiProc == -1)
            {
                closed = false;

                if (nextPoint == checkEdge[0])
                {
                    Foam::Swap(parPts[0], parPts[nPoints]);
                }
            }

            // Increment point index
            nPoints++;
        }

        // Now look through processors, and add their points
        forAll(procIndices_, pI)
        {
            label proc = procIndices_[pI];

            // Fetch reference to subMesh
            const coupleMap& cMap = recvMeshes_[pI].map();
            const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

            label sI = -1, sP = -1;

            if ((sI = cMap.findSlave(coupleMap::EDGE, eIndex)) == -1)
            {
                continue;
            }

            // If this is a new edge, bail out for now.
            // This will be handled at the next time-step.
            if (sI >= mesh.nOldEdges_)
            {
                return false;
            }

            const edge& slaveEdge = mesh.edges_[sI];
            const labelList& seFaces = mesh.edgeFaces_[sI];

            // Determine the point index that corresponds to checkEdge[0]
            edge cE
            (
                cMap.findMaster(coupleMap::POINT, slaveEdge[0]),
                cMap.findMaster(coupleMap::POINT, slaveEdge[1])
            );

            if (cE[0] == checkEdge[0])
            {
                sP = slaveEdge[0];
            }
            else
            if (cE[1] == checkEdge[0])
            {
                sP = slaveEdge[1];
            }
            else
            {
                label meP = whichEdgePatch(eIndex);
                word mN(meP < 0 ? "Internal" : boundaryMesh()[meP].name());

                label seP = mesh.whichEdgePatch(sI);
                word sN(seP < 0 ? "Internal" : mesh.boundaryMesh()[seP].name());

                Pout<< " Can't find master point: " << checkEdge[0] << nl
                    << " on master: " << eIndex << "::" << checkEdge
                    << " Patch: " << mN << nl
                    << " with slave: " << sI << "::" << slaveEdge
                    << " Patch: " << sN << nl
                    << " cE: " << cE << " on proc: " << proc << nl
                    << abort(FatalError);
            }

            forAll(seFaces, faceI)
            {
                label slavePatch = mesh.whichPatch(seFaces[faceI]);

                // If this is talking to a lower-ranked processor,
                // skip the insertion step.
                label neiProc = mesh.getNeighbourProcessor(slavePatch);

                if
                (
                    neiProc > -1
                 && priority(neiProc, lessOp<label>(), proc)
                )
                {
                    continue;
                }

                // This face is either internal or on a physical boundary
                const face& slaveFace = mesh.faces_[seFaces[faceI]];

                // Find the isolated point
                meshOps::findIsolatedPoint
                (
                    slaveFace,
                    slaveEdge,
                    otherPoint,
                    nextPoint
                );

                // Insert point and index
                parPts.append(mesh.points_[otherPoint]);
                parVtx.append(nPoints);

                // Physical patch: Is this an appropriate start face?
                //  - If yes, swap with first index
                if (slavePatch > -1 && neiProc == -1)
                {
                    closed = false;

                    if (nextPoint == sP)
                    {
                        Foam::Swap(parPts[0], parPts[nPoints]);
                    }
                }

                // Increment point index
                nPoints++;
            }

            nProcs++;
        }

        // Sort points / indices in counter-clockwise order
        SortableList<scalar> angles(nPoints, 0.0);

        // Fetch 2 * pi
        scalar twoPi = mathematicalConstant::twoPi;

        // Define a base direction
        // from the start point
        vector dir = (parPts[0] - xe);
        dir -= ((dir & te) * te);
        dir /= mag(dir) + VSMALL;

        // Local coordinate system
        coordinateSystem cs("cs", xe, te, dir);

        for (label i = 1; i < nPoints; i++)
        {
            // Convert to local csys and determine angle
            vector local = cs.localPosition(parPts[i]);
            scalar angle = atan2(local.y(), local.x());

            // Account for 3rd and 4th quadrants
            angles[i] = (angle < 0.0 ? angle + twoPi : angle);
        }

        // Sort by angle
        angles.sort();

        // Reorder points and transfer
        List<point> sortedParPts(nPoints);

        const labelList& indices = angles.indices();

        forAll(sortedParPts, pointI)
        {
            sortedParPts[pointI] = parPts[indices[pointI]];
        }

        parPts.transfer(sortedParPts);

        // Fill the last two points for the edge
        edge parEdge(-1, -1);

        parPts.append(lpr.start());
        parEdge[0] = nPoints++;

        parPts.append(lpr.end());
        parEdge[1] = nPoints++;

        // Compute minQuality with this loop
        minQuality = computeMinQuality(parEdge, parVtx, parPts, closed);

        if (debug > 4 || minQuality < 0.0)
        {
            // Write out edge connectivity
            writeEdgeConnectivity(eIndex);

            meshOps::writeVTK
            (
                (*this),
                "parPts_" + Foam::name(eIndex),
                parPts.size(),
                parPts.size(),
                parPts.size(),
                pointField(parPts)
            );

            if (minQuality < 0.0)
            {
                Switch isClosed(closed);

                Pout<< " * * * Error in fillTables * * * " << nl
                    << " Edge: " << eIndex << " :: " << checkEdge << nl
                    << " minQuality: " << minQuality << nl
                    << " Closed: " << isClosed << nl
                    << abort(FatalError);
            }
        }

        // Fill in the size
        m[0] = parVtx.size();

        // Check if a table-resize is necessary
        if (m[0] > maxTetsPerEdge_)
        {
            if (allowTableResize_)
            {
                // Resize the tables to account for
                // more tets per edge
                label& mtpe = const_cast<label&>(maxTetsPerEdge_);

                mtpe = m[0];

                // Clear tables for this index.
                Q[0].clear();
                K[0].clear();
                triangulations[0].clear();

                // Resize for this index.
                initTables(m, Q, K, triangulations);
            }
            else
            {
                // Can't resize. Bail out.
                return false;
            }
        }

        // Fill dynamic programming tables
        fillTables
        (
            parEdge,
            minQuality,
            m[0],
            parVtx,
            parPts,
            Q[0],
            K[0],
            triangulations[0]
        );

        success = true;
    }

    return success;
}


// Initialize coupled weights calculation
void dynamicTopoFvMesh::initCoupledWeights()
{
    if (!Pstream::parRun())
    {
        return;
    }

    forAll(procIndices_, pI)
    {
        // Fetch reference to subMesh
        const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

        mesh.cells();
        mesh.cellCells();
        mesh.cellPoints();

        const polyBoundaryMesh& sBoundary = mesh.boundaryMesh();

        forAll(sBoundary, patchI)
        {
            sBoundary[patchI].faceFaces();
        }
    }
}


// Additional mapping contributions for coupled entities
void dynamicTopoFvMesh::computeCoupledWeights
(
    const label index,
    const label dimension,
    labelList& parents,
    scalarField& weights,
    vectorField& centres,
    bool output
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    // Fetch offsets from mapper
    const labelList& cStarts = mapper_->cellStarts();
    const labelListList& pSizes = mapper_->patchSizes();

    if (dimension == 2)
    {
        dynamicLabelList faceParents(10);

        label patchIndex = whichPatch(index);

        forAll(procIndices_, pI)
        {
            // Ensure that the patch is physical
            if (patchIndex < 0 || patchIndex >= pSizes[pI].size())
            {
                Pout<< " Face: " << index
                    << " Patch: " << patchIndex
                    << " does not belong to a physical patch." << nl
                    << " nPhysicalPatches: " << pSizes[pI].size()
                    << abort(FatalError);
            }

            // Fetch reference to subMesh
            const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

            // Prepare lists
            labelList coupleObjects;
            scalarField coupleWeights;
            vectorField coupleCentres;

            // Convex-set algorithm for faces
            faceSetAlgorithm faceAlgorithm
            (
                mesh,
                oldPoints_,
                edges_,
                faces_,
                cells_,
                owner_,
                neighbour_
            );

            // Initialize the bounding box
            faceAlgorithm.computeNormFactor(index);

            // Loop through all subMesh faces, and check for bounds
            const polyBoundaryMesh& boundary = mesh.boundaryMesh();

            label pSize = boundary[patchIndex].size();
            label pStart = boundary[patchIndex].start();

            for (label faceI = pStart; faceI < (pStart + pSize); faceI++)
            {
                if (faceAlgorithm.contains(faceI))
                {
                    faceParents.append(faceI);
                }
            }

            // Obtain weighting factors for this face.
            faceAlgorithm.computeWeights
            (
                index,
                pStart,
                faceParents,
                boundary[patchIndex].faceFaces(),
                coupleObjects,
                coupleWeights,
                coupleCentres,
                output
            );

            // Add contributions with offsets
            if (coupleObjects.size())
            {
                // Fetch patch size on master
                label patchSize = boundaryMesh()[patchIndex].size();

                // Resize lists
                label oldSize = parents.size();

                parents.setSize(oldSize + coupleObjects.size());
                weights.setSize(oldSize + coupleObjects.size());
                centres.setSize(oldSize + coupleObjects.size());

                forAll(coupleObjects, indexI)
                {
                    parents[indexI + oldSize] =
                    (
                        patchSize + coupleObjects[indexI]
                    );

                    weights[indexI + oldSize] = coupleWeights[indexI];
                    centres[indexI + oldSize] = coupleCentres[indexI];
                }
            }

            // Clear list
            faceParents.clear();
        }
    }
    else
    if (dimension == 3)
    {
        dynamicLabelList cellParents(10);

        forAll(procIndices_, pI)
        {
            // Fetch reference to subMesh
            const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

            // Prepare lists
            labelList coupleObjects;
            scalarField coupleWeights;
            vectorField coupleCentres;

            // Convex-set algorithm for cells
            cellSetAlgorithm cellAlgorithm
            (
                mesh,
                oldPoints_,
                edges_,
                faces_,
                cells_,
                owner_,
                neighbour_
            );

            // Initialize the bounding box
            cellAlgorithm.computeNormFactor(index);

            // Loop through all subMesh cells, and check for bounds
            const cellList& meshCells = mesh.cells();

            forAll(meshCells, cellI)
            {
                if (cellAlgorithm.contains(cellI))
                {
                    cellParents.append(cellI);
                }
            }

            // Obtain weighting factors for this cell.
            cellAlgorithm.computeWeights
            (
                index,
                0,
                cellParents,
                mesh.polyMesh::cellCells(),
                coupleObjects,
                coupleWeights,
                coupleCentres,
                output
            );

            // Add contributions with offsets
            if (coupleObjects.size())
            {
                label cellStart = cStarts[pI];

                // Resize lists
                label oldSize = parents.size();

                parents.setSize(oldSize + coupleObjects.size());
                weights.setSize(oldSize + coupleObjects.size());
                centres.setSize(oldSize + coupleObjects.size());

                forAll(coupleObjects, indexI)
                {
                    parents[indexI + oldSize] =
                    (
                        cellStart + coupleObjects[indexI]
                    );

                    weights[indexI + oldSize] = coupleWeights[indexI];
                    centres[indexI + oldSize] = coupleCentres[indexI];
                }
            }

            // Clear list
            cellParents.clear();
        }
    }
    else
    {
        FatalErrorIn("scalar dynamicTopoFvMesh::computeCoupledWeights()")
            << " Incorrect dimension: " << dimension << nl
            << abort(FatalError);
    }
}


// Fetch length-scale info for processor entities
scalar dynamicTopoFvMesh::processorLengthScale(const label index) const
{
    scalar procScale = 0.0;

    if (is2D())
    {
        // First check the master processor
        procScale += lengthScale_[owner_[index]];

        // Next, check the slave processor
        bool foundSlave = false;

        forAll(procIndices_, pI)
        {
            // Fetch non-const reference to subMeshes
            const label faceEnum = coupleMap::FACE;
            const coupleMap& cMap = recvMeshes_[pI].map();
            const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

            label sIndex = -1;

            if ((sIndex = cMap.findSlave(faceEnum, index)) > -1)
            {
                procScale += mesh.lengthScale_[mesh.owner_[sIndex]];

                foundSlave = true;
                break;
            }
        }

        // Should have found at least one slave
        if (!foundSlave)
        {
            FatalErrorIn
            (
                "scalar dynamicTopoFvMesh::processorLengthScale"
                "(const label index) const"
            )
                << "Processor lengthScale lookup failed: " << nl
                << " Master face: " << index
                << " :: " << faces_[index] << nl
                << abort(FatalError);
        }

        // Average the scale
        procScale *= 0.5;
    }
    else
    {
        // Check if this is a 'pure' processor edge
        bool pure = processorCoupledEntity(index, false, true, true);

        const label edgeEnum = coupleMap::EDGE;
        const labelList& eFaces = edgeFaces_[index];

        bool foundSlave = false;

        if (pure)
        {
            // First check the master processor
            label nC = 0;

            forAll(eFaces, faceI)
            {
                label own = owner_[eFaces[faceI]];
                label nei = neighbour_[eFaces[faceI]];

                procScale += lengthScale_[own];
                nC++;

                if (nei > -1)
                {
                    procScale += lengthScale_[nei];
                    nC++;
                }
            }

            // Next check slaves
            forAll(procIndices_, pI)
            {
                const coupleMap& cMap = recvMeshes_[pI].map();
                const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

                label sIndex = -1;

                if ((sIndex = cMap.findSlave(edgeEnum, index)) > -1)
                {
                    // Fetch connectivity from patchSubMesh
                    const labelList& peFaces = mesh.edgeFaces_[sIndex];

                    foundSlave = true;

                    forAll(peFaces, faceI)
                    {
                        label own = mesh.owner_[peFaces[faceI]];
                        label nei = mesh.neighbour_[peFaces[faceI]];

                        procScale += mesh.lengthScale_[own];
                        nC++;

                        if (nei > -1)
                        {
                            procScale += mesh.lengthScale_[nei];
                            nC++;
                        }
                    }
                }
            }

            // Average the final scale
            procScale /= nC;
        }
        else
        {
            // Processor is adjacent to physical patch types.
            // Search for boundary faces, and average their scale

            // First check the master processor
            label nBoundary = 0;

            forAll(eFaces, faceI)
            {
                label facePatch = whichPatch(eFaces[faceI]);

                // Skip internal faces
                if (facePatch == -1)
                {
                    continue;
                }

                // Skip processor patches
                if (getNeighbourProcessor(facePatch) > -1)
                {
                    continue;
                }

                // If this is a floating face, pick the owner length-scale
                if (lengthEstimator().isFreePatch(facePatch))
                {
                    procScale += lengthScale_[owner_[eFaces[faceI]]];
                }
                else
                {
                    // Fetch fixed length-scale
                    procScale +=
                    (
                        lengthEstimator().fixedLengthScale
                        (
                            eFaces[faceI],
                            facePatch
                        )
                    );
                }

                nBoundary++;
            }

            forAll(procIndices_, pI)
            {
                const coupleMap& cMap = recvMeshes_[pI].map();
                const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

                label sIndex = -1;

                if ((sIndex = cMap.findSlave(edgeEnum, index)) > -1)
                {
                    // Fetch connectivity from patchSubMesh
                    const labelList& peFaces = mesh.edgeFaces_[sIndex];

                    foundSlave = true;

                    forAll(peFaces, faceI)
                    {
                        label facePatch = mesh.whichPatch(peFaces[faceI]);

                        // Skip internal faces
                        if (facePatch == -1)
                        {
                            continue;
                        }

                        // Skip processor patches
                        if (mesh.getNeighbourProcessor(facePatch) > -1)
                        {
                            continue;
                        }

                        // If this is a floating face,
                        // pick the owner length-scale
                        if (lengthEstimator().isFreePatch(facePatch))
                        {
                            procScale +=
                            (
                                mesh.lengthScale_[mesh.owner_[peFaces[faceI]]]
                            );
                        }
                        else
                        {
                            // Fetch fixed length-scale
                            procScale +=
                            (
                                lengthEstimator().fixedLengthScale
                                (
                                    peFaces[faceI],
                                    facePatch
                                )
                            );
                        }

                        nBoundary++;
                    }
                }
            }

            if (nBoundary != 2)
            {
                // Write out for post-processing
                writeEdgeConnectivity(index);

                FatalErrorIn
                (
                    "scalar dynamicTopoFvMesh::processorLengthScale"
                    "(const label index) const"
                )
                    << " Expected two physical boundary patches: " << nl
                    << " nBoundary: " << nBoundary
                    << " Master edge: " << index
                    << " :: " << edges_[index]
                    << " Patch: "
                    << boundaryMesh()[whichEdgePatch(index)].name() << nl
                    << abort(FatalError);
            }

            procScale *= 0.5;
        }

        // Should have found at least one slave
        if (!foundSlave)
        {
            FatalErrorIn
            (
                "scalar dynamicTopoFvMesh::processorLengthScale"
                "(const label index) const"
            )
                << "Processor lengthScale lookup failed: " << nl
                << " Master edge: " << index
                << " :: " << edges_[index] << nl
                << abort(FatalError);
        }
    }

    return procScale;
}


// Method to determine whether the master face is locally coupled
bool dynamicTopoFvMesh::locallyCoupledEntity
(
    const label index,
    bool checkSlaves,
    bool checkProcs,
    bool checkFace
) const
{
    // Bail out if no patchCoupling is present
    if (patchCoupling_.empty())
    {
        return false;
    }

    if (is2D() || checkFace)
    {
        label patch = whichPatch(index);

        if (patch == -1)
        {
            return false;
        }

        // Processor checks receive priority.
        if (checkProcs)
        {
            if (getNeighbourProcessor(patch) > -1)
            {
                return false;
            }
        }

        // Check coupled master patches.
        if (patchCoupling_(patch))
        {
            return true;
        }
        else
        if (checkSlaves)
        {
            // Check on slave patches as well.
            forAll(patchCoupling_, pI)
            {
                if (patchCoupling_(pI))
                {
                    const coupleMap& cMap = patchCoupling_[pI].map();

                    if (cMap.slaveIndex() == patch)
                    {
                        return true;
                    }
                }
            }
        }
    }
    else
    {
        const labelList& eFaces = edgeFaces_[index];

        // Search for boundary faces, and determine boundary type.
        forAll(eFaces, faceI)
        {
            if (neighbour_[eFaces[faceI]] == -1)
            {
                label patch = whichPatch(eFaces[faceI]);

                // Processor checks receive priority.
                if (checkProcs)
                {
                    if (getNeighbourProcessor(patch) > -1)
                    {
                        return false;
                    }
                }

                // Check coupled master patches.
                if (patchCoupling_(patch))
                {
                    return true;
                }

                if (checkSlaves)
                {
                    // Check on slave patches as well.
                    forAll(patchCoupling_, pI)
                    {
                        if (patchCoupling_(pI))
                        {
                            const coupleMap& cMap =
                            (
                                patchCoupling_[pI].map()
                            );

                            if (cMap.slaveIndex() == patch)
                            {
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }

    // Could not find any faces on locally coupled patches.
    return false;
}


// Method to determine the locally coupled patch index
label dynamicTopoFvMesh::locallyCoupledEdgePatch(const label eIndex) const
{
    const labelList& eFaces = edgeFaces_[eIndex];

    // Search for boundary faces, and determine boundary type.
    forAll(eFaces, faceI)
    {
        if (neighbour_[eFaces[faceI]] == -1)
        {
            label patch = whichPatch(eFaces[faceI]);

            // Check coupled master patches.
            if (patchCoupling_(patch))
            {
                return patch;
            }

            // Check on slave patches as well.
            forAll(patchCoupling_, pI)
            {
                if (patchCoupling_(pI))
                {
                    const coupleMap& cMap = patchCoupling_[pI].map();

                    if (cMap.slaveIndex() == patch)
                    {
                        return patch;
                    }
                }
            }
        }
    }

    // Could not find any faces on locally coupled patches.
    FatalErrorIn
    (
        "label dynamicTopoFvMesh::locallyCoupledEdgePatch"
        "(const label cIndex) const"
    )
        << "Edge: " << eIndex << ":: " << edges_[eIndex]
        << " does not lie on any coupled patches."
        << abort(FatalError);

    return -1;
}


// Method to determine if the entity is on a processor boundary
//  - Also provides an additional check for 'pure' processor edges
//    i.e., edges that do not abut a physical patch. This is necessary
//    while deciding on collapse cases towards bounding curves.
bool dynamicTopoFvMesh::processorCoupledEntity
(
    const label index,
    bool checkFace,
    bool checkEdge,
    bool checkPure,
    FixedList<label, 2>* patchLabels,
    FixedList<vector, 2>* patchNormals
) const
{
    // Skip check for serial runs
    if (!Pstream::parRun())
    {
        return false;
    }

    label patch = -2;

    if ((is2D() || checkFace) && !checkEdge)
    {
        patch = whichPatch(index);

        if (patch == -1)
        {
            return false;
        }

        if (getNeighbourProcessor(patch) > -1)
        {
            return true;
        }
    }
    else
    {
        const labelList& eFaces = edgeFaces_[index];

        label nPhysical = 0, nProcessor = 0;

        // Search for boundary faces, and determine boundary type.
        forAll(eFaces, faceI)
        {
            label patch = whichPatch(eFaces[faceI]);

            if (patch == -1)
            {
                continue;
            }

            if (getNeighbourProcessor(patch) > -1)
            {
                // Increment the processor patch count
                nProcessor++;

                if (!checkPure)
                {
                    // We don't have to validate if this
                    // is a 'pure' processor edge, so bail out.
                    return true;
                }
            }
            else
            {
                // Physical patch.
                if (patchLabels)
                {
                    (*patchLabels)[nPhysical] = patch;
                }

                if (patchNormals)
                {
                    (*patchNormals)[nPhysical] =
                    (
                        faces_[eFaces[faceI]].normal(points_)
                    );
                }

                nPhysical++;
            }
        }

        if (patchLabels && patchNormals)
        {
            // Check other coupled-edges as well
            forAll(procIndices_, pI)
            {
                // Fetch reference to subMesh
                const coupleMap& cMap = recvMeshes_[pI].map();
                const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

                label sI = -1;

                if ((sI = cMap.findSlave(coupleMap::EDGE, index)) == -1)
                {
                    continue;
                }

                const labelList& seFaces = mesh.edgeFaces_[sI];

                forAll(seFaces, faceI)
                {
                    label sPatch = mesh.whichPatch(seFaces[faceI]);
                    label neiProc = mesh.getNeighbourProcessor(sPatch);

                    // Skip internal / processor faces
                    if (sPatch == -1 || neiProc > -1)
                    {
                        continue;
                    }

                    const face& sFace = mesh.faces_[seFaces[faceI]];

                    // Fill patch / normal info
                    (*patchLabels)[nPhysical] = sPatch;
                    (*patchNormals)[nPhysical] = sFace.normal(mesh.points_);

                    nPhysical++;
                }
            }
        }

        // Purity check
        if (checkPure)
        {
            if (nProcessor && !nPhysical)
            {
                return true;
            }
        }
    }

    // Could not find any faces on processor patches.
    return false;
}


// Build a list of entities that need to be avoided
// by regular topo-changes.
void dynamicTopoFvMesh::buildEntitiesToAvoid
(
    labelHashSet& entities,
    bool checkSubMesh
)
{
    entities.clear();

    // Build a set of entities to avoid during regular modifications,
    // and build a master stack for coupled modifications.

    // Determine locally coupled slave patches.
    labelHashSet localMasterPatches, localSlavePatches;

    forAll(patchCoupling_, patchI)
    {
        if (patchCoupling_(patchI))
        {
            const coupleMap& cMap = patchCoupling_[patchI].map();

            localMasterPatches.insert(cMap.masterIndex());
            localSlavePatches.insert(cMap.slaveIndex());
        }
    }

    // Loop through boundary faces and check whether
    // they belong to master/slave coupled patches.
    for (label faceI = nOldInternalFaces_; faceI < faces_.size(); faceI++)
    {
        // Add only valid faces
        if (faces_[faceI].empty())
        {
            continue;
        }

        label pIndex = whichPatch(faceI);

        if (pIndex == -1)
        {
            continue;
        }

        // Check if this is a coupled face.
        if
        (
            localMasterPatches.found(pIndex) ||
            localSlavePatches.found(pIndex) ||
            getNeighbourProcessor(pIndex) > -1
        )
        {
            if (is2D())
            {
                // Avoid this face during regular modification.
                if (!entities.found(faceI))
                {
                    entities.insert(faceI);
                }
            }
            else
            {
                const labelList& fEdges = faceEdges_[faceI];

                forAll(fEdges, edgeI)
                {
                    // Avoid this edge during regular modification.
                    if (!entities.found(fEdges[edgeI]))
                    {
                        entities.insert(fEdges[edgeI]);
                    }
                }
            }
        }
    }

    // Loop through entities contained in patchSubMeshes, if requested
    if (checkSubMesh)
    {
        forAll(procIndices_, pI)
        {
            const coupleMap& cMap = sendMeshes_[pI].map();
            const Map<label> rEdgeMap = cMap.reverseEntityMap(coupleMap::EDGE);

            if (cMap.slaveIndex() == Pstream::myProcNo())
            {
                forAllConstIter(Map<label>, rEdgeMap, eIter)
                {
                    if (is2D())
                    {
                        const labelList& eFaces = edgeFaces_[eIter.key()];

                        forAll(eFaces, faceI)
                        {
                            if (faces_[eFaces[faceI]].size() == 4)
                            {
                                if (!entities.found(eFaces[faceI]))
                                {
                                    entities.insert(eFaces[faceI]);
                                }
                            }
                        }
                    }
                    else
                    {
                        const edge& check = edges_[eIter.key()];
                        const labelList& eFaces = edgeFaces_[eIter.key()];

                        // Skip deleted edges
                        if (eFaces.size())
                        {
                            forAll(check, pI)
                            {
                                const labelList& pE = pointEdges_[check[pI]];

                                forAll(pE, edgeI)
                                {
                                    if (!entities.found(pE[edgeI]))
                                    {
                                        entities.insert(pE[edgeI]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (debug > 3)
    {
        Pout<< nl << "nEntitiesToAvoid: " << entities.size() << endl;

        if (debug > 4)
        {
            // Write out entities
            label elemType = is2D() ? 2 : 1;

            writeVTK
            (
                "entitiesToAvoid_"
              + Foam::name(Pstream::myProcNo()),
                entities.toc(),
                elemType
            );
        }
    }
}


// Check whether the specified edge is a coupled master edge.
bool dynamicTopoFvMesh::isCoupledMaster
(
    const label eIndex
) const
{
    if (!coupledModification_)
    {
        return true;
    }

    return locallyCoupledEntity(eIndex);
}


} // End namespace Foam

// ************************************************************************* //

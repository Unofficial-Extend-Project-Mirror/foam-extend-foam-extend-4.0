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
    coupleMap

Description
    Implementation of the coupleMap class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "coupleMap.H"
#include "boolList.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* coupleMap::names[coupleMap::INVALID + 1] =
{
    "BISECTION",
    "COLLAPSE_FIRST",
    "COLLAPSE_SECOND",
    "COLLAPSE_MIDPOINT",
    "REMOVE_CELL",
    "MOVE_POINT",
    "CONVERT_PATCH",
    "CONVERT_PHYSICAL",
    "INVALID"
};

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

const char* coupleMap::asText(const opType oType)
{
    return names[oType];
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupleMap::coupleMap
(
    const IOobject& io,
    const bool twoDMesh,
    const bool isLocal,
    const bool isSend,
    const label patchIndex,
    const label masterIndex,
    const label slaveIndex
)
:
    regIOobject(io),
    twoDMesh_(twoDMesh),
    isLocal_(isLocal),
    isSend_(isSend),
    patchIndex_(patchIndex),
    masterIndex_(masterIndex),
    slaveIndex_(slaveIndex),
    nEntities_(-1),
    edgesPtr_(NULL),
    facesPtr_(NULL),
    faceEdgesPtr_(NULL)
{
    if
    (
        (io.readOpt() == IOobject::MUST_READ)
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        // Construct an Istream and read from disk.
        readData(readStream(typeName));
        close();
    }
}


// Construct as copy
coupleMap::coupleMap(const coupleMap& cm)
:
    regIOobject(cm, true),
    twoDMesh_(cm.twoDMesh_),
    isLocal_(cm.isLocal_),
    isSend_(cm.isSend_),
    patchIndex_(cm.patchIndex_),
    masterIndex_(cm.masterIndex_),
    slaveIndex_(cm.slaveIndex_),
    nEntities_(cm.nEntities_),
    edgesPtr_(NULL),
    facesPtr_(NULL),
    faceEdgesPtr_(NULL)
{
    if
    (
        (cm.readOpt() == IOobject::MUST_READ)
     || (cm.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        // Construct an Istream and read from disk.
        readData(readStream(typeName));
        close();
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

coupleMap::~coupleMap()
{
    clearMaps();
    clearBuffers();
    clearAddressing();
}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void coupleMap::clearAddressing() const
{
    deleteDemandDrivenData(edgesPtr_);
    deleteDemandDrivenData(facesPtr_);
    deleteDemandDrivenData(faceEdgesPtr_);
}


void coupleMap::makeEdges() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (edgesPtr_)
    {
        FatalErrorIn("void coupleMap::makeEdges() const")
            << "Edges have already been calculated."
            << abort(FatalError);
    }

    label nEdges = nEntities(coupleMap::EDGE);
    const labelList& eBuffer = entityBuffer(coupleMap::EDGE);

    edgesPtr_ = new edgeList(nEdges);

    edgeList& edges = *edgesPtr_;

    forAll(edges, edgeI)
    {
        edges[edgeI][0] = eBuffer[(2*edgeI)+0];
        edges[edgeI][1] = eBuffer[(2*edgeI)+1];
    }
}


void coupleMap::makeFaces() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (facesPtr_ || faceEdgesPtr_)
    {
        FatalErrorIn("void coupleMap::makeFaces() const")
            << "Faces have already been calculated."
            << abort(FatalError);
    }

    label nFaces = nEntities(coupleMap::FACE);

    const labelList& fBuffer = entityBuffer(coupleMap::FACE);
    const labelList& feBuffer = entityBuffer(coupleMap::FACE_EDGE);
    const labelList& nfeBuffer = entityBuffer(coupleMap::NFE_BUFFER);

    facesPtr_ = new faceList(nFaces);
    faceEdgesPtr_ = new labelListList(nFaces);

    faceList& faces = *facesPtr_;
    labelListList& faceEdges = *faceEdgesPtr_;

    label sumNFE = 0;

    forAll(faces, faceI)
    {
        face& f = faces[faceI];
        labelList& fe = faceEdges[faceI];

        // Fetch the buffer value
        label nfe = nfeBuffer[faceI];

        // Size up the lists
        f.setSize(nfe, -1);
        fe.setSize(nfe, -1);

        for (label p = 0; p < nfe; p++)
        {
            f[p] = fBuffer[sumNFE + p];
            fe[p] = feBuffer[sumNFE + p];
        }

        sumNFE += nfe;
    }

    if (sumNFE != nEntities(coupleMap::NFE_SIZE))
    {
        FatalErrorIn("void coupleMap::makeFaces() const")
            << " Mismatched buffer." << nl
            << " sumNFE: " << sumNFE << nl
            << " NFE_SIZE: " << nEntities(coupleMap::NFE_SIZE) << nl
            << abort(FatalError);
    }
}


void coupleMap::makeFaceMap() const
{
    // It is an error to attempt to recalculate
    // if the map is already calculated
    if (faceMap_.size())
    {
        FatalErrorIn("void coupleMap::makeFaceMap() const")
            << "faceMap has already been calculated."
            << abort(FatalError);
    }

    const Map<label>& mapFaces = entityMap(coupleMap::FACE);

    // Size the list
    faceMap_.setSize(mapFaces.size(), -1);

    // Fill-in entries
    forAllConstIter(Map<label>, mapFaces, fIter)
    {
        faceMap_[fIter.key()] = fIter();
    }
}


void coupleMap::makeCellMap() const
{
    // It is an error to attempt to recalculate
    // if the map is already calculated
    if (cellMap_.size())
    {
        FatalErrorIn("void coupleMap::makeCellMap() const")
            << "cellMap has already been calculated."
            << abort(FatalError);
    }

    const Map<label>& mapCells = entityMap(coupleMap::CELL);

    // Size the list
    cellMap_.setSize(mapCells.size(), -1);

    // Fill-in entries
    forAllConstIter(Map<label>, mapCells, cIter)
    {
        cellMap_[cIter.key()] = cIter();
    }
}


void coupleMap::makeInternalFaceMap() const
{
    // It is an error to attempt to recalculate
    // if the map is already calculated
    if (internalFaceMap_.size())
    {
        FatalErrorIn("void coupleMap::makeInternalFaceMap() const")
            << "internal faceMap has already been calculated."
            << abort(FatalError);
    }

    // Slice for internal faces
    internalFaceMap_ = SubList<label>(faceMap(), nEntities(INTERNAL_FACE));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

pointField& coupleMap::pointBuffer() const
{
    return pointBuffer_;
}


pointField& coupleMap::oldPointBuffer() const
{
    return oldPointBuffer_;
}


labelList& coupleMap::subMeshPoints() const
{
    return subMeshPoints_;
}


List<labelPair>& coupleMap::globalProcPoints() const
{
    return globalProcPoints_;
}


void coupleMap::allocateBuffers() const
{
    forAll(nEntities_, entityI)
    {
        if (nEntities_[entityI] < 0)
        {
            FatalErrorIn("void coupleMap::allocateBuffers() const")
                << " Entity sizes are not valid." << nl
                << " nEntities: " << nEntities_
                << abort(FatalError);
        }
    }

    // Size up point buffers
    pointBuffer().setSize(nEntities(POINT));
    oldPointBuffer().setSize(nEntities(POINT));

    // Size up connectivity buffers
    entityBuffer(POINT).setSize
    (
        nEntities(GLOBAL_POINT)
      - nEntities(SHARED_POINT)
    );

    // Set edge buffer
    entityBuffer(EDGE).setSize(2 * nEntities(EDGE));

    // Set addressing buffers
    entityBuffer(OWNER).setSize(nEntities(FACE));
    entityBuffer(NEIGHBOUR).setSize(nEntities(INTERNAL_FACE));

    // Size up boundary buffers
    entityBuffer(FACE_STARTS).setSize(nEntities(NBDY));
    entityBuffer(FACE_SIZES).setSize(nEntities(NBDY));
    entityBuffer(EDGE_STARTS).setSize(nEntities(NBDY));
    entityBuffer(EDGE_SIZES).setSize(nEntities(NBDY));
    entityBuffer(PATCH_ID).setSize(nEntities(NBDY));

    // Set face-sizes
    entityBuffer(NFE_BUFFER).setSize(nEntities(FACE));

    // Allocate for variable size face-lists
    entityBuffer(FACE).setSize(nEntities(NFE_SIZE));
    entityBuffer(FACE_EDGE).setSize(nEntities(NFE_SIZE));
}


label coupleMap::findSlave
(
    const label eType,
    const label Index
) const
{
    Map<label>::const_iterator it = entityMap_[eType].find(Index);

    if (it == entityMap_[eType].end())
    {
        return -1;
    }
    else
    {
        return it();
    }
}


label coupleMap::findMaster
(
    const label eType,
    const label Index
) const
{
    Map<label>::const_iterator it = reverseEntityMap_[eType].find(Index);

    if (it == reverseEntityMap_[eType].end())
    {
        return -1;
    }
    else
    {
        return it();
    }
}


void coupleMap::removeSlave
(
    const label eType,
    const label Index
) const
{
    Map<label>::iterator it = reverseEntityMap_[eType].find(Index);

    if (it != reverseEntityMap_[eType].end())
    {
        reverseEntityMap_[eType].erase(it);
    }
}


void coupleMap::removeMaster
(
    const label eType,
    const label Index
) const
{
    Map<label>::iterator it = entityMap_[eType].find(Index);

    if (it != entityMap_[eType].end())
    {
        entityMap_[eType].erase(it);
    }
}


void coupleMap::mapSlave
(
    const label eType,
    const label master,
    const label slave
) const
{
    entityMap_[eType].set(master, slave);
}


void coupleMap::mapMaster
(
    const label eType,
    const label slave,
    const label master
) const
{
    reverseEntityMap_[eType].set(slave, master);
}


void coupleMap::pushOperation
(
    const label index,
    const opType oType
) const
{
    entityIndices_.setSize(entityIndices_.size() + 1, index);
    entityOperations_.setSize(entityOperations_.size() + 1, oType);
}


void coupleMap::pushOperation
(
    const label index,
    const opType oType,
    const label pIndex
) const
{
    if (oType == coupleMap::CONVERT_PHYSICAL)
    {
        entityIndices_.setSize(entityIndices_.size() + 1, index);
        entityOperations_.setSize(entityOperations_.size() + 1, oType);

        patchIndices_.setSize(patchIndices_.size() + 1, pIndex);
    }
    else
    {
        opType t = (oType < coupleMap::INVALID ? oType : coupleMap::INVALID);

        FatalErrorIn("void coupleMap::pushOperation() const")
            << " Expected CONVERT_PHYSICAL" << nl
            << " Found: " << coupleMap::names[t]
            << abort(FatalError);
    }
}


void coupleMap::pushOperation
(
    const label index,
    const opType oType,
    const point& newPoint,
    const point& oldPoint
) const
{
    entityIndices_.setSize(entityIndices_.size() + 1, index);
    entityOperations_.setSize(entityOperations_.size() + 1, oType);

    if
    (
        oType == coupleMap::MOVE_POINT ||
        oType == coupleMap::CONVERT_PATCH
    )
    {
        moveNewPoints_.setSize(moveNewPoints_.size() + 1, newPoint);
        moveOldPoints_.setSize(moveOldPoints_.size() + 1, oldPoint);
    }
    else
    {
        opType t = (oType < coupleMap::INVALID ? oType : coupleMap::INVALID);

        FatalErrorIn("void coupleMap::pushOperation() const")
            << " Expected either MOVE_POINT or CONVERT_PATCH" << nl
            << " Found: " << coupleMap::names[t]
            << abort(FatalError);
    }
}


void coupleMap::transferMaps
(
    const label eType,
    Map<label>& newEntityMap,
    Map<label>& newReverseEntityMap
) const
{
    entityMap_[eType].transfer(newEntityMap);
    reverseEntityMap_[eType].transfer(newReverseEntityMap);
}


void coupleMap::clearMaps() const
{
    faceMap_.clear();
    cellMap_.clear();
    internalFaceMap_.clear();

    subMeshPointMap_.clear();
    subMeshEdgeMap_.clear();

    forAll(entityMap_, mapI)
    {
        entityMap_[mapI].clear();
        reverseEntityMap_[mapI].clear();
    }
}


void coupleMap::clearBuffers() const
{
    pointBuffer_.clear();
    oldPointBuffer_.clear();

    subMeshPoints_.clear();
    globalProcPoints_.clear();

    forAll(entityBuffer_, bufferI)
    {
        entityBuffer_[bufferI].clear();
    }

    entityIndices_.clear();
    entityOperations_.clear();

    patchIndices_.clear();
    moveNewPoints_.clear();
    moveOldPoints_.clear();
}


label coupleMap::nInternalFaces() const
{
    return nEntities(INTERNAL_FACE);
}


const labelList& coupleMap::owner() const
{
    return entityBuffer(OWNER);
}


const labelList& coupleMap::neighbour() const
{
    return entityBuffer(NEIGHBOUR);
}


const edgeList& coupleMap::edges() const
{
    if (!edgesPtr_)
    {
        makeEdges();
    }

    return *edgesPtr_;
}


const faceList& coupleMap::faces() const
{
    if (!facesPtr_)
    {
        makeFaces();
    }

    return *facesPtr_;
}


const labelListList& coupleMap::faceEdges() const
{
    if (!faceEdgesPtr_)
    {
        makeFaces();
    }

    return *faceEdgesPtr_;
}


const labelList& coupleMap::faceMap() const
{
    if (faceMap_.empty())
    {
        makeFaceMap();
    }

    return faceMap_;
}


const labelList& coupleMap::cellMap() const
{
    if (cellMap_.empty())
    {
        makeCellMap();
    }

    return cellMap_;
}


const labelList& coupleMap::internalFaceMap() const
{
    if (internalFaceMap_.empty())
    {
        makeInternalFaceMap();
    }

    return internalFaceMap_;
}


bool coupleMap::readData(Istream& is)
{
    Map<label> tmpMap(is);

    entityMap(coupleMap::POINT).transfer(tmpMap);

    // Prepare the reversePointMap as well.
    const Map<label>& pMap = entityMap(coupleMap::POINT);
    Map<label>& rpMap = reverseEntityMap(coupleMap::POINT);

    forAllConstIter(Map<label>, pMap, pIter)
    {
        rpMap.set(pIter(), pIter.key());
    }

    return !is.bad();
}


bool coupleMap::writeData(Ostream& os) const
{
    // Only write-out point-map information
    // to avoid geometric checking.
    // The rest can be constructed topologically.
    return (os << entityMap(coupleMap::POINT)).good();;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void coupleMap::operator=(const coupleMap& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("coupleMap::operator=(const coupleMap&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


} // End namespace Foam

// ************************************************************************* //

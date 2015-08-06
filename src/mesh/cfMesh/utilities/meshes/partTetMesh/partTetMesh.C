/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "demandDrivenData.H"
#include "partTetMesh.H"
#include "polyMeshGenModifier.H"
#include "VRWGraphList.H"
#include "polyMeshGenAddressing.H"
#include "helperFunctions.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

partTetMesh::partTetMesh(polyMeshGen& mesh, const labelLongList& lockedPoints)
:
    origMesh_(mesh),
    points_(),
    tets_(),
    nodeLabelInOrigMesh_(),
    smoothVertex_(),
    pointTets_(),
    internalPointsOrderPtr_(NULL),
    boundaryPointsOrderPtr_(NULL),
    globalPointLabelPtr_(NULL),
    pAtProcsPtr_(NULL),
    globalToLocalPointAddressingPtr_(NULL),
    neiProcsPtr_(NULL),
    pAtParallelBoundariesPtr_(NULL),
    pAtBufferLayersPtr_(NULL)
{
    List<direction> useCell(mesh.cells().size(), direction(1));

    boolList lockedPoint(mesh.points().size(), false);
    forAll(lockedPoints, i)
        lockedPoint[lockedPoints[i]] = true;

    createPointsAndTets(useCell, lockedPoint);
}

partTetMesh::partTetMesh
(
    polyMeshGen& mesh,
    const labelLongList& lockedPoints,
    const direction nLayers
)
:
    origMesh_(mesh),
    points_(),
    tets_(),
    nodeLabelInOrigMesh_(),
    smoothVertex_(),
    pointTets_(),
    internalPointsOrderPtr_(NULL),
    boundaryPointsOrderPtr_(NULL),
    globalPointLabelPtr_(NULL),
    pAtProcsPtr_(NULL),
    globalToLocalPointAddressingPtr_(NULL),
    neiProcsPtr_(NULL),
    pAtParallelBoundariesPtr_(NULL),
    pAtBufferLayersPtr_(NULL)
{
    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();
    const labelList& owner = mesh.owner();
    const PtrList<boundaryPatch>& boundaries = mesh.boundaries();
    const VRWGraph& pointCells = mesh.addressingData().pointCells();

    List<direction> useCell(cells.size(), direction(0));

    //- select cells containing at least one vertex of the bad faces
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label size = boundaries[patchI].patchSize();

        for(label fI=0;fI<size;++fI)
        {
            useCell[owner[start+fI]] = 1;
        }
    }

    //- add additional layer of cells
    for(direction layerI=1;layerI<(nLayers+1);++layerI)
    {
        forAll(useCell, cI)
            if( useCell[cI] == layerI )
            {
                const cell& c = cells[cI];

                forAll(c, fI)
                {
                    const face& f = faces[c[fI]];

                    forAll(f, pI)
                    {
                        forAllRow(pointCells, f[pI], pcI)
                        {
                            const label cLabel = pointCells(f[pI], pcI);
                            if( !useCell[cLabel] )
                                useCell[cLabel] = layerI + 1;
                        }
                    }
                }
            }

        if( Pstream::parRun() )
        {
            const labelLongList& globalPointLabel =
                mesh.addressingData().globalPointLabel();
            const VRWGraph& pProcs = mesh.addressingData().pointAtProcs();
            const Map<label>& globalToLocal =
                mesh.addressingData().globalToLocalPointAddressing();

            std::map<label, LongList<label> > eData;
            forAllConstIter(Map<label>, globalToLocal, iter)
            {
                const label pointI = iter();

                forAllRow(pProcs, pointI, procI)
                {
                    const label neiProc = pProcs(pointI, procI);
                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    if( eData.find(neiProc) == eData.end() )
                    {
                        eData.insert
                        (
                            std::make_pair(neiProc, LongList<label>())
                        );
                    }

                    forAllRow(pointCells, pointI, pcI)
                        if( useCell[pointCells(pointI, pcI)] == layerI )
                        {
                            eData[neiProc].append(globalPointLabel[pointI]);
                            break;
                        }
                }
            }

            //- exchange data with other processors
            labelLongList receivedData;
            help::exchangeMap(eData, receivedData);

            forAll(receivedData, i)
            {
                const label pointI = globalToLocal[receivedData[i]];

                forAllRow(pointCells, pointI, pcI)
                {
                    const label cLabel = pointCells(pointI, pcI);
                    if( !useCell[cLabel] )
                        useCell[cLabel] = layerI + 1;
                }
            }
        }
    }

    boolList lockedPoint(mesh.points().size(), false);
    forAll(lockedPoints, i)
        lockedPoint[lockedPoints[i]] = true;

    createPointsAndTets(useCell, lockedPoint);
}

partTetMesh::partTetMesh
(
    polyMeshGen& mesh,
    const labelLongList& lockedPoints,
    labelHashSet& badFaces,
    const direction additionalLayers
)
:
    origMesh_(mesh),
    points_(),
    tets_(),
    nodeLabelInOrigMesh_(),
    smoothVertex_(),
    pointTets_(),
    internalPointsOrderPtr_(NULL),
    boundaryPointsOrderPtr_(NULL),
    globalPointLabelPtr_(NULL),
    pAtProcsPtr_(NULL),
    globalToLocalPointAddressingPtr_(NULL),
    neiProcsPtr_(NULL),
    pAtParallelBoundariesPtr_(NULL),
    pAtBufferLayersPtr_(NULL)
{
    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();
    const VRWGraph& pointCells = mesh.addressingData().pointCells();

    List<direction> useCell(cells.size(), direction(0));

    //- select cells containing at least one vertex of the bad faces
    forAll(faces, faceI)
        if( badFaces.found(faceI) )
        {
            const face& f = faces[faceI];

            forAll(f, pI)
            {
                forAllRow(pointCells, f[pI], pcI)
                    useCell[pointCells(f[pI], pcI)] = 1;
            }
        }

    //- add additional layer of cells
    for(direction layerI=1;layerI<(additionalLayers+1);++layerI)
    {
        forAll(useCell, cI)
            if( useCell[cI] == layerI )
            {
                const cell& c = cells[cI];

                forAll(c, fI)
                {
                    const face& f = faces[c[fI]];

                    forAll(f, pI)
                    {
                        forAllRow(pointCells, f[pI], pcI)
                        {
                            const label cLabel = pointCells(f[pI], pcI);
                            if( !useCell[cLabel] )
                                useCell[cLabel] = layerI + 1;
                        }
                    }
                }
            }

        if( Pstream::parRun() )
        {
            const labelLongList& globalPointLabel =
                mesh.addressingData().globalPointLabel();
            const VRWGraph& pProcs = mesh.addressingData().pointAtProcs();
            const Map<label>& globalToLocal =
                mesh.addressingData().globalToLocalPointAddressing();

            std::map<label, LongList<label> > eData;
            forAllConstIter(Map<label>, globalToLocal, iter)
            {
                const label pointI = iter();

                forAllRow(pProcs, pointI, procI)
                {
                    const label neiProc = pProcs(pointI, procI);
                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    if( eData.find(neiProc) == eData.end() )
                    {
                        eData.insert
                        (
                            std::make_pair(neiProc, LongList<label>())
                        );
                    }

                    forAllRow(pointCells, pointI, pcI)
                        if( useCell[pointCells(pointI, pcI)] == layerI )
                        {
                            eData[neiProc].append(globalPointLabel[pointI]);
                            break;
                        }
                }
            }

            //- exchange data with other processors
            labelLongList receivedData;
            help::exchangeMap(eData, receivedData);

            forAll(receivedData, i)
            {
                const label pointI = globalToLocal[receivedData[i]];

                forAllRow(pointCells, pointI, pcI)
                {
                    const label cLabel = pointCells(pointI, pcI);
                    if( !useCell[cLabel] )
                        useCell[cLabel] = layerI + 1;
                }
            }
        }
    }

    boolList lockedPoint(mesh.points().size(), false);
    forAll(lockedPoints, i)
        lockedPoint[lockedPoints[i]] = true;

    createPointsAndTets(useCell, lockedPoint);
}

partTetMesh::~partTetMesh()
{
    deleteDemandDrivenData(internalPointsOrderPtr_);
    deleteDemandDrivenData(boundaryPointsOrderPtr_);
    deleteDemandDrivenData(globalPointLabelPtr_);
    deleteDemandDrivenData(pAtProcsPtr_);
    deleteDemandDrivenData(globalToLocalPointAddressingPtr_);
    deleteDemandDrivenData(neiProcsPtr_);
    deleteDemandDrivenData(pAtParallelBoundariesPtr_);
    deleteDemandDrivenData(pAtBufferLayersPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const VRWGraph& partTetMesh::internalPointOrdering() const
{
    # ifdef USE_OMP
    if( omp_in_parallel() )
    {
        FatalErrorIn
        (
            "const VRWGraph& partTetMesh::internalPointOrdering() const"
        ) << "Calculating addressing inside a parallel region."
          << " This is not thread safe" << exit(FatalError);
    }
    # endif

    if( !internalPointsOrderPtr_ )
        createSMOOTHPointsOrdering();

    return *internalPointsOrderPtr_;
}

const VRWGraph& partTetMesh::boundaryPointOrdering() const
{
    # ifdef USE_OMP
    if( omp_in_parallel() )
    {
        FatalErrorIn
        (
            "const VRWGraph& partTetMesh::boundaryPointOrdering() const"
        ) << "Calculating addressing inside a parallel region."
          << " This is not thread safe" << exit(FatalError);
    }
    # endif

    if( !boundaryPointsOrderPtr_ )
        createBOUNDARYPointsOrdering();

    return *boundaryPointsOrderPtr_;
}
void partTetMesh::updateVertex(const label pointI, const point& newP)
{
    points_[pointI] = newP;

    if( smoothVertex_[pointI] & (FACECENTRE+CELLCENTRE) )
    {
        Warning << "Smoothing auxiliary vertex."
            << " This has no effect on the original mesh" << endl;
        return;
    }

    //- find face centres attached
    DynList<label, 64> helper;
    forAllRow(pointTets_, pointI, ptI)
    {
        const label centreI = tets_[pointTets_(pointI, ptI)][2];
        if( smoothVertex_[centreI] & FACECENTRE )
            helper.appendIfNotIn(centreI);
    }

    //- update coordinates of FACECENTRE vertices
    forAll(helper, i)
    {
        const label centreI = helper[i];

        point centre(vector::zero);
        scalar faceArea(0.0);
        forAllRow(pointTets_, centreI, ptI)
        {
            const partTet& tet = tets_[pointTets_(centreI, ptI)];
            point c(vector::zero);
            for(label i=0;i<3;++i)
                c += points_[tet[i]];
            c /= 3;
            const scalar area = Foam::mag(tet.Sd(points_)) + VSMALL;

            centre += c * area;
            faceArea += area;
        }

        points_[centreI] = centre / faceArea;
    }

    //- find cell centres attached
    helper.clear();
    forAllRow(pointTets_, pointI, ptI)
    {
        const label centreI = tets_[pointTets_(pointI, ptI)][3];
        if( smoothVertex_[centreI] & CELLCENTRE )
            helper.appendIfNotIn(centreI);
    }

    //- update coordinates of CELLCENTRE vertices
    forAll(helper, i)
    {
        const label centreI = helper[i];

        point centre(vector::zero);
        scalar cellVol(0.0);
        forAllRow(pointTets_, centreI, ptI)
        {
            const partTet& tet = tets_[pointTets_(centreI, ptI)];
            const point c = tet.centroid(points_);
            const scalar vol = Foam::mag(tet.mag(points_)) + VSMALL;

            centre += c * vol;
            cellVol += vol;
        }

        points_[centreI] = centre / cellVol;
    }
}

void partTetMesh::updateVerticesSMP(const List<LongList<labelledPoint> >& np)
{
    List<direction> updateType(points_.size(), direction(0));

    # ifdef USE_OMP
    # pragma omp parallel num_threads(np.size())
    # endif
    {
        # ifdef USE_OMP
        const LongList<labelledPoint>& newPoints = np[omp_get_thread_num()];
        # else
        const LongList<labelledPoint>& newPoints = np[0];
        # endif

        forAll(newPoints, i)
        {
            const labelledPoint& lp = newPoints[i];
            const label pointI = lp.pointLabel();

            points_[pointI] = lp.coordinates();

            forAllRow(pointTets_, pointI, ptI)
            {
                const partTet& pt = tets_[pointTets_(pointI, ptI)];

                if( smoothVertex_[pt[3]] & CELLCENTRE )
                    updateType[pt[3]] |= CELLCENTRE;
                if( smoothVertex_[pt[2]] & FACECENTRE )
                    updateType[pt[2]] |= FACECENTRE;
            }
        }

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp flush(updateType)

        //- update coordinates of CELLCENTRE and FACECENTRE vertices
        # pragma omp for schedule(dynamic, 20)
        # endif
        forAll(updateType, pI)
        {
            if( updateType[pI] & CELLCENTRE )
            {
                point centre(vector::zero);
                scalar cellVol(0.0);
                forAllRow(pointTets_, pI, ptI)
                {
                    const partTet& tet = tets_[pointTets_(pI, ptI)];
                    const point c = tet.centroid(points_);
                    const scalar vol = Foam::mag(tet.mag(points_)) + VSMALL;

                    centre += c * vol;
                    cellVol += vol;
                }

                points_[pI] = centre / cellVol;
            }
            else if( updateType[pI] & FACECENTRE )
            {
                point centre(vector::zero);
                scalar faceArea(0.0);
                forAllRow(pointTets_, pI, ptI)
                {
                    const partTet& tet = tets_[pointTets_(pI, ptI)];
                    point c(vector::zero);
                    for(label i=0;i<3;++i)
                        c += points_[tet[i]];
                    c /= 3;
                    const scalar area = Foam::mag(tet.Sd(points_)) + VSMALL;

                    centre += c * area;
                    faceArea += area;
                }

                points_[pI] = centre / faceArea;
            }
        }
    }
}

void partTetMesh::updateOrigMesh(boolList* changedFacePtr)
{
    pointFieldPMG& pts = origMesh_.points();

    boolList changedNode(pts.size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for if( pts.size() > 1000 ) \
    schedule(guided, 10)
    # endif
    forAll(nodeLabelInOrigMesh_, pI)
        if( nodeLabelInOrigMesh_[pI] != -1 )
        {
            changedNode[nodeLabelInOrigMesh_[pI]] = true;
            pts[nodeLabelInOrigMesh_[pI]] = points_[pI];
        }

    if( changedFacePtr )
    {
        boolList& chF = *changedFacePtr;
        chF = false;

        const cellListPMG& cells = origMesh_.cells();
        const VRWGraph& pointCells = origMesh_.addressingData().pointCells();

        # ifdef USE_OMP
        # pragma omp parallel for if( pointCells.size() > 100 ) \
        schedule(dynamic, 20)
        # endif
        forAll(pointCells, pointI)
        {
            if( changedNode[pointI] )
            {
                forAllRow(pointCells, pointI, pcI)
                {
                    const cell& c = cells[pointCells(pointI, pcI)];

                    forAll(c, fI)
                        chF[c[fI]] = true;
                }
            }
        }

        //- make sure that neighbouring processors get the same information
        const PtrList<processorBoundaryPatch>& pBnd = origMesh_.procBoundaries();
        forAll(pBnd, patchI)
        {
            const label start = pBnd[patchI].patchStart();
            const label size = pBnd[patchI].patchSize();

            labelLongList sendData;
            for(label faceI=0;faceI<size;++faceI)
            {
                if( chF[start+faceI] )
                    sendData.append(faceI);
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                pBnd[patchI].neiProcNo(),
                sendData.byteSize()
            );

            toOtherProc << sendData;
        }

        forAll(pBnd, patchI)
        {
            labelList receivedData;

            IPstream fromOtherProc
            (
                Pstream::blocking,
                pBnd[patchI].neiProcNo()
            );

            fromOtherProc >> receivedData;

            const label start = pBnd[patchI].patchStart();
            forAll(receivedData, i)
                chF[start+receivedData[i]] = true;
        }

        //- update geometry information
        const_cast<polyMeshGenAddressing&>
        (
            origMesh_.addressingData()
        ).updateGeometry(chF);
    }
    else
    {
        const_cast<polyMeshGenAddressing&>
        (
            origMesh_.addressingData()
        ).clearGeom();
    }
}

void partTetMesh::createPolyMesh(polyMeshGen& pmg) const
{
    polyMeshGenModifier meshModifier(pmg);

    pointFieldPMG& pAccess = meshModifier.pointsAccess();
    pAccess.setSize(points_.size());
    forAll(points_, pI)
        pAccess[pI] = points_[pI];

    VRWGraphList cellFaces;

    forAll(tets_, tetI)
    {
        const partTet& tet = tets_[tetI];

        FixedList<FixedList<label, 3>, 4> tetFaces;

        //- face 0
        tetFaces[0][0] = tet[0];
        tetFaces[0][1] = tet[2];
        tetFaces[0][2] = tet[1];

        //- face 1
        tetFaces[1][0] = tet[0];
        tetFaces[1][1] = tet[1];
        tetFaces[1][2] = tet[3];

        //- face 2
        tetFaces[2][0] = tet[0];
        tetFaces[2][1] = tet[3];
        tetFaces[2][2] = tet[2];

        //- face 3
        tetFaces[3][0] = tet[1];
        tetFaces[3][1] = tet[2];
        tetFaces[3][2] = tet[3];

        cellFaces.appendGraph(tetFaces);
    }

    meshModifier.addCells(cellFaces);
    meshModifier.reorderBoundaryFaces();

    //- store points into subsets
    const label bndPointID = pmg.addPointSubset("boundaryPoints");
    const label smoothPointID = pmg.addPointSubset("smoothPoints");
    const label faceCentreID = pmg.addPointSubset("faceCentres");
    const label cellCentreID = pmg.addPointSubset("cellCentres");

    forAll(smoothVertex_, pointI)
    {
        if( smoothVertex_[pointI] & SMOOTH )
            pmg.addPointToSubset(smoothPointID, pointI);
        if( smoothVertex_[pointI] & BOUNDARY )
            pmg.addPointToSubset(bndPointID, pointI);
        if( smoothVertex_[pointI] & FACECENTRE )
            pmg.addPointToSubset(faceCentreID, pointI);
        if( smoothVertex_[pointI] & CELLCENTRE )
            pmg.addPointToSubset(cellCentreID, pointI);
    }

    const VRWGraph& internalOrdering = internalPointOrdering();
    const VRWGraph& boundaryOrdering = boundaryPointOrdering();

    forAll(internalOrdering, i)
    {
        const word name = "smoothPoints_"+help::scalarToText(i);
        const label orderID = pmg.addPointSubset(name);

        forAllRow(internalOrdering, i, nI)
            pmg.addPointToSubset(orderID, internalOrdering(i, nI));
    }

    forAll(boundaryOrdering, i)
    {
        const word name = "boundaryPoints_"+help::scalarToText(i);
        const label orderID = pmg.addPointSubset(name);

        forAllRow(boundaryOrdering, i, nI)
            pmg.addPointToSubset(orderID, boundaryOrdering(i, nI));
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

#include "extrudeLayer.H"
#include "helperFunctions.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "labelledPointScalar.H"

# ifdef USE_OMP
#include <omp.h>
# endif

#ifdef DEBUGExtrudeLayer
#include "polyMeshGenChecks.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Nested class definition

extrudeLayer::addressingCalculator::addressingCalculator
(
    const faceListPMG& faces,
    const LongList<labelPair>& extrudedFaces,
    const LongList<bool>& pairOrientation,
    const VRWGraph& pointExtruded
)
:
    faces_(faces),
    extrudedFaces_(extrudedFaces),
    pairOrientation_(pairOrientation),
    pointExtruded_(pointExtruded)
{}

extrudeLayer::addressingCalculator::~addressingCalculator()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void extrudeLayer::createDuplicateFrontFaces(const LongList<labelPair>& front)
{
    polyMeshGenModifier meshModifier(mesh_);

    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();
    faceListPMG& faces = meshModifier.facesAccess();
    cellListPMG& cells = meshModifier.cellsAccess();

    labelList faceInFront(faces.size(), -1);
    LongList<labelPair> newFaceLabels;
    label counter(0);
    forAll(front, ffI)
    {
        const labelPair& lp = front[ffI];

        if( faceInFront[lp.first()] == -1 )
        {
            faceInFront[lp.first()] = newFaceLabels.size();
            newFaceLabels.append(labelPair(-1, -1));
        }

        labelPair& cf = newFaceLabels[faceInFront[lp.first()]];

        if( (lp.second() == owner[lp.first()]) && (cf.first() == -1) )
        {
            cf.first() = counter;
            ++counter;
        }
        else if( (lp.second() == neighbour[lp.first()]) && (cf.second() == -1) )
        {
            cf.second() = counter;
            ++counter;
        }
    }

    //- create a copy of the faces in the extrusion front
    //- to create space for the layer
    faces.setSize(nOrigFaces_+counter);
    extrudedFaces_.setSize(counter);
    pairOrientation_.setSize(counter);

    # ifdef USE_OMP
    # pragma omp parallel for if( faceInFront.size() > 100 ) schedule(guided)
    # endif
    forAll(faceInFront, faceI)
    {
        if( faceInFront[faceI] < 0 )
            continue;

        const label fOwn = newFaceLabels[faceInFront[faceI]].first();
        const label fNei = newFaceLabels[faceInFront[faceI]].second();

        if( neighbour[faceI] < 0 )
        {
            //- boundary faces
            if( fOwn != -1 )
            {
                faces[nOrigFaces_+fOwn] = faces[faceI];
                extrudedFaces_[fOwn] = labelPair(nOrigFaces_+fOwn, faceI);
                pairOrientation_[fOwn] = true;
            }
        }
        else
        {
            //- internal face
            if( fOwn != -1 && fNei != -1 )
            {
                if( fOwn < fNei )
                {
                    faces[nOrigFaces_+fOwn] = faces[faceI];
                    extrudedFaces_[fOwn] = labelPair(nOrigFaces_+fOwn, faceI);
                    pairOrientation_[fOwn] = true;

                    faces[nOrigFaces_+fNei] = faces[faceI].reverseFace();
                    extrudedFaces_[fNei] = labelPair(nOrigFaces_+fNei, faceI);
                    pairOrientation_[fNei] = false;
                }
                else
                {
                    faces[nOrigFaces_+fOwn].transfer(faces[faceI]);
                    extrudedFaces_[fOwn] = labelPair(nOrigFaces_+fOwn, faceI);
                    pairOrientation_[fOwn] = false;

                    faces[faceI] = faces[nOrigFaces_+fOwn].reverseFace();

                    faces[nOrigFaces_+fNei] = faces[faceI];
                    extrudedFaces_[fNei] = labelPair(nOrigFaces_+fNei, faceI);
                    pairOrientation_[fNei] = true;
                }
            }
            else if( fOwn != -1 )
            {
                faces[nOrigFaces_+fOwn].transfer(faces[faceI]);
                faces[faceI] = faces[nOrigFaces_+fOwn].reverseFace();
                extrudedFaces_[fOwn] = labelPair(nOrigFaces_+fOwn, faceI);
                pairOrientation_[fOwn] = false;
            }
            else if( fNei != -1 )
            {
                faces[nOrigFaces_+fNei] = faces[faceI].reverseFace();
                extrudedFaces_[fNei] = labelPair(nOrigFaces_+fNei, faceI);
                pairOrientation_[fNei] = false;
            }
        }
    }

    //- renumber the cells
    # ifdef USE_OMP
    # pragma omp parallel for if( faceInFront.size() > 100 ) schedule(guided)
    # endif
    forAll(faceInFront, faceI)
    {
        if( faceInFront[faceI] < 0 )
            continue;

        const labelPair& lp = newFaceLabels[faceInFront[faceI]];

        if( lp.first() != -1 )
        {
            cell& c = cells[owner[faceI]];

            forAll(c, fI)
                if( c[fI] == faceI )
                    c[fI] = nOrigFaces_ + lp.first();
        }
        if( lp.second() != -1 )
        {
            cell& c = cells[neighbour[faceI]];

            forAll(c, fI)
                if( c[fI] == faceI )
                    c[fI] = nOrigFaces_ + lp.second();
        }
    }

    meshModifier.clearAll();
}

void extrudeLayer::createNewVertices()
{
    polyMeshGenModifier meshModifier(mesh_);
    meshModifier.clearAll();
    pointFieldPMG& points = meshModifier.pointsAccess();
    faceListPMG& faces = meshModifier.facesAccess();

    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    //- find the points in the marked front
    List<direction> frontPoints(points.size(), NONE);

    # ifdef USE_OMP
    # pragma omp parallel for if( points.size() > 1000 ) schedule(guided)
    # endif
    forAll(extrudedFaces_, efI)
    {
        const face& f = faces[extrudedFaces_[efI].first()];

        forAll(f, pI)
            frontPoints[f[pI]] |= FRONTVERTEX;
    }

    //- propagate this information to other processors in case of a parallel run
    if( Pstream::parRun() )
    {
        const polyMeshGenAddressing& addr = mesh_.addressingData();
        const VRWGraph& pAtProcs = addr.pointAtProcs();
        const Map<label>& globalToLocal = addr.globalToLocalPointAddressing();
        const DynList<label>& pProcs = addr.pointNeiProcs();

        //- allocate the map
        std::map<label, labelLongList> exchangeData;
        forAll(pProcs, i)
            exchangeData.insert(std::make_pair(pProcs[i], labelLongList()));

        //- collect the information about markes points at processor boundaries
        forAllConstIter(Map<label>, globalToLocal, it)
            if( frontPoints[it()] & FRONTVERTEX )
            {
                //- mark points at processor boundaries
                frontPoints[it()] |= FRONTVERTEXPROCBND;

                forAllRow(pAtProcs, it(), i)
                {
                    const label neiProc = pAtProcs(it(), i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append(it.key());
                }
            }

        LongList<label> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        # ifdef USE_OMP
        # pragma omp parallel for if( receivedData.size() > 1000 ) \
        schedule(guided)
        # endif
        forAll(receivedData, i)
        {
            frontPoints[globalToLocal[receivedData[i]]] =
                FRONTVERTEX+FRONTVERTEXPROCBND;
        }
    }

    //- create a new vertex for each group of faces
    //- faces are grouped such that the faces belonging to the same group
    //- can be visited over cells without crossing any front faces
    VRWGraph pointFaces;
    pointFaces.reverseAddressing(points.size(), faces);

    if( Pstream::parRun() )
    {
        Pout << "Creating new points at processor boundaries" << endl;
        for(label procI=0;procI<Pstream::nProcs();++procI)
        {
            if( Pstream::myProcNo() == procI )
            {
               Pout << "Front points are " << frontPoints << endl;
            }

            returnReduce(1, sumOp<label>());
        }

        //- create new vertices at processor boundaries
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();
        const polyMeshGenAddressing& addr = mesh_.addressingData();
        const VRWGraph& pAtProcs = addr.pointAtProcs();
        const labelLongList& globalPointLabel = addr.globalPointLabel();
        const Map<label>& globalToLocal = addr.globalToLocalPointAddressing();
        const DynList<label>& pProcs = addr.pointNeiProcs();
        const labelLongList& globalCellLabel = addr.globalCellLabel();

        //- create the information which faces are attached to points
        //- at parallel boundaries in dual form where each edge represents
        //- global labels of cells sharing a face
        typedef std::map<label, DynList<edge> > dualEdgesMap;
        dualEdgesMap procPointsDual;

        //- fill in local data
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            if( frontPoints[it()] & FRONTVERTEXPROCBND )
            {
                const label pointI = it();
                DynList<edge>& pEdges = procPointsDual[pointI];

                forAllRow(pointFaces, pointI, pfI)
                {
                    const label faceI = pointFaces(pointI, pfI);

                    //- do not store faces at processor boundaries
                    //- and newly created faces
                    if( (neighbour[faceI] < 0) || (owner[faceI] < 0) )
                        continue;

                    const edge e
                    (
                        globalCellLabel[owner[faceI]],
                        globalCellLabel[neighbour[faceI]]
                    );

                    pEdges.append(e);
                }
            }
        }

        for(label procI=0;procI<Pstream::nProcs();++procI)
        {
            if( Pstream::myProcNo() == procI )
            {
               forAllConstIter(dualEdgesMap, procPointsDual, it)
                    Pout << "Point " << it->first << " local dual edges " << it->second << endl;
            }

            returnReduce(1, sumOp<label>());
        }

        //- fill-in with data at processor boundaries. Store edges
        //- on the processor with the lower label not to duplicate the data
        returnReduce(1, sumOp<label>());
        Pout << "Adding data from processor boundaries" << endl;
        forAll(procBoundaries, patchI)
        {
            if( procBoundaries[patchI].owner() )
                continue;

            const label start = procBoundaries[patchI].patchStart();
            labelList globalLabels(procBoundaries[patchI].patchSize());
            forAll(globalLabels, fI)
            {
                if( owner[start+fI] < 0 )
                {
                    globalLabels[fI] = -1;
                }
                else
                {
                    globalLabels[fI] = globalCellLabel[owner[start+fI]];
                }
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                globalLabels.byteSize()
            );

            toOtherProc << globalLabels;
        }

        forAll(procBoundaries, patchI)
        {
            if( !procBoundaries[patchI].owner() )
                continue;

            labelList receivedData;
            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtherProc >> receivedData;

            const label start = procBoundaries[patchI].patchStart();

            forAll(receivedData, i)
            {
                if( owner[start+i] < 0 )
                    continue;

                const face& f = faces[start+i];

                forAll(f, pI)
                {
                    if( frontPoints[f[pI]] & FRONTVERTEXPROCBND )
                    {
                        dualEdgesMap::iterator dIter =
                            procPointsDual.find(f[pI]);
                        const label cOwn = globalCellLabel[owner[start+i]];
                        const label cNei = receivedData[i];
                        dIter->second.append(edge(cOwn, cNei));
                    }
                }
            }
        }

        //- exchange this information with neighbouring processors
        returnReduce(1, sumOp<label>());
        Pout << "Exchanging data with other processors" << endl;

        std::map<label, labelLongList> exchangeData;
        forAll(pProcs, i)
            exchangeData.insert(std::make_pair(pProcs[i], labelLongList()));

        //- fill in the exchangeData map
        forAllConstIter(dualEdgesMap, procPointsDual, dIter)
        {
            const label pointI = dIter->first;

            forAllRow(pAtProcs, pointI, i)
            {
                const label neiProc = pAtProcs(pointI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                labelLongList& dts = exchangeData[neiProc];

                dts.append(globalPointLabel[pointI]);

                const DynList<edge>& dualEdges = dIter->second;
                dts.append(dualEdges.size());
                forAll(dualEdges, edgeI)
                {
                    const edge& e = dualEdges[edgeI];
                    dts.append(e.start());
                    dts.append(e.end());
                }
            }
        }

        //- exchange data with other processors
        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        //- update local data
        label counter(0);
        while( counter < receivedData.size() )
        {
            const label pointI = globalToLocal[receivedData[counter++]];

            DynList<edge>& dualEdges = procPointsDual[pointI];

            const label nDualEdges = receivedData[counter++];
            for(label eI=0;eI<nDualEdges;++eI)
            {
                edge e;
                e.start() = receivedData[counter++];
                e.end() = receivedData[counter++];

                dualEdges.append(e);
            }
        }

        for(label procI=0;procI<Pstream::nProcs();++procI)
        {
            if( Pstream::myProcNo() == procI )
            {
               forAllConstIter(dualEdgesMap, procPointsDual, it)
                    Pout << "Point " << it->first << " dual edges " << it->second << endl;
            }

            returnReduce(1, sumOp<label>());
        }

        //- Finally, find groups of faces and create new vertices
        returnReduce(1, sumOp<label>());
        Pout << "Finding groups of edges at vertex" << endl;
        forAllConstIter(dualEdgesMap, procPointsDual, dIter)
        {
            const label pointI = dIter->first;
            const DynList<edge>& dEdges = dIter->second;

            DynList<label> edgeGroup;
            edgeGroup.setSize(dEdges.size());
            edgeGroup = -1;

            //- check edge connections and store all edges which can be reached
            //- over other edges into the same group
            label group(0);

            forAll(dEdges, eI)
            {
                if( edgeGroup[eI] != -1 )
                    continue;

                edgeGroup[eI] = group;
                DynList<label> front;
                front.append(eI);

                while( front.size() != 0 )
                {
                    const label eLabel = front.removeLastElement();

                    forAll(dEdges, eJ)
                    {
                        if( edgeGroup[eJ] != -1 )
                            continue;
                        if( eJ == eLabel )
                            continue;

                        if( dEdges[eLabel].commonVertex(dEdges[eJ]) != -1 )
                        {
                            front.append(eJ);
                            edgeGroup[eJ] = group;
                        }
                    }
                }

                ++group;
            }

            //- find face groups from the groups assigned to dual edges
            DynList<label> faceGroups;
            faceGroups.setSize(pointFaces.sizeOfRow(pointI));
            faceGroups = -1;

            forAllRow(pointFaces, pointI, pfI)
            {
                const label faceI = pointFaces(pointI, pfI);

                if( owner[faceI] < 0 )
                    continue;

                const label own = globalCellLabel[owner[faceI]];

                forAll(dEdges, eI)
                {
                    const edge& dualEdge = dEdges[eI];

                    if( (own == dualEdge.start()) || (own == dualEdge.end()) )
                    {
                        faceGroups[pfI] = edgeGroup[eI];
                        break;
                    }
                }
            }

            Info << "Edge groups for point " << pointI << " are " << edgeGroup << endl;
            Info << "Face groups at point " << pointI << " are " << faceGroups
                << " point faces " << pointFaces[pointI] << endl;

            //- stop in case there is only one group
            //- of faces attached to this point
            if( group == 1 )
            {
                bool problem(true);
                forAll(faceGroups, i)
                    if( faceGroups[i] != 0 )
                    {
                        problem = false;
                        break;
                    }

                if( problem )
                    FatalErrorIn("void extrudeLayer::createNewVertices()")
                        << "Front is not defined at point " << pointI
                        << ". Cannot continue.." << exit(FatalError);
            }

            //- create new vertices
            const label currentNumPoints = points.size();
            for(label i=0;i<group;++i)
            {
                const point p = points[pointI];
                origPointLabel_.append(pointI);
                points.append(p);
            }

            //- renumber faces
            forAllRow(pointFaces, pointI, pfI)
            {
                if( faceGroups[pfI] == -1 )
                    continue;

                face& f = faces[pointFaces(pointI, pfI)];

                const label pos = f.which(pointI);
                f[pos] = currentNumPoints + faceGroups[pfI];
            }
        }

        Pout << "Finishing creating new vertices at paralle boundaries" << endl;
        returnReduce(1, sumOp<label>());
    }

    //- treat local points
    forAll(pointFaces, pointI)
    {
        if( frontPoints[pointI] != FRONTVERTEX )
            continue;

        //- assign groups to faces and cells
        DynList<label> faceGroup;
        faceGroup.setSize(pointFaces.sizeOfRow(pointI));
        faceGroup = -1;

        label group(0);

        forAllRow(pointFaces, pointI, pfI)
        {
            if( pointFaces(pointI, pfI) < nOrigFaces_ )
                continue;
            if( faceGroup[pfI] != -1 )
                continue;

            DynList<label> front;
            front.append(pfI);
            faceGroup[pfI] = group;

            while( front.size() )
            {
                const label fLabel = front.removeLastElement();

                const label own = owner[pointFaces(pointI, fLabel)];
                const label nei = neighbour[pointFaces(pointI, fLabel)];

                forAllRow(pointFaces, pointI, pfJ)
                {
                    if( faceGroup[pfJ] != -1 )
                        continue;

                    const label faceJ = pointFaces(pointI, pfJ);

                    if( owner[faceJ] < 0 )
                        continue;

                    if( owner[faceJ] == own || owner[faceJ] == nei )
                    {
                        faceGroup[pfJ] = group;
                        front.append(pfJ);
                    }

                    if( neighbour[faceJ] < 0 )
                        continue;

                    if( neighbour[faceJ] == own || neighbour[faceJ] == nei )
                    {
                        faceGroup[pfJ] = group;
                        front.append(pfJ);
                    }
                }
            }

            ++group;
        }

        //- stop in case there is only one group of faces attached to this point
        if( group == 1 )
        {
            bool problem(true);
            forAll(faceGroup, i)
                if( faceGroup[i] != 0 )
                {
                    problem = false;
                    break;
                }

            if( problem )
                FatalErrorIn("void extrudeLayer::createNewVertices()")
                    << "Front is not defined at point " << pointI
                    << ". Cannot continue.." << exit(FatalError);
        }

        //- create new vertices
        const label currentNumPoints = points.size();
        for(label i=0;i<group;++i)
        {
            const point p = points[pointI];
            origPointLabel_.append(pointI);
            points.append(p);
        }

        //- renumber faces
        forAllRow(pointFaces, pointI, pfI)
        {
            if( faceGroup[pfI] == -1 )
                continue;

            face& f = faces[pointFaces(pointI, pfI)];

            const label pos = f.which(pointI);
            f[pos] = currentNumPoints + faceGroup[pfI];
        }
    }

    polyMeshGenModifier(mesh_).clearOut();

    # ifdef DEBUGExtrudeLayer
    const label procID = mesh_.addPointSubset("parPoints");
    const label frontID = mesh_.addPointSubset("frontPoints");
    forAll(frontPoints, pI)
    {
        if( frontPoints[pI] & FRONTVERTEXPROCBND )
            mesh_.addPointToSubset(procID, pI);
        if( frontPoints[pI] & FRONTVERTEX )
            mesh_.addPointToSubset(frontID, pI);
    }

    returnReduce(1, sumOp<label>());
    //::exit(1);
    # endif
}

void extrudeLayer::movePoints()
{
    pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    vectorField displacements(points.size()-nOrigPoints_);
    boolList pointAtProcBnd(displacements.size(), false);

    VRWGraph pointFaces;
    pointFaces.reverseAddressing(points.size(), mesh_.faces());

    if( Pstream::parRun() )
    {
        const polyMeshGenAddressing& addr = mesh_.addressingData();
        const Map<label>& globalToLocal = addr.globalToLocalPointAddressing();
        const VRWGraph& pAtProcs = addr.pointAtProcs();
        const DynList<label>& pProcs = addr.pointNeiProcs();

        std::map<label, vector> normals;
        std::map<label, scalar> distances;

        //- create a map for exchanging data
        std::map<label, LongList<labelledPointScalar> > exchangeData;
        forAll(pProcs, i)
            exchangeData.insert
            (
                std::make_pair(pProcs[i], LongList<labelledPointScalar>())
            );

        //- create displacements from local data
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            if( it() >= nOrigPoints_ )
            {
                const label pointI = it();

                pointAtProcBnd[pointI-nOrigPoints_] = true;

                //- create local part of the displacement vector
                labelledPointScalar lps;
                lps.pointLabel() = it.key();
                lps.coordinates() = vector::zero;
                lps.scalarValue() = VGREAT;

                forAllRow(pointFaces, pointI, pfI)
                {
                    const label faceI = pointFaces(pointI, pfI);

                    if( faceI >= nOrigFaces_ )
                    {
                        const face& f = faces[faceI];

                        lps.coordinates() -= f.normal(points);

                        if( thickness_ < 0.0 )
                        {
                            const vector c = f.centre(points);
                            scalar d(VGREAT);

                            forAll(f, pI)
                                d = Foam::min(d, Foam::mag(points[f[pI]] - c));

                            d *= 0.4;
                            lps.scalarValue() = Foam::min(lps.scalarValue(), d);
                        }
                        else
                        {
                            lps.scalarValue() = thickness_;
                        }
                    }
                }

                normals[pointI] = lps.coordinates();
                distances[pointI] = lps.scalarValue();

                //- store data in the exchangeData map
                forAllRow(pAtProcs, pointI, i)
                {
                    const label neiProc = pAtProcs(pointI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append(lps);
                }
            }
        }

        LongList<labelledPointScalar> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const labelledPointScalar& lps = receivedData[i];
            const label pointI = globalToLocal[lps.pointLabel()];

            normals[pointI] += lps.coordinates();
            scalar& d = distances[pointI];
            d = Foam::min(d, lps.scalarValue());
        }

        //- calculate displacements of vertices at processor boundaries
        for
        (
            std::map<label, vector>::const_iterator it=normals.begin();
            it!=normals.end();
            ++it
        )
        {
            vector n = it->second;
            if( mag(n) > VSMALL )
                n /= mag(n);
            displacements[it->first-nOrigPoints_] = n * distances[it->first];
        }
    }

    # ifdef USE_OMP
    # pragma omp parallel if( displacements.size() > 100 )
    # endif
    {
        //- find displacement vectors
        # ifdef USE_OMP
        # pragma omp for schedule(guided)
        # endif
        forAll(displacements, pI)
        {
            if( pointAtProcBnd[pI] )
                continue;

            const label pointI = nOrigPoints_ + pI;

            vector normal(vector::zero);
            scalar thickness(VGREAT);

            forAllRow(pointFaces, pointI, pfI)
            {
                const label faceI = pointFaces(pointI, pfI);

                if( faceI >= nOrigFaces_ )
                {
                    const face& f = faces[faceI];

                    normal -= f.normal(points);

                    if( thickness_ < 0.0 )
                    {
                        const vector c = f.centre(points);
                        scalar d(VGREAT);

                        forAll(f, pI)
                            d = Foam::min(d, Foam::mag(points[f[pI]] - c));

                        thickness = Foam::min(thickness, d);
                    }
                }
            }

            if( thickness_ >= 0.0 )
            {
                thickness = thickness_;
            }
            else
            {
                thickness *= 0.4;
            }

            const scalar d = mag(normal);
            if( d > VSMALL )
                normal /= d;

            displacements[pI] = normal * thickness;
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        # ifdef USE_OMP
        # pragma omp for schedule(guided)
        # endif
        forAll(displacements, pI)
            points[nOrigPoints_+pI] += displacements[pI];
    }

    # ifdef DEBUGExtrudeLayer
    mesh_.write();
    returnReduce(1, sumOp<label>());
    //::exit(1);
    # endif
}

void extrudeLayer::createNewFacesParallel()
{
    if( !Pstream::parRun() )
        return;

    VRWGraph newProcFaces;
    labelLongList faceProcPatch;

    //- add faces into the mesh
    polyMeshGenModifier(mesh_).addProcessorFaces(newProcFaces, faceProcPatch);
}

void extrudeLayer::createLayerCells()
{
    const faceListPMG& faces = mesh_.faces();

    VRWGraphList newCells;

    //- create cells from corners and edges
    faceList::subList newFaces(faces, faces.size()-nOrigFaces_, nOrigFaces_);
    VRWGraph pointFaces;
    pointFaces.reverseAddressing(mesh_.points().size(), newFaces);

    //- create new cells extruded from the faces in the selected front
    forAll(extrudedFaces_, extrudedI)
    {
        const face& f = faces[extrudedFaces_[extrudedI].first()];
        const face& of = faces[extrudedFaces_[extrudedI].second()];

        //- create new cell from the front pair
        DynList<DynList<label, 4> > newCell;
        newCell.setSize(f.size()+2);

        newCell[0] = f;
        newCell[1] = of;

        if( pairOrientation_[extrudedI] )
        {
            //- create a new cell. Faces are of the same orientation
            forAll(f, pI)
            {
                DynList<label, 4>& ef = newCell[pI+2];

                ef.setSize(4);
                ef[0] = f[pI];
                ef[1] = f.nextLabel(pI);
                ef[2] = of[f.fcIndex(pI)];
                ef[3] = of[pI];
            }
        }
        else
        {
            //- create a new cell. Faces are of the opposite orientation
            forAll(f, pI)
            {
                DynList<label, 4>& ef = newCell[pI+2];

                ef.setSize(4);
                ef[0] = f[pI];
                ef[1] = f.nextLabel(pI);
                ef[2] = of[(f.size()-f.fcIndex(pI))%f.size()];
                ef[3] = of[(f.size()-pI)%f.size()];
            }
        }

        newCells.appendGraph(newCell);
    }

    addressingCalculator adc
    (
        faces,
        extrudedFaces_,
        pairOrientation_,
        pointFaces
    );

    //- create cells created as a consequence of self-intersection
    //- over edges. And edge is transformed into a hex cell
    forAll(extrudedFaces_, extrudedI)
    {
        const face& f = faces[extrudedFaces_[extrudedI].first()];

        forAll(f, eI)
        {
            const label pointI = f[eI];
            if( pointFaces.sizeOfRow(pointI) == 0 )
                continue;
            const label nextI = f.nextLabel(eI);
            if( pointFaces.sizeOfRow(nextI) == 0 )
                continue;
            if( pointI < nextI )
                continue;

            //- find point labels of the edge which is the origin
            //- of the edge (pointI, nextI) with respect to face extrudedI
            const label origFacePointI = adc.origPoint(extrudedI, pointI);
            const label origFaceNextI = adc.origPoint(extrudedI, nextI);

            //- points of the current edge and of the original edge must have
            //- the same point as their origin
            if( origPointLabel_[pointI] != origPointLabel_[origFacePointI] )
                continue;
            if( origPointLabel_[nextI] != origPointLabel_[origFaceNextI] )
                continue;

            //- Finally, start creating a cell from this edge
            //- find the other extruded face which shares this edge
            const label otherExtI = adc.faceSharingEdge(extrudedI, eI);

            //- find point labels of the edge which is the origin
            //- of the edge (pointI, nextI)
            //- with respect to face otherExtrudedFace
            const label otherOrigFacePointI = adc.origPoint(otherExtI, pointI);
            const label otherOrigFaceNextI = adc.origPoint(otherExtI, nextI);

            //- find points of the edge opposite to the edge (pointI, nextI)
            //- these points make the original edge which is blown into
            //- a hex cell
            label origPointI(-1), origNextI(-1);

            if(
                (origFacePointI >= nOrigPoints_) &&
                (origFaceNextI >= nOrigPoints_)
            )
            {
                //- find an extruded face attached to edge
                //- (origFacePointI, origFaceNextI). The must exist only one
                //- such face
                DynList<label> fc;
                adc.facesSharingEdge(origFacePointI, origFaceNextI, fc);

                if( fc.size() != 1 )
                    FatalErrorIn("void extrudeLayer::createLayerCells()")
                        << "Cannot find searched face" << abort(FatalError);

                origPointI = adc.origPoint(fc[0], origFacePointI);
                origNextI = adc.origPoint(fc[0], origFaceNextI);
            }
            else
            {
                FatalErrorIn("void extrudeLayer::createLayerCells()")
                    << "Cannot find points" << abort(FatalError);
            }

            //- create new cell and add it to the list
            FixedList<FixedList<label, 4>, 6> newCell;

            //- face 0
            newCell[0][0] = pointI;
            newCell[0][1] = origFacePointI;
            newCell[0][2] = origFaceNextI;
            newCell[0][3] = nextI;

            //- face 1
            newCell[1][0] = pointI;
            newCell[1][1] = nextI;
            newCell[1][2] = otherOrigFaceNextI;
            newCell[1][3] = otherOrigFacePointI;

            //- face 2
            newCell[2][0] = otherOrigFacePointI;
            newCell[2][1] = otherOrigFaceNextI;
            newCell[2][2] = origNextI;
            newCell[2][3] = origPointI;

            //- face 3
            newCell[3][0] = origFacePointI;
            newCell[3][1] = origPointI;
            newCell[3][2] = origNextI;
            newCell[3][3] = origFaceNextI;

            //- face 4
            newCell[4][0] = pointI;
            newCell[4][1] = otherOrigFacePointI;
            newCell[4][2] = origPointI;
            newCell[4][3] = origFacePointI;

            //- face 5
            newCell[5][0] = nextI;
            newCell[5][1] = origFaceNextI;
            newCell[5][2] = origNextI;
            newCell[5][3] = otherOrigFaceNextI;

            newCells.appendGraph(newCell);
        }
    }

    //- create cells at points where three or more self-intersecting layers meet
    forAll(pointFaces, pointI)
    {
        if( pointFaces.sizeOfRow(pointI) < 3 )
            continue;

        //- find labels of points
        DynList<label> origFacePoints;
        origFacePoints.setSize(pointFaces.sizeOfRow(pointI));
        origFacePoints = -1;

        forAllRow(pointFaces, pointI, pfI)
        {
            const label extI = pointFaces(pointI, pfI);

            origFacePoints[pfI] = adc.origPoint(extI, pointI);
        }

        //- cell shall be created only if all points are different
        bool createCell(true);
        for(label i=origFacePoints.size()-1;i>0;--i)
        {
            for(label j=i-1;j>=0;--j)
                if( origFacePoints[i] == origFacePoints[j] )
                {
                    createCell = false;
                    break;
                }

            if( !createCell )
                break;
        }

        if( !createCell )
            continue;

        DynList<FixedList<label, 4>, 6> newCell;
        newCell.setSize(2*origFacePoints.size());

        //- start creating faces attached to the pointI
        //- this is performed by finding original points with respect to
        //- the face extI and the face which shares edge eI of the face extI
        //- this is needed to ensure correct normal orientation
        forAllRow(pointFaces, pointI, pfI)
        {
            const label extI = pointFaces(pointI, pfI);

            //- find original labels of points making a forward circular edge
            //- with respect to pointI
            const face& f = faces[extrudedFaces_[extI].first()];
            const label eI = f.which(pointI);
            const label origFacePointI = origFacePoints[pfI];
            const label origFaceNextI = adc.origPointLabel(extI, f.fcIndex(eI));

            //- find another face attached to the edge eI
            const label fFace = adc.faceSharingEdge(extI, eI);
            const label pos = pointFaces.containsAtPosition(pointI, fFace);

            //- find an extrusion face attached to the edge consisting of points
            //- origFacePointI and origFaceNextI
            DynList<label> fe;
            adc.facesSharingEdge(origFacePointI, origFaceNextI, fe);
            const label origPointI = adc.origPoint(fe[0], origFacePointI);

            //- create a face attached to pointI
            FixedList<label, 4>& cf = newCell[pfI];
            cf[0] = pointI;
            cf[1] = origFacePoints[pfI];
            cf[2] = origPointI;
            cf[3] = origFacePoints[pos];
        }

        //- close the cell by creating new faces from the existing
        //- faces which obey pre-determined order. If a face contains
        //- a point in the origFacePoints list at the second position, then
        //- the point after shall be the previous point of the face. If a face
        //- contains such point at the last position then the point before it
        //- shall be the next point in the face. The last point of the face
        //- is the original points from which all these points were generated.
        forAll(origFacePoints, opI)
        {
            FixedList<label, 4>& cf = newCell[opI+origFacePoints.size()];

            //- find previous and next point
            label prev(-1), next(-1);
            forAll(origFacePoints, fJ)
            {
                if( newCell[fJ][1] == origFacePoints[opI] )
                    prev = newCell[fJ][2];
                if( newCell[fJ][3] == origFacePoints[opI] )
                    next = newCell[fJ][2];
            }

            if( (prev < 0) || (next < 0) )
                FatalErrorIn("void extrudeLayer::createLayerCells()")
                    << "Corner cell is invalid" << abort(FatalError);

            cf[0] = prev;
            cf[1] = origFacePoints[opI];
            cf[2] = next;
            cf[3] = origPointLabel_[pointI];
        }

        newCells.appendGraph(newCell);
    }

    //- create new faces at parallel boundaries
    createNewFacesParallel();

    //- add cells into the mesh
    polyMeshGenModifier(mesh_).addCells(newCells);
}

void extrudeLayer::updateBoundary()
{
    wordList patchNames(mesh_.boundaries().size());
    wordList patchTypes(mesh_.boundaries().size());
    forAll(patchNames, patchI)
    {
        patchNames[patchI] = mesh_.boundaries()[patchI].patchName();
        patchTypes[patchI] = mesh_.boundaries()[patchI].patchType();
    }

    VRWGraph newBoundaryFaces;
    labelLongList newBoundaryOwners;
    labelLongList newBoundaryPatches;

    meshSurfaceEngine mse(mesh_);
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bfOwner = mse.faceOwners();
    const labelList& facePatch = mse.boundaryFacePatches();
    const labelList& bp = mse.bp();

    meshSurfacePartitioner mPart(mse);
    const VRWGraph& pointPatches = mPart.pointPatches();

    //- store existing boundary faces. They remain in the mesh
    forAll(bFaces, bfI)
    {
        newBoundaryFaces.appendList(bFaces[bfI]);
        newBoundaryOwners.append(bfOwner[bfI]);
        newBoundaryPatches.append(facePatch[bfI]);
    }

    //- find and store boundary faces which have been generated as a consequence
    //- of layer insertion
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    for(label faceI=nOrigFaces_;faceI<faces.size();++faceI)
    {
        if( neighbour[faceI] >= 0 )
            continue;

        const face& f = faces[faceI];

        if( f.size() != 4 )
            continue;

        DynList<label> origPts, newPts;
        forAll(f, pI)
            if( origPointLabel_[f[pI]] < 0 )
            {
                origPts.appendIfNotIn(f[pI]);
            }
            else
            {
                origPts.appendIfNotIn(origPointLabel_[f[pI]]);
                newPts.append(f[pI]);
            }

        if( origPts.size() > 2 )
            continue;
        if( newPts.size() == 0 )
            continue;

        //- this face belongs to the extruded layer
        //- find patches of boundary faces attached to the newly created points
        DynList<label> patches;
        DynList<label> nAppearances;
        forAll(newPts, npI)
        {
            const label bpI = bp[newPts[npI]];

            if( bpI < 0 )
                continue;

            forAllRow(pointPatches, bpI, ppI)
            {
                const label patchI = pointPatches(bpI, ppI);
                if( patches.contains(patchI) )
                {
                    ++nAppearances[patches.containsAtPosition(patchI)];
                }
                else
                {
                    patches.append(patchI);
                    nAppearances.append(1);
                }
            }
        }

        label patch(-1);

        forAll(patches, i)
        {
            if( nAppearances[i] == newPts.size() )
            {
                if( patch != -1 )
                    FatalErrorIn("void extrudeLayer::updateBoundary()")
                        << "Found more than one patch!" << abort(FatalError);

                patch = patches[i];
            }
        }

        if( patch != -1 )
        {
            //- append boundary face
            newBoundaryFaces.appendList(f);
            newBoundaryOwners.append(owner[faceI]);
            newBoundaryPatches.append(patch);
        }
    }

    polyMeshGenModifier meshModifier(mesh_);
    meshModifier.reorderBoundaryFaces();
    meshModifier.replaceBoundary
    (
        patchNames,
        newBoundaryFaces,
        newBoundaryOwners,
        newBoundaryPatches
    );

    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
        boundaries[patchI].patchType() = patchTypes[patchI];

    meshModifier.clearAll();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh reference
extrudeLayer::extrudeLayer
(
    polyMeshGen& mesh,
    const LongList<labelPair>& extrusionFront,
    const scalar thickness
)
:
    mesh_(mesh),
    thickness_(thickness),
    nOrigPoints_(mesh.points().size()),
    nOrigFaces_(mesh.faces().size()),
    nOrigCells_(mesh.cells().size()),
    extrudedFaces_(),
    pairOrientation_(),
    origPointLabel_(nOrigPoints_, -1)
{
    createDuplicateFrontFaces(extrusionFront);

    createNewVertices();

    movePoints();

    createLayerCells();

    updateBoundary();

    mesh_.clearAddressingData();

    # ifdef DEBUGExtrudeLayer
    polyMeshGenChecks::checkMesh(mesh_, true);
    # endif
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

extrudeLayer::~extrudeLayer()
{
    mesh_.clearAddressingData();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

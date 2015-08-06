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
#include "boundaryLayerOptimisation.H"
#include "meshSurfacePartitioner.H"
#include "meshSurfaceEngine.H"
#include "detectBoundaryLayers.H"
//#include "helperFunctions.H"
//#include "labelledScalar.H"
//#include "refLabelledPoint.H"
//#include "refLabelledPointScalar.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngineModifier.H"
//#include "partTetMeshSimplex.H"
//#include "volumeOptimizer.H"
#include "OFstream.H"

//#define DEBUGLayer

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayerOptimisation::writeVTK
(
    const fileName& fName,
    const pointField& origin,
    const vectorField& vecs
)
{
    if( origin.size() != vecs.size() )
        FatalErrorIn
        (
            "void boundaryLayerOptimisation::writeVTK(const fileName&,"
            " const pointField&, const vectorField&)"
        ) << "Sizes do not match" << abort(FatalError);

    OFstream file(fName);

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "POINTS " << 2*origin.size() << " float\n";
    forAll(origin, pI)
    {
        const point& p = origin[pI];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;

        const point op = p + vecs[pI];

        file << op.x() << ' ' << op.y() << ' ' << op.z() << nl;
    }

    //- write lines
    file << "\nLINES " << vecs.size()
         << " " << 3*vecs.size() << nl;
    forAll(vecs, eI)
    {
        file << 2 << " " << 2*eI << " " << (2*eI+1) << nl;
    }

    file << "\n";
}

void boundaryLayerOptimisation::writeHairEdges
(
    const fileName& fName,
    const direction eType,
    const vectorField& vecs
) const
{
    if( vecs.size() != hairEdges_.size() )
        FatalErrorIn
        (
            "void boundaryLayerOptimisation::writeHairEdges"
            "(const fileName&, const direction, const vectorField&) const"
        ) << "Sizes do not match" << abort(FatalError);

    //- count the number of hair edges matching this criteria
    label counter(0);

    forAll(hairEdgeType_, heI)
        if( hairEdgeType_[heI] & eType )
            ++counter;

    //- copy edge vector
    vectorField copyVecs(counter);
    pointField pts(counter);

    counter = 0;

    const pointFieldPMG& points = mesh_.points();

    forAll(hairEdgeType_, heI)
    {
        if( hairEdgeType_[heI] & eType )
        {
            const edge& he = hairEdges_[heI];

            pts[counter] = points[he.start()];
            copyVecs[counter] = vecs[heI] * he.mag(points);

            ++counter;
        }
    }

    //- write data to file
    writeVTK(fName, pts, copyVecs);
}

void boundaryLayerOptimisation::writeHairEdges
(
    const fileName& fName,
    const direction eType
) const
{
    //- count the number of hair edges matching this criteria
    label counter(0);

    forAll(hairEdgeType_, heI)
        if( hairEdgeType_[heI] & eType )
            ++counter;

    //- copy edge vector
    vectorField vecs(counter);
    pointField pts(counter);

    counter = 0;

    const pointFieldPMG& points = mesh_.points();

    forAll(hairEdgeType_, heI)
    {
        if( hairEdgeType_[heI] & eType )
        {
            const edge& he = hairEdges_[heI];

            pts[counter] = points[he.start()];
            vecs[counter] = he.vec(points);

            ++counter;
        }
    }

    //- write data to file
    writeVTK(fName,pts, vecs);
}

const meshSurfaceEngine& boundaryLayerOptimisation::meshSurface() const
{
    if( !meshSurfacePtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const meshSurfaceEngine&"
                " boundaryLayerOptimisation::meshSurface()"
            ) << "Cannot generate meshSurfaceEngine" << abort(FatalError);
        # endif

        meshSurfacePtr_ = new meshSurfaceEngine(mesh_);
    }

    return *meshSurfacePtr_;
}

const meshSurfacePartitioner&
boundaryLayerOptimisation::surfacePartitioner() const
{
    if( !partitionerPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const meshSurfacePartitioner& "
                "boundaryLayerOptimisation::surfacePartitioner()"
            ) << "Cannot generate meshSurfacePartitioner" << abort(FatalError);
        # endif

        partitionerPtr_ = new meshSurfacePartitioner(meshSurface());
    }

    return *partitionerPtr_;
}

void boundaryLayerOptimisation::calculateHairEdges()
{
    const meshSurfaceEngine& mse = meshSurface();
    const edgeList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const labelList& faceOwner = mse.faceOwners();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();

    const meshSurfacePartitioner& mPart = surfacePartitioner();

    //- detect layers in the mesh
    const detectBoundaryLayers detectLayers(mPart);

    hairEdges_ = detectLayers.hairEdges();
    hairEdgesAtBndPoint_ = detectLayers.hairEdgesAtBndPoint();

    //- mark boundary faces which are base face for the boundary layer
    const labelList& layerAtBndFace = detectLayers.faceInLayer();
    isBndLayerBase_.setSize(bFaces.size());
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(layerAtBndFace, bfI)
    {
        if( layerAtBndFace[bfI] < 0 )
        {
            isBndLayerBase_[bfI] = false;
        }
        else
        {
            isBndLayerBase_[bfI] = true;
        }
    }

    # ifdef DEBUGLayer
    const label bndLayerFaceId = mesh_.addFaceSubset("bndLayerFaces");
    const label startBndFaceI = mesh_.boundaries()[0].patchStart();
    forAll(isBndLayerBase_, bfI)
        if( isBndLayerBase_[bfI] )
            mesh_.addFaceToSubset(bndLayerFaceId, startBndFaceI+bfI);
    # endif

    //- check if a face is an exiting face for a bnd layer
    isExitFace_.setSize(isBndLayerBase_.size());
    isExitFace_ = false;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(edgeFaces, edgeI)
    {
        //- avoid edges at inter-processor boundaries
        if( edgeFaces.sizeOfRow(edgeI) != 2 )
            continue;

        const label f0 = edgeFaces(edgeI, 0);
        const label f1 = edgeFaces(edgeI, 1);

        //- both faces have to be part of the same cell
        if( faceOwner[f0] != faceOwner[f1] )
            continue;

        //- check if the feature edge is attachd to exittin faces
        if
        (
            (isBndLayerBase_[f0] && (bFaces[f1].size() == 4)) &&
            (isBndLayerBase_[f1] && (bFaces[f0].size() == 4))
        )
        {
            isExitFace_[f0] = true;
            isExitFace_[f1] = true;
        }
    }

    # ifdef DEBUGLayer
    const label exittingFaceId = mesh_.addFaceSubset("exittingFaces");
    forAll(isExitFace_, bfI)
        if( isExitFace_[bfI] )
            mesh_.addFaceToSubset(exittingFaceId, startBndFaceI+bfI);
    # endif

    //- classify hair edges as they require different tretment
    //- in the smoothing procedure
    hairEdgeType_.setSize(hairEdges_.size());

    const labelHashSet& corners = mPart.corners();
    const labelHashSet& edgePoints = mPart.edgePoints();
    const labelHashSet& featureEdges = mPart.featureEdges();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        hairEdgeType_[hairEdgeI] = NONE;

        const edge& e = hairEdges_[hairEdgeI];
        const label bpI = bp[e.start()];

        //- check if it is a boundary edge
        forAllRow(bpEdges, bpI, peI)
        {
            const label beI = bpEdges(bpI, peI);

            const edge& be = edges[bpEdges(bpI, peI)];

            if( be == e )
            {
                hairEdgeType_[hairEdgeI] |= BOUNDARY;

                if( featureEdges.found(beI) )
                    hairEdgeType_[hairEdgeI] |= FEATUREEDGE;
            }

            if( corners.found(bpI) )
            {
                hairEdgeType_[hairEdgeI] |= ATCORNER;
            }
            else if( edgePoints.found(bpI) )
            {
                hairEdgeType_[hairEdgeI] |= ATEDGE;
            }
        }

        if( !(hairEdgeType_[hairEdgeI] & BOUNDARY) )
            hairEdgeType_[hairEdgeI] |= INSIDE;
    }

    thinnedHairEdge_.setSize(hairEdges_.size());

    //- calculate which other hair edges influence a hair edges
    //- and store it in a graph
    hairEdgesNearHairEdge_.setSize(hairEdges_.size());

    const cellListPMG& cells = mesh_.cells();
    const faceList& faces = mesh_.faces();

    VRWGraph bpFacesHelper(bpEdges.size());
    forAll(faceOwner, bfI)
    {
        const label cellI = faceOwner[bfI];

        const cell& c = cells[cellI];

        forAll(c, fI)
        {
            const face& f = faces[c[fI]];

            forAll(f, pI)
            {
                const label bpI = bp[f[pI]];
                if( bpI < 0 )
                    continue;

                bpFacesHelper.appendIfNotIn(bpI, c[fI]);
            }
        }
    }

    forAll(hairEdges_, hairEdgeI)
    {
        const edge& e = hairEdges_[hairEdgeI];
        const label bpI = bp[e.start()];

        DynList<label> neiHairEdges;

        //- find mesh faces comprising of the current hair edge
        forAllRow(bpFacesHelper, bpI, pfI)
        {
            const face& f = faces[bpFacesHelper(bpI, pfI)];

            //- face must be a quad
            if( f.size() != 4 )
                continue;

            //- check if the current face comprises of the hair edge
            label faceEdge(-1);
            forAll(f, eI)
                if( f.faceEdge(eI) == e )
                {
                    faceEdge = eI;
                    break;
                }

            if( faceEdge != -1 )
            {
                //- check if the opposite edge is also a hair edge
                const label eJ = (faceEdge+2) % 4;

                const edge fe = f.faceEdge(eJ);

                for(label i=0;i<2;++i)
                {
                    const label bpJ = bp[fe[i]];

                    if( bpJ >= 0 )
                    {
                        forAllRow(hairEdgesAtBndPoint_, bpJ, pI)
                        {
                            const label heJ = hairEdgesAtBndPoint_(bpJ, pI);
                            if( hairEdges_[heJ] == fe )
                                neiHairEdges.append(heJ);
                        }
                    }
                }
            }
        }

        hairEdgesNearHairEdge_.setRow(hairEdgeI, neiHairEdges);
    }

    # ifdef DEBUGLayer
    const label hairEdgesId = mesh_.addPointSubset("hairEdgePoints");
    const label bndHairEdgeId = mesh_.addPointSubset("bndHairEdgePoints");
    const label featureHairEdgeId = mesh_.addPointSubset("featureEdgePoints");
    const label cornerHairEdgeId = mesh_.addPointSubset("cornerHairEdgePoints");
    const label hairEdgeAtEdgeId = mesh_.addPointSubset("hairEdgeAtEdgePoints");

    forAll(hairEdgeType_, heI)
    {
        const edge& e = hairEdges_[heI];

        mesh_.addPointToSubset(hairEdgesId, e.start());
        mesh_.addPointToSubset(hairEdgesId, e.end());
        if( hairEdgeType_[heI] & FEATUREEDGE)
        {
            mesh_.addPointToSubset(featureHairEdgeId, e.start());
            mesh_.addPointToSubset(featureHairEdgeId, e.end());
        }

        if( hairEdgeType_[heI] & BOUNDARY)
        {
            mesh_.addPointToSubset(bndHairEdgeId, e.start());
            mesh_.addPointToSubset(bndHairEdgeId, e.end());
        }

        if( hairEdgeType_[heI] & ATCORNER)
        {
            mesh_.addPointToSubset(cornerHairEdgeId, e.start());
            mesh_.addPointToSubset(cornerHairEdgeId, e.end());
        }

        if( hairEdgeType_[heI] & ATEDGE)
        {
            mesh_.addPointToSubset(hairEdgeAtEdgeId, e.start());
            mesh_.addPointToSubset(hairEdgeAtEdgeId, e.end());
        }
    }

    mesh_.write();
    # endif
}

bool boundaryLayerOptimisation::optimiseLayersAtExittingFaces()
{
    bool modified(false);

    //- find edge points inside the mesh with more than one hair edge
    //- attached to it
    labelList nEdgesAtPoint(mesh_.points().size(), 0);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges_, heI)
    {
        # ifdef USE_OMP
        # pragma omp atomic
        # endif
        ++nEdgesAtPoint[hairEdges_[heI].end()];
    }

    //- find the points with more than one hair edge which was modified
    //- in the previous procedure
    boolList thinnedPoints(mesh_.points().size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(thinnedHairEdge_, heI)
    {
        if
        (
            thinnedHairEdge_[heI] &&
            (nEdgesAtPoint[hairEdges_[heI].end()] > 2)
        )
        {
            modified = true;
            thinnedPoints[hairEdges_[heI].end()] = true;
        }
    }

    reduce(modified, maxOp<bool>());

    if( !modified )
        return false;

    Info << "Hair edges at exitting faces shall "
         << "be modified due to inner constraints" << endl;

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayerOptimisation::optimiseLayer()
{
    //- create surface smoother
    meshSurfaceOptimizer surfOpt(meshSurface());

    //- lock exitting faces and feature edges
    labelLongList lockedFaces;
    forAll(isExitFace_, bfI)
        if( isExitFace_[bfI] )
            lockedFaces.append(bfI);
    surfOpt.lockBoundaryFaces(lockedFaces);
    surfOpt.lockFeatureEdges();

    label nIter(0);
    do
    {
        thinnedHairEdge_ = false;

        //- calculate normals at the boundary
        optimiseHairNormalsAtTheBoundary();

        //- smoothing thickness variation of boundary hairs
        optimiseThicknessVariation(BOUNDARY);

        if( true )
        {
            meshSurfaceEngineModifier bMod(meshSurface());
            bMod.updateGeometry();

            surfOpt.optimizeSurface(2);
            bMod.updateGeometry();
        }

        # ifdef DEBUGLayer
        label counter(0);
        forAll(thinnedHairEdge_, heI)
            if( thinnedHairEdge_[heI] )
                ++counter;
        reduce(counter, sumOp<label>());
        Info << "Thinned " << counter << " bnd hair edges" << endl;
        # endif

        //- optimise normals inside the mesh
        optimiseHairNormalsInside();

        //- optimise thickness variation inside the mesh
        optimiseThicknessVariation(INSIDE);

        # ifdef DEBUGLayer
        label intCounter = 0;
        forAll(thinnedHairEdge_, heI)
            if( thinnedHairEdge_[heI] )
                ++intCounter;
        Info << "Thinned " << (intCounter - counter)
             << " inner hair edges" << endl;
        # endif
    } while( optimiseLayersAtExittingFaces() && (++nIter < maxNumIterations_) );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

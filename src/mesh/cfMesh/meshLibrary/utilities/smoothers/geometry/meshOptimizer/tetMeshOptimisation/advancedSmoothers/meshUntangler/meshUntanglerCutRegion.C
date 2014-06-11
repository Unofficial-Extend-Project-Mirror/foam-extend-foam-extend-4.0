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
#include "meshUntangler.H"
#include "plane.H"
#include "primitiveMesh.H"
#include "polyMeshGenModifier.H"
#include "sortEdgesIntoChains.H"

//#define DEBUGSmooth

#ifdef DEBUGSmooth
#include "Time.H"
#include "objectRegistry.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshUntangler::cutRegion::createInitialConfiguration
(
    const boundBox& bb
)
{
    pointsPtr_ = new DynList<point, 64>();
    DynList<point, 64>& bVertices = *pointsPtr_;
    edgesPtr_ = new DynList<edge, 128>();
    DynList<edge, 128>& bEdges = *edgesPtr_;
    facesPtr_ = new DynList<DynList<label, 8>, 64>();
    DynList<DynList<label, 8>, 64>& bFaces = *facesPtr_;

    //- set vertices
    const point c = (bb.max() + bb.min()) / 2.0;
    const point vec = (bb.max() - bb.min()) / 2.0;
    bVertices.append
    (
        point(c.x() - vec.x(), c.y() - vec.y(), c.z() - vec.z())
    );
    bVertices.append
    (
        point(c.x() + vec.x(), c.y() - vec.y(), c.z() - vec.z())
    );
    bVertices.append
    (
        point(c.x() + vec.x(), c.y() + vec.y(), c.z() - vec.z())
    );
    bVertices.append
    (
        point(c.x() - vec.x(), c.y() + vec.y(), c.z() - vec.z())
    );
    bVertices.append
    (
        point(c.x() - vec.x(), c.y() - vec.y(), c.z() + vec.z())
    );
    bVertices.append
    (
        point(c.x() + vec.x(), c.y() - vec.y(), c.z() + vec.z())
    );
    bVertices.append
    (
        point(c.x() + vec.x(), c.y() + vec.y(), c.z() + vec.z())
    );
    bVertices.append
    (
        point(c.x() - vec.x(), c.y() + vec.y(), c.z() + vec.z())
    );

    //- set edges

    //- edges in x direction
    bEdges.append(edge(0, 1));
    bEdges.append(edge(3, 2));
    bEdges.append(edge(7, 6));
    bEdges.append(edge(4, 5));

    //- edges in y direction
    bEdges.append(edge(1, 2));
    bEdges.append(edge(0, 3));
    bEdges.append(edge(4, 7));
    bEdges.append(edge(5, 6));

    //- edges in z direction
    bEdges.append(edge(0, 4));
    bEdges.append(edge(1, 5));
    bEdges.append(edge(2, 6));
    bEdges.append(edge(3, 7));

    //- set faces
    DynList<label, 8> f;
    f.setSize(4);

    //- faces in x direction
    f[0] = 5;
    f[1] = 11;
    f[2] = 6;
    f[3] = 8;
    bFaces.append(f);
    f[0] = 4;
    f[1] = 10;
    f[2] = 7;
    f[3] = 9;
    bFaces.append(f);
    //- faces in y direction
    f[0] = 0;
    f[1] = 8;
    f[2] = 3;
    f[3] = 9;
    bFaces.append(f);
    f[0] = 1;
    f[1] = 11;
    f[2] = 2;
    f[3] = 10;
    bFaces.append(f);
    //- faces in z direction
    f[0] = 0;
    f[1] = 4;
    f[2] = 1;
    f[3] = 5;
    bFaces.append(f);
    f[0] = 3;
    f[1] = 7;
    f[2] = 2;
    f[3] = 6;
    bFaces.append(f);

    # ifdef DEBUGSmooth
    Info << "Original vertices " << *pointsPtr_ << endl;
    Info << "Original edges " << *edgesPtr_ << endl;
    Info << "Original faces " << *facesPtr_ << endl;
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

meshUntangler::cutRegion::cutRegion(const boundBox& bb)
:
    pointsPtr_(NULL),
    edgesPtr_(NULL),
    facesPtr_(NULL),
    cPtsPtr_(NULL),
    cEdgesPtr_(NULL),
    cFacesPtr_(NULL),
    newVertexLabel_(),
    vertexDistance_(),
    vertexTypes_(),
    newEdgeLabel_(),
    origNumVertices_(),
    tol_(SMALL * bb.mag()),
    valid_(true)
{
    createInitialConfiguration(bb);
}

meshUntangler::cutRegion::~cutRegion()
{
    deleteDemandDrivenData(pointsPtr_);
    deleteDemandDrivenData(edgesPtr_);
    deleteDemandDrivenData(facesPtr_);
    deleteDemandDrivenData(cPtsPtr_);
    deleteDemandDrivenData(cEdgesPtr_);
    deleteDemandDrivenData(cFacesPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshUntangler::cutRegion::planeCut(const plane& plane)
{
    if( !valid_ )
        return;

    # ifdef DEBUGSmooth
    if( cFacesPtr_ || cPtsPtr_ || cEdgesPtr_ )
    {
        FatalErrorIn
        (
            "void meshUntangler::"
            "cutRegion::planeCut(const plane& plane)"
        ) << "Pointers should not be allocated!" << abort(FatalError);
    }

    Foam::Time runTime
    (
        Foam::Time::controlDictName,
        "../.",
        "testSmoothing"
    );

    objectRegistry oR(runTime);

    polyMeshGen pmg
    (
        oR
    );
    this->createPolyMeshFromRegion(pmg);
    # endif

    if( findNewVertices(plane) )
    {
        findNewEdges();

        findNewFaces();

        if( !valid_ ) return;

        deleteDemandDrivenData(pointsPtr_);
        pointsPtr_ = cPtsPtr_;
        cPtsPtr_ = NULL;

        deleteDemandDrivenData(edgesPtr_);
        edgesPtr_ = cEdgesPtr_;
        cEdgesPtr_ = NULL;

        deleteDemandDrivenData(facesPtr_);
        facesPtr_ = cFacesPtr_;
        cFacesPtr_ = NULL;
    }
}

void meshUntangler::cutRegion::createPolyMeshFromRegion
(
    polyMeshGen& mesh
) const
{
    polyMeshGenModifier meshModifier(mesh);
    pointFieldPMG& points = meshModifier.pointsAccess();
    points.setSize(pointsPtr_->size());
    forAll(points, pI)
        points[pI] = (*pointsPtr_)[pI];

    faceListPMG& faces = meshModifier.facesAccess();
    cellListPMG& cells = meshModifier.cellsAccess();
    cells.setSize(1);
    cells[0].setSize(facesPtr_->size());
    faces.setSize(facesPtr_->size());

    const DynList<edge, 128>& edges = *edgesPtr_;
    const DynList<DynList<label, 8>, 64>& fcs = *facesPtr_;
    forAll(faces, fI)
    {
        DynList<edge> fEdges;
        const DynList<label, 8>& f = fcs[fI];
        forAll(f, eI)
            fEdges.append(edges[f[eI]]);

        Info << "Edges forming face " << fI << " are " << fEdges << endl;
        labelListList sf = sortEdgesIntoChains(fEdges).sortedChains();
        if( sf.size() != 1 )
            FatalErrorIn
            (
                "void meshOptimizer::meshUntangler::"
                "cutRegion::createPolyMeshFromRegion(polyMesgGen&)"
            ) << "More than one face created!" << abort(FatalError);

        faces[fI] = face(sf[0]);
        cells[0][fI] = fI;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

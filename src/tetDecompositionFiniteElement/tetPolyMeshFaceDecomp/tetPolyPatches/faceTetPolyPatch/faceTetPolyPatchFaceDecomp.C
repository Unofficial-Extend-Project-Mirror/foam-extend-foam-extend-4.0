/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "faceTetPolyPatchFaceDecomp.H"
#include "tetPolyBoundaryMeshFaceDecomp.H"
#include "tetPolyMeshFaceDecomp.H"
#include "demandDrivenData.H"
#include "boolList.H"
#include "SubField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(faceTetPolyPatchFaceDecomp, 0);
defineRunTimeSelectionTable(faceTetPolyPatchFaceDecomp, polyPatch);

addToRunTimeSelectionTable
(
    faceTetPolyPatchFaceDecomp,
    faceTetPolyPatchFaceDecomp,
    polyPatch
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

labelList faceTetPolyPatchFaceDecomp::calcLocalEdgesIndices
(
    const primitivePatch& p
) const
{
    if (debug)
    {
        Info<< "labelList faceTetPolyPatchFaceDecomp::calcLocalEdgesIndices("
            << "const primitivePatch& p ) const : " << endl
            << "calculating local edge indices"
            << endl;
    }

    // Count number of edges in the patch

    // Get reference to the mesh
    const tetPolyMeshFaceDecomp& mesh = boundaryMesh().mesh();

    // Get reference to edges of the primitive patch
    const edgeList& patchEdges = p.edges();

    // Get reference to faces of the primitive patch
    const faceList& patchFaces = p;

    // Get reference to mesh points of the primitive patch
    const labelList& patchMeshPoints = p.meshPoints();

    // edges of the polyPatch
    label nEdgesInPatch = patchEdges.size();

    // diagonal edges across faces
    forAll (patchFaces, faceI)
    {
        nEdgesInPatch += patchFaces[faceI].size();
    }

    // Prepare result
    labelList localEdgeInd(nEdgesInPatch, -1);

    label nEdges = 0;

    const lduAddressing& lduAddr = mesh.lduAddr();

    // First do the edges of the primitive patch
    forAll (patchEdges, edgeI)
    {
        localEdgeInd[nEdges] =
            lduAddr.triIndex
            (
                patchMeshPoints[patchEdges[edgeI].start()],
                patchMeshPoints[patchEdges[edgeI].end()]
            );

        nEdges++;
    }

    // Now do the face-internal edges

    // Calculate the offset to the first face centre
    label offset = mesh.faceOffset() + patch().start();

    forAll (patchFaces, faceI)
    {
        const face& curFace = patchFaces[faceI];

        forAll (curFace, pointI)
        {
            localEdgeInd[nEdges] =
                lduAddr.triIndex(curFace[pointI], offset + faceI);

            nEdges++;
        }
    }

#   ifdef DEBUGtetFemMatrix
    if (min(localEdgeInd) < 0)
    {
        FatalErrorIn
        (
            "void faceTetPolyPatchFaceDecomp::calcLocalEdgesIndices() const"
        )   << "Problem in local edge addressing"
            << abort(FatalError);
    }
#   endif

    if (debug)
    {
        Info<< "labelList faceTetPolyPatchFaceDecomp::calcLocalEdgesIndices("
            << "const primitivePatch& p ) const : " << endl
            << "finished calculating local edge indices"
            << endl;
    }

    return localEdgeInd;
}


labelList faceTetPolyPatchFaceDecomp::calcCutEdgeIndices
(
    const primitivePatch& p
) const
{
    // Make a list over all edges in the mesh.  Mark the ones that are local
    // to the patch and then collect the rest

    if (debug)
    {
        Info<< "void faceTetPolyPatchFaceDecomp::calcCutEdgeIndices() const : "
            << endl << "calculating cut edge indices"
            << endl;
    }

    boolList isLocal(boundaryMesh().mesh().nEdges(), false);

    // get reference to local edge indices
    const labelList& localEdges = localEdgeIndices();

    forAll (localEdges, edgeI)
    {
        isLocal[localEdges[edgeI]] = true;
    }

    // Count the maximum number of edges coming from the patch
    label maxEdgesOnPatch = 0;

    const tetPolyMeshFaceDecomp& mesh = boundaryMesh().mesh();

    const labelList& mp = meshPoints();

    forAll (mp, pointI)
    {
        maxEdgesOnPatch += mesh.nEdgesForPoint(mp[pointI]);
    }

    // Prepare the result
    labelList cutEdgeInd(maxEdgesOnPatch, -1);

    label nCutEdgeInd = 0;

    // Go through all the local points and get all the edges coming
    // from that point.  Check if the edge has been marked as local;
    // if not, add it to the list of cut edges.
    forAll (mp, pointI)
    {
        labelList curEdges = mesh.edgesForPoint(mp[pointI]);

        forAll (curEdges, edgeI)
        {
            if (!isLocal[curEdges[edgeI]])
            {
                cutEdgeInd[nCutEdgeInd] = curEdges[edgeI];
                nCutEdgeInd++;
            }
        }
    }

    // Reset the size of the edge list
    cutEdgeInd.setSize(nCutEdgeInd);

    if (debug)
    {
        Info<< "void faceTetPolyPatchFaceDecomp::calcCutEdgeIndices() const : "
            << endl << "finished calculating cut edge indices"
            << endl;
    }

    return cutEdgeInd;
}


labelList faceTetPolyPatchFaceDecomp::calcMeshPoints
(
    const primitivePatch& p
) const
{
    if (debug)
    {
        Info<< "faceTetPolyPatchFaceDecomp::calcMeshPoints() : " << endl
            << "calculating mesh points"
            << endl;
    }

    // Prepare the result
    labelList mp(size(), -1);

    // insert the vertices first
    const labelList& polyPatchMeshPoints = p.meshPoints();

    label nPoints = 0;

    forAll (polyPatchMeshPoints, pointI)
    {
        mp[nPoints] = polyPatchMeshPoints[pointI];

        nPoints++;
    }

    // insert faces using the offset
    const label faceOffset = boundaryMesh().mesh().faceOffset();

    const label polyPatchStart = patch().start() + faceOffset;
    const label polyPatchStartPlusSize = polyPatchStart + patch().size();

    for (label i = polyPatchStart; i < polyPatchStartPlusSize; i++)
    {
        mp[nPoints] = i;

        nPoints++;
    }

    if (debug)
    {
        Info<< "faceTetPolyPatchFaceDecomp::calcMeshPoints() : " << endl
            << "finished calculating mesh points"
            << endl;
    }

    return mp;
}


void faceTetPolyPatchFaceDecomp::calcLocalPoints() const
{
    if (debug)
    {
        Info<< "faceTetPolyPatchFaceDecomp::calcLocalPoints() : " << endl
            << "calculating local points"
            << endl;
    }

    if (localPointsPtr_)
    {
        FatalErrorIn
        (
            "void faceTetPolyPatchFaceDecomp::calcLocalPoints() const"
        )   << "localPointsPtr_ already allocated"
            << abort(FatalError);
    }

    // Prepare the result
    localPointsPtr_ = new pointField(size());
    pointField& lp = *localPointsPtr_;

    // insert the vertices first
    const pointField& polyPatchLocalPoints = patch().localPoints();

    label nPoints = 0;

    forAll (polyPatchLocalPoints, pointI)
    {
        lp[nPoints] = polyPatchLocalPoints[pointI];

        nPoints++;
    }

    // insert faces
    const vectorField::subField polyPatchFaceCentres = patch().faceCentres();

    forAll (polyPatchFaceCentres, faceI)
    {
        lp[nPoints] = polyPatchFaceCentres[faceI];

        nPoints++;
    }

    if (debug)
    {
        Info<< "faceTetPolyPatchFaceDecomp::calcLocalPoints() : " << endl
            << "finished calculating local points"
            << endl;
    }
}


void faceTetPolyPatchFaceDecomp::calcPointNormals() const
{
    if (debug)
    {
        Info<< "faceTetPolyPatchFaceDecomp::calcPointNormals() : " << endl
            << "calculating point normals"
            << endl;
    }

    if (pointNormalsPtr_)
    {
        FatalErrorIn
        (
            "void faceTetPolyPatchFaceDecomp::calcPointNormals() const"
        )   << "pointNormalsPtr_ already allocated"
            << abort(FatalError);
    }

    // Prepare the result
    pointNormalsPtr_ = new vectorField(size());
    vectorField& pn = *pointNormalsPtr_;

    // insert the vertices first
    const vectorField& polyPatchPointNormals = patch().pointNormals();

    label nPoints = 0;

    forAll (polyPatchPointNormals, pointI)
    {
        pn[nPoints] = polyPatchPointNormals[pointI];

        nPoints++;
    }

    // insert faces
    vectorField polyPatchFaceNormals(patch().faceAreas());

    polyPatchFaceNormals /= mag(polyPatchFaceNormals);

    forAll (polyPatchFaceNormals, faceI)
    {
        pn[nPoints] = polyPatchFaceNormals[faceI];

        nPoints++;
    }

    if (debug)
    {
        Info<< "faceTetPolyPatchFaceDecomp::calcPointNormals() : " << endl
            << "finished calculating point normals"
            << endl;
    }
}


void faceTetPolyPatchFaceDecomp::clearOut()
{
    deleteDemandDrivenData(meshPointsPtr_);
    deleteDemandDrivenData(localPointsPtr_);
    deleteDemandDrivenData(pointNormalsPtr_);
    deleteDemandDrivenData(localEdgeIndicesPtr_);
    deleteDemandDrivenData(cutEdgeIndicesPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyPatch
faceTetPolyPatchFaceDecomp::faceTetPolyPatchFaceDecomp
(
    const polyPatch& p,
    const tetPolyBoundaryMeshFaceDecomp& bm
)
:
    tetPolyPatchFaceDecomp(bm),
    boundaryIndex_(p.index()),
    size_(p.meshPoints().size() + p.size()),
    meshPointsPtr_(NULL),
    localPointsPtr_(NULL),
    pointNormalsPtr_(NULL),
    localEdgeIndicesPtr_(NULL),
    cutEdgeIndicesPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

faceTetPolyPatchFaceDecomp::~faceTetPolyPatchFaceDecomp()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const polyPatch& faceTetPolyPatchFaceDecomp::patch() const
{
    return boundaryMesh().mesh()().polyMesh::boundaryMesh()[index()];
}


const labelList& faceTetPolyPatchFaceDecomp::meshPoints() const
{
    if (!meshPointsPtr_)
    {
        meshPointsPtr_ = new labelList(calcMeshPoints(patch()));
    }

    return *meshPointsPtr_;
}


const pointField& faceTetPolyPatchFaceDecomp::localPoints() const
{
    if (!localPointsPtr_)
    {
        calcLocalPoints();
    }
    else
    {
        // Check if the mesh has moved since the last access
        const pointField& polyPatchLocalPoints = patch().localPoints();

        vectorField::subField tetPatchPointSlice
        (
            *localPointsPtr_,
            polyPatchLocalPoints.size()
        );

        if
        (
            sum(mag(polyPatchLocalPoints - tetPatchPointSlice))
          > polyPatchLocalPoints.size()*SMALL
        )
        {
            deleteDemandDrivenData(localPointsPtr_);

            calcLocalPoints();
        }
    }

    return *localPointsPtr_;
}


const vectorField& faceTetPolyPatchFaceDecomp::pointNormals() const
{
    if (!pointNormalsPtr_)
    {
        calcPointNormals();
    }
    else
    {
        // Check if the mesh has moved since the last access
        const vectorField& polyPatchPointNormals = patch().pointNormals();

        vectorField::subField tetPatchNormalSlice
        (
            *pointNormalsPtr_,
            polyPatchPointNormals.size()
        );

        if
        (
            sum(mag(polyPatchPointNormals - tetPatchNormalSlice))
          > polyPatchPointNormals.size()*SMALL
        )
        {
            deleteDemandDrivenData(pointNormalsPtr_);

           calcPointNormals();
        }
    }

    return *pointNormalsPtr_;
}


triFaceList faceTetPolyPatchFaceDecomp::faceTriangles
(
    const label faceID
) const
{
    const face& f =
        boundaryMesh().mesh()()
            .boundaryMesh()[index()].localFaces()[faceID];

    // Create a list of triangles to keep the triangles that
    // have already been added
    triFaceList result(f.size());

    const label offset = patch().meshPoints().size();

    forAll (f, triI)
    {
        result[triI] =
            triFace
            (
                f[triI],
                f.nextLabel(triI),
                offset + faceID
            );
    }

    return result;
}


faceList faceTetPolyPatchFaceDecomp::triFaces() const
{
    const faceList& mf =
        boundaryMesh().mesh()().boundaryMesh()[index()];

    // Count the number of faces
    label nTriFaces = 0;

    forAll (mf, faceI)
    {
        nTriFaces += mf[faceI].size();
    }

    faceList result(nTriFaces);

    // Reset the counter for re-use
    nTriFaces = 0;

    face triangleFace(3);
    const label faceOffset =
        boundaryMesh().mesh().faceOffset() + patch().start();

    forAll (mf, faceI)
    {
        const face& f = mf[faceI];

        forAll (f, triI)
        {
            triangleFace[0] = f[triI];
            triangleFace[1] = f.nextLabel(triI);
            triangleFace[2] = faceOffset + faceI;

            result[nTriFaces] = triangleFace;
            nTriFaces++;
        }
    }

    return result;
}


const labelList& faceTetPolyPatchFaceDecomp::localEdgeIndices() const
{
    if (!localEdgeIndicesPtr_)
    {
        localEdgeIndicesPtr_ = new labelList(calcLocalEdgesIndices(patch()));
    }

    return *localEdgeIndicesPtr_;
}


const labelList& faceTetPolyPatchFaceDecomp::cutEdgeIndices() const
{
    if (!cutEdgeIndicesPtr_)
    {
        cutEdgeIndicesPtr_ = new labelList(calcCutEdgeIndices(patch()));
    }

    return *cutEdgeIndicesPtr_;
}


void faceTetPolyPatchFaceDecomp::updateMesh()
{
    clearOut();

    size_ = patch().meshPoints().size() + patch().size();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

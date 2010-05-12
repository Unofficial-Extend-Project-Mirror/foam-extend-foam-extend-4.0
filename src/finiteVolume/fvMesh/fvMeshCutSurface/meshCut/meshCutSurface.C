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

Description

    Cuts tets resulting from tet decomposition. See Userguide,
    'Cell shapes in the traditional mesh description'

        f = face
        fp = vertex of face
        f[fp] = mesh vertex
        fc = face centre
        cc = cell centre

    Caters for both types of decomposition:

    cellDecomp: introduction of cellCentre. Faces get triangulated by
    introducing diagonals.

    faceDecomp: introduction of cellCentre and faceCentres.


    Tet given by its 4 vertices:

            faceDecomp                          cellDecomp
            ----------                          ----------
        0:  f[fp]                               f[fp]
        1:  fc                                  f[0]
        2:  f[fp+1]                             f[fp+1]
        3:  cc                                  cc

    Tet given by its 6 edges:

            faceDecomp                          cellDecomp
            ----------                          ----------
        0:  fp-fc        (face edge)            fp-f0   (diagonal edge)
        1:  fp-fp+1      ('real', mesh edge)    fp-fp+1 (mesh edge)
        2:  fp-cc        (pyramid edge)         fp-cc   (pyramid edge)
        3:  fc-cc        (centre edge)          f0-cc   (pyramid edge)
        4:  fc-fp+1      (face edge)            f0-fp+1 (diagonal edge)
        5:  fp+1 - cc    (pyramid edge)         fp+1-cc (pyramid edge)

\*---------------------------------------------------------------------------*/

#include "meshCutSurface.H"
#include "faceDecompCuts.H"
#include "cellDecompCuts.H"
#include "primitiveMesh.H"
#include "cellAddressing.H"
#include "Map.H"
#include "boolList.H"
#include "cellModeller.H"
#include "orientedSurface.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Symbolic names of tet vertices in faceDecomposition (replace FC by f[0] for
// cellDecomposition)
Foam::label Foam::meshCutSurface::FP0 = 0;
Foam::label Foam::meshCutSurface::FC  = 1;
Foam::label Foam::meshCutSurface::FP1 = 2;
Foam::label Foam::meshCutSurface::CC  = 3;

// Symbolic names of tet edges in faceDecomposition
Foam::label Foam::meshCutSurface::FP0_FC  = 0;
Foam::label Foam::meshCutSurface::FP0_FP1 = 1;
Foam::label Foam::meshCutSurface::FP0_CC  = 2;
Foam::label Foam::meshCutSurface::FC_CC   = 3;
Foam::label Foam::meshCutSurface::FC_FP1  = 4;
Foam::label Foam::meshCutSurface::FP1_CC  = 5;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class Hash>
Foam::label Foam::meshCutSurface::find
(
    const HashTable<label, T, Hash>& table,
    const T& key
)
{
    typename HashTable<label, T, Hash>::const_iterator iter = table.find(key);

    if (iter == table.end())
    {
        return -1;
    }
    else
    {
        return iter();
    }
}


Foam::label Foam::meshCutSurface::findEdge
(
    const edgeList& edges,
    const labelList& edgeLabels,
    const edge& e
)
{
    forAll(edgeLabels, edgeLabelI)
    {
        label edgeI = edgeLabels[edgeLabelI];

        if (edges[edgeI] == e)
        {
            return edgeI;
        }
    }

    forAll(edgeLabels, edgeLabelI)
    {
        label edgeI = edgeLabels[edgeLabelI];

        Pout<< "edge:" << edgeI << " vertices:" << edges[edgeI] << endl;
    }

    FatalErrorIn
    (
        "meshCutSurface::findEdge(const edgeList&, const labelList&, const edge&)"
    )   << "Can not find edge " << e << " among edge labels " << edgeLabels
        << exit(FatalError);

    return -1;
}


Foam::label Foam::meshCutSurface::findCutEdge
(
    const cellAddressing& model,
    const labelList& tetEdgeCuts,
    const label faceI,
    const label cutEdgeI
)
{
    const labelList& fEdges = model.faceEdges()[faceI];

    forAll(fEdges, fEdgeI)
    {
       label edgeI = fEdges[fEdgeI];

       if ((edgeI != cutEdgeI) && (tetEdgeCuts[edgeI] != -1))
       {
            return edgeI;
       }
    }
    
    FatalErrorIn
    (
        "meshCutSurface::findCutEdge(const cellAddressing&, const labelList&"
        ", const label, const label)"
    )   << "Problem: face " << faceI << " consists of edges " << fEdges
        << " of which only one is cut"
        << exit(FatalError);

    return -1;
}


Foam::triSurface Foam::meshCutSurface::removeDuplicates(const triSurface& surf)
{
    boolList includeTri(surf.size(), true);

    forAll(surf, triI)
    {
        const labelledTri& t = surf[triI];

        const labelList& nbs = surf.faceFaces()[triI];

        forAll(nbs, nbI)
        {
            if (nbs[nbI] <= triI)
            {
                // lower numbered faces already checked
                continue;
            }

            const labelledTri& nb = surf[nbs[nbI]];

            if
            (
                ((t[0] == nb[0]) || (t[0] == nb[1]) || (t[0] == nb[2]))
             && ((t[1] == nb[0]) || (t[1] == nb[1]) || (t[1] == nb[2]))
             && ((t[2] == nb[0]) || (t[2] == nb[1]) || (t[2] == nb[2]))
            )
            {
                includeTri[nbs[nbI]] = false;
            }
        }
    }

    labelList pointMap, faceMap;

    return surf.subsetMesh(includeTri, pointMap, faceMap);
}


void Foam::meshCutSurface::cutTet
(
    const cellAddressing& model,
    const labelList& tetVertCuts,
    const labelList& tetEdgeCuts,
    const label cellI,

    List<labelledTri>& tris,
    label& triI
)
{
    label nVertCuts = 0;
    forAll(tetVertCuts, tetVertI)
    {
        if (tetVertCuts[tetVertI] != -1)
        {
            nVertCuts++;
        }
    }

    label nEdgeCuts = 0;
    forAll(tetEdgeCuts, tetEdgeI)
    {
        if (tetEdgeCuts[tetEdgeI] != -1)
        {
            nEdgeCuts++;
        }
    }

    if (nVertCuts != 0)
    {
        // Cuts through vertices as well as edges. Handle separately.
        cutTetThroughVerts
        (
            model,
            tetVertCuts,
            nVertCuts,
            tetEdgeCuts,
            nEdgeCuts,
            cellI,

            tris,
            triI
        );
    }
    else
    {
        // Cuts through edges. 'Simple' case.
        cutTetThroughEdges
        (
            model,
            tetEdgeCuts,
            nEdgeCuts,
            cellI,

            tris,
            triI
        );
    }
}


void Foam::meshCutSurface::cutTetThroughEdges
(
    const cellAddressing& model,
    const labelList& tetEdgeCuts,
    const label nEdgeCuts,
    const label cellI,

    List<labelledTri>& tris,
    label& triI
)
{
    if (nEdgeCuts < 3)
    {
        //Pout<< "Tet only cut by " << nEdgeCuts << endl;
    }
    else if (nEdgeCuts == 3)
    {
        label triVertI = 0;

        labelledTri& tri = tris[triI++];

        forAll(tetEdgeCuts, tetEdgeI)
        {
            if (tetEdgeCuts[tetEdgeI] != -1)
            {
                tri[triVertI++] = tetEdgeCuts[tetEdgeI];
            }
        }

        //Pout<< "Triangle formed by vertices:" << tri << endl;
        tri.region() = cellI;
    }
    else if (nEdgeCuts == 4)
    {
        // Find ordering of vertices

        label edge0 = -1;
        forAll(tetEdgeCuts, tetEdgeI)
        {
            if (tetEdgeCuts[tetEdgeI] != -1)
            {
                edge0 = tetEdgeI;

                break;
            }
        }

        const labelList& nbFaces = model.edgeFaces()[edge0];

        // Find cut edges among faces using edge0
        label edge1 = findCutEdge(model, tetEdgeCuts, nbFaces[0], edge0);
        label edge2 = findCutEdge(model, tetEdgeCuts, nbFaces[1], edge0);

        // Find the other edge
        label edge3 = -1;
        forAll(tetEdgeCuts, tetEdgeI)
        {
            if
            (
                (tetEdgeCuts[tetEdgeI] != -1)
             && (tetEdgeI != edge0)
             && (tetEdgeI != edge1)
             && (tetEdgeI != edge2)
            )
            {
                // Fourth cut edge
                edge3 = tetEdgeI;

                break;
            }
        }

        //Pout<< "Triangle 1 formed by vertices:" << tri1 << endl;
        tris[triI++] =
            labelledTri
            (
                tetEdgeCuts[edge0],
                tetEdgeCuts[edge1],
                tetEdgeCuts[edge2],
                cellI
            );

        //Pout<< "Triangle 2 formed by vertices:" << tri2 << endl;
        tris[triI++] = 
            labelledTri
            (
                tetEdgeCuts[edge1],
                tetEdgeCuts[edge2],
                tetEdgeCuts[edge3],
                cellI
            );

    }
    else
    {
        //FatalErrorIn("meshCutSurface::cutTet")
        //    << "Cannot handle tets with " << nEdgeCuts << " edges cut"
        //    << exit(FatalError);
    }
}


void Foam::meshCutSurface::cutTetThroughVerts
(
    const cellAddressing& model,
    const labelList& tetVertCuts,
    const label nVertCuts,
    const labelList& tetEdgeCuts,
    const label nEdgeCuts,
    const label cellI,

    List<labelledTri>& tris,
    label& triI
)
{
    if (nVertCuts == 1)
    {
        if (nEdgeCuts == 1)
        {
            return;
        }
        else if (nEdgeCuts != 2)
        {
            // Illegal
            return;
        }
        else    // nEdgeCuts == 2
        {
            // Find cut-edges not connected to cut-vertex.

            // Triangle to fill.
            labelledTri tri;
            label triVertI = 0;

            // Get label of vertex cut
            label cutVertI = -1;

            forAll(tetVertCuts, tetVertI)
            {
                if (tetVertCuts[tetVertI] != -1)
                {
                    tri[triVertI++] = tetVertCuts[tetVertI];

                    cutVertI = tetVertI;

                    break; 
                }
            }

            // Mark all edges connected to the cut vertex
            boolList connectedEdge(6, false);

            const labelList& vEdges = model.pointEdges()[cutVertI];

            forAll(vEdges, vEdgeI)
            {
                connectedEdge[vEdges[vEdgeI]] = true;
            }

            // Go through all edges and collect vertices from cut edges.
            forAll(tetEdgeCuts, tetEdgeI)
            {
                if (tetEdgeCuts[tetEdgeI] != -1)
                {
                    // Cut edge.

                    if (connectedEdge[tetEdgeI])
                    {
                        // Illegal
                        return;
                    }
                    else
                    {
                        // Store edge cut vertex
                        tri[triVertI++] = tetEdgeCuts[tetEdgeI];
                    }
                }
            }


            if (triVertI == 3)
            {
                // Collected all triVerts.
                tri.region() = cellI;
                tris[triI++] = tri;
            }
            else
            {
                // Illegal.
            }
        }
    }
    else if (nVertCuts == 2)
    {
        if (nEdgeCuts == 1)
        {
            // Only interesting if the cut edge is the one opposite the
            // edge connecting the two cut vertices.

            // Triangle to fill.
            labelledTri tri;
            label triVertI = 0;

            // Get label of edge cut
            label edgeI = -1;

            forAll(tetEdgeCuts, tetEdgeI)
            {
                if (tetEdgeCuts[tetEdgeI] != -1)
                {
                    tri[triVertI++] = tetEdgeCuts[tetEdgeI];

                    edgeI = tetEdgeI;

                    break; 
                }
            }


            // Find cut vertices which are not on edge edgeI.
            const edge& e = model.edges()[edgeI];

            forAll(tetVertCuts, tetVertI)
            {
                if
                (
                    (tetVertCuts[tetVertI] != -1)
                 && (tetVertI != e.start())
                 && (tetVertI != e.end())
                )
                {
                    tri[triVertI++] = tetVertCuts[tetVertI];
                }
            }


            if (triVertI == 3)
            {
                // Collected all triVerts.
                //Pout<< "Triangle formed by vertices:" << tri << endl;
                tri.region() = cellI;
                tris[triI++] = tri;
            }
            else
            {
                // Illegal.
            }
        }
        else    // nEdgeCuts != 1
        {
            // Illegal if 2 or more edges cut.
        }   
    }
    else if (nVertCuts == 3)
    {
        if (nEdgeCuts == 0)
        {
            // Cut coincides with face. Form triangle out of cut vertices.

            // Triangle to fill.
            labelledTri tri;
            label triVertI = 0;

            forAll(tetVertCuts, tetVertI)
            {
                if (tetVertCuts[tetVertI] != -1)
                {
                    tri[triVertI++] = tetVertCuts[tetVertI];
                }
            }
            tri.region() = cellI;
            tris[triI++] = tri;
        }
        else    // nEdgeCuts != 0
        {
            // Illegal.
        }
    }
    else // nVertCuts <= 0 || > 3
    {
        // Illegal.
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from edge cuts. Is given edges and position on edge.
//  - cuts.cells(): list of cells affected in any way by cutting.
//  - cuts.meshVerts(): labels of mesh vertices that are cut
//    (cut exactly through point
//  - cuts.faceCentres(): face labels whose faceCentres are cut
//  - cuts.cellCentres(): cell labels whose cellcentres are cut
//  - cuts.meshEdges(): labels of mesh edges that are cut
//  - cuts.pyrEdges(): pyramid edges cut
//  - cuts.centreEdges(): face-cellCentre edges cut
//  - cuts.faceEdges(): face-face edges cut
Foam::meshCutSurface::meshCutSurface(const faceDecompCuts& cuts)
:
    triSurface()
{
    const primitiveMesh& mesh = cuts.mesh();

    // Collect points and create reverse map from edge to new vertex.
    pointField points(cuts.size());
    label pointI = 0;

    //
    // Cuts through existing points
    //

    // Mesh vertices cut
    Map<label> meshVertToVert(cuts.meshVerts().size()*2);

    forAll(cuts.meshVerts(), cutVertI)
    {
        points[pointI++] = mesh.points()[cuts.meshVerts()[cutVertI]];

        meshVertToVert.insert(cuts.meshVerts()[cutVertI], pointI - 1);
    }

    // Face centres cut
    Map<label> faceCentreToVert(cuts.meshFaceCentres().size()*2);

    forAll(cuts.meshFaceCentres(), cutFaceCentreI)
    {
        label faceI = cuts.meshFaceCentres()[cutFaceCentreI];

        points[pointI++] = mesh.faceCentres()[faceI];

        faceCentreToVert.insert(faceI, pointI - 1);
    }

    // Cell centres cut
    Map<label> cellCentreToVert(cuts.meshCellCentres().size()*2);

    forAll(cuts.meshCellCentres(), cutCellCentreI)
    {
        label cellI = cuts.meshCellCentres()[cutCellCentreI];

        points[pointI++] = mesh.cellCentres()[cellI];

        cellCentreToVert.insert(cellI, pointI - 1);
    }


    //
    // Cuts through edges
    //

    // Cuts through existing (mesh) edges
    Map<label> meshEdgeToVert(cuts.meshEdges().size()*2);

    forAll(cuts.meshEdges(), meshCutSurfaceEdgeI)
    {
        label edgeI = cuts.meshEdges()[meshCutSurfaceEdgeI];

        const edge& e = mesh.edges()[edgeI];

        scalar weight = cuts.meshEdgeWeights()[meshCutSurfaceEdgeI];

        points[pointI++] =
            (1-weight)*mesh.points()[e.start()]
          + weight*mesh.points()[e.end()];

        meshEdgeToVert.insert(edgeI, pointI - 1);
    }

    // Cuts on tetDecomposition: Pyramid edges
    HashTable<label, pyramidEdge, pyramidEdge::pyramidEdgeHash>
        pyrEdgeToVert(cuts.pyrEdges().size()*2);

    forAll(cuts.pyrEdges(), edgeI)
    {
        const pyramidEdge& e = cuts.pyrEdges()[edgeI];

        points[pointI++] = e.coord(mesh, cuts.pyrEdgeWeights()[edgeI]);

        pyrEdgeToVert.insert(e, pointI - 1);
    }

    // Cuts on tetDecomposition: Face-Cell centre edges
    HashTable<label, centreEdge, centreEdge::centreEdgeHash>
        centreEdgeToVert(cuts.centreEdges().size()*2);

    forAll(cuts.centreEdges(), edgeI)
    {
        const centreEdge& e = cuts.centreEdges()[edgeI];

        points[pointI++] = e.coord(mesh, cuts.centreEdgeWeights()[edgeI]);

        centreEdgeToVert.insert(e, pointI - 1);
    }

    // Cuts on tetDecomposition: Face edges
    HashTable<label, faceEdge, faceEdge::faceEdgeHash>
        faceEdgeToVert(cuts.faceEdges().size()*2);

    forAll(cuts.faceEdges(), edgeI)
    {
        const faceEdge& e = cuts.faceEdges()[edgeI];

        points[pointI++] = e.coord(mesh, cuts.faceEdgeWeights()[edgeI]);

        faceEdgeToVert.insert(e, pointI - 1);
    }

    if (pointI != points.size())
    {
        FatalErrorIn
        (
            "meshCutSurface::meshCutSurface(const faceDecompCuts& cuts)"
        )   << "pointI:" << pointI << "  points.size():" << points.size()
            << abort(FatalError);
    }


    // Generate tet model + additional addressing
    const cellModel& tet = *(cellModeller::lookup("tet"));

    cellAddressing tetPlus(tet);

    // Tet edges cut. -1 if edge not cut, vertex label otherwise.
    labelList tetEdgeCuts(6);

    // Tet vertices cut. -1 if edge not cut, vertex label otherwise.
    labelList tetVertCuts(4);

    List<labelledTri> tris(2*cuts.size());
    label triI = 0;

    forAll(cuts.cells(), cutCellI)
    {
        label cellI = cuts.cells()[cutCellI];

        tetVertCuts[CC] = find(cellCentreToVert, cellI);

        const cell& cellFaces = mesh.cells()[cellI];

        forAll(cellFaces, cellFaceI)
        {
            label faceI = cellFaces[cellFaceI];

            tetVertCuts[FC] = find(faceCentreToVert, faceI);

            bool isOwner = (mesh.faceOwner()[faceI] == cellI);

            tetEdgeCuts[FC_CC] = find
            (
                centreEdgeToVert,
                centreEdge(faceI, isOwner)
            );

            const face& f = mesh.faces()[faceI];

            forAll(f, fp)
            {
                label fp1 = (fp+1) % f.size();

                tetVertCuts[FP0] = find(meshVertToVert, f[fp]);
                tetVertCuts[FP1] = find(meshVertToVert, f[fp1]);


                tetEdgeCuts[FP0_FC] = find
                (
                    faceEdgeToVert,
                    faceEdge(faceI, fp)
                );


                label edgeI = findEdge
                (
                    mesh.edges(),
                    mesh.faceEdges()[faceI],
                    edge(f[fp], f[fp1])
                );
                tetEdgeCuts[FP0_FP1] = find
                (
                    meshEdgeToVert,
                    edgeI
                );

                tetEdgeCuts[FP0_CC] = find
                (
                    pyrEdgeToVert,
                    pyramidEdge(f[fp], cellI)
                );

                tetEdgeCuts[FC_FP1] = find
                (
                    faceEdgeToVert,
                    faceEdge(faceI, fp1)
                );

                tetEdgeCuts[FP1_CC] = find
                (
                    pyrEdgeToVert,
                    pyramidEdge(f[fp1], cellI)
                );


                cutTet
                (
                    tetPlus,
                    tetVertCuts,
                    tetEdgeCuts,
                    cellI,

                    tris,
                    triI
                );
            }
        }
    }

    tris.setSize(triI);

    // Get a point outside all of the surface (over all processors)
    boundBox bb(points, true);

    point outsidePt(1.5*bb.max() - 0.5*bb.min());

    // Orient the surface after removing duplicate faces.
    triSurface::operator=
    (
        orientedSurface
        (
            removeDuplicates
            (
                triSurface(tris, points)
            ),
            outsidePt
        )
    );
}


// Construct from edge cuts. Is given edges and position on edge.
//  - cuts.cells(): list of cells affected in any way by cutting.
//  - cuts.meshVerts(): labels of mesh vertices that are cut (cut exactly through
//    point
//  - cuts.cellCentres(): cell labels whose cellcentres are cut
//  - cuts.meshEdges(): labels of mesh edges that are cut
//  - cuts.pyrEdges(): pyramid edges cut
//  - cuts.diagEdges(): diagonal edges cut
Foam::meshCutSurface::meshCutSurface(const cellDecompCuts& cuts)
:
    triSurface()
{
    const primitiveMesh& mesh = cuts.mesh();

    // Collect points and create reverse map from edge to new vertex.
    pointField points(cuts.size());
    label pointI = 0;

    //
    // Cuts through existing points
    //

    // Mesh vertices cut
    Map<label> meshVertToVert(cuts.meshVerts().size()*2);

    forAll(cuts.meshVerts(), cutVertI)
    {
        points[pointI++] = mesh.points()[cuts.meshVerts()[cutVertI]];

        meshVertToVert.insert(cuts.meshVerts()[cutVertI], pointI - 1);
    }

    // Cell centres cut
    Map<label> cellCentreToVert(cuts.meshCellCentres().size()*2);

    forAll(cuts.meshCellCentres(), cutCellCentreI)
    {
        label cellI = cuts.meshCellCentres()[cutCellCentreI];

        points[pointI++] = mesh.cellCentres()[cellI];

        cellCentreToVert.insert(cellI, pointI - 1);
    }


    //
    // Cuts through edges
    //

    // Cuts through existing (mesh) edges
    Map<label> meshEdgeToVert(cuts.meshEdges().size()*2);

    forAll(cuts.meshEdges(), meshCutSurfaceEdgeI)
    {
        label edgeI = cuts.meshEdges()[meshCutSurfaceEdgeI];

        const edge& e = mesh.edges()[edgeI];

        scalar weight = cuts.meshEdgeWeights()[meshCutSurfaceEdgeI];

        points[pointI++] =
            (1-weight)*mesh.points()[e.start()]
          + weight*mesh.points()[e.end()];

        meshEdgeToVert.insert(edgeI, pointI - 1);
    }

    // Cuts on tetDecomposition: Pyramid edges
    HashTable<label, pyramidEdge, pyramidEdge::pyramidEdgeHash>
        pyrEdgeToVert(cuts.pyrEdges().size()*2);

    forAll(cuts.pyrEdges(), edgeI)
    {
        const pyramidEdge& e = cuts.pyrEdges()[edgeI];

        points[pointI++] = e.coord(mesh, cuts.pyrEdgeWeights()[edgeI]);

        pyrEdgeToVert.insert(e, pointI - 1);
    }

    // Cuts on tetDecomposition: diagonal edges
    HashTable<label, diagonalEdge, diagonalEdge::diagonalEdgeHash>
        diagEdgeToVert(cuts.diagEdges().size()*2);

    forAll(cuts.diagEdges(), edgeI)
    {
        const diagonalEdge& e = cuts.diagEdges()[edgeI];

        points[pointI++] = e.coord(mesh, cuts.diagEdgeWeights()[edgeI]);

        diagEdgeToVert.insert(e, pointI - 1);
    }



    // Generate tet model + additional addressing
    const cellModel& tet = *(cellModeller::lookup("tet"));

    cellAddressing tetPlus(tet);

    // Tet edges cut. -1 if edge not cut, vertex label otherwise.
    labelList tetEdgeCuts(6);

    // Tet vertices cut. -1 if edge not cut, vertex label otherwise.
    labelList tetVertCuts(4);

    List<labelledTri> tris(2*cuts.size());
    label triI = 0;

    forAll(cuts.cells(), cutCellI)
    {
        label cellI = cuts.cells()[cutCellI];

        tetVertCuts[CC] = find(cellCentreToVert, cellI);

        const cell& cellFaces = mesh.cells()[cellI];

        forAll(cellFaces, cellFaceI)
        {
            label faceI = cellFaces[cellFaceI];

            const face& f = mesh.faces()[faceI];

            // Handle f0 separately
            tetVertCuts[FC] = find(meshVertToVert, f[0]);
            tetEdgeCuts[FC_CC] = find
            (
                pyrEdgeToVert,
                pyramidEdge(f[0], cellI)
            );

            for (label fp = 1; fp < f.size()-1; fp++)
            {
                label fp1 = fp+1;

                tetVertCuts[FP0] = find(meshVertToVert, f[fp]);
                tetVertCuts[FP1] = find(meshVertToVert, f[fp1]);

                tetEdgeCuts[FP0_FP1] = find
                (
                    meshEdgeToVert,
                    findEdge
                    (
                        mesh.edges(),
                        mesh.faceEdges()[faceI],
                        edge(f[fp], f[fp1])
                    )
                );

                if (fp == 1)
                {
                    // f[1] to f[0] is normal mesh edge
                    tetEdgeCuts[FP0_FC] = find
                    (
                        meshEdgeToVert,
                        findEdge
                        (
                            mesh.edges(),
                            mesh.faceEdges()[faceI],
                            edge(f[0], f[fp])
                        )
                    );
                }
                else
                {
                    // Diagonal cut
                    tetEdgeCuts[FP0_FC] = find
                    (
                        diagEdgeToVert,
                        diagonalEdge(faceI, 0, fp)
                    );
                }

                if (fp1 == f.size()-1)
                {
                    // Last face vertex to f[0] -> normal mesh edge
                    tetEdgeCuts[FC_FP1] = find
                    (
                        meshEdgeToVert,
                        findEdge
                        (
                            mesh.edges(),
                            mesh.faceEdges()[faceI],
                            edge(f[0], f[fp1])
                        )
                    );
                }
                else
                {
                    tetEdgeCuts[FC_FP1] = find
                    (
                        diagEdgeToVert,
                        diagonalEdge(faceI, 0, fp1)
                    );
                }

                tetEdgeCuts[FP0_CC] = find
                (
                    pyrEdgeToVert,
                    pyramidEdge(f[fp], cellI)
                );

                tetEdgeCuts[FP1_CC] = find
                (
                    pyrEdgeToVert,
                    pyramidEdge(f[fp1], cellI)
                );


                cutTet
                (
                    tetPlus,
                    tetVertCuts,
                    tetEdgeCuts,
                    cellI,

                    tris,
                    triI
                );
            }
        }
    }

    tris.setSize(triI);

    // Get a point outside all of the surface (over all processors)
    boundBox bb(points, true);

    point outsidePt(1.5*bb.max() - 0.5*bb.min());

    // Orient the surface after removing duplicate faces.
    triSurface::operator=
    (
        orientedSurface
        (
            removeDuplicates
            (
                triSurface(tris, points)
            ),
            outsidePt
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //

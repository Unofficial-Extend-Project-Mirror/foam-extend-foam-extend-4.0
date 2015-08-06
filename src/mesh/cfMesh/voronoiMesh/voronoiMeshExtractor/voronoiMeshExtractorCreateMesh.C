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

#include "voronoiMeshExtractor.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGVoronoi

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void voronoiMeshExtractor::createPoints()
{
    const LongList<point>& tetPoints = tetCreator_.tetPoints();
    const LongList<partTet>& tets = tetCreator_.tets();

    pointFieldPMG& points = mesh_.points();
    points.setSize(tets.size());

    # ifdef DEBUGVoronoi
    Info << "Number of tets " << tets.size() << endl;
    # endif

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(tets, tetI)
    {
        points[tetI] = tets[tetI].centroid(tetPoints);

        # ifdef DEBUGVoronoi
        Info << "Point " << tetI << " has coordinates "
            << points[tetI] << endl;
        Info << "Tet of origin " << tetI << " has nodes "
            << tets[tetI] << endl;
        # endif
    }
}

void voronoiMeshExtractor::createPolyMesh()
{
    const VRWGraph& pointEdges = this->pointEdges();
    const VRWGraph& edgeTets = this->edgeTets();
    const boolList& boundaryEdge = this->boundaryEdge();
    const LongList<edge>& edges = this->edges();

    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();
    cellListPMG& cells = meshModifier.cellsAccess();

    //- count the number of cells
    label nCells(0);
    labelList cellLabel(pointEdges.size(), -1);

    forAll(pointEdges, pointI)
    {
        bool create(true);
        forAllRow(pointEdges, pointI, eI)
            if( boundaryEdge[pointEdges(pointI, eI)] )
            {
                create = false;
                break;
            }

        if( !create || (pointEdges.sizeOfRow(pointI) == 0) )
            continue;

        cellLabel[pointI] = nCells;
        ++nCells;
    }

    # ifdef DEBUGVoronoi
    Info << "Number of cells " << nCells << endl;
    # endif

    //- count the number of faces
    label nFaces(0);
    labelList faceLabel(edges.size(), -1);

    forAll(boundaryEdge, edgeI)
    {
        if( boundaryEdge[edgeI] )
            continue;

        const edge& e = edges[edgeI];
        if( cellLabel[e[0]] < 0 && cellLabel[e[1]] < 0 )
            continue;

        faceLabel[edgeI] = nFaces;
        ++nFaces;
    }

    # ifdef DEBUGVoronoi
    Info << "Number of mesh faces " << nFaces << endl;
    # endif

    //- create faces
    Info << "Creating faces " << endl;
    faces.setSize(nFaces);
    labelList nFacesInCell(nCells, 0);

    forAll(edges, edgeI)
    {
        if( faceLabel[edgeI] < 0 )
            continue;

        const label faceI = faceLabel[edgeI];

        face& f = faces[faceI];
        f.setSize(edgeTets.sizeOfRow(edgeI));

        //- fill the faces with the node labels
        forAllRow(edgeTets, edgeI, pI)
            f[pI] = edgeTets(edgeI, pI);

        const edge& e = edges[edgeI];
        const label cOwn = cellLabel[e[0]];
        const label cNei = cellLabel[e[1]];

        if( cOwn < 0 || ((cOwn > cNei) && (cNei != -1)) )
        {
            # ifdef DEBUGVoronoi
            Info << "Reversing face " << cOwn << " " << cNei << endl;
            # endif

            f = f.reverseFace();
        }

        if( cOwn >= 0 )
            ++nFacesInCell[cOwn];
        if( cNei >= 0 )
            ++nFacesInCell[cNei];
    }

    //- create cells
    # ifdef DEBUGVoronoi
    Info << "Setting cell sizes" << endl;
    # endif

    cells.setSize(nCells);
    forAll(nFacesInCell, cellI)
    {
        cells[cellI].setSize(nFacesInCell[cellI]);
        nFacesInCell[cellI] = 0;
    }

    # ifdef DEBUGVoronoi
    Info << "Filling cells" << endl;
    # endif

    forAll(edges, edgeI)
    {
        if( faceLabel[edgeI] < 0 )
            continue;

        const label faceI = faceLabel[edgeI];
        const edge& e = edges[edgeI];
        const label cOwn = cellLabel[e[0]];
        const label cNei = cellLabel[e[1]];

        if( cOwn >= 0 )
            cells[cOwn][nFacesInCell[cOwn]++] = faceI;
        if( cNei >= 0 )
            cells[cNei][nFacesInCell[cNei]++] = faceI;
    }

    # ifdef DEBUGVoronoi
    Info << "Finished generating cells" << endl;
    # endif

    mesh_.clearAddressingData();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

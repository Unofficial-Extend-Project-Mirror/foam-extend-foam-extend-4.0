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

\*----------------------p-----------------------------------------------------*/

#include "decomposeCells.H"
#include "helperFunctions.H"
#include "triFace.H"

//#define DEBUGDecompose

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void decomposeCells::findAddressingForCell
(
    const label cellI,
    DynList<label, 32>& vrt,
    DynList<edge, 64>& edges,
    DynList<DynList<label, 8> >& faceEdges,
    DynList<DynList<label, 2>, 64>& edgeFaces
) const
{
    const cell& c = mesh_.cells()[cellI];

    vrt.clear();
    edges.clear();
    edgeFaces.clear();
    faceEdges.setSize(c.size());

    const faceListPMG& faces = mesh_.faces();
    forAll(faceEdges, feI)
    {
        faceEdges[feI].setSize(faces[c[feI]].size());
        faceEdges[feI] = -1;
    }

    forAll(c, fI)
    {
        const face& f = faces[c[fI]];

        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            bool store(true);
            forAll(vrt, vI)
                if( vrt[vI] == f[eI] )
                {
                    store = false;
                    break;
                }
            if( store )
            {
                vrt.append(f[eI]);
            }

            //- check if the edge alreready exists
            store = true;

            forAll(edges, eJ)
                if( e == edges[eJ] )
                {
                    store = false;
                    faceEdges[fI][eI] = eJ;
                    edgeFaces[eJ].append(fI);
                    break;
                }

            if( store )
            {
                faceEdges[fI][eI] = edges.size();
                DynList<label, 2> ef;
                ef.append(fI);
                edgeFaces.append(ef);
                edges.append(e);
            }
        }
    }

    //- check if the cell is topologically closed
    forAll(edgeFaces, efI)
        if( edgeFaces[efI].size() != 2 )
        {
            forAll(c, fI)
                Info << "Face " << c[fI] << " is " << faces[c[fI]] << endl;

            Info << "Edges " << edges << endl;
            Info << "faceEdges " << faceEdges << endl;
            Info << "edgeFaces " << edgeFaces << endl;
            mesh_.write();
            FatalErrorIn
            (
                "void decomposeCells::findAddressingForCell"
                "(const label, DynList<label, 32>&, DynList<edge, 64>&"
                ", DynList<DynList<label, 8> >&"
                ", DynList<DynList<label, 2>, 64>&) const"
            ) << "Cell " << cellI << " is not topologically closed!"
                << abort(FatalError);
        }
}

label decomposeCells::findTopVertex
(
    const label cellI,
    const DynList<label, 32>& vrt,
    const DynList<edge, 64>& edges,
    const DynList<DynList<label, 2>, 64>& edgeFaces
)
{
    const cell& c = mesh_.cells()[cellI];
    const faceListPMG& faces = mesh_.faces();

    pointFieldPMG& pointsAccess = mesh_.points();

    //- there is no vertex in 3 or more patches
    //- find boundary faces
    label topVertex(-1);

    const labelList cp = c.labels(faces);
    point p(vector::zero);
    forAll(cp, cpI)
        p += pointsAccess[cp[cpI]];
    p /= cp.size();

    topVertex = pointsAccess.size();
    pointsAccess.append(p);

    # ifdef DEBUGDecompose
    Info << "Top vertex is " << topVertex << endl;
    # endif

    return topVertex;
}

void decomposeCells::decomposeCellIntoPyramids(const label cellI)
{
    const cellListPMG& cells = mesh_.cells();
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();

    const cell& c = cells[cellI];

    # ifdef DEBUGDecompose
    Info << "Starting decomposing cell " << cellI << endl;
    Info << "Cell consists of faces " << c << endl;
    forAll(c, fI)
        Info << "Face " << c[fI] << " is " << faces[c[fI]] << endl;
    # endif

    //- calculate edges, faceEdges and edgeFaces addressing
    DynList<label, 32> vrt;
    DynList<edge, 64> edges;
    DynList<DynList<label, 8> > faceEdges;
    faceEdges.setSize(c.size());
    DynList<DynList<label, 2>, 64> edgeFaces;
    findAddressingForCell(cellI, vrt, edges, faceEdges, edgeFaces);

    // find a vertex which will be the top of the pyramids
    //- if there exist a corner vertex which is in 3 or more patches then
    //- it is selected as the top vertex
    label topVertex = findTopVertex(cellI, vrt, edges, edgeFaces);

    //- start generating pyramids
    forAll(c, fI)
    {
        # ifdef DEBUGDecompose
        Info << "Face " << faces[c[fI]] << " is a base face" << endl;
        #endif
        const face& f = faces[c[fI]];
        DynList<DynList<label, 8> > cellFaces;
        cellFaces.setSize(f.size() + 1);

        DynList<triFace> triFaces;
        triFaces.setSize(f.size());
        forAll(triFaces, pI)
        {
            triFaces[pI][0] = f.nextLabel(pI);
            triFaces[pI][1] = f[pI];
            triFaces[pI][2] = topVertex;
        }

        label cfI(0);
        if( owner[c[fI]] == cellI )
        {
            cellFaces[cfI++] = faces[c[fI]];

            forAll(triFaces, tfI)
            {
                cellFaces[cfI++] = triFaces[tfI];
            }
        }
        else
        {
            cellFaces[cfI++] = faces[c[fI]].reverseFace();

            forAll(triFaces, tfI)
            {
                triFace rf;
                rf[0] = triFaces[tfI][0];
                rf[1] = triFaces[tfI][2];
                rf[2] = triFaces[tfI][1];
                cellFaces[cfI++] = rf;
            }
        }

        # ifdef DEBUGDecompose
        Info << "Cell for face is " << cellFaces << endl;

        DynList<edge, 64> cEdges;
        DynList<DynList<label, 2>, 64> eFaces;
        forAll(cellFaces, fI)
        {
            const DynList<label, 8>& f = cellFaces[fI];
            forAll(f, eI)
            {
                const edge e(f[eI], f[(eI+1)%f.size()]);

                const label pos = cEdges.contains(e);

                if( pos < 0 )
                {
                    cEdges.append(e);
                    DynList<label, 2> ef;
                    ef.append(fI);
                    eFaces.append(ef);
                }
                else
                {
                    eFaces[pos].append(fI);
                }
            }
        }

        forAll(eFaces, eI)
            if( eFaces[eI].size() != 2 )
                Pout << "This pyrmid is not topologically closed" << endl;
        # endif

        facesOfNewCells_.appendGraph(cellFaces);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //

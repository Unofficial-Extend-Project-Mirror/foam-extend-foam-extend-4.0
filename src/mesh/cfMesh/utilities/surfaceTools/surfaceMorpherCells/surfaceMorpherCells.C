/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "surfaceMorpherCells.H"
#include "demandDrivenData.H"

//#define DEBUGMorph

# ifdef DEBUGMorph
#include <sstream>
#include "polyMeshGenAddressing.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void surfaceMorpherCells::replaceMeshBoundary()
{
    wordList patchNames(1);
    patchNames[0] = "defaultFaces";

    polyMeshGenModifier(mesh_).replaceBoundary
    (
        patchNames,
        newBoundaryFaces_,
        newBoundaryOwners_,
        newBoundaryPatches_
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

surfaceMorpherCells::surfaceMorpherCells
(
    polyMeshGen& mesh
)
    :
    mesh_(mesh),
    nIntFaces_(0),
    boundaryVertex_(mesh.points().size()),
    cellFlags_(mesh.cells().size()),
    newBoundaryFaces_(),
    newBoundaryOwners_(),
    newBoundaryPatches_()
{
}

surfaceMorpherCells::~surfaceMorpherCells()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void surfaceMorpherCells::morphMesh()
{
    //- perform surface morphing
    bool changed;

    # ifdef DEBUGMorph
    label iter(0);
    # endif

    do
    {
        changed = false;
        # ifdef DEBUGMorph
        Info << "Iteration " << ++iter << endl;
        # endif

        findBoundaryVertices();

        findBoundaryCells();

        if( removeCellsWithAllVerticesAtTheBoundary() )
        {
            changed = true;
            continue;
        }

        if( morphInternalFaces() )
        {
            changed = true;
            continue;
        }

        if( morphBoundaryFaces() )
        {
            changed = true;
            continue;
        }

        # ifdef DEBUGMorph
        mesh_.write();
        mesh_.addressingData().checkMesh(true);
        fileName name("morphedMesh");
        std::ostringstream ss;
        ss << iter;
        name += ss.str();
        Info << "name " << name << endl;
        ++iter;

        const cellListPMG& cells = mesh_.cells();
        const faceListPMG& faces = mesh_.faces();
        forAll(cells, cellI)
        {
            const cell& c = cells[cellI];

            const edgeList edges = c.edges(faces);
            List<direction> nAppearances(edges.size(), direction(0));

            forAll(c, fI)
            {
                const edgeList fEdges = faces[c[fI]].edges();

                forAll(fEdges, eI)
                    forAll(edges, eJ)
                        if( fEdges[eI] == edges[eJ] )
                        {
                            ++nAppearances[eJ];
                            break;
                        }
            }

            bool closed(true);
            forAll(nAppearances, eI)
                if( nAppearances[eI] != 2 )
                {
                    closed = false;
                    Info << "Edge " << edges[eI] << " appears "
                        << label(nAppearances[eI]) << " times in cell "
                        << cellI << endl;
                }

            if( !closed )
            {
                Info << "Cell " << cellI << " consists of faces " << c << endl;
                forAll(c, fI)
                    Info << "Face " << c[fI] << " is " << faces[c[fI]] << endl;
                FatalErrorIn
                (
                    "void surfaceMorpherCells::morphMesh()"
                ) << "Cell " << cellI << " is not topologically closed"
                    << abort(FatalError);
            }
        }
        # endif

    } while( changed );

    polyMeshGenModifier(mesh_).removeUnusedVertices();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

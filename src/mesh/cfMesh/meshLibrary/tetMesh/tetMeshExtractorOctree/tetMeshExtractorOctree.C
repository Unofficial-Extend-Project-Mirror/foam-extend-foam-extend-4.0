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

#include "tetMeshExtractorOctree.H"
#include "meshOctree.H"
#include "triSurface.H"
#include "polyMeshGenModifierAddCellByCell.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// #define DEBUGTets

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetMeshExtractorOctree::createPoints()
{
    polyMeshGenModifier meshModifier ( mesh_ );
    pointFieldPMG& points = meshModifier.pointsAccess();

    const LongList<point>& tetPoints = tetCreator_.tetPoints();

    points.setSize(tetPoints.size());

    # ifdef USE_OMP
    # pragma omp parallel for
    # endif
    forAll(tetPoints, pointI)
        points[pointI] = tetPoints[pointI];
}

void tetMeshExtractorOctree::createPolyMesh()
{
    polyMeshGenModifier meshModifier ( mesh_ );

    faceListPMG& faces = meshModifier.facesAccess();
    cellListPMG& cells = meshModifier.cellsAccess();
    meshModifier.boundariesAccess().setSize ( 0 );
    meshModifier.procBoundariesAccess().setSize ( 0 );

    const LongList<partTet>& tets = tetCreator_.tets();

    VRWGraph pTets;
    pTets.reverseAddressing(mesh_.points().size(), tets);

    //- set the number of cells
    cells.setSize(tets.size());

    //- all faces of tetrahedral cells
    faces.setSize(4*tets.size());
    boolList removeFace(faces.size());

    # ifdef USE_OMP
    # pragma omp parallel if( tets.size() > 1000 )
    # endif
    {
        //- set face labels
        # ifdef USE_OMP
        # pragma omp for
        # endif
        forAll(removeFace, faceI)
            removeFace[faceI] = false;

        //- set sizes of cells and create all faces
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 20)
        # endif
        forAll(tets, elmtI)
        {
            cells[elmtI].setSize(4);

            const partTet& elmt = tets[elmtI];

            //tessellationElement telmt(elmt[0], elmt[1], elmt[2], elmt[3]);

            label faceI = 4 * elmtI;

            //- first face
            cells[elmtI][0] = faceI;

            face& f0 = faces[faceI];
            f0.setSize(3);

            f0[0] = elmt.a();
            f0[1] = elmt.c();
            f0[2] = elmt.b();

            ++faceI;

            //- second face
            cells[elmtI][1] = faceI;

            face& f1 = faces[faceI];
            f1.setSize(3);

            f1[0] = elmt.a();
            f1[1] = elmt.b();
            f1[2] = elmt.d();

            ++faceI;

            //- third face
            cells[elmtI][2] = faceI;

            face& f2 = faces[faceI];
            f2.setSize ( 3 );

            f2[0] = elmt.b();
            f2[1] = elmt.c();
            f2[2] = elmt.d();

            ++faceI;

            //- fourth face
            cells[elmtI][3] = faceI;

            face& f3 = faces[faceI];
            f3.setSize ( 3 );

            f3[0] = elmt.c();
            f3[1] = elmt.a();
            f3[2] = elmt.d();
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- find duplicate faces
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 20)
        # endif
        forAll(cells, cellI)
        {
            cell& c = cells[cellI];

            forAll(c, fI)
            {
                const face& f = faces[c[fI]];
                const label pointI = f[0];

                forAllRow(pTets, pointI, ptI)
                {
                    //- do not check cells with greater labels
                    //- they cannot be face owners
                    if( pTets(pointI, ptI) >= cellI )
                        continue;

                    const cell& otherTet = cells[pTets(pointI, ptI)];

                    //- check faces created from a tet
                    forAll(otherTet, ofI)
                    {
                        //- do not compare faces with greater labels
                        //- they shall not be removed here
                        if( otherTet[ofI] >= c[fI] )
                            continue;

                        //- check if the faces are equal
                        if( f == faces[otherTet[ofI]] )
                        {
                            removeFace[c[fI]] = true;
                            c[fI] = otherTet[ofI];
                        }
                    }
                }
            }
        }
    }

    //- remove duplicate faces
    label nFaces(0);
    labelLongList newFaceLabel(faces.size(), -1);

    forAll(faces, faceI)
    {
        if( !removeFace[faceI] )
        {
            if( nFaces < faceI )
                faces[nFaces].transfer(faces[faceI]);

            newFaceLabel[faceI] = nFaces;
            ++nFaces;
        }
    }

    //- set the size of faces
    faces.setSize(nFaces);

    //- change cells
    # ifdef USE_OMP
    # pragma omp for schedule(dynamic, 40)
    # endif
    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        DynList<label> newC;

        forAll(c, fI)
        {
            if( newFaceLabel[c[fI]] != -1 )
                newC.append(newFaceLabel[c[fI]]);
        }

        c.setSize(newC.size());
        forAll(c, fI)
        c[fI] = newC[fI];
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from octree and mesh data
tetMeshExtractorOctree::tetMeshExtractorOctree
(
    const meshOctree& octree,
    const IOdictionary& meshDict,
    polyMeshGen& mesh
)
:
    tetCreator_(octree, meshDict),
    mesh_(mesh)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tetMeshExtractorOctree::~tetMeshExtractorOctree()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetMeshExtractorOctree::createMesh()
{
    Info << "Extracting tetMesh" << endl;

    //- copy tet points into the mesh
    createPoints();

    //- create the mesh
    createPolyMesh();

    polyMeshGenModifier(mesh_).reorderBoundaryFaces();
    polyMeshGenModifier(mesh_).removeUnusedVertices();

    Info << "Mesh has :" << nl
    << mesh_.points().size() << " vertices " << nl
    << mesh_.faces().size() << " faces" << nl
    << mesh_.cells().size() << " cells" << endl;

    Info << "Finished extracting tetMesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

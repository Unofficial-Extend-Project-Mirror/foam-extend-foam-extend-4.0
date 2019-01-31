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


#include "decomposeFaces.H"
#include "boolList.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGDec

# ifdef DEBUGDec
#include "polyMeshGenAddressing.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//- Constructor
decomposeFaces::decomposeFaces(polyMeshGen& mesh)
:
    mesh_(mesh),
    newFacesForFace_(),
    done_(false)
{}

//- Destructor
decomposeFaces::~decomposeFaces()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void decomposeFaces::decomposeMeshFaces(const boolList& decomposeFace)
{
    done_ = false;

    newFacesForFace_.setSize(mesh_.faces().size());
    forAll(newFacesForFace_, fI)
        newFacesForFace_.setRowSize(fI, 0);

    const label nIntFaces = mesh_.nInternalFaces();
    label nFaces(0), nPoints = mesh_.points().size();
    point p;

    pointFieldPMG& points = mesh_.points();
    forAll(decomposeFace, fI)
        if( decomposeFace[fI] )
            ++nFaces;

    points.setSize(nPoints + nFaces);

    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();

    if( decomposeFace.size() != faces.size() )
        FatalErrorIn
        (
            "void decomposeFaces::decomposeMeshFaces(const boolList&)"
        ) << "Incorrect size of the decomposeFace list!" << abort(FatalError);

    nFaces = 0;
    VRWGraph newFaces;

    //- decompose internal faces
    for(label faceI=0;faceI<nIntFaces;++faceI)
    {
        const face& f = faces[faceI];

        if( decomposeFace[faceI] )
        {
            # ifdef DEBUGDec
            Info << "Decomposing internal face " << faceI << " with nodes "
                << f << endl;
            # endif

            FixedList<label, 3> newF;

            forAll(f, pI)
            {
                newF[0] = f[pI];
                newF[1] = f.nextLabel(pI);
                newF[2] = nPoints;

                # ifdef DEBUGDec
                Info << "Storing face " << newF << " with label "
                    << nFaces << endl;
                # endif

                newFaces.appendList(newF);
                newFacesForFace_.append(faceI, nFaces++);
            }

            p = f.centre(points);
            points[nPoints] = p;
            ++nPoints;
        }
        else
        {
            # ifdef DEBUGDec
            Info << "Storing internal face " << faceI << " with nodes "
                << f << " as new face " << faceI << endl;
            # endif

            newFaces.appendList(f);
            newFacesForFace_.append(faceI, nFaces++);
        }
    }

    //- decompose boundary faces
    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        boundaries[patchI].patchStart() = nFaces;

        for(label bfI=start;bfI<end;++bfI)
        {
            const face& f = faces[bfI];
            if( decomposeFace[bfI] )
            {
                # ifdef DEBUGDec
                Info << "Decomposing boundary face " << bfI
                    << " with nodes " << f << endl;
                # endif

                FixedList<label, 3> newF;

                forAll(f, pI)
                {
                    newF[0] = f[pI];
                    newF[1] = f.nextLabel(pI);
                    newF[2] = nPoints;

                    # ifdef DEBUGDec
                    Info << "Storing face " << newF << " with label "
                        << nFaces << endl;
                    # endif

                    newFaces.appendList(newF);
                    newFacesForFace_.append(bfI, nFaces++);
                }

                p = f.centre(points);
                points[nPoints++] = p;
            }
            else
            {
                # ifdef DEBUGDec
                Info << "Storing boundary face " << bfI << " in patch"
                    << patchI << " as new face " << bfI << endl;
                # endif

                newFaces.appendList(f);
                newFacesForFace_.append(bfI, nFaces++);
            }
        }

        boundaries[patchI].patchSize() =
            nFaces - boundaries[patchI].patchStart();
    }

    //- decompose processor faces
    if( Pstream::parRun() )
    {
        PtrList<processorBoundaryPatch>& procBoundaries =
            meshModifier.procBoundariesAccess();

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label end = start + procBoundaries[patchI].patchSize();

            const bool own = procBoundaries[patchI].owner();

            procBoundaries[patchI].patchStart() = nFaces;

            for(label bfI=start;bfI<end;++bfI)
            {
                face f;
                if( own )
                {
                    f = faces[bfI];
                }
                else
                {
                    f = faces[bfI].reverseFace();
                }

                if( decomposeFace[bfI] )
                {
                    # ifdef DEBUGDec
                    Info << "Decomposing processor boundary face " << bfI
                        << " with nodes " << f << endl;
                    # endif

                    face newF(3);

                    forAll(f, pI)
                    {
                        newF[0] = f[pI];
                        newF[1] = f.nextLabel(pI);
                        newF[2] = nPoints;

                        # ifdef DEBUGDec
                        Info << "Storing face " << newF << " with label "
                            << nFaces << endl;
                        # endif

                        if( own )
                        {
                            newFaces.appendList(newF);
                        }
                        else
                        {
                            newFaces.appendList(newF.reverseFace());
                        }
                        newFacesForFace_.append(bfI, nFaces++);
                    }

                    p = f.centre(points);
                    points[nPoints++] = p;
                }
                else
                {
                    # ifdef DEBUGDec
                    Info << "Storing boundary face " << bfI << " in patch"
                        << patchI << " as new face " << bfI << endl;
                    # endif

                    if( own )
                    {
                        newFaces.appendList(f);
                    }
                    else
                    {
                        newFaces.appendList(f.reverseFace());
                    }
                    newFacesForFace_.append(bfI, nFaces++);
                }
            }

            procBoundaries[patchI].patchSize() =
                nFaces - procBoundaries[patchI].patchStart();
        }
    }

    //- store the faces back into their list
    faces.setSize(nFaces);
    forAll(faces, faceI)
    {
        face& f = faces[faceI];
        f.setSize(newFaces.sizeOfRow(faceI));

        forAll(f, pI)
            f[pI] = newFaces(faceI, pI);
    }
    newFaces.setSize(0);

    //- update subsets
    mesh_.updateFaceSubsets(newFacesForFace_);

    //- change the mesh
    cellListPMG& cells = meshModifier.cellsAccess();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        DynList<label> newC;

        forAll(c, fJ)
        {
            const label faceI = c[fJ];
            forAllRow(newFacesForFace_, faceI, nfI)
            newC.append(newFacesForFace_(faceI, nfI));
        }

        # ifdef DEBUGDec
        Info << "Cell " << cellI << " with faces " << c
            << " is changed into " << newC << endl;
        # endif

        c.setSize(newC.size());
        forAll(newC, fJ)
        c[fJ] = newC[fJ];
    }

    meshModifier.clearAll();

    done_ = true;

    # ifdef DEBUGDec
    Info << "New cells are " << cells << endl;
    mesh_.write();
    # endif
}

void decomposeFaces::decomposeConcaveInternalFaces
(
    const boolList& concaveVertex
)
{
    if( Pstream::parRun() )
    {
        FatalErrorIn
        (
            "void decomposeFaces::decomposeConcaveInternalFaces"
            "(const boolList& concaveVertex)"
        ) << "This procedure is not parallelised!" << exit(FatalError);
    }

    done_ = false;

    newFacesForFace_.setSize(mesh_.faces().size());
    forAll(newFacesForFace_, fI)
        newFacesForFace_.setRowSize(fI, 0);

    const label nIntFaces = mesh_.nInternalFaces();

    polyMeshGenModifier meshModifier(mesh_);
    pointFieldPMG& points = meshModifier.pointsAccess();
    faceListPMG& faces = meshModifier.facesAccess();

    if( concaveVertex.size() != mesh_.points().size() )
        FatalErrorIn
        (
            "void decomposeFaces::decomposeMeshFaces(const boolList&)"
        ) << "Incorrect size of the concaveVertex list!" << abort(FatalError);

    VRWGraph newFaces;
    DynList<label> newF;
    newF.setSize(3);

    # ifdef DEBUGDec
    const label id = mesh_.addFaceSubset("decomposedFaces");
    # endif

    //- decompose internal faces
    for(label faceI=0;faceI<nIntFaces;++faceI)
    {
        const face& f = faces[faceI];

        DynList<label> concavePos;
        forAll(f, pI)
            if( concaveVertex[f[pI]] )
            {
                concavePos.append(pI);
            }

        if( concavePos.size() == 1 )
        {
            # ifdef DEBUGDec
            Info << "1. Decomposing internal face " << faceI << " with nodes "
                << f << endl;
            mesh_.addFaceToSubset(id, faceI);
            # endif

            newF[0] = f[concavePos[0]];

            for(label pI=1;pI<(f.size()-1);++pI)
            {
                const label pJ = (concavePos[0] + pI) % f.size();
                newF[1] = f[pJ];
                newF[2] = f.nextLabel(pJ);

                # ifdef DEBUGDec
                Info << "Storing face " << newF << " with label "
                    << newFaces.size() << endl;
                # endif

                newFacesForFace_.append(faceI, newFaces.size());
                newFaces.appendList(newF);
            }
        }
        else if( concavePos.size() > 1 )
        {
            # ifdef DEBUGDec
            Info << "2. Decomposing internal face " << faceI << " with nodes "
                << f << endl;
            mesh_.addFaceToSubset(id, faceI);
            # endif

            newF[0] = points.size();
            forAll(f, pI)
            {
                newF[1] = f[pI];
                newF[2] = f.nextLabel(pI);

                # ifdef DEBUGDec
                Info << "2. Storing face " << newF << " with label "
                    << newFaces.size() << endl;
                # endif

                newFacesForFace_.append(faceI, newFaces.size());
                newFaces.appendList(newF);
            }

            const point fCent = f.centre(points);
            points.append(fCent);
        }
        else
        {
            # ifdef DEBUGDec
            Info << "Storing internal face " << faceI << " with nodes "
                << f << " as new face " << newFaces.size() << endl;
            # endif

            newFacesForFace_.append(faceI, newFaces.size());
            newFaces.appendList(f);
        }
    }

    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        //- set new patch start
        boundaries[patchI].patchStart() = newFaces.size();

        //- store faces into newFaces
        for(label bfI=start;bfI<end;++bfI)
        {
            newFacesForFace_.append(bfI, newFaces.size());
            newFaces.appendList(faces[bfI]);
        }
    }

    //- copy new faces into the faceListPMG
    faces.setSize(newFaces.size());
    forAll(newFaces, faceI)
    {
        faces[faceI].setSize(newFaces.sizeOfRow(faceI));

        forAllRow(newFaces, faceI, pI)
            faces[faceI][pI] = newFaces(faceI, pI);
    }

    newFaces.setSize(0);

    //- update cells
    cellListPMG& cells = meshModifier.cellsAccess();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        DynList<label, 24> newC;

        forAll(c, fJ)
        {
            const label faceI = c[fJ];
            forAllRow(newFacesForFace_, faceI, nfI)
            newC.append(newFacesForFace_(faceI, nfI));
        }

        # ifdef DEBUGDec
        Info << "Cell " << cellI << " with faces " << c
            << " is changed into " << newC << endl;
        # endif

        c.setSize(newC.size());
        forAll(newC, fJ)
            c[fJ] = newC[fJ];
    }

    meshModifier.clearAll();

    //- update subsets
    mesh_.updateFaceSubsets(newFacesForFace_);

    # ifdef DEBUGDec
    Info << "New cells are " << cells << endl;
    mesh_.write();
    ::exit(1);
    # endif

    meshModifier.removeUnusedVertices();

    done_ = true;
}

const VRWGraph& decomposeFaces::newFacesForFace() const
{
    if( !done_ )
        WarningIn
        (
            "const VRWGraph& decomposeFaces::newFacesForFace() const"
        ) << "Decomposition is not yet performed!" << endl;

    return newFacesForFace_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //

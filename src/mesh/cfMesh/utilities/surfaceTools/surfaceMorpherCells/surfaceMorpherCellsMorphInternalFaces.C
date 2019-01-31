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

#include "surfaceMorpherCells.H"
#include "demandDrivenData.H"
#include "helperFunctions.H"
#include "HashSet.H"

//#define DEBUGMorph

#ifdef DEBUGMorph
#include "polyMeshGenAddressing.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void surfaceMorpherCells::findBoundaryVertices()
{
    const faceListPMG& faces = mesh_.faces();

    boundaryVertex_.setSize(mesh_.points().size());
    boundaryVertex_ = false;

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        for(label faceI=start;faceI<end;++faceI)
        {
            const face& f = faces[faceI];

            forAll(f, pI)
            {
                # ifdef DEBUGMorph
                if( (f[pI] >= boundaryVertex_.size()) || (f[pI] < 0) )
                {
                    Info << f << endl;
                    FatalError << "Wrong label " << f[pI] << " in face "
                        << faceI << abort(FatalError);
                }
                # endif

                boundaryVertex_[f[pI]] = true;
            }
        }
    }

    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();

        bool changed;
        do
        {
            changed = false;

            //- send data about boundary vertices to other processors
            forAll(procBoundaries, patchI)
            {
                const label start = procBoundaries[patchI].patchStart();
                const label end = start + procBoundaries[patchI].patchSize();

                //- create information about bnd nodes which must be exchanged
                //- with other processors
                labelHashSet addToSend;
                labelLongList dts;
                for(label faceI=start;faceI<end;++faceI)
                {
                    const face& f = faces[faceI];

                    forAll(f, pI)
                        if( boundaryVertex_[f[pI]] && !addToSend.found(f[pI]) )
                        {
                            addToSend.insert(f[pI]);
                            dts.append(faceI-start);
                            dts.append((f.size()-pI)%f.size());
                        }
                }

                labelList bndVertsToSend(dts.size());
                forAll(dts, i)
                    bndVertsToSend[i] = dts[i];

                //- send the list of other processor
                OPstream toOtherProc
                (
                    Pstream::commsTypes::blocking,
                    procBoundaries[patchI].neiProcNo(),
                    bndVertsToSend.byteSize()
                );
                toOtherProc << bndVertsToSend;
            }

            //- receive boundary vertices from other processor and
            //- continue sending and receiving as long as it makes some change
            forAll(procBoundaries, patchI)
            {
                labelList receivedBndNodes;
                IPstream fromOtherProc
                (
                    Pstream::commsTypes::blocking,
                    procBoundaries[patchI].neiProcNo()
                );
                fromOtherProc >> receivedBndNodes;
                const label start = procBoundaries[patchI].patchStart();

                label entryI(0);
                while( entryI < receivedBndNodes.size() )
                {
                    const label fI = receivedBndNodes[entryI++];
                    const label pI = receivedBndNodes[entryI++];

                    const face& f = faces[start+fI];

                    if( !boundaryVertex_[f[pI]] )
                    {
                        boundaryVertex_[f[pI]] = true;
                        changed = true;
                    }
                }
            }

            reduce(changed, maxOp<bool>());
        } while( changed );
    }
}

void surfaceMorpherCells::findBoundaryCells()
{
    const labelList& owner = mesh_.owner();

    cellFlags_.setSize(mesh_.cells().size());
    cellFlags_ = NONE;

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        for(label faceI=start;faceI<end;++faceI)
            cellFlags_[owner[faceI]] = BOUNDARY;
    }
}

bool surfaceMorpherCells::morphInternalFaces()
{
    Info << "Morphing internal faces" << endl;

    # ifdef DEBUGMorph
    Serr << "Zip check 1" << endl;
    labelHashSet zipCellsBefore;
    mesh_.addressingData().checkCellsZipUp(true, &zipCellsBefore);
    if( zipCellsBefore.size() )
    {
        Serr << Pstream::myProcNo() << "Open cells found!" << endl;
        Serr << "Cells " << zipCellsBefore << " are not zipped!!" << endl;
    }
    mesh_.clearAddressingData();
    # endif

    bool changed(false);

    newBoundaryFaces_.setSize(0);
    newBoundaryOwners_.setSize(0);
    newBoundaryPatches_.setSize(0);

    const label nIntFaces = mesh_.nInternalFaces();
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    //- copy boundary faces
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        for(label faceI=start;faceI<end;++faceI)
        {
            newBoundaryFaces_.appendList(faces[faceI]);
            newBoundaryOwners_.append(owner[faceI]);
            newBoundaryPatches_.append(0);
        }
    }

    //- start morphing internal faces
    for(label faceI=0;faceI<nIntFaces;++faceI)
    {
        # ifdef DEBUGMorph
        Info << "Morphing internal face " << faceI << endl;
        # endif

        const face& f = faces[faceI];

        DynList<bool> removeFaceVertex(f.size(), false);

        face newF(f.size());
        label i(0);

        forAll(f, pI)
            if(
                boundaryVertex_[f.prevLabel(pI)] &&
                boundaryVertex_[f[pI]] &&
                boundaryVertex_[f.nextLabel(pI)]
            )
            {
                removeFaceVertex[pI] = true;

                # ifdef DEBUGMorph
                Info << "Removing vertex " << f[pI] << " from face "
                    << f << endl;
                # endif
            }
            else
            {
                newF[i++] = f[pI];
            }

        if( i < f.size() )
        {
            changed = true;

            //- store shrinked face
            newF.setSize(i);

            # ifdef DEBUGMorph
            Info << "Removing edges " << removeEdge << endl;
            Info << "Face label " << faceI << " owner " << owner[faceI]
                << " neighbour " << neighbour[faceI] << endl;
            Info << "Internal face " << f << " has been shrinked into "
                << newF << endl;
            # endif

            //- create new boundary faces from the removed part
            label mat(1);
            DynList<direction> nodeMaterial(f.size(), direction(0));
            DynList<DynList<edge>, 2> edgeMats;
            forAll(nodeMaterial, nI)
                if( !nodeMaterial[nI] && removeFaceVertex[nI] )
                {
                    edgeMats.append(DynList<edge>());
                    DynList<label> front;
                    front.append(nI);

                    do
                    {
                        DynList<label> newFront;
                        forAll(front, pI)
                        {
                            const label fLabel = front[pI];
                            if( nodeMaterial[fLabel] )
                                continue;
                            nodeMaterial[fLabel] = mat;
                            edgeMats[mat-1].appendIfNotIn(f.faceEdge(fLabel));
                            edgeMats[mat-1].appendIfNotIn
                            (
                                f.faceEdge(f.rcIndex(fLabel))
                            );

                            if( removeFaceVertex[f.rcIndex(fLabel)] )
                                newFront.append(f.rcIndex(fLabel));
                            if( removeFaceVertex[f.fcIndex(fLabel)] )
                                newFront.append(f.fcIndex(fLabel));
                        }
                        front = newFront;
                    }
                    while( front.size() != 0 );

                    ++mat;
                }

            forAll(edgeMats, matI)
            {
                edgeMats[matI].shrink();
                const face rf = help::removeEdgesFromFace(f, edgeMats[matI]);
                const face newBf = help::createFaceFromRemovedPart(f, rf);

                # ifdef DEBUGMorph
                Info << "Adding face " << newBf << " as boundary faces" << endl;
                Info << "Owner and neighbour are " << owner[faceI] << " and "
                    << neighbour[faceI] << endl;
                # endif

                newBoundaryFaces_.appendList(newBf);
                newBoundaryOwners_.append(owner[faceI]);
                newBoundaryPatches_.append(0);

                newBoundaryFaces_.appendList(newBf.reverseFace());
                newBoundaryOwners_.append(neighbour[faceI]);
                newBoundaryPatches_.append(0);
            }

            const_cast<face&>(faces[faceI]) = newF;
        }
    }

    //- treat processor boundaries
    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label end = start + procBoundaries[patchI].patchSize();
            const bool isOwner = procBoundaries[patchI].owner();

            for(label faceI=start;faceI<end;++faceI)
            {
                face copy;
                if( isOwner )
                {
                    copy = faces[faceI];
                }
                else
                {
                    copy = faces[faceI].reverseFace();
                }

                const face& f = copy;

                DynList<bool> removeFaceVertex(f.size(), false);

                face newF(f.size());
                label i(0);

                forAll(f, pI)
                    if(
                        boundaryVertex_[f.prevLabel(pI)] &&
                        boundaryVertex_[f[pI]] &&
                        boundaryVertex_[f.nextLabel(pI)]
                    )
                    {
                        removeFaceVertex[pI] = true;

                        # ifdef DEBUGMorph
                        Info << "Removing vertex " << f[pI] << " from face "
                            << f << endl;
                        # endif
                    }
                    else
                    {
                        newF[i++] = f[pI];
                    }

                if( i < f.size() )
                {
                    changed = true;

                    //- store shrinked face
                    newF.setSize(i);

                    # ifdef DEBUGMorph
                    Info << "Removing edges " << removeEdge << endl;
                    Info << "Face label " << faceI << " owner " << owner[faceI]
                        << " neighbour " << neighbour[faceI] << endl;
                    Info << "Internal face " << f << " has been shrinked into "
                        << newF << endl;
                    # endif

                    //- create new boundary faces from the removed part
                    label mat(1);
                    DynList<direction> nodeMaterial(f.size(), direction(0));
                    DynList< DynList<edge> > edgeMats;
                    forAll(nodeMaterial, nI)
                        if( !nodeMaterial[nI] && removeFaceVertex[nI] )
                        {
                            edgeMats.append(DynList<edge>());
                            DynList<label> front;
                            front.append(nI);

                            do
                            {
                                DynList<label> newFront;
                                forAll(front, pI)
                                {
                                    const label fLabel = front[pI];
                                    if( nodeMaterial[fLabel] ) continue;
                                    nodeMaterial[fLabel] = mat;
                                    edgeMats[mat-1].appendIfNotIn
                                    (
                                        f.faceEdge(fLabel)
                                    );
                                    edgeMats[mat-1].appendIfNotIn
                                    (
                                        f.faceEdge(f.rcIndex(fLabel))
                                    );

                                    if( removeFaceVertex[f.rcIndex(fLabel)] )
                                        newFront.append(f.rcIndex(fLabel));
                                    if( removeFaceVertex[f.fcIndex(fLabel)] )
                                        newFront.append(f.fcIndex(fLabel));
                                }
                                front = newFront;
                            }
                            while( front.size() != 0 );

                            ++mat;
                        }

                    forAll(edgeMats, matI)
                    {
                        edgeMats[matI].shrink();
                        const face rf =
                            help::removeEdgesFromFace(f, edgeMats[matI]);
                        const face newBf =
                            help::createFaceFromRemovedPart(f, rf);

                        # ifdef DEBUGMorph
                        Info << "Adding face " << newBf
                        << " as boundary faces" << endl;
                        Info << "Owner and neighbour are "
                            << owner[faceI] << endl;
                        # endif

                        if( isOwner )
                        {
                            newBoundaryFaces_.appendList(newBf);
                        }
                        else
                        {
                            newBoundaryFaces_.appendList(newBf.reverseFace());
                        }
                        newBoundaryOwners_.append(owner[faceI]);
                        newBoundaryPatches_.append(0);
                    }

                    if( isOwner || (newF.size() < 3) )
                    {
                        const_cast<face&>(faces[faceI]) = newF;
                    }
                    else
                    {
                        const_cast<face&>(faces[faceI]) = newF.reverseFace();
                    }
                }
            }
        }
    }

    polyMeshGenModifier meshModifier(mesh_);

    if( Pstream::parRun() )
    {
        reduce(changed, maxOp<bool>());
    }

    if( changed )
    {
        //- replace boundary of the mesh
        replaceMeshBoundary();

        # ifdef DEBUGMorph
        Serr << "Zip check 2" << endl;
        const cellListPMG& cells = mesh_.cells();
        forAll(cells, cI)
        {
            const cell& c = cells[cI];

            DynList<edge> edges;
            DynList<direction> nAppearances;

            forAll(c, fI)
            {
                const face& f = faces[c[fI]];
                forAll(f, eI)
                {
                    const label pos = edges.containsAtPosition(f.faceEdge(eI));

                    if( pos == -1 )
                    {
                        edges.append(f.faceEdge(eI));
                        nAppearances.append(1);
                    }
                    else
                    {
                        ++nAppearances[pos];
                    }
                }
            }

            forAll(nAppearances, eI)
                if( nAppearances[eI] != 2 )
                {
                    Warning << "Cell " << cI << " is not closed" << endl;
                    Serr << "Cell faces " << c << endl;
                    forAll(c, fI)
                        Serr << "Face " << c[fI] << " is "
                            << faces[c[fI]] << endl;
                }
        }
        # endif

        //- remove faces which do not exist any more
        boolList removeFace(faces.size(), false);
        bool removeFaces(false);

        for(label faceI=0;faceI<nIntFaces;++faceI)
            if( faces[faceI].size() < 3 )
            {
                removeFace[faceI] = true;
                removeFaces = true;
            }

        if( Pstream::parRun() )
        {
            const PtrList<processorBoundaryPatch>& procBoundaries =
                mesh_.procBoundaries();

            forAll(procBoundaries, patchI)
            {
                const label start = procBoundaries[patchI].patchStart();
                const label end = start + procBoundaries[patchI].patchSize();
                for(label faceI=start;faceI<end;++faceI)
                    if( faces[faceI].size() < 3 )
                    {
                        removeFace[faceI] = true;
                        removeFaces = true;
                    }
            }

            reduce(removeFaces, maxOp<bool>());
        }

        if( removeFaces )
            meshModifier.removeFaces(removeFace);
    }

    # ifdef DEBUGMorph
    Serr << "Zip check 3" << endl;
    labelHashSet zipCells;
    mesh_.addressingData().checkCellsZipUp(true, &zipCells);
    if( zipCells.size() )
    {
        Serr << Pstream::myProcNo() << "Open cells appeared!" << endl;
        //mesh_.write();
        Serr << "Cells " << zipCells << " are not zipped!!" << endl;
    }
    mesh_.clearAddressingData();
    # endif

    Info << "Finished morphing internal faces" << endl;

    return changed;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

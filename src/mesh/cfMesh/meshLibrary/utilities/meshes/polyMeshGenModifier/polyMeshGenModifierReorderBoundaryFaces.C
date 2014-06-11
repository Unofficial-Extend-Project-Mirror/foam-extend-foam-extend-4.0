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

#include "polyMeshGenModifier.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::reorderBoundaryFaces()
{
    Info << "Reordering boundary faces " << endl;

    if( Pstream::parRun() )
        reorderProcBoundaryFaces();

    faceListPMG& faces = mesh_.faces_;
    cellListPMG& cells = mesh_.cells_;

    const labelList& neighbour = mesh_.neighbour();
    const label nInternalFaces = mesh_.nInternalFaces();

    //- count internal and boundary faces
    const label numBFaces = faces.size() - nInternalFaces;

    labelLongList newFaceLabel(faces.size(), -1);

    //- find faces which should be repositioned
    label nReplaced(0);
    labelList internalToChange;
    labelList boundaryToChange;

    # ifdef USE_OMP
    const label nThreads = 3 * omp_get_num_procs();
    # else
    const label nThreads(1);
    # endif
    labelList nInternalToChangeThread(nThreads);
    labelList nBoundaryToChangeThread(nThreads);

    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI = 0;
        # endif

        label& nItc = nInternalToChangeThread[threadI];
        label& nBtc = nBoundaryToChangeThread[threadI];

        labelLongList internalToChangeLocal, boundaryToChangeLocal;

        //- find the boundary faces within the range of internal faces
        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        for(label faceI=0;faceI<nInternalFaces;++faceI)
        {
            if( neighbour[faceI] == -1 )
                internalToChangeLocal.append(faceI);
        }

        nItc = internalToChangeLocal.size();

        //- find the internal faces within the range of boundary faces
        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        for(label faceI=nInternalFaces;faceI<faces.size();++faceI)
        {
            if( neighbour[faceI] != -1 )
                boundaryToChangeLocal.append(faceI);
        }

        nBtc = boundaryToChangeLocal.size();

        //- perform reduction such that all threads know how many faces
        //- need to be swapped
        # ifdef USE_OMP
        # pragma omp critical
        # endif
        nReplaced += nBtc;

        # ifdef USE_OMP
        # pragma omp barrier

        # pragma omp master
        # endif
        {
            internalToChange.setSize(nReplaced);
            boundaryToChange.setSize(nReplaced);
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        label localStart(0);
        for(label i=0;i<threadI;++i)
            localStart += nInternalToChangeThread[i];

        forAll(internalToChangeLocal, i)
            internalToChange[localStart++] = internalToChangeLocal[i];

        localStart = 0;
        for(label i=0;i<threadI;++i)
            localStart += nBoundaryToChangeThread[i];

        forAll(boundaryToChangeLocal, i)
            boundaryToChange[localStart++] = boundaryToChangeLocal[i];

        # ifdef USE_OMP
        # pragma omp barrier

        //- start moving positions of faces
        # pragma omp for schedule(static)
        # endif
        forAll(internalToChange, fI)
        {
            //- swap with the face at the location the face should be
            face f;
            f.transfer(faces[internalToChange[fI]]);
            faces[internalToChange[fI]].transfer(faces[boundaryToChange[fI]]);
            faces[boundaryToChange[fI]].transfer(f);
            newFaceLabel[internalToChange[fI]] = boundaryToChange[fI];
            newFaceLabel[boundaryToChange[fI]] = internalToChange[fI];
        }

        # ifdef USE_OMP
        # pragma omp barrier

        //- renumber cells
        # pragma omp for schedule(dynamic, 40)
        # endif
        forAll(cells, cellI)
        {
            cell& c = cells[cellI];

            forAll(c, fI)
                if( newFaceLabel[c[fI]] != -1 )
                    c[fI] = newFaceLabel[c[fI]];
        }
    }

    //- re-create boundary data
    PtrList<boundaryPatch>& boundaries = mesh_.boundaries_;
    if( boundaries.size() != 1 )
    {
        boundaries.clear();
        boundaries.setSize(1);
        boundaries.set
        (
            0,
            new boundaryPatch
            (
                "defaultFaces",
                "patch",
                numBFaces,
                nInternalFaces
            )
        );
    }
    else
    {
        boundaries[0].patchStart() = nInternalFaces;
        boundaries[0].patchSize() = numBFaces;
    }

    if( Pstream::parRun() )
    {
        //- processor boundary faces must be contained at the end
        label nProcFaces(0);
        forAll(mesh_.procBoundaries_, procPatchI)
            nProcFaces += mesh_.procBoundaries_[procPatchI].patchSize();

        boundaries[0].patchSize() -= nProcFaces;
    }

    //- update face subsets
    mesh_.updateFaceSubsets(newFaceLabel);

    //- delete invalid data
    mesh_.clearOut();
    this->clearOut();

    Info << "Finished reordering boundary faces" << endl;
}

void polyMeshGenModifier::reorderProcBoundaryFaces()
{
    PtrList<processorBoundaryPatch>& procBoundaries = mesh_.procBoundaries_;
    if( procBoundaries.size() == 0 )
    {
        Warning << "Processor " << Pstream::myProcNo() << " has no "
            << "processor boundaries!" << endl;
        return;
    }

    //- check if there exist any internal or ordinary bnd faces
    //- which appear after processor bnd faces. Move those faces before
    //- the processor boundary
    const label origProcStart = procBoundaries[0].patchStart();
    label nProcFaces(0);
    forAll(procBoundaries, patchI)
        nProcFaces += procBoundaries[patchI].patchSize();

    faceListPMG& faces = mesh_.faces_;
    cellListPMG& cells = mesh_.cells_;

    const label shift = faces.size() - (origProcStart + nProcFaces);
    if( shift == 0 )
        return;
    if( shift < 0 )
        FatalErrorIn
        (
            "void polyMeshGenModifier::reorderProcBoundaryFaces()"
        ) << "Missing some faces!" << abort(FatalError);

    labelLongList newFaceLabel(faces.size(), -1);

    //- faces added after processor boundaries should be moved up front
    faceList facesAtEnd(shift);
    label counter(0);
    for(label faceI=(origProcStart + nProcFaces);faceI<faces.size();++faceI)
    {
        facesAtEnd[counter].transfer(faces[faceI]);
        newFaceLabel[faceI] = origProcStart + counter;
        ++counter;
    }

    //- shift proc faces
    forAllReverse(procBoundaries, patchI)
    {
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + procBoundaries[patchI].patchSize();

        //- set patch start to the new value
        procBoundaries[patchI].patchStart() += shift;

        for(label faceI=end-1;faceI>=start;--faceI)
        {
            faces[faceI+shift].transfer(faces[faceI]);
            newFaceLabel[faceI] = faceI + shift;
        }
    }

    //- store faces taken from the end
    forAll(facesAtEnd, fI)
    {
        faces[origProcStart+fI].transfer(facesAtEnd[fI]);
    }

    //- set correct patch size
    PtrList<boundaryPatch>& boundaries = mesh_.boundaries_;
    if( boundaries.size() == 1 )
    {
        boundaries[0].patchSize() =
            procBoundaries[0].patchStart() - boundaries[0].patchStart();
    }
    else
    {
        const label start = boundaries[0].patchStart();

        boundaries.clear();
        boundaries.setSize(1);
        boundaries.set
        (
            0,
            new boundaryPatch
            (
                "defaultFaces",
                "patch",
                procBoundaries[0].patchStart() - start,
                start
            )
        );
    }

    //- renumber cells
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        forAll(c, fI)
            if( newFaceLabel[c[fI]] != -1 )
                c[fI] = newFaceLabel[c[fI]];
    }

    //- update face subsets
    mesh_.updateFaceSubsets(newFaceLabel);

    //- delete invalid data
    mesh_.clearOut();
    this->clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

#include "meshOptimizer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class labelListType>
inline void meshOptimizer::lockCells(const labelListType& l)
{
    boolList lockedFace(mesh_.faces().size(), false);
    const cellListPMG& cells = mesh_.cells();
    forAll(l, lcI)
    {
        const cell& c = cells[l[lcI]];

        forAll(c, fI)
            lockedFace[c[fI]] = true;
    }

    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();

        forAll(procBoundaries, patchI)
        {
            labelLongList dataToSend;

            const label start = procBoundaries[patchI].patchStart();
            const label end = start+procBoundaries[patchI].patchSize();

            for(label faceI=start;faceI<end;++faceI)
                if( lockedFace[faceI] )
                    dataToSend.append(faceI-start);

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dataToSend.byteSize()
            );

            toOtherProc << dataToSend;
        }

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            labelList receivedData;
            fromOtherProc >> receivedData;

            forAll(receivedData, i)
                lockedFace[start+receivedData[i]];
        }
    }

    //- Finally, mark locked points and faces
    const faceListPMG& faces = mesh_.faces();
    forAll(lockedFace, faceI)
    {
        if( lockedFace[faceI] )
        {
            lockedFaces_.append(faceI);

            const face& f = faces[faceI];

            forAll(f, pI)
                vertexLocation_[f[pI]] |= LOCKED;
        }
    }

    # ifdef DEBUGSmoothing
    const label lockedFacesI = mesh_.addFaceSubset("lockedFaces");
    forAll(lockedFaces_, lfI)
        mesh_.addFaceToSubset(lockedFacesI, lockedFaces_[lfI]);

    const label lockPointsI = mesh_.addPointSubset("lockedPoints");
    forAll(vertexLocation_, pointI)
        if( vertexLocation_[pointI] & LOCKED )
            mesh_.addPointToSubset(lockPointsI, pointI);
    # endif
}

template<class labelListType>
inline void meshOptimizer::lockFaces(const labelListType& lf)
{
    boolList lockedFace(mesh_.faces().size(), false);
    forAll(lf, lfI)
    {
        lockedFace[lf[lfI]] = true;
    }

    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();

        forAll(procBoundaries, patchI)
        {
            labelLongList dataToSend;

            const label start = procBoundaries[patchI].patchStart();
            const label end = start+procBoundaries[patchI].patchSize();

            for(label faceI=start;faceI<end;++faceI)
                if( lockedFace[faceI] )
                    dataToSend.append(faceI-start);

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dataToSend.byteSize()
            );

            toOtherProc << dataToSend;
        }

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            labelList receivedData;
            fromOtherProc >> receivedData;

            forAll(receivedData, i)
                lockedFace[start+receivedData[i]];
        }
    }

    //- Finally, mark locked points and faces
    const faceListPMG& faces = mesh_.faces();
    forAll(lockedFace, faceI)
    {
        if( lockedFace[faceI] )
        {
            lockedFaces_.append(faceI);

            const face& f = faces[faceI];

            forAll(f, pI)
                vertexLocation_[f[pI]] |= LOCKED;
        }
    }

    # ifdef DEBUGSmoothing
    const label lockedFacesI = mesh_.addFaceSubset("lockedFaces");
    forAll(lockedFaces_, lfI)
        mesh_.addFaceToSubset(lockedFacesI, lockedFaces_[lfI]);

    const label lockPointsI = mesh_.addPointSubset("lockedPoints");
    forAll(vertexLocation_, pointI)
        if( vertexLocation_[pointI] & LOCKED )
            mesh_.addPointToSubset(lockPointsI, pointI);
    # endif
}

template<class labelListType>
inline void meshOptimizer::lockPoints(const labelListType& lp)
{
    forAll(lp, lpI)
        vertexLocation_[lp[lpI]] |= LOCKED;

    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();
        const faceListPMG& faces = mesh_.faces();

        forAll(procBoundaries, patchI)
        {
            labelLongList dataToSend;

            const label start = procBoundaries[patchI].patchStart();
            const label end = start+procBoundaries[patchI].patchSize();

            for(label faceI=start;faceI<end;++faceI)
            {
                const face& f = faces[faceI];

                forAll(f, pI)
                {
                    if( vertexLocation_[pI] & LOCKED )
                    {
                        // send the face lable and the location in the face
                        dataToSend.append(faceI-start);
                        dataToSend.append((f.size()-pI)%f.size());
                    }
                }
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dataToSend.byteSize()
            );

            toOtherProc << dataToSend;
        }

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            labelList receivedData;
            fromOtherProc >> receivedData;

            label counter(0);
            while( counter < receivedData.size() )
            {
                const face& f = faces[start+receivedData[counter++]];
                vertexLocation_[f[receivedData[counter++]]] |= LOCKED;
            }
        }
    }

    # ifdef DEBUGSmoothing
    const label lockPointsI = mesh_.addPointSubset("lockedPoints");
    forAll(vertexLocation_, pointI)
        if( vertexLocation_[pointI] & LOCKED )
            mesh_.addPointToSubset(lockPointsI, pointI);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

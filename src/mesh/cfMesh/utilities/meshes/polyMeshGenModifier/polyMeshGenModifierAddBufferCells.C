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
#include "polyMeshGenAddressing.H"
#include "demandDrivenData.H"
#include "labelledPoint.H"
#include "helperFunctions.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::addBufferCells()
{
    if( !Pstream::parRun() )
        return;

    Info << "Adding buffer layers" << endl;

    const labelList& owner = mesh_.owner();
    pointFieldPMG& points = mesh_.points();
    faceListPMG& faces = facesAccess();
    const cellListPMG& cells = mesh_.cells();
    const PtrList<processorBoundaryPatch>& procBoundaries = mesh_.procBoundaries();
    const polyMeshGenAddressing& addressing = mesh_.addressingData();
    const labelLongList& globalPointLabel = addressing.globalPointLabel();
    const Map<label>& globalToLocal = addressing.globalToLocalPointAddressing();

    //- receive vertices
    forAll(procBoundaries, patchI)
    {
        labelHashSet pointsToSend;

        label faceI = procBoundaries[patchI].patchStart();
        const label end = faceI + procBoundaries[patchI].patchSize();
        for(;faceI<end;++faceI)
        {
            const cell& c = cells[owner[faceI]];
            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                forAll(f, pI)
                    pointsToSend.insert(f[pI]);
            }
        }

        faceI = 0;
        List<labelledPoint> ptsToSend(pointsToSend.size());
        forAllConstIter(labelHashSet, pointsToSend, it)
            ptsToSend[faceI++] =
            labelledPoint
            (
                globalPointLabel[it.key()],
                points[it.key()]
            );

        OPstream toOtherProc
        (
            Pstream::commsTypes::blocking,
            procBoundaries[patchI].neiProcNo()
        );

        toOtherProc << ptsToSend;
    }

    Map<label> globalToLocalReceived;
    forAll(procBoundaries, patchI)
    {
        List<labelledPoint> receivedPoints;
        IPstream fromOtherProc
        (
            Pstream::commsTypes::blocking,
            procBoundaries[patchI].neiProcNo()
        );

        fromOtherProc >> receivedPoints;

        forAll(receivedPoints, i)
        {
            if( globalToLocal.found(receivedPoints[i].pointLabel()) )
                continue;
            if( globalToLocalReceived.found(receivedPoints[i].pointLabel()) )
                continue;

            globalToLocalReceived.insert
            (
                receivedPoints[i].pointLabel(),
                points.size()
            );
            points.append(receivedPoints[i].coordinates());
        }
    }

    //- send cells to other processors
    forAll(procBoundaries, patchI)
    {
        labelHashSet cellsToSend;

        label faceI = procBoundaries[patchI].patchStart();
        const label end = faceI + procBoundaries[patchI].patchSize();
        for(;faceI<end;++faceI)
            cellsToSend.insert(owner[faceI]);

        labelLongList flattenedCells;
        forAllConstIter(labelHashSet, cellsToSend, it)
        {
            const cell& c = cells[it.key()];

            //- the number of faces in the cell
            flattenedCells.append(c.size());

            //- add faces
            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                //- the number of vertices in the face
                flattenedCells.append(f.size());
                forAll(f, pI)
                    flattenedCells.append(globalPointLabel[f[pI]]);
            }
        }

        OPstream toOtherProc
        (
            Pstream::commsTypes::blocking,
            procBoundaries[patchI].neiProcNo()
        );

        toOtherProc << flattenedCells;
    }

    forAll(procBoundaries, patchI)
    {
        word subsetName = "processor_";
        subsetName += help::scalarToText(procBoundaries[patchI].neiProcNo());
        const label subsetID = mesh_.addCellSubset(subsetName);

        labelList receivedCells;

        IPstream fromOtherProc
        (
            Pstream::commsTypes::blocking,
            procBoundaries[patchI].neiProcNo()
        );

        fromOtherProc >> receivedCells;

        label counter(0);
        while( counter < receivedCells.size() )
        {
            faceList c(receivedCells[counter++]);
            forAll(c, fI)
            {
                c[fI].setSize(receivedCells[counter++]);
                forAll(c[fI], pI)
                {
                    const label gpI = receivedCells[counter++];

                    if( globalToLocal.found(gpI) )
                    {
                        c[fI][pI] = globalToLocal[gpI];
                    }
                    else
                    {
                        c[fI][pI] = globalToLocalReceived[gpI];
                    }
                }
            }

            mesh_.addCellToSubset(subsetID, cells.size());
            addCell(c);
        }
    }

    mesh_.clearOut();

    Info << "Finished adding buffer layers" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

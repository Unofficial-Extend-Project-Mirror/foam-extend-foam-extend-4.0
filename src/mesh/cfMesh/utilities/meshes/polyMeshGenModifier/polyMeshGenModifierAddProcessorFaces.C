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

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::addProcessorFaces
(
    const VRWGraph& procFaces,
    const labelLongList& facePatches
)
{
    Info << "Adding processor faces" << endl;

    PtrList<processorBoundaryPatch>& procBoundaries = mesh_.procBoundaries_;

    labelList nAddedFaces(procBoundaries.size(), 0);
    forAll(facePatches, fI)
        ++nAddedFaces[facePatches[fI]];

    labelList newPatchStart(procBoundaries.size());
    newPatchStart[0] = procBoundaries[0].patchStart();
    for(label i=1;i<procBoundaries.size();++i)
        newPatchStart[i] =
            newPatchStart[i-1] +
            procBoundaries[i-1].patchSize() + nAddedFaces[i-1];

    //- set new size of the faceListPMG
    faceListPMG& faces = mesh_.faces_;
    const label nFaces = faces.size();
    faces.setSize(nFaces+procFaces.size());

    label endProcFaces(0);
    forAllReverse(procBoundaries, patchI)
    {
        const processorBoundaryPatch& wp = procBoundaries[patchI];
        endProcFaces = Foam::max(endProcFaces, wp.patchStart()+wp.patchSize());
    }

    //- move faces to their new positions
    labelLongList newFaceLabel(nFaces, -1);

    if( endProcFaces != nFaces )
    {
        for(label faceI=nFaces-1;faceI>=endProcFaces;--faceI)
        {
            newFaceLabel[faceI] = faceI+facePatches.size();
            faces[faceI+facePatches.size()].transfer(faces[faceI]);
        }
    }

    labelList faceIndex(procBoundaries.size());
    for(label patchI=procBoundaries.size()-1;patchI>=0;--patchI)
    {
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + procBoundaries[patchI].patchSize();
        const label shift = newPatchStart[patchI] - start;

        if( shift != 0 )
        {
            for(label faceI=end-1;faceI>=start;--faceI)
            {
                faces[faceI+shift].transfer(faces[faceI]);
                newFaceLabel[faceI] = faceI+shift;
            }
        }

        //- set new start for the given patch
        procBoundaries[patchI].patchStart() = newPatchStart[patchI];
        faceIndex[patchI] =
            newPatchStart[patchI] + procBoundaries[patchI].patchSize();
        procBoundaries[patchI].patchSize() += nAddedFaces[patchI];
    }

    //- add new faces into patches
    forAll(procFaces, fI)
    {
        face f(procFaces.sizeOfRow(fI));
        forAll(f, pI)
            f[pI] = procFaces(fI, pI);

        faces[faceIndex[facePatches[fI]]++].transfer(f);
    }

    //- renumber cells
    cellListPMG& cells = mesh_.cells_;
    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided)
    # endif
    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        forAll(c, fI)
            if( newFaceLabel[c[fI]] != -1 )
                c[fI] = newFaceLabel[c[fI]];
    }

    this->clearOut();
    mesh_.clearOut();
    mesh_.updateFaceSubsets(newFaceLabel);

    Info << "Finished adding processor faces" << endl;
}

label polyMeshGenModifier::addProcessorPatch(const label otherProcLabel)
{
    const label nProcPatches = mesh_.procBoundaries().size();

    PtrList<processorBoundaryPatch>& procBoundaries =
        this->procBoundariesAccess();

    procBoundaries.setSize(nProcPatches + 1);

    std::ostringstream ss;
    ss << Pstream::myProcNo();
    std::ostringstream ssNei;
    ssNei << otherProcLabel;
    const word name("processor"+ss.str()+"to"+ssNei.str());

    procBoundaries.set
    (
        nProcPatches,
        new processorBoundaryPatch
        (
            name,
            "processor",
            0,
            0,
            Pstream::myProcNo(),
            otherProcLabel
        )
    );

    return nProcPatches;
}

bool polyMeshGenModifier::removeEmptyProcessorPatches()
{
    PtrList<processorBoundaryPatch>& procBoundaries =
        this->procBoundariesAccess();

    label nValidPatches(0);
    forAll(procBoundaries, patchI)
    {
        if( procBoundaries[patchI].patchSize() != 0 )
            ++nValidPatches;
    }

    if( nValidPatches == procBoundaries.size() )
        return false;

    PtrList<processorBoundaryPatch> newProcBoundaries(nValidPatches);

    nValidPatches = 0;
    forAll(procBoundaries, patchI)
    {
        if( procBoundaries[patchI].patchSize() != 0 )
        {
            newProcBoundaries.set
            (
                nValidPatches++,
                new processorBoundaryPatch(procBoundaries[patchI])
            );
        }
    }

    procBoundaries.transfer(newProcBoundaries);

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

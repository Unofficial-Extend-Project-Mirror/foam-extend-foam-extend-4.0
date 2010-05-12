/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    boundary faces
    - use pointCells when searching for connectivity
    - initialize the cell connectivity with '-1'
    - find both cell faces corresponding to the baffles and mark them
      to prevent a connection
    - standard connectivity checks

\*---------------------------------------------------------------------------*/

#include "meshReader.H"
#include "Time.H"
#include "polyPatch.H"
#include "emptyPolyPatch.H"
#include "preservePatchTypes.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshReader::addPolyBoundaryFace
(
    const label cellId,
    const label cellFaceId,
    const label nCreatedFaces
)
{
#ifdef DEBUG_BOUNDARY
    Info<< nCreatedFaces
    << " add bnd for cell " << cellId
    << " face " << cellFaceId
    << " (original cell " << origCellId_[cellId] << ")"
    << endl;
#endif

    if (cellFaceId < 0)
    {
        FatalErrorIn
        (
            "void meshReader::addPolyBoundaryFace\n"
            "(\n"
            "    const label cellId,\n"
            "    const label cellFaceId,\n"
            "    const label nCreatedFaces\n"
            ")"
        )   << "Trying to add non-existing boundary face: " << cellFaceId
            << " for cell " << cellId
            << abort(FatalError);
    }

    // standard case: volume cells
    const face& thisFace = cellFaces_[cellId][cellFaceId];

    // Debugging
    if (cellPolys_[cellId][cellFaceId] > nInternalFaces_)
    {
        Info<< "meshReader::createPolyBoundary(): "
            << "Problem with face: " << thisFace << endl
            << "Probably multiple definitions "
            << "of a single boundary face." << endl
            << endl;
    }
    else if (cellPolys_[cellId][cellFaceId] >= 0)
    {
        Info<< "meshReader::createPolyBoundary(): "
            << "Problem with face: " << thisFace << endl
            << "Probably trying to define a boundary face "
            << "on a previously matched internal face." << endl
            << "Internal face: "
            << meshFaces_[cellPolys_[cellId][cellFaceId]]
            << endl;
    }

    meshFaces_[nCreatedFaces] = thisFace;
    cellPolys_[cellId][cellFaceId] = nCreatedFaces;
}


void meshReader::createPolyBoundary()
{
    label nBoundaryFaces = 0;
    label nMissingFaces = 0;
    label nInterfaces = 0;

    const faceListList & cFaces = cellFaces();

    // determine number of non-patched faces:
    forAll(cellPolys_, cellI)
    {
        cell & curCell = cellPolys_[cellI];

        forAll(curCell, fI)
        {
            if (curCell[fI] < 0)
            {
                nMissingFaces++;
            }
        }
    }

    forAll(boundaryCells_, patchI)
    {
        nBoundaryFaces += boundaryCells_[patchI].size();
    }

    Info<< nl << "There are " << nMissingFaces
        << " faces to be patched and " << nBoundaryFaces
        << " specified - collect missed boundaries to final patch" << endl;

    patchStarts_.setSize(boundaryCells_.size());
    patchSizes_.setSize(boundaryCells_.size());

    label nCreatedFaces = nInternalFaces_;
    label baffleOffset  = cFaces.size();
    interfaces_.setSize(baffleCellIds_.size());
    nBoundaryFaces = 0;

    forAll(boundaryCells_, patchI)
    {
        const labelList& cList = boundaryCells_[patchI];
        const labelList& fList = boundaryFaces_[patchI];

        patchStarts_[patchI] = nCreatedFaces;

        // Write each baffle side separately
        if (patchPhysicalTypes_[patchI] == "baffle")
        {
            label count = 0;

            for (label side = 0; side < 2; ++side)
            {
                label position = nInterfaces;

                forAll(cList, bndI)
                {
                    label baffleI = cList[bndI] - baffleOffset;

                    if
                    (
                        baffleI >= 0
                     && baffleI < baffleFaces_.size()
                     && baffleCellIds_[baffleI].size()
                    )
                    {
                        addPolyBoundaryFace
                        (
                            baffleCellIds_[baffleI][side],
                            baffleFaceIds_[baffleI][side],
                            nCreatedFaces
                        );

                        // Remove applied boundaries
                        if (side == 1)
                        {
                            baffleCellIds_[baffleI].clear();
                            baffleFaceIds_[baffleI].clear();
                        }

                        interfaces_[position][side] = nCreatedFaces;

                        nBoundaryFaces++;
                        nCreatedFaces++;
                        position++;
                        count++;
                    }
                }
            }

            nInterfaces += (count - (count % 2)) / 2;
        }
        else
        {
            forAll(cList, bndI)
            {
                label cellId = cList[bndI];

                // standard case: volume cells
                if (cellId < baffleOffset)
                {
                    addPolyBoundaryFace
                    (
                        cellId,
                        fList[bndI],
                        nCreatedFaces
                    );

                    nBoundaryFaces++;
                    nCreatedFaces++;
                }
            }
        }

        patchSizes_[patchI] = nCreatedFaces - patchStarts_[patchI];
    }

    // add in missing faces
    Info<< "Missing faces added to patch after face "
        << nCreatedFaces << ":" <<endl;
    nMissingFaces = 0;

    // look for baffles first - keep them together at the start of the patch
    for (label side = 0; side < 2; ++side)
    {
        label position = nInterfaces;

        forAll(baffleCellIds_, baffleI)
        {
            if (baffleCellIds_[baffleI].size())
            {
                // Add each side for each baffle
                addPolyBoundaryFace
                (
                    baffleCellIds_[baffleI][side],
                    baffleFaceIds_[baffleI][side],
                    nCreatedFaces
                );

                interfaces_[position][side] = nCreatedFaces;

                // Remove applied boundaries
                if (side == 1)
                {
                    baffleCellIds_[baffleI].clear();
                    baffleFaceIds_[baffleI].clear();
                }

                nMissingFaces++;
                nCreatedFaces++;
                position++;
            }
        }
    }

    nInterfaces += (nMissingFaces - (nMissingFaces % 2)) / 2;

    // scan for any other missing faces
    forAll(cellPolys_, cellI)
    {
        const labelList& curFaces = cellPolys_[cellI];

        forAll(curFaces, cellFaceI)
        {
            if (curFaces[cellFaceI] < 0)
            {
                // Just report the first 10 or so
                if (nMissingFaces < 10)
                {
                    const face & thisFace = cFaces[cellI][cellFaceI];

                    Info<< "  cell " << cellI << " face " << cellFaceI
                        << " (original cell " << origCellId_[cellI] << ")"
                        << " - Face: " << thisFace
                        << endl;
                }
                else if (nMissingFaces == 10)
                {
                    Info<< "  ..." << nl << endl;
                }

                addPolyBoundaryFace(cellI, cellFaceI, nCreatedFaces);
                nMissingFaces++;
                nCreatedFaces++;
            }
        }
    }

    Info << "Added " << nMissingFaces << " unmatched faces" << endl;

    if (nMissingFaces > 0)
    {
        patchSizes_[patchSizes_.size() - 1] = nMissingFaces;
    }
    else
    {
        patchStarts_.setSize(patchStarts_.size() - 1);
    }

    // reset the size of the face list
    meshFaces_.setSize(nCreatedFaces);

    // check the mesh for face mismatch
    // (faces addressed once or more than twice)
    labelList markupFaces(meshFaces_.size(), 0);

    forAll(cellPolys_, cellI)
    {
        const labelList& curFaces = cellPolys_[cellI];

        forAll(curFaces, faceI)
        {
            markupFaces[curFaces[faceI]]++;
        }
    }

    for (label i = nInternalFaces_; i < markupFaces.size(); i++)
    {
        markupFaces[i]++;
    }

    label nProblemFaces = 0;

    forAll(markupFaces, faceI)
    {
        if (markupFaces[faceI] != 2)
        {
            const face& problemFace = meshFaces_[faceI];

            Info << "meshReader::createPolyBoundary() : "
                << "problem with face " << faceI << ": addressed "
                << markupFaces[faceI] << " times (should be 2!). Face: "
                << problemFace << endl;

            nProblemFaces++;
        }
    }

    if (nProblemFaces > 0)
    {
        Info<< "Number of incorrectly matched faces: "
            << nProblemFaces << endl;
    }

    // adjust for missing members
    if (nInterfaces < interfaces_.size())
    {
        interfaces_.setSize(nInterfaces);
    }

    Info << "Number of boundary faces: " << nBoundaryFaces << endl;
    Info << "Total number of faces: " << nCreatedFaces << endl;
    Info << "Number of interfaces: " << nInterfaces << endl;
}


List<polyPatch*> meshReader::polyBoundaryPatches(const polyMesh& pMesh)
{
    List<polyPatch*> p(patchStarts_.size());

    // Default boundary patch types
    word defaultFacesType(emptyPolyPatch::typeName);

    preservePatchTypes
    (
        runTime_,
        runTime_.constant(),
        polyMesh::defaultRegion,
        patchNames_,
        patchTypes_,
        defaultFacesType,
        patchPhysicalTypes_
    );

    forAll(patchStarts_, patchI)
    {
        p[patchI] = polyPatch::New
        (
            patchTypes_[patchI],
            patchNames_[patchI],
            patchSizes_[patchI],
            patchStarts_[patchI],
            patchI,
            pMesh.boundaryMesh()
        ).ptr();
    }

    return p;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

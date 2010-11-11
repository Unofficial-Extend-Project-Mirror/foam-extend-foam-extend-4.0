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

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMesh::initMesh()
{
    if (debug)
    {
        Info<< "void polyMesh::initMesh() : "
            << "initialising primitiveMesh" << endl;
    }

    // For backward compatibility check if the owner array is the same
    // length as the number of allFaces and shrink to remove the -1s padding
    if (min(owner_) < 0)
    {
        label nActiveFaces = 0;

        forAll(owner_, faceI)
        {
            if (owner_[faceI] == -1)
            {
                break;
            }
            else
            {
                nActiveFaces++;
            }
        }

        InfoIn("void polyMesh::initMesh()")
            << "Truncating owner list at " << nActiveFaces
            << " for backward compatibility" << endl;

        owner_.setSize(nActiveFaces);
    }

    // For backward compatibility check if the neighbour array is the same
    // length as the owner and shrink to remove the -1s padding
    if (min(neighbour_) < 0)
    {
        label nIntFaces = 0;

        forAll(neighbour_, faceI)
        {
            if (neighbour_[faceI] == -1)
            {
                break;
            }
            else
            {
                nIntFaces++;
            }
        }

        InfoIn("void polyMesh::initMesh()")
            << "Truncating neighbour list at " << nIntFaces
            << " for backward compatibility" << endl;

        neighbour_.setSize(nIntFaces);
    }

    label nCells = -1;

    forAll(owner_, facei)
    {
        nCells = max(nCells, owner_[facei]);
    }

    // The neighbour array may or may not be the same length as the owner
    forAll(neighbour_, facei)
    {
        nCells = max(nCells, neighbour_[facei]);
    }

    nCells++;

    label nUsedFaces = 0;
    label nIntFaces = 0;

    // Use patch info if provided, use all faces otherwise
    if (boundary_.size())
    {
        nUsedFaces =
            boundary_[boundary_.size() - 1].start()
          + boundary_[boundary_.size() - 1].size();
        nIntFaces = boundary_[0].start();
    }
    else
    {
        // No patch info. Assume all faces are used.
        nUsedFaces = owner_.size();
        nIntFaces = neighbour_.size();
    }


    label nUsedPoints = allPoints_.size();

    // If not all faces are being used, check that all unused
    // faces are at the back of the list and check the number of
    // used points.
    if (nUsedFaces < allFaces_.size())
    {
        if (debug)
        {
            Info<< "void polyMesh::initMesh() : "
                << "unused faces detected.  "
                << "Number of used faces: " << nUsedFaces
                << ".  Total number of faces: " << allFaces_.size() << endl;
        }

        // Mark the points used by live faces.
        boolList usedPoints(allPoints_.size(), false);

        for (label faceI = 0; faceI < nUsedFaces; faceI++)
        {
            const face& curFace = allFaces_[faceI];

            forAll (curFace, pointI)
            {
                usedPoints[curFace[pointI]] = true;
            }
        }

        forAll (usedPoints, pointI)
        {
            if (!usedPoints[pointI])
            {
                nUsedPoints = pointI;
                break;
            }
        }

        if (nUsedPoints < allPoints_.size())
        {
            if (debug)
            {
                Info<< "void polyMesh::initMesh() : unused points "
                    << "detected.  Number of used points: "
                    << nUsedPoints << ". Total number of points: "
                    << allPoints_.size() << endl;
            }

            for (label i = nUsedPoints; i < allPoints_.size(); i++)
            {
                if (usedPoints[i])
                {
                    FatalErrorIn("void polyMesh::initMesh()")
                        << "Error in point ordering: mixed used and unused "
                        << "points at the end of point list." << nl
                        << "Last used point: " << nUsedPoints
                        << " " << allPoints_[nUsedPoints] << nl
                        << "First unused point: " << nUsedPoints + 1
                        << " " << allPoints_[nUsedPoints + 1] << nl
                        << "and point " << i << " " << allPoints_[i]
                        << " is used by a live face." << endl;

                    // Find out which face uses this point
                    // Problem point = i
                    for (label faceI = 0; faceI < nUsedFaces; faceI++)
                    {
                        const face& curFace = allFaces_[faceI];

                        forAll (curFace, pointI)
                        {
                            if (curFace[pointI] == i)
                            {
                                FatalError
                                    << "Face " << faceI << " " << curFace
                                    << " with points "
                                    << curFace.points(allPoints_)
                                    << endl;
                                break;
                            }
                        }
                    }

                    allPoints_.write();
                    allFaces_.write();
                    owner_.write();
                    neighbour_.write();
                    boundary_.write();
                    FatalError << "Done. " << abort(FatalError);
                }
            }
        }
    }

    // Set sliced lists
    points_.reset(allPoints_, nUsedPoints);
    faces_.reset(allFaces_, owner_.size());

    // Reset the primitiveMesh with the sizes of the primitive arrays
    primitiveMesh::reset
    (
        nUsedPoints,
        neighbour_.size(),
        owner_.size(),
        nCells
    );

    string meshInfo =
        "nPoints: " + Foam::name(nPoints())
      + " nCells: " + Foam::name(this->nCells())
      + " nFaces: " + Foam::name(nFaces())
      + " nInternalFaces: " + Foam::name(nInternalFaces());

    owner_.note() = meshInfo;
    neighbour_.note() = meshInfo;
}


void Foam::polyMesh::initMesh(cellList& c)
{
    if (debug)
    {
        Info<< "void polyMesh::initMesh(cellList& c) : "
            << "calculating owner-neighbour arrays" << endl;
    }

    owner_.setSize(allFaces_.size(), -1);
    neighbour_.setSize(allFaces_.size(), -1);

    boolList markedFaces(allFaces_.size(), false);

    label nInternalFaces = 0;

    forAll(c, cellI)
    {
        // get reference to face labels for current cell
        const labelList& cellfaces = c[cellI];

        forAll (cellfaces, faceI)
        {
            if (!markedFaces[cellfaces[faceI]])
            {
                // First visit: owner
                owner_[cellfaces[faceI]] = cellI;
                markedFaces[cellfaces[faceI]] = true;
            }
            else
            {
                // Second visit: neighbour
                neighbour_[cellfaces[faceI]] = cellI;
                nInternalFaces++;
            }
        }
    }

    // The neighbour array is initialised with the same length as the owner
    // padded with -1s and here it is truncated to the correct size of the
    // number of internal faces.
    neighbour_.setSize(nInternalFaces);

    label nUsedPoints = allPoints_.size();
    label nUsedFaces = allFaces_.size();

    forAll (owner_, faceI)
    {
        if (owner_[faceI] < 0)
        {
            nUsedFaces = faceI;
            break;
        }
    }

    // If not all faces are being used, check that all unused
    // faces are at the back of the list and check the number of
    // used points.
    if (nUsedFaces < owner_.size())
    {
        if (debug)
        {
            Info<< "void polyMesh::initMesh(cellList& c) : "
                << "unused faces detected.  "
                << "Number of used faces: " << nUsedFaces
                << ".  Total number of faces: " << owner_.size() << endl;
        }

        for (label i = nUsedFaces; i < owner_.size(); i++)
        {
            if (owner_[i] >= 0)
            {
                FatalErrorIn("void polyMesh::initMesh(cellList& c)")
                    << "Error in face ordering: mixed used and unused "
                    << "faces at the end of face list." << nl
                    << "Number of used faces: " << nUsedFaces
                    << "  and face " << i << " is owned by cell " << owner_[i]
                    << abort(FatalError);
            }
        }

        // Reset the size of owner array to the number of live faces
        owner_.setSize(nUsedFaces);

        // Mark the points used by live faces.
        boolList usedPoints(allPoints_.size(), false);

        for (label faceI = 0; faceI < nUsedFaces; faceI++)
        {
            const face& curFace = allFaces_[faceI];

            forAll (curFace, pointI)
            {
                usedPoints[curFace[pointI]] = true;
            }
        }

        forAll (usedPoints, pointI)
        {
            if (!usedPoints[pointI])
            {
                nUsedPoints = pointI;
                break;
            }
        }

        if (nUsedPoints < allPoints_.size())
        {
            if (debug)
            {
                Info<< "void polyMesh::initMesh(cellList& c) : unused points "
                    << "detected.  Number of used points: "
                    << nUsedPoints << ". Total number of points: "
                    << allPoints_.size() << endl;
            }

            for (label i = nUsedPoints; i < allPoints_.size(); i++)
            {
                if (usedPoints[i])
                {
                    FatalErrorIn("void polyMesh::initMesh(cellList& c)")
                        << "Error in point ordering: mixed used and unused "
                        << "points at the end of point list." << nl
                        << "Number of used points: " << nUsedPoints
                        << "  and point " << i
                        << " is used by a live face." << endl;

                    // Find out which face uses this point
                    // Problem point = i
                    for (label faceI = 0; faceI < nUsedFaces; faceI++)
                    {
                        const face& curFace = allFaces_[faceI];

                        forAll (curFace, pointI)
                        {
                            if (curFace[pointI] == i)
                            {
                                FatalError
                                    << "Face " << faceI << " " << curFace
                                    << endl;
                                break;
                            }
                        }
                    }

                    FatalError << "Done. " << abort(FatalError);
                }
            }
        }
    }


    // Reset the primitiveMesh
    primitiveMesh::reset
    (
        nUsedPoints,
        neighbour_.size(),
        owner_.size(),
        c.size(),
        c
    );

    string meshInfo =
        "nPoints: " + Foam::name(nPoints())
      + " nCells: " + Foam::name(nCells())
      + " nFaces: " + Foam::name(nFaces())
      + " nInternalFaces: " + Foam::name(this->nInternalFaces());

    owner_.note() = meshInfo;
    neighbour_.note() = meshInfo;
}


// ************************************************************************* //

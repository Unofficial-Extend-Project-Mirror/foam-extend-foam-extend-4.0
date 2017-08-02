/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "nearWallDist.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "cellDistFuncs.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nearWallDist::doAll()
{
    // Correct all cells with face on wall
    const volVectorField& cellCentres = mesh_.C();

    // HR 12.02.18: Use hashSet to determine nbs
    // This should removes a possible error due to wrong sizing since
    // getPointNeighbours may still run over AND should be faster
    // since the linear search (again getPointNeighbours) is removed.
    labelHashSet nbs(20);

    // Correct all cells with face on wall
    forAll(mesh_.boundary(), patchI)
    {
        fvPatchScalarField& ypatch = operator[](patchI);

        const fvPatch& patch = mesh_.boundary()[patchI];

        if (patch.isWall())
        {
            const polyPatch& pPatch = patch.patch();
            const pointField& points = pPatch.points();
            const unallocLabelList& faceCells = pPatch.faceCells();

            // Check cells with face on wall
            forAll(patch, patchFaceI)
            {
                const face& f = pPatch.localFaces()[patchFaceI];

                scalar minDist = GREAT;

                // Loop over points
                forAll(f, fI)
                {
                    const labelList& pointNbs = pPatch.pointFaces()[f[fI]];

                    // Loop over faces sharing current point
                    // This will include the face itself
                    forAll(pointNbs, pointNbsI)
                    {
                        const label nbr = pointNbs[pointNbsI];
                        if (nbs.insert(nbr))
                        {
                            const pointHit curHit = pPatch[nbr].nearestPoint
                            (
                                cellCentres[faceCells[nbr]],
                                points
                            );

                            if (curHit.distance() < minDist)
                            {
                                minDist = curHit.distance();
                            }
                        }
                    }
                }

                ypatch[patchFaceI] = minDist;

                nbs.clear();
            }
        }
        else
        {
            ypatch = 0.0;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearWallDist::nearWallDist(const Foam::fvMesh& mesh)
:
    volScalarField::GeometricBoundaryField
    (
        mesh.boundary(),
        mesh.V(),           // Dummy internal field,
        calculatedFvPatchScalarField::typeName
    ),
    mesh_(mesh)
{
    doAll();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nearWallDist::~nearWallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nearWallDist::correct()
{
    if (mesh_.changing())
    {
        // Update size of GeometricBoundaryField
        forAll(mesh_.boundary(), patchI)
        {
            operator[](patchI).setSize(mesh_.boundary()[patchI].size());
        }
    }

    doAll();
}


// ************************************************************************* //

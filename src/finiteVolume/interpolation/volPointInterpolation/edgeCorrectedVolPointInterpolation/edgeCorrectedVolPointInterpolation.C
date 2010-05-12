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

\*---------------------------------------------------------------------------*/

#include "edgeCorrectedVolPointInterpolation.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::edgeCorrectedVolPointInterpolation::makeExtrapolationVectors() const
{
    // Interpolation:
    // For each point to be corrected calculate the extrapolated value
    // for each face around it and combine them using the inverse
    // distance weighting factors

    extrapolationVectors

    if (extrapolationVectorsPtr_)
    {
        FatalErrorIn
        (
            "void edgeCorrectedVolPointInterpolation::"
            "makeExtrapolationVectors() const"
        )   << "extrapolation vectors already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        Info<< "void edgeCorrectedVolPointInterpolation::"
            << "makeExtrapolationVectors() const : "
            << "constructing extrapolation vectors"
            << endl;
    }

    const labelList& ptc = boundaryPoints();

    // Calculate the correction vectors

    extrapolationVectors_.clear();

    extrapolationVectors_.setSize(ptc.size());

    const labelListList& pf = vMesh().pointFaces();

    const volVectorField& centres = vMesh().C();

    const fvBoundaryMesh& bm = vMesh().boundary();

    forAll (ptc, pointI)
    {
        const label curPoint = ptc[pointI];

        const labelList& curFaces = pf[curPoint];

        extraVecs.hook(new vectorField(curFaces.size()));
        vectorField& curExtraVectors = extraVecs[pointI];

        label nFacesAroundPoint = 0;

        const vector& pointLoc = vMesh().points()[curPoint];

        // Go through all the faces
        forAll (curFaces, faceI)
        {
            if (!vMesh().isInternalFace(curFaces[faceI]))
            {
                // This is a boundary face.  If not in the empty patch
                // or coupled calculate the extrapolation vector
                label patchID =
                    vMesh().boundaryMesh().whichPatch(curFaces[faceI]);

                if
                (
                    !isA<emptyFvPatch>(bm[patchID])
                 && !bm[patchID].coupled()
                )
                {
                    // Found a face for extrapolation
                    curExtraVectors[nFacesAroundPoint] =
                        pointLoc
                      - centres.boundaryField()[patchID]
                            [bm[patchID].patch().whichFace(curFaces[faceI])];

                    nFacesAroundPoint++;
                }
            }
        }

        curExtraVectors.setSize(nFacesAroundPoint);
    }

    if (debug)
    {
        Info<< "void edgeCorrectedVolPointInterpolation::"
            << "makeExtrapolationVectors() const : "
            << "finished constructing extrapolation vectors"
            << endl;
    }
}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::edgeCorrectedVolPointInterpolation::edgeCorrectedVolPointInterpolation
(
    const fvMesh& vm,
    const pointMesh& pm
)
:
    volPointInterpolation(vm, pm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::edgeCorrectedVolPointInterpolation::~edgeCorrectedVolPointInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::edgeCorrectedVolPointInterpolation::updateMesh()
{
    makeExtrapolationVectors();
    return volPointInterpolation::movePoints();
}


// Do what is neccessary if the mesh has moved
bool Foam::edgeCorrectedVolPointInterpolation::movePoints()
{
    makeExtrapolationVectors();
    return volPointInterpolation::movePoints();
}


// ************************************************************************* //

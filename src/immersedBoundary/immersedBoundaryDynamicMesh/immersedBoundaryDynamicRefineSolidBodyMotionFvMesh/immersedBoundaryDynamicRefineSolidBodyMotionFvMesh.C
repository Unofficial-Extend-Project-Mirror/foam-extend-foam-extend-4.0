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

#include "immersedBoundaryDynamicRefineSolidBodyMotionFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvCFD.H"
#include "syncTools.H"
#include "pointFields.H"
#include "directTopoChange.H"
#include "immersedBoundaryPolyPatch.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug
(
    immersedBoundaryDynamicRefineSolidBodyMotionFvMesh,
    0
);

addToRunTimeSelectionTable
(
    dynamicFvMesh,
    immersedBoundaryDynamicRefineSolidBodyMotionFvMesh,
    IOobject
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryDynamicRefineSolidBodyMotionFvMesh::
immersedBoundaryDynamicRefineSolidBodyMotionFvMesh(const IOobject& io)
:
    dynamicPolyRefinementFvMesh(io),
    ibMotions_()
{
    // Read motion function for all regions
    dictionary dynamicMeshCoeffs
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    );

    PtrList<entry> motionDicts(refinementDict().lookup("motionFunctions"));

    ibMotions_.setSize(motionDicts.size());

    forAll (motionDicts, mI)
    {
        ibMotions_.set
        (
            mI,
            new movingImmersedBoundary
            (
                motionDicts[mI].keyword(),
                *this,
                motionDicts[mI].dict()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

immersedBoundaryDynamicRefineSolidBodyMotionFvMesh::
~immersedBoundaryDynamicRefineSolidBodyMotionFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool immersedBoundaryDynamicRefineSolidBodyMotionFvMesh::update()
{
    // Handling multiple calls in a single time step
    if (!firstUpdate())
    {
        // This is not the first call to update, simply return false
        return false;
    }
    else
    {
        // Grab old volumes before moving the mesh
        // This MUST be followed by mesh motion.  HJ, 29/Dec/2017
        setV0();
    }

    forAll (ibMotions_, ibI)
    {
        ibMotions_[ibI].movePoints();
    }

    bool hasChanged = dynamicPolyRefinementFvMesh::update();

    // If the background mesh has not changed, execute dummy mesh sync
    // and mesh motion to re-calculate immersed boundary parameters,
    // since the immersed boundary has been moved
    // HJ, 17/May/2018
    if (!hasChanged)
    {
        fvMesh::syncUpdateMesh();
    }

    // Execute dummy mesh motion for the background mesh
    const pointField oldPoints = allPoints();
    fvMesh::movePoints(oldPoints);

    return true;
}


// ************************************************************************* //

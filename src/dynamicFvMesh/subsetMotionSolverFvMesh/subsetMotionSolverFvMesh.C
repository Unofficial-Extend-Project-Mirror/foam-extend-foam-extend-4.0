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

#include "subsetMotionSolverFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "cellSet.H"
#include "motionSolver.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(subsetMotionSolverFvMesh, 0);

    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        subsetMotionSolverFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subsetMotionSolverFvMesh::subsetMotionSolverFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    movingMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    subsetMesh_
    (
        IOobject
        (
            "motion",
            io.time().constant(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    motionPtr_(NULL),
    alpha_(readScalar(movingMeshCoeffs_.lookup("alpha")))
{
    // Create subset
    word setName = movingMeshCoeffs_.lookup("set");

    cellSet currentSet
    (
        IOobject
        (
            setName,
            io.time().constant(),
            polyMesh::meshSubDir/"sets",
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    subsetMesh_.setLargeCellSubset(currentSet, -1);

    // Create motion solver on the subset
    motionPtr_ = motionSolver::New(subsetMesh_.subMesh());

    // Read motion under-relaxation
    if (alpha_ < 0 || alpha_ > 1.0)
    {
        FatalErrorIn
        (
            "subsetMotionSolverFvMesh::subsetMotionSolverFvMesh"
            "(const IOobject& io)"
        )   << "Ill-defined motion under-relaxation: should be between 0 and 1."
            << "  Alpha = " << alpha_ << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::subsetMotionSolverFvMesh::~subsetMotionSolverFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvMeshSubset& Foam::subsetMotionSolverFvMesh::subsetMesh() const
{
    return subsetMesh_;
}


bool Foam::subsetMotionSolverFvMesh::update()
{
    // Get the points from the moving part
    pointField subsetPoints = motionPtr_->newPoints();

    //- Get copy of mesh points
    pointField p = points();

    //- Map the moving part
    const labelList& subsetPointAddr = subsetMesh_.pointMap();

    forAll (subsetPoints, subsetI)
    {
        p[subsetPointAddr[subsetI]] = subsetPoints[subsetI];
    }

    subsetMesh_.subMesh().movePoints(subsetPoints);

    // Under-relax mesh motion
    p = alpha_*p + (1 - alpha_)*points();

    fvMesh::movePoints(p);

    // Mesh motion only - return false
    return false;
}


// ************************************************************************* //

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
    Virtual base class for mesh motion solver.

\*---------------------------------------------------------------------------*/

#include "tetDecompositionMotionSolver.H"
#include "tetFec.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tetDecompositionMotionSolver, 0);
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::tetDecompositionMotionSolver::applyConstraints
(
    tetFemVectorMatrix& matrix
)
{
    forAll (fixedPoints_, i)
    {
        matrix.addConstraint(fixedPoints_[i], fixedVelocity_[i]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetDecompositionMotionSolver::tetDecompositionMotionSolver
(
    const polyMesh& mesh
)
:
    motionSolver(mesh),
    tetMesh_(mesh),
    motionU_
    (
        IOobject
        (
            "motionU",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        tetMesh_
    ),
    fixedPoints_(),
    fixedVelocity_(),
    totDisplacementPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tetDecompositionMotionSolver::~tetDecompositionMotionSolver()
{
    deleteDemandDrivenData(totDisplacementPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::tetDecompositionMotionSolver::curPoints() const
{
    // Process current point positions

    // Grab all point locations
    tmp<pointField> tcurPoints(new pointField(tetMesh_().allPoints()));
    pointField& cp = tcurPoints();

    // Move live points from mesh motion
    vectorField mp
    (
        vectorField::subField
        (
            motionU_.internalField(),
            tetMesh_().nPoints()
        )*tetMesh_.time().deltaT().value()
    );

    // Note: mp is the size of live points, while cp contains all points.
    //       Looping over the smaller array.  HJ, 5/Sep/2009
    forAll (mp, i)
    {
        cp[i] += mp[i];
    }

    twoDCorrectPoints(tcurPoints());

    return tcurPoints;
}


void Foam::tetDecompositionMotionSolver::setConstraint
(
    const label pointID,
    const vector& fixedVel
)
{
    // Add constraint to the list
    fixedPoints_.append(pointID);
    fixedVelocity_.append(fixedVel);
}


void Foam::tetDecompositionMotionSolver::clearConstraints()
{
    fixedPoints_.clear();
    fixedVelocity_.clear();
}


void Foam::tetDecompositionMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    tetPolyMeshMapper mapper(tetMesh_, mpm);
    tetMesh_.updateMesh(mapper);

    // Reset the internal field of the motion variable to avoid mapping isses
    motionU_.internalField() = vector::zero;

    motionSolver::updateMesh(mpm);
}


Foam::tmp<Foam::elementScalarField>
Foam::tetDecompositionMotionSolver::distortionEnergy() const
{
    tmp<elementScalarField> tUd
    (
        new elementScalarField
        (
            IOobject
            (
                "distortionEnergy",
                tetMesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh(),
            dimensionedScalar("0.0", dimless, 0.0)
        )
    );

    elementScalarField& Ud = tUd();

    elementTensorField gradU =
        tetFec::elementGrad(motionU_)*tetMesh().time().deltaT();

    Ud = (1.0/2.0)*((gradU && gradU) + (gradU && gradU.T()))
       - (1.0/3.0)*tr(gradU)*tr(gradU);

    return tUd;
}


Foam::tmp<Foam::elementScalarField>
Foam::tetDecompositionMotionSolver::deformationEnergy() const
{
    tmp<elementScalarField> tUd
    (
        new elementScalarField
        (
            IOobject
            (
                "deforamationEnergy",
                tetMesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh(),
            dimensionedScalar("0.0", dimless, 0.0)
        )
    );

    elementScalarField& Ud = tUd();

    elementTensorField gradU =
        tetFec::elementGrad(motionU_)*tetMesh().time().deltaT();

    scalar nu = 0.25;

    Ud = (1.0/2.0)*((gradU&&gradU) + (gradU&&gradU.T()));
       + (nu/(1.0+2.0*nu))*tr(gradU)*tr(gradU);

    return tUd;
}



Foam::tmp<Foam::elementScalarField>
Foam::tetDecompositionMotionSolver::totDistortionEnergy() const
{
    tmp<elementScalarField> tUd
    (
        new elementScalarField
        (
            IOobject
            (
                "totDeformationEnergy",
                tetMesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh(),
            dimensionedScalar("0.0", dimless, 0.0)
        )
    );

    elementScalarField& Ud = tUd();

    if (!needTotDisplacement())
    {
        FatalErrorIn
        (
            "tetDecompositionMotionSolver::totDeformationEnergy()"
        )   << "Total displacement field is not stored "
            << "in tetDecompositionMotionSolver object." << endl
            << exit(FatalError);
    }

    elementTensorField gradU = tetFec::elementGrad(totDisplacement());

    Ud = (1.0/2.0)*((gradU&&gradU) + (gradU&&gradU.T()))
       - (1.0/3.0)*tr(gradU)*tr(gradU);

    return tUd;
}



Foam::tmp<Foam::elementScalarField>
Foam::tetDecompositionMotionSolver::totDeformationEnergy() const
{
    tmp<elementScalarField> tUd
    (
        new elementScalarField
        (
            IOobject
            (
                "totDistortionEnergy",
                tetMesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh(),
            dimensionedScalar("0.0", dimless, 0.0)
        )
    );

    elementScalarField& Ud = tUd();

    if (!needTotDisplacement())
    {
        FatalErrorIn
        (
            "tetDecompositionMotionSolver::totDistortionEnergy()"
        )   << "Total displacement field is not stored." << endl
            << exit(FatalError);
    }

    elementTensorField gradU = tetFec::elementGrad(totDisplacement());

    scalar nu = 0.25;

    Ud = (1.0/2.0)*((gradU && gradU) + (gradU && gradU.T()));
       + (nu/(1.0+2.0*nu))*tr(gradU)*tr(gradU);

    return tUd;
}


// ************************************************************************* //

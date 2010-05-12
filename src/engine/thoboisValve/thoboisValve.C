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

#include "thoboisValve.H"
#include "engineTime.H"
#include "polyMesh.H"
#include "interpolateXY.H"
#include "IFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::thoboisValve::adjustCrankAngle(const scalar theta) const
{
    if (theta < liftProfileStart_)
    {
        scalar adjustedTheta = theta;

        while (adjustedTheta < liftProfileStart_)
        {
            adjustedTheta += liftProfileEnd_ - liftProfileStart_;
        }

        return adjustedTheta;
    }
    else if (theta > liftProfileEnd_)
    {
        scalar adjustedTheta = theta;

        while (adjustedTheta > liftProfileEnd_)
        {
            adjustedTheta -= liftProfileEnd_ - liftProfileStart_;
        }

        return adjustedTheta;
    }
    else
    {
        return theta;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::thoboisValve::thoboisValve
(
    const word& name,
    const polyMesh& mesh,
    const autoPtr<coordinateSystem>& valveCS,
    const word& bottomPatchName,
    const word& poppetPatchName,
    const word& sidePatchName,
    const word& stemPatchName,
    const word& detachInCylinderPatchName,
    const word& detachInPortPatchName,
    const word& detachFacesName,
    const graph& liftProfile,
    const scalar minLift,
    const scalar diameter,
    const word& staticPointsName,
    const word& movingPointsName,
    const word& movingInternalPointsName,
    const word& staticCellsName,
    const word& movingCellsName
)
:
    name_(name),
    mesh_(mesh),
    engineDB_(refCast<const engineTime>(mesh.time())),
    csPtr_(valveCS),
    bottomPatch_(bottomPatchName, mesh.boundaryMesh()),
    poppetPatch_(poppetPatchName, mesh.boundaryMesh()),
    sidePatch_(sidePatchName, mesh.boundaryMesh()),
    stemPatch_(stemPatchName, mesh.boundaryMesh()),
    detachInCylinderPatch_(detachInCylinderPatchName, mesh.boundaryMesh()),
    detachInPortPatch_(detachInPortPatchName, mesh.boundaryMesh()),
    detachFacesName_(detachFacesName),
    liftProfile_(liftProfile),
    liftProfileStart_(min(liftProfile_.x())),
    liftProfileEnd_(max(liftProfile_.x())),
    minLift_(minLift),
    diameter_(diameter),
    staticPointsName_(staticPointsName),
    movingPointsName_(movingPointsName),
    movingInternalPointsName_(movingInternalPointsName),
    staticCellsName_(staticCellsName),
    movingCellsName_(movingCellsName)
{}


// Construct from dictionary
Foam::thoboisValve::thoboisValve
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    engineDB_(refCast<const engineTime>(mesh_.time())),
    csPtr_
    (
        coordinateSystem::New
        (
            "coordinateSystem",
            dict.subDict("coordinateSystem")
        )
    ),
    bottomPatch_(dict.lookup("bottomPatch"), mesh.boundaryMesh()),
    poppetPatch_(dict.lookup("poppetPatch"), mesh.boundaryMesh()),
    sidePatch_(dict.lookup("sidePatch"), mesh.boundaryMesh()),
    stemPatch_(dict.lookup("stemPatch"), mesh.boundaryMesh()),
    detachInCylinderPatch_
    (
        dict.lookup("detachInCylinderPatch"),
        mesh.boundaryMesh()
    ),
    detachInPortPatch_
    (
        dict.lookup("detachInPortPatch"),
        mesh.boundaryMesh()
    ),
    detachFacesName_(dict.lookup("detachFaces")),
    liftProfile_
    (
        "theta",
        "lift",
        name_,
        IFstream
        (
            mesh.time().path()/mesh.time().constant()/
            word(dict.lookup("liftProfileFile"))
        )()
    ),
    liftProfileStart_(min(liftProfile_.x())),
    liftProfileEnd_(max(liftProfile_.x())),
    minLift_(readScalar(dict.lookup("minLift"))),
    diameter_(readScalar(dict.lookup("diameter"))),
    staticPointsName_(dict.lookup("staticPoints")),
    movingPointsName_(dict.lookup("movingPoints")),
    movingInternalPointsName_(dict.lookup("movingInternalPoints")),
    staticCellsName_(dict.lookup("staticCells")),
    movingCellsName_(dict.lookup("movingCells"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::thoboisValve::lift(const scalar theta) const
{
    return interpolateXY
    (
        adjustCrankAngle(theta),
        liftProfile_.x(),
        liftProfile_.y()
    );
}


bool Foam::thoboisValve::isOpen() const
{
    return lift(engineDB_.theta()) >= minLift_;
}


Foam::scalar Foam::thoboisValve::curLift() const
{
    return max
    (
        lift(engineDB_.theta()),
        minLift_
    );
}


Foam::scalar Foam::thoboisValve::curVelocity() const
{
    return
       -(
             curLift()
           - max
             (
                 lift(engineDB_.theta() - engineDB_.deltaTheta()),
                 minLift_
             )
        )/(engineDB_.deltaT().value() + VSMALL);
}


Foam::labelList Foam::thoboisValve::movingPatchIDs() const
{
    labelList mpIDs(2);
    label nMpIDs = 0;

    if (bottomPatch_.active())
    {
        mpIDs[nMpIDs] = bottomPatch_.index();
        nMpIDs++;
    }

    if (poppetPatch_.active())
    {
        mpIDs[nMpIDs] = poppetPatch_.index();
        nMpIDs++;
    }

    mpIDs.setSize(nMpIDs);

    return mpIDs;
}


void Foam::thoboisValve::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK;

    cs().writeDict(os);

    os  << "bottomPatch " << bottomPatch_.name() << token::END_STATEMENT << nl
        << "poppetPatch " << poppetPatch_.name() << token::END_STATEMENT << nl
        << "sidePatch " << sidePatch_.name() << token::END_STATEMENT << nl
        << "stemPatch " << stemPatch_.name() << token::END_STATEMENT << nl
        << "detachInCylinderPatch " << detachInCylinderPatch_.name()
        << token::END_STATEMENT << nl
        << "detachInPortPatch " << detachInPortPatch_.name()
        << token::END_STATEMENT << nl
        << "detachFaces " << detachFacesName_ << token::END_STATEMENT << nl
        << "liftProfile " << nl << token::BEGIN_BLOCK
        << liftProfile_ << token::END_BLOCK << token::END_STATEMENT << nl
        << "minLift " << minLift_ << token::END_STATEMENT << nl
        << "diameter " << diameter_ << token::END_STATEMENT << nl
        << "staticCells " << staticCellsName_ << token::END_STATEMENT << nl
        << "movingCells " << movingCellsName_ << token::END_STATEMENT << nl
        << "staticPoints " << staticPointsName_ << token::END_STATEMENT << nl
        << "movingPoints " << movingPointsName_ << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// ************************************************************************* //

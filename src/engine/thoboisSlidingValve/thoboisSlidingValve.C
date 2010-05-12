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

#include "thoboisSlidingValve.H"
#include "engineTime.H"
#include "polyMesh.H"
#include "interpolateXY.H"
#include "IFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::thoboisSlidingValve::adjustCrankAngle(const scalar theta) const
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
Foam::thoboisSlidingValve::thoboisSlidingValve
(
    const word& name,
    const polyMesh& mesh,
    const autoPtr<coordinateSystem>& valveCS,
    const word& bottomPatchName,
    const word& poppetPatchName,
    const word& sidePatchName,
    const word& stemPatchName,
    const word& curtainInPortPatchName,
    const word& curtainInCylinderPatchName,
    const word& detachInCylinderPatchName,
    const word& detachInPortPatchName,
    const labelList& detachFaces,
    const graph& liftProfile,
    const scalar minLift,
    const scalar diameter,
    const scalar deformationLift,
    const word& layeringFacesTopName,
    const word& layeringFacesBottomName,
    const word& movingCellsTopName,
    const word& movingCellsBottomName,
    const word& movingPointsTopName,
    const word& movingPointsBottomName,
    const scalar minTopLayer,
    const scalar maxTopLayer,
    const scalar minBottomLayer,
    const scalar maxBottomLayer,
    const word& staticPointsName,
    const word& movingPointsName,
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
    curtainInCylinderPatch_(curtainInPortPatchName, mesh.boundaryMesh()),
    curtainInPortPatch_(curtainInPortPatchName, mesh.boundaryMesh()),
    detachInCylinderPatch_(detachInCylinderPatchName, mesh.boundaryMesh()),
    detachInPortPatch_(detachInPortPatchName, mesh.boundaryMesh()),
    detachFaces_(detachFaces),
    liftProfile_(liftProfile),
    liftProfileStart_(min(liftProfile_.x())),
    liftProfileEnd_(max(liftProfile_.x())),
    minLift_(minLift),
    diameter_(diameter),
    deformationLift_(deformationLift),
    layeringFacesTopName_(layeringFacesTopName),
    layeringFacesBottomName_(layeringFacesBottomName),
    movingCellsTopName_(movingCellsTopName),
    movingCellsBottomName_(movingCellsBottomName),
    movingPointsTopName_(movingPointsTopName),
    movingPointsBottomName_(movingPointsBottomName),
    minTopLayer_(minTopLayer),
    maxTopLayer_(maxTopLayer),
    minBottomLayer_(minBottomLayer),
    maxBottomLayer_(maxBottomLayer),
    staticPointsName_(staticPointsName),
    movingPointsName_(movingPointsName),
    staticCellsName_(staticCellsName),
    movingCellsName_(movingCellsName)
{}


// Construct from dictionary
Foam::thoboisSlidingValve::thoboisSlidingValve
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
    curtainInCylinderPatch_
    (
        dict.lookup("curtainInCylinderPatch"),
        mesh.boundaryMesh()
    ),
    curtainInPortPatch_
    (
        dict.lookup("curtainInPortPatch"),
        mesh.boundaryMesh()
    ),
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
    detachFaces_(dict.lookup("detachFaces")),
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
    deformationLift_(readScalar(dict.lookup("deformationLift"))),    
    layeringFacesTopName_(dict.lookup("layeringFacesTop")),
    layeringFacesBottomName_(dict.lookup("layeringFacesBottom")),
    movingCellsTopName_(dict.lookup("movingCellsTop")),
    movingCellsBottomName_(dict.lookup("movingCellsBottom")),
    movingPointsTopName_(dict.lookup("movingPointsTop")),
    movingPointsBottomName_(dict.lookup("movingPointsBottom")),
    minTopLayer_(readScalar(dict.lookup("minTopLayer"))),
    maxTopLayer_(readScalar(dict.lookup("maxTopLayer"))),
    minBottomLayer_(readScalar(dict.lookup("minBottomLayer"))),
    maxBottomLayer_(readScalar(dict.lookup("maxBottomLayer"))),
    staticPointsName_(dict.lookup("staticPoints")),
    movingPointsName_(dict.lookup("movingPoints")),
    staticCellsName_(dict.lookup("staticCells")),
    movingCellsName_(dict.lookup("movingCells"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::thoboisSlidingValve::lift(const scalar theta) const
{
    return interpolateXY
    (
        adjustCrankAngle(theta),
        liftProfile_.x(),
        liftProfile_.y()
    );
}


bool Foam::thoboisSlidingValve::isOpen() const
{
    return lift(engineDB_.theta()) >= minLift_;
}


Foam::scalar Foam::thoboisSlidingValve::curLift() const
{
    return max
    (
        lift(engineDB_.theta()),
        minLift_
    );
}


Foam::scalar Foam::thoboisSlidingValve::curVelocity() const
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


Foam::labelList Foam::thoboisSlidingValve::movingPatchIDs() const
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


void Foam::thoboisSlidingValve::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK;

    cs().writeDict(os);

    os  << "bottomPatch " << bottomPatch_.name() << token::END_STATEMENT << nl
        << "poppetPatch " << poppetPatch_.name() << token::END_STATEMENT << nl
        << "sidePatch " << sidePatch_.name() << token::END_STATEMENT << nl
        << "stemPatch " << stemPatch_.name() << token::END_STATEMENT << nl
        << "curtainInCylinderPatch " << curtainInCylinderPatch_.name()
        << token::END_STATEMENT << nl
        << "curtainInPortPatch " << curtainInPortPatch_.name()
        << token::END_STATEMENT << nl
//        << "detachInCylinderPatch " << detachInCylinderPatch_.name()
//        << token::END_STATEMENT << nl
        << "liftProfile " << nl << token::BEGIN_BLOCK
        << liftProfile_ << token::END_BLOCK << token::END_STATEMENT << nl
        << "minLift " << minLift_ << token::END_STATEMENT << nl
        << "diameter " << diameter_ << token::END_STATEMENT << nl
         << token::END_BLOCK << endl;
}


// ************************************************************************* //

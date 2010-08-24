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

#include "engineVerticalValve.H"
#include "engineTime.H"
#include "polyMesh.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::engineVerticalValve::engineVerticalValve
(
    const word& name,
    const polyMesh& mesh,
    const autoPtr<coordinateSystem>& valveCS,
    const word& bottomPatchName,
    const word& poppetPatchName,
    const word& stemPatchName,
    const word& curtainInPortPatchName,
    const word& curtainInCylinderPatchName,
    const word& detachInCylinderPatchName,
    const word& detachInPortPatchName,
    const labelList& detachFaces,
    const graph& liftProfile,
    const scalar minLift,
    const scalar minTopLayer,
    const scalar maxTopLayer,
    const scalar minBottomLayer,
    const scalar maxBottomLayer,
    const scalar diameter,
    const word& valveHeadPatchName,
    const scalar topLayerOffset,
    const scalar topLayerTol,
    const scalar bottomLayerOffset,
    const scalar bottomLayerTol,
    const scalar detachDistance,
    const scalar detachTol,
    const scalar deformationLift
)
:
    engineValve
    (
        name,
        mesh,
        valveCS,
        bottomPatchName,
        poppetPatchName,
        stemPatchName,
        curtainInPortPatchName,
        curtainInCylinderPatchName,
        detachInCylinderPatchName,
        detachInPortPatchName,
        detachFaces,
        liftProfile,
        minLift,
        minTopLayer,
        maxTopLayer,
        minBottomLayer,
        maxBottomLayer,
        diameter
    ),
    valveHeadPatch_(valveHeadPatchName, mesh.boundaryMesh()),
    topLayerOffset_(topLayerOffset),
    topLayerTol_(topLayerTol),
    bottomLayerOffset_(bottomLayerOffset),
    bottomLayerTol_(bottomLayerTol),
    detachDistance_(detachDistance),
    detachTol_(detachTol),
    deformationLift_(deformationLift)
{}


// Construct from dictionary
Foam::engineVerticalValve::engineVerticalValve
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    engineValve
    (
        name,
        mesh,
        dict
    ),
    valveHeadPatch_
    (
        dict.lookup("valveHeadPatch"),
        mesh.boundaryMesh()
    ),
    topLayerOffset_(readScalar(dict.lookup("topLayerOffset"))),
    topLayerTol_(readScalar(dict.lookup("topLayerTol"))),
    bottomLayerOffset_(readScalar(dict.lookup("bottomLayerOffset"))),
    bottomLayerTol_(readScalar(dict.lookup("bottomLayerTol"))),
    detachDistance_(readScalar(dict.lookup("detachDistance"))),
    detachTol_(readScalar(dict.lookup("detachTol"))),
    deformationLift_(readScalar(dict.lookup("deformationLift")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


void Foam::engineVerticalValve::writeDict(Ostream& os) const
{
    engineValve::writeDict(os);


    os  << "headPatch " << valveHeadPatch_.name() << token::END_STATEMENT << nl
        << "bottomLayerOffset " << bottomLayerOffset_
        << token::END_STATEMENT << nl
        << "bottomLayerOffsetayerTol " << bottomLayerTol_
       << token::END_STATEMENT << nl
        << "topLayerOffset " << topLayerOffset_ << token::END_STATEMENT << nl
        << "topLayerTol " << topLayerTol_ << token::END_STATEMENT << nl
        << "detachDistance " << detachDistance_ << token::END_STATEMENT << nl
        << "detachTol " << detachTol_ << token::END_STATEMENT << nl
        << "deformationLift " << deformationLift_ << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "contactPatchPair.H"
#include "contactProblem.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::contactPatchPair::contactPatchPair
(
    const word& name,
    const contactProblem& cp,
    const word& masterPatchName,
    const word& slavePatchName,
    const dimensionedScalar& frictionCoeff,
    const scalar contactTol,
    const intersection::algorithm alg,
    const intersection::direction dir
)
:
    name_(name),
    cp_(cp),
    masterPatch_(masterPatchName, cp.mesh().boundaryMesh()),
    slavePatch_(slavePatchName, cp.mesh().boundaryMesh()),
    frictionCoeff_(frictionCoeff),
    contactTol_(contactTol),
    masterInterpolate_
    (
        cp.mesh().boundaryMesh()[masterPatch_.index()]
    ),
    slaveInterpolate_
    (
        cp.mesh().boundaryMesh()[slavePatch_.index()]
    ),
    masterToSlaveInterpolate_
    (
        cp.mesh().boundaryMesh()[masterPatch_.index()],   // from patch
        cp.mesh().boundaryMesh()[slavePatch_.index()],    // to patch
        alg,
        dir
    ),
    slaveToMasterInterpolate_
    (
        cp.mesh().boundaryMesh()[slavePatch_.index()],    // from patch
        cp.mesh().boundaryMesh()[masterPatch_.index()],   // to patch
        alg,
        dir
    )
{}


// Construct from dictionary
Foam::contactPatchPair::contactPatchPair
(
    const word& name,
    const contactProblem& cp,
    const dictionary& dict
)
:
    name_(name),
    cp_(cp),
    masterPatch_(dict.lookup("masterPatch"), cp.mesh().boundaryMesh()),
    slavePatch_(dict.lookup("slavePatch"), cp.mesh().boundaryMesh()),
    frictionCoeff_(dict.lookup("frictionCoeff")),
    contactTol_(readScalar(dict.lookup("contactTol"))),
    masterInterpolate_
    (
        cp.mesh().boundaryMesh()[masterPatch_.index()]
    ),
    slaveInterpolate_
    (
        cp.mesh().boundaryMesh()[slavePatch_.index()]
    ),
    masterToSlaveInterpolate_
    (
        cp.mesh().boundaryMesh()[masterPatch_.index()],    // from patch
        cp.mesh().boundaryMesh()[slavePatch_.index()],     // to patch
        intersection::algorithmNames_.read(dict.lookup("projectionAlgo")),
        intersection::directionNames_.read(dict.lookup("projectionDir"))
        
    ),
    slaveToMasterInterpolate_
    (
        cp.mesh().boundaryMesh()[slavePatch_.index()],     // from patch
        cp.mesh().boundaryMesh()[masterPatch_.index()],    // to patch
        intersection::algorithmNames_.read(dict.lookup("projectionAlgo")),
        intersection::directionNames_.read(dict.lookup("projectionDir"))
        
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::contactPatchPair::masterTouchFraction() const
{
    // Get reference to displacement field and mesh
    const volVectorField& U = cp_.U();
    const fvMesh& mesh = cp_.mesh();

    // Interpolate slave displacement into master vertices
    vectorField masterVertexU =
        slaveToMasterInterpolate_.pointInterpolate<vector>
        (
            slaveInterpolate_.faceToPointInterpolate
            (
                U.boundaryField()[slavePatch_.index()]
            )
        );

    const vectorField& projectionDir =
        mesh.boundaryMesh()[masterPatch_.index()].pointNormals();


    // Calculate master gap function
    scalarField vertexMasterGap =
    (
        (
            masterVertexU
          - masterInterpolate_.faceToPointInterpolate
            (
                U.boundaryField()[masterPatch_.index()]
            )
        )
        & projectionDir
    ) + slaveToMasterInterpolate_.pointDistanceToIntersection() - contactTol_;

    // Calculate area in contact

    const faceList& masterPatchLocalFaces =
        mesh.boundaryMesh()[masterPatch_.index()].localFaces();

    const pointField& masterPatchLocalPoints =
        mesh.boundaryMesh()[masterPatch_.index()].localPoints();

    tmp<scalarField> ttouchFrac
    (
        new scalarField(masterPatchLocalFaces.size(), 0)
    );
    scalarField& touchFrac = ttouchFrac();

    forAll (masterPatchLocalFaces, faceI)
    {
        touchFrac[faceI] =
            masterPatchLocalFaces[faceI].areaInContact
            (
                masterPatchLocalPoints,
                vertexMasterGap
            );
    }

    return ttouchFrac;
}


Foam::tmp<Foam::scalarField>
Foam::contactPatchPair::slaveTouchFraction() const
{
    // Get reference to displacement field and mesh
    const volVectorField& U = cp_.U();
    const fvMesh& mesh = cp_.mesh();

    // Interpolate master displacement into slave vertices
    vectorField slaveVertexU =
        masterToSlaveInterpolate_.pointInterpolate<vector>
        (
            masterInterpolate_.faceToPointInterpolate
            (
                U.boundaryField()[masterPatch_.index()]
            )
        );

    const vectorField& projectionDir =
        mesh.boundaryMesh()[slavePatch_.index()].pointNormals();


    // Calculate slave gap function
    scalarField vertexSlaveGap =
    (
        (
            slaveVertexU
          - slaveInterpolate_.faceToPointInterpolate
            (
                U.boundaryField()[slavePatch_.index()]
            )
        )
        & projectionDir
    ) + masterToSlaveInterpolate_.pointDistanceToIntersection() - contactTol_;

    // Calculate area in contact

    const faceList& slavePatchLocalFaces =
        mesh.boundaryMesh()[slavePatch_.index()].localFaces();

    const pointField& slavePatchLocalPoints =
        mesh.boundaryMesh()[slavePatch_.index()].localPoints();

    tmp<scalarField> ttouchFrac
    (
        new scalarField(slavePatchLocalFaces.size(), 0)
    );
    scalarField& touchFrac = ttouchFrac();

    forAll (slavePatchLocalFaces, faceI)
    {
        touchFrac[faceI] =
            slavePatchLocalFaces[faceI].areaInContact
            (
                slavePatchLocalPoints,
                vertexSlaveGap
            );
    }

    return ttouchFrac;
}


void Foam::contactPatchPair::correct
(
    const FieldField<Field, vector>& curTraction,
    FieldField<Field, vector>& newTraction,
    FieldField<Field, vector>& refValue,
    FieldField<Field, scalar>& valueFraction
)
{
    // Get reference to displacement field and mesh
    const volVectorField::GeometricBoundaryField& Upatches =
        cp_.U().boundaryField();

    const fvMesh& mesh = cp_.mesh();
    const surfaceVectorField::GeometricBoundaryField& Apatches =
        mesh.Sf().boundaryField();
    const surfaceScalarField::GeometricBoundaryField& magApatches =
        mesh.magSf().boundaryField();

    // Get patch indices
    const label masterIndex = masterPatch_.index();
    const label slaveIndex = slavePatch_.index();

    // Calculate patch normals
    vectorField nMasterPatch = Apatches[masterIndex]/magApatches[masterIndex];

    vectorField nSlavePatch = Apatches[slaveIndex]/magApatches[slaveIndex];


    // Calculate slave pressure and tangential force

    scalarField slavePressure = -( nSlavePatch & curTraction[slaveIndex]);

    // Enforce gradient condition on the master patch

    // Calculate relative tangential velocity for master patch
    vectorField relUmaster =
        slaveToMasterInterpolate_.faceInterpolate<vector>
        (
            Upatches[slaveIndex]
        )
      - Upatches[masterIndex];

    relUmaster -= nMasterPatch*(nMasterPatch & relUmaster);
    relUmaster /= mag(relUmaster) + VSMALL;

    // Calculate tangential master traction
    scalarField magMasterTangential =
        Foam::mag((I - nMasterPatch*nMasterPatch) & curTraction[masterIndex]);

    // Calculate master pressure
    scalarField masterPressure =
        max
        (
            slaveToMasterInterpolate_.faceInterpolate<scalar>
            (
                slavePressure
            ),
            0.0
        );

    // Calculate master traction, using the positive part of
    // slave pressure and tangential fricton
    // Mind the signs: pressure = negative gradient (minus master normal)
    //                 friction = positive pressure
    newTraction[masterIndex] +=
        masterTouchFraction()*
        (
           -nMasterPatch*masterPressure
          + relUmaster*
            min
            (
                frictionCoeff_.value()*masterPressure,
                magMasterTangential
            )
        );

    // Enforce direction mixed condition on the slave patch

    // Calculate slave fraction.  Correct for negative pressure
    // (if the pressure is negative, the contact is released)
    //HJ, fiddle pos pressure!!!
    scalarField slaveFrac = slaveTouchFraction();

    // Calculate slave displacement
    vectorField slaveVertexU =
        masterToSlaveInterpolate_.pointInterpolate<vector>
        (
            masterInterpolate_.faceToPointInterpolate
            (
                Upatches[masterIndex]
            )
        );

    const vectorField& projectionDir =
        mesh.boundaryMesh()[slaveIndex].pointNormals();

    // Calculate slave displacement
    vectorField slaveDisp =
        slaveInterpolate_.pointToFaceInterpolate
        (
            slaveVertexU
          + masterToSlaveInterpolate_.pointDistanceToIntersection()
            *projectionDir
        );

    // Accumulate normal of slave displacement
    refValue[slaveIndex] +=
        nSlavePatch*
        min
        (
            pos(slaveFrac)*
            (
                (nSlavePatch & Upatches[slaveIndex])
              + slaveFrac*contactTol_
            ),
            (nSlavePatch & slaveDisp)
        );


    // Accumulate slave friction

    // Calculate relative tangential velocity for slave patch
    vectorField relUslave =
        masterToSlaveInterpolate_.faceInterpolate<vector>
        (
            Upatches[masterIndex]
        )
      - Upatches[slaveIndex];

    relUslave -= nSlavePatch*(nSlavePatch & relUslave);
    relUslave /= mag(relUslave) + VSMALL;

    // Take out normal component out of slave traction and find the
    // magnitude of the tangential traction.
    scalarField magSlaveTangential =
        Foam::mag((I - nSlavePatch*nSlavePatch) & curTraction[slaveIndex]);

    // Calculate slave traction
    newTraction[slaveIndex] +=
        slaveFrac*relUslave*
        min
        (
            frictionCoeff_.value()*max(slavePressure, scalar(0)),
            magSlaveTangential
        );

    // Accumulate slave touch fraction
    valueFraction[slaveIndex] += slaveFrac;

/*
    Info << "slavePressure: " << slavePressure << nl
         << "slaveTouchFrac: " << slaveTouchFraction() << nl
//          << "slaveFrac: " << slaveFrac << nl
         << "refValueSlave: " << refValue[slaveIndex].component(vector::Y) << nl
//          << "slaveTraction: " << newTraction[slaveIndex] << nl
         << "masterTouchFrac: " << masterTouchFraction() << nl
//          << "interpolated slave pressure: "
//          << slaveToMasterInterpolate_.faceInterpolate<scalar>
//             (
//                 slavePressure
//             )
//          << nl
//          << "masterTraction: "
//          << newTraction[masterIndex].component(vector::Y)
         << endl;
*/
}


void Foam::contactPatchPair::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK;

    os  << "masterPatch " << masterPatch_.name() << token::END_STATEMENT << nl
        << "slavePatch " << slavePatch_.name() << token::END_STATEMENT << nl
        << "frictionCoeff " << frictionCoeff_ << token::END_STATEMENT << nl
        << "contactTol " << contactTol_ << token::END_STATEMENT << nl
        << "projectionAlgo "
        << intersection::algorithmNames_
               [masterToSlaveInterpolate_.projectionAlgo()]
        << token::END_STATEMENT << nl
        << "projectionDir "
        << intersection::directionNames_
               [masterToSlaveInterpolate_.projectionDir()]
        << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //

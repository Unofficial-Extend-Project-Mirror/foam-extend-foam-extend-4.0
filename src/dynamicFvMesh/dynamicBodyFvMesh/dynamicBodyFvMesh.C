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

#include "dynamicBodyFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "tetDecompositionMotionSolver.H"
#include "laplaceTetDecompositionMotionSolver.H"
#include "fixedValueTetPolyPatchFields.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicBodyFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicBodyFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicBodyFvMesh::dynamicBodyFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
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
    motionPtr_(motionSolver::New(*this)),
    bodyPatchName_
    (
        dynamicMeshCoeffs_.lookup("bodyPatchName")
    ),
    bodyPatchID_(-1),
    translationDirection_
    (
        dynamicMeshCoeffs_.lookup("translationDirection")
    ),
    translationAmplitude_
    (
        readScalar(dynamicMeshCoeffs_.lookup("translationAmplitude"))
    ),
    translationFrequency_
    (
        readScalar(dynamicMeshCoeffs_.lookup("translationFrequency"))
    ),
    initialRotationOrigin_
    (
        dynamicMeshCoeffs_.lookup("initialRotationOrigin")
    ),
    rotationAxis_
    (
        dynamicMeshCoeffs_.lookup("rotationAxis")
    ),
    rotationAmplitude_
    (
        readScalar(dynamicMeshCoeffs_.lookup("rotationAmplitude"))
    ),
    rotationFrequency_
    (
        readScalar(dynamicMeshCoeffs_.lookup("rotationFrequency"))
    )
{
    bodyPatchID_ = boundaryMesh().findPatchID(bodyPatchName_);

    if(bodyPatchID_<0)
    {
        FatalErrorIn
        (
            "dynamicBodyFvMesh::dynamicBodyFvMesh(const IOobject& io)"
        )   
            << "Can't find patch: " << bodyPatchName_
                << exit(FatalError);
    }

    translationDirection_ /= mag(translationDirection_) + SMALL;

    rotationAxis_ /= mag(rotationAxis_) + SMALL;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicBodyFvMesh::~dynamicBodyFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicBodyFvMesh::update()
{
    scalar curTime = time().value();
    scalar oldTime = curTime - time().deltaT().value();

    if
    (
        motionPtr_->type()
     == laplaceTetDecompositionMotionSolver::typeName
    )
    {
        tetDecompositionMotionSolver& mSolver = 
            dynamic_cast<tetDecompositionMotionSolver&>
            (
                motionPtr_()
            );

        vector trans =
            translationAmplitude_
           *(
                sin(2*mathematicalConstant::pi*translationFrequency_*curTime)
              - sin(2*mathematicalConstant::pi*translationFrequency_*oldTime)
            )
           *translationDirection_;


        scalar rotAngle =
            rotationAmplitude_
           *(
                sin(2*mathematicalConstant::pi*rotationFrequency_*curTime)
              - sin(2*mathematicalConstant::pi*rotationFrequency_*oldTime)
            );

        vector curRotationOrigin = 
            initialRotationOrigin_
          + translationDirection_
           *translationAmplitude_
           *sin(2*mathematicalConstant::pi*translationFrequency_*oldTime);

        const pointField& oldPoints =
            mSolver.tetMesh().boundary()[bodyPatchID_].localPoints();

        vector r0(1, 1, 1);
        r0 -= rotationAxis_*(rotationAxis_ & r0);
        r0 /= mag(r0);

        // http://mathworld.wolfram.com/RotationFormula.html
        vector r1 = 
            r0*cos(rotAngle)
          + rotationAxis_*(rotationAxis_ & r0)*(1 - cos(rotAngle))
          + (r0 ^ rotationAxis_)*sin(rotAngle);

        tensor T = rotationTensor(r0, r1);

        vectorField rot = 
            transform(T, oldPoints - curRotationOrigin) 
          + curRotationOrigin
          - oldPoints;


        tetPointVectorField& motionU = mSolver.motionU();

        if
        (
            motionU.boundaryField()[bodyPatchID_].type()
         == fixedValueTetPolyPatchVectorField::typeName
        )
        {
            fixedValueTetPolyPatchVectorField& motionUBodyPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[bodyPatchID_]
                );

            motionUBodyPatch == (trans + rot)/time().deltaT().value();
        }
        else
        {
            FatalErrorIn("dynamicBodyFvMesh::update()")   
                << "Bounary condition on " << motionU.name() 
                    <<  " for " << bodyPatchName_ << " patch is " 
                    << motionU.boundaryField()[bodyPatchID_].type() 
                    << ", instead " 
                    << fixedValueTetPolyPatchVectorField::typeName
                    << exit(FatalError);
        }
    }
    else
    {
        FatalErrorIn("dynamicBodyFvMesh::update()")   
            << "Selected mesh motion solver is "
                << motionPtr_->type()
                << ", instead " 
                << tetDecompositionMotionSolver::typeName
                << exit(FatalError);
    }

    fvMesh::movePoints(motionPtr_->newPoints());

    // Mesh motion only - return false
    return false;
}


// ************************************************************************* //

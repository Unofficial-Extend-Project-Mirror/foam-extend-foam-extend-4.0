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

#include "unsTotalLagrangianSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGradf.H"
#include "fvcCellLimitedGrad.H"

#include "tractionDisplacementFvPatchVectorField.H"
#include "velocityTractionDisplacementFvPatchVectorField.H"
#include "skewCorrectionVectors.H"
#include "multiMaterial.H"
#include "twoDPointCorrector.H"
#include "thermalModel.H"
#include "componentReferenceList.H"

#include "eig3Field.H"
#include "nonLinearGeometry.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidSolvers
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(unsTotalLagrangianSolid, 0);
addToRunTimeSelectionTable(solidSolver, unsTotalLagrangianSolid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsTotalLagrangianSolid::unsTotalLagrangianSolid(const fvMesh& mesh)
:
    solidSolver(typeName, mesh),
    D_
    (
        IOobject
        (
            "D",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    U_
    (
        IOobject
        (
//             "U",
            "ddt(" + D_.name() + ")",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
//         fvc::ddt(D_)
        mesh,
        dimensionedVector("0", dimVelocity, vector::zero)
    ),
    pMesh_(mesh),
    pointD_
    (
        IOobject
        (
            "pointD",
            runTime().timeName(),
            mesh,
//             IOobject::READ_IF_PRESENT,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_
//         dimensionedVector("0", dimLength, vector::zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    rheology_(sigma_, D_),
    volToPoint_(mesh),
    gradDf_
    (
        IOobject
        (
            "gradDf",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    sigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    gradD_
    (
        IOobject
        (
            "grad(" + D_.name() + ")",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    rho_("rho", rheology_.rho()),
    mu_(rheology_.mu()),
    muf_("muf", fvc::interpolate(mu_)),
    lambda_(rheology_.lambda()),
    lambdaf_("lambdaf", fvc::interpolate(lambda_)),
    interface_(NULL),
    threeKPtr_(NULL),
    alphaPtr_(NULL),
    threeKfPtr_(NULL),
    alphafPtr_(NULL)
{
    pointD_.oldTime();
//     D_.oldTime();
//     U_.oldTime();

    if (rheology_.law().type() == multiMaterial::typeName)
    {
        interface_.set(new TLMaterialInterface(D_, pointD_));
    }

    if (interface().valid())
    {
        muf_ = interface()->interpolate(mu_);
        lambdaf_ = interface()->interpolate(lambda_);

//         interface()->modifyProperty(muf_);
//         interface()->modifyProperty(lambdaf_);
    }

    if (interface().valid())
    {
//         interface()->updateDisplacement(pointD_); // restart issue
        interface()->updateDisplacementGradient(gradD_, gradDf_);
    }
    else
    {
//         volToPoint_.interpolate(D_, pointD_);  // restart issue
        gradD_ = fvc::cellLimitedGrad(D_, pointD_, 0);
        gradDf_ = fvc::fGrad(D_, pointD_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector unsTotalLagrangianSolid::pointU(label pointID) const
{
    pointVectorField pointU
    (
        IOobject
        (
            "pointU",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimVelocity, vector::zero)
    );

    volToPoint_.interpolate(U_, pointU);

    return pointU.internalField()[pointID];
}

//- Patch point displacement
tmp<vectorField> unsTotalLagrangianSolid::patchPointDisplacementIncrement
(
    const label patchID
) const
{
    tmp<vectorField> tPointDisplacement
    (
        new vectorField
        (
            mesh().boundaryMesh()[patchID].localPoints().size(),
            vector::zero
        )
    );

    tPointDisplacement() =
        vectorField
        (
            pointD_.internalField() - pointD_.oldTime().internalField(),
            mesh().boundaryMesh()[patchID].meshPoints()
        );

    return tPointDisplacement;
}

//- Face zone point displacement
tmp<vectorField> unsTotalLagrangianSolid::faceZonePointDisplacementIncrement
(
    const label zoneID
) const
{
    tmp<vectorField> tPointDisplacement
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().localPoints().size(),
            vector::zero
        )
    );
    vectorField& pointDisplacement = tPointDisplacement();

    const vectorField& pointDI = pointD_.internalField();
    const vectorField& oldPointDI = pointD_.oldTime().internalField();

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone

        const labelList& curPointMap =
            globalToLocalFaceZonePointMap()[globalZoneIndex];

        const labelList& zoneMeshPoints =
            mesh().faceZones()[zoneID]().meshPoints();

        vectorField zonePointsDisplGlobal
        (
            zoneMeshPoints.size(),
            vector::zero
        );

        //- Inter-proc points are shared by multiple procs
        //  pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(zoneMeshPoints.size(), 0);

        forAll(zonePointsDisplGlobal, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            if(zoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = zoneMeshPoints[localPoint];

                zonePointsDisplGlobal[globalPointI] =
                    pointDI[procPoint] - oldPointDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }

        if (Pstream::parRun())
        {
            reduce(zonePointsDisplGlobal, sumOp<vectorField>());
            reduce(pointNumProcs, sumOp<scalarField>());

            //- now average the displacement between all procs
            zonePointsDisplGlobal /= pointNumProcs;
        }

        forAll(pointDisplacement, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            pointDisplacement[localPoint] =
                zonePointsDisplGlobal[globalPointI];
        }
    }
    else
    {
        tPointDisplacement() =
            vectorField
            (
                pointDI - oldPointDI,
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    return tPointDisplacement;
}


//- Patch point displacement
tmp<vectorField> unsTotalLagrangianSolid::patchPointDisplacement
(
    const label patchID
) const
{
    tmp<vectorField> tPointDisplacement
    (
        new vectorField
        (
            mesh().boundaryMesh()[patchID].localPoints().size(),
            vector::zero
        )
    );

    tPointDisplacement() =
        vectorField
        (
            pointD_.oldTime().internalField(),
            mesh().boundaryMesh()[patchID].meshPoints()
        );

    return tPointDisplacement;
}

//- Face zone point displacement
tmp<vectorField> unsTotalLagrangianSolid::faceZonePointDisplacement
(
    const label zoneID
) const
{
    tmp<vectorField> tPointDisplacement
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().localPoints().size(),
            vector::zero
        )
    );
    vectorField& pointDisplacement = tPointDisplacement();

    const vectorField& oldPointDI = pointD_.oldTime().internalField();

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone

        const labelList& curPointMap =
            globalToLocalFaceZonePointMap()[globalZoneIndex];

        const labelList& zoneMeshPoints =
            mesh().faceZones()[zoneID]().meshPoints();

        vectorField zonePointsDisplGlobal
        (
            zoneMeshPoints.size(),
            vector::zero
        );

        //- Inter-proc points are shared by multiple procs
        //  pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(zoneMeshPoints.size(), 0);

        forAll(zonePointsDisplGlobal, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            if(zoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = zoneMeshPoints[localPoint];

                zonePointsDisplGlobal[globalPointI] =
                    oldPointDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }

        if (Pstream::parRun())
        {
            reduce(zonePointsDisplGlobal, sumOp<vectorField>());
            reduce(pointNumProcs, sumOp<scalarField>());

            //- now average the displacement between all procs
            zonePointsDisplGlobal /= pointNumProcs;
        }

        forAll(pointDisplacement, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            pointDisplacement[localPoint] =
                zonePointsDisplGlobal[globalPointI];
        }
    }
    else
    {
        tPointDisplacement() =
            vectorField
            (
                oldPointDI,
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    return tPointDisplacement;
}

//- Patch face acceleration
tmp<Foam::vectorField> unsTotalLagrangianSolid::patchFaceAcceleration
(
    const label patchID
) const
{
    tmp<vectorField> tAcceleration
    (
        new vectorField
        (
            mesh().boundary()[patchID].size(),
            vector::zero
        )
    );

     volVectorField a = fvc::d2dt2(D_);
//      volVectorField a = fvc::ddt(U_);

    tAcceleration() = a.boundaryField()[patchID];

    return tAcceleration;
}

//- Patch face acceleration
tmp<vectorField> unsTotalLagrangianSolid::faceZoneAcceleration
(
    const label zoneID,
    const label patchID
) const
{
    tmp<vectorField> tAcceleration
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().size(),
            vector::zero
        )
    );
    vectorField& acceleration = tAcceleration();

    volVectorField a = fvc::d2dt2(D_);
//     volVectorField a = fvc::ddt(U_);

    vectorField patchAcceleration = a.boundaryField()[patchID];

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(patchAcceleration, i)
    {
        acceleration[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
            patchAcceleration[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(acceleration, sumOp<vectorField>());

    return tAcceleration;
}


//- Patch face acceleration
tmp<vectorField> unsTotalLagrangianSolid::faceZoneVelocity
(
    const label zoneID,
    const label patchID
) const
{
    tmp<vectorField> tVelocity
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().size(),
            vector::zero
        )
    );
    vectorField& velocity = tVelocity();

    vectorField patchVelocity = U_.boundaryField()[patchID];

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(patchVelocity, i)
    {
        velocity[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
            patchVelocity[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(velocity, sumOp<vectorField>());


//     const vectorField& patchVelocity = U_.boundaryField()[patchID];

//     label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

//     if (globalZoneIndex != -1)
//     {
//         // global face zone

//         const label patchStart =
//             mesh().boundaryMesh()[patchID].start();

//         forAll(patchVelocity, i)
//         {
//             velocity
//             [
//                 mesh().faceZones()[zoneID].whichFace(patchStart + i)
//             ] =
//                 patchVelocity[i];
//         }

//         // Parallel data exchange: collect field on all processors
//         reduce(velocity, sumOp<vectorField>());
//     }
//     else
//     {
//         velocity = patchVelocity;
//     }

    return tVelocity;
}


//- Face zone point displacement
tmp<tensorField> unsTotalLagrangianSolid::faceZoneSurfaceGradientOfVelocity
(
    const label zoneID,
    const label patchID
) const
{
    tmp<tensorField> tVelocityGradient
    (
        new tensorField
        (
            mesh().faceZones()[zoneID]().size(),
            tensor::zero
        )
    );
    tensorField& velocityGradient = tVelocityGradient();

    vectorField pPointU =
        volToPoint_.interpolate(mesh().boundaryMesh()[patchID], U_);

    const faceList& localFaces =
        mesh().boundaryMesh()[patchID].localFaces();

    vectorField localPoints =
        mesh().boundaryMesh()[patchID].localPoints();
    localPoints += pointD_.boundaryField()[patchID].patchInternalField();

    PrimitivePatch<face, List, const pointField&> patch
    (
        localFaces,
        localPoints
    );

    tensorField patchGradU = fvc::fGrad(patch, pPointU);

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone

        const label patchStart =
            mesh().boundaryMesh()[patchID].start();

        forAll(patchGradU, i)
        {
            velocityGradient
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ] =
                patchGradU[i];
        }

        // Parallel data exchange: collect field on all processors
        reduce(velocityGradient, sumOp<tensorField>());
    }
    else
    {
        velocityGradient = patchGradU;
    }

    return tVelocityGradient;
}


tmp<vectorField>
unsTotalLagrangianSolid::currentFaceZonePoints(const label zoneID) const
{
    vectorField pointDisplacement
    (
        mesh().faceZones()[zoneID]().localPoints().size(),
        vector::zero
    );

    const vectorField& pointDI = pointD_.internalField();

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone
        const labelList& curPointMap =
            globalToLocalFaceZonePointMap()[globalZoneIndex];

        const labelList& zoneMeshPoints =
            mesh().faceZones()[zoneID]().meshPoints();

        vectorField zonePointsDisplGlobal
        (
            zoneMeshPoints.size(),
            vector::zero
        );

        //- Inter-proc points are shared by multiple procs
        //  pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(zoneMeshPoints.size(), 0);

        forAll(zonePointsDisplGlobal, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            if(zoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = zoneMeshPoints[localPoint];

                zonePointsDisplGlobal[globalPointI] =
                    pointDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }

        if (Pstream::parRun())
        {
            reduce(zonePointsDisplGlobal, sumOp<vectorField>());
            reduce(pointNumProcs, sumOp<scalarField>());

            //- now average the displacement between all procs
            zonePointsDisplGlobal /= pointNumProcs;
        }

        forAll(pointDisplacement, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            pointDisplacement[localPoint] =
                zonePointsDisplGlobal[globalPointI];
        }
    }
    else
    {
        pointDisplacement =
            vectorField
            (
                pointDI,
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    tmp<vectorField> tCurrentPoints
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().localPoints()
          + pointDisplacement
        )
    );

    return tCurrentPoints;
}


//- Face zone point displacement
tmp<vectorField> unsTotalLagrangianSolid::faceZoneNormal
(
    const label zoneID,
    const label patchID
) const
{
    tmp<vectorField> tNormals
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().size(),
            vector::zero
        )
    );
    vectorField& normals = tNormals();

    const faceList& localFaces =
        mesh().boundaryMesh()[patchID].localFaces();

    vectorField localPoints =
        mesh().boundaryMesh()[patchID].localPoints();
    localPoints += pointD_.boundaryField()[patchID].patchInternalField();

    PrimitivePatch<face, List, const pointField&> patch
    (
        localFaces,
        localPoints
    );

    vectorField patchNormals(patch.size(), vector::zero);

    forAll(patchNormals, faceI)
    {
        patchNormals[faceI] =
            localFaces[faceI].normal(localPoints);
    }

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone

        const label patchStart =
            mesh().boundaryMesh()[patchID].start();

        forAll(patchNormals, i)
        {
            normals
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ] =
                patchNormals[i];
        }

        // Parallel data exchange: collect field on all processors
        reduce(normals, sumOp<vectorField>());
    }
    else
    {
        normals = patchNormals;
    }

    return tNormals;
}

void unsTotalLagrangianSolid::setTraction
(
    const label patchID,
    const vectorField& traction
)
{
    if
    (
        D_.boundaryField()[patchID].type()
     != tractionDisplacementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void unsTotalLagrangianSolid::setTraction(...)")
            << "Bounary condition on " << D_.name()
                <<  " is "
                << D_.boundaryField()[patchID].type()
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead "
                << tractionDisplacementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementFvPatchVectorField& patchU =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            D_.boundaryField()[patchID]
        );

    patchU.traction() = traction;
}

void unsTotalLagrangianSolid::setPressure
(
    const label patchID,
    const scalarField& pressure
)
{
    if
    (
        D_.boundaryField()[patchID].type()
     != tractionDisplacementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void unsTotalLagrangianSolid::setTraction(...)")
            << "Bounary condition on " << D_.name()
                <<  " is "
                << D_.boundaryField()[patchID].type()
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead "
                << tractionDisplacementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementFvPatchVectorField& patchU =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            D_.boundaryField()[patchID]
        );

    patchU.pressure() = pressure;
}

void unsTotalLagrangianSolid::setTraction
(
    const label patchID,
    const label zoneID,
    const vectorField& faceZoneTraction
)
{
    vectorField patchTraction(mesh().boundary()[patchID].size(), vector::zero);

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(patchTraction, i)
    {
        patchTraction[i] =
            faceZoneTraction
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
    }

    setTraction(patchID, patchTraction);
}

void unsTotalLagrangianSolid::setPressure
(
    const label patchID,
    const label zoneID,
    const scalarField& faceZonePressure
)
{
    scalarField patchPressure(mesh().boundary()[patchID].size(), 0.0);

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(patchPressure, i)
    {
        patchPressure[i] =
            faceZonePressure
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
    }

    setPressure(patchID, patchPressure);
}


//- Set traction at specified patch
void unsTotalLagrangianSolid::setVelocityAndTraction
(
    const label patchID,
    const vectorField& traction,
    const vectorField& velocity,
    const vectorField& normal
)
{
    if
    (
        D_.boundaryField()[patchID].type()
     != velocityTractionDisplacementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn
        (
            "void unsTotalLagrangianSolid::"
            "setVelocityAndTraction(...)"
        )
            << "Bounary condition on " << D_.name()
                <<  " is "
                << D_.boundaryField()[patchID].type()
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead "
                << velocityTractionDisplacementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    velocityTractionDisplacementFvPatchVectorField& patchU =
        refCast<velocityTractionDisplacementFvPatchVectorField>
        (
            D_.boundaryField()[patchID]
        );

    patchU.traction() = traction;
    patchU.velocity() = velocity;
    patchU.normal() = normal;
}


//- Set traction at specified patch
void unsTotalLagrangianSolid::setVelocityAndTraction
(
    const label patchID,
    const label zoneID,
    const vectorField& faceZoneTraction,
    const vectorField& faceZoneVelocity,
    const vectorField& faceZoneNormal
)
{
    vectorField patchTraction(mesh().boundary()[patchID].size(), vector::zero);
    vectorField patchVelocity(mesh().boundary()[patchID].size(), vector::zero);
    vectorField patchNormal(mesh().boundary()[patchID].size(), vector::zero);

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(patchTraction, i)
    {
        patchTraction[i] =
            faceZoneTraction
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
        patchVelocity[i] =
            faceZoneVelocity
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
        patchNormal[i] =
            faceZoneNormal
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
    }

    setVelocityAndTraction
    (
        patchID,
        patchTraction,
        patchVelocity,
        patchNormal
    );
}


tmp<vectorField> unsTotalLagrangianSolid::predictTraction
(
    const label patchID,
    const label zoneID
)
{
    // Predict traction on patch
    if
    (
        D_.boundaryField()[patchID].type()
     != tractionDisplacementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void unsTotalLagrangianSolid::setTraction(...)")
            << "Bounary condition on " << D_.name()
                <<  " is "
                << D_.boundaryField()[patchID].type()
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead "
                << tractionDisplacementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementFvPatchVectorField& patchUo =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            D_.oldTime().boundaryField()[patchID]
        );

    tractionDisplacementFvPatchVectorField& patchUoo =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            D_.oldTime().oldTime().boundaryField()[patchID]
        );


    vectorField ptF = 2*patchUo.traction() - patchUoo.traction();

    tmp<vectorField> ttF
    (
        new vectorField(mesh().faceZones()[zoneID].size(), vector::zero)
    );
    vectorField& tF = ttF();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(ptF, i)
    {
        tF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] = ptF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(tF, sumOp<vectorField>());

    return ttF;
}


tmp<scalarField> unsTotalLagrangianSolid::predictPressure
(
    const label patchID,
    const label zoneID
)
{
    // Predict pressure field on patch
    if
    (
        D_.boundaryField()[patchID].type()
     != tractionDisplacementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void unsTotalLagrangianSolid::setTraction(...)")
            << "Bounary condition on " << D_.name()
                <<  " is "
                << D_.boundaryField()[patchID].type()
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead "
                << tractionDisplacementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementFvPatchVectorField& patchUo =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            D_.oldTime().boundaryField()[patchID]
        );

    tractionDisplacementFvPatchVectorField& patchUoo =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            D_.oldTime().oldTime().boundaryField()[patchID]
        );


    scalarField pPF = 2*patchUo.pressure() - patchUoo.pressure();

    tmp<scalarField> tpF
    (
        new scalarField(mesh().faceZones()[zoneID].size(), 0)
    );
    scalarField& pF = tpF();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(pPF, i)
    {
        pF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] = pPF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(pF, sumOp<scalarField>());

    return tpF;
}


bool unsTotalLagrangianSolid::evolve()
{
    Info << "Evolving solid solver: "
        << unsTotalLagrangianSolid::typeName << endl;

    int nCorr
    (
        readInt(solidProperties().lookup("nCorrectors"))
    );

    scalar convergenceTolerance
    (
        readScalar(solidProperties().lookup("convergenceTolerance"))
    );

    scalar relConvergenceTolerance = 0;
    if (solidProperties().found("relConvergenceTolerance"))
    {
        relConvergenceTolerance =
            readScalar(solidProperties().lookup("relConvergenceTolerance"));
    }

    // componentReferenceList cr
    // (
    //     solidProperties().lookup("componentReference"),
    //     componentReference::iNew(mesh())
    // );


    // Non-linear
    nonLinearGeometry::nonLinearType nonLinear =
        nonLinearGeometry::nonLinearNames_.read
        (
            solidProperties().lookup("nonLinear")
        );
    // Update solidMechanics dictionary
    const_cast<dictionary&>
    (
        mesh().solutionDict().subDict("solidMechanics")
    ).set("nonLinear", nonLinearGeometry::nonLinearNames_[nonLinear]);



    Switch debug(solidProperties().lookup("debug"));

    dimensionedScalar K("K", dimless/dimTime, 0);
    if (solidProperties().found("K"))
    {
        K = dimensionedScalar(solidProperties().lookup("K"));
        Info << K << endl;
    }

//     dimensionedVector g("g", dimVelocity/dimTime, vector::zero);
    if (solidProperties().found("g"))
    {
        bodyForce() = dimensionedVector(solidProperties().lookup("g"));
    }

    bool enforceLinear = false;
    solidProperties().set("enforceLinear", enforceLinear);

    if (thermalStress())
    {
        threeK();
        alpha();
        threeKf();
        alphaf();
    }

    int iCorr = 0;
    scalar initialResidual = 0;
    lduSolverPerformance solverPerf;
    scalar res = 1;
    scalar maxRes = 0;
    scalar curConvergenceTolerance = convergenceTolerance;

    lduMatrix::debug = debug;

//     label DRefCell = 0;
//     vector DRefValue = vector::zero;
//     setRefCellVector(D_, solidProperties(), DRefCell, DRefValue);

    OFstream* resFilePtr = NULL;
    if (lduMatrix::debug && runTime().outputTime())
    {
        Info << "Creating residual output file" << endl;
        mkDir(runTime().timePath());
        resFilePtr = new OFstream(runTime().timePath()/"residual.dat");
    }

    do
    {
        if (lduMatrix::debug)
        {
            Info<< "Time: " << runTime().timeName()
                << ", outer iteration: " << iCorr << endl;
        }

        D_.storePrevIter();

//         // Demirdzic-Muzaferija
//         fvVectorMatrix DEqn
//         (
//             rho_*fvm::d2dt2(D_)
//          == fvm::laplacian(muf_, D_, "laplacian(DD,D)")
//           + fvc::div
//             (
//                 mesh().Sf()
//               & (
//                     muf_*gradDf_.T()
//                   + lambdaf_*(I*tr(gradDf_))
//                 )
//             )
//           + rho_*bodyForce()
//         );

        fvVectorMatrix DEqn
        (
            rho_*fvm::d2dt2(D_)
         == fvm::laplacian(2*muf_ + lambdaf_, D_, "laplacian(DD,D)")
          + fvc::div
            (
                mesh().Sf()
              & (
                  - (muf_ + lambdaf_)*gradDf_
                  + muf_*gradDf_.T()
                  + lambdaf_*(I*tr(gradDf_))
                )
            )
          + rho_*bodyForce()
        );

        // Add damping
        if (K.value() > SMALL)
        {
            DEqn += K*rho_*fvm::ddt(D_);
        }

        if (nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            surfaceSymmTensorField Ef =
                symm(gradDf_) + 0.5*symm(gradDf_ & gradDf_.T());

            sigmaf_ = 2*muf_*Ef  + I*(lambdaf_*tr(Ef));

            DEqn -=
                fvc::div
                (
                    muf_*(mesh().Sf() & (gradDf_ & gradDf_.T()))
                  + 0.5*lambdaf_*tr(gradDf_ & gradDf_.T())*mesh().Sf()
                )
              + fvc::div(mesh().Sf() & (sigmaf_ & gradDf_));
        }

        if (thermalStress())
        {
            DEqn +=
                fvc::div
                (
                    mesh().Sf()*threeKf()*alphaf()*DTf()
                );
        }

        if (interface().valid())
        {
            interface()->correct(DEqn);
        }

        // forAll (cr, crI)
        // {
        //     DEqn.setComponentReference
        //     (
        //         cr[crI].patchIndex(),
        //         cr[crI].faceIndex(),
        //         cr[crI].dir(),
        //         cr[crI].value()
        //     );
        // }

        if (mesh().solutionDict().relaxEquation(D_.name()))
        {
            DEqn.relax(mesh().solutionDict().fieldRelaxationFactor(D_.name()));
        }
        // else
        // {
        //     DEqn.relax(1);
        // }

        solverPerf = DEqn.solve();

        D_.relax();

        if(iCorr == 0)
        {
            initialResidual = solverPerf.initialResidual();
        }

        if (interface().valid())
        {
            interface()->updateDisplacement(pointD_);
            interface()->updateDisplacementGradient(gradD_, gradDf_);
        }
        else
        {
            volToPoint_.interpolate(D_, pointD_);
            gradD_ = fvc::cellLimitedGrad(D_, pointD_, lsCellGrad(), 0);
            gradDf_ = fvc::fGrad(D_, pointD_, interpFaceGrad());
        }

        if (nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            surfaceScalarField Det = det(I+gradDf_);

            scalar minDetFf = min(Det).value();
            reduce(minDetFf, minOp<scalar>());

            scalar maxDetFf = max(Det).value();
            reduce(maxDetFf, maxOp<scalar>());

            if ( (minDetFf<0.01) || (maxDetFf>100) )
            {
                Info << "det: " << minDetFf << ", " << maxDetFf << endl;

                const_cast<dictionary&>
                (
                    mesh().solutionDict().subDict("solidMechanics")
                ).set("nonLinear", word("off"));

                nonLinear =
                    nonLinearGeometry::nonLinearNames_.read
                    (
                        mesh().solutionDict().subDict
                        (
                            "solidMechanics"
                        ).lookup("nonLinear")
                    );
            }
        }

        // Calculate relative momentum residual
        res = residual();

        if (res > maxRes)
        {
            maxRes = res;
        }

        curConvergenceTolerance = maxRes*relConvergenceTolerance;
        if (curConvergenceTolerance < convergenceTolerance)
        {
            curConvergenceTolerance = convergenceTolerance;
        }

        if (lduMatrix::debug)
        {
            Info << "Relative residual = " << res << endl;

            if (resFilePtr)
            {
                (*resFilePtr) << iCorr+1 << " " << res << endl;
            }
        }

        if (iCorr==0)
        {
            res = 1;
        }
    }
    while
    (
        res > curConvergenceTolerance
     && ++iCorr < nCorr
    );

    // Update number of outter correctors
    nCorrectors() = iCorr;

    U_.oldTime();
    U_ = fvc::ddt(D_);

    // Calculate second Piola-Kirchhoff stress
    {
        epsilon_ = symm(gradD_);
        surfaceSymmTensorField epsilonf = symm(gradDf_);

        if(nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            epsilon_ += 0.5*symm(gradD_ & gradD_.T());
            epsilonf += 0.5*symm(gradDf_ & gradDf_.T());
        }

        sigma_ = 2*mu_*epsilon_ + I*(lambda_*tr(epsilon_));
        sigmaf_ = 2*muf_*epsilonf + I*(lambdaf_*tr(epsilonf));

        if (thermalStress())
        {
            sigma_ -= symmTensor(1,0,0,1,0,1)*threeK()*alpha()*DT();
            sigmaf_ -= symmTensor(1,0,0,1,0,1)*threeKf()*alphaf()*DTf();
        }
    }

    Info << solverPerf.solverName() << ": Solving for " << D_.name()
        << ", Initial residula = " << initialResidual
        << ", Final residual = " << solverPerf.initialResidual()
        << ", No outer iterations = " << iCorr
        << "\nMax relative residual = " << maxRes
        << ", Relative residual = " << res
        << ", enforceLinear = " << enforceLinear << endl;

    lduMatrix::debug = 1;

    nonLinearGeometry::nonLinearType nonLinearOriginal =
        nonLinearGeometry::nonLinearNames_.read
        (
            solidProperties().lookup("nonLinear")
        );

    if (nonLinear != nonLinearOriginal)
    {
        return false;
    }

    return true;
}


tmp<volScalarField> unsTotalLagrangianSolid::hydPressure() const
{
    Switch nonLinear(solidProperties().lookup("nonLinear"));

    if (nonLinear)
    {
        volTensorField F = I + gradD_.T();
        volScalarField J = det(F);
        volSymmTensorField sigmaCauchy
        (
            "sigmaCauchy",
            (1/J) * symm(F & sigma_ & F.T())
        );

        tmp<volScalarField> tHydPressure
        (
            new volScalarField
            (
                tr(sigmaCauchy)/3
            )
        );

        return tHydPressure;
    }

    tmp<volScalarField> tHydPressure
    (
        new volScalarField
        (
            tr(sigma_)/3
        )
    );

    return tHydPressure;
}


void unsTotalLagrangianSolid::predict()
{
    Info << "Predicting solid" << endl;

    D_ = D_ + U_*runTime().deltaT();

    if (interface().valid())
    {
        interface()->updateDisplacement(pointD_);
        interface()->updateDisplacementGradient(gradD_, gradDf_);
    }
    else
    {
        volToPoint_.interpolate(D_, pointD_);
        gradD_ = fvc::grad(D_, pointD_);
        gradDf_ = fvc::fGrad(D_, pointD_);
    }

//     D_.correctBoundaryConditions();
}


void unsTotalLagrangianSolid::updateFields()
{
    if (interface().valid())
    {
        interface()->updateDisplacement(pointD_);
        interface()->updateDisplacementGradient(gradD_, gradDf_);
    }
    else
    {
        volToPoint_.interpolate(D_, pointD_);
        gradD_ = fvc::grad(D_, pointD_);
        gradDf_ = fvc::fGrad(D_, pointD_);
    }

    rho_ = rheology_.rho();
    mu_ = rheology_.mu();
    lambda_ = rheology_.lambda();

    muf_ = surfaceScalarField("muf", fvc::interpolate(mu_));
    lambdaf_ = surfaceScalarField("lambdaf", fvc::interpolate(lambda_));

    if (interface().valid())
    {
        muf_ = interface()->interpolate(mu_);
        lambdaf_ = interface()->interpolate(lambda_);

//         interface()->modifyProperty(muf_);
//         interface()->modifyProperty(lambdaf_);
    }

    nonLinearGeometry::nonLinearType nonLinear =
        nonLinearGeometry::nonLinearNames_.read
        (
            solidProperties().lookup("nonLinear")
        );
//     Switch nonLinear(solidProperties().lookup("nonLinear"));

    // Calculate second Piola-Kirchhoff stress
    {
        epsilon_ = symm(gradD_);
        surfaceSymmTensorField epsilonf = symm(gradDf_);

        if(nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            epsilon_ += 0.5*symm(gradD_ & gradD_.T());
            epsilonf += 0.5*symm(gradDf_ & gradDf_.T());
        }

        sigma_ = 2*mu_*epsilon_ + I*(lambda_*tr(epsilon_));
        sigmaf_ = 2*muf_*epsilonf + I*(lambdaf_*tr(epsilonf));

        if (thermalStress())
        {
            sigma_ -= symmTensor(1,0,0,1,0,1)*threeK()*alpha()*DT();
            sigmaf_ -= symmTensor(1,0,0,1,0,1)*threeKf()*alphaf()*DTf();
        }
    }
}


tmp<surfaceVectorField> unsTotalLagrangianSolid::traction() const
{
    tmp<surfaceVectorField> tTraction
    (
        new surfaceVectorField
        (
            (mesh().Sf() & sigmaf_)/mesh().magSf()
        )
    );

    if (interface().valid())
    {
        interface()->correct(tTraction());
    }

    return tTraction;
}


bool unsTotalLagrangianSolid::writeObject
(
    IOstream::streamFormat,
    IOstream::versionNumber,
    IOstream::compressionType
) const
{
    // Switch moveMesh(solidProperties().lookup("moveMesh"));
    Switch moveMesh(false);

    nonLinearGeometry::nonLinearType nonLinear =
        nonLinearGeometry::nonLinearNames_.read
        (
            solidProperties().lookup("nonLinear")
        );
//     Switch nonLinear(solidProperties().lookup("nonLinear"));

    if (moveMesh)
    {
        pointIOField curPoints
        (
            IOobject
            (
                "points",
                runTime().timeName(),
                polyMesh::meshSubDir,
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh().allPoints()
        );

        const vectorField& pointDI = pointD_.internalField();

        forAll (pointDI, pointI)
        {
            curPoints[pointI] += pointDI[pointI];
        }

        // Unused points (procedure developed by Philip Cardiff, UCD)
        forAll(globalFaceZones(), zoneI)
        {
            const label curZoneID = globalFaceZones()[zoneI];

            const labelList& curMap =
                globalToLocalFaceZonePointMap()[zoneI];

            const labelList& curZoneMeshPoints =
                mesh().faceZones()[curZoneID]().meshPoints();

            vectorField curGlobalZonePointDispl
            (
                curZoneMeshPoints.size(),
                vector::zero
            );

            //-Inter-proc points are shared by multiple procs
            // pointNumProc is the number of procs which a point lies on
            scalarField pointNumProcs(curZoneMeshPoints.size(), 0);

            forAll(curGlobalZonePointDispl, globalPointI)
            {
                label localPoint = curMap[globalPointI];

                if(curZoneMeshPoints[localPoint] < mesh().nPoints())
                {
                    label procPoint = curZoneMeshPoints[localPoint];

                    curGlobalZonePointDispl[globalPointI] = pointDI[procPoint];

                    pointNumProcs[globalPointI] = 1;
                }
            }

            if (Pstream::parRun())
            {
                reduce(curGlobalZonePointDispl, sumOp<vectorField>());
                reduce(pointNumProcs, sumOp<scalarField>());

                //- now average the displacement between all procs
                curGlobalZonePointDispl /= pointNumProcs;
            }

            //- The curZonePointsDisplGlobal now contains the correct
            //  face zone displacement in a global master processor order,
            //  now convert them back into the local proc order

            vectorField curZonePointDispl
            (
                curZoneMeshPoints.size(),
                vector::zero
            );

            forAll(curGlobalZonePointDispl, globalPointI)
            {
                label localPoint = curMap[globalPointI];

                curZonePointDispl[localPoint] =
                    curGlobalZonePointDispl[globalPointI];
            }

            forAll(curZonePointDispl, pointI)
            {
                // unused points
                if (curZoneMeshPoints[pointI] >= mesh().nPoints())
                {
                    curPoints[curZoneMeshPoints[pointI]] +=
                        curZonePointDispl[pointI];
                }
            }
        }

        twoDPointCorrector twoDCorrector(mesh());
        twoDCorrector.correctPoints(curPoints);

        curPoints.write();

        // meshPhi must be present in order to reconstruction procedure works
        surfaceScalarField meshPhi
        (
            IOobject
            (
                "meshPhi",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("0", dimVolume/dimTime, 0.0)
        );
        meshPhi.write();
    }

    // Calculate equivalent stress
    if (nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        volTensorField F = I + gradD_.T();
        volScalarField J = det(F);
        volSymmTensorField sigmaCauchy
        (
            "sigmaCauchy",
            (1/J) * symm(F & sigma_ & F.T())
        );
        sigmaCauchy.write();

        volScalarField sigmaCauchyEq
        (
            IOobject
            (
                "sigmaCauchyEq",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqrt((3.0/2.0)*magSqr(dev(sigmaCauchy)))
        );
        sigmaCauchyEq.write();

        Info<< "Max sigmaCauchyEq = " << max(sigmaCauchyEq).value()
            << endl;

        Info<< "SigmaCauchyEq, max: " << gMax(sigmaCauchyEq.internalField())
            << ", avg: " << gAverage(sigmaCauchyEq.internalField())
            << ", min: " << gMin(sigmaCauchyEq.internalField()) << endl;
    }

    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma_)))
    );
    sigmaEq.write();

    Info<< "Max sigmaEq = " << max(sigmaEq).value()
        << endl;

    Info<< "SigmaEq, max: " << gMax(sigmaEq.internalField())
        << ", avg: " << gAverage(sigmaEq.internalField())
        << ", min: " << gMin(sigmaEq.internalField()) << endl;

    // Write point sigma field
    if (false)
    {
        pointSymmTensorField pointSigma
        (
            IOobject
            (
                "pointSigam",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh_,
            dimensioned<symmTensor>("0", sigma_.dimensions(), symmTensor::zero)
        );

        for (direction cmpt = 0; cmpt < symmTensor::nComponents; cmpt++)
        {
            volScalarField cmptSigma = sigma_.component(cmpt);

            pointScalarField cmptPointSigma
            (
                IOobject
                (
                    "cmptPointSigma",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ
                ),
                pMesh_,
                dimensioned<scalar>
                (
                    "0",
                    sigma_.dimensions(),
                    0
                )
            );

            volToPoint_.interpolate(cmptSigma, cmptPointSigma);

            pointSigma.internalField().replace
            (
                cmpt,
                cmptPointSigma.internalField()
            );
        }

//         pointTensorField pointEigenVector
//         (
//             IOobject
//             (
//                 "pointEigenVector",
//                 runTime().timeName(),
//                 mesh(),
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             pMesh_,
//             dimensioned<tensor>("0", dimless, tensor::zero)
//         );

//         pointVectorField pointEigenValue
//         (
//             IOobject
//             (
//                 "pointEigenValue",
//                 runTime().timeName(),
//                 mesh(),
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             pMesh_,
//             dimensioned<vector>("0", sigma_.dimensions(), vector::zero)
//         );

//         eig3Field(pointSigma, pointEigenVector, pointEigenValue);

//         pointVectorField pointEigenVector0
//         (
//             IOobject
//             (
//                 "pointEigenVector0",
//                 runTime().timeName(),
//                 mesh(),
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             pMesh_,
//             dimensioned<vector>("0", dimless, vector::zero)
//         );

//         pointEigenVector0.internalField().replace
//         (
//             0, //vector::X,
//             pointEigenVector.internalField().component(0)
//         );
//         pointEigenVector0.internalField().replace
//         (
//             1, //vector::Y,
//             pointEigenVector.internalField().component(1)
//         );
//         pointEigenVector0.internalField().replace
//         (
//             2, //vector::Z,
//             pointEigenVector.internalField().component(2)
//         );

//         pointVectorField pointEigenVector1
//         (
//             IOobject
//             (
//                 "pointEigenVector1",
//                 runTime().timeName(),
//                 mesh(),
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             pMesh_,
//             dimensioned<vector>("0", dimless, vector::zero)
//         );

//         pointEigenVector1.internalField().replace
//         (
//             0, //vector::X,
//             pointEigenVector.internalField().component(3)
//         );
//         pointEigenVector1.internalField().replace
//         (
//             1, //vector::Y,
//             pointEigenVector.internalField().component(4)
//         );
//         pointEigenVector1.internalField().replace
//         (
//             2, //vector::Z,
//             pointEigenVector.internalField().component(5)
//         );

//         pointVectorField pointEigenVector2
//         (
//             IOobject
//             (
//                 "pointEigenVector2",
//                 runTime().timeName(),
//                 mesh(),
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             pMesh_,
//             dimensioned<vector>("0", dimless, vector::zero)
//         );

//         pointEigenVector2.internalField().replace
//         (
//             0, //vector::X,
//             pointEigenVector.internalField().component(6)
//         );
//         pointEigenVector2.internalField().replace
//         (
//             1, //vector::Y,
//             pointEigenVector.internalField().component(7)
//         );
//         pointEigenVector2.internalField().replace
//         (
//             2, //vector::Z,
//             pointEigenVector.internalField().component(8)
//         );

        pointSigma.write();
//         pointEigenValue.write();
//         pointEigenVector0.write();
//         pointEigenVector1.write();
//         pointEigenVector2.write();
    }

    // Write cell and point hydrostatic pressure field
    if (false)
    {
        volSymmTensorField E = symm(gradD_);

        if(nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            E += 0.5*symm(gradD_ & gradD_.T());
        }

        volScalarField p("p", (lambda_+(2.0/3.0)*mu_)*tr(E));

        if (thermalStress())
        {
            p -= threeK()*alpha()*DT();
        }

        p.write();
    }

    return true;
}


const surfaceScalarField& unsTotalLagrangianSolid::DTf() const
{
//     if (!DTfPtr_)
//     {
//         DTfPtr_ =
//             const_cast<surfaceScalarField*>
//             (
//                 &mesh().lookupObject<surfaceScalarField>("DTf")
//             );
//     }

//     return *DTfPtr_;

    return mesh().lookupObject<surfaceScalarField>("DTf");
}


const volScalarField& unsTotalLagrangianSolid::DT() const
{
//     if (!DTPtr_)
//     {
//         DTPtr_ =
//             const_cast<volScalarField*>
//             (
//                 &mesh().lookupObject<volScalarField>("DT")
//             );
//     }

//     return *DTPtr_;

    return mesh().lookupObject<volScalarField>("DT");
}

const volScalarField& unsTotalLagrangianSolid::threeK() const
{
    if (!threeKPtr_)
    {
        threeKPtr_ =
            new volScalarField
            (
                "threeK",
                rheology_.threeK()*rheology_.rho()
            );
    }

    return *threeKPtr_;
}

const volScalarField& unsTotalLagrangianSolid::alpha() const
{
    if (!alphaPtr_)
    {
        const thermalModel& thermal =
            mesh().objectRegistry::lookupObject<thermalModel>
            (
                "thermalProperties"
            );

        alphaPtr_ = new volScalarField(thermal.alpha());
    }

    return *alphaPtr_;
}

const surfaceScalarField& unsTotalLagrangianSolid::threeKf() const
{
    if (!threeKfPtr_)
    {
        volScalarField threeK
        (
            "threeK",
            rheology_.threeK()*rheology_.rho()
        );

        if (interface().valid())
        {
            threeKfPtr_ =
                new surfaceScalarField
                (
                    "threeKf",
                    interface()->interpolate(threeK)
                );
        }
        else
        {
            threeKfPtr_ =
                new surfaceScalarField
                (
                    "threeKf",
                    fvc::interpolate(threeK)
                );
        }
    }

    return *threeKfPtr_;
}

const surfaceScalarField& unsTotalLagrangianSolid::alphaf() const
{
    if (!alphafPtr_)
    {
        const thermalModel& thermal =
            mesh().objectRegistry::lookupObject<thermalModel>
            (
                "thermalProperties"
            );

        if (interface().valid())
        {
            volScalarField alpha = thermal.alpha();

            alphafPtr_ =
                new surfaceScalarField
                (
                    "alphaf",
                    interface()->interpolate(alpha)
                );
        }
        else
        {

            alphafPtr_ =
                new surfaceScalarField
                (
                    "alphaf",
                    fvc::interpolate(thermal.alpha())
                );
        }
    }

    return *alphafPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidSolvers
} // End namespace Foam

// ************************************************************************* //

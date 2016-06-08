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

#include "unsIncrTotalLagrangianSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGradf.H"
#include "tractionDisplacementIncrementFvPatchVectorField.H"
#include "skewCorrectionVectors.H"
#include "multiMaterial.H"
#include "twoDPointCorrector.H"

// #include "componentReferenceList.H"
#include "nonLinearGeometry.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidSolvers
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(unsIncrTotalLagrangianSolid, 0);
addToRunTimeSelectionTable
(
    solidSolver,
    unsIncrTotalLagrangianSolid,
    dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unsIncrTotalLagrangianSolid::unsIncrTotalLagrangianSolid(const fvMesh& mesh)
:
    solidSolver(typeName, mesh),
    DD_
    (
        IOobject
        (
            "DD",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    D_
    (
        IOobject
        (
            "D",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength, vector::zero)
    ),
//     DU_
//     (
//         IOobject
//         (
//             "DU",
//             runTime().timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::NO_WRITE
//         ),
//         mesh,
//         dimensionedVector("0", dimVelocity, vector::zero)
//     ),
    U_
    (
        IOobject
        (
            "U",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimVelocity, vector::zero)
    ),
    pMesh_(mesh),
    pointDD_
    (
        IOobject
        (
            "pointDD",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
//             IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    pointD_
    (
        IOobject
        (
            "pointD",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    DSigma_
    (
        IOobject
        (
            "DSigma",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
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
    DEpsilon_
    (
        IOobject
        (
            "DEpsilon",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonf_
    (
        IOobject
        (
            "DEpsilonf",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
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
    epsilonf_
    (
        IOobject
        (
            "epsilonf",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
//     epsilonP_
//     (
//         IOobject
//         (
//             "epsilonP",
//             runTime().timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::NO_WRITE
//         ),
//         mesh,
//         dimensionedSymmTensor("zero", dimless, symmTensor::zero)
//     ),
    epsilonPf_
    (
        IOobject
        (
            "epsilonPf",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DSigmaf_
    (
        IOobject
        (
            "DSigmaf",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
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
    rheology_(DSigma_, DD_),
    volToPoint_(mesh),
    gradDDf_
    (
        IOobject
        (
            "gradDDf",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
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
    gradDD_
    (
        IOobject
        (
            "grad(" + DD_.name() + ")",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    gradD_
    (
        IOobject
        (
            "grad(" + D_.name() + ")",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    rho_(rheology_.rho()),
    mu_(rheology_.mu()),
    muf_("muf", fvc::interpolate(mu_)),
    lambda_(rheology_.lambda()),
    lambdaf_("lambdaf", fvc::interpolate(lambda_)),
    interface_(NULL),
    curTimeIndex_(runTime().timeIndex())
{
    pointDD_.oldTime();
//     DD_.oldTime();
//     D_.oldTime();

    if (rheology_.law().type() == multiMaterial::typeName)
    {
        interface_.set(new ITLMaterialInterface(DD_, pointDD_, D_));
    }

    if (interface().valid())
    {
        muf_ = interface()->interpolate(mu_);
        lambdaf_ = interface()->interpolate(lambda_);
//         interface()->modifyProperty(muf_);
//         interface()->modifyProperty(lambdaf_);
    }

//     bool enforceLinear = false;
//     solidProperties().set("enforceLinear", enforceLinear);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector unsIncrTotalLagrangianSolid::pointU(label pointID) const
{
//     volVectorField totU = U_ + DU_;

    return volToPoint_.interpolate(pointID, U_);
}

//- Patch point displacement
tmp<vectorField> unsIncrTotalLagrangianSolid::patchPointDisplacementIncrement
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
            pointDD_.internalField(),
            mesh().boundaryMesh()[patchID].meshPoints()
        );

    return tPointDisplacement;
}

//- Face zone point displacement
tmp<vectorField> unsIncrTotalLagrangianSolid
::faceZonePointDisplacementIncrement
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

    const vectorField& pointDDI = pointDD_.internalField();

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

                zonePointsDisplGlobal[globalPointI] = pointDDI[procPoint];

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
                pointDDI,
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    return tPointDisplacement;
}


//- Patch point displacement
tmp<vectorField> unsIncrTotalLagrangianSolid::patchPointDisplacement
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
            pointD_.internalField(),
            mesh().boundaryMesh()[patchID].meshPoints()
        );

    return tPointDisplacement;
}

//- Face zone point displacement
tmp<vectorField> unsIncrTotalLagrangianSolid
::faceZonePointDisplacement
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

                zonePointsDisplGlobal[globalPointI] = pointDI[procPoint];

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
                pointDI,
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    return tPointDisplacement;
}


//- Face zone point displacement
tmp<tensorField> unsIncrTotalLagrangianSolid
::faceZoneSurfaceGradientOfVelocity
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
unsIncrTotalLagrangianSolid::currentFaceZonePoints(const label zoneID) const
{
    vectorField pointDisplacement
    (
        mesh().faceZones()[zoneID]().localPoints().size(),
        vector::zero
    );

    const vectorField& pointDI = pointD_.internalField();
    const vectorField& pointDDI = pointDD_.internalField();

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
                    pointDI[procPoint] + pointDDI[procPoint];

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
                pointDI + pointDDI,
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
tmp<vectorField> unsIncrTotalLagrangianSolid::faceZoneNormal
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
    localPoints +=
        pointD_.boundaryField()[patchID].patchInternalField()
      + pointDD_.boundaryField()[patchID].patchInternalField();

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

void unsIncrTotalLagrangianSolid::setTraction
(
    const label patchID,
    const vectorField& traction
)
{
    if
    (
        DD_.boundaryField()[patchID].type()
     != tractionDisplacementIncrementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void unsIncrTotalLagrangianSolid::setTraction(...)")
            << "Bounary condition on " << DD_.name()
                <<  " is "
                << DD_.boundaryField()[patchID].type()
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead "
                << tractionDisplacementIncrementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementIncrementFvPatchVectorField& patchDD =
        refCast<tractionDisplacementIncrementFvPatchVectorField>
        (
            DD_.boundaryField()[patchID]
        );

    patchDD.traction() = traction;
}

void unsIncrTotalLagrangianSolid::setPressure
(
    const label patchID,
    const scalarField& pressure
)
{
    if
    (
        DD_.boundaryField()[patchID].type()
     != tractionDisplacementIncrementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void unsIncrTotalLagrangianSolid::setTraction(...)")
            << "Bounary condition on " << DD_.name()
                <<  " is "
                << DD_.boundaryField()[patchID].type()
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead "
                << tractionDisplacementIncrementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementIncrementFvPatchVectorField& patchDD =
        refCast<tractionDisplacementIncrementFvPatchVectorField>
        (
            DD_.boundaryField()[patchID]
        );

    patchDD.pressure() = pressure;
}

void unsIncrTotalLagrangianSolid::setTraction
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

void unsIncrTotalLagrangianSolid::setPressure
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


tmp<vectorField> unsIncrTotalLagrangianSolid::predictTraction
(
    const label patchID,
    const label zoneID
)
{
    // Predict traction on patch
    if
    (
        DD_.boundaryField()[patchID].type()
     != tractionDisplacementIncrementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void unsIncrTotalLagrangianSolid::predictTraction(...)")
            << "Bounary condition on " << DD_.name()
                <<  " is "
                << DD_.boundaryField()[patchID].type()
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead "
                << tractionDisplacementIncrementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementIncrementFvPatchVectorField& patchDUo =
        refCast<tractionDisplacementIncrementFvPatchVectorField>
        (
            DD_.oldTime().boundaryField()[patchID]
        );

    tractionDisplacementIncrementFvPatchVectorField& patchDUoo =
        refCast<tractionDisplacementIncrementFvPatchVectorField>
        (
            DD_.oldTime().oldTime().boundaryField()[patchID]
        );

    vectorField ptF = 2*patchDUo.traction() - patchDUoo.traction();

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


tmp<scalarField> unsIncrTotalLagrangianSolid::predictPressure
(
    const label patchID,
    const label zoneID
)
{
    // Predict pressure field on patch
    if
    (
        DD_.boundaryField()[patchID].type()
     != tractionDisplacementIncrementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn("void unsIncrTotalLagrangianSolid::predictTraction(...)")
            << "Bounary condition on " << DD_.name()
                <<  " is "
                << DD_.boundaryField()[patchID].type()
                << "for patch" << mesh().boundary()[patchID].name()
                << ", instead "
                << tractionDisplacementIncrementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementIncrementFvPatchVectorField& patchDUo =
        refCast<tractionDisplacementIncrementFvPatchVectorField>
        (
            DD_.oldTime().boundaryField()[patchID]
        );

    tractionDisplacementIncrementFvPatchVectorField& patchDUoo =
        refCast<tractionDisplacementIncrementFvPatchVectorField>
        (
            DD_.oldTime().oldTime().boundaryField()[patchID]
        );


    scalarField pPF = 2*patchDUo.pressure() - patchDUoo.pressure();

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


bool unsIncrTotalLagrangianSolid::evolve()
{
    Info << "Evolving solid solver: "
        << unsIncrTotalLagrangianSolid::typeName << endl;

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
    // Switch nonLinear(solidProperties().lookup("nonLinear"));


    Switch debug(solidProperties().lookup("debug"));

    // componentReferenceList cr
    // (
    //     solidProperties().lookup("componentReference"),
    //     componentReference::iNew(mesh())
    // );

    dimensionedScalar K("K", dimless/dimTime, 0);
    if (solidProperties().found("K"))
    {
        K = dimensionedScalar(solidProperties().lookup("K"));
    }

    int iCorr = 0;
    scalar initialResidual = 0;
    lduMatrix::solverPerformance solverPerf;
    scalar res = 1;
    scalar maxRes = 0;
    scalar curConvergenceTolerance = convergenceTolerance;

    lduMatrix::debug = debug;

//     bool enforceLinear = false;
//     solidProperties().set("enforceLinear", enforceLinear);

    // Activate old time fields
    fvc::ddt(D_);

    surfaceVectorField n = mesh().Sf()/mesh().magSf();

    do
    {
        if (lduMatrix::debug)
        {
            Info<< "Time: " << runTime().timeName()
                << ", outer iteration: " << iCorr << endl;
        }

        DD_.storePrevIter();
        DSigmaf_.storePrevIter();

        fvVectorMatrix DDEqn
        (
            rho_*fvm::d2dt2(DD_)
          - fvm::laplacian(2*muf_ + lambdaf_, DD_, "laplacian(DDD,DD)")
         == fvc::div
            (
                mesh().Sf()
              & (
                  - (muf_ + lambdaf_)*gradDDf_
                  + muf_*gradDDf_.T() + lambdaf_*(I*tr(gradDDf_))
                )
            )
        );

        // Add damping if not zero
        if (K.value() > SMALL)
        {
            DDEqn += K*rho_*fvm::ddt(DD_);
        }

        // Update strain increment
        DEpsilonf_ = symm(gradDDf_);
//         if (nonLinear && !enforceLinear)
        if (nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            DEpsilonf_ += 0.5*symm(gradDDf_ & gradDDf_.T());
            DEpsilonf_ += 0.5*symm(gradDDf_ & gradDf_.T());
            DEpsilonf_ += 0.5*symm(gradDf_ & gradDDf_.T());
        }

        DSigmaf_ = 2*muf_*DEpsilonf_ + I*(lambdaf_*tr(DEpsilonf_));

        if (rheology_.plasticityActive())
        {
            DSigmaf_ -= 2*muf_*fvc::interpolate(rheology_.DEpsilonP());
//             DSigmaf_ -= 2*muf_*rheology_.DEpsilonPf();
        }

        if (nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN)
//         if (nonLinear && !enforceLinear)
        {
            DDEqn -=
                fvc::div
                (
                    muf_*(mesh().Sf() & (gradDDf_ & gradDf_.T()))
                  + muf_*(mesh().Sf() & (gradDf_ & gradDDf_.T()))
                  + muf_*(mesh().Sf() & (gradDDf_ & gradDDf_.T()))
                  + 0.5*lambdaf_*tr(gradDDf_ & gradDf_.T())*mesh().Sf()
                  + 0.5*lambdaf_*tr(gradDf_ & gradDDf_.T())*mesh().Sf()
                  + 0.5*lambdaf_*tr(gradDDf_ & gradDDf_.T())*mesh().Sf()
                )
              + fvc::div(mesh().Sf() & (DSigmaf_ & gradDf_))
              + fvc::div(mesh().Sf() & ((sigmaf_ + DSigmaf_) & gradDDf_));
        }

        if (rheology_.plasticityActive())
        {
            DDEqn +=
                fvc::div
                (
                    2*muf_
                   *(
                        mesh().Sf()
                      & fvc::interpolate(rheology_.DEpsilonP())
                    )
                );

//             DDEqn +=
//                 fvc::div(2*muf_*(mesh().Sf() & rheology_.DEpsilonPf()));
        }

        if (interface().valid())
        {
            interface()->correct(DDEqn);
        }

        // forAll (cr, crI)
        // {
        //     DDEqn.setComponentReference
        //     (
        //         cr[crI].patchIndex(),
        //         cr[crI].faceIndex(),
        //         cr[crI].dir(),
        //         cr[crI].value()
        //     );
        // }

        solverPerf = DDEqn.solve();

        if(iCorr == 0)
        {
            initialResidual = solverPerf.initialResidual();
        }

        DD_.relax();

        if (interface().valid())
        {
            interface()->updateDisplacement(pointDD_);
            interface()->updateDisplacementGradient
            (
                gradDD_,
                gradDDf_
            );
        }
        else
        {
            volToPoint_.interpolate(DD_, pointDD_);
            gradDD_ = fvc::grad(DD_, pointDD_);
            gradDDf_ = fvc::fGrad(DD_, pointDD_);
        }


        // Calculate momentu residual
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
        }

        // Calculate strain increment
        {
            DEpsilon_ = symm(gradDD_);

            if(nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN)
//             if(nonLinear && !enforceLinear)
            {
                DEpsilon_ += 0.5*symm(gradDD_ & gradDD_.T());
                DEpsilon_ += 0.5*symm(gradDD_ & gradD_.T());
                DEpsilon_ += 0.5*symm(gradD_ & gradDD_.T());
            }
        }

        // Correct plasticity term
        rheology_.correct();

        if (iCorr==0)
        {
            res = 1;
        }
    }
    while
    (
        (res > curConvergenceTolerance)
     && (++iCorr < nCorr)
    );

//     fvc::ddt(D_);

    // Calculate second Piola-Kirchhoff stress increment
    {
        DSigma_ = 2*mu_*DEpsilon_ + I*(lambda_*tr(DEpsilon_));

        if (rheology_.plasticityActive())
        {
            DSigma_ -= 2*mu_*rheology_.DEpsilonP();
        }
    }

    Info << solverPerf.solverName() << ": Solving for " << DD_.name()
        << ", Initial residula = " << initialResidual
        << ", Final residual = " << solverPerf.initialResidual()
        << ", No outer iterations = " << iCorr
        << "\nMax relative residual = " << maxRes
        << ", Relative momentum residual = " << res << endl;
//         << ", enforceLinear = " << enforceLinear << endl;

    lduMatrix::debug = 1;

//     if (nonLinear && enforceLinear)
//     {
//         return false;
//     }

    return true;
}


void unsIncrTotalLagrangianSolid::predict()
{
    Info << "Predicting solid solver" << endl;

    DD_ = U_*runTime().deltaT();

    if (interface().valid())
    {
        interface()->updateDisplacement(pointDD_);
        interface()->updateDisplacementGradient(gradDD_, gradDDf_);
    }
    else
    {
        volToPoint_.interpolate(DD_, pointDD_);
        gradDD_ = fvc::grad(DD_, pointDD_);
        gradDDf_ = fvc::fGrad(DD_, pointDD_);
    }

    DD_.boundaryField().updateCoeffs();
}


tmp<surfaceVectorField> unsIncrTotalLagrangianSolid::traction() const
{
    tmp<surfaceVectorField> tTraction
    (
        new surfaceVectorField
        (
            (mesh().Sf() & (sigmaf_ + DSigmaf_))/mesh().magSf()
        )
    );

    if (interface().valid())
    {
        interface()->correct(tTraction());
    }

    return tTraction;
}


void unsIncrTotalLagrangianSolid::updateTotalFields()
{
    Info << "Update total fields" << endl;

    D_ += DD_;
    pointD_ += pointDD_;
    gradD_ += gradDD_;
    gradDf_ += gradDDf_;
    sigma_ += DSigma_;
    sigmaf_ += DSigmaf_;
    epsilon_ += DEpsilon_;
    epsilonf_ += DEpsilonf_;

    U_ = fvc::ddt(D_);

    if (rheology_.plasticityActive())
    {
//         epsilonP_ += rheology_.DEpsilonP();
//         epsilonPf_ += rheology_.DEpsilonPf();
        rheology_.updateYieldStress();
    }
}


bool unsIncrTotalLagrangianSolid::writeObject
(
    IOstream::streamFormat,
    IOstream::versionNumber,
    IOstream::compressionType
) const
{
    Switch moveMesh(false);
    // Switch moveMesh(solidProperties().lookup("moveMesh"));

    nonLinearGeometry::nonLinearType nonLinear =
        nonLinearGeometry::nonLinearNames_.read
        (
            solidProperties().lookup("nonLinear")
        );
   // Switch nonLinear(solidProperties().lookup("nonLinear"));

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

                    curGlobalZonePointDispl[globalPointI] =
                        pointDI[procPoint];

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


    // Calculate second Piola-Kirchhoff equivalent stress

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

    Info<< "Max sigmaEq = " << max(sigmaEq).value() << endl;

    Info<< "SigmaEq, max: " << gMax(sigmaEq.internalField())
        << ", avg: " << gAverage(sigmaEq.internalField())
        << ", min: " << gMin(sigmaEq.internalField()) << endl;

//     volTensorField F = I + gradD_;
//     volScalarField J = det(F);
//     J.write();


    // Calculate Cauchy stress (rotate stress field)
    if (nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        volTensorField F = I + gradD_;
        volScalarField J = det(F);
        J.write();
        volSymmTensorField sigmaCauchy
        (
            "sigmaCauchy",
            1/J * symm(F.T() & sigma_ & F)
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

        Info<< "Max sigmaCauchyEq = " << max(sigmaCauchyEq).value() << endl;

        Info<< "SigmaCauchyEq, max: " << gMax(sigmaCauchyEq.internalField())
            << ", avg: " << gAverage(sigmaCauchyEq.internalField())
            << ", min: " << gMin(sigmaCauchyEq.internalField()) << endl;
    }

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

        pointSigma.write();
    }


    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidSolvers
} // End namespace Foam

// ************************************************************************* //

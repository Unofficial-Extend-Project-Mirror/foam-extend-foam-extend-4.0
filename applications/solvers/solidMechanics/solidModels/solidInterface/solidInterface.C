/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

#include "solidInterface.H"
#include "fvc.H"
#include "processorFvPatchFields.H"
#include "fvMatrices.H"
#include "skewCorrectionVectors.H"
#include "leastSquaresGrad.H"
#include "gaussGrad.H"
#include "faceSet.H"
#include "faMesh.H"
#include "faCFD.H"
#include "zeroGradientFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidInterface, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void solidInterface::makeSubMesh() const
{
    if (debug)
    {
        Info<< "void solidInterface::makeSubMesh() const : "
            << "creating sum-mesh near interface"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (subMeshPtr_)
    {
        FatalErrorIn("solidInterface::makeSubMesh() const")
            << "sub-mesh already exist"
            << abort(FatalError);
    }

    const volScalarField& materials =
        mesh_.lookupObject<volScalarField>("materials");

    const scalarField& materialsI = materials.internalField();

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    labelHashSet interfaceCellSet;

    forAll(neighbour, faceI)
    {
        if
        (
            mag(materialsI[neighbour[faceI]] - materialsI[owner[faceI]])
          > SMALL
        )
        {
            interfaceCellSet.insert(neighbour[faceI]);
            interfaceCellSet.insert(owner[faceI]);
        }
    }

    labelList firstRowCells(interfaceCellSet.toc());

    const labelListList& cellCells = mesh_.cellCells();

    forAll(firstRowCells, cellI)
    {
        const labelList& curCells = cellCells[firstRowCells[cellI]];

        forAll(curCells, cellI)
        {
            if(!interfaceCellSet.found(curCells[cellI]))
            {
                interfaceCellSet.insert(curCells[cellI]);
            }
        }
    }

    subMeshPtr_ = new fvMeshSubset
    (
        IOobject
        (
            "interfaceSubMesh",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    subMeshPtr_->setLargeCellSubset(interfaceCellSet);
}


void solidInterface::makeGlobalInterFaces() const
{
    if (debug)
    {
        Info<< "void solidInterface::makeGlobalInterFaces() const : "
            << "creating global inter-faces addressing near interface"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (globalInterFacesPtr_)
    {
        FatalErrorIn("solidInterface::makeGlobalInterFaces() const")
            << "global inter-faces addressing already exist"
            << abort(FatalError);
    }

    if (mesh_.foundObject<volScalarField>("materials"))
    {
        const volScalarField& materials =
            mesh_.lookupObject<volScalarField>("materials");

        const scalarField& materialsI = materials.internalField();

        const unallocLabelList& owner = mesh_.owner();
        const unallocLabelList& neighbour = mesh_.neighbour();

        labelHashSet interFacesSet;

        forAll(neighbour, faceI)
        {
            if
            (
                mag(materialsI[neighbour[faceI]] - materialsI[owner[faceI]])
              > SMALL
            )
            {
                interFacesSet.insert(faceI);
            }
        }

        globalInterFacesPtr_ = new labelList(interFacesSet.toc());

        faceSet faMeshSet(mesh_, "faMeshSet", interFacesSet);
        faMeshSet.write();
    }
    else
    {
        globalInterFacesPtr_ = new labelList(0);
    }
}


void solidInterface::makeLocalInterFaces() const
{
    if (debug)
    {
        Info<< "void solidInterface::makeLocalInterFaces() const : "
            << "creating local inter-faces addressing near interface"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (localInterFacesPtr_)
    {
        FatalErrorIn("solidInterface::makeLocalInterFaces() const")
            << "local inter-faces addressing already exist"
            << abort(FatalError);
    }

    const volScalarField& materials =
        mesh_.lookupObject<volScalarField>("materials");

    volScalarField subMaterials = subMesh().interpolate(materials);

    const scalarField& materialsI = subMaterials.internalField();

    const fvMesh& sMesh = subMesh().subMesh();

    const unallocLabelList& owner = sMesh.owner();
    const unallocLabelList& neighbour = sMesh.neighbour();

    labelHashSet interFacesSet;

    forAll(neighbour, faceI)
    {
        if
        (
            mag(materialsI[neighbour[faceI]] - materialsI[owner[faceI]])
          > SMALL
        )
        {
            interFacesSet.insert(faceI);
        }
    }

    localInterFacesPtr_ = new labelList(interFacesSet.toc());
}


void solidInterface::makeInterfaceDisplacement() const
{
    if (debug)
    {
        Info<< "void solidInterface::makeInterfaceDisplacement() const : "
            << "creating interface displacement field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (interfaceUPtr_)
    {
        FatalErrorIn("solidInterface::makeInterfaceDisplacement() const")
            << "interface displacement field already exist"
            << abort(FatalError);
    }

    interfaceUPtr_ = new vectorField(globalInterFaces().size(), vector::zero);
}


void solidInterface::makeProcessorPatches() const
{
    if (debug)
    {
        Info<< "void solidInterface::makeProcessorPatches() const : "
            << "creating list of processor patches at the interface"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (processorPatchesPtr_)
    {
        FatalErrorIn("solidInterface::makeProcessorPatches() const")
            << "list of processor patches already exist"
            << abort(FatalError);
    }

    if (mesh_.foundObject<volScalarField>("materials"))
    {
        const volScalarField& materials =
            mesh_.lookupObject<volScalarField>("materials");

        labelHashSet processorPatches;

        forAll(materials.boundaryField(), patchI)
        {
            if (mesh_.boundary()[patchI].type() == processorFvPatch::typeName)
            {
                scalarField ownMat =
                    materials.boundaryField()[patchI].patchInternalField();

                scalarField ngbMat =
                    materials.boundaryField()[patchI].patchNeighbourField();

                forAll(ownMat, faceI)
                {
                    if (mag(ownMat[faceI] - ngbMat[faceI]) > SMALL)
                    {
                        processorPatches.insert(patchI);
                        break;
                    }
                }
            }
        }

        processorPatchesPtr_ = new labelList(processorPatches.toc());
    }
    else
    {
        processorPatchesPtr_ = new labelList(0);
    }
}


void solidInterface::makeProcessorPatchFaces() const
{
    if (debug)
    {
        Info<< "void solidInterface::makeProcessorPatchFaces() const : "
            << "creating list of processor patch faces at the interface"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (processorPatchFacesPtr_)
    {
        FatalErrorIn("solidInterface::makeProcessorPatchFaces() const")
            << "list of processor patch faces already exist"
            << abort(FatalError);
    }

    const labelList& procPatches = processorPatches();

    processorPatchFacesPtr_ = new labelListList(procPatches.size());
    labelListList& processorPatchFaces = *processorPatchFacesPtr_;

    forAll(procPatches, patchI)
    {
        const volScalarField& materials =
            mesh_.lookupObject<volScalarField>("materials");

        label curProcPatch = procPatches[patchI];

        scalarField ownMat =
            materials.boundaryField()[curProcPatch].patchInternalField();

        scalarField ngbMat =
            materials.boundaryField()[curProcPatch].patchNeighbourField();

        labelHashSet curProcPatchFaces;

        forAll(ownMat, faceI)
        {
            if (mag(ownMat[faceI] - ngbMat[faceI]) > SMALL)
            {
                curProcPatchFaces.insert(faceI);
            }
        }

        processorPatchFaces[patchI] = labelList(curProcPatchFaces.toc());
    }
}


void solidInterface::makeProcessorInterfaceDisplacement() const
{
    if (debug)
    {
        Info<< "void solidInterface::makeProcessorInterfaceDisplacement() const : "
            << "creating processor inter-faces displacement"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (processorInterfaceUPtr_)
    {
        FatalErrorIn("solidInterface::makeProcessorInterfaceDisplacement() const")
            << "processor interface displacement already exist"
            << abort(FatalError);
    }

    processorInterfaceUPtr_ =
        new FieldField<Field, vector>(processorPatches().size());
    FieldField<Field, vector>& processorInterfaceU = *processorInterfaceUPtr_;

    forAll(processorInterfaceU, patchI)
    {
        processorInterfaceU.set
        (
            patchI,
            new vectorField(processorPatchFaces()[patchI].size(), vector::zero)
        );
    }
}


void solidInterface::makeIndicator() const
{
    if (debug)
    {
        Info<< "void solidInterface::makeIndicator() const : "
            << "creating interface indicator"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (indicatorPtr_)
    {
        FatalErrorIn("solidInterface::makeIndicator() const")
            << "interface indicator already exist"
            << abort(FatalError);
    }

    indicatorPtr_ =
        new List<labelPair>(globalInterFaces().size(), labelPair(0, 0));

    List<labelPair>& indicator = *indicatorPtr_;

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    if (mesh_.foundObject<volScalarField>("materials"))
    {
        const volScalarField& materials =
            mesh_.lookupObject<volScalarField>("materials");

        forAll(globalInterFaces(), faceI)
        {
            label curFace = globalInterFaces()[faceI];

            labelPair& curPair = indicator[faceI];

            if (materials[owner[curFace]] < materials[neighbour[curFace]])
            {
                curPair.first() = int(materials[owner[curFace]]);
                curPair.second() = int(materials[neighbour[curFace]]);
            }
            else
            {
                curPair.second() = int(materials[owner[curFace]]);
                curPair.first() = int(materials[neighbour[curFace]]);
            }
        }
    }

//     Info << indicator << endl;
}


void solidInterface::clearOut()
{
    deleteDemandDrivenData(subMeshPtr_);
    deleteDemandDrivenData(globalInterFacesPtr_);
    deleteDemandDrivenData(localInterFacesPtr_);
    deleteDemandDrivenData(interfaceUPtr_);
    deleteDemandDrivenData(muPtr_);
    deleteDemandDrivenData(lambdaPtr_);
    deleteDemandDrivenData(processorPatchesPtr_);
    deleteDemandDrivenData(processorPatchFacesPtr_);
    deleteDemandDrivenData(processorInterfaceUPtr_);
    deleteDemandDrivenData(indicatorPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


solidInterface::solidInterface
(
    const fvMesh& mesh,
    const rheologyModel& rheology
)
:
    mesh_(mesh),
    rheology_(rheology),
    subMeshPtr_(NULL),
    globalInterFacesPtr_(NULL),
    localInterFacesPtr_(NULL),
    interfaceUPtr_(NULL),
    muPtr_(NULL),
    lambdaPtr_(NULL),
    processorPatchesPtr_(NULL),
    processorPatchFacesPtr_(NULL),
    processorInterfaceUPtr_(NULL),
    indicatorPtr_(NULL)
{
    muPtr_ = new volScalarField(rheology_.mu());
    lambdaPtr_ = new volScalarField(rheology_.lambda());
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

solidInterface::~solidInterface()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


const fvMeshSubset& solidInterface::subMesh() const
{
    if (!subMeshPtr_)
    {
        makeSubMesh();
    }

    return *subMeshPtr_;
}


const labelList& solidInterface::globalInterFaces() const
{
    if (!globalInterFacesPtr_)
    {
        makeGlobalInterFaces();
    }

    return *globalInterFacesPtr_;
}


const labelList& solidInterface::localInterFaces() const
{
    if (!localInterFacesPtr_)
    {
        makeLocalInterFaces();
    }

    return *localInterFacesPtr_;
}


vectorField& solidInterface::interfaceDisplacement()
{
    if (!interfaceUPtr_)
    {
        makeInterfaceDisplacement();
    }

    return *interfaceUPtr_;
}


const vectorField& solidInterface::interfaceDisplacement() const
{
    if (!interfaceUPtr_)
    {
        makeInterfaceDisplacement();
    }

    return *interfaceUPtr_;
}


const labelList& solidInterface::processorPatches() const
{
    if (!processorPatchesPtr_)
    {
        makeProcessorPatches();
    }

    return *processorPatchesPtr_;
}


const labelListList& solidInterface::processorPatchFaces() const
{
    if (!processorPatchFacesPtr_)
    {
        makeProcessorPatchFaces();
    }

    return *processorPatchFacesPtr_;
}


const FieldField<Field, vector>& solidInterface
::processorInterfaceDisplacement() const
{
    if (!processorInterfaceUPtr_)
    {
        makeProcessorInterfaceDisplacement();
    }

    return *processorInterfaceUPtr_;
}


FieldField<Field, vector>& solidInterface::processorInterfaceDisplacement()
{
    if (!processorInterfaceUPtr_)
    {
        makeProcessorInterfaceDisplacement();
    }

    return *processorInterfaceUPtr_;
}


void solidInterface::correct(fvVectorMatrix& UEqn)
{
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const volVectorField& U = UEqn.psi();
    const vectorField& UI = U.internalField();

    const volTensorField& gradU =
        mesh_.lookupObject<volTensorField>("grad(" + U.name() + ')');
    const tensorField& gradUI = gradU.internalField();

    const volScalarField& mu = *muPtr_;
    const volScalarField& lambda = *lambdaPtr_;

    const vectorField& SI  = mesh_.Sf().internalField();
    const scalarField& magSI  = mesh_.magSf().internalField();
    const scalarField& deltaCoeffs = mesh_.deltaCoeffs().internalField();
    const scalarField& w = mesh_.weights().internalField();

    scalarField& diag = UEqn.diag();
    scalarField& upper = UEqn.upper();
    vectorField& source = UEqn.source();
    FieldField<Field, vector>& boundaryCoeffs = UEqn.boundaryCoeffs();

    vectorField& interU = interfaceDisplacement();


//     // Calc surface gradient using FAM
//     tensorField interSGradU(interU.size(), tensor::zero);
//     {
//         faMesh aMesh(mesh_);

//         const labelList& faceLabels = aMesh.faceLabels();

//         wordList patchFieldTypes
//         (
//             aMesh.boundary().size(),
//             zeroGradientFaPatchVectorField::typeName
//         );

//         areaVectorField Us
//         (
//             IOobject
//             (
//                 "Us",
//                 mesh_.time().timeName(),
//                 mesh_,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             aMesh,
//             dimensioned<vector>("Us", dimLength, vector::zero),
//             patchFieldTypes
//         );

//         forAll(globalInterFaces(), faceI)
//         {
//             label curFace = globalInterFaces()[faceI];
//             label index = findIndex(faceLabels, curFace);

//             Us.internalField()[index] = interU[faceI];
//         }
//         Us.correctBoundaryConditions();

//         tensorField sGradU = fac::grad(Us)().internalField();
//         forAll(globalInterFaces(), faceI)
//         {
//             label curFace = globalInterFaces()[faceI];
//             label index = findIndex(faceLabels, curFace);

//             interSGradU[faceI] = sGradU[index];
//         }
//     }

    // Internal faces
    forAll(globalInterFaces(), faceI)
    {
        label curGlobalFace = globalInterFaces()[faceI];

        label curOwner = owner[curGlobalFace];
        label curNeighbour = neighbour[curGlobalFace];

        vector ownN = SI[curGlobalFace]/magSI[curGlobalFace];
        vector ngbN = -ownN;

        scalar magS = magSI[curGlobalFace];

        tensor ownSGradU = ((I-ownN*ownN)&gradUI[curOwner]);
        tensor ngbSGradU = ((I-ngbN*ngbN)&gradUI[curNeighbour]);

//         ownSGradU = interSGradU[faceI];
//         ngbSGradU = interSGradU[faceI];

        scalar ownTrSGradUt = tr(ownSGradU&(I-ownN*ownN));
        scalar ngbTrSGradUt = tr(ngbSGradU&(I-ownN*ownN));

//         ownTrSGradUt = 0.5*(ownTrSGradUt + ngbTrSGradUt);
//         ngbTrSGradUt = ownTrSGradUt;

        vector ownSGradUn = (ownSGradU&ownN);
        vector ngbSGradUn = (ngbSGradU&ownN);

//         ownSGradUn = 0.5*(ownSGradUn + ngbSGradUn);
//         ngbSGradUn = ownSGradUn;

        scalar ngbDn = w[curGlobalFace]*(1.0/deltaCoeffs[curGlobalFace]);
        scalar ownDn = (1.0/deltaCoeffs[curGlobalFace]) - ngbDn;

        scalar ownMu  = mu.internalField()[curOwner];
        scalar ngbMu  = mu.internalField()[curNeighbour];

        scalar ownLambda  = lambda.internalField()[curOwner];
        scalar ngbLambda  = lambda.internalField()[curNeighbour];

        vector ownUt = ((I-ownN*ownN)&UI[curOwner]);
        vector ngbUt = ((I-ngbN*ngbN)&UI[curNeighbour]);

        vector ownUn = ownN*(ownN&UI[curOwner]);
        vector ngbUn = ngbN*(ngbN&UI[curNeighbour]);


        // Interface displacement

        vector curInterUt =
            (
                ownMu*ownUt*ngbDn + ngbMu*ngbUt*ownDn
              + ownDn*ngbDn*(ngbMu*ngbSGradUn - ownMu*ownSGradUn)
            )
           /(ownMu*ngbDn + ngbMu*ownDn);

        vector curInterUn =
            (
                (2*ownMu + ownLambda)*ownUn*ngbDn
              + (2*ngbMu + ngbLambda)*ngbUn*ownDn
              + ownDn*ngbDn
               *(ngbLambda*ngbTrSGradUt - ownLambda*ownTrSGradUt)*ownN
            )
           /((2*ownMu + ownLambda)*ngbDn + (2*ngbMu + ngbLambda)*ownDn);

        interU[faceI] = curInterUn + curInterUt;


        // Implicit coupling

        scalar wRevLin = 1.0 - w[curGlobalFace];

        scalar ownK = (2*ownMu + ownLambda);
        scalar ngbK = (2*ngbMu + ngbLambda);

        scalar Kf = 1.0/(wRevLin/ownK + (1.0-wRevLin)/ngbK);
        scalar muf = 1.0/(wRevLin/ownMu + (1.0-wRevLin)/ngbMu);

        scalar Dnf = 1.0/deltaCoeffs[curGlobalFace];

        // Owner
        diag[curOwner] += Kf*magS/Dnf;

        upper[curGlobalFace] -= Kf*magS/Dnf;

        source[curOwner] +=
          - (Kf - muf)
           *((UI[curNeighbour] - UI[curOwner])&(I-ownN*ownN))*magS/Dnf;

        source[curOwner] +=
            (
                ownK*ngbDn*ngbLambda*ngbTrSGradUt
              + ngbK*ownDn*ownLambda*ownTrSGradUt
            )*ownN*magS
           /(ownK*ngbDn + ngbK*ownDn)
          + (
              ownMu*ngbMu*ngbDn*ngbSGradUn
            + ownMu*ngbMu*ownDn*ownSGradUn
            )*magS
           /(ownMu*ngbDn + ngbMu*ownDn);

        // Neighbour
        diag[curNeighbour] += Kf*magS/Dnf;

        source[curNeighbour] -=
          - (Kf - muf)
           *((UI[curNeighbour] - UI[curOwner])&(I-ownN*ownN))*magS/Dnf;

        source[curNeighbour] -=
            (
                ownK*ngbDn*ngbLambda*ngbTrSGradUt
              + ngbK*ownDn*ownLambda*ownTrSGradUt
            )*ownN*magS
           /(ownK*ngbDn + ngbK*ownDn)
          + (
              ownMu*ngbMu*ngbDn*ngbSGradUn
            + ownMu*ngbMu*ownDn*ownSGradUn
            )*magS
           /(ownMu*ngbDn + ngbMu*ownDn);
    }



    // Processor faces

    forAll(processorPatchFaces(), patchI)
    {
        label curPatch = processorPatches()[patchI];

        vectorField& curProcInterU =
            processorInterfaceDisplacement()[patchI];

        const vectorField curProcOwnU =
            U.boundaryField()[curPatch].patchInternalField();
        const vectorField curProcNgbU =
            U.boundaryField()[curPatch].patchNeighbourField();

        const tensorField curProcOwnGradU =
            gradU.boundaryField()[curPatch].patchInternalField();
        const tensorField curProcNgbGradU =
            gradU.boundaryField()[curPatch].patchNeighbourField();

        const vectorField& curProcS =
            mesh_.Sf().boundaryField()[curPatch];
        const scalarField& curProcMagS =
            mesh_.magSf().boundaryField()[curPatch];
        const scalarField& curProcDeltaCoeffs =
            mesh_.deltaCoeffs().boundaryField()[curPatch];
        const scalarField& curProcW =
            mesh_.weights().boundaryField()[curPatch];

        const scalarField curProcOwnMu =
            mu.boundaryField()[curPatch].patchInternalField();
        const scalarField curProcNgbMu =
            mu.boundaryField()[curPatch].patchNeighbourField();

        const scalarField curProcOwnLambda =
            lambda.boundaryField()[curPatch].patchInternalField();
        const scalarField curProcNgbLambda =
            lambda.boundaryField()[curPatch].patchNeighbourField();

        const unallocLabelList& curProcFaceCells =
            mesh_.boundary()[curPatch].faceCells();

        forAll(processorPatchFaces()[patchI], faceI)
        {
            label curFace = processorPatchFaces()[patchI][faceI];

            scalar ngbDn = curProcW[curFace]*(1.0/curProcDeltaCoeffs[curFace]);
            scalar ownDn = (1.0/curProcDeltaCoeffs[curFace]) - ngbDn;

            scalar magS = curProcMagS[curFace];

            scalar ownMu  = curProcOwnMu[curFace];
            scalar ngbMu  = curProcNgbMu[curFace];

            scalar ownLambda  = curProcOwnLambda[curFace];
            scalar ngbLambda  = curProcNgbLambda[curFace];

            vector ownN = curProcS[curFace]/curProcMagS[curFace];
            vector ngbN = -ownN;

            vector ownUt = ((I-ownN*ownN)&curProcOwnU[curFace]);
            vector ngbUt = ((I-ngbN*ngbN)&curProcNgbU[curFace]);

            vector ownUn = ownN*(ownN&curProcOwnU[curFace]);
            vector ngbUn = ngbN*(ngbN&curProcNgbU[curFace]);

            tensor ownSGradU = ((I-ownN*ownN)&curProcOwnGradU[curFace]);
            tensor ngbSGradU = ((I-ngbN*ngbN)&curProcNgbGradU[curFace]);

            scalar ownTrSGradUt = tr(ownSGradU&(I-ownN*ownN));
            scalar ngbTrSGradUt = tr(ngbSGradU&(I-ngbN*ngbN));

            vector ownSGradUn = (ownSGradU&ownN);
            vector ngbSGradUn = (ngbSGradU&ngbN);


            // Interface displacement

            vector curInterUt =
            (
                ownMu*ownUt/ownDn + ngbMu*ngbUt/ngbDn
              - (ngbMu*ngbSGradUn + ownMu*ownSGradUn)
            )
           /(ownMu/ownDn + ngbMu/ngbDn);

            vector curInterUn =
            (
                (2*ownMu + ownLambda)*ownUn/ownDn
              + (2*ngbMu + ngbLambda)*ngbUn/ngbDn
              + (ngbLambda*ngbTrSGradUt - ownLambda*ownTrSGradUt)*ownN
            )
           /((2*ownMu + ownLambda)/ownDn + (2*ngbMu + ngbLambda)/ngbDn);

            curProcInterU[faceI] = curInterUn + curInterUt;


            // Implicit coupling

            scalar wRevLin = 1.0 - curProcW[curFace];

            scalar ownK = (2*ownMu + ownLambda);
            scalar ngbK = (2*ngbMu + ngbLambda);

            scalar Kf = 1.0/(wRevLin/ownK + (1.0-wRevLin)/ngbK);
            scalar muf = 1.0/(wRevLin/ownMu + (1.0-wRevLin)/ngbMu);

            scalar Dnf = 1.0/curProcDeltaCoeffs[curFace];

            // Owner
            diag[curProcFaceCells[curFace]] += Kf*magS/Dnf;


            boundaryCoeffs[curPatch][curFace] += Kf*magS*vector::one/Dnf;
//                 upper[curGlobalFace] -= Kf*magS/Dnf;

            source[curProcFaceCells[curFace]] +=
              - (Kf - muf)*((ngbUt - ownUt)&(I-ownN*ownN))*magS/Dnf
              + ownLambda*ownTrSGradUt*ownN*magS
              + ownMu*ownSGradUn*magS;

            source[curProcFaceCells[curFace]] +=
                (ownK*ngbDn/(ownK*ngbDn + ngbK*ownDn))
               *(ngbLambda*ngbTrSGradUt - ownLambda*ownTrSGradUt)*ownN*magS
              - (ownMu*ngbDn/(ownMu*ngbDn + ngbMu*ownDn))
               *(ngbMu*ngbSGradUn + ownMu*ownSGradUn)*magS;
        }
    }
}


void solidInterface::modifyProperties
(
    surfaceScalarField& muf,
    surfaceScalarField& lambdaf
) const
{
    forAll(globalInterFaces(), faceI)
    {
        label curGlobalFace = globalInterFaces()[faceI];

        muf.internalField()[curGlobalFace] = 0;
        lambdaf.internalField()[curGlobalFace] = 0;
    }

    forAll(processorPatchFaces(), patchI)
    {
        label curPatch = processorPatches()[patchI];

        forAll(processorPatchFaces()[patchI], faceI)
        {
            label curFace = processorPatchFaces()[patchI][faceI];

            muf.boundaryField()[curPatch][curFace] = 0;
            lambdaf.boundaryField()[curPatch][curFace] = 0;
        }
    }
}


tmp<volTensorField> solidInterface::grad(volVectorField& U) const
{
    tmp<volTensorField> tGradU
    (
        new volTensorField
        (
            IOobject
            (
                "grad(" + U.name() + ')',
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor("zero", dimless, tensor::zero)
        )
    );
    volTensorField& gradU = tGradU();

    surfaceVectorField Uf = fvc::interpolate(U);


    // Skew-correction

    if (skewCorrectionVectors::New(this->mesh_).skew())
    {
        const skewCorrectionVectors& scv = skewCorrectionVectors::New(mesh_);

        const volTensorField& gradU =
            mesh_.lookupObject<volTensorField>("grad(" + U.name() + ')');

        Uf +=
        (
            scv()
          & linear<tensor>(mesh_).interpolate
            (
                gradU
            )
        );

//         Uf +=
//         (
//             scv()
//           & linear<tensor>(mesh_).interpolate
//             (
//                 fv::leastSquaresGrad<vector>(mesh_).grad(U)
//             )
//         );
    }


    // Interface correction

    const vectorField& interU = interfaceDisplacement();

    forAll(globalInterFaces(), faceI)
    {
        label curGlobalFace = globalInterFaces()[faceI];

        Uf.internalField()[curGlobalFace] = interU[faceI];
    }

    forAll(processorPatchFaces(), patchI)
    {
        label curPatch = processorPatches()[patchI];

        const vectorField& curProcInterU =
            processorInterfaceDisplacement()[patchI];

        forAll(processorPatchFaces()[patchI], faceI)
        {
            label curFace = processorPatchFaces()[patchI][faceI];

            Uf.boundaryField()[curPatch][curFace] = curProcInterU[faceI];
        }
    }


    // Gradient calculation using Gauss method

    gradU = fv::gaussGrad<vector>(mesh_).grad(Uf);
    fv::gaussGrad<vector>(mesh_).correctBoundaryConditions(U, gradU);


    return tGradU;
}


tmp<symmTensorField> solidInterface::sigmaA() const
{
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const volVectorField& U =
        mesh_.lookupObject<volVectorField>("U");
    const vectorField& UI = U.internalField();

    const volTensorField& gradU =
        mesh_.lookupObject<volTensorField>("grad(" + U.name() + ')');
    const tensorField& gradUI = gradU.internalField();

    const volScalarField& mu = *muPtr_;
    const volScalarField& lambda = *lambdaPtr_;

    const vectorField& interU = interfaceDisplacement();

    const vectorField& SI  = mesh_.Sf().internalField();
    const scalarField& magSI  = mesh_.magSf().internalField();
    const scalarField& deltaCoeffs = mesh_.deltaCoeffs().internalField();
    const scalarField& w = mesh_.weights().internalField();

    tmp<symmTensorField> tSigmaA
    (
        new symmTensorField(globalInterFaces().size(), symmTensor::zero)
    );
    symmTensorField& sigmaA = tSigmaA();

    if (mesh_.foundObject<volScalarField>("materials"))
    {
        const volScalarField& materials =
            mesh_.lookupObject<volScalarField>("materials");

        tensorField gradUA(sigmaA.size(), tensor::zero);
        scalarField muA(sigmaA.size(), 0);
        scalarField lambdaA(sigmaA.size(), 0);

//         // Calc surface gradient using FAM
//         tensorField interSGradU(interU.size(), tensor::zero);
//         {
//             faMesh aMesh(mesh_);

//             const labelList& faceLabels = aMesh.faceLabels();

//             wordList patchFieldTypes
//             (
//                 aMesh.boundary().size(),
//                 zeroGradientFaPatchVectorField::typeName
//             );

//             areaVectorField Us
//             (
//                 IOobject
//                 (
//                     "Us",
//                     mesh_.time().timeName(),
//                     mesh_,
//                     IOobject::NO_READ,
//                     IOobject::AUTO_WRITE
//                 ),
//                 aMesh,
//                 dimensioned<vector>("Us", dimLength, vector::zero),
//                 patchFieldTypes
//             );

//             forAll(globalInterFaces(), faceI)
//             {
//                 label curFace = globalInterFaces()[faceI];
//                 label index = findIndex(faceLabels, curFace);

//                 Us.internalField()[index] = interU[faceI];
//             }
//             Us.correctBoundaryConditions();

//             tensorField sGradU = fac::grad(Us)().internalField();
//             forAll(globalInterFaces(), faceI)
//             {
//                 label curFace = globalInterFaces()[faceI];
//                 label index = findIndex(faceLabels, curFace);

//                 interSGradU[faceI] = sGradU[index];
//             }
//         }

        forAll(globalInterFaces(), faceI)
        {
            label curFace = globalInterFaces()[faceI];

            scalar ngbDn = w[curFace]*(1.0/deltaCoeffs[curFace]);
            scalar ownDn = (1.0/deltaCoeffs[curFace]) - ngbDn;

            label curCell = owner[curFace];
            vector n = SI[curFace]/magSI[curFace];
            scalar dn = ownDn;

            if (materials[owner[curFace]] < materials[neighbour[curFace]])
            {
                curCell = neighbour[curFace];
                n *= -1;
                dn = ngbDn;
            }

            muA[faceI] = mu[curCell];
            lambdaA[faceI] = lambda[curCell];

            gradUA[faceI] = n*(interU[faceI] - UI[curCell])/dn;

            gradUA[faceI] += ((I-n*n)&gradUI[curCell]);
//             gradUA[faceI] += interSGradU[faceI];
        }

        sigmaA = 2*muA*symm(gradUA)  + lambdaA*(I*tr(gradUA));
    }

    return tSigmaA;
}


tmp<symmTensorField> solidInterface::sigmaB() const
{
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const volVectorField& U =
        mesh_.lookupObject<volVectorField>("U");
    const vectorField& UI = U.internalField();

    const volTensorField& gradU =
        mesh_.lookupObject<volTensorField>("grad(" + U.name() + ')');
    const tensorField& gradUI = gradU.internalField();

    const volScalarField& mu = *muPtr_;
    const volScalarField& lambda = *lambdaPtr_;

    const vectorField& interU = interfaceDisplacement();

    const vectorField& SI  = mesh_.Sf().internalField();
    const scalarField& magSI  = mesh_.magSf().internalField();
    const scalarField& deltaCoeffs = mesh_.deltaCoeffs().internalField();
    const scalarField& w = mesh_.weights().internalField();

    tmp<symmTensorField> tSigmaB
    (
        new symmTensorField(globalInterFaces().size(), symmTensor::zero)
    );
    symmTensorField& sigmaB = tSigmaB();

    if (mesh_.foundObject<volScalarField>("materials"))
    {
        const volScalarField& materials =
            mesh_.lookupObject<volScalarField>("materials");

        tensorField gradUB(sigmaB.size(), tensor::zero);
        scalarField muB(sigmaB.size(), 0);
        scalarField lambdaB(sigmaB.size(), 0);

//         // Calc surface gradient using FAM
//         tensorField interSGradU(interU.size(), tensor::zero);
//         {
//             faMesh aMesh(mesh_);

//             const labelList& faceLabels = aMesh.faceLabels();

//             wordList patchFieldTypes
//             (
//                 aMesh.boundary().size(),
//                 zeroGradientFaPatchVectorField::typeName
//             );

//             areaVectorField Us
//             (
//                 IOobject
//                 (
//                     "Us",
//                     mesh_.time().timeName(),
//                     mesh_,
//                     IOobject::NO_READ,
//                     IOobject::AUTO_WRITE
//                 ),
//                 aMesh,
//                 dimensioned<vector>("Us", dimLength, vector::zero),
//                 patchFieldTypes
//             );

//             forAll(globalInterFaces(), faceI)
//             {
//                 label curFace = globalInterFaces()[faceI];
//                 label index = findIndex(faceLabels, curFace);

//                 Us.internalField()[index] = interU[faceI];
//             }
//             Us.correctBoundaryConditions();

//             tensorField sGradU = fac::grad(Us)().internalField();
//             forAll(globalInterFaces(), faceI)
//             {
//                 label curFace = globalInterFaces()[faceI];
//                 label index = findIndex(faceLabels, curFace);

//                 interSGradU[faceI] = sGradU[index];
//             }
//         }

        forAll(globalInterFaces(), faceI)
        {
            label curFace = globalInterFaces()[faceI];

            scalar ngbDn = w[curFace]*(1.0/deltaCoeffs[curFace]);
            scalar ownDn = (1.0/deltaCoeffs[curFace]) - ngbDn;

            label curCell = neighbour[curFace];
            vector n = -SI[curFace]/magSI[curFace];
            scalar dn = ngbDn;

            if (materials[owner[curFace]] < materials[neighbour[curFace]])
            {
                curCell = owner[curFace];
                n *= -1;
                dn = ownDn;
            }

            muB[faceI] = mu[curCell];
            lambdaB[faceI] = lambda[curCell];

            gradUB[faceI] = n*(interU[faceI] - UI[curCell])/dn;

            gradUB[faceI] += ((I-n*n)&gradUI[curCell]);
//             gradUB[faceI] += interSGradU[faceI];
        }

        sigmaB = 2*muB*symm(gradUB)  + lambdaB*(I*tr(gradUB));
    }

    return tSigmaB;
}


const List<labelPair>& solidInterface::indicator() const
{
    if (!indicatorPtr_)
    {
        makeIndicator();
    }

    return *indicatorPtr_;
}


void solidInterface::correctGrad
(
    const volVectorField& U,
    volTensorField& gradU
) const
{
    const fvMesh& sMesh = subMesh().subMesh();

    const unallocLabelList& owner = sMesh.owner();
    const unallocLabelList& neighbour = sMesh.neighbour();

    volVectorField subU = subMesh().interpolate(U);

    surfaceVectorField Us = linearInterpolate(subU);

    const vectorField& interU = interfaceDisplacement();

    forAll(localInterFaces(), faceI)
    {
        label curFace = localInterFaces()[faceI];

        Us.internalField()[curFace] = interU[faceI];
    }

    volTensorField gaussGradU =
        fv::gaussGrad<vector>(sMesh).grad(Us);
    fv::gaussGrad<vector>(sMesh).correctBoundaryConditions
    (
        subU,
        gaussGradU
    );

    forAll(localInterFaces(), faceI)
    {
        label curFace = localInterFaces()[faceI];

        gradU.internalField()
            [subMesh().cellMap()[owner[curFace]]] =
            gaussGradU.internalField()[owner[curFace]];

        gradU.internalField()
            [subMesh().cellMap()[neighbour[curFace]]] =
            gaussGradU.internalField()[neighbour[curFace]];
    }

    fv::gaussGrad<vector>(mesh_).correctBoundaryConditions(U, gradU);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

Description

\*---------------------------------------------------------------------------*/

#include "ULLSMaterialInterface.H"
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
#include "OStringStream.H"
#include "IStringStream.H"
#include "solidSolver.H"
#include "directTopoChange.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ULLSMaterialInterface, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void ULLSMaterialInterface::makeDisplacementIncrement() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface::"
            << "makeInterfaceDisplacementIncrement() const : "
            << "creating interface displacement field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (displacementIncrementPtr_)
    {
        FatalErrorIn("ULLSMaterialInterface::makeDisplacementIncrement() const")
            << "interface displacement increment field already exist"
            << abort(FatalError);
    }

    displacementIncrementPtr_ =
        new vectorField(faces().size(), vector::zero);

    // Initialize displacement increment
    surfaceVectorField DDf = fvc::interpolate(DD_);

    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh().nInternalFaces())
        {
            (*displacementIncrementPtr_)[faceI] = DDf.internalField()[curFace];
        }
        else
        {
            label curPatch = mesh().boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh().boundaryMesh()[curPatch].start();

            (*displacementIncrementPtr_)[faceI] =
                DDf.boundaryField()[curPatch][curPatchFace];
        }
    }
}


void ULLSMaterialInterface::makeDisplacement() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface"
            << "::makeInterfaceDisplacement() const : "
            << "creating interface displacement field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (displacementPtr_)
    {
        FatalErrorIn("ULLSMaterialInterface::makeDisplacement() const")
            << "interface displacement field already exist"
            << abort(FatalError);
    }

    displacementPtr_ = new vectorField(faces().size(), vector::zero);

    // Initialize displacement
    surfaceVectorField Df = fvc::interpolate(D_);

    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh().nInternalFaces())
        {
            (*displacementPtr_)[faceI] = Df.internalField()[curFace];
        }
        else
        {
            label curPatch = mesh().boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh().boundaryMesh()[curPatch].start();

            (*displacementPtr_)[faceI] =
                Df.boundaryField()[curPatch][curPatchFace];
        }
    }
}


// void ULLSMaterialInterface::makeTractionIncrement() const
// {
//     if (debug)
//     {
//         Info<< "void ULLSMaterialInterface::makeTractionIncrement() const : "
//             << "creating interface traction increment field"
//             << endl;
//     }

//     // It is an error to attempt to recalculate
//     // if the pointer is already set
//     if (tractionIncrementPtr_)
//     {
//         FatalErrorIn("ULLSMaterialInterface::makeTractionIncrement() const")
//             << "interface traction increment field already exist"
//             << abort(FatalError);
//     }

//     tractionIncrementPtr_ =
//         new vectorField(faces().size(), vector::zero);
// }


void ULLSMaterialInterface::makeTraction() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface::makeTraction() const : "
            << "creating interface traction field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (tractionPtr_)
    {
        FatalErrorIn("ULLSMaterialInterface::makeTraction() const")
            << "interface traction field already exist"
            << abort(FatalError);
    }

    tractionPtr_ = new vectorField(faces().size(), vector::zero);
}


void ULLSMaterialInterface::makeRelF() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface::makeRelF() const : "
            << "creating interface relative deformation gradient field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ngbRelFPtr_ || ownRelFPtr_)
    {
        FatalErrorIn("ULLSMaterialInterface::makeRelF() const")
            << "interface relative deformation gradient field "
                << "already exists"
                << abort(FatalError);
    }

    ownRelFPtr_ = new tensorField(faces().size(), I);
    ngbRelFPtr_ = new tensorField(faces().size(), I);
}


void ULLSMaterialInterface::makeF() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface::makeF() const : "
            << "creating interface deformation gradient field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ngbFPtr_ || ownFPtr_)
    {
        FatalErrorIn("ULLSMaterialInterface::makeF() const")
            << "interface old deformation gradient field "
                << "already exists"
                << abort(FatalError);
    }

    ownFPtr_ = new tensorField(faces().size(), I);
    ngbFPtr_ = new tensorField(faces().size(), I);
}


void ULLSMaterialInterface::makeBbar() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface::makeF() const : "
            << "creating interface deformation gradient field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ngbBbarPtr_ || ownBbarPtr_)
    {
        FatalErrorIn("ULLSMaterialInterface::makeF() const")
            << "interface old deformation gradient field "
                << "already exists"
                << abort(FatalError);
    }

    ownBbarPtr_ = new symmTensorField(faces().size(), I);
    ngbBbarPtr_ = new symmTensorField(faces().size(), I);
}


void ULLSMaterialInterface::makeDEpsilonP() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface::makeDEpsilonP() const : "
            << "creating interface plastic strain increment field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ngbDEpsilonPPtr_ || ownDEpsilonPPtr_)
    {
        FatalErrorIn("ULLSMaterialInterface::makeDEpsilonP() const")
            << "interface plastic strain increment field "
                << "already exists"
                << abort(FatalError);
    }

    ownDEpsilonPPtr_ = new symmTensorField(faces().size(), symmTensor::zero);
    ngbDEpsilonPPtr_ = new symmTensorField(faces().size(), symmTensor::zero);
}


void ULLSMaterialInterface::makeSigmaY() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface::makeSigmaY() const : "
            << "creating yield stress field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ngbSigmaYPtr_ || ownSigmaYPtr_)
    {
        FatalErrorIn("ULLSMaterialInterface::makeSigmaY() const")
            << "interface yield stress field "
                << "already exists"
                << abort(FatalError);
    }

    ownSigmaYPtr_ = new scalarField(faces().size(), 0);
    ngbSigmaYPtr_ = new scalarField(faces().size(), 0);
}


void ULLSMaterialInterface::makeDSigmaY() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface::makeDSigmaY() const : "
            << "creating yield stress increment field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ngbDSigmaYPtr_ || ownDSigmaYPtr_)
    {
        FatalErrorIn("ULLSMaterialInterface::makeDSigmaY() const")
            << "interface yield stress increment field "
                << "already exists"
                << abort(FatalError);
    }

    ownDSigmaYPtr_ = new scalarField(faces().size(), 0);
    ngbDSigmaYPtr_ = new scalarField(faces().size(), 0);
}


void ULLSMaterialInterface::makeEpsilonPEq() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface::makeEpsilonPEq() const : "
            << "creating equivalent plastic strain field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ngbEpsilonPEqPtr_ || ownEpsilonPEqPtr_)
    {
        FatalErrorIn("ULLSMaterialInterface::makeEpsilonPEq() const")
            << "interface equivalent plastic strain field "
                << "already exists"
                << abort(FatalError);
    }

    ownEpsilonPEqPtr_ = new scalarField(faces().size(), 0);
    ngbEpsilonPEqPtr_ = new scalarField(faces().size(), 0);
}


void ULLSMaterialInterface::makeDEpsilonPEq() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface::makeDEpsilonPEq() const : "
            << "creating equivalent plastic strain increment field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ngbDEpsilonPEqPtr_ || ownDEpsilonPEqPtr_)
    {
        FatalErrorIn("ULLSMaterialInterface::makeDEpsilonPEq() const")
            << "interface equivalent plastic strain increment field "
                << "already exists"
                << abort(FatalError);
    }

    ownDEpsilonPEqPtr_ = new scalarField(faces().size(), 0);
    ngbDEpsilonPEqPtr_ = new scalarField(faces().size(), 0);
}


void ULLSMaterialInterface::makeSubMeshPointDD() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface::makeSubMeshPointDD() const : "
            << "creating point displacements increment fields"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshPointDD_.empty())
    {
        FatalErrorIn("ULLSMaterialInterface::makeSubMeshPointDD() const")
            << "Point displacement increment fields already exist"
            << abort(FatalError);
    }

    subMeshPointDD_.setSize(subMeshes().size());

    forAll(subMeshPointDD_, meshI)
    {
        subMeshPointDD_.set
        (
            meshI,
            new pointVectorField
            (
                subMeshes()[meshI].interpolate(pointDD_)
            )
        );
    }
}


void ULLSMaterialInterface::makeSubMeshDD() const
{
    if (debug)
    {
        Info<< "void ULLSMaterialInterface::makeSubMeshDD() const : "
            << "creating displacements fields"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshDD_.empty())
    {
        FatalErrorIn("ULLSMaterialInterface::makeSubMeshDD() const")
            << "Displacement increment fields already exist"
            << abort(FatalError);
    }

    subMeshDD_.setSize(subMeshes().size());

    forAll(subMeshDD_, meshI)
    {
        OStringStream SubsetName;
        SubsetName() << Pstream::myProcNo() << '_' << meshI << '_'
            << DD_.name();

        subMeshDD_.set
        (
            meshI,
            new volVectorField
            (
                IOobject
                (
                    word(SubsetName.str()),
                    subMeshes()[meshI].subMesh().time().timeName(),
                    subMeshes()[meshI].subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMeshes()[meshI].subMesh(),
                dimensionedVector("0", DD_.dimensions(), vector::zero)
            )
        );

        subMeshDD_[meshI] = subMeshes()[meshI].interpolate(DD_);

        // Correct displacement field at the interface
        const labelList& patchMap = subMeshes()[meshI].patchMap();
        label interfacePatchIndex = -1;
        forAll(patchMap, patchI)
        {
            if (patchMap[patchI] == -1)
            {
                interfacePatchIndex = patchI;
                break;
            }
        }

        if (interfacePatchIndex != -1)
        {
            vectorField& interfaceDD =
                subMeshDD_[meshI].boundaryField()[interfacePatchIndex];

            const labelList& fm = subMeshes()[meshI].faceMap();

            label interfacePatchStart =
                subMeshes()[meshI].subMesh().boundaryMesh()
                [
                    interfacePatchIndex
                ].start();

            forAll(interfaceDD, faceI)
            {
                label curInterFace =
                    findIndex(faces(), fm[interfacePatchStart + faceI]);

                interfaceDD[faceI] =
                    displacementIncrement()[curInterFace];
            }
        }
    }
}


void ULLSMaterialInterface::clearOut()
{
    deleteDemandDrivenData(displacementIncrementPtr_);
    deleteDemandDrivenData(displacementPtr_);
//     deleteDemandDrivenData(tractionIncrementPtr_);
    deleteDemandDrivenData(tractionPtr_);

//     deleteDemandDrivenData(ownDSigmafPtr_);
//     deleteDemandDrivenData(ngbDSigmafPtr_);

//     deleteDemandDrivenData(ownSigmafPtr_);
//     deleteDemandDrivenData(ngbSigmafPtr_);

    deleteDemandDrivenData(ownRelFPtr_);
    deleteDemandDrivenData(ngbRelFPtr_);
    deleteDemandDrivenData(ownFPtr_);
    deleteDemandDrivenData(ngbFPtr_);
    deleteDemandDrivenData(ownBbarPtr_);
    deleteDemandDrivenData(ngbBbarPtr_);

    deleteDemandDrivenData(ownDEpsilonPPtr_);
    deleteDemandDrivenData(ngbDEpsilonPPtr_);
    deleteDemandDrivenData(ownSigmaYPtr_);
    deleteDemandDrivenData(ngbSigmaYPtr_);
    deleteDemandDrivenData(ownDSigmaYPtr_);
    deleteDemandDrivenData(ngbDSigmaYPtr_);

    deleteDemandDrivenData(ownEpsilonPEqPtr_);
    deleteDemandDrivenData(ngbEpsilonPEqPtr_);
    deleteDemandDrivenData(ownDEpsilonPEqPtr_);
    deleteDemandDrivenData(ngbDEpsilonPEqPtr_);

    subMeshDD_.clear();
    subMeshPointDD_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


ULLSMaterialInterface::ULLSMaterialInterface
(
    const volVectorField& DD,
    const pointVectorField& pointDD,
    const volVectorField& D
)
:
    materialInterface(DD.mesh()),
    DD_(DD),
    pointDD_(pointDD),
    D_(D),
    curTimeIndex_(-1),
    displacementIncrementPtr_(NULL),
    displacementPtr_(NULL),
//     tractionIncrementPtr_(NULL),
    tractionPtr_(NULL),
//     ownDSigmafPtr_(NULL),
//     ngbDSigmafPtr_(NULL),
//     ownSigmafPtr_(NULL),
//     ngbSigmafPtr_(NULL),
    ownRelFPtr_(NULL),
    ngbRelFPtr_(NULL),
    ownFPtr_(NULL),
    ngbFPtr_(NULL),
    ownBbarPtr_(NULL),
    ngbBbarPtr_(NULL),
    ownDEpsilonPPtr_(NULL),
    ngbDEpsilonPPtr_(NULL),
    ownSigmaYPtr_(NULL),
    ngbSigmaYPtr_(NULL),
    ownDSigmaYPtr_(NULL),
    ngbDSigmaYPtr_(NULL),
    ownEpsilonPEqPtr_(NULL),
    ngbEpsilonPEqPtr_(NULL),
    ownDEpsilonPEqPtr_(NULL),
    ngbDEpsilonPEqPtr_(NULL),
    subMeshDD_(0),
    subMeshPointDD_(0)
{
    // Looking up solid solver
    const solidSolver& solid =
        mesh().lookupObject<solidSolver>
        (
            "solidProperties"
        );

    // If plasticity active
    if (solid.rheology().plasticityActive())
    {
        const unallocLabelList& owner = mesh().owner();
        const unallocLabelList& neighbour = mesh().neighbour();

        const volScalarField& sigmaY =
            mesh().lookupObject<volScalarField>("sigmaY");

        forAll(faces(), faceI)
        {
            label curFace = faces()[faceI];

            if (curFace < mesh().nInternalFaces())
            {
                ownSigmaY()[faceI] = sigmaY[owner[curFace]];
                ngbSigmaY()[faceI] = sigmaY[neighbour[curFace]];
            }
            else
            {
                label curPatch = mesh().boundaryMesh().whichPatch(curFace);
                label curPatchFace =
                    curFace - mesh().boundaryMesh()[curPatch].start();

                const unallocLabelList& faceCells =
                    mesh().boundary()[curPatch].faceCells();

                ownSigmaY()[faceI] =
                    sigmaY[faceCells[curPatchFace]];
                ngbSigmaY()[faceI] =
                    sigmaY.boundaryField()[curPatch][curPatchFace];
            }
        }
    }
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

ULLSMaterialInterface::~ULLSMaterialInterface()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

vectorField& ULLSMaterialInterface::displacementIncrement()
{
    if (!displacementIncrementPtr_)
    {
        makeDisplacementIncrement();
    }

    return *displacementIncrementPtr_;
}

const vectorField& ULLSMaterialInterface::displacementIncrement() const
{
    if (!displacementIncrementPtr_)
    {
        makeDisplacementIncrement();
    }

    return *displacementIncrementPtr_;
}

vectorField& ULLSMaterialInterface::displacement()
{
    if (!displacementPtr_)
    {
        makeDisplacement();
    }

    return *displacementPtr_;
}

const vectorField& ULLSMaterialInterface::displacement() const
{
    if (!displacementPtr_)
    {
        makeDisplacement();
    }

    return *displacementPtr_;
}

// vectorField& ULLSMaterialInterface::tractionIncrement()
// {
//     if (!tractionIncrementPtr_)
//     {
//         makeTractionIncrement();
//     }

//     return *tractionIncrementPtr_;
// }

// const vectorField& ULLSMaterialInterface::tractionIncrement() const
// {
//     if (!tractionIncrementPtr_)
//     {
//         makeTractionIncrement();
//     }

//     return *tractionIncrementPtr_;
// }

vectorField& ULLSMaterialInterface::traction()
{
    if (!tractionPtr_)
    {
        makeTraction();
    }

    return *tractionPtr_;
}

const vectorField& ULLSMaterialInterface::traction() const
{
    if (!tractionPtr_)
    {
        makeTraction();
    }

    return *tractionPtr_;
}

// symmTensorField& ULLSMaterialInterface::ownSigmaf()
// {
//     if (!ownSigmafPtr_)
//     {
//         makeSigmaf();
//     }

//     return *ownSigmafPtr_;
// }

// const symmTensorField& ULLSMaterialInterface::ownSigmaf() const
// {
//     if (!ownSigmafPtr_)
//     {
//         makeSigmaf();
//     }

//     return *ownSigmafPtr_;
// }

// symmTensorField& ULLSMaterialInterface::ngbSigmaf()
// {
//     if (!ngbSigmafPtr_)
//     {
//         makeSigmaf();
//     }

//     return *ngbSigmafPtr_;
// }

// const symmTensorField& ULLSMaterialInterface::ngbSigmaf() const
// {
//     if (!ngbSigmafPtr_)
//     {
//         makeSigmaf();
//     }

//     return *ngbSigmafPtr_;
// }

// symmTensorField& ULLSMaterialInterface::ownDSigmaf()
// {
//     if (!ownDSigmafPtr_)
//     {
//         makeDSigmaf();
//     }

//     return *ownDSigmafPtr_;
// }

// const symmTensorField& ULLSMaterialInterface::ownDSigmaf() const
// {
//     if (!ownDSigmafPtr_)
//     {
//         makeDSigmaf();
//     }

//     return *ownDSigmafPtr_;
// }

// symmTensorField& ULLSMaterialInterface::ngbDSigmaf()
// {
//     if (!ngbDSigmafPtr_)
//     {
//         makeDSigmaf();
//     }

//     return *ngbDSigmafPtr_;
// }

// const symmTensorField& ULLSMaterialInterface::ngbDSigmaf() const
// {
//     if (!ngbDSigmafPtr_)
//     {
//         makeDSigmaf();
//     }

//     return *ngbDSigmafPtr_;
// }

tensorField& ULLSMaterialInterface::ownRelF()
{
    if (!ownRelFPtr_)
    {
        makeRelF();
    }

    return *ownRelFPtr_;
}

tensorField& ULLSMaterialInterface::ngbRelF()
{
    if (!ngbRelFPtr_)
    {
        makeRelF();
    }

    return *ngbRelFPtr_;
}

const tensorField& ULLSMaterialInterface::ownRelF() const
{
    if (!ownRelFPtr_)
    {
        makeRelF();
    }

    return *ownRelFPtr_;
}

const tensorField& ULLSMaterialInterface::ngbRelF() const
{
    if (!ngbRelFPtr_)
    {
        makeRelF();
    }

    return *ngbRelFPtr_;
}

tensorField& ULLSMaterialInterface::ownF()
{
    if (!ownFPtr_)
    {
        makeF();
    }

    return *ownFPtr_;
}

tensorField& ULLSMaterialInterface::ngbF()
{
    if (!ngbFPtr_)
    {
        makeF();
    }

    return *ngbFPtr_;
}

const tensorField& ULLSMaterialInterface::ownF() const
{
    if (!ownFPtr_)
    {
        makeF();
    }

    return *ownFPtr_;
}

const tensorField& ULLSMaterialInterface::ngbF() const
{
    if (!ngbFPtr_)
    {
        makeF();
    }

    return *ngbFPtr_;
}

symmTensorField& ULLSMaterialInterface::ownBbar()
{
    if (!ownBbarPtr_)
    {
        makeBbar();
    }

    return *ownBbarPtr_;
}

symmTensorField& ULLSMaterialInterface::ngbBbar()
{
    if (!ngbBbarPtr_)
    {
        makeBbar();
    }

    return *ngbBbarPtr_;
}

const symmTensorField& ULLSMaterialInterface::ownBbar() const
{
    if (!ownBbarPtr_)
    {
        makeBbar();
    }

    return *ownBbarPtr_;
}

const symmTensorField& ULLSMaterialInterface::ngbBbar() const
{
    if (!ngbBbarPtr_)
    {
        makeBbar();
    }

    return *ngbBbarPtr_;
}



symmTensorField& ULLSMaterialInterface::ownDEpsilonP()
{
    if (!ownDEpsilonPPtr_)
    {
        makeDEpsilonP();
    }

    return *ownDEpsilonPPtr_;
}

symmTensorField& ULLSMaterialInterface::ngbDEpsilonP()
{
    if (!ngbDEpsilonPPtr_)
    {
        makeDEpsilonP();
    }

    return *ngbDEpsilonPPtr_;
}

const symmTensorField& ULLSMaterialInterface::ownDEpsilonP() const
{
    if (!ownDEpsilonPPtr_)
    {
        makeDEpsilonP();
    }

    return *ownDEpsilonPPtr_;
}

const symmTensorField& ULLSMaterialInterface::ngbDEpsilonP() const
{
    if (!ngbDEpsilonPPtr_)
    {
        makeDEpsilonP();
    }

    return *ngbDEpsilonPPtr_;
}



scalarField& ULLSMaterialInterface::ownSigmaY()
{
    if (!ownSigmaYPtr_)
    {
        makeSigmaY();
    }

    return *ownSigmaYPtr_;
}

scalarField& ULLSMaterialInterface::ngbSigmaY()
{
    if (!ngbSigmaYPtr_)
    {
        makeSigmaY();
    }

    return *ngbSigmaYPtr_;
}

const scalarField& ULLSMaterialInterface::ownSigmaY() const
{
    if (!ownSigmaYPtr_)
    {
        makeSigmaY();
    }

    return *ownSigmaYPtr_;
}

const scalarField& ULLSMaterialInterface::ngbSigmaY() const
{
    if (!ngbSigmaYPtr_)
    {
        makeSigmaY();
    }

    return *ngbSigmaYPtr_;
}

scalarField& ULLSMaterialInterface::ownDSigmaY()
{
    if (!ownDSigmaYPtr_)
    {
        makeDSigmaY();
    }

    return *ownDSigmaYPtr_;
}

scalarField& ULLSMaterialInterface::ngbDSigmaY()
{
    if (!ngbDSigmaYPtr_)
    {
        makeDSigmaY();
    }

    return *ngbDSigmaYPtr_;
}

const scalarField& ULLSMaterialInterface::ownDSigmaY() const
{
    if (!ownDSigmaYPtr_)
    {
        makeDSigmaY();
    }

    return *ownDSigmaYPtr_;
}

const scalarField& ULLSMaterialInterface::ngbDSigmaY() const
{
    if (!ngbDSigmaYPtr_)
    {
        makeDSigmaY();
    }

    return *ngbDSigmaYPtr_;
}



scalarField& ULLSMaterialInterface::ownEpsilonPEq()
{
    if (!ownEpsilonPEqPtr_)
    {
        makeEpsilonPEq();
    }

    return *ownEpsilonPEqPtr_;
}

scalarField& ULLSMaterialInterface::ngbEpsilonPEq()
{
    if (!ngbEpsilonPEqPtr_)
    {
        makeEpsilonPEq();
    }

    return *ngbEpsilonPEqPtr_;
}

const scalarField& ULLSMaterialInterface::ownEpsilonPEq() const
{
    if (!ownEpsilonPEqPtr_)
    {
        makeEpsilonPEq();
    }

    return *ownEpsilonPEqPtr_;
}

const scalarField& ULLSMaterialInterface::ngbEpsilonPEq() const
{
    if (!ngbEpsilonPEqPtr_)
    {
        makeEpsilonPEq();
    }

    return *ngbEpsilonPEqPtr_;
}



scalarField& ULLSMaterialInterface::ownDEpsilonPEq()
{
    if (!ownDEpsilonPEqPtr_)
    {
        makeDEpsilonPEq();
    }

    return *ownDEpsilonPEqPtr_;
}

scalarField& ULLSMaterialInterface::ngbDEpsilonPEq()
{
    if (!ngbDEpsilonPEqPtr_)
    {
        makeDEpsilonPEq();
    }

    return *ngbDEpsilonPEqPtr_;
}

const scalarField& ULLSMaterialInterface::ownDEpsilonPEq() const
{
    if (!ownDEpsilonPEqPtr_)
    {
        makeDEpsilonPEq();
    }

    return *ownDEpsilonPEqPtr_;
}

const scalarField& ULLSMaterialInterface::ngbDEpsilonPEq() const
{
    if (!ngbDEpsilonPEqPtr_)
    {
        makeDEpsilonPEq();
    }

    return *ngbDEpsilonPEqPtr_;
}

const PtrList<volVectorField>& ULLSMaterialInterface::subMeshDD() const
{
    if (subMeshDD_.empty())
    {
        makeSubMeshDD();
    }

    return subMeshDD_;
}

PtrList<volVectorField>& ULLSMaterialInterface::subMeshDD()
{
    if (subMeshDD_.empty())
    {
        makeSubMeshDD();
    }

    return subMeshDD_;
}

const PtrList<pointVectorField>& ULLSMaterialInterface::subMeshPointDD() const
{
    if (subMeshPointDD_.empty())
    {
        makeSubMeshPointDD();
    }

    return subMeshPointDD_;
}

PtrList<pointVectorField>& ULLSMaterialInterface::subMeshPointDD()
{
    if (subMeshPointDD_.empty())
    {
        makeSubMeshPointDD();
    }

    return subMeshPointDD_;
}

void ULLSMaterialInterface::correct(fvVectorMatrix& DEqn)
{}


void ULLSMaterialInterface::correct(surfaceVectorField& trac) const
{
    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh().nInternalFaces())
        {
            trac.internalField()[curFace] = traction()[faceI];
        }
        else
        {
            label curPatch = mesh().boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh().boundaryMesh()[curPatch].start();

            trac.boundaryField()[curPatch][curPatchFace] = traction()[faceI];
        }
    }
}


tmp<symmTensorField> ULLSMaterialInterface::trialBbar
(
    const label faceI,
    const volTensorField& gradDD,
    const surfaceTensorField& gradDDf
) const
{
    tmp<symmTensorField> tTrialBbar(new symmTensorField(2, symmTensor::zero));

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const vectorField& C  = mesh().C().internalField();
    const vectorField& Cf  = mesh().Cf().internalField();
    const vectorField& SI = mesh().Sf().internalField();
    const scalarField& magSI  = mesh().magSf().internalField();
    const scalarField& deltaCoeffs = mesh().deltaCoeffs().internalField();
    const scalarField& w = mesh().weights().internalField();

    const vectorField& DDI = DD_.internalField();

    const tensorField& gradDDI = gradDD.internalField();
    const tensorField& gradDDfI = gradDDf.internalField();

    label curFace = faces()[faceI];

    const vectorField& interDD = displacementIncrement();

    // Internal faces
    if (curFace < mesh().nInternalFaces())
    {
        label curOwner = owner[curFace];
        label curNeighbour = neighbour[curFace];

        vector ownN = SI[curFace]/magSI[curFace];
        vector ngbN = -ownN;

        vector ownCorrVec = Cf[curFace] - C[curOwner];
        ownCorrVec -= ownN*(ownN & ownCorrVec);
        vector ngbCorrVec = Cf[curFace] - C[curNeighbour];
        ngbCorrVec -= ngbN*(ngbN&ngbCorrVec);

        vector ownDD = DDI[curOwner];
        ownDD += (ownCorrVec & gradDDI[curOwner]);
        vector ngbDD = DDI[curNeighbour];
        ngbDD += (ngbCorrVec & gradDDI[curNeighbour]);

        scalar ngbDeltaN = w[curFace]*(1.0/deltaCoeffs[curFace]);
        scalar ownDeltaN = (1.0/deltaCoeffs[curFace]) - ngbDeltaN;

        tensor ownGradDD = gradDDfI[curFace];
        ownGradDD +=
            ownN*(interDD[faceI] - ownDD)/ownDeltaN
          - ownN*(ownN & ownGradDD);

        tensor ngbGradDD = gradDDfI[curFace];
        ngbGradDD +=
            ngbN*(interDD[faceI] - ngbDD)/ngbDeltaN
          - ngbN*(ngbN & ngbGradDD);

        tensor ownRelFcur = I + ownGradDD.T();
        tensor ngbRelFcur = I + ngbGradDD.T();

        scalar ownRelJ = det(ownRelFcur);
        scalar ngbRelJ = det(ngbRelFcur);

        if (ownRelJ < 0 || ngbRelJ < 0)
        {
            Warning << "Determinant of relative deformation gradient "
                << "is less then zero at the interface. "
                << "Performing stabilisation procedure." << endl;

            ownRelFcur = I;
            ngbRelFcur = I;

            ownRelJ = 1;
            ngbRelJ = 1;
        }

//         tensor ownInvRelF = hinv(ownRelFcur);
//         tensor ngbInvRelF = hinv(ngbRelFcur);

//         vector ownS = (ownInvRelF.T() & SI[curFace])*ownRelJ;
//         vector ngbS = (ngbInvRelF.T() & SI[curFace])*ngbRelJ;
//         vector S = 0.5*(ownS + ngbS);
//         vector N = S/(mag(S) + SMALL);

        tensor ownRelFbar = pow(ownRelJ, -1.0/3.0)*ownRelFcur;
        tensor ngbRelFbar = pow(ngbRelJ, -1.0/3.0)*ngbRelFcur;

//         tensor ownFcur = (ownRelFcur & ownF()[faceI]);
//         tensor ngbFcur = (ngbRelFcur & ngbF()[faceI]);

//         scalar ownJ = det(ownFcur);
//         scalar ngbJ = det(ngbFcur);

        symmTensor ownBbarCur =
            symm(ownRelFbar & ownBbar()[faceI] & ownRelFbar.T());
        symmTensor ngbBbarCur =
            symm(ngbRelFbar & ngbBbar()[faceI] & ngbRelFbar.T());

        tTrialBbar()[0] = ownBbarCur;
        tTrialBbar()[1] = ngbBbarCur;
    }
    else
    {
        label curPatch = mesh().boundaryMesh().whichPatch(curFace);
        label curPatchFace =
            curFace - mesh().boundaryMesh()[curPatch].start();

        const unallocLabelList& faceCells =
            mesh().boundary()[curPatch].faceCells();

        label curOwner = faceCells[curPatchFace];

        vector ownN =
            mesh().Sf().boundaryField()[curPatch][curPatchFace]
           /mesh().magSf().boundaryField()[curPatch][curPatchFace];
        vector ngbN = -ownN;

        scalar ngbDeltaN =
            mesh().weights().boundaryField()[curPatch][curPatchFace]
           *(
                1.0
               /mesh().deltaCoeffs().boundaryField()[curPatch][curPatchFace]
            );
        scalar ownDeltaN =
            (
                1.0
               /mesh().deltaCoeffs().boundaryField()[curPatch][curPatchFace]
            )
          - ngbDeltaN;

        vector ownDD = DD_.internalField()[curOwner];
        vector ngbDD = DD_.boundaryField()[curPatch][curPatchFace];

        vector ownCorrVec =
            mesh().Cf().boundaryField()[curPatch][curPatchFace]
          - mesh().C().internalField()[curOwner];
        ownCorrVec -= ownN*(ownN & ownCorrVec);

        vector ngbCorrVec =
            mesh().Cf().boundaryField()[curPatch][curPatchFace]
          - mesh().C().boundaryField()[curPatch][curPatchFace];
        ngbCorrVec -= ngbN*(ngbN & ngbCorrVec);

        vector ownDDCorr =
            (
                ownCorrVec
              & gradDD.internalField()[curOwner]
            );
        vector ngbDDCorr =
            (
                ngbCorrVec
              & gradDD.boundaryField()[curPatch][curPatchFace]
            );

        ownDD += ownDDCorr;
        ngbDD += ngbDDCorr;

        tensor ownGradDD = gradDDf.boundaryField()[curPatch][curPatchFace];
        ownGradDD +=
            ownN*(interDD[faceI] - ownDD)/ownDeltaN
          - ownN*(ownN & ownGradDD);

        tensor ngbGradDD = gradDDf.boundaryField()[curPatch][curPatchFace];
        ngbGradDD +=
            ngbN*(interDD[faceI] - ngbDD)/ngbDeltaN
          - ngbN*(ngbN & ngbGradDD);

        tensor ownRelFcur = I + ownGradDD.T();
        tensor ngbRelFcur = I + ngbGradDD.T();

        scalar ownRelJ = det(ownRelFcur);
        scalar ngbRelJ = det(ngbRelFcur);

        if (ownRelJ < 0 || ngbRelJ < 0)
        {
            Warning << "Determinant of relative deformation gradient "
                << "is less then zero at the interface. "
                << "Performing stabilisation procedure." << endl;

            ownRelFcur = I;
            ngbRelFcur = I;

            ownRelJ = 1;
            ngbRelJ = 1;
        }

//         tensor ownInvRelF = hinv(ownRelFcur);
//         tensor ngbInvRelF = hinv(ngbRelFcur);

//         vector ownS = (ownInvRelF.T() & SI[curFace])*ownRelJ;
//         vector ngbS = (ngbInvRelF.T() & SI[curFace])*ngbRelJ;
//         vector S = 0.5*(ownS + ngbS);
//         vector N = S/(mag(S) + SMALL);

        tensor ownRelFbar = pow(ownRelJ, -1.0/3.0)*ownRelFcur;
        tensor ngbRelFbar = pow(ngbRelJ, -1.0/3.0)*ngbRelFcur;

//         tensor ownFcur = (ownRelFcur & ownF()[faceI]);
//         tensor ngbFcur = (ngbRelFcur & ngbF()[faceI]);

//         scalar ownJ = det(ownFcur);
//         scalar ngbJ = det(ngbFcur);

        symmTensor ownBbarCur =
            symm(ownRelFbar & ownBbar()[faceI] & ownRelFbar.T());
        symmTensor ngbBbarCur =
            symm(ngbRelFbar & ngbBbar()[faceI] & ngbRelFbar.T());

        tTrialBbar()[0] = ownBbarCur;
        tTrialBbar()[1] = ngbBbarCur;
    }

    return tTrialBbar;
}


tmp<scalarField> ULLSMaterialInterface::J
(
    const label faceI,
    const volTensorField& gradDD,
    const surfaceTensorField& gradDDf
) const
{
    tmp<scalarField> tJ(new scalarField(2, 0));

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const vectorField& C  = mesh().C().internalField();
    const vectorField& Cf  = mesh().Cf().internalField();
    const vectorField& SI = mesh().Sf().internalField();
    const scalarField& magSI  = mesh().magSf().internalField();
    const scalarField& deltaCoeffs = mesh().deltaCoeffs().internalField();
    const scalarField& w = mesh().weights().internalField();

    const vectorField& DDI = DD_.internalField();

    const tensorField& gradDDI = gradDD.internalField();
    const tensorField& gradDDfI = gradDDf.internalField();

    label curFace = faces()[faceI];

    const vectorField& interDD = displacementIncrement();

    // Internal faces
    if (curFace < mesh().nInternalFaces())
    {
        label curOwner = owner[curFace];
        label curNeighbour = neighbour[curFace];

        vector ownN = SI[curFace]/magSI[curFace];
        vector ngbN = -ownN;

        vector ownCorrVec = Cf[curFace] - C[curOwner];
        ownCorrVec -= ownN*(ownN & ownCorrVec);
        vector ngbCorrVec = Cf[curFace] - C[curNeighbour];
        ngbCorrVec -= ngbN*(ngbN&ngbCorrVec);

        vector ownDD = DDI[curOwner];
        ownDD += (ownCorrVec & gradDDI[curOwner]);
        vector ngbDD = DDI[curNeighbour];
        ngbDD += (ngbCorrVec & gradDDI[curNeighbour]);

        scalar ngbDeltaN = w[curFace]*(1.0/deltaCoeffs[curFace]);
        scalar ownDeltaN = (1.0/deltaCoeffs[curFace]) - ngbDeltaN;

        tensor ownGradDD = gradDDfI[curFace];
        ownGradDD +=
            ownN*(interDD[faceI] - ownDD)/ownDeltaN
          - ownN*(ownN & ownGradDD);

        tensor ngbGradDD = gradDDfI[curFace];
        ngbGradDD +=
            ngbN*(interDD[faceI] - ngbDD)/ngbDeltaN
          - ngbN*(ngbN & ngbGradDD);

        tensor ownRelFcur = I + ownGradDD.T();
        tensor ngbRelFcur = I + ngbGradDD.T();

        scalar ownRelJ = det(ownRelFcur);
        scalar ngbRelJ = det(ngbRelFcur);

        if (ownRelJ < 0 || ngbRelJ < 0)
        {
            Warning << "Determinant of relative deformation gradient "
                << "is less then zero at the interface. "
                << "Performing stabilisation procedure." << endl;

            ownRelFcur = I;
            ngbRelFcur = I;

            ownRelJ = 1;
            ngbRelJ = 1;
        }

//         tensor ownRelFbar = pow(ownRelJ, -1.0/3.0)*ownRelFcur;
//         tensor ngbRelFbar = pow(ngbRelJ, -1.0/3.0)*ngbRelFcur;

        tensor ownFcur = (ownRelFcur & ownF()[faceI]);
        tensor ngbFcur = (ngbRelFcur & ngbF()[faceI]);

        scalar ownJ = det(ownFcur);
        scalar ngbJ = det(ngbFcur);

        tJ()[0] = ownJ;
        tJ()[1] = ngbJ;
    }
    else
    {
        label curPatch = mesh().boundaryMesh().whichPatch(curFace);
        label curPatchFace =
            curFace - mesh().boundaryMesh()[curPatch].start();

        const unallocLabelList& faceCells =
            mesh().boundary()[curPatch].faceCells();

        label curOwner = faceCells[curPatchFace];

        vector ownN =
            mesh().Sf().boundaryField()[curPatch][curPatchFace]
           /mesh().magSf().boundaryField()[curPatch][curPatchFace];
        vector ngbN = -ownN;

        scalar ngbDeltaN =
            mesh().weights().boundaryField()[curPatch][curPatchFace]
           *(
                1.0
               /mesh().deltaCoeffs().boundaryField()[curPatch][curPatchFace]
            );
        scalar ownDeltaN =
            (
                1.0
               /mesh().deltaCoeffs().boundaryField()[curPatch][curPatchFace]
            )
          - ngbDeltaN;

        vector ownDD = DD_.internalField()[curOwner];
        vector ngbDD = DD_.boundaryField()[curPatch][curPatchFace];

        vector ownCorrVec =
            mesh().Cf().boundaryField()[curPatch][curPatchFace]
          - mesh().C().internalField()[curOwner];
        ownCorrVec -= ownN*(ownN & ownCorrVec);

        vector ngbCorrVec =
            mesh().Cf().boundaryField()[curPatch][curPatchFace]
          - mesh().C().boundaryField()[curPatch][curPatchFace];
        ngbCorrVec -= ngbN*(ngbN & ngbCorrVec);

        vector ownDDCorr =
            (
                ownCorrVec
              & gradDD.internalField()[curOwner]
            );
        vector ngbDDCorr =
            (
                ngbCorrVec
              & gradDD.boundaryField()[curPatch][curPatchFace]
            );

        ownDD += ownDDCorr;
        ngbDD += ngbDDCorr;

        tensor ownGradDD = gradDDf.boundaryField()[curPatch][curPatchFace];
        ownGradDD +=
            ownN*(interDD[faceI] - ownDD)/ownDeltaN
          - ownN*(ownN & ownGradDD);

        tensor ngbGradDD = gradDDf.boundaryField()[curPatch][curPatchFace];
        ngbGradDD +=
            ngbN*(interDD[faceI] - ngbDD)/ngbDeltaN
          - ngbN*(ngbN & ngbGradDD);

        tensor ownRelFcur = I + ownGradDD.T();
        tensor ngbRelFcur = I + ngbGradDD.T();

        scalar ownRelJ = det(ownRelFcur);
        scalar ngbRelJ = det(ngbRelFcur);

        if (ownRelJ < 0 || ngbRelJ < 0)
        {
            Warning << "Determinant of relative deformation gradient "
                << "is less then zero at the interface. "
                << "Performing stabilisation procedure." << endl;

            ownRelFcur = I;
            ngbRelFcur = I;

            ownRelJ = 1;
            ngbRelJ = 1;
        }

//         tensor ownRelFbar = pow(ownRelJ, -1.0/3.0)*ownRelFcur;
//         tensor ngbRelFbar = pow(ngbRelJ, -1.0/3.0)*ngbRelFcur;

        tensor ownFcur = (ownRelFcur & ownF()[faceI]);
        tensor ngbFcur = (ngbRelFcur & ngbF()[faceI]);

        scalar ownJ = det(ownFcur);
        scalar ngbJ = det(ngbFcur);

        tJ()[0] = ownJ;
        tJ()[1] = ngbJ;
    }

    return tJ;
}


tmp<vectorField> ULLSMaterialInterface::cauchyTraction
(
    const vector& interDD,
    const label faceI,
    const volTensorField& gradDD,
    const surfaceTensorField& gradDDf,
    const volScalarField& mu,
    const volScalarField& lambda,
    const bool plasticityActive
) const
{
    tmp<vectorField> tCauchyTraction(new vectorField(2, vector::zero));

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const vectorField& C  = mesh().C().internalField();
    const vectorField& Cf  = mesh().Cf().internalField();
    const vectorField& SI = mesh().Sf().internalField();
    const scalarField& magSI  = mesh().magSf().internalField();
    const scalarField& deltaCoeffs = mesh().deltaCoeffs().internalField();
    const scalarField& w = mesh().weights().internalField();

    const volVectorField& DD = DD_.prevIter();
    const vectorField& DDI = DD.internalField();

    const tensorField& gradDDI = gradDD.internalField();
    const tensorField& gradDDfI = gradDDf.internalField();

    label curFace = faces()[faceI];

    // Internal faces
    if (curFace < mesh().nInternalFaces())
    {
        label curOwner = owner[curFace];
        label curNeighbour = neighbour[curFace];

        vector ownN = SI[curFace]/magSI[curFace];
        vector ngbN = -ownN;

        vector ownCorrVec = Cf[curFace] - C[curOwner];
        ownCorrVec -= ownN*(ownN & ownCorrVec);
        vector ngbCorrVec = Cf[curFace] - C[curNeighbour];
        ngbCorrVec -= ngbN*(ngbN&ngbCorrVec);

        vector ownDD = DDI[curOwner];
        ownDD += (ownCorrVec & gradDDI[curOwner]);
        vector ngbDD = DDI[curNeighbour];
        ngbDD += (ngbCorrVec & gradDDI[curNeighbour]);

        scalar ngbDeltaN = w[curFace]*(1.0/deltaCoeffs[curFace]);
        scalar ownDeltaN = (1.0/deltaCoeffs[curFace]) - ngbDeltaN;

        tensor ownGradDD = gradDDfI[curFace];
        ownGradDD +=
            ownN*(interDD - ownDD)/ownDeltaN
          - ownN*(ownN & ownGradDD);

        tensor ngbGradDD = gradDDfI[curFace];
        ngbGradDD +=
            ngbN*(interDD - ngbDD)/ngbDeltaN
          - ngbN*(ngbN & ngbGradDD);

        tensor ownRelFcur = I + ownGradDD.T();
        tensor ngbRelFcur = I + ngbGradDD.T();

        scalar ownRelJ = det(ownRelFcur);
        scalar ngbRelJ = det(ngbRelFcur);

        if (ownRelJ < 0 || ngbRelJ < 0)
        {
            Warning << "Determinant of relative deformation gradient "
                << "is less then zero at the interface. "
                << "Performing stabilisation procedure." << endl;

            ownRelFcur = I;
            ngbRelFcur = I;

            ownRelJ = 1;
            ngbRelJ = 1;
        }

        tensor ownInvRelF = hinv(ownRelFcur);
        tensor ngbInvRelF = hinv(ngbRelFcur);

        vector ownS = (ownInvRelF.T() & SI[curFace])*ownRelJ;
        vector ngbS = (ngbInvRelF.T() & SI[curFace])*ngbRelJ;
        vector S = 0.5*(ownS + ngbS);
        vector N = S/(mag(S) + SMALL);

        tensor ownRelFbar = pow(ownRelJ, -1.0/3.0)*ownRelFcur;
        tensor ngbRelFbar = pow(ngbRelJ, -1.0/3.0)*ngbRelFcur;

        tensor ownFcur = (ownRelFcur & ownF()[faceI]);
        tensor ngbFcur = (ngbRelFcur & ngbF()[faceI]);

        scalar ownJ = det(ownFcur);
        scalar ngbJ = det(ngbFcur);

        symmTensor ownBbarCur =
            symm(ownRelFbar & ownBbar()[faceI] & ownRelFbar.T());
        symmTensor ngbBbarCur =
            symm(ngbRelFbar & ngbBbar()[faceI] & ngbRelFbar.T());

        scalar ownMu = mu.internalField()[curOwner];
        scalar ngbMu = mu.internalField()[curNeighbour];

        scalar ownLambda = lambda.internalField()[curOwner];
        scalar ngbLambda = lambda.internalField()[curNeighbour];

        scalar ownK  = (ownLambda + (2.0/3.0)*ownMu);
        scalar ngbK  = (ngbLambda + (2.0/3.0)*ngbMu);

        symmTensor ownTau =
            ownMu*dev(ownBbarCur)
          + 0.5*ownK*(sqr(ownJ) - 1.0)*I;
        symmTensor ngbTau =
            ngbMu*dev(ngbBbarCur)
          + 0.5*ngbK*(sqr(ngbJ) - 1.0)*I;

        if (plasticityActive)
        {
            ownTau -=  2*ownMu*ownDEpsilonP()[faceI];
            ngbTau -=  2*ngbMu*ngbDEpsilonP()[faceI];
        }

        tCauchyTraction()[0] = (N & ownTau)/ownJ;
        tCauchyTraction()[1] = (N & ngbTau)/ngbJ;
    }
    else
    {
        label curPatch = mesh().boundaryMesh().whichPatch(curFace);
        label curPatchFace =
            curFace - mesh().boundaryMesh()[curPatch].start();

        const unallocLabelList& faceCells =
            mesh().boundary()[curPatch].faceCells();

        label curOwner = faceCells[curPatchFace];

        vector ownN =
            mesh().Sf().boundaryField()[curPatch][curPatchFace]
           /mesh().magSf().boundaryField()[curPatch][curPatchFace];
        vector ngbN = -ownN;

        scalar ngbDeltaN =
            mesh().weights().boundaryField()[curPatch][curPatchFace]
           *(
                1.0
               /mesh().deltaCoeffs().boundaryField()[curPatch][curPatchFace]
            );
        scalar ownDeltaN =
            (
                1.0
               /mesh().deltaCoeffs().boundaryField()[curPatch][curPatchFace]
            )
          - ngbDeltaN;

        vector ownDD = DD.internalField()[curOwner];
        vector ngbDD = DD.boundaryField()[curPatch][curPatchFace];

        vector ownCorrVec =
            mesh().Cf().boundaryField()[curPatch][curPatchFace]
          - mesh().C().internalField()[curOwner];
        ownCorrVec -= ownN*(ownN & ownCorrVec);

        vector ngbCorrVec =
            mesh().Cf().boundaryField()[curPatch][curPatchFace]
          - mesh().C().boundaryField()[curPatch][curPatchFace];
        ngbCorrVec -= ngbN*(ngbN & ngbCorrVec);

        vector ownDDCorr =
            (
                ownCorrVec
              & gradDD.internalField()[curOwner]
            );
        vector ngbDDCorr =
            (
                ngbCorrVec
              & gradDD.boundaryField()[curPatch][curPatchFace]
            );

        ownDD += ownDDCorr;
        ngbDD += ngbDDCorr;

        tensor ownGradDD = gradDDf.boundaryField()[curPatch][curPatchFace];
        ownGradDD +=
            ownN*(interDD - ownDD)/ownDeltaN
          - ownN*(ownN & ownGradDD);

        tensor ngbGradDD = gradDDf.boundaryField()[curPatch][curPatchFace];
        ngbGradDD +=
            ngbN*(interDD - ngbDD)/ngbDeltaN
          - ngbN*(ngbN & ngbGradDD);

        tensor ownRelFcur = I + ownGradDD.T();
        tensor ngbRelFcur = I + ngbGradDD.T();

        scalar ownRelJ = det(ownRelFcur);
        scalar ngbRelJ = det(ngbRelFcur);

        if (ownRelJ < 0 || ngbRelJ < 0)
        {
            Warning << "Determinant of relative deformation gradient "
                << "is less then zero at the interface. "
                << "Performing stabilisation procedure." << endl;

            ownRelFcur = I;
            ngbRelFcur = I;

            ownRelJ = 1;
            ngbRelJ = 1;
        }

        tensor ownInvRelF = hinv(ownRelFcur);
        tensor ngbInvRelF = hinv(ngbRelFcur);

        vector ownS = (ownInvRelF.T() & SI[curFace])*ownRelJ;
        vector ngbS = (ngbInvRelF.T() & SI[curFace])*ngbRelJ;
        vector S = 0.5*(ownS + ngbS);
        vector N = S/(mag(S) + SMALL);

        tensor ownRelFbar = pow(ownRelJ, -1.0/3.0)*ownRelFcur;
        tensor ngbRelFbar = pow(ngbRelJ, -1.0/3.0)*ngbRelFcur;

        tensor ownFcur = (ownRelFcur & ownF()[faceI]);
        tensor ngbFcur = (ngbRelFcur & ngbF()[faceI]);

        scalar ownJ = det(ownFcur);
        scalar ngbJ = det(ngbFcur);

        symmTensor ownBbarCur =
            symm(ownRelFbar & ownBbar()[faceI] & ownRelFbar.T());
        symmTensor ngbBbarCur =
            symm(ngbRelFbar & ngbBbar()[faceI] & ngbRelFbar.T());

        scalar ownMu  = mu.internalField()[curOwner];
        scalar ngbMu  = mu.boundaryField()[curPatch][curPatchFace];

        scalar ownLambda = lambda.internalField()[curOwner];
        scalar ngbLambda = lambda.boundaryField()[curPatch][curPatchFace];

        scalar ownK  = (ownLambda + (2.0/3.0)*ownMu);
        scalar ngbK  = (ngbLambda + (2.0/3.0)*ngbMu);

        symmTensor ownTau =
            ownMu*dev(ownBbarCur)
          + 0.5*ownK*(sqr(ownJ) - 1.0)*I;
        symmTensor ngbTau =
            ngbMu*dev(ngbBbarCur)
          + 0.5*ngbK*(sqr(ngbJ) - 1.0)*I;

        if (plasticityActive)
        {
            ownTau -=  2*ownMu*ownDEpsilonP()[faceI];
            ngbTau -=  2*ngbMu*ngbDEpsilonP()[faceI];
        }

        tCauchyTraction()[0] = (N & ownTau)/ownJ;
        tCauchyTraction()[1] = (N & ngbTau)/ngbJ;
    }

    return tCauchyTraction;
}


void ULLSMaterialInterface::updateTotalFields()
{
//     const unallocLabelList& owner = mesh().owner();
//     const unallocLabelList& neighbour = mesh().neighbour();

//     const vectorField& DDI = DD_.internalField();

//     const volTensorField& gradDD =
//         mesh().lookupObject<volTensorField>("grad(" + DD_.name() + ')');
//     const tensorField& gradDDI = gradDD.internalField();

//     const surfaceTensorField& gradDDf =
//         mesh().lookupObject<surfaceTensorField>("grad" + DD_.name() + 'f');
//     const tensorField& gradDDfI = gradDDf.internalField();

//     const vectorField& C  = mesh().C().internalField();
//     const vectorField& Cf  = mesh().Cf().internalField();

//     const vectorField& SI = mesh().Sf().internalField();
//     const scalarField& magSI  = mesh().magSf().internalField();
//     const scalarField& deltaCoeffs = mesh().deltaCoeffs().internalField();
//     const scalarField& w = mesh().weights().internalField();

//     vectorField& interDD = displacementIncrement();

//     // Interface faces
//     tensorField ownGradDD(faces().size(), tensor::zero);
//     tensorField ngbGradDD(faces().size(), tensor::zero);
//     forAll(faces(), faceI)
//     {
//         label curFace = faces()[faceI];

//         // Internal faces
//         if (curFace < mesh().nInternalFaces())
//         {
//             label curOwner = owner[curFace];
//             label curNeighbour = neighbour[curFace];

//             vector ownDD = DDI[curOwner];
//             vector ngbDD = DDI[curNeighbour];

//             vector ownN = SI[curFace]/magSI[curFace];
//             vector ngbN = -ownN;

//             scalar ngbDeltaN = w[curFace]*(1.0/deltaCoeffs[curFace]);
//             scalar ownDeltaN = (1.0/deltaCoeffs[curFace]) - ngbDeltaN;

//             vector ownCorrVec = Cf[curFace] - C[curOwner];
//             ownCorrVec -= ownN*(ownN&ownCorrVec);
//             vector ngbCorrVec = Cf[curFace] - C[curNeighbour];
//             ngbCorrVec -= ngbN*(ngbN&ngbCorrVec);

//             ownDD += (ownCorrVec & gradDDI[curOwner]);
//             ngbDD += (ngbCorrVec & gradDDI[curNeighbour]);

//             ownGradDD[faceI] = gradDDfI[curFace];
//             ownGradDD[faceI] +=
//                 ownN*(interDD[faceI] - ownDD)/ownDeltaN
//               - ownN*(ownN & ownGradDD[faceI]);

//             ngbGradDD[faceI] = gradDDfI[curFace];
//             ngbGradDD[faceI] +=
//                 ngbN*(interDD[faceI] - ngbDD)/ngbDeltaN
//               - ngbN*(ngbN & ngbGradDD[faceI]);
//         }
//         else
//         {
//             label curPatch = mesh().boundaryMesh().whichPatch(curFace);
//             label curPatchFace =
//                 curFace - mesh().boundaryMesh()[curPatch].start();

//             const unallocLabelList& faceCells =
//                 mesh().boundary()[curPatch].faceCells();

//             label curOwner = faceCells[curPatchFace];

//             vector ownN =
//                 mesh().Sf().boundaryField()[curPatch][curPatchFace]
//                /mesh().magSf().boundaryField()[curPatch][curPatchFace];
//             vector ngbN = -ownN;

//             scalar ngbDeltaN =
//                 mesh().weights().boundaryField()[curPatch][curPatchFace]
//                *(
//                    1.0
//                   /mesh().deltaCoeffs().boundaryField()[curPatch][curPatchFace]
//                 );
//             scalar ownDeltaN =
//                 (
//                     1.0
//                    /mesh().deltaCoeffs()
//                    .boundaryField()[curPatch][curPatchFace]
//                 )
//               - ngbDeltaN;

//             vector ownDD = DD_.internalField()[curOwner];
//             vector ngbDD = DD_.boundaryField()[curPatch][curPatchFace];

//             vector ownCorrVec =
//                 mesh().Cf().boundaryField()[curPatch][curPatchFace]
//               - mesh().C().internalField()[curOwner];
//             ownCorrVec -= ownN*(ownN & ownCorrVec);

//             vector ngbCorrVec =
//                 mesh().Cf().boundaryField()[curPatch][curPatchFace]
//               - mesh().C().boundaryField()[curPatch][curPatchFace];
//             ngbCorrVec -= ngbN*(ngbN & ngbCorrVec);

//             vector ownDDCorr =
//                 (
//                     ownCorrVec
//                   & gradDD.internalField()[curOwner]
//                 );
//             vector ngbDDCorr =
//                 (
//                     ngbCorrVec
//                   & gradDD.boundaryField()[curPatch][curPatchFace]
//                 );

//             ownDD += ownDDCorr;
//             ngbDD += ngbDDCorr;

//             ownGradDD[faceI] = gradDDf.boundaryField()[curPatch][curPatchFace];
//             ownGradDD[faceI] +=
//                 ownN*(interDD[faceI] - ownDD)/ownDeltaN
//               - ownN*(ownN & ownGradDD[faceI]);

//             ngbGradDD[faceI] = gradDDf.boundaryField()[curPatch][curPatchFace];
//             ngbGradDD[faceI] +=
//                 ngbN*(interDD[faceI] - ngbDD)/ngbDeltaN
//               - ngbN*(ngbN & ngbGradDD[faceI]);
//         }
//     }

    displacement() += displacementIncrement();

//     ownRelF() = I + ownGradDD.T();
//     ngbRelF() = I + ngbGradDD.T();

    scalarField ownRelJ = det(ownRelF());
    scalarField ngbRelJ = det(ngbRelF());

    ownF() = (ownRelF() & ownF());
    ngbF() = (ngbRelF() & ngbF());

    tensorField ownRelFbar = pow(ownRelJ, -1.0/3.0)*ownRelF();
    tensorField ngbRelFbar = pow(ngbRelJ, -1.0/3.0)*ngbRelF();

    ownBbar() = symm(ownRelFbar & ownBbar() & ownRelFbar);
    ngbBbar() = symm(ngbRelFbar & ngbBbar() & ngbRelFbar);


    // Looking up solid solver
    const solidSolver& solid =
        mesh().lookupObject<solidSolver>
        (
            "solidProperties"
        );

    // If plasticity active
    if (solid.rheology().plasticityActive())
    {
        ownSigmaY() += ownDSigmaY();
        ngbSigmaY() += ngbDSigmaY();

        ownEpsilonPEq() += ownDEpsilonPEq();
        ngbEpsilonPEq() += ngbDEpsilonPEq();

        ownBbar() -= 2*ownDEpsilonP();
        ngbBbar() -= 2*ngbDEpsilonP();
    }
}


void ULLSMaterialInterface::updateInterfaceDisplacementIncrement()
{
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const vectorField& DDI = DD_.internalField();
    const volVectorField& prevDD = DD_.prevIter();
    const vectorField& prevDDI = prevDD.internalField();

    const volTensorField& gradDD =
        mesh().lookupObject<volTensorField>("grad(" + DD_.name() + ')');
    const tensorField& gradDDI = gradDD.internalField();

    const surfaceTensorField& gradDDf =
        mesh().lookupObject<surfaceTensorField>("grad" + DD_.name() + 'f');
//     const tensorField& gradDDfI = gradDDf.internalField();

    const volScalarField& mu =
        mesh().lookupObject<volScalarField>("mu");
    const volScalarField& lambda =
        mesh().lookupObject<volScalarField>("lambda");

    const vectorField& C  = mesh().C().internalField();
    const vectorField& Cf  = mesh().Cf().internalField();

    const vectorField& SI = mesh().Sf().internalField();
    const scalarField& magSI  = mesh().magSf().internalField();
    const scalarField& deltaCoeffs = mesh().deltaCoeffs().internalField();
    const scalarField& w = mesh().weights().internalField();

    vectorField& interDD = displacementIncrement();
    vectorField& interTraction = traction();

    // Interface faces
    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        // Internal faces
        if (curFace < mesh().nInternalFaces())
        {
            label curOwner = owner[curFace];
            label curNeighbour = neighbour[curFace];

            vector ownDD = DDI[curOwner];
            vector ngbDD = DDI[curNeighbour];

            vector ownPrevDD = prevDDI[curOwner];
            vector ngbPrevDD = prevDDI[curNeighbour];

            // Calculate initial interface displacement increment
            // using extrapolation

//             vector ownDelta = Cf[curFace] - C[curOwner];
//             vector ngbDelta = Cf[curFace] - C[curNeighbour];

//             interDD[faceI] =
//             (
//                 ownDD + (ownDelta & gradDDI[curOwner])
//               + ngbDD + (ngbDelta & gradDDI[curNeighbour])
//             )/2;


            // Calculate interface displacement using equilibrium equation
            vector ownN = SI[curFace]/magSI[curFace];
            vector ngbN = -ownN;

            scalar ngbDeltaN = w[curFace]*(1.0/deltaCoeffs[curFace]);
            scalar ownDeltaN = (1.0/deltaCoeffs[curFace]) - ngbDeltaN;

            vector ownCorrVec = Cf[curFace] - C[curOwner];
            ownCorrVec -= ownN*(ownN&ownCorrVec);
            vector ngbCorrVec = Cf[curFace] - C[curNeighbour];
            ngbCorrVec -= ngbN*(ngbN&ngbCorrVec);

            ownDD += (ownCorrVec & gradDDI[curOwner]);
            ngbDD += (ngbCorrVec & gradDDI[curNeighbour]);

            ownPrevDD += (ownCorrVec & gradDDI[curOwner]);
            ngbPrevDD += (ngbCorrVec & gradDDI[curNeighbour]);

            scalar ownMu  = mu.internalField()[curOwner];
            scalar ngbMu  = mu.internalField()[curNeighbour];

            scalar ownLambda  = lambda.internalField()[curOwner];
            scalar ngbLambda  = lambda.internalField()[curNeighbour];

            vectorField curCauchyTraction =
                cauchyTraction
                (
                    interDD[faceI],
                    faceI,
                    gradDD,
                    gradDDf,
                    mu,
                    lambda
                );

            vector ngbNonLinTrac =
                curCauchyTraction[1]
              - (2*ngbMu + ngbLambda)
               *(ngbPrevDD - interDD[faceI])/ngbDeltaN;

            vector ownNonLinTrac =
                curCauchyTraction[0]
              - (2*ownMu + ownLambda)
               *(interDD[faceI] - ownPrevDD)/ownDeltaN;

            interDD[faceI] =
                (
                    (2*ownMu + ownLambda)*ownDD*ngbDeltaN
                  + (2*ngbMu + ngbLambda)*ngbDD*ownDeltaN
                  + ownDeltaN*ngbDeltaN*(ngbNonLinTrac - ownNonLinTrac)
                )
               /(
                    (2*ownMu + ownLambda)*ngbDeltaN
                  + (2*ngbMu + ngbLambda)*ownDeltaN
                );

            interTraction[faceI] = curCauchyTraction[0];

//             vectorField curCauchyTraction_ =
//                 cauchyTraction
//                 (
//                     interDD[faceI],
//                     faceI,
//                     gradDD,
//                     gradDDf,
//                     mu,
//                     lambda
//                 );
//             Info << curCauchyTraction_[0] << ", "
//                 << curCauchyTraction_[1] << endl;
        }
        else
        {
            label curPatch = mesh().boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh().boundaryMesh()[curPatch].start();

            const unallocLabelList& faceCells =
                mesh().boundary()[curPatch].faceCells();

            label curOwner = faceCells[curPatchFace];

            vector ownN =
                mesh().Sf().boundaryField()[curPatch][curPatchFace]
               /mesh().magSf().boundaryField()[curPatch][curPatchFace];
            vector ngbN = -ownN;

            scalar ngbDeltaN =
                mesh().weights().boundaryField()[curPatch][curPatchFace]
               *(
                   1.0
                  /mesh().deltaCoeffs().boundaryField()[curPatch][curPatchFace]
                );
            scalar ownDeltaN =
                (
                    1.0
                   /mesh().deltaCoeffs()
                   .boundaryField()[curPatch][curPatchFace]
                )
              - ngbDeltaN;

            vector ownDD = DD_.internalField()[curOwner];
            vector ngbDD = DD_.boundaryField()[curPatch][curPatchFace];

            vector ownPrevDD = prevDD.internalField()[curOwner];
            vector ngbPrevDD = prevDD.boundaryField()[curPatch][curPatchFace];

            vector ownCorrVec =
                mesh().Cf().boundaryField()[curPatch][curPatchFace]
              - mesh().C().internalField()[curOwner];
            ownCorrVec -= ownN*(ownN & ownCorrVec);

            vector ngbCorrVec =
                mesh().Cf().boundaryField()[curPatch][curPatchFace]
              - mesh().C().boundaryField()[curPatch][curPatchFace];
            ngbCorrVec -= ngbN*(ngbN & ngbCorrVec);

            vector ownDDCorr =
                (
                    ownCorrVec
                  & gradDD.internalField()[curOwner]
                );
            vector ngbDDCorr =
                (
                    ngbCorrVec
                  & gradDD.boundaryField()[curPatch][curPatchFace]
                );

            ownDD += ownDDCorr;
            ngbDD += ngbDDCorr;

            ownPrevDD += ownDDCorr;
            ngbPrevDD += ngbDDCorr;

            scalar ownMu  = mu.internalField()[curOwner];
            scalar ngbMu  = mu.boundaryField()[curPatch][curPatchFace];

            scalar ownLambda =
                lambda.internalField()[curOwner];
            scalar ngbLambda = lambda.boundaryField()[curPatch][curPatchFace];

            vectorField curCauchyTraction =
                cauchyTraction
                (
                    interDD[faceI],
                    faceI,
                    gradDD,
                    gradDDf,
                    mu,
                    lambda
                );

            vector ngbNonLinTrac =
                curCauchyTraction[1]
              - (2*ngbMu + ngbLambda)
               *(ngbPrevDD - interDD[faceI])/ngbDeltaN;

            vector ownNonLinTrac =
                curCauchyTraction[0]
              - (2*ownMu + ownLambda)
               *(interDD[faceI] - ownPrevDD)/ownDeltaN;

            interDD[faceI] =
                (
                    (2*ownMu + ownLambda)*ownDD*ngbDeltaN
                  + (2*ngbMu + ngbLambda)*ngbDD*ownDeltaN
                  + ownDeltaN*ngbDeltaN*(ngbNonLinTrac - ownNonLinTrac)
                )
               /(
                    (2*ownMu + ownLambda)*ngbDeltaN
                  + (2*ngbMu + ngbLambda)*ownDeltaN
                );

            interTraction[faceI] = curCauchyTraction[0];
        }
    }
}


void ULLSMaterialInterface::correctInterfaceGradient
(
    surfaceTensorField& gradDDf
)
{
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const vectorField& DDI = DD_.internalField();

    const volTensorField& gradDD =
        mesh().lookupObject<volTensorField>("grad(" + DD_.name() + ')');
    const tensorField& gradDDI = gradDD.internalField();

//     const volScalarField& mu =
//         mesh().lookupObject<volScalarField>("mu");
//     const volScalarField& lambda =
//         mesh().lookupObject<volScalarField>("lambda");

    const vectorField& SI  = mesh().Sf().internalField();
    const scalarField& magSI  = mesh().magSf().internalField();
    const scalarField& deltaCoeffs = mesh().deltaCoeffs().internalField();
    const scalarField& w = mesh().weights().internalField();

    const vectorField& C  = mesh().C().internalField();
    const vectorField& Cf  = mesh().Cf().internalField();

    const vectorField& interDD = displacementIncrement();

    tensorField& gradDDfI = gradDDf.internalField();

    // Interface faces
    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        // Internal faces
        if (curFace < mesh().nInternalFaces())
        {
            label curOwner = owner[curFace];
            label curNeighbour = neighbour[curFace];

            vector ownN = SI[curFace]/magSI[curFace];
            vector ngbN = -ownN;

            scalar ngbDeltaN = w[curFace]*(1.0/deltaCoeffs[curFace]);
            scalar ownDeltaN = (1.0/deltaCoeffs[curFace]) - ngbDeltaN;

//             scalar ownMu = mu.internalField()[curOwner];
//             scalar ngbMu = mu.internalField()[curNeighbour];

//             scalar ownLambda = lambda.internalField()[curOwner];
//             scalar ngbLambda = lambda.internalField()[curNeighbour];

//             scalar ownK  = (ownLambda + (2.0/3.0)*ownMu);
//             scalar ngbK  = (ngbLambda + (2.0/3.0)*ngbMu);

            vector ownDD = DDI[curOwner];
            vector ngbDD = DDI[curNeighbour];

            vector ownCorrVec = Cf[curFace] - C[curOwner];
            ownCorrVec -= ownN*(ownN & ownCorrVec);

            vector ngbCorrVec = Cf[curFace] - C[curNeighbour];
            ngbCorrVec -= ngbN*(ngbN&ngbCorrVec);

            vector ownDDCorr = (ownCorrVec & gradDDI[curOwner]);
            vector ngbDDCorr = (ngbCorrVec & gradDDI[curNeighbour]);

            ownDD += ownDDCorr;
            ngbDD += ngbDDCorr;

            tensor ownGradDD = gradDDfI[curFace];
            ownGradDD +=
                ownN*(interDD[faceI] - ownDD)/ownDeltaN
              - ownN*(ownN & ownGradDD);

            tensor ngbGradDD = gradDDfI[curFace];
            ngbGradDD +=
                ngbN*(interDD[faceI] - ngbDD)/ngbDeltaN
              - ngbN*(ngbN & ngbGradDD);


            gradDDfI[curFace] = ownGradDD;
//             gradDDfI[curFace] -= ownN*(ownN & gradDDfI[curFace]);
//             gradDDfI[curFace] +=
//                 ownN*(interDD[faceI] - (ownDD + ownDDCorr))/ownDeltaN;

            // Update relative deformation gradient at the interface
            ownRelF()[faceI] = I + ownGradDD.T();
            ngbRelF()[faceI] = I + ngbGradDD.T();

//             if (ngbK > ownK)
//             {
//                 gradDDfI[curFace] -= ngbN*(ngbN & gradDDfI[curFace]);
//                 gradDDfI[curFace] +=
//                     ngbN*(interDD[faceI] - (ngbDD + ngbDDCorr))/ngbDeltaN;
//             }
//             else
//             {
//                 gradDDfI[curFace] -= ownN*(ownN & gradDDfI[curFace]);
//                 gradDDfI[curFace] +=
//                     ownN*(interDD[faceI] - (ownDD + ownDDCorr))/ownDeltaN;
//             }
        }
        else
        {
            label curPatch = mesh().boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh().boundaryMesh()[curPatch].start();

            const unallocLabelList& faceCells =
                mesh().boundary()[curPatch].faceCells();

            label curOwner = faceCells[curPatchFace];
//             label curNeighbour = neighbour[curFace];

            vector ownN =
                mesh().Sf().boundaryField()[curPatch][curPatchFace]
               /mesh().magSf().boundaryField()[curPatch][curPatchFace];
            vector ngbN = -ownN;

            scalar ngbDeltaN =
                mesh().weights().boundaryField()[curPatch][curPatchFace]
               *(
                   1.0
                  /mesh().deltaCoeffs().boundaryField()[curPatch][curPatchFace]
                );
            scalar ownDeltaN =
                (
                    1.0
                   /mesh().deltaCoeffs()
                   .boundaryField()[curPatch][curPatchFace]
                )
              - ngbDeltaN;


//             scalar ownMu  = mu.internalField()[curOwner];
//             scalar ngbMu  = mu.boundaryField()[curPatch][curPatchFace];

//             scalar ownLambda =
//                 lambda.internalField()[curOwner];
//             scalar ngbLambda = lambda.boundaryField()[curPatch][curPatchFace];

//             scalar ownK  = (ownLambda + (2.0/3.0)*ownMu);
//             scalar ngbK  = (ngbLambda + (2.0/3.0)*ngbMu);

            vector ownDD = DD_.internalField()[curOwner];
            vector ngbDD = DD_.boundaryField()[curPatch][curPatchFace];

            vector ownCorrVec =
                mesh().Cf().boundaryField()[curPatch][curPatchFace]
              - mesh().C().internalField()[curOwner];
            ownCorrVec -= ownN*(ownN & ownCorrVec);

            vector ngbCorrVec =
                mesh().Cf().boundaryField()[curPatch][curPatchFace]
              - mesh().C().boundaryField()[curPatch][curPatchFace];
            ngbCorrVec -= ngbN*(ngbN & ngbCorrVec);

            vector ownDDCorr =
                (
                    ownCorrVec
                  & gradDD.internalField()[curOwner]
                );
            vector ngbDDCorr =
                (
                    ngbCorrVec
                  & gradDD.boundaryField()[curPatch][curPatchFace]
                );

            ownDD += ownDDCorr;
            ngbDD += ngbDDCorr;

            tensor ownGradDD = gradDDf.boundaryField()[curPatch][curPatchFace];
            ownGradDD +=
                ownN*(interDD[faceI] - ownDD)/ownDeltaN
              - ownN*(ownN & ownGradDD);

            tensor ngbGradDD = gradDDf.boundaryField()[curPatch][curPatchFace];
            ngbGradDD +=
                ngbN*(interDD[faceI] - ngbDD)/ngbDeltaN
              - ngbN*(ngbN & ngbGradDD);

//             gradDDf.boundaryField()[curPatch][curPatchFace] -=
//                 ownN
//                *(
//                     ownN & gradDDf.boundaryField()[curPatch][curPatchFace]
//                 );
//             gradDDf.boundaryField()[curPatch][curPatchFace] +=
//                 ownN*(interDD[faceI] - (ownDD + ownDDCorr))/ownDeltaN;

            const processorPolyPatch & procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh().boundaryMesh()[curPatch]
                );

            if (procPatch.neighbour())
//             if (ngbK > ownK)
            {
                gradDDf.boundaryField()[curPatch][curPatchFace] = ngbGradDD;

//                 gradDDf.boundaryField()[curPatch][curPatchFace] -=
//                     ngbN
//                    *(
//                         ngbN & gradDDf.boundaryField()[curPatch][curPatchFace]
//                     );
//                 gradDDf.boundaryField()[curPatch][curPatchFace] +=
//                     ngbN
//                    *(
//                        interDD[faceI] - (ngbDD + ngbDDCorr)
//                     )/ngbDeltaN;
            }
            else
            {
                gradDDf.boundaryField()[curPatch][curPatchFace] = ownGradDD;

//                 gradDDf.boundaryField()[curPatch][curPatchFace] -=
//                     ownN
//                    *(
//                        ownN & gradDDf.boundaryField()[curPatch][curPatchFace]
//                     );
//                 gradDDf.boundaryField()[curPatch][curPatchFace] +=
//                     ownN*(interDD[faceI] - (ownDD + ownDDCorr))/ownDeltaN;
            }

            // Update relative deformation gradient at the interface
            ownRelF()[faceI] = I + ownGradDD.T();
            ngbRelF()[faceI] = I + ngbGradDD.T();
        }
    }
}


void ULLSMaterialInterface::updateDisplacementIncrement
(
    pointVectorField& pointDD
)
{
    if (debug)
    {
        Info<< "ULLSMaterialInterface::updateDisplacementIncrement("
            << "pointVectorField&)"
            << "interpolating displacement incr field from cells to points"
            << endl;
    }

    updateInterfaceDisplacementIncrement();

    subMeshPointDD_.clear();

    materialInterface::updateDisplacement
    (
        DD_,
        displacementIncrement(),
        pointDD,
        subMeshDD(),
        subMeshPointDD()
    );
}


void ULLSMaterialInterface::updateDisplacementIncrementGradient
(
    volTensorField& gradDD,
    surfaceTensorField& gradDDf
)
{
    materialInterface::updateDisplacementGradient
    (
        DD_,
        subMeshDD(),
        subMeshPointDD(),
        gradDD,
        gradDDf
    );

    correctInterfaceGradient(gradDDf);
}


void ULLSMaterialInterface::modifyProperty
(
    surfaceScalarField& muf,
    const volScalarField& mu
) const
{
    Info << "ULLSMaterialInterface::modifyProperty" << endl;

    const unallocLabelList& owner = mesh().owner();
//     const unallocLabelList& neighbour = mesh().neighbour();

//     const volScalarField& Mu =
//         mesh().lookupObject<volScalarField>("mu");
//     const volScalarField& Lambda =
//         mesh().lookupObject<volScalarField>("lambda");

    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh().nInternalFaces())
        {
//             scalar ownK =
//                 Lambda[owner[curFace]]
//               + (2.0/3.0)*Mu[owner[curFace]];
//             scalar ngbK =
//                 Lambda[neighbour[curFace]]
//               + (2.0/3.0)*Mu[neighbour[curFace]];

//             if (ngbK > ownK)
//             {
//                 muf.internalField()[curFace] = mu[neighbour[curFace]];
//             }
//             else
//             {
//                 muf.internalField()[curFace] = mu[owner[curFace]];
//             }

            muf.internalField()[curFace] = mu[owner[curFace]];
        }
        else
        {
            label curPatch = mesh().boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh().boundaryMesh()[curPatch].start();

            const unallocLabelList& faceCells =
                mesh().boundary()[curPatch].faceCells();

//             scalar ownK =
//                 Lambda[faceCells[curPatchFace]]
//               + (2.0/3.0)*Mu[faceCells[curPatchFace]];
//             scalar ngbK =
//                 Lambda.boundaryField()[curPatch][curPatchFace]
//               + (2.0/3.0)*Mu.boundaryField()[curPatch][curPatchFace];


            const processorPolyPatch & procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh().boundaryMesh()[curPatch]
                );

            if (procPatch.neighbour())
//             if (ngbK > ownK)
            {
                muf.boundaryField()[curPatch][curPatchFace] =
                    mu.boundaryField()[curPatch][curPatchFace];
            }
            else
            {
                muf.boundaryField()[curPatch][curPatchFace] =
                    mu[faceCells[curPatchFace]];
            }

//             muf.boundaryField()[curPatch][curPatchFace] =
//                 mu[faceCells[curPatchFace]];

//             muf.boundaryField()[curPatch][curPatchFace] = 0;
        }
    }
}


bool ULLSMaterialInterface::updateMesh(const mapPolyMesh& map)
{
    labelList oldFaces = faces();

    subMeshDD_.clear();
    subMeshPointDD_.clear();

    materialInterface::updateMesh(map);

    const labelList& faceMap = map.faceMap();

    vectorField oldDisplacementIncrement = displacementIncrement();
    displacementIncrement().setSize(faces().size(), vector::zero);
    vectorField oldDisplacement = displacement();
    displacement().setSize(faces().size(), vector::zero);

//     vectorField oldTractionIncrement = tractionIncrement();
//     tractionIncrement().setSize(faces().size(), vector::zero);
    vectorField oldTraction = traction();
    traction().setSize(faces().size(), vector::zero);

    tensorField oldOwnRelF = ownRelF();
    tensorField oldNgbRelF = ngbRelF();
    ownRelF().setSize(faces().size(), I);
    ngbRelF().setSize(faces().size(), I);

    tensorField oldOwnF = ownF();
    tensorField oldNgbF = ngbF();
    ownF().setSize(faces().size(), I);
    ngbF().setSize(faces().size(), I);

    symmTensorField oldOwnBbar = ownBbar();
    symmTensorField oldNgbBbar = ngbBbar();
    ownBbar().setSize(faces().size(), I);
    ngbBbar().setSize(faces().size(), I);

    symmTensorField oldOwnDEpsilonP = ownDEpsilonP();
    symmTensorField oldNgbDEpsilonP = ngbDEpsilonP();
    ownDEpsilonP().setSize(faces().size(), I);
    ngbDEpsilonP().setSize(faces().size(), I);

    scalarField oldOwnSigmaY = ownSigmaY();
    scalarField oldNgbSigmaY = ngbSigmaY();
    ownSigmaY().setSize(faces().size(), 0);
    ngbSigmaY().setSize(faces().size(), 0);

    scalarField oldOwnDSigmaY = ownDSigmaY();
    scalarField oldNgbDSigmaY = ngbDSigmaY();
    ownDSigmaY().setSize(faces().size(), 0);
    ngbDSigmaY().setSize(faces().size(), 0);

    scalarField oldOwnEpsilonPEq = ownEpsilonPEq();
    scalarField oldNgbEpsilonPEq = ngbEpsilonPEq();
    ownEpsilonPEq().setSize(faces().size(), 0);
    ngbEpsilonPEq().setSize(faces().size(), 0);

    scalarField oldOwnDEpsilonPEq = ownDEpsilonPEq();
    scalarField oldNgbDEpsilonPEq = ngbDEpsilonPEq();
    ownDEpsilonPEq().setSize(faces().size(), 0);
    ngbDEpsilonPEq().setSize(faces().size(), 0);

    // Mam local fields
    forAll(faces(), faceI)
    {
        label oldFace = faceMap[faces()[faceI]];
        label oldInterFace = findIndex(oldFaces, oldFace);

        if (oldInterFace != -1)
        {
            // Map previous values
            displacementIncrement()[faceI] =
                oldDisplacementIncrement[oldInterFace];
            displacement()[faceI] = oldDisplacement[oldInterFace];

//             tractionIncrement()[faceI] = oldTractionIncrement[oldInterFace];
            traction()[faceI] = oldTraction[oldInterFace];

            ownRelF()[faceI] = oldOwnRelF[oldInterFace];
            ngbRelF()[faceI] = oldNgbRelF[oldInterFace];

            ownF()[faceI] = oldOwnF[oldInterFace];
            ngbF()[faceI] = oldNgbF[oldInterFace];

            ownBbar()[faceI] = oldOwnBbar[oldInterFace];
            ngbBbar()[faceI] = oldNgbBbar[oldInterFace];

            ownDEpsilonP()[faceI] = oldOwnDEpsilonP[oldInterFace];
            ngbDEpsilonP()[faceI] = oldOwnDEpsilonP[oldInterFace];

            ownSigmaY()[faceI] = oldOwnSigmaY[oldInterFace];
            ngbSigmaY()[faceI] = oldNgbSigmaY[oldInterFace];

            ownDSigmaY()[faceI] = oldOwnDSigmaY[oldInterFace];
            ngbDSigmaY()[faceI] = oldNgbDSigmaY[oldInterFace];

            ownEpsilonPEq()[faceI] = oldOwnEpsilonPEq[oldInterFace];
            ngbEpsilonPEq()[faceI] = oldNgbEpsilonPEq[oldInterFace];

            ownDEpsilonPEq()[faceI] = oldOwnDEpsilonPEq[oldInterFace];
            ngbDEpsilonPEq()[faceI] = oldNgbDEpsilonPEq[oldInterFace];
        }
        else
        {
            // No new faces
        }
    }

    return true;
}


bool ULLSMaterialInterface::update()
{
    labelList oldFaces = faces();

    subMeshDD_.clear();
    subMeshPointDD_.clear();

    materialInterface::update();

//     const labelList& faceMap = map.faceMap();

    labelList faceMap(mesh().nFaces(), -1);
    forAll(faceMap, faceI)
    {
        faceMap[faceI] = faceI;
    }

    vectorField oldDisplacementIncrement = displacementIncrement();
    displacementIncrement().setSize(faces().size(), vector::zero);
    vectorField oldDisplacement = displacement();
    displacement().setSize(faces().size(), vector::zero);

//     vectorField oldTractionIncrement = tractionIncrement();
//     tractionIncrement().setSize(faces().size(), vector::zero);
    vectorField oldTraction = traction();
    traction().setSize(faces().size(), vector::zero);

    tensorField oldOwnRelF = ownRelF();
    tensorField oldNgbRelF = ngbRelF();
    ownRelF().setSize(faces().size(), I);
    ngbRelF().setSize(faces().size(), I);

    tensorField oldOwnF = ownF();
    tensorField oldNgbF = ngbF();
    ownF().setSize(faces().size(), I);
    ngbF().setSize(faces().size(), I);

    symmTensorField oldOwnBbar = ownBbar();
    symmTensorField oldNgbBbar = ngbBbar();
    ownBbar().setSize(faces().size(), I);
    ngbBbar().setSize(faces().size(), I);

    symmTensorField oldOwnDEpsilonP = ownDEpsilonP();
    symmTensorField oldNgbDEpsilonP = ngbDEpsilonP();
    ownDEpsilonP().setSize(faces().size(), I);
    ngbDEpsilonP().setSize(faces().size(), I);

    scalarField oldOwnSigmaY = ownSigmaY();
    scalarField oldNgbSigmaY = ngbSigmaY();
    ownSigmaY().setSize(faces().size(), 0);
    ngbSigmaY().setSize(faces().size(), 0);

    scalarField oldOwnDSigmaY = ownDSigmaY();
    scalarField oldNgbDSigmaY = ngbDSigmaY();
    ownDSigmaY().setSize(faces().size(), 0);
    ngbDSigmaY().setSize(faces().size(), 0);

    scalarField oldOwnEpsilonPEq = ownEpsilonPEq();
    scalarField oldNgbEpsilonPEq = ngbEpsilonPEq();
    ownEpsilonPEq().setSize(faces().size(), 0);
    ngbEpsilonPEq().setSize(faces().size(), 0);

    scalarField oldOwnDEpsilonPEq = ownDEpsilonPEq();
    scalarField oldNgbDEpsilonPEq = ngbDEpsilonPEq();
    ownDEpsilonPEq().setSize(faces().size(), 0);
    ngbDEpsilonPEq().setSize(faces().size(), 0);


    // Initialize displacement
    surfaceVectorField DDf = fvc::interpolate(DD_);
    surfaceVectorField Df = fvc::interpolate(D_);

    const surfaceTensorField& gradDDf =
        mesh().lookupObject<surfaceTensorField>("gradDDf");

    const surfaceTensorField& Ff =
        mesh().lookupObject<surfaceTensorField>("Ff");

    const surfaceSymmTensorField& bBarf =
        mesh().lookupObject<surfaceSymmTensorField>("bBarf");

    // Looking up solid solver
    const solidSolver& solid =
        mesh().lookupObject<solidSolver>
        (
            "solidProperties"
        );


    // Mam local fields
    forAll(faces(), faceI)
    {
        label oldFace = faceMap[faces()[faceI]];
        label oldInterFace = findIndex(oldFaces, oldFace);

        if (oldInterFace != -1)
        {
            // Map previous values
            displacementIncrement()[faceI] =
                oldDisplacementIncrement[oldInterFace];
            displacement()[faceI] = oldDisplacement[oldInterFace];

//             tractionIncrement()[faceI] = oldTractionIncrement[oldInterFace];
            traction()[faceI] = oldTraction[oldInterFace];

            ownRelF()[faceI] = oldOwnRelF[oldInterFace];
            ngbRelF()[faceI] = oldNgbRelF[oldInterFace];

            ownF()[faceI] = oldOwnF[oldInterFace];
            ngbF()[faceI] = oldNgbF[oldInterFace];

            ownBbar()[faceI] = oldOwnBbar[oldInterFace];
            ngbBbar()[faceI] = oldNgbBbar[oldInterFace];

            ownDEpsilonP()[faceI] = oldOwnDEpsilonP[oldInterFace];
            ngbDEpsilonP()[faceI] = oldOwnDEpsilonP[oldInterFace];

            ownSigmaY()[faceI] = oldOwnSigmaY[oldInterFace];
            ngbSigmaY()[faceI] = oldNgbSigmaY[oldInterFace];

            ownDSigmaY()[faceI] = oldOwnDSigmaY[oldInterFace];
            ngbDSigmaY()[faceI] = oldNgbDSigmaY[oldInterFace];

            ownEpsilonPEq()[faceI] = oldOwnEpsilonPEq[oldInterFace];
            ngbEpsilonPEq()[faceI] = oldNgbEpsilonPEq[oldInterFace];

            ownDEpsilonPEq()[faceI] = oldOwnDEpsilonPEq[oldInterFace];
            ngbDEpsilonPEq()[faceI] = oldNgbDEpsilonPEq[oldInterFace];
        }
        else
        {
            // New faces
            displacementIncrement()[faceI] = DDf[faces()[faceI]];
            displacement()[faceI] = Df[faces()[faceI]];

//             tractionIncrement()[faceI];
//             traction()[faceI];

            ownRelF()[faceI] = I + gradDDf[faces()[faceI]].T();
            ngbRelF()[faceI] = I + gradDDf[faces()[faceI]].T();

            ownF()[faceI] = Ff[faces()[faceI]];
            ngbF()[faceI] = Ff[faces()[faceI]];

            ownBbar()[faceI] = bBarf[faces()[faceI]];
            ngbBbar()[faceI] = bBarf[faces()[faceI]];

            // If plasticity active
            if (solid.rheology().plasticityActive())
            {
                ownDEpsilonP()[faceI];
                ngbDEpsilonP()[faceI];

                ownSigmaY()[faceI];
                ngbSigmaY()[faceI];

                ownDSigmaY()[faceI];
                ngbDSigmaY()[faceI];

                ownEpsilonPEq()[faceI];
                ngbEpsilonPEq()[faceI];

                ownDEpsilonPEq()[faceI];
                ngbDEpsilonPEq()[faceI];
            }
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

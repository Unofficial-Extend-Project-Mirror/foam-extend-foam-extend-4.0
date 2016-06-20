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

#include "ITLMaterialInterface.H"
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ITLMaterialInterface, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void ITLMaterialInterface::makeDisplacementIncrement() const
{
    if (debug)
    {
        Info<< "void ITLMaterialInterface::"
            << "makeInterfaceDisplacementIncrement() const : "
            << "creating interface displacement field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (displacementIncrementPtr_)
    {
        FatalErrorIn("ITLMaterialInterface::makeDisplacementIncrement() const")
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


void ITLMaterialInterface::makeDisplacement() const
{
    if (debug)
    {
        Info<< "void ITLMaterialInterface"
            << "::makeInterfaceDisplacement() const : "
            << "creating interface displacement field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (displacementPtr_)
    {
        FatalErrorIn("ITLMaterialInterface::makeDisplacement() const")
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


void ITLMaterialInterface::makeTractionIncrement() const
{
    if (debug)
    {
        Info<< "void ITLMaterialInterface::makeTractionIncrement() const : "
            << "creating interface traction increment field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (tractionIncrementPtr_)
    {
        FatalErrorIn("ITLMaterialInterface::makeTractionIncrement() const")
            << "interface traction increment field already exist"
            << abort(FatalError);
    }

    tractionIncrementPtr_ =
        new vectorField(faces().size(), vector::zero);
}


void ITLMaterialInterface::makeTraction() const
{
    if (debug)
    {
        Info<< "void ITLMaterialInterface::makeTraction() const : "
            << "creating interface traction field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (tractionPtr_)
    {
        FatalErrorIn("ITLMaterialInterface::makeTraction() const")
            << "interface traction field already exist"
            << abort(FatalError);
    }

    tractionPtr_ = new vectorField(faces().size(), vector::zero);
}


void ITLMaterialInterface::makeDSigmaf() const
{
    if (debug)
    {
        Info<< "void ITLMaterialInterface::makeDSigmaf() const : "
            << "creating interface second Piola-Kirchhoff stress incr field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ngbDSigmafPtr_ || ownDSigmafPtr_)
    {
        FatalErrorIn("ITLMaterialInterface::makeDSigmaf() const")
            << "interface second Piola-Kirchhoff stress incr field "
                << "already exists"
                << abort(FatalError);
    }

    ownDSigmafPtr_ = new symmTensorField(faces().size(), symmTensor::zero);
    ngbDSigmafPtr_ = new symmTensorField(faces().size(), symmTensor::zero);
}


void ITLMaterialInterface::makeSigmaf() const
{
    if (debug)
    {
        Info<< "void ITLMaterialInterface::makeSigmaf() const : "
            << "creating interface old second Piola-Kirchhoff stress field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ngbSigmafPtr_ || ownSigmafPtr_)
    {
        FatalErrorIn("ITLMaterialInterface::makeSigmaf() const")
            << "interface old second Piola-Kirchhoff stress field "
                << "already exists"
                << abort(FatalError);
    }

    ownSigmafPtr_ = new symmTensorField(faces().size(), symmTensor::zero);
    ngbSigmafPtr_ = new symmTensorField(faces().size(), symmTensor::zero);
}


void ITLMaterialInterface::makeSubMeshPointDD() const
{
    if (debug)
    {
        Info<< "void ITLMaterialInterface::makeSubMeshPointDD() const : "
            << "creating point displacements increment fields"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshPointDD_.empty())
    {
        FatalErrorIn("ITLMaterialInterface::makeSubMeshPointDD() const")
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


void ITLMaterialInterface::makeSubMeshDD() const
{
    if (debug)
    {
        Info<< "void ITLMaterialInterface::makeSubMeshDD() const : "
            << "creating displacements fields"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshDD_.empty())
    {
        FatalErrorIn("ITLMaterialInterface::makeSubMeshDD() const")
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


void ITLMaterialInterface::clearOut()
{
    deleteDemandDrivenData(displacementIncrementPtr_);
    deleteDemandDrivenData(displacementPtr_);
    deleteDemandDrivenData(tractionIncrementPtr_);
    deleteDemandDrivenData(tractionPtr_);

    deleteDemandDrivenData(ownDSigmafPtr_);
    deleteDemandDrivenData(ngbDSigmafPtr_);

    deleteDemandDrivenData(ownSigmafPtr_);
    deleteDemandDrivenData(ngbSigmafPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


ITLMaterialInterface::ITLMaterialInterface
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
    tractionIncrementPtr_(NULL),
    tractionPtr_(NULL),
    ownDSigmafPtr_(NULL),
    ngbDSigmafPtr_(NULL),
    ownSigmafPtr_(NULL),
    ngbSigmafPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

ITLMaterialInterface::~ITLMaterialInterface()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

vectorField& ITLMaterialInterface::displacementIncrement()
{
    if (!displacementIncrementPtr_)
    {
        makeDisplacementIncrement();
    }

    return *displacementIncrementPtr_;
}

const vectorField& ITLMaterialInterface::displacementIncrement() const
{
    if (!displacementIncrementPtr_)
    {
        makeDisplacementIncrement();
    }

    return *displacementIncrementPtr_;
}

vectorField& ITLMaterialInterface::displacement()
{
    if (!displacementPtr_)
    {
        makeDisplacement();
    }

    return *displacementPtr_;
}

const vectorField& ITLMaterialInterface::displacement() const
{
    if (!displacementPtr_)
    {
        makeDisplacement();
    }

    return *displacementPtr_;
}

vectorField& ITLMaterialInterface::tractionIncrement()
{
    if (!tractionIncrementPtr_)
    {
        makeTractionIncrement();
    }

    return *tractionIncrementPtr_;
}

const vectorField& ITLMaterialInterface::tractionIncrement() const
{
    if (!tractionIncrementPtr_)
    {
        makeTractionIncrement();
    }

    return *tractionIncrementPtr_;
}

vectorField& ITLMaterialInterface::traction()
{
    if (!tractionPtr_)
    {
        makeTraction();
    }

    return *tractionPtr_;
}

const vectorField& ITLMaterialInterface::traction() const
{
    if (!tractionPtr_)
    {
        makeTraction();
    }

    return *tractionPtr_;
}

symmTensorField& ITLMaterialInterface::ownSigmaf()
{
    if (!ownSigmafPtr_)
    {
        makeSigmaf();
    }

    return *ownSigmafPtr_;
}

const symmTensorField& ITLMaterialInterface::ownSigmaf() const
{
    if (!ownSigmafPtr_)
    {
        makeSigmaf();
    }

    return *ownSigmafPtr_;
}

symmTensorField& ITLMaterialInterface::ngbSigmaf()
{
    if (!ngbSigmafPtr_)
    {
        makeSigmaf();
    }

    return *ngbSigmafPtr_;
}

const symmTensorField& ITLMaterialInterface::ngbSigmaf() const
{
    if (!ngbSigmafPtr_)
    {
        makeSigmaf();
    }

    return *ngbSigmafPtr_;
}

symmTensorField& ITLMaterialInterface::ownDSigmaf()
{
    if (!ownDSigmafPtr_)
    {
        makeDSigmaf();
    }

    return *ownDSigmafPtr_;
}

const symmTensorField& ITLMaterialInterface::ownDSigmaf() const
{
    if (!ownDSigmafPtr_)
    {
        makeDSigmaf();
    }

    return *ownDSigmafPtr_;
}

symmTensorField& ITLMaterialInterface::ngbDSigmaf()
{
    if (!ngbDSigmafPtr_)
    {
        makeDSigmaf();
    }

    return *ngbDSigmafPtr_;
}

const symmTensorField& ITLMaterialInterface::ngbDSigmaf() const
{
    if (!ngbDSigmafPtr_)
    {
        makeDSigmaf();
    }

    return *ngbDSigmafPtr_;
}


const PtrList<volVectorField>& ITLMaterialInterface::subMeshDD() const
{
    if (subMeshDD_.empty())
    {
        makeSubMeshDD();
    }

    return subMeshDD_;
}


PtrList<volVectorField>& ITLMaterialInterface::subMeshDD()
{
    if (subMeshDD_.empty())
    {
        makeSubMeshDD();
    }

    return subMeshDD_;
}


const PtrList<pointVectorField>& ITLMaterialInterface::subMeshPointDD() const
{
    if (subMeshPointDD_.empty())
    {
        makeSubMeshPointDD();
    }

    return subMeshPointDD_;
}


PtrList<pointVectorField>& ITLMaterialInterface::subMeshPointDD()
{
    if (subMeshPointDD_.empty())
    {
        makeSubMeshPointDD();
    }

    return subMeshPointDD_;
}


void ITLMaterialInterface::correct(fvVectorMatrix& DEqn)
{
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const vectorField& DDI = DD_.internalField();

    const vectorField& DI = D_.internalField();

    const volTensorField& gradDD =
        mesh().lookupObject<volTensorField>("grad(" + DD_.name() + ')');
    const tensorField& gradDDI = gradDD.internalField();

    const volTensorField& gradD =
        mesh().lookupObject<volTensorField>("grad(" + D_.name() + ')');
    const tensorField& gradDI = gradD.internalField();

    const surfaceTensorField& gradDDf =
        mesh().lookupObject<surfaceTensorField>("grad" + DD_.name() + 'f');
    const tensorField& gradDDfI = gradDDf.internalField();

    const surfaceTensorField& gradDf =
        mesh().lookupObject<surfaceTensorField>("grad" + D_.name() + 'f');
    const tensorField& gradDfI = gradDf.internalField();

    const volScalarField& mu =
        mesh().lookupObject<volScalarField>("mu");
    const volScalarField& lambda =
        mesh().lookupObject<volScalarField>("lambda");

    const vectorField& SI  = mesh().Sf().internalField();
    const scalarField& magSI  = mesh().magSf().internalField();
    const scalarField& deltaCoeffs = mesh().deltaCoeffs().internalField();
    const scalarField& w = mesh().weights().internalField();

    const vectorField& C  = mesh().C().internalField();
    const vectorField& Cf  = mesh().Cf().internalField();

    scalarField& diag = DEqn.diag();
    scalarField& upper = DEqn.upper();
    vectorField& source = DEqn.source();
    FieldField<Field, vector>& boundaryCoeffs = DEqn.boundaryCoeffs();

    vectorField& interDD = displacementIncrement();
    vectorField& interDT = tractionIncrement();

    if (curTimeIndex_ != mesh().time().timeIndex())
    {
        // Update total fields
        ownSigmaf() += ownDSigmaf();
        ngbSigmaf() += ngbDSigmaf();

        displacement() += displacementIncrement();
        traction() += tractionIncrement();

        curTimeIndex_ = mesh().time().timeIndex();
    }

    const vectorField& interD = displacement();
    const symmTensorField& ownSigmaf_ = ownSigmaf();
    const symmTensorField& ngbSigmaf_ = ngbSigmaf();

    symmTensorField& ownDSigmaf_ = ownDSigmaf();
    symmTensorField& ngbDSigmaf_ = ngbDSigmaf();

    // Looking up solid solver
    const solidSolver& stress =
        mesh().objectRegistry::lookupObject<solidSolver>
        (
            "solidProperties"
        );

    Switch nonLinear
    (
        stress.solidProperties().lookup("nonLinear")
    );

    Switch enforceLinear
    (
        stress.solidProperties().lookup("enforceLinear")
    );

    // Internal faces
    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh().nInternalFaces())
        {
            label curOwner = owner[curFace];
            label curNeighbour = neighbour[curFace];

            vector ownN = SI[curFace]/magSI[curFace];
            vector ngbN = -ownN;

            scalar magS = magSI[curFace];

            tensor ownSGradDD = ((I-ownN*ownN)&gradDDfI[curFace]);
            tensor ngbSGradDD = ((I-ngbN*ngbN)&gradDDfI[curFace]);

            scalar ownTrSGradDDt = tr(ownSGradDD&(I-ownN*ownN));
            scalar ngbTrSGradDDt = tr(ngbSGradDD&(I-ownN*ownN));

            vector ownSGradDDn = (ownSGradDD&ownN);
            vector ngbSGradDDn = (ngbSGradDD&ownN);

            scalar ngbDeltaN = w[curFace]*(1.0/deltaCoeffs[curFace]);
            scalar ownDeltaN = (1.0/deltaCoeffs[curFace]) - ngbDeltaN;

            scalar ownMu  = mu.internalField()[curOwner];
            scalar ngbMu  = mu.internalField()[curNeighbour];

            scalar ownLambda  = lambda.internalField()[curOwner];
            scalar ngbLambda  = lambda.internalField()[curNeighbour];

            vector ownDD = DDI[curOwner];
            vector ngbDD = DDI[curNeighbour];

            vector ownDDt = ((I-ownN*ownN)&DDI[curOwner]);
            vector ngbDDt = ((I-ngbN*ngbN)&DDI[curNeighbour]);

            vector ownDDn = ownN*(ownN&DDI[curOwner]);
            vector ngbDDn = ngbN*(ngbN&DDI[curNeighbour]);

            vector ownCorrVec = Cf[curFace] - C[curOwner];
            ownCorrVec -= ownN*(ownN&ownCorrVec);

            vector ngbCorrVec = Cf[curFace] - C[curNeighbour];
            ngbCorrVec -= ngbN*(ngbN&ngbCorrVec);

            vector ownDDCorr = (ownCorrVec&gradDDI[curOwner]);
            vector ngbDDCorr = (ngbCorrVec&gradDDI[curNeighbour]);

            vector ownDDnCorr = ownN*(ownN&ownDDCorr);
            vector ngbDDnCorr = ngbN*(ngbN&ngbDDCorr);

            vector ownDDtCorr = ((I-ownN*ownN)&ownDDCorr);
            vector ngbDDtCorr = ((I-ngbN*ngbN)&ngbDDCorr);

            vector ownNonLinTracIncr = vector::zero;
            vector ngbNonLinTracIncr = vector::zero;

            if (nonLinear && (!enforceLinear))
            {
                // Gradient of displacement increment at the interface
                tensor ownGradDDf =
                    ownSGradDD
                  + ownN*(interDD[faceI] - (DDI[curOwner] + ownDDCorr))
                   /ownDeltaN;
                tensor ngbGradDDf =
                    ngbSGradDD
                  + ownN*((DDI[curNeighbour] + ngbDDCorr) - interDD[faceI])
                   /ngbDeltaN;

                // Gradient of displacement at the interface
                tensor ownSGradD = ((I-ownN*ownN)&gradDfI[curFace]);
                tensor ngbSGradD = ((I-ngbN*ngbN)&gradDfI[curFace]);

                vector ownDCorr = (ownCorrVec&gradDI[curOwner]);
                vector ngbDCorr = (ngbCorrVec&gradDI[curNeighbour]);

                tensor ownGradDf =
                    ownSGradD
                  + ownN*(interD[faceI] - (DI[curOwner] + ownDCorr))
                   /ownDeltaN;
                tensor ngbGradDf =
                    ngbSGradD
                  + ownN*((DI[curNeighbour] + ngbDCorr) - interD[faceI])
                   /ngbDeltaN;

                // Strain increment at the interface
                symmTensor ownDEf =
                    symm(ownGradDDf)
                  + 0.5*symm(ownGradDDf & ownGradDDf.T())
                  + 0.5*symm(ownGradDf & ownGradDDf.T())
                  + 0.5*symm(ownGradDDf & ownGradDf.T());
                symmTensor ngbDEf =
                    symm(ngbGradDDf)
                  + 0.5*symm(ngbGradDDf & ngbGradDDf.T())
                  + 0.5*symm(ngbGradDf & ngbGradDDf.T())
                  + 0.5*symm(ngbGradDDf & ngbGradDf.T());

                // Stress increment at the interface
                ownDSigmaf_[faceI] =
                    2*ownMu*ownDEf  + I*(ownLambda*tr(ownDEf));
                ngbDSigmaf_[faceI] =
                    2*ngbMu*ngbDEf  + I*(ngbLambda*tr(ngbDEf));


                // Nonlinear part of the traction increment
                tensor ownNonLinDSigmaf =
                    ownMu*(ownGradDDf & ownGradDDf.T())
                  + ownMu*(ownGradDf & ownGradDDf.T())
                  + ownMu*(ownGradDDf & ownGradDf.T())
                  + 0.5*ownLambda*tr(ownGradDDf & ownGradDDf.T())*I
                  + 0.5*ownLambda*tr(ownGradDf & ownGradDDf.T())*I
                  + 0.5*ownLambda*tr(ownGradDDf & ownGradDf.T())*I
                  + (ownDSigmaf_[faceI]&ownGradDf)
                  + ((ownSigmaf_[faceI]+ownDSigmaf_[faceI])&ownGradDDf);

                tensor ngbNonLinDSigmaf =
                    ngbMu*(ngbGradDDf & ngbGradDDf.T())
                  + ngbMu*(ngbGradDf & ngbGradDDf.T())
                  + ngbMu*(ngbGradDDf & ngbGradDf.T())
                  + 0.5*ngbLambda*tr(ngbGradDDf & ngbGradDDf.T())*I
                  + 0.5*ngbLambda*tr(ngbGradDf & ngbGradDDf.T())*I
                  + 0.5*ngbLambda*tr(ngbGradDDf & ngbGradDf.T())*I
                  + (ngbDSigmaf_[faceI]&ngbGradDf)
                  + ((ngbSigmaf_[faceI]+ngbDSigmaf_[faceI])&ngbGradDDf);

                ownNonLinTracIncr = (ownN&ownNonLinDSigmaf);
                ngbNonLinTracIncr = (ownN&ngbNonLinDSigmaf);
            }

            vector ownNonLinTracIncrN = ownN*(ownN&ownNonLinTracIncr);
            vector ngbNonLinTracIncrN = ngbN*(ngbN&ngbNonLinTracIncr);

            vector ownNonLinTracIncrT = ((I-ownN*ownN)&ownNonLinTracIncr);
            vector ngbNonLinTracIncrT = ((I-ngbN*ngbN)&ngbNonLinTracIncr);


            // Interface displacement

            vector curInterDDt =
                (
                    ownMu*(ownDDt+ownDDtCorr)*ngbDeltaN
                  + ngbMu*(ngbDDt+ngbDDtCorr)*ownDeltaN
                  + ownDeltaN*ngbDeltaN*(ngbMu*ngbSGradDDn - ownMu*ownSGradDDn)
                  + ownDeltaN*ngbDeltaN
                   *(ngbNonLinTracIncrT - ownNonLinTracIncrT)
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            vector curInterDDn =
                (
                    (2*ownMu + ownLambda)*(ownDDn+ownDDnCorr)*ngbDeltaN
                  + (2*ngbMu + ngbLambda)*(ngbDDn+ngbDDnCorr)*ownDeltaN
                  + ownDeltaN*ngbDeltaN
                   *(ngbLambda*ngbTrSGradDDt - ownLambda*ownTrSGradDDt)*ownN
                  + ownDeltaN*ngbDeltaN
                   *(ngbNonLinTracIncrN - ownNonLinTracIncrN)
                )
               /(
                   (2*ownMu + ownLambda)*ngbDeltaN
                 + (2*ngbMu + ngbLambda)*ownDeltaN
                );

            interDD[faceI] = curInterDDn + curInterDDt;


            // Implicit coupling

            scalar wRevLin = 1.0 - w[curFace];

            scalar ownK = (2*ownMu + ownLambda);
            scalar ngbK = (2*ngbMu + ngbLambda);

            scalar Kf = 1.0/(wRevLin/ownK + (1.0-wRevLin)/ngbK);
            scalar muf = 1.0/(wRevLin/ownMu + (1.0-wRevLin)/ngbMu);

            scalar DeltaNf = ownDeltaN + ngbDeltaN;

            // Owner
            diag[curOwner] += Kf*magS/DeltaNf;

            upper[curFace] -= Kf*magS/DeltaNf;

            source[curOwner] +=
              - (Kf - muf)
               *((DDI[curNeighbour] - DDI[curOwner])&(I-ownN*ownN))
               *magS/DeltaNf;

            source[curOwner] +=
                (
                    ownK*ngbDeltaN*ngbLambda*ngbTrSGradDDt
                  + ngbK*ownDeltaN*ownLambda*ownTrSGradDDt
                )*ownN*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownMu*ngbMu*ngbDeltaN*ngbSGradDDn
                  + ownMu*ngbMu*ownDeltaN*ownSGradDDn
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            source[curOwner] +=
                (
                    ownK*ngbDeltaN*ngbNonLinTracIncrN
                  + ngbK*ownDeltaN*ownNonLinTracIncrN
                )*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbNonLinTracIncrT
                  + ngbMu*ownDeltaN*ownNonLinTracIncrT
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            source[curOwner] +=
                Kf*(ngbDDCorr - ownDDCorr)*magS/DeltaNf
              - (Kf - muf)*(ngbDDtCorr - ownDDtCorr)*magS/DeltaNf;

            // Neighbour
            diag[curNeighbour] += Kf*magS/DeltaNf;

            source[curNeighbour] -=
              - (Kf - muf)
               *((DDI[curNeighbour] - DDI[curOwner])&(I-ownN*ownN))
               *magS/DeltaNf;

            source[curNeighbour] -=
                (
                    ownK*ngbDeltaN*ngbLambda*ngbTrSGradDDt
                  + ngbK*ownDeltaN*ownLambda*ownTrSGradDDt
                )*ownN*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownMu*ngbMu*ngbDeltaN*ngbSGradDDn
                  + ownMu*ngbMu*ownDeltaN*ownSGradDDn
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            source[curNeighbour] -=
                (
                    ownK*ngbDeltaN*ngbNonLinTracIncrN
                  + ngbK*ownDeltaN*ownNonLinTracIncrN
                )*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbNonLinTracIncrT
                  + ngbMu*ownDeltaN*ownNonLinTracIncrT
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            source[curNeighbour] -=
                Kf*(ngbDDCorr - ownDDCorr)*magS/DeltaNf
              - (Kf - muf)*(ngbDDtCorr - ownDDtCorr)*magS/DeltaNf;


            // Interface Cauchy traction increment

            vector curInterDTt =
                muf*((ngbDD+ngbDDCorr) - (ownDD+ownDDCorr))/DeltaNf
              + (
                    ownMu*ngbMu*ngbDeltaN*ngbSGradDDn
                  + ownMu*ngbMu*ownDeltaN*ownSGradDDn
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbNonLinTracIncrT
                  + ngbMu*ownDeltaN*ownNonLinTracIncrT
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            vector curInterDTn =
                Kf*((ngbDDn+ngbDDnCorr) - (ownDDn+ownDDnCorr))/DeltaNf
              + (
                    ownK*ngbDeltaN*ngbLambda*ngbTrSGradDDt
                  + ngbK*ownDeltaN*ownLambda*ownTrSGradDDt
                )*ownN
              + (
                    ownK*ngbDeltaN*ngbNonLinTracIncrN
                  + ngbK*ownDeltaN*ownNonLinTracIncrN
                )
               /(ownK*ngbDeltaN + ngbK*ownDeltaN);

            interDT[faceI] = curInterDTt + curInterDTn;
        }
        else
        {
            // Processor faces
            label curPatch = mesh().boundaryMesh().whichPatch(curFace);

            label curPatchFace =
                curFace - mesh().boundaryMesh()[curPatch].start();

            const labelList& curPatchCells =
                mesh().boundaryMesh()[curPatch].faceCells();

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

            scalar magS =
                mesh().magSf().boundaryField()[curPatch][curPatchFace];

            scalar ownMu  = mu.internalField()[curPatchCells[curPatchFace]];
            scalar ngbMu  = mu.boundaryField()[curPatch][curPatchFace];

            scalar ownLambda =
                lambda.internalField()[curPatchCells[curPatchFace]];
            scalar ngbLambda = lambda.boundaryField()[curPatch][curPatchFace];

            vector ownN =
                mesh().Sf().boundaryField()[curPatch][curPatchFace]
               /mesh().magSf().boundaryField()[curPatch][curPatchFace];
            vector ngbN = -ownN;

            vector ownDD = DD_.internalField()[curPatchCells[curPatchFace]];
            vector ngbDD = DD_.boundaryField()[curPatch][curPatchFace];

            vector ownDDt = ((I-ownN*ownN)&ownDD);
            vector ngbDDt = ((I-ngbN*ngbN)&ngbDD);

            vector ownDDn = ownN*(ownN&ownDD);
            vector ngbDDn = ngbN*(ngbN&ngbDD);

            tensor ownSGradDD =
                (
                    (I-ownN*ownN)
                   &gradDDf.boundaryField()[curPatch][curPatchFace]
                );
            tensor ngbSGradDD = ownSGradDD;

            scalar ownTrSGradDDt = tr(ownSGradDD&(I-ownN*ownN));
            scalar ngbTrSGradDDt = tr(ngbSGradDD&(I-ngbN*ngbN));

            vector ownSGradDDn = (ownSGradDD&ownN);
            vector ngbSGradDDn = (ngbSGradDD&ownN);

            vector ownCorrVec =
                mesh().Cf().boundaryField()[curPatch][curPatchFace]
              - mesh().C().internalField()[curPatchCells[curPatchFace]];
            ownCorrVec -= ownN*(ownN&ownCorrVec);

            vector ngbCorrVec =
                mesh().Cf().boundaryField()[curPatch][curPatchFace]
              - mesh().C().boundaryField()[curPatch][curPatchFace];
            ngbCorrVec -= ngbN*(ngbN&ngbCorrVec);

            vector ownDDCorr =
                (
                    ownCorrVec
                  & gradDD.internalField()[curPatchCells[curPatchFace]]
                );
            vector ngbDDCorr =
                (ngbCorrVec&gradDD.boundaryField()[curPatch][curPatchFace]);

            vector ownDDnCorr = ownN*(ownN&ownDDCorr);
            vector ngbDDnCorr = ownN*(ownN&ngbDDCorr);

            vector ownDDtCorr = ((I-ownN*ownN)&ownDDCorr);
            vector ngbDDtCorr = ((I-ngbN*ngbN)&ngbDDCorr);


            vector ownNonLinTracIncr = vector::zero;
            vector ngbNonLinTracIncr = vector::zero;

            if (nonLinear && (!enforceLinear))
            {
                tensor ownGradDDf =
                    ownSGradDD
                  + ownN*(interDD[faceI] - (ownDD + ownDDCorr))
                   /ownDeltaN;

                tensor ngbGradDDf =
                    ngbSGradDD
                  + ownN*((ngbDD + ngbDDCorr) - interDD[faceI])
                   /ngbDeltaN;

                symmTensor ownDEf =
                    symm(ownGradDDf) + 0.5*symm(ownGradDDf & ownGradDDf.T());
                symmTensor ngbDEf =
                    symm(ngbGradDDf) + 0.5*symm(ngbGradDDf & ngbGradDDf.T());

                symmTensor ownDSigmaf =
                    2*ownMu*ownDEf  + I*(ownLambda*tr(ownDEf));
                symmTensor ngbDSigmaf =
                    2*ngbMu*ngbDEf  + I*(ngbLambda*tr(ngbDEf));

                tensor ownNonLinDSigmaf =
                    ownMu*(ownGradDDf & ownGradDDf.T())
                  + 0.5*ownLambda*tr(ownGradDDf & ownGradDDf.T())*I
                  + (ownDSigmaf&ownGradDDf);
                tensor ngbNonLinDSigmaf =
                    ngbMu*(ngbGradDDf & ngbGradDDf.T())
                  + 0.5*ngbLambda*tr(ngbGradDDf & ngbGradDDf.T())*I
                  + (ngbDSigmaf&ngbGradDDf);

                ownNonLinTracIncr = (ownN&ownNonLinDSigmaf);
                ngbNonLinTracIncr = (ownN&ngbNonLinDSigmaf);
            }

            vector ownNonLinTracIncrN = ownN*(ownN&ownNonLinTracIncr);
            vector ngbNonLinTracIncrN = ngbN*(ngbN&ngbNonLinTracIncr);

            vector ownNonLinTracIncrT = ((I-ownN*ownN)&ownNonLinTracIncr);
            vector ngbNonLinTracIncrT = ((I-ngbN*ngbN)&ngbNonLinTracIncr);

            // Interface displacement

            vector curInterDDt =
                (
                    ownMu*(ownDDt+ownDDtCorr)*ngbDeltaN
                  + ngbMu*(ngbDDt+ngbDDtCorr)*ownDeltaN
                  + ownDeltaN*ngbDeltaN*(ngbMu*ngbSGradDDn - ownMu*ownSGradDDn)
                  + ownDeltaN*ngbDeltaN
                   *(ngbNonLinTracIncrT - ownNonLinTracIncrT)
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            vector curInterDDn =
                (
                    (2*ownMu + ownLambda)*(ownDDn+ownDDnCorr)*ngbDeltaN
                  + (2*ngbMu + ngbLambda)*(ngbDDn+ngbDDnCorr)*ownDeltaN
                  + ownDeltaN*ngbDeltaN
                   *(ngbLambda*ngbTrSGradDDt - ownLambda*ownTrSGradDDt)*ownN
                  + ownDeltaN*ngbDeltaN
                   *(ngbNonLinTracIncrN - ownNonLinTracIncrN)
                )
               /(
                   (2*ownMu + ownLambda)*ngbDeltaN
                 + (2*ngbMu + ngbLambda)*ownDeltaN
                );

            interDD[faceI] = curInterDDn + curInterDDt;


            // Implicit coupling

            scalar wRevLin =
                1.0
              - mesh().weights().boundaryField()[curPatch][curPatchFace];

            scalar ownK = (2*ownMu + ownLambda);
            scalar ngbK = (2*ngbMu + ngbLambda);

            scalar Kf = 1.0/(wRevLin/ownK + (1.0-wRevLin)/ngbK);
            scalar muf = 1.0/(wRevLin/ownMu + (1.0-wRevLin)/ngbMu);

            scalar DeltaNf =
                1.0
               /mesh().deltaCoeffs().boundaryField()[curPatch][curPatchFace];

            // Owner
            diag[curPatchCells[curPatchFace]] += Kf*magS/DeltaNf;

            boundaryCoeffs[curPatch][curPatchFace] +=
                Kf*magS*vector::one/DeltaNf;

            source[curPatchCells[curPatchFace]] +=
              - (Kf - muf)*(ngbDDt - ownDDt)*magS/DeltaNf;

            source[curPatchCells[curPatchFace]] +=
                (
                    ownK*ngbDeltaN*ngbLambda*ngbTrSGradDDt
                  + ngbK*ownDeltaN*ownLambda*ownTrSGradDDt
                )*ownN*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownMu*ngbMu*ngbDeltaN*ngbSGradDDn
                  + ownMu*ngbMu*ownDeltaN*ownSGradDDn
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            source[curPatchCells[curPatchFace]] +=
                (
                    ownK*ngbDeltaN*ngbNonLinTracIncrN
                  + ngbK*ownDeltaN*ownNonLinTracIncrN
                )*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbNonLinTracIncrT
                  + ngbMu*ownDeltaN*ownNonLinTracIncrT
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            source[curPatchCells[curPatchFace]] +=
                Kf*(ngbDDCorr - ownDDCorr)*magS/DeltaNf
              - (Kf - muf)*(ngbDDtCorr - ownDDtCorr)*magS/DeltaNf;


            // Interface traction increment

            vector curInterDTt =
                muf*((ngbDD+ngbDDCorr) - (ownDD+ownDDCorr))/DeltaNf
              + (
                    ownMu*ngbMu*ngbDeltaN*ngbSGradDDn
                  + ownMu*ngbMu*ownDeltaN*ownSGradDDn
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbNonLinTracIncrT
                  + ngbMu*ownDeltaN*ownNonLinTracIncrT
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            vector curInterDTn =
                Kf*((ngbDDn+ngbDDnCorr) - (ownDDn+ownDDnCorr))/DeltaNf
              + (
                    ownK*ngbDeltaN*ngbLambda*ngbTrSGradDDt
                  + ngbK*ownDeltaN*ownLambda*ownTrSGradDDt
                )*ownN
              + (
                    ownK*ngbDeltaN*ngbNonLinTracIncrN
                  + ngbK*ownDeltaN*ownNonLinTracIncrN
                )
               /(ownK*ngbDeltaN + ngbK*ownDeltaN);

            interDT[faceI] = curInterDTt + curInterDTn;
        }
    }
}


void ITLMaterialInterface::correct(surfaceVectorField& trac) const
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


void ITLMaterialInterface::updateDisplacementIncrement
(
    pointVectorField& pointDD
)
{
    if (debug)
    {
        Info<< "ITLMaterialInterface::updateDisplacementIncrement("
            << "pointVectorField&)"
            << "interpolating displacement incr field from cells to points"
            << endl;
    }

    materialInterface::updateDisplacement
    (
        DD_,
        displacementIncrement(),
        pointDD,
        subMeshDD(),
        subMeshPointDD()
    );
}


void ITLMaterialInterface::updateDisplacementIncrementGradient
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
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

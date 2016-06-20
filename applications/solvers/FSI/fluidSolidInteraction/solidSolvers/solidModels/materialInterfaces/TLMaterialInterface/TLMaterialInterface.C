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

#include "TLMaterialInterface.H"
#include "fvc.H"
#include "processorFvPatchFields.H"
#include "fvMatrices.H"
#include "skewCorrectionVectors.H"
#include "leastSquaresGrad.H"
#include "gaussGrad.H"
#include "faceSet.H"

#include "processorFvsPatchFields.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "solidSolver.H"
#include "fvcGradf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(TLMaterialInterface, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void TLMaterialInterface::makeDisplacement() const
{
    if (debug)
    {
        Info<< "void TLMaterialInterface::makeInterfaceDisplacement() const : "
            << "creating interface displacement field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (displacementPtr_)
    {
        FatalErrorIn("TLMaterialInterface::makeDisplacement() const")
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


void TLMaterialInterface::makeTraction() const
{
    if (debug)
    {
        Info<< "void TLMaterialInterface::makeTraction() const : "
            << "creating interface traction field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (tractionPtr_)
    {
        FatalErrorIn("TLMaterialInterface::makeTraction() const")
            << "interface traction field already exist"
            << abort(FatalError);
    }

    tractionPtr_ = new vectorField(faces().size(), vector::zero);
}


void TLMaterialInterface::makeSubMeshPointD() const
{
    if (debug)
    {
        Info<< "void TLMaterialInterface::makeSubMeshPointD() const : "
            << "creating point displacements fields"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshPointD_.empty())
    {
        FatalErrorIn("TLMaterialInterface::makeSubMeshPointD() const")
            << "Point displacement fields already exist"
            << abort(FatalError);
    }

    subMeshPointD_.setSize(subMeshes().size());

    forAll(subMeshPointD_, meshI)
    {
        subMeshPointD_.set
        (
            meshI,
            new pointVectorField
            (
                subMeshes()[meshI].interpolate(pointD_)
            )
        );
    }
}


void TLMaterialInterface::makeSubMeshD() const
{
    if (debug)
    {
        Info<< "void TLMaterialInterface::makeSubMeshD() const : "
            << "creating displacements fields"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshD_.empty())
    {
        FatalErrorIn("TLMaterialInterface::makeSubMeshD() const")
            << "Displacement fields already exist"
            << abort(FatalError);
    }

    subMeshD_.setSize(subMeshes().size());

    forAll(subMeshD_, meshI)
    {
        OStringStream SubsetName;
        SubsetName() << Pstream::myProcNo() << '_' << meshI << '_'
            << D_.name();

        subMeshD_.set
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
                dimensionedVector("0", D_.dimensions(), vector::zero)
            )
        );

        subMeshD_[meshI] = subMeshes()[meshI].interpolate(D_);

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
            vectorField& interfaceD =
                subMeshD_[meshI].boundaryField()[interfacePatchIndex];

            const labelList& fm = subMeshes()[meshI].faceMap();

            label interfacePatchStart =
                subMeshes()[meshI].subMesh().boundaryMesh()
                [
                    interfacePatchIndex
                ].start();

            forAll(interfaceD, faceI)
            {
                label curInterFace =
                    findIndex(faces(), fm[interfacePatchStart + faceI]);

                interfaceD[faceI] = displacement()[curInterFace];
            }
        }
    }
}


void TLMaterialInterface::clearOut()
{
    deleteDemandDrivenData(displacementPtr_);
    deleteDemandDrivenData(tractionPtr_);

//     subMeshes_.clear();
//     subMeshVolToPoint_.clear();

//     deleteDemandDrivenData(pointNumOfMaterialsPtr_);
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


TLMaterialInterface::TLMaterialInterface
(
    const volVectorField& D,
    const pointVectorField& pointD
)
:
    materialInterface(D.mesh()),
    D_(D),
    pointD_(pointD),
    displacementPtr_(NULL),
    tractionPtr_(NULL),
    subMeshD_(0),
    subMeshPointD_(0)
{}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

TLMaterialInterface::~TLMaterialInterface()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


vectorField& TLMaterialInterface::displacement()
{
    if (!displacementPtr_)
    {
        makeDisplacement();
    }

    return *displacementPtr_;
}

const vectorField& TLMaterialInterface::displacement() const
{
    if (!displacementPtr_)
    {
        makeDisplacement();
    }

    return *displacementPtr_;
}

vectorField& TLMaterialInterface::traction()
{
    if (!tractionPtr_)
    {
        makeTraction();
    }

    return *tractionPtr_;
}

const vectorField& TLMaterialInterface::traction() const
{
    if (!tractionPtr_)
    {
        makeTraction();
    }

    return *tractionPtr_;
}


const PtrList<volVectorField>& TLMaterialInterface::subMeshD() const
{
    if (subMeshD_.empty())
    {
        makeSubMeshD();
    }

    return subMeshD_;
}


PtrList<volVectorField>& TLMaterialInterface::subMeshD()
{
    if (subMeshD_.empty())
    {
        makeSubMeshD();
    }

    return subMeshD_;
}


const PtrList<pointVectorField>& TLMaterialInterface::subMeshPointD() const
{
    if (subMeshPointD_.empty())
    {
        makeSubMeshPointD();
    }

    return subMeshPointD_;
}


PtrList<pointVectorField>& TLMaterialInterface::subMeshPointD()
{
    if (subMeshPointD_.empty())
    {
        makeSubMeshPointD();
    }

    return subMeshPointD_;
}


void TLMaterialInterface::correct(fvVectorMatrix& DEqn)
{
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const vectorField& DI = D_.internalField();

    const volTensorField& gradD =
        mesh().lookupObject<volTensorField>("grad(" + D_.name() + ')');
    const tensorField& gradDI = gradD.internalField();

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

    vectorField& interD = displacement();
    vectorField& interT = traction();

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

    Switch thermalStress = stress.thermalStress();

    surfaceScalarField* DTfPtr = NULL;
    volScalarField* threeKPtr = NULL;
    volScalarField* alphaPtr = NULL;

    if (thermalStress)
    {
        DTfPtr =
            const_cast<surfaceScalarField*>
            (
                &mesh().lookupObject<surfaceScalarField>("DTf")
            );
        threeKPtr =
            const_cast<volScalarField*>
            (
                &mesh().lookupObject<volScalarField>("threeK")
            );
        alphaPtr =
            const_cast<volScalarField*>
            (
                &mesh().lookupObject<volScalarField>("alpha")
            );
    }
    const surfaceScalarField& DTf = *DTfPtr;
    const volScalarField& threeK  = *threeKPtr;
    const volScalarField& alpha  = *alphaPtr;

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

//             tensor ownSGradD = ((I-ownN*ownN)&gradDI[curOwner]);
//             tensor ngbSGradD = ((I-ngbN*ngbN)&gradDI[curNeighbour]);
            tensor ownSGradD = ((I-ownN*ownN)&gradDfI[curFace]);
            tensor ngbSGradD = ((I-ngbN*ngbN)&gradDfI[curFace]);


            scalar ownTrSGradDt = tr(ownSGradD&(I-ownN*ownN));
            scalar ngbTrSGradDt = tr(ngbSGradD&(I-ownN*ownN));

            vector ownSGradDn = (ownSGradD&ownN);
            vector ngbSGradDn = (ngbSGradD&ownN);

            scalar ngbDeltaN = w[curFace]*(1.0/deltaCoeffs[curFace]);
            scalar ownDeltaN = (1.0/deltaCoeffs[curFace]) - ngbDeltaN;

            scalar ownMu  = mu.internalField()[curOwner];
            scalar ngbMu  = mu.internalField()[curNeighbour];

            scalar ownLambda  = lambda.internalField()[curOwner];
            scalar ngbLambda  = lambda.internalField()[curNeighbour];

            vector ownD = DI[curOwner];
            vector ngbD = DI[curNeighbour];

            vector ownDt = ((I-ownN*ownN)&DI[curOwner]);
            vector ngbDt = ((I-ngbN*ngbN)&DI[curNeighbour]);

            vector ownDn = ownN*(ownN&DI[curOwner]);
            vector ngbDn = ngbN*(ngbN&DI[curNeighbour]);

            vector ownCorrVec = Cf[curFace] - C[curOwner];
            ownCorrVec -= ownN*(ownN&ownCorrVec);

            vector ngbCorrVec = Cf[curFace] - C[curNeighbour];
            ngbCorrVec -= ngbN*(ngbN&ngbCorrVec);

            vector ownDCorr = (ownCorrVec&gradDI[curOwner]);
            vector ngbDCorr = (ngbCorrVec&gradDI[curNeighbour]);

            vector ownDnCorr = ownN*(ownN&ownDCorr);
            vector ngbDnCorr = ngbN*(ngbN&ngbDCorr);

            vector ownDtCorr = ((I-ownN*ownN)&ownDCorr);
            vector ngbDtCorr = ((I-ngbN*ngbN)&ngbDCorr);


            vector ownThermalTraction = vector::zero;
            vector ngbThermalTraction = vector::zero;

            if (thermalStress)
            {
                scalar ownThreeK = threeK.internalField()[curOwner];
                scalar ngbThreeK = threeK.internalField()[curNeighbour];
                scalar ownAlpha = alpha.internalField()[curOwner];
                scalar ngbAlpha = alpha.internalField()[curNeighbour];
                scalar interfaceDT = DTf.internalField()[curFace];
                ownThermalTraction = -ownThreeK*ownAlpha*interfaceDT*ownN;
                ngbThermalTraction = -ngbThreeK*ngbAlpha*interfaceDT*ownN;
            }


            vector ownNonLinTrac = vector::zero;
            vector ngbNonLinTrac = vector::zero;

            if (nonLinear && (!enforceLinear))
            {
                tensor ownGradDf =
                    ownSGradD
                  + ownN*(interD[faceI] - (DI[curOwner] + ownDCorr))
                   /ownDeltaN;

                tensor ngbGradDf =
                    ngbSGradD
                  + ownN*((DI[curNeighbour] + ngbDCorr) - interD[faceI])
                   /ngbDeltaN;

                symmTensor ownEf =
                    symm(ownGradDf) + 0.5*symm(ownGradDf & ownGradDf.T());
                symmTensor ngbEf =
                    symm(ngbGradDf) + 0.5*symm(ngbGradDf & ngbGradDf.T());

                symmTensor ownSigmaf =
                    2*ownMu*ownEf  + I*(ownLambda*tr(ownEf));
                symmTensor ngbSigmaf =
                    2*ngbMu*ngbEf  + I*(ngbLambda*tr(ngbEf));

                tensor ownNonLinSigmaf =
                    ownMu*(ownGradDf & ownGradDf.T())
                  + 0.5*ownLambda*tr(ownGradDf & ownGradDf.T())*I
                  + (ownSigmaf&ownGradDf);
                tensor ngbNonLinSigmaf =
                    ngbMu*(ngbGradDf & ngbGradDf.T())
                  + 0.5*ngbLambda*tr(ngbGradDf & ngbGradDf.T())*I
                  + (ngbSigmaf&ngbGradDf);

                ownNonLinTrac = (ownN&ownNonLinSigmaf);
                ngbNonLinTrac = (ownN&ngbNonLinSigmaf);

                if (thermalStress)
                {
                    tensor ownF = I+ownGradDf;
                    tensor ownInvF = inv(ownF);
                    scalar ownJ = det(ownF);

                    ownThermalTraction = ownJ*(ownInvF&ownThermalTraction);

                    tensor ngbF = I+ngbGradDf;
                    tensor ngbInvF = inv(ngbF);
                    scalar ngbJ = det(ngbF);

                    ngbThermalTraction = ngbJ*(ngbInvF&ngbThermalTraction);
                }
            }

            vector ownNonLinTracN = ownN*(ownN&ownNonLinTrac);
            vector ngbNonLinTracN = ngbN*(ngbN&ngbNonLinTrac);

            vector ownNonLinTracT = ((I-ownN*ownN)&ownNonLinTrac);
            vector ngbNonLinTracT = ((I-ngbN*ngbN)&ngbNonLinTrac);

            vector ownThermalTractionN = ownN*(ownN&ownThermalTraction);
            vector ngbThermalTractionN = ngbN*(ngbN&ngbThermalTraction);

            vector ownThermalTractionT = ((I-ownN*ownN)&ownThermalTraction);
            vector ngbThermalTractionT = ((I-ngbN*ngbN)&ngbThermalTraction);

            // Interface displacement

            vector curInterDt =
                (
                    ownMu*(ownDt+ownDtCorr)*ngbDeltaN
                  + ngbMu*(ngbDt+ngbDtCorr)*ownDeltaN
                  + ownDeltaN*ngbDeltaN*(ngbMu*ngbSGradDn - ownMu*ownSGradDn)
                  + ownDeltaN*ngbDeltaN*(ngbNonLinTracT - ownNonLinTracT)
                  + ownDeltaN*ngbDeltaN
                   *(ngbThermalTractionT - ownThermalTractionT)
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            vector curInterDn =
                (
                    (2*ownMu + ownLambda)*(ownDn+ownDnCorr)*ngbDeltaN
                  + (2*ngbMu + ngbLambda)*(ngbDn+ngbDnCorr)*ownDeltaN
                  + ownDeltaN*ngbDeltaN
                   *(ngbLambda*ngbTrSGradDt - ownLambda*ownTrSGradDt)*ownN
                  + ownDeltaN*ngbDeltaN
                   *(ngbNonLinTracN - ownNonLinTracN)
                  + ownDeltaN*ngbDeltaN
                   *(ngbThermalTractionN - ownThermalTractionN)
                )
               /(
                   (2*ownMu + ownLambda)*ngbDeltaN
                 + (2*ngbMu + ngbLambda)*ownDeltaN
                );

            interD[faceI] = curInterDn + curInterDt;


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
               *((DI[curNeighbour] - DI[curOwner])&(I-ownN*ownN))*magS/DeltaNf;

            source[curOwner] +=
                (
                    ownK*ngbDeltaN*ngbLambda*ngbTrSGradDt
                  + ngbK*ownDeltaN*ownLambda*ownTrSGradDt
                )*ownN*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownMu*ngbMu*ngbDeltaN*ngbSGradDn
                  + ownMu*ngbMu*ownDeltaN*ownSGradDn
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            source[curOwner] +=
                (
                    ownK*ngbDeltaN*ngbNonLinTracN
                  + ngbK*ownDeltaN*ownNonLinTracN
                )*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownK*ngbDeltaN*ngbThermalTractionN
                  + ngbK*ownDeltaN*ownThermalTractionN
                )*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbNonLinTracT
                  + ngbMu*ownDeltaN*ownNonLinTracT
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbThermalTractionT
                  + ngbMu*ownDeltaN*ownThermalTractionT
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            source[curOwner] +=
                Kf*(ngbDCorr - ownDCorr)*magS/DeltaNf
              - (Kf - muf)*(ngbDtCorr - ownDtCorr)*magS/DeltaNf;

            // Neighbour
            diag[curNeighbour] += Kf*magS/DeltaNf;

            source[curNeighbour] -=
              - (Kf - muf)
               *((DI[curNeighbour] - DI[curOwner])&(I-ownN*ownN))*magS/DeltaNf;

            source[curNeighbour] -=
                (
                    ownK*ngbDeltaN*ngbLambda*ngbTrSGradDt
                  + ngbK*ownDeltaN*ownLambda*ownTrSGradDt
                )*ownN*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownMu*ngbMu*ngbDeltaN*ngbSGradDn
                  + ownMu*ngbMu*ownDeltaN*ownSGradDn
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            source[curNeighbour] -=
                (
                    ownK*ngbDeltaN*ngbNonLinTracN
                  + ngbK*ownDeltaN*ownNonLinTracN
                )*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownK*ngbDeltaN*ngbThermalTractionN
                  + ngbK*ownDeltaN*ownThermalTractionN
                )*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbNonLinTracT
                  + ngbMu*ownDeltaN*ownNonLinTracT
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbThermalTractionT
                  + ngbMu*ownDeltaN*ownThermalTractionT
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            source[curNeighbour] -=
                Kf*(ngbDCorr - ownDCorr)*magS/DeltaNf
              - (Kf - muf)*(ngbDtCorr - ownDtCorr)*magS/DeltaNf;


            // Interface traction

            vector curInterTt =
                muf*((ngbD+ngbDCorr) - (ownD+ownDCorr))/DeltaNf
              + (
                    ownMu*ngbMu*ngbDeltaN*ngbSGradDn
                  + ownMu*ngbMu*ownDeltaN*ownSGradDn
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbNonLinTracT
                  + ngbMu*ownDeltaN*ownNonLinTracT
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbThermalTractionT
                  + ngbMu*ownDeltaN*ownThermalTractionT
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            vector curInterTn =
                Kf*((ngbDn+ngbDnCorr) - (ownDn+ownDnCorr))/DeltaNf
              + (
                    ownK*ngbDeltaN*ngbLambda*ngbTrSGradDt
                  + ngbK*ownDeltaN*ownLambda*ownTrSGradDt
                )*ownN
              + (
                    ownK*ngbDeltaN*ngbNonLinTracN
                  + ngbK*ownDeltaN*ownNonLinTracN
                )
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownK*ngbDeltaN*ngbThermalTractionN
                  + ngbK*ownDeltaN*ownThermalTractionN
                )
               /(ownK*ngbDeltaN + ngbK*ownDeltaN);

            interT[faceI] = curInterTt + curInterTn;
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

            vector ownD = D_.internalField()[curPatchCells[curPatchFace]];
            vector ngbD = D_.boundaryField()[curPatch][curPatchFace];

            vector ownDt = ((I-ownN*ownN)&ownD);
            vector ngbDt = ((I-ngbN*ngbN)&ngbD);

            vector ownDn = ownN*(ownN&ownD);
            vector ngbDn = ngbN*(ngbN&ngbD);

            tensor ownSGradD =
                ((I-ownN*ownN)&gradDf.boundaryField()[curPatch][curPatchFace]);
            tensor ngbSGradD = ownSGradD;

            scalar ownTrSGradDt = tr(ownSGradD&(I-ownN*ownN));
            scalar ngbTrSGradDt = tr(ngbSGradD&(I-ngbN*ngbN));

            vector ownSGradDn = (ownSGradD&ownN);
            vector ngbSGradDn = (ngbSGradD&ownN);

            vector ownCorrVec =
                mesh().Cf().boundaryField()[curPatch][curPatchFace]
              - mesh().C().internalField()[curPatchCells[curPatchFace]];
            ownCorrVec -= ownN*(ownN&ownCorrVec);

            vector ngbCorrVec =
                mesh().Cf().boundaryField()[curPatch][curPatchFace]
              - mesh().C().boundaryField()[curPatch][curPatchFace];
            ngbCorrVec -= ngbN*(ngbN&ngbCorrVec);

            vector ownDCorr =
                (
                    ownCorrVec
                  & gradD.internalField()[curPatchCells[curPatchFace]]
                );
            vector ngbDCorr =
                (ngbCorrVec&gradD.boundaryField()[curPatch][curPatchFace]);

            vector ownDnCorr = ownN*(ownN&ownDCorr);
            vector ngbDnCorr = ownN*(ownN&ngbDCorr);

            vector ownDtCorr = ((I-ownN*ownN)&ownDCorr);
            vector ngbDtCorr = ((I-ngbN*ngbN)&ngbDCorr);


            vector ownThermalTraction = vector::zero;
            vector ngbThermalTraction = vector::zero;

            if (thermalStress)
            {
                scalar ownThreeK =
                    threeK.internalField()[curPatchCells[curPatchFace]];
                scalar ngbThreeK =
                    threeK.boundaryField()[curPatch][curPatchFace];
                scalar ownAlpha =
                    alpha.internalField()[curPatchCells[curPatchFace]];
                scalar ngbAlpha =
                    alpha.boundaryField()[curPatch][curPatchFace];
                scalar interfaceDT =
                    DTf.boundaryField()[curPatch][curPatchFace];
                ownThermalTraction = -ownThreeK*ownAlpha*interfaceDT*ownN;
                ngbThermalTraction = -ngbThreeK*ngbAlpha*interfaceDT*ownN;
            }


            vector ownNonLinTrac = vector::zero;
            vector ngbNonLinTrac = vector::zero;

            if (nonLinear && (!enforceLinear))
            {
                tensor ownGradDf =
                    ownSGradD
                  + ownN*(interD[faceI] - (ownD + ownDCorr))
                   /ownDeltaN;

                tensor ngbGradDf =
                    ngbSGradD
                  + ownN*((ngbD + ngbDCorr) - interD[faceI])
                   /ngbDeltaN;

                symmTensor ownEf =
                    symm(ownGradDf) + 0.5*symm(ownGradDf & ownGradDf.T());
                symmTensor ngbEf =
                    symm(ngbGradDf) + 0.5*symm(ngbGradDf & ngbGradDf.T());

                symmTensor ownSigmaf =
                    2*ownMu*ownEf  + I*(ownLambda*tr(ownEf));
                symmTensor ngbSigmaf =
                    2*ngbMu*ngbEf  + I*(ngbLambda*tr(ngbEf));

                tensor ownNonLinSigmaf =
                    ownMu*(ownGradDf & ownGradDf.T())
                  + 0.5*ownLambda*tr(ownGradDf & ownGradDf.T())*I
                  + (ownSigmaf&ownGradDf);
                tensor ngbNonLinSigmaf =
                    ngbMu*(ngbGradDf & ngbGradDf.T())
                  + 0.5*ngbLambda*tr(ngbGradDf & ngbGradDf.T())*I
                  + (ngbSigmaf&ngbGradDf);

                ownNonLinTrac = (ownN&ownNonLinSigmaf);
                ngbNonLinTrac = (ownN&ngbNonLinSigmaf);

                if (thermalStress)
                {
                    tensor ownF = I+ownGradDf;
                    tensor ownInvF = inv(ownF);
                    scalar ownJ = det(ownF);

                    ownThermalTraction = ownJ*(ownInvF&ownThermalTraction);

                    tensor ngbF = I+ngbGradDf;
                    tensor ngbInvF = inv(ngbF);
                    scalar ngbJ = det(ngbF);

                    ngbThermalTraction = ngbJ*(ngbInvF&ngbThermalTraction);
                }
            }

            vector ownNonLinTracN = ownN*(ownN&ownNonLinTrac);
            vector ngbNonLinTracN = ngbN*(ngbN&ngbNonLinTrac);

            vector ownNonLinTracT = ((I-ownN*ownN)&ownNonLinTrac);
            vector ngbNonLinTracT = ((I-ngbN*ngbN)&ngbNonLinTrac);

            vector ownThermalTractionN = ownN*(ownN&ownThermalTraction);
            vector ngbThermalTractionN = ngbN*(ngbN&ngbThermalTraction);

            vector ownThermalTractionT = ((I-ownN*ownN)&ownThermalTraction);
            vector ngbThermalTractionT = ((I-ngbN*ngbN)&ngbThermalTraction);


            // Interface displacement

            vector curInterDt =
                (
                    ownMu*(ownDt+ownDtCorr)*ngbDeltaN
                  + ngbMu*(ngbDt+ngbDtCorr)*ownDeltaN
                  + ownDeltaN*ngbDeltaN*(ngbMu*ngbSGradDn - ownMu*ownSGradDn)
                  + ownDeltaN*ngbDeltaN*(ngbNonLinTracT - ownNonLinTracT)
                  + ownDeltaN*ngbDeltaN
                   *(ngbThermalTractionT - ownThermalTractionT)
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            vector curInterDn =
                (
                    (2*ownMu + ownLambda)*(ownDn+ownDnCorr)*ngbDeltaN
                  + (2*ngbMu + ngbLambda)*(ngbDn+ngbDnCorr)*ownDeltaN
                  + ownDeltaN*ngbDeltaN
                   *(ngbLambda*ngbTrSGradDt - ownLambda*ownTrSGradDt)*ownN
                  + ownDeltaN*ngbDeltaN
                   *(ngbNonLinTracN - ownNonLinTracN)
                  + ownDeltaN*ngbDeltaN
                   *(ngbThermalTractionN - ownThermalTractionN)
                )
               /(
                   (2*ownMu + ownLambda)*ngbDeltaN
                 + (2*ngbMu + ngbLambda)*ownDeltaN
                );

            interD[faceI] = curInterDn + curInterDt;


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
              - (Kf - muf)*(ngbDt - ownDt)*magS/DeltaNf;

            source[curPatchCells[curPatchFace]] +=
                (
                    ownK*ngbDeltaN*ngbLambda*ngbTrSGradDt
                  + ngbK*ownDeltaN*ownLambda*ownTrSGradDt
                )*ownN*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownMu*ngbMu*ngbDeltaN*ngbSGradDn
                  + ownMu*ngbMu*ownDeltaN*ownSGradDn
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            source[curPatchCells[curPatchFace]] +=
                (
                    ownK*ngbDeltaN*ngbNonLinTracN
                  + ngbK*ownDeltaN*ownNonLinTracN
                )*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownK*ngbDeltaN*ngbThermalTractionN
                  + ngbK*ownDeltaN*ownThermalTractionN
                )*magS
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbNonLinTracT
                  + ngbMu*ownDeltaN*ownNonLinTracT
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbThermalTractionT
                  + ngbMu*ownDeltaN*ownThermalTractionT
                )*magS
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            source[curPatchCells[curPatchFace]] +=
                Kf*(ngbDCorr - ownDCorr)*magS/DeltaNf
              - (Kf - muf)*(ngbDtCorr - ownDtCorr)*magS/DeltaNf;


            // Interface traction

            vector curInterTt =
                muf*((ngbD+ngbDCorr) - (ownD+ownDCorr))/DeltaNf
              + (
                    ownMu*ngbMu*ngbDeltaN*ngbSGradDn
                  + ownMu*ngbMu*ownDeltaN*ownSGradDn
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbNonLinTracT
                  + ngbMu*ownDeltaN*ownNonLinTracT
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN)
              + (
                    ownMu*ngbDeltaN*ngbThermalTractionT
                  + ngbMu*ownDeltaN*ownThermalTractionT
                )
               /(ownMu*ngbDeltaN + ngbMu*ownDeltaN);

            vector curInterTn =
                Kf*((ngbDn+ngbDnCorr) - (ownDn+ownDnCorr))/DeltaNf
              + (
                    ownK*ngbDeltaN*ngbLambda*ngbTrSGradDt
                  + ngbK*ownDeltaN*ownLambda*ownTrSGradDt
                )*ownN
              + (
                    ownK*ngbDeltaN*ngbNonLinTracN
                  + ngbK*ownDeltaN*ownNonLinTracN
                )
               /(ownK*ngbDeltaN + ngbK*ownDeltaN)
              + (
                    ownK*ngbDeltaN*ngbThermalTractionN
                  + ngbK*ownDeltaN*ownThermalTractionN
                )
               /(ownK*ngbDeltaN + ngbK*ownDeltaN);

            interT[faceI] = curInterTt + curInterTn;
        }
    }
}


void TLMaterialInterface::correct(surfaceVectorField& trac) const
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


void TLMaterialInterface::updateDisplacement(pointVectorField& pointD)
{
    if (debug)
    {
        Info<< "TLMaterialInterface::updateDisplacement("
            << "pointVectorField&)"
            << "interpolating fields from cells to points"
            << endl;
    }

    subMeshPointD_.clear();

//     subMeshD();

//     subMeshPointD();

    materialInterface::updateDisplacement
    (
        D_,
        displacement(),
        pointD,
        subMeshD(),
        subMeshPointD()
    );
}


void TLMaterialInterface::updateDisplacementGradient
(
    volTensorField& gradD,
    surfaceTensorField& gradDf
)
{
    materialInterface::updateDisplacementGradient
    (
        D_,
        subMeshD(),
        subMeshPointD(),
        gradD,
        gradDf
    );
}


bool TLMaterialInterface::update()
{
    clearOut();

    subMeshD_.clear();
    subMeshPointD_.clear();

    materialInterface::update();

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

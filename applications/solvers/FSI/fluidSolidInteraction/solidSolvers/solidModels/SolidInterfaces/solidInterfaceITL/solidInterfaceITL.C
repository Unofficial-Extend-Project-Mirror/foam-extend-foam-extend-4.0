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

#include "solidInterfaceITL.H"
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

defineTypeNameAndDebug(solidInterfaceITL, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void solidInterfaceITL::makeFaces() const
{
    if (debug)
    {
        Info<< "void solidInterfaceITL::makeFaces() const : "
            << "creating list of interface faces"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (facesPtr_)
    {
        FatalErrorIn("solidInterfaceITL::makeFaces() const")
            << "list of interface faces already exists"
            << abort(FatalError);
    }

    const fvMesh& mesh_ = D_.mesh();

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
                        label globalFaceID =
                            mesh_.boundaryMesh()[patchI].start() + faceI;

                        interFacesSet.insert(globalFaceID);
                    }
                }
            }
        }

        facesPtr_ = new labelList(interFacesSet.toc());
    }
    else
    {
        facesPtr_ = new labelList(0);
    }
}

void solidInterfaceITL::makeDisplacementIncrement() const
{
    if (debug)
    {
        Info<< "void solidInterfaceITL::"
            << "makeInterfaceDisplacementIncrement() const : "
            << "creating interface displacement field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (displacementIncrementPtr_)
    {
        FatalErrorIn("solidInterfaceITL::makeDisplacementIncrement() const")
            << "interface displacement increment field already exist"
            << abort(FatalError);
    }

    displacementIncrementPtr_ =
        new vectorField(faces().size(), vector::zero);

    // Initialize displacement increment
    surfaceVectorField DDf = fvc::interpolate(DD_);

    const fvMesh& mesh_ = DD_.mesh();

    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh_.nInternalFaces())
        {
            (*displacementIncrementPtr_)[faceI] = DDf.internalField()[curFace];
        }
        else
        {
            label curPatch = mesh_.boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh_.boundaryMesh()[curPatch].start();

            (*displacementIncrementPtr_)[faceI] =
                DDf.boundaryField()[curPatch][curPatchFace];
        }
    }
}

void solidInterfaceITL::makeDisplacement() const
{
    if (debug)
    {
        Info<< "void solidInterfaceITL::makeInterfaceDisplacement() const : "
            << "creating interface displacement field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (displacementPtr_)
    {
        FatalErrorIn("solidInterfaceITL::makeDisplacement() const")
            << "interface displacement field already exist"
            << abort(FatalError);
    }

    displacementPtr_ = new vectorField(faces().size(), vector::zero);

    // Initialize displacement
    surfaceVectorField Df = fvc::interpolate(D_);

    const fvMesh& mesh_ = D_.mesh();

    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh_.nInternalFaces())
        {
            (*displacementPtr_)[faceI] = Df.internalField()[curFace];
        }
        else
        {
            label curPatch = mesh_.boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh_.boundaryMesh()[curPatch].start();

            (*displacementPtr_)[faceI] =
                Df.boundaryField()[curPatch][curPatchFace];
        }
    }
}

void solidInterfaceITL::makeTractionIncrement() const
{
    if (debug)
    {
        Info<< "void solidInterfaceITL::makeTractionIncrement() const : "
            << "creating interface traction increment field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (tractionIncrementPtr_)
    {
        FatalErrorIn("solidInterfaceITL::makeTractionIncrement() const")
            << "interface traction increment field already exist"
            << abort(FatalError);
    }

    tractionIncrementPtr_ =
        new vectorField(faces().size(), vector::zero);
}

void solidInterfaceITL::makeTraction() const
{
    if (debug)
    {
        Info<< "void solidInterfaceITL::makeTraction() const : "
            << "creating interface traction field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (tractionPtr_)
    {
        FatalErrorIn("solidInterfaceITL::makeTraction() const")
            << "interface traction field already exist"
            << abort(FatalError);
    }

    tractionPtr_ = new vectorField(faces().size(), vector::zero);
}

void solidInterfaceITL::makeDSigmaf() const
{
    if (debug)
    {
        Info<< "void solidInterfaceITL::makeDSigmaf() const : "
            << "creating interface second Piola-Kirchhoff stress incr field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ngbDSigmafPtr_ || ownDSigmafPtr_)
    {
        FatalErrorIn("solidInterfaceITL::makeDSigmaf() const")
            << "interface second Piola-Kirchhoff stress incr field "
                << "already exists"
                << abort(FatalError);
    }

    ownDSigmafPtr_ = new symmTensorField(faces().size(), symmTensor::zero);
    ngbDSigmafPtr_ = new symmTensorField(faces().size(), symmTensor::zero);
}

void solidInterfaceITL::makeSigmaf() const
{
    if (debug)
    {
        Info<< "void solidInterfaceITL::makeSigmaf() const : "
            << "creating interface old second Piola-Kirchhoff stress field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ngbSigmafPtr_ || ownSigmafPtr_)
    {
        FatalErrorIn("solidInterfaceITL::makeSigmaf() const")
            << "interface old second Piola-Kirchhoff stress field "
                << "already exists"
                << abort(FatalError);
    }

    ownSigmafPtr_ = new symmTensorField(faces().size(), symmTensor::zero);
    ngbSigmafPtr_ = new symmTensorField(faces().size(), symmTensor::zero);
}

void solidInterfaceITL::makeSubMeshes() const
{
    if (debug)
    {
        Info<< "void solidInterfaceITL::makeSubMeshes() const : "
            << "creating material sub-meshes"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshes_.empty())
    {
        FatalErrorIn("solidInterfaceITL::makeSubMeshes() const")
            << "material sub-meshes already exist"
            << abort(FatalError);
    }

    const fvMesh& mesh_ = D_.mesh();

    const volScalarField& materials =
        mesh_.lookupObject<volScalarField>("materials");
    const scalarField& materialsI = materials.internalField();

    labelList region(materialsI.size(), -1);
    forAll(region, cellI)
    {
        region[cellI] = materialsI[cellI];
    }

    label nMat = gMax(materialsI) + 1;

    subMeshes_.setSize(nMat);

    forAll(subMeshes_, matI)
    {
        OStringStream SubsetName;
        SubsetName() << Pstream::myProcNo() << '_' << matI << '_';

        subMeshes_.set
        (
            matI,
            new fvMeshSubset
            (
                IOobject
                (
                    word(SubsetName.str()),
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        subMeshes_[matI].setLargeCellSubset(region, matI);
    }
}


void solidInterfaceITL::makeSubMeshVolToPoint() const
{
    if (debug)
    {
        Info<< "void solidInterfaceITL::makeVolToPointInterpolators() const : "
            << "creating cell-to-point interpolators"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshVolToPoint_.empty())
    {
        FatalErrorIn("solidInterfaceITL::makeVolToPointInterpolators() const")
            << "Cell-to-point intrpolators already exist"
            << abort(FatalError);
    }

    subMeshVolToPoint_.setSize(subMeshes().size());

    forAll(subMeshVolToPoint_, meshI)
    {
        subMeshVolToPoint_.set
        (
            meshI,
            new leastSquaresVolPointInterpolation
            (
                subMeshes()[meshI].subMesh()
            )
        );
    }
}

void solidInterfaceITL::makePointNumOfMaterials() const
{
    if (debug)
    {
        Info<< "void solidInterfaceITL::makePointNoMaterials() const : "
            << "creating number of materials for each point"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (pointNumOfMaterialsPtr_)
    {
        FatalErrorIn("solidInterfaceITL::makePointNumOfMaterials() const")
            << "Point numbero of materials already exist"
            << abort(FatalError);
    }

    const fvMesh& mesh_ = D_.mesh();

    pointNumOfMaterialsPtr_ = new labelList(mesh_.nPoints(), 0);
    labelList& pointNumOfMaterials = *pointNumOfMaterialsPtr_;

    const labelListList& pointCells = mesh_.pointCells();

    const volScalarField& materials =
        mesh_.lookupObject<volScalarField>("materials");
    const scalarField& materialsI = materials.internalField();

    forAll(pointNumOfMaterials, pointI)
    {
        const labelList& curCells = pointCells[pointI];

        labelHashSet matSet;

        forAll(curCells, cellI)
        {
            if (!matSet.found(materialsI[curCells[cellI]]))
            {
                matSet.insert(materialsI[curCells[cellI]]);
            }
        }

        pointNumOfMaterials[pointI] = matSet.toc().size();
    }
}


void solidInterfaceITL::clearOut()
{
    deleteDemandDrivenData(facesPtr_);
    deleteDemandDrivenData(displacementIncrementPtr_);
    deleteDemandDrivenData(displacementPtr_);
    deleteDemandDrivenData(tractionIncrementPtr_);
    deleteDemandDrivenData(tractionPtr_);

    deleteDemandDrivenData(ownDSigmafPtr_);
    deleteDemandDrivenData(ngbDSigmafPtr_);

    deleteDemandDrivenData(ownSigmafPtr_);
    deleteDemandDrivenData(ngbSigmafPtr_);

//     subMeshes_.clear();
//     subMeshVolToPoint_.clear();

    deleteDemandDrivenData(pointNumOfMaterialsPtr_);
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


solidInterfaceITL::solidInterfaceITL
(
    const volVectorField& DD,
    const volVectorField& D
)
:
    regIOobject
    (
        IOobject
        (
            "solidInterfaceITL",
            D.mesh().time().constant(),
            D.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    DD_(DD),
    D_(D),
    curTimeIndex_(-1),
    facesPtr_(NULL),
    displacementIncrementPtr_(NULL),
    displacementPtr_(NULL),
    tractionIncrementPtr_(NULL),
    tractionPtr_(NULL),
    ownDSigmafPtr_(NULL),
    ngbDSigmafPtr_(NULL),
    ownSigmafPtr_(NULL),
    ngbSigmafPtr_(NULL),
    subMeshes_(0),
    subMeshVolToPoint_(0),
    pointNumOfMaterialsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

solidInterfaceITL::~solidInterfaceITL()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


const labelList& solidInterfaceITL::faces() const
{
    if (!facesPtr_)
    {
        makeFaces();
    }

    return *facesPtr_;
}

vectorField& solidInterfaceITL::displacementIncrement()
{
    if (!displacementIncrementPtr_)
    {
        makeDisplacementIncrement();
    }

    return *displacementIncrementPtr_;
}

const vectorField& solidInterfaceITL::displacementIncrement() const
{
    if (!displacementIncrementPtr_)
    {
        makeDisplacementIncrement();
    }

    return *displacementIncrementPtr_;
}

vectorField& solidInterfaceITL::displacement()
{
    if (!displacementPtr_)
    {
        makeDisplacement();
    }

    return *displacementPtr_;
}

const vectorField& solidInterfaceITL::displacement() const
{
    if (!displacementPtr_)
    {
        makeDisplacement();
    }

    return *displacementPtr_;
}

vectorField& solidInterfaceITL::tractionIncrement()
{
    if (!tractionIncrementPtr_)
    {
        makeTractionIncrement();
    }

    return *tractionIncrementPtr_;
}

const vectorField& solidInterfaceITL::tractionIncrement() const
{
    if (!tractionIncrementPtr_)
    {
        makeTractionIncrement();
    }

    return *tractionIncrementPtr_;
}

vectorField& solidInterfaceITL::traction()
{
    if (!tractionPtr_)
    {
        makeTraction();
    }

    return *tractionPtr_;
}

const vectorField& solidInterfaceITL::traction() const
{
    if (!tractionPtr_)
    {
        makeTraction();
    }

    return *tractionPtr_;
}

symmTensorField& solidInterfaceITL::ownSigmaf()
{
    if (!ownSigmafPtr_)
    {
        makeSigmaf();
    }

    return *ownSigmafPtr_;
}

const symmTensorField& solidInterfaceITL::ownSigmaf() const
{
    if (!ownSigmafPtr_)
    {
        makeSigmaf();
    }

    return *ownSigmafPtr_;
}

symmTensorField& solidInterfaceITL::ngbSigmaf()
{
    if (!ngbSigmafPtr_)
    {
        makeSigmaf();
    }

    return *ngbSigmafPtr_;
}

const symmTensorField& solidInterfaceITL::ngbSigmaf() const
{
    if (!ngbSigmafPtr_)
    {
        makeSigmaf();
    }

    return *ngbSigmafPtr_;
}

symmTensorField& solidInterfaceITL::ownDSigmaf()
{
    if (!ownDSigmafPtr_)
    {
        makeDSigmaf();
    }

    return *ownDSigmafPtr_;
}

const symmTensorField& solidInterfaceITL::ownDSigmaf() const
{
    if (!ownDSigmafPtr_)
    {
        makeDSigmaf();
    }

    return *ownDSigmafPtr_;
}

symmTensorField& solidInterfaceITL::ngbDSigmaf()
{
    if (!ngbDSigmafPtr_)
    {
        makeDSigmaf();
    }

    return *ngbDSigmafPtr_;
}

const symmTensorField& solidInterfaceITL::ngbDSigmaf() const
{
    if (!ngbDSigmafPtr_)
    {
        makeDSigmaf();
    }

    return *ngbDSigmafPtr_;
}

const PtrList<fvMeshSubset>& solidInterfaceITL::subMeshes() const
{
    if (subMeshes_.empty())
    {
        makeSubMeshes();
    }

    return subMeshes_;
}

const PtrList<leastSquaresVolPointInterpolation>& solidInterfaceITL
::subMeshVolToPoint() const
{
    if (subMeshVolToPoint_.empty())
    {
        makeSubMeshVolToPoint();
    }

    return subMeshVolToPoint_;
}


const labelList& solidInterfaceITL::pointNumOfMaterials() const
{
    if (!pointNumOfMaterialsPtr_)
    {
        makePointNumOfMaterials();
    }

    return *pointNumOfMaterialsPtr_;
}


void solidInterfaceITL::correct(fvVectorMatrix& DEqn)
{
    const fvMesh& mesh_ = D_.mesh();

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const vectorField& DDI = DD_.internalField();

    const vectorField& DI = D_.internalField();

    const volTensorField& gradDD =
        mesh_.lookupObject<volTensorField>("grad(" + DD_.name() + ')');
    const tensorField& gradDDI = gradDD.internalField();

    const volTensorField& gradD =
        mesh_.lookupObject<volTensorField>("grad(" + D_.name() + ')');
    const tensorField& gradDI = gradD.internalField();

    const surfaceTensorField& gradDDf =
        mesh_.lookupObject<surfaceTensorField>("grad" + DD_.name() + 'f');
    const tensorField& gradDDfI = gradDDf.internalField();

    const surfaceTensorField& gradDf =
        mesh_.lookupObject<surfaceTensorField>("grad" + D_.name() + 'f');
    const tensorField& gradDfI = gradDf.internalField();

    const volScalarField& mu =
        mesh_.lookupObject<volScalarField>("mu");
    const volScalarField& lambda =
        mesh_.lookupObject<volScalarField>("lambda");

    const vectorField& SI  = mesh_.Sf().internalField();
    const scalarField& magSI  = mesh_.magSf().internalField();
    const scalarField& deltaCoeffs = mesh_.deltaCoeffs().internalField();
    const scalarField& w = mesh_.weights().internalField();

    const vectorField& C  = mesh_.C().internalField();
    const vectorField& Cf  = mesh_.Cf().internalField();

    scalarField& diag = DEqn.diag();
    scalarField& upper = DEqn.upper();
    vectorField& source = DEqn.source();
    FieldField<Field, vector>& boundaryCoeffs = DEqn.boundaryCoeffs();

    vectorField& interDD = displacementIncrement();
    vectorField& interDT = tractionIncrement();

    if (curTimeIndex_ != mesh_.time().timeIndex())
    {
        // Update total fields
        ownSigmaf() += ownDSigmaf();
        ngbSigmaf() += ngbDSigmaf();

        displacement() += displacementIncrement();
        traction() += tractionIncrement();

        curTimeIndex_ = mesh_.time().timeIndex();
    }

    const vectorField& interD = displacement();
    const symmTensorField& ownSigmaf_ = ownSigmaf();
    const symmTensorField& ngbSigmaf_ = ngbSigmaf();

    symmTensorField& ownDSigmaf_ = ownDSigmaf();
    symmTensorField& ngbDSigmaf_ = ngbDSigmaf();

    // Looking up solid solver
    const solidSolver& stress =
        mesh_.objectRegistry::lookupObject<solidSolver>
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

        if (curFace < mesh_.nInternalFaces())
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
            label curPatch = mesh_.boundaryMesh().whichPatch(curFace);

            label curPatchFace =
                curFace - mesh_.boundaryMesh()[curPatch].start();

            const labelList& curPatchCells =
                mesh_.boundaryMesh()[curPatch].faceCells();

            scalar ngbDeltaN =
                mesh_.weights().boundaryField()[curPatch][curPatchFace]
               *(
                   1.0
                  /mesh_.deltaCoeffs().boundaryField()[curPatch][curPatchFace]
                );
            scalar ownDeltaN =
                (
                    1.0
                   /mesh_.deltaCoeffs().boundaryField()[curPatch][curPatchFace]
                )
              - ngbDeltaN;

            scalar magS =
                mesh_.magSf().boundaryField()[curPatch][curPatchFace];

            scalar ownMu  = mu.internalField()[curPatchCells[curPatchFace]];
            scalar ngbMu  = mu.boundaryField()[curPatch][curPatchFace];

            scalar ownLambda =
                lambda.internalField()[curPatchCells[curPatchFace]];
            scalar ngbLambda = lambda.boundaryField()[curPatch][curPatchFace];

            vector ownN =
                mesh_.Sf().boundaryField()[curPatch][curPatchFace]
               /mesh_.magSf().boundaryField()[curPatch][curPatchFace];
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
                mesh_.Cf().boundaryField()[curPatch][curPatchFace]
              - mesh_.C().internalField()[curPatchCells[curPatchFace]];
            ownCorrVec -= ownN*(ownN&ownCorrVec);

            vector ngbCorrVec =
                mesh_.Cf().boundaryField()[curPatch][curPatchFace]
              - mesh_.C().boundaryField()[curPatch][curPatchFace];
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
              - mesh_.weights().boundaryField()[curPatch][curPatchFace];

            scalar ownK = (2*ownMu + ownLambda);
            scalar ngbK = (2*ngbMu + ngbLambda);

            scalar Kf = 1.0/(wRevLin/ownK + (1.0-wRevLin)/ngbK);
            scalar muf = 1.0/(wRevLin/ownMu + (1.0-wRevLin)/ngbMu);

            scalar DeltaNf =
                1.0
               /mesh_.deltaCoeffs().boundaryField()[curPatch][curPatchFace];

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


void solidInterfaceITL::correct(surfaceVectorField& trac)
{
    const fvMesh& mesh_ = D_.mesh();

    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh_.nInternalFaces())
        {
            trac.internalField()[curFace] = traction()[faceI];
        }
        else
        {
            label curPatch = mesh_.boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh_.boundaryMesh()[curPatch].start();

            trac.boundaryField()[curPatch][curPatchFace] = traction()[faceI];
        }
    }
}


void solidInterfaceITL::modifyProperties
(
    surfaceScalarField& muf,
    surfaceScalarField& lambdaf
) const
{
    const fvMesh& mesh_ = D_.mesh();

    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh_.nInternalFaces())
        {
            muf.internalField()[curFace] = 0;
            lambdaf.internalField()[curFace] = 0;
        }
        else
        {
            label curPatch = mesh_.boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh_.boundaryMesh()[curPatch].start();

            muf.boundaryField()[curPatch][curPatchFace] = 0;
            lambdaf.boundaryField()[curPatch][curPatchFace] = 0;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

#include "solidInterfaceTL.H"
#include "fvc.H"
#include "processorFvPatchFields.H"
#include "fvMatrices.H"
#include "skewCorrectionVectors.H"
#include "leastSquaresGrad.H"
#include "gaussGrad.H"
#include "faceSet.H"
// #include "faMesh.H"
// #include "faCFD.H"
// #include "zeroGradientFaPatchFields.H"
#include "processorFvsPatchFields.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "solidSolver.H"
#include "fvcGradf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidInterfaceTL, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void solidInterfaceTL::makeFaces() const
{
    if (debug)
    {
        Info<< "void solidInterfaceTL::makeFaces() const : "
            << "creating list of interface faces"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (facesPtr_)
    {
        FatalErrorIn("solidInterfaceTL::makeFaces() const")
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

void solidInterfaceTL::makeDisplacement() const
{
    if (debug)
    {
        Info<< "void solidInterfaceTL::makeInterfaceDisplacement() const : "
            << "creating interface displacement field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (displacementPtr_)
    {
        FatalErrorIn("solidInterfaceTL::makeDisplacement() const")
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


void solidInterfaceTL::makeTraction() const
{
    if (debug)
    {
        Info<< "void solidInterfaceTL::makeTraction() const : "
            << "creating interface traction field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (tractionPtr_)
    {
        FatalErrorIn("solidInterfaceTL::makeTraction() const")
            << "interface traction field already exist"
            << abort(FatalError);
    }

    tractionPtr_ = new vectorField(faces().size(), vector::zero);
}


void solidInterfaceTL::makeSubMeshes() const
{
    if (debug)
    {
        Info<< "void solidInterfaceTL::makeSubMeshes() const : "
            << "creating material sub-meshes"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshes_.empty())
    {
        FatalErrorIn("solidInterfaceTL::makeSubMeshes() const")
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


void solidInterfaceTL::makeSubMeshVolToPoint() const
{
    if (debug)
    {
        Info<< "void solidInterfaceTL::makeVolToPointInterpolators() const : "
            << "creating cell-to-point interpolators"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshVolToPoint_.empty())
    {
        FatalErrorIn("solidInterfaceTL::makeVolToPointInterpolators() const")
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


void solidInterfaceTL::makeSubMeshPointD() const
{
    if (debug)
    {
        Info<< "void solidInterfaceTL::makeSubMeshPointD() const : "
            << "creating point displacements fields"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshPointD_.empty())
    {
        FatalErrorIn("solidInterfaceTL::makeSubMeshPointD() const")
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


void solidInterfaceTL::makeSubMeshD() const
{
    if (debug)
    {
        Info<< "void solidInterfaceTL::makeSubMeshD() const : "
            << "creating displacements fields"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshD_.empty())
    {
        FatalErrorIn("solidInterfaceTL::makeSubMeshD() const")
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


void solidInterfaceTL::makePointNumOfMaterials() const
{
    if (debug)
    {
        Info<< "void solidInterfaceTL::makePointNoMaterials() const : "
            << "creating number of materials for each point"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (pointNumOfMaterialsPtr_)
    {
        FatalErrorIn("solidInterfaceTL::makePointNumOfMaterials() const")
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


void solidInterfaceTL::clearOut()
{
    deleteDemandDrivenData(facesPtr_);
    deleteDemandDrivenData(displacementPtr_);
    deleteDemandDrivenData(tractionPtr_);
//     deleteDemandDrivenData(muPtr_);
//     deleteDemandDrivenData(lambdaPtr_);

//     subMeshes_.clear();
//     subMeshVolToPoint_.clear();

    deleteDemandDrivenData(pointNumOfMaterialsPtr_);
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


solidInterfaceTL::solidInterfaceTL
(
    const volVectorField& D,
    const pointVectorField& pointD
//     const constitutiveModel& rheology
)
:
    regIOobject
    (
        IOobject
        (
            "solidInterfaceTL",
            D.mesh().time().constant(),
            D.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    D_(D),
    pointD_(pointD),
//     rheology_(rheology),
    facesPtr_(NULL),
    displacementPtr_(NULL),
    tractionPtr_(NULL),
//     muPtr_(NULL),
//     lambdaPtr_(NULL),
    subMeshes_(0),
    subMeshVolToPoint_(0),
    subMeshD_(0),
    subMeshPointD_(0),
    pointNumOfMaterialsPtr_(NULL)
{
//     muPtr_ = new volScalarField(rheology_.mu());
//     lambdaPtr_ = new volScalarField(rheology_.lambda());
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

solidInterfaceTL::~solidInterfaceTL()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


const labelList& solidInterfaceTL::faces() const
{
    if (!facesPtr_)
    {
        makeFaces();
    }

    return *facesPtr_;
}

vectorField& solidInterfaceTL::displacement()
{
    if (!displacementPtr_)
    {
        makeDisplacement();
    }

    return *displacementPtr_;
}

const vectorField& solidInterfaceTL::displacement() const
{
    if (!displacementPtr_)
    {
        makeDisplacement();
    }

    return *displacementPtr_;
}

vectorField& solidInterfaceTL::traction()
{
    if (!tractionPtr_)
    {
        makeTraction();
    }

    return *tractionPtr_;
}

const vectorField& solidInterfaceTL::traction() const
{
    if (!tractionPtr_)
    {
        makeTraction();
    }

    return *tractionPtr_;
}

const PtrList<fvMeshSubset>& solidInterfaceTL::subMeshes() const
{
    if (subMeshes_.empty())
    {
        makeSubMeshes();
    }

    return subMeshes_;
}


const PtrList<leastSquaresVolPointInterpolation>& solidInterfaceTL
::subMeshVolToPoint() const
{
    if (subMeshVolToPoint_.empty())
    {
        makeSubMeshVolToPoint();
    }

    return subMeshVolToPoint_;
}


const PtrList<volVectorField>& solidInterfaceTL::subMeshD() const
{
    if (subMeshD_.empty())
    {
        makeSubMeshD();
    }

    return subMeshD_;
}


PtrList<volVectorField>& solidInterfaceTL::subMeshD()
{
    if (subMeshD_.empty())
    {
        makeSubMeshD();
    }

    return subMeshD_;
}


const PtrList<pointVectorField>& solidInterfaceTL::subMeshPointD() const
{
    if (subMeshPointD_.empty())
    {
        makeSubMeshPointD();
    }

    return subMeshPointD_;
}


PtrList<pointVectorField>& solidInterfaceTL::subMeshPointD()
{
    if (subMeshPointD_.empty())
    {
        makeSubMeshPointD();
    }

    return subMeshPointD_;
}


const labelList& solidInterfaceTL::pointNumOfMaterials() const
{
    if (!pointNumOfMaterialsPtr_)
    {
        makePointNumOfMaterials();
    }

    return *pointNumOfMaterialsPtr_;
}


void solidInterfaceTL::correct(fvVectorMatrix& DEqn)
{
    const fvMesh& mesh_ = D_.mesh();

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const vectorField& DI = D_.internalField();

    const volTensorField& gradD =
        mesh_.lookupObject<volTensorField>("grad(" + D_.name() + ')');
    const tensorField& gradDI = gradD.internalField();

    const surfaceTensorField& gradDf =
        mesh_.lookupObject<surfaceTensorField>("grad" + D_.name() + 'f');
    const tensorField& gradDfI = gradDf.internalField();

//     const volScalarField& mu = *muPtr_;
//     const volScalarField& lambda = *lambdaPtr_;
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

    vectorField& interD = displacement();
    vectorField& interT = traction();

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

    Switch thermalStress = stress.thermalStress();

    surfaceScalarField* DTfPtr = NULL;
    volScalarField* threeKPtr = NULL;
    volScalarField* alphaPtr = NULL;
    if (thermalStress)
    {
        DTfPtr =
            const_cast<surfaceScalarField*>
            (
                &mesh_.lookupObject<surfaceScalarField>("DTf")
            );
        threeKPtr =
            const_cast<volScalarField*>
            (
                &mesh_.lookupObject<volScalarField>("threeK")
            );
        alphaPtr =
            const_cast<volScalarField*>
            (
                &mesh_.lookupObject<volScalarField>("alpha")
            );
    }
    const surfaceScalarField& DTf = *DTfPtr;
    const volScalarField& threeK  = *threeKPtr;
    const volScalarField& alpha  = *alphaPtr;

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
                mesh_.Cf().boundaryField()[curPatch][curPatchFace]
              - mesh_.C().internalField()[curPatchCells[curPatchFace]];
            ownCorrVec -= ownN*(ownN&ownCorrVec);

            vector ngbCorrVec =
                mesh_.Cf().boundaryField()[curPatch][curPatchFace]
              - mesh_.C().boundaryField()[curPatch][curPatchFace];
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


void solidInterfaceTL::correct(surfaceVectorField& trac)
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


void solidInterfaceTL::updateDisplacement(pointVectorField& pointD)
{
    if (debug)
    {
        Info<< "solidInterfaceTL::volToPointInterpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    vectorField& pointDI = pointD.internalField();
    pointDI = vector::zero;

    const labelList& noMat = pointNumOfMaterials();

    const fvMesh& mesh_ = D_.mesh();

    const labelList& spLabels =
        mesh_.globalData().sharedPointLabels();

    const labelList& spAddressing =
        mesh_.globalData().sharedPointAddr();

    List<List<Map<vector> > > glData(Pstream::nProcs());
    forAll(glData, procI)
    {
        glData[procI] =
            List<Map<vector> >
            (
                mesh_.globalData().nGlobalPoints(),
                Map<vector>()
            );
    }

    forAll(subMeshes(), meshI)
    {
        // Update sub-mesh cell-centre displacement
        subMeshD()[meshI] = subMeshes()[meshI].interpolate(D_);

        // Correct displacement at the interface
        const labelList& patchMap = subMeshes()[meshI].patchMap();
        forAll(patchMap, patchI)
        {
            if (patchMap[patchI] == -1)
            {
                label interfacePatchIndex = patchI;

                vectorField& interfacePatchD =
                    subMeshD()[meshI].boundaryField()[interfacePatchIndex];

                const labelList& fm = subMeshes()[meshI].faceMap();

                label interfacePatchStart =
                    subMeshes()[meshI].subMesh().boundaryMesh()
                    [
                        interfacePatchIndex
                    ].start();

                forAll(interfacePatchD, faceI)
                {
                    label curInterFace =
                        findIndex(faces(), fm[interfacePatchStart + faceI]);

                    interfacePatchD[faceI] = displacement()[curInterFace];
                }
            }
        }

        // Calc point sub-mesh point displacement
        subMeshVolToPoint()[meshI].interpolate
            (
                subMeshD()[meshI],
                subMeshPointD()[meshI]
            );

        const vectorField& subMeshPointDI =
            subMeshPointD()[meshI].internalField();

        // Map point field from sub-mesh to global mesh
        const labelList& pointMap = subMeshes()[meshI].pointMap();
        forAll(pointMap, pointI)
        {
            label curMeshPoint = pointMap[pointI];

            bool sharedPoint(findIndex(spLabels, curMeshPoint) != -1);

            if (sharedPoint)
            {
                label k = findIndex(spLabels, curMeshPoint);
                label curSpIndex = spAddressing[k];
                glData[Pstream::myProcNo()][curSpIndex].insert
                    (
                        meshI,
                        subMeshPointDI[pointI]
                    );
            }
            else
            {
                pointDI[curMeshPoint] +=
                    subMeshPointDI[pointI]/noMat[curMeshPoint];
            }
        }
    }

    Pstream::gatherList(glData);
    Pstream::scatterList(glData);

    // Gloabal points
    if (mesh_.globalData().nGlobalPoints())
    {
        for (label k=0; k<mesh_.globalData().nGlobalPoints(); k++)
        {
            label curSpIndex = findIndex(spAddressing, k);

            if (curSpIndex != -1)
            {
                List<label> matN(subMeshes().size(), 0);
                List<vector> matAvg(subMeshes().size(), vector::zero);

                forAll(glData, procI)
                {
                    const Map<vector>& curProcGlData = glData[procI][k];
                    for (label i=0; i<subMeshes().size(); i++)
                    {
                        if (curProcGlData.found(i))
                        {
                            matAvg[i] += curProcGlData[i];
                            matN[i]++;
                        }
                    }
                }

                label nMat = 0;
                vector avg = vector::zero;
                forAll(matAvg, matI)
                {
                    if (matN[matI])
                    {
                        matAvg[matI] /= matN[matI];
                        avg += matAvg[matI];
                        nMat++;
                    }
                }
                avg /= nMat;

                label curMeshPoint = spLabels[curSpIndex];
                pointDI[curMeshPoint] = avg;
            }
        }
    }

    pointD.correctBoundaryConditions();
}

void solidInterfaceTL::updateDisplacementGradient
(
    volTensorField& gradD,
    surfaceTensorField& gradDf
)
{
    tensorField& gradDI = gradD.internalField();
    tensorField& gradDfI = gradDf.internalField();
    gradDfI = tensor::zero;

    const fvMesh& mesh = D_.mesh();

    // Calculate displacement gradient for sub-meshes
    forAll(subMeshes(), meshI)
    {
        volTensorField subMeshGradD =
            fvc::grad(subMeshD()[meshI], subMeshPointD()[meshI]);

        // Face gradient in tangential direction
        surfaceTensorField subMeshGradDf =
            fvc::fsGrad(subMeshD()[meshI], subMeshPointD()[meshI]);

        const tensorField& subMeshGradDI = subMeshGradD.internalField();
        const tensorField& subMeshGradDfI = subMeshGradDf.internalField();

        // Map iternal field from sub-mesh to global mesh
        const labelList& cellMap = subMeshes()[meshI].cellMap();
        forAll(subMeshGradDI, cellI)
        {
            label curGlobalMeshCell = cellMap[cellI];

            gradDI[curGlobalMeshCell] = subMeshGradDI[cellI];
        }

        const labelList& faceMap = subMeshes()[meshI].faceMap();
        forAll(subMeshGradDfI, faceI)
        {
            label curGlobalMeshFace = faceMap[faceI];

            gradDfI[curGlobalMeshFace] = subMeshGradDfI[faceI];
        }

        // Map boundary field
        const labelList& patchMap = subMeshes()[meshI].patchMap();
        forAll(subMeshGradD.boundaryField(), patchI)
        {
            const fvPatchField<tensor>& subMeshPatchGradD =
                subMeshGradD.boundaryField()[patchI];

            const fvsPatchField<tensor>& subMeshPatchGradDf =
                subMeshGradDf.boundaryField()[patchI];

            label start = subMeshPatchGradD.patch().patch().start();

            if (patchMap[patchI] != -1)
            {
                fvPatchField<tensor>& patchGradD =
                    gradD.boundaryField()[patchMap[patchI]];

                if (!patchGradD.coupled())
                {
                    forAll(subMeshPatchGradD, faceI)
                    {
                        label globalGlobalMeshFace = faceMap[start+faceI];

                        label curGlobalMeshPatchFace =
                            globalGlobalMeshFace
                          - mesh.boundaryMesh()[patchMap[patchI]].start();

                        patchGradD[curGlobalMeshPatchFace] =
                            subMeshPatchGradD[faceI];
                    }
                }

                fvsPatchField<tensor>& patchGradDf =
                    gradDf.boundaryField()[patchMap[patchI]];

                forAll(subMeshPatchGradDf, faceI)
                {
                    label globalGlobalMeshFace = faceMap[start+faceI];

                    label curGlobalMeshPatchFace =
                        globalGlobalMeshFace
                      - mesh.boundaryMesh()[patchMap[patchI]].start();

                    patchGradDf[curGlobalMeshPatchFace] =
                        subMeshPatchGradDf[faceI];
                }
            }
            else // interface faces
            {
                forAll(subMeshPatchGradDf, faceI)
                {
                    label globalGlobalMeshFace = faceMap[start+faceI];

                    if (globalGlobalMeshFace < mesh.nInternalFaces())
                    {
                        gradDfI[globalGlobalMeshFace] +=
                            0.5*subMeshPatchGradDf[faceI];
                    }
                    else
                    {
                        label curPatch =
                            mesh.boundaryMesh().whichPatch
                            (
                                globalGlobalMeshFace
                            );
                        label curPatchFace =
                            globalGlobalMeshFace
                          - mesh.boundaryMesh()[curPatch].start();

                        gradDf.boundaryField()[curPatch][curPatchFace] =
                            subMeshPatchGradDf[faceI];
                    }
                }
            }
        }
    }

    // Correct cell gradient at the boundary
    fv::gaussGrad<vector>(mesh).correctBoundaryConditions(D_, gradD);

    // Make sure face gradient is consistent accross processor patches
    forAll(gradDf.boundaryField(), patchI)
    {
        fvsPatchField<tensor>& patchGradDf =
            gradDf.boundaryField()[patchI];

        if (patchGradDf.type() == processorFvsPatchTensorField::typeName)
        {
            const tensorField& patchField = patchGradDf;

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh.boundaryMesh()[patchI]
                );

            OPstream::write
            (
                Pstream::blocking,
                procPatch.neighbProcNo(),
                reinterpret_cast<const char*>(patchField.begin()),
                patchField.byteSize()
            );
        }
    }
    forAll(gradDf.boundaryField(), patchI)
    {
        fvsPatchField<tensor>& patchGradDf =
            gradDf.boundaryField()[patchI];

        if (patchGradDf.type() == processorFvsPatchTensorField::typeName)
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh.boundaryMesh()[patchI]
                );

            tensorField ngbPatchField(procPatch.size(), tensor::zero);

            IPstream::read
            (
                Pstream::blocking,
                procPatch.neighbProcNo(),
                reinterpret_cast<char*>(ngbPatchField.begin()),
                ngbPatchField.byteSize()
            );

            tensorField& patchField = patchGradDf;

            patchField = 0.5*(patchField + ngbPatchField);
        }
    }

    // Add normal gradient
    gradDf += mesh.Sf()*fvc::snGrad(D_)/mesh.magSf();
}

void solidInterfaceTL::modifyProperties
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


void solidInterfaceTL::modifyProperty
(
    surfaceScalarField& muf
) const
{
    const fvMesh& mesh = D_.mesh();

    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh.nInternalFaces())
        {
            muf.internalField()[curFace] = 0;
        }
        else
        {
            label curPatch = mesh.boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh.boundaryMesh()[curPatch].start();

            muf.boundaryField()[curPatch][curPatchFace] = 0;
        }
    }
}


tmp<volTensorField> solidInterfaceTL::grad(volVectorField& D) const
{
    const fvMesh& mesh_ = D_.mesh();

    tmp<volTensorField> tGradD
    (
        new volTensorField
        (
            IOobject
            (
                "grad(" + D.name() + ')',
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor("zero", dimless, tensor::zero)
        )
    );
    volTensorField& gradD = tGradD();

    surfaceVectorField Df = fvc::interpolate(D);


    // Skew-correction

    if (skewCorrectionVectors::New(mesh_).skew())
    {
        const skewCorrectionVectors& scv = skewCorrectionVectors::New(mesh_);

        const volTensorField& gradD =
            mesh_.lookupObject<volTensorField>("grad(" + D.name() + ')');

        Df +=
        (
            scv()
          & linear<tensor>(mesh_).interpolate
            (
                gradD
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

    const vectorField& interD = displacement();

    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh_.nInternalFaces())
        {
            Df.internalField()[curFace] = interD[faceI];
        }
        else
        {
            label curPatch = mesh_.boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh_.boundaryMesh()[curPatch].start();

            Df.boundaryField()[curPatch][curPatchFace] = interD[faceI];
        }
    }

//     forAll(processorPatchFaces(), patchI)
//     {
//         label curPatch = processorPatches()[patchI];

//         const vectorField& curProcInterU =
//             processorInterfaceDisplacement()[patchI];

//         forAll(processorPatchFaces()[patchI], faceI)
//         {
//             label curFace = processorPatchFaces()[patchI][faceI];

//             Uf.boundaryField()[curPatch][curFace] = curProcInterU[faceI];
//         }
//     }


    // Gradient calculation using Gauss method

    gradD = fv::gaussGrad<vector>(mesh_).gradf(Df, "grad(D)");
    fv::gaussGrad<vector>(mesh_).correctBoundaryConditions(D, gradD);

    return tGradD;
}


tmp<symmTensorField> solidInterfaceTL::sigmaA() const
{
    const fvMesh& mesh_ = D_.mesh();

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const vectorField& DI = D_.internalField();

    const volTensorField& gradD =
        mesh_.lookupObject<volTensorField>("grad(" + D_.name() + ')');
    const tensorField& gradDI = gradD.internalField();

    const surfaceTensorField& gradDf =
        mesh_.lookupObject<surfaceTensorField>("grad" + D_.name() + 'f');
    const tensorField& gradDfI = gradDf.internalField();

    const volScalarField& mu =
        mesh_.lookupObject<volScalarField>("mu");
    const volScalarField& lambda =
        mesh_.lookupObject<volScalarField>("lambda");
//     const volScalarField& mu = *muPtr_;
//     const volScalarField& lambda = *lambdaPtr_;

    const vectorField& interD = displacement();

    const vectorField& SI  = mesh_.Sf().internalField();
    const scalarField& magSI  = mesh_.magSf().internalField();
    const scalarField& deltaCoeffs = mesh_.deltaCoeffs().internalField();
    const scalarField& w = mesh_.weights().internalField();

    const vectorField& C  = mesh_.C().internalField();
    const vectorField& Cf  = mesh_.Cf().internalField();

    // Looking up solid solver
    const solidSolver& stress =
        mesh_.objectRegistry::lookupObject<solidSolver>
        (
            "solidProperties"
        );

    Switch nonLinear
    (
        stress.lookup("nonLinear")
    );

    Switch enforceLinear
    (
        stress.lookup("enforceLinear")
    );

    tmp<symmTensorField> tSigmaA
    (
        new symmTensorField(faces().size(), symmTensor::zero)
    );
    symmTensorField& sigmaA = tSigmaA();

    if (mesh_.foundObject<volScalarField>("materials"))
    {
        const volScalarField& materials =
            mesh_.lookupObject<volScalarField>("materials");

        symmTensorField EA(sigmaA.size(), symmTensor::zero);
        scalarField muA(sigmaA.size(), 0);
        scalarField lambdaA(sigmaA.size(), 0);

        forAll(faces(), faceI)
        {
            label curFace = faces()[faceI];

            if (curFace < mesh_.nInternalFaces())
            {
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

                vector corrVec = Cf[curFace] - C[curCell];
                corrVec -= n*(n&corrVec);
                vector DCorr = (corrVec&gradDI[curCell]);

                tensor gradDA = n*(interD[faceI] - (DI[curCell]+DCorr))/dn;
                gradDA += ((I-n*n)&gradDfI[curFace]);

                EA[faceI] = symm(gradDA);
                if (nonLinear && (!enforceLinear))
                {
                    EA[faceI] += 0.5*symm(gradDA & gradDA.T());
                }
            }
        }

        sigmaA = 2*muA*EA + lambdaA*(I*tr(EA));
    }

    return tSigmaA;
}


tmp<symmTensorField> solidInterfaceTL::sigmaB() const
{
    const fvMesh& mesh_ = D_.mesh();

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const vectorField& DI = D_.internalField();

    const volTensorField& gradD =
        mesh_.lookupObject<volTensorField>("grad(" + D_.name() + ')');
    const tensorField& gradDI = gradD.internalField();

    const surfaceTensorField& gradDf =
        mesh_.lookupObject<surfaceTensorField>("grad" + D_.name() + 'f');
    const tensorField& gradDfI = gradDf.internalField();

    const volScalarField& mu =
        mesh_.lookupObject<volScalarField>("mu");
    const volScalarField& lambda =
        mesh_.lookupObject<volScalarField>("lambda");
//     const volScalarField& mu = *muPtr_;
//     const volScalarField& lambda = *lambdaPtr_;

    const vectorField& interD = displacement();

    const vectorField& SI  = mesh_.Sf().internalField();
    const scalarField& magSI  = mesh_.magSf().internalField();
    const scalarField& deltaCoeffs = mesh_.deltaCoeffs().internalField();
    const scalarField& w = mesh_.weights().internalField();

    const vectorField& C  = mesh_.C().internalField();
    const vectorField& Cf  = mesh_.Cf().internalField();

    // Looking up solid solver
    const solidSolver& stress =
        mesh_.objectRegistry::lookupObject<solidSolver>
        (
            "solidProperties"
        );

    Switch nonLinear
    (
        stress.lookup("nonLinear")
    );

    Switch enforceLinear
    (
        stress.lookup("enforceLinear")
    );

    tmp<symmTensorField> tSigmaB
    (
        new symmTensorField(faces().size(), symmTensor::zero)
    );
    symmTensorField& sigmaB = tSigmaB();

    if (mesh_.foundObject<volScalarField>("materials"))
    {
        const volScalarField& materials =
            mesh_.lookupObject<volScalarField>("materials");

        symmTensorField EB(sigmaB.size(), symmTensor::zero);
        scalarField muB(sigmaB.size(), 0);
        scalarField lambdaB(sigmaB.size(), 0);

        forAll(faces(), faceI)
        {
            label curFace = faces()[faceI];

            if (curFace < mesh_.nInternalFaces())
            {
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

                vector corrVec = Cf[curFace] - C[curCell];
                corrVec -= n*(n&corrVec);
                vector DCorr = (corrVec&gradDI[curCell]);

                tensor gradDB = n*(interD[faceI] - (DI[curCell]+DCorr))/dn;
                gradDB += ((I-n*n)&gradDfI[curFace]);

                EB[faceI] = symm(gradDB);
                if (nonLinear && (!enforceLinear))
                {
                    EB[faceI] += 0.5*symm(gradDB & gradDB.T());
                }
            }
        }

        sigmaB = 2*muB*EB + lambdaB*(I*tr(EB));
    }

    return tSigmaB;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

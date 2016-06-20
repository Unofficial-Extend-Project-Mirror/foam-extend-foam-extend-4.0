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

Class
    solidInterface

Description
    Material rheology for solids.

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
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidInterface, 0);
defineRunTimeSelectionTable(solidInterface, dictionary);


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
            if (!interfaceCellSet.found(curCells[cellI]))
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
        //faMeshSet.write();
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
        Info<< "void solidInterface::"
            << "makeProcessorInterfaceDisplacement() const : "
            << "creating processor inter-faces displacement"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (processorInterfaceUPtr_)
    {
        FatalErrorIn
            ("solidInterface::makeProcessorInterfaceDisplacement() const")
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

void solidInterface::makeIndicatorFieldMap() const
{
  if (indicatorFieldMapPtr_)
    {
      FatalErrorIn("solidInterface::makeIndicatorFieldMap() const")
    << "interface indicator field map already exists"
    << abort(FatalError);
    }

    indicatorFieldMapPtr_  =
      new labelList
      (
       mesh_.owner().size(),
       -1
       );
    labelList& indicatorFieldMap = *indicatorFieldMapPtr_;

    forAll(globalInterFaces(), faceI)
      {
        label curGlobalFace = globalInterFaces()[faceI];

    indicatorFieldMap[curGlobalFace] = faceI;
      }

}

void solidInterface::makeProcessorPatchFacesMap() const
{
  if (processorPatchFacesMapPtr_)
    {
      FatalErrorIn("solidInterface::makeProcessorPatchFacesMap() const")
    << "processor interface indicator field map already exists"
    << abort(FatalError);
    }

    processorPatchFacesMapPtr_  =
      new labelListList
      (
       mesh_.boundary().size(),
       labelList(0,-1)
       );
    processorPatchMapPtr_  =
      new labelList
      (
       mesh_.boundary().size(),
       -1
       );
    labelListList& processorPatchFacesMap = *processorPatchFacesMapPtr_;
    labelList& processorPatchMap = *processorPatchMapPtr_;
    forAll(mesh_.boundary(), patchI)
      {
    processorPatchFacesMap[patchI].setSize
      (
       mesh_.boundary()[patchI].size(),
       -1
       );
      }

    forAll(processorPatchFaces(), patchI)
      {
    label curPatch = processorPatches()[patchI];

    processorPatchMap[curPatch] = patchI;

        forAll(processorPatchFaces()[patchI], faceI)
      {
            label curFace = processorPatchFaces()[patchI][faceI];
        processorPatchFacesMap[curPatch][curFace] = faceI;
      }
    }
}

  // this is a public member now to allow data to be cleared
  // when there are topological changes
void solidInterface::clearOut()
{
    deleteDemandDrivenData(subMeshPtr_);
    deleteDemandDrivenData(globalInterFacesPtr_);
    deleteDemandDrivenData(localInterFacesPtr_);
    deleteDemandDrivenData(interfaceUPtr_);
    deleteDemandDrivenData(processorPatchesPtr_);
    deleteDemandDrivenData(processorPatchFacesPtr_);
    deleteDemandDrivenData(processorInterfaceUPtr_);
    deleteDemandDrivenData(indicatorPtr_);
    deleteDemandDrivenData(indicatorFieldMapPtr_);
    deleteDemandDrivenData(processorPatchMapPtr_);
    deleteDemandDrivenData(processorPatchFacesMapPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidInterface::solidInterface
(
    const word& name,
    const fvMesh& mesh,
    const constitutiveModel& rheology
)
:
  name_(name),
  mesh_(mesh),
  rheology_(rheology),
  subMeshPtr_(NULL),
  globalInterFacesPtr_(NULL),
  localInterFacesPtr_(NULL),
  interfaceUPtr_(NULL),
  processorPatchesPtr_(NULL),
  processorPatchFacesPtr_(NULL),
  processorInterfaceUPtr_(NULL),
  indicatorPtr_(NULL),
  indicatorFieldMapPtr_(NULL),
  processorPatchMapPtr_(NULL),
  processorPatchFacesMapPtr_(NULL)
{}


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

void solidInterface::modifyProperties
(
    surfaceScalarField& s
) const
{
    forAll(globalInterFaces(), faceI)
    {
        label curGlobalFace = globalInterFaces()[faceI];

        s.internalField()[curGlobalFace] = 0;
    }

    forAll(processorPatchFaces(), patchI)
    {
        label curPatch = processorPatches()[patchI];

        forAll(processorPatchFaces()[patchI], faceI)
        {
            label curFace = processorPatchFaces()[patchI][faceI];

            s.boundaryField()[curPatch][curFace] = 0;
        }
    }
}

void solidInterface::modifyProperties
(
 surfaceSymmTensor4thOrderField& st
 ) const
{
    forAll(globalInterFaces(), faceI)
    {
        label curGlobalFace = globalInterFaces()[faceI];

        st.internalField()[curGlobalFace] = symmTensor4thOrder::zero;
    }

    forAll(processorPatchFaces(), patchI)
    {
        label curPatch = processorPatches()[patchI];

        forAll(processorPatchFaces()[patchI], faceI)
        {
            label curFace = processorPatchFaces()[patchI][faceI];

            st.boundaryField()[curPatch][curFace] = symmTensor4thOrder::zero;
        }
    }
}

void solidInterface::modifyProperties
(
 surfaceDiagTensorField& dt
 ) const
{
    forAll(globalInterFaces(), faceI)
    {
        label curGlobalFace = globalInterFaces()[faceI];

        dt.internalField()[curGlobalFace] = diagTensor::zero;
    }

    forAll(processorPatchFaces(), patchI)
    {
        label curPatch = processorPatches()[patchI];

        forAll(processorPatchFaces()[patchI], faceI)
        {
            label curFace = processorPatchFaces()[patchI][faceI];

            dt.boundaryField()[curPatch][curFace] = diagTensor::zero;
        }
    }
}

void solidInterface::modifyProperties
(
    surfaceScalarField& mu,
    surfaceScalarField& lambda
) const
{
    forAll(globalInterFaces(), faceI)
    {
        label curGlobalFace = globalInterFaces()[faceI];

        mu.internalField()[curGlobalFace] = 0;
        lambda.internalField()[curGlobalFace] = 0;
    }

    forAll(processorPatchFaces(), patchI)
    {
        label curPatch = processorPatches()[patchI];

        forAll(processorPatchFaces()[patchI], faceI)
        {
            label curFace = processorPatchFaces()[patchI][faceI];

            mu.boundaryField()[curPatch][curFace] = 0;
            lambda.boundaryField()[curPatch][curFace] = 0;
        }
    }
}

void solidInterface::modifyProperties
(
    surfaceScalarField& mu,
    surfaceScalarField& lambda,
    surfaceScalarField& threeKalpha
) const
{
    forAll(globalInterFaces(), faceI)
    {
        label curGlobalFace = globalInterFaces()[faceI];

        mu.internalField()[curGlobalFace] = 0;
        lambda.internalField()[curGlobalFace] = 0;
        threeKalpha.internalField()[curGlobalFace] = 0;
    }

    forAll(processorPatchFaces(), patchI)
    {
        label curPatch = processorPatches()[patchI];

        forAll(processorPatchFaces()[patchI], faceI)
        {
            label curFace = processorPatchFaces()[patchI][faceI];

            mu.boundaryField()[curPatch][curFace] = 0;
            lambda.boundaryField()[curPatch][curFace] = 0;
        threeKalpha.boundaryField()[curPatch][curFace] = 0;
        }
    }
}


void solidInterface::modifyProperties
(
    surfaceSymmTensor4thOrderField& C,
    surfaceDiagTensorField& K
) const
{
    forAll(globalInterFaces(), faceI)
    {
        label curGlobalFace = globalInterFaces()[faceI];

        C.internalField()[curGlobalFace] = symmTensor4thOrder::zero;
        K.internalField()[curGlobalFace] = diagTensor::zero;
    }

    forAll(processorPatchFaces(), patchI)
    {
        label curPatch = processorPatches()[patchI];

        forAll(processorPatchFaces()[patchI], faceI)
        {
            label curFace = processorPatchFaces()[patchI][faceI];

        C.boundaryField()[curPatch][curFace] = symmTensor4thOrder::zero;
        K.boundaryField()[curPatch][curFace] = diagTensor::zero;
        }
    }
}


// tmp<volTensorField> solidInterface::grad(volVectorField& U) const
// {
//     tmp<volTensorField> tGradU
//     (
//         new volTensorField
//         (
//             IOobject
//             (
//                 "grad(" + U.name() + ')',
//                 mesh_.time().timeName(),
//                 mesh_,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             mesh_,
//             dimensionedTensor("zero", dimless, tensor::zero)
//         )
//     );
//     volTensorField& gradU = tGradU();

//     surfaceVectorField Uf = fvc::interpolate(U);


//     // Skew-correction

//     if (skewCorrectionVectors::New(this->mesh_).skew())
//     {
//         const skewCorrectionVectors& scv = skewCorrectionVectors::New(mesh_);

//         const volTensorField& gradU =
//             mesh_.lookupObject<volTensorField>("grad(" + U.name() + ')');

//         Uf +=
//         (
//             scv()
//           & linear<tensor>(mesh_).interpolate
//             (
//                 gradU
//             )
//         );

// //         Uf +=
// //         (
// //             scv()
// //           & linear<tensor>(mesh_).interpolate
// //             (
// //                 fv::leastSquaresGrad<vector>(mesh_).grad(U)
// //             )
// //         );
//     }


//     // Interface correction

//     const vectorField& interU = interfaceDisplacement();

//     forAll(globalInterFaces(), faceI)
//     {
//         label curGlobalFace = globalInterFaces()[faceI];

//         Uf.internalField()[curGlobalFace] = interU[faceI];
//     }

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


//     // Gradient calculation using Gauss method

//     gradU = fv::gaussGrad<vector>(mesh_).grad(Uf);
//     fv::gaussGrad<vector>(mesh_).correctBoundaryConditions(U, gradU);


//     return tGradU;
// }

const List<labelPair>& solidInterface::indicator() const
{
    if (!indicatorPtr_)
    {
        makeIndicator();
    }

    return *indicatorPtr_;
}

const labelList& solidInterface::indicatorFieldMap() const
{
    if (!indicatorFieldMapPtr_)
    {
        makeIndicatorFieldMap();
    }

    return *indicatorFieldMapPtr_;
}

const labelList& solidInterface::processorPatchMap() const
{
    if (!processorPatchMapPtr_)
    {
        makeProcessorPatchFacesMap();
    }

    return *processorPatchMapPtr_;
}

const labelListList& solidInterface::processorPatchFacesMap() const
{
    if (!processorPatchFacesMapPtr_)
    {
        makeProcessorPatchFacesMap();
    }

    return *processorPatchFacesMapPtr_;
}


} // End namespace Foam

// ************************************************************************* //

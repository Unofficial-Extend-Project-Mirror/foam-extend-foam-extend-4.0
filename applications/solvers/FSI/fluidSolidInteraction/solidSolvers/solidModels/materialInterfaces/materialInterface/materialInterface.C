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

#include "materialInterface.H"
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
#include "fvcCellLimitedGrad.H"
#include "componentMixedPointPatchFields.H"
#include "processorPointPatchFields.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(materialInterface, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void materialInterface::makeFaces() const
{
    if (debug)
    {
        Info<< "void materialInterface::makeFaces() const : "
            << "creating list of interface faces"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (facesPtr_)
    {
        FatalErrorIn("materialInterface::makeFaces() const")
            << "list of interface faces already exists"
            << abort(FatalError);
    }

    if (mesh().foundObject<volScalarField>("materials"))
    {
        const volScalarField& materials =
            mesh().lookupObject<volScalarField>("materials");

        const scalarField& materialsI = materials.internalField();

        const unallocLabelList& owner = mesh().owner();
        const unallocLabelList& neighbour = mesh().neighbour();

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
            if (mesh().boundary()[patchI].type() == processorFvPatch::typeName)
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
                            mesh().boundaryMesh()[patchI].start() + faceI;

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


void materialInterface::makeSubMeshes() const
{
    if (debug)
    {
        Info<< "void materialInterface::makeSubMeshes() const : "
            << "creating material sub-meshes"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshes_.empty())
    {
        FatalErrorIn("materialInterface::makeSubMeshes() const")
            << "material sub-meshes already exist"
            << abort(FatalError);
    }

    const volScalarField& materials =
        mesh().lookupObject<volScalarField>("materials");
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
                    mesh().time().constant(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh()
            )
        );

        subMeshes_[matI].setLargeCellSubset(region, matI);
    }
}


void materialInterface::makeSubMeshVolToPoint() const
{
    if (debug)
    {
        Info<< "void materialInterface::makeVolToPointInterpolators() const : "
            << "creating cell-to-point interpolators"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!subMeshVolToPoint_.empty())
    {
        FatalErrorIn("materialInterface::makeVolToPointInterpolators() const")
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


void materialInterface::makePointNumOfMaterials() const
{
    if (debug)
    {
        Info<< "void materialInterface::makePointNoMaterials() const : "
            << "creating number of materials for each point"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (pointNumOfMaterialsPtr_)
    {
        FatalErrorIn("materialInterface::makePointNumOfMaterials() const")
            << "Point numbero of materials already exist"
            << abort(FatalError);
    }

    pointNumOfMaterialsPtr_ = new labelList(mesh().nPoints(), 0);
    labelList& pointNumOfMaterials = *pointNumOfMaterialsPtr_;

    const labelListList& pointCells = mesh().pointCells();

    const volScalarField& materials =
        mesh().lookupObject<volScalarField>("materials");
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


void materialInterface::makeIsolatedInterfacePoints() const
{
    if (debug)
    {
        Pout<< "void materialInterface::makeIsolatedInterfacePoints() const : "
            << "creating isolated interface points"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (isolatedInterfacePointsPtr_)
    {
        FatalErrorIn("materialInterface::makeIsolatedInterfacePoints() const")
            << "Isolated interface points already exist"
            << abort(FatalError);
    }

    const volScalarField& materials =
        mesh().lookupObject<volScalarField>("materials");
    const scalarField& materialsI = materials.internalField();

    pointMesh pMesh(mesh());

    pointScalarField pointMaterials
    (
        IOobject
        (
            "pointMaterials",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh,
        dimensionedScalar("0", dimless, 0)
    );
    scalarField& pointMaterialsI = pointMaterials.internalField();

    const labelListList& pointCells = mesh().pointCells();

    forAll(pointMaterialsI, pointI)
    {
        const labelList& curPointCells = pointCells[pointI];

        forAll(curPointCells, cellI)
        {
            label curCell = curPointCells[cellI];
            pointMaterialsI[pointI] += materialsI[curCell];
        }

        pointMaterialsI[pointI] /= curPointCells.size() + SMALL;
    }

    pointMaterials.correctBoundaryConditions();

    scalarField matInter(pointMaterialsI.size(), 0);

//     label count = 0;
    forAll(pointMaterialsI, pointI)
    {
        const labelList& curPointCells = pointCells[pointI];

        if
        (
            mag(pointMaterialsI[pointI] - materialsI[curPointCells[0]])
          > SMALL
        )
        {
            matInter[pointI] = 1;
//             count++;
        }
    }

    labelHashSet isolatedPointsSet;

    const labelList& noMat = pointNumOfMaterials();

    const labelList& spLabels =
        mesh().globalData().sharedPointLabels();

    const labelListList& pointFaces = mesh().pointFaces();
    forAll(matInter, pointI)
    {
        if (matInter[pointI] && (noMat[pointI] == 1))
//         if (materialInterface_[pointI] && (noMat[pointI] == 1))
        {
            bool sharedPoint(findIndex(spLabels, pointI) != -1);

            if (!sharedPoint)
            {
                bool hasProcessorFace = false;

                const labelList& curPointFaces = pointFaces[pointI];
                forAll(curPointFaces, faceI)
                {
                    label faceID = curPointFaces[faceI];
                    label patchID =
                        mesh().boundaryMesh().whichPatch(faceID);

                    if (patchID != -1)
                    {
                        if
                        (
                            isA<processorPolyPatch>
                            (
                                mesh().boundaryMesh()[patchID]
                            )
                        )
                        {
                            if (findIndex(faces(), faceID) == -1)
                            {
                                hasProcessorFace = true;
                                break;
                            }
                        }
                    }
//                     if (findIndex(faces(), curPointFaces[faceI]) != -1)
//                     {
//                         noInterfaceFaces = false;
//                         break;
//                     }
                }

                if (hasProcessorFace)
                {
                    isolatedPointsSet.insert(pointI);
                }
            }
        }
    }

    isolatedInterfacePointsPtr_ = new labelList(isolatedPointsSet.toc());

//     Pout << "Num of material interface points: "
//         << sum(materialInterface_) << "/"
//         << materialInterface_.size() << endl;

//     Pout << "Isolated interface points: "
//         << isolatedPointsSet.toc() << endl;
}


void materialInterface::clearOut()
{
    deleteDemandDrivenData(facesPtr_);

    subMeshVolToPoint_.clear();
    subMeshes_.clear();

    deleteDemandDrivenData(pointNumOfMaterialsPtr_);
    deleteDemandDrivenData(isolatedInterfacePointsPtr_);
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


materialInterface::materialInterface
(
    const fvMesh& mesh
)
:
    regIOobject
    (
        IOobject
        (
            "materialInterface",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
//     pMesh_(mesh),
//     materialInterface_
//     (
//         IOobject
//         (
//             "materialInterface",
//             mesh.time().timeName(),
//             mesh,
//             IOobject::MUST_READ,
//             IOobject::AUTO_WRITE
//         ),
//         pMesh_
//     ),
    facesPtr_(NULL),
    subMeshes_(0),
    subMeshVolToPoint_(0),
    pointNumOfMaterialsPtr_(NULL),
    isolatedInterfacePointsPtr_(NULL)
{
    isolatedInterfacePoints();
// Pout << isolatedInterfacePoints() << endl;
//     Pout << vectorField(mesh_.points(), isolatedInterfacePoints()) << endl;
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

materialInterface::~materialInterface()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


const labelList& materialInterface::faces() const
{
    if (!facesPtr_)
    {
        makeFaces();
    }

    return *facesPtr_;
}


const PtrList<fvMeshSubset>& materialInterface::subMeshes() const
{
    if (subMeshes_.empty())
    {
        makeSubMeshes();
    }

    return subMeshes_;
}


const PtrList<leastSquaresVolPointInterpolation>& materialInterface
::subMeshVolToPoint() const
{
    if (subMeshVolToPoint_.empty())
    {
        makeSubMeshVolToPoint();
    }

    return subMeshVolToPoint_;
}


const labelList& materialInterface::pointNumOfMaterials() const
{
    if (!pointNumOfMaterialsPtr_)
    {
        makePointNumOfMaterials();
    }

    return *pointNumOfMaterialsPtr_;
}


const labelList& materialInterface::isolatedInterfacePoints() const
{
    if (!isolatedInterfacePointsPtr_)
    {
        makeIsolatedInterfacePoints();
    }

    return *isolatedInterfacePointsPtr_;
}


void materialInterface::updateDisplacement
(
    const volVectorField& displacement,
    const vectorField& interfaceDisplacement,
    pointVectorField& pointDisplacement,
    PtrList<volVectorField>& subMeshDisplacement,
    PtrList<pointVectorField>& subMeshPointDisplacement
)
{
    if (debug)
    {
        Info<< "materialInterface::upateDisplacement("
            << "PtrList<volVectorField>& subMeshDisplacement, "
            << "PtrList<pointVectorField>& subMeshPointDisplacement, "
            << "volVectorField& displacement, "
            << "pointVectorField& pointDisplacement, "
            << "const vectorField& interfaceDisplacement"
            << ") : "
            << "interpolating fields from cells to points"
            << endl;
    }

    vectorField& pointDI = pointDisplacement.internalField();
    pointDI = vector::zero;

    const labelList& noMat = pointNumOfMaterials();

    const labelList& spLabels =
        mesh().globalData().sharedPointLabels();

    const labelList& spAddressing =
        mesh().globalData().sharedPointAddr();

    List<List<Map<vector> > > glData(Pstream::nProcs());
    forAll(glData, procI)
    {
        glData[procI] =
            List<Map<vector> >
            (
                mesh().globalData().nGlobalPoints(),
                Map<vector>()
            );
    }

    forAll(subMeshes(), meshI)
    {
        // Update sub-mesh cell-centre displacement
        subMeshDisplacement[meshI] ==
            subMeshes()[meshI].interpolate(displacement);

        // Correct displacement at the interface
        const labelList& patchMap = subMeshes()[meshI].patchMap();
        forAll(patchMap, patchI)
        {
            if (patchMap[patchI] == -1)
            {
                label interfacePatchIndex = patchI;

                vectorField& interfacePatchD =
                    subMeshDisplacement[meshI]
                   .boundaryField()[interfacePatchIndex];

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

                    interfacePatchD[faceI] =
                        interfaceDisplacement[curInterFace];
                }
            }
        }

        // Calc point sub-mesh point displacement
        subMeshVolToPoint()[meshI].interpolate
            (
                subMeshDisplacement[meshI],
                subMeshPointDisplacement[meshI]
            );

        const vectorField& subMeshPointDI =
            subMeshPointDisplacement[meshI].internalField();

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
    if (mesh().globalData().nGlobalPoints())
    {
        for (label k=0; k<mesh().globalData().nGlobalPoints(); k++)
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


    label size = isolatedInterfacePoints().size();

    reduce(size, sumOp<label>());


    if (size)
    {
        // Correct isolated points
        forAll(pointDisplacement.boundaryField(), patchI)
        {
            if
            (
                pointDisplacement.boundaryField()[patchI].type()
             == processorPointPatchVectorField::typeName
            )
            {
                if (Pstream::parRun())
                {
                    const processorPointPatchVectorField& procPatchDispl =
                        refCast<processorPointPatchVectorField>
                        (
                            pointDisplacement.boundaryField()[patchI]
                        );

                    const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>
                        (
                            mesh().boundaryMesh()[patchI]
                        );

                    vectorField pif =
                        procPatchDispl.patchInternalField
                        (
                            pointDisplacement.internalField()
                        );

                    OPstream::write
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo(),
                        reinterpret_cast<const char*>(pif.begin()),
                        pif.byteSize()
                    );
                }
            }
        }

        forAll(pointDisplacement.boundaryField(), patchI)
        {
            if
            (
                pointDisplacement.boundaryField()[patchI].type()
             == processorPointPatchVectorField::typeName
            )
            {
                if (Pstream::parRun())
                {
                    const processorPointPatchVectorField& procPatchDispl =
                        refCast<processorPointPatchVectorField>
                        (
                            pointDisplacement.boundaryField()[patchI]
                        );

                    const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>
                        (
                            mesh().boundaryMesh()[patchI]
                        );

                    tmp<vectorField> tNgbProcPatchDispl
                    (
                        new vectorField(procPatchDispl.size(), vector::zero)
                    );
                    vectorField& ngbProcPatchDispl = tNgbProcPatchDispl();

                    IPstream::read
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo(),
                        reinterpret_cast<char*>(ngbProcPatchDispl.begin()),
                        ngbProcPatchDispl.byteSize()
                    );

                    forAll(isolatedInterfacePoints(), pI)
                    {
                        label curPoint =
                            isolatedInterfacePoints()[pI];

                        const labelList& curPointFaces =
                            mesh().pointFaces()[curPoint];

                        forAll(curPointFaces, fI)
                        {
                            label faceID = curPointFaces[fI];
                            label patchID =
                                mesh().boundaryMesh().whichPatch(faceID);

                            if (patchID == patchI)
                            {
                                label curPatchPoint =
                                    mesh().boundaryMesh()[patchI]
                                    .meshPointMap()[curPoint];

                                pointDI[curPoint] =
                                    ngbProcPatchDispl[curPatchPoint];

                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    pointDisplacement.correctBoundaryConditions();

    // Re-evaluate componentMixed patches
    forAll(pointDisplacement.boundaryField(), patchI)
    {
        if
        (
            pointDisplacement.boundaryField()[patchI].type()
         == componentMixedPointPatchVectorField::typeName
        )
        {
            pointDisplacement.boundaryField()[patchI].evaluate();
        }
    }
}


void materialInterface::updateDisplacementGradient
(
    const volVectorField& displacement,
    const PtrList<volVectorField>& subMeshDisplacement,
    const PtrList<pointVectorField>& subMeshPointDisplacement,
    volTensorField& cellDisplacementGradient,
    surfaceTensorField& faceDisplacementGradient
)
{
    tensorField& gradDI = cellDisplacementGradient.internalField();
    tensorField& gradDfI = faceDisplacementGradient.internalField();
    gradDfI = tensor::zero;

    // Calculate displacement gradient for sub-meshes
    forAll(subMeshes(), meshI)
    {
        // Cell gradient
        volTensorField subMeshGradD =
            fvc::cellLimitedGrad
            (
                subMeshDisplacement[meshI],
                subMeshPointDisplacement[meshI],
                0
            );

        // Face gradient in tangential direction
        surfaceTensorField subMeshGradDf =
            fvc::fsGrad
            (
                subMeshDisplacement[meshI],
                subMeshPointDisplacement[meshI]
            );

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
                    cellDisplacementGradient.boundaryField()[patchMap[patchI]];

                if (!patchGradD.coupled())
                {
                    forAll(subMeshPatchGradD, faceI)
                    {
                        label globalGlobalMeshFace = faceMap[start+faceI];

                        label curGlobalMeshPatchFace =
                            globalGlobalMeshFace
                          - mesh().boundaryMesh()[patchMap[patchI]].start();

                        patchGradD[curGlobalMeshPatchFace] =
                            subMeshPatchGradD[faceI];
                    }
                }

                fvsPatchField<tensor>& patchGradDf =
                    faceDisplacementGradient.boundaryField()[patchMap[patchI]];

                forAll(subMeshPatchGradDf, faceI)
                {
                    label globalGlobalMeshFace = faceMap[start+faceI];

                    label curGlobalMeshPatchFace =
                        globalGlobalMeshFace
                      - mesh().boundaryMesh()[patchMap[patchI]].start();

                    patchGradDf[curGlobalMeshPatchFace] =
                        subMeshPatchGradDf[faceI];
                }
            }
            else // interface faces
            {
                forAll(subMeshPatchGradDf, faceI)
                {
                    label globalGlobalMeshFace = faceMap[start+faceI];

                    if (globalGlobalMeshFace < mesh().nInternalFaces())
                    {
                        gradDfI[globalGlobalMeshFace] +=
                            0.5*subMeshPatchGradDf[faceI];
                    }
                    else
                    {
                        label curPatch =
                            mesh().boundaryMesh().whichPatch
                            (
                                globalGlobalMeshFace
                            );
                        label curPatchFace =
                            globalGlobalMeshFace
                          - mesh().boundaryMesh()[curPatch].start();

                        faceDisplacementGradient
                       .boundaryField()[curPatch][curPatchFace] =
                            subMeshPatchGradDf[faceI];
                    }
                }
            }
        }
    }

    // Correct cell gradient at the boundary
    fv::gaussGrad<vector>(mesh()).correctBoundaryConditions
        (
            displacement,
            cellDisplacementGradient
        );

    // Make sure face gradient is consistent accross processor patches
    forAll(faceDisplacementGradient.boundaryField(), patchI)
    {
        fvsPatchField<tensor>& patchGradDf =
            faceDisplacementGradient.boundaryField()[patchI];

        if (patchGradDf.type() == processorFvsPatchTensorField::typeName)
        {
            const tensorField& patchField = patchGradDf;

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
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
    forAll(faceDisplacementGradient.boundaryField(), patchI)
    {
        fvsPatchField<tensor>& patchGradDf =
            faceDisplacementGradient.boundaryField()[patchI];

        if (patchGradDf.type() == processorFvsPatchTensorField::typeName)
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
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
    faceDisplacementGradient +=
        mesh().Sf()*fvc::snGrad(displacement)/mesh().magSf();
}


void materialInterface::modifyProperty
(
    surfaceScalarField& muf,
    const volScalarField& mu
) const
{
    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh().nInternalFaces())
        {
            muf.internalField()[curFace] = 0;
        }
        else
        {
            label curPatch = mesh().boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh().boundaryMesh()[curPatch].start();

            muf.boundaryField()[curPatch][curPatchFace] = 0;
        }
    }
}


tmp<surfaceScalarField> materialInterface::interpolate
(
    const volScalarField& mu
) const
{
    tmp<surfaceScalarField> tMuf
    (
        new surfaceScalarField
        (
            IOobject
            (
                mu.name() + 'f',
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", mu.dimensions(), 0)
        )
    );
    surfaceScalarField& muf = tMuf();

    muf = fvc::interpolate(mu);

    this->modifyProperty(muf, mu);

    return tMuf;
}


tmp<volTensorField> materialInterface::grad
(
    const volVectorField& displacement,
    const vectorField& interfaceDisplacement
) const
{
    tmp<volTensorField> tGradD
    (
        new volTensorField
        (
            IOobject
            (
                "grad(" + displacement.name() + ')',
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("zero", dimless, tensor::zero)
        )
    );
    volTensorField& gradD = tGradD();

    surfaceVectorField Df = fvc::interpolate(displacement);


    // Skew-correction

    if (skewCorrectionVectors::New(mesh()).skew())
    {
        const skewCorrectionVectors& scv = skewCorrectionVectors::New(mesh());

        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>
            (
                "grad(" + displacement.name() + ')'
            );

        Df +=
        (
            scv()
          & linear<tensor>(mesh()).interpolate
            (
                gradD
            )
        );
    }


    // Interface correction

    const vectorField& interD = interfaceDisplacement;

    forAll(faces(), faceI)
    {
        label curFace = faces()[faceI];

        if (curFace < mesh().nInternalFaces())
        {
            Df.internalField()[curFace] = interD[faceI];
        }
        else
        {
            label curPatch = mesh().boundaryMesh().whichPatch(curFace);
            label curPatchFace =
                curFace - mesh().boundaryMesh()[curPatch].start();

            Df.boundaryField()[curPatch][curPatchFace] = interD[faceI];
        }
    }

    // Gradient calculation using Gauss method

    gradD = fv::gaussGrad<vector>(mesh()).gradf(Df, "grad(D)");
    fv::gaussGrad<vector>(mesh()).correctBoundaryConditions
    (
        displacement,
        gradD
    );

    return tGradD;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

Description

\*---------------------------------------------------------------------------*/

#include "meshReader.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "emptyPolyPatch.H"
#include "cellModeller.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Pointers to cell shape models
const cellModel* meshReader::unknownModel = cellModeller::lookup("unknown");
const cellModel* meshReader::tetModel = cellModeller::lookup("tet");
const cellModel* meshReader::pyrModel = cellModeller::lookup("pyr");
const cellModel* meshReader::prismModel = cellModeller::lookup("prism");
const cellModel* meshReader::hexModel = cellModeller::lookup("hex");

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// write cellTableId as volume field
// this is crucial for later conversion back to ccm
void meshReader::volFieldCellTableId(fvMesh & mesh)
{
    const word & propertyName = "cellTableId";

    volScalarField calculatedField
    (
        IOobject
        (
            propertyName,
            "0",
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(propertyName, dimless, 0.0)
    );

    forAll(cellTableId_, cellI)
    {
        calculatedField[cellI] = cellTableId_[cellI];
    }

    Info<< "Writing " << propertyName << " to "
        << calculatedField.objectPath() << endl;

    // NOTE:
    // the cellTableId is an integer and almost always < 1000, thus ASCII
    // will be compacter than binary and makes external scripting easier
    // Fixed by HJ, 5/Jan/2007
    calculatedField.write();
}

void meshReader::writeMesh()
{
    // Make polyhedral mesh data (packing)
    Info << "Creating a polyMesh" << endl;

    createPolyCells();

    Info<< "Number of internal faces: " << nInternalFaces_ << endl;

    createPolyBoundary();

    clearExtraStorage();

    fvMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime_.constant(),
            runTime_
        ),
        points(),
        meshFaces_,
        cellPolys_
    );

    // adding patches also checks the mesh
    mesh.addFvPatches(polyBoundaryPatches(mesh));
    Info << "Writing polyMesh" << endl;
    mesh.write();

    volFieldCellTableId(mesh);    
    writeAux(mesh);
}

// Clear extra storage before creation of the mesh to reduce the memory usage
void meshReader::clearExtraStorage()
{
    Info << "Clearing extra storage" << endl;

    origCellId_.clear();
    cellFaces_.clear();
    baffleFaces_.clear();
    boundaryCells_.clear();
    boundaryFaces_.clear();
    baffleCellIds_.clear();
    baffleFaceIds_.clear();

    deleteDemandDrivenData(pointCellsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// Construct from components
meshReader::meshReader
(
    const fileName& file_or_prefix,
    const Time& runtime,
    const scalar scaleFactor
)
:
    runTime_(runtime),
    nInternalFaces_(0),
    patchStarts_(0),
    patchSizes_(0),
    interfaces_(0),
    baffleCellIds_(0),
    baffleFaceIds_(0),
    meshFaces_(0),
    cellPolys_(0),
    origCellId_(0),
    boundaryCells_(0),
    boundaryFaces_(0),
    patchTypes_(0),
    patchNames_(0),
    patchPhysicalTypes_(0),
    cellFaces_(0),
    baffleFaces_(0),
    cellTableId_(0),
    pointCellsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshReader::~meshReader()
{
    deleteDemandDrivenData(pointCellsPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

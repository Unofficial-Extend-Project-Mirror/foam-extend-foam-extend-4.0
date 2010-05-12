/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006 Mark Olesen
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
    Writes out the FOAM mesh in pro-STAR (v4) bnd/cel/vrt format.

    Alternatively, extracts the surface of the FOAM mesh into
    pro-STAR (v4) .cel/.vrt/ format.
    This can be useful, for example, for surface morphing in an external
    package.

    The cellTableId and cellTable information are used (if available).
    Otherwise the cellZones are used (if available).

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "volFields.H"
#include "cellModeller.H"
#include "SortableList.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Cell shape models
static const cellModel* unknownPtr_ = cellModeller::lookup("unknown");
static const cellModel* tetPtr_ = cellModeller::lookup("tet");
static const cellModel* pyrPtr_ = cellModeller::lookup("pyr");
static const cellModel* prismPtr_ = cellModeller::lookup("prism");
static const cellModel* hexPtr_ = cellModeller::lookup("hex");

// face addressing from foam faces -> pro-STAR faces for primitive shapes
static const label foamToStarFaceAddressing[4][6] =
{
    { 4, 5, 2, 3, 0, 1 },     // 11 = hex
    { 0, 1, 4, 5, 2, -1 },    // 12 = prism
    { 5, 4, 2, 0, -1, -1 },   // 13 = tetra
    { 0, 4, 3, 5, 2, -1 }     // 14 = pyramid
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// use file globals until we make a class
static labelList cellTableId_;

// a very stripped-down cell table
// - map the cell type id -> material type id (1 = fluid, 2 = solid)
static Map<label> cellTableMap_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Prostar 4+ header format
void prostarHeader(Ostream &os, const char * filetype)
{
    os  << "PROSTAR_" << filetype << endl
    << 4000 << " "
    << 0 << " "
    << 0 << " "
    << 0 << " "
    << 0 << " "
    << 0 << " "
    << 0 << " "
    << 0 << endl;
}

void getCellTable(const fvMesh & mesh)
{
    cellTableMap_.clear();
    cellTableId_.setSize(mesh.nCells(), 1);

    IOdictionary cellTableDict
    (
    IOobject
    (
        "cellTable",
        "constant",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    )
    );

    volScalarField volField
    (
    IOobject
    (
        "cellTableId",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedScalar("cellTableId", dimless, 1.0)
    );

    // get cellTableId information from the volScalarField if possible
    if (volField.headerOk())
    {
    const scalarField & field = volField.internalField();

    forAll(field, cellI)
    {
        cellTableId_[cellI] = static_cast<int>(field[cellI]);
    }

    if (cellTableDict.headerOk())
    {
        // convert dictionary to map
        wordList toc = cellTableDict.toc();

        forAll(toc, i)
        {
        word keyword = toc[i];
        if (!cellTableDict.isDict(keyword)) continue;

        const dictionary & dict = cellTableDict.subDict(keyword);

        if (dict.found("Id") && dict.found("MaterialType"))
        {
            label Id;
            dict["Id"] >> Id;
            dict["MaterialType"] >> keyword;

            if (keyword == "fluid")
            {
            cellTableMap_.insert(Id, 1);
            }
            else if (keyword == "solid")
            {
            cellTableMap_.insert(Id, 2);
            }
        }
        }
    }
    else
    {
        Info<< "No cellTable information available" << endl;
    }
    }
    else
    {
    Info<< "No cellTableId information available - using cellZones (if available)" << endl;

    const Map<label> & zoneMap = mesh.cellZones().zoneMap();

    // start zoned cells at cell type 1
    label typeOffset = 1;

    // fewer zoned cells than total cells
    // - leave unzoned cells as type 1 and start zoned cells at cell type 2
    if (zoneMap.size() < mesh.nCells())
    {
        typeOffset = 2;
    }

    forAllConstIter(Map<label>, zoneMap, iter)
    {
        cellTableId_[iter.key()] = iter() + typeOffset;
    }
    }
}

void writePoints
(
    const polyMesh& mesh,
    const fileName& timeName,
    const scalar scaleFactor
)
{
    fileName name(mesh.time().path()/"meshExport_" + timeName + ".vrt");

    OFstream outputFile(name);
    prostarHeader(outputFile, "VERTEX");

    // Set the precision of the points data to 10
    outputFile.precision(10);

    // force decimal point for Fortran input
    outputFile.setf(std::ios::showpoint);

    const pointField& points = mesh.points();

    Info<< "Writing " << name << " : "
    << points.size() << " points" << endl;

    forAll(points, ptI)
    {
    // convert [m] -> [mm]
    outputFile
        << ptI + 1 << " "
        << scaleFactor * points[ptI].x() << " "
        << scaleFactor * points[ptI].y() << " "
        << scaleFactor * points[ptI].z() << endl;
    }
}

void writeCells(const polyMesh& mesh, const fileName& timeName)
{
    fileName name(mesh.time().path()/"meshExport_" + timeName + ".cel");

    OFstream outputFile(name);
    prostarHeader(outputFile, "CELL");

    // this is what we seem to need
    // map foam cellModeller index -> star shape
    Map<label> shapeLookupIndex;
    shapeLookupIndex.insert(hexPtr_->index(), 11);
    shapeLookupIndex.insert(prismPtr_->index(), 12);
    shapeLookupIndex.insert(tetPtr_->index(), 13);
    shapeLookupIndex.insert(pyrPtr_->index(), 14);

    const cellShapeList& shapes = mesh.cellShapes();
    const cellList & cells  = mesh.cells();
    const faceList & faces  = mesh.faces();
    const labelList & owner = mesh.faceOwner();

    Info<< "Writing " << name << " : "
    << cells.size() << " cells" << endl;

    forAll(cells, cellId)
    {
    label tableId = cellTableId_[cellId];
    label materialType  = 1;    // 1(fluid)
    if (cellTableMap_.found(tableId))
    {
        materialType = cellTableMap_[tableId];
    }

    const cellShape & shape = shapes[cellId];
    label mapIndex = shape.model().index();

    // a registered primitive type
    if (shapeLookupIndex.found(mapIndex))
    {
        label shapeId = shapeLookupIndex[mapIndex];
        const labelList & vrtList = shapes[cellId];

        outputFile
        << cellId + 1 << " "
        << shapeId << " "
        << vrtList.size() << " "
        << tableId << " "
        << materialType;

        // primitives have <= 8 vertices, but prevent overrun anyhow
        label count = 0;
        forAll(vrtList, i)
        {
        if ((count % 8) == 0)
        {
            outputFile << endl;
            outputFile << cellId + 1;
        }
        outputFile << " " << vrtList[i] + 1;
        count++;
        }
        outputFile << endl;

    }
    else
    {
        label shapeId = 255;    // treat as general polyhedral
        const labelList & cFaces  = cells[cellId];

        // create (beg,end) indices
        List<label> indices(cFaces.size() + 1);
        indices[0] = indices.size();

        label count = indices.size();
        // determine the total number of vertices
        forAll(cFaces, faceI)
        {
        count += faces[cFaces[faceI]].size();
        indices[faceI+1] = count;
        }

        outputFile
        << cellId + 1 << " "
        << shapeId << " "
        << count << " "
        << tableId << " "
        << materialType;

        // write indices - max 8 per line
        count = 0;
        forAll(indices, i)
        {
        if ((count % 8) == 0)
        {
            outputFile << endl;
            outputFile << cellId + 1;
        }
        outputFile << " " << indices[i];
        count++;
        }

        // write faces - max 8 per line
        forAll(cFaces, faceI)
        {
        label meshFace = cFaces[faceI];
        face f;

        if (owner[meshFace] == cellId)
        {
            f = faces[meshFace];
        }
        else
        {
            f = faces[meshFace].reverseFace();
        }

        forAll(f, i)
        {
            if ((count % 8) == 0)
            {
            outputFile << endl;
            outputFile << cellId + 1;
            }

            outputFile << " " << f[i] + 1;
            count++;
        }
        }

        outputFile << endl;
    }
    }
}

void writeBoundary(const polyMesh& mesh, const fileName& timeName)
{
    fileName name(mesh.time().path()/"meshExport_" + timeName + ".bnd");

    OFstream outputFile(name);
    prostarHeader(outputFile, "BOUNDARY");

    const cellShapeList& shapes = mesh.cellShapes();
    const cellList & cells  = mesh.cells();
    const faceList & faces  = mesh.faces();
    const labelList & owner = mesh.faceOwner();
    const polyBoundaryMesh & patches = mesh.boundaryMesh();

    // this is what we seem to need
    // these MUST correspond to foamToStarFaceAddressing
    //
    Map<label> faceLookupIndex;
    faceLookupIndex.insert(hexPtr_->index(), 0);
    faceLookupIndex.insert(prismPtr_->index(), 1);
    faceLookupIndex.insert(tetPtr_->index(), 2);
    faceLookupIndex.insert(pyrPtr_->index(), 3);

    Info<< "Writing " << name << " : "
    << (mesh.nFaces() - patches[0].start()) << " boundaries" << endl;

    label boundId = 0;
    // Write boundary faces
    //
    forAll(patches, patchI)
    {
    label patchStart = patches[patchI].start();
    label patchSize  = patches[patchI].size();

    for
    (
        label faceI = patchStart;
        faceI < (patchStart + patchSize);
        ++faceI
    )
    {
        label cellId = owner[faceI];
        const labelList & cFaces  = cells[cellId];
        const cellShape & shape = shapes[cellId];
        label cellFaceId = findIndex(cFaces, faceI);

        //         Info<< "cell " << cellId + 1 << " face " << faceI
        //         << " == " << faces[faceI]
        //         << " is index " << cellFaceId << " from " << cFaces;

        // Unfortunately, the order of faces returned by
        //   primitiveMesh::cells() is not necessarily the same
        //   as defined by primitiveMesh::cellShapes()
        // Thus, for registered primitive types, do the lookup ourselves.
        // Finally, the cellModel face number is re-mapped to the
        // Star-CD local face number

        label mapIndex = shape.model().index();

        // a registered primitive type
        if (faceLookupIndex.found(mapIndex))
        {
        const faceList sFaces = shape.faces();
        forAll(sFaces, sFaceI)
        {
            if (faces[faceI] == sFaces[sFaceI])
            {
            cellFaceId = sFaceI;
            break;
            }
        }

        mapIndex = faceLookupIndex[mapIndex];
        cellFaceId = foamToStarFaceAddressing[mapIndex][cellFaceId];
        }
        // Info<< endl;

        boundId++;

        outputFile
        << boundId << " "
        << cellId + 1 << " "
        << cellFaceId + 1 << " "
        << patchI + 1 << " "
        << 0 << " "
        << "PATCH" << endl;

    }
    }
}

void writeVolumeMesh
(
    const fvMesh& mesh,
    const fileName& timeName,
    const scalar scaleFactor
)
{
    getCellTable(mesh);
    writePoints(mesh, timeName, scaleFactor);
    writeCells(mesh, timeName);
    writeBoundary(mesh, timeName);
}

void writeSurfaceMesh
(
    const fvMesh& mesh,
    const fileName& timeName,
    const scalar scaleFactor
)
{
    getCellTable(mesh);

    word prefix("surfaceExport_" + timeName);
    fileName name(mesh.time().path()/prefix + ".cel");

    Info << "Writing " << name << endl;

    OFstream celFile(name);
    prostarHeader(celFile, "CELL");

    // mesh and patch info
    const pointField & points = mesh.points();
    const labelList & owner = mesh.faceOwner();
    const faceList  & meshFaces = mesh.faces();
    const polyBoundaryMesh & patches = mesh.boundaryMesh();

    label shapeId = 3;  // shell/baffle element
    label typeId  = 4;    // 4(shell)

    // remember which points need to be written
    labelHashSet pointHash;

    // Write boundary faces as normal Star-CD mesh
    // use the face Id as the cell Id,
    // use the cell table id of the face owner - allows separation of parts
    forAll(patches, patchI)
    {
    label patchStart = patches[patchI].start();
    label patchSize  = patches[patchI].size();

    // use face id as cell id
    for
    (
        label faceI = patchStart;
        faceI < (patchStart + patchSize);
        ++faceI
    )
    {
        const labelList & vrtList = meshFaces[faceI];
        label cellId = faceI;

        celFile
        << cellId + 1 << " "
        << shapeId << " "
        << vrtList.size() << " "
        << cellTableId_[owner[faceI]] << " "
        << typeId;

        // likely <= 8 vertices, but prevent overrun anyhow
        label count = 0;
        forAll(vrtList, i)
        {
        if ((count % 8) == 0)
        {
            celFile << endl;
            celFile << cellId + 1;
        }
        // remember which points we'll need to write
        pointHash.insert(vrtList[i]);
        celFile << " " << vrtList[i] + 1;
        count++;
        }
        celFile << endl;
    }
    }

    name = (mesh.time().path()/prefix + ".vrt");

    Info << "Writing " << name << endl;
    OFstream vrtFile(name);
    prostarHeader(vrtFile, "VERTEX");

    vrtFile.precision(10);
    vrtFile.setf(std::ios::showpoint);    // force decimal point for Fortran

    // build sorted table of contents
    SortableList<label> toc(pointHash.size());
    {
    label i = 0;
    forAllConstIter(labelHashSet, pointHash, iter)
    {
        toc[i++] = iter.key();
    }
    }
    toc.sort();
    pointHash.clear();

    // write points in sorted order
    forAll(toc, i)
    {
    label vrtId = toc[i];
    vrtFile
        << vrtId + 1 << " "
        << scaleFactor * points[vrtId].x() << " "
        << scaleFactor * points[vrtId].y() << " "
        << scaleFactor * points[vrtId].z() << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("noscale", "");
    argList::validOptions.insert("surface", "");

#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

    // rescale from [m] to [mm] by default
    scalar scaleFactor = 1000.0;

    if (args.options().found("noscale"))
    {
        scaleFactor = 1.0;
    }

    bool surfaceOnly = false;

    if (args.options().found("surface"))
    {
        surfaceOnly = true;
    }

#   include "createMesh.H"

    bool firstCheck = true;

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        polyMesh::readUpdateState state = mesh.readUpdate();

        if (firstCheck || state != polyMesh::UNCHANGED)
        {
            if (surfaceOnly)
            {
                writeSurfaceMesh(mesh, runTime.timeName(), scaleFactor);
            }
            else
            {
                writeVolumeMesh(mesh, runTime.timeName(), scaleFactor);
            }
        }
        else
        {
            Info << "No mesh." << endl;
        }

        firstCheck = false;

        Info << endl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //

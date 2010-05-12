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
    Read pro-STAR vrt/cel/bnd files.

 * Starting with Pro-STAR version 4, the files have become easier to read.
 * Vertices are space-delimited.
 * The cell format is logical.
 * Trimmed and degenerate cells are saved as polyhedral.
 * The boundaries corresponds to cells and their faces.
 *
 *--------------------------------------------------------------------------
 * changes:
 *
 * * re-write for polyhedral cells and new parsing formats
 *   copyright (C) 2006 Mark Olesen
\*---------------------------------------------------------------------------*/

#include "starMeshReader.H"
#include "cyclicPolyPatch.H"
#include "emptyPolyPatch.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "cellModeller.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// face addressing from pro-STAR faces -> foam faces
static const label starToFoamFaceAddressing[4][6] =
{
    { 4, 5, 2, 3, 0, 1 },     // 11 = hex
    { 0, 1, 4, -1, 2, 3 },    // 12 = prism
    { 3, -1, 2, -1, 1, 0 },   // 13 = tetra
    { 0, -1, 4, 2, 1, 3 }     // 14 = pyramid
};

// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

// useful enough to wander to meshes/meshShapes/face ?

// make face canonical - lowest vertex number first
static void face_canonical(face& thisFace)
{
    label pivot = findMin(thisFace);

    if (pivot > 0)
    {
        label count = thisFace.size();
        labelList oldToNew(count);

        forAll(oldToNew, i)
        {
            oldToNew[i] = i - pivot;

            // wrap
            if (oldToNew[i] < 0)
            {
                oldToNew[i] += count;
            }
        }

        inplaceReorder(oldToNew, thisFace);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ignore the rest of the line
void readToNewline(ISstream& is)
{
    char ch = '\n';
    do
    {
        (is).get(ch);
    }
    while ((is) && ch != '\n');
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static bool readProstarHeader(IFstream &is)
{
    word header;
    label majorVersion;

    is >> header;
    is >> majorVersion;

    // add checks ...

    // skip the rest of the line
    readToNewline(is);

    return true;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// read in the points from the .vrt file
//
/*---------------------------------------------------------------------------*\
Line 1:
  PROSTAR_VERTEX [newline]

Line 2:
  <version> 0 0 0 0 0 0 0 [newline]

Body:
  <vertexId>  <x>  <y>  <z> [newline]

\*---------------------------------------------------------------------------*/
void starMeshReader::readPoints
(
    const fileName & inputFileName,
    const scalar scaleFactor
)
{
    label nPoints = 0;
    label maxId = 0;

    // first pass
    // count nPoints and maxId
    {
        IFstream inputFile(inputFileName);

        // Pass 1: get # points and maximum vertex label

        if (inputFile.good())
        {
            label lineLabel;
            scalar x, y, z;

            readProstarHeader(inputFile);

            while ((inputFile >> lineLabel).good())
            {
                nPoints++;
                maxId = max(maxId, lineLabel);
                inputFile >> x >> y >> z;
            }
        }
        else
        {
            FatalErrorIn("starMeshReader::readPoints()")
                << "cannot read file " << inputFileName
                << abort(FatalError);
        }
    }

    Info<< "Number of points  = " << nPoints << endl;

    // Original Point number for a given vertex
    labelList origPointId(nPoints);

    points_.setSize(nPoints);
    mapToFoamPointId_.setSize(maxId+1);

    // reset to invalid values
    origPointId = -1;
    mapToFoamPointId_ = -1;

    // Pass 2: construct pointlist and conversion table
    // from Star vertex numbers to Foam point labels
    if (nPoints > 0)
    {
        IFstream inputFile(inputFileName);
        label lineLabel;

        readProstarHeader(inputFile);

        for (label pointI = 0; (inputFile >> lineLabel).good(); ++pointI)
        {
            inputFile
                >> points_[pointI].x()
                >> points_[pointI].y()
                >> points_[pointI].z();

            origPointId[pointI] = lineLabel;
            mapToFoamPointId_[lineLabel] = pointI;
        }

        if (scaleFactor > 1.0 + SMALL || scaleFactor < 1.0 - SMALL)
        {
            points_ *= scaleFactor;
        }
    }
    else
    {
        FatalErrorIn("starMeshReader::readPoints()")
            << "no points in file " << inputFileName
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// read in the cells from the .cel file
//
/*---------------------------------------------------------------------------*\
Line 1:
  PROSTAR_CELL [newline]

Line 2:
  <version> 0 0 0 0 0 0 0 [newline]

Body:
  <cellId>  <shapeId>  <nVertices>  <cellTableId>  <typeId> [newline]
  <cellId>  <int1> .. <int8>
  <cellId>  <int9> .. <int16>

 with shapeId:
 *   1 = point
 *   2 = line
 *   3 = shell
 *  11 = hexa
 *  12 = prism
 *  13 = tetra
 *  14 = pyramid
 * 255 = polyhedron

 with typeId
 *   1 = fluid
 *   2 = solid
 *   3 = baffle
 *   4 = shell
 *   5 = line
 *   6 = point

For primitive cell shapes, the number of vertices will never exceed 8 (hexa).
For polyhedral, the number of vertices includes the index table comprising
beg/end pairs for each cell face.

Strictly speaking, we only need the cellModeller for adding boundaries
\*---------------------------------------------------------------------------*/

void starMeshReader::readCells
(
    const fileName & inputFileName
)
{
    label nCells = 0;
    label nBaffles = 0;
    label maxId = 0;

    // save cellTable Id MaterialType information
    Map<label> cellTableMap;

    // first pass
    // count nCells, nBaffles and maxId
    // also see if polyhedral cells were used
    {
        IFstream inputFile(inputFileName);

        if (inputFile.good())
        {
            label lastLabel = -1;
            label lineLabel, shapeId, nVertices, cellTableId, typeId;

            readProstarHeader(inputFile);

            while ((inputFile >> lineLabel).good())
            {
                if (lineLabel != lastLabel)
                {
                    lastLabel = lineLabel;

                    inputFile
                        >> shapeId >> nVertices >> cellTableId >> typeId;

                    // fluid/solid cells
                    if (typeId == 1 || typeId == 2)
                    {
                        if (!cellTableMap.found(cellTableId))
                        {
                            cellTableMap.insert(cellTableId, typeId);
                        }

                        maxId = max(maxId, lineLabel);
                        nCells++;
                    }
                    else if (typeId == 3)
                    {
                        // baffles
                        maxId = max(maxId, lineLabel);
                        nBaffles++;
                    }
                }

                // skip the rest of the line
                readToNewline(inputFile);
            }
        }
        else
        {
            FatalErrorIn("starMeshReader::readCells()")
                << "cannot read file " << inputFileName
                << abort(FatalError);
        }
    }

    Info<< "Number of cells   = " << nCells  << endl
        << "Number of baffles = " << nBaffles << endl << endl;

    cellTable_.clear();

    forAllIter(Map<label>, cellTableMap, iter)
    {
        label id = iter.key();
        label materialType = cellTableMap[id];

        dictionary dict;
        dict.add("Id", id);
        if (materialType == 1)
        {
            dict.add("MaterialType", "fluid");
            cellTable_.insert(id, dict);
        }
        else if (materialType == 2)
        {
            dict.add("MaterialType", "solid");
            cellTable_.insert(id, dict);
        }
    }

    cellFaces_.setSize(nCells);
    cellShapes_.setSize(nCells);
    cellTableId_.setSize(nCells);

    // information for the interfaces
    baffleFaces_.setSize(nBaffles);

    // extra space for baffles
    origCellId_.setSize(nCells+nBaffles);
    mapToFoamCellId_.setSize(maxId+1);
    mapToFoamCellId_ = -1;

    // Pass 2: construct cellFaces_ and possibly cellShapes_
    if (nCells > 0)
    {
        IFstream inputFile(inputFileName);

        readProstarHeader(inputFile);

        labelList starLabels(64);
        label lineLabel, shapeId, nVertices, cellTableId, typeId;

        label cellI = 0;
        label baffleI = 0;

        while ((inputFile >> lineLabel).good())
        {
            label starCellId = lineLabel;
            inputFile
                >> shapeId >> nVertices >> cellTableId >> typeId;

            label nLabels = nVertices;
            if (nLabels > starLabels.size())
            {
                starLabels.setSize(nLabels);
            }
            starLabels = -1;

            // read indices - max 8 per line
            for (label i = 0; i < nLabels; ++i)
            {
                if ((i % 8) == 0)
                {
                    inputFile >> lineLabel;
                }
                inputFile >> starLabels[i];
            }

            // determine the foam cell shape
            const cellModel* curModelPtr = NULL;

            // fluid/solid cells
            switch (shapeId)
            {
                case 11:
                    curModelPtr = hexModel;
                    break;
                case 12:
                    curModelPtr = prismModel;
                    break;
                case 13:
                    curModelPtr = tetModel;
                    break;
                case 14:
                    curModelPtr = pyrModel;
                    break;
            }

            if (curModelPtr)
            {
                // primitive cell - use shapes

                // Record original Star cell number and lookup
                origCellId_[cellI] = starCellId;
                mapToFoamCellId_[starCellId] = cellI;

                // convert Star vertex number to point label
                for (label i=0; i < nVertices; ++i)
                {
                    label pointId = mapToFoamPointId_[starLabels[i]];

                    if (pointId < 0)
                    {
                        Info<< "Cells inconsistent with vertex file. "
                            << "Star vertex " << starLabels[i]
                            << " does not exist" << endl;
                    }

                    starLabels[i] = pointId;
                }

                cellTableId_[cellI] = cellTableId;
                cellShapes_[cellI] = cellShape
                (
                    *curModelPtr,
                    SubList<label>(starLabels, nVertices)
                );

                cellFaces_[cellI] = cellShapes_[cellI].faces();
                cellI++;
            }
            else if (shapeId == 255)
            {
                // poly cell - create cell faces directly
                label nFaces = starLabels[0] - 1;

                // Record original Star cell number and lookup
                origCellId_[cellI] = starCellId;
                mapToFoamCellId_[lineLabel] = cellI;

                // convert Star vertex number to point label
                for (label i=starLabels[0]; i < nVertices; ++i)
                {
                    label pointId = mapToFoamPointId_[starLabels[i]];

                    if (pointId < 0)
                    {
                        Info<< "Cells inconsistent with vertex file. "
                            << "Star vertex " << starLabels[i]
                            << " does not exist" << endl;
                    }

                    starLabels[i] = pointId;
                }

                // avoid empty faces and/or strange indices
                label skipFaces = 0;
                for (label i = 0; i < nFaces; ++i)
                {
                    if (starLabels[i + 1] <= starLabels[i])
                    {
                        ++skipFaces;
                    }
                }

                // ignore these faces
                nFaces -= skipFaces;

                if (skipFaces)
                {
                    Info<< "star cell " << starCellId << " has " << skipFaces
                        << " empty faces - could cause boundary "
                        << "addressing problems"
                        << endl;
                }

                cellFaces_[cellI] = faceList(nFaces);
                faceList & faces = cellFaces_[cellI];

                for (label i=0, faceI=0; i < nFaces+skipFaces; ++i)
                {
                    // valid faces only
                    if (starLabels[i+1] > starLabels[i])
                    {
                        faces[faceI] = face
                        (
                            SubList<label>
                            (
                                starLabels,
                                starLabels[i+1] - starLabels[i],
                                starLabels[i]
                            )
                        );

                        // this would ideally not be required
                        // - ie., should be added to the face() constructor
                        face_canonical(faces[faceI]);
                        ++faceI;
                    }
                }

                cellTableId_[cellI] = cellTableId;
                cellShapes_[cellI] = cellShape(*unknownModel, labelList(0));
                cellI++;
            }
            else if (typeId == 3)
            {
                // baffle cells

                // convert Star vertex number to point label
                for (label i=0; i < nVertices; ++i)
                {
                    label pointId = mapToFoamPointId_[starLabels[i]];

                    if (pointId < 0)
                    {
                        Info<< "Cells inconsistent with vertex file. "
                            << "Star vertex " << starLabels[i]
                            << " does not exist" << endl;
                    }

                    starLabels[i] = pointId;
                }

                // create face indices
                baffleFaces_[baffleI] = face
                (
                    SubList<label>(starLabels, nVertices)
                );

                // this is required
                face_canonical(baffleFaces_[baffleI]);

                // insert lookup addressing in normal list
                mapToFoamCellId_[starCellId] = nCells + baffleI;
                origCellId_[nCells + baffleI] = starCellId;
                baffleI++;
            }
        }
    }
    else
    {
        FatalErrorIn("starMeshReader::readCells()")
            << "no cells in file " << inputFileName
            << abort(FatalError);
    }

#ifdef DEBUG_READING
    Info << "READ CELLS" << endl;
#endif

    // cleanup
    mapToFoamPointId_.clear();
}


// read in the boundaries from the .bnd file
//
/*---------------------------------------------------------------------------*\
Line 1:
  PROSTAR_BOUNDARY [newline]

Line 2:
  <version> 0 0 0 0 0 0 0 [newline]

Body:
  <boundId>  <cellId>  <cellFace>  <regionId>  0  <boundaryType> [newline]

where boundaryType is truncated to 4 characters from one of the following:
INLET
PRESSSURE
OUTLET
BAFFLE
etc,
\*---------------------------------------------------------------------------*/

void starMeshReader::readBoundary
(
    const fileName & inputFileName
)
{
    label nPatches = 0, nFaces = 0, nBafflePatches = 0;
    label maxId = 0;

    label lineLabel, starCellId, cellFaceId, starRegion, configNumber;
    word patchType;

    labelList mapToFoamPatchId(1000, -1);
    labelList nPatchFaces(1000, 0);
    labelList origRegion(1000, 0);
    patchTypes_.setSize(1000);

    // this is what we seem to need
    // these MUST correspond to starToFoamFaceAddressing
    //
    Map<label> faceLookupIndex;

    faceLookupIndex.insert(hexModel->index(), 0);
    faceLookupIndex.insert(prismModel->index(), 1);
    faceLookupIndex.insert(tetModel->index(), 2);
    faceLookupIndex.insert(pyrModel->index(), 3);

    // Pass 1:
    // collect
    // no. of faces (nFaces), no. of patches (nPatches)
    // and for each of these patches the number of faces
    // (nPatchFaces[patchLabel])
    //
    // and a conversion table from Star regions to (Foam) patchLabels
    //
    // additionally note the no. of baffle patches (nBafflePatches)
    // so that we sort these to the end of the patch list
    // - this makes it easier to transfer them to an adjacent patch if reqd
    {
        IFstream inputFile(inputFileName);

        if (inputFile.good())
        {
            readProstarHeader(inputFile);

            while ((inputFile >> lineLabel).good())
            {
                nFaces++;
                inputFile
                    >> starCellId
                    >> cellFaceId
                    >> starRegion
                    >> configNumber
                    >> patchType;

                // Build translation table to convert star patch to foam patch
                label patchLabel = mapToFoamPatchId[starRegion];
                if (patchLabel == -1)
                {
                    patchLabel = nPatches;
                    mapToFoamPatchId[starRegion] = patchLabel;
                    origRegion[patchLabel] = starRegion;
                    patchTypes_[patchLabel] = patchType;

                    maxId = max(maxId, starRegion);

                    if (patchType == "BAFF")
                    {
                        nBafflePatches++;
                    }
                    nPatches++;
                }

                nPatchFaces[patchLabel]++;
            }

            if (nPatches == 0)
            {
                Info<< "No boundary faces in file " << inputFileName << endl;
            }
        }
        else
        {
            Info<< "Could not read boundary file " << inputFileName << endl;
        }
    }

    // keep empty patch region in reserve
    nPatches++;
    Info<< "Number of patches = " << nPatches
        << " (including extra for missing)" << endl;

    // resize
    origRegion.setSize(nPatches);
    patchTypes_.setSize(nPatches);
    patchNames_.setSize(nPatches);
    nPatchFaces.setSize(nPatches);

    // add our empty patch
    origRegion[nPatches - 1] = 0;
    nPatchFaces[nPatches - 1] = 0;
    patchTypes_[nPatches - 1] = "NONE";

    // create names
    forAll(patchTypes_, patchI)
    {
        patchNames_[patchI] = patchTypes_[patchI] + name(origRegion[patchI]);
    }
    patchNames_[nPatches - 1] = "Default_Boundary_Region";

    // re-sort to have baffles near the end
    if (nBafflePatches)
    {
        labelList oldToNew = identity(nPatches);
        label newIndex = 0;
        label baffleIndex = (nPatches-1 - nBafflePatches);

        for(label i = 0; i < oldToNew.size()-1; ++i)
        {
            if (patchTypes_[i] == "BAFF")
            {
                oldToNew[i] = baffleIndex++;
            }
            else
            {
                oldToNew[i] = newIndex++;
            }
        }

        inplaceReorder(oldToNew, origRegion);
        inplaceReorder(oldToNew, patchTypes_);
        inplaceReorder(oldToNew, patchNames_);
        inplaceReorder(oldToNew, nPatchFaces);
    }

    mapToFoamPatchId.setSize(maxId+1, -1);

    forAll(origRegion, patchI)
    {
        mapToFoamPatchId[origRegion[patchI]] = patchI;
    }

    boundaryCells_.setSize(nPatches);
    boundaryFaces_.setSize(nPatches);

    // size the lists and reset the counters to be used again
    forAll(boundaryCells_, patchI)
    {
        boundaryCells_[patchI].setSize(nPatchFaces[patchI], -1);
        boundaryFaces_[patchI].setSize(nPatchFaces[patchI], -1);
        nPatchFaces[patchI] = 0;
    }

    if (nPatches > 1)
    {
        IFstream inputFile(inputFileName);
        readProstarHeader(inputFile);

        while ((inputFile >> lineLabel).good())
        {
            inputFile
                >> starCellId
                >> cellFaceId
                >> starRegion
                >> configNumber
                >> patchType;

            label patchI = mapToFoamPatchId[starRegion];

            // zero-based indexing
            cellFaceId--;

            label cellId = -1;

            // convert to foam cell number
            if (starCellId < mapToFoamCellId_.size())
            {
                cellId = mapToFoamCellId_[starCellId];
            }

            if (cellId < 0)
            {
                Info<< "Boundaries inconsistent with cell file. "
                    << "Star cell " << starCellId
                    << " does not exist" << endl;
            }
            else
            {
                // restrict lookup to volume cells (no baffles)
                if (cellId < cellShapes_.size())
                {
                    label index = cellShapes_[cellId].model().index();
                    
                    if (faceLookupIndex.found(index))
                    {
                        index = faceLookupIndex[index];
                        cellFaceId =
                            starToFoamFaceAddressing[index][cellFaceId];
                    }
                }
                else
                {
                    // we currently use cellId >= nCells to tag baffles,
                    // we can also use a negative face number
                    cellFaceId = -1;
                }

                boundaryCells_[patchI][nPatchFaces[patchI]] = cellId;
                boundaryFaces_[patchI][nPatchFaces[patchI]] = cellFaceId;

#ifdef DEBUG_BOUNDARY
                Info << "bnd " << cellId << " " << cellFaceId << endl;
#endif
                // increment counter of faces in current patch
                nPatchFaces[patchI]++;
            }
        }
    }

    // retain original information ih patchPhysicalTypes_ - overwrite latter
    patchPhysicalTypes_.setSize(patchTypes_.size());

    forAll(boundaryCells_, patchI)
    {
        // resize - avoid invalid boundaries
        if (nPatchFaces[patchI] < boundaryCells_[patchI].size())
        {
            boundaryCells_[patchI].setSize(nPatchFaces[patchI]);
            boundaryFaces_[patchI].setSize(nPatchFaces[patchI]);
        }

        word origType = patchTypes_[patchI];
        patchPhysicalTypes_[patchI] = origType;

        if (origType == "SYMP")
        {
            patchTypes_[patchI] = symmetryPolyPatch::typeName;
            patchPhysicalTypes_[patchI] = patchTypes_[patchI];
        }
        else if (origType == "WALL")
        {
            patchTypes_[patchI] = wallPolyPatch::typeName;
            patchPhysicalTypes_[patchI] = patchTypes_[patchI];
        }
        else if (origType == "CYCL")
        {
            // incorrect. should be cyclicPatch but this
            // requires info on connected faces.
            patchTypes_[patchI] = cyclicPolyPatch::typeName;
            patchPhysicalTypes_[patchI] = patchTypes_[patchI];
        }
        else if (origType == "BAFF")
        {
            // incorrect. tag the patch until we get proper support.
            // set physical type to a canonical "baffle"
            patchTypes_[patchI] = emptyPolyPatch::typeName;
            patchPhysicalTypes_[patchI] = "baffle";
        }
        else
        {
            patchTypes_[patchI] = polyPatch::typeName;
        }

        Info<< "patch " << patchI
            << " (orig region " << origRegion[patchI]
            << " " << origType << ") is type '" << patchTypes_[patchI]
            << "' with name: " << patchNames_[patchI] << endl;
    }

    // cleanup
    mapToFoamCellId_.clear();
    cellShapes_.clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void starMeshReader::readMesh
(
    const fileName &prefix,
    const scalar scaleFactor
)
{
    readPoints(fileName(prefix + ".vrt"), scaleFactor);
    readCells(fileName(prefix + ".cel"));
    readBoundary(fileName(prefix + ".bnd"));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from reading a file
starMeshReader::starMeshReader
(
    const fileName& prefix,
    const Time& runtime,
    const scalar scaleFactor
)
:
    meshReader(prefix, runtime, scaleFactor),
    points_(0),
    cellShapes_(0),
    mapToFoamPointId_(0),
    mapToFoamCellId_(0)
{
    readMesh(prefix, scaleFactor);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

starMeshReader::~starMeshReader()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

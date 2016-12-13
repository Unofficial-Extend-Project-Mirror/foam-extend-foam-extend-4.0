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

\*---------------------------------------------------------------------------*/

#include "ElmerMeshWriter.H"

#include "foamTime.H"
#include "SortableList.H"
#include "OFstream.H"

#include "hexMatcher.H"
#include "wedgeMatcher.H"
#include "prismMatcher.H"
#include "pyrMatcher.H"
#include "tetWedgeMatcher.H"
#include "tetMatcher.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshWriters::Elmer::getCellTable()
{
    // read constant/polyMesh/propertyName
    IOList<label> ioList
    (
        IOobject
        (
            "cellTableId",
            "constant",
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    bool useCellZones = false;
    cellTableId_.setSize(mesh_.nCells(), -1);

    // get information from constant/polyMesh/cellTableId if possible
    if (ioList.headerOk())
    {
        if (ioList.size() == mesh_.nCells())
        {
            cellTableId_.transfer(ioList);

            if (cellTable_.empty())
            {
                Info<< "No cellTable information available" << endl;
            }
        }
        else
        {
            WarningIn("Elmer::getCellTable()")
                << ioList.objectPath() << " has incorrect number of cells "
                << " - use cellZone information"
                << endl;

            ioList.clear();
            useCellZones = true;
        }
    }
    else
    {
        useCellZones = true;
    }


    if (useCellZones)
    {
        if (cellTable_.empty())
        {
            Info<< "Created cellTable from cellZones" << endl;
            cellTable_ = mesh_;
        }

        // track if there are unzoned cells
        label nUnzoned = mesh_.nCells();

        // get the cellZone <-> cellTable correspondence
        Info<< "Matching cellZones to cellTable" << endl;

        forAll (mesh_.cellZones(), zoneI)
        {
            const cellZone& cZone = mesh_.cellZones()[zoneI];
            if (cZone.size())
            {
                nUnzoned -= cZone.size();

                label tableId = cellTable_.findIndex(cZone.name());
                if (tableId < 0)
                {
                    dictionary dict;

                    dict.add("Label", cZone.name());
                    dict.add("MaterialType", "fluid");
                    tableId = cellTable_.append(dict);
                }

                forAll (cZone, i)
                {
                    cellTableId_[cZone[i]] = tableId;
                }
            }
        }

        if (nUnzoned)
        {
            dictionary dict;

            dict.add("Label", "__unZonedCells__");
            dict.add("MaterialType", "fluid");
            label tableId = cellTable_.append(dict);

            forAll (cellTableId_, i)
            {
                if (cellTableId_[i] < 0)
                {
                    cellTableId_[i] = tableId;
                }
            }
        }
    }
}


Foam::label Foam::meshWriters::Elmer::getFaceType
(
    const label nvert,
    const word& zname
) const
{
    switch(nvert)
    {
        case 3:
            return ELMER_ETYPE_TRIA;

        case 4:
            return ELMER_ETYPE_QUAD;

        default:
            Info<< "Found a bad face with " << nvert
                << " vertices on zone " << zname
                << "." << endl;
            return ELMER_ETYPE_BAD;
    }
}


void Foam::meshWriters::Elmer::writeNames() const
{
    OFstream os("mesh.names");

    const cellZoneMesh& czones = mesh_.cellZones();
    const faceZoneMesh& fzones = mesh_.faceZones();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    label boundaryID = 0;


    Info<< "Writing " << os.name() << "." << endl;

    os << "! ----- Names for bodies -----" << nl;
    forAll (czones, zoneI)
    {
        os  << "$ " << czones[zoneI].name() << " = " << zoneI+1 << nl;
    }

    os << "! ----- Names for exterior boundaries -----" << nl;
    forAll(patches, patchI)
    {
        os  << "$ " << patches[patchI].name() << " = " << ++boundaryID << nl;
    }

    os << "! ----- Names for interior boundaries -----" << nl;
    forAll(fzones, fzoneI)
    {
        if (!faceZoneExcludePattern.match(fzones[fzoneI].name()))
        {
            os  << "$ " << fzones[fzoneI].name() << " = " << ++boundaryID << nl;
        }
    }
}


bool Foam::meshWriters::Elmer::writeHeader() const
{
    const pointField& points = mesh_.points();
    const cellList& cells = mesh_.cells();
    const faceList& faces = mesh_.faces();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const faceZoneMesh& fzones = mesh_.faceZones();
    label extZones = 0; // Count exterior boundary zones
    label extFaces = 0; // Count exterior boundary faces
    label intZones = 0; // Used interior boundary zones
    label intFaces = 0; // Used interior boundary faces

    // Count the different cell types
    hexMatcher hex;
    prismMatcher pri;
    pyrMatcher pyr;
    tetMatcher tet;

    label nHex = 0;
    label nPri = 0;
    label nPyr = 0;
    label nTet = 0;
    label nTri = 0;
    label nQua = 0;
    label nBad = 0;


    // Count all exterior boundary zone + faces + the types
    forAll(patches, patchI)
    {
        const label patchStart = patches[patchI].start();
        const label patchEnd = patchStart + patches[patchI].size();

        extZones++;
        extFaces += patches[patchI].size();

        // Count triangles and quads
        for(label faceI = patchStart; faceI < patchEnd; ++faceI)
        {
            if (faces[faceI].size() == 4)
            {
               nQua++;
            }
            else
            {
               nTri++;
            }
        }
    }

    // Count all interior boundary zones + faces with names not matching
    // the exclude pattern
    forAll(fzones, fzoneI)
    {
        if (faceZoneExcludePattern.match(fzones[fzoneI].name()))
        {
           continue;
        }

        const faceZone& fzone = fzones[fzoneI];

        intZones++;
        intFaces += fzone.size();

        // Count triangles and quads
        forAll(fzone, i)
        {
            if (faces[fzone[i]].size() == 4)
            {
                nQua++;
            }
            else
            {
                nTri++;
            }
        }
    }

    // Count the volume element types
    for(label cellI = 0; cellI < mesh_.nCells(); cellI++)
    {
        if (hex.isA(mesh_, cellI))
        {
            nHex++;
        }
        else if (tet.isA(mesh_, cellI))
        {
            nTet++;
        }
        else if (pyr.isA(mesh_, cellI))
        {
            nPyr++;
        }
        else if (pri.isA(mesh_, cellI))
        {
            nPri++;
        }
        else
        {
            nBad++; // Not a valid element type for Elmer
        }
    }

    OFstream os("mesh.header");

    os << points.size() << " "
       << cells.size()  << " "
       << extFaces+intFaces << nl;

    Info<< "Writing " << os.name() << "." << endl;

    // Save the number of different element types used
    os  << (nTet>0) + (nPri>0) + (nPyr>0) +
           (nHex>0) + (nQua>0) + (nTri>0) + (nBad>0) << endl;

    // Save type and count if count > 0
    if (nTet > 0) os << ELMER_ETYPE_TET   << " " << nTet << endl;
    if (nPri > 0) os << ELMER_ETYPE_PRISM << " " << nPri << endl;
    if (nPyr > 0) os << ELMER_ETYPE_PYRAM << " " << nPyr << endl;
    if (nHex > 0) os << ELMER_ETYPE_HEX   << " " << nHex << endl;
    if (nQua > 0) os << ELMER_ETYPE_QUAD  << " " << nQua << endl;
    if (nTri > 0) os << ELMER_ETYPE_TRIA  << " " << nTri << endl;
    if (nBad > 0) os << ELMER_ETYPE_BAD   << " " << nBad << endl;

    Info<< "Mesh statistics:" << nl
        << "    Nodes          : " << points.size()  << nl
        << "    Elements       : " << cells.size() << nl
        << "        Hexahedra  : " << nHex << nl
        << "        Wedges     : " << nPri << nl
        << "        Pyramids   : " << nPyr << nl
        << "        Tetrahedra : " << nTet << nl
        << "        Polyhedra  : " << nBad << nl
        << "    Regions        : " << mesh_.cellZones().size() << nl
        << "    Ext. boundaries: " << extZones << ", "
        << extFaces << " faces" << nl
        << "    Int. boundaries: " << intZones << ", "
        << intFaces << " faces" << nl
        << "        Tri faces  : " << nTri << nl
        << "        Quad faces : " << nQua << nl
        << endl;

    return !nBad;
}


void Foam::meshWriters::Elmer::writeNodes() const
{
    OFstream os("mesh.nodes");

    // Set the precision of the points data to 10
    os.precision(10);

    // Force decimal point for Fortran90 input style
    os.setf(std::ios::showpoint);

    const pointField& points = mesh_.points();

    Info<< "Writing " << os.name() << ": scale=" << scaleFactor_ << endl;

    forAll(points, ptI)
    {
        os << ptI + 1 << " -1 "
           << scaleFactor_ * points[ptI].x() << " "
           << scaleFactor_ * points[ptI].y() << " "
           << scaleFactor_ * points[ptI].z() << nl;
    }
}


void Foam::meshWriters::Elmer::writeElements() const
{
    OFstream os("mesh.elements");

    // map foam cellModeller index -> Elmer element types
    Map<label> shapeLookupIndex;
    shapeLookupIndex.insert(tetModel->index(), ELMER_ETYPE_TET  );
    shapeLookupIndex.insert(prismModel->index(), ELMER_ETYPE_PRISM);
    shapeLookupIndex.insert(pyrModel->index(), ELMER_ETYPE_PYRAM);
    shapeLookupIndex.insert(hexModel->index(), ELMER_ETYPE_HEX  );

    const cellShapeList& shapes = mesh_.cellShapes();
    const cellList& cells  = mesh_.cells();

    Info<< "Writing " << os.name() << "." << endl;

    forAll(cells, cellId)
    {
        const cellShape& shape = shapes[cellId];
        label mapIndex = shape.model().index();

        // a registered primitive type
        if (shapeLookupIndex.found(mapIndex))
        {
            os << cellId+1 << " "
               << cellTableId_[cellId] << " "
               << shapeLookupIndex[mapIndex];

            const labelList& vrtList = shapes[cellId];
            forAll(vrtList, i)
            {
                os << " " << vrtList[i] + 1;
            }
            os << endl;
        }
        else
        {
            Info<< "***WARNING: polyhedron " << cellId << " ignored." << endl;
        }
    }
}


void Foam::meshWriters::Elmer::writeBoundary() const
{
    OFstream os("mesh.boundary");

    const faceList&  faces = mesh_.faces();
    const labelList& owner = mesh_.faceOwner();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    label bfaceID = 0;
    label boundaryID = 0;


    Info<< "Writing " << os.name() << "." << endl;

    // Write all patches == exterior boundary faces
    forAll(patches, patchI)
    {
        const word& pname = patches[patchI].name();
        const label patchStart = patches[patchI].start();
        const label patchEnd = patchStart + patches[patchI].size();

        boundaryID++;
        for(label faceI = patchStart; faceI < patchEnd; ++faceI)
        {
            const labelList& vrtList = faces[faceI];

            os << ++bfaceID
               << " " << boundaryID
               << " " << owner[faceI] + 1  // ID of parent cell
               << " 0 "                    // 0 == exterior boundary face
               << getFaceType(vrtList.size(),pname);

            forAll(vrtList, i)
            {
                os << " " << vrtList[i] + 1;
            }
            os << endl;
        }
    }


    // Write face zones (== interior boundary faces) with names not matching
    // the exclude pattern
    const faceZoneMesh& fzones = mesh_.faceZones();
    forAll(fzones, fzoneI)
    {
        if(faceZoneExcludePattern.match(fzones[fzoneI].name()))
        {
           continue;
        }

        const faceZone& fzone = fzones[fzoneI];
        const labelList& master = fzone.masterCells();
        const labelList& slave = fzone.slaveCells();
        const word& zname = fzone.name();

        boundaryID++;
        forAll(fzone, i)
        {
            const labelList& vrtList = faces[fzone[i]];

            os << ++bfaceID
               << " " << boundaryID
               << " " << master[i] + 1 // ID of parent cell
               << " " << slave[i]  + 1 // ID of neigbour cell
               << " " << getFaceType(vrtList.size(),zname);

            forAll(vrtList, j)
            {
                os << " " << vrtList[j] + 1;
            }
            os << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshWriters::Elmer::Elmer
(
    const polyMesh& mesh,
    const wordRe& excludePattern,
    const scalar scaleFactor
)
:
    meshWriter(mesh, scaleFactor),
    faceZoneExcludePattern(excludePattern)
{
    boundaryRegion_.readDict(mesh_);
    cellTable_.readDict(mesh_);
    getCellTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshWriters::Elmer::~Elmer()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::meshWriters::Elmer::write(const fileName& dirName) const
{
    fileName baseName(dirName);

    if (baseName.empty())
    {
        baseName = meshWriter::defaultMeshName;

        if
        (
            mesh_.time().timeName() != "0"
         && mesh_.time().timeName() != "constant"
        )
        {
            baseName += "_" + mesh_.time().timeName();
        }
    }

    // Create the mesh directory (elemerMesh), chdir into and cleanup
    mkDir(baseName);
    chDir(baseName);
    rm("mesh.header");
    rm("mesh.nodes");
    rm("mesh.elements");
    rm("mesh.boundary");
    rm("mesh.names");

    bool success = writeHeader();
    if (success)
    {
        writeNodes();
        writeElements();
        writeBoundary();
        writeNames();
    }
    else
    {
        rm("mesh.header");
    }

    chDir("..");

    return success;
}

// ************************************************************************* //

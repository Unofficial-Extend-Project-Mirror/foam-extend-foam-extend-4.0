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
    Transcribe cellZones to 0/cellTableId volume field

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "volFields.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static Map<dictionary> cellTable_;

// writing the cell table as a dictionary of dictionaries eases reading later
// - generate the keys ourselves, since the 'Label' entries are not
//   guaranteed to be unique
void writeCellTable(const polyMesh& mesh)
{
    IOdictionary ioObj
    (
        IOobject
        (
            "cellTable",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        )
    );

    ioObj.note() = "persistent data for pro-STAR <-> foam translation";

    Info << "Writing cellTable to " << ioObj.objectPath() << endl;

    OFstream os(ioObj.objectPath());
    ioObj.writeHeader(os);

    labelList toc = cellTable_.toc();
    forAll(toc, i)
    {
        label key = toc[i];

        os << word("cellTableId_" + Foam::name(key));
        cellTable_[key].write(os, true);
        os << endl;
    }

    os  << "// ************************************************************************* //"
    << endl;
}

// write cellZones as volume field cellTableId
void cellZoneToCellTableId
(
    const fvMesh & mesh
)
{
    cellTable_.clear();
    const word & propertyName = "cellTableId";

    volScalarField calculatedField
    (
        IOobject
        (
            propertyName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(propertyName, dimless, 0.0)
    );

    const Map<label> & zoneMap = mesh.cellZones().zoneMap();
    wordList zoneNames = mesh.cellZones().names();

    // start zoned cells at cell type 1
    label typeOffset = 1;

    // fewer zoned cells than total cells
    // - leave unzoned cells as type 1 and start zoned cells at cell type 2
    if (zoneMap.size() < mesh.nCells())
    {
        typeOffset = 2;
        label id = 1;

        dictionary dict;
        dict.add("Id", id);
        dict.add("Label", "__unzonedCells__");
        cellTable_.insert(id, dict);
    }

    forAll(zoneNames, zoneI) 
    {
        label id = zoneI + typeOffset;

        dictionary dict;
        dict.add("Id", id);
        dict.add("Label", zoneNames[zoneI]);
        cellTable_.insert(id, dict);
    }

    writeCellTable(mesh);

    forAllConstIter(Map<label>, zoneMap, iter)
    {
        calculatedField[iter.key()] = iter() + typeOffset;
    }
    
    Info<<"Writing the " << mesh.cellZones().size() << " cellZones as "
        << propertyName << " to " << calculatedField.objectPath()
        << endl;

    // NOTE:
    // the cellTableId is an integer and almost always < 1000, thus ASCII
    // will be compacter than binary and makes external scripting easier
    //
    calculatedField.write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"

    bool firstCheck = true;

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        polyMesh::readUpdateState state = mesh.readUpdate();

        if (firstCheck || state != polyMesh::UNCHANGED)
        {
            cellZoneToCellTableId(mesh);
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

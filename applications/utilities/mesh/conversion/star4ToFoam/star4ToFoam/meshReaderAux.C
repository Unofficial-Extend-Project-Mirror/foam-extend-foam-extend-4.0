/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) ...
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

Class
    meshReader

Description

\*----------------------------------------------------------------------------*/
# include "meshReader.H"
# include "Time.H"
# include "OFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// writing the cell table as a dictionary of dictionaries eases reading later
// - generate the keys ourselves, since the 'Label' entries are not
//   guaranteed to be unique
void Foam::meshReader::writeCellTable(const polyMesh& mesh)
{
    IOdictionary ioObj
    (
        IOobject
        (
            "cellTable",
            runTime_.constant(),
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


// write interface (baffle) mapping
void Foam::meshReader::writeInterfaces(const polyMesh& mesh)
{
    IOList<labelList> ioObj
    (
        IOobject
        (
            "interfaces",
            runTime_.constant(),
            fvMesh::meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        )
    );

    ioObj.note() = "as yet unsupported interfaces (baffles)";

    Info << "Writing interfaces to " << ioObj.objectPath() << endl;

    OFstream os(ioObj.objectPath());
    ioObj.writeHeader(os);

    os  << interfaces_
	<< "// ************************************************************************* //"
	<< endl;

}

// ************************************************************************* //

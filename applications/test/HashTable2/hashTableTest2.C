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

#include "HashTable.H"
#include "HashPtrTable.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    HashTable<label, Foam::string> testTable(0);

    testTable.insert("kjhk", 10);
    testTable.insert("kjhk2", 12);

    Info<< testTable << endl;
    Info<< testTable.toc() << endl;

    HashTable<label, label, Hash<label> > testTable2(10);

    testTable2.insert(3, 10);
    testTable2.insert(5, 12);
    testTable2.insert(7, 16);

    Info<< testTable2 << endl;
    Info<< testTable2.toc() << endl;

    HashTable<label, label, Hash<label> > testTable3(1);
    testTable3.transfer(testTable2);

    Info<< testTable2 << endl;
    Info<< testTable2.toc() << endl;

    Info<< testTable3 << endl;
    Info<< testTable3.toc() << endl;

    Foam::HashPtrTable<label, Foam::string> testPtrTable(0);
    testPtrTable.insert("kjhkjh", new label(10));

    Info<< testPtrTable.toc() << endl;

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

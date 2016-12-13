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

#include "dynamicFvMesh.H"
#include "foamTime.H"
#include "dlLibraryTable.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dynamicFvMesh> Foam::dynamicFvMesh::New(const IOobject& io)
{
    wordList libNames(1);
    libNames[0]=word("mesquiteMotionSolver");

    forAll(libNames,i) {
        const word libName("lib"+libNames[i]+".so");

        bool libLoaded = false;

        dlLibraryTable& ll = dlLibraryTable::loadedLibraries;

        forAllConstIter(dlLibraryTable, ll, llI)
        {
            if (ll(llI.key()) == libName)
            {
                libLoaded = true;
                break;
            }
        }

        if (!libLoaded)
        {
            bool ok=dlLibraryTable::open(libName);
            if(!ok) {
                WarningIn("dynamicFvMesh::New(const IOobject& io)")
                    << "Loading of dynamic mesh library " << libName
                    << " unsuccesful. Some dynamic mesh  methods may not be "
                    << " available"
                    << endl;
            }
        }
    }

    // Enclose the creation of the dynamicMesh to ensure it is
    // deleted before the dynamicFvMesh is created otherwise the dictionary
    // is entered in the database twice
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            io.time().constant(),
            (io.name() == dynamicFvMesh::defaultRegion ? "" : io.name() ),
            io.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word dynamicFvMeshTypeName(dynamicMeshDict.lookup("dynamicFvMesh"));

    Info<< "Selecting dynamicFvMesh " << dynamicFvMeshTypeName << endl;

    dlLibraryTable::open
    (
        dynamicMeshDict,
        "dynamicFvMeshLibs",
        IOobjectConstructorTablePtr_
    );

    if (!IOobjectConstructorTablePtr_)
    {
        FatalErrorIn
        (
            "dynamicFvMesh::New(const IOobject&)"
        )   << "dynamicFvMesh table is empty"
            << exit(FatalError);
    }

    IOobjectConstructorTable::iterator cstrIter =
        IOobjectConstructorTablePtr_->find(dynamicFvMeshTypeName);

    if (cstrIter == IOobjectConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "dynamicFvMesh::New(const IOobject&)"
        )   << "Unknown dynamicFvMesh type " << dynamicFvMeshTypeName
            << endl << endl
            << "Valid dynamicFvMesh types are :" << endl
            << IOobjectConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<dynamicFvMesh>(cstrIter()(io));
}


// ************************************************************************* //

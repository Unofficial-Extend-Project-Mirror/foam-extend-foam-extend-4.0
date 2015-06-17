/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "engineTopoChangerMesh.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::engineTopoChangerMesh> Foam::engineTopoChangerMesh::New
(
    const Foam::IOobject& io
)
{
    word engineTopoChangerMeshTypeName;

    // Enclose the creation of the engineGeometry to ensure it is
    // deleted before the engineTopoChangerMesh is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary engineGeometryDict
        (
            IOobject
            (
                "engineGeometry",
                io.time().constant(),
                io.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        engineGeometryDict.lookup("engineTopoChangerMesh") >> engineTopoChangerMeshTypeName;
    }

    Info<< "Selecting engineTopoChangerMesh " << engineTopoChangerMeshTypeName << endl;

    IOobjectConstructorTable::iterator cstrIter =
        IOobjectConstructorTablePtr_->find(engineTopoChangerMeshTypeName);

    if (cstrIter == IOobjectConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "engineTopoChangerMesh::New(const IOobject&)"
        )   << "Unknown engineTopoChangerMesh type " << engineTopoChangerMeshTypeName
            << endl << endl
            << "Valid engineTopoChangerMesh types are :" << endl
            << IOobjectConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<engineTopoChangerMesh>(cstrIter()(io));
}


// ************************************************************************* //

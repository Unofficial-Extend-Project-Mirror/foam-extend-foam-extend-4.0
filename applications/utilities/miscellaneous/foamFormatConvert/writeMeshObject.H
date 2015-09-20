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

Global
    writeMeshObject

Description
    Write a mesh object in the format specified in controlDict.

\*---------------------------------------------------------------------------*/

#ifndef writeMeshObject_H
#define writeMeshObject_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
inline bool writeMeshObject
(
    const word& name,
    const fileName& meshDir,
    Time& runTime
)
{
    IOobject io
    (
        name,
        runTime.timeName(),
        meshDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    bool writeOk = false;

    if (io.headerOk())
    {
        Info<< "        Reading " << io.headerClassName()
            << " : " << name << endl;

        // Switch off type checking (for reading e.g. faceZones as
        // generic list of dictionaries).
        const word oldTypeName = T::typeName;
        const_cast<word&>(T::typeName) = word::null;

        T meshObject(io);

        const_cast<word&>(T::typeName) = oldTypeName;
        // Fake type back to what was in field
        const_cast<word&>(meshObject.type()) = io.headerClassName();

        Info<< "        Writing " << name << endl;
        writeOk = meshObject.regIOobject::write();
    }

    return writeOk;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

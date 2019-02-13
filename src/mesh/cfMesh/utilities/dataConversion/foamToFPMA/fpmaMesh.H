/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    fpmaMesh

Description

SourceFiles
    fpmaMesh.C

\*---------------------------------------------------------------------------*/

#ifndef fpmaMesh_H
#define fpmaMesh_H

#include "polyMeshGen.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fpmaMesh Declaration
\*---------------------------------------------------------------------------*/

class fpmaMesh
{
    // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fpmaMesh(const fpmaMesh&);

        //- Disallow default bitwise assignment
        void operator=(const fpmaMesh&);

        void writePoints(OFstream& fpmaGeometryFile) const;

        void writeFaces(OFstream& fpmaGeometryFile) const;

        void writeCells(OFstream& fpmaGeometryFile) const;

        void writeSubsets(OFstream& fpmaGeometryFile) const;

public:

    // Public data

        const polyMeshGen& mesh_;

    // Constructors

        //- Construct from polyMeshGen
        fpmaMesh(const polyMeshGen& mesh);


    // Destructor

        ~fpmaMesh();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write

            void write(OFstream& fpmaGeometryFile) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

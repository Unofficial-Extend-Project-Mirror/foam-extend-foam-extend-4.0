/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

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

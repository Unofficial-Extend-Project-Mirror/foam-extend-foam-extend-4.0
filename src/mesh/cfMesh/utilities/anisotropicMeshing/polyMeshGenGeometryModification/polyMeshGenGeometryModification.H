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
    polyMeshGenGeometryModification

Description
    Modifies coordinates of points in the mesh to achieve anistropy

SourceFiles
    polyMeshGenGeometryModification.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenGeometryModification_H
#define polyMeshGenGeometryModification_H

#include "coordinateModification.H"
#include "polyMeshGenModifier.H"
#include "typeInfo.H"
#include "coordinateModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

/*---------------------------------------------------------------------------*\
                Class polyMeshGenGeometryModification Declaration
\*---------------------------------------------------------------------------*/

class polyMeshGenGeometryModification
{
    // Private data
        //- reference to polyMeshGen
        polyMeshGen& mesh_;

        //- length of box sides
        const dictionary& meshDict_;

        //- contruct coordinate modification
        coordinateModifier* coordinateModifierPtr_;

        //- is mofdification active
        bool modificationActive_;


    // Private member functions
        //- check existence of geometry modifiers
        void checkModification();

        //- disable bitwise copy construct
        polyMeshGenGeometryModification(const polyMeshGenGeometryModification&);

        //- disable bitwise assignment
        void operator=(const polyMeshGenGeometryModification&);

public:

    // Constructors

        //- Construct from polyMeshGen and dictionary
        polyMeshGenGeometryModification(polyMeshGen&, const dictionary&);

    // Destructor
        ~polyMeshGenGeometryModification();

    // Member Functions

        //- is geometry modification active
        bool activeModification() const;

        //- modify coordinates
        void modifyGeometry();

        //- revert modification of coorinates
        void revertGeometryModification();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

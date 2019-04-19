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
    partTriMeshSimplex

Description
    A simplex which is used for smoothing purposes

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef partTriMeshSimplex_H
#define partTriMeshSimplex_H

#include "partTriMesh.H"
#include "triFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class partTriMesh;

/*---------------------------------------------------------------------------*\
                           Class partTriMeshSimplex Declaration
\*---------------------------------------------------------------------------*/

class partTriMeshSimplex
{
    // Private data
        //- points making the simplex
        DynList<point, 32> pts_;

        //- triangles making the simplex
        DynList<triFace, 32> trias_;

    public:

    // Constructors
        //- Construct from partTriMesh and point label
        partTriMeshSimplex(const partTriMesh& tm, const label pI);

    // Destructor
        ~partTriMeshSimplex();

    // Member functions
        //- return points
        inline DynList<point, 32>& pts()
        {
            return pts_;
        }

        //- return points
        inline const DynList<point, 32>& pts() const
        {
            return pts_;
        }

        //- return triangles
        inline const DynList<triFace, 32>& triangles() const
        {
            return trias_;
        }

        //- return centre point coordinates
        inline const point& centrePoint() const
        {
            return pts_[trias_[0][0]];
        }

        //- return the normal of the simplex
        vector normal() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

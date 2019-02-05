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
    partTetMeshSimplex

Description
    A simplex which is used for smoothing purposes

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef partTetMeshSimplex_H
#define partTetMeshSimplex_H

#include "partTetMesh.H"
#include "parPartTet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class partTetMesh;

/*---------------------------------------------------------------------------*\
                           Class partTetMeshSimplex Declaration
\*---------------------------------------------------------------------------*/

class partTetMeshSimplex
{
    // Private data
        //- points making the simplex
        DynList<point, 128> pts_;

        //- tets making the simplex
        DynList<partTet, 128> tets_;

    public:

    // Constructors
        //- Construct from partTetMeshSimplex and point label
        partTetMeshSimplex(const partTetMesh& tm, const label pI);

        //- Construct from the list of parPartTet and point label
        partTetMeshSimplex
        (
            const DynList<parPartTet>& pt,
            const label gpI
        );

        //- contruct from a list of points and a list of partTets
        //- ad a point label
        partTetMeshSimplex
        (
            const DynList<point, 128>& pts,
            const DynList<partTet, 128>& tets,
            const label pointI
        );

    // Destructor
        ~partTetMeshSimplex();

    // Member functions
        //- return points
        inline DynList<point, 128>& pts()
        {
            return pts_;
        }

        //- return points
        inline const DynList<point, 128>& pts() const
        {
            return pts_;
        }

        //- return tets
        inline const DynList<partTet, 128>& tets() const
        {
            return tets_;
        }

        //- return centre point coordinates
        inline const point& centrePoint() const
        {
            return pts_[tets_[0][3]];
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "partTetMeshSimplexI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    simplexSmoother

Description
    Mesh smoothing without any topological changes. The vertex is put into
    location which minimises the objective function.

SourceFiles
    simplexSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef simplexSmoother_H
#define simplexSmoother_H

#include "point.H"
#include "DynList.H"
#include "partTet.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class partTetMeshSimplex;

/*---------------------------------------------------------------------------*\
                    Class simplexSmoother Declaration
\*---------------------------------------------------------------------------*/

//- class for volume optimizer
class simplexSmoother
{
    protected:

    // Protected data
        //- mesh points
        DynList<point, 128>& points_;

        //- list tets around the given vertex
        const DynList<partTet, 128>& tets_;

        //- label of the point
        const label pointI_;

        //- bound box
        boundBox bb_;

private:

    // Private member functions
        //- Disallow default bitwise copy construct
        simplexSmoother(const simplexSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const simplexSmoother&);

    public:

    // Constructor
        //- construct from partTetMeshSimplex
        simplexSmoother(partTetMeshSimplex& simplex);

    // Destructor
        virtual ~simplexSmoother();

    // Member functions
        //- improve the position of the centre node
        virtual void optimizeNodePosition(const scalar tol = 0.001) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

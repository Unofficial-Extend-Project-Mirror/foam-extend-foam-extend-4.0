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
    triSurface2DCheck

Description


SourceFiles
    triSurface2DCheck.C

\*---------------------------------------------------------------------------*/

#ifndef triSurface2DCheck_H
#define triSurface2DCheck_H

#include "triSurf.H"
#include "symmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class triSurface2DCheck Declaration
\*---------------------------------------------------------------------------*/

class triSurface2DCheck
{
    // Private data
        //- reference to triSurf
        const triSurf& surf_;

        //- covariance matrix
        symmTensor covarianceMatrix_;

    // Private member functions
        //- create covariance matrix
        void createCovarianceMatrix();

        //- Disallow default bitwise copy construct
        triSurface2DCheck(const triSurface2DCheck&);

        //- Disallow default bitwise assignment
        void operator=(const triSurface2DCheck&);

public:

    // Constructors

        //- Construct from octree
        triSurface2DCheck(const triSurf& surface);

    // Destructor

        ~triSurface2DCheck();


    // Member Functions

        //- checks if the surface is a 2D triangulation
        bool is2DSurface() const;

        //- create subset containing invalid facets
        void createSubsets();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

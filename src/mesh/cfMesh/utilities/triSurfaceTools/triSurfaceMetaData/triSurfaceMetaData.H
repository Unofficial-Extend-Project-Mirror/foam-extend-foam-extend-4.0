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
    triSurfaceMetaData

Description
    Generates meta data of the surface mesh. It contains:
    The number of points
    The number of triangles
    The number of patches
    The number of feature edges
    The number of subsets, names, and the number of elements in each subset

SourceFiles
    triSurfaceMetaData.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceMetaData_H
#define triSurfaceMetaData_H

#include "triSurf.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class triSurfaceMetaData Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceMetaData
{
    // Private data
        //- const reference to triSurf
        const triSurf& surf_;

        //- patches/subsets for removal
        dictionary metaDict_;

    // Private member functions
        //- create dictionary with meta data
        void createMetaData();

        //- Disallow default bitwise copy construct
        triSurfaceMetaData(const triSurfaceMetaData&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceMetaData&);

public:

    // Constructors

        //- Construct from triSurf
        triSurfaceMetaData(const triSurf& surface);

    // Destructor

        ~triSurfaceMetaData();


    // Member Functions

        //- return a constant reference to meta data
        const dictionary& metaData() const
        {
            return metaDict_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

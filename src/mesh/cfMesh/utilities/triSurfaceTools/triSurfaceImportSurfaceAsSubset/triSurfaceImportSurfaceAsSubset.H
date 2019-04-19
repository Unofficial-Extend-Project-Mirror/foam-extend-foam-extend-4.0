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
    triSurfaceImportSurfaceAsSubset

Description
    Creates a subset in the master surface consisting of facets which are
    near the other surface

SourceFiles
    triSurfaceImportSurfaceAsSubset.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceImportSurfaceAsSubset_H
#define triSurfaceImportSurfaceAsSubset_H

#include "triSurf.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class meshOctree;

/*---------------------------------------------------------------------------*\
             Class triSurfaceImportSurfaceAsSubset Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceImportSurfaceAsSubset
{
    // Private data
        //- reference to triSurf
        triSurf& surf_;

        //- pointer to meshOctree, needed for searching on the master surface
        meshOctree* octreePtr_;

    // Private member functions
        void createOctree(const triSurf&, meshOctree&);

        //- Disallow default bitwise copy construct
        triSurfaceImportSurfaceAsSubset(const triSurfaceImportSurfaceAsSubset&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceImportSurfaceAsSubset&);

public:

    // Constructors

        //- Construct from triSurf
        triSurfaceImportSurfaceAsSubset(triSurf& surface);

    // Destructor

        ~triSurfaceImportSurfaceAsSubset();


    // Member Functions

        //- finds the nearest faces in the surface to the import surf
        //- and creates a subset
        void addSurfaceAsSubset
        (
            const triSurf& importSurf,
            const word& subsetName,
            const scalar angleTol = 5.*M_PI/180.
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

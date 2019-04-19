/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Class
    Foam::movingImmersedBoundary

Description
    Moving immesed boundary.  Motion is prescribed using a solid body motion
    function.

SourceFiles
    movingImmersedBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef movingImmersedBoundary_H
#define movingImmersedBoundary_H

#include "polyMesh.H"
#include "solidBodyMotionFunction.H"
#include "triSurfaceMesh.H"
#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class movingImmersedBoundary Declaration
\*---------------------------------------------------------------------------*/

class movingImmersedBoundary
{
    // Private data

        //- Patch name
        const word name_;

        //- Mesh reference
        const polyMesh& mesh_;

        //- Overset region motion control function
        autoPtr<solidBodyMotionFunction> sbmfPtr_;

        //- Reference tri surface mesh position
        triSurfaceMesh refIbSurface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        movingImmersedBoundary(const movingImmersedBoundary&);

        //- Disallow default bitwise assignment
        void operator=(const movingImmersedBoundary&);


public:

    // Constructors

        //- Construct from dictionary
        movingImmersedBoundary
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~movingImmersedBoundary();


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Move immersed boundary patch points
        void movePoints() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

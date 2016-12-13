/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::movingBodyTopoFvMesh

Description
    Sample topoChangerFvMesh that moves a cell zone in a given direction,
    introducing/removing layers around it.

    Cells in the movingCells zone shall be moved given the prescribed velocity
    and will be bounded in "front" and "back" by other cell zones.
    Layer addition/removal interfaces are inserted at boundaries between the
    moving zone and front and back, pointing outside of the moving cell zone

Author and rewrite
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    movingBodyTopoFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef movingBodyTopoFvMesh_H
#define movingBodyTopoFvMesh_H

#include "topoChangerFvMesh.H"
#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class movingBodyTopoFvMesh Declaration
\*---------------------------------------------------------------------------*/

class movingBodyTopoFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Dynamic mesh dictionary
        dictionary dict_;

        //- Moving cell zone name
        word movingCellsName_;

        //- Front face zone name
        word frontFacesName_;

        //- Back face zone name
        word backFacesName_;

        //- Motion control function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;

        //- Vertex motion mask
        scalarField motionMask_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        movingBodyTopoFvMesh(const movingBodyTopoFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const movingBodyTopoFvMesh&);


        //- Add mixer zones and modifiers
        void addZonesAndModifiers();

        //- Mark motion vertices
        tmp<scalarField> calcMotionMask() const;


public:

    //- Runtime type information
    TypeName("movingBodyTopoFvMesh");


    // Constructors

        //- Construct from database
        explicit movingBodyTopoFvMesh(const IOobject& io);


    // Destructor

        virtual ~movingBodyTopoFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

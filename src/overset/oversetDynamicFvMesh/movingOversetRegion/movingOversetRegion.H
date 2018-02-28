/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Foam::movingOversetRegion

Description
    Moving region in overset motion hierarchy.  Motion is prescribed using a
    solid body motion function and moving cells are chosen as a set of
    cell zones.

SourceFiles
    movingOversetRegion.C

\*---------------------------------------------------------------------------*/

#ifndef movingOversetRegion_H
#define movingOversetRegion_H

#include "fvMesh.H"
#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class movingOversetRegion Declaration
\*---------------------------------------------------------------------------*/

class movingOversetRegion
{
    // Private data

        //- Name
        const word name_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Overset region motion control function
        autoPtr<solidBodyMotionFunction> sbmfPtr_;

        //- Moving cell zone names
        wordList movingZoneNames_;

        //- Region motion mask
        mutable scalarField* motionMaskPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        movingOversetRegion(const movingOversetRegion&);

        //- Disallow default bitwise assignment
        void operator=(const movingOversetRegion&);


        //- Calculate motion mask
        void calcMotionMask() const;

        // Clear storage
        void clearOut();


public:

    // Constructors

        //- Construct from dictionary
        movingOversetRegion
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~movingOversetRegion();

    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return motion mask
        const scalarField& motionMask() const;

        //- Return motion increment
        const tmp<pointField> motionIncrement
        (
            const pointField& undisplacedPoints
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

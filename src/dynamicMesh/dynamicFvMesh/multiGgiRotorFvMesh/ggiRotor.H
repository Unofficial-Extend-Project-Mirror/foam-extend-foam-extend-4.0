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
    ggiRotor

Description
    Single rotating region with GGI interfaces

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    ggiRotor.C

\*---------------------------------------------------------------------------*/

#ifndef ggiRotor_H
#define ggiRotor_H

#include "cylindricalCS.H"
#include "Switch.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class regionSplit;
class polyTopoChanger;

/*---------------------------------------------------------------------------*\
                        Class ggiRotor Declaration
\*---------------------------------------------------------------------------*/

class ggiRotor
{
    // Private data

        //- Name
        const word name_;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Coordinate system
        cylindricalCS cs_;

        // - Rotational speed in rotations per minute (rpm)
        const scalar rpm_;

        //- Moving cells zone
        const word movingCellsZoneName_;

        //- Moving patches: for zones
        const wordList movingPatches_;

        //- Static patches: for zones
        const wordList staticPatches_;

        //- Markup field for points.  Moving points marked with 1
        mutable scalarField* movingPointsMaskPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ggiRotor(const ggiRotor&);

        //- Disallow default bitwise assignment
        void operator=(const ggiRotor&);


        //- Return coordinate system
        const cylindricalCS& cs() const
        {
            return cs_;
        }

        //- Calculate moving mask
        void calcMovingMask() const;

        //- Return moving points mask
        const scalarField& movingPointsMask() const;

        //- Clear moving points mask
        void clearPointMask();


public:

    // Constructors

        //- Construct from  dictionary
        ggiRotor
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ggiRotor();


    // Member Functions

        //- Return accumulative point motion
        virtual tmp<vectorField> pointMotion() const;

        //- Update topology
        void updateTopology();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::trajectoryCollisionModel

Description
    Collision model by P.J. O'Rourke.
    Implemented as described in the KIVA manual

SourceFiles
    trajectoryModel.C

\*---------------------------------------------------------------------------*/

#ifndef trajectoryCollisionModel_H
#define trajectoryCollisionModel_H

#include "collisionModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class trajectoryCollisionModel Declaration
\*---------------------------------------------------------------------------*/

class trajectoryCollisionModel
:
    public collisionModel
{

private:
    // Reference to the cell volumes
        dictionary coeffsDict_;

        scalar cSpace_;
        scalar cTime_;
        Switch coalescence_;

public:

    //- Runtime type information
    TypeName("trajectory");


    // Constructors

        //- Construct from components
        trajectoryCollisionModel
        (
            const dictionary& dict,
            spray& sm,
            Random& rndGen
        );

    // Destructor

        ~trajectoryCollisionModel();


    // Member Functions

        void collideParcels(const scalar dt) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

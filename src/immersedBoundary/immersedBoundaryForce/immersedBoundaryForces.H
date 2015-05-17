/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Foam::immersedBoundaryForces

Description
    Calculates the forces and moments by integrating the pressure and
    skin-friction forces over a given list of immersed boundary patches.

    Member function calcForcesMoment() calculates and returns the forces and
    moments on immersed boundary patches.

Author
    Hrvoje Jasak.  All rights reserved.

SourceFiles
    immersedBoundaryForces.C
    IOimmersedBoundaryForces.H

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryForces_H
#define immersedBoundaryForces_H

#include "forces.H"
#include "fvPatchFields.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class immersedBoundaryForces Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryForces
:
    public forces
{
public:

    //- Runtime type information
    TypeName("immersedBoundaryForces");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        immersedBoundaryForces
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~immersedBoundaryForces();


    // Member Functions

        //- Calculate and return forces and moment
        virtual forcesMoments calcForcesMoment() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

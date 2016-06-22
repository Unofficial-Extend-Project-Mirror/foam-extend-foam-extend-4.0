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
    Foam::reflectParcel

Description
    The parcels are reflected back from the wall in an inelastic collision.

\*---------------------------------------------------------------------------*/

#ifndef reflectParcel_H
#define reflectParcel_H

#include "wallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class reflectParcel Declaration
\*---------------------------------------------------------------------------*/

class reflectParcel
:
    public wallModel
{

private:
    // Private data

        const volVectorField& U_;

        dictionary coeffsDict_;

    // model constants

    //- Elasticity constant v^2_out = elasticity * v^2_in
        scalar elasticity_;

public:

    //- Runtime type information
    TypeName("reflect");


    // Constructors

        //- Construct from components
        reflectParcel
        (
            const dictionary& dict,
            const volVectorField& U,
            spray& sm
        );


    // Destructor

        virtual ~reflectParcel();


    // Member Functions

        //- Return true if parcel is to be kept, and false if it is to be
        // removed
        bool wallTreatment
        (
            parcel& parcel,
            const label facei
        ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

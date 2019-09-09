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
    inverseDistanceInterpolation

Description
    Inverse distance overset interpolation algorithm
    For each acceptor, the value is interpolated based on inverse distance from
    a set of donors (master donor and neighbouring donors).

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

SourceFiles
    inverseDistanceInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef inverseDistanceInterpolation_H
#define inverseDistanceInterpolation_H

#include "oversetInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class inverseDistanceInterpolation Declaration
\*---------------------------------------------------------------------------*/

class inverseDistanceInterpolation
:
    public oversetInterpolation
{
    // Private data

        // Demand-driven data

            //- Weights for all donors of all acceptors in all regions
            // Note: for each acceptor in each region, we have a list of scalars
            // for all donors: master + neighbouring. Master weight is the first
            // [index 0], followed by neighbouring weights.
            mutable scalarFieldFieldList* weightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        inverseDistanceInterpolation(const inverseDistanceInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const inverseDistanceInterpolation&);


        // Functions for calculating/clearing weights

            //- Calculate weights
            void calcWeights() const;

            //- Clear weights
            void clearWeights() const;


public:

    //- Runtime type information
    TypeName("inverseDistance");


    // Constructors

        //- Construct from dictionary
        inverseDistanceInterpolation
        (
            const oversetMesh& overset,
            const word& name
        );


    //-Destructor
    virtual ~inverseDistanceInterpolation();


    // Member Functions

        // Access

            //- Return weights
            virtual const scalarFieldFieldList& weights() const;


        //- Update the interpolation
        virtual void update() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    injectionInterpolation

Description
    Injection overset interpolation algorithm
    For each acceptor, the donor data is extracted from the cell where
    the donor centroid falls. Acceptor value is equal to the master donor value,
    while the neighbouring donors do not have any influence at all.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

SourceFiles
    injectionInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef injectionInterpolation_H
#define injectionInterpolation_H

#include "oversetInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class injectionInterpolation Declaration
\*---------------------------------------------------------------------------*/

class injectionInterpolation
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
        injectionInterpolation(const injectionInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const injectionInterpolation&);


        // Functions for calculating/clearing weights

            //- Calculate weights
            void calcWeights() const;

            //- Clear weights
            void clearWeights() const;


public:

    //- Runtime type information
    TypeName("injection");


    // Constructors

        //- Construct from dictionary
        injectionInterpolation
        (
            const oversetMesh& overset,
            const word& name
        );


    //- Destructor
    virtual ~injectionInterpolation();


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

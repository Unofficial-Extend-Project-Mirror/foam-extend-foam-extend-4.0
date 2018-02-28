/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

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
            mutable ListScalarFieldField* weightsPtr_;


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


    // Destructor

        virtual ~inverseDistanceInterpolation();


    // Member Functions

        // Access

            //- Return weights
            virtual const ListScalarFieldField& weights() const;


        //- Update the interpolation
        virtual void update() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

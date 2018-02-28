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
    leastSquaresGradientInterpolation

Description
    Least squares gradient overset interpolation aglorithm
    For each acceptor, the value is extrapolated from master donor cell based on
    least squares evaluation of a gradient using neighbouring donors.

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

SourceFiles
    leastSquaresGradientInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresGradientInterpolation_H
#define leastSquaresGradientInterpolation_H

#include "oversetInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class leastSquaresGradientInterpolation Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresGradientInterpolation
:
    public oversetInterpolation
{
    // Private data

        // Demand-driven data

            // Note: for each master donor (acceptor), we have a list of scalars
            // for all donors: master + neighbouring. Master weight is the first
            // [index 0], followed by neighbouring weights.

            //- Weights for local donor/acceptor pairs (all donors and acceptor
            //  are on the same processor)
            mutable ScalarFieldField* localWeightsPtr_;

            // Master processor data (only accessible on master processor)

                //- Weights for remote donor/acceptor pairs (all donors are on a
                //  remote processor). List for all processors is stored on
                //  master processor (trading off memory for CPU performance)
                mutable ListScalarFieldField* remoteWeightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        leastSquaresGradientInterpolation
        (
            const leastSquaresGradientInterpolation&
        );

        //- Disallow default bitwise assignment
        void operator=(const leastSquaresGradientInterpolation&);


        //- Calculate weights
        void calcWeights() const;

        //- Clear weights
        void clearWeights() const;


public:

    //- Runtime type information
    TypeName("leastSquaresGradient");


    // Constructors

        //- Construct from dictionary
        leastSquaresGradientInterpolation
        (
            const oversetMesh& overset,
            const dictionary& dict
        );


    // Destructor

        virtual ~leastSquaresGradientInterpolation();


    // Member Functions

        //- Return local donor weights
        virtual const ScalarFieldField& localWeights() const;

        //- Return remote donor weights
        virtual const ListScalarFieldField& remoteWeights() const;

        //- Update the interpolation
        virtual void update() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

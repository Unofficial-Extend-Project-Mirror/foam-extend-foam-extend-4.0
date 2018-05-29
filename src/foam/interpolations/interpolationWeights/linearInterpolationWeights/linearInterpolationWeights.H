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
    Foam::linearInterpolationWeights

Description

SourceFiles
    linearInterpolationWeights.C

\*---------------------------------------------------------------------------*/

#ifndef linearInterpolationWeights_H
#define linearInterpolationWeights_H

#include "interpolationWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class linearInterpolationWeights Declaration
\*---------------------------------------------------------------------------*/

class linearInterpolationWeights
:
    public interpolationWeights
{

private:

    // Private data

        //- Cached index in samples from previous invocation
        mutable label index_;

    // Private Member Functions

        //- Get weights of i and i+1 to calculate integration from t to
        //  samples_[i+1]
        Pair<scalar> integrationWeights
        (
            const label i,
            const scalar t
        ) const;

public:

    //- Runtime type information
    TypeName("linear");

    // Constructors

        //- Construct from components
        linearInterpolationWeights
        (
            const scalarField& samples
        );


    //- Destructor
    virtual ~linearInterpolationWeights()
    {}


    // Member Functions

        //- Calculate weights and indices to calculate t from samples.
        //  Returns true if indices changed.
        virtual bool valueWeights
        (
            const scalar t,
            labelList& indices,
            scalarField& weights
        ) const;

        //- Calculate weights and indices to calculate integrand of t1..t2
        //  from samples. Returns true if indices changed.
        virtual bool integrationWeights
        (
            const scalar t1,
            const scalar t2,
            labelList& indices,
            scalarField& weights
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

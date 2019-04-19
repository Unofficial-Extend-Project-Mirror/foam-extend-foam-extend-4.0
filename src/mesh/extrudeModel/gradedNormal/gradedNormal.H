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
    Foam::extrudeModels::gradedNormal

Description
    Extrudes by transforming points normal to the surface. Input parameters are:
    1. Extrusion thickness (in meters),
    2. Initial delta (in meters),
    3. Number of layers.

    Expansion factor is calculated numerically using bisection algorithm.

Author
    Vuko Vukcevic. FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef gradedNormal_H
#define gradedNormal_H

#include "point.H"
#include "extrudeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

/*---------------------------------------------------------------------------*\
                    Class gradedNormal Declaration
\*---------------------------------------------------------------------------*/

class gradedNormal
:
    public extrudeModel
{
    // Private data

        //- Layer thickness
        scalar thickness_;

        //- Initial delta (cell size at the extruded patch)
        scalar delta0_;

        //- Expansion factor
        scalar expansionFactor_;


    // Expansion factor equation (functor class used by BisectionRoot)
    class expansionFactorEqn
    {
        // Private data

            //- Const reference to underlying gradedNormal model
            const gradedNormal& gnm_;


    public:

        //- Construct given gradedNormal model
        explicit expansionFactorEqn(const gradedNormal& gnm)
        :
            gnm_(gnm)
        {}


        //- Function call operator
        scalar operator()(const scalar& x) const
        {
            scalar result = gnm_.thickness();

            for (label i = 0; i <= gnm_.nLayers(); ++i)
            {
                result -= gnm_.delta0()*pow(x, i);
            }

            return result;
        }
    };


public:

    //- Runtime type information
    TypeName("gradedNormal");

    // Constructors

        //- Construct from components
        gradedNormal(const dictionary& dict);


    //- Destructor
    ~gradedNormal();


    // Access functions

        //- Return const reference to thickness
        const scalar& thickness() const
        {
            return thickness_;
        }


        //- Return const reference to initial delta
        const scalar& delta0() const
        {
            return delta0_;
        }


    // Member Operators

        point operator()
        (
            const point& surfacePoint,
            const vector& surfaceNormal,
            const label layer
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


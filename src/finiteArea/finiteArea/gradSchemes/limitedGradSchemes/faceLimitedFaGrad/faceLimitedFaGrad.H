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
    Foam::fa::faceLimitedGrad

Description
    faceLimitedGrad gradient scheme applied to a runTime selected base gradient
    scheme.

    The scalar limiter based on limiting the extrapolated face values
    between the maximum and minumum cell and cell neighbour values and is
    applied to all components of the gradient.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    faceLimitedFaGrads.C

\*---------------------------------------------------------------------------*/

#ifndef faceLimitedFaGrad_H
#define faceLimitedFaGrad_H

#include "faGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                       Class faceLimitedGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faceLimitedGrad
:
    public fa::gradScheme<Type>
{
    // Private Data

        tmp<fa::gradScheme<Type> > basicGradScheme_;

        //- Limiter coefficient
        scalar k_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faceLimitedGrad(const faceLimitedGrad&);

        //- Disallow default bitwise assignment
        void operator=(const faceLimitedGrad&);


public:

    //- RunTime type information
    TypeName("faceLimited");


    // Constructors

        //- Construct from mesh and schemeData
        faceLimitedGrad(const faMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            basicGradScheme_(fa::gradScheme<Type>::New(mesh, schemeData)),
            k_(readScalar(schemeData))
        {
            if (k_ < 0 || k_ > 1)
            {
                FatalIOErrorIn
                (
                    "faceLimitedGrad(const faMesh& mesh, Istream& schemeData)",
                    schemeData
                )   << "coefficient = " << k_
                    << " should be >= 0 and <= 1"
                    << exit(FatalIOError);
            }
        }


    // Member Functions

        static inline void limitEdge
        (
            Type& limiter,
            const Type& maxDelta,
            const Type& minDelta,
            const Type& extrapolate
        );


        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>
        > grad
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

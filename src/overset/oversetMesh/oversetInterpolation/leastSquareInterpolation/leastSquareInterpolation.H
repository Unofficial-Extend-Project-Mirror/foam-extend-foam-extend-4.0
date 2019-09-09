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
    leastSquareInterpolation

Description
    Least squares overset interpolation algorithm

    Donor cells are collected directly

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    leastSquareInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquareInterpolation_H
#define leastSquareInterpolation_H

#include "oversetInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class leastSquareInterpolation Declaration
\*---------------------------------------------------------------------------*/

class leastSquareInterpolation
:
    public oversetInterpolation
{
    // Private data

        // Demand-driven data

            //- Addressing
            mutable labelListList* addressingPtr_;

            //- Weights
            mutable FieldField<Field, scalar>* weightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        leastSquareInterpolation(const leastSquareInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const leastSquareInterpolation&);


        //- Calculate addressing
        void calcAddressing() const;

        //- Clear addressing
        void clearAddressing() const;


public:

    //- Runtime type information
    TypeName("leastSquares");


    // Constructors

        //- Construct from dictionary
        leastSquareInterpolation
        (
            const oversetMesh& overset,
            const dictionary& dict
        );


    //- Destructor
    virtual ~leastSquareInterpolation();


    // Member Functions

            //- Return interpolation local donor list
            virtual const labelList& localDonors() const;

            //- Return interpolation donor list
            virtual const labelList& remoteDonors() const;

            //- Return acceptor addressing
            virtual const labelListList& addressing() const;

            //- Return acceptor weights
            virtual const FieldField<Field, scalar>& weights() const;


        //- Update the interpolation
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

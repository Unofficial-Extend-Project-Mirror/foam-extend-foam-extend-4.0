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


    // Destructor

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

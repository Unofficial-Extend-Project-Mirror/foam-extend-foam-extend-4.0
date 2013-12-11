/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    tractionBoundaryGradient

Description
    Calculates the implicit boundary gradient corresponding to the user
    specified boundary traction.
    This class, which is a function class, can be used in traction boundary
    conditions to calculate the gradient to be applied

SourceFiles
    tractionBoundaryGradient.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef tractionBoundaryGradient_H
#define tractionBoundaryGradient_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "rheologyLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class tractionBoundaryGradient Declaration
\*---------------------------------------------------------------------------*/

class tractionBoundaryGradient
{

public:

    //- Runtime type information
    TypeName("tractionBoundaryGradient");

    // Constructors

        tractionBoundaryGradient()
        {}


    // Destructor

        virtual ~tractionBoundaryGradient()
        {}


    // Member functions

        //- Return the boundary Cauchy traction corresponding to
        //  the given gradient
        tmp<vectorField> traction
        (
            const tensorField& gradField,
            const word fieldName,
            const fvPatch& patch,
            Switch orthotropic,
            word nonLinear
        ) const;


    // Operators

        tmp<vectorField> operator()
        (
            const vectorField& traction,
            const scalarField& pressure,
            const word fieldName,
            const fvPatch& patch,
            Switch orthotropic,
            word nonLinear
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

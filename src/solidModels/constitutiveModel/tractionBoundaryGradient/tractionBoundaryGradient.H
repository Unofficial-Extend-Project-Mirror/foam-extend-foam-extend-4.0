/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Clean-up and re-factoring Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef tractionBoundaryGradient_H
#define tractionBoundaryGradient_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "rheologyLaw.H"
#include "nonLinearGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class tractionBoundaryGradient Declaration
\*---------------------------------------------------------------------------*/

class tractionBoundaryGradient
{
public:

    // Static member functions

        //- Return the boundary Cauchy traction corresponding to
        //  the given gradient
        static tmp<vectorField> traction
        (
            const tensorField& gradField,
            const word& workingFieldName,   // Working variable
            const word& integralFieldName,  // Integrated displacement
            const fvPatch& patch,
            const bool orthotropic,
            const nonLinearGeometry::nonLinearType& nonLinear,
            const bool incremental
        );


    // Operators

        //- Return surface-normal gradient given traction and pressure
        static tmp<vectorField> snGrad
        (
            const vectorField& traction,
            const scalarField& pressure,
            const word& workingFieldName,   // Working variable
            const word& integralFieldName,  // Integrated displacement
            const fvPatch& patch,
            const bool orthotropic,
            const nonLinearGeometry::nonLinearType& nonLinear,
            const bool incremental
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

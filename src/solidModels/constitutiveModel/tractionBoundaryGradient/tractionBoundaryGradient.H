/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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

        // Return the boundary Cauchy traction corresponding to the given gradient
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

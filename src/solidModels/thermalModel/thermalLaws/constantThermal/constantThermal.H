/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008 Hrvoje Jasak
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
    constantThermal

Description
    Constant thermal properties

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    constantThermal.C

\*---------------------------------------------------------------------------*/

#ifndef constantThermal_H
#define constantThermal_H

#include "thermalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantThermal Declaration
\*---------------------------------------------------------------------------*/

class constantThermal
:
    public thermalLaw
{
    // Private data

        //- Specific heat capacity
        dimensionedScalar C_;

        //- Thermal conductivity
        dimensionedScalar k_;

        //- Thermal expansion coefficient
        dimensionedScalar alpha_;

        //- Reference temperature
        dimensionedScalar T0_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantThermal(const constantThermal&);

        //- Disallow default bitwise assignment
        void operator=(const constantThermal&);


public:

    //- Runtime type information
    TypeName("constant");

    // Static data members


    // Constructors

        //- Construct from dictionary
        constantThermal
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~constantThermal();


    // Member Functions

        //- Return specific heat capacity
        virtual tmp<volScalarField> C() const;

        //- Return thermal conductivity
        virtual tmp<volScalarField> k() const;

        //- Return thermal expansion coefficient
        virtual tmp<volScalarField> alpha() const;

        //- Return reference temperature
        virtual tmp<volScalarField> T0() const;

        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

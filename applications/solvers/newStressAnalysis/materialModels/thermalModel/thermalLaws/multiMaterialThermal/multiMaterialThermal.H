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
    multiMaterialThermal

Description
    Zoned multi-material thermal model controlled by an indicator field.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    multiMaterialThermal.C

\*---------------------------------------------------------------------------*/

#ifndef multiMaterialThermal_H
#define multiMaterialThermal_H

#include "thermalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiMaterialThermal Declaration
\*---------------------------------------------------------------------------*/

class multiMaterialThermal
:
    public thermalLaw,
    public PtrList<thermalLaw>
{
    // Private data

        //- Material indicator field
        volScalarField materials_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiMaterialThermal(const multiMaterialThermal&);

        //- Disallow default bitwise assignment
        void operator=(const multiMaterialThermal&);


        //- Calculate indicator field given index
        tmp<scalarField> indicator(const label index) const;

public:

    //- Runtime type information
    TypeName("multiMaterial");

    // Static data members


    // Constructors

        //- Construct from dictionary
        multiMaterialThermal
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~multiMaterialThermal();


    // Member Functions

        //- Return specific heat capacity
        virtual tmp<volScalarField> C() const;

        //- Return thermal conductivity
        virtual tmp<volScalarField> k() const;

        //- Return thermal expansion coefficient
        virtual tmp<volScalarField> alpha() const;

        //- Return reference temperature
        virtual tmp<volScalarField> T0() const;

        //- Correct the rheological model
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

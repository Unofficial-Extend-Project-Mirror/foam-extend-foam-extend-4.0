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
    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiMaterialThermal(const multiMaterialThermal&);

        //- Disallow default bitwise assignment
        void operator=(const multiMaterialThermal&);


        //- Calculate indicator field given index
        tmp<volScalarField> indicator(const label index) const;

protected:

    // Protected data

        //- Material indicator field
        volScalarField materials_;

    // Protected member functions

        //- Check that the laws and the material field are okay
        void readLaws
        (
            const volScalarField& T,
            const dictionary& dict
        );

        //- Check that the laws and the material field are okay
        void checkLaws() const;


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

        //- Construct from dictionary and create default material field
        multiMaterialThermal
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict,
            const scalar defaultMaterial
        );

    // Destructor

        virtual ~multiMaterialThermal();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return specific heat capacity
        virtual tmp<volScalarField> C() const;

        //- Return thermal conductivity
        virtual tmp<volScalarField> k() const;

        //- Correct the rheological model
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

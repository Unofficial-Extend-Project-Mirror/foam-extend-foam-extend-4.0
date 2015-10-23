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
    Foam::solid

Description
    The thermophysical properties of a solid

SourceFiles
    solid.C

\*---------------------------------------------------------------------------*/

#ifndef solid_H
#define solid_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solid;

Ostream& operator<<
(
     Ostream&,
     const solid&
);


/*---------------------------------------------------------------------------*\
                            Class solid Declaration
\*---------------------------------------------------------------------------*/

class solid
{

    // Private data

        //- Density [kg/m3]
        scalar rho_;

        //- Specific heat capacity [J/(kg.K)]
        scalar cp_;

        //- Thermal conductivity [W/(m.K)]
        scalar K_;

        //- Heat of formation [J/kg]
        scalar Hf_;

        //- Emissivity
        scalar emissivity_;


public:

    //- Runtime type information
    TypeName("solid");


    // Declare run-time constructor selection tables

    declareRunTimeSelectionTable
    (
        autoPtr,
        solid,
        ,
        (),
        ()
    );

    declareRunTimeSelectionTable
    (
        autoPtr,
        solid,
        Istream,
        (Istream& is),
        (is)
    );


    // Constructors

        //- Construct from components
        solid
        (
            scalar rho,
            scalar cp,
            scalar K,
            scalar Hf,
            scalar emissivity
        );

        //- Construct from Istream
        solid(Istream& is);

        //- Return a pointer to a new solid created from input
        static autoPtr<solid> New(Istream& is);


    // Destructor

        virtual ~solid()
        {}


    // Member Functions

        // Phisical constants which define the solid

            //- Density [kg/m3]
            inline scalar rho() const;

            //- Specific heat capacity [J/(kg.K)]
            inline scalar cp() const;

            //- Thermal conductivity [W/(m.K)]
            inline scalar K() const;

            //- Heat of formation [J/kg]
            inline scalar Hf() const;

            //- Emissivity []
            inline scalar emissivity() const;


        // I-O

            //- Write the solid properties
            virtual void writeData(Ostream& os) const;


        // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const solid& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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
        )
        :
            rho_(rho),
            cp_(cp),
            K_(K),
            Hf_(Hf),
            emissivity_(emissivity)
        {}

        //- Construct from Istream
        solid(Istream& is)
        :
            rho_(readScalar(is)),
            cp_(readScalar(is)),
            K_(readScalar(is)),
            Hf_(readScalar(is)),
            emissivity_(readScalar(is))
        {}

        //- Return a pointer to a new solid created from input
        static autoPtr<solid> New(Istream& is);


    // Destructor

        virtual ~solid()
        {}


    // Member Functions

        // Phisical constants which define the solid

            //- Density [kg/m3]
            scalar rho() const
            {
                return rho_;
            }

            //- Specific heat capacity [J/(kg.K)]
            scalar cp() const
            {
                return cp_;
            }

            //- Thermal conductivity [W/(m.K)]
            scalar K() const
            {
                return K_;
            }

            //- Heat of formation [J/kg]
            scalar Hf() const
            {
                return Hf_;
            }

            //- emissivity
            scalar emissivity() const
            {
                return emissivity_;
            }


        // I-O

            //- Write the solid properties
            virtual void writeData(Ostream& os) const
            {
                os<< rho_ << token::SPACE
                  << cp_ << token::SPACE
                  << K_ << token::SPACE
                  << Hf_ << token::SPACE
                  << emissivity_;
            }


        // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const solid& s)
        {
            s.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

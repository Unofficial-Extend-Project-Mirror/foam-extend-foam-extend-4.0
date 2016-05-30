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
    Foam::basicThermo

Description
    Basic thermodynamic properties

SourceFiles
    basicThermo.C
    newBasicThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicThermo_H
#define basicThermo_H

#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class basicThermo Declaration
\*---------------------------------------------------------------------------*/

class basicThermo
:
    public IOdictionary
{

protected:

    // Protected data

        // Fields

            //- Pressure [Pa]
            volScalarField p_;

            //- Compressibility [s^2/m^2]
            volScalarField psi_;

            //- Temperature [K]
            volScalarField T_;

            //- Laminar dynamic viscosity [kg/m/s]
            volScalarField mu_;

            //- Laminar thermal diffusuvity [kg/m/s]
            volScalarField alpha_;


    // Protected member functions

        // Enthalpy

            //- Return the enthalpy field boundary types by interrogating the
            //  temperature field boundary types
            wordList hBoundaryTypes();

            //- Correct the enthalpy field boundaries
            void hBoundaryCorrection(volScalarField& h);


        // Internal energy

            //- Return the internal energy  field boundary types by
            //  interrogating the temperature field boundary types
            wordList eBoundaryTypes();

            //- Correct the internal energy field boundaries
            void eBoundaryCorrection(volScalarField& e);


        //- Construct as copy (not implemented)
        basicThermo(const basicThermo&);


public:

    //- Runtime type information
    TypeName("basicThermo");


    // Constructors

        //- Construct from mesh and object registry
        basicThermo(const fvMesh&, const objectRegistry&);


    //- Destructor
    virtual ~basicThermo();


    // Member functions

        //- Update properties
        virtual void correct() = 0;


        // Access to thermodynamic state variables

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p();

            //- Pressure [Pa]
            virtual const volScalarField& p() const;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const;

            //- Total enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& h();

            //- Total enthalpy [J/kg]
            virtual const volScalarField& h() const;

            //- Total enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Total enthalpy for patch [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Sensible enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& hs();

            //- Sensible enthalpy [J/kg]
            virtual const volScalarField& hs() const;

            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Sensible enthalpy for patch [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const;

            //- Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& e();

            //- Internal energy [J/kg]
            virtual const volScalarField& e() const;

            //- Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> e
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Internal energy for patch [J/kg]
            virtual tmp<scalarField> e
            (
                const scalarField& T,
                const label patchi
            ) const;


        // Fields derived from thermodynamic state variables

            //- Temperature [K]
            virtual const volScalarField& T() const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for cell-set [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual const volScalarField& mu() const;

            //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
            virtual const volScalarField& alpha() const;


        //- Read thermophysicalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

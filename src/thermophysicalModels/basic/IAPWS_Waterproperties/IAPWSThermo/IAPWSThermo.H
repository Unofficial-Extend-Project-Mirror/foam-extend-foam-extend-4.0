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
    Foam::IAPWSThermo

Description:

    Waterproperties based on the IAPWS 97 tables
    The water properties are caluclated using freeSteam (http://freesteam.sourceforge.net/example.php)

    General paper decribing the water tables:

    "Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam"

SourceFiles
    IAPWSThermo.C

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig
Germany

\*---------------------------------------------------------------------------*/

#ifndef IAPWSThermo_H
#define IAPWSThermo_H

#include "basicPsiThermo.H"
#include "IAPWS-IF97.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class IAPWSTHERMO Declaration
\*---------------------------------------------------------------------------*/

class IAPWSThermo
:
    public basicPsiThermo
{
    // Private data

        //- Enthalpy field
        volScalarField h_;

        //- DensityField
        volScalarField rho_;

        // CL:drhodh Field
        // CL:needed for pressure equation
        volScalarField drhodh_;


    // Private member functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        IAPWSThermo(const IAPWSThermo&);


public:

    //- Runtime type information
    TypeName("IAPWSThermo");

    // Constructors

        //- Construct from mesh
        IAPWSThermo(const fvMesh&, const objectRegistry& obj);


    //- Destructor
    virtual ~IAPWSThermo();

    // Member functions



        //- Update properties
        virtual void correct();


        // Access to thermodynamic state variables

            //- Enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& h()
            {
                return h_;
            }

            //- Enthalpy [J/kg]
            virtual const volScalarField& h() const
            {
                return h_;
            }

        // Fields derived from thermodynamic state variables

            //- Enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const labelList& cells
            ) const;


            //- Enthalpy for patch [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const label patchi
            ) const;


            //- Density for patch [J/kg]
            virtual tmp<scalarField> rho
            (
                const scalarField& p,
                const scalarField& h,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            // dummy function needed for BC
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
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

            //- Gradient drhodh @ constant pressure
            virtual const volScalarField& drhodh() const
            {
                return drhodh_;
            }

             //- Density [kg/m^3] - uses current value of pressure
             virtual tmp<volScalarField> rho() const;


        //- Read thermophysicalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IAPWSThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

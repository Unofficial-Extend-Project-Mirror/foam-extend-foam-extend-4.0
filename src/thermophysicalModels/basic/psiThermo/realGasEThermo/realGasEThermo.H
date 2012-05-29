/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | .
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
    Foam::realGasEThermo

Description
    Internal energy for a real gas fluid libary


SourceFiles
    realGasEThermo.C

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany


\*---------------------------------------------------------------------------*/

#ifndef realGasEThermo_H
#define realGasEThermo_H

#include "basicPsiThermo.H"
#include "basicMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class realGasEThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class realGasEThermo
:
    public basicPsiThermo,
    public MixtureType
{
    // Private data

        //- Enthalpy field
        volScalarField e_;
	
        //- DensityField
        volScalarField rho_;

        //- drhode_Field
        volScalarField drhode_;

    // Private member functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        realGasEThermo(const realGasEThermo<MixtureType>&);


public:

    //- Runtime type information
    TypeName("realGasEThermo");


    // Constructors

        //- Construct from mesh
        realGasEThermo(const fvMesh&);

    //- Destructor
    virtual ~realGasEThermo();


    // Member functions

        //- Return the compostion of the mixture
        virtual basicMixture& composition()
        {
            return *this;
        }

        //- Return the compostion of the mixture
        virtual const basicMixture& composition() const
        {
            return *this;
        }

        //- Update properties
        virtual void correct();


        // Access to thermodynamic state variables

            //- Enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& e()
            {
                return e_;
            }

            //- Enthalpy [J/kg]
            virtual const volScalarField& e() const
            {
                return e_;
            }

            //CL: drhode needed for pressure equation of the real gas solver 
            virtual const volScalarField& drhode() const
            {
                return drhode_;
            }

        // Fields derived from thermodynamic state variables

            //- Enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> e
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy for patch [J/kg]
            virtual tmp<scalarField> e
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Density for patch [J/kg]
            virtual tmp<scalarField> rho
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
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

            //- Density [kg/m^3] - uses current value of pressure
            virtual tmp<volScalarField> rho() const;


        //- Read thermophysicalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "realGasEThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

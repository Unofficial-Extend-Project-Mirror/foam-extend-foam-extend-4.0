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
#include "runTimeSelectionTables.H"

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

        volScalarField p_;
        volScalarField T_;
        volScalarField psi_;

        volScalarField mu_;
        volScalarField alpha_;


    // Protected member functions

        wordList hBoundaryTypes();
        void hBoundaryCorrection(volScalarField& h);

        //- Construct as copy (not implemented)
        basicThermo(const basicThermo&);


public:

    //- Runtime type information
    TypeName("basicThermo");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            basicThermo,
            fvMesh,
            (const fvMesh& mesh),
            (mesh)
        );


    // Constructors

        //- Construct from mesh
        basicThermo(const fvMesh&);


    // Selectors

        static autoPtr<basicThermo> New(const fvMesh&);


    // Destructor

        virtual ~basicThermo();


    // Member functions

        //- Update properties
        virtual void correct() = 0;


        // Access to thermodynamic state variables.

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            volScalarField& p()
            {
                return p_;
            }

            //- Pressure [Pa]
            const volScalarField& p() const
            {
                return p_;
            }

            //- Enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& h()
            {
                notImplemented("basicThermo::h()");
                return const_cast<volScalarField&>(volScalarField::null());
            }

            //- Enthalpy [J/kg]
            virtual const volScalarField& h() const
            {
                notImplemented("basicThermo::h() const");
                return volScalarField::null();
            }

            //- Enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const labelList& cells
            ) const
            {
                notImplemented
                (
                    "basicThermo::h"
                    "(const scalarField& T, const labelList& cells) const"
                );
                return tmp<scalarField>(NULL);
            }

            //- Enthalpy for patch [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const label patchi
            ) const
            {
                notImplemented
                (
                    "basicThermo::h"
                    "(const scalarField& T, const label patchi) const"
                );
                return tmp<scalarField>(NULL);
            }

            //- Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& e()
            {
                notImplemented("basicThermo::e()");
                return const_cast<volScalarField&>(volScalarField::null());
            }

            //- Internal energy [J/kg]
            virtual const volScalarField& e() const
            {
                notImplemented("basicThermo::e()");
                return volScalarField::null();
            }


        // Fields derived from thermodynamic state variables

            //- Temperature [K]
            const volScalarField& T() const
            {
                return T_;
            }

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const
            {
                return p_*psi();
            }

            //- Compressibility [s^2/m^2]
            const volScalarField& psi() const
            {
                return psi_;
            }

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const
            {
                notImplemented
                (
                    "basicThermo::Cp"
                    "(const scalarField& T, const label patchi) const"
                );
                return tmp<scalarField>(NULL);
            }

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const
            {
                notImplemented("basicThermo::Cp() const");
                return volScalarField::null();
            }

            //- Heat capacity at constant voluve for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const
            {
                notImplemented
                (
                    "basicThermo::Cv"
                    "(const scalarField& T, const label patchi) const"
                );
                return tmp<scalarField>(NULL);
            }

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const
            {
                notImplemented("basicThermo::Cv() const");
                return volScalarField::null();
            }


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/ms]
            const volScalarField& mu() const
            {
                return mu_;
            }

            //- Thermal diffusivity for enthalpy of mixture [kg/ms]
            const volScalarField& alpha() const
            {
                return alpha_;
            }


        //- Read thermophysicalProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

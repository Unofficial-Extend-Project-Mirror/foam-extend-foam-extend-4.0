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
    Foam::hhuCombustionThermo

Description
    Foam::hhuCombustionThermo

SourceFiles
    hhuCombustionThermo.C
    newhhuCombustionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hhuCombustionThermo_H
#define hhuCombustionThermo_H

#include "hCombustionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class hhuCombustionThermo Declaration
\*---------------------------------------------------------------------------*/

class hhuCombustionThermo
:
    public hCombustionThermo
{

protected:

    // Protected data

        volScalarField Tu_;
        volScalarField hu_;


    // Protected member functions

        wordList huBoundaryTypes();
        void huBoundaryCorrection(volScalarField& hu);


public:

    //- Runtime type information
    TypeName("hhuCombustionThermo");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            hhuCombustionThermo,
            fvMesh,
            (const fvMesh& mesh, const objectRegistry& obj),
            (mesh, obj)
        );


    // Constructors

        //- Construct from dictionary and mesh
        hhuCombustionThermo(const fvMesh&, const objectRegistry&);


    // Selectors

        static autoPtr<hhuCombustionThermo> New(const fvMesh&, const objectRegistry&);

        static autoPtr<hhuCombustionThermo> New(const fvMesh& mesh)
        {
            return New(mesh, mesh);
        }


    //- Destructor
    virtual ~hhuCombustionThermo();


    // Member functions

        //- Update properties
        virtual void correct() = 0;


        // Access to thermodynamic state variables.

            //- Unburnt gas enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& hu()
            {
                return hu_;
            }

            //- Unburnt gas enthalpy [J/kg]
            virtual const volScalarField& hu() const
            {
                return hu_;
            }


        // Fields derived from thermodynamic state variables

            //- Unburnt gas enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hu
            (
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Unburnt gas enthalpy for patch [J/kg]
            virtual tmp<scalarField> hu
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Unburnt gas temperature [K]
            virtual const volScalarField& Tu() const
            {
                return Tu_;
            }

            //- Burnt gas temperature [K]
            virtual tmp<volScalarField> Tb() const = 0;

            //- Unburnt gas density [kg/m^3]
            virtual tmp<volScalarField> rhou() const
            {
                return p_*psiu();
            }

            //- Burnt gas density [kg/m^3]
            virtual tmp<volScalarField> rhob() const
            {
                return p_*psib();
            }

            //- Unburnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psiu() const = 0;

            //- Burnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psib() const = 0;

            //- Dynamic viscosity of unburnt gas [kg/ms]
            virtual tmp<volScalarField> muu() const = 0;

            //- Dynamic viscosity of burnt gas [kg/ms]
            virtual tmp<volScalarField> mub() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

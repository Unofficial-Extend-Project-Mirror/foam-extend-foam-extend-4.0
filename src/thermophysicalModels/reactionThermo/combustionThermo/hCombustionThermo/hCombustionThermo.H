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
    Foam::hCombustionThermo

Description
    Foam::hCombustionThermo

SourceFiles
    hCombustionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hCombustionThermo_H
#define hCombustionThermo_H

#include "basicPsiThermo.H"
#include "basicMultiComponentMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class hCombustionThermo Declaration
\*---------------------------------------------------------------------------*/

class hCombustionThermo
:
    public basicPsiThermo
{

protected:

    // Protected data

        //- Enthalpy field
        volScalarField h_;


public:

    //- Runtime type information
    TypeName("hCombustionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        hCombustionThermo,
        fvMesh,
        (const fvMesh& mesh, const objectRegistry& obj),
        (mesh, obj)
    );


    // Constructors

        //- Construct from dictionary and mesh
        hCombustionThermo(const fvMesh&, const objectRegistry&);


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<hCombustionThermo> New(const fvMesh&, const objectRegistry&);

        static autoPtr<hCombustionThermo> New(const fvMesh& mesh)
        {
            return New(mesh, mesh);
        }

        //- Select and check that package contains 'thermoType'
        static autoPtr<hCombustionThermo> NewType
        (
            const fvMesh&,
            const objectRegistry&,
            const word& thermoType
        );


    //- Destructor
    virtual ~hCombustionThermo();


    // Member functions

        //- Return the composition of the multi-component mixture
        virtual basicMultiComponentMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicMultiComponentMixture& composition() const = 0;


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


        //- Update properties
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

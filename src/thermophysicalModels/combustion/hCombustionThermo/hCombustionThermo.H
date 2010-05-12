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
    Foam::hCombustionThermo

Description
    Foam::hCombustionThermo

SourceFiles
    hCombustionThermo.C
    newhCombustionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hCombustionThermo_H
#define hCombustionThermo_H

#include "basicThermo.H"
#include "combustionMixture.H"
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
    public basicThermo
{

protected:

    // Protected data

        volScalarField h_;


public:

    //- Runtime type information
    TypeName("hCombustionThermo");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            hCombustionThermo,
            fvMesh,
            (const fvMesh& mesh),
            (mesh)
        );


    // Constructors

        //- Construct from dictionary and mesh
        hCombustionThermo(const fvMesh&);


    // Selectors

        static autoPtr<hCombustionThermo> New(const fvMesh&);


    // Destructor

        virtual ~hCombustionThermo();


    // Member functions

        //- Return the composition of the combustion mixture
        virtual combustionMixture& composition() = 0;

        //- Return the composition of the combustion mixture
        virtual const combustionMixture& composition() const = 0;


        // Access to thermodynamic state variables

            //- Enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            volScalarField& h()
            {
                return h_;
            }

            //- Enthalpy [J/kg]
            const volScalarField& h() const
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

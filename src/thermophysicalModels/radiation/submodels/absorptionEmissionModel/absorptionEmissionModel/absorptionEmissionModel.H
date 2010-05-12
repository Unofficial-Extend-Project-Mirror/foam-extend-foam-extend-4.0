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
    Foam::radiation::absorptionEmissionModel

Description
    Model to supply absorption and emission coefficients for radiation
    modelling

\*---------------------------------------------------------------------------*/

#ifndef radiationAbsorptionEmissionModel_H
#define radiationAbsorptionEmissionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                   Class absorptionEmissionModel Declaration
\*---------------------------------------------------------------------------*/

class absorptionEmissionModel
{

protected:

    // Protected data

        //- Radiation model dictionary
        const dictionary dict_;

        //- Reference to the fvMesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("absorptionEmissionModel");

    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            absorptionEmissionModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Constructors

        //- Construct from components
        absorptionEmissionModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Selector

        static autoPtr<absorptionEmissionModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor

        virtual ~absorptionEmissionModel();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient (net)
                virtual tmp<volScalarField> a() const;

                //- Absorption coefficient for continuous phase
                virtual tmp<volScalarField> aCont() const;

                //- Absorption coefficient for dispersed phase
                virtual tmp<volScalarField> aDisp() const;


            // Emission coefficient

                //- Emission coefficient (net)
                virtual tmp<volScalarField> e() const;

                //- Return emission coefficient for continuous phase
                virtual tmp<volScalarField> eCont() const;

                //- Return emission coefficient for dispersed phase
                virtual tmp<volScalarField> eDisp() const;


            // Emission contribution

                //- Emission contribution (net)
                virtual tmp<volScalarField> E() const;

                //- Emission contribution for continuous phase
                virtual tmp<volScalarField> ECont() const;

                //- Emission contribution for dispersed phase
                virtual tmp<volScalarField> EDisp() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


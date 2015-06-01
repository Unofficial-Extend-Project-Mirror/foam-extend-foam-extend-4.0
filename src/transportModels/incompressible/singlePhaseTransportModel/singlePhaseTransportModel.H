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
    Foam::singlePhaseTransportModel

Description
    A simple single-phase transport model based on viscosityModel.

    Used by the incompressible single-phase solvers like simpleFoam,
    turbFoam etc.

SourceFiles
    singlePhaseTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseTransportModel_H
#define singlePhaseTransportModel_H

#include "incompressible/transportModel/transportModel.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class viscosityModel;

/*---------------------------------------------------------------------------*\
                Class singlePhaseTransportModel Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseTransportModel
:
    public transportModel
{
    // Private Data

        autoPtr<viscosityModel> viscosityModelPtr_;


    // Private Member Functions

        //- Disallow copy construct
        singlePhaseTransportModel(const singlePhaseTransportModel&);

        //- Disallow default bitwise assignment
        void operator=(const singlePhaseTransportModel&);


public:

    // Constructors

        //- Construct from components
        singlePhaseTransportModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~singlePhaseTransportModel();


    // Member Functions

        //- Return the laminar viscosity
        virtual const volScalarField& nu() const;

        //- Correct the laminar viscosity
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

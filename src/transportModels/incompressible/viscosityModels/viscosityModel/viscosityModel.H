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

Namespace
    Foam::viscosityModels

Description
    A namespace for various incompressible viscosityModel implementations.

Class
    Foam::viscosityModel

Description
    An abstract base class for incompressible viscosityModels

SourceFiles
    viscosityModel.C
    newViscosityModel.C

\*---------------------------------------------------------------------------*/

#ifndef viscosityModel_H
#define viscosityModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class viscosityModel Declaration
\*---------------------------------------------------------------------------*/

class viscosityModel
{

protected:

    // Protected data

        word name_;
        dictionary viscosityProperties_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;


    // Private Member Functions

        //- Disallow copy construct
        viscosityModel(const viscosityModel&);

        //- Disallow default bitwise assignment
        void operator=(const viscosityModel&);


public:

    //- Runtime type information
    TypeName("viscosityModel");


    // Declare run-time constructor selection table

#ifndef SWIG
        declareRunTimeSelectionTable
        (
            autoPtr,
            viscosityModel,
            dictionary,
            (
                const word& name,
                const dictionary& viscosityProperties,
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (name, viscosityProperties, U, phi)
        );
#endif


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<viscosityModel> New
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        viscosityModel
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~viscosityModel()
        {}


    // Member Functions

        //- Return the phase transport properties dictionary
        const dictionary& viscosityProperties() const
        {
            return viscosityProperties_;
        }

        //- Return the strain rate mag(grad(U))
        tmp<volScalarField> strainRate() const;

        //- Return the laminar viscosity
        virtual const volScalarField& nu() const = 0;

        //- Correct the laminar viscosity
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

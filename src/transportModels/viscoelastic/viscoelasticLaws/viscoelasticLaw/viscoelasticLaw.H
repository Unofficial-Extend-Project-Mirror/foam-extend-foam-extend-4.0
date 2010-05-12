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
    viscoelasticLaw

Description
    Abstract base class for viscoelastic flow models.

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    viscoelasticLaw.C
    newViscoelasticLaw.C

\*---------------------------------------------------------------------------*/

#ifndef viscoelasticLaw_H
#define viscoelasticLaw_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class viscoelasticLaw Declaration
\*---------------------------------------------------------------------------*/

class viscoelasticLaw
{
private:

    // Private data

        //- Name
        const word name_;

        //- Reference to velocity field
        const volVectorField& U_;

        //- Reference to face flux field
        const surfaceScalarField& phi_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        viscoelasticLaw(const viscoelasticLaw&);

        //- Disallow default bitwise assignment
        void operator=(const viscoelasticLaw&);


public:

    //- Runtime type information
    TypeName("viscoelasticLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            viscoelasticLaw,
            dictionary,
            (
                const word& name,
                const volVectorField& U,
                const surfaceScalarField& phi,
                const dictionary& dict
            ),
            (name, U, phi, dict)
        );


    // Constructors

        //- Construct from components
        viscoelasticLaw
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Selectors

        //- Return a reference to the selected viscoelastic law
        static autoPtr<viscoelasticLaw> New
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~viscoelasticLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return the velocity field
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the face flux field
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const = 0;

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const = 0;

        //- Correct the viscoelastic stress
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

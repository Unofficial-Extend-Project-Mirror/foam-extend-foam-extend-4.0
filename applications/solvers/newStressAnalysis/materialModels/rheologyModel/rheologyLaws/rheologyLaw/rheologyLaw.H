/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    rheologyLaw

Description
    Material rheology for solids.

SourceFiles
    rheologyLaw.C
    newRheologyLaw.C

\*---------------------------------------------------------------------------*/

#ifndef rheologyLaw_H
#define rheologyLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rheologyLaw Declaration
\*---------------------------------------------------------------------------*/

class rheologyLaw
{
    // Private data

        //- Name
        const word name_;

        //- Reference to stress field
        const volSymmTensorField& sigma_;


    // Private Member Functions

        //- Disallow copy construct
        rheologyLaw(const rheologyLaw&);

        //- Disallow default bitwise assignment
        void operator=(const rheologyLaw&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return sigma_.mesh();
        }

        //- Return reference to stress field
        const volSymmTensorField& sigma() const
        {
            return sigma_;
        }


public:

    //- Runtime type information
    TypeName("rheologyLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            rheologyLaw,
            dictionary,
            (
                const word name,
                const volSymmTensorField& sigma,
                const dictionary& dict
            ),
            (name, sigma, dict)
        );


    // Selectors

        //- Return a reference to the selected rheology model
        static autoPtr<rheologyLaw> New
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        rheologyLaw
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Destructor

        virtual ~rheologyLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }


        //- Return density
        virtual tmp<volScalarField> rho() const = 0;

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E() const = 0;

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu() const = 0;

        //- Return density
        virtual tmp<volScalarField> rho(scalar t) const = 0;

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E(scalar t) const = 0;

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu(scalar t) const = 0;

        //- Return creep compliance
        virtual tmp<volScalarField> J(scalar t) const = 0;

        //- Correct the rheological model
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

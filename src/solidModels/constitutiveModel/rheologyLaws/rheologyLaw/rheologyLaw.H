/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E(const volScalarField& epsEq) const = 0;

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu(scalar t) const = 0;

        //- Return creep compliance
        virtual tmp<volScalarField> J(scalar t) const = 0;

        //- Return yield stress
        virtual tmp<volScalarField> sigmaY() const = 0;

        //- Return yield stress for given cellID
        // cellID needed for multiple materials
        virtual scalar sigmaY
        (const scalar epsilonPEq, const label cellID) const = 0;

        //- If the rheologyLaw needs plasticity procedures
        virtual bool plasticityModelNeeded() const
        {
            return false;
        };

        //- If effects are active
        virtual bool viscoActive() const
        {
            return false;
        };

        //- Return plastic modulus
        virtual tmp<volScalarField> Ep() const = 0;

        //- Return plastic modulus
        virtual tmp<volScalarField> Ep(const volScalarField& epsEq) const = 0;

        //- orthotropic tensorial diffusivity
        virtual tmp<volDiagTensorField> K() const;

        //- fourth order elastic stiffness tensor
        virtual tmp<volSymmTensor4thOrderField> C() const;

        //- Correct the rheological model
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    frictionLaw

Description
    friction law e.g. Coulomb, Orowan, Shaw, etc.
    Given the pressure (or other variable for more complex laws),
    the friction law returns the shear traction.

SourceFiles
    frictionLaw.C
    newFrictionLaw.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef frictionLaw_H
#define frictionLaw_H

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
                         Class frictionLaw Declaration
\*---------------------------------------------------------------------------*/

class frictionLaw
{
    // Private data

        // friction law dict
        dictionary frictionLawDict_;

    // Private Member Functions

        //- Disallow copy construct
        frictionLaw(const frictionLaw&);

        //- Disallow default bitwise assignment
        void operator=(const frictionLaw&);


public:

    //- Runtime type information
    TypeName("frictionLaw");

    // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            frictionLaw,
            dictionary,
            (
                const word name,
                const dictionary& dict
            ),
            (name, dict)
        );


    // Selectors

        //- Return a reference to the selected rheology model
        static autoPtr<frictionLaw> New
        (
            const word& name,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        frictionLaw
        (
            const word& name,
            const dictionary& dict
        );


    // Destructor

        virtual ~frictionLaw()
        {}


    // Member Functions

        //- Return friction law dict
        const dictionary& frictionLawDict() const
        {
            return frictionLawDict_;
        }

        //- Return slip traction field
        // the traction that would cause slipping
        // for the given pressure
        virtual scalar slipTraction(const scalar pressure) = 0;

        //- Return slip traction field
        // the traction that would cause slipping
        // for the given pressure
        virtual tmp<scalarField> slipTraction(const scalarField& pressure);

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::rotationalConstraint

Description
    Abstract base class containing interface for rotation constraints used
    within sixDOFODE classes.

    The constraint is implicitly defined as:

    g(omegaDot, t) = f(t)*omegaDot + a(t) = 0,

    where omegaDot is angular acceleration.

    Interface provides all the necessary data for inserting the constraint into
    the resulting linear system via Lagrangian multipliers:
    1. matrixContribution() - corresponding to f(t) (prefactor multiplying
                              omegaDot), to be inserted into the matrix.
    2. sourceContribution() - corresponding to a(t), to be inserted into right
                              hand side vector.

Author
    Viktor Pandza, FSB Zagreb.  All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    rotationalConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef rotationalConstraint_H
#define rotationalConstraint_H

#include "sixDOFODE.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class rotationalConstraint Declaration
\*---------------------------------------------------------------------------*/

class rotationalConstraint
{
    // Private data

        //- Name of the constraint
        const word name_;

        //- Reference to underlying sixDOFODE
        const sixDOFODE& sixDOFODE_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        rotationalConstraint(const rotationalConstraint&);

        //- Disallow default bitwise assignment
        void operator=(const rotationalConstraint&);


public:

    //- Runtime type information
    TypeName("rotationalConstraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            rotationalConstraint,
            word,
            (
                const word& name,
                const sixDOFODE& sixDOFODE,
                const dictionary& dict
            ),
            (name, sixDOFODE, dict)
        );


        //- Class used for the read-construction of
        //  PtrLists of rotationalConstraint
        class iNew
        {
            const sixDOFODE& sixDOFODE_;

        public:

            iNew(const sixDOFODE& sixDOFODE)
            :
                sixDOFODE_(sixDOFODE)
            {}

            autoPtr<rotationalConstraint> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);
                return rotationalConstraint::New(name, sixDOFODE_, dict);
            }
        };


    // Constructors

        //- Construct from dictionary
        rotationalConstraint
        (
            const word& name,
            const sixDOFODE& sixDOFODE,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected rotationalConstraint
        static autoPtr<rotationalConstraint> New
        (
            const word& name,
            const sixDOFODE& sixDOFODE,
            const dictionary& dict
        );


    // Destructor

        virtual ~rotationalConstraint();


    // Member Functions

        // Access functions

            //- Return const reference to name of the constraint
            const word& name() const
            {
                return name_;
            }

            //- Return const reference to underlying sixDOFODE object
            const sixDOFODE& sixDOF() const
            {
                return sixDOFODE_;
            }


        // Constraint specific functions

            //- Return matrix contribution defined by constraint, f(t)
            virtual vector matrixContribution
            (
                const scalar t,
                const tensor& toRelative
            ) const = 0;

            //- Return source contribution defined by constraint, a(t)
            virtual scalar sourceContribution(const scalar t) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

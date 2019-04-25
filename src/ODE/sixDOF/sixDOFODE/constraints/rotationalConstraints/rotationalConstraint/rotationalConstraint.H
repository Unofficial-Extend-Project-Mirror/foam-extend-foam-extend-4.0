/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

    Notes:
    1. Constraints are usually used alongside appropriate initial
       conditions (rotational rate in a given direction),
    2. According to DAE theory, a stabilisation on orientation/angular velocity
       is necessary based on "differentiation index" used to obtain the final
       form of the constraint.

Author
    Viktor Pandza, FSB Zagreb.  All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    rotationalConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef rotationalConstraint_H
#define rotationalConstraint_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "dimensionedTypes.H"
#include "dictionary.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class sixDOFODE;

// Forward declaration of friend functions and operators
class rotationalConstraint;

Ostream& operator<<(Ostream&, const rotationalConstraint&);


/*---------------------------------------------------------------------------*\
                   Class rotationalConstraint Declaration
\*---------------------------------------------------------------------------*/

class rotationalConstraint
{
    // Private Data

        //- Name of the constraint
        word name_;

        //- Reference to underlying sixDOFODE
        const sixDOFODE& sixDOF_;

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
                const dictionary& dict,
                const sixDOFODE& sixDOF
            ),
            (name, dict, sixDOF)
        );


        //- Class used for the read-construction of
        //  PtrLists of rotationalConstraint
        class iNew
        {
            const sixDOFODE& sixDOF_;

        public:

            iNew(const sixDOFODE& sixDOF)
            :
                sixDOF_(sixDOF)
            {}

            autoPtr<rotationalConstraint> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);
                return rotationalConstraint::New(name, dict, sixDOF_);
            }
        };


    // Constructors

        //- Construct from dictionary
        rotationalConstraint
        (
            const word& name,
            const dictionary& dict,
            const sixDOFODE& sixDOF
        );

        //- Construct and return a clone
        virtual autoPtr<rotationalConstraint> clone() const = 0;


    // Selectors

        //- Return a reference to the selected rotationalConstraint
        static autoPtr<rotationalConstraint> New
        (
            const word& name,
            const dictionary& dict,
            const sixDOFODE& sixDOF
        );


    // Destructor

        virtual ~rotationalConstraint();


    // Member Functions

        // Access functions

            //- Return name
            const word& name() const;

            //- Return underlying sixDOFODE object
            const sixDOFODE& sixDOF() const;


        // Constraint specific functions

            //- Return matrix contribution defined by constraint, f(t)
            virtual vector matrixContribution
            (
                const scalar t,
                const tensor& toRelative,
                const vector& omega
            ) const = 0;

            //- Return source contribution defined by constraint, a(t)
            virtual scalar sourceContribution
            (
                const scalar t,
                const tensor& toRelative,
                const vector& omega
            ) const = 0;

            //- Stabilise the constraint (necessary for constraints with
            //  differentiation index higher than zero)
            virtual void stabilise
            (
                const scalar t,
                vector& omega
            ) const = 0;


    // I-O Functions and Operators

        //- Virtual write function
        virtual void write(Ostream& os) const = 0;

        //- Ostream operator implemented in terms of write operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const rotationalConstraint& rc
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

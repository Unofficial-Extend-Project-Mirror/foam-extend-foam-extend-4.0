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
    Foam::rotationalRestraint

Description
    Abstract base class containing interface for rotational restraints used
    within sixDOFODE classes.

    Interface provides restraining moment to the rotational equations of
    motion via moment() member function.

Author
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    rotationalRestraint.C

\*---------------------------------------------------------------------------*/

#ifndef rotationalRestraint_H
#define rotationalRestraint_H

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
class rotationalRestraint;

Ostream& operator<<(Ostream&, const rotationalRestraint&);


/*---------------------------------------------------------------------------*\
                   Class rotationalRestraint Declaration
\*---------------------------------------------------------------------------*/

class rotationalRestraint
{
    // Private Data

        //- Name of the restraint
        word name_;

        //- Reference to underlying sixDOFODE
        const sixDOFODE& sixDOF_;

public:

    //- Runtime type information
    TypeName("rotationalRestraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            rotationalRestraint,
            word,
            (
                const word& name,
                const dictionary& dict,
                const sixDOFODE& sixDOF
            ),
            (name, dict, sixDOF)
        );


        //- Class used for the read-construction of
        //  PtrLists of rotationalRestraint
        class iNew
        {
            const sixDOFODE& sixDOF_;

        public:

            iNew(const sixDOFODE& sixDOF)
            :
                sixDOF_(sixDOF)
            {}

            autoPtr<rotationalRestraint> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);
                return rotationalRestraint::New(name, dict, sixDOF_);
            }
        };


    // Constructors

        //- Construct from dictionary
        rotationalRestraint
        (
            const word& name,
            const dictionary& dict,
            const sixDOFODE& sixDOF
        );

        //- Construct and return a clone
        virtual autoPtr<rotationalRestraint> clone() const = 0;


    // Selectors

        //- Return a reference to the selected rotationalRestraint
        static autoPtr<rotationalRestraint> New
        (
            const word& name,
            const dictionary& dict,
            const sixDOFODE& sixDOF
        );


    // Destructor

        virtual ~rotationalRestraint();


    // Member Functions

        // Access functions

            //- Return name
            const word& name() const;

            //- Return underlying sixDOFODE object
            const sixDOFODE& sixDOF() const;


        // Restraint specific functions

            //- Return restraining moment (in the local coordinate system)
            virtual vector restrainingMoment
            (
                const scalar t,
                const tensor& toRelative,
                const vector& omega
            ) const = 0;


    // I-O Functions and Operators

        //- Virtual write function
        virtual void write(Ostream& os) const = 0;

        //- Ostream operator implemented in terms of write operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const rotationalRestraint& rr
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

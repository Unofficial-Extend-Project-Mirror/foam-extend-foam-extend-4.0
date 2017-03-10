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
    Foam::translationalRestraint

Description
    Abstract base class containing interface for translational restraints used
    within sixDOFODE classes.

    Interface provides restraining forces to the translational equations of
    motion via force() member function.

Author
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    translationalRestraint.C

\*---------------------------------------------------------------------------*/

#ifndef translationalRestraint_H
#define translationalRestraint_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "dimensionedTypes.H"
#include "dictionary.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class translationalRestraint Declaration
\*---------------------------------------------------------------------------*/

class translationalRestraint
{
    // Private Data

        //- Name of the restraint
        word name_;


public:

    //- Runtime type information
    TypeName("translationalRestraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            translationalRestraint,
            word,
            (
                const word& name,
                const dictionary& dict
            ),
            (name, dict)
        );


        //- Class used for the read-construction of
        //  PtrLists of translationalRestraint
        class iNew
        {
        public:

            iNew()
            {}

            autoPtr<translationalRestraint> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);
                return translationalRestraint::New(name, dict);
            }
        };


    // Constructors

        //- Construct from dictionary
        translationalRestraint
        (
            const word& name,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<translationalRestraint> clone() const = 0;


    // Selectors

        //- Return a reference to the selected translationalRestraint
        static autoPtr<translationalRestraint> New
        (
            const word& name,
            const dictionary& dict
        );


    // Destructor

        virtual ~translationalRestraint();


    // Member Functions

        // Restraint specific functions

            //- Return restraining force (in the global coordinate system)
            virtual vector restrainingForce
            (
                const scalar t,
                const tensor& toRelative,
                const vector& x,
                const vector& u
            ) const = 0;


    // I-O Functions and Operators

        //- Virtual write function
        virtual void write(Ostream& os) const = 0;

        //- Ostream operator implemented in terms of write operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const translationalRestraint& tr
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

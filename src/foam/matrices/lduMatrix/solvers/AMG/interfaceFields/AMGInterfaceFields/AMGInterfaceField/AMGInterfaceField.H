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
    Foam::AMGInterfaceField

Description
     Abstract base class for AMG agglomerated interface fields.

SourceFiles
    AMGInterfaceField.C
    newAmgInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef AMGInterfaceField_H
#define AMGInterfaceField_H

#include "lduInterfaceField.H"
#include "AMGInterface.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class AMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class AMGInterfaceField
:
    public lduInterfaceField
{
    // Private data

        //- Local reference cast into the interface
        const AMGInterface& interface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        AMGInterfaceField(const AMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const AMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("AMGInterfaceField");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            AMGInterfaceField,
            lduInterface,
            (
                const AMGInterface& AMGCp,
                const lduInterfaceField& fineInterface
            ),
            (AMGCp, fineInterface)
        );


    // Selectors

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<AMGInterfaceField> New
        (
            const AMGInterface& AMGCp,
            const lduInterfaceField& fineInterface
        );


    // Constructors

        //- Construct from AMG interface and fine level interface field
        AMGInterfaceField
        (
            const AMGInterface& AMGCp,
            const lduInterfaceField&
        )
        :
            lduInterfaceField(AMGCp),
            interface_(AMGCp)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

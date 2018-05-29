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
    Foam::SAMGInterfaceField

Description
     Abstract base class for AMG selected interface fields.

SourceFiles
    SAMGInterfaceField.C
    newAmgInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef SAMGInterfaceField_H
#define SAMGInterfaceField_H

#include "lduInterfaceField.H"
#include "SAMGInterface.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class SAMGInterfaceField
:
    public lduInterfaceField
{
    // Private data

        //- Local reference cast into the interface
        const SAMGInterface& interface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        SAMGInterfaceField(const SAMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const SAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("SAMGInterfaceField");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            SAMGInterfaceField,
            lduInterface,
            (
                const SAMGInterface& SAMGCp,
                const lduInterfaceField& fineInterface
            ),
            (SAMGCp, fineInterface)
        );


    // Selectors

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<SAMGInterfaceField> New
        (
            const SAMGInterface& SAMGCp,
            const lduInterfaceField& fineInterface
        );


    // Constructors

        //- Construct from AMG interface and fine level interface field
        SAMGInterfaceField
        (
            const SAMGInterface& SAMGCp,
            const lduInterfaceField&
        )
        :
            lduInterfaceField(SAMGCp),
            interface_(SAMGCp)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

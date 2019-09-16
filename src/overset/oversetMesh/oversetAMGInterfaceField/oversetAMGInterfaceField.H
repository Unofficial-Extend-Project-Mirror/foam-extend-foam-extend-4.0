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
    Foam::oversetAMGInterfaceField

Description
    AMG agglomerated overset interface field.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    oversetAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef oversetAMGInterfaceField_H
#define oversetAMGInterfaceField_H

#include "AMGInterfaceField.H"
#include "oversetAMGInterface.H"
#include "oversetLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class oversetAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class oversetAMGInterfaceField
:
    public AMGInterfaceField,
    virtual public oversetLduInterfaceField
{
    // Private Data Members

        //- Local reference cast into the overset interface
        const oversetAMGInterface& oversetInterface_;


    // Private Member Functions

        // Copy control

            //- Delete default bitwise copy construct
            oversetAMGInterfaceField(const oversetAMGInterfaceField&) = delete;

            //- Delete default bitwise assignment
            void operator=(const oversetAMGInterfaceField&) = delete;


public:

    //- Runtime type information
    TypeName("overset");


    // Constructors

        //- Construct from AMG interface and fine level interface field
        oversetAMGInterfaceField
        (
            const AMGInterface& AMGCp,
            const lduInterfaceField& fineInterfaceField
        );


    //- Destructor
    virtual ~oversetAMGInterfaceField() = default;


    // Member Functions

        // Access functions


        // Coupled interface matrix update

            //- Transform neighbour field
            virtual void transformCoupleField
            (
                scalarField& pnf,
                const direction cmpt
            ) const
            {
                // Does nothing
            }

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField&,
                scalarField&,
                const lduMatrix&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField&,
                scalarField&,
                const lduMatrix&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

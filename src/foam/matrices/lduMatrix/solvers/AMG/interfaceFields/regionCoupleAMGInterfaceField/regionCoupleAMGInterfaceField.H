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
    Foam::regionCoupleAMGInterfaceField

Description
    AMG agglomerated regionCouple interface field.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    regionCoupleAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupleAMGInterfaceField_H
#define regionCoupleAMGInterfaceField_H

#include "regionCoupleAMGInterface.H"
#include "ggiAMGInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class regionCoupleAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupleAMGInterfaceField
:
    public ggiAMGInterfaceField
{
    // Private data

        //- Local reference cast into the regionCouple interface
        const regionCoupleAMGInterface& regionCoupleInterface_;

public:

    //- Runtime type information
    TypeName("regionCouple");


    // Constructors

        //- Construct from AMG interface and fine level interface field
        regionCoupleAMGInterfaceField
        (
            const AMGInterface& AMGCp,
            const lduInterfaceField& fineInterfaceField
        );


    // Destructor

        virtual ~regionCoupleAMGInterfaceField();


    // Member Functions

        // Coupled interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

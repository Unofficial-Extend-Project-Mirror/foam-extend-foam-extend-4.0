/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::regionCoupleGAMGInterfaceField

Description
    GAMG agglomerated regionCouple interface field.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    regionCoupleGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupleGAMGInterfaceField_H
#define regionCoupleGAMGInterfaceField_H

#include "regionCoupleGAMGInterface.H"
#include "ggiGAMGInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class regionCoupleGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupleGAMGInterfaceField
:
    public ggiGAMGInterfaceField
{
    // Private data

        //- Local reference cast into the regionCouple interface
        const regionCoupleGAMGInterface& regionCoupleInterface_;

public:

    //- Runtime type information
    TypeName("regionCouple");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        regionCoupleGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );


    // Destructor

        virtual ~regionCoupleGAMGInterfaceField();


    // Member Functions

        // Interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

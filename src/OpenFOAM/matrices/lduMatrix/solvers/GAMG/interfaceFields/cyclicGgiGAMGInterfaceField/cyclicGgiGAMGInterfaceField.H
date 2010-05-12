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
    Foam::cyclicGgiGAMGInterfaceField

Description
    GAMG agglomerated cyclicGgi interface field.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    cyclicGgiGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicGgiGAMGInterfaceField_H
#define cyclicGgiGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "ggiGAMGInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicGgiGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class cyclicGgiGAMGInterfaceField
:
    public ggiGAMGInterfaceField
{
public:

    //- Runtime type information
    TypeName("cyclicGgi");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        cyclicGgiGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );


    // Destructor

        virtual ~cyclicGgiGAMGInterfaceField();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    overlapGgiLduInterfaceField

Description
    Abstract base class for overlapGgi coupled interface fields.

Author:
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    overlapGgiLduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef overlapGgiLduInterfaceField_H
#define overlapGgiLduInterfaceField_H

#include "primitiveFieldsFwd.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class overlapGgiLduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class overlapGgiLduInterfaceField
{

public:

    //- Runtime type information
    TypeName("overlapGgiLduInterfaceField");


    // Constructors

        //- Construct given coupled patch
        overlapGgiLduInterfaceField()
        {}


    // Destructor

        virtual ~overlapGgiLduInterfaceField();


    // Member Functions

        // Access

            //- Is the transform required
            virtual bool doTransform() const = 0;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const = 0;

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const = 0;

            //- Return rank of component for transform
            virtual int rank() const = 0;


        //- Transform given patch internal field
        void transformCoupleField
        (
            scalarField& psiInternal,
            const direction cmpt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

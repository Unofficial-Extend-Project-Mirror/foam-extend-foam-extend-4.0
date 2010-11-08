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
    BlockLduInterface

Description
    An abstract base class for implicitly-coupled interface fields
    e.g. processor and cyclic patch fields.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    BlockLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduInterface_H
#define BlockLduInterface_H

#include "primitiveFieldsFwd.H"
#include "FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class BlockLduMatrix;

/*---------------------------------------------------------------------------*\
                     Class BlockLduInterface Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduInterface
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockLduInterface(const BlockLduInterface&);

        //- Disallow default bitwise assignment
        void operator=(const BlockLduInterface&);


public:

    //- Runtime type information
    TypeName("BlockLduInterface");


    // Destructor

        virtual ~BlockLduInterface()
        {}


    // Member Functions

        // Coupled interface matrix update

            //- Initialise matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const BlockLduMatrix<Type>& matrix,
                Field<Type>& Ax,
                const Field<Type>& x
            ) const
            {}

            //- Update result based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const BlockLduMatrix<Type>& matrix,
                const FieldField<CoeffField, Type>& coeffs,
                Field<Type>& Ax,
                const Field<Type>& x
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

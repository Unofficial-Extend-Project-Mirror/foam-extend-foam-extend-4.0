/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
\    /   O peration     |
 \  /    A nd           | Copyright held by original author
  \/     M anipulation  |
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
    BlockLduInterfaceField

Description
    Abstract base class for interface fields with block coefficients

Author
    Ivor Clifford, ivor.clifford@gmail.com

SourceFiles
    BlockLduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduInterfaceField_H
#define BlockLduInterfaceField_H

#include "lduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class BlockLduMatrix;

template<class Type>
class CoeffField;

/*---------------------------------------------------------------------------*\
                           Class BlockLduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduInterfaceField
:
    public lduInterfaceField
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockLduInterfaceField(const BlockLduInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const BlockLduInterfaceField&);


public:

    //- Runtime Type information
    TypeName("BlockLduInterfaceField");

    // Constructors

        //- Construct given coupled patch
        BlockLduInterfaceField(const lduInterface& patch)
        :
            lduInterfaceField(patch)
        {}


    // Destructor

        virtual ~BlockLduInterfaceField()
        {}



    // Member Functions

        // Coupled interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
            ) const
            {}

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

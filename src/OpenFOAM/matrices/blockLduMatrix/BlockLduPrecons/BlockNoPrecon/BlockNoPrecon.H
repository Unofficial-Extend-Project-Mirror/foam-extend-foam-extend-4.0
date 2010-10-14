/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    BlockNoPrecon

Description
    No preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    BlockNoPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef BlockNoPrecon_H
#define BlockNoPrecon_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class BlockLduPrecon;

/*---------------------------------------------------------------------------*\
                       Class BlockNoPrecon Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockNoPrecon
:
    public BlockLduPrecon<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockNoPrecon(const BlockNoPrecon&);

        //- Disallow default bitwise assignment
        void operator=(const BlockNoPrecon&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        BlockNoPrecon
        (
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        )
        :
            BlockLduPrecon<Type>(matrix)
        {}


    // Destructor

        virtual ~BlockNoPrecon()
        {}


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            Field<Type>& x,
            const Field<Type>& b
        ) const
        {
            x = b;
        }

        //- Execute preconditioning
        virtual void preconditionT
        (
            Field<Type>& xT,
            const Field<Type>& bT
        ) const
        {
            precondition(xT, bT);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

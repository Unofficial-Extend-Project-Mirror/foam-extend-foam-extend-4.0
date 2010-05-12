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
    symGaussSeidelSmoother

Description
    Symmetric Gauss-Seidel smoother

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    symGaussSeidelSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef symGaussSeidelSmoother_H
#define symGaussSeidelSmoother_H

#include "symGaussSeidelPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class symGaussSeidelSmoother Declaration
\*---------------------------------------------------------------------------*/

class symGaussSeidelSmoother
:
    public lduSmoother
{
    // Private Data

        //- Gauss-Seidel preconditioner
        symGaussSeidelPrecon gs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        symGaussSeidelSmoother(const symGaussSeidelSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const symGaussSeidelSmoother&);


public:

    //- Runtime type information
    TypeName("symGaussSeidel");


    // Constructors

        //- Construct from components
        symGaussSeidelSmoother
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );


    // Destructor

        virtual ~symGaussSeidelSmoother()
        {}


    // Member Functions

        //- Execute smoothing
        virtual void smooth
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt,
            const label nSweeps
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

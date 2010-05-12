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
    coupledGaussSeidelSmoother

Description
    GaussSeidel smoother for coupled matrix

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    coupledGaussSeidelSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef coupledGaussSeidelSmoother_H
#define coupledGaussSeidelSmoother_H

#include "coupledLduSmoother.H"
#include "coupledGaussSeidelPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class coupledGaussSeidelSmoother Declaration
\*---------------------------------------------------------------------------*/

class coupledGaussSeidelSmoother
:
    public coupledLduSmoother
{
    // Private Data

        // GaussSeidel preconditioner
        coupledGaussSeidelPrecon gs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledGaussSeidelSmoother(const coupledGaussSeidelSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const coupledGaussSeidelSmoother&);


public:

    //- Runtime type information
    TypeName("GaussSeidel");


    // Constructors

        //- Construct from matrix
        coupledGaussSeidelSmoother
        (
            const coupledLduMatrix& matrix,
            const PtrList<FieldField<Field, scalar> >& bouCoeffs,
            const PtrList<FieldField<Field, scalar> >& intCoeffs,
            const lduInterfaceFieldPtrsListList& interfaces
        );


    // Destructor

        virtual ~coupledGaussSeidelSmoother()
        {}


    // Member Functions

        //- Execute smoothing.
        virtual void smooth
        (
            FieldField<Field, scalar>& x,
            const FieldField<Field, scalar>& b,
            const direction cmpt,
            const label nSweeps
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    crMatrix

Description
    Sparse matrix in compressed row format

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    crMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef crMatrix_H
#define crMatrix_H

#include "crAddressing.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class crMatrix Declaration
\*---------------------------------------------------------------------------*/

class crMatrix
:
    public refCount
{
    // Private data

        //- Addressing
        crAddressing crAddr_;

        //- Coefficients
        scalarField coeffs_;


public:

    // Constructors

        //- Construct from addressing
        crMatrix(const crAddressing& addr);

        //- Construct from row and column size and row count
        crMatrix
        (
            const label nRows,
            const label nCols,
            const labelList& count
        );

        //- Construct from components of addressing
        crMatrix
        (
            const label nRows,
            const label nCols,
            const labelList& row,
            const labelList& col
        );

        //- Construct as copy
        crMatrix(const crMatrix& m);

        //- Construct as copy of tmp<crMatrix> deleting argument
        crMatrix(const tmp<crMatrix>& tm);

        //- Construct from Istream
        crMatrix(Istream&);


    // Destructor - default


    // Member Functions

        // Access

            //- Return addressing
            const crAddressing& crAddr() const
            {
                return crAddr_;
            }

            //- Return coefficients
            const scalarField& coeffs() const
            {
                return coeffs_;
            }

            //- Return coefficients to be set
            scalarField& coeffs()
            {
                return coeffs_;
            }

            //- Return column array to be set
            labelList& col()
            {
                return crAddr_.col();
            }


        // Matrix operations

            //- Return transpose matrix
            tmp<crMatrix> T() const;

        // Matrix-vector operations

            //- Calculate b += A*x
            void dotPlus(scalarField& b, const scalarField& x) const;


    // Member operators

        void operator=(const crMatrix&);
        void operator=(const tmp<crMatrix>&);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const crMatrix&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

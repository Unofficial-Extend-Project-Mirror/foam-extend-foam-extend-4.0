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
    Foam::DICGaussSeidelSmoother

Description
    Combined DIC/GaussSeidel smoother for symmetric matrices in which DIC
    smoothing is followed by GaussSeidel to ensure that any "spikes" created
    by the DIC sweeps are smoothed-out.

SourceFiles
    DICGaussSeidelSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef DICGaussSeidelSmoother_H
#define DICGaussSeidelSmoother_H

#include "DICSmoother.H"
#include "GaussSeidelSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class DICGaussSeidelSmoother Declaration
\*---------------------------------------------------------------------------*/

class DICGaussSeidelSmoother
:
    public lduSmoother
{
    // Private data

        DICSmoother dicSmoother_;

        GaussSeidelSmoother gsSmoother_;


public:

    //- Runtime type information
    TypeName("DICGaussSeidel");


    // Constructors

        //- Construct from matrix components
        DICGaussSeidelSmoother
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );


    // Member Functions

        //- Execute smoothing
        virtual void smooth
        (
            scalarField& psi,
            const scalarField& Source,
            const direction cmpt,
            const label nSweeps
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

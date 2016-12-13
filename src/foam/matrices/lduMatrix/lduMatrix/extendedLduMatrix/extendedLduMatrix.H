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
    Foam::extendedLduMatrix

Description
    Extension of lower and upper coeffs of lduMatrix. Uses faceMap from
    extendedLduAddressing to define extended lower and upper coeffs from
    ordinary lower and upper coeffs given by ordinary lduMatrix. Used in ILUCp
    preconditioner.

SourceFiles
    extendedLduMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef extendedLduMatrix_H
#define extendedLduMatrix_H

#include "extendedLduAddressing.H"
#include "lduMatrix.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                      Class extendedLduMatrix Declaration
\*---------------------------------------------------------------------------*/

class extendedLduMatrix
{
    // Private data

        //- Reference to underlying lduMatrix
        const lduMatrix& basicLduMatrix_;

        //- Reference to extended ldu addressing
        const extendedLduAddressing& extLduAddr_;

        //- Extended lower coeffs
        scalarField* extendedLowerPtr_;

        //- Extended upper coeffs
        scalarField* extendedUpperPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        extendedLduMatrix(const extendedLduMatrix&);

        //- Disallow default bitwise assignement
        void operator=(const extendedLduMatrix&);


        // Helper functions

            //- Clear all demand driven data - helper function
            void clearOut();


public:


    ClassName("extendedLduMatrix");


    // Constructors

        //- Construct given lduMatrix, extension level and polyMesh
        extendedLduMatrix
        (
            const lduMatrix&,
            const extendedLduAddressing&
        );


    //- Destructor
    ~extendedLduMatrix();


    // Member functions

        // Access

            //- Const access to underlying basic lduMatrix
            const lduMatrix& basicLduMatrix() const
            {
                return basicLduMatrix_;
            };

            //- Const access to extendedLduAddressing
            const extendedLduAddressing& extendedLduAddr() const
            {
                return extLduAddr_;
            };

            //- Non-const access to extended lower coeffs
            scalarField& extendedLower();

            //- Non-const access to extended upper coeffs
            scalarField& extendedUpper();

            //- Const access to extended lower coeffs
            const scalarField& extendedLower() const;

            //- Const access to extended upper coeffs
            const scalarField& extendedUpper() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::extendedBlockLduMatrix

Description
    Extension of lower and upper coeffs of BlockLduMatrix. Uses faceMap from
    extendedLduAddressing to define extended lower and upper coeffs from
    ordinary lower and upper coeffs given by ordinary BlockLduMatrix. Used in
    ILUCp preconditioner.

SourceFiles
    extendedBlockLduMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef extendedBlockLduMatrix_H
#define extendedBlockLduMatrix_H

#include "extendedLduAddressing.H"
#include "BlockLduMatrix.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                   Class extendedBlockLduMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class extendedBlockLduMatrix
{
    // Public data type

        typedef CoeffField<Type> TypeCoeffField;


private:

    // Private data

        //- Reference to underlying BlockLduMatrix
        const BlockLduMatrix<Type>& basicBlockLduMatrix_;

        //- Reference to extended ldu addressing
        const extendedLduAddressing& extLduAddr_;

        //- Extended lower coeffs
        CoeffField<Type>* extendedLowerPtr_;

        //- Extended upper coeffs
        CoeffField<Type>* extendedUpperPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        extendedBlockLduMatrix(const extendedBlockLduMatrix<Type>&);

        //- Disallow default bitwise assignement
        void operator=(const extendedBlockLduMatrix<Type>&);


        // Helper functions

            //- Clear all demand driven data - helper function
            void clearOut();

            //- Map upper and lower coeffs from ordinary block matrix to
            //  extended block matrix
            void mapOffDiagCoeffs(const BlockLduMatrix<Type>&);


public:


    ClassName("extendedBlockLduMatrix");


    // Constructors

        //- Construct given BlockduMatrix, extension level and polyMesh
        extendedBlockLduMatrix
        (
            const BlockLduMatrix<Type>&,
            const extendedLduAddressing&
        );


    //- Destructor
    ~extendedBlockLduMatrix();


    // Member functions

        // Access

            //- Const access to underlying basic BlockLduMatrix
            const lduMatrix& basicBlockLduMatrix() const
            {
                return basicBlockLduMatrix_;
            };

            //- Const access to extendedLduAddressing
            const extendedLduAddressing& extendedLduAddr() const
            {
                return extLduAddr_;
            };

            //- Non-const access to extended lower coeffs
            TypeCoeffField& extendedLower();

            //- Non-const access to extended upper coeffs
            TypeCoeffField& extendedUpper();

            //- Const access to extended lower coeffs
            const TypeCoeffField& extendedLower() const;

            //- Const access to extended upper coeffs
            const TypeCoeffField& extendedUpper() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "extendedBlockLduMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

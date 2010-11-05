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
    TensorN

Description

SourceFiles
    TensorNI.H

\*---------------------------------------------------------------------------*/

#ifndef TensorN_H
#define TensorN_H

#include "VectorN.H"
#include "SphericalTensorN.H"
#include "DiagTensorN.H"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

template <class Cmpt, int length>
class DiagTensorN;

template <class Cmpt, int length>
class SphericalTensorN;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class TensorN Declaration
\*---------------------------------------------------------------------------*/

template <class Cmpt, int length>
class TensorN
:
    public VectorSpace<TensorN<Cmpt, length>, Cmpt, length*length>
{

public:

    // Member constants

        enum
        {
            rank = 2,            // Rank of TensorN is 2
            rowLength = length   // Number of components in a row
        };


    // Static data members

        static const char* const typeName;

        static const TensorN zero;
        static const TensorN one;


    // Constructors

        //- Construct null
        inline TensorN();

        //- Construct given VectorSpace
        inline TensorN
        (
            const VectorSpace<TensorN<Cmpt, length>, Cmpt, length*length>&
        );

        //- Construct given component value.  Special use only!
        explicit inline TensorN(const Cmpt& tx);

        //- Construct from Istream
        TensorN(Istream&);


    // Member Functions

        //- Return (i, j) component
        inline const Cmpt& operator()
        (
            const direction i,
            const direction j
        ) const;

        //- Return access to (i, j) component
        inline Cmpt& operator()
        (
            const direction i,
            const direction j
        );

        //- Diagonal
        inline DiagTensorN<Cmpt, length> diag() const;

        //- Transpose
        inline TensorN<Cmpt, length> T() const;

        //- Negative sum the vertical off-diagonal components
        inline TensorN<Cmpt, length> negSumDiag() const;

    // Member Operators

        //- Assign to a SphericalTensorN
        inline void operator=(const SphericalTensorN<Cmpt, length>&);

        //- Assign to a DiagTensorN
        inline void operator=(const DiagTensorN<Cmpt, length>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "TensorNI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::SymmTensor2D

Description
    Templated 2D symmetric tensor derived from VectorSpace adding construction
    from 4 components, element access using xx(), xy() etc. member functions
    and the inner-product (dot-product) and outer-product of two Vectors
    (tensor-product) operators.

SourceFiles
    SymmTensor2DTemplateI.H

\*---------------------------------------------------------------------------*/

#ifndef SymmTensor2D_H
#define SymmTensor2D_H

#include "VectorTemplate.H"
#include "SphericalTensor2DTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SymmTensor2D Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class SymmTensor2D
:
    public VectorSpace<SymmTensor2D<Cmpt>, Cmpt, 3>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef SymmTensor2D<label> labelType;


    // Member constants

        enum
        {
            rank = 2 // Rank of SymmTensor2D is 2
        };


    // Static data members

        static const char* const typeName;
        static const char* componentNames[];

        static const SymmTensor2D zero;
        static const SymmTensor2D one;
        static const SymmTensor2D max;
        static const SymmTensor2D min;
        static const SymmTensor2D I;


    //- Component labeling enumeration
    enum components { XX, XY, YY };


    // Constructors

        //- Construct null
        inline SymmTensor2D();

        //- Construct given VectorSpace
        inline SymmTensor2D(const VectorSpace<SymmTensor2D<Cmpt>, Cmpt, 3>&);

        //- Construct given SphericalTensor
        inline SymmTensor2D(const SphericalTensor2D<Cmpt>&);

        //- Construct given the three components
        inline SymmTensor2D
        (
            const Cmpt txx, const Cmpt txy,
                            const Cmpt tyy
        );

        //- Construct from Istream
        SymmTensor2D(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& xx() const;
            inline const Cmpt& xy() const;
            inline const Cmpt& yy() const;

            inline Cmpt& xx();
            inline Cmpt& xy();
            inline Cmpt& yy();

        //- Transpose
        inline const SymmTensor2D<Cmpt>& T() const;


    // Member Operators

        //- Construct given SphericalTensor2D
        inline void operator=(const SphericalTensor2D<Cmpt>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "SymmTensor2DTemplateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

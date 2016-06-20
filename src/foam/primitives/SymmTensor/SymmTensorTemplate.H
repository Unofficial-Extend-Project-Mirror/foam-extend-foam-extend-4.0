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
    Foam::SymmTensor

Description
    Templated 3D symmetric tensor derived from VectorSpace adding construction
    from 6 components, element access using xx(), xy() etc. member functions
    and the inner-product (dot-product) and outer-product of two Vectors
    (tensor-product) operators.

SourceFiles
    SymmTensorITemplate.H

\*---------------------------------------------------------------------------*/

#ifndef SymmTensorTemplate_H
#define SymmTensorTemplate_H

#include "VectorSpace.H"
#include "SphericalTensorTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SymmTensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class SymmTensor
:
    public VectorSpace<SymmTensor<Cmpt>, Cmpt, 6>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef SymmTensor<label> labelType;


    // Member constants

        enum
        {
            rank = 2 // Rank of SymmTensor is 2
        };


    // Static data members

        static const char* const typeName;
        static const char* componentNames[];

        static const SymmTensor zero;
        static const SymmTensor one;
        static const SymmTensor max;
        static const SymmTensor min;
        static const SymmTensor I;


    //- Component labeling enumeration
    enum components { XX, XY, XZ, YY, YZ, ZZ };


    // Constructors

        //- Construct null
        inline SymmTensor();

        //- Construct given VectorSpace
        inline SymmTensor(const VectorSpace<SymmTensor<Cmpt>, Cmpt, 6>&);

        //- Construct given SphericalTensor
        inline SymmTensor(const SphericalTensor<Cmpt>&);

        //- Construct given the six components
        inline SymmTensor
        (
            const Cmpt txx, const Cmpt txy, const Cmpt txz,
                            const Cmpt tyy, const Cmpt tyz,
                                            const Cmpt tzz
        );

        //- Construct from Istream
        SymmTensor(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& xx() const;
            inline const Cmpt& xy() const;
            inline const Cmpt& xz() const;
            inline const Cmpt& yy() const;
            inline const Cmpt& yz() const;
            inline const Cmpt& zz() const;

            inline Cmpt& xx();
            inline Cmpt& xy();
            inline Cmpt& xz();
            inline Cmpt& yy();
            inline Cmpt& yz();
            inline Cmpt& zz();

        //- Transpose
        inline const SymmTensor<Cmpt>& T() const;


    // Member Operators

        //- Construct given SphericalTensor
        inline void operator=(const SphericalTensor<Cmpt>&);
};


template<class Cmpt>
class symmTypeOfRank<Cmpt, 2>
{
public:

    typedef SymmTensor<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "SymmTensorTemplateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

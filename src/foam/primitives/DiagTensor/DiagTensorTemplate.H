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
    Foam::DiagTensor

Description
    Templated 3D DiagTensor derived from VectorSpace.

    Adding construction from 3 components, element access using xx(), yy()
    and zz() member functions and the inner-product (dot-product) and
    outer-product operators.

SourceFiles
    DiagTensorTemplateI.H

\*---------------------------------------------------------------------------*/

#ifndef DiagTensorTemplate_H
#define DiagTensorTemplate_H

#include "TensorTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DiagTensor Declaration
\*---------------------------------------------------------------------------*/

template <class Cmpt>
class DiagTensor
:
    public VectorSpace<DiagTensor<Cmpt>, Cmpt, 3>
{

public:

    // Member constants

        enum
        {
            rank = 2 // Rank of DiagTensor is 2
        };


    // Static data members

        static const char* const typeName;
        static const char* componentNames[];
        static const DiagTensor zero;
        static const DiagTensor one;
        static const DiagTensor max;
        static const DiagTensor min;


    //- Component labeling enumeration
    enum components { XX, YY, ZZ };


    // Constructors

        //- Construct null
        inline DiagTensor();

        //- Construct given VectorSpace
        inline DiagTensor(const VectorSpace<DiagTensor<Cmpt>, Cmpt, 3>&);

        //- Construct given three components
        inline DiagTensor(const Cmpt& txx, const Cmpt& tyy, const Cmpt& tzz);

        //- Construct from Istream
        inline DiagTensor(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& xx() const;
            inline const Cmpt& yy() const;
            inline const Cmpt& zz() const;

            inline Cmpt& xx();
            inline Cmpt& yy();
            inline Cmpt& zz();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "DiagTensorTemplateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    CoeffField<T>

Description
    Template specialisation for symmTensor4thOrder coefficients

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    symmTensor4thOrderCoeffField.C

\*---------------------------------------------------------------------------*/

#ifndef symmTensor4thOrderCoeffField_H
#define symmTensor4thOrderCoeffField_H

#include "CoeffField.H"
#include "DecoupledCoeffField.H"
#include "symmTensor4thOrder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
class CoeffField<symmTensor4thOrder>
:
    public DecoupledCoeffField<symmTensor4thOrder>
{
public:

    // Constructors

        //- Construct given size
        explicit CoeffField(const label);

        //- Construct as copy
        CoeffField(const CoeffField<symmTensor4thOrder>&);

        //- Construct as copy of base
        explicit CoeffField(const DecoupledCoeffField<symmTensor4thOrder>&);

        //- Construct as copy of base
        explicit CoeffField
        (
            const tmp<DecoupledCoeffField<symmTensor4thOrder> >&
        );

        //- Construct from Istream
        explicit CoeffField(Istream&);


    // Member operators

        void operator=(const CoeffField<symmTensor4thOrder>&);
        void operator=(const tmp<CoeffField<symmTensor4thOrder> >&);

        void operator=(const scalarTypeField&);
        void operator=(const tmp<scalarTypeField>&);
        void operator=(const linearTypeField&);
        void operator=(const tmp<linearTypeField>&);

    // IOstream operators

        friend Ostream& operator<<
        (
            Ostream&,
            const CoeffField<symmTensor4thOrder>&
        );

        friend Ostream& operator<<
        (
            Ostream&,
            const tmp<CoeffField<symmTensor4thOrder> >&
        );
};


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

tmp<CoeffField<symmTensor4thOrder> >
inv(const CoeffField<symmTensor4thOrder>& f);

template<>
void multiply
(
    Field<symmTensor4thOrder>& f,
    const CoeffField<symmTensor4thOrder>& f1,
    const symmTensor4thOrder& f2
);

template<>
void multiply
(
    Field<symmTensor4thOrder>& f,
    const CoeffField<symmTensor4thOrder>& f1,
    const Field<symmTensor4thOrder>& f2
);

template<>
void multiply
(
    Field<symmTensor4thOrder>& f,
    const Field<symmTensor4thOrder>& f1,
    const CoeffField<symmTensor4thOrder>& f2
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

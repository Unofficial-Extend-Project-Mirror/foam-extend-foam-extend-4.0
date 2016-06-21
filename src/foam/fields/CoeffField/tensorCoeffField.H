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
    Template specialisation for tensor coefficients

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    tensorCoeffField.C

\*---------------------------------------------------------------------------*/

#ifndef tensorCoeffField_H
#define tensorCoeffField_H

#include "CoeffField.H"
#include "DecoupledCoeffField.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
class CoeffField<tensor>
:
    public DecoupledCoeffField<tensor>
{
public:

    // Constructors

        //- Construct given size
        explicit CoeffField(const label);

        //- Construct as copy
        CoeffField(const CoeffField<tensor>&);

        //- Construct as copy of base
        explicit CoeffField(const DecoupledCoeffField<tensor>&);

        //- Construct as copy of base
        explicit CoeffField(const tmp<DecoupledCoeffField<tensor> >&);

        //- Construct from Istream
        explicit CoeffField(Istream&);


    // Member operators

        void operator=(const CoeffField<tensor>&);
        void operator=(const tmp<CoeffField<tensor> >&);

        void operator=(const scalarTypeField&);
        void operator=(const tmp<scalarTypeField>&);
        void operator=(const linearTypeField&);
        void operator=(const tmp<linearTypeField>&);


    // IOstream operators

        friend Ostream& operator<<
        (
            Ostream&,
            const CoeffField<tensor>&
        );

        friend Ostream& operator<<
        (
            Ostream&,
            const tmp<CoeffField<tensor> >&
        );
};


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

template<>
tmp<CoeffField<tensor> > inv(const CoeffField<tensor>& f);

template<>
void multiply
(
    Field<tensor>& f,
    const CoeffField<tensor>& f1,
    const tensor& f2
);

template<>
void multiply
(
    Field<tensor>& f,
    const CoeffField<tensor>& f1,
    const Field<tensor>& f2
);

template<>
void multiply
(
    Field<tensor>& f,
    const Field<tensor>& f1,
    const CoeffField<tensor>& f2
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Template specialisation for diagTensor coefficients

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    diagTensorCoeffField.C

\*---------------------------------------------------------------------------*/

#ifndef diagTensorCoeffField_H
#define diagTensorCoeffField_H

#include "CoeffField.H"
#include "DecoupledCoeffField.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
class CoeffField<diagTensor>
:
    public DecoupledCoeffField<diagTensor>
{
public:

    // Constructors

        //- Construct given size
        explicit CoeffField(const label);

        //- Construct as copy
        CoeffField(const CoeffField<diagTensor>&);

        //- Construct as copy of base
        explicit CoeffField(const DecoupledCoeffField<diagTensor>&);

        //- Construct as copy of base
        explicit CoeffField(const tmp<DecoupledCoeffField<diagTensor> >&);

        //- Construct from Istream
        explicit CoeffField(Istream&);


    // Member operators

        void operator=(const CoeffField<diagTensor>&);
        void operator=(const tmp<CoeffField<diagTensor> >&);

        void operator=(const scalarTypeField&);
        void operator=(const tmp<scalarTypeField>&);
        void operator=(const linearTypeField&);
        void operator=(const tmp<linearTypeField>&);

    // IOstream operators

        friend Ostream& operator<<
        (
            Ostream&,
            const CoeffField<diagTensor>&
        );

        friend Ostream& operator<<
        (
            Ostream&,
            const tmp<CoeffField<diagTensor> >&
        );
};


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

tmp<CoeffField<diagTensor> > inv(const CoeffField<diagTensor>& f);

template<>
void multiply
(
    Field<diagTensor>& f,
    const CoeffField<diagTensor>& f1,
    const diagTensor& f2
);

template<>
void multiply
(
    Field<diagTensor>& f,
    const CoeffField<diagTensor>& f1,
    const Field<diagTensor>& f2
);

template<>
void multiply
(
    Field<diagTensor>& f,
    const Field<diagTensor>& f1,
    const CoeffField<diagTensor>& f2
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

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

};


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

tmp<CoeffField<tensor> > inv(const CoeffField<tensor>& f);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

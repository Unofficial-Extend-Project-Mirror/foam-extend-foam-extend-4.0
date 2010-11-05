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
    Template specialisation for sphericalTensor coefficients

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    sphericalTensorCoeffField.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalTensorCoeffField_H
#define sphericalTensorCoeffField_H

#include "CoeffField.H"
#include "DecoupledCoeffField.H"
#include "sphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
class CoeffField<sphericalTensor>
:
    public DecoupledCoeffField<sphericalTensor>
{
public:

    // Constructors

        //- Construct given size
        explicit CoeffField(const label);

        //- Construct as copy
        CoeffField(const CoeffField<sphericalTensor>&);

        //- Construct as copy of base
        explicit CoeffField(const DecoupledCoeffField<sphericalTensor>&);

        //- Construct as copy of base
        explicit CoeffField(const tmp<DecoupledCoeffField<sphericalTensor> >&);

        //- Construct from Istream
        explicit CoeffField(Istream&);


    // Member operators

        void operator=(const CoeffField<sphericalTensor>&);
        void operator=(const tmp<CoeffField<sphericalTensor> >&);

        void operator=(const scalarTypeField&);
        void operator=(const tmp<scalarTypeField>&);
        void operator=(const linearTypeField&);
        void operator=(const tmp<linearTypeField>&);

};


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

tmp<CoeffField<sphericalTensor> > inv(const CoeffField<sphericalTensor>& f);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

Description
    Function to calculate eigen values and eigen vectors of volSymmTensorField
    Using the main procedure/code from here:
    http://barnesc.blogspot.ie/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
    and code here:
    http://www.connellybarnes.com/code/c/eig3-1.0.0.zip
    Note: built-in OpenFOAM functions mess-up on a number of different tensors.

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef eig3Field_H
#define eig3Field_H

//#include "fvCFD.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


void eig3Field
(
    const volTensorField& A, volTensorField& V, volVectorField& d
);


void eig3Field
(
    const volSymmTensorField& A, volTensorField& V, volVectorField& d
);


void eig3Field
(
    const volSymmTensorField& A, volTensorField& V, volDiagTensorField& d
);


void eig3Field
(
    const pointSymmTensorField& A,
    pointTensorField& V,
    pointVectorField& d
);


tmp<volSymmTensorField> log(const volSymmTensorField& vf);

symmTensor log(const symmTensor& t);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

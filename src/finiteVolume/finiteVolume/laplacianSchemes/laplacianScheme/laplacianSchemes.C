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

\*---------------------------------------------------------------------------*/

#include "laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define the constructor function hash tables

#define makeLaplacianGTypeScheme(Type, GType)                                 \
    typedef laplacianScheme<Type, GType> laplacianScheme##Type##GType;        \
    defineTemplateRunTimeSelectionTable(laplacianScheme##Type##GType, Istream);

#define makeLaplacianScheme(Type)                                             \
    makeLaplacianGTypeScheme(Type, scalar);                                   \
    makeLaplacianGTypeScheme(Type, symmTensor);                               \
    makeLaplacianGTypeScheme(Type, diagTensor);                               \
    makeLaplacianGTypeScheme(Type, tensor);

makeLaplacianScheme(scalar);
makeLaplacianScheme(vector);
makeLaplacianScheme(sphericalTensor);
makeLaplacianScheme(symmTensor);
makeLaplacianScheme(diagTensor);
makeLaplacianScheme(tensor);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "equationSource.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#ifdef __INTEL_COMPILER
#   define defineEquationSourceTypeNames(type, Type)                          \
    const char* const Foam::equationSource<type>::typeName =                  \
        "equation" #Type "Source"
#else
#   define defineEquationSourceTypeNames(type, Type)                          \
    template<>                                                                \
    const char* const Foam::equationSource<type>::typeName =                  \
        "equation" #Type "Source"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineEquationSourceTypeNames(scalar, Scalar);
    defineEquationSourceTypeNames(vector, Vector);
    defineEquationSourceTypeNames(tensor, Tensor);
    defineEquationSourceTypeNames(diagTensor, DiagTensor);
    defineEquationSourceTypeNames(symmTensor, SymmTensor);
    defineEquationSourceTypeNames(sphericalTensor, SphericalTensor);

} // End namespace Foam

#include "equationScalarSource.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

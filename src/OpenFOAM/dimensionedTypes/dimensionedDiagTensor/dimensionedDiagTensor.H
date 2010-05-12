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

Class
    dimensionedDiagTensor

Description
    Dimensioned diagonal tensor obtained from generic dimensioned type.

Author
    Dubravko Matijasevic
    Hrvoje Jasak

SourceFiles
    dimensionedDiagTensor.C

\*---------------------------------------------------------------------------*/

#ifndef dimensionedDiagTensor_H
#define dimensionedDiagTensor_H

#include "dimensionedVector.H"
#include "dimensionedTensor.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef dimensioned<diagTensor> dimensionedDiagTensor;


// friend functions

dimensionedScalar tr(const dimensionedDiagTensor&);
dimensionedScalar det(const dimensionedDiagTensor&);
dimensionedDiagTensor inv(const dimensionedDiagTensor&);
dimensionedDiagTensor diag(const dimensionedTensor&);

// friend operators


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


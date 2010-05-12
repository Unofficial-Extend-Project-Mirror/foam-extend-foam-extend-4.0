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

Description
    Dimensioned diagonal tensor obtained from generic dimensioned type.

\*---------------------------------------------------------------------------*/

#include "dimensionedDiagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

dimensionedScalar tr(const dimensionedDiagTensor& ddt)
{
    return dimensionedScalar
    (
        "tr(" + ddt.name() + ')',
        ddt.dimensions(),
        tr(ddt.value())
    );
}


dimensionedScalar det(const dimensionedDiagTensor& ddt)
{
    return dimensionedScalar
    (
        "det(" + ddt.name() + ')',
        pow(ddt.dimensions(), tensor::dim),
        det(ddt.value())
    );
}


dimensionedDiagTensor inv(const dimensionedDiagTensor& ddt)
{
    return dimensionedDiagTensor
    (
        "inv(" + ddt.name() + ')',
        dimless/ddt.dimensions(),
        inv(ddt.value())
    );
}


dimensionedDiagTensor diag(const dimensionedTensor& dt)
{
    return dimensionedDiagTensor
    (
        "diag(" + dt.name() + ')',
        dt.dimensions(),
        diag(dt.value())
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

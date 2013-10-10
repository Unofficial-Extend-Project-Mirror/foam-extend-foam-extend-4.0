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

//#include "equationScalarSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
label Foam::equationSource<scalar>::lookupComponentIndex
(
    const word& componentName
) const
{
    // scalar specialization: also returns 0 if word::null is given
    if ((componentName == word::null) || componentName == "x")
    {
        return 0;
    }

    return -1;
}


template<>
const scalar& equationSource<scalar>::singleValue
(
    label sourceIndex,
    label componentIndex
) const
{
    return singles_[sourceIndex];
}

template<>
const scalar& equationSource<scalar>::fieldValue
(
    label sourceIndex,
    label componentIndex,
    label cellIndex,
    label geoIndex
) const
{
    return fields_[sourceIndex][geoIndex][cellIndex];
}

template<>
void equationSource<scalar>::fullFieldValue
(
    scalarField& result,
    label sourceIndex,
    label componentIndex,
    label geoIndex
) const
{
    //result.setSize(fields_[sourceIndex][geoIndex].size());
    result = fields_[sourceIndex][geoIndex];
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

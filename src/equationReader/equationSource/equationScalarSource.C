/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

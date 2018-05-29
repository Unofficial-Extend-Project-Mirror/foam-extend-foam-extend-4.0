/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "regionSizeDistribution.H"
#include "regionSplit.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::Map<Type> Foam::regionSizeDistribution::regionSum
(
    const regionSplit& regions,
    const Field<Type>& fld
) const
{
    // Per region the sum of fld
    Map<Type> regionToSum(regions.nRegions()/Pstream::nProcs());

    forAll(fld, cellI)
    {
        label regionI = regions[cellI];

        typename Map<Type>::iterator fnd = regionToSum.find(regionI);
        if (fnd == regionToSum.end())
        {
            regionToSum.insert(regionI, fld[cellI]);
        }
        else
        {
            fnd() += fld[cellI];
        }
    }
    Pstream::mapCombineGather(regionToSum, plusEqOp<Type>());
    Pstream::mapCombineScatter(regionToSum);

    return regionToSum;
}


// Get data in sortedToc order
template<class Type>
Foam::List<Type> Foam::regionSizeDistribution::extractData
(
    const UList<label>& keys,
    const Map<Type>& regionData
) const
{
    List<Type> sortedData(keys.size());

    forAll(keys, i)
    {
        sortedData[i] = regionData[keys[i]];
    }
    return sortedData;
}


// ************************************************************************* //

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

#include "channelIndex.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::Field<T> Foam::channelIndex::regionSum(const Field<T>& cellField) const
{
    Field<T> regionField(cellRegion_().nRegions(), pTraits<T>::zero);

    forAll(cellRegion_(), cellI)
    {
        regionField[cellRegion_()[cellI]] += cellField[cellI];
    }

    // Global sum
    Pstream::listCombineGather(regionField, plusEqOp<T>());
    Pstream::listCombineScatter(regionField);

    return regionField;
}


template<class T>
Foam::Field<T> Foam::channelIndex::collapse
(
    const Field<T>& cellField,
    const bool asymmetric
) const
{
    // Average and order
    const Field<T> summedField(regionSum(cellField));

    Field<T> regionField
    (
        summedField
      / regionCount_,
        sortMap_
    );

    // Symmetry?
    if (symmetric_)
    {
        label nlb2 = cellRegion_().nRegions()/2;

        if (asymmetric)
        {
            for (label j=0; j<nlb2; j++)
            {
                regionField[j] =
                    0.5
                  * (
                        regionField[j]
                      - regionField[cellRegion_().nRegions() - j - 1]
                    );
            }
        }
        else
        {
            for (label j=0; j<nlb2; j++)
            {
                regionField[j] =
                    0.5
                  * (
                        regionField[j]
                      + regionField[cellRegion_().nRegions() - j - 1]
                    );
            }
        }

        regionField.setSize(nlb2);
    }

    return regionField;
}


// ************************************************************************* //

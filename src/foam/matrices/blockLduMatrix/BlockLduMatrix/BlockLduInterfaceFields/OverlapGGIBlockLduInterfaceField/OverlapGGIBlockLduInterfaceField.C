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

#include "OverlapGGIBlockLduInterfaceField.H"
#include "diagTensorField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::OverlapGGIBlockLduInterfaceField<Type>::
~OverlapGGIBlockLduInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::OverlapGGIBlockLduInterfaceField<Type>::transformCoupleField
(
    scalarField& f,
    const direction cmpt
) const
{
    // KRJ: 2013-02-08: Transform not tested
    if (doTransform())
    {
        if (forwardT().size() == 1)
        {
            f *= pow(diag(forwardT()[0]).component(cmpt), pTraits<Type>::rank);
        }
        else
        {
            f *= pow(diag(forwardT())().component(cmpt), pTraits<Type>::rank);
        }
    }
}


template<class Type>
void Foam::OverlapGGIBlockLduInterfaceField<Type>::transformCoupleField
(
    Field<Type>& f
) const
{
    // KRJ: 2013-02-08: Transform not tested
    if (doTransform())
    {
        if (forwardT().size() == 1)
        {
            transform(f, forwardT()[0], f);
        }
        else
        {
            transform(f, forwardT(), f);
        }
    }
}


// ************************************************************************* //

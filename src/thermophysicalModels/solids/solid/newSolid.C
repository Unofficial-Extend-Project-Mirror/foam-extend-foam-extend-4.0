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

#include "solid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solid> Foam::solid::New(Istream& is)
{
    if (debug)
    {
        Info<< "solid::New(Istream&): "
            << "constructing solid"
            << endl;
    }

    word solidType(is);

    word coeffs(is);

    if (coeffs == "defaultCoeffs")
    {
        ConstructorTable::iterator cstrIter =
            ConstructorTablePtr_->find(solidType);

        if (cstrIter == ConstructorTablePtr_->end())
        {
            FatalErrorIn("solid::New(Istream&)")
                << "Unknown solid type " << solidType << nl << nl
                << "Valid solid types are:" << endl
                << ConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<solid>(cstrIter()());
    }
    else if (coeffs == "coeffs")
    {
        IstreamConstructorTable::iterator cstrIter =
            IstreamConstructorTablePtr_->find(solidType);

        if (cstrIter == IstreamConstructorTablePtr_->end())
        {
            FatalErrorIn("solid::New(Istream&)")
                << "Unknown solid type " << solidType << nl << nl
                << "Valid solid types are:" << endl
                << IstreamConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<solid>(cstrIter()(is));
    }
    else
    {
        FatalErrorIn("solid::New(Istream&)")
            << "solid type " << solidType
            << ", option " << coeffs << " given"
            << ", should be coeffs or defaultCoeffs"
            << exit(FatalError);

        return autoPtr<solid>(NULL);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

#include "absorptionCoeffs.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::absorptionCoeffs::absorptionCoeffs(Istream& is)
:
   Tcommon_(readScalar(is)),
   Tlow_(readScalar(is)),
   Thigh_(readScalar(is)),
   invTemp_(readBool(is))
{
    for (label coefLabel=0; absorptionCoeffs::nCoeffs_; coefLabel++)
    {
        is >> highACoeffs_[coefLabel];
    }

    for (label coefLabel=0; absorptionCoeffs::nCoeffs_; coefLabel++)
    {
        is >> lowACoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::radiation::absorptionCoeffs::~absorptionCoeffs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::absorptionCoeffs::checkT(const scalar T) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        FatalErrorIn
        (
            "absorptionCoeffs::checkT(const scalar T) const"
        )   << "attempt to use absCoeff out of temperature range:" << nl
            << "    " << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << nl << abort(FatalError);
    }
}


const Foam::radiation::absorptionCoeffs::coeffArray&
Foam::radiation::absorptionCoeffs::coeffs
(
    const scalar T
) const
{
    checkT(T);

    if (T < Tcommon_)
    {
        return lowACoeffs_;
    }
    else
    {
        return highACoeffs_;
    }
}


void Foam::radiation::absorptionCoeffs::initialise(Istream&)
{
    absorptionCoeffs(Istream);
}


void Foam::radiation::absorptionCoeffs::initialise(const dictionary& dict)
{
    dict.lookup("Tcommon") >> Tcommon_;
    dict.lookup("Tlow") >> Tlow_;
    dict.lookup("Thigh") >> Thigh_;
    dict.lookup("invTemp") >> invTemp_;

    dict.lookup("loTcoeffs") >> lowACoeffs_;
    dict.lookup("hiTcoeffs") >> highACoeffs_;
}


// ************************************************************************* //

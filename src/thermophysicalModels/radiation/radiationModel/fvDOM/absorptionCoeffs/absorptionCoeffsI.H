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

inline bool Foam::radiation::absorptionCoeffs::invTemp() const
{
    return  invTemp_;
}


inline Foam::scalar Foam::radiation::absorptionCoeffs::Tcommon() const
{
    return  Tcommon_;
}


inline Foam::scalar Foam::radiation::absorptionCoeffs::Tlow() const
{
    return  Tlow_;
}


inline Foam::scalar Foam::radiation::absorptionCoeffs::Thigh() const
{
    return  Thigh_;
}


inline const Foam::radiation::absorptionCoeffs::coeffArray&
Foam::radiation::absorptionCoeffs::highACoeffs() const
{
    return  highACoeffs_;
}


inline const Foam::radiation::absorptionCoeffs::coeffArray&
Foam::radiation::absorptionCoeffs::lowACoeffs() const
{
    return  lowACoeffs_;
}


// ************************************************************************* //

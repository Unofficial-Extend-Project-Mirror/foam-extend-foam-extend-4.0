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

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const T Foam::SVD::sign(const T& a, const T& b)
{
    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::scalarRectangularMatrix& Foam::SVD::U() const
{
    return U_;
}

inline const Foam::scalarRectangularMatrix& Foam::SVD::V() const
{
    return V_;
}

inline const Foam::scalarDiagonalMatrix& Foam::SVD::S() const
{
    return S_;
}

inline const Foam::scalarRectangularMatrix& Foam::SVD::VSinvUt() const
{
    return VSinvUt_;
}

inline Foam::label Foam::SVD::nZeros() const
{
    return nZeros_;
}

inline Foam::scalar Foam::SVD::minNonZeroS() const
{
    scalar minS = S_[0];
    for(label i = 1; i < S_.size(); i++)
    {
        scalar s = S_[i];
        if (s > VSMALL && s < minS) minS = s;
    }
    return minS;
}


// ************************************************************************* //

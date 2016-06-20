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

Class
    FadOne

Description

\*---------------------------------------------------------------------------*/

#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int nVars>
inline Foam::FadOne<nVars>::FadOne()
{
    for (int i = 0; i < Length; i++)
    {
        data_[i] = 0;
    }
}


template<int nVars>
inline Foam::FadOne<nVars>::FadOne(const double& v)
{
    this->value() = v;

    for (int i = 0; i < nVars; i++)
    {
        this->deriv(i) = 0;
    }
}


template<int nVars>
inline Foam::FadOne<nVars>::FadOne(Istream& is)
{
    is >> *this;
}


template<int nVars>
inline Foam::FadOne<nVars>::FadOne(const FadOne<nVars>& f)
{
    for (int i = 0; i < Length; i++)
    {
        data_[i] = f.data_[i];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int nVars>
inline const double& Foam::FadOne<nVars>::value() const
{
    return data_[0];
}


template<int nVars>
inline double& Foam::FadOne<nVars>::value()
{
    return data_[0];
}


template<int nVars>
inline const double& Foam::FadOne<nVars>::deriv(const label i) const
{
#   ifdef FULLDEBUG
    if (i >= nVars)
    {
        FatalErrorIn
        (
            "double& FadOne<nVars>::deriv(const label i)"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    return data_[i + 1];
}


template<int nVars>
inline double& Foam::FadOne<nVars>::deriv(const label i)
{
#   ifdef FULLDEBUG
    if (i >= nVars)
    {
        FatalErrorIn
        (
            "double& FadOne<nVars>::deriv(const label i)"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    return data_[i + 1];
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<int nVars>
inline void Foam::FadOne<nVars>::operator=(const double& d)
{
    for (int i = 0; i < nVars; i++)
    {
        this->deriv(i) = 0;
    }

    this->value() = d;
}


template<int nVars>
inline void Foam::FadOne<nVars>::operator=(const FadOne<nVars>& v)
{
    for (int i = 0; i < Length; i++)
    {
        data_[i] = v.data_[i];
    }
}


template<int nVars>
inline void Foam::FadOne<nVars>::operator+=(const double& d)
{
    this->value() += d;
}


template<int nVars>
inline void Foam::FadOne<nVars>::operator+=(const FadOne<nVars>& v)
{
    for (int i = 0; i < Length; i++)
    {
        data_[i] += v.data_[i];
    }
}


template<int nVars>
inline void Foam::FadOne<nVars>::operator-=(const double& d)
{
    this->value() -= d;
}


template<int nVars>
inline void Foam::FadOne<nVars>::operator-=(const FadOne<nVars>& v)
{
    for (int i = 0; i < Length; i++)
    {
        data_[i] -= v.data_[i];
    }
}


template<int nVars>
inline void Foam::FadOne<nVars>::operator*=(const double& d)
{
    for (int i = 0; i < Length; i++)
    {
        data_[i] *= d;
    }
}


template<int nVars>
inline void Foam::FadOne<nVars>::operator*=(const FadOne<nVars>& v)
{
    for (int i = 0; i < nVars; i++)
    {
        this->deriv(i) = this->deriv(i)*v.value() + this->value()*v.deriv(i);
    }

    this->value() *= v.value();
}


template<int nVars>
inline void Foam::FadOne<nVars>::operator/=(const double& d)
{
    for (int i = 0; i < Length; i++)
    {
        data_[i] /= d;
    }
}


template<int nVars>
inline void Foam::FadOne<nVars>::operator/=(const FadOne<nVars>& v)
{
    double& lv = this->value();
    const double& vv = v.value();

    for (int i = 0; i < nVars; i++)
    {
        this->deriv(i) = (this->deriv(i)*vv - lv*v.deriv(i))/(vv*vv);
    }

    lv /= vv;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<int nVars>
inline Foam::FadOne<nVars> Foam::operator+
(
    const FadOne<nVars>& u,
    const FadOne<nVars>& v
)
{
    FadOne<nVars> r(u);
    r += v;
    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::operator+
(
    const FadOne<nVars>& u,
    const double v
)
{
    FadOne<nVars> r(u);
    r += v;
    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::operator+
(
    const double u,
    const FadOne<nVars>& v
)
{
    FadOne<nVars> r(v);
    r += u;
    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::operator-
(
    const FadOne<nVars>& u
)
{
    FadOne<nVars> r;
    for (int i = 0; i < FadOne<nVars>::Length; i++)
    {
        r.data_[i] = -u.data_[i];
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::operator-
(
    const FadOne<nVars>& u,
    const FadOne<nVars>& v
)
{
    FadOne<nVars> r(u);
    r -= v;
    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::operator-
(
    const FadOne<nVars>& u,
    const double v
)
{
    FadOne<nVars> r(u);
    r -= v;
    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::operator-
(
    const double u,
    const FadOne<nVars>& v
)
{
    FadOne<nVars> r(v);
    r -= u;
    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::operator*
(
    const FadOne<nVars>& u,
    const FadOne<nVars>& v
)
{
    FadOne<nVars> r(u);
    r *= v;
    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::operator*
(
    const FadOne<nVars>& u,
    const double v
)
{
    FadOne<nVars> r(u);
    r *= v;
    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::operator*
(
    const double u,
    const FadOne<nVars>& v
)
{
    FadOne<nVars> r(v);
    r *= u;
    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::operator/
(
    const FadOne<nVars>& u,
    const FadOne<nVars>& v
)
{
    FadOne<nVars> r(u);
    r /= v;
    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::operator/
(
    const FadOne<nVars>& u,
    const double v
)
{
    FadOne<nVars> r(u);
    r /= v;
    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::operator/
(
    const double u,
    const FadOne<nVars>& v
)
{
    FadOne<nVars> r(v);
    r /= u;
    return r;
}


template<int nVars>
inline bool Foam::operator<
(
    const FadOne<nVars>& u,
    const FadOne<nVars>& v
)
{
    return u.value() < v.value();
}


template<int nVars>
inline bool Foam::operator<=
(
    const FadOne<nVars>& u,
    const FadOne<nVars>& v
)
{
    return u.value() <= v.value();
}


template<int nVars>
inline bool Foam::operator>
(
    const FadOne<nVars>& u,
    const FadOne<nVars>& v
)
{
    return u.value() > v.value();
}


template<int nVars>
inline bool Foam::operator>=
(
    const FadOne<nVars>& u,
    const FadOne<nVars>& v
)
{
    return u.value() >= v.value();
}


template<int nVars>
inline bool Foam::operator==
(
    const FadOne<nVars>& u,
    const FadOne<nVars>& v
)
{
    return u.value() == v.value();
}


template<int nVars>
inline bool Foam::operator!=
(
    const FadOne<nVars>& u,
    const FadOne<nVars>& v
)
{
    return u.value() != v.value();
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


template<int nVars>
Foam::Istream& Foam::operator>>(Istream& is, FadOne<nVars>& f)
{
    // Read begining of FadOne
    if (token(is) != token::BEGIN_SQR)
    {
        Info<< "expected a " << token::BEGIN_SQR << " in FadOne<nVars>"
            << endl << "in stream " << is.info() << endl;
    }

    for (int i = 0; i < FadOne<nVars>::Length ; i++)
    {
        is >> f.data_[i];
    }

    // Check end of FadOne
    if (token(is) != token::END_SQR)
    {
        Info<< "expected a " << token::END_SQR << " in FadOne<nVars>"
            << endl << "in stream " << is.info() << endl;
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, FadOne<nVars>&)");

    return is;
}


template<int nVars>
Foam::Ostream& Foam::operator<<(Ostream& os, const FadOne<nVars>& f)
{
    os << token::BEGIN_SQR;

    for (int i = 0; i < FadOne<nVars>::Length - 1; i++)
    {
        os << f.data_[i] << token::SPACE;
    }
    os << f.data_[FadOne<nVars>::Length -1 ] << token::END_SQR;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const FadOne<nVars>&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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

Class
    FadOne

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<int nVars>
inline Foam::FadOne<nVars> Foam::mag(const FadOne<nVars>& u)
{
    FadOne<nVars> r(u.value());

    for (int i = 0; i < nVars; i++)
    {
        if (u.deriv(i) == 0)
        {
            r.deriv(i) = 0;
        }
        else
        {
            if (u.value() > 0)
            {
                r.deriv(i) = u.deriv(i);
            }
            else
            {
                r.deriv(i) = -u.deriv(i);
            }
        }
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::sin(const FadOne<nVars>& u)
{
    FadOne<nVars> r(Foam::sin(u.value()));

    for (int i = 0; i < nVars; i++)
    {
        r.deriv(i) = u.deriv(i)*Foam::cos(u.value());
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::asin(const FadOne<nVars>& u)
{
    const double& x = u.value();

    FadOne<nVars> r(Foam::asin(x));

    for (int i = 0; i < nVars; i++)
    {
        const double& d = u.deriv(i);
        r.deriv(i) = d/Foam::sqrt(1 - Foam::sqr(x));
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::cos(const FadOne<nVars>& u)
{
    FadOne<nVars> r(Foam::cos(u.value()));

    for (int i = 0; i < nVars; i++)
    {
        r.deriv(i) = -u.deriv(i)*Foam::sin(u.value());
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::acos(const FadOne<nVars>& u)
{
    const double& x = u.value();

    FadOne<nVars> r(Foam::acos(x));

    for (int i = 0; i < nVars; i++)
    {
        const double& d = u.deriv(i);
        r.deriv(i) = -d/Foam::sqrt(1 - Foam::sqr(x));
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::tan(const FadOne<nVars>& u)
{
    FadOne<nVars> r(Foam::tan(u.value()));

    for (int i = 0; i < nVars; i++)
    {
        r.deriv(i) = u.deriv(i)/(Foam::sqr(Foam::cos(u.value())));
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::atan(const FadOne<nVars>& u)
{
    const double& x = u.value();

    FadOne<nVars> r(Foam::atan(x));

    for (int i = 0; i < nVars; i++)
    {
        const double& d = u.deriv(i);
        r.deriv(i) = d/Foam::sqrt(1 + Foam::sqr(x));
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::tanh(const FadOne<nVars>& u)
{
    FadOne<nVars> r(Foam::tanh(u.value()));

    for (int i = 0; i < nVars; i++)
    {
        r.deriv(i) = u.deriv(i)*(1 - Foam::sqr(Foam::tanh(u.value())));
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::exp(const FadOne<nVars>& u)
{
    FadOne<nVars> r(Foam::exp(u.value()));

    for (int i = 0; i < nVars; i++)
    {
        r.deriv(i) = u.deriv(i)*Foam::exp(u.value());
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::log(const FadOne<nVars>& u)
{
    FadOne<nVars> r(Foam::log(u.value()));

    for (int i = 0; i < nVars; i++)
    {
        r.deriv(i) = u.deriv(i)/u.value();
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::sqr(const FadOne<nVars>& u)
{
    FadOne<nVars> r(Foam::sqr(u.value()));

    for (int i = 0; i < nVars; i++)
    {
        r.deriv(i) = 2*u.deriv(i)*u.value();
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::sqrt(const FadOne<nVars>& u)
{
    FadOne<nVars> r(Foam::sqrt(u.value()));

    for (int i = 0; i < nVars; i++)
    {
        if (u.deriv(i) > 0)
        {
            r.deriv(i) = 0.5*u.deriv(i)/Foam::sqrt(u.value());
        }
        else
        {
            r.deriv(i) = 0;
        }
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::pow
(
    const FadOne<nVars>& u,
    const FadOne<nVars>& v
)
{
    FadOne<nVars> r(Foam::pow(u.value(), v.value()));

    for (int i = 0; i < nVars; i++)
    {
        r.deriv(i) =
            v.deriv(i)*Foam::log(u.value())*Foam::pow(u.value(), v.value())
          + v.value()*u.deriv(i)*Foam::pow(u.value(), v.value() - 1);
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::pow
(
    const FadOne<nVars>& u,
    const double v
)
{
    FadOne<nVars> r(Foam::pow(u.value(), v));

    for (int i = 0; i < nVars; i++)
    {
        r.deriv(i) = v*u.deriv(i)*Foam::pow(u.value(), v - 1);
    }

    return r;
}


template<int nVars>
inline Foam::FadOne<nVars> Foam::pow
(
    const double& u,
    const FadOne<nVars>& v
)
{
    FadOne<nVars> r(Foam::pow(u, v.value()));

    for (int i = 0; i < nVars; i++)
    {
        r.deriv(i) = v.deriv(i)*Foam::log(u)*Foam::pow(u, v.value());
    }

    return r;
}


template<int nVars>
inline Foam::label Foam::sign(const FadOne<nVars>& u)
{
    return (u.value() >= 0) ? 1: -1;
}


template<int nVars>
inline Foam::label Foam::pos(const FadOne<nVars>& u)
{
    return (u.value() >= 0) ? 1: 0;
}


template<int nVars>
inline Foam::label Foam::neg(const FadOne<nVars>& u)
{
    return (u.value() < 0) ? 1: 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //

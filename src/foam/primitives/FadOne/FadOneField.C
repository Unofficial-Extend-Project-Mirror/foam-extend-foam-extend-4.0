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
    FadOneFields

\*----------------------------------------------------------------------------*/

#include "FadOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<int nVars>
Foam::tmp<Foam::Field<Foam::FadOne<nVars> > >
Foam::ValueFadOneField(const UList<scalar>& u)
{
    tmp<Foam::Field<FadOne<nVars> > > tr
    (
        new Foam::Field<FadOne<nVars> > (u.size())
    );

    Field<FadOne<nVars> >& r = tr();

    forAll (r, i)
    {
        r[i].value() = u[i];
    }

    return tr;
}


template<int nVars>
Foam::tmp<Foam::scalarField> Foam::FadOneValue(const Field<FadOne<nVars> >& u)
{
    tmp<scalarField> tr(new scalarField(u.size()));
    scalarField& r = tr();

    forAll (r, i)
    {
        r[i] = u[i].value();
    }

    return tr;
}


template<int nVars>
void Foam::FadOneSetValue
(
    Field<FadOne<nVars> >& u,
    const scalarField& val
)
{
    forAll (u, i)
    {
        u[i].value() = val[i];
    }
}


template<int nVars>
Foam::tmp<Foam::scalarField> Foam::FadOneDeriv
(
    const Field<FadOne<nVars> >& u,
    const direction d
)
{
    tmp<scalarField> tr(new scalarField(u.size()));
    scalarField& r = tr();

    forAll (r, i)
    {
        r[i] = u[i].deriv(d);
    }

    return tr;
}


template<int nVars>
void Foam::FadOneSetDeriv
(
    Field<FadOne<nVars> >& u,
    const direction d,
    const scalarField& der
)
{
    forAll (u, i)
    {
        u[i].deriv(d) = der[i];
    }
}


// ************************************************************************* //

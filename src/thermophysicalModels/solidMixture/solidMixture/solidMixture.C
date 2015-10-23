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

#include "solidMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::solidMixture::solidMixture
(
    const dictionary& thermophysicalProperties
)
:
    components_(thermophysicalProperties.lookup("solidComponents")),
    properties_(components_.size())
{

    forAll(components_, i)
    {
        properties_.set
        (
            i,
            solid::New(thermophysicalProperties.lookup(components_[i]))
        );
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidMixture> Foam::solidMixture::New
(
    const dictionary& thermophysicalProperties
)
{
    return autoPtr<solidMixture>(new solidMixture(thermophysicalProperties));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField Foam::solidMixture::X
(
    const scalarField& Y
) const
{
    scalarField X(Y.size());
    scalar rhoInv = 0.0;
    forAll(X, i)
    {
        rhoInv += Y[i]/properties_[i].rho();
        X[i] = Y[i]/properties_[i].rho();
    }

    return X/rhoInv;
}


Foam::scalar Foam::solidMixture::rho
(
    const scalarField& X
) const
{
    scalar tmp = 0.0;
    forAll(properties_, i)
    {
        tmp += properties_[i].rho()*X[i];
    }
    return tmp;
}


Foam::scalar Foam::solidMixture::cp
(
    const scalarField& Y
) const
{
    scalar tmp = 0.0;
    forAll(properties_, i)
    {
        tmp += properties_[i].cp()*Y[i];
    }
    return tmp;
}


// ************************************************************************* //

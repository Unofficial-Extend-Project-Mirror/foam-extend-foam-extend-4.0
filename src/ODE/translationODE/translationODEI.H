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
    translationODE

Description
    Ordinary differential equation for three degrees of freedom
    solid body translation

Author
    Hrvoje Jasak
    Dubravko Matijasevic

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar& Foam::translationODE::mass() const
{
    return mass_;
}


const Foam::dimensionedVector& Foam::translationODE::Xrel() const
{
    return Xrel_;
}


Foam::dimensionedVector Foam::translationODE::X() const
{
    return xEquilibrium_ + Xrel();
}


const Foam::dimensionedVector& Foam::translationODE::U() const
{
    return U_;
}


const Foam::dimensionedVector& Foam::translationODE::Uold() const
{
    return Uold_;
}


Foam::dimensionedVector Foam::translationODE::A() const
{
    return A(Xrel(), U());
}


Foam::dimensionedVector Foam::translationODE::Uaverage() const
{
    return 0.5*(U_ + Uold_);
}


const Foam::dimensionedVector& Foam::translationODE::force() const
{
    return force_;
}


Foam::dimensionedVector& Foam::translationODE::force()
{
    return force_;
}


// ************************************************************************* //

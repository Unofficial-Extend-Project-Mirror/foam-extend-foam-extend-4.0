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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar& Foam::sixDOFODE::mass() const
{
    return mass_;
}


Foam::dimensionedScalar& Foam::sixDOFODE::mass()
{
    return mass_;
}


const Foam::dimensionedDiagTensor& Foam::sixDOFODE::momentOfInertia() const
{
    return momentOfInertia_;
}


Foam::dimensionedDiagTensor& Foam::sixDOFODE::momentOfInertia()
{
    return momentOfInertia_;
}


const Foam::dimensionedVector& Foam::sixDOFODE::Xequilibrium() const
{
    return Xequilibrium_;
}


Foam::dimensionedVector& Foam::sixDOFODE::Xequilibrium()
{
    return Xequilibrium_;
}


const Foam::dimensionedVector& Foam::sixDOFODE::force() const
{
    return force_;
}


Foam::dimensionedVector& Foam::sixDOFODE::force()
{
    return force_;
}


const Foam::dimensionedVector& Foam::sixDOFODE::moment() const
{
    return moment_;
}


Foam::dimensionedVector& Foam::sixDOFODE::moment()
{
    return moment_;
}


const Foam::dimensionedVector& Foam::sixDOFODE::forceRelative() const
{
    return forceRelative_;
}


Foam::dimensionedVector& Foam::sixDOFODE::forceRelative()
{
    return forceRelative_;
}


const Foam::dimensionedVector& Foam::sixDOFODE::momentRelative() const
{
    return momentRelative_;
}


Foam::dimensionedVector& Foam::sixDOFODE::momentRelative()
{
    return momentRelative_;
}


Foam::dimensionedVector Foam::sixDOFODE::forceTotal() const
{
    return force_ + (this->toAbsolute() & forceRelative_);
}


Foam::dimensionedVector Foam::sixDOFODE::momentTotal() const
{
    return moment_ + (this->toAbsolute() & momentRelative_);
}


// ************************************************************************* //

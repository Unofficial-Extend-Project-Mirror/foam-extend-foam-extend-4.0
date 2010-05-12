/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

Description
    6-DOF solver using quaternions

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Updates by Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteRotation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar& Foam::sixDOFqODE::mass() const
{
    return mass_;
}


const Foam::dimensionedDiagTensor& Foam::sixDOFqODE::momentOfInertia() const
{
    return momentOfInertia_;
}


const Foam::dimensionedVector& Foam::sixDOFqODE::Xequilibrium() const
{
    return Xequilibrium_;
}


const Foam::dimensionedVector& Foam::sixDOFqODE::Xrel() const
{
    return Xrel_;
}


const Foam::dimensionedVector& Foam::sixDOFqODE::omega() const
{
    return omega_;
}


Foam::dimensionedVector Foam::sixDOFqODE::omegaDot() const
{
    return OmegaDot(rotation_.eCurrent(), omega_);
}


Foam::dimensionedVector Foam::sixDOFqODE::X() const
{
    return Xequilibrium_ + Xrel_;
}


const Foam::dimensionedVector& Foam::sixDOFqODE::U() const
{
    return U_;
}


Foam::dimensionedVector Foam::sixDOFqODE::A() const
{
    return A(Xrel_, U_, rotation_.eCurrent());
}


Foam::vector Foam::sixDOFqODE::rotVector() const
{
    return rotation_.rotVector();
}


Foam::dimensionedScalar Foam::sixDOFqODE::rotAngle() const
{
    return dimensionedScalar("rotAngle", dimless, rotation_.rotAngle());
}


Foam::dimensionedVector& Foam::sixDOFqODE::U()
{
    return U_;
}


Foam::dimensionedVector& Foam::sixDOFqODE::omega()
{
    return omega_;
}


const Foam::dimensionedVector& Foam::sixDOFqODE::Uaverage() const
{
    return Uaverage_;
}


Foam::vector Foam::sixDOFqODE::rotVectorAverage() const
{
    return rotation_.rotVectorAverage();
}


const Foam::dimensionedVector& Foam::sixDOFqODE::omegaAverage() const
{
    return omegaAverage_;
}


const Foam::dimensionedVector& Foam::sixDOFqODE::omegaAverageAbsolute() const
{
    return omegaAverageAbsolute_;
}


const Foam::dimensionedVector& Foam::sixDOFqODE::force() const
{
    return force_;
}


Foam::dimensionedVector& Foam::sixDOFqODE::force()
{
    return force_;
}


const Foam::dimensionedVector& Foam::sixDOFqODE::moment() const
{
    return moment_;
}


Foam::dimensionedVector& Foam::sixDOFqODE::moment()
{
    return moment_;
}


const Foam::dimensionedVector& Foam::sixDOFqODE::forceRelative() const
{
    return forceRelative_;
}


Foam::dimensionedVector& Foam::sixDOFqODE::forceRelative()
{
    return forceRelative_;
}


const Foam::dimensionedVector& Foam::sixDOFqODE::momentRelative() const
{
    return momentRelative_;
}


Foam::dimensionedVector& Foam::sixDOFqODE::momentRelative()
{
    return momentRelative_;
}


Foam::dimensionedVector Foam::sixDOFqODE::forceTotal() const
{
    return force() + (toAbsolute() & forceRelative());
}


Foam::dimensionedVector Foam::sixDOFqODE::momentTotal() const
{
    return moment() + (toAbsolute() & momentRelative());
}


Foam::tensor Foam::sixDOFqODE::toRelative() const
{
    return rotation_.eCurrent().R();
}


Foam::tensor Foam::sixDOFqODE::toAbsolute() const
{
    return rotation_.eCurrent().invR();
}


const Foam::tensor& Foam::sixDOFqODE::rotIncrementTensor() const
{
    return rotation_.rotIncrementTensor();
}


// ************************************************************************* //

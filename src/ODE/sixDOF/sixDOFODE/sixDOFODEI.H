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
    sixDOFODE

Description
    Abstract base class for six-degrees-of-freedom (6DOF) ordinary differential
    equations with necessary interface for two-way coupling with CFD solvers.

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Hrvoje Jasak, FSB Zagreb. All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

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


const Foam::dimensionedDiagTensor& Foam::sixDOFODE::linSpringCoeffs() const
{
    return linSpringCoeffs_;
}


const Foam::dimensionedDiagTensor& Foam::sixDOFODE::linDampingCoeffs() const
{
    return linDampingCoeffs_;
}


const Foam::dimensionedVector& Foam::sixDOFODE::force() const
{
    return force_;
}


const Foam::dimensionedVector& Foam::sixDOFODE::moment() const
{
    return moment_;
}


void Foam::sixDOFODE::setExternalForceAndMoment
(
    const vector& externalForce,
    const vector& externalMoment
)
{
    // Initialise the state before setting external forces and moments
    initState();

    // Set forces and moments
    force_ = externalForce;
    moment_ = externalMoment;

    // Relax acceleration if the Aitkens relaxation is turned on
    if (aitkensRelaxation_)
    {
        relaxAcceleration();
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

Foam::scalar Foam::sixDOFODE::odeStepFraction(const scalar odeTime) const
{
    // Get current global time. Note: assuming that the ODESolver for a given
    // time step is integrating from t - deltaT to t.
    const scalar globalTime = dict().time().value();
    const scalar globalDeltaT = dict().time().deltaT().value();

    // Calculate current ODE step time step size
    const scalar odeDeltaT = globalDeltaT - (globalTime - odeTime);

    // Calculate and return ODE step size fraction
    return odeDeltaT/globalDeltaT;
}


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


const Foam::dimensionedVector& Foam::sixDOFODE::moment() const
{
    return moment_;
}


void Foam::sixDOFODE::setExternalForceAndMoment
(
    const dimensionedVector& externalForce,
    const dimensionedVector& externalMoment
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


void Foam::sixDOFODE::initExternalForceAndMoment
(
    const dimensionedVector& externalForce,
    const dimensionedVector& externalMoment
)
{
    // Initialise force and moment only for the first time step
    if (curTimeIndex_ == -1)
    {
        force_ = externalForce;
        moment_ = externalMoment;
    }
}


const Foam::PtrList<Foam::translationalConstraint>&
Foam::sixDOFODE::translationalConstraints() const
{
    return translationalConstraints_;
}


const Foam::PtrList<Foam::rotationalConstraint>&
Foam::sixDOFODE::rotationalConstraints() const
{
    return rotationalConstraints_;
}


const Foam::PtrList<Foam::translationalRestraint>&
Foam::sixDOFODE::translationalRestraints() const
{
    return translationalRestraints_;
}


const Foam::PtrList<Foam::rotationalRestraint>&
Foam::sixDOFODE::rotationalRestraints() const
{
    return rotationalRestraints_;
}

// ************************************************************************* //

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


Foam::dimensionedScalar& Foam::sixDOFqODE::mass()
{
    return mass_;
}


const Foam::dimensionedDiagTensor& Foam::sixDOFqODE::momentOfInertia() const
{
    return momentOfInertia_;
}


Foam::dimensionedDiagTensor& Foam::sixDOFqODE::momentOfInertia()
{
    return momentOfInertia_;
}


const Foam::dimensionedVector& Foam::sixDOFqODE::Xequilibrium() const
{
    return Xequilibrium_;
}


Foam::dimensionedVector& Foam::sixDOFqODE::Xequilibrium()
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


void Foam::sixDOFqODE::setXrel(const vector& x)
{
    Xrel_.value() = x;

    // Set X in coefficients
    coeffs_[0] = x.x();
    coeffs_[1] = x.y();
    coeffs_[2] = x.z();
}


void Foam::sixDOFqODE::setU(const vector& U)
{
    U_.value() = U;
    Uaverage_ = U_;

    // Set U in coefficients
    coeffs_[3] = U.x();
    coeffs_[4] = U.y();
    coeffs_[5] = U.z();
}


void Foam::sixDOFqODE::setRotation(const HamiltonRodriguezRot& rot)
{
    // Cannot call update rotation: simply set up coefficients
    // HJ, 23/Mar/2015

    coeffs_[9] = rot.e0();
    coeffs_[10] = rot.e1();
    coeffs_[11] = rot.e2();
    coeffs_[12] = rot.e3();
}


void Foam::sixDOFqODE::setOmega(const vector& omega)
{
    omega_.value() = omega;
    omegaAverage_ = omega_;
    omegaAverageAbsolute_ = omega_;

    // Set omega in coefficients
    coeffs_[6] = omega.x();
    coeffs_[7] = omega.y();
    coeffs_[8] = omega.z();
}


void Foam::sixDOFqODE::setReferentRotation(const HamiltonRodriguezRot& rot)
{
    referentRotation_ = rot;

    referentMotionConstraints_ = true;
}


const Foam::dimensionedVector& Foam::sixDOFqODE::Uaverage() const
{
    return Uaverage_;
}


const Foam::finiteRotation& Foam::sixDOFqODE::rotation() const
{
    return rotation_;
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

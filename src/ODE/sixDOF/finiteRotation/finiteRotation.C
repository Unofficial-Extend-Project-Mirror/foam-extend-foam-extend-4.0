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

Class
    finiteRotation

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Update by Hrvoje Jasak

\*----------------------------------------------------------------------------*/

#include "finiteRotation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::finiteRotation::rotVector(const tensor& rotT)
{
    vector ur = - *( inv(I + rotT) & (I - rotT) );

    // Scaling to a unit vector.  HJ, problems with round-off
    // HJ, 4/Aug/2008
    return ur/(mag(ur) + SMALL);
}


Foam::scalar Foam::finiteRotation::rotAngle(const tensor& rotT)
{
    // Alternative formulation: Daniel Schmode, 15/Feb/2009
    scalar x = rotT.zy() - rotT.yz();
    scalar y = rotT.xz() - rotT.zx();
    scalar z = rotT.yx() - rotT.xy();

    scalar r = hypot(x, hypot(y, z));
    scalar t = tr(rotT);

    return atan2(r, t - 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteRotation::finiteRotation(const HamiltonRodriguezRot& rot)
:
    eInitial_(rot),
    rotTensor_(rot.R()),
    rotIncrementTensor_(tensor::zero)
{}


Foam::finiteRotation::finiteRotation
(
    const vector& r,
    const scalar& angle
)
:
    eInitial_(HamiltonRodriguezRot(r, angle)),
    rotTensor_(eInitial_.R()),
    rotIncrementTensor_(tensor::zero)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::finiteRotation::~finiteRotation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::finiteRotation::updateRotation(const HamiltonRodriguezRot& rot)
{
    tensor rotR = rot.R();

    rotIncrementTensor_ = (rotR & (rotTensor_.T()));
    rotTensor_ = rotR;
}


const Foam::HamiltonRodriguezRot& Foam::finiteRotation::eInitial() const
{
    return eInitial_;
}


Foam::HamiltonRodriguezRot Foam::finiteRotation::eCurrent() const
{
    return HamiltonRodriguezRot(rotVector(), rotAngle());
}


const Foam::tensor& Foam::finiteRotation::rotTensor() const
{
    return rotTensor_;
}


Foam::vector Foam::finiteRotation::rotVector() const
{
    return rotVector(rotTensor_);
}


Foam::scalar Foam::finiteRotation::rotAngle() const
{
    return rotAngle(rotTensor_);
}


const Foam::tensor& Foam::finiteRotation::rotIncrementTensor() const
{
    return rotIncrementTensor_;
}


Foam::vector Foam::finiteRotation::rotIncrementVector() const
{
    return rotVector(rotIncrementTensor_);
}


Foam::scalar Foam::finiteRotation::rotIncrementAngle() const
{
    return rotAngle(rotIncrementTensor_);
}


Foam::vector Foam::finiteRotation::omegaAverage(const scalar deltaT) const
{
    return (rotIncrementAngle()/deltaT)*rotIncrementVector();
}


Foam::tensor Foam::finiteRotation::rotTensorAverage() const
{
    return
    (
        HamiltonRodriguezRot
        (
            rotIncrementVector(), 
            0.5*rotIncrementAngle()
        ).R().T()
      & rotTensor_
    );
}


Foam::vector Foam::finiteRotation::rotVectorAverage() const
{
    return rotVector(rotTensorAverage());
}


Foam::vector Foam::finiteRotation::omegaAverageAbsolute
(
    const scalar deltaT
) const
{
    return (rotTensorAverage().T() & omegaAverage(deltaT));
}


// ************************************************************************* //

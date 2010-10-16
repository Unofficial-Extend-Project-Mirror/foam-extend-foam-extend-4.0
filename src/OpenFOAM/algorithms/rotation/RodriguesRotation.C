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

\*---------------------------------------------------------------------------*/

#include "RodriguesRotation.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tensor Foam::RodriguesRotation
(
    const vector& rotationAxis,
    const scalar& rotationAngle,
    const bool inDegrees
)
{
    tensor rotTensor;
    scalar theta = rotationAngle;

    if (inDegrees)
    {
        theta *= mathematicalConstant::pi/180.0;
    }

    scalar sinTheta = sin(theta);
    scalar cosTheta = cos(theta);
    scalar oneMinusCosTheta = 1.0 - cosTheta;

    scalar magRotAxis = mag(rotationAxis);

    if (magRotAxis < SMALL)
    {
        FatalErrorIn
        (
            "tensor RodriguesRotation\n"
            "(\n"
            "    const vector& rotationAxis,\n"
            "    const scalar& rotationAngle\n"
            ")"
        )   << "Incorrectly defined axis: " << rotationAxis
            << abort(FatalError);
    }

    vector unitVector = rotationAxis/magRotAxis;

    scalar wx = unitVector.x();
    scalar wy = unitVector.y();
    scalar wz = unitVector.z();

    rotTensor.xx() = cosTheta + sqr(wx)*oneMinusCosTheta;
    rotTensor.yy() = cosTheta + sqr(wy)*oneMinusCosTheta;
    rotTensor.zz() = cosTheta + sqr(wz)*oneMinusCosTheta;

    rotTensor.xy() = wx*wy*oneMinusCosTheta - wz*sinTheta;
    rotTensor.xz() = wy*sinTheta + wx*wz*oneMinusCosTheta;

    rotTensor.yx() =  wz*sinTheta + wx*wy*oneMinusCosTheta;
    rotTensor.yz() = -wx*sinTheta + wy*wz*oneMinusCosTheta;

    rotTensor.zx() = -wy*sinTheta + wx*wz*oneMinusCosTheta;
    rotTensor.zy() =  wx*sinTheta + wy*wz*oneMinusCosTheta;

    return rotTensor;
}


// ************************************************************************* //

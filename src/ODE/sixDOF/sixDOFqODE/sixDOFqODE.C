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
    sixDOFqODE

Description
    6-DOF solver using quaternions

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "sixDOFqODE.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Runtime type information
// Not possible because of I/O error: incorrect type, expecting dictionary
// HJ, 11/Feb/2008
// namespace Foam
// {
//     defineTypeNameAndDebug(sixDOFqODE, 0);
// }


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sixDOFqODE::setCoeffs()
{
    // Set ODE coefficients from position and rotation

    // Linear displacement relative to spring equilibrium
    const vector& Xval = Xrel_.value();
    coeffs_[0] = Xval.x();
    coeffs_[1] = Xval.y();
    coeffs_[2] = Xval.z();

    // Linear velocity
    const vector& Uval = U_.value();
    coeffs_[3] = Uval.x();
    coeffs_[4] = Uval.y();
    coeffs_[5] = Uval.z();

    // Rotational velocity in non - inertial coordinate system
    const vector& omegaVal = omega_.value();
    coeffs_[6] = omegaVal.x();
    coeffs_[7] = omegaVal.y();
    coeffs_[8] = omegaVal.z();

    // Quaternions
    coeffs_[9] = rotation_.eInitial().e0();
    coeffs_[10] = rotation_.eInitial().e1();
    coeffs_[11] = rotation_.eInitial().e2();
    coeffs_[12] = rotation_.eInitial().e3();
}


Foam::dimensionedVector Foam::sixDOFqODE::A
(
    const dimensionedVector& xR,
    const dimensionedVector& uR,
    const HamiltonRodriguezRot& rotation
) const
{
    return
    (
       - (linSpringCoeffs_ & xR)    // spring
       - (linDampingCoeffs_ & uR)   // damping
       + force()
         // To absolute
       + (rotation.invR() & forceRelative())
    )/mass_;
}


Foam::dimensionedVector Foam::sixDOFqODE::OmegaDot
(
    const HamiltonRodriguezRot& rotation,
    const dimensionedVector& omega
) const
{
    return 
        inv(momentOfInertia_)
      & (
            E(omega)
            // To relative
          + (rotation.R() & moment())
          + momentRelative()
        );
}


Foam::dimensionedVector Foam::sixDOFqODE::E
(
    const dimensionedVector& omega
) const
{
    return (*(momentOfInertia_ & omega) & omega);
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::sixDOFqODE::sixDOFqODE(const IOobject& io)
:
    IOdictionary(io),
    mass_(lookup("mass")),
    momentOfInertia_(lookup("momentOfInertia")),
    Xequilibrium_(lookup("equilibriumPosition")),
    linSpringCoeffs_(lookup("linearSpring")),
    linDampingCoeffs_(lookup("linearDamping")),
    Xrel_(lookup("Xrel")),
    U_(lookup("U")),
    Uaverage_(U_),
    rotation_
    (
        vector(lookup("rotationVector")),
        dimensionedScalar(lookup("rotationAngle")).value()
    ),
    omega_(lookup("omega")),
    omegaAverage_(omega_),
    omegaAverageAbsolute_(omega_),
    force_(lookup("force")),
    moment_(lookup("moment")),
    forceRelative_(lookup("forceRelative")),
    momentRelative_(lookup("momentRelative")),
    coeffs_(13, 0.0)
{
    setCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDOFqODE::~sixDOFqODE()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDOFqODE::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    // Set the derivatives for displacement
    dydx[0] = y[3];
    dydx[1] = y[4];
    dydx[2] = y[5];

    dimensionedVector curX("curX", dimLength, vector(y[0], y[1], y[2]));
    dimensionedVector curU("curU", dimVelocity, vector(y[3], y[4], y[5]));
    const HamiltonRodriguezRot curRotation
    (
        y[9],
        y[10],
        y[11],
        y[12]
    );

    const vector accel = A(curX, curU, curRotation).value();

    dydx[3] = accel.x();
    dydx[4] = accel.y();
    dydx[5] = accel.z();

    // Set the derivatives for rotation
    dimensionedVector curOmega
    (
        "curOmega",
        dimless/dimTime,
        vector(y[6], y[7], y[8])
    );

    const vector omegaDot = OmegaDot(curRotation, curOmega).value();

    dydx[6]  = omegaDot.x();
    dydx[7] = omegaDot.y();
    dydx[8] = omegaDot.z();

    dydx[9] = curRotation.eDot(curOmega.value(), 0);
    dydx[10] = curRotation.eDot(curOmega.value(), 1);
    dydx[11] = curRotation.eDot(curOmega.value(), 2);
    dydx[12] = curRotation.eDot(curOmega.value(), 3);
}


void Foam::sixDOFqODE::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    Info << "jacobian(...)" << endl;
    notImplemented("sixDOFqODE::jacobian(...) const");
}


void Foam::sixDOFqODE::update(const scalar delta)
{
    // Update position
    vector Xold = Xrel_.value();

    vector& Xval = Xrel_.value();

    Xval.x() = coeffs_[0];
    Xval.y() = coeffs_[1];
    Xval.z() = coeffs_[2];

    // Update velocity
    Uaverage_.value() = (Xval - Xold)/delta;

    vector& Uval = U_.value();

    Uval.x() = coeffs_[3];
    Uval.y() = coeffs_[4];
    Uval.z() = coeffs_[5];

    // Update omega
    vector& omegaVal = omega_.value();

    omegaVal.x() = coeffs_[6];
    omegaVal.y() = coeffs_[7];
    omegaVal.z() = coeffs_[8];

    rotation_.updateRotation
    (
        HamiltonRodriguezRot
        (
            coeffs_[9],
            coeffs_[10],
            coeffs_[11],
            coeffs_[12]
        )
    );

    omegaAverage_.value() = rotation_.omegaAverage(delta);
    omegaAverageAbsolute_.value() = rotation_.omegaAverageAbsolute(delta);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::sixDOFqODE::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const sixDOFqODE& sds)
{
    os.writeKeyword("mass") << tab << sds.mass_ << token::END_STATEMENT << nl;
    os.writeKeyword("momentOfInertia") << tab << sds.momentOfInertia_
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("equilibriumPosition") << tab << sds.Xequilibrium_
        << token::END_STATEMENT << nl;
    os.writeKeyword("linearSpring") << tab << sds.linSpringCoeffs_
        << token::END_STATEMENT << nl;
    os.writeKeyword("linearDamping") << tab << sds.linDampingCoeffs_
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("Xrel") << tab << sds.Xrel() << token::END_STATEMENT << nl;
    os.writeKeyword("U") << tab << sds.U() << token::END_STATEMENT << nl;
    os.writeKeyword("rotationVector") << tab << sds.rotVector()
        << token::END_STATEMENT << nl;
    os.writeKeyword("rotationAngle") << tab << sds.rotAngle()
        << token::END_STATEMENT << nl;
    os.writeKeyword("omega") << tab << sds.omega()
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("force") << tab << sds.force()
        << token::END_STATEMENT << nl;
    os.writeKeyword("moment") << tab << sds.moment()
        << token::END_STATEMENT << nl;
    os.writeKeyword("forceRelative") << tab << sds.forceRelative()
        << token::END_STATEMENT << nl;
    os.writeKeyword("momentRelative") << tab << sds.momentRelative()
        << token::END_STATEMENT << endl;

    return os;
}


// ************************************************************************* //

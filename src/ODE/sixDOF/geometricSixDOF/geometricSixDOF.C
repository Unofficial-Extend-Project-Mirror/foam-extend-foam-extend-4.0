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
    geometricSixDOF

Description
    6-DOF solver using a geometric method for integration of rotations.

Author
    Viktor Pandza, FSB Zagreb.  All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    geometricSixDOF.C

\*---------------------------------------------------------------------------*/

#include "geometricSixDOF.H"
#include "OutputControlDictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(geometricSixDOF, 0);
addToRunTimeSelectionTable(sixDOFODE, geometricSixDOF, dictionary);

}


const Foam::debug::tolerancesSwitch Foam::geometricSixDOF::rotIncTensorTol_
(
    "geometrixSixDOFRotIncTensorTol",
    1e-10
);


const Foam::debug::tolerancesSwitch Foam::geometricSixDOF::rotIncRateTol_
(
    "geometrixSixDOFRotIncRateTol",
    1e-6
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dimensionedVector Foam::geometricSixDOF::A
(
    const dimensionedVector& xR,
    const dimensionedVector& uR,
    const tensor& R
) const
{
    // Fix the total force in global coordinate system
    dimensionedVector fAbs =
        // Force in global coordinate system
        force()
        // Force in local coordinate system
      + (R.T() & forceRelative())
        // Spring force in global coordinate system
      - (linSpringCoeffs() & xR)
        // Damping force in global coordinate system
      - (linDampingCoeffs() & uR);

    // Constrain translation simply by setting the total force to zero
    constrainTranslation(fAbs.value());

    return fAbs/mass();
}


Foam::dimensionedVector Foam::geometricSixDOF::OmegaDot
(
    const tensor& R,
    const dimensionedVector& omega
) const
{
    // External moment (torque) in local coordinate system
    dimensionedVector mRel =
        // Moment in global coordinate system
        (dimensionedTensor("R", dimless, R) & moment())
        // Moment in local coordinate system
      + momentRelative();

    // Note: constraints not implemented at the moment. They shall be
    // implemented in terms of Lagrange multipliers.

    return
        inv(momentOfInertia())
      & (
            E(omega)
          + mRel
        );
}


Foam::dimensionedVector Foam::geometricSixDOF::E
(
    const dimensionedVector& omega
) const
{
    return (*(momentOfInertia() & omega) & omega);
}


void Foam::geometricSixDOF::constrainTranslation(vector& vec) const
{
    // Constrain the vector with respect to referent or global coordinate system
    if (referentMotionConstraints_)
    {
        vector consVec(referentRotation_.R() & vec);

        if (fixedSurge_)
        {
            consVec.x() = 0;
        }
        if (fixedSway_)
        {
            consVec.y() = 0;
        }
        if (fixedHeave_)
        {
            consVec.z() = 0;
        }

        vec = referentRotation_.invR() & consVec;
    }
    else
    {
        if (fixedSurge_)
        {
            vec.x() = 0;
        }
        if (fixedSway_)
        {
            vec.y() = 0;
        }
        if (fixedHeave_)
        {
            vec.z() = 0;
        }
    }
}


Foam::tensor Foam::geometricSixDOF::expMap(const vector& rotInc) const
{
    tensor R;

    // Calculate the magnitude of the rotational increment vector
    const scalar magRotInc = mag(rotInc);

    if (magRotInc < rotIncTensorTol_)
    {
        // No rotational increment
        R = I;
    }
    else
    {
        // Calculate rotational increment tensor using the exponential map

        // Skew-symmetric tensor corresponding to the rotation increment
        const tensor skewRotInc(*rotInc);

        R = I
          + skewRotInc*sin(magRotInc)/magRotInc
          + (skewRotInc & skewRotInc)*(1.0 - cos(magRotInc))/sqr(magRotInc);
    }

    return R.T();
}


Foam::vector Foam::geometricSixDOF::dexpMap
(
    const vector& rotInc,
    const vector& omega
) const
{
    vector rotIncDot;

    // Calculate the magnitude of the rotational increment vector
    const scalar magRotInc = mag(rotInc);

    if (magRotInc < rotIncRateTol_)
    {
        // Stabilised calculation of rotation increment to avoid small
        // denominators
        rotIncDot = omega;
    }
    else
    {
        // Calculate rate of the rotational increment vector using the
        // differential of the exponential map

        // Skew-symmetric tensor corresponding to the rotation increment
        const tensor skewRotInc(*rotInc);

        rotIncDot =
        (
            I
          + 0.5*skewRotInc
          - (skewRotInc & skewRotInc)*
            (magRotInc/tan(magRotInc/2.0) - 2.0)/(2.0*sqr(magRotInc))
        )
      & omega;
    }

    return rotIncDot;
}


// * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * * * //

void Foam::geometricSixDOF::setState(const sixDOFODE& sd)
{
    // First set the state in base class subobject
    sixDOFODE::setState(sd);

    // Cast sixDOFODE& to geometricSixDOF&
    const geometricSixDOF& gsd = refCast<const geometricSixDOF>(sd);

    // Set state variables for this class
    Xrel_ = gsd.Xrel_;
    U_ = gsd.U_;
    Uaverage_ = gsd.Uaverage_;
    rotation_ = gsd.rotation_;
    omega_ = gsd.omega_;
    omegaAverage_ = gsd.omegaAverage_;

    // Copy ODE coefficients: this carries actual ODE state
    // HJ, 23/Mar/2015
    coeffs_ = gsd.coeffs_;

    fixedSurge_ = gsd.fixedSurge_;
    fixedSway_ = gsd.fixedSway_;
    fixedHeave_ = gsd.fixedHeave_;
    fixedRoll_ = gsd.fixedRoll_;
    fixedPitch_ = gsd.fixedPitch_;
    fixedYaw_ = gsd.fixedYaw_;
    referentMotionConstraints_ = gsd.referentMotionConstraints_;
    referentRotation_ = gsd.referentRotation_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometricSixDOF::geometricSixDOF(const IOobject& io)
:
    sixDOFODE(io),

    Xrel_(dict().lookup("Xrel")),
    U_(dict().lookup("U")),
    Uaverage_("Uaverage", U_),
    rotation_(tensor(dict().lookup("rotationTensor"))),
    rotIncrement_
    (
        dict().lookupOrDefault<tensor>("rotationIncrementTensor", tensor::zero)
    ),
    omega_(dict().lookup("omega")),
    omegaAverage_("omegaAverage", omega_),

    nEqns_(),
    coeffs_(),

    fixedSurge_(dict().lookup("fixedSurge")),
    fixedSway_(dict().lookup("fixedSway")),
    fixedHeave_(dict().lookup("fixedHeave")),
    fixedRoll_(dict().lookup("fixedRoll")),
    fixedPitch_(dict().lookup("fixedPitch")),
    fixedYaw_(dict().lookup("fixedYaw")),
    referentMotionConstraints_
    (
        dict().lookupOrDefault<Switch>
        (
            "referentMotionConstraints",
            false
        )
    ),
    referentRotation_(vector(1, 0, 0), 0)
{
    // Missing constraints. Count how many rotational constraints we have in
    // order to count the number of equations.

    nEqns_ = 12;

    // Set size for ODE coefficients depending on number of equations
    coeffs_.setSize(nEqns_);

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

    // Increment of the rotation vector (zero for initial condition)
    coeffs_[9] = 0;
    coeffs_[10] = 0;
    coeffs_[11] = 0;
}


Foam::geometricSixDOF::geometricSixDOF
(
    const word& name,
    const geometricSixDOF& gsd
)
:
    sixDOFODE(name, gsd),

    Xrel_(gsd.Xrel_.name(), gsd.Xrel_),
    U_(gsd.U_.name(), gsd.U_),
    Uaverage_(gsd.Uaverage_.name(), gsd.Uaverage_),
    rotation_(gsd.rotation_),
    omega_(gsd.omega_.name(), gsd.omega_),
    omegaAverage_(gsd.omegaAverage_.name(), gsd.omegaAverage_),

    nEqns_(gsd.nEqns_),
    coeffs_(gsd.coeffs_),

    fixedSurge_(gsd.fixedSurge_),
    fixedSway_(gsd.fixedSway_),
    fixedHeave_(gsd.fixedHeave_),
    fixedRoll_(gsd.fixedRoll_),
    fixedPitch_(gsd.fixedPitch_),
    fixedYaw_(gsd.fixedYaw_),
    referentMotionConstraints_(gsd.referentMotionConstraints_),
    referentRotation_(gsd.referentRotation_)
{}


Foam::autoPtr<Foam::sixDOFODE> Foam::geometricSixDOF::clone
(
    const word& name
) const
{
    // Create and return an autoPtr to the new geometricSixDOF object with a
    // new name
    return autoPtr<sixDOFODE>
    (
        new geometricSixDOF
        (
            name,
            *this
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::geometricSixDOF::~geometricSixDOF()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedVector& Foam::geometricSixDOF::Xrel() const
{
    return Xrel_;
}


const Foam::dimensionedVector& Foam::geometricSixDOF::omega() const
{
    return omega_;
}


Foam::dimensionedVector Foam::geometricSixDOF::X() const
{
    return Xequilibrium() + Xrel_;
}


const Foam::dimensionedVector& Foam::geometricSixDOF::U() const
{
    return U_;
}

const Foam::dimensionedVector& Foam::geometricSixDOF::Uaverage() const
{
    return Uaverage_;
}


const Foam::dimensionedVector& Foam::geometricSixDOF::omegaAverage() const
{
    return omegaAverage_;
}


Foam::tensor Foam::geometricSixDOF::toRelative() const
{
    return rotation_;
}


Foam::tensor Foam::geometricSixDOF::toAbsolute() const
{
    return rotation_.T();
}


const Foam::tensor& Foam::geometricSixDOF::rotIncrementTensor() const
{
    return rotIncrement_;
}


void Foam::geometricSixDOF::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    // Translation

    // Set the derivatives for displacement
    dydx[0] = y[3];
    dydx[1] = y[4];
    dydx[2] = y[5];

    dimensionedVector curX("curX", dimLength, vector(y[0], y[1], y[2]));
    dimensionedVector curU("curU", dimVelocity, vector(y[3], y[4], y[5]));

    // Get rotational increment vector (u)
    const vector rotIncrementVector(y[9], y[10], y[11]);

    // Calculate current rotation tensor obtained with exponential map
    const tensor curRot = (expMap(rotIncrementVector) & rotation_);

    // Calculate translational acceleration using current rotation
    const vector accel = A(curX, curU, curRot).value();

    // Set the derivatives for velocity
    dydx[3] = accel.x();
    dydx[4] = accel.y();
    dydx[5] = accel.z();

    // Rotation

    dimensionedVector curOmega
    (
        "curOmega",
        dimless/dimTime,
        vector(y[6], y[7], y[8])
    );

    // Calculate rotational acceleration using current rotation
    const vector omegaDot = OmegaDot(curRot, curOmega).value();

    dydx[6] = omegaDot.x();
    dydx[7] = omegaDot.y();
    dydx[8] = omegaDot.z();

    // Calculate derivative of rotIncrementVector using current rotation
    // information
    const vector rotIncrementVectorDot =
        dexpMap(rotIncrementVector, curOmega.value());

    // Set the derivatives for rotation
    dydx[9] = rotIncrementVectorDot.x();
    dydx[10] = rotIncrementVectorDot.y();
    dydx[11] = rotIncrementVectorDot.z();
}


void Foam::geometricSixDOF::update(const scalar delta)
{
    // Translation

    // Update displacement
    const vector Xold = Xrel_.value();

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

    // Constrain velocity and re-set coefficients
    constrainTranslation(Uval);
    coeffs_[3] = Uval.x();
    coeffs_[4] = Uval.y();
    coeffs_[5] = Uval.z();

    // Rotation

    // Update omega
    const vector omegaOld = omega_.value();

    vector& omegaVal = omega_.value();

    omegaVal.x() = coeffs_[6];
    omegaVal.y() = coeffs_[7];
    omegaVal.z() = coeffs_[8];

    // Update rotational increment tensor
    rotIncrement_ = expMap(vector(coeffs_[9], coeffs_[10], coeffs_[11]));

    // Update rotational tensor
    rotation_ = (rotIncrement_ & rotation_);

    // Reset increment vector in coefficients for the next step
    coeffs_[9] = 0;
    coeffs_[10] = 0;
    coeffs_[11] = 0;

    // Consider calculating average omega using rotational tensor and rotational
    // increment tensors
    omegaAverage_.value() = 0.5*(omegaVal + omegaOld);
}


bool Foam::geometricSixDOF::writeData(Ostream& os) const
{
    // First write the part related to base class subobject
    sixDOFODE::writeData(os);

    // Write type name
    os.writeKeyword("type") << tab << type() << token::END_STATEMENT << endl;

    // Write data
    os.writeKeyword("Xrel") << tab << Xrel_
        << token::END_STATEMENT << nl;
    os.writeKeyword("U") << tab << U_ << token::END_STATEMENT << nl;
    os.writeKeyword("rotationTensor") << tab << rotation_
        << token::END_STATEMENT << nl;
    os.writeKeyword("rotationIncrementTensor") << tab << rotIncrement_
        << token::END_STATEMENT << nl;
    os.writeKeyword("omega") << tab << omega_
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("fixedSurge") << tab << fixedSurge_ <<
        token::END_STATEMENT << nl;
    os.writeKeyword("fixedSway") << tab << fixedSway_ <<
        token::END_STATEMENT << nl;
    os.writeKeyword("fixedHeave") << tab << fixedHeave_ <<
        token::END_STATEMENT << nl;
    os.writeKeyword("fixedRoll") << tab << fixedRoll_ <<
        token::END_STATEMENT << nl;
    os.writeKeyword("fixedPitch") << tab << fixedPitch_ <<
        token::END_STATEMENT << nl;
    os.writeKeyword("fixedYaw") << tab << fixedYaw_ <<
        token::END_STATEMENT << nl << endl;

    return os.good();
}


// ************************************************************************* //

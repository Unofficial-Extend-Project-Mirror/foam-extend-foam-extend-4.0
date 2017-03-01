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
    quaternionSixDOF

Description
    6-DOF solver using quaternions

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Hrvoje Jasak, FSB Zagreb.  All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    quaternionSixDOF.C

\*---------------------------------------------------------------------------*/

#include "quaternionSixDOF.H"
#include "OutputControlDictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(quaternionSixDOF, 0);
addToRunTimeSelectionTable(sixDOFODE, quaternionSixDOF, dictionary);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dimensionedVector Foam::quaternionSixDOF::A
(
    const dimensionedVector& xR,
    const dimensionedVector& uR,
    const HamiltonRodriguezRot& rotation
) const
{
    // Fix the total force in global coordinate system
    dimensionedVector fAbs =
        // Force in global coordinate system
        force()
        // Force in local coordinate system
      + (rotation.invR() & forceRelative())
        // Spring force in global coordinate system
      - (linSpringCoeffs() & xR)
        // Damping force in global coordinate system
      - (linDampingCoeffs() & uR);

    // Constrain translation
    constrainTranslation(fAbs.value());

    return fAbs/mass();
}


Foam::dimensionedVector Foam::quaternionSixDOF::OmegaDot
(
    const HamiltonRodriguezRot& rotation,
    const dimensionedVector& omega
) const
{
    // Fix the global moment for global rotation constraints
    dimensionedVector mAbs = moment();

    // Constrain rotation
    constrainRotation(mAbs.value());

    return
        inv(momentOfInertia())
      & (
            E(omega)
            // To relative
          + (rotation.R() & mAbs)
          + momentRelative()
        );
}


Foam::dimensionedVector Foam::quaternionSixDOF::E
(
    const dimensionedVector& omega
) const
{
    return (*(momentOfInertia() & omega) & omega);
}


void Foam::quaternionSixDOF::constrainRotation(vector& vec) const
{
    // Constrain the vector with respect to referent or global coordinate system
    if (referentMotionConstraints_)
    {
        vector consVec(referentRotation_.R() & vec);

        if (fixedRoll_)
        {
            consVec.x() = 0;
        }
        if (fixedPitch_)
        {
            consVec.y() = 0;
        }
        if (fixedYaw_)
        {
            consVec.z() = 0;
        }

        vec = referentRotation_.invR() & consVec;
    }
    else
    {
        if (fixedRoll_)
        {
            vec.x() = 0;
        }
        if (fixedPitch_)
        {
            vec.y() = 0;
        }
        if (fixedYaw_)
        {
            vec.z() = 0;
        }
    }
}


void Foam::quaternionSixDOF::constrainTranslation(vector& vec) const
{
    // Constrain the vector in respect to referent or global coordinate system
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::quaternionSixDOF::quaternionSixDOF(const IOobject& io)
:
    sixDOFODE(io),

    Xrel_(dict().lookup("Xrel")),
    U_(dict().lookup("U")),
    Uaverage_("Uaverage", U_),
    rotation_
    (
        vector(dict().lookup("rotationVector")),
        dimensionedScalar(dict().lookup("rotationAngle")).value()
    ),
    omega_(dict().lookup("omega")),
    omegaAverage_("omegaAverage", omega_),
    omegaAverageAbsolute_("omegaAverageAbsolute", omega_),

    coeffs_(13, 0.0),

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


Foam::quaternionSixDOF::quaternionSixDOF
(
    const word& name,
    const quaternionSixDOF& qsd
)
:
    sixDOFODE
    (
        IOobject
        (
            name,
            qsd.dict().instance(),
            qsd.dict().local(),
            qsd.dict().db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),

    Xrel_(qsd.Xrel_.name(), qsd.Xrel_),
    U_(qsd.U_.name(), qsd.U_),
    Uaverage_(qsd.Uaverage_.name(), qsd.Uaverage_),
    rotation_(qsd.rotation_),
    omega_(qsd.omega_.name(), qsd.omega_),
    omegaAverage_(qsd.omegaAverage_.name(), qsd.omegaAverage_),
    omegaAverageAbsolute_
    (
        qsd.omegaAverageAbsolute_.name(),
        qsd.omegaAverageAbsolute_
    ),

    coeffs_(qsd.coeffs_),

    fixedSurge_(qsd.fixedSurge_),
    fixedSway_(qsd.fixedSway_),
    fixedHeave_(qsd.fixedHeave_),
    fixedRoll_(qsd.fixedRoll_),
    fixedPitch_(qsd.fixedPitch_),
    fixedYaw_(qsd.fixedYaw_),
    referentMotionConstraints_(qsd.referentMotionConstraints_),
    referentRotation_(qsd.referentRotation_)
{}


Foam::autoPtr<Foam::sixDOFODE> Foam::quaternionSixDOF::clone
(
    const word& name
) const
{
    // Create and return an autoPtr to the new quaternionSixDOF object with a
    // new name
    return autoPtr<sixDOFODE>
    (
        new quaternionSixDOF
        (
            name,
            *this
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::quaternionSixDOF::~quaternionSixDOF()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedVector& Foam::quaternionSixDOF::Xrel() const
{
    return Xrel_;
}


const Foam::dimensionedVector& Foam::quaternionSixDOF::omega() const
{
    return omega_;
}


Foam::dimensionedVector Foam::quaternionSixDOF::X() const
{
    return Xequilibrium() + Xrel_;
}


const Foam::dimensionedVector& Foam::quaternionSixDOF::U() const
{
    return U_;
}


const Foam::dimensionedVector& Foam::quaternionSixDOF::Uaverage() const
{
    return Uaverage_;
}


void Foam::quaternionSixDOF::setState(const sixDOFODE& sd)
{
    // First set the state in base class subobject
    sixDOFODE::setState(sd);

    // Cast sixDOFODE& to quaternionSixDOF&
    const quaternionSixDOF& qsd = refCast<const quaternionSixDOF>(sd);

    // Set state variables for this class
    Xrel_ = qsd.Xrel_;
    U_ = qsd.U_;
    Uaverage_ = qsd.Uaverage_;
    rotation_ = qsd.rotation_;
    omega_ = qsd.omega_;
    omegaAverage_ = qsd.omegaAverage_;
    omegaAverageAbsolute_ = qsd.omegaAverageAbsolute_;

    // Copy ODE coefficients: this carries actual ODE state
    // HJ, 23/Mar/2015
    coeffs_ = qsd.coeffs_;

    fixedSurge_ = qsd.fixedSurge_;
    fixedSway_ = qsd.fixedSway_;
    fixedHeave_ = qsd.fixedHeave_;
    fixedRoll_ = qsd.fixedRoll_;
    fixedPitch_ = qsd.fixedPitch_;
    fixedYaw_ = qsd.fixedYaw_;
    referentMotionConstraints_ = qsd.referentMotionConstraints_;
    referentRotation_ = qsd.referentRotation_;
}


const Foam::dimensionedVector& Foam::quaternionSixDOF::omegaAverage() const
{
    return omegaAverage_;
}


Foam::tensor Foam::quaternionSixDOF::toRelative() const
{
    return rotation_.eCurrent().R();
}


Foam::tensor Foam::quaternionSixDOF::toAbsolute() const
{
    return rotation_.eCurrent().invR();
}


const Foam::tensor& Foam::quaternionSixDOF::rotIncrementTensor() const
{
    return rotation_.rotIncrementTensor();
}


void Foam::quaternionSixDOF::derivatives
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

    dydx[6] = omegaDot.x();
    dydx[7] = omegaDot.y();
    dydx[8] = omegaDot.z();

    dydx[9] = curRotation.eDot(curOmega.value(), 0);
    dydx[10] = curRotation.eDot(curOmega.value(), 1);
    dydx[11] = curRotation.eDot(curOmega.value(), 2);
    dydx[12] = curRotation.eDot(curOmega.value(), 3);

    // Add rotational constraints by setting RHS of given components to zero
    if (fixedRoll_)
    {
        dydx[10] = 0; // Roll axis (roll quaternion evolution RHS)
    }
    if (fixedPitch_)
    {
        dydx[11] = 0; // Pitch axis (pitch quaternion evolution RHS)
    }
    if (fixedYaw_)
    {
        dydx[12] = 0; // Yaw axis (yaw quaternion evolution RHS)
    }
}


void Foam::quaternionSixDOF::update(const scalar delta)
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

    // Constrain velocity
    constrainTranslation(Uval);
    coeffs_[3] = Uval.x();
    coeffs_[4] = Uval.y();
    coeffs_[5] = Uval.z();

    // Update omega
    vector& omegaVal = omega_.value();

    omegaVal.x() = coeffs_[6];
    omegaVal.y() = coeffs_[7];
    omegaVal.z() = coeffs_[8];

    // Constrain omega
    constrainRotation(omegaVal);
    coeffs_[6] = omegaVal.x();
    coeffs_[7] = omegaVal.y();
    coeffs_[8] = omegaVal.z();

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

    // Calculate and constrain omegaAverageAbsolute appropriately
    vector& omegaAverageAbsoluteValue = omegaAverageAbsolute_.value();
    omegaAverageAbsoluteValue = rotation_.omegaAverageAbsolute(delta);

    if (fixedRoll_)
    {
        omegaAverageAbsoluteValue.x() = 0;
    }
    if (fixedPitch_)
    {
        omegaAverageAbsoluteValue.y() = 0;
    }
    if (fixedYaw_)
    {
        omegaAverageAbsoluteValue.z() = 0;
    }
}


bool Foam::quaternionSixDOF::writeData(Ostream& os) const
{
    // First write the part related to base class subobject
    sixDOFODE::writeData(os);

    // Write type name
    os.writeKeyword("type") << tab << type() << token::END_STATEMENT << endl;

    // Write data
    os.writeKeyword("Xrel") << tab << Xrel_
        << token::END_STATEMENT << nl;
    os.writeKeyword("U") << tab << U_ << token::END_STATEMENT << nl;
    os.writeKeyword("rotationVector") << tab << rotation_.rotVector()
        << token::END_STATEMENT << nl;
    os.writeKeyword("rotationAngle") << tab
        << dimensionedScalar("rotAngle", dimless, rotation_.rotAngle())
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

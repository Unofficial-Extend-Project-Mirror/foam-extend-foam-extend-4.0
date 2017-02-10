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
    sixDOFqODE

Description
    6-DOF solver using quaternions

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "sixDOFqODE.H"

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
    // Fix the global force for global rotation constraints
    dimensionedVector fAbs = force();

    // Constrain translation
    constrainTranslation(fAbs.value());

    return
    (
       - (linSpringCoeffs_ & xR)    // spring
       - (linDampingCoeffs_ & uR)   // damping
       + fAbs
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
    // Fix the global moment for global rotation constraints
    dimensionedVector mAbs = moment();

    // Constrain rotation
    constrainRotation(mAbs.value());

    return
        inv(momentOfInertia_)
      & (
            E(omega)
            // To relative
          + (rotation.R() & mAbs)
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


void Foam::sixDOFqODE::constrainRotation(vector& vec) const
{
    vector consVec(vector::zero);

    // Constrain the vector in respect to referent or global coordinate system
    if (referentMotionConstraints_)
    {
        consVec = referentRotation_.R() & vec;

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


void Foam::sixDOFqODE::constrainTranslation(vector& vec) const
{
    vector consVec(vector::zero);

    // Constrain the vector in respect to referent or global coordinate system
    if (referentMotionConstraints_)
    {
        consVec = referentRotation_.R() & vec;

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


void Foam::sixDOFqODE::aitkensRelaxation
(
    const scalar min,
    const scalar max
)
{
    // Calculate translational relax factor
    const scalar saveOldRelFacT = oldRelaxFactorT_;
    oldRelaxFactorT_ = relaxFactorT_;

    if(magSqr(A_[0] - An_[1] - A_[1] - An_[2]) > SMALL)
    {
        relaxFactorT_ = saveOldRelFacT + (saveOldRelFacT - 1)*
        (
            (A_[1] - An_[2])
          & (A_[0] - An_[1] - A_[1] - An_[2])
        )/
        magSqr(A_[0] - An_[1] - A_[1] - An_[2]);
    }
    else
    {
        relaxFactorT_ = min;
    }

    // Bound the relaxation factor for stability
    if (relaxFactorT_ > max)
    {
        relaxFactorT_ = max;
    }
    else if (relaxFactorT_ < min)
    {
        relaxFactorT_ = min;
    }

    // Calculate rotational relax factor
    const scalar saveOldRelFacR = oldRelaxFactorR_;
    oldRelaxFactorR_ = relaxFactorR_;

    if
    (
        magSqr(OmegaDot_[0] - OmegaDotn_[1] - OmegaDot_[1] - OmegaDotn_[2])
      > SMALL
    )
    {
        relaxFactorR_ =
        saveOldRelFacR + (saveOldRelFacR - 1)*
        ((OmegaDot_[1] - OmegaDotn_[2]) &
        (OmegaDot_[0] - OmegaDotn_[1] - OmegaDot_[1] - OmegaDotn_[2]))/
        magSqr(OmegaDot_[0] - OmegaDotn_[1] - OmegaDot_[1] - OmegaDotn_[2]);
    }
    else
    {
        relaxFactorR_ = min;
    }

    // Bound the relaxation factor for stability
    if(relaxFactorR_ > max)
    {
        relaxFactorR_ = max;
    }
    else if(relaxFactorR_ < min)
    {
        relaxFactorR_ = min;
    }
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
    relaxFactorT_(1.0),
    relaxFactorR_(1.0),
    oldRelaxFactorT_(1.0),
    oldRelaxFactorR_(1.0),

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

    A_(3, vector::zero),
    OmegaDot_(3, vector::zero),
    An_(3, vector::zero),
    OmegaDotn_(3, vector::zero),
    force_(lookup("force")),
    moment_(lookup("moment")),
    forceRelative_(lookup("forceRelative")),
    momentRelative_(lookup("momentRelative")),

    coeffs_(13, 0.0),

    fixedSurge_(lookup("fixedSurge")),
    fixedSway_(lookup("fixedSway")),
    fixedHeave_(lookup("fixedHeave")),
    fixedRoll_(lookup("fixedRoll")),
    fixedPitch_(lookup("fixedPitch")),
    fixedYaw_(lookup("fixedYaw")),
    referentMotionConstraints_
    (
        lookupOrDefault<Switch>
        (
            "referentMotionConstraints",
            false
        )
    ),
    referentRotation_(vector(1, 0, 0), 0)
{
    setCoeffs();
}


// Construct as copy
Foam::sixDOFqODE::sixDOFqODE
(
    const word& name,
    const sixDOFqODE& sd
)
:
    IOdictionary
    (
        IOobject
        (
            name,
            sd.instance(),
            sd.local(),
            sd.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),

    mass_(sd.mass_.name(), sd.mass_),
    momentOfInertia_(sd.momentOfInertia_.name(), sd.momentOfInertia_),

    Xequilibrium_(sd.Xequilibrium_.name(), sd.Xequilibrium_),
    linSpringCoeffs_(sd.linSpringCoeffs_.name(), sd.linSpringCoeffs_),
    linDampingCoeffs_(sd.linDampingCoeffs_.name(), sd.linDampingCoeffs_),
    relaxFactorT_(sd.relaxFactorT_),
    relaxFactorR_(sd.relaxFactorR_),
    oldRelaxFactorT_(sd.oldRelaxFactorT_),
    oldRelaxFactorR_(sd.oldRelaxFactorR_),

    Xrel_(sd.Xrel_.name(), sd.Xrel_),
    U_(sd.U_.name(), sd.U_),
    Uaverage_(sd.Uaverage_.name(), sd.Uaverage_),
    rotation_(sd.rotation_),
    omega_(sd.omega_.name(), sd.omega_),
    omegaAverage_(sd.omegaAverage_.name(), sd.omegaAverage_),
    omegaAverageAbsolute_
    (
        sd.omegaAverageAbsolute_.name(),
        sd.omegaAverageAbsolute_
    ),

    A_(sd.A_),
    OmegaDot_(sd.OmegaDot_),
    An_(sd.An_),
    OmegaDotn_(sd.OmegaDotn_),
    force_(sd.force_.name(), sd.force_),
    moment_(sd.moment_.name(), sd.moment_),
    forceRelative_(sd.forceRelative_.name(), sd.forceRelative_),
    momentRelative_(sd.momentRelative_.name(), sd.momentRelative_),

    coeffs_(sd.coeffs_),

    fixedSurge_(sd.fixedSurge_),
    fixedSway_(sd.fixedSway_),
    fixedHeave_(sd.fixedHeave_),
    fixedRoll_(sd.fixedRoll_),
    fixedPitch_(sd.fixedPitch_),
    fixedYaw_(sd.fixedYaw_),
    referentMotionConstraints_(sd.referentMotionConstraints_),
    referentRotation_(sd.referentRotation_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDOFqODE::~sixDOFqODE()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDOFqODE::setState(const sixDOFqODE& sd)
{
    // Set state does not copy AList_, AOld_, relaxFactor_ and relaxFactorOld_ to
    // allow Aitkens relaxation (IG 5/May/2016)
    mass_ = sd.mass_;
    momentOfInertia_ = sd.momentOfInertia_;

    Xequilibrium_ = sd.Xequilibrium_;
    linSpringCoeffs_ = sd.linSpringCoeffs_;
    linDampingCoeffs_ = sd.linDampingCoeffs_;

    Xrel_ = sd.Xrel_;
    U_ = sd.U_;
    Uaverage_ = sd.Uaverage_;
    rotation_ = sd.rotation_;
    omega_ = sd.omega_;
    omegaAverage_ = sd.omegaAverage_;
    omegaAverageAbsolute_ = sd.omegaAverageAbsolute_;

    force_ = sd.force_;
    moment_ = sd.moment_;
    forceRelative_ = sd.forceRelative_;
    momentRelative_ = sd.momentRelative_;

    // Copy ODE coefficients: this carries actual ODE state
    // HJ, 23/Mar/2015
    coeffs_ = sd.coeffs_;

    fixedSurge_ = sd.fixedSurge_;
    fixedSway_ = sd.fixedSway_;
    fixedHeave_ = sd.fixedHeave_;
    fixedRoll_ = sd.fixedRoll_;
    fixedPitch_ = sd.fixedPitch_;
    fixedYaw_ = sd.fixedYaw_;
}


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


void Foam::sixDOFqODE::relaxAcceleration
(
    const scalar minRelFactor,
    const scalar maxRelFactor
)
{
    if (mag(minRelFactor - maxRelFactor) < SMALL)
    {
        // Fixed relaxation
       relaxFactorT_ = minRelFactor;
       relaxFactorR_ = minRelFactor;
    }
    else
    {
        // Use Aitkens relaxation

        // Update Aitkens relaxation factor
        aitkensRelaxation(minRelFactor, maxRelFactor);

        // Update non relaxed accelerations
        An_[1] = An_[2];
        An_[2] = (force()/mass_).value();
        OmegaDotn_[1] = OmegaDotn_[2];
        OmegaDotn_[2] = (inv(momentOfInertia_) & moment()).value();
    }

    const vector Aold = A_[2];

    const vector OmegaDotold = OmegaDot_[2];

    force().value() =
        Aold*mass_.value()
      + relaxFactorT_*(force().value() - Aold*mass_.value());

    moment().value() =
        (momentOfInertia_.value() & OmegaDotold)
      + relaxFactorR_*
        (
            moment().value()
          - (momentOfInertia_.value() & OmegaDotold)
        );

    // Update relaxed old accelerations
    A_[0] = A_[1];
    A_[1] = A_[2];
    A_[2] = (force()/mass_).value();
    OmegaDot_[0] = OmegaDot_[1];
    OmegaDot_[1] = OmegaDot_[2];
    OmegaDot_[2] = (inv(momentOfInertia_) & moment()).value();
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

    os.writeKeyword("fixedSurge") << tab << sds.fixedSurge_ <<
        token::END_STATEMENT << endl;
    os.writeKeyword("fixedSway") << tab << sds.fixedSway_ <<
        token::END_STATEMENT << endl;
    os.writeKeyword("fixedHeave") << tab << sds.fixedHeave_ <<
        token::END_STATEMENT << endl;
    os.writeKeyword("fixedRoll") << tab << sds.fixedRoll_ <<
        token::END_STATEMENT << endl;
    os.writeKeyword("fixedPitch") << tab << sds.fixedPitch_ <<
        token::END_STATEMENT << endl;
    os.writeKeyword("fixedYaw") << tab << sds.fixedYaw_ <<
        token::END_STATEMENT << endl;

    return os;
}


// ************************************************************************* //

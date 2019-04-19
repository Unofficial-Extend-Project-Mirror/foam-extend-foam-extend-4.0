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
    quaternionSixDOF

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
    const HamiltonRodriguezRot& rotation,
    const scalar t
) const
{
    // Create a scalar square matrix representing Newton equations with
    // constraints and the corresponding source (right hand side vector).
    // Note: size of the matrix is 3 + number of constraints
    scalarField rhs(translationalConstraints().size() + 3, 0.0);
    scalarSquareMatrix M(rhs.size(), 0.0);

    // Insert mass and explicit forcing into the system. Note: translations are
    // solved in the global coordinate system and the explicit forcing contains
    // restraining forces
    const dimensionedVector explicitForcing =
        force
        (
            t,
            rotation.R(),
            xR.value(),
            uR.value()
        );

    const vector& efVal = explicitForcing.value();
    const scalar& m = mass().value();

    forAll(efVal, dirI)
    {
        M[dirI][dirI] = m;

        rhs[dirI] = efVal[dirI];
    }

    // Insert contributions from the constraints
    forAll(translationalConstraints(), tcI)
    {
        // Get reference to current constraint
        const translationalConstraint& curTc = translationalConstraints()[tcI];

        // Get matrix contribution from constraint
        const vector mc =
            curTc.matrixContribution
            (
                t,
                rotation.R(),
                xR.value(),
                uR.value()
            );

        // Get matrix index
        const label index = tcI + 3;

        // Insert contributions into the matrix
        forAll(mc, dirI)
        {
            M[dirI][index] = mc[dirI];
            M[index][dirI] = mc[dirI];
        }

        // Insert source contribution (remainder of the constraint function)
        rhs[index] =
            curTc.sourceContribution
           (
               t,
               rotation.R(),
               xR.value(),
               uR.value()
           );
    }

    // Solve the matrix using LU decomposition. Note: solution is in the rhs and
    // it contains accelerations in the first three entries and corresponding
    // Lagrangian multipliers in other entries.
    scalarSquareMatrix::LUsolve(M, rhs);

    return
        dimensionedVector
        (
            "A",
            force().dimensions()/mass().dimensions(),
            vector(rhs[0], rhs[1], rhs[2])
        );
}


Foam::dimensionedVector Foam::quaternionSixDOF::OmegaDot
(
    const HamiltonRodriguezRot& rotation,
    const dimensionedVector& omega,
    const scalar t
) const
{
    // Create a scalar square matrix representing Euler equations with
    // constraints and the corresponding source (right hand side vector).
    // Note: size of the matrix is 3 + number of constraints
    scalarField rhs(rotationalConstraints().size() + 3, 0.0);
    scalarSquareMatrix J(rhs.size(), 0.0);

    // Get current inertial-to-local transformation. Note: different convention
    // (R represents coordinate transformation from global to local)
    const dimensionedTensor R("R", dimless, rotation.R());

    // Insert moment of inertia and explicit forcing into the system
    const dimensionedVector explicitForcing
    (
        E(omega) // Euler part
      + (
          R.value()
        & moment
          (
              t,
              R.value(),
              omega.value()
          )
        ) // External torque with restraints
    );
    const vector& efVal = explicitForcing.value();
    const diagTensor& I = momentOfInertia().value();

    forAll(efVal, dirI)
    {
        J[dirI][dirI] = I[dirI];

        rhs[dirI] = efVal[dirI];
    }

    // Insert contributions from the constraints
    forAll(rotationalConstraints(), rcI)
    {
        // Get reference to current constraint
        const rotationalConstraint& curRc = rotationalConstraints()[rcI];

        // Get matrix contribution from the constraint
        const vector mc =
            curRc.matrixContribution
            (
                t,
                R.value(),
                omega.value()
            );

        // Get matrix index
        const label index = rcI + 3;

        // Insert contributions into the matrix
        forAll(mc, dirI)
        {
            J[dirI][index] = mc[dirI];
            J[index][dirI] = mc[dirI];
        }

        // Insert source contribution (remainder of the constraint function)
        rhs[index] = curRc.sourceContribution(t, R.value(), omega.value());
    }

    // Solve the matrix using LU decomposition. Note: solution is in the rhs and
    // it contains OmegaDot's in the first three entries and corresponding
    // Lagrangian multipliers in other entries.
    scalarSquareMatrix::LUsolve(J, rhs);

    return
        dimensionedVector
        (
            "OmegaDot",
            moment().dimensions()/momentOfInertia().dimensions(),
            vector(rhs[0], rhs[1], rhs[2])
        );
}


Foam::dimensionedVector Foam::quaternionSixDOF::E
(
    const dimensionedVector& omega
) const
{
    return (*(momentOfInertia() & omega) & omega);
}


// * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * * * //

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

    // Copy ODE coefficients: this carries actual ODE state
    // HJ, 23/Mar/2015
    coeffs_ = qsd.coeffs_;
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

    coeffs_(13, 0.0)
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
    sixDOFODE(name, qsd),

    Xrel_(qsd.Xrel_.name(), qsd.Xrel_),
    U_(qsd.U_.name(), qsd.U_),
    Uaverage_(qsd.Uaverage_.name(), qsd.Uaverage_),
    rotation_(qsd.rotation_),
    omega_(qsd.omega_.name(), qsd.omega_),
    omegaAverage_(qsd.omegaAverage_.name(), qsd.omegaAverage_),

    coeffs_(qsd.coeffs_)
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


const Foam::dimensionedVector& Foam::quaternionSixDOF::omegaAverage() const
{
    return omegaAverage_;
}


Foam::dimensionedVector
Foam::quaternionSixDOF::translationalAcceleration() const
{
    // Calculate and return translational acceleration in global c. s.
    return A(Xrel(), Uaverage(), rotation_.eCurrent(), dict().time().value());
}


Foam::dimensionedVector
Foam::quaternionSixDOF::rotationalAcceleration() const
{
    // Calculate and return rotational acceleration in relative c. s.
    return
        OmegaDot(rotation_.eCurrent(), omegaAverage(), dict().time().value());
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

    const vector accel = A(curX, curU, curRotation, x).value();

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

    const vector omegaDot = OmegaDot(curRotation, curOmega, x).value();

    dydx[6] = omegaDot.x();
    dydx[7] = omegaDot.y();
    dydx[8] = omegaDot.z();

    dydx[9] = curRotation.eDot(curOmega.value(), 0);
    dydx[10] = curRotation.eDot(curOmega.value(), 1);
    dydx[11] = curRotation.eDot(curOmega.value(), 2);
    dydx[12] = curRotation.eDot(curOmega.value(), 3);
}


void Foam::quaternionSixDOF::update(const scalar delta)
{
    // Translation

    // Get displacement
    const vector Xold = Xrel_.value();

    vector& Xval = Xrel_.value();
    Xval.x() = coeffs_[0];
    Xval.y() = coeffs_[1];
    Xval.z() = coeffs_[2];

    // Get velocity
    vector& Uval = U_.value();
    Uval.x() = coeffs_[3];
    Uval.y() = coeffs_[4];
    Uval.z() = coeffs_[5];

    // Stabilise translational constraints if necessary
    forAll(translationalConstraints(), tcI)
    {
        // Note: get end time for this ODE step from mesh, assuming that the top
        // level solves this ode from [t - deltaT, t], yielding solution at
        // t. Done this way to preserve the interface of ODE class.
        // VV, 10/Mar/2017.
        const scalar t = dict().time().value();

        translationalConstraints()[tcI].stabilise(t, Xval, Uval);
    }

    // Update (possibly constrained) displacement
    coeffs_[0] = Xval.x();
    coeffs_[1] = Xval.y();
    coeffs_[2] = Xval.z();

    // Update (possibly constrained) velocity
    coeffs_[3] = Uval.x();
    coeffs_[4] = Uval.y();
    coeffs_[5] = Uval.z();

    // Update average velocity
    Uaverage_.value() = (Xval - Xold)/delta;


    // Rotation

    // Get angular velocity
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

    // Stabilise rotational constraints if necessary
    forAll(rotationalConstraints(), rcI)
    {
        // Note: get end time for this ODE step from mesh, assuming that the top
        // level solves this ode from [t - deltaT, t], yielding solution at
        // t. Done this way to preserve the interface of ODE class.
        // VV, 10/Mar/2017.
        const scalar t = dict().time().value();

        rotationalConstraints()[rcI].stabilise(t, omegaVal);
    }

    // Update (possibly constrained) omega
    coeffs_[6] = omegaVal.x();
    coeffs_[7] = omegaVal.y();
    coeffs_[8] = omegaVal.z();

    // Update average omega
    omegaAverage_.value() = rotation_.omegaAverage(delta);
}


bool Foam::quaternionSixDOF::writeData(Ostream& os) const
{
    // First write the part related to base class subobject
    sixDOFODE::writeData(os);

    // Write type name
    os.writeKeyword("type") << tab << type()
        << token::END_STATEMENT << nl << nl;

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

    return os.good();
}


// ************************************************************************* //

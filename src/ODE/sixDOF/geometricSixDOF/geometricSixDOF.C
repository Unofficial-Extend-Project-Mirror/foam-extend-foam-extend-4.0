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
    geometricSixDOF

\*---------------------------------------------------------------------------*/

#include "geometricSixDOF.H"
#include "scalarSquareMatrix.H"
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
    "geometricSixDOFRotIncTensorTol",
    1e-10
);


const Foam::debug::tolerancesSwitch Foam::geometricSixDOF::rotIncRateTol_
(
    "geometricSixDOFRotIncRateTol",
    1e-6
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dimensionedVector Foam::geometricSixDOF::A
(
    const dimensionedVector& xR,
    const dimensionedVector& uR,
    const tensor& R,
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
            R.T(),
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
                R.T(),
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
        rhs[index] = curTc.sourceContribution(t, R.T(), xR.value(), uR.value());
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


Foam::dimensionedVector Foam::geometricSixDOF::OmegaDot
(
    const tensor& R,
    const dimensionedVector& omega,
    const scalar t
) const
{
    // Create a scalar square matrix representing Euler equations with
    // constraints and the corresponding source (right hand side vector).
    // Note: size of the matrix is 3 + number of constraints
    scalarField rhs(rotationalConstraints().size() + 3, 0.0);
    scalarSquareMatrix J(rhs.size(), 0.0);

    // Get current inertial-to-local transformation
    const dimensionedTensor RT("RT", dimless, R.T());

    // Insert moment of inertia and explicit forcing into the system
    const dimensionedVector explicitForcing
    (
        E(omega) // Euler part
      + (
          RT
        & moment
          (
              t,
              RT.value(),
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
                RT.value(),
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
        rhs[index] = curRc.sourceContribution(t, RT.value(), omega.value());
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


Foam::dimensionedVector Foam::geometricSixDOF::E
(
    const dimensionedVector& omega
) const
{
    return (*(momentOfInertia() & omega) & omega);
}


Foam::tensor Foam::geometricSixDOF::expMap(const vector& rotInc) const
{
    tensor R;

    // Calculate the magnitude of the rotational increment vector
    const scalar magRotInc = mag(rotInc);

    if (magRotInc < rotIncTensorTol_())
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

    return R;
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

    if (magRotInc < rotIncRateTol_())
    {
        // Stabilised calculation of rotation increment to avoid small
        // denominators
        const tensor lbRotIncOmega(lieBracket(rotInc, omega));

        rotIncDot =
        (
            I
          + lbRotIncOmega/2.0
          + lieBracket(rotInc, *lbRotIncOmega)/12.0
        ) & omega;
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
        ) & omega;
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

    coeffs_(12, 0.0)
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

    coeffs_(gsd.coeffs_)
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


Foam::dimensionedVector Foam::geometricSixDOF::translationalAcceleration() const
{
    // Calculate and return translational acceleration in global c. s.
    return A(Xrel(), Uaverage(), rotation_, dict().time().value());
}


Foam::dimensionedVector Foam::geometricSixDOF::rotationalAcceleration() const
{
    // Calculate and return rotational acceleration in relative c. s.
    return OmegaDot(rotation_, omegaAverage(), dict().time().value());
}


Foam::tensor Foam::geometricSixDOF::toRelative() const
{
    return rotation_.T();
}


Foam::tensor Foam::geometricSixDOF::toAbsolute() const
{
    return rotation_;
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
    const tensor curRot = (rotation_ & expMap(rotIncrementVector));

    // Calculate translational acceleration using current rotation
    const vector accel = A(curX, curU, curRot, x).value();

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
    const vector omegaDot = OmegaDot(curRot, curOmega, x).value();

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
    const vector omegaOld = omega_.value();

    vector& omegaVal = omega_.value();
    omegaVal.x() = coeffs_[6];
    omegaVal.y() = coeffs_[7];
    omegaVal.z() = coeffs_[8];

    // Update rotational increment tensor
    rotIncrement_ = expMap(vector(coeffs_[9], coeffs_[10], coeffs_[11]));

    // Update rotational tensor
    rotation_ = (rotation_ & rotIncrement_);

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

    // Reset increment vector in coefficients for the next step
    coeffs_[9] = 0;
    coeffs_[10] = 0;
    coeffs_[11] = 0;

    // Update average omega
    omegaAverage_.value() = 0.5*(omegaVal + omegaOld);
}


bool Foam::geometricSixDOF::writeData(Ostream& os) const
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
    os.writeKeyword("rotationTensor") << tab << rotation_
        << token::END_STATEMENT << nl;
    os.writeKeyword("rotationIncrementTensor") << tab << rotIncrement_
        << token::END_STATEMENT << nl;
    os.writeKeyword("omega") << tab << omega_
        << token::END_STATEMENT << nl << endl;

    return os.good();
}


// ************************************************************************* //

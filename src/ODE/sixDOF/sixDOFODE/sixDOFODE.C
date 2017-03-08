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
    sixDOFODE

Description
    Abstract base class for six-degrees-of-freedom (6DOF) ordinary differential
    equations with necessary interface for two-way coupling with CFD solvers.

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Hrvoje Jasak, FSB Zagreb.  All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "sixDOFODE.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(sixDOFODE, 0);
defineRunTimeSelectionTable(sixDOFODE, dictionary);

}

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

void Foam::sixDOFODE::updateRelaxFactors()
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
        relaxFactorT_ = minRelaxFactor_;
    }

    // Bound the relaxation factor for stability
    if (relaxFactorT_ > maxRelaxFactor_)
    {
        relaxFactorT_ = maxRelaxFactor_;
    }
    else if (relaxFactorT_ < minRelaxFactor_)
    {
        relaxFactorT_ = minRelaxFactor_;
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
        relaxFactorR_ = saveOldRelFacR
          + (saveOldRelFacR - 1)*
            (
                (OmegaDot_[1] - OmegaDotn_[2])
              & (OmegaDot_[0] - OmegaDotn_[1] - OmegaDot_[1] - OmegaDotn_[2])
            )/
            magSqr(OmegaDot_[0] - OmegaDotn_[1] - OmegaDot_[1] - OmegaDotn_[2]);
    }
    else
    {
        relaxFactorR_ = minRelaxFactor_;
    }

    // Bound the relaxation factor for stability
    if (relaxFactorR_ > maxRelaxFactor_)
    {
        relaxFactorR_ = maxRelaxFactor_;
    }
    else if (relaxFactorR_ < minRelaxFactor_)
    {
        relaxFactorR_ = minRelaxFactor_;
    }
}


void Foam::sixDOFODE::relaxAcceleration()
{
    if (mag(minRelaxFactor_ - maxRelaxFactor_) < SMALL)
    {
       // Fixed relaxation
       relaxFactorT_ = minRelaxFactor_;
       relaxFactorR_ = minRelaxFactor_;
    }
    else
    {
        // Use Aitkens relaxation

        // Update Aitkens relaxation factor
        updateRelaxFactors();

        // Update non relaxed accelerations
        An_[1] = An_[2];
        An_[2] = (force_/mass_).value();
        OmegaDotn_[1] = OmegaDotn_[2];
        OmegaDotn_[2] = (inv(momentOfInertia_) & moment_).value();
    }

    const vector Aold = A_[2];
    const vector OmegaDotold = OmegaDot_[2];

    force_.value() =
        Aold*mass_.value()
      + relaxFactorT_*(force_.value() - Aold*mass_.value());

    moment_.value() =
        (momentOfInertia_.value() & OmegaDotold)
      + relaxFactorR_*
        (
            moment_.value()
          - (momentOfInertia_.value() & OmegaDotold)
        );

    // Update relaxed old accelerations
    A_[0] = A_[1];
    A_[1] = A_[2];
    A_[2] = (force_/mass_).value();
    OmegaDot_[0] = OmegaDot_[1];
    OmegaDot_[1] = OmegaDot_[2];
    OmegaDot_[2] = (inv(momentOfInertia_) & moment_).value();
}


void Foam::sixDOFODE::initState()
{
    // Get time index
    const label timeIndex = dict().time().timeIndex();

    if (curTimeIndex_ < timeIndex)
    {
        // First call in this time index, store data
        oldStatePtr_ = this->clone(dict().name() + "_0");

        Info<< "First 6DOF solution within a time step, storing old data..."
            << endl;
    }
    else
    {
        // Multiple calls in this time index, reset this data
        this->setState(oldStatePtr_());

        Info<< "Repeated 6DOF solution within a time step, restoring data..."
            << endl;
    }

    // Update local time index
    curTimeIndex_ = timeIndex;
}


// * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * * * //

void Foam::sixDOFODE::setState(const sixDOFODE& sd)
{
    // Set state does not copy AList_, AOld_, relaxFactor_ and relaxFactorOld_.
    // In case of multiple updates, overwriting Aitkens relaxation parameters
    // would invalidate the underrelaxation.  IG, 5/May/2016
    mass_ = sd.mass_;
    momentOfInertia_ = sd.momentOfInertia_;

    Xequilibrium_ = sd.Xequilibrium_;
    linSpringCoeffs_ = sd.linSpringCoeffs_;
    linDampingCoeffs_ = sd.linDampingCoeffs_;

    force_ = sd.force_;
    moment_ = sd.moment_;
}


Foam::dimensionedVector Foam::sixDOFODE::force(const scalar t) const
{
    // Get ODE step fraction
    const scalar alpha = odeStepFraction(t);

    // Return linearly interpolated external force
    return (alpha*oldStatePtr_->force() + (1 - alpha)*force());
}


Foam::dimensionedVector Foam::sixDOFODE::moment(const scalar t) const
{
    // Get ODE step fraction
    const scalar alpha = odeStepFraction(t);

    // Return linearly interpolated external moment
    return (alpha*oldStatePtr_->moment() + (1 - alpha)*moment());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDOFODE::sixDOFODE(const IOobject& io)
:
    ODE(),
    dict_(io, *this),

    mass_(dict_.lookup("mass")),
    momentOfInertia_(dict_.lookup("momentOfInertia")),

    Xequilibrium_(dict_.lookup("equilibriumPosition")),
    linSpringCoeffs_(dict_.lookup("linearSpring")),
    linDampingCoeffs_(dict_.lookup("linearDamping")),

    aitkensRelaxation_
    (
        dict_.lookupOrDefault<Switch>("useAitkensRelaxation", false)
    ),
    minRelaxFactor_(dict_.lookupOrDefault<scalar>("minRelaxFactor", 0.1)),
    maxRelaxFactor_(dict_.lookupOrDefault<scalar>("maxRelaxFactor", 0.5)),

    relaxFactorT_(1.0),
    relaxFactorR_(1.0),
    oldRelaxFactorT_(1.0),
    oldRelaxFactorR_(1.0),

    A_(3, vector::zero),
    OmegaDot_(3, vector::zero),
    An_(3, vector::zero),
    OmegaDotn_(3, vector::zero),

    force_(dict_.lookup("force")),
    moment_(dict_.lookup("moment")),

    curTimeIndex_(-1),
    oldStatePtr_()
{
    // Sanity checks
    if (mass_.value() < SMALL)
    {
        FatalErrorIn("sixDOFODE::sixDOFODE(const IOobject& io)")
            << "Zero or negative mass detected: " << mass_.value()
            << nl << "Please check " << dict_.name() << "dictionary."
            << exit(FatalError);
    }

    if (cmptMin(momentOfInertia_.value()) < SMALL)
    {
        FatalErrorIn("sixDOFODE::sixDOFODE(const IOobject& io)")
            << "Zero or negative moment of inertia detected: "
            << momentOfInertia_.value()
            << nl << "Please check " << dict_.name() << "dictionary."
            << exit(FatalError);
    }
}


Foam::sixDOFODE::sixDOFODE(const word& name, const sixDOFODE& sd)
:
    ODE(),
    dict_(sd.dict_),

    mass_(sd.mass_),
    momentOfInertia_(sd.momentOfInertia_),

    Xequilibrium_(sd.Xequilibrium_),
    linSpringCoeffs_(sd.linSpringCoeffs_),
    linDampingCoeffs_(sd.linDampingCoeffs_),

    aitkensRelaxation_(sd.aitkensRelaxation_),
    minRelaxFactor_(sd.minRelaxFactor_),
    maxRelaxFactor_(sd.maxRelaxFactor_),

    relaxFactorT_(1.0),
    relaxFactorR_(1.0),
    oldRelaxFactorT_(1.0),
    oldRelaxFactorR_(1.0),

    A_(3, vector::zero),
    OmegaDot_(3, vector::zero),
    An_(3, vector::zero),
    OmegaDotn_(3, vector::zero),

    force_(sd.force_),
    moment_(sd.moment_),

    curTimeIndex_(sd.curTimeIndex_),
    oldStatePtr_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDOFODE::~sixDOFODE()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::OutputControlDictionary<Foam::sixDOFODE>&
Foam::sixDOFODE::dict() const
{
    return dict_;
}


bool Foam::sixDOFODE::writeData(Ostream& os) const
{
    os.writeKeyword("mass") << tab << mass_ << token::END_STATEMENT << nl;
    os.writeKeyword("momentOfInertia") << tab << momentOfInertia_
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("equilibriumPosition") << tab << Xequilibrium_
        << token::END_STATEMENT << nl;
    os.writeKeyword("linearSpring") << tab << linSpringCoeffs_
        << token::END_STATEMENT << nl;
    os.writeKeyword("linearDamping") << tab << linDampingCoeffs_
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("useAitkensRelaxation") << tab << aitkensRelaxation_
        << token::END_STATEMENT << nl;
    os.writeKeyword("minRelaxFactor") << tab << minRelaxFactor_
        << token::END_STATEMENT << nl;
    os.writeKeyword("maxRelaxFactor") << tab << maxRelaxFactor_
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("force") << tab << force_
        << token::END_STATEMENT << nl;
    os.writeKeyword("moment") << tab << moment_
        << token::END_STATEMENT << nl << nl;

    return os.good();
}


// ************************************************************************* //

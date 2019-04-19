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
    sixDOFODE

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

    force_ = sd.force_;
    moment_ = sd.moment_;


    // Copy constraints

    // Translational constraints
    translationalConstraints_.setSize(sd.translationalConstraints_.size());
    forAll(sd.translationalConstraints_, trI)
    {
        translationalConstraints_.set
        (
            trI,
            sd.translationalConstraints_[trI].clone()
        );
    }

    // Rotational constraints
    rotationalConstraints_.setSize(sd.rotationalConstraints_.size());
    forAll(sd.rotationalConstraints_, rrI)
    {
        rotationalConstraints_.set
        (
            rrI,
            sd.rotationalConstraints_[rrI].clone()
        );
    }


    // Copy restraints

    // Translational restraints
    translationalRestraints_.setSize(sd.translationalRestraints_.size());
    forAll(sd.translationalRestraints_, trI)
    {
        translationalRestraints_.set
        (
            trI,
            sd.translationalRestraints_[trI].clone()
        );
    }

    // Rotational restraints
    rotationalRestraints_.setSize(sd.rotationalRestraints_.size());
    forAll(sd.rotationalRestraints_, rrI)
    {
        rotationalRestraints_.set
        (
            rrI,
            sd.rotationalRestraints_[rrI].clone()
        );
    }

    // Combined restraints
    combinedRestraints_.setSize(sd.combinedRestraints_.size());
    forAll(sd.combinedRestraints_, rrI)
    {
        combinedRestraints_.set
        (
            rrI,
            sd.combinedRestraints_[rrI].clone()
        );
    }
}


Foam::dimensionedVector Foam::sixDOFODE::force
(
    const scalar t,
    const tensor& toRelative,
    const vector& x,
    const vector& u
) const
{
    // Get ODE step fraction
    const scalar alpha = odeStepFraction(t);

    // Calculate restraining force
    dimensionedVector rForce("zero", dimForce, vector::zero);

    forAll(translationalRestraints_, trI)
    {
        rForce.value() += translationalRestraints_[trI].restrainingForce
        (
            t,          // Time
            toRelative, // Transformation tensor
            x,          // Position in the global c.s.
            u           // Velocity in the global c.s.
        );
    }

    forAll(combinedRestraints_, crI)
    {
        rForce.value() += combinedRestraints_[crI].restrainingForce
        (
            t,          // Time
            toRelative, // Transformation tensor
            x,          // Position in the global c.s.
            u           // Velocity in the global c.s.
        );
    }

    // Return linearly interpolated external force with restraining force
    return (alpha*oldStatePtr_->force() + (1 - alpha)*force()) + rForce;
}


Foam::dimensionedVector Foam::sixDOFODE::moment
(
    const scalar t,
    const tensor& toRelative,
    const vector& omega
) const
{
    // Get ODE step fraction
    const scalar alpha = odeStepFraction(t);

    // Calculate restraining moment
    dimensionedVector rMoment("zero", dimForce*dimLength, vector::zero);

    forAll(rotationalRestraints_, rrI)
    {
        rMoment.value() += rotationalRestraints_[rrI].restrainingMoment
        (
            t,          // Time
            toRelative, // Transformation tensor
            omega       // Angular velocity in local c.s.
        );
    }

    forAll(combinedRestraints_, crI)
    {
        rMoment.value() += combinedRestraints_[crI].restrainingMoment
        (
            t,          // Time
            toRelative, // Transformation tensor
            omega       // Angular velocity in local c.s.
        );
    }

    // Return linearly interpolated external moment with restraining moment
    return (alpha*oldStatePtr_->moment() + (1 - alpha)*moment() + rMoment);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDOFODE::sixDOFODE(const IOobject& io)
:
    ODE(),
    dict_(io, *this),

    mass_(dict_.lookup("mass")),
    momentOfInertia_(dict_.lookup("momentOfInertia")),

    Xequilibrium_(dict_.lookup("equilibriumPosition")),

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
    oldStatePtr_(),

    translationalConstraints_(),
    rotationalConstraints_(),

    translationalRestraints_(),
    rotationalRestraints_(),
    combinedRestraints_()
{
    // Sanity checks
    if (mass_.value() < SMALL)
    {
        FatalIOErrorIn
        (
            "sixDOFODE::sixDOFODE(const IOobject& io)",
            dict_
        )   << "Zero or negative mass detected: " << mass_.value()
            << nl << "Please check " << dict_.name() << "dictionary."
            << exit(FatalIOError);
    }

    if (cmptMin(momentOfInertia_.value()) < SMALL)
    {
        FatalIOErrorIn
        (
            "sixDOFODE::sixDOFODE(const IOobject& io)",
            dict_
        )   << "Zero or negative moment of inertia detected: "
            << momentOfInertia_.value()
            << nl << "Please check " << dict_.name() << "dictionary."
            << exit(FatalIOError);
    }

    if
    (
        (minRelaxFactor_ < SMALL)
     || (maxRelaxFactor_ > 1.0)
     || ((maxRelaxFactor_ - minRelaxFactor_) < 0)
    )
    {
        FatalIOErrorIn
        (
            "sixDOFODE::sixDOFODE(const IOobject& io)",
            dict_
        )   << "Invalid minRelaxFactor and maxRelaxFactor specified."
            << nl << "Please use values within 0 and 1."
            << exit(FatalIOError);
    }

    // Read and construct constraints and restraints

    // Read translation constraints if they are present
    if (dict().found("translationalConstraints"))
    {
        PtrList<translationalConstraint> tcList
        (
            dict().lookup("translationalConstraints"),
            translationalConstraint::iNew(*this)
        );
        translationalConstraints_.transfer(tcList);
    }

    // Read rotation constraints if they are present
    if (dict().found("rotationalConstraints"))
    {
        PtrList<rotationalConstraint> rcList
        (
            dict().lookup("rotationalConstraints"),
            rotationalConstraint::iNew(*this)
        );
        rotationalConstraints_.transfer(rcList);
    }

    // Read translation restraints if they are present
    if (dict().found("translationalRestraints"))
    {
        PtrList<translationalRestraint> tcList
        (
            dict().lookup("translationalRestraints"),
            translationalRestraint::iNew(*this)
        );
        translationalRestraints_.transfer(tcList);
    }

    // Read rotation restraints if they are present
    if (dict().found("rotationalRestraints"))
    {
        PtrList<rotationalRestraint> rcList
        (
            dict().lookup("rotationalRestraints"),
            rotationalRestraint::iNew(*this)
        );
        rotationalRestraints_.transfer(rcList);
    }

    // Read rotation restraints if they are present
    if (dict().found("combinedRestraints"))
    {
        PtrList<combinedRestraint> rcList
        (
            dict().lookup("combinedRestraints"),
            combinedRestraint::iNew(*this)
        );
        combinedRestraints_.transfer(rcList);
    }
}


Foam::sixDOFODE::sixDOFODE(const word& name, const sixDOFODE& sd)
:
    ODE(),
    dict_(sd.dict_),

    mass_(sd.mass_),
    momentOfInertia_(sd.momentOfInertia_),

    Xequilibrium_(sd.Xequilibrium_),

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
    oldStatePtr_(),

    translationalConstraints_(sd.translationalConstraints_),
    rotationalConstraints_(sd.rotationalConstraints_),

    translationalRestraints_(sd.translationalRestraints_),
    rotationalRestraints_(sd.rotationalRestraints_),
    combinedRestraints_(sd.combinedRestraints_)
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

    if (!translationalConstraints_.empty())
    {
        os.writeKeyword("translationalConstraints")
            << translationalConstraints_
            << token::END_STATEMENT << nl << endl;
    }

    if (!rotationalConstraints_.empty())
    {
        os.writeKeyword("rotationalConstraints")
            << rotationalConstraints_
            << token::END_STATEMENT << endl;
    }

    if (!translationalRestraints_.empty())
    {
        os.writeKeyword("translationalRestraints")
            << translationalRestraints_
            << token::END_STATEMENT << nl << endl;
    }

    if (!rotationalRestraints_.empty())
    {
        os.writeKeyword("rotationalRestraints")
            << rotationalRestraints_
            << token::END_STATEMENT << endl;
    }

    if (!combinedRestraints_.empty())
    {
        os.writeKeyword("combinedRestraints")
            << combinedRestraints_
            << token::END_STATEMENT << endl;
    }

    return os.good();
}


// ************************************************************************* //

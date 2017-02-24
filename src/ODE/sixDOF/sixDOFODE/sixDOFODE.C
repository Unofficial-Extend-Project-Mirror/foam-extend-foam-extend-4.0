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
    equations

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Hrvoje Jasak, FSB Zagreb. All rights reserved.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::sixDOFODE::aitkensRelaxation
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

Foam::sixDOFODE::sixDOFODE(const IOobject& io)
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

    A_(3, vector::zero),
    OmegaDot_(3, vector::zero),
    An_(3, vector::zero),
    OmegaDotn_(3, vector::zero),

    force_(lookup("force")),
    moment_(lookup("moment")),
    forceRelative_(lookup("forceRelative")),
    momentRelative_(lookup("momentRelative"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDOFODE::~sixDOFODE()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::sixDOFODE::relaxAcceleration
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


void Foam::sixDOFODE::setState(const sixDOFODE& sd)
{
    // Set state does not copy AList_, AOld_, relaxFactor_ and
    // relaxFactorOld_. In case of multiple updates, overwriting Aitkens
    // relaxation parameters would invalidate the underrelaxation.
    // IG, 5/May/2016
    mass_ = sd.mass_;
    momentOfInertia_ = sd.momentOfInertia_;

    Xequilibrium_ = sd.Xequilibrium_;
    linSpringCoeffs_ = sd.linSpringCoeffs_;
    linDampingCoeffs_ = sd.linDampingCoeffs_;

    force_ = sd.force_;
    moment_ = sd.moment_;
    forceRelative_ = sd.forceRelative_;
    momentRelative_ = sd.momentRelative_;
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

    os.writeKeyword("force") << tab << force_
        << token::END_STATEMENT << nl;
    os.writeKeyword("moment") << tab << moment_
        << token::END_STATEMENT << nl;
    os.writeKeyword("forceRelative") << tab << forceRelative_
        << token::END_STATEMENT << nl;
    os.writeKeyword("momentRelative") << tab << momentRelative_
        << token::END_STATEMENT << endl;

    return os.good();
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const sixDOFODE& sds)
{
    sds.writeData(os);

    os.check("Ostream& operator<<(Ostream&, const sixDOFODE");

    return os;
}


// ************************************************************************* //

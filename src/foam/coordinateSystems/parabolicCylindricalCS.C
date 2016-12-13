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

\*---------------------------------------------------------------------------*/

#include "parabolicCylindricalCS.H"
#include "mathematicalConstants.H"
#include "boundBox.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(parabolicCylindricalCS, 0);
    addToRunTimeSelectionTable
    (
        coordinateSystem,
        parabolicCylindricalCS,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parabolicCylindricalCS::parabolicCylindricalCS(const bool inDegrees)
:
    coordinateSystem(),
    inDegrees_(inDegrees)
{}


Foam::parabolicCylindricalCS::parabolicCylindricalCS
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr,
    const bool inDegrees
)
:
    coordinateSystem(name, origin, cr),
    inDegrees_(inDegrees)
{}


Foam::parabolicCylindricalCS::parabolicCylindricalCS
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict),
    inDegrees_(dict.lookupOrDefault<Switch>("degrees", true))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::coordinateSystem::spanInfo
Foam::parabolicCylindricalCS::spanLimited() const
{
    spanInfo b(Pair<bool>(true, true));

    // Upper bound or r is unlimited
    b[0] = Pair<bool>(true, false);

    // z is unlimited
    b[2] = Pair<bool>(false, false);

    return b;
}


Foam::boundBox Foam::parabolicCylindricalCS::spanBounds() const
{
    return boundBox
    (
        vector
        (
            0,
            ( inDegrees_ ? -180.0 : -mathematicalConstant::pi ),
            -GREAT
        ),
        vector
        (
            GREAT,
            ( inDegrees_ ? 180.0 : mathematicalConstant::pi ),
            GREAT
        )
    );
}


bool Foam::parabolicCylindricalCS::inDegrees() const
{
    return inDegrees_;
}


Foam::Switch& Foam::parabolicCylindricalCS::inDegrees()
{
    return inDegrees_;
}


Foam::vector Foam::parabolicCylindricalCS::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    // Notation: u = local.x() v = local.y() z = local.z();
    if (local.y() < 0.0)
    {
        FatalErrorIn
        (
            "parabolicCylindricalCS::localToGlobal(const vector&, bool) const"
        )
            << "parabolic cylindrical coordinates v < 0"
            << abort(FatalError);
    }

    return coordinateSystem::localToGlobal
    (
        vector
        (
            0.5*(sqr(local.x()) - sqr(local.y())),
            local.x()*local.y(),
            local.z()
        ),
        translate
    );
}


Foam::tmp<Foam::vectorField> Foam::parabolicCylindricalCS::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    if (min(local.component(vector::Y)) < 0.0)
    {
        FatalErrorIn
        (
            "parabolicCylindricalCS::localToGlobal"
            "(const vectorField&, bool) const"
        )   << "parabolic cylindrical coordinates v < 0"
            << abort(FatalError);
    }

    vectorField lc(local.size());
    lc.replace
    (
        vector::X,
        0.5*
        (
            sqr(local.component(vector::X))
          - sqr(local.component(vector::Y))
        )
    );

    lc.replace
    (
        vector::Y,
        local.component(vector::X) * local.component(vector::Y)
    );

    lc.replace
    (
        vector::Z,
        local.component(vector::Z)
    );

    return coordinateSystem::localToGlobal(lc, translate);
}


Foam::vector Foam::parabolicCylindricalCS::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    notImplemented
    (
        "parabolicCylindricalCS::globalToLocal(const vector&, bool) const"
    );

    return vector::zero;
}

Foam::tmp<Foam::vectorField> Foam::parabolicCylindricalCS::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    notImplemented
    (
        "parabolicCylindricalCS::globalToLocal(const vectorField&, bool) const"
    );

    return tmp<vectorField>(vectorField::null());
}


void Foam::parabolicCylindricalCS::write(Ostream& os) const
{
    coordinateSystem::write(os);
    os << " inDegrees: " << inDegrees() << endl;
}


void Foam::parabolicCylindricalCS::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    coordinateSystem::writeDict(os, false);
    os.writeKeyword("inDegrees") << inDegrees_ << token::END_STATEMENT << nl;

    if (subDict)
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}


// ************************************************************************* //

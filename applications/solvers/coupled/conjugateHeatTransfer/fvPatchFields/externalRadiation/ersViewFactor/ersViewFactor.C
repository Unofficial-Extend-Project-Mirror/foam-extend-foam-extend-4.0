/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "ersViewFactor.H"
#include "addToRunTimeSelectionTable.H"
#include "radiationConstants.H"

namespace Foam
{

defineTypeNameAndDebug(ersViewFactor, 0);
addToRunTimeSelectionTable(externalRadiationSource, ersViewFactor, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ersViewFactor::ersViewFactor
(
    const word& name,
    const dictionary& dict,
    const fvPatch& p
)
:
    externalRadiationSource(name),
    Tinf_(readScalar(dict.lookup("Tinf"))),
    F_("F", dict, p.size()),
    epsilon_(readScalar(dict.lookup("epsilon")))
{}

ersViewFactor::ersViewFactor
(
    const word& name,
    const dictionary& dict
)
:
    externalRadiationSource(name),
    Tinf_(readScalar(dict.lookup("Tinf"))),
    epsilon_(readScalar(dict.lookup("epsilon")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> ersViewFactor::q
(
    const scalarField& Tw
) const
{
    return epsilon_*F_*radiation::sigmaSB.value()*(pow4(Tw) - pow4(Tinf_));
}


void ersViewFactor::write(Ostream& os) const
{
    externalRadiationSource::write(os);

    os.writeKeyword("Tinf") << Tinf_ << token::END_STATEMENT << nl;
    F_.writeEntry("F", os);
    os.writeKeyword("epsilon") << epsilon_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

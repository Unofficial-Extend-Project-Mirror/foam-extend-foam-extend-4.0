/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "viewFactor.H"
#include "addToRunTimeSelectionTable.H"
#include "radiationConstants.H"

namespace Foam
{

defineTypeNameAndDebug(viewFactor, 0);
addToRunTimeSelectionTable(externalRadiationSource, viewFactor, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

viewFactor::viewFactor
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

viewFactor::viewFactor
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

tmp<scalarField> viewFactor::q
(
    const scalarField& Tw
) const
{
    return epsilon_*F_*radiation::sigmaSB.value()*(pow4(Tw) - pow4(Tinf_));
}


void viewFactor::write(Ostream& os) const
{
    externalRadiationSource::write(os);

    os.writeKeyword("Tinf") << Tinf_ << token::END_STATEMENT << nl;
    F_.writeEntry("F", os);
    os.writeKeyword("epsilon") << epsilon_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

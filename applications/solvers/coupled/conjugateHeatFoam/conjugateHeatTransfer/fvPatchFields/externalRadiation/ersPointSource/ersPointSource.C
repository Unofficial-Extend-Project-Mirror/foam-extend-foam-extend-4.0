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

#include "ersPointSource.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(ersPointSource, 0);
addToRunTimeSelectionTable(externalRadiationSource, ersPointSource, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ersPointSource::ersPointSource
(
    const word& name,
    const dictionary& dict,
    const fvPatch& p
)
:
    ersConstantFlux(name),
    qmax_(readScalar(dict.lookup("qmax"))),
    alpha_(readScalar(dict.lookup("alpha"))),
    direction_(dict.lookup("direction"))
{
    q() = -alpha_*qmax_*min(direction_ & p.Sf()/p.magSf(), 0.0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ersPointSource::write(Ostream& os) const
{
    ersConstantFlux::write(os);
    os.writeKeyword("qmax") << qmax_ << token::END_STATEMENT << nl;
    os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;
    os.writeKeyword("direction") << direction_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

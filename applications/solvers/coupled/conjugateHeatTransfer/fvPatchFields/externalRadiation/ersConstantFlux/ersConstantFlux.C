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

#include "ersConstantFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(ersConstantFlux, 0);
addToRunTimeSelectionTable(externalRadiationSource, ersConstantFlux, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ersConstantFlux::ersConstantFlux
(
    const word& name,
    const dictionary& dict,
    const fvPatch& p
)
:
    externalRadiationSource(name),
    q_("q", dict, p.size())
{}

ersConstantFlux::ersConstantFlux
(
    const word& name
)
:
    externalRadiationSource(name)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ersConstantFlux::write(Ostream& os) const
{
    externalRadiationSource::write(os);

    q_.writeEntry("q", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
